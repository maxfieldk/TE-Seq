module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
library(tidyverse)
library(magrittr)
library(flextable)
library(magrittr)
library(modelr)
library(tidybayes)
library(tidybayes.rethinking)
library(rstan)
library(rethinking)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

outputdir <- "bayes"
dir.create(outputdir)
# Load data and prepare for model
df <- read_csv("meth_for_bayes.csv")

# Prepare data with integer indices for hierarchical model
dat <- df %>%
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos)),
        sample = as.integer(factor(sample)),
        condition = factor(condition),
        gene_id = factor(gene_id)
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id)

dat_index <- df %>%
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = as.integer(factor(as.character(consensus_pos))),
        sample = as.integer(factor(sample)),
        condition = as.integer(factor(condition)),
        gene_id = as.integer(factor(gene_id))
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id)
dat %$% UTR_status
# Define and fit the model
m_var1 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a[consensus_pos] + b[consensus_pos] * condition - 1,
        a[consensus_pos] ~ dnorm(a_bar, a_sigma),
        b[consensus_pos] ~ dnorm(b_bar, b_sigma),
        a_bar ~ dnorm(0, 0.5),
        a_sigma ~ dexp(0.5),
        b_bar ~ dnorm(0, 0.5),
        b_sigma ~ dexp(0.5)
    ),
    data = dat, chains = 4, cores = 4, log_lik = FALSE
)

m_var2 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a_c[consensus_pos] + a_g[gene_id] + b * condition - 1,
        a_c[consensus_pos] ~ dnorm(a_c_bar, a_c_sigma),
        a_g[consensus_pos] ~ dnorm(a_g_bar, a_g_sigma),
        b ~ dnorm(b_bar, b_sigma),
        a_c_bar ~ dnorm(0, 0.5),
        a_c_sigma ~ dexp(0.5),
        a_g_bar ~ dnorm(0, 0.5),
        a_g_sigma ~ dexp(0.5),
        b ~ dnorm(0, 0.5),
    ),
    data = dat, chains = 4, cores = 4, log_lik = FALSE
)

m_var2_uncentered <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a_c_bar + z_c[consensus_pos] * a_c_sigma +
            a_g_bar + z_g[gene_id] * a_g_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),

        # Hyperpriors
        a_c_bar ~ dnorm(0, 1),
        a_c_sigma ~ dexp(1),
        a_g_bar ~ dnorm(0, 1),
        a_g_sigma ~ dexp(1),

        # Prior for b
        b ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- a_c_bar + z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- a_g_bar + z_g * a_g_sigma
    ),
    data = dat_index, chains = 4, cores = 4, log_lik = FALSE
)

saveRDS(m_var2_uncentered, "m_var2_uncentered.rds", compress = FALSE)

m_var2_uncentered <- readRDS("m_var2_uncentered.rds")

m_var2_stratified <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),

        # The logistic regression for p, now with a varying b by gene_id and consensus_pos
        logit(p) <- a_c_bar + z_c[consensus_pos] * a_c_sigma +
            a_g_bar + z_g[gene_id] * a_g_sigma +
            b_bar + z_b[gene_id, consensus_pos] * b_sigma,

        # Non-centered parameters for consensus_pos and gene_id
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),

        # Non-centered varying effect for b by gene_id and consensus_pos
        matrix[gene_id, consensus_pos]:z_b ~ dnorm(0, 1),

        # Hyperpriors for varying intercepts
        a_c_bar ~ dnorm(0, 1),
        a_c_sigma ~ dexp(1),
        a_g_bar ~ dnorm(0, 1),
        a_g_sigma ~ dexp(1),

        # Hyperpriors for b stratified by gene_id and consensus_pos
        b_bar ~ dnorm(0, 1),
        b_sigma ~ dexp(1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- a_c_bar + z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- a_g_bar + z_g * a_g_sigma,
        gq > matrix[gene_id, consensus_pos]:b_stratified <<- b_bar + z_b * b_sigma
    ),
    data = dat_index, chains = 4, cores = 4, log_lik = FALSE
)

saveRDS(m_var2_stratified, "m_var2_stratified.rds", compress = FALSE)

# m2.3 <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a_position[consensus_pos_index] + a_gene_id[gene_id_index] + a_condition[condition_index],
#         a_position[consensus_pos_index] ~ dnorm(a_position_bar, sigma_position),
#         a_position_bar ~ dnorm(0, 2),
#         sigma_position ~ dexp(0.5),
#         a_gene_id[gene_id_index] ~ dnorm(a_gene_id_bar, sigma_gene_id),
#         a_gene_id_bar ~ dnorm(0, 2),
#         sigma_gene_id ~ dexp(0.5),
#         a_condition[condition_index] ~ dnorm(a_condition_bar, sigma_condition),
#         a_condition_bar ~ dnorm(0, 2),
#         sigma_condition ~ dexp(0.5)
#     ),
#     data = df_m %>% dplyr::select(-UTR_status, -sample_index), chains = 4, cores = 4, log_lik = FALSE
# )

summary(m_var)
m_var1 <- m_var1 %>%
    recover_types(dat)
model_name <- "m_var1"
p <- m_var1 %>%
    spread_draws(a[consensus_pos], b[consensus_pos], a_bar, b_bar) %>%
    mutate(consensus_pos = as.integer(consensus_pos)) %>%
    ggplot(aes(y = a, x = consensus_pos)) +
    stat_pointinterval() +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/%s/a_summary.pdf", outputdir, model_name), 10, 4)

p <- m_var1 %>%
    spread_draws(a[consensus_pos], b[consensus_pos], a_bar, b_bar) %>%
    mutate(consensus_pos = as.integer(consensus_pos)) %>%
    ggplot(aes(y = b, x = consensus_pos)) +
    stat_pointinterval() +
    mtclosed
mysaveandstore(sprintf("%s/%s/b_summary.pdf", outputdir, model_name), 10, 4)

p <- m_var1 %>%
    spread_draws(a[consensus_pos], b[consensus_pos], a_bar, b_bar) %>%
    ggplot(aes(y = b_bar, x = "b_bar")) +
    stat_halfeye() +
    mtclosed
mysaveandstore(sprintf("%s/%s/b_summary.pdf", outputdir, model_name), 4, 4)

p <- m_var1 %>%
    spread_draws(a[consensus_pos], b[consensus_pos], a_bar, b_bar) %>%
    ggplot(aes(y = a_bar, x = "a_bar")) +
    stat_halfeye() +
    mtclosed
mysaveandstore(sprintf("%s/%s/a_summary.pdf", outputdir, model_name), 4, 4)

# lin pred
lin_pred_draws <- dat %>%
    mutate(condition = as.integer(condition) - 1) %>%
    data_grid(consensus_pos, condition) %>%
    add_linpred_draws(m_var1)
p <- lin_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(y = .linpred, x = consensus_pos)) +
    stat_pointinterval(.width = c(.66, .95)) +
    labs(y = "P(CpG is Meth)") +
    mtclosed
mysaveandstore(sprintf("%s/%s/linpred_a.pdf", outputdir, model_name), 10, 4)


p <- lin_pred_draws %>%
    mutate(condition = as.factor(condition)) %>% # Ensure condition is a factor
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(y = .linpred, x = consensus_pos, group = condition, color = condition)) +
    stat_pointinterval(.width = c(.66, .95), position = position_dodge2(width = 0.5)) + # Adjust position dodging here
    scale_color_brewer(palette = "Set2") +
    labs(y = "P(CpG is Meth)") +
    mtclosed
mysaveandstore(sprintf("%s/%s/linpred_by_condition.pdf", outputdir, model_name), 10, 4)

# posterior predictions
post_pred_draws <- dat %>%
    mutate(condition = as.integer(condition) - 1) %>%
    data_grid(consensus_pos, condition) %>%
    mutate(cov = 30) %>%
    add_predicted_draws(m_var1) %>%
    mutate(.prediction = .prediction / 30)


p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    scale_color_brewer() +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred.pdf", outputdir, model_name), 10, 4)

p <- post_pred_draws %>%
    group_by(consensus_pos, condition) %>%
    slice_sample(n = 2) %>%
    ungroup() %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    geom_point() +
    labs(y = "Predicted Pct Meth") +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred.pdf", outputdir, model_name), 10, 4)


p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    geom_point(aes(y = pct_meth, alpha = 1 / 4),
        data = dat %>%
            mutate(pct_meth = Nmeth / cov) %>%
            mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
            group_by(consensus_pos, condition) %>%
            slice_sample(n = 5)
    ) +
    scale_color_brewer() +
    labs(y = "Predicted Pct Meth") +
    mtclosed

mysaveandstore(sprintf("%s/%s/postpred_with10rawdata.pdf", outputdir, model_name), 20, 4)




##### EXAMPLES
set.seed(5)
n <- 10
n_condition <- 5
ABC <-
    tibble(
        condition = factor(rep(c("A", "B", "C", "D", "E"), n)),
        response = rnorm(n * 5, c(0, 1, 2, 1, -1), 0.5)
    )

ABC %>%
    ggplot(aes(y = fct_rev(condition), x = response)) +
    geom_point()

m <- ulam(
    alist(
        response ~ normal(mu, sigma),

        # submodel for conditional mean
        mu <- intercept[condition],
        intercept[condition] ~ normal(mu_condition, tau_condition),
        mu_condition ~ normal(0, 5),
        tau_condition ~ exponential(1),

        # submodel for conditional standard deviation
        log(sigma) <- sigma_intercept[condition],
        sigma_intercept[condition] ~ normal(0, 1)
    ),
    data = ABC,
    chains = 1,
    cores = parallel::detectCores(),
    iter = 2000
)

summary(m)

m %>%
    spread_draws(intercept[condition]) %>%
    head(10)

m %>%
    recover_types(ABC) %>%
    spread_draws(intercept[condition]) %>%
    head(10)

m %<>% recover_types(ABC)


m %>%
    spread_draws(mu_condition, tau_condition) %>%
    head(10)

m %>%
    spread_draws(intercept[condition]) %>%
    median_qi() %>%
    ggplot(aes(y = fct_rev(condition), x = intercept, xmin = .lower, xmax = .upper)) +
    geom_pointinterval()

m %>%
    spread_draws(intercept[condition]) %>%
    ggplot(aes(y = fct_rev(condition), x = intercept)) +
    stat_pointinterval()
plot(1, 1)

m %>%
    spread_draws(mu_condition, intercept[condition]) %>%
    head(10)

m %>%
    spread_draws(mu_condition, intercept[condition]) %>%
    mutate(condition_offset = intercept - mu_condition) %>%
    median_qi(condition_offset)


ABC %>%
    data_grid(condition) %>%
    add_linpred_draws(m) %>%
    head(10)

ABC %>%
    data_grid(condition) %>%
    add_linpred_draws(m) %>%
    ggplot(aes(x = .linpred, y = fct_rev(condition))) +
    stat_pointinterval(.width = c(.66, .95))


ABC %>%
    data_grid(condition) %>%
    add_predicted_draws(m) %>%
    ggplot(aes(x = .prediction, y = condition)) +
    stat_slab()

grid <- ABC %>%
    data_grid(condition)

fits <- grid %>%
    add_linpred_draws(m)

preds <- grid %>%
    add_predicted_draws(m)

ABC %>%
    ggplot(aes(y = condition, x = response)) +
    stat_interval(aes(x = .prediction), data = preds) +
    stat_pointinterval(aes(x = .linpred), data = fits, .width = c(.66, .95), position = position_nudge(y = -0.3)) +
    geom_point() +
    scale_color_brewer()


####
mtcars_clean <- mtcars %>%
    mutate(cyl = factor(cyl))

m_mpg <- ulam(
    alist(
        mpg ~ normal(mu, sigma),
        mu <- intercept[cyl] + slope[cyl] * hp,
        intercept[cyl] ~ normal(20, 10),
        slope[cyl] ~ normal(0, 10),
        sigma ~ exponential(1)
    ),
    data = mtcars_clean,
    chains = 1,
    cores = parallel::detectCores(),
    iter = 2000
)

mtcars_clean %>%
    group_by(cyl) %>%
    data_grid(hp = seq_range(hp, n = 51)) %>%
    add_linpred_draws(m_mpg) %>%
    ggplot(aes(x = hp, y = mpg, color = cyl)) +
    stat_lineribbon(aes(y = .linpred)) +
    geom_point(data = mtcars_clean) +
    scale_fill_brewer(palette = "Greys") +
    scale_color_brewer(palette = "Set2")


mtcars_clean %>%
    group_by(cyl) %>%
    data_grid(hp = seq_range(hp, n = 101)) %>%
    add_predicted_draws(m_mpg) %>%
    ggplot(aes(x = hp, y = mpg, color = cyl, fill = cyl)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1 / 4) +
    geom_point(data = mtcars_clean) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2")
