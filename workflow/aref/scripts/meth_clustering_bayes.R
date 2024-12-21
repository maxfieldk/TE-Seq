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
library(ggridges)

set.seed(123)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

outputdir <- "bayes"
dir.create(outputdir)
stan_workdir <- sprintf("%s/stan_workdir", outputdir)
dir.create(stan_workdir, recursive = TRUE)

# Load data and prepare for model
source("conf/sample_table_source.R")
sample_table$age_scaled <- scale(sample_table$age)
df <- read_csv("meth_for_bayes.csv")

df <- df %>% left_join(sample_table %>% dplyr::select(sample_name, sex, age_scaled, braak, apoe) %>% dplyr::rename(sample = sample_name))

dat_alpha <- df %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos), levels = as.character(df$consensus_pos %>% unique() %>% sort())),
        sample = factor(sample),
        condition = factor(condition, levels = c("CTRL", "AD")),
        gene_id = factor(gene_id),
        sex = factor(sex, levels = c("F", "M")),
        braak = factor(braak, levels = c())
    )

dat_alpha$age
restore_sample <- tibble(sample_val = dat_alpha$sample %>% unique()) %>% mutate(sample = sample_val %>% as.integer())
restore_condition <- tibble(condition_val = dat_alpha$condition %>% unique()) %>% mutate(condition = condition_val %>% unique() %>% as.integer() - 1)
restore_consensus_pos <- tibble(consensus_pos_val = dat_alpha$consensus_pos %>% unique()) %>% mutate(consensus_pos = consensus_pos_val %>% as.integer())
restore_gene_id <- tibble(gene_id_val = dat_alpha$gene_id %>% unique()) %>% mutate(gene_id = gene_id_val %>% as.integer())

dat_index <-  df %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos), levels = as.character(df$consensus_pos %>% unique() %>% sort())),
        sample = factor(sample),
        condition = factor(condition, levels = c("CTRL", "AD")),
        gene_id = factor(gene_id)
    ) %>%
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        consensus_pos = as.integer(consensus_pos),
        sample = as.integer(sample),
        condition = as.integer(as.integer(condition) - 1),
        gene_id = as.integer(gene_id)
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id)


dat_index_328 <-  df %>%
    filter(consensus_pos <= 328) %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos), levels = as.character(df$consensus_pos %>% unique() %>% sort())),
        sample = factor(sample),
        condition = factor(condition, levels = c("CTRL", "AD")),
        gene_id = factor(gene_id),
        sex = factor(sex, levels = c("F", "M"))) %>%
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        consensus_pos = as.integer(consensus_pos),
        sample = as.integer(sample),
        condition = as.integer(as.integer(condition) - 1),
        gene_id = as.integer(gene_id),
        sex = as.integer(as.integer(sex) - 1),
        braak = as.integer(braak)
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id, sex, braak, age_scaled)

dat_index_500 <- df %>%
    filter(consensus_pos <= 500) %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos), levels = as.character(df$consensus_pos %>% unique() %>% sort())),
        sample = factor(sample),
        condition = factor(condition, levels = c("CTRL", "AD")),
        gene_id = factor(gene_id)
    ) %>%
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        consensus_pos = as.integer(consensus_pos),
        sample = as.integer(sample),
        condition = as.integer(as.integer(condition) - 1),
        gene_id = as.integer(gene_id)
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id)
dat_index_909 <- df %>%
    filter(consensus_pos <= 909) %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(pctM / 100 * cov),
        UTR_status = as.integer(factor(ifelse(consensus_pos < 909, "UTR", "Body"))),
        consensus_pos = factor(as.character(consensus_pos), levels = as.character(df$consensus_pos %>% unique() %>% sort())),
        sample = factor(sample),
        condition = factor(condition, levels = c("CTRL", "AD")),
        gene_id = factor(gene_id)
    ) %>% 
    mutate(
        Nmeth = as.integer(pctM / 100 * cov),
        consensus_pos = as.integer(consensus_pos),
        sample = as.integer(sample),
        condition = as.integer(as.integer(condition) - 1),
        gene_id = as.integer(gene_id)
    ) %>%
    dplyr::select(Nmeth, cov, UTR_status, consensus_pos, sample, condition, gene_id)


# models
model_name <- "mv2_328"
mv2_328 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_328, chains = 4, cores = 4, threads = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model_name <- "mv2_328_prior"
mv2_328__prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_328, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)

model_name <- "mv2_500"
mv2_500 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_500, chains = 4, cores = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir, threads = 4
)
model_name <- "mv2_500_prior"
mv2_500_prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_500, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)

model_name <- "mv2_909"
mv2_909 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_909, chains = 4, cores = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir, threads = 4
)
model_name <- "mv2_909_prior"
mv2_909_prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_909, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)

model_name <- "mv2_full"
mv2 <- ulam(
  alist(
    Nmeth ~ dbinom(cov, p),
    logit(p) <- a +
      z_c[consensus_pos] * a_c_sigma +
      z_g[gene_id] * a_g_sigma +
      z_s[sample] * a_s_sigma +
      b[UTR_status] * condition,  # Slope for condition varies by UTR_status

    # Non-centered parameters
    vector[consensus_pos]:z_c ~ dnorm(0, 1),
    vector[gene_id]:z_g ~ dnorm(0, 1),
    vector[sample]:z_s ~ dnorm(0, 1),
    vector[UTR_status]:b ~ dnorm(0, b_sigma),  # Prior for varying slopes of condition

    # Hyperpriors
    a_c_sigma ~ dexp(1),
    a_g_sigma ~ dexp(1),
    a_s_sigma ~ dexp(1),
    b_sigma ~ dexp(1),  # Hyperprior for variability of b across UTR_status
    a ~ dnorm(0, 1)
  ),
  data = dat_index, chains = 4, cores = 4, log_lik = FALSE
)
model_name <- "mv2_full_prior"
mv2_full_prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b[UTR_status] * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)


####### with covariates

model_name <- "mv2_328_with_sex_age"
mv2_328_with_sex_age <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition +
            beta_sex * sex +           # Effect of sex
            beta_age * age_scaled,            # Effect of age

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),

        # Priors for coefficients
        b ~ dnorm(0, 0.5),             # Prior for condition slope
        beta_sex ~ dnorm(0, 0.5),      # Prior for sex effect
        beta_age ~ dnorm(0, 0.5),      # Prior for age effect
        a ~ dnorm(0, 1),               # Prior for intercept

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_328, chains = 4, cores = 4, threads = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)

model_name <- "mv2_328_with_sex_age_prior"
mv2_328_with_sex_age_prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition +
            beta_sex * sex +           # Effect of sex
            beta_age * age_scaled,            # Effect of age

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),

        # Priors for coefficients
        b ~ dnorm(0, 0.5),             # Prior for condition slope
        beta_sex ~ dnorm(0, 0.5),      # Prior for sex effect
        beta_age ~ dnorm(0, 0.5),      # Prior for age effect
        a ~ dnorm(0, 1),               # Prior for intercept

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dat_index_328, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)


dat <- list(
    Nmeth = dat_index_328$Nmeth,
    cov = dat_index_328$cov,
    consensus_pos = dat_index_328$consensus_pos,
    gene_id = dat_index_328$gene_id,
    sample = dat_index_328$sample,
    braak = dat_index_328$braak, # BRAAK as an ordered index
    alpha = rep(2, max(dat_index_328$braak) - 1) # Dirichlet prior
)

model_name <- "mv2_braak_binomial"
mv2_braak_binomial <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            bB * sum(delta_j[1:braak]),
        
        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),
        
        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        
        # Coefficient priors
        bB ~ normal(0, 1),
        a ~ normal(0, 1),
        
        # Delta parameters for BRAAK stage
        vector[8]:delta_j <<- append_row(0, delta), # delta_j includes a fixed 0
        simplex[7]:delta ~ dirichlet(alpha)     # Dirichlet prior on delta
    ),
    data = dat,
    chains = 4, cores = 4, threads = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)



# model_name <- "mv2"
# mv2 <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a +
#             z_c[consensus_pos] * a_c_sigma +
#             z_g[gene_id] * a_g_sigma +
#             z_s[sample] * a_s_sigma +
#             b * condition,

#         # Non-centered parameters
#         vector[consensus_pos]:z_c ~ dnorm(0, 1),
#         vector[gene_id]:z_g ~ dnorm(0, 1),
#         vector[sample]:z_s ~ dnorm(0, 1),

#         # Hyperpriors
#         a_c_sigma ~ dexp(1),
#         a_g_sigma ~ dexp(1),
#         a_s_sigma ~ dexp(1),
#         # Prior for b
#         b ~ dnorm(0, 0.5),
#         a ~ dnorm(0, 1),

#         # Generated quantities for easier interpretation
#         gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
#         gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
#         gq > vector[sample]:a_s <<- z_s * a_s_sigma
#     ),
#     data = dat_index, chains = 4, cores = 4, log_lik = FALSE,
#     file = sprintf("%s/%s", outputdir, model_name),
#     output_dir = stan_workdir
# )


# mf1 <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a + a_c[consensus_pos] +
#             a_g[gene_id] +
#             a_s[sample] +
#             b * condition,
#         b ~ dnorm(0, 1),
#         a ~ dnorm(0, 1),
#         a_c[consensus_pos] ~ dnorm(0, 1),
#         a_g[gene_id] ~ dnorm(0, 1),
#         a_s[sample] ~ dnorm(0, 1)
#     ),
#     data = datindexs, chains = 4, cores = 4, log_lik = FALSE
# )


# mv1 <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a +
#             a_s_bar + z_s[sample] * a_s_sigma +
#             b * condition,

#         # Non-centered parameters
#         vector[sample]:z_s ~ dnorm(0, 1),

#         # Hyperpriors
#         a_s_bar ~ dnorm(0, 1),
#         a_s_sigma ~ dexp(1),
#         # Prior for b
#         b ~ dnorm(0, 1),
#         a ~ dnorm(0, 1),

#         # Generated quantities for easier interpretation
#         gq > vector[sample]:a_s <<- a_s_bar + z_s * a_s_sigma
#     ),
#     data = datindexs, chains = 4, cores = 4, log_lik = FALSE
# )

# mv2 <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a + a_c_bar + z_c[consensus_pos] * a_c_sigma +
#             a_g_bar + z_g[gene_id] * a_g_sigma +
#             a_s_bar + z_s[sample] * a_s_sigma +
#             b * condition,

#         # Non-centered parameters
#         vector[consensus_pos]:z_c ~ dnorm(0, 1),
#         vector[gene_id]:z_g ~ dnorm(0, 1),
#         vector[sample]:z_s ~ dnorm(0, 1),

#         # Hyperpriors
#         a_c_bar ~ dnorm(0, 1),
#         a_c_sigma ~ dexp(1),
#         a_g_bar ~ dnorm(0, 1),
#         a_g_sigma ~ dexp(1),
#         a_s_bar ~ dnorm(0, 1),
#         a_s_sigma ~ dexp(1),
#         # Prior for b
#         b ~ dnorm(0, 1),
#         a ~ dnorm(0, 1),

#         # Generated quantities for easier interpretation
#         gq > vector[consensus_pos]:a_c <<- a_c_bar + z_c * a_c_sigma,
#         gq > vector[gene_id]:a_g <<- a_g_bar + z_g * a_g_sigma,
#         gq > vector[sample]:a_s <<- a_s_bar + z_s * a_s_sigma
#     ),
#     data = datindexs, chains = 4, cores = 4, log_lik = FALSE
# )

################




model <- mv2_braak_binomial
modelprior <- mv2_328_with_sex_age_prior
model_name <- "mv2_braak_binomial"
model_data <- dat
precis(model %>%
    spread_draws(a, bB))
precis(model_no_confounds %>%
    spread_draws(a, b, a_c[..], a_g[..], a_s[..]))
precis(model %>%
    spread_draws(b, beta_age))
precis(model_no_confounds %>%
    spread_draws(b))


draws_spread <- model %>%
    spread_draws(a, b, a_c_sigma, a_c[..], a_g_sigma, a_g[..], a_s_sigma, a_s[..])

draws_gathered <- model %>%
    gather_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    left_join(restore_consensus_pos) %>%
    left_join(restore_sample) %>%
    left_join(restore_gene_id) %>% 
    ungroup()


draws_spread_prior <- model %>%
    spread_draws(a, b, a_c_sigma, a_c[..], a_g_sigma, a_g[..], a_s_sigma, a_s[..])

draws_gathered_prior <- modelprior %>%
    gather_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    left_join(restore_consensus_pos) %>%
    left_join(restore_sample) %>%
    left_join(restore_gene_id) %>%
    ungroup()

p <- draws_gathered %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/%s/summary.pdf", outputdir, model_name), 4, 4)

# these are simulation plots showing real effect and estimate
p <- draws_gathered %>%
    filter(.variable == "a_c") %>%
    filter(!is.na(consensus_pos_val)) %>%
    ggplot(aes(y = consensus_pos_val, x = .value)) +
    stat_halfeye() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(
        data =
            draws_gathered_prior %>%
                filter(.variable == "a_c") %>%
                filter(!is.na(consensus_pos_val)),
        alpha = 0.2, color = "blue"
    ) +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_ac.pdf", outputdir, model_name), 4, 4)


p <- draws_gathered %>%
    filter(.variable == "a_s") %>%
    filter(!is.na(sample)) %>%
    ggplot(aes(y = sample, x = .value)) +
    stat_halfeye() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(
        data = draws_gathered_prior %>%
            filter(.variable == "a_s") %>%
            filter(!is.na(sample)),
     alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_as.pdf", outputdir, model_name), 4, 4)

p <- draws_gathered %>%
    filter(.variable == "a_g") %>%
    filter(!is.na(gene_id)) %>%
    ggplot(aes(y = gene_id, x = .value)) +
    stat_halfeye() +
    stat_halfeye(
        data =
            draws_gathered_prior %>%
                filter(.variable == "a_g") %>%
                filter(!is.na(gene_id)), alpha = 0.2, color = "blue"
    ) +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_ag.pdf", outputdir, model_name), 4, 40)

p <- draws_gathered %>%
    filter(.variable == "b") %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye() +
    stat_halfeye(
        data =
            draws_gathered_prior %>%
                filter(.variable == "b"),
        alpha = 0.2, color = "blue"
    ) +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_b.pdf", outputdir, model_name), 4, 4)

b_samples <- draws_spread %>% pull(b)
left_tail <- mean(b_samples < 0)
right_tail <- mean(b_samples > 0)
one_tailed_pval <- min(left_tail, right_tail)
two_tailed_pval <- 2 * min(left_tail, right_tail)
statsframe <- tibble(parameter = c("b"), one_tailed_pval = one_tailed_pval, two_tailed_pval = two_tailed_pval)
write_delim(statsframe, sprintf("%s/%s/stats_b.tsv", outputdir, model_name))

# lin pred
# note that I am only using first 5 gene_ids and first 10 consensus_pos else the draws frame is 200 M rows long...
lin_pred_draws <- model_data %>%
    data_grid(condition, sample, gene_id = 1:5, consensus_pos = 1:10) %>%
    add_linpred_draws(model)
p <- lin_pred_draws %>% ggplot(aes(x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred.pdf", outputdir, model_name), 10, 4)

p <- lin_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred_bycondition.pdf", outputdir, model_name), 10, 4)

# posterior predictions
post_pred_draws <- model_data %>%
    data_grid(condition, sample, gene_id = 1:5, consensus_pos = 1:10) %>%
    mutate(cov = 30) %>%
    add_predicted_draws(model) %>%
    mutate(.prediction = .prediction / 30)
p <- post_pred_draws %>% ggplot(aes(x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred.pdf", outputdir, model_name), 10, 4)

p <- post_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred_condition.pdf", outputdir, model_name), 10, 4)


p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    scale_color_brewer() +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred_by_cpos.pdf", outputdir, model_name), 10, 4)
}



################

summary(mf1)
summary(mv1)
summary(mv2)

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


########### SIMULATION
mean_cov <- 30
Nconsensus_pos <- 5
Ngene_id <- 10
Nsample <- 12

sf <- crossing(consensus_pos = 1:Nconsensus_pos, sample = 1:Nsample, gene_id = 1:Ngene_id) %>%
    mutate(condition = ifelse(sample <= Nsample / 2, 0, 1))

effect_consensus_pos <- tibble(consensus_pos = 1:Nconsensus_pos, consensus_pos_effect = rnorm(n = Nconsensus_pos, 0, sd = 0.2)) %>%
    mutate(consensus_pos_effect = consensus_pos_effect - mean(consensus_pos_effect))
effect_sample <- tibble(sample = 1:Nsample, sample_effect = rnorm(n = Nsample, 0, sd = 0.2)) %>%
    mutate(sample_effect = sample_effect - mean(sample_effect))

effect_gene_id <- tibble(gene_id = 1:Ngene_id, gene_id_effect = rnorm(n = Ngene_id, 0, sd = 0.2)) %>%
    mutate(gene_id_effect = gene_id_effect - mean(gene_id_effect))

effect_condition <- tibble(condition = c(0, 1), condition_effect = c(0.15, -0.15))

dsim_space <- sf %>%
    left_join(effect_consensus_pos) %>%
    left_join(effect_sample) %>%
    left_join(effect_gene_id) %>%
    left_join(effect_condition) %>%
    mutate(condition = as.integer(condition)) %>%
    mutate(cov = rpois(length(.$sample), mean_cov)) %>%
    mutate(alpha = 1) %>%
    mutate(p = logistic(alpha + sample_effect + consensus_pos_effect + gene_id_effect + condition_effect)) %>%
    rowwise() %>%
    ungroup()

reapply_bind_rows <- function(data, times) {
    result <- data
    for (i in seq_len(times)) {
        result <- bind_rows(result, data)
    }
    result
}



sampled_dsim_space <- dsim_space %>% sample_n(500, replace = TRUE)
# Apply `bind_rows` to `dsim_space` 3 times
sampled_dsim_space <- reapply_bind_rows(dsim_space, 1)
dsim <- sampled_dsim_space %>%
    rowwise() %>%
    mutate(Nmeth = sum(rbinom(cov, 1, p))) %>%
    ungroup()
dsim_for_model <- dsim %>% dplyr::select(Nmeth, cov, gene_id, sample, consensus_pos, condition)
dsim_for_model %>%
    mutate(pctM = Nmeth / 30) %>%
    group_by(condition) %>%
    summarise(meanmeth = mean(pctM))


mv2sim_prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a + a_c_bar + z_c[consensus_pos] * a_c_sigma +
            a_g_bar + z_g[gene_id] * a_g_sigma +
            a_s_bar + z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_bar ~ dnorm(0, 0.5),
        a_c_sigma ~ dexp(1),
        a_g_bar ~ dnorm(0, 0.5),
        a_g_sigma ~ dexp(1),
        a_s_bar ~ dnorm(0, 0.5),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- a_c_bar + z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- a_g_bar + z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- a_s_bar + z_s * a_s_sigma
    ),
    data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE
)

model <- mv2sim1
model_name <- "mv2sim_prior"
precis(model %>%
    spread_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]))

p <- model %>%
    gather_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_density_ridges() +
    lims(x = c(-2, 2)) +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/%s/summary.pdf", outputdir, model_name), 10, 4)


# lin pred
lin_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    add_linpred_draws(model)
p <- lin_pred_draws %>% ggplot(aes(x = .linpred)) +
    geom_density()
mysaveandstore(sprintf("%s/%s/linpred.pdf", outputdir, model_name), 10, 4)



# posterior predictions
post_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    mutate(cov = 30) %>%
    add_predicted_draws(model) %>%
    mutate(.prediction = .prediction / 30)
p <- post_pred_draws %>% ggplot(aes(x = .prediction)) +
    geom_density()
mysaveandstore(sprintf("%s/%s/pred.pdf", outputdir, model_name), 10, 4)

p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    scale_color_brewer() +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred.pdf", outputdir, model_name), 10, 4)

### prior looks good now sample
mfsimprior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            a_c[consensus_pos] +
            a_g[gene_id] +
            a_s[sample] +
            b * condition,
        vector[consensus_pos]:a_c ~ dnorm(0, 0.5),
        vector[gene_id]:a_g ~ dnorm(0, 0.5),
        vector[sample]:a_s ~ dnorm(0, 0.5),
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1)
    ),
    data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE
)
mfsim <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            a_c[consensus_pos] +
            a_g[gene_id] +
            a_s[sample] +
            b * condition,
        vector[consensus_pos]:a_c ~ dnorm(0, 0.5),
        vector[gene_id]:a_g ~ dnorm(0, 0.5),
        vector[sample]:a_s ~ dnorm(0, 0.5),
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1)
    ),
    data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE
)



model <- mv2sim1
modelprior <- mv2sim1
model_name <- "mv2"
precis(model %>%
    spread_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]))

model %<>% recover_types(dsim_for_model)

ttestdf <- dsim_for_model %>%
    mutate(pmeth = Nmeth / cov) %>%
    group_by(sample, condition) %>%
    summarize(pmeth = mean(pmeth)) %>%
    mutate(logodsmeth = logit(pmeth))
t.test(pmeth ~ condition, data = ttestdf, var.equal = TRUE)
t.test(logodsmeth ~ condition, data = ttestdf, var.equal = TRUE)



draws_spread <- model %>%
    spread_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
draws_gathered <- model %>%
    gather_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))


draws_spread_prior <- modelprior %>%
    spread_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
draws_gathered_prior <- modelprior %>%
    gather_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
p <- draws_gathered %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye() +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/%s/summary.pdf", outputdir, model_name), 10, 4)

p <- draws_spread %>%
    ggplot(aes(y = consensus_pos, x = a_c)) +
    stat_halfeye() +
    geom_point(aes(x = consensus_pos_effect), data = dsim %>% mutate(consesus_pos = as.factor(consensus_pos)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1))
mysaveandstore(sprintf("%s/%s/summary_ac.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = sample, x = a_s)) +
    stat_halfeye() +
    geom_point(aes(x = sample_effect), data = dsim %>% mutate(sample = as.factor(sample)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_as.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = gene_id, x = a_g)) +
    stat_halfeye() +
    geom_point(aes(x = gene_id_effect), data = dsim %>% mutate(gene_id = as.factor(gene_id)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_ag.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = "condition", x = b)) +
    stat_halfeye() +
    geom_point(aes(x = condition_effect), data = dsim %>% filter(condition == 1), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_b.pdf", outputdir, model_name), 4, 4)


# lin pred
lin_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    add_linpred_draws(model)
p <- lin_pred_draws %>% ggplot(aes(x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred.pdf", outputdir, model_name), 10, 4)

p <- lin_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred_bycondition.pdf", outputdir, model_name), 10, 4)

# posterior predictions
post_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    mutate(cov = 30) %>%
    add_predicted_draws(model) %>%
    mutate(.prediction = .prediction / 30)
p <- post_pred_draws %>% ggplot(aes(x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred.pdf", outputdir, model_name), 10, 4)

p <- post_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred_condition.pdf", outputdir, model_name), 10, 4)


p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    scale_color_brewer() +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred.pdf", outputdir, model_name), 10, 4)



####
mv2sim1prior <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a + z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE, sample_prior = TRUE
)

model_name <- "mv2sim1"
mv2sim1 <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a + z_c[consensus_pos] * a_c_sigma +
            z_g[gene_id] * a_g_sigma +
            z_s[sample] * a_s_sigma +
            b * condition,

        # Non-centered parameters
        vector[consensus_pos]:z_c ~ dnorm(0, 1),
        vector[gene_id]:z_g ~ dnorm(0, 1),
        vector[sample]:z_s ~ dnorm(0, 1),

        # Hyperpriors
        a_c_sigma ~ dexp(1),
        a_g_sigma ~ dexp(1),
        a_s_sigma ~ dexp(1),
        # Prior for b
        b ~ dnorm(0, 0.5),
        a ~ dnorm(0, 1),

        # Generated quantities for easier interpretation
        gq > vector[consensus_pos]:a_c <<- z_c * a_c_sigma,
        gq > vector[gene_id]:a_g <<- z_g * a_g_sigma,
        gq > vector[sample]:a_s <<- z_s * a_s_sigma
    ),
    data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)


# THIS MODEL EFFECTIVELY INCLUDED BOTH LEFT AND RIGHT LEG PREDICTORS - CAN'T HAVE MULTIPLE A BARS.
# mv2sim <- ulam(
#     alist(
#         Nmeth ~ dbinom(cov, p),
#         logit(p) <- a + a_c_bar + z_c[consensus_pos] * a_c_sigma +
#             a_g_bar + z_g[gene_id] * a_g_sigma +
#             a_s_bar + z_s[sample] * a_s_sigma +
#             b * condition,

#         # Non-centered parameters
#         vector[consensus_pos]:z_c ~ dnorm(0, 1),
#         vector[gene_id]:z_g ~ dnorm(0, 1),
#         vector[sample]:z_s ~ dnorm(0, 1),

#         # Hyperpriors
#         a_c_bar ~ dnorm(0, 0.5),
#         a_c_sigma ~ dexp(1),
#         a_g_bar ~ dnorm(0, 0.5),
#         a_g_sigma ~ dexp(1),
#         a_s_bar ~ dnorm(0, 0.5),
#         a_s_sigma ~ dexp(1),
#         # Prior for b
#         b ~ dnorm(0, 0.5),
#         a ~ dnorm(0, 1),

#         # Generated quantities for easier interpretation
#         gq > vector[consensus_pos]:a_c <<- a_c_bar + z_c * a_c_sigma,
#         gq > vector[gene_id]:a_g <<- a_g_bar + z_g * a_g_sigma,
#         gq > vector[sample]:a_s <<- a_s_bar + z_s * a_s_sigma
#     ),
#     data = dsim_for_model, chains = 4, cores = 4, log_lik = FALSE
# )


model <- mv2sim1
modelprior <- mv2sim1prior
model_name <- "mv2sim1"
precis(model %>%
    spread_draws(a, b, a_c[consensus_pos], a_g[gene_id], a_s[sample]))

model %<>% recover_types(dsim_for_model)

draws_spread <- model %>%
    spread_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
draws_gathered <- model %>%
    gather_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))


draws_spread_prior <- modelprior %>%
    spread_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
draws_gathered_prior <- modelprior %>%
    gather_draws(a, b, a_c_sigma, a_c[consensus_pos], a_g_sigma, a_g[gene_id], a_s_sigma, a_s[sample]) %>%
    mutate(consensus_pos = as.factor(consensus_pos)) %>%
    mutate(gene_id = as.factor(gene_id)) %>%
    mutate(sample = as.factor(sample))
p <- draws_gathered %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye() +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/%s/summary.pdf", outputdir, model_name), 10, 4)

p <- draws_spread %>%
    ggplot(aes(y = consensus_pos, x = a_c)) +
    stat_halfeye() +
    geom_point(aes(x = consensus_pos_effect), data = dsim %>% mutate(consesus_pos = as.factor(consensus_pos)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1))
mysaveandstore(sprintf("%s/%s/summary_ac.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = sample, x = a_s)) +
    stat_halfeye() +
    geom_point(aes(x = sample_effect), data = dsim %>% mutate(sample = as.factor(sample)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_as.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = gene_id, x = a_g)) +
    stat_halfeye() +
    geom_point(aes(x = gene_id_effect), data = dsim %>% mutate(gene_id = as.factor(gene_id)), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_ag.pdf", outputdir, model_name), 4, 4)

p <- draws_spread %>%
    ggplot(aes(y = "condition", x = b)) +
    stat_halfeye() +
    geom_point(aes(x = condition_effect), data = dsim %>% filter(condition == 1), color = "red") +
    stat_halfeye(data = draws_spread_prior, alpha = 0.2, color = "blue") +
    lims(x = c(-1, 1)) +
    mtclosed
mysaveandstore(sprintf("%s/%s/summary_b.pdf", outputdir, model_name), 4, 4)


b_samples <- draws_spread %>% pull(b)
left_tail <- mean(b_samples < 0)
right_tail <- mean(b_samples > 0)
two_tailed_p_value_alt <- 2 * min(left_tail, right_tail)

# lin pred
lin_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    add_linpred_draws(model)
p <- lin_pred_draws %>% ggplot(aes(x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred.pdf", outputdir, model_name), 10, 4)

p <- lin_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .linpred)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/linpred_bycondition.pdf", outputdir, model_name), 10, 4)

# posterior predictions
post_pred_draws <- dsim_for_model %>%
    data_grid(consensus_pos, condition, gene_id, sample) %>%
    mutate(cov = 30) %>%
    add_predicted_draws(model) %>%
    mutate(.prediction = .prediction / 30)
p <- post_pred_draws %>% ggplot(aes(x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred.pdf", outputdir, model_name), 10, 4)

p <- post_pred_draws %>% ggplot(aes(y = as.factor(condition), x = .prediction)) +
    stat_halfeye()
mysaveandstore(sprintf("%s/%s/pred_condition.pdf", outputdir, model_name), 10, 4)


p <- post_pred_draws %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    ggplot(aes(x = consensus_pos, y = .prediction)) +
    stat_interval(.width = c(.50, .80, .95, .99)) +
    scale_color_brewer() +
    mtclosed
mysaveandstore(sprintf("%s/%s/postpred.pdf", outputdir, model_name), 10, 4)


############# READ MODEL

region = "L1HS_intactness_req_ALL"
mod_code_var = "m"
distance_from_start_to_probe = 909
consider_reads_spanning_fraction = 0.7
fraction_meth_threshold = 0.5
context = "CpG"

rdf <- read_delim(sprintf("bayes/data/reads_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.tsv", 
 region, distance_from_start_to_probe, consider_reads_spanning_fraction, 
 fraction_meth_threshold, mod_code_var, context), col_names = TRUE
)
rdf <- rdf %>% mutate(highly_methylated = ifelse(fraction_meth > 0.5, 1, 0)) %>% 
    group_by(sample) %>%
    summarise(Nmeth = sum(highly_methylated), cov = n())

rdf$Nmeth / rdf$cov
rdf <- rdf %>% left_join(sample_table %>% dplyr::select(sample_name, condition, sex, age_scaled, braak, apoe) %>% dplyr::rename(sample = sample_name))

d1 <- rdf %>%
    mutate(
        cov = as.integer(cov),
        Nmeth = as.integer(Nmeth),
        sample = as.integer(factor(sample)),
        condition = as.integer(as.integer(factor(condition, levels = c("CTRL", "AD"))) -1),
        sex = as.integer(as.integer(factor(sex, levels = c("F", "M")))-1),
        braak = as.integer(braak),
        apoe = as.integer(apoe)
    )


dat <- list(
    Nmeth = d1$Nmeth,
    cov = d1$cov,
    sample = d1$sample,
    braak = d1$braak, # BRAAK as an ordered index
    condition = d1$condition,
    apoe = d1$apoe,
    sex = d1$sex,
    age_scaled = d1$age_scaled,
    alpha = rep(2, max(d1$braak) - 1) # Dirichlet prior
)

model_name <- "r_braak"
r_braak <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_s[sample] * a_s_sigma +
            bB * sum(delta_j[1:braak]),
        
        # Non-centered parameters
        vector[sample]:z_s ~ dnorm(0, 1),
        
        # Hyperpriors
        a_s_sigma ~ dexp(1),
        
        # Coefficient priors
        bB ~ normal(0, 1),
        a ~ normal(0, 1),
        
        # Delta parameters for BRAAK stage
        vector[8]:delta_j <<- append_row(0, delta), # delta_j includes a fixed 0
        simplex[7]:delta ~ dirichlet(alpha)     # Dirichlet prior on delta
    ),
    data = dat,
    chains = 4, cores = 4, threads = 1, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model <- r_braak
model_name = "r_braak"
precis(model %>%
    spread_draws(a, bB))


model_name <- "r_braak_cov"
r_braak_cov <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_s[sample] * a_s_sigma +
            bB * sum(delta_j[1:braak]) +
            beta_sex * sex +           # Effect of sex
            beta_age * age_scaled,    
        # Non-centered parameters
        vector[sample]:z_s ~ dnorm(0, 1),
        
        # Hyperpriors
        a_s_sigma ~ dexp(1),
        
        # Coefficient priors
        bB ~ normal(0, 1),
        a ~ normal(0, 1),
        beta_sex ~ dnorm(0, 0.5),      # Prior for sex effect
        beta_age ~ dnorm(0, 0.5),      # Prior for age effect
        # Delta parameters for BRAAK stage
        vector[8]:delta_j <<- append_row(0, delta), # delta_j includes a fixed 0
        simplex[7]:delta ~ dirichlet(alpha)     # Dirichlet prior on delta
    ),
    data = dat,
    chains = 4, cores = 4, threads = 1, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model <- r_braak_cov
model_name = "r_braak_cov"
precis(model %>%
    spread_draws(a, bB))

model_name <- "r_ad_cov"
r_ad_cov <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_s[sample] * a_s_sigma +
            z_e[apoe] * a_e_sigma +
            b * condition +
            beta_sex * sex +           # Effect of sex
            beta_age * age_scaled,    
        # Non-centered parameters
        vector[sample]:z_s ~ dnorm(0, 1),
        vector[apoe]:z_e ~ dnorm(0, 1),

        # Hyperpriors
        a_s_sigma ~ dexp(1),
        a_e_sigma ~ dexp(1),

        # Coefficient priors
        b ~ normal(0, 0.5),
        a ~ normal(0, 1),
        beta_sex ~ dnorm(0, 0.5),      # Prior for sex effect
        beta_age ~ dnorm(0, 0.5)      # Prior for age effect
    ),
    data = dat,
    chains = 4, cores = 4, threads = 1, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model <- r_ad_cov
model_name = "r_ad_cov"
precis(model %>%
    spread_draws(a, b, beta_sex, beta_age, z_s[..], z_e[..]))

b_samples <- model %>%
    spread_draws(a, b) %>% pull(b)
left_tail <- mean(b_samples < 0)
right_tail <- mean(b_samples > 0)
one_tailed_pval <- min(left_tail, right_tail)
two_tailed_pval <- 2 * min(left_tail, right_tail)
statsframe <- tibble(parameter = c("b"), one_tailed_pval = one_tailed_pval, two_tailed_pval = two_tailed_pval)
path <- sprintf("%s/%s/stats_b.tsv", outputdir, model_name)
dir.create(dirname(path), recursive = TRUE)
write_delim(statsframe, path)

model_name <- "r_ad"
r_ad <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            z_s[sample] * a_s_sigma +
            b * condition,
        # Non-centered parameters
        vector[sample]:z_s ~ dnorm(0, 1),
        
        # Hyperpriors
        a_s_sigma ~ dexp(1),
        
        # Coefficient priors
        b ~ normal(0, 0.5),
        a ~ normal(0, 1)
    ),
    data = dat,
    chains = 4, cores = 4, threads = 1, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model <- r_ad
model_name = "r_ad"
precis(model %>%
    spread_draws(a, b, z_s[..]))

b_samples <- model %>%
    spread_draws(a, b) %>% pull(b)
left_tail <- mean(b_samples < 0)
right_tail <- mean(b_samples > 0)
one_tailed_pval <- min(left_tail, right_tail)
two_tailed_pval <- 2 * min(left_tail, right_tail)
statsframe <- tibble(parameter = c("b"), one_tailed_pval = one_tailed_pval, two_tailed_pval = two_tailed_pval)
path <- sprintf("%s/%s/stats_b.tsv", outputdir, model_name)
dir.create(dirname(path), recursive = TRUE)
write_delim(statsframe, path)



model_name <- "r_ad_nosample"
r_ad_nosample <- ulam(
    alist(
        Nmeth ~ dbinom(cov, p),
        logit(p) <- a +
            b * condition +
            beta_sex * sex +           # Effect of sex
            beta_age * age_scaled,    
        
        # Coefficient priors
        b ~ normal(0, 0.5),
        a ~ normal(0, 1),
        beta_sex ~ dnorm(0, 0.5),      # Prior for sex effect
        beta_age ~ dnorm(0, 0.5)      # Prior for age effect
    ),
    data = dat,
    chains = 4, cores = 4, threads = 1, log_lik = FALSE,
    file = sprintf("%s/%s", outputdir, model_name),
    output_dir = stan_workdir
)
model <- r_ad_nosample
model_name = "r_ad_nosample"
precis(model %>%
    spread_draws(a, b, beta_sex, beta_age))

b_samples <- model %>%
    spread_draws(a, b) %>% pull(b)
left_tail <- mean(b_samples < 0)
right_tail <- mean(b_samples > 0)
one_tailed_pval <- min(left_tail, right_tail)
two_tailed_pval <- 2 * min(left_tail, right_tail)
statsframe <- tibble(parameter = c("b"), one_tailed_pval = one_tailed_pval, two_tailed_pval = two_tailed_pval)
path <- sprintf("%s/%s/stats_b.tsv", outputdir, model_name)
dir.create(dirname(path), recursive = TRUE)
write_delim(statsframe, path)


ttt <- rdf %>% mutate(pct = Nmeth/cov)
summary(lm(pct ~ condition + sex + age_scaled, data = ttt))
