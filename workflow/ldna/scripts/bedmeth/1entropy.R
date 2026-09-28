module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")

cell_fractions <- read_delim("ldna/results/m/tables/scMD_cell_type_fractions_all.csv")
cell_fractions_scaled <- cell_fractions %>%
    mutate(Neuron = Inh + Exc) %>%
    mutate(across(all_of(c(cell_types_cols, "Neuron")), ~ as.numeric(scale(.)), .names = "{.col}_z"))

sample_table <- sample_table %>% left_join(cell_fractions_scaled)
sample_table <- sample_table %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = sample_name)
samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(zoo)
library(pryr)
library(circlize)
library(rGREAT)
library(reactome.db)
library(msigdb)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)
# library(ReMapEnrich)
library(msigdbr)
library(Biostrings)
library(ggpubr)
library(PCAtools)
library(patchwork)
library(ggh4x)
library(ggnewscale)
library(betareg)
library(scales)
library(ggnewscale)
library(glmmTMB)
library(broom.mixed)

conditions <- conf$levels
contrasts <- conf$contrasts

# Helper to parse a contrast string into condition1 (reference) and condition2 (test)
parse_contrast <- function(contrast_string) {
    parts <- str_match(contrast_string, "^condition_(.+)_vs_(.+)$")
    list(condition2 = parts[1, 2], condition1 = parts[1, 3])
}

# For backwards compatibility, set condition1/condition2 from first contrast
first_contrast <- parse_contrast(contrasts[[1]])
condition1 <- first_contrast$condition1
condition2 <- first_contrast$condition2
condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
enough_samples_per_condition_for_stats <- ifelse(length(condition1samples) < 3 | length(condition2samples) < 3, FALSE, TRUE)

adjustment_set <- c(conf$linear_model_adjustment_set)
lm_right_hand_side <- ifelse(is.null(adjustment_set), "condition", paste0(c("condition", adjustment_set), collapse = " + "))
adjustment_set_categorical <- if (is.null(adjustment_set)) {
    NULL
} else {
    sample_table %>%
        dplyr::select(all_of(adjustment_set)) %>%
        map(class) %>%
        unlist() %>%
        grep(pattern = "character|factor", value = TRUE) %>%
        names()
}

asc <- ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1])


{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")
    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
    if ("chrY" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, grep("_chrX_|_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes)
        } else {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, grep("_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX)
        }
    } else if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrY, grep("_chrX_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrY)
    } else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
    }
}
#################### functions and themes
named_group_split <- function(.tbl, ...) {
    grouped <- group_by(.tbl, ...)
    names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

    grouped %>%
        group_split() %>%
        rlang::set_names(names)
}
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", sample_table$sample_name, sample_table$sample_name),
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/m/tables/dmrs.tsv",
            dmls = "ldna/results/m/tables/dmls.tsv"
        ), env = globalenv())
        assign("params", list(
            contrasts = conf$contrasts,
            mod_code = "m"
        ), env = globalenv())
        assign("outputs", list(
            promoters_bed = "ldna/Rintermediates/m/promoters_t05.bed",
            dmrpromoterhyper_bed = "ldna/Rintermediates/m/promoters_dmhyperregions_t05.bed",
            dmrpromoterhypo_bed = "ldna/Rintermediates/m/promoters_dmhyporegions_t05.bed"
        ), env = globalenv())
    }
)


merge_with_grs <- function(grs, rte_frame) {
    mbo <- mergeByOverlaps(grs, rte_frame)
    methdf <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    rte_only_frame <- mbo$rte_frame %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
    rtedf_promoters <- bind_cols(methdf, rte_only_frame)
    return(rtedf_promoters)
}

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

# rmannextended <- get_repeat_annotations(
#     default_or_extended = "default",
#     keep_non_central = FALSE
# )

# rmannextended %>% filter(rte_subfamily == "L1HS") %>% filter(refstatus == "Ref") %>% filter(intactness_req == "Intact")
rmannextended <- get_repeat_annotations(
    default_or_extended = "extended",
    keep_non_central = FALSE
)






##########################
# covariate / confounder analysis for per-sample mean entropy
# Shannon entropy estimated from reads is biased downward at low depth (fewer reads -> fewer distinct
# methylation patterns observed), so a condition with lower coverage will look less entropic. These
# helpers check whether condition differences survive adjustment for mean_num_reads and for the
# sample_table covariates (sex, age, ancestry, apoe, cell type fractions), and how balanced those
# covariates are across conditions in the first place.
library(ggrepel)
entropy_covariates <- c("mean_num_reads", "age", "sex", "ancestry", "apoe", paste0(c(cell_types_cols, "Neuron"), "_z"))
# covariates for the jointly adjusted model. age is left out because it is (nearly) collinear with condition in
# any contrast involving YNG; Inh/Exc are left out because they sum to Neuron.
entropy_full_model_covariates <- c("mean_num_reads", "sex", "Neuron_z", "Oligo_z")

entropy_covariate_analysis <- function(df, out_dir, title) {
    # df: one row per sample with sample_name, condition, me, mean_num_reads
    df <- df %>%
        dplyr::select(sample_name, condition, me, mean_num_reads) %>%
        mutate(sample_name = as.character(sample_name)) %>%
        left_join(sample_table %>% mutate(sample_name = as.character(sample_name)) %>% dplyr::select(-any_of(c("condition", "mean_num_reads", "sample"))), by = "sample_name") %>%
        mutate(condition = fct_drop(condition))
    ref_condition <- levels(df$condition)[1]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # only covariates that actually vary among these samples
    covs <- entropy_covariates[map_lgl(entropy_covariates, ~ .x %in% names(df) && n_distinct(df[[.x]], na.rm = TRUE) > 1)]
    num_covs <- covs[map_lgl(covs, ~ is.numeric(df[[.x]]))]
    cat_covs <- setdiff(covs, num_covs)

    ## 1. are the covariates balanced across conditions?
    covariate_balance <- map_dfr(covs, function(cv) {
        if (cv %in% num_covs) {
            fit <- lm(reformulate("condition", cv), df)
            tibble(covariate = cv, test = "anova", p.value = anova(fit)$`Pr(>F)`[1], cor_with_me = cor(df$me, df[[cv]], use = "complete.obs"))
        } else {
            tibble(covariate = cv, test = "fisher", p.value = fisher.test(table(df$condition, df[[cv]]))$p.value, cor_with_me = NA_real_)
        }
    }) %>% arrange(p.value)
    covariate_means <- df %>%
        group_by(condition) %>%
        summarise(n = n(), across(all_of(num_covs), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    write_delim(covariate_balance, file.path(out_dir, "covariate_balance.tsv"), delim = "\t")
    write_delim(covariate_means, file.path(out_dir, "covariate_means_by_condition.tsv"), delim = "\t")

    ## 2. condition effect under a series of adjustments
    model_rhs <- c(
        list(unadjusted = "condition"),
        setNames(map(covs, ~ c("condition", .x)), paste0("+", covs)),
        setNames(map(setdiff(covs, "mean_num_reads"), ~ c("condition", "mean_num_reads", .x)), paste0("+mean_num_reads+", setdiff(covs, "mean_num_reads"))),
        list(full = c("condition", intersect(entropy_full_model_covariates, covs)))
    )
    coef_tbl <- imap_dfr(model_rhs, function(rhs, name) {
        fit <- lm(reformulate(rhs, "me"), df)
        broom::tidy(fit) %>%
            mutate(model = name, formula = paste("me ~", paste(rhs, collapse = " + ")), n = nobs(fit), df.residual = fit$df.residual, adj.r.squared = summary(fit)$adj.r.squared, .before = 1)
    })
    write_delim(coef_tbl, file.path(out_dir, "covariate_models.tsv"), delim = "\t")

    cond_terms <- coef_tbl %>%
        filter(grepl("^condition", term)) %>%
        mutate(term = paste0(gsub("^condition", "", term), " vs ", ref_condition)) %>%
        mutate(model = factor(model, levels = rev(names(model_rhs))))
    p <- cond_terms %>%
        ggplot(aes(x = estimate, y = model)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_pointrange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error, color = p.value < 0.05)) +
        scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), name = "p < 0.05") +
        facet_wrap(~term, scales = "free_x") +
        labs(x = "Condition effect on mean entropy (95% CI)", y = "", title = title) +
        mtclosed
    mysaveandstore(file.path(out_dir, "condition_effect_by_adjustment.pdf"), 4 + 2.5 * n_distinct(cond_terms$term), 1 + 0.25 * length(model_rhs), pl = p)

    ## 3. entropy vs read depth, with per-condition fits, and the group comparison at equal depth
    p <- df %>%
        ggplot(aes(x = mean_num_reads, y = me)) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed") +
        geom_smooth(aes(color = condition), method = "lm", formula = y ~ x, se = FALSE) +
        geom_point(aes(color = condition), size = 2) +
        geom_text_repel(aes(label = sample_name), size = 2.5, max.overlaps = Inf) +
        scale_conditions +
        labs(x = "Mean reads per window", y = "Mean entropy", title = title) +
        mtclosed
    mysaveandstore(file.path(out_dir, "entropy_vs_reads.pdf"), 5, 4, pl = p)

    fit_reads <- lm(me ~ condition + mean_num_reads, df)
    df <- df %>% mutate(me_reads_adj = me - coef(fit_reads)[["mean_num_reads"]] * (mean_num_reads - mean(mean_num_reads)))
    p <- df %>%
        pivot_longer(c(me, me_reads_adj), names_to = "measure", values_to = "value") %>%
        mutate(measure = factor(recode(measure, me = "Raw", me_reads_adj = "Read-depth adjusted"), levels = c("Raw", "Read-depth adjusted"))) %>%
        ggplot(aes(x = condition, y = value)) +
        geom_violin(aes(fill = condition), color = "black", show.legend = FALSE) +
        geom_quasirandom(width = 0.15, size = 1.2) +
        stat_summary(fun.data = mean_se, geom = "pointrange", color = "black", size = 0.3) +
        facet_wrap(~measure) +
        scale_conditions +
        labs(x = "", y = "Mean entropy", title = title) +
        mtclosed
    mysaveandstore(file.path(out_dir, "entropy_raw_vs_reads_adjusted.pdf"), 5, 4, pl = p)

    ## 4. entropy vs every numeric covariate
    p <- df %>%
        pivot_longer(all_of(num_covs), names_to = "covariate", values_to = "value") %>%
        ggplot(aes(x = value, y = me)) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed") +
        geom_point(aes(color = condition), size = 1.5) +
        facet_wrap(~covariate, scales = "free_x") +
        scale_conditions +
        labs(x = "", y = "Mean entropy", title = title) +
        mtclosed
    mysaveandstore(file.path(out_dir, "entropy_vs_covariates.pdf"), 9, 7, pl = p)

    invisible(list(coefficients = coef_tbl, balance = covariate_balance))
}

###########################




# global entropy
entropydfs <- list()
for (sample in sample_table$sample_name) {
    enttemp <- read_delim(str_glue("ldna/results/m/tables/entropy/{sample}/regions.bed"))
    enttemp$sample_name <- sample
    entropydfs[[sample]] <- enttemp
}
entdfglobal <- entropydfs %>%
    purrr::reduce(bind_rows) %>%
    distinct(sample_name, region_name, .keep_all = TRUE)
n_entropy_samples_global <- n_distinct(entdfglobal$sample_name)

pf <- entdfglobal %>%
    group_by(region_name) %>%
    mutate(n_represented = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_represented == n_entropy_samples_global) %>%
    left_join(sample_table) %>%
    group_by(sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    ungroup()
write_delim(pf, str_glue("ldna/results/{params$mod_code}/tables/entropy_global_per_sample.tsv"), delim = "\t")


entropy_covariate_results <- list()
entropy_covariate_results[["global"]] <- entropy_covariate_analysis(pf, str_glue("ldna/results/{params$mod_code}/plots/entropy_global/covariates"), "Genome-wide entropy")


# dmr entropy
entropydfs <- list()
for (contrast1 in contrasts) {
    for (subset in c("t05", "hyper_t05", "hypo_t05", "t01", "hyper_t01", "hypo_t01")) {
        for (sample in sample_table$sample_name) {
            enttemp <- read_delim(str_glue("ldna/results/{params$mod_code}/tables/entropy_ds15_{contrast1}/{subset}/{sample}/regions.bed"))
            enttemp$sample_name <- sample
            enttemp$contrast <- contrast1
            enttemp$type <- subset
            entropydfs[[paste0(sample, contrast1, subset)]] <- enttemp
        }
    }
}
# some entropy files contain every region twice (duplicate runs appended, differing only in float rounding),
# so keep one row per region per sample
entdf <- entropydfs %>%
    purrr::reduce(bind_rows) %>%
    distinct(contrast, type, sample_name, region_name, .keep_all = TRUE)
n_entropy_samples <- n_distinct(entdf$sample_name)


entdf %>%
    group_by(contrast, type, region_name) %>%
    summarise(n_represented = n_distinct(sample_name), .groups = "drop") %>%
    count(contrast, type, n_represented)

# per-sample mean entropy over the DMRs of each contrast, keeping only DMRs represented in every sample and
# only the two conditions that belong to the contrast
contrast_conditions <- lapply(contrasts, function(x) unlist(parse_contrast(x))) %>% setNames(contrasts)
dmr_entropy_per_sample <- entdf %>%
    mutate(
        threshold = str_extract(type, "t0[15]$"),
        direction = case_when(grepl("^hyper", type) ~ "Hyper", grepl("^hypo", type) ~ "Hypo", TRUE ~ "All")
    ) %>%
    group_by(contrast, type, region_name) %>%
    mutate(n_represented = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_represented == n_entropy_samples) %>%
    left_join(sample_table) %>%
    group_by(contrast, threshold, direction, sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads), n_dmrs = n(), .groups = "drop") %>%
    rowwise() %>%
    filter(condition %in% contrast_conditions[[contrast]]) %>%
    ungroup() %>%
    mutate(condition = fct_drop(condition), contrast = factor(contrast, levels = contrasts))

plot_dmr_entropy <- function(df, title) {
    ggplot(df, aes(x = condition, y = me)) +
        geom_violin(aes(fill = condition), color = "black", show.legend = FALSE) +
        geom_quasirandom(width = 0.15, size = 1.2) +
        stat_summary(fun.data = mean_se, geom = "pointrange", color = "black", size = 0.3) +
        labs(x = "", y = "Mean DMR entropy") +
        ggtitle(title) +
        scale_conditions +
        mtclosed
}

for (threshold in c("t05", "t01")) {
    for (c1 in contrasts) {
        df_c1 <- dmr_entropy_per_sample %>% filter(threshold == !!threshold, contrast == c1)
        if (nrow(df_c1) == 0) next
        cp <- parse_contrast(c1)
        df_c1 <- df_c1 %>% mutate(condition = factor(condition, levels = c(cp$condition1, cp$condition2)))

        p <- plot_dmr_entropy(df_c1 %>% filter(direction == "All"), sprintf("DMR entropy (%s)", c1))
        mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{threshold}/{c1}/mean_entropy.pdf"), 3, 4, pl = p)

        p <- plot_dmr_entropy(df_c1 %>% filter(direction != "All"), sprintf("DMR entropy (%s)", c1)) +
            facet_wrap(~direction)
        mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{threshold}/{c1}/mean_entropy_by_direction.pdf"), 5, 4, pl = p)
    }

    # all contrasts side by side; free x so each panel only shows its own two conditions
    df_all <- dmr_entropy_per_sample %>% filter(threshold == !!threshold)
    if (nrow(df_all) == 0) next
    p <- plot_dmr_entropy(df_all %>% filter(direction == "All"), "DMR entropy") +
        facet_wrap(~contrast, scales = "free_x")
    mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{threshold}/mean_entropy_per_contrast.pdf"), 2 + 2 * length(contrasts), 4, pl = p)

    p <- plot_dmr_entropy(df_all %>% filter(direction != "All"), "DMR entropy") +
        facet_grid(direction ~ contrast, scales = "free_x")
    mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{threshold}/mean_entropy_per_contrast_by_direction.pdf"), 2 + 2 * length(contrasts), 6.5, pl = p)
}

# does the DMR entropy difference survive adjustment for read depth / covariates?
for (threshold in c("t05", "t01")) {
    for (c1 in contrasts) {
        for (dir1 in c("All", "Hyper", "Hypo")) {
            df_c1 <- dmr_entropy_per_sample %>% filter(threshold == !!threshold, contrast == c1, direction == dir1)
            if (nrow(df_c1) < 4) next
            cp <- parse_contrast(c1)
            df_c1 <- df_c1 %>% mutate(condition = factor(condition, levels = c(cp$condition1, cp$condition2)))
            out_dir <- str_glue("ldna/results/{params$mod_code}/plots/entropy_{threshold}/{c1}/covariates", if (dir1 == "All") "" else paste0("_", tolower(dir1)))
            res <- entropy_covariate_analysis(df_c1, out_dir, sprintf("DMR entropy (%s, %s, %s)", c1, threshold, dir1))
            entropy_covariate_results[[paste(threshold, c1, dir1)]] <- res$coefficients %>% mutate(threshold = threshold, contrast = c1, direction = dir1, .before = 1)
        }
    }
}

entropy_covariate_results[["global"]] <- entropy_covariate_results[["global"]]$coefficients %>% mutate(threshold = "global", contrast = "all", direction = "All", .before = 1)
write_delim(bind_rows(entropy_covariate_results), str_glue("ldna/results/{params$mod_code}/tables/entropy_covariate_models.tsv"), delim = "\t")


sflist <- list()
for (subset in c("t05", "hyper_t05", "hypo_t05", "t01", "hyper_t01", "hypo_t01")) {
    tempdf <- entdf %>% filter(type == subset)

    pf <- tempdf %>%
        group_by(contrast, region_name) %>%
        mutate(n_represented = n_distinct(sample_name)) %>%
        ungroup() %>%
        filter(n_represented == n_entropy_samples) %>%
        left_join(sample_table) %>%
        group_by(contrast, sample_name, condition, sex, age) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ungroup()
    p <- pf %>%
        ggviolin(x = "condition", y = "me", fill = "condition", add = c("mean_se", "dotplot"), facet.by = "contrast") +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{subset}/mean_entropy.pdf"), 5, 5, pl = p)

    for (c1 in contrasts) {
        pf_c1 <- pf %>% filter(contrast == c1)
        sftemp <- broom::tidy(summary(lm(me ~ condition + age + sex, pf_c1)))
        sftemp$type <- subset
        sftemp$contrast <- c1
        sflist[[paste0(subset, c1)]] <- sftemp
    }

    p <- entdf %>%
        group_by(contrast, sample_name) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ggplot() +
        geom_point(aes(x = mean_num_reads, y = me)) +
        facet_wrap(~contrast) +
        mtclosed
    mysaveandstore("zztemp.324111111111111111111111111111111111.pdf")

    pfregions <- tempdf %>%
        group_by(contrast, region_name) %>%
        mutate(n_represented = n_distinct(sample_name)) %>%
        ungroup() %>%
        filter(n_represented == n_entropy_samples) %>%
        left_join(sample_table) %>%
        group_by(contrast, region_name, sample_name, condition, sex, age) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ungroup()

    library(ggridges)

    pfmeans <- pf %>%
        group_by(contrast, condition) %>%
        summarise(me = mean(me))
    p <- pfregions %>%
        group_by(contrast, region_name, sample_name, condition) %>%
        summarise(me = mean(me)) %>%
        ggplot(aes(y = sample_name, x = me, fill = sample_name)) +
        geom_density_ridges2(quantile_lines = TRUE, quantiles = c(0.5)) +
        scale_samples_unique +
        facet_wrap(~contrast) +
        ggtitle("DMR Entropy") +
        mtclosed
    mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{subset}/entropy_dif_density_by_sample_unique111.pdf"), 5, 5, pl = p)


    dif_dfs <- list()
    for (c1 in unique(pfregions$contrast)) {
        cc <- parse_contrast(c1)
        dif_temp <- pfregions %>%
            filter(contrast == c1) %>%
            group_by(region_name, condition) %>%
            summarise(me = mean(me)) %>%
            pivot_wider(names_from = condition, values_from = me) %>%
            mutate(dif = !!sym(cc$condition2) - !!sym(cc$condition1)) %>%
            mutate(contrast = c1)
        dif_dfs[[c1]] <- dif_temp
    }
    dif_df <- bind_rows(dif_dfs)
    p <- dif_df %>%
        ggplot() +
        geom_histogram(aes(x = dif)) +
        geom_vline(xintercept = 0) +
        facet_wrap(~contrast) +
        ggtitle("DMR Entropy Difference") +
        mtclosed
    mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_{subset}/entropy_dif_hist.pdf"), 5, 3.5, pl = p)
}

tempdf <- entdf %>%
    filter(type != "t05") %>%
    filter(type != "t01")
sf <- purrr::reduce(sflist, bind_rows)
pf <- tempdf %>%
    group_by(contrast, type, region_name) %>%
    mutate(n_represented = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_represented == n_entropy_samples) %>%
    left_join(sample_table) %>%
    group_by(contrast, type, sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    mutate(stat_thresh = gsub(".*_", "", type)) %>%
    mutate(direction = gsub("_.*", "", type)) %>%
    ungroup()
p <- pf %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", facet.by = c("contrast", "direction"), add = c("dotplot")) +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_both/entropy_by_direction.pdf"), 6, 4, pl = p, sf = sf)
p <- pf %>%
    filter(direction == "hyper") %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", facet.by = "contrast", add = c("dotplot", "mean_se")) +
    labs(y = "Mean Entropy", x = "Condition") +
    ggtitle("Hyper") +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_both/entropy_hyper.pdf"), 3.5, 4, pl = p, sf = sf)
p <- pf %>%
    filter(direction == "hypo") %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", facet.by = "contrast", add = c("dotplot", "mean_se")) +
    labs(y = "Mean Entropy", x = "Condition") +
    ggtitle("Hypo") +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_both/entropy_hypo.pdf"), 3.5, 4, pl = p, sf = sf)

pfregions <- tempdf %>%
    group_by(contrast, type, region_name) %>%
    mutate(n_represented = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_represented == n_entropy_samples) %>%
    left_join(sample_table) %>%
    group_by(contrast, type, region_name, sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    ungroup()

dif_dfs <- list()
for (c1 in unique(pfregions$contrast)) {
    cc <- parse_contrast(c1)
    dif_temp <- pfregions %>%
        filter(contrast == c1) %>%
        group_by(type, region_name, condition) %>%
        summarise(me = mean(me)) %>%
        pivot_wider(names_from = condition, values_from = me) %>%
        mutate(dif = !!sym(cc$condition2) - !!sym(cc$condition1)) %>%
        mutate(contrast = c1)
    dif_dfs[[c1]] <- dif_temp
}
dif_df <- bind_rows(dif_dfs)
p <- dif_df %>%
    ggplot() +
    geom_histogram(aes(x = dif)) +
    geom_vline(xintercept = 0) +
    facet_wrap(~ contrast + type, scales = "free_y") +
    ggtitle("DMR Entropy Difference") +
    mtclosed
mysaveandstore(str_glue("ldna/results/{params$mod_code}/plots/entropy_t05/entropy_dif_hist1.pdf"), 5, 5, pl = p)

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "1entropy", params$mod_code))
