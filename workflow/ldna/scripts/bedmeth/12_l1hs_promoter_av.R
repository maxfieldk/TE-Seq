module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = sample_name)
samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")

cell_fractions <- read_delim("ldna/results/m/tables/scMD_cell_type_fractions_all.csv")
cell_fractions_scaled <- cell_fractions %>%
    mutate(Neuron = Inh + Exc) %>%
    mutate(across(all_of(c(cell_types_cols, "Neuron")), ~ as.numeric(scale(.)), .names = "{.col}_z"))

sample_table <- sample_table %>% left_join(cell_fractions_scaled)

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





##################

perl1hs_5utr_region <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE) %>% mutate(region = ordered(region, levels = c("110", "328", "500", "909", "ASP")))
perl1hs_5utr_region$sample <- factor(perl1hs_5utr_region$sample, levels = conf$samples)
perl1hs_5utr_region$condition <- factor(perl1hs_5utr_region$condition, levels = conf$levels)

l1hs_intrautr <- read_delim(sprintf("ldna/Rintermediates/%s/l1hs_intrautr.tsv", params$mod_code), col_names = TRUE)
l1hs_intrautr$sample <- factor(l1hs_intrautr$sample, levels = conf$samples)
l1hs_intrautr$condition <- factor(l1hs_intrautr$condition, levels = conf$levels)

###########################


{
    pf <- perl1hs_5utr_region

    p <- pf %>%
        group_by(gene_id, region) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ggplot(aes(x = mean_meth)) +
        geom_histogram(color = "black") +
        facet_wrap(~region, nrow = 2) +
        mtclosedgrid
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_histogram_5utr_regions.pdf", params$mod_code), 5, 5)


    # Define function to compute ECDF percentages
    compute_ecdf_df <- function(df, region, group_vars, breakpoints, by_sample) {
        df |>
            filter(region == region) |>
            group_by(across(all_of(group_vars))) |>
            summarise(mean_meth = mean(mean_meth), .groups = "drop") |>
            pull(mean_meth) |>
            (\(x) {
                pct_below <- ecdf(x)(breakpoints) * 100
                tibble(threshold = breakpoints, percent_below = pct_below)
            })() |>
            mutate(region = region, by_sample = by_sample)
    }

    # Define regions and breakpoints
    regions <- c("909", "500", "328")
    breakpoints <- seq(0, 100, 5)

    # Compute ECDF for both grouping methods
    ecdf_gene <- map_dfr(regions, ~ compute_ecdf_df(pf, .x, "gene_id", breakpoints, FALSE))
    ecdf_sample <- map_dfr(regions, ~ compute_ecdf_df(pf, .x, c("gene_id", "sample"), breakpoints, TRUE))

    # Combine both results
    ecdf_df <- bind_rows(ecdf_gene, ecdf_sample)
    ecdf_df %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_meth_ecdf.csv", params$mod_code))


    # quantiles
    # Define function to compute quantiles
    compute_quantile_df <- function(df, region, group_vars, probs, by_sample) {
        df |>
            filter(region == region) |>
            group_by(across(all_of(group_vars))) |>
            summarise(mean_meth = mean(mean_meth), .groups = "drop") |>
            pull(mean_meth) |>
            (\(x) {
                quant_values <- quantile(x, probs)
                tibble(quantile = probs, value = quant_values)
            })() |>
            mutate(region = region, by_sample = by_sample)
    }

    # Define regions and quantile probabilities
    regions <- c("909", "500", "328")
    quantile_probs <- seq(0, 1, 0.05) # 5% increments

    # Compute quantiles for both grouping methods
    quantile_gene <- map_dfr(regions, ~ compute_quantile_df(pf, .x, "gene_id", quantile_probs, FALSE))
    quantile_sample <- map_dfr(regions, ~ compute_quantile_df(pf, .x, c("gene_id", "sample"), quantile_probs, TRUE))

    # Combine both results
    quantile_df <- bind_rows(quantile_gene, quantile_sample)
    quantile_df %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_meth_quantiles.csv", params$mod_code))

    p <- pf %>%
        filter(region == 909) %>%
        group_by(sample, condition) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        mutate(condition = factor(condition, levels = conf$levels)) %>%
        ggboxplot(x = "condition", y = "mean_meth", fill = "condition", add = c("mean_se", "dotplot")) +
        scale_conditions +
        mtopen
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/mean_meth_point_fll1hs_box.pdf", params$mod_code), 4.5, 3.75)


    p <- pf %>%
        ggplot(aes(x = region, y = mean_meth, color = condition)) +
        ggbeeswarm::geom_quasirandom(dodge.width = 0.75) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        xlab("") +
        ylab("Average CpG Methylation Per Element") +
        ggtitle("L1HS CpG Methylation") +
        mtopen +
        scale_conditions
    tryCatch(
        {
            mdf <- l1hs_intrautr %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

            senseelement <- mdf %>%
                filter(rte_strand == "+") %$% gene_id %>%
                pluck(1)

            antisenseelement <- mdf %>%
                filter(rte_strand == "-") %$% gene_id %>%
                pluck(1)

            cpgmapping_check <- cg_positions_df %>%
                filter(gene_id == senseelement) %$% sequence_pos %>%
                unique() %>%
                sort()
            methdf_check <- mdf %>%
                filter(gene_id == senseelement) %>%
                relocate(sequence_pos) %$% sequence_pos %>%
                unique() %>%
                sort()
            print(methdf_check)
            print(cpgmapping_check)
            cpgmapping_check <- cg_positions_df %>%
                filter(gene_id == antisenseelement) %$% sequence_pos %>%
                unique() %>%
                sort()
            methdf_check <- mdf %>%
                filter(gene_id == antisenseelement) %>%
                relocate(sequence_pos) %$% sequence_pos %>%
                unique() %>%
                sort()
            print(methdf_check)
            print(cpgmapping_check)

            merged <- left_join(mdf, cg_positions_df, by = c("gene_id", "sequence_pos"))
            cpg_order <- merged %$% consensus_pos %>%
                unique() %>%
                sort()
            merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))

            dat <- merged %>%
                filter(!is.na(pctM))
            dat1 <- merged %>%
                filter(!is.na(consensus_pos))
            merged %>%
                filter(is.na(pctM)) %>%
                pw()
            # write_csv(data_model, "meth_for_sarah.csv")

            data_model <- dat %>%
                mutate(consensus_pos = as.character(consensus_pos)) %>%
                filter(region == "500") %>%
                dplyr::rename(sample_name = sample) %>%
                left_join(sample_table) %>%
                mutate(
                    total_sites = cov,
                    methylated_sites = round(cov * pctM / 100)
                ) %>%
                mutate(age_z = as.numeric(scale(age))) %>%
                mutate(condition = factor(condition, levels = conf$levels))
            # dplyr::select(sample_name, condition, braak, ancestry, cell_types_z_cols, sex, age_z, seqnames, start, gene_id, consensus_pos, total_sites, methylated_sites)

            library(glmmTMB)
            library(broom.mixed)

            global_model <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + ancestry + (1 | sample_name),
                data = data_model,
                family = binomial()
            )

            global_model3 <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name),
                data = data_model,
                family = binomial()
            )

            global_model4 <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name) + (1 | gene_id) + (1 | consensus_pos),
                data = data_model,
                family = binomial()
            )
            global_model5 <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name) + (1 | gene_id) + (1 | consensus_pos),
                data = data_model,
                family = betabinomial()
            )

            broom::tidy(global_model)
            broom::tidy(global_model2) %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909.csv", params$mod_code))
            broom::tidy(global_model3)
            broom::tidy(global_model4) %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909.csv", params$mod_code))
            broom::tidy(global_model5)
            # //ANCHOR - global l1hs stats


            # now data for bayes
            data_model %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909_data.csv", params$mod_code))


            global_model_with_interaction <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition * sex + age_z + (1 | sample_name) + (1 | gene_id) + (1 | consensus_pos),
                data = data_model,
                family = binomial()
            )
            broom::tidy(global_model_with_interaction) %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909_withsexinteraction.csv", params$mod_code))

            gene_model_results <- data_model %>%
                group_by(gene_id) %>%
                group_split() %>%
                map_df(~ {
                    df <- .x
                    if (n_distinct(df$condition) < 2) {
                        return(NULL)
                    } # skip genes with only one condition

                    model <- tryCatch(
                        {
                            glmmTMB(
                                cbind(methylated_sites, total_sites - methylated_sites) ~
                                    condition + sex + age_z + (1 | sample_name) + (1 | consensus_pos),
                                data = df,
                                family = binomial()
                            )
                        },
                        error = function(e) {
                            return(NULL)
                        }
                    )

                    if (!is.null(model)) {
                        est <- summary(model)$coefficients$cond
                        data.frame(
                            gene_id = unique(df$gene_id),
                            condition_effect = est[paste0("condition", condition2), "Estimate"],
                            std_error = est[paste0("condition", condition2), "Std. Error"],
                            z_value = est[paste0("condition", condition2), "z value"],
                            p_value = est[paste0("condition", condition2), "Pr(>|z|)"]
                        )
                    } else {
                        NULL
                    }
                }) %>%
                mutate(p_adj = p.adjust(p_value, method = "fdr")) # adjust for multiple testing

            gene_model_results %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_gene_level_hierarchical_model_909.csv", params$mod_code))

            gene_model_results %>%
                mutate(sig = ifelse(p_adj <= 0.05, TRUE, FALSE)) %>%
                mutate(sig_dir = case_when(
                    condition_effect > 0 & sig == TRUE ~ "sig_up",
                    condition_effect <= 0 & sig == TRUE ~ "sig_down",
                    condition_effect > 0 & sig == FALSE ~ "notsig_up",
                    condition_effect <= 0 & sig == FALSE ~ "notsig_down",
                )) %$% sig_dir %>%
                table()


            ####


            res909 <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 909) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
                broom::tidy() %>%
                mutate(region = "909") %>%
                mutate(model_type = "no_interaction")

            res500 <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 500) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
                broom::tidy() %>%
                mutate(region = "500") %>%
                mutate(model_type = "no_interaction")
            res328 <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 328) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
                broom::tidy() %>%
                mutate(region = "328") %>%
                mutate(model_type = "no_interaction")
            resASP <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == "ASP") %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
                broom::tidy() %>%
                mutate(region = "ASP") %>%
                mutate(model_type = "no_interaction")

            res909interaction <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 909) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
                broom::tidy() %>%
                mutate(region = "909") %>%
                mutate(model_type = "interaction")
            res500interaction <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 500) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
                broom::tidy() %>%
                mutate(region = "500") %>%
                mutate(model_type = "interaction")
            res328interaction <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == 328) %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
                broom::tidy() %>%
                mutate(region = "328") %>%
                mutate(model_type = "interaction")
            resASPinteraction <- pf %>%
                group_by(sample, region) %>%
                summarise(pctM = mean(mean_meth) / 100) %>%
                filter(region == "ASP") %>%
                left_join(sample_table) %>%
                betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
                broom::tidy() %>%
                mutate(region = "ASP") %>%
                mutate(model_type = "interaction")

            stats <- bind_rows(res909, res500, res328, resASP, res909interaction, res500interaction, res328interaction, resASPinteraction)
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4, sf = stats)
        },
        error = function(e) {
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4)
        }
    )


    p <- pf %>%
        left_join(rmannextended) %>%
        filter(region == "909") %>%
        ggplot(aes(x = loc_lowres_integrative_stranded, y = mean_meth, color = condition)) +
        ggbeeswarm::geom_quasirandom(dodge.width = 0.75) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        xlab("") +
        ylab("Average CpG Methylation Per Element") +
        ggtitle("L1HS CpG Methylation") +
        mtopen +
        scale_conditions
    tryCatch(
        {
            res <- pf %>%
                left_join(rmannextended) %>% # Join with rmannextended
                filter(region == "909") %>% # Filter for region 909
                group_by(sample, loc_lowres_integrative_stranded) %>%
                summarise(pctM = mean(mean_meth), .groups = "drop") %>% # Summarize mean methylation
                left_join(sample_table) %>% # Join sample table
                nest(data = -loc_lowres_integrative_stranded) %>% # Nest data for each locus
                mutate(
                    model_results = map(
                        data,
                        ~ lm(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .x) %>%
                            summary() %>%
                            broom::tidy()
                    )
                ) %>%
                dplyr::select(-data) %>% # Drop nested data
                unnest(model_results)
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_gene_loc.pdf", params$mod_code), 5, 4, sf = res)
        },
        error = function(e) {
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_gene_loc.pdf", params$mod_code), 5, 4)
        }
    )
}

########## PCA
if ((conf$single_condition == "no")) {
    library(PCAtools)
    pcaframe <- perl1hs_5utr_region %>%
        filter(region == "909") %>%
        dplyr::select(sample, mean_meth, gene_id) %>%
        mutate(sample = factor(sample, levels = conf$samples)) %>%
        pivot_wider(names_from = gene_id, values_from = mean_meth) %>%
        arrange(sample) %>%
        column_to_rownames(var = "sample") %>%
        as.matrix() %>%
        t()
    pcaframe <- pcaframe[complete.cases(pcaframe), ]

    pcaObj <- pca(pcaframe, center = TRUE, scale = FALSE, metadata = sample_table %>% column_to_rownames(var = "sample_name"))

    p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/screeplot.pdf", params$mod_code), 5, 4)
    p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/screeplot.pdf", params$mod_code), 5, 4)


    # p <- plotloadings(pcaObj,
    #     components = getComponents(pcaObj, seq_len(3)),
    #     rangeRetain = 0.045, labSize = 2
    # ) + mtopen
    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/loadings.pdf", params$mod_code), 5, 4)
    # p <- plotloadings(pcaObj,
    #     components = getComponents(pcaObj, seq_len(3)),
    #     rangeRetain = 0.045, labSize = 2
    # ) + mtopen
    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/loadings.pdf", params$mod_code), 5, 4)

    if (is.null(adjustment_set_categorical)) {
        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", colby = "condition",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen + scale_conditions
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
    } else {
        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1]), colby = "condition",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen + scale_conditions
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
    }
    if (is.null(adjustment_set_categorical)) {
        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", colby = "condition",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen + scale_conditions
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
    } else {
        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1]), colby = "condition",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen + scale_conditions
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
    }





    pf <- pcaframe %>%
        colMeans() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(sample_name = "rowname", mean_meth = ".") %>%
        tibble() %>%
        left_join(sample_table) %>%
        mutate(sample_name = fct_reorder(sample_name, mean_meth))
    p <- pf %>%
        ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
        geom_col() +
        scale_conditions +
        new_scale_fill() +
        geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
        scale_palette +
        mtopen
    # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar1.pdf", params$mod_code), 5, 4)

    tryCatch(
        {
            p <- pf %>%
                mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
                ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
                geom_col() +
                scale_conditions +
                new_scale_fill() +
                geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
                scale_palette +
                mtopen
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage.pdf", params$mod_code), 5, 4)

            p <- pf %>%
                mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
                ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
                geom_col() +
                scale_conditions +
                new_scale_fill() +
                geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
                scale_palette +
                mtopen
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage_ordered.pdf", params$mod_code), 5, 4)


            p <- pf %>%
                mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
                ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
                geom_point(size = 3) +
                scale_conditions +
                geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
                geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
                geom_text_repel(aes(label = apoe)) +
                # new_scale_fill() +
                # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
                scale_conditions +
                mtopen
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_ordered.pdf", params$mod_code), 5, 4)

            p <- pf %>%
                mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
                ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
                geom_point(size = 3) +
                scale_conditions +
                geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
                geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
                geom_text_repel(aes(label = apoe)) +
                # new_scale_fill() +
                # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
                scale_conditions +
                mtopen
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_orderedmeth.pdf", params$mod_code), 5, 4)
        },
        error = function(e) {

        }
    )
}

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "12_l1hs_promoter_av", params$mod_code))
