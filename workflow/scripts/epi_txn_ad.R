module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
set.seed(123)

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(purrr)
library(GenomicRanges)
library(AnnotationDbi)
library(zoo)
library(rtracklayer)
library(Biostrings)
library(ggpubr)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "outputdir" = "integrated/genomebrowserplots/dorado",
            "regions_of_interest" = "conf/integrated_regions_of_interest.bed",
            "r_annotation_fragmentsjoined" = conf[["srna"]]$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf[["srna"]]$r_repeatmasker_annotation,
            mod_code = "m",
            "txdb" = conf[["srna"]]$txdb
        ), env = globalenv())
        assign("inputs", list(
            "srna_results" = "srna/results/agg/deseq/resultsdf.tsv",
            "lrna_results" = "lrna/results/agg/deseq/resultsdf.tsv",
            "ldna_methylation" = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", conf[["ldna"]][["samples"]], conf[["ldna"]][["samples"]]),
            "rte_promoter_meth_by_gene" = "ldna/Rintermediates/m/perelementdf_promoters.tsv",
            "dmrs" = "ldna/results/m/tables/dmrs.CG_m.tsv",
            "dmls" = "ldna/results/m/tables/dmrs.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "plots" = "integrated/epigenome_transcriptome_correlation/plots.rda"
        ), env = globalenv())
    }
)


rmannShared <- read_csv(conf$rmann)
rmannSamples <- list()
for (sample in sample_table$sample_name) {
    df <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
    df$sample_name <- sample
    rmannSamples[[sample]] <- df
}
rmannnonref <- do.call(rbind, rmannSamples) %>% tibble()
RMdf <- bind_rows(rmannShared, rmannnonref) %>%
    filter(refstatus != "NonCentral") %>%
    mutate(gene_id = case_when(
        grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
        TRUE ~ gene_id
    ))

RM <- GRanges(RMdf)

modules_run <- confALL$pipelines_to_deploy

ldna_sample_table <- read_csv("conf/sample_table_ldna.csv")
if ("srna" %in% modules_run) {
    srna_sample_table <- read_csv("conf/sample_table_srna.csv")
    srna_df <- read_delim(inputs$srna_results, delim = "\t") %>% filter(counttype == "telescope_multi")
    colnames(srna_df) <- paste0("srna_", colnames(srna_df)) %>% gsub("srna_gene_id", "gene_id", .)
}
if ("lrna" %in% modules_run) {
    lrna_sample_table <- read_csv("conf/sample_table_lrna.csv")
    lrna_df <- read_delim(inputs$lrna_results, delim = "\t") %>% filter(counttype == "telescope_multi")
    colnames(lrna_df) <- paste0("lrna_", colnames(lrna_df)) %>% gsub("lrna_gene_id", "gene_id", .)
}
if (confALL$ldna$single_condition == "no") {
    if (("lrna" %in% modules_run) & ("srna" %in% modules_run)) {
        lrna_conditions <- lrna_sample_table$condition
        srna_conditions <- srna_sample_table$condition
        ldna_conditions <- ldna_sample_table$condition
        shared_conditions <- intersect(intersect(srna_conditions, lrna_conditions), ldna_conditions)

        srna_samples <- srna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
        lrna_samples <- lrna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
    } else if (("srna" %in% modules_run)) {
        srna_conditions <- srna_sample_table$condition
        ldna_conditions <- ldna_sample_table$condition
        shared_conditions <- intersect(srna_conditions, ldna_conditions)

        srna_samples <- srna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
    } else if (("lrna" %in% modules_run)) {
        lrna_conditions <- lrna_sample_table$condition
        ldna_conditions <- ldna_sample_table$condition

        shared_conditions <- intersect(lrna_conditions, ldna_conditions)
        shared_contrasts <- intersect(lrna_contrasts, ldna_contrasts)
        lrna_samples <- lrna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
    }

    ldna_samples <- ldna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
} else {
    integrated_sample_table <- read_csv("conf/sample_table_integrated.csv")
    ldna_samples <- integrated_sample_table %$% ldna
    if (("lrna" %in% modules_run) & ("srna" %in% modules_run)) {
        srna_samples <- integrated_sample_table %$% srna
        lrna_samples <- integrated_sample_table %$% lrna
    } else if (("srna" %in% modules_run)) {
        srna_samples <- integrated_sample_table %$% srna
    } else if (("lrna" %in% modules_run)) {
        lrna_samples <- integrated_sample_table %$% lrna
    }
}

outputdir <- "integrated/epigenome_transcriptome_correlation"
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

promoter_methdf_by_cpg <- read_delim("ldna/Rintermediates/m/refseq_gene_promoter_methylation.tsv", col_names = TRUE)
promoter_meth_by_gene <- promoter_methdf_by_cpg %>%
    group_by(gene_id, sample) %>%
    summarise(meanmeth = mean(pctM)) %>%
    pivot_wider(names_from = sample, values_from = meanmeth) %>%
    ungroup()
colnames(promoter_meth_by_gene) <- paste0("ldna_", colnames(promoter_meth_by_gene)) %>% gsub("ldna_gene_id", "gene_id", .)

promoter_meth_by_gene_tidy <- pivot_longer(promoter_meth_by_gene, cols = -gene_id, names_to = "sample_name", values_to = "pctM") %>%
    mutate(sample_name = gsub("ldna_", "", sample_name)) %>%
    filter(sample_name %in% ldna_samples)

if (promoter_meth_by_gene_tidy %$% sample_name %>% unique() %>% length() == 1) {
    promoter_meth_by_gene_tidy <- promoter_meth_by_gene_tidy %>% dplyr::select(-sample_name)
}


rte_promoter_meth_by_gene_tidy <- read_delim(inputs$rteprommeth) %>%
    dplyr::rename(sample_name = sample) %>%
    mutate(gene_id = case_when(
        grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
        TRUE ~ gene_id
    ))
rte_promoter_meth_by_gene_tidy <- rte_promoter_meth_by_gene_tidy %>% dplyr::select(gene_id, sample_name, mean_meth)

rte_promoter_meth_by_condition <- rte_promoter_meth_by_gene %>%
    group_by(gene_id, condition) %>%
    summarise(meanmeth = mean(mean_meth)) %>%
    pivot_wider(names_from = condition, values_from = meanmeth) %>%
    ungroup()
# colnames(perrepeatdf) <- paste0("ldna_", colnames(perrepeatdf)) %>% gsub("ldna_gene_id", "gene_id", .)

l1hs_promoterregion_meth_by_gene_tidy <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE) %>%
    dplyr::rename(sample_name = sample) %>%
    mutate(region = ordered(region, levels = c("328", "500", "909", "ASP"))) %>%
    mutate(gene_id = case_when(
        grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
        TRUE ~ gene_id
    )) %>%
    left_join(RMdf %>% dplyr::rename(NI_sample = sample_name))

l1hs_promoterregion_meth_by_gene_tidy_strictly_intergenic <- l1hs_promoterregion_meth_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic")


l1hs_promoterregion_meth_by_sample_tidy <- l1hs_promoterregion_meth_by_gene_tidy %>%
    group_by(sample_name, condition, region) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup()

l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic <- l1hs_promoterregion_meth_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic") %>%
    group_by(sample_name, condition, region) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup()

l1hs_promoterregion_hdmr_by_gene_tidy <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, "L1HS_intactness_req_ALL")) %>%
    dplyr::rename(sample_name = sample) %>%
    mutate(gene_id = case_when(
        grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
        TRUE ~ gene_id
    )) %>%
    left_join(RMdf %>% dplyr::rename(NI_sample = sample_name))

l1hs_promoterregion_hdmr_by_gene_tidy_strictly_intergenic <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic")

l1hs_promoterregion_hdmr_by_sample_tidy <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    group_by(sample_name, condition, subset_threshold) %>%
    summarise(propUnmeth = mean(propUnmeth)) %>%
    ungroup()

l1hs_promoterregion_hdmr_by_sample_tidy_strictly_intergenic <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic") %>%
    group_by(sample_name, condition, subset_threshold) %>%
    summarise(propUnmeth = mean(propUnmeth)) %>%
    ungroup()

srna_gene_expression_tidy <- srna_df %>%
    dplyr::select(gene_id, paste0("srna_", srna_samples)) %>%
    pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "srna_expression") %>%
    mutate(sample_name = gsub("srna_", "", sample_name)) %>%
    mutate(gene_id = case_when(
        grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
        TRUE ~ gene_id
    ))

srna_df

fll1hs_gene_expression_tidy <- srna_gene_expression_tidy %>%
    filter(gene_id %in% (RMdf %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %$% gene_id)) %>%
    left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name))
# srna_gene_expression_tidy %>% filter(grepl("__AS$", gene_id))
fll1hs_gene_expression_tidy_strictly_intergenic <- srna_gene_expression_tidy %>%
    filter(gene_id %in% (RMdf %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %$% gene_id)) %>%
    left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name)) %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic")
# srna_gene_expression_tidy %>% filter(grepl("__AS$", gene_id))

fll1hs_gene_expression_tidy_ASP <- srna_gene_expression_tidy %>%
    filter(grepl("__AS$", gene_id)) %>%
    dplyr::select(sample_name, srna_expression, gene_id) %>%
    dplyr::rename(ASP_counts = srna_expression) %>%
    mutate(gene_id = case_when(
        TRUE ~ gsub("__AS$", "", gene_id)
    )) %>%
    left_join(srna_gene_expression_tidy)
fll1hs_gene_expression_by_sample_tidy_ASP <- fll1hs_gene_expression_tidy_ASP %>%
    group_by(sample_name) %>%
    summarise(srna_expression = mean(srna_expression))

fll1hs_gene_expression_tidy_ASP_strictly_intergenic <- srna_gene_expression_tidy %>%
    filter(grepl("__AS$", gene_id)) %>%
    dplyr::select(sample_name, srna_expression, gene_id) %>%
    dplyr::rename(ASP_counts = srna_expression) %>%
    mutate(gene_id = case_when(
        TRUE ~ gsub("__AS$", "", gene_id)
    )) %>%
    left_join(srna_gene_expression_tidy) %>%
    left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name)) %>%
    filter(loc_superlowres_integrative_stranded_incl_pred == "Intergenic")

fll1hs_gene_expression_by_sample_tidy_ASP_strictly_intergenic <- fll1hs_gene_expression_tidy_ASP_strictly_intergenic %>%
    group_by(sample_name) %>%
    summarise(srna_expression = mean(srna_expression))


fll1hs_gene_expression_by_sample_tidy <- fll1hs_gene_expression_tidy %>%
    group_by(sample_name) %>%
    summarise(srna_expression = mean(srna_expression))
fll1hs_gene_expression_tidy %>% filter(is.na(loc_lowres_integrative_stranded))

fll1hs_gene_expression_by_sample_tidy_strictly_intergenic <- fll1hs_gene_expression_tidy_strictly_intergenic %>%
    group_by(sample_name) %>%
    summarise(srna_expression = mean(srna_expression))
fll1hs_gene_expression_tidy %>% filter(is.na(loc_lowres_integrative_stranded))

# genes
fll1hs_gene_expression_tidy %$% loc_lowres_integrative_stranded %>% unique()

# rtes



# L1HS
mrg_l1hs_by_sample_meth <- l1hs_promoterregion_meth_by_sample_tidy %>% left_join(fll1hs_gene_expression_by_sample_tidy)
mrg_l1hs_by_sample_hdmr <- l1hs_promoterregion_hdmr_by_sample_tidy %>% left_join(fll1hs_gene_expression_by_sample_tidy)

mrg_l1hs_by_sample_meth_strictly_intergenic <- l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic %>% left_join(fll1hs_gene_expression_by_sample_tidy_strictly_intergenic)
mrg_l1hs_by_sample_hdmr_strictly_intergenic <- l1hs_promoterregion_hdmr_by_sample_tidy_strictly_intergenic %>% left_join(fll1hs_gene_expression_by_sample_tidy_strictly_intergenic)


mrg_l1hs_by_gene_meth <- l1hs_promoterregion_meth_by_gene_tidy %>% left_join(fll1hs_gene_expression_tidy)
mrg_l1hs_by_gene_hdmr <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    left_join(fll1hs_gene_expression_tidy) %>%
    mutate(srna_expression = case_when(
        is.na(srna_expression) ~ 0,
        TRUE ~ srna_expression
    ))

mrg_l1hs_by_gene_meth_strictly_intergenic <- l1hs_promoterregion_meth_by_gene_tidy_strictly_intergenic %>% left_join(fll1hs_gene_expression_tidy_strictly_intergenic)
mrg_l1hs_by_gene_hdmr_strictly_intergenic <- l1hs_promoterregion_hdmr_by_gene_tidy_strictly_intergenic %>%
    left_join(fll1hs_gene_expression_tidy_strictly_intergenic) %>%
    mutate(srna_expression = case_when(
        is.na(srna_expression) ~ 0,
        TRUE ~ srna_expression
    ))


##########
### By Sample Correlations
## Methylation
# all elements
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_sample_meth$region %>% unique()) {
    tdf <- mrg_l1hs_by_sample_meth %>%
        filter(region == sthresh) %>%
        left_join(sample_table)
    model <- lm(srna_expression ~ mean_meth, data = tdf)
    model_with_RIN <- lm(srna_expression ~ mean_meth + RIN, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "no RIN")
    stats_with_RIN <- model_with_RIN %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "with RIN")
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    model_stats_with_RIN <- model_with_RIN %>%
        broom::glance() %>%
        mutate(call = model_with_RIN %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    statsframe <- bind_rows(stats, stats_with_RIN)
    modelframe <- bind_rows(model_stats, model_stats_with_RIN)
    p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = statsframe, sfm = modelframe)

    p <- tdf %>% ggplot(aes(x = RIN, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/RIN_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = statsframe, sfm = modelframe)

    model_stats_list[[sthresh]] <- model_stats
}

mf <- model_stats_list %>%
    imap_dfr(~ .x %>%
        dplyr::select(r.squared, p.value) %>%
        mutate(region = .y))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_all_rsq.pdf"), w = 3, h = 4, sf = mf)

# now strictly intergenic
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_sample_meth_strictly_intergenic$region %>% unique()) {
    tdf <- mrg_l1hs_by_sample_meth_strictly_intergenic %>%
        filter(region == sthresh) %>%
        left_join(sample_table)
    model <- lm(srna_expression ~ mean_meth, data = tdf)
    model_with_RIN <- lm(srna_expression ~ mean_meth + RIN, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "no RIN")
    stats_with_RIN <- model_with_RIN %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "with RIN")
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    model_stats_with_RIN <- model_with_RIN %>%
        broom::glance() %>%
        mutate(call = model_with_RIN %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    statsframe <- bind_rows(stats, stats_with_RIN)
    modelframe <- bind_rows(model_stats, model_stats_with_RIN)
    p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_sample_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = statsframe, sfm = modelframe)

    p <- tdf %>% ggplot(aes(x = RIN, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/RIN_by_sample_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

    model_stats_list[[sthresh]] <- model_stats
}

mf <- model_stats_list %>%
    imap_dfr(~ .x %>%
        dplyr::select(r.squared, p.value) %>%
        mutate(region = .y))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_strictly_intergenic_all_rsq.pdf"), w = 3, h = 4, sf = mf)
###########
## HDRs
# all elements
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_sample_hdmr$subset_threshold %>% unique()) {
    tdf <- mrg_l1hs_by_sample_hdmr %>% filter(subset_threshold == sthresh)
    model <- lm(srna_expression ~ propUnmeth, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy()
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

    model_stats_list[[sthresh]] <- model_stats
}

# strictly intergenic
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_sample_hdmr_strictly_intergenic$subset_threshold %>% unique()) {
    tdf <- mrg_l1hs_by_sample_hdmr_strictly_intergenic %>% filter(subset_threshold == sthresh)
    model <- lm(srna_expression ~ propUnmeth, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy()
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

    model_stats_list[[sthresh]] <- model_stats
}


model_stats_list
mf <- model_stats_list %>%
    imap_dfr(~ .x %>%
        dplyr::select(r.squared, p.value) %>%
        mutate(subset = .y)) %>%
    separate(subset, sep = "_", into = c("region", "threshold"))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
    facet_wrap(~region) +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_all_rsq.pdf"), w = 5, h = 4, sf = mf)
##########
### By Gene
### By Sample Correlations
## Methylation
# all elements
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_gene_meth$region %>% unique()) {
    tdf <- mrg_l1hs_by_gene_meth %>% filter(region == sthresh)
    model <- lm(srna_expression ~ mean_meth, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy()
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_gene_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

    model_stats_list[[sthresh]] <- model_stats
}

mf <- model_stats_list %>%
    imap_dfr(~ .x %>%
        dplyr::select(r.squared, p.value) %>%
        mutate(region = .y))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_gene_all_rsq.pdf"), w = 3, h = 4, sf = mf)

# now strictly intergenic
model_stats_list <- list()
for (sthresh in mrg_l1hs_by_gene_meth_strictly_intergenic$region %>% unique()) {
    tdf <- mrg_l1hs_by_gene_meth_strictly_intergenic %>% filter(region == sthresh)
    model <- lm(srna_expression ~ mean_meth, data = tdf)
    stats <- model %>%
        summary() %>%
        broom::tidy()
    model_stats <- model %>%
        broom::glance() %>%
        mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition), size = 3) +
        mtclosed +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

    model_stats_list[[sthresh]] <- model_stats
}

mf <- model_stats_list %>%
    imap_dfr(~ .x %>%
        dplyr::select(r.squared, p.value) %>%
        mutate(region = .y))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_all_rsq.pdf"), w = 3, h = 4, sf = mf)


## HDRs
# All elements

model_stats_list <- list()
model_stats_list_filtered <- list()
for (sthresh in mrg_l1hs_by_gene_hdmr$subset_threshold %>% unique()) {
    for (element_type in mrg_l1hs_by_gene_hdmr %$% loc_superlowres_integrative_stranded_incl_pred %>%
        unique() %>%
        discard(is.na)) {
        tdf <- mrg_l1hs_by_gene_hdmr %>%
            filter(subset_threshold == sthresh) %>%
            filter(loc_lowres_integrative_stranded == element_type)
        tdf %>%
            group_by(gene_id) %>%
            summarise(mean_expr = mean(srna_expression)) %$% mean_expr %>%
            quantile(seq(0, 1, 0.05))
        # idea is to remove those intergenic eleemnts whose epression is likely being driven by adjacent genes - huge inflection point for elements above 90th percentile expression for expression values
        expression_90th_percentile <- tdf %>%
            group_by(gene_id) %>%
            summarise(mean_expr = mean(srna_expression)) %$% mean_expr %>%
            quantile(0.9)
        elements_to_keep <- tdf %>%
            group_by(gene_id) %>%
            summarise(mean_expr = mean(srna_expression)) %>%
            filter(mean_expr < expression_90th_percentile) %$% gene_id

        tdf_filtered <- tdf %>% filter(gene_id %in% elements_to_keep)

        tryCatch(
            {
                model <- lm(srna_expression ~ propUnmeth, data = tdf)
                stats <- model %>%
                    summary() %>%
                    broom::tidy()
                model_stats <- model %>%
                    broom::glance() %>%
                    mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
                p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
                    stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
                    ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                        formula = y ~ x, parse = TRUE, label.y.npc = "top"
                    ) +
                    geom_point(aes(color = condition), size = 3) +
                    mtclosed +
                    scale_conditions +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_gene_{sthresh}_{element_type}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

                model_stats_list[[sthresh]][[element_type]] <- model_stats
            },
            error = function(e) {

            }
        )

        # now for filtered elements
        tryCatch(
            {
                model <- lm(srna_expression ~ propUnmeth, data = tdf_filtered)
                stats <- model %>%
                    summary() %>%
                    broom::tidy()
                model_stats <- model %>%
                    broom::glance() %>%
                    mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
                p <- tdf_filtered %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
                    stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
                    ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                        formula = y ~ x, parse = TRUE, label.y.npc = "top"
                    ) +
                    geom_point(aes(color = condition), size = 3) +
                    mtclosed +
                    scale_conditions +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outputdir}/L1HS_90thpercentilefiltered/hdr_by_gene_{sthresh}_{element_type}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

                model_stats_list_filtered[[sthresh]][[element_type]] <- model_stats
            },
            error = function(e) {

            }
        )
    }
}

mf <- model_stats_list %>%
    map_dfr(
        ~ map_dfr(.x, ~ .x %>%
            dplyr::select(r.squared, p.value), .id = "element_type"),
        .id = "subset"
    ) %>%
    separate(subset, sep = "_", into = c("region", "threshold"))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
    facet_grid(cols = vars(region), rows = vars(element_type)) +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_sample_all_rsq.pdf"), w = 8, h = 8, sf = mf)

mf <- model_stats_list_filtered %>%
    map_dfr(
        ~ map_dfr(.x, ~ .x %>%
            dplyr::select(r.squared, p.value), .id = "element_type"),
        .id = "subset"
    ) %>%
    separate(subset, sep = "_", into = c("region", "threshold"))
sf <- mf
p <- mf %>% ggplot() +
    geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
    facet_grid(cols = vars(region), rows = vars(element_type)) +
    mtclosed +
    anchorbar
mysaveandstore(pl = p, str_glue("{outputdir}/L1HS_90thpercentilefiltered/hdr_by_sample_all_rsq.pdf"), w = 8, h = 8, sf = mf)
#############



model_stats_list <- list()
for (sthresh in mrg_l1hs_by_gene_hdmr$subset_threshold %>% unique()) {
    tdf <- mrg_l1hs_by_gene_hdmr %>% filter(subset_threshold == sthresh)
    gene_ids <- tdf %>%
        group_by(gene_id) %>%
        group_keys() %>%
        pull(gene_id)


    models <- tdf %>%
        group_split(gene_id) %>%
        map2_dfr(gene_ids, ~ {
            model <- lm(srna_expression ~ propUnmeth, data = .x)
            glance_df <- broom::glance(model)
            coef_df <- broom::tidy(model) %>%
                filter(term == "propUnmeth") %>%
                dplyr::select(estimate)
            glance_df %>%
                mutate(
                    coef_estimate = coef_df$estimate,
                    call = deparse(model$call) %>% paste(collapse = " "),
                    gene_id = .y
                )
        }) %>%
        mutate(subset_threshold = sthresh) %>%
        separate(subset_threshold, sep = "_", into = c("region", "threshold"))

    p <- models %>%
        mutate(direction = ifelse(is.na(coef_estimate), "NA", ifelse(coef_estimate > 0, "positive", "negative"))) %>%
        mutate(sig = ifelse(is.na(p.value), "NA", ifelse(p.value <= 0.05, "Sig", "NS"))) %>%
        mutate(direction = ifelse(sig == "NA", "NA", direction)) %>%
        ggplot() +
        geom_bar(aes(x = sig, fill = direction), color = "black") +
        labs(y = "Count", title = str_glue("FL L1HS Expr ~ HDR {sthresh}")) +
        mtclosed +
        anchorbar +
        scale_palette
    mysaveandstore(str_glue("{outputdir}/L1HS/unique_loci/aaa_hdr_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = models)


    for (gene in mrg_l1hs_by_gene_hdmr %>%
        filter(subset_threshold == sthresh) %$% gene_id %>%
        unique()) {
        tdf <- mrg_l1hs_by_gene_hdmr %>%
            filter(subset_threshold == sthresh) %>%
            filter(gene_id == gene)
        model <- lm(srna_expression ~ propUnmeth, data = tdf)
        stats <- model %>%
            summary() %>%
            broom::tidy()
        model_stats <- model %>%
            broom::glance() %>%
            mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
        p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/unique_loci/{gene}_hdr_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)
    }
}
