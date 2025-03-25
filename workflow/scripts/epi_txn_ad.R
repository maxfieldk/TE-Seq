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




ldna_sample_table <- read_csv("conf/sample_table_ldna.csv")
if ("srna" %in% modules_run) {
    srna_sample_table <- read_csv("conf/sample_table_srna.csv") %>% mutate(condition = factor(condition, levels = confALL$srna$levels))
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
    )) %>%
    mutate(loc_superlowres_integrative_stranded_incl_pred2 = case_when(
        loc_superlowres_integrative_stranded_incl_pred == "Genic_Sense" ~ "Genic_Sense",
        loc_superlowres_integrative_stranded_incl_pred == "Genic_Antisense" ~ "Genic_Antisense",
        loc_superlowres_integrative_stranded_incl_pred == "GenicPred_Sense" ~ "GenicPred_Sense",
        loc_superlowres_integrative_stranded_incl_pred == "GenicPred_Antisense" ~ "GenicPred_Antisense",
        loc_integrative == "CdgTxAdj" & coding_tx_adjacent_orientation_loc == "Sense" ~ "TxAdj_Sense",
        loc_integrative == "NoncdgTxAdj" & noncoding_tx_adjacent_orientation_loc == "Sense" ~ "TxAdj_Sense",
        loc_integrative == "CdgTxAdj" & coding_tx_adjacent_orientation_loc == "Antisense" ~ "TxAdj_Antisense",
        loc_integrative == "NoncdgTxAdj" & noncoding_tx_adjacent_orientation_loc == "Antisense" ~ "TxAdj_Antisense",
        loc_superlowres_integrative_stranded_incl_pred == "Intergenic" ~ "Intergenic",
        TRUE ~ "Other"
    ))

RM <- GRanges(RMdf)

gencode <- import("/users/mkelsey/data/Nanopore/alz/gencodeV35hs1.bed", format = "BED")
mcols(gencode) <- mcols(gencode)[, 1:11]
gencodedf <- as.data.frame(gencode) %>% tibble()
fll1hs <- RMdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL")

fll1hs_t2tcoord <- fll1hs %>%
    mutate(
        seq_split = str_split(seqnames, "_"), # Split once and reuse for efficiency
        seqnames = case_when(
            grepl("_NI_", gene_id) ~ map_chr(seq_split, ~ .[3]),
            TRUE ~ seqnames
        ),
        start = case_when(
            grepl("_NI_", gene_id) ~ as.numeric(map_chr(seq_split, ~ .[4])),
            TRUE ~ start
        ),
        end = case_when(
            grepl("_NI_", gene_id) ~ as.numeric(map_chr(seq_split, ~ .[5])),
            TRUE ~ end
        )
    ) %>%
    dplyr::select(-seq_split)

fll1hs_t2tcoordgr <- GRanges(fll1hs_t2tcoord)
fll1hs_genic_sense_gr <- fll1hs_t2tcoordgr %>% subsetByOverlaps(gencode, ignore.strand = FALSE)
fll1hs_genic_gr <- fll1hs_t2tcoordgr %>% subsetByOverlaps(gencode, ignore.strand = TRUE)
fll1hs_genic_antisense_gr <- fll1hs_genic_gr %>% subsetByOverlaps(fll1hs_genic_sense_gr, ignore.strand = FALSE, invert = TRUE)
fll1hs_intergenic_gr <- fll1hs_t2tcoordgr %>% subsetByOverlaps(gencode, ignore.strand = TRUE, invert = TRUE)
fll1hs_nonsense <- c(fll1hs_intergenic_gr, fll1hs_genic_antisense_gr)
fll1hs_nonsense_asp <- c(fll1hs_intergenic_gr, fll1hs_genic_sense_gr)
fll1hs_nonsense_gencode <- mcols(fll1hs_nonsense)$gene_id # note that some are ruled out by refseq, so will have to double filter
fll1hs_nonsense_asp_gencode <- mcols(fll1hs_nonsense_asp)$gene_id


modules_run <- confALL$pipelines_to_deploy


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


# rte_promoter_meth_by_gene_tidy <- read_delim(inputs$rteprommeth) %>%
#     dplyr::rename(sample_name = sample) %>%
#     mutate(gene_id = case_when(
#         grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
#         TRUE ~ gene_id
#     ))
# rte_promoter_meth_by_gene_tidy <- rte_promoter_meth_by_gene_tidy %>% dplyr::select(gene_id, sample_name, mean_meth)

# rte_promoter_meth_by_condition <- rte_promoter_meth_by_gene %>%
#     group_by(gene_id, condition) %>%
#     summarise(meanmeth = mean(mean_meth)) %>%
#     pivot_wider(names_from = condition, values_from = meanmeth) %>%
#     ungroup()
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
    filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Antisense", "GenicPred_Antisense", "TxAdj_Antisense")) %>%
    filter(gene_id %in% fll1hs_nonsense_gencode)
l1hs_promoterregion_meth_by_gene_tidy_strictly_intergenic %$% loc_superlowres_integrative_stranded_incl_pred2 %>% table()

l1hs_promoterregion_meth_by_sample_tidy <- l1hs_promoterregion_meth_by_gene_tidy %>%
    group_by(sample_name, condition, region) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup()

l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic <- l1hs_promoterregion_meth_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Antisense", "GenicPred_Antisense", "TxAdj_Antisense")) %>%
    filter(gene_id %in% fll1hs_nonsense_gencode) %>%
    group_by(sample_name, condition, region) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup()
l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic_sense <- l1hs_promoterregion_meth_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Sense", "GenicPred_Sense", "TxAdj_Sense")) %>%
    filter(gene_id %in% fll1hs_nonsense_asp_gencode) %>%
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
    filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Antisense", "GenicPred_Antisense", "TxAdj_Antisense")) %>%
    filter(gene_id %in% fll1hs_nonsense_gencode)

l1hs_promoterregion_hdmr_by_sample_tidy <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    group_by(sample_name, condition, subset_threshold) %>%
    summarise(propUnmeth = mean(propUnmeth)) %>%
    ungroup()

l1hs_promoterregion_hdmr_by_sample_tidy_strictly_intergenic <- l1hs_promoterregion_hdmr_by_gene_tidy %>%
    filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Antisense", "GenicPred_Antisense", "TxAdj_Antisense")) %>%
    filter(gene_id %in% fll1hs_nonsense_gencode) %>%
    group_by(sample_name, condition, subset_threshold) %>%
    summarise(propUnmeth = mean(propUnmeth)) %>%
    ungroup()

for (tecounttype in c("telescope_unique")) {
    srna_df <- read_delim(inputs$srna_results, delim = "\t") %>% filter(counttype == tecounttype)
    colnames(srna_df) <- paste0("srna_", colnames(srna_df)) %>% gsub("srna_gene_id", "gene_id", .)


    outputdir <- str_glue("integrated/epigenome_transcriptome_correlation/{tecounttype}")
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)


    # prep data


    srna_gene_expression_tidy <- srna_df %>%
        dplyr::select(gene_id, paste0("srna_", srna_samples)) %>%
        pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "srna_expression") %>%
        mutate(sample_name = gsub("srna_", "", sample_name)) %>%
        mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        ))


    fll1hs_gene_expression_tidy <- srna_gene_expression_tidy %>%
        filter(gene_id %in% (RMdf %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %$% gene_id)) %>%
        left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name))
    # srna_gene_expression_tidy %>% filter(grepl("__AS$", gene_id))
    fll1hs_gene_expression_tidy_strictly_intergenic <- srna_gene_expression_tidy %>%
        filter(gene_id %in% (RMdf %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %$% gene_id)) %>%
        left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name)) %>%
        filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Antisense", "GenicPred_Antisense", "TxAdj_Antisense")) %>%
        filter(gene_id %in% fll1hs_nonsense_gencode)


    fll1hs_stringly_intergenic_gene_exp_ranked <- fll1hs_gene_expression_tidy_strictly_intergenic %>%
        group_by(gene_id, seqnames, start, end) %>%
        summarise(srna_expression = mean(srna_expression)) %>%
        ungroup() %>%
        filter(!grepl("_NI_", gene_id)) %>%
        arrange(-srna_expression) %>%
        mutate(rank = row_number())
    write_csv(fll1hs_stringly_intergenic_gene_exp_ranked, str_glue("{outputdir}/fll1hs_gene_expression_tidy_strictly_intergenic_expression_ranked.csv"))

    nelements <- fll1hs_stringly_intergenic_gene_exp_ranked %>% nrow()
    fll1hs_stringly_intergenic_gene_exp_ranked %>%
        filter(gene_id == "L1HS_2p13.2_1") %>%
        mutate(percentile = rank / !!nelements)
    fll1hs_stringly_intergenic_gene_exp_ranked %>%
        filter(gene_id == "L1HS_2q21.1_2") %>%
        mutate(percentile = rank / !!nelements)

    { # needed plot
        pf <- fll1hs_gene_expression_tidy_strictly_intergenic %>%
            group_by(sample_name) %>%
            summarise(sumcounts = sum(srna_expression)) %>%
            left_join(srna_sample_table)
        model <- lm(sumcounts ~ condition, data = pf)
        model_with_RIN <- lm(sumcounts ~ condition + RIN, data = pf)
        model_with_RIN_with_cov <- lm(sumcounts ~ condition + RIN + age + sex, data = pf)
        stats <- model %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "no cov no RIN")
        stats_with_RIN <- model_with_RIN %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "RIN")
        stats_with_RIN_with_cov <- model_with_RIN_with_cov %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "RIN with Cov")
        statsframe <- bind_rows(stats, stats_with_RIN, stats_with_RIN_with_cov)
        p <- pf %>%
            ggplot(aes(x = condition, y = sumcounts)) +
            geom_boxplot(aes(fill = condition)) +
            geom_point() +
            scale_conditions +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HSFLexpression_box.pdf"), w = 4, h = 4, sf = statsframe)
        p <- pf %>%
            ggplot(aes(x = condition, y = sumcounts)) +
            geom_col(data = pf %>% group_by(condition) %>% summarise(sumcounts = mean(sumcounts)), aes(fill = condition), color = "black") +
            geom_point() +
            scale_conditions +
            anchorbar +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HSFLexpression_bar.pdf"), w = 4, h = 4, sf = statsframe)
    }

    # srna_gene_expression_tidy %>% filter(grepl("__AS$", gene_id))


    fll1hs_gene_expression_tidy_ASP <- srna_gene_expression_tidy %>%
        filter(grepl("__AS$", gene_id)) %>%
        dplyr::select(sample_name, srna_expression, gene_id) %>%
        dplyr::rename(ASP_counts = srna_expression) %>%
        mutate(gene_id = case_when(
            TRUE ~ gsub("__AS$", "", gene_id)
        )) %>%
        left_join(srna_gene_expression_tidy)

    fll1hs_gene_expression_tidy_ASP_strictly_intergenic <- srna_gene_expression_tidy %>%
        filter(grepl("__AS$", gene_id)) %>%
        dplyr::select(sample_name, srna_expression, gene_id) %>%
        dplyr::rename(ASP_counts = srna_expression) %>%
        mutate(gene_id = case_when(
            TRUE ~ gsub("__AS$", "", gene_id)
        )) %>%
        left_join(srna_gene_expression_tidy) %>%
        left_join(RMdf %>% dplyr::rename(NI_sample_name = sample_name)) %>%
        filter(loc_superlowres_integrative_stranded_incl_pred2 %in% c("Intergenic", "Genic_Sense", "GenicPred_Sense", "TxAdj_Sense")) %>%
        filter(gene_id %in% fll1hs_nonsense_asp_gencode)

    { # needed plot
        pf <- fll1hs_gene_expression_tidy_ASP_strictly_intergenic %>%
            group_by(sample_name) %>%
            summarise(sumcounts = sum(ASP_counts)) %>%
            left_join(srna_sample_table)
        model <- lm(sumcounts ~ condition, data = pf)
        model_with_RIN <- lm(sumcounts ~ condition + RIN, data = pf)
        model_with_RIN_with_cov <- lm(sumcounts ~ condition + RIN + age + sex, data = pf)
        stats <- model %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "no cov no RIN")
        stats_with_RIN <- model_with_RIN %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "RIN")
        stats_with_RIN_with_cov <- model_with_RIN_with_cov %>%
            summary() %>%
            broom::tidy() %>%
            mutate(type = "RIN with Cov")
        statsframe <- bind_rows(stats, stats_with_RIN, stats_with_RIN_with_cov)
        p <- pf %>%
            ggplot(aes(x = condition, y = sumcounts)) +
            geom_boxplot(aes(fill = condition)) +
            geom_point() +
            scale_conditions +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HSFLASPexpression_box.pdf"), w = 4, h = 4, sf = statsframe)
        p <- pf %>%
            ggplot(aes(x = condition, y = sumcounts)) +
            geom_col(data = pf %>% group_by(condition) %>% summarise(sumcounts = mean(sumcounts)), aes(fill = condition), color = "black") +
            geom_point() +
            scale_conditions +
            anchorbar +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HSFLASPexpression_bar.pdf"), w = 4, h = 4, sf = statsframe)
    }


    fll1hs_gene_expression_by_sample_tidy_ASP <- fll1hs_gene_expression_tidy_ASP %>%
        group_by(sample_name) %>%
        summarise(srna_expression = mean(ASP_counts))



    fll1hs_gene_expression_by_sample_tidy_ASP_strictly_intergenic <- fll1hs_gene_expression_tidy_ASP_strictly_intergenic %>%
        group_by(sample_name) %>%
        summarise(srna_expression = sum(ASP_counts))


    fll1hs_gene_expression_by_sample_tidy <- fll1hs_gene_expression_tidy %>%
        group_by(sample_name) %>%
        summarise(srna_expression = sum(srna_expression))

    fll1hs_gene_expression_by_sample_tidy_strictly_intergenic <- fll1hs_gene_expression_tidy_strictly_intergenic %>%
        group_by(sample_name) %>%
        summarise(srna_expression = sum(srna_expression))

    # genes

    # rtes

    # L1HS
    mrg_l1hs_by_sample_meth <- l1hs_promoterregion_meth_by_sample_tidy %>% left_join(fll1hs_gene_expression_by_sample_tidy)
    mrg_l1hs_by_sample_hdmr <- l1hs_promoterregion_hdmr_by_sample_tidy %>% left_join(fll1hs_gene_expression_by_sample_tidy)

    mrg_l1hs_by_sample_meth_strictly_intergenic <- l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic %>% left_join(fll1hs_gene_expression_by_sample_tidy_strictly_intergenic)
    mrg_l1hs_by_sample_meth_strictly_intergenic_ASP <- l1hs_promoterregion_meth_by_sample_tidy_strictly_intergenic_sense %>% left_join(fll1hs_gene_expression_by_sample_tidy_ASP_strictly_intergenic)

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
            left_join(srna_sample_table)



        model_mean_meth_vs_RIN <- lm(mean_meth ~ RIN, data = tdf) %>%
            summary() %>%
            broom::tidy()
        model_srna_expression_vs_RIN <- lm(srna_expression ~ RIN, data = tdf) %>%
            summary() %>%
            broom::tidy()

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
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_all_rsq.pdf"), w = 4, h = 4, sf = mf)

    # now strictly intergenic
    model_stats_list <- list()
    for (sthresh in mrg_l1hs_by_sample_meth_strictly_intergenic$region %>% unique()) {
        tdf <- mrg_l1hs_by_sample_meth_strictly_intergenic %>%
            filter(region == sthresh) %>%
            left_join(srna_sample_table)


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
        mysaveandstore(str_glue("{outputdir}/L1HS/RIN_by_sample_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4)

        model_stats_list[[sthresh]] <- model_stats
    }

    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))

    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_strictly_intergenic_all_rsq.pdf"), w = 4, h = 4, sf = mf)

    # now strictly intergenic ASP
    model_stats_list <- list()
    for (sthresh in mrg_l1hs_by_sample_meth_strictly_intergenic_ASP$region %>% unique()) {
        tdf <- mrg_l1hs_by_sample_meth_strictly_intergenic_ASP %>%
            filter(region == sthresh) %>%
            left_join(srna_sample_table)


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
        mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_sample_ASP_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = statsframe, sfm = modelframe)

        p <- tdf %>% ggplot(aes(x = RIN, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/RIN_by_sample_ASP_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4)

        model_stats_list[[sthresh]] <- model_stats
    }

    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        anchorbar +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_ASP_strictly_intergenic_all_rsq.pdf"), w = 4, h = 4, sf = mf)


    # now strictly intergenic without ctrl1 and 2
    model_stats_list <- list()
    for (sthresh in mrg_l1hs_by_sample_meth_strictly_intergenic$region %>% unique()) {
        tdf <- mrg_l1hs_by_sample_meth_strictly_intergenic %>%
            filter(region == sthresh) %>%
            left_join(srna_sample_table) %>%
            filter(!(sample_name %in% c("CTRL1", "CTRL2")))
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
        mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_sample_strictly_intergenic_{sthresh}_noctrl12.pdf"), w = 4, h = 4, sf = statsframe, sfm = modelframe)

        p <- tdf %>% ggplot(aes(x = RIN, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/RIN_by_sample_strictly_intergenic_{sthresh}_noctrl12.pdf"), w = 4, h = 4)

        model_stats_list[[sthresh]] <- model_stats
    }

    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))

    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_sample_strictly_intergenic_all_rsq_noctrl12.pdf"), w = 3, h = 4, sf = mf)

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
        p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_{sthresh}_norot.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

        model_stats_list[[sthresh]] <- model_stats
    }

    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(subset = .y)) %>%
        separate(subset, sep = "_", into = c("region", "threshold"))
    p <- mf %>% ggplot() +
        geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
        facet_wrap(~region, nrow = 1) +
        mtclosed +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))

    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_all_rsq.pdf"), w = 6, h = 4, sf = mf)


    p <- mf %>%
        filter(region != "400to600") %>%
        ggplot() +
        geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
        facet_wrap(~region, nrow = 1) +
        mtclosed +
        scale_fill_distiller(palette = "Blues", direction = 1, name = "p-value") + # Use distiller for continuous
        scale_y_continuous(expand = expansion(mult = c(0, .075)))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_noasp_rsq.pdf"), w = 6, h = 4, sf = mf)
    p <- mf %>%
        filter(region != "400to600") %>%
        ggplot() +
        geom_col(aes(x = threshold, y = r.squared, fill = p.value), color = "black") +
        facet_wrap(~region, nrow = 1) +
        mtclosed +
        scale_fill_distiller(palette = "Blues", direction = 1, name = "p-value") + # Use distiller for continuous
        scale_y_continuous(expand = expansion(mult = c(0, .075)), limits = c(0, 1))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_sample_strictly_intergenic_noasp_rsq_fullrange.pdf"), w = 6, h = 4, sf = mf)


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
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        anchorbar
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_gene_all_rsq.pdf"), w = 4, h = 4, sf = mf)

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
        p <- tdf %>% ggplot(aes(x = log(mean_meth + 1), y = log(srna_expression + 1))) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

        p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = loc_superlowres_integrative_stranded_incl_pred2), size = 3, alpha = 0.3) +
            mtclosed +
            scale_palette +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_{sthresh}_colorregion.pdf"), w = 8, h = 4, sf = stats, sfm = model_stats)

        p <- tdf %>% ggplot(aes(x = mean_meth, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = loc_superlowres_integrative_stranded_incl_pred2), size = 3, alpha = 0.3) +
            mtclosed +
            scale_palette +
            lims(x = c(50, 100), y = c(0, 10)) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_{sthresh}_colorregion_zoom.pdf"), w = 8, h = 4, sf = stats, sfm = model_stats)

        model_stats_list[[sthresh]] <- model_stats
    }

    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        anchorbar
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/meth_by_gene_strictly_intergenic_all_rsq.pdf"), w = 4, h = 4, sf = mf)


    ## HDRs
    # Strictly Intergenic

    model_stats_list <- list()
    model_stats_list_filtered <- list()
    for (sthresh in mrg_l1hs_by_gene_hdmr_strictly_intergenic$subset_threshold %>% unique()) {
        tdf <- mrg_l1hs_by_gene_hdmr_strictly_intergenic %>%
            filter(subset_threshold == sthresh)
        tdf %>%
            group_by(gene_id) %>%
            summarise(mean_expr = mean(srna_expression)) %$% mean_expr %>%
            quantile(seq(0, 1, 0.05))
        # idea is to remove those intergenic eleemnts whose epression is likely being driven by adjacent genes - huge inflection point for elements above 90th percentile expression for expression values
        expression_90th_percentile <- tdf %>%
            group_by(gene_id) %>%
            summarise(mean_expr = mean(srna_expression)) %$% mean_expr %>%
            quantile(0.95)
        elements_to_remove <- tdf %>%
            filter(srna_expression > expression_90th_percentile) %$% gene_id %>%
            unique()

        tdf_filtered <- tdf %>% filter(!(gene_id %in% elements_to_remove))

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
                    geom_point(aes(color = condition), size = 3, alpha = 0.4) +
                    mtclosed +
                    scale_conditions +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_gene_strictly_intergenic_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

                p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
                    stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
                    ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                        formula = y ~ x, parse = TRUE, label.y.npc = "top"
                    ) +
                    geom_point(aes(color = loc_superlowres_integrative_stranded_incl_pred), size = 3, alpha = 0.4) +
                    mtclosed +
                    scale_palette +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outputdir}/L1HS/hdr_by_gene_strictly_intergenic_{sthresh}_colorbyregion.pdf"), w = 8, h = 4, sf = stats, sfm = model_stats)

                model_stats_list[[sthresh]] <- model_stats
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
                mysaveandstore(str_glue("{outputdir}/L1HS_95thpercentilefiltered/hdr_by_gene_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)

                model_stats_list_filtered[[sthresh]] <- model_stats
            },
            error = function(e) {

            }
        )
    }


    mf <- model_stats_list %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS/hdr_by_gene_strictly_intergenic_all_rsq.pdf"), w = 4, h = 4, sf = mf)

    mf <- model_stats_list_filtered %>%
        imap_dfr(~ .x %>%
            dplyr::select(r.squared, p.value) %>%
            mutate(region = .y))
    p <- mf %>% ggplot() +
        geom_col(aes(x = region, y = r.squared, fill = p.value), color = "black") +
        mtclosed +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, str_glue("{outputdir}/L1HS_95thpercentile_filtered/hdr_by_gene_strictly_intergenic_all_rsq.pdf"), w = 4, h = 4, sf = mf)
    #############



    #     model_stats_list <- list()
    #     for (sthresh in mrg_l1hs_by_gene_hdmr$subset_threshold %>% unique()) {
    #         tdf <- mrg_l1hs_by_gene_hdmr %>% filter(subset_threshold == sthresh)
    #         gene_ids <- tdf %>%
    #             group_by(gene_id) %>%
    #             group_keys() %>%
    #             pull(gene_id)


    #         models <- tdf %>%
    #             group_split(gene_id) %>%
    #             map2_dfr(gene_ids, ~ {
    #                 model <- lm(srna_expression ~ propUnmeth, data = .x)
    #                 glance_df <- broom::glance(model)
    #                 coef_df <- broom::tidy(model) %>%
    #                     filter(term == "propUnmeth") %>%
    #                     dplyr::select(estimate)
    #                 glance_df %>%
    #                     mutate(
    #                         coef_estimate = coef_df$estimate,
    #                         call = deparse(model$call) %>% paste(collapse = " "),
    #                         gene_id = .y
    #                     )
    #             }) %>%
    #             mutate(subset_threshold = sthresh) %>%
    #             separate(subset_threshold, sep = "_", into = c("region", "threshold"))

    #         p <- models %>%
    #             mutate(direction = ifelse(is.na(coef_estimate), "NA", ifelse(coef_estimate > 0, "positive", "negative"))) %>%
    #             mutate(sig = ifelse(is.na(p.value), "NA", ifelse(p.value <= 0.05, "Sig", "NS"))) %>%
    #             mutate(direction = ifelse(sig == "NA", "NA", direction)) %>%
    #             ggplot() +
    #             geom_bar(aes(x = sig, fill = direction), color = "black") +
    #             labs(y = "Count", title = str_glue("FL L1HS Expr ~ HDR {sthresh}")) +
    #             mtclosed +
    #             anchorbar +
    #             scale_palette
    #         mysaveandstore(str_glue("{outputdir}/L1HS/unique_loci/aaa_hdr_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = models)


    #         for (gene in mrg_l1hs_by_gene_hdmr %>%
    #             filter(subset_threshold == sthresh) %$% gene_id %>%
    #             unique()) {
    #             tdf <- mrg_l1hs_by_gene_hdmr %>%
    #                 filter(subset_threshold == sthresh) %>%
    #                 filter(gene_id == gene)
    #             model <- lm(srna_expression ~ propUnmeth, data = tdf)
    #             stats <- model %>%
    #                 summary() %>%
    #                 broom::tidy()
    #             model_stats <- model %>%
    #                 broom::glance() %>%
    #                 mutate(call = model %>% pluck("call") %>% deparse() %>% paste(collapse = " "))
    #             p <- tdf %>% ggplot(aes(x = propUnmeth, y = srna_expression)) +
    #                 stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
    #                 ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
    #                     formula = y ~ x, parse = TRUE, label.y.npc = "top"
    #                 ) +
    #                 geom_point(aes(color = condition), size = 3) +
    #                 mtclosed +
    #                 scale_conditions +
    #                 theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #             mysaveandstore(str_glue("{outputdir}/L1HS/unique_loci/{gene}_hdr_by_sample_{sthresh}.pdf"), w = 4, h = 4, sf = stats, sfm = model_stats)
    #         }
    #     }
    # }
}


required_fraction_of_total_cg <- "0.75"
region <- "L1HS_FL"
stats <- read_csv(sprintf("ldna/results/%s/tables/reads_new/%s_%s/by_gene.csv", params$mod_code, region, required_fraction_of_total_cg))


sigDM <- stats %>%
    filter(subset == "0to909") %>%
    filter(meth_threshold == 0.5) %>%
    filter(term == "conditionAD") %>%
    filter(p.value <= 0.05) %>%
    mutate(direction = case_when(
        statistic > 0 ~ "Hypo",
        statistic < 0 ~ "Hyper",
        TRUE ~ "Other"
    ))


{ # needed plot
    pf <- fll1hs_gene_expression_tidy_strictly_intergenic %>%
        filter(gene_id %in% (sigDM %>% filter(direction == "Hypo") %$% gene_id)) %>%
        group_by(sample_name) %>%
        summarise(sumcounts = sum(srna_expression)) %>%
        left_join(srna_sample_table)
    model <- lm(sumcounts ~ condition, data = pf)
    model_with_RIN <- lm(sumcounts ~ condition + RIN, data = pf)
    model_with_RIN_with_cov <- lm(sumcounts ~ condition + RIN + age + sex, data = pf)
    stats <- model %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "no cov no RIN")
    stats_with_RIN <- model_with_RIN %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "RIN")
    stats_with_RIN_with_cov <- model_with_RIN_with_cov %>%
        summary() %>%
        broom::tidy() %>%
        mutate(type = "RIN with Cov")
    statsframe <- bind_rows(stats, stats_with_RIN, stats_with_RIN_with_cov)
    p <- pf %>%
        ggplot(aes(x = condition, y = sumcounts)) +
        geom_boxplot(aes(fill = condition)) +
        geom_point() +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/L1HSFLexpression_box_DM_hypo.pdf"), w = 4, h = 4, sf = statsframe)
    p <- pf %>%
        ggplot(aes(x = condition, y = sumcounts)) +
        geom_col(data = pf %>% group_by(condition) %>% summarise(sumcounts = mean(sumcounts)), aes(fill = condition), color = "black") +
        geom_point() +
        scale_conditions +
        anchorbar +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/L1HSFLexpression_bar_DM_hypo.pdf"), w = 4, h = 4, sf = statsframe)
}


ea <- read_delim("srna/results/agg/enrichment_analysis/results_table_unbiased.tsv", delim = "\t")



great_hypo <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/ldna/results/m/tables/great/GO:BP/Hypogreat_enrichment.tsv", delim = " ")
great_hyper <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/ldna/results/m/tables/great/GO:BP/Hypergreat_enrichment.tsv", delim = " ")

great_hyper %>% filter(p_adjust <= 0.05)

load("ldna/results/m/tables/great/t05/gs_subcat/et_noextension/tablesMsigdb.rds")
great_t05 <- tablesMsigdb
load("ldna/results/m/tables/great/t01/gs_subcat/et_noextension/tablesMsigdb.rds")
great_t01 <- tablesMsigdb

great_df05 <- map_df(names(great_t05), function(ontology) {
    map_df(names(great_t05[[ontology]]), function(direction) {
        great_t05[[ontology]][[direction]] %>%
            mutate(ontology = ontology, direction = direction)
    })
}) %>%
    mutate(threhold = "t05") %>%
    tibble()

great_df01 <- map_df(names(great_t01), function(ontology) {
    map_df(names(great_t01[[ontology]]), function(direction) {
        great_t01[[ontology]][[direction]] %>%
            mutate(ontology = ontology, direction = direction)
    })
}) %>%
    mutate(threhold = "t01") %>%
    tibble()

# View the result
head(great_df)

directions <- c("Hypo", "Hyper", "Dif")
for (direc in directions) {
    meth_ids <- great_df05 %>%
        filter(direction == direc) %>%
        filter(p_adjust <= 0.05) %$% id
    if (direc == "Hypo") {
        rna_ids <<- ea %>%
            filter(NES > 0) %>%
            filter(p.adjust <= 0.05) %$% ID
    } else if (direc == "Hyper") {
        rna_ids <<- ea %>%
            filter(NES < 0) %>%
            filter(p.adjust <= 0.05) %$% ID
    } else {
        rna_ids <<- ea %>% filter(p.adjust <= 0.05) %$% ID
    }

    rna_ids_all <- ea %$% ID

    insetn <- meth_ids %in% rna_ids %>% sum()
    total_n <- meth_ids %in% rna_ids_all %>% sum()

    DMinDEfrac <- insetn / total_n
    DEfrac <- length(rna_ids) / length(rna_ids_all)
    print(str_glue("{direc}: dm in de {DMinDEfrac} versus DE frac {DEfrac}"))

    concordant <- meth_ids[meth_ids %in% rna_ids]
}
