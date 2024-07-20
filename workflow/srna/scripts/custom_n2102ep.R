module_name <- "srna"
# module_name <- snakemake@params$module_name
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
set.seed(123)
# whether or not to store plots in list for figure generation at script end
store_var <- "yes"

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggVennDiagram")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)
library(rtracklayer)
library(ComplexUpset)
library(patchwork)
library(scales)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        print(module_name)
        if (module_name == "srna") {
            assign("params", list(
                "inputdir" = "srna/results/agg/deseq",
                "outputdir" = "srna/results/agg/repeatanalysis",
                "counttype" = "telescope_multi"
            ), env = globalenv())
            assign("inputs", list(
                "resultsdf" = "srna/results/agg/deseq/resultsdf.tsv",
                "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
                "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
            ), env = globalenv())
            assign("outputs", list(
                "environment" = "srna/results/agg/repeatanalysis/telescope_multi/repeatanalysisplots_environment.RData"
            ), env = globalenv())
        } else if (module_name == "lrna") {
            assign("params", list(
                "inputdir" = "lrna/results/agg/deseq",
                "outputdir" = "lrna/results/agg/repeatanalysis/relaxed",
                "counttype" = "relaxed"
            ), env = globalenv())
            assign("inputs", list(
                "resultsdf" = "lrna/results/agg/deseq/resultsdf.tsv",
                "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
                "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
            ), env = globalenv())
            assign("outputs", list(
                "environment" = "lrna/results/agg/repeatanalysis/relaxed/repeatanalysisplots_environment.RData"
            ), env = globalenv())
        } else {
            print("varible assignment failed")
        }
    }
)

outputdir <- params$outputdir
contrasts <- conf$contrasts
sample_table <- read.csv(conf$sample_table)
counttype <- params$counttype

## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t") %>% filter(counttype == !!counttype)
r_annotation_fragmentsjoined <- read_csv(inputs$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(inputs$r_repeatmasker_annotation) %>%
    mutate(req_integrative = factor(req_integrative, levels = c("Old Trnc", "Old FL", "Yng Trnc", "Yng FL", "Yng Intact"))) %>%
    mutate(ltr_viral_status = factor(ltr_viral_status, levels = c("Int (Has 5LTR)", "Int (No 5'LTR)", "5'LTR (FL Int)", "3'LTR (FL Int)", "5'LTR (Trnc Int)", "3'LTR (Trnc Int)", "LTR (Solo)", "Other")))
resultsdfwithgenes <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)
resultsdf <- resultsdfwithgenes %>% filter(gene_or_te != "gene")

gres <- resultsdfwithgenes %>% filter(gene_or_te == "gene")
colnames(gres)
gres %>%
    arrange(padj_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    print(n = 100)

gres %>%
    arrange(-log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    filter(padj_condition_N21_ORF1p_500_vs_N21_rIgG_500 < 0.05) %>%
    write_delim("top_genes.tsv", delim = "\t")

gres %>%
    filter(grepl("^RN", gene_id)) %>%
    print(n = 200)

gres %>%
    filter(gene_id == "LOC105373777") %>%
    print(width = Inf)


gene_ids <- gres %>%
    arrange(-log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    filter(log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500 > 0) %>%
    filter(log2FoldChange_condition_N21_ORF1p_300_vs_N21_rIgG_300 > 0) %>%
    filter(N21_ORF1p_500_b1_l1 > 25) %$% gene_id %>%
    head(20)

i <- 1
for (id in gene_ids) {
    pf <- gres %>%
        filter(gene_id == id) %>%
        dplyr::select(sample_table$sample_name) %>%
        pivot_longer(cols = everything()) %>%
        mutate(name_copy = name) %>%
        separate_wider_delim(cols = name_copy, delim = "_", names = c("cell_type", "fraction", "stringency", "batch", "lane"))
    p <- pf %>% ggplot(aes(x = name, y = value, fill = fraction)) +
        geom_col() +
        facet_wrap(~cell_type) +
        labs(title = id) +
        coord_flip() +
        mtclosed
    mysaveandstore(sprintf("srna/results/special_analysis/log2fc_large_n2102/%s_%s_orf1p_associated_gene.pdf", i, id), 7, 10)
    i <- i + 1
}

gene_ids <- gres %>%
    filter(grepl("^LOC", gene_id)) %>%
    arrange(-log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    filter(log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500 > 0) %>%
    filter(log2FoldChange_condition_N21_ORF1p_300_vs_N21_rIgG_300 > 0) %>%
    filter(N21_ORF1p_500_b1_l1 > 25) %$% gene_id %>%
    head(20)

i <- 1
for (id in gene_ids) {
    pf <- gres %>%
        filter(gene_id == id) %>%
        dplyr::select(sample_table$sample_name) %>%
        pivot_longer(cols = everything()) %>%
        mutate(name_copy = name) %>%
        separate_wider_delim(cols = name_copy, delim = "_", names = c("cell_type", "fraction", "stringency", "batch", "lane"))
    p <- pf %>% ggplot(aes(x = name, y = value, fill = fraction)) +
        geom_col() +
        facet_wrap(~cell_type) +
        labs(title = id) +
        coord_flip() +
        mtclosed
    mysaveandstore(sprintf("srna/results/special_analysis/log2fc_large_n2102_LOC/%s_%s_orf1p_associated_gene.pdf", i, id), 7, 10)
    i <- i + 1
}

gene_ids <- gres %>%
    arrange(-log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    filter(log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500 > 0) %>%
    filter(log2FoldChange_condition_N21_ORF1p_300_vs_N21_rIgG_300 > 0) %>%
    filter(padj_condition_N21_ORF1p_500_vs_N21_rIgG_500 < 0.05) %>%
    filter(padj_condition_N21_ORF1p_300_vs_N21_rIgG_300 < 0.05) %>%
    filter(padj_condition_MT6_ORF1p_500_vs_MT6_rIgG_500 < 0.05) %>%
    filter(log2FoldChange_condition_MT6_ORF1p_500_vs_MT6_rIgG_500 > 0) %>%
    filter(N21_ORF1p_500_b1_l1 > 25) %$% gene_id %>%
    head(20)

i <- 1
for (id in gene_ids) {
    pf <- gres %>%
        filter(gene_id == id) %>%
        dplyr::select(sample_table$sample_name) %>%
        pivot_longer(cols = everything()) %>%
        mutate(name_copy = name) %>%
        separate_wider_delim(cols = name_copy, delim = "_", names = c("cell_type", "fraction", "stringency", "batch", "lane"))
    p <- pf %>% ggplot(aes(x = name, y = value, fill = fraction)) +
        geom_col() +
        facet_wrap(~cell_type) +
        labs(title = id) +
        coord_flip() +
        mtclosed
    mysaveandstore(sprintf("srna/results/special_analysis/all_critieria/%s_%s_orf1p_associated_gene.pdf", i, id), 7, 10)
    i <- i + 1
}

gene_ids <- gres %>%
    filter(grepl("^LOC", gene_id)) %>%
    arrange(-log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500) %>%
    filter(log2FoldChange_condition_N21_ORF1p_500_vs_N21_rIgG_500 > 0) %>%
    filter(log2FoldChange_condition_N21_ORF1p_300_vs_N21_rIgG_300 > 0) %>%
    filter(padj_condition_N21_ORF1p_500_vs_N21_rIgG_500 < 0.05) %>%
    filter(padj_condition_N21_ORF1p_300_vs_N21_rIgG_300 < 0.05) %>%
    filter(padj_condition_MT6_ORF1p_500_vs_MT6_rIgG_500 < 0.05) %>%
    filter(log2FoldChange_condition_MT6_ORF1p_500_vs_MT6_rIgG_500 > 0) %>%
    filter(N21_ORF1p_500_b1_l1 > 25) %$% gene_id %>%
    head(20)

i <- 1
for (id in gene_ids) {
    pf <- gres %>%
        filter(gene_id == id) %>%
        dplyr::select(sample_table$sample_name) %>%
        pivot_longer(cols = everything()) %>%
        mutate(name_copy = name) %>%
        separate_wider_delim(cols = name_copy, delim = "_", names = c("cell_type", "fraction", "stringency", "batch", "lane"))
    p <- pf %>% ggplot(aes(x = name, y = value, fill = fraction)) +
        geom_col() +
        facet_wrap(~cell_type) +
        labs(title = id) +
        coord_flip() +
        mtclosed
    mysaveandstore(sprintf("srna/results/special_analysis/all_critieria_LOC/%s_%s_orf1p_associated_gene.pdf", i, id), 7, 10)
    i <- i + 1
}
