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
library(GGally)
library(ggpubr)




conf <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")

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
            "txdb" = conf[["srna"]]$txdb
        ), env = globalenv())
        assign("inputs", list(
            "srna_results" = "srna/results/agg/deseq/resultsdf.tsv",
            "lrna_results" = "lrna/results/agg/deseq/resultsdf.tsv",
            "ldna_methylation" = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", conf[["ldna"]][["samples"]], conf[["ldna"]][["samples"]]),
            "rteprommeth" = "ldna/Rintermediates/perelementdf_promoters.tsv",
            "dmrs" = "ldna/results/tables/dmrs.CG_m.tsv",
            "dmls" = "ldna/results/tables/dmrs.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "plots" = "integrated/epigenome_transcriptome_correlation/objects/plots.rda"
        ), env = globalenv())
    }
)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
RMdf <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
RM <- GRanges(RMdf)

modules_run <- conf$pipelines_to_deploy


if (conf$ldna$single_condition == "no") {
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

    ldna_sample_table <- read_csv("conf/sample_table_ldna.csv")
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


### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "_.*family")]
    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        big_ontology_groups <- c(big_ontology_groups, RMdf %>%
            pull(!!sym(ontology)) %>%
            unique())
    }
    big_ontology_groups <- big_ontology_groups %>% unique()

    modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
    region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
    element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
}



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




outputdir <- dirname(outputs$plots)
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)



promoter_methdf <- read_delim("ldna/Rintermediates/refseq_gene_promoter_methylation.tsv", col_names = TRUE)
pergenedf <- promoter_methdf %>%
    group_by(gene_id, sample) %>%
    summarise(meanmeth = mean(pctM)) %>%
    pivot_wider(names_from = sample, values_from = meanmeth) %>%
    ungroup()
colnames(pergenedf) <- paste0("ldna_", colnames(pergenedf)) %>% gsub("ldna_gene_id", "gene_id", .)

pergenedf_tidy <- pivot_longer(pergenedf, cols = -gene_id, names_to = "sample_name", values_to = "pctM") %>%
    mutate(sample_name = gsub("ldna_", "", sample_name)) %>%
    filter(sample_name %in% ldna_samples)

if (pergenedf_tidy %$% sample_name %>% unique() %>% length() == 1) {
    pergenedf_tidy <- pergenedf_tidy %>% dplyr::select(-sample_name)
}
rteprommeth %>% filter(rte_subfamily == "L1HS") %$% rte_length_req
rteprommeth <- read_delim(inputs$rteprommeth)
perrepeatdf <- rteprommeth %>%
    group_by(gene_id, condition) %>%
    summarise(meanmeth = mean(mean_meth)) %>%
    pivot_wider(names_from = condition, values_from = meanmeth) %>%
    ungroup()
colnames(perrepeatdf) <- paste0("ldna_", colnames(perrepeatdf)) %>% gsub("ldna_gene_id", "gene_id", .)

perrepeatdf_tidy <- pivot_longer(perrepeatdf, cols = -gene_id, names_to = "sample_name", values_to = "pctM") %>%
    mutate(sample_name = gsub("ldna_", "", sample_name)) %>%
    filter(sample_name %in% ldna_samples)
if (perrepeatdf_tidy %$% sample_name %>% unique() %>% length() == 1) {
    perrepeatdf_tidy <- perrepeatdf_tidy %>% dplyr::select(-sample_name)
}
methdf_tidy <- bind_rows(pergenedf_tidy %>% mutate(gene_or_te = "gene"), perrepeatdf_tidy %>% mutate(gene_or_te = "repeat"))

if ("srna" %in% modules_run) {
    srna_gene_expression_tidy <- srna_df %>%
        dplyr::select(gene_id, paste0("srna_", srna_samples)) %>%
        pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "srna_expression") %>%
        mutate(sample_name = gsub("srna_", "", sample_name))
    # if (conf$ldna$single_condition == "yes") {
    #     srna_gene_expression_tidy <- srna_gene_expression_tidy %>%
    #         left_join(integrated_sample_table %>% dplyr::select(ldna, srna) %>% dplyr::rename(sample_name = srna)) %>%
    #         mutate(sample_name = ldna) %>%
    #         dplyr::select(-ldna)
    # }

    srnapctMdf <- full_join(methdf_tidy, srna_gene_expression_tidy) %>%
        mutate(srna_expressionlog10 = log10(srna_expression + 1)) %>%
        drop_na()


    if (conf$ldna$single_condition == "yes") {
        srnapctMdfmean <- full_join(methdf_tidy, srna_gene_expression_tidy) %>%
            group_by(gene_id, gene_or_te) %>%
            summarise(pctM = mean(pctM), srna_expression = mean(srna_expression)) %>%
            mutate(srna_expressionlog10 = log10(srna_expression + 1)) %>%
            ungroup() %>%
            drop_na()
    }

    p <- ggscatter(srnapctMdf %>% filter(gene_or_te == "gene"),
        x = "pctM", y = "srna_expressionlog10", color = "sample_name", add = "reg.line", alpha = 0.25,
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
    ) +
        mtclosed
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes.pdf", outputdir), w = 6, h = 6)
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes.pdf", outputdir), raster = TRUE, w = 6, h = 6)

    p <- ggscatter(srnapctMdf %>% filter(gene_or_te == "gene"),
        x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
    ) +
        mtclosed
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes_no_fill.pdf", outputdir), w = 6, h = 6)
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes_no_fill.pdf", outputdir), raster = TRUE, w = 6, h = 6)


    p <- ggscatter(srnapctMdfmean %>% filter(gene_or_te == "gene"),
        x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
    ) + mtclosed
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes_mean.pdf", outputdir), w = 6, h = 6)
    mysaveandstore(sprintf("%s/srna_vs_pctM_genes_mean.pdf", outputdir), raster = TRUE, w = 6, h = 6)

    testrnapctMdf <- srnapctMdf %>%
        filter(gene_or_te == "repeat") %>%
        left_join(r_repeatmasker_annotation)
    testrnapctMdfmean <- srnapctMdfmean %>%
        filter(gene_or_te == "repeat") %>%
        left_join(r_repeatmasker_annotation)
    for (repeat_family in testrnapctMdf %$% rte_subfamily %>% unique()) {
        p <- ggscatter(testrnapctMdf %>% filter(rte_subfamily == repeat_family),
            x = "pctM", y = "srna_expressionlog10", color = "sample_name", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(testrnapctMdf %>% filter(rte_subfamily == repeat_family),
            x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(testrnapctMdfmean %>% filter(rte_subfamily == repeat_family),
            x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)
    
        if (repeat_family == "L1HS") {
            p <- ggscatter(testrnapctMdfmean %>% filter(rte_subfamily == repeat_family),
                x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25, color = "intactness_req",
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
            ) +
                mtclosed
            mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), w = 4, h = 4)
            mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), raster = TRUE, w = 4, h = 4)
        
        }
    }
}


if ("lrna" %in% modules_run) {
    lrna_gene_expression_tidy <- lrna_df %>%
        select(gene_id, paste0("lrna_", lrna_samples)) %>%
        pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "lrna_expression") %>%
        mutate(sample_name = gsub("lrna_", "", sample_name))
    if (conf$ldna$single_condition == "yes") {
        lrna_gene_expression_tidy <- lrna_gene_expression_tidy %>%
            left_join(integrated_sample_table %>% dplyr::select(ldna, lrna) %>% dplyr::rename(sample_name = lrna)) %>%
            mutate(sample_name = ldna) %>%
            dplyr::select(-ldna)
    }
    lrnapctMdf <- full_join(methdf_tidy, lrna_gene_expression_tidy) %>%
        mutate(lrna_expressionlog10 = log10(lrna_expression + 1)) %>%
        drop_na()


    p <- ggscatter(lrnapctMdf %>% filter(gene_or_te == "gene"),
        x = "pctM", y = "lrna_expressionlog10", add = "reg.line", alpha = 0.25,
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
    ) +
        mtclosed
    mysaveandstore(sprintf("%s/lrna_vs_pctM_genes.pdf", outputdir), w = 6, h = 6)
    mysaveandstore(sprintf("%s/lrna_vs_pctM_genes.pdf", outputdir), raster = TRUE, w = 6, h = 6)

    testrnapctMdf <- lrnapctMdf %>%
        filter(gene_or_te == "repeat") %>%
        left_join(r_repeatmasker_annotation)
    for (repeat_family in testrnapctMdf %$% rte_subfamily %>% unique()) {
        p <- ggscatter(testrnapctMdf %>% filter(rte_subfamily == repeat_family),
            x = "pctM", y = "lrna_expressionlog10", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/lrna_vs_pctM_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/lrna_vs_pctM_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)
    }
}


dir.create(dirname(outputs$plots), showWarnings = FALSE, recursive = TRUE)
save(plots, file = outputs$plots)
