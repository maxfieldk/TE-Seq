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
            "rteprommeth" = "ldna/Rintermediates/m/perelementdf_promoters.tsv",
            "dmrs" = "ldna/results/m/tables/dmrs.CG_m.tsv",
            "dmls" = "ldna/results/m/tables/dmrs.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "plots" = "integrated/epigenome_transcriptome_correlation/objects/plots.rda"
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
    filter(refstatus != "NonCentral")

RM <- GRanges(RMdf)

modules_run <- confALL$pipelines_to_deploy
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

ldna_sample_table <- read_csv("conf/sample_table_ldna.csv")

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



outputdir <- dirname(outputs$plots)
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)



promoter_methdf <- read_delim("ldna/Rintermediates/m/refseq_gene_promoter_methylation.tsv", col_names = TRUE)
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


rteprommeth <- read_delim(inputs$rteprommeth)
perrepeatdf_tidy <- rteprommeth %>% dplyr::select(gene_id, sample, mean_meth)
perrepeatdfbycondition <- rteprommeth %>%
    group_by(gene_id, condition) %>%
    summarise(meanmeth = mean(mean_meth)) %>%
    pivot_wider(names_from = condition, values_from = meanmeth) %>%
    ungroup()
colnames(perrepeatdf) <- paste0("ldna_", colnames(perrepeatdf)) %>% gsub("ldna_gene_id", "gene_id", .)



perl1hs_5utr_region <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE) %>% mutate(region = ordered(region, levels = c("328", "500", "909")))
# perl1hs_5utr_region$sample <- factor(perl1hs_5utr_region$sample, levels = conf$samples)
# perl1hs_5utr_region$condition <- factor(perl1hs_5utr_region$condition, levels = conf$levels)

# perrepeatdf_tidy <- pivot_longer(perrepeatdf, cols = -gene_id, names_to = "sample_name", values_to = "pctM") %>%
#     mutate(sample_name = gsub("ldna_", "", sample_name)) %>%
#     filter(sample_name %in% ldna_samples)
# if (perrepeatdf_tidy %$% sample_name %>% unique() %>% length() == 1) {
#     perrepeatdf_tidy <- perrepeatdf_tidy %>% dplyr::select(-sample_name)
# }

methdf_tidy <- bind_rows(pergenedf_tidy %>% mutate(gene_or_te = "gene"), perrepeatdf_tidy %>% mutate(gene_or_te = "repeat") %>% dplyr::rename(pctM = mean_meth, sample_name = sample))

if ("srna" %in% modules_run) {
    srna_gene_expression_tidy <- srna_df %>%
        dplyr::select(gene_id, paste0("srna_", srna_samples)) %>%
        pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "srna_expression") %>%
        mutate(sample_name = gsub("srna_", "", sample_name)) %>%
        mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        ))

    srnapctMdf <- full_join(methdf_tidy, srna_gene_expression_tidy) %>%
        mutate(srna_expressionlog10 = log10(srna_expression + 1)) %>%
        drop_na()

    testrnapctMdf <- srnapctMdf %>%
        filter(gene_or_te == "repeat") %>%
        mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
        inner_join(RMdf %>% mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
            dplyr::select(-sample_name))

    srnapctMdf_regions <- left_join(perl1hs_5utr_region %>% dplyr::rename(sample_name = sample), srna_gene_expression_tidy) %>%
        mutate(srna_expressionlog10 = log10(srna_expression + 1)) %>%
        drop_na()
    testrnapctMdf_regions <- srnapctMdf_regions %>%
        mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
        inner_join(RMdf %>% filter(rte_subfamily == "L1HS") %>% mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
            dplyr::select(-sample_name), by = c("gene_id"))


    if (confALL$ldna$single_condition == "yes") {
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
        mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
        inner_join(RMdf %>% mutate(gene_id = case_when(
            grepl("NI_", gene_id) ~ paste0(sample_name, "_", gene_id),
            TRUE ~ gene_id
        )) %>%
            dplyr::select(-sample_name))
    # testrnapctMdfmean <- srnapctMdfmean %>%
    #     filter(gene_or_te == "repeat") %>%
    #     left_join(r_repeatmasker_annotation)
    region =  "L1HS_intactness_req_ALL"
    by_gene_id <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region))

    by_sample <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region)) %>% dplyr::rename(region1 = region)

    for (subset_threshold in by_sample %$% subset_threshold %>% unique()) {
        tdf <- srna_gene_expression_tidy %>% filter(rte_subfamily == "")
        
    }


    for (regionstring in testrnapctMdf_regions %$% region %>% unique()) {
        tdf <- testrnapctMdf_regions %>% filter(region == regionstring)

        group_level <- tdf %>%
            group_by(sample_name, condition, rte_subfamily) %>%
            summarise(av_pctM = mean(mean_meth), av_srna = mean(srna_expression)) %>%
            left_join(by_sample)

        group_level <- group_level %>%
            ungroup() %>%
            dplyr::select(-condition) %>%
            left_join(sample_table)
        stats <- lm(av_pctM ~ av_srna, data = group_level) %>%
            summary() %>%
            broom::tidy()
        p <- group_level %>% ggplot(aes(x = av_pctM, y = av_srna)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions
        mysaveandstore(sprintf("%s/group_level_srna_vs_pctM_L1HS_%s.pdf", outputdir, regionstring), w = 4, h = 4, sf = stats)
    }






    for (repeat_family in testrnapctMdf %$% rte_subfamily %>% unique()) {
        tdf <- testrnapctMdf %>% filter(rte_subfamily == repeat_family)

        group_level <- tdf %>%
            group_by(sample_name, condition, rte_subfamily) %>%
            summarise(av_pctM = mean(pctM), av_srna = mean(srna_expression))

        group_level <- group_level %>%
            ungroup() %>%
            dplyr::select(-condition) %>%
            left_join(sample_table)
        stats <- lm(av_pctM ~ av_srna, data = group_level) %>%
            summary() %>%
            broom::tidy()
        p <- group_level %>% ggplot(aes(x = av_pctM, y = av_srna)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions
        mysaveandstore(sprintf("%s/group_level_srna_vs_pctM_%s.pdf", outputdir, repeat_family), w = 4, h = 4, sf = stats)
    }


    conditions <- conf$levels
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

    top_l1hs_movers <- rteprommeth %>%
        filter(rte_subfamily == "L1HS") %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        group_by(rte_subfamily) %>%
        mutate(rank_change = row_number()) %>%
        mutate(top_change = ifelse(rank_change <= 10, "Top", "NotTop")) %>%
        ungroup() %>%
        filter(rte_subfamily == "L1HS") %$% gene_id %>%
        head(n = 10)

    for (element in top_l1hs_movers) {
        tdf <- testrnapctMdf %>% filter(gene_id == element)
        # TODO need to fix the condition color

        stats <- lm(pctM ~ srna_expression, data = tdf) %>%
            summary() %>%
            broom::tidy()
        p <- tdf %>% ggplot(aes(x = pctM, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions
        mysaveandstore(sprintf("%s/unique_loci/srna_vs_pctM_%s.pdf", outputdir, element), w = 4, h = 4, sf = stats)
    }

    for (element in rteprommeth %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(t05 == "Hypo") %$% gene_id %>%
        unique()) {
        tdf <- testrnapctMdf %>% filter(gene_id == element)

        stats <- lm(pctM ~ srna_expression, data = tdf) %>%
            summary() %>%
            broom::tidy()
        p <- tdf %>% ggplot(aes(x = pctM, y = srna_expression)) +
            stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
            ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
                formula = y ~ x, parse = TRUE, label.y.npc = "top"
            ) +
            geom_point(aes(color = condition), size = 3) +
            mtclosed +
            scale_conditions
        mysaveandstore(sprintf("%s/unique_loci/srna_vs_pctM_%s.pdf", outputdir, element), w = 4, h = 4, sf = stats)
    }




    for (repeat_family in testrnapctMdf %$% rte_subfamily %>% unique()) {
        tdf <- testrnapctMdf %>% filter(rte_subfamily == repeat_family)

        group_level <- tdf %>%
            group_by(sample_name, condition, rte_subfamily) %>%
            summarise(av_pctM = mean(pctM), av_srna = mean(srna_expression))

        group_level <- group_level %>%
            ungroup() %>%
            dplyr::select(-condition) %>%
            left_join(sample_table)
        stats <- lm(av_pctM ~ av_srna, data = group_level) %>%
            summary() %>%
            broom::tidy()
        p <- group_level %>% ggplot(aes(x = av_pctM, y = av_srna, color = condition)) +
            geom_point(size = 3) +
            stat_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) + # Least squares line
            mtclosed +
            scale_conditions
        p <- ggscatter(group_level,
            x = "av_pctM", y = "av_srna", color = "condition", add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed +
            scale_conditions
        mysaveandstore(sprintf("%s/group_level_srna_vs_pctM1_%s.pdf", outputdir, repeat_family), w = 6, h = 6)


        stats <- lm(srna_expression ~ pctM, data = tdf) %>%
            summary() %>%
            broom::tidy()
        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expressionlog10", color = "sample_name", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expression", color = "sample_name", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_linear_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_linear_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)


        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expression", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_linear_no_fill.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_linear_no_fill.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        tdf <-
            p <- ggscatter(tdf,
                x = "pctM", y = "srna_expressionlog10", color = "sample_name", add = "reg.line", alpha = 0.25,
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
            ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expression", color = "sample_name", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_linear_%s.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_linear_%s.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)


        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_no_fill.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        p <- ggscatter(tdf,
            x = "pctM", y = "srna_expression", add = "reg.line", alpha = 0.25,
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        ) +
            mtclosed
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_linear_no_fill.pdf", outputdir, repeat_family), w = 6, h = 6)
        mysaveandstore(sprintf("%s/srna_vs_pctM_%s_linear_no_fill.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)


        # p <- ggscatter(testrnapctMdfmean %>% filter(rte_subfamily == repeat_family),
        #     x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25,
        #     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        #     conf.int = TRUE, # Add confidence interval
        #     cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        #     cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        # ) +
        #     mtclosed
        # mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), w = 6, h = 6)
        # mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), raster = TRUE, w = 6, h = 6)

        # if (repeat_family == "L1HS") {
        #     p <- ggscatter(testrnapctMdfmean %>% filter(rte_subfamily == repeat_family),
        #         x = "pctM", y = "srna_expressionlog10", add = "reg.line", alpha = 0.25, color = "intactness_req",
        #         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        #         conf.int = TRUE, # Add confidence interval
        #         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        #         cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n", color = "red")
        #     ) +
        #         mtclosed
        #     mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), w = 4, h = 4)
        #     mysaveandstore(sprintf("%s/srna_vs_pctM_%s_mean.pdf", outputdir, repeat_family), raster = TRUE, w = 4, h = 4)
        # }
    }
}


if ("lrna" %in% modules_run) {
    lrna_gene_expression_tidy <- lrna_df %>%
        select(gene_id, paste0("lrna_", lrna_samples)) %>%
        pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "lrna_expression") %>%
        mutate(sample_name = gsub("lrna_", "", sample_name))
    if (confALL$ldna$single_condition == "yes") {
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
