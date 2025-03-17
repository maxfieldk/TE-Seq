module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
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
# library(ggrastr)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "ref" = "aref/default/A.REF.fa",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "annotation_genes" = conf$annotation_genes
        ), env = globalenv())
        assign("outputs", list(
            environment = "srna/results/agg/tpm_sources/telescope_multi/characterize_tpm_environment.RData"
        ), env = globalenv())
        assign("inputs", list(
            tpm = "srna/outs/agg/tpm/telescope_multi/tpmdf.tsv"
        ), env = globalenv())
    }
)

counttype <- params[["counttype"]]
counttype_label <- gsub("telescope_", "", counttype) %>%
    gsub("counttype_", "", .) %>%
    str_to_title()
outputdir <- dirname(outputs$environment)

rmfragments <- read_csv(params$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(params$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

tpmdf <- read_delim(inputs$tpm)
repeats <- tpmdf %>% filter(gene_id %in% c(rmann %$% gene_id))
genes <- tpmdf %>% filter(!(gene_id %in% c(rmann %$% gene_id)))
trepeats <- repeats %>%
    pivot_longer(cols = -gene_id) %>%
    dplyr::rename(sample = name, tpm = value) %>%
    mutate(gene_or_te = "repeat") %>%
    left_join(rmann)
wrepeats <- repeats %>% left_join(rmann)
tgenes <- genes %>%
    pivot_longer(cols = -gene_id) %>%
    dplyr::rename(sample = name, tpm = value) %>%
    mutate(gene_or_te = "gene")
wgenes <- genes


resultsdf <- wrepeats
tidydf <- trepeats
library(ggpubr)



# pan contrast
pancontrastbarplots <- function(tdf = tidydf, ontology_column = "rte_subfamily", ontology_column_value = "L1HS", ontology_column_value_modifier = "", facetvars = c("req_integrative", "genic_loc"), refstatus_to_include = c("Ref", "NonRef")) {
    # Generated many variants of a simple grouped barplot with and without various statistics
    facetvarsstring <- paste(facetvars, collapse = "_")
    refstatusstring <- paste(refstatus_to_include, collapse = "_")
    nconditions <- length(conf$levels)
    nhorizontalfacets <- tdf[[facetvars[2]]] %>%
        unique() %>%
        length()
    width <- 4 * 1 / 3 * nconditions * 1 / 2 * nhorizontalfacets
    height <- 8
    # Apply filters and transformations
    df <- tdf %>%
        filter(!!sym(ontology_column) == ontology_column_value) %>%
        mutate(refstatus = as.character(refstatus)) # Convert refstatus to character if it's a factor

    ontology_column_value <- paste0(ontology_column_value, ontology_column_value_modifier)
    # Perform filtering
    pf <- df %>%
        filter(refstatus %in% refstatus_to_include) %>%
        group_by(sample, condition, across(all_of(facetvars))) %>%
        summarise(sample_sum = sum(tpm), condition = dplyr::first(condition), n = n()) %>%
        ungroup() %>%
        filter(if_all(all_of(facetvars), ~ . != "Other")) %>%
        arrange(across(all_of(facetvars)))

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Sum TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_sum/%s/%s_bar_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Sum TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, ref.group = conf$levels[1], bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_sum/%s/%s_bar_stats_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)
    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Sum TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_sum/%s/%s_bar_stats_allcomps_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Sum TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = FALSE, ref.group = conf$levels[1], bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_sum/%s/%s_bar_stats_allsigannot_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)
    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Sum TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_sum/%s/%s_bar_stats_allsigannot_allcomps_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)

    # now mean across elements
    pf <- df %>%
        filter(refstatus %in% refstatus_to_include) %>%
        group_by(sample, condition, across(all_of(facetvars))) %>%
        summarise(sample_mean = mean(tpm), condition = dplyr::first(condition), n = n()) %>%
        ungroup() %>%
        filter(if_all(all_of(facetvars), ~ . != "Other")) %>%
        arrange(across(all_of(facetvars)))

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_mean", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Mean TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_mean/%s/%s_bar_%s_%s4.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_mean", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Mean TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, ref.group = conf$levels[1], bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_mean/%s/%s_bar_stats_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)
    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_mean", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Mean TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_mean/%s/%s_bar_stats_allcomps_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)

    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_mean", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Mean TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = FALSE, ref.group = conf$levels[1], bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_mean/%s/%s_bar_stats_allsigannot_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)
    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_mean", fill = "condition", facet.by = facetvars, add = c("mean_se", "dotplot"), scales = "free_y") +
        geom_text(aes(x = -Inf, y = Inf, label = paste0("N = ", n)),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        labs(x = "", y = "Mean TPM", subtitle = counttype_label, title = ontology_column_value) +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, bracket.nudge.y = -0.1, step.increase = .1)
    mysaveandstore(pl = p, fn = sprintf("%s/%s/pan_contrast/bar_mean/%s/%s_bar_stats_allsigannot_allcomps_%s_%s.pdf", outputdir, counttype, ontology_column_value, ontology_column_value, facetvarsstring, refstatusstring), width, height)
}

tryCatch(
    {
        tidy_ASP <- tidydf %>%
            filter(grepl("__AS$", gene_id)) %>%
            dplyr::select(sample, tpm, gene_id) %>%
            dplyr::rename(ASP_tpm = tpm) %>%
            mutate(gene_id = case_when(
                TRUE ~ gsub("__AS$", "", gene_id)
            )) %>%
            left_join(tidydf)


        p <- tidy_ASP %>%
            filter(rte_subfamily == "L1HS") %>%
            ggplot(aes(x = tpm, y = ASP_tpm)) +
            geom_point() +
            facet_wrap(~loc_lowres_integrative_stranded) +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HS_SP_ASP_cor.pdf"), w = 10, h = 6)

        p <- tidy_ASP %>%
            filter(rte_subfamily == "L1HS") %>%
            ggplot(aes(x = log(tpm + 1), y = log(ASP_tpm + 1))) +
            geom_point() +
            facet_wrap(~loc_lowres_integrative_stranded) +
            mtclosed
        mysaveandstore(str_glue("{outputdir}/L1HS_SP_ASP_log_cor.pdf"), w = 10, h = 6)


        tASP <- tidy_ASP %>%
            dplyr::select(-tpm) %>%
            dplyr::rename(tpm = ASP_tpm) %>%
            left_join(sample_table %>% dplyr::rename(sample = sample_name))
        pancontrastbarplots(tdf = tASP, ontology_column = "rte_subfamily", ontology_column_value = "L1HS", ontology_column_value_modifier = "_ASP", facetvars = c("req_integrative", "loc_lowres_integrative_stranded"), refstatus_to_include = c("Ref", "NonRef"))
    },
    error = function(e) {

    }
)

# pan contrast
for (g_var in c("rte_family", "rte_subfamily")) {
    groups <- tidydf[[g_var]] %>%
        unique() %>%
        na.omit()
    groups <- groups[groups != "Other"]
    for (group in groups) {
        width <- 8
        height <- 8
        df <- tidydf %>% filter(!!sym(g_var) == group)

        pf <- df %>%
            group_by(sample, req_integrative, genic_loc) %>%
            summarise(sample_sum = sum(tpm)) %>%
            ungroup()
        nsamples <- pf %$% sample %>%
            unique() %>%
            length()

        p <- pf %>%
            mutate(req_integrative = factor(req_integrative, levels = c("Old Trnc", "Old FL", "Yng Trnc", "Yng FL", "Yng Intact"))) %>%
            mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "sample", facet.by = c("req_integrative", "genic_loc"), scales = "free_y") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_samples +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar.pdf", outputdir, group), 1 + nsamples / 2.5, height)


        p <- pf %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc"), scales = "free_y") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_palette +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative.pdf", outputdir, group), 1 + nsamples / 2.5, 5)

        p <- pf %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc")) +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_palette +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative_same_scale.pdf", outputdir, group), 1 + nsamples / 2.5, 5)

        pff <- df %>%
            group_by(sample, req_integrative) %>%
            summarise(sample_sum = sum(tpm)) %>%
            ungroup()
        p <- pff %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_palette +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative_nofacet.pdf", outputdir, group), 1 + nsamples / 2.5, 5)


        pff <- df %>%
            group_by(sample, req_integrative) %>%
            summarise(sample_sum = sum(tpm)) %>%
            ungroup() %>%
            filter(grepl("total", sample))
        p <- pff %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_palette +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative_nofacet_fl.pdf", outputdir, group), 1 + 15 / 2.5, 5)


        if (df %$% rte_superfamily %>% unique() == "LTR") {
            pf <- df %>%
                group_by(sample, ltr_viral_status, genic_loc) %>%
                summarise(sample_sum = sum(tpm)) %>%
                ungroup() %>%
                filter(ltr_viral_status != "Other")
            p <- pf %>%
                arrange(ltr_viral_status) %>%
                mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
                ggbarplot(x = "sample", y = "sample_sum", fill = "sample", facet.by = c("ltr_viral_status", "genic_loc"), scales = "free_y") +
                labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
                mtclosedgridh +
                scale_samples +
                scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            mysaveandstore(sprintf("%s/single_group/%s_bar_ltr_viral_status_.pdf", outputdir, group), width, height + 4)
        }
    }
}

joined <- full_join(trepeats, tgenes)

p <- joined %>%
    group_by(sample, gene_or_te) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "gene_or_te") +
    scale_fill_brewer("paired") +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/multiple_groups/gene_or_repeat.pdf", outputdir), w = 8, h = 5)

pf <- joined %>% mutate(Class = ifelse(gene_or_te == "gene", "Gene", repeat_superfamily))
p <- pf %>%
    mutate(Class = factor(Class, levels = c("Gene", "Other", "LowComp", "SimpleRep", "Retroposon", "SAT", "DNA", "LTR", "LINE", "SINE"))) %>%
    group_by(sample, Class) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "Class") +
    xlab("") +
    ylab("TPM") +
    scale_palette +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/multiple_groups/gene_or_repeat_type.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)

p_rep_only <- pf %>%
    mutate(Class = factor(Class, levels = c("Gene", "LowComp", "SimpleRep", "Retroposon", "SAT", "DNA", "LTR", "LINE", "SINE", "Other"))) %>%
    group_by(sample, Class) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    filter(Class != "Gene") %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "Class") +
    xlab("") +
    ylab("TPM") +
    scale_palette +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(pl = p_rep_only, sprintf("%s/multiple_groups/repeat_only_type.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)


refseq <- import(params$annotation_genes)
refseqdf <- refseq %>%
    as.data.frame() %>%
    tibble()
biotypes <- refseqdf %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_biotype)

biotype_levels <- c("repeat", "Other", "lncRNA", "miRNA", "snRNA", "tRNA", "rRNA", "protein_coding")
pf <- joined %>%
    left_join(biotypes) %>%
    mutate(Class = ifelse(gene_or_te == "repeat", "repeat", gene_biotype))

pff <- pf %>% mutate(Class = case_when(
    Class %in% biotype_levels ~ Class,
    TRUE ~ "Other"
))
pff %$% Class %>% table()
p <- pff %>%
    mutate(Class = factor(Class, levels = biotype_levels)) %>%
    group_by(sample, Class) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "Class") +
    xlab("") +
    ylab("TPM") +
    scale_palette +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/multiple_groups/gene_biotype_or_repeat.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)
p_gene_only <- pff %>%
    mutate(Class = factor(Class, levels = biotype_levels)) %>%
    group_by(sample, Class) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    filter(Class != "repeat") %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "Class") +
    xlab("") +
    ylab("TPM") +
    scale_palette +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(pl = p_gene_only, sprintf("%s/multiple_groups/gene_biotype_only.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)


ptch <- p_gene_only + p_rep_only + plot_layout(nrow = 2)
mysaveandstore(pl = ptch, sprintf("%s/multiple_groups/gene_repeat_split.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 8)

x <- tibble(OUT = "")
write_tsv(x, file = outputs$environment)
