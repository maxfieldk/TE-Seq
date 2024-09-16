module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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




p <- tgenes %>%
    group_by(sample) %>%
    arrange(-tpm) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(sample == conf$samples[[1]]) %>%
    ggplot(aes(x = rank, y = log10(tpm))) +
    geom_point()
mysaveandstore()

resultsdf <- wrepeats
tidydf <- trepeats
library(ggpubr)
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
        p <- pf %>%
            mutate(req_integrative = factor(req_integrative, levels = c("Old Trnc", "Old FL", "Yng Trnc", "Yng FL", "Yng Intact"))) %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "sample", facet.by = c("req_integrative", "genic_loc"), scales = "free_y") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_fill_viridis_d(option = "viridis") +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar.pdf", outputdir, group), width, height)

        p <- pf %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc"), scales = "free_y") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_fill_viridis_d() +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative.pdf", outputdir, group), 10, 5)

        p <- pf %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc")) +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_fill_viridis_d() +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative_same_scale.pdf", outputdir, group), 10, 5)

        pff <- df %>%
            group_by(sample, req_integrative) %>%
            summarise(sample_sum = sum(tpm)) %>%
            ungroup()
        p <- pff %>%
            arrange(req_integrative) %>%
            ggbarplot(x = "sample", y = "sample_sum", fill = "req_integrative") +
            labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
            mtclosedgridh +
            scale_fill_viridis_d() +
            scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, .075))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/single_group/%s_bar_reqintegrative_nofacet.pdf", outputdir, group), 6, 5)

        if (df %$% rte_superfamily %>% unique() == "LTR") {
            pf <- df %>%
                group_by(sample, ltr_viral_status, genic_loc) %>%
                summarise(sample_sum = sum(tpm)) %>%
                ungroup() %>%
                filter(ltr_viral_status != "Other")
            p <- pf %>%
                arrange(ltr_viral_status) %>%
                ggbarplot(x = "sample", y = "sample_sum", fill = "sample", facet.by = c("ltr_viral_status", "genic_loc"), scales = "free_y") +
                labs(x = "", y = "TPM", subtitle = counttype_label, title = group) +
                mtclosedgridh +
                scale_fill_viridis_d() +
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
    mutate(Class = factor(Class, levels = c("Gene", "Other", "Retroposon", "SAT", "DNA", "LTR", "LINE", "SINE"))) %>%
    group_by(sample, Class) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample", y = "tpm", fill = "Class") +
    xlab("") +
    ylab("TPM") +
    scale_fill_viridis_d() +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/multiple_groups/gene_or_repeat_type.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)

refseq <- import(params$annotation_genes)
refseqdf <- refseq %>%
    as.data.frame() %>%
    tibble()
biotypes <- refseqdf %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_biotype)

biotype_levels <- c("repeat", "other", "lncRNA", "miRNA", "snRNA", "tRNA", "rRNA", "protein_coding")
pf <- joined %>%
    left_join(biotypes) %>%
    mutate(Class = ifelse(gene_or_te == "repeat", "repeat", gene_biotype))

pff <- pf %>% mutate(Class = case_when(
    Class %in% biotype_levels ~ Class,
    TRUE ~ "other"
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
    scale_fill_viridis_d() +
    mtclosed +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(sprintf("%s/multiple_groups/gene_biotype_or_repeat.pdf", outputdir), w = 1 + 1.5 * length(conf$samples) / 2.4, h = 4)



x <- tibble(OUT = "")
write_tsv(x, file = outputs$environment)
