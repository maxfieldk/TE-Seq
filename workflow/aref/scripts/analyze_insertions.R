module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(ggpubr)
library(ggh4x)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = "aref/A.REF_tldr/A.REF.table.txt",
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/default/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/default/ref_pre_ins_filtering.fa",
            contigs_to_keep = "aref/default/contigs_to_keep.txt",
            filtered_tldr = "aref/default/A.REF.table.kept_in_updated_ref.txt"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/default/A.REF_Analysis/tldr_plots/tldr_plots.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)
df_filtered <- read_delim(inputs$filtered_tldr) %>% dplyr::rename(nonref_UUID = UUID)

nrdf <- rmann %>%
    filter(refstatus == "NonRef") %>%
    left_join(df_filtered, by = "nonref_UUID") %>%
    mutate(homozygosity = factor(ifelse(fraction_reads_count >= 0.95, "homozygous", "heterozygous"), levels = c("homozygous", "heterozygous"))) %>%
    mutate(known_nonref = factor(ifelse(is.na(NonRef), "novel", "known"), levels = c("novel", "known")))


p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative") +
    facet_grid2(rows = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))

mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, known_nonref) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "known_nonref") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_nofacet_known.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, homozygosity) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "homozygosity") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_nofacet_zygosity.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_nofacet.pdf", outputdir), 6, 4)


novel_frac_df <- nrdf %>%
    group_by(rte_subfamily, known_nonref, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(known_nonref == "novel") %>%
    mutate(frac_novel = frac_total, frac_known = 1 - frac_total)

homozyg_frac_df <- nrdf %>%
    group_by(rte_subfamily, homozygosity, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(homozygosity == "heterozygous") %>%
    mutate(frac_heterozygous = frac_total, frac_homozygous = 1 - frac_total)


hp <- homozyg_frac_df %>%
    ggbarplot(x = "rte_subfamily", y = "frac_heterozygous") +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Subfamily")
np <- novel_frac_df %>%
    ggbarplot(x = "rte_subfamily", y = "frac_novel") +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Subfamily")
ptch <- p1 / np / hp + plot_layout(heights = c(1, 0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_nofacet1.pdf", outputdir), 6, 6)

for (rte in homozyg_frac_df %$% rte_subfamily) {
    pf <- homozyg_frac_df %>%
        filter(rte_subfamily == rte) %>%
        pivot_longer(cols = c(frac_heterozygous, frac_homozygous), names_to = "label")

    p <- ggpie(pf, "value", "label", fill = "label") +
        scale_fill_manual(values = c("yellow", "blue"))
    mysaveandstore(sprintf("%s/%s_zygosity_pie.pdf", outputdir, rte), 3, 3)
}



p2 <- nrdf %>%
    group_by(loc_integrative) %>%
    summarise(count = n()) %>%
    mutate(loc_integrative = fct_reorder(loc_integrative, count)) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "loc_integrative", y = "counts", fill = "loc_integrative") +
    ggtitle("Non-reference RTE Insertions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Genomic Context", y = "Count") +
    scale_palette +
    coord_flip() +
    mtclosed + anchorbar
mysaveandstore(pl = p2, sprintf("%s/insertions_genomic_context.pdf", outputdir), 6, 4)
library(patchwork)
p <- p1 + p2 + plot_layout(widths = c(1, 1), guides = "collect")
mysaveandstore(sprintf("%s/insertions_sub_fig.pdf", outputdir), 12, 5)


p <- nrdf %>%
    filter(rte_subfamily == "L1HS") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mtopen +
    scale_palette +
    anchorbar
mysaveandstore(sprintf("%s/l1hs_length_in_updated_ref.pdf", outputdir), 5, 5)


x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)
