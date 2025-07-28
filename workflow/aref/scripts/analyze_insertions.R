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


rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE
)

df_filtered <- read_delim(inputs$filtered_tldr) %>% dplyr::rename(nonref_UUID = UUID)

nrdf <- rmann %>%
    filter(refstatus == "NonRef") %>%
    filter(!grepl("__AS$", gene_id)) %>%
    left_join(df_filtered, by = "nonref_UUID") %>%
    mutate(zygosity = factor(ifelse(fraction_reads_count >= 0.95, "homozygous", "heterozygous"), levels = c("homozygous", "heterozygous"))) %>%
    mutate(known_nonref = factor(
        case_when(
            conf$update_ref_with_tldr$known_nonref$response == "no" ~ "NotCalled",
            is.na(NonRef) ~ "novel",
            TRUE ~ "known"
        ),
        levels = c("novel", "known", "NotCalled")
    ))

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
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative") +
    facet_grid2(rows = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_tall.pdf", outputdir), 6, 5)

p1 <- nrdf %>%
    filter(rte_subfamily != "Other") %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative", label = TRUE) +
    facet_grid2(rows = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_tall_noother.pdf", outputdir), 6, 5)

p1 <- nrdf %>%
    filter(rte_subfamily != "Other") %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative", label = TRUE) +
    facet_grid2(cols = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_tall_noother_wide.pdf", outputdir), 7, 4)


nrdf %>%
    filter(rte_subfamily == "Other") %>%
    pw()

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
    group_by(rte_family, rte_subfamily, zygosity) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "zygosity") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_nofacet_zygosity.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative", label = TRUE) +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily_nofacet.pdf", outputdir), 6, 4)


novel_frac_df <- nrdf %>%
    filter(rte_subfamily != "Other") %>%
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
    filter(rte_subfamily != "Other") %>%
    group_by(rte_subfamily, zygosity, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(zygosity == "heterozygous") %>%
    mutate(frac_heterozygous = frac_total, frac_homozygous = 1 - frac_total)

hp <- homozyg_frac_df %>%
    pivot_longer(cols = c(frac_heterozygous, frac_homozygous)) %>%
    ggbarplot(x = "rte_subfamily", y = "value", fill = "name") +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("tan", "lightblue")) +
    labs(x = "Subfamily")
np <- novel_frac_df %>%
    pivot_longer(cols = c(frac_novel, frac_known)) %>%
    ggbarplot(x = "rte_subfamily", y = "value", fill = "name") +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("brown", "grey")) +
    labs(x = "Subfamily")
ptch <- p1 / np / hp + plot_layout(heights = c(0.4, 0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch.pdf", outputdir), 6, 6)

ptch <- p1 / hp + plot_layout(heights = c(0.4, 0.2), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch_2.pdf", outputdir), 6, 5)

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
mysaveandstore(pl = p2, sprintf("%s/insertions_genomic_context.pdf", outputdir), 4.5, 3.5)

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
mysaveandstore(sprintf("%s/l1hs_length_in_updated_ref.pdf", outputdir), 3, 4)


x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)



tryCatch(
    {
        rmann %$% masked_in_ref %>% unique()
        tt <- rmann %>%
            filter(refstatus != "NonCentral") %>%
            filter(rte_subfamily != "Other") %>%
            filter(asp != "asp") %>%
            filter(rte_family %in% c("L1", "Alu", "SVA"))
        p <- tt %>%
            filter(masked_in_ref == "yes") %>%
            filter(fraction_masked == 1) %>%
            group_by(rte_subfamily) %>%
            summarise(n = n()) %>%
            ggplot() +
            geom_bar(aes(x = rte_subfamily, y = n), color = "black", fill = "lightgrey", stat = "identity") +
            geom_text(aes(x = rte_subfamily, y = n, label = n),
                position = position_stack(vjust = 1.1), size = 3
            ) +
            labs(title = "NonReference Deletions", x = "Subfamily", y = "Count") +
            coord_flip() +
            mtclosed +
            anchorbar +
            scale_palette
        mysaveandstore(sprintf("%s/homozygous_deletions1111.pdf", outputdir), w = 3, h = 4.35)

        p <- rmann %>%
            filter(rte_subfamily == "L1HS") %>%
            filter(refstatus != "NonCentral") %>%
            mutate(length = ifelse(str_detect(rte_length_req, "Trnc"), "Trnc", "FL")) %>%
            mutate(req_integrative = case_when(
                fraction_masked == 1 ~ "Masked",
                rte_length_req == "Trnc" ~ "Trnc",
                req_integrative == "Yng FL" ~ "FL",
                req_integrative == "Yng Intact" ~ "Intact",
                TRUE ~ req_integrative
            )) %>%
            group_by(refstatus, req_integrative, length) %>%
            summarise(n = n()) %>%
            mutate(req_integrative = factor(req_integrative, levels = c("Masked", "Trnc", "FL", "Intact"))) %>%
            ungroup() %>%
            ggplot() +
            geom_bar(aes(x = refstatus, y = n, fill = req_integrative), color = "black", stat = "identity") +
            geom_text(aes(x = refstatus, y = n, label = n, group = req_integrative),
                position = position_stack(vjust = 0.5), size = 3
            ) +
            labs(title = "FL L1HS", x = "Subfamily", y = "Count", fill = "") +
            facet_wrap(~length, scales = "free_y") +
            mtclosed +
            anchorbar +
            scale_palette +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/L1HS_FL_masked.pdf", outputdir), w = 4, h = 5)
    },
    error = function(e) {

    }
)
