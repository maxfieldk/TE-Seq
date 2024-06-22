module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

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
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/ref_pre_ins_filtering.fa",
            contigs_to_keep = "aref/contigs_to_keep.txt",
            filtered_tldr = "aref/A.REF_tldr/A.REF.table.kept_in_updated_ref.txt"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/A.REF_Analysis/tldr_plots/tldr_plots.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

nrdf <- rmann %>% filter(refstatus == "NonRef")
nrdf %$% family %>% table()

p1 <- nrdf %>% group_by(rte_family, rte_subfamily, req_integrative) %>% summarise(count = n()) %>% 
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative") +
    facet_grid2(rows = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))

mysaveandstore(pl = p1, sprintf("%s/insertions_subfamily.pdf", outputdir), 6, 4)



p2 <- nrdf %>% group_by(loc_integrative) %>% summarise(count = n()) %>% 
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


df <- read_delim(inputs$tldroutput) %>%
    mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS")

df_filtered <- read_delim(inputs$filtered_tldr)


df %$% Subfamily %>% table()
df %$% Filter %>% table()
p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    mtopen + anchorbar
mysaveandstore(sprintf("%s/insertions.pdf", outputdir), 5, 5)


p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_grid(~Family, scales = "free", space = "free_x") +
    geom_bar(color = "black") +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    mtclosed + anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/insertions_subfamily.pdf", outputdir), 5, 5)

p <- df %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mtopen +
    anchorbar
mysaveandstore(sprintf("%s/l1hs_length.pdf", outputdir), 5, 5)


# only insertions that are in the updated reference
df_filtered %$% Subfamily %>% table()
df_filtered %$% Filter %>% table()
p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    mtopen + anchorbar

mysaveandstore(sprintf("%s/insertions_in_updated_ref.pdf", outputdir), 3, 4)


p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_grid(~Family, scales = "free", space = "free_x") +
    geom_bar(color = "black") +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    mtclosed + anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/insertions_subfamily_in_updated_ref.pdf", outputdir), 6, 4)

p <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mtopen + paletteer::scale_fill_paletteer_d(conf$default_palette) +
    anchorbar
mysaveandstore(sprintf("%s/l1hs_length_in_updated_ref.pdf", outputdir), 5, 5)

###############

trsd <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd %>% head(n = 4) %$% Transduction_3p

trsd %$% TransductionLen %>% summary()

save(mysaveandstoreplots, file = outputs$plots)
