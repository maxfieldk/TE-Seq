source("~/data/common/myDefaults.r")
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

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")[["aref"]]
)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = "aref/tldr/tldr.table.txt",
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/annotations/repeatmasker_annotation.csv",
            ref = "aref/ref_pre_ins_filtering.fa"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/RefAnalysis/tldr_plots/tldr_plots.rds",
            contigs_to_keep = "aref/contigs_to_keep.txt",
            filtered_tldr = "aref/tldr/tldr.table.kept_in_updated_ref.txt"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


# rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
# rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
# rmann <- left_join(rmfragments, rmfamilies)

df <- read_delim(inputs$tldroutput) %>%
    mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS")

df_filtered <- read_delim("aref/tldr/tldr.table.kept_in_updated_ref.txt")


df %$% Subfamily %>% table()
df %$% Filter %>% table()
p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysaveandstore(sprintf("%s/insertions.png", outputdir), 5, 5)


p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_wrap(~Family, scales = "free") +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysaveandstore(sprintf("%s/insertions_subfamily.png", outputdir), 5, 5)

p <- df %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mythemeconditions
mysaveandstore(sprintf("%s/l1hs_length.png", outputdir), 5, 5)


# only insertions that are in the updated reference
df_filtered %$% Subfamily %>% table()
df_filtered %$% Filter %>% table()
p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysaveandstore(sprintf("%s/insertions_in_updated_ref.png", outputdir), 5, 5)


p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_wrap(~Family, scales = "free") +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysaveandstore(sprintf("%s/insertions_subfamily_in_updated_ref.png", outputdir), 6, 4)

p <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mythemeconditions
mysaveandstore(sprintf("%s/l1hs_length_in_updated_ref.png", outputdir), 5, 5)

###############

trsd <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd %>% head(n = 4) %$% Transduction_3p

trsd %$% TransductionLen %>% summary()

save(mysaveandstoreplots, file = outputs$plots)
