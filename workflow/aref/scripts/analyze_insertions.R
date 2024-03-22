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
            tldroutput = "aref/updated_ref/tldr.table.txt"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/RefAnalysis/tldr_plots/tldr_plots.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
plots <- list()

dff <- read_delim(inputs$tldroutput)
df <- dff %>% filter(Filter == "PASS")

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
mysave(sprintf("%s/insertions.png", outputdir), 5, 5)
plots[["l1hs_length"]] <- p


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
mysave(sprintf("%s/insertions_subfamily.png", outputdir), 5, 5)
plots[["l1hs_length"]] <- p

p <- df %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 5850, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mythemeconditions
mysave(sprintf("%s/l1hs_length.png", outputdir), 5, 5)
plots[["l1hs_length"]] <- p



trsd <- df %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p))

trsd %$% TransductionLen %>% summary()

save(plots, file = outputs$plots)
