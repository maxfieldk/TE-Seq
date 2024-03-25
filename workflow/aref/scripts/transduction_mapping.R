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
            filtered_tldr = "aref/tldr/tldr.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/annotations/repeatmasker_annotation.csv",
            ref = "aref/ref.fa",
            blast_njs = "aref/ref.njs"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/RefAnalysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/RefAnalysis/transduction_df.csv"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

df_filtered <- read_delim("aref/tldr/tldr.table.kept_in_updated_ref.txt")

trsd <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd %>% head(n = 4) %$% Transduction_3p

trsd %$% TransductionLen %>% summary()

save(mysaveandstoreplots, file = outputs$plots)
