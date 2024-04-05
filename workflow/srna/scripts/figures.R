source("workflow/scripts/defaults.R")
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
library(plotly)
library(DT)
library(ggExtra)
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)


conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]


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
            "inputdir" = "srna/results/agg/deseq_telescope",
            "outputdir" = "srna/results/agg/repeatanalysis_telescope",
            "tecounttypes" = c("telescope_multi"),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
        repeatanalysis_plots <- "srna/results/agg/repeatanalysis_telescope/repeatanalysisplots_plots.RData",
        enrichment_analysis_repeats_plots = "srna/results/agg/enrichment_analysis_repeats/{tecounttype}/enrichment_analysis_repeats_plots.RData",
        deseq_plots = "srna/results/agg/deseq_telescope/{tecounttype}/deseq_plots.RData",
        enrichment_analysis_plots = "srna/results/agg/enrichment_analysis/enrichment_analysis_plots.RData"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "srna/results/agg/figures/outfile.txt"
        ), env = globalenv())
    }
)










x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)