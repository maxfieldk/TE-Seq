source("workflow/scripts/defaults.R")
module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")


library("tibble")
library("tidyr")
library(DT)
library(ggpubr)
library(GenomicRanges)
library(paletteer)
library(ggplot2)
library(scales)


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
        repeatanalysis_plots = "srna/results/agg/repeatanalysis_telescope/repeatanalysisplots_plots.RData",
        enrichment_analysis_repeats_plots = sprintf("srna/results/agg/enrichment_analysis_repeats/%s/enrichment_analysis_repeats_plots.RData", params$tecounttypes),
        deseq_plots = sprintf("srna/results/agg/deseq_telescope/%s/deseq_plots.RData", params$tecounttypes),
        enrichment_analysis_plots = "srna/results/agg/enrichment_analysis/enrichment_analysis_plots.RData"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "srna/results/agg/figures/outfile.txt"
        ), env = globalenv())
    }
)

load(inputs$repeatanalysis_plots)
repeatanalysis_plots <- mysaveandstoreplots
names(repeatanalysis_plots)

p <- repeatanalysis_plots[["srna/results/agg/repeatanalysis_telescope/telescope_multi/condition_RS_0w_low_O2_vs_PRO_low_O2/pvp/L1HS_l1_intactness_req_genic_loc_log2labels.png"]]
mysave()





x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)