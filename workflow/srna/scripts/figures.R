module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")

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
            "inputdir" = "srna/results/agg/deseq",
            "outputdir" = "srna/results/agg/repeatanalysis",
            "counttypes" = c("telescope_multi"),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
        repeatanalysis_plots = "srna/results/agg/repeatanalysis/telescope_multi/repeatanalysisplots_plots.RData",
        enrichment_analysis_repeats_plots = sprintf("srna/results/agg/enrichment_analysis_repeats/%s/enrichment_analysis_repeats_plots.RData", params$counttypes),
        deseq_plots = sprintf("srna/results/agg/deseq/%s/deseq_plots.RData", params$counttypes),
        enrichment_analysis_plots = "srna/results/agg/enrichment_analysis/enrichment_analysis_plots.RData"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "srna/results/agg/figures/outfile.txt"
        ), env = globalenv())
    }
)




#Repeats
load(inputs$repeatanalysis_plots)
repeatanalysis_plots <- mysaveandstoreplots
names(repeatanalysis_plots)

p <- repeatanalysis_plots[["srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_RS_0w_low_O2_vs_PRO_low_O2/pvp/L1HS_intactness_req_genic_loc_log2labels.pdf"]]
mysave()

p1<-repeatanalysis_plots[["srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_LSEN_vs_PRO/stripp/AluY_rte_length_req_ALL_stats.pdf"]]
p2<-repeatanalysis_plots[["srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_LSEN_vs_PRO/stripp/HERVK_INT_rte_length_req_ALL.pdf"]]
p3<-repeatanalysis_plots[["srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_LSEN_vs_PRO/stripp/L1HS_rte_length_req_ALL_stats.pdf"]]
mysaveandstore(pl = p1)
print(p1)
ptch <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect")


p1<-"srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_LSEN_vs_PRO/pvp/L1HS_rte_length_req_genic_loc_labels.pdf"
p1<-"srna/results/agg/repeatanalysis/telescope_multi/telescope_multi/condition_LSEN_vs_PRO/myheatmap/L1HS_ALL_genic_loc_DE_scaled.pdf"


#EA


load(inputs$enrichment_analysis_plots)
enrichment_analysis_plots <- mysaveandstoreplots
names(enrichment_analysis_plots)

"srna/results/agg/enrichment_analysis/heatmap_Senescence.pdf"

"srna/results/agg/enrichment_analysis/condition_LSEN_vs_PRO/gsea/msigdbH/nes15.pdf"
"srna/results/agg/enrichment_analysis/condition_ESEN_vs_PRO/gsea/msigdbH/nes15.pdf"
"srna/results/agg/enrichment_analysis/condition_LSEN_vs_ESEN/gsea/msigdbH/nes15.pdf"
"srna/results/agg/enrichment_analysis/condition_LSEN_vs_ESEN/gsea/msigdbH/all_genes/scatter_HALLMARK_INFLAMMATORY_RESPONSE.pdf"
x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)