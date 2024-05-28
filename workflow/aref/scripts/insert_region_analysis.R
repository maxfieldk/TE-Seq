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
library(msigdbr)
library(regioneR)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            filtered_tldr = "aref/A.REF_tldr/A.REF.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/A.REF.fa",
            blast_njs = "aref/A.REF.njs"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/A.REF_Analysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/A.REF_Analysis/transduction_df.csv"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "aref/A.REF"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "aref/A.REF"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)
refseq <- import(conf$ref_refseq_gtf)
temptxs <- refseq[!is.na(mcols(refseq)$tag)]
refseq_select <- temptxs[mcols(temptxs)$tag == "RefSeq Select"]


gene_sets = msigdbr(species = "human")
sets_of_interest <- c("REACTOME_ONCOGENIC_MAPK_SIGNALING", "REACTOME_ONCOGENE_INDUCED_SENESCENCE", "KEGG_PATHWAYS_IN_CANCER", "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_P53_SIGNALING_PATHWAY", "KEGG_CITRATE_CYCLE_TCA_CYCLE")
all.genes <- refseq_select[refseq_select$type == "transcript"]

rmann %>% filter(rte_subfamily == "L1HS") %$% req_integrative %>% unique()
intact_l1hs <- rmann %>% filter(rte_subfamily == "L1HS") %>% filter(grepl("Intact", req_integrative)) %>% GRanges()
nonintact_l1hs_fl <- rmann %>% filter(rte_subfamily == "L1HS") %>% filter(grepl("Full Length", req_integrative)) %>% GRanges()
l1hs_fl <- rmann %>% filter(rte_subfamily == "L1HS") %>% filter(grepl(">", rte_length_req)) %>% GRanges()
(intact_l1hs %$% gene_id) %in% (nonintact_l1hs_fl %$% gene_id)

#I'll run two test per set of interest
# 1. Are intact L1hs more likely to be adjacent to members of this set vs any gene
# 2. Are intact L1hs more likely to be adjacent to members of this set vs nonintact L1hs
sets_of_interest <- gene_sets %>% filter(gs_cat == "H") %$% gs_name %>% unique()
library(BiocParallel)
register(MulticoreParam(12))


{
test1_pval <- list()
test1_zscore <- list()
test2_pval <- list()
test2_zscore <- list()
for (set_name in sets_of_interest) {
soi <- gene_sets %>% filter(gs_name %in% set_name) %$% gene_symbol
soigr <- refseq_select[refseq_select$gene_id %in% soi & refseq_select$type == "transcript"]
#add 5000 extra bp on both sides of teh regions
soigr <- resize(soigr, width=width(soigr)+10000, fix="center")

pt1 <- permTest(A=soigr, ntimes=20, randomize.function=resampleRegions, universe=all.genes,
evaluate.function=numOverlaps, B=intact_l1hs, verbose=TRUE)

pt2 <- permTest(A=intact_l1hs, ntimes=20, randomize.function=resampleRegions, universe=nonintact_l1hs_fl,
evaluate.function=numOverlaps, B=soigr, verbose=TRUE)

test1_pval[[set_name]] <- pt1$numOverlaps$pval
test1_zscore[[set_name]] <- pt1$numOverlaps$zscore
test2_pval[[set_name]] <- pt2$numOverlaps$pval
test2_zscore[[set_name]] <- pt2$numOverlaps$zscore
}

test1_padj <- p.adjust(unlist(test1_pval), method = "fdr")
test2_padj <- p.adjust(unlist(test2_pval), method = "fdr")

resdf <- tibble(
    set_name = names(test1_pval),
    test1_pval = unlist(test1_pval),
    test1_padj = test1_padj,
    test1_zscore = unlist(test1_zscore),
    test2_pval = unlist(test2_pval),
    test2_padj = test2_padj,
    test2_zscore = unlist(test2_zscore)
)
resdf %>% arrange(-test1_zscore)
hits <- resdf %>% filter(test1_zscore > 3) %$% set_name

hit_test1_pval <- list()
hit_test1_zscore <- list()
hit_test2_pval <- list()
hit_test2_zscore <- list()
for (set_name in hits) {
soi <- gene_sets %>% filter(gs_name %in% set_name) %$% gene_symbol
soigr <- refseq_select[refseq_select$gene_id %in% soi & refseq_select$type == "transcript"]
#add 5000 extra bp on both sides of teh regions
soigr <- resize(soigr, width=width(soigr)+10000, fix="center")

pt1 <- permTest(A=soigr, ntimes=50, randomize.function=resampleRegions, universe=all.genes,
evaluate.function=numOverlaps, B=intact_l1hs, verbose=TRUE)

pt2 <- permTest(A=intact_l1hs, ntimes=50, randomize.function=resampleRegions, universe=nonintact_l1hs_fl,
evaluate.function=numOverlaps, B=soigr, verbose=TRUE)

hit_test1_pval[[set_name]] <- pt1$numOverlaps$pval
hit_test1_zscore[[set_name]] <- pt1$numOverlaps$zscore
hit_test2_pval[[set_name]] <- pt2$numOverlaps$pval
hit_test2_zscore[[set_name]] <- pt2$numOverlaps$zscore
} 

hit_test1_padj <- p.adjust(unlist(hit_test1_pval), method = "fdr")
hit_test2_padj <- p.adjust(unlist(hit_test2_pval), method = "fdr")

hit_resdf <- tibble(
    set_name = names(hit_test1_pval),
    test1_pval = unlist(hit_test1_pval),
    test1_padj = hit_test1_padj,
    test1_zscore = unlist(hit_test1_zscore),
    test2_pval = unlist(hit_test2_pval),
    test2_padj = hit_test2_padj,
    test2_zscore = unlist(hit_test2_zscore)
)
}



pt <- permTest(A=my.regions, B=repeats, randomize.function=randomizeRegions,
evaluate.function=numOverlaps)