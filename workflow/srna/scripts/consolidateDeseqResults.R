library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")

conf <- configr::read.config(file = "conf/config.yaml")["srna"]

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "inputdir" = "results/agg/deseq_telescope",
            "outputdir" = "results/agg/repeatanalysis_telescope"
        ), env = globalenv())
        assign("outputs", list(
            "resultsdf" = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
        ), env = globalenv())
    }
)



outputdir <- params$outputdir
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
# order matters for the colors!
contrast_colors <- conf$contrast_colors
condition_colors <- conf$condition_colors
contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap
tecounttypes <- conf$tecounttypes
tecounttypes <- conf$tecounttypes
lengthreq <- conf$lengthreq
maincontrast <- contrasts[1]



contrast_colors <- unname(unlist(contrast_colors))
condition_colors <- unname(unlist(condition_colors))
peptable <- read.csv(conf$peptable)




##########################



# build dds frames
RESLIST <- list()
for (tecounttype in tecounttypes) {
    contrastl <- list()
    for (contrast in contrasts) {
        ddsres_rtes <- read_csv(paste(params$inputdir, tecounttype, contrast, "results_rtes.csv", sep = "/"))
        ddsres_rtes$gene_or_te <- "repeat"
        ddsres_genes <- read_csv(paste(params$inputdir, tecounttype, contrast, "results_genes.csv", sep = "/"))
        ddsres_genes$gene_or_te <- "gene"
        ddsres <- bind_rows(ddsres_rtes, ddsres_genes)
        ddsresmod <- ddsres %>%
            dplyr::rename(gene_id = ...1) %>%
            mutate(Significance = ifelse(padj < 0.05, ifelse(padj < 0.001, "< 0.001", "< 0.05"), "> 0.05")) %>%
            dplyr::select(c(gene_id, log2FoldChange, stat, padj, Significance, gene_or_te)) %>%
            dplyr::rename(!!paste0("log2FoldChange_", contrast) := log2FoldChange) %>%
            dplyr::rename(!!paste0("stat_", contrast) := stat) %>%
            dplyr::rename(!!paste0("padj_", contrast) := padj) %>%
            dplyr::rename(!!paste0("Significance_", contrast) := Significance)
        contrastl[[contrast]] <- ddsresmod
    }
    ddsrestetype <- Reduce(function(x, y) merge(x, y, by = c("gene_id", "gene_or_te")), contrastl, accumulate = FALSE)
    ddsrestetype <- ddsrestetype %>% mutate(tecounttype = tecounttype)
    RESLIST[[tecounttype]] <- ddsrestetype
}
ddsfinal <- bind_rows(RESLIST)

COUNTLIST <- list()
for (tecounttype in tecounttypes) {
    conditions <- peptable$condition %>% unique()
    ddscounts <- read_csv(paste(params$inputdir, tecounttype, "counttablesizenormed.csv", sep = "/")) %>%
        rename(gene_id = ...1) %>%
        mutate(tecounttype = tecounttype)
    COUNTLIST[[tecounttype]] <- ddscounts
}
countsfinal <- bind_rows(COUNTLIST)
# merge counts and dds, then add TE annotations


# merge counts and dds, then add TE annotations
resultsdf <- full_join(ddsfinal, countsfinal, by = c("gene_id", "tecounttype"))


write_tsv(resultsdf, outputs$resultsdf)
