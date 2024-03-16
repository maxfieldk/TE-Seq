library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("org.Hs.eg.db")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")

conf <- configr::read.config(file = "conf/config.yaml")["lrna"]


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "inputdir" = "results/agg/deseq",
        ), env = globalenv())
        assign("outputs", list(
            "resultsdf" = "results/agg/deseq/resultsdf.tsv"
        ), env = globalenv())
    }
)


samples <- conf$samples
sample_table <- read_csv("conf/sample_table.csv")
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap
tecounttypes <- conf$tecounttypes
lengthreq <- conf$lengthreq
maincontrast <- contrasts[1]
##########################



# build dds frames
RESLIST <- list()
contrastl <- list()
for (contrast in contrasts) {
    ddsres_rtes <- read_csv(paste(params$inputdir, contrast, "results_rtes.csv", sep = "/"))
    ddsres_rtes$gene_or_te <- "repeat"
    ddsres_genes <- read_csv(paste(params$inputdir, contrast, "results_genes.csv", sep = "/"))
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
ddsrestetype <- Reduce(function(x, y) merge(x, y, by = "gene_id"), contrastl, accumulate = FALSE)

ddsfinal <- ddsrestetype

ddscounts <- read_csv(paste(params$inputdir, "counttablesizenormed.csv", sep = "/")) %>%
    rename(gene_id = ...1)

countsfinal <- ddscounts

resultsdf <- full_join(ddsfinal, countsfinal, by = c("gene_id"))
resultsdf$gene_id <- resultsdf$gene_id %>% gsub(pattern = "gene-", "", x = .)

write_tsv(resultsdf, outputs$resultsdf)
