set.seed(123)
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
# conf <- configr::read.config(file = "conf/config.yaml")[["lrna"]]
module_name <- snakemake@params$module_name
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "inputdir" = "results/agg/deseq",
            "outputdir" = "results/agg/repeatanalysis"
        ), env = globalenv())
        assign("outputs", list(
            "resultsdf" = "results/agg/repeatanalysis/resultsdf.tsv"
        ), env = globalenv())
    }
)



outputdir <- params$outputdir
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
contrasts <- conf$contrasts
counttypes <- conf$counttypes

ifelse(module_name == "srna", assign("counttypes", conf$counttypes, env = globalenv()), 
ifelse(module_name == "lrna", assign("counttypes", conf$counttypes, env = globalenv()), NA))


lengthreq <- conf$lengthreq
maincontrast <- contrasts[1]



sample_table <- read.csv(conf$sample_table)




##########################



# build dds frames
RESLIST <- list()
for (counttype in counttypes) {
    contrastl <- list()
    for (contrast in contrasts) {
        ddsres_rtes <- read_csv(paste(params$inputdir, counttype, contrast, "results_rtes.csv", sep = "/"))
        ddsres_rtes$gene_or_te <- "repeat"
        ddsres_genes <- read_csv(paste(params$inputdir, counttype, contrast, "results_genes.csv", sep = "/"))
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
    ddsrestetype <- ddsrestetype %>% mutate(counttype = counttype)
    RESLIST[[counttype]] <- ddsrestetype
}
ddsfinal <- bind_rows(RESLIST)

COUNTLIST <- list()
for (counttype in counttypes) {
    conditions <- sample_table$condition %>% unique()
    ddscounts <- read_csv(paste(params$inputdir, counttype, "counttablesizenormed.csv", sep = "/")) %>%
        rename(gene_id = ...1) %>%
        mutate(counttype = counttype)
    COUNTLIST[[counttype]] <- ddscounts
}
countsfinal <- bind_rows(COUNTLIST)
# merge counts and dds, then add TE annotations


# merge counts and dds, then add TE annotations
resultsdf <- full_join(ddsfinal, countsfinal, by = c("gene_id", "counttype"))


write_tsv(resultsdf, outputs$resultsdf)
