module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library("DESeq2")
library("readr")
library("pheatmap")
library("ggplot2")
library("tibble")
library("genefilter")
library("RColorBrewer")
library("cowplot")
library("PCAtools")
library("GGally")
library("tidyr")
library("bcbioRNASeq")
library("DESeqAnalysis")
library(AnnotationHub)
library(KEGGREST)
library(clusterProfiler)
library(AnnotationDbi)
library("biomaRt")
library(stringr)
library("dplyr")
library(EnhancedVolcano)
library(limma)
library(ComplexHeatmap)
library(patchwork)




tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = sprintf("%s/results/agg/deseq", conf$prefix),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())
        assign("outputs", list(
            plots = sprintf("srna/results/agg/deseq/telescope_multi/deseq_plots.RData"),
            sizefactors = "srna/results/agg/deseq/telescope_multi/sizefactors.csv"
        ), env = globalenv())
        assign("inputs", list(
            counts = sprintf("%s/outs/agg/featurecounts_genes/counts.txt", conf$prefix),
            rte_counts = sprintf("%s/outs/%s/telescope/telescope-run_stats.tsv", conf$prefix, conf$samples)
        ), env = globalenv())
    }
)



counttype <- params[["counttype"]]

print(counttype)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
countspath <- outputs$counts_normed
countsbatchnotremovedpath <- paste(outputdir, counttype, "counttablesizenormedbatchnotremoved.csv", sep = "/")
print("countspath")
print(countspath)
coldata <- read.csv(params[["sample_table"]])
samples <- conf$samples
coldata <- coldata[match(conf$samples, coldata$sample_name), ]

if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

df <- read_delim(inputs[["rte_counts"]][1], comment = "#", col_names = FALSE)

if (counttype == "telescope_multi") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (path in inputs$rte_counts) {
        bounddf <- full_join(bounddf, read_delim(path, comment = "#", col_names = FALSE) %>% dplyr::select(X1, X3) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

if (counttype == "telescope_unique") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (path in inputs$rte_counts) {
        bounddf <- full_join(bounddf, read_delim(path, comment = "#", col_names = FALSE) %>% dplyr::select(X1, X6) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]

gene_cts <- read.delim(inputs$counts)
str(gene_cts)
length(gene_cts$gene_id)
colnames(gene_cts) <- c("gene_id", conf$samples)


cts <- rbind(gene_cts, as.data.frame(bounddf1 %>% replace(is.na(.), 0)))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# keep only genes with counts in at least one sample
# cts <- cts[rowSums(cts > 0) != 0, ]
# rounding since genes are not allowed fractional counts
cts <- cts %>% mutate(across(everything(), ~ as.integer(round(.))))
if (any(grepl("batch", colnames(coldata)))) {
    dds <- DESeqDataSetFromMatrix(
        countData = cts,
        colData = coldata,
        design = formula(paste0("~", paste0(grep("batch", colnames(coldata), value = TRUE), collapse = " + "), " + condition"))
    )
} else {
    dds <- DESeqDataSetFromMatrix(
        countData = cts,
        colData = coldata,
        design = ~condition
    )
}

colData(dds)
# sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

# I estimate the size factors using genes, and not RTEs, since there are 5M repeats and most have very very low counts
dds1 <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)
sf <- as.data.frame(sizeFactors(dds1)) %>%
    as_tibble(rownames = "sample_name") %>%
    dplyr::rename(sizefactor = `sizeFactors(dds1)`)
dir.create(dirname(outputs$sizefactors), recursive = TRUE, showWarnings = FALSE)

dds_all_sf <- estimateSizeFactors(dds)
sf_all <- as.data.frame(sizeFactors(dds_all_sf)) %>%
    as_tibble(rownames = "sample_name") %>%
    dplyr::rename(sizefactor = `sizeFactors(dds_all_sf)`)
dir.create(dirname(outputs$sizefactors), recursive = TRUE, showWarnings = FALSE)

is.na(colnames(dds1)) %>% sum()
colData(dds1)
ddsrtes <- dds1[!(rownames(dds1) %in% gene_cts$gene_id), ]
ddsgenes <- dds1[rownames(dds1) %in% gene_cts$gene_id, ]

ddsallrtes <- dds_all_sf[!(rownames(dds_all_sf) %in% gene_cts$gene_id), ]
ddsallgenes <- dds_all_sf[rownames(dds_all_sf) %in% gene_cts$gene_id, ]



####
ddsrteslist <- list()
ddsgeneslist <- list()
ddsallrteslist <- list()
ddsallgeneslist <- list()
ddstogetherlist <- list()

# determine all the DESeq calls that will need to be run, a different one for each base level that is used in contrasts
baselevels <- contrasts %>%
    str_extract("vs_.*") %>%
    str_remove("vs_") %>%
    unique()
for (baselevel in baselevels[[1]]) {
    levels_temp <- c(baselevel, levels[levels != baselevel])
    # this sets the reference level since its first in the vector
    dds$condition <- factor(dds$condition, levels = levels_temp)
    ddsrtes$condition <- factor(ddsrtes$condition, levels = levels_temp)
    ddsgenes$condition <- factor(ddsgenes$condition, levels = levels_temp)
    ddsallrtes$condition <- factor(ddsallrtes$condition, levels = levels_temp)
    ddsallgenes$condition <- factor(ddsallgenes$condition, levels = levels_temp)
    dds_all_sf$condition <- factor(dds_all_sf$condition, levels = levels_temp)
    if (params$paralellize_bioc) {
        # dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
        ddsrtes <- DESeq(ddsrtes, parallel = TRUE, BPPARAM = MulticoreParam(16))
        ddsgenes <- DESeq(ddsgenes, parallel = TRUE, BPPARAM = MulticoreParam(16))
        ddsrteslist[[baselevel]] <- ddsrtes
        ddsgeneslist[[baselevel]] <- ddsgenes

        ddsallrtes <- DESeq(ddsallrtes, parallel = TRUE, BPPARAM = MulticoreParam(16))
        ddsallgenes <- DESeq(ddsallgenes, parallel = TRUE, BPPARAM = MulticoreParam(16))
        ddsallrteslist[[baselevel]] <- ddsallrtes
        ddsallgeneslist[[baselevel]] <- ddsallgenes

        ddstogether <- DESeq(dds_all_sf, parallel = TRUE, BPPARAM = MulticoreParam(16))
        ddstogetherlist[[baselevel]] <- ddstogether

    } else {
        ddsrtes <- DESeq(ddsrtes)
        ddsgenes <- DESeq(ddsgenes)
        ddsrteslist[[baselevel]] <- ddsrtes
        ddsgeneslist[[baselevel]] <- ddsgenes
    }
}

ctsrtes <- counts(ddsrteslist[[1]], normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble()
ctsgenes <- counts(ddsgeneslist[[1]], normalized = TRUE)  %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble()
ctsrtesall <- counts(ddsallrteslist[[1]], normalized = TRUE)  %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble()
ctsgenesall <- counts(ddsallgeneslist[[1]], normalized = TRUE)  %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble()
ctstogether <- counts(ddstogetherlist[[1]], normalized = TRUE)  %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble()
ctsrtestogether <- ctstogether %>% filter(!(gene_id %in% gene_cts$gene_id))
ctsgenestogether <- ctstogether %>% filter((gene_id %in% gene_cts$gene_id))

tctsrtes <- ctsrtes %>% pivot_longer(-gene_id) %>% mutate(norm_type = "only_genes")
tctsrtesall <- ctsrtesall %>% pivot_longer(-gene_id) %>% mutate(norm_type = "all")
tctsrtestogether <- ctsrtestogether %>% pivot_longer(-gene_id) %>% mutate(norm_type = "together")
trtes <- bind_rows(tctsrtes, tctsrtesall, tctsrtestogether)
trtes <- trtes %>% dplyr::rename(sample_name = name, counts = value) %>% left_join(sample_table)

aa <- trtes %>% filter(grepl("L1HS", gene_id)) %>% group_by(norm_type, sample_name) %>% summarise(mean_counts = mean(counts))
aa %>% print(n = Inf)

contrast = contrasts[1]
resrtes <- results(ddsrteslist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "only_genes")
resrtesall <- results(ddsallrteslist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "all")
resrtestogether <- results(ddstogetherlist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% filter(!(gene_id %in% gene_cts$gene_id)) %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "together")

rrtes <- bind_rows(resrtes, resrtesall, resrtestogether)
ar <- rrtes %>% filter(grepl("^L1", gene_id)) 
ar %>% filter(padj < 0.05) %>% group_by(norm_type) %>% summarise(n = n())
rrtes %>% filter(padj < 0.05) %>% group_by(norm_type) %>% summarise(n = n())

resgenes <- results(ddsgeneslist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "only_genes")
resgenesall <- results(ddsallgeneslist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "all")
resgenestogether <- results(ddstogetherlist[[1]], name = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% tibble() %>% filter(!(gene_id %in% gene_cts$gene_id)) %>% dplyr::select(gene_id, padj)%>% mutate(norm_type = "together")

rgenes <- bind_rows(resgenes, resgenesall, resgenestogether)
rgenes %>% filter(padj < 0.05) %>% group_by(norm_type) %>% summarise(n = n())

