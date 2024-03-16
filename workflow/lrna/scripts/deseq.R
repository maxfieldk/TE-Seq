source("~/data/common/myDefaults.r")

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

conf <- configr::read.config(file = "conf/config.yaml")["lrna"]


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = "results/agg/deseq/dorado/relaxed",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())

        assign("inputs", list(
            counts = sprintf("intermediates/%s/counts/genome/dorado/relaxed/%srtesandgenes.counts.txt", conf$samples, conf$samples)
        ), env = globalenv())
    }
)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
coldata <- read.csv(params[["sample_table"]])
samples <- conf$samples
coldata <- coldata[match(conf$samples, coldata$sample_name), ]

if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

df <- read_delim(inputs[["counts"]][1], comment = "#", col_names = TRUE)
bounddf <- tibble(df[, 1]) %>% rename(gene_id = Geneid)
for (sample in conf$samples) {
    tempdf <- read_delim(grep(sprintf("/%s/", sample), inputs$counts, value = TRUE), comment = "#", col_names = TRUE)
    colnames(tempdf) <- c("gene_id", sample)
    bounddf <- full_join(bounddf, tempdf, by = "gene_id")
}

cts <- as.data.frame(bounddf)
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# keep only genes with counts in at least one sample
# cts <- cts[rowSums(cts > 0) != 0, ]
# rounding since genes are not allowed fractional counts
cts <- cts %>% mutate(across(everything(), ~ as.integer(round(.))))
cts[rownames(cts) == "gene-CDKN1A", ]
cts[rownames(cts) == "gene-IL6", ]
condition <- coldata$condition
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition
)
# sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

# I estimate the size factors using genes, and not RTEs, since there are 5M repeats and most have very very low counts
genes <- grep(pattern = "gene-", x = rownames(cts), value = TRUE)
dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% genes)

ddsrtes <- dds[!(rownames(dds) %in% genes), ]
ddsgenes <- dds[rownames(dds) %in% genes, ]

# this sets prol as the reference level since its first in the vector
dds$condition <- factor(dds$condition, levels = levels)
ddsrtes$condition <- factor(ddsrtes$condition, levels = levels)
ddsgenes$condition <- factor(ddsgenes$condition, levels = levels)
####
if (params$paralellize_bioc) {
    # dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
    ddsrtes <- DESeq(ddsrtes, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
    ddsgenes <- DESeq(ddsgenes, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
} else {
    ddsrtes <- DESeq(ddsrtes)
    ddsgenes <- DESeq(ddsgenes)
}


####
counttablesizenormedrtes <- counts(ddsrtes, normalized = TRUE)
counttablesizenormedgenes <- counts(ddsgenes, normalized = TRUE)
counttablesizenormed <- rbind(as.data.frame(counttablesizenormedrtes), as.data.frame(counttablesizenormedgenes))
countspath <- paste(outputdir, "counttablesizenormed.csv", sep = "/")
dir.create(dirname(countspath), recursive = TRUE, showWarnings = FALSE)
write.csv(counttablesizenormed, file = countspath)
# write.csv(as.data.frame(assay(vst_assaydf)), file = paste(outputdir, "vstcounts.csv", sep = "/"))

# tag PLOTS
deseq_plots <- list()
for (subset in c("genes", "rtes")) {
    if (subset == "rtes") {
        ddstemp <- ddsrtes
    } else {
        ddstemp <- ddsgenes
    }
    for (contrast in contrasts) {
        res <- results(ddstemp, name = contrast)
        res <- res[order(res$pvalue), ]
        respath <- paste(outputdir, contrast, sprintf("results_%s.csv", subset), sep = "/")
        dir.create(dirname(respath), recursive = TRUE, showWarnings = FALSE)
        write.csv(as.data.frame(res), file = respath)

        if (subset == "genes") {
            res <- res[rownames(res) %in% genes, ]
        }

        p <- EnhancedVolcano(res,
            lab = rownames(res),
            selectLab = c(""),
            title = contrast,
            drawConnectors = TRUE,
            x = "log2FoldChange",
            y = "padj",
            ylab = expression(-Log[10] ~ P["adj"]),
            legendPosition = "none",
            pCutoff = 0.05,
            xlim = c(-10, 10),
            ylim = c(-1, 15)
        ) + theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))
        mysave(paste(outputdir, subset, contrast, "deplot.png", sep = "/"), 8, 8)
        deseq_plots[[subset]][["volcano"]][[contrast]] <- p
    }
    keep <- rowSums(counts(ddstemp) >= 2) >= 2
    ddstemp <- ddstemp[keep]
    vst <- varianceStabilizingTransformation(ddstemp, blind = FALSE)
    vst_assay <- assay(vst)
    tryCatch(
        {
            sampleDists <- dist(t(vst_assay))

            ## PCA plots
            pcaObj <- pca(vst_assay, metadata = colData(ddstemp), removeVar = 0.1)

            p <- screeplot(pcaObj, title = "") +
                theme_cowplot() +
                mytheme
            mysave(paste(outputdir, subset, "screeplot.png", sep = "/"), 4, 4)
            deseq_plots[[subset]][["scree"]] <- p


            p <- plotloadings(pcaObj,
                components = getComponents(pcaObj, seq_len(3)),
                rangeRetain = 0.045, labSize = 4
            ) +
                theme(legend.position = "none") +
                mytheme
            mysave(paste(outputdir, subset, "loadings.png", sep = "/"), 4, 4)
            deseq_plots[[subset]][["loadings"]] <- p


            p <- biplot(pcaObj,
                showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "condition", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) +
                theme_gray() +
                mytheme
            mysave(paste(outputdir, subset, "pca.png", sep = "/"), 4, 4)
            deseq_plots[[subset]][["pca"]] <- p


            p <- biplot(pcaObj,
                x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "condition", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mytheme
            mysave(paste(outputdir, subset, "pca34.png", sep = "/"), 4, 4)
            deseq_plots[[subset]][["pca34"]] <- p



            sampleDistMatrix <- as.matrix(sampleDists)
            rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep = "-")
            colnames(sampleDistMatrix) <- NULL
            colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

            p <- pheatmap::pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors
            )
            mysave(paste(outputdir, subset, "pheatmap.png", sep = "/"), 4, 4)
            deseq_plots[[subset]][["dist_heatmap"]] <- p
        },
        error = function(e) {
            print(e)
        }
    )
}

save(deseq_plots, file = paste(outputdir, "deseq_plots.RData", sep = "/"))
save(ddsrtes, file = paste(outputdir, "dds_rtes.RData", sep = "/"))
save(ddsgenes, file = paste(outputdir, "dds_genes.RData", sep = "/"))
