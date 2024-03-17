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


### functions
plotSave <- function(path, plot, width = 6, height = 6) {
    dir.create(dirname(path), recursive = TRUE)
    png(path, width = width, height = height, units = "in", res = 300)
    print(plot)
    dev.off()
}


conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]



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
            "outputdir" = "results/agg/deseq2"
        ), env = globalenv())
        assign("inputs", list(
            counts = "outs/agg/featurecounts_genes.txt",
        ), env = globalenv())
        assign("outputs", list(outfile = "results/agg/deseq2/featurecounts_genesoutfile.txt"), env = globalenv())
    }
)

if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}


### inputs

counttype <- params[["counttype"]]
coldata <- read.csv(params[["sample_table"]])
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]

print(counttype)
cts <- read.delim(inputs$counts)
rownames(cts) <- cts$Geneid
cts <- dplyr::select(cts, -Geneid)
cnames <- colnames(cts)
###

condition <- coldata$condition
colnames(cts) <- coldata$sample_name
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition
)
# used to filter here
# this sets prol as the reference level since its first in the vector
dds$condition <- factor(dds$condition, levels = levels)

####
if (params$paralellize_bioc) {
    dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
} else {
    dds <- DESeq(dds)
}

keep <- rowSums(counts(dds)) >= 10
ddsfiltered <- dds[keep, ]
####
resultsNames(dds) # lists the coefficients
counttablesizenormed <- counts(dds, normalized = T)

############
for (contrast in contrasts) {
    res <- results(dds, name = contrast)
    # or to shrink log fold changes association with condition:
    if (params$paralellize_bioc) {
        resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm", parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
    } else {
        resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")
    }
    # Save the plot as a png file
    library(EnhancedVolcano)
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
    plotSave(paste(outputdir, counttype, contrast, "deplot.png", sep = "/"), p, 8, 8)


    resOrdered <- res[order(res$pvalue), ]
    resSig <- subset(resOrdered, padj < 0.1)
    rld <- vst(dds, blind = FALSE)

    write.csv(as.data.frame(resOrdered), file = paste(outputdir, counttype, contrast, "results.csv", sep = "/"))
    write.csv(as.data.frame(resSig), file = paste(outputdir, counttype, contrast, "resultsSig.csv", sep = "/"))
    write.csv(as.data.frame(counttablesizenormed), file = paste(outputdir, counttype, contrast, "counttablesizenormed.csv", sep = "/"))
    write.csv(as.data.frame(assay(rld)), file = paste(outputdir, counttype, contrast, "rlogcounts.csv", sep = "/"))


    active_conditions <- str_split_1(gsub("condition_", "", contrast), "_vs_")
    tmp <- coldata %>% filter(str_detect(condition, active_conditions[1]) | str_detect(condition, active_conditions[2]))
    df <- data.frame(condition = tmp$condition)
    rownames(df) <- tmp$sample_name
    topexpressed <- order(rowMeans(counttablesizenormed[, rownames(df)]),
        decreasing = TRUE
    )[1:40]
    topSig <- resOrdered %>%
        head(n = 40) %>%
        rownames()
    topVarGenes <- head(order(rowVars(assay(rld)[, rownames(df)]), decreasing = TRUE), 500)
    topVarGenes40 <- head(order(rowVars(assay(rld)[, rownames(df)]), decreasing = TRUE), 40)
    interestingsets <- list("Top 40 Variable Genes" = topVarGenes40, "Top Highly Expressed Genes" = topexpressed, "Top Differentially Expressed Genes" = topSig, "Top Variable Genes" = topVarGenes)


    for (name in names(interestingsets)) {
        for (scale in c("none", "row")) {
            for (binary in c(TRUE, FALSE)) {
                if (binary) {
                    blabel <- "rownames"
                } else {
                    blabel <- "norownames"
                }
                set <- interestingsets[[name]]
                dirname <- file.path(outputdir, counttype, contrast)
                dir.create(dirname, recursive = TRUE)
                filename <- paste0(name, blabel, scale, ".png")
                path <- file.path(dirname, filename)
                png(path, width = 10, height = 8, units = "in", res = 300)
                pheatmap(assay(rld)[set, rownames(df)], annotation_col = df, scale = scale, show_rownames = binary, main = name)
                dev.off()
            }
        }
    }
}

rld <- vst(dds, blind = FALSE)

tmp <- coldata
df <- data.frame(condition = tmp$condition)
rownames(df) <- tmp$sample_name
topexpressed <- order(rowMeans(counttablesizenormed),
    decreasing = TRUE
)[1:40]
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 500)
topVarGenes40 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 40)
interestingsets <- list("Top 40 Variable Genes" = topVarGenes40, "Top Highly Expressed Genes" = topexpressed, "Top Variable Genes" = topVarGenes)


for (name in names(interestingsets)) {
    for (scale in c("none", "row")) {
        for (binary in c(TRUE, FALSE)) {
            if (binary) {
                blabel <- "rownames"
            } else {
                blabel <- "norownames"
            }
            set <- interestingsets[[name]]
            dirname <- file.path(outputdir, counttype, "plots")
            dir.create(dirname, recursive = TRUE)
            filename <- paste0(name, blabel, scale, ".png")
            path <- file.path(dirname, filename)
            png(path, width = 10, height = 8, units = "in", res = 300)
            pheatmap(assay(rld)[set, ], annotation_col = df, scale = scale, show_rownames = binary, main = name)
            dev.off()
        }
    }
}
####
# for pca ill use the ddsfiltered but will keep the unfiltered for repeat analysis purposes
# most repeats just have a couple counts - these should not be thrown away, they can be aggregated
vst <- assay(vst(ddsfiltered))
vst_full <- vst(ddsfiltered)
sampleDists <- dist(t(vst))

## PCA plots
pcaObj <- pca(vst, metadata = colData(dds), removeVar = 0.1)

screep <- screeplot(pcaObj, title = "") +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, counttype, "plots", "screeplot.png", sep = "/"), width = 5, height = 6, units = "in", res = 300)
print(screep)
dev.off()

loadingsp <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 4
) +
    theme(legend.position = "none") +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, counttype, "plots", "loadings.png", sep = "/"), width = 6, height = 6, units = "in", res = 300)
print(loadingsp)
dev.off()

pcaplot <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
pcap <- pcaplot +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, counttype, "plots", "pcaplot.png", sep = "/"), width = 7, height = 6, units = "in", res = 300)
print(pcap)
dev.off()

pcaplot34 <- biplot(pcaObj,
    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
pcap34 <- pcaplot34 +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, counttype, "plots", "pcaplotPC34.png", sep = "/"), width = 7, height = 6, units = "in", res = 300)
print(pcap)
dev.off()

legend <- get_legend(
    # create some space to the left of the legend
    pcap + theme(legend.box.margin = margin(0, 0, 0, 1), legend.position = "right")
)

prow <- plot_grid(screep + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA)),
    loadingsp + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA)),
    pcap + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 0.5, fill = NA)),
    nrow = 1,
    rel_widths = c(1, 1, 1), labels = "AUTO",
    align = "vh",
    axis = "bt"
) + theme(legend.position = "none")
p <- plot_grid(prow, legend, nrow = 1, rel_widths = c(3, 0.4))

png(paste(outputdir, counttype, "plots", "PCAgrid.png", sep = "/"), width = 16, height = 6, units = "in", res = 300)
print(p)
dev.off()

# png("eigencorr.png", width=10, height=8)
# eigencorplot(p, metavars = c('condition', 'sample_name'))
# dev.off()


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst_full$condition, vst_full$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(paste(outputdir, counttype, "plots", "heatmapplot.png", sep = "/"), width = 10, height = 8, units = "in", res = 300)
pheatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
)
dev.off()

####

x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)
