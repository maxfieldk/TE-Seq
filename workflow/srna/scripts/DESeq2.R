source("workflow/scripts/defaults.R")

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

{
    #build color palettes:
    sample_palette <- setNames(c("darkgrey", as.character(paletteer_d(conf$default_palette, direction = 1, n = length(conf$samples)-1))), conf$samples)
    contrast_palette <- setNames(paletteer_d(conf$default_palette, direction = 1, n = length(conf$levels)), conf$levels)
    direction_palette <- setNames(c("red", "blue"), c("up", "down"))

    scale_contrasts <- list(scale_fill_manual(values = contrast_palette), scale_color_manual(values = contrast_palette))
    scale_samples <- list(scale_fill_manual(values = samples_palette), scale_color_manual(values = samples_palette))
    scale_directions <- list(scale_fill_manual(values = directions_palette), scale_color_manual(values = directions_palette))
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
    mysaveandstore(paste(outputdir, counttype, contrast, "deplot.png", sep = "/"), w = 8, h = 8, res = 300)

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
                filename <- paste0(name, blabel, scale, ".png")
                path <- file.path(dirname, filename)
                p <- pheatmap(assay(rld)[set, rownames(df)], annotation_col = df, scale = scale, show_rownames = binary, main = name)
                mysaveandstore(path, w = 10, h = 8, res = 300)
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
            filename <- paste0(name, blabel, scale, ".png")
            path <- file.path(dirname, filename)
            p <- pheatmap(assay(rld)[set, ], annotation_col = df, scale = scale, show_rownames = binary, main = name)
            mysaveandstore(path, w = 10, h = 8, res = 300)
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

p <- screeplot(pcaObj, title = "") +
    mtopen
mysaveandstore(paste(outputdir, counttype, "plots", "screeplot.png", sep = "/"), w = 5, h = 6, res = 300)

p <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 4
) +
    theme(legend.position = "none") +
    mtopen
mysaveandstore(paste(outputdir, counttype, "plots", "loadings.png", sep = "/"), w = 6, h = 6, res = 300)

pcaplot <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
p <- pcaplot +
    mtopen
mysaveandstore(paste(outputdir, counttype, "plots", "pcaplot.png", sep = "/"), w = 7, h = 6, res = 300)

pcaplot34 <- biplot(pcaObj,
    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
p <- pcaplot34 +
    mtopen
mysaveandstore(paste(outputdir, counttype, "plots", "pcaplotPC34.png", sep = "/"), w = 7, h = 6, res = 300)


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst_full$condition, vst_full$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

p <- pheatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
)
mysaveandstore(paste(outputdir, counttype, "plots", "heatmapplot.png", sep = "/"), w = 10, h = 8, res = 300)
####

save(mysaveandstoreplots, file = outputs$plots)
x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)
