source("~/data/common/myDefaults.r")
module_name <- "srna"
load(sprintf("conf/colors_%s.RData", module_name))
conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]

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



tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "tecounttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = sprintf("%s/results/agg/deseq_telescope", conf$prefix),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())

        assign("inputs", list(
            counts = sprintf("%s/outs/agg/featurecounts_genes/counts.txt", conf$prefix),
            rte_counts = sprintf("%s/outs/%s/telescope/telescope-run_stats.tsv", conf$prefix, conf$samples)
        ), env = globalenv())
    }
)



tecounttype <- params[["tecounttype"]]

print(tecounttype)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
countspath <- paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/")
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

if (tecounttype == "telescope_multi") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (path in inputs$rte_counts) {
        bounddf <- full_join(bounddf, read_delim(path, comment = "#", col_names = FALSE) %>% dplyr::select(X1, X3) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

if (tecounttype == "telescope_unique") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (path in inputs$rte_counts) {
        bounddf <- full_join(bounddf, read_delim(path, comment = "#", col_names = FALSE) %>% dplyr::select(X1, X6) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]

gene_cts <- read.delim(inputs$counts)

length(gene_cts$gene_id)
colnames(gene_cts) <- c("gene_id", conf$samples)


cts <- rbind(gene_cts, as.data.frame(bounddf1 %>% replace(is.na(.), 0)))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# keep only genes with counts in at least one sample
# cts <- cts[rowSums(cts > 0) != 0, ]
# rounding since genes are allowed fractional counts
cts <- cts %>% mutate(across(everything(), ~ as.integer(round(.))))
cts[rownames(cts) == "CDKN1A", ]
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~batch + condition
)

colData(dds)
# sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

# I estimate the size factors using genes, and not RTEs, since there are 5M repeats and most have very very low counts
dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)
is.na(colnames(dds)) %>% sum()
colData(dds)
ddsrtes <- dds[!(rownames(dds) %in% gene_cts$gene_id), ]
ddsgenes <- dds[rownames(dds) %in% gene_cts$gene_id, ]



####
ddsrteslist <- list()
ddsgeneslist <- list()
# determine all the DESeq calls that will need to be run, a different one for each base level that is used in contrasts
baselevels <- contrasts %>%
    str_extract("vs_\\w+") %>%
    str_remove("vs_") %>%
    unique()
for (baselevel in baselevels) {
    levels_temp <- c(baselevel, levels[levels != baselevel])
    # this sets the reference level since its first in the vector
    dds$condition <- factor(dds$condition, levels = levels_temp)
    ddsrtes$condition <- factor(ddsrtes$condition, levels = levels_temp)
    ddsgenes$condition <- factor(ddsgenes$condition, levels = levels_temp)
    if (params$paralellize_bioc) {
        # dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
        ddsrtes <- DESeq(ddsrtes, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
        ddsgenes <- DESeq(ddsgenes, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
        ddsrteslist[[baselevel]] <- ddsrtes
        ddsgeneslist[[baselevel]] <- ddsgenes
    } else {
        ddsrtes <- DESeq(ddsrtes)
        ddsgenes <- DESeq(ddsgenes)
        ddsrteslist[[baselevel]] <- ddsrtes
        ddsgeneslist[[baselevel]] <- ddsgenes
    }
}

####
counttablesizenormedrtes <- counts(ddsrteslist[[1]], normalized = TRUE)
counttablesizenormedgenes <- counts(ddsgeneslist[[1]], normalized = TRUE)
counttablesizenormed <- rbind(as.data.frame(counttablesizenormedrtes), as.data.frame(counttablesizenormedgenes))
countspath <- paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/")
dir.create(dirname(countspath), recursive = TRUE, showWarnings = FALSE)
write.csv(counttablesizenormed, file = countspath)
# write.csv(as.data.frame(assay(vst_assaydf)), file = paste(outputdir, tecounttype, "vstcounts.csv", sep = "/"))

# tag PLOTS
deseq_plots <- list()
for (subset in c("rtes", "genes")) {
    if (subset == "rtes") {
        ddstemplist <- ddsrteslist
    } else {
        ddstemplist <- ddsgeneslist
    }
    for (contrast in contrasts) {
        baselevel <- str_extract(contrast, "vs_\\w+") %>% str_remove("vs_")
        ddstemp <- ddstemplist[[baselevel]]
        colData(ddstemp)$condition
        res <- results(ddstemp, name = contrast)
        res <- res[order(res$pvalue), ]
        respath <- paste(outputdir, tecounttype, contrast, sprintf("results_%s.csv", subset), sep = "/")
        dir.create(dirname(respath), recursive = TRUE, showWarnings = FALSE)
        write.csv(as.data.frame(res), file = respath)

        if (subset == "genes") {
            res <- res[rownames(res) %in% gene_cts$gene_id, ]
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
        mysaveandstore(paste(outputdir, tecounttype, subset, contrast, "deplot.png", sep = "/"), 8, 8)
        deseq_plots[[tecounttype]][[subset]][["volcano"]][[contrast]] <- p
    }

    vst <- varianceStabilizingTransformation(ddstemp, blind = FALSE)
    vst_assay <- assay(vst)
    sampleDists <- dist(t(vst_assay))

    ## PCA plots
    pcaObj <- pca(vst_assay, metadata = colData(ddstemp), removeVar = 0.1)

    p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
    mysaveandstore(paste(outputdir, tecounttype, subset, "screeplot.png", sep = "/"), 4, 4)


    p <- plotloadings(pcaObj,
        components = getComponents(pcaObj, seq_len(3)),
        rangeRetain = 0.045, labSize = 4
    ) + mtopen
    mysaveandstore(paste(outputdir, tecounttype, subset, "loadings.png", sep = "/"), 10, 7)


    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
        colby = "condition", legendPosition = "right",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(paste(outputdir, tecounttype, subset, "pca.png", sep = "/"), 4, 4)


    p <- biplot(pcaObj,
        x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
        colby = "condition", legendPosition = "right",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(paste(outputdir, tecounttype, subset, "pca34.png", sep = "/"), 4, 4)



    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep = "-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    p <- pheatmap::pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        col = colors
    )
    mysaveandstore(paste(outputdir, tecounttype, subset, "pheatmap.png", sep = "/"), 4, 4)
}

save(ddsrteslist, file = paste(outputdir, tecounttype, "dds_rtes.RData", sep = "/"))
save(ddsgeneslist, file = paste(outputdir, tecounttype, "dds_genes.RData", sep = "/"))

save(mysaveandstoreplots, file = outputs$plots)
