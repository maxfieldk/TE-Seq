source("workflow/scripts/defaults.R")
module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

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




tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "tecounttype" = "relaxed",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = "lrna/results/agg/deseq/relaxed",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())
        assign("outputs", list(
        counts_normed = "lrna/results/agg/deseq/relaxed/counttablesizenormed.csv",
        sizefactors = "lrna/results/agg/deseq/relaxed/sizefactors.csv",
        plots = "lrna/results/agg/deseq/relaxed/deseq_plots.RData", env = globalenv()))
        
        assign("inputs", list(
            counts = sprintf("lrna/intermediates/%s/counts/relaxed/%srtesandgenes.counts.txt", conf$samples, conf$samples)
        ), env = globalenv())
    }
)



tecounttype <- params[["tecounttype"]]

print(tecounttype)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
countspath <- paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/")
countsbatchnotremovedpath <- paste(outputdir, tecounttype, "counttablesizenormedbatchnotremoved.csv", sep = "/")
print("countspath")
print(countspath)

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
gene_cts <- cts %>% filter(grepl("gene-", gene_id))
gene_cts$gene_id <- gene_cts$gene_id %>% gsub("gene-", "", .)

rownames(cts) <- cts$gene_id %>% gsub("gene-", "", .)
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# keep only genes with counts in at least one sample
# cts <- cts[rowSums(cts > 0) != 0, ]
# rounding since genes are not allowed fractional counts
cts <- cts %>% mutate(across(everything(), ~ as.integer(round(.))))
#remove rows with NA
cts <- cts %>% na.omit()

if ("batch" %in% colnames(coldata)) {
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~batch + condition)
} else {
    dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition)
}

colData(dds)
# sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

# I estimate the size factors using genes, and not RTEs, since there are 5M repeats and most have very very low counts
dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)
sf <- as.data.frame(sizeFactors(dds)) %>% as_tibble(rownames = "sample_name") %>% dplyr::rename(sizefactor = `sizeFactors(dds)`)
dir.create(dirname(outputs$sizefactors), recursive = TRUE, showWarnings = FALSE)
write_csv(sf, outputs$sizefactors)

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
counttablesizenormedbatchnotremoved <- rbind(as.data.frame(counttablesizenormedrtes), as.data.frame(counttablesizenormedgenes))
colnames(counttablesizenormedbatchnotremoved) == coldata$sample_name

if ("batch" %in% colnames(coldata)) {
counttablesizenormed <- removeBatchEffect(counttablesizenormedbatchnotremoved, batch=coldata$batch, design=model.matrix(~coldata$condition))
countsbatchnotremovedpath <- paste(outputdir, tecounttype, "counttablesizenormedbatchnotremoved.csv", sep = "/")
dir.create(dirname(countsbatchnotremovedpath), recursive = TRUE, showWarnings = FALSE)
write.csv(counttablesizenormedbatchnotremoved, file = countsbatchnotremovedpath)
} else {
    counttablesizenormed <- counttablesizenormedbatchnotremoved
}

countspath <- paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/")
dir.create(dirname(countspath), recursive = TRUE, showWarnings = FALSE)
write.csv(counttablesizenormed, file = countspath)

# write.csv(as.data.frame(assay(vst_assaydf)), file = paste(outputdir, tecounttype, "vstcounts.csv", sep = "/"))

# tag PLOTS

for (batchnormed in c("yes", "no")) {

if (batchnormed == "yes" & !("batch" %in% colnames(coldata))) { next }

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

        DE_UP <- res %>% as.data.frame() %>% tibble() %>% filter(log2FoldChange > 0) %>% filter(padj < 0.05) %>% nrow()
        DE_DOWN <- res %>% as.data.frame() %>% tibble() %>% filter(log2FoldChange < 0) %>% filter(padj < 0.05) %>% nrow()
        TOTAL <- res %>% as.data.frame() %>% tibble() %>% nrow()
        p <- EnhancedVolcano(res,
        lab = rownames(res),
        selectLab = c(""),
        title = contrast,
        drawConnectors = TRUE,
        x = "log2FoldChange",
        y = "padj",
        legendPosition = 'none',
        ylab = expression(-Log[10] ~ P["adj"]),
                    pCutoff = 0.05,
        ) + mtopen + labs(subtitle = NULL, caption = sprintf("DE UP: %s\nDE DOWN: %s\nTOTAL: %s", DE_UP, DE_DOWN, TOTAL)) + theme(legend.position = "none")
        mysaveandstore(paste(outputdir, tecounttype, subset, contrast, "deplot.png", sep = "/"), 6, 6)

        p <- DESeq2::plotMA(res, alpha = 0.05) + mtclosed
        mysaveandstore(paste(outputdir, tecounttype, subset, contrast, "maplot.png", sep = "/"), 6, 6)   
        }

        vst <- varianceStabilizingTransformation(ddstemp, blind = FALSE)
        vst_assay <- assay(vst)
        if (batchnormed == "yes") {
        vst_assay <- removeBatchEffect(vst_assay, batch = colData(ddstemp)$batch, design = model.matrix(~colData(ddstemp)$condition))
        }

        p <- vsn::meanSdPlot(vst_assay)
        mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "vstmeansdplot.png", sep = "/"), 6, 6)

        sampleDists <- dist(t(vst_assay))

        ## PCA plots
        pcaObj <- pca(vst_assay, metadata = colData(ddstemp), removeVar = 0.1)

        p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
        mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "screeplot.png", sep = "/"), 4, 4)


        p <- plotloadings(pcaObj,
            components = getComponents(pcaObj, seq_len(3)),
            rangeRetain = 0.045, labSize = 4
        ) + mtopen
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "loadings.png", sep = "/"), 10, 7)

        if ("batch" %in% colnames(coldata)) {

            p <- biplot(pcaObj,
                showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "condition", legendPosition = "right", shape = "batch",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mtopen + scale_conditions
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca.png", sep = "/"), 5, 5)

                p <- biplot(pcaObj,
                showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "batch", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mtopen
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_batch.png", sep = "/"), 5, 5)

            p <- biplot(pcaObj,
                showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "batch", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mtopen
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_batch_large.png", sep = "/"), 16, 16)

            p <- pairsplot(pcaObj, colby = 'batch', title = 'Batch', legendPosition = "right")
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_pairs_batch.png", sep = "/"), 15, 15)

            p <- eigencorplot(pcaObj,metavars = c("batch", "condition"))
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_pairs.png", sep = "/"), 8, 4)

        } else {
            print("no batch")
        }

            p <- biplot(pcaObj,
                showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "condition", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mtopen + scale_conditions
    mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca.png", sep = "/"), 5, 5)

    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
        colby = "condition", legendPosition = "right",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_large.png", sep = "/"), 16, 16)


    # p <- pairsplot(pcaObj, colby = "condition", title = 'Condition', legendPosition = "right")
    # mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_pairs_condition.png", sep = "/"), 15, 15)


            p <- biplot(pcaObj,
                x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                colby = "condition", legendPosition = "right",
                labSize = 5, pointSize = 5, sizeLoadingsNames = 5
            ) + mtopen + scale_conditions
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca34.png", sep = "/"), 4, 4)



            sampleDistMatrix <- as.matrix(sampleDists)
            rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep = "-")
            colnames(sampleDistMatrix) <- NULL
            colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

            p <- pheatmap::pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors
            )
            mysaveandstore(paste(outputdir, tecounttype, subset,sprintf("batchRemoved_%s", batchnormed), "pheatmap.png", sep = "/"), 4, 4)
        }
}

save(ddsrteslist, file = paste(outputdir, tecounttype, "dds_rtes.RData", sep = "/"))
save(ddsgeneslist, file = paste(outputdir, tecounttype, "dds_genes.RData", sep = "/"))

save(mysaveandstoreplots, file = outputs$plots)
