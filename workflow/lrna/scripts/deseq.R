module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
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
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = "lrna/results/agg/deseq",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "normalization_method" = conf$normalization_method,
            "paralellize_bioc" = 8
        ), env = globalenv())
        assign("outputs", list(
            counts_normed = "lrna/results/agg/deseq/telescope_multi/counttablesizenormed.csv",
            sizefactors = "lrna/results/agg/deseq/telescope_multi/sizefactors.csv",
            environment = "lrna/results/agg/deseq/telescope_multi/deseq_environment.RData"
        ), env = globalenv())
        sample_table_temp <- read.csv(conf$sample_table)
        sample_table_temp <- sample_table_temp[match(conf$samples, sample_table_temp$sample_name), ]
        assign("inputs", list(
            counts = "lrna/outs/agg/featurecounts_genes/counts.txt",
            rte_counts = sprintf("lrna/outs/%s/telescope/telescope-run_stats.tsv", conf$samples),
            alignment_stats = sprintf("lrna/qc/%s/%s/%s.%s.%s.sorted.bam.stats.txt",
                sample_table_temp$sample_name, sample_table_temp$rate,
                sample_table_temp$sample_name, sample_table_temp$type,
                sample_table_temp$modification_string)
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
continuous_cov <- colnames(coldata)[grepl("batchCon", colnames(coldata))]
coldata <- coldata %>%
    mutate(across(all_of(continuous_cov), ~ as.numeric(scale(.)))) # center and scale continuous covariates

if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

if (counttype == "telescope_multi") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- gsub("run_stats", "TE_counts", grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE))
        tdf <- read_delim(path, comment = "#", col_names = TRUE) %>%
            dplyr::rename(gene_id = transcript)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

if (counttype == "telescope_unique") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE)
        tdf <- read_delim(path, comment = "#", col_names = FALSE) %>%
            dplyr::select(X1, X6) %>%
            dplyr::rename(gene_id = X1)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

# Gene counts from featurecounts (relaxed), RTE counts from telescope
bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]

gene_cts <- read.delim(inputs$counts)
str(gene_cts)
length(gene_cts$gene_id)
colnames(gene_cts) <- colnames(gene_cts) %>%
    gsub("Geneid", "gene_id", .) %>%
    gsub("lrna.intermediates.", "", .) %>%
    gsub(".alignments.*", "", .)
gene_cts <- gene_cts %>% dplyr::select(gene_id, conf$samples)

cts <- bind_rows(gene_cts, as.data.frame(bounddf1 %>% replace(is.na(.), 0)))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# remove rows with NA, round to integer
cts <- cts %>% na.omit()
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

normalization_method <- params[["normalization_method"]]
if (normalization_method == "library_size") {
    # Parse mapped reads from individual samtools stats files
    norm_by_aligned_reads <- tibble(sample_name = character(), reads_mapped = numeric())
    for (stats_file in inputs$alignment_stats) {
        stats_lines <- readLines(stats_file)
        mapped_line <- grep("^SN\treads mapped:", stats_lines, value = TRUE)
        reads_mapped <- as.numeric(gsub(".*:\t(\\d+).*", "\\1", mapped_line))
        sname <- basename(dirname(dirname(stats_file)))
        norm_by_aligned_reads <- bind_rows(norm_by_aligned_reads, tibble(sample_name = sname, reads_mapped = reads_mapped))
    }
    norm_by_aligned_reads <- norm_by_aligned_reads %>%
        mutate(mean_mapped = mean(reads_mapped)) %>%
        mutate(scale_factor = reads_mapped / mean_mapped)
    sizeFactors(dds) <- setNames(norm_by_aligned_reads$scale_factor, norm_by_aligned_reads$sample_name)[colnames(dds)]
} else {
    # median of ratios using genes as control (default)
    dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)
}

sf <- as.data.frame(sizeFactors(dds)) %>%
    as_tibble(rownames = "sample_name") %>%
    dplyr::rename(sizefactor = `sizeFactors(dds)`)
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
    str_extract("vs_.*") %>%
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

if (any(grepl("batch", colnames(coldata)))) {
    if (sum(grepl("batch", colnames(coldata))) == 2) {
        batches <- grep("batch", colnames(coldata), value = TRUE)
        batch_vector <- coldata[[batches[1]]]
        batch2_vector <- coldata[[batches[2]]]
        counttablesizenormed <- removeBatchEffect(counttablesizenormedbatchnotremoved, batch = batch_vector, batch2 = batch2_vector, design = model.matrix(~ coldata$condition))
    } else {
        batches <- grep("batch", colnames(coldata), value = TRUE)
        counttablesizenormed <- removeBatchEffect(counttablesizenormedbatchnotremoved, batch = coldata[[batches[1]]], design = model.matrix(~ coldata$condition))
    }
    countsbatchnotremovedpath <- paste(outputdir, counttype, "counttablesizenormedbatchnotremoved.csv", sep = "/")
    dir.create(dirname(countsbatchnotremovedpath), recursive = TRUE, showWarnings = FALSE)
    write.csv(counttablesizenormedbatchnotremoved, file = countsbatchnotremovedpath)
} else {
    counttablesizenormed <- counttablesizenormedbatchnotremoved
}

dir.create(dirname(countspath), recursive = TRUE, showWarnings = FALSE)
write.csv(counttablesizenormed, file = countspath)

# write.csv(as.data.frame(assay(vst_assaydf)), file = paste(outputdir, counttype, "vstcounts.csv", sep = "/"))

# tag PLOTS

for (batchnormed in c("yes", "no")) {
    if (batchnormed == "yes" & !(any(grepl("batch", colnames(coldata))))) {
        next
    }

    for (subset in c("rtes", "genes")) {
        if (subset == "rtes") {
            ddstemplist <- ddsrteslist
        } else {
            ddstemplist <- ddsgeneslist
        }
        for (contrast in contrasts) {
            baselevel <- str_extract(contrast, "vs_.*") %>% str_remove("vs_")
            ddstemp <- ddstemplist[[baselevel]]
            colData(ddstemp)$condition
            res <- results(ddstemp, name = contrast)
            res <- res[order(res$pvalue), ]


            if (subset == "genes") {
                res <- res[rownames(res) %in% gene_cts$gene_id, ]
            } else {
                res <- res[!(rownames(res) %in% gene_cts$gene_id), ]
            }

            respath <- paste(outputdir, counttype, contrast, sprintf("results_%s.csv", subset), sep = "/")
            dir.create(dirname(respath), recursive = TRUE, showWarnings = FALSE)
            write.csv(as.data.frame(res), file = respath)

            DE_UP <- res %>%
                as.data.frame() %>%
                tibble() %>%
                filter(log2FoldChange > 0) %>%
                filter(padj < 0.05) %>%
                nrow()
            DE_DOWN <- res %>%
                as.data.frame() %>%
                tibble() %>%
                filter(log2FoldChange < 0) %>%
                filter(padj < 0.05) %>%
                nrow()
            TOTAL <- res %>%
                as.data.frame() %>%
                tibble() %>%
                nrow()
            p <- EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = c(""),
                title = contrast,
                drawConnectors = TRUE,
                x = "log2FoldChange",
                y = "padj",
                legendPosition = "none",
                ylab = expression(-Log[10] ~ P["adj"]),
                pCutoff = 0.05,
            ) + mtopen + labs(subtitle = NULL, caption = sprintf("DE UP: %s\nDE DOWN: %s\nTOTAL: %s", DE_UP, DE_DOWN, TOTAL)) + theme(legend.position = "none")
            mysaveandstore(paste(outputdir, counttype, subset, contrast, "deplot.pdf", sep = "/"), 6, 6)
            mysaveandstore(paste(outputdir, counttype, subset, contrast, "deplot.pdf", sep = "/"), raster = TRUE, 6, 6)

            p <- DESeq2::plotMA(res, alpha = 0.05) + mtclosed
            mysaveandstore(paste(outputdir, counttype, subset, contrast, "maplot.pdf", sep = "/"), 6, 6)
            mysaveandstore(paste(outputdir, counttype, subset, contrast, "maplot.pdf", sep = "/"), raster = TRUE, 6, 6)
        }
        tryCatch(
            {
                vst <- varianceStabilizingTransformation(ddstemp, blind = TRUE)
                vst_assay <- assay(vst)

                if (batchnormed == "yes") {
                    if (sum(grepl("batch", colnames(coldata))) == 2) {
                        batches <- grep("batch", colnames(coldata), value = TRUE)
                        batch_vector <- coldata[[batches[1]]]
                        batch2_vector <- coldata[[batches[2]]]
                        vst_assay <- removeBatchEffect(vst_assay, batch = batch_vector, batch2 = batch2_vector, design = model.matrix(~ colData(ddstemp)$condition))
                    } else {
                        batches <- grep("batch", colnames(coldata), value = TRUE)
                        vst_assay <- removeBatchEffect(vst_assay, batch = coldata[[batches[1]]], design = model.matrix(~ colData(ddstemp)$condition))
                    }
                }

                p <- vsn::meanSdPlot(vst_assay)
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "vstmeansdplot.pdf", sep = "/"), 6, 6)

                sampleDists <- dist(t(vst_assay))

                ## PCA plots
                pcaObj <- pca(vst_assay, metadata = colData(ddstemp), removeVar = 0.1)

                p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "screeplot.pdf", sep = "/"), 4, 4)


                p <- plotloadings(pcaObj,
                    components = getComponents(pcaObj, seq_len(3)),
                    rangeRetain = 0.045, labSize = 4
                ) + mtopen
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "loadings.pdf", sep = "/"), 10, 7)

                if (any(grepl("batch", colnames(coldata)))) {
                    p <- biplot(pcaObj,
                        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                        colby = "condition", legendPosition = "right", shape = "batch",
                        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                    ) + mtopen + scale_conditions
                    mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca.pdf", sep = "/"), 5, 5)

                    p <- biplot(pcaObj,
                        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                        colby = "batch", legendPosition = "right",
                        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                    ) + mtopen
                    mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_batch.pdf", sep = "/"), 5, 5)

                    p <- biplot(pcaObj,
                        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                        colby = "batch", legendPosition = "right",
                        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                    ) + mtopen
                    mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_batch_large.pdf", sep = "/"), 16, 16)

                    p <- pairsplot(pcaObj, colby = "batch", title = "Batch", legendPosition = "right")
                    mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_pairs_batch.pdf", sep = "/"), 15, 15)

                    p <- eigencorplot(pcaObj, metavars = c("batch", "condition"))
                    mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_pairs.pdf", sep = "/"), 8, 4)
                } else {
                    print("no batch")
                }

                p <- biplot(pcaObj,
                    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                    colby = "condition", legendPosition = "right",
                    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                ) + mtopen + scale_conditions
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca.pdf", sep = "/"), 5, 5)

                p <- biplot(pcaObj,
                    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                    colby = "condition", legendPosition = "right",
                    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                ) + mtopen + scale_conditions
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_large.pdf", sep = "/"), 16, 16)


                # p <- pairsplot(pcaObj, colby = "condition", title = 'Condition', legendPosition = "right")
                # mysaveandstore(paste(outputdir, counttype, subset,sprintf("batchRemoved_%s", batchnormed), "pca_pairs_condition.pdf", sep = "/"), 15, 15)


                p <- biplot(pcaObj,
                    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
                    colby = "condition", legendPosition = "right",
                    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
                ) + mtopen + scale_conditions
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca34.pdf", sep = "/"), 4, 4)



                sampleDistMatrix <- as.matrix(sampleDists)
                rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep = "-")
                colnames(sampleDistMatrix) <- NULL
                colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

                hm <- sampleDistMatrix %>%
                    Heatmap(
                        name = "Sample Distance",
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_names_rot = 90,
                        col = c("red", "white"),
                        border_gp = gpar(col = "black")
                    )
                p <- wrap_elements(grid.grabExpr(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right")))
                mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "sample_dist_heatmap.pdf", sep = "/"), 4, 4)
            },
            error = function(e) {

            }
        )
    }
}

save(ddsrteslist, file = paste(outputdir, counttype, "dds_rtes.RData", sep = "/"))
save(ddsgeneslist, file = paste(outputdir, counttype, "dds_genes.RData", sep = "/"))

x <- tibble(OUT = "")
write_tsv(x, file = outputs$environment)

# figures: modify plot compositions at will!
# load(outputs$environment)

names(mysaveandstoreplots)




tryCatch(
    {
        library(patchwork)
        names(mysaveandstoreplots)
        p1 <- mysaveandstoreplots[[paste(outputdir, counttype, "genes/batchRemoved_no/pca.pdf", sep = "/")]]
        p2 <- mysaveandstoreplots[[paste(outputdir, counttype, "genes/batchRemoved_no/screeplot.pdf", sep = "/")]]
        p3 <- mysaveandstoreplots[[paste(outputdir, counttype, "genes/batchRemoved_no/sample_dist_heatmap.pdf", sep = "/")]]
        ptch <- ((p1 / p2) | p3) + plot_layout(guides = "collect")
        mysaveandstore(pl = ptch, fn = paste(outputdir, counttype, "figs/f32.pdf", sep = "/"), w = 8, h = 6)

        p1 <- mysaveandstoreplots[[paste(outputdir, counttype, "genes/batchRemoved_no/pca.pdf", sep = "/")]]
        p2 <- mysaveandstoreplots[[paste(outputdir, counttype, "genes/batchRemoved_no/screeplot.pdf", sep = "/")]]
        ptch <- (p1 | p2) + plot_layout(guides = "collect")
        mysaveandstore(pl = ptch, fn = paste(outputdir, counttype, "figs/pca_scree.pdf", sep = "/"), w = 9, h = 3.4)

        ptch <- wrap_plots(mysaveandstoreplots[names(mysaveandstoreplots) %>%
            grep("genes", ., value = TRUE) %>%
            grep("maplot", ., value = TRUE) %>%
            grep(paste(contrasts, collapse = "|"), ., value = TRUE)])
        mysaveandstore(pl = ptch, fn = paste(outputdir, counttype, "figs/maplots_genes.pdf", sep = "/"), raster = TRUE, w = 14, h = 6)

        ptch <- wrap_plots(mysaveandstoreplots[names(mysaveandstoreplots) %>%
            grep("genes", ., value = TRUE) %>%
            grep("deplot", ., value = TRUE) %>%
            grep(paste(contrasts, collapse = "|"), ., value = TRUE)])
        mysaveandstore(pl = ptch, fn = paste(outputdir, counttype, "figs/volcano_genes.pdf", sep = "/"), raster = TRUE, w = 14, h = 6)
    },
    error = function(e) {

    }
)
