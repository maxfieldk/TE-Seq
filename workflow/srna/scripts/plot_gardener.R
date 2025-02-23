if (interactive()) {
    module_name <<- "srna"
} else {
    module_name <<- snakemake@params$module_name
}
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)

set.seed(123)


library(igvR)
library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggVennDiagram")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(plotly)
library(DT)
library(ggExtra)
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(plotgardener)
library(AnnotationDbi)
library(zoo)
library(rtracklayer)

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "outputdir" = "srna/results/agg/genomebrowserplots",
            "regions_of_interest" = "conf/integrated_regions_of_interest.bed",
            "rmann" = conf$rmann,
            "txdbrefseq" = "aref/default/A.REF_annotations/refseq.sqlite",
            "txdbrep" = "aref/default/A.REF_annotations/A.REF_repeatmasker.complete.sqlite",
            "txdb" = "aref/default/A.REF_annotations/A.REF_repeatmasker_refseq.complete.sqlite",
            "counttype" = "telescope_multi"
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "srna/results/agg/deseq/resultsdf.tsv",
            bwF = sprintf("srna/outs/%s/star_output/%s.F.bw", conf$samples, conf$samples),
            bwR = sprintf("srna/outs/%s/star_output/%s.R.bw", conf$samples, conf$samples)
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "outfiles/genomebrowserplots.out"
        ), env = globalenv())
    }
)

tryCatch(
    {
        assign("norm_by_aligned_reads", read_delim("srna/qc/multiqc/multiqc_data/samtools-stats-dp.txt") %>%
            filter(Sample %in% sample_table$sample_name) %>%
            mutate(meanmandp = mean(`Mapped &amp; paired`)) %>%
            mutate(scale_factor = `Mapped &amp; paired` / meanmandp) %>%
            dplyr::select(Sample, scale_factor, `Mapped &amp; paired`, meanmandp) %>%
            dplyr::rename(sample_name = Sample), env = globalenv())
    },
    error = function(e) {
        assign("norm_by_aligned_reads", read_delim("srna/qc/multiqc/multiqc_data/multiqc_samtools_stats.txt") %>%
            filter(Sample %in% sample_table$sample_name) %>%
            mutate(meanmandp = mean(reads_mapped_and_paired)) %>%
            mutate(scale_factor = reads_mapped_and_paired / meanmandp) %>%
            dplyr::select(Sample, scale_factor, reads_mapped_and_paired, meanmandp) %>%
            dplyr::rename(sample_name = Sample), env = globalenv())
    }
)

paths_bwF <- inputs$bwF
paths_bwR <- inputs$bwR


grs_list <- list()
grs_total_score <- list()
for (sample in sample_table$sample_name) {
    bwF <- import(grep(sprintf("/%s/", sample), paths_bwF, value = TRUE))
    strand(bwF) <- "+"
    bwR <- import(grep(sprintf("/%s/", sample), inputs$bwR, value = TRUE))
    strand(bwR) <- "-"
    grstemp <- c(bwF, bwR)
    grstemp <- grstemp[!grepl("^NI", seqnames(grstemp))]
    seqlevels(grstemp, pruning.mode = "coarse") <- seqlevelsInUse(grstemp)
    mcols(grstemp)$sample_name <- sample
    score <- grstemp$score %>% sum()
    grs_total_score[[sample]] <- score
    mcols(grstemp)$score <- mcols(grstemp)$score / score
    # mcols(grstemp)$score <- mcols(grstemp)$score / norm_by_aligned_reads$scale_factor[norm_by_aligned_reads$sample_name == sample]
    grs_list[[sample]] <- grstemp
}
grs <- Reduce(c, grs_list)


counttype <- params$counttype
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)

resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
resultsdf <- resultsdf1 %>%
    filter(counttype == !!counttype) %>%
    full_join(rmann)

txdbrefseq <- loadDb(params$txdbrefseq)
txdbrep <- loadDb(params$txdbrep)
seqlevels(txdbrep)
# columns(txdb)
# keys <- keys(txdb) %>% head()
# AnnotationDbi::select(txdb, keys = keys, columns = "TXNAME", keytype = "GENEID")
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
assembly <- assembly(Genome = "custom", TxDb = txdbrefseq, OrgDb = "org.Hs.eg.db", BSgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0, gene.id.column = "GENEID", display.column = "GENEID")
assemblyrep <- assembly(Genome = "custom", TxDb = txdbrep, OrgDb = "org.Hs.eg.db", BSgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0, gene.id.column = "GENEID", display.column = "GENEID")

samples_to_plot <- conf$samples

rnasignallistF <- list()
for (sample in samples_to_plot) {
    path <- grep(sprintf("/%s/", sample), inputs$bwF, value = TRUE)
    temp <- readBigwig(path)
    temp$strand <- "+"
    rnasignallistF[[sample]] <- temp
}
rnasignallistR <- list()
for (sample in samples_to_plot) {
    path <- grep(sprintf("/%s/", sample), inputs$bwR, value = TRUE)
    temp <- readBigwig(path)
    temp$strand <- "-"
    rnasignallistR[[sample]] <- temp
}

gois <- resultsdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    filter(if_any(starts_with("padj_"), ~ . <= 0.05)) %$%
    gene_id

resultsdf_unique <- resultsdf1 %>%
    filter(counttype == "telescope_unique") %>%
    full_join(rmann)
p <- resultsdf_unique %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    filter(if_any(starts_with("padj_"), ~ . <= 0.05))
p <- resultsdf_unique %>%
    filter(gene_id %in% gois) %>%
    pw()
p <- resultsdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    filter(if_any(starts_with("padj_"), ~ . <= 0.05)) %>%
    dplyr::select(gene_id, conf$samples, starts_with("Significance"), loc_integrative_stranded, nearest_coding_tx, dist_to_nearest_coding_tx, nearest_noncoding_tx, dist_to_nearest_noncoding_tx) %>%
    ggtexttable()
mysaveandstore(sprintf("%s/plotted_insert_df2.pdf", params$outputdir), h = 5, w = 40)

p <- resultsdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    filter(if_any(starts_with("padj_"), ~ . <= 0.05)) %>%
    dplyr::select(gene_id, seqnames, start, end, strand, loc_integrative_stranded, nearest_coding_tx, dist_to_nearest_coding_tx, nearest_noncoding_tx, dist_to_nearest_noncoding_tx) %>%
    ggtexttable()
mysaveandstore(sprintf("%s/plotted_insert_df_geneinfo2.pdf", params$outputdir), h = 5, w = 40)


##### plotting config

plotrm <- function(rnasignallistF, rnasignallistR, row, max_sig = TRUE, flank = 10000, samples_to_plot = conf$samples, levels = conf$levels) {
    margin <- 0.2
    labelright <- 1
    xval <- margin
    yval <- margin
    pagewidth <- 6
    trackwidth <- pagewidth - 2 * xval - labelright
    height <- 0.3
    y_padding <- 0.05
    genesheight <- 2
    total_height <- (length(conf$levels) + length(samples_to_plot)) * (height + y_padding) + 2 * yval + genesheight + genesheight + 0.1

    chr <- row$seqnames
    roi_not_ext_grs <- GRanges(row)
    strand <- row$strand
    start <- row$start - flank
    end <- row$end + flank
    roi <- GRanges(paste0(chr, ":", start, "-", end))
    # genes_in_region <- genes(txdbrep, filter = list(tx_chrom = "chr14")) %>%
    #     subsetByOverlaps(roi) %>%
    #     as.data.frame() %>%
    #     tibble() %>%
    #     filter(gene_id != element)
    # genes_in_region_sense <- genes_in_region %>%
    #     filter(strand == "+") %>%
    #     mutate(color = "blue")
    # genes_in_region_antisense <- genes_in_region %>%
    #     filter(strand == "-") %>%
    #     mutate(color = "green")
    # highlight_df <- tibble(gene = c(row$old_id, genes_in_region_sense$gene_id, genes_in_region_antisense$gene_id), color = c("red", genes_in_region_sense$color, genes_in_region_antisense$color))

    highlight_df <- tibble(gene = c(row$old_id), color = c("red"))
    if (max_sig == TRUE) {
        roi_to_norm <- roi
    } else {
        roi_to_norm <- roi_not_ext_grs
    }
    maxscore <- 0
    for (ii in c(1:length(samples_to_plot))) {
        scoreF <- GRanges(rnasignallistF[[sample]]) %>%
            subsetByOverlaps(roi_to_norm) %>%
            as.data.frame() %>%
            tibble() %$% score %>%
            max()
        scoreR <- GRanges(rnasignallistR[[sample]]) %>%
            subsetByOverlaps(roi_to_norm) %>%
            as.data.frame() %>%
            tibble() %$% score %>%
            max()
        maxscore <- max(scoreF, scoreR, maxscore)
        print(maxscore)
    }

    {
        path <- sprintf("%s/%s.png", params$outputdir, element)
        dir.create(dirname(path), recursive = TRUE)
        pdf(sprintf("%s/%s_flank_%s_maxsig_%s.pdf", params$outputdir, element, flank, ifelse(max_sig == TRUE, "T", "F")), width = pagewidth, height = total_height)
        pageCreate(width = pagewidth, height = total_height, default.units = "inches")
        i <- 0
        for (ii in c(1:length(samples_to_plot))) {
            sample <- samples_to_plot[ii]
            print(sample)
            condition <- sample_table %>% filter(sample_name == sample) %$% condition
            plotSignal(
                data = rnasignallistF[[sample]],
                assembly = assembly,
                chrom = chr, chromstart = start, chromend = end,
                x = xval, y = yval + (height + y_padding) * (i + ii - 1),
                width = trackwidth, height = height / 2,
                baseline = FALSE,
                scale = TRUE,
                range = c(0, ifelse(maxscore == 0, 10, ceiling(1.1 * maxscore))),
                linecolor = "blue",
                default.units = "inches"
            )
            plotText(
                label = paste0(sample, " +"), fonsize = 10, fontcolor = "black",
                x = xval + trackwidth, y = yval + (height + y_padding) * (i + ii - 1) + height / 4,
                just = c("left", "center"),
                default.units = "inches"
            )
            plotSignal(
                data = rnasignallistR[[sample]],
                assembly = assembly,
                chrom = chr, chromstart = start, chromend = end,
                x = xval, y = yval + (height + y_padding) * (i + ii - 1) + height / 2,
                width = trackwidth, height = height / 2,
                baseline = FALSE,
                scale = TRUE,
                range = c(0, ifelse(maxscore == 0, 10, ceiling(1.1 * maxscore))),
                linecolor = "green",
                default.units = "inches"
            )
            plotText(
                label = paste0(sample, " -"), fonsize = 10, fontcolor = "black",
                x = xval + trackwidth, y = yval + (height + y_padding) * (i + ii - 1) + height / 2 + height / 4, just = c("left", "center"),
                default.units = "inches"
            )
        }
        tryCatch(
            {
                plotGenes(
                    chrom = chr, chromstart = start, chromend = end,
                    assembly = assembly,
                    x = xval, y = yval + (height + y_padding) * (i + ii), width = trackwidth, height = genesheight, just = c("left", "top"),
                    default.units = "inches",
                    fontsize = 6
                )
            },
            error = function(e) {
                print(e)
            }
        )
        tryCatch(
            {
                plotGenes(
                    chrom = chr, chromstart = start, chromend = end,
                    assembly = assemblyrep,
                    x = xval, y = yval + (height + y_padding) * (i + ii) + genesheight, width = trackwidth, height = genesheight, just = c("left", "top"),
                    default.units = "inches",
                    fontsize = 6
                )
            },
            error = function(e) {
                print(e)
            }
        )
        plotGenomeLabel(
            chrom = chr,
            chromstart = start, chromend = end,
            assembly = assembly,
            x = xval, y = yval + (height + y_padding) * (i + ii) + genesheight + genesheight, length = trackwidth,
            default.units = "inches"
        )
        pageGuideHide()
        dev.off()
    }
}

for (element in gois) {
    row <- rmann %>% filter(gene_id == element)
    plotrm(rnasignallistF,
        rnasignallistR,
        row,
        flank = 10000,
        max_sig = FALSE,
        samples_to_plot = conf$samples,
        levels = conf$levels
    )
}



x <- data.frame()
write_delim(x, outputs$outfile, delim = "\t")
