source("~/data/common/myDefaults.r")
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

conf <- configr::read.config(file = "conf/config.yaml")[["lrna"]]


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
            "outputdir" = "results/agg/genomebrowserplots/dorado",
            "regions_of_interest" = "conf/regions_of_interest.bed",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "txdb" = conf$txdb
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "results/agg/deseq/dorado/relaxed/resultsdf.tsv",
            "rnasignalsF" = list(
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro1/alignments/genome/guppy/pro1.F.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro2/alignments/genome/guppy/pro2.F.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro3/alignments/genome/guppy/pro3.F.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen1/alignments/genome/guppy/sen1.F.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen2/alignments/genome/guppy/sen2.F.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen3/alignments/genome/guppy/sen3.F.bw"
            ),
            "rnasignalsR" = list(
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro1/alignments/genome/guppy/pro1.R.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro2/alignments/genome/guppy/pro2.R.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/pro3/alignments/genome/guppy/pro3.R.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen1/alignments/genome/guppy/sen1.R.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen2/alignments/genome/guppy/sen2.R.bw",
                "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen3/alignments/genome/guppy/sen3.R.bw"
            ),
            "dnamethylation" = list(
                "/users/mkelsey/data/Nanopore/p2_1/intermediates/PRO1/methylation/PRO1_CG_m_dss.tsv",
                "/users/mkelsey/data/Nanopore/p2_1/intermediates/SEN1/methylation/SEN1_CG_m_dss.tsv"
            )
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "outfiles/genomebrowserplots.out"
        ), env = globalenv())
    }
)
samples <- conf$samples
sample_table <- read_csv("conf/sample_table.csv")
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

outputdir <- params$outputdir
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap


roi <- import("conf/regions_of_interest.bed")
dertes <- list.files(path = "./results/agg/repeatanalysis/dorado/relaxed/tables/differentially_expressed_elements/condition_SEN_vs_PRO", pattern = "\\.bed$", full.names = TRUE)

grs <- GRanges()
grs <- c(grs, roi)
for (file in dertes) {
    print(file)
    tryCatch(
        {
            gr <- import(file)
            print(gr)
            print(head(gr))
            print(length(gr))
            grs <- c(grs, gr)
        },
        error = function(e) {
            print(e)
        }
    )
}
grsdf <- as.data.frame(grs) %>% tibble()


txdb <- loadDb(params$txdb)
# columns(txdb)
# keys <- keys(txdb) %>% head()
# AnnotationDbi::select(txdb, keys = keys, columns = "TXNAME", keytype = "GENEID")
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
assembly <- assembly(Genome = "custom", TxDb = txdb, OrgDb = "org.Hs.eg.db", BSgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0, gene.id.column = "GENEID", display.column = "GENEID")

condition_colors <- conf$condition_colors
conditions_to_plot <- conf$levels
dnamethylationlist <- list()
for (condition in conditions_to_plot) {
    path <- grep(sprintf("%s", condition), inputs$dnamethylation, value = TRUE)
    df <- read_delim(path, col_names = TRUE, delim = "\t") %>%
        mutate(start = pos) %>%
        mutate(end = pos + 1) %>%
        mutate(pctM = X / N) %>%
        dplyr::select(chr, start, end, pctM)
    dnamethylationlist[[condition]] <- df
}

condition1 <- conditions_to_plot[1]
condition2 <- conditions_to_plot[2]

methdif1 <- left_join(dnamethylationlist[[condition1]] %>% dplyr::rename(condition1 = pctM), dnamethylationlist[[condition2]] %>% dplyr::rename(condition2 = pctM))
methdif <- methdif1 %>%
    mutate(dif = condition2 - condition1)

samples_to_plot <- conf$samples

rnasignallistF <- list()
for (sample in samples_to_plot) {
    path <- grep(sprintf("/%s/", sample), inputs$rnasignalsF, value = TRUE)
    rnasignallistF[[sample]] <- readBigwig(path)
}
rnasignallistR <- list()
for (sample in samples_to_plot) {
    path <- grep(sprintf("/%s/", sample), inputs$rnasignalsR, value = TRUE)
    rnasignallistR[[sample]] <- readBigwig(path)
}

##### plotting config
margin <- 0.2
labelright <- 1
xval <- margin
yval <- margin
pagewidth <- 6
trackwidth <- pagewidth - 2 * xval - labelright
height <- 0.3
y_padding <- 0.05
genesheight <- 2
total_height <- (length(conditions_to_plot) + length(samples_to_plot)) * (height + y_padding) + 2 * yval + genesheight + 0.1

for (element in grsdf$name) {
    row <- grsdf[grsdf$name == element, ]
    chr <- row$seqnames
    strand <- row$strand
    if (strand == "+") {
        start <- row$start - 5000
        end <- row$end + 2000
    } else {
        start <- row$start - 2000
        end <- row$end + 5000
    }

    roi <- GRanges(paste0(chr, ":", start, "-", end))

    maxscore <- 0
    for (ii in c(1:length(samples_to_plot))) {
        score <- GRanges(rnasignallist[[sample]]) %>%
            subsetByOverlaps(roi) %>%
            as.data.frame() %>%
            tibble() %$% score %>%
            max()
        if (score > maxscore) {
            maxscore <- score
        }
        print(maxscore)
    }
    {
        png(sprintf("%s/%s.png", outputdir, element), width = pagewidth, height = total_height, units = "in", res = 300)
        pageCreate(width = pagewidth, height = total_height, default.units = "inches")
        i <- 0
        for (i in c(1:length(conditions_to_plot))) {
            condition <- conditions_to_plot[i]
            print(condition)
            df <- GRanges(dnamethylationlist[[condition]]) %>%
                subsetByOverlaps(roi) %>%
                as.data.frame() %>%
                tibble()
            df <- df %>%
                mutate(rM = rollmean(pctM, 5, na.pad = TRUE, align = "center")) %>%
                filter(!is.na(rM)) %>%
                dplyr::rename(score = rM) %>%
                dplyr::select(seqnames, start, end, score)

            p <- plotSignal(
                data = df,
                assembly = assembly,
                chrom = chr, chromstart = start, chromend = end,
                x = xval, y = yval + (height + y_padding) * (i - 1),
                width = trackwidth, height = height,
                linecolor = condition_colors[[condition]],
                range = c(0, 1),
                scale = TRUE,
                baseline = FALSE,
                default.units = "inches"
            )
            # condition_colors[[condition]]
            plotText(
                label = paste0(condition, " mCpG %"), fontsize = 10, fontcolor = "black",
                x = xval + trackwidth, y = yval + (height + y_padding) * (i - 1) + height / 2, just = c("left", "center"),
                default.units = "inches"
            )
        }
        if (strand == "+") {
            rnasignallist <- rnasignallistF
        } else {
            rnasignallist <- rnasignallistR
        }
        for (ii in c(1:length(samples_to_plot))) {
            sample <- samples_to_plot[ii]
            print(sample)
            condition <- sample_table %>% filter(sample_name == sample) %$% condition
            plotSignal(
                data = rnasignallist[[sample]],
                assembly = assembly,
                chrom = chr, chromstart = start, chromend = end,
                x = xval, y = yval + (height + y_padding) * (i + ii - 1),
                width = trackwidth, height = height,
                baseline = FALSE,
                scale = TRUE,
                range = c(0, ifelse(maxscore == 0, 10, ceiling(1.1 * maxscore))),
                linecolor = condition_colors[[condition]],
                default.units = "inches"
            )
            plotText(
                label = paste0(sample, " RNA CPM"), fonsize = 10, fontcolor = "black",
                x = xval + trackwidth, y = yval + (height + y_padding) * (i + ii - 1) + height / 2, just = c("left", "center"),
                default.units = "inches"
            )
        }
        tryCatch(
            {
                plotTranscripts(
                    chrom = chr, chromstart = start, chromend = end,
                    assembly = assembly,
                    x = xval, y = yval + (height + y_padding) * (i + ii), width = trackwidth, height = genesheight, just = c("left", "top"),
                    default.units = "inches",
                    labels = "gene",
                    fontsize = 6,
                    strandSplit = TRUE
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
            x = xval, y = yval + (height + y_padding) * (i + ii) + genesheight, length = trackwidth,
            default.units = "inches"
        )
        pageGuideHide()
        dev.off()
    }

    difscalefactor <- 3
    total_height <- (length(samples_to_plot) + difscalefactor) * (height + y_padding) + 2 * yval + genesheight + 0.1
    {
        png(sprintf("%s/%s_methdif.png", outputdir, element), width = pagewidth, height = total_height, units = "in", res = 300)
        pageCreate(width = pagewidth, height = total_height, default.units = "inches")
        i <- 0
        df <- GRanges(methdif) %>%
            subsetByOverlaps(roi) %>%
            as.data.frame() %>%
            tibble()
        df <- df %>%
            mutate(rM = rollmean(dif, 5, na.pad = TRUE, align = "center")) %>%
            filter(!is.na(rM)) %>%
            dplyr::rename(score = rM) %>%
            dplyr::select(seqnames, start, end, score)
        plotSignal(
            data = df,
            assembly = assembly,
            chrom = chr, chromstart = start, chromend = end,
            x = xval, y = yval,
            width = trackwidth, height = difscalefactor * height,
            linecolor = "black",
            scale = TRUE,
            baseline = TRUE,
            negData = TRUE,
            default.units = "inches"
        )
        plotText(
            label = paste0(condition2, " - ", condition1, " mCpG"), fontsize = 10, fontcolor = "black",
            x = xval + trackwidth, y = yval + difscalefactor * height / 2, just = c("left", "center"),
            default.units = "inches"
        )
        i <- i + difscalefactor
        if (strand == "+") {
            rnasignallist <- rnasignallistF
        } else {
            rnasignallist <- rnasignallistR
        }
        for (ii in c(1:length(samples_to_plot))) {
            sample <- samples_to_plot[ii]
            print(sample)
            condition <- sample_table %>% filter(sample_name == sample) %$% condition
            plotSignal(
                data = rnasignallist[[sample]],
                assembly = assembly,
                chrom = chr, chromstart = start, chromend = end,
                x = xval, y = yval + (height + y_padding) * (i + ii - 1),
                width = trackwidth, height = height,
                baseline = FALSE,
                scale = TRUE,
                range = c(0, ifelse(maxscore == 0, 10, ceiling(1.1 * maxscore))),
                linecolor = condition_colors[[condition]],
                default.units = "inches"
            )
            plotText(
                label = paste0(sample, " RNA CPM"), fonsize = 10, fontcolor = "black",
                x = xval + trackwidth, y = yval + (height + y_padding) * (i + ii - 1) + height / 2, just = c("left", "center"),
                default.units = "inches"
            )
        }
        tryCatch(
            {
                plotTranscripts(
                    chrom = chr, chromstart = start, chromend = end,
                    assembly = assembly,
                    x = xval, y = yval + (height + y_padding) * (i + ii), width = trackwidth, height = genesheight, just = c("left", "top"),
                    default.units = "inches",
                    labels = "gene",
                    fontsize = 6,
                    strandSplit = TRUE
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
            x = xval, y = yval + (height + y_padding) * (i + ii) + genesheight, length = trackwidth,
            default.units = "inches"
        )
        pageGuideHide()
        dev.off()
    }
}


x <- data.frame()
write_delim(x, outputs$outfile, delim = "\t")

# png("temp.png", width = 6, height = 6, units = "in", res = 300)
# plotMultiSignal(
#     data = inputs$rnasignalsF,
#     binSize = 10,
#     binCap = FALSE,
#     assembly = assembly,
#     chrom = chr, chromstart = start, chromend = end
# )

# dev.off()

# plotGenes(
#     assembly = assembly,
#     chrom = "chr15", chromstart = 210398, chromend = 247535
# )
