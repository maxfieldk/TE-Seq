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
library(GenomicFeatures)
library(rtracklayer)

tryCatch(
    {
        inputs <- inputs
        outputs <- outputs
    },
    error = function(e) {
        assign("inputs", list(
            refseq = "aref/annotations/refseq.gff3",
            repeatmasker = "aref/annotations/repeatmasker.complete.gff3",
            genome2bit = "aref/ref.2bit"
        ), env = globalenv())
        assign("outputs", list(
            txdb = "aref/annotations/repeatmasker_refseq.complete.sqlite",
            txdbrefseq = "aref/annotations/refseq.sqlite",
            txdbrepeatmasker = "aref/annotations/repeatmasker.complete.sqlite"
        ), env = globalenv())
    }
)


grs_refseq <- import(inputs$refseq)
grs_repeatmasker <- import(inputs$repeatmasker)
grs <- c(grs_refseq, grs_repeatmasker)
txdb <- makeTxDbFromGRanges(grs)
txdbrepeatmasker <- makeTxDbFromGRanges(grs_repeatmasker)
txdbrefseq <- makeTxDbFromGRanges(grs_refseq)
saveDb(txdb, file = outputs$txdb)
saveDb(txdbrefseq, file = outputs$txdbrefseq)
saveDb(txdbrepeatmasker, file = outputs$txdbrepeatmasker)


# library(Biostrings)
# library(BSgenome)


# seed_files <- system.file("extdata", "GentlemanLab", package = "BSgenome")
# t2tseed <- grep("T2T.*seed", list.files(seed_files, full.names = TRUE), value = TRUE)
# t2tseedlines <- cat(readLines(t2tseed), sep = "\n") %>% as.character()
# gsub("Package: ", "Package: CustomBSgenome Ignore details below")

# keylist <- list(
#     "Package" = "BSgenome.Hsapiens."
# )
