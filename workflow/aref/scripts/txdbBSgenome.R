module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(GenomicFeatures)
library(rtracklayer)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            refseq = "aref/default/A.REF_annotations/refseq.gff3",
            repeatmasker = "aref/default/A.REF_annotations/A.REF_repeatmasker.complete.gff3",
            genome2bit = "aref/default/A.REF.2bit"
        ), env = globalenv())
        assign("outputs", list(
            txdb = "aref/default/A.REF_annotations/A.REF_repeatmasker_refseq.complete.sqlite",
            txdbrefseq = "aref/default/A.REF_annotations/refseq.sqlite",
            txdbrepeatmasker = "aref/default/A.REF_annotations/A.REF_repeatmasker.complete.sqlite"
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
