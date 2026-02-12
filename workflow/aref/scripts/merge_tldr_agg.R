module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)


tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldr = sprintf("aref/A.REF_tldr/%s.table.txt", c(paste0("chr", 1:22), "chrX", "chrY"))
        ), env = globalenv())
        assign("outputs", list(
            df = "aref/A.REF_tldr/A.REF.table.txt"
        ), env = globalenv())
    }
)

dflist <- list()
for (df in inputs$tldr) {
    tt <- read_delim(df)
    dflist[[df]] <- tt
}
dff <- Reduce(bind_rows, dflist)
write_delim(dff, outputs$df)

list.dirs(dirname(outputs$df))
for (dir in gsub(".table.txt", "", inputs$tldr)) {
    dir.create(paste0(dirname(dir), "/A.REF"))
    # Move all files from dir one level up
    system(paste0("mv ", file.path(dir, "*"), " ", paste0(dirname(dir), "/A.REF")))
}
