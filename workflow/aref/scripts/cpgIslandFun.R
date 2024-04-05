source("workflow/scripts/defaults.R")
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(zoo)
library(pryr)
library(circlize)
library(rGREAT)
library(reactome.db)
library(msigdb)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)

tryCatch(
    {
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "ref_cpgislands" = "annotations/cpgIslands.txt"
        ), env = globalenv())
        assign("outputs", list(
            "cpg_islands_fullinfo" = "annotations/cpg_islands.tsv",
            "cpg_islands" = "annotations/cpg_islands.bed",
            "cpgi_shores" = "annotations/cpgi_shores.bed",
            "cpgi_shelves" = "annotations/cpgi_shelves.bed"
        ), env = globalenv())
    }
)

# island analysis
cpgislands <- read_delim(params$ref_cpgislands, col_names = TRUE)
colnames(cpgislands) <- colnames(cpgislands) %>% str_remove("#")
cpg_islands <- GRanges(cpgislands)

##### CODE FOR CPG SHORE/SHELF FROM https://www.r-bloggers.com/2014/03/cpg-island-shelves/
# extract the shore defined by 2000 bp upstream of cpg islands
shore1 <- flank(cpg_islands, 2000)
# extract the shore defined by 2000 bp downstream of cpg islands
shore2 <- flank(cpg_islands, 2000, FALSE)
# perform intersection and combine the shores where they overlap
shore1_2 <- reduce(c(shore1, shore2))
# extract the features (ranges) that are present in shores only and not in cpg_islands (ie., shores not overlapping islands)
cpgi_shores <- setdiff(shore1_2, cpg_islands)
cpgi_shores$name <- "shore"
# extract the shore defined by 4000 bp upstream of cpg islands
shelves1 <- flank(cpg_islands, 4000)
# extract the shore defined by 2000 bp downstream of cpg islands
shelves2 <- flank(cpg_islands, 4000, FALSE)
# perform intersection and combine the shelves where they overlap
shelves1_2 <- reduce(c(shelves1, shelves2))
# create a set of ranges consisting CpG Islands, Shores
island_shores <- c(cpg_islands, cpgi_shores)
# extract the features (ranges) that are present in shelves only and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
cpgi_shelves <- setdiff(shelves1_2, island_shores)
cpgi_shelves$name <- "shelf"
#############
write_delim(cpgislands, outputs$cpg_islands_fullinfo, delim = "\t")

export(cpg_islands, outputs$cpg_islands)
export(cpgi_shores, outputs$cpgi_shores)
export(cpgi_shelves, outputs$cpgi_shelves)
