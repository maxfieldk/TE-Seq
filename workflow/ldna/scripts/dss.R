library("stringr")
library("dplyr")
library("ggplot2")
library("readr")
library("magrittr")
library("purrr")
library("tibble")


module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]

library(DSS)
library(BiocParallel)
library(readr)
require(bsseq)
library(Biostrings)


tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("inputs", list(
            "data" = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            "r_annotation_fragmentsjoined" = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            "r_repeatmasker_annotation" = "annotations/repeatmasker_annotation.csv"
        ), env = globalenv())
        assign("outputs", list(
            dmls_unfiltered = "ldna/results/tables/dmls.CG_m.unfiltered.tsv",
            dmrs_unfiltered = "ldna/results/tables/dmrs.CG_m.unfiltered.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv",
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmrs_bed = "ldna/results/tables/dmrs.CG_m.bed",
            dmrs_hypo_bed = "ldna/results/tables/dmrs_hypo.CG_m.bed",
            dmrs_hyper_bed = "ldna/results/tables/dmrs_hyper.CG_m.bed"
        ), env = globalenv())
        assign("params", list(
            chromosome = "chr10"
        ), env = globalenv())
    }
)



print(sprintf("chromosome is %s", params$chromosome))

sample_dfs <- list()
for (sample in sample_table$sample_name) {
    sample_dfs[[sample]] <- read_delim(grep(sprintf("/%s/", sample), inputs$data, value = TRUE), col_names = TRUE) %>% filter(chr == params$chromosome)
}

BSobj <- makeBSseqData(sample_dfs, names(sample_dfs))
mParam <- MulticoreParam(workers = 12, progressbar = TRUE)

conditions <- conf$levels
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name
# condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name[1]
# condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name[1]
# need to adjust numbering given idiosyncracies of dmltest
dmlTest <- DMLtest(BSobj, group1 = condition2samples, group2 = condition1samples, smoothing = TRUE)

dmlTestNONAMECHANGE <- dmlTest
tryCatch(
    {
        head(dmlTest)
        dmls <- callDML(dmlTest)
        head(dmls)
        dmlcols <- colnames(dmls)
        dmlcols1 <- gsub("1", "_ctwo", dmlcols)
        dmlcols2 <- gsub("2", "_cone", dmlcols1)
        dmlcols3 <- gsub("two", "2", dmlcols2)
        dmlcols4 <- gsub("one", "1", dmlcols3)
        dmlcols5 <- gsub("diff", "diff_c2_minus_c1", dmlcols4)
        colnames(dmls) <- dmlcols5
        head(dmls)
    },
    error = function(e) {
        print("no DMLs")
        dmls <- data.frame()
    }
)

tryCatch(
    {
        dmrs <- callDMR(dmlTest)
        dmrs <- as.data.frame(dmrs)
        head(dmrs)
        dmrcols <- colnames(dmrs)
        dmrcols1 <- gsub("1", "_ctwo", dmrcols)
        dmrcols2 <- gsub("2", "_cone", dmrcols1)
        dmrcols3 <- gsub("two", "2", dmrcols2)
        dmrcols4 <- gsub("one", "1", dmrcols3)
        dmrcols5 <- gsub("diff.Methy", "diff_c2_minus_c1", dmrcols4)
        colnames(dmrs) <- dmrcols5
    },
    error = function(e) {
        print("no DMRs")
        dmrs <- data.frame()
    }
)
# save results
options(scipen = 500)
dir.create(dirname(outputs$dmls), recursive = TRUE, showWarnings = FALSE)
write_delim(dmls, outputs$dmls_unfiltered, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs_unfiltered, delim = "\t", col_names = TRUE)

{
genome_lengths <- fasta.seqlengths(conf$reference)
chromosomesAll <- names(genome_lengths)
nonrefchromosomes <- grep("nonref", chromosomesAll, value = TRUE)
refchromosomes <- grep("^chr", chromosomesAll, value = TRUE)
autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE)
chrX <- c("chrX")
chrY <- c("chrY")

MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
if ("chrY" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
    if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, grep("_chrX_|_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes)
    } else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, grep("_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX)
    }
} else if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
    CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrY, grep("_chrX_", nonrefchromosomes, invert = TRUE, value = TRUE))
    CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrY)
} else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
}
}

dmls <- dmls %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
dmrs <- dmrs %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
dmrs <- dmrs %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, "Hyper", "Hypo"))
dmrs$direction <- factor(dmrs$direction, levels = c("Hyper", "Hypo"))

dmls <- dmls %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, "Hyper", "Hypo"))
dmls$direction <- factor(dmls$direction, levels = c("Hyper", "Hypo"))




write_delim(dmls, outputs$dmls, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs, delim = "\t", col_names = TRUE)