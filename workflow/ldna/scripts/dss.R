library("stringr")
library("dplyr")
library("ggplot2")
library("readr")
library("magrittr")
library("purrr")
library("tibble")
set.seed(123)

module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("conf/sample_table_source.R")



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
            "data" = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
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
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name %>% as.character()
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name %>% as.character()


design <- data.frame(sample_table)
X <- model.matrix(~condition, design)

if (length(condition1samples) < 3) {
    tryCatch(
        {
            dmlFit <- DMLtest(BSobj, group1 = condition2samples, group2 = condition1samples, smoothing = TRUE)
            head(dmlFit)
            dmls <- callDML(dmlFit)
            dmls <- dmls %>%
                dplyr::rename(pvals = pval, fdrs = fdr) %>%
                dplyr::select(chr, pos, stat, pvals, fdrs)

            dmrs_f1 <- callDMR(dmlFit, p.threshold = 0.05)
            dmrs_f1$dmr_type <- "t05"
            dmrs_f1 <- dmrs_f1 %>%
                dplyr::rename() %>%
                dplyr::select(chr, start, end, length, nCG, areaStat, dmr_type)

            dmrs_f2 <- callDMR(dmlFit, p.threshold = 0.05, minCG = 10)
            dmrs_f2$dmr_type <- "t05CG10"
            dmrs_f2 <- dmrs_f2 %>%
                dplyr::rename() %>%
                dplyr::select(chr, start, end, length, nCG, areaStat, dmr_type)

            dmrs_f3 <- callDMR(dmlFit, p.threshold = 0.01)
            dmrs_f3$dmr_type <- "t01"
            dmrs_f3 <- dmrs_f3 %>%
                dplyr::rename() %>%
                dplyr::select(chr, start, end, length, nCG, areaStat, dmr_type)

            dmrs_f4 <- callDMR(dmlFit, p.threshold = 0.001)
            dmrs_f4$dmr_type <- "t001"
            dmrs_f4 <- dmrs_f4 %>%
                dplyr::rename() %>%
                dplyr::select(chr, start, end, length, nCG, areaStat, dmr_type)
        },
        error = function(e) {
            print("no DMLs")
            dmls <- data.frame()
        }
    )
} else {
    sampleNames(BSobj)
    DMLfit <- DMLfit.multiFactor(BSobj, design = design, smoothing = TRUE, formula = ~condition)
    colnames(DMLfit$X)
    dmls <- dmlFit.multiFactor(DMLfit, term = "condition")
    dmrs_f1 <- callDMR(dmls, p.threshold = 0.05)
    dmrs_f1$dmr_type <- "t05"
    dmrs_f2 <- callDMR(dmls, p.threshold = 0.05, minCG = 10)
    dmrs_f2$dmr_type <- "t05CG10"
    dmrs_f3 <- callDMR(dmls, p.threshold = 0.01)
    dmrs_f3$dmr_type <- "t01"
    dmrs_f4 <- callDMR(dmls, p.threshold = 0.001)
    dmrs_f4$dmr_type <- "t001"
}


dmrs <- bind_rows(dmrs_f1, dmrs_f2, dmrs_f3, dmrs_f4)

options(scipen = 500)
dir.create(dirname(outputs$dmls), recursive = TRUE, showWarnings = FALSE)
write_delim(dmls, outputs$dmls_unfiltered, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs_unfiltered, delim = "\t", col_names = TRUE)

{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE)
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
dmrs <- dmrs %>% mutate(direction = ifelse(areaStat > 0, "Hyper", "Hypo"))
dmrs$direction <- factor(dmrs$direction, levels = c("Hyper", "Hypo"))

dmls <- dmls %>% mutate(direction = ifelse(stat > 0, "Hyper", "Hypo"))
dmls$direction <- factor(dmls$direction, levels = c("Hyper", "Hypo"))


write_delim(dmls, outputs$dmls, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs, delim = "\t", col_names = TRUE)
