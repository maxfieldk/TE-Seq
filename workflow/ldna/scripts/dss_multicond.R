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
            "data" = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(
            dmls_unfiltered = "ldna/results/m/tables/condition_AD_vs_YNG/dmls.chr10.unfiltered.tsv",
            dmrs_unfiltered = "ldna/results/m/tables/condition_AD_vs_YNG/dmrs.chr10.unfiltered.tsv",
            dmls = "ldna/results/m/tables/condition_AD_vs_YNG/dmls.chr10.tsv",
            dmrs = "ldna/results/m/tables/condition_AD_vs_YNG/dmrs.chr10.tsv"
        ), env = globalenv())
        assign("params", list(
            chromosome = "chr10",
            contrast = "condition_AD_vs_YNG"
        ), env = globalenv())
    }
)


# Parse the contrast string to get the test condition and reference
# Format: "condition_TEST_vs_REFERENCE" e.g. "condition_AD_vs_YNG"
contrast_string <- params$contrast
contrast_parts <- str_match(contrast_string, "^condition_(.+)_vs_(.+)$")
condition_test <- contrast_parts[1, 2]
condition_ref <- contrast_parts[1, 3]
print(sprintf("Contrast: %s vs %s (test vs reference)", condition_test, condition_ref))
print(sprintf("Chromosome: %s", params$chromosome))

# Load cell fractions and z-score them
cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")
cell_fractions <- read_csv(inputs$cell_typefractions)
cell_fractions_scaled <- cell_fractions %>%
    mutate(across(all_of(cell_types_cols), ~ as.numeric(scale(.)), .names = "{.col}_z"))

# Build design table with all samples and covariates
design <- sample_table %>%
    left_join(cell_fractions_scaled %>% dplyr::select(sample_name, ends_with("_z")), by = "sample_name")

# Set reference level for condition
design$condition <- relevel(factor(design$condition), ref = condition_ref)
design$sex <- factor(design$sex)
design$ancestry <- factor(design$ancestry)

print(sprintf("All samples: %s", paste(design$sample_name, collapse = ", ")))
print(sprintf("Conditions: %s", paste(levels(design$condition), collapse = ", ")))

# Load DSS data for all samples
sample_dfs <- list()
for (sample in design$sample_name) {
    sample_dfs[[sample]] <- read_delim(grep(sprintf("/%s/", sample), inputs$data, value = TRUE), col_names = TRUE) %>%
        filter(chr == params$chromosome)
}

BSobj <- makeBSseqData(sample_dfs, names(sample_dfs))

# The coefficient to test corresponds to the contrast
# e.g. for condition_AD_vs_YNG with ref=YNG, the coef is "conditionAD"
coef_name <- paste0("condition", condition_test)
print(sprintf("Testing coefficient: %s", coef_name))

# Fit the full GLM with all covariates
# Formula: ~ condition + sex + ancestry + cell type fractions (z-scored)
dss_design <- data.frame(design)

tryCatch(
    {
        DMLfit <- DMLfit.multiFactor(BSobj,
            design = dss_design,
            smoothing = TRUE,
            formula = ~ condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z
        )
        print("Design matrix columns:")
        print(colnames(DMLfit$X))

        # Test the specific contrast coefficient
        dmls <- DMLtest.multiFactor(DMLfit, coef = coef_name)

        dmrs_f1 <- callDMR(dmls, p.threshold = 0.05)
        dmrs_f1$dmr_type <- "t05"
        dmrs_f2 <- callDMR(dmls, p.threshold = 0.05, minCG = 10)
        dmrs_f2$dmr_type <- "t05CG10"
        dmrs_f3 <- callDMR(dmls, p.threshold = 0.01)
        dmrs_f3$dmr_type <- "t01"
        dmrs_f4 <- callDMR(dmls, p.threshold = 0.001)
        dmrs_f4$dmr_type <- "t001"
    },
    error = function(e) {
        print(paste("DMLfit.multiFactor error:", e$message))
        assign("dmls", data.frame(chr = character(), pos = integer(), stat = numeric(), pvals = numeric(), fdrs = numeric()), envir = parent.env(environment()))
        empty_dmrs <- data.frame(chr = character(), start = integer(), end = integer(), length = integer(), nCG = integer(), areaStat = numeric(), dmr_type = character())
        assign("dmrs_f1", empty_dmrs, envir = parent.env(environment()))
        assign("dmrs_f2", empty_dmrs, envir = parent.env(environment()))
        assign("dmrs_f3", empty_dmrs, envir = parent.env(environment()))
        assign("dmrs_f4", empty_dmrs, envir = parent.env(environment()))
    }
)


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
