source("workflow/scripts/defaults.R")
module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(circlize)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(Biostrings)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            sniffles = sprintf("ldna/intermediates/%s/sniffles/sniffles.vcf", sample_table$sample_name),
            snpeff = sprintf("ldna/intermediates/%s/snpeff/snpeff.pass.vcf", sample_table$sample_name),
            filtered_tldr = sprintf("aref/%s_tldr/%s.table.kept_in_updated_ref.txt", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(plots = "ldna/results/variant_analysis/variant_analysis.rds"
), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
