library("stringr")
library("dplyr")
library("ggplot2")
library("readr")
library("magrittr")
library("purrr")
library("tibble")
library(readr)
set.seed(123)

module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]

conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS
chrom_to_keep <- c(paste0("chr", 1:22), "chrX", "chrY") %>% setdiff(conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS)
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
            dmrs = sprintf("ldna/results/m/tables/condition_AD_vs_CTRL/dmrs.%s.tsv", chrom_to_keep),
            dmls = sprintf("ldna/results/m/tables/condition_AD_vs_CTRL/dmls.%s.tsv", chrom_to_keep)
        ), env = globalenv())
        assign("outputs", list(
            dmls = "ldna/results/m/tables/condition_AD_vs_CTRL/dmls.tsv",
            dmrs = "ldna/results/m/tables/condition_AD_vs_CTRL/dmrs.tsv",
            dmrs_bed = "ldna/results/m/tables/condition_AD_vs_CTRL/dmrs.bed",
            dmrs_hypo_bed = "ldna/results/m/tables/condition_AD_vs_CTRL/dmrs_hypo.bed",
            dmrs_hyper_bed = "ldna/results/m/tables/condition_AD_vs_CTRL/dmrs_hyper.bed"
        ), env = globalenv())
    }
)



chr_dmrs <- list()
for (filepath in inputs$dmrs) {
    chr_dmrs[[filepath]] <- read_delim(filepath, col_names = TRUE)
}

chr_dmls <- list()
for (filepath in inputs$dmls) {
    chr_dmls[[filepath]] <- read_delim(filepath, col_names = TRUE)
}
dmrs <- Reduce(bind_rows, chr_dmrs)
dmls <- Reduce(bind_rows, chr_dmls)
dmrs %$% direction %>% table()
dmls %$% direction %>% table()


write_delim(dmls, outputs$dmls, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs, delim = "\t", col_names = TRUE)

write_delim(dmrs %>% dplyr::select(chr, start, end), outputs$dmrs_bed, delim = "\t", col_names = FALSE)
write_delim(dmrs %>% filter(direction == grep("Hypo", dmrs %$% direction %>% unique(), value = TRUE)) %>% dplyr::select(chr, start, end), outputs$dmrs_hypo_bed, delim = "\t", col_names = FALSE)
write_delim(dmrs %>% filter(direction == grep("Hyper", dmrs %$% direction %>% unique(), value = TRUE)) %>% dplyr::select(chr, start, end), outputs$dmrs_hyper_bed, delim = "\t", col_names = FALSE)
