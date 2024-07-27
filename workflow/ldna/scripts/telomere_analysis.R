module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

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
library(ggpubr)

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
            telometer_out = sprintf("ldna/results/tables/telometer/%s_telometer_4000.tsv", samples)
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/outfiles/bedmethylanalysis.txt"
        ), env = globalenv())
    }
)

dflist <- list()
for (sample_name in samples) {
    df <- read_delim(
        grep(sprintf("/%s", sample_name),
            inputs$telometer_out,
            value = TRUE
        )
    )
    df$sample_name <- sample_name
    dflist[[sample_name]] <- df
}
df <- Reduce(rbind, dflist)

df1 <- df %>%
    mutate(chr_arm = paste0(chromosome, arm)) %>%
    filter(!grepl("^NI", chromosome)) %>%
    mutate(merged_chr = paste0(gsub("_.*", "", gsub("^alt", "", gsub("^chr", "", chromosome)))))

df %>% arrange(-telomere_length)

p <- df1 %>%
    ggboxplot(x = "merged_chr", y = "telomere_length", fill = "sample_name")
mysaveandstore(pl = p, w = 30, h = 4)
