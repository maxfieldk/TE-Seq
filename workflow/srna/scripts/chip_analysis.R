if (interactive()) {
    module_name <<- "srna"
} else {
    module_name <<- snakemake@params$module_name
}
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)

set.seed(123)

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggVennDiagram")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)
library(rtracklayer)
library(ComplexUpset)
library(patchwork)
library(scales)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "inputdir" = "srna/results/agg/deseq",
            "outputdir" = "srna/chip/results/agg/",
            "counttype" = "telescope_multi",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "srna/outs/agg/featurecounts_rtes_chip/counts.txt"
        ), env = globalenv())
        assign("outputs", list(
            "environment" = "srna/results/agg/repeatanalysis/telescope_multi/repeatanalysisplots_environment.RData"
        ), env = globalenv())
    }
)

outputdir <- params$outputdir
contrasts <- conf$contrasts

rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE,
    append_NI_samplename_modifier = FALSE
)
resultsmito <- read_delim("srna/outs/agg/featurecounts_genes_chip3/counts.txt")
colnames(resultsmito) <- colnames(resultsmito) %>%
    gsub("Geneid", "gene_id", .) %>%
    gsub("srna.outs.", "", .) %>%
    gsub(".bowtie2.*", "", .)

resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
colnames(resultsdf1) <- colnames(resultsdf1) %>%
    gsub("Geneid", "gene_id", .) %>%
    gsub("srna.outs.", "", .) %>%
    gsub(".bowtie2.*", "", .)
resultsdf <- resultsdf1 %>% dplyr::select(gene_id, conf$samples)
rm(resultsdf1)

norm_factors <- read_delim("/users/mkelsey/data/cgasnes/TE-Seq/srna/outs/qc/chip/reads_prop_paired.tsv") %>%
    mutate(sample_name = gsub("Geneid", "gene_id", file) %>%
        gsub("srna.outs.", "", .) %>%
        gsub(".bowtie2.*", "", .)) %>%
    mutate(norm_factor = reads_properly_paired / 1000000)

tidydf <- resultsdf %>% pivot_longer(cols = conf$samples, names_to = "sample_name", values_to = "counts")
tidydf <- tidydf %>% left_join(norm_factors %>% dplyr::select(sample_name, norm_factor))
tidydf <- tidydf %>% mutate(normed_counts = counts / norm_factor)





tdf <- tidydf %>% left_join(rmann)

tdf %>%
    group_by(sample_name) %>%
    summarise(sc = sum(counts))
tdf %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))

tdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(counts))


tdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))


tdf %>%
    filter(rte_subfamily == "L1PA5") %>%
    filter(rte_length_req == "Trnc") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))

tdf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "Trnc") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))

tdf %>%
    filter(rte_subfamily == "AluY") %>%
    filter(rte_length_req == "FL") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))


tdf %>%
    filter(rte_subfamily == "L1PA6") %>%
    group_by(sample_name) %>%
    summarise(sc = sum(normed_counts))


# Load required package
library(glmTMB)

# Count values
counts <- c(100, 120, 110, 150, 160, 170)

# Condition labels
condition <- factor(rep(c("input", "treatment"), each = 3))

# Create data frame
df <- data.frame(
    counts = counts,
    condition = condition
)
# Fit Poisson model
poisson_model <- glmTMB(counts ~ condition, data = df, family = poisson)

summary(poisson_model)
