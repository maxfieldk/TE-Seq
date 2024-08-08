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
source("conf/sample_table_source.R")


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            deconv = "ldna/results/wgbs/hg38/tables/deconv.csv"
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/outfiles/deconv.outfile"
        ), env = globalenv())
    }
)

dc <- read_csv(inputs$deconv)
dc1 <- dc %>%
    pivot_longer(cols = -CellType) %>%
    dplyr::rename(sample_name = name) %>%
    mutate(sample_name = gsub("\\..*", "", sample_name)) %>%
    left_join(sample_table)
p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "sample_name", color = "sample_name", position = position_dodge()) +
    scale_samples_unique +
    coord_flip() +
    mtclosedgridv
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all_nofitler.pdf"), pl = p, w = 5, h = 10)


dc1 <- dc %>%
    pivot_longer(cols = -CellType) %>%
    dplyr::rename(sample_name = name) %>%
    mutate(sample_name = gsub("\\..*", "", sample_name)) %>%
    filter(value > 0.05) %>%
    left_join(sample_table)
p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "sample_name", position = position_dodge()) +
    scale_samples_unique +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all_5pct.pdf"), pl = p, w = 12, h = 4)
p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "condition", position = position_dodge(), add = c("mean_se")) +
    geom_pwc(aes(group = condition)) +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all_condition_5pct.pdf"), pl = p, w = 4, h = 4)
p <- dc1 %>%
    filter(CellType == "Neuron") %>%
    ggbarplot(x = "condition", y = "value", fill = "condition", add = c("mean_se", "dotplot")) +
    stat_compare_means() +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all_5pct.pdf"), pl = p, w = 4, h = 4)


dc1 <- dc %>%
    pivot_longer(cols = -CellType) %>%
    dplyr::rename(sample_name = name) %>%
    mutate(sample_name = gsub("\\..*", "", sample_name)) %>%
    filter(value > 0.10) %>%
    left_join(sample_table)
p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "sample_name", position = position_dodge()) +
    scale_samples_unique +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all.pdf"), pl = p, w = 12, h = 4)
p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "condition", position = position_dodge(), add = c("mean_se")) +
    geom_pwc(aes(group = condition)) +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all_condition.pdf"), pl = p, w = 4, h = 4)
p <- dc1 %>%
    filter(CellType == "Neuron") %>%
    ggbarplot(x = "condition", y = "value", fill = "condition", add = c("mean_se", "dotplot")) +
    stat_compare_means() +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all.pdf"), pl = p, w = 4, h = 4)


dc <- read_csv("ldna/results/wgbs/hg38/tables/deconv.brainsubset.csv")

dc1 <- dc %>%
    pivot_longer(cols = -CellType) %>%
    dplyr::rename(sample_name = name) %>%
    mutate(sample_name = gsub("\\..*", "", sample_name)) %>%
    filter(value > 0.10) %>%
    left_join(sample_table)

p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "sample_name", position = position_dodge()) +
    scale_samples_unique +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/brain_sample.pdf"), pl = p, w = 12, h = 4)

p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "condition", position = position_dodge(), add = c("mean_se")) +
    geom_pwc(aes(group = condition)) +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/brain_condition.pdf"), pl = p, w = 4, h = 4)



p <- dc1 %>%
    filter(CellType == "Neuron") %>%
    ggbarplot(x = "condition", y = "value", fill = "condition", add = c("mean_se", "dotplot")) +
    stat_compare_means() +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/brain1.pdf"), pl = p, w = 4, h = 4)
