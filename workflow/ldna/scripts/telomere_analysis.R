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
            telometer_out = sprintf("ldna/results/tables/telometer/%s_telometer_4000.tsv", samples),
            telometer_out_split1 = sprintf("ldna/results/tables/telometer/%s_1_telometer_4000.tsv", samples),
            telometer_out_split2 = sprintf("ldna/results/tables/telometer/%s_2_telometer_4000.tsv", samples),
            telometer_out_split3 = sprintf("ldna/results/tables/telometer/%s_3_telometer_4000.tsv", samples)
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/outfiles/bedmethylanalysis.txt"
        ), env = globalenv())
    }
)

{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
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

mqc <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/aref/qc/multiqc_data/multiqc_general_stats.txt") %>%
    filter(grepl("pyco", Sample)) %>%
    mutate(sample_name = gsub("pycoQC", "", Sample)) %>%
    dplyr::rename(median_read_length = `pycoQC_mqc-generalstats-pycoqc-passed_median_read_length`) %>%
    dplyr::select(sample_name, median_read_length)


dflist <- list()
for (sample_name in samples) {
    tryCatch(
        {
            df <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dflist[[sample_name]] <- df
        },
        error = function(e) {

        }
    )
}

dfsplitlist <- list()
for (sample_name in samples) {
    tryCatch(
        {
            df <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out_split1,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dfsplitlist[[paste0(sample_name, "1")]] <- df

            df <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out_split2,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dfsplitlist[[paste0(sample_name, "2")]] <- df

            df <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out_split3,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dfsplitlist[[paste0(sample_name, "3")]] <- df
        },
        error = function(e) {

        }
    )
}

df <- Reduce(rbind, dflist)
dfsplitlist <- Reduce(rbind, dfsplitlist)

df2 <- rbind(df, dfsplitlist)

duped_reads <- df2 %$% read_id %>% table()
df3 <- df2 %>%
    group_by(read_id) %>%
    arrange(-telomere_length) %>%
    filter(row_number() == 1) %>%
    ungroup()


df %$% sample_name
df1 <- df3 %>%
    mutate(chr_arm = paste0(chromosome, arm)) %>%
    filter(!grepl("^NI", chromosome)) %>%
    mutate(merged_chr = paste0("chr", gsub("_.*", "", gsub("^alt", "", gsub("^chr", "", chromosome))))) %>%
    mutate(merged_chr = factor(merged_chr, levels = refchromosomes, ordered = TRUE))
df1 %>% print(width = Inf)


stat_frame <- df1 %>%
    group_by(merged_chr, sample_name, condition) %>%
    summarize(telomere_length = mean(telomere_length)) %>%
    ungroup()
p <- df1 %>%
    ggboxplot(x = "merged_chr", y = "telomere_length", fill = "condition") +
    geom_pwc(data = stat_frame, aes(group = condition)) +
    labs(x = "", y = "Telomere Length (bp)") +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/all_chr.pdf"), pl = p, w = 30, h = 4)

p <- df1 %>%
    ggboxplot(x = "condition", y = "telomere_length", fill = "condition") +
    labs(x = "", y = "Telomere Length (bp)") +
    stat_compare_means() +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/all.pdf"), pl = p, w = 4, h = 4)

stat_frame2 <- stat_frame %>%
    group_by(merged_chr) %>%
    mutate(chr_mean = mean(telomere_length)) %>%
    ungroup() %>%
    mutate(telomere_length_chrnormed = telomere_length - chr_mean) %>%
    group_by(sample_name, condition) %>%
    summarize(sample_adjusted_mean_telo = mean(telomere_length_chrnormed)) %>%
    left_join(sample_table)
summary(lm(sample_adjusted_mean_telo ~ braak, stat_frame2))

p <- df1 %>%
    group_by(sample_name) %>%
    summarise(ml = mean(telomere_length)) %>%
    arrange(-ml) %>%
    left_join(mqc) %>%
    ggscatter(x = "median_read_length", y = "ml", color = "sample_name") +
    scale_samples_unique
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/telo_vs_medianreadlength.pdf"), pl = p, w = 4, h = 4)





####### Deconv
dc <- read_csv("ldna/results/wgbs/hg38/tables/deconv.csv")
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
