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

conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]

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
            df.1 <- read_delim(
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
            df.1 <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out_split1,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dfsplitlist[[paste0(sample_name, "1")]] <- df

            df.1 <- read_delim(
                grep(sprintf("/%s", sample_name),
                    inputs$telometer_out_split2,
                    value = TRUE
                )
            )
            df$sample_name <- sample_name
            df$condition <- sample_table %>% filter(sample_name == !!sample_name) %$% condition
            dfsplitlist[[paste0(sample_name, "2")]] <- df

            df.1 <- read_delim(
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

df.1 <- Reduce(rbind, dflist)
dfsplitlist <- Reduce(rbind, dfsplitlist)

df2 <- rbind(df, dfsplitlist)

duped_reads <- df2 %$% read_id %>% table()
df3 <- df2 %>%
    group_by(read_id) %>%
    arrange(-telomere_length) %>%
    filter(row_number() == 1) %>%
    ungroup()


df.1 %$% sample_name
df <- df3 %>%
    filter(!grepl("^NI", chromosome)) %>%
    filter(mapping_quality > 50) %>%
    mutate(merged_chr = paste0("chr", gsub("_.*", "", gsub("^alt", "", gsub("^chr", "", chromosome))))) %>%
    filter(!(grepl("chrX|chrY", merged_chr))) %>%
    mutate(merged_chr = factor(merged_chr, levels = refchromosomes, ordered = TRUE)) %>%
    mutate(chr_arm = paste0(merged_chr, arm))

df %>% print(width = Inf)
df %$% read_length %>% median()


read_length_filter <- c(4000, 7000, 9000, 11000)
for (filter in read_length_filter) {
    dftemp <- df %>% filter(read_length >= filter)
    stat_frame <- dftemp %>%
        group_by(merged_chr, sample_name, condition) %>%
        summarize(telomere_length = mean(telomere_length)) %>%
        ungroup()
    p <- dftemp %>%
        ggboxplot(x = "merged_chr", y = "telomere_length", fill = "condition") +
        geom_pwc(data = stat_frame, aes(group = condition)) +
        labs(x = "", y = "Telomere Length (bp)") +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_chr.pdf", filter), pl = p, w = 30, h = 4)

    p <- dftemp %>%
        ggboxplot(x = "condition", y = "telomere_length", fill = "condition") +
        labs(x = "", y = "Telomere Length (bp)") +
        stat_compare_means() +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all.pdf", filter), pl = p, w = 4, h = 4)

    stat_frame2 <- stat_frame %>%
        group_by(merged_chr) %>%
        mutate(chr_mean = mean(telomere_length)) %>%
        ungroup() %>%
        mutate(telomere_length_chrnormed = telomere_length - chr_mean) %>%
        group_by(sample_name, condition) %>%
        summarize(sample_adjusted_mean_telo = mean(telomere_length_chrnormed)) %>%
        left_join(sample_table)
    summary(lm(sample_adjusted_mean_telo ~ braak, stat_frame2))

    p <- dftemp %>%
        group_by(sample_name) %>%
        summarise(ml = mean(telomere_length)) %>%
        arrange(-ml) %>%
        left_join(mqc) %>%
        ggscatter(x = "median_read_length", y = "ml", color = "sample_name") +
        scale_samples_unique
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/telo_vs_medianreadlength.pdf", filter), pl = p, w = 4, h = 4)

    dftemp %$% chr_arm %>% table()

    df4 <- dftemp %>%
        group_by(chr_arm) %>%
        mutate(n = n()) %>%
        filter(n > max(length(sample_table$sample_name), 10)) %>%
        mutate(chr_arm_n = paste0(chr_arm, "\nn=", n())) %>%
        ungroup() %>%
        mutate(chr_arm_n = fct_reorder(chr_arm_n, .x = telomere_length, .fun = median))

    p <- df4 %>%
        ggviolin(x = "chr_arm_n", y = "telomere_length", fill = "condition", draw_quantiles = c(.50)) +
        geom_pwc(aes(group = condition)) +
        labs(x = "", y = "Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_chr_arm_violin.pdf", filter), pl = p, w = 30, h = 4)

    p <- df4 %>%
        ggstripchart(x = "chr_arm_n", y = "telomere_length", color = "condition", add = c("median"), add.params = list(color = "black", group = "condition"), position = position_jitterdodge(dodge.width = 1), alpha = 0.85) +
        labs(x = "", y = "Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        lims(y = c(0, 15000)) +
        geom_pwc(aes(group = condition)) +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_chr_arm_strip.pdf", filter), pl = p, w = 30, h = 4)


    p1 <- df4 %>%
        ggstripchart(x = "chr_arm_n", y = "telomere_length", color = "condition", add = c("median"), add.params = list(color = "black", group = "condition"), position = position_jitterdodge(dodge.width = 1), alpha = 0.85) +
        labs(x = "", y = "Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_pwc(aes(group = condition)) +
        scale_conditions +
        mtclosedgridh
    pf <- df4 %>%
        group_by(chr_arm_n, condition) %>%
        summarise(mv = mean(telomere_length)) %>%
        pivot_wider(names_from = condition, values_from = mv) %>%
        mutate(dif = !!sym(condition2) - !!sym(condition1))
    p2 <- ggplot(pf) +
        geom_col(aes(x = chr_arm_n, y = dif, fill = ifelse(dif > 0, TRUE, FALSE)), color = "black") +
        mtclosedgridh
    library(patchwork)
    p <- (p1 / p2) + plot_layout(axes = "collect")
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_chr_arm_strip_stats.pdf", filter), pl = p, w = 30, h = 8)



    df5 <- df4 %>%
        group_by(chr_arm) %>%
        mutate(chr_arm_mean = mean(telomere_length)) %>%
        ungroup() %>%
        mutate(normed_telolength = telomere_length - chr_arm_mean)

    p <- df5 %>%
        ggstripchart(x = "condition", y = "normed_telolength", color = "condition", add = c("median"), add.params = list(color = "black", group = "condition"), position = position_jitterdodge(dodge.width = 1), alpha = 0.85) +
        labs(x = "", y = "Relative Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_pwc(aes(group = condition)) +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_strip_stats.pdf", filter), pl = p, w = 4, h = 4)

    p <- df5 %>%
        ggboxplot(x = "condition", y = "normed_telolength", color = "condition") +
        labs(x = "", y = "Relative Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_pwc(aes(group = condition)) +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_boxplot_stats.pdf", filter), pl = p, w = 4, h = 4)

    p <- df5 %>%
        group_by(sample_name, condition) %>%
        summarise(mean_relative_telo = mean(normed_telolength)) %>%
        ggboxplot(x = "condition", y = "mean_relative_telo", color = "condition") +
        labs(x = "", y = "Relative Telomere Length (bp)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_pwc(aes(group = condition)) +
        scale_conditions +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/plots/telometer/read_filter_%s/all_boxplot_sample_stats.pdf", filter), pl = p, w = 4, h = 4)
    df5
}

df5 %$% telomere_length %>% quantile()
