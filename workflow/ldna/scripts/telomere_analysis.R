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
    tryCatch({
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
    })

}
df <- Reduce(rbind, dflist)

df1 <- df %>%
    mutate(chr_arm = paste0(chromosome, arm)) %>%
    filter(!grepl("^NI", chromosome)) %>%
    mutate(merged_chr = paste0("chr",gsub("_.*", "", gsub("^alt", "", gsub("^chr", "", chromosome))))) %>%
    mutate(merged_chr = factor(merged_chr, levels = refchromosomes, ordered = TRUE))
df1 %>% print(width = Inf)

p <- df1 %>%
    ggboxplot(x = "merged_chr", y = "telomere_length", fill = "condition") +
    geom_pwc(aes(group = condition)) +
    labs(x = "", y = "Telomere Length (bp)") +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/all_chr.pdf"), pl = p, w = 30, h = 4)

p <- df1 %>%
    ggboxplot(x = "condition", y = "telomere_length", fill = "condition") +
    labs(x = "", y = "Telomere Length (bp)") +
    stat_compare_means()+
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/all.pdf"), pl = p, w = 4, h = 4)


p <- df1 %>% group_by(sample_name) %>% summarise(ml = mean(telomere_length)) %>% arrange(-ml) %>% left_join(mqc) %>% 
    ggscatter(x = "median_read_length", y = "ml", color = "sample_name") +
    scale_samples_unique
mysaveandstore(fn = sprintf("ldna/results/plots/telometer/telo_vs_medianreadlength.pdf"), pl = p, w = 4, h = 4)





####### Deconv
dc <- read_csv("ldna/results/wgbs/hg38/tables/deconv.csv")

dc1 <- dc %>% pivot_longer(cols = -CellType) %>% dplyr::rename(sample_name = name) %>%
    mutate(sample_name = gsub("\\..*", "", sample_name)) %>% filter(value > 0.10) %>% 
    left_join(sample_table)

p <- dc1 %>% ggbarplot(x = "CellType", y = "value", fill = "sample_name", position = position_dodge()) +
    scale_samples_unique +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all.pdf"), pl = p, w = 12, h = 4)


p <- dc1 %>% filter(CellType == "Neuron") %>%
    ggbarplot(x = "condition", y = "value", fill = "condition", add = c("mean_se", "dotplot")) +
    stat_compare_means() +
    scale_conditions +
    mtclosedgridh
mysaveandstore(fn = sprintf("ldna/results/plots/deconvolve/all.pdf"), pl = p, w = 4, h = 4)


dc1 %>% stat_compare_means()



