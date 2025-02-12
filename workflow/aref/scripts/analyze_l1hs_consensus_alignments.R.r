library(BSgenome)

module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(jsonlite)
library(ggpubr)
library(ggh4x)
library(seqinr)
library(rBLAST)
library(GenomicAlignments)


tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = if (conf$update_ref_with_tldr$per_sample == "yes") {
                sprintf("aref/%s_tldr/%s.table.txt", conf$samples, conf$samples)
            } else {
                rep(sprintf("aref/%s_tldr/%s.table.txt", "A.REF", "A.REF"), length(sample_table$sample_name))
            },
            blast_njs = if (conf$update_ref_with_tldr$per_sample == "yes") {
                sprintf("aref/default/blastdb/%s.njs", conf$samples)
            } else {
                sprintf("aref/default/blastdb/%s.njs", "A.REF")
            },
            json = sprintf("aref/qc/%s/%spycoQC.json", conf$samples, conf$samples),
            filtered_tldr = if (conf$update_ref_with_tldr$per_sample == "yes") {
                paste0("aref/default/", sample_table$sample_name, ".table.kept_in_updated_ref.txt")
            } else {
                sprintf("aref/default/%s.table.kept_in_updated_ref.txt", "A.REF")
            },
            bam = sprintf("ldna/intermediates/%s/alignments/consensus_seq_alignments/l1hs_%s.bam", sample_table$sample_name, sample_table$sample_name),
            l1hs_aln = sprintf("ldna/intermediates/%s/alignments/consensus_seq_alignments/l1hs_pulled_from_genome_aln_%s.bam", sample_table$sample_name, sample_table$sample_name),
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/default/A.REF_annotations/A.REF_repeatmasker_annotation.csv"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/results/somatic_insertions/analyze_nongermline_insertions.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

library(Rsamtools)
library(GenomicAlignments)

# load sequencing meta data
rm(sample_sequencing_data)
for (sample in sample_table$sample_name) {
    json <- fromJSON(grep(sprintf("/%s/", sample), inputs$json, value = TRUE))
    reads_number <- json[["All Reads"]][["basecall"]]$reads_number
    N50 <- json[["All Reads"]][["basecall"]]$N50
    bases_number <- json[["All Reads"]][["basecall"]]$bases_number
    row <- tibble(sample_name = sample, reads_number = reads_number, N50 = N50, bases_number = bases_number)
    if (!exists("sample_sequencing_data")) {
        sample_sequencing_data <- row
    } else {
        sample_sequencing_data <- rbind(sample_sequencing_data, row)
    }
}


fa <- Rsamtools::FaFile("resources/sequences/L1HS_consensus_polyA_stripped.fa")

flag <- scanBamFlag(
    isUnmappedQuery = NA,
    isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
    isDuplicate = NA, isSupplementaryAlignment = NA
)


for (sample in conf$samples) {
    print(sample)


    bampath <- grep(sample, inputs$bam, value = TRUE)

    aln1 <- readGAlignments(bampath, param = ScanBamParam(what = scanBamWhat(), tag = c("SA"), flag = flag))
    alndf <- as.data.frame(aln1) %>%
        tibble()

    sequenced_l1hs_width <- alndf %$% width %>% sum()
    l1hs_grs <- rtracklayer::import("aref/default/A.REF_annotations/A.REF_rte_beds/L1HS_ALL_ALL.bed")
    total_l1hs_width <- width(l1hs_grs) %>% sum()
    sample_cov <- sample_sequencing_data %>%
        mutate(cov = bases_number / (3.1 * 10**9)) %>%
        filter(sample_name == sample) %$% cov

    sequenced_l1hs_width / (total_l1hs_width * sample_cov)
    excess_bp_in_fll1hs_equivalents <- (sequenced_l1hs_width - total_l1hs_width) / 6000
    numberfl_nonref_needed_to_account <- excess_bp_in_fll1hs_equivalents / sample_cov


    df <- alndf %>%
        filter(width > 1000) %>%
        filter(mapq > 20) %>%
        mutate(last_cigar_element = str_extract(cigar, "[0-9]+[IDHS]$")) %>%
        mutate(last_cigar_element_num = as.numeric(str_extract(str_extract(cigar, "[0-9]+[IDHS]$"), "[0-9]+"))) %>%
        mutate(last_cigar_element_char = str_extract(cigar, "[IDHS]$")) %>%
        mutate(ends_within_100bp = ifelse(last_cigar_element_char == "S" | last_cigar_element_char == "H", ifelse(last_cigar_element_num > 100, 0, 1), 1)) %>%
        mutate(ends_within_50bp = ifelse(last_cigar_element_char == "S" | last_cigar_element_char == "H", ifelse(last_cigar_element_num > 50, 0, 1), 1)) %>%
        mutate(ends_within_10bp = ifelse(last_cigar_element_char == "S" | last_cigar_element_char == "H", ifelse(last_cigar_element_num > 10, 0, 1), 1))

    pf <- df %>%
        mutate(end_bin = cut(end, 20)) %>%
        group_by(end_bin) %>%
        summarise(
            mean_ends_within_10bp = mean(ends_within_10bp, na.rm = TRUE),
            mean_ends_within_50bp = mean(ends_within_50bp, na.rm = TRUE),
            mean_ends_within_100bp = mean(ends_within_100bp, na.rm = TRUE),
        )

    p <- pf %>% ggplot() +
        geom_point(aes(x = end_bin, y = mean_ends_within_50bp))
    mysaveandstore()
}
