module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
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
library(rBLAST)
library(msa)
library(seqinr)
library(ape)
library(ggtree)
library(treeio)
library(ggmsa)
library(gggenes)
library(ggnewscale)
library(RColorBrewer)
library(ORFik)


tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/default/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/default/A.REF.fa"
        ), env = globalenv())
        assign("params", list(l13 = conf$l13fasta), env = globalenv())
        assign("outputs", list(
            outfile = "aref/default/A.REF_Analysis/l1element_analysis.outfile",
            plots = "aref/default/A.REF_Analysis/l1element_analysis.rds"
        ), env = globalenv())
    }
)


library(Rsamtools)
fa <- Rsamtools::FaFile(inputs$ref)
rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies) %>%
    filter(refstatus != "NonCentral")

rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)

##############################
outputdir <- "/users/mkelsey/data/LF1/RTE/ldna/results/plots/l1_alignment_meth"
outputdir <- "ldna/results/plots/l1_alignment_meth"

subfam <- "L1HS"
grs_fl <- rmann %>%
    filter(rte_length_req == "FL") %>%
    filter(rte_subfamily == subfam) %>%
    group_by(rte_subfamily) %>%
    mutate(length99 = quantile(length, 0.99)) %>%
    ungroup() %>%
    filter(length < length99 * 1.05) %>%
    GRanges()
grs_fl_ss <- getSeq(fa, grs_fl)
names(grs_fl_ss) <- mcols(grs_fl)$gene_id
path <- sprintf("%s/alignments/%s_fl.fa", outputdir, subfam)
dir.create(dirname(path), recursive = TRUE)
writeXStringSet(grs_fl_ss, path)
system(sprintf("echo $(pwd); mafft --auto %s/alignments/%s_fl.fa > %s/alignments/%s_fl.aln.fa", outputdir, subfam, outputdir, subfam))


aln <- readDNAMultipleAlignment(sprintf("%s/alignments/%s_fl.aln.fa", outputdir, subfam))
alnss <- readDNAStringSet(sprintf("%s/alignments/%s_fl.aln.fa", outputdir, subfam))

# Function to get CG positions in alignment and sequence-specific positions
get_cg_positions_df <- function(seq, seq_name) {
    # Alignment sequence as a character string (with gaps)
    seq_str <- as.character(seq)

    # Sequence without gaps for sequence-specific positions
    seq_no_gaps <- gsub("-", "", seq_str)

    # Keep a running index for sequence-specific position (ignores gaps)
    seq_specific_index <- 0

    # Prepare to store positions
    alignment_positions <- numeric()
    sequence_specific_positions <- numeric()

    # Loop through the alignment, tracking both alignment and sequence-specific positions
    for (i in seq(1, nchar(seq_str))) {
        if (substring(seq_str, i, i) != "-") {
            seq_specific_index <- seq_specific_index + 1 # Update sequence-specific index (ignore gaps)
        }

        # Check for CG dinucleotides in the alignment (ignore gaps)
        if (i < nchar(seq_str) && substring(seq_str, i, i + 1) == "CG") {
            alignment_positions <- c(alignment_positions, i)
            sequence_specific_positions <- c(sequence_specific_positions, seq_specific_index)
        }
    }

    # Create a data frame with the alignment and sequence-specific positions
    df <- data.frame(
        Sequence = seq_name,
        Alignment_Position = alignment_positions,
        Sequence_Specific_Position = sequence_specific_positions
    )

    return(df)
}


cg_positions_dfs <- lapply(seq_along(alnss), function(i) {
    get_cg_positions_df(alnss[[i]], names(alnss)[i])
})

# Combine the individual data frames into one tidy data frame
cg_positions_df <- bind_rows(cg_positions_dfs)
write_csv(cg_positions_df, sprintf("%s/%s_fl_cpg_mapping_table.csv", outputdir, subfam))
cg_positions_df %>%
    tibble() %$% Sequence %>%
    table()
df <- cg_positions_df %>%
    tibble() %>%
    mutate(cpgID = paste0("cpg", Alignment_Position))
num_elements <- df %$% Sequence %>%
    unique() %>%
    length()
alnpos_counts <- df %$% Alignment_Position %>% table()

alnpos_keep <- alnpos_counts[alnpos_counts > num_elements / 3] %>% names()
df <- df %>% filter(Alignment_Position %in% alnpos_keep)


methdf <- rtedf %>% filter(rte_subfamily == subfam)
mdf <- methdf %>% mutate(cpg_rel_start = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))


senseelement <- mdf %>%
    filter(rte_strand == "+") %$% gene_id %>%
    pluck(1)
antisenseelement <- mdf %>%
    filter(rte_strand == "-") %$% gene_id %>%
    pluck(1)

cpgmapping_check <- df %>%
    filter(Sequence == senseelement) %$% Sequence_Specific_Position %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == senseelement) %>%
    relocate(cpg_rel_start) %$% cpg_rel_start %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)
cpgmapping_check <- df %>%
    filter(Sequence == antisenseelement) %$% Sequence_Specific_Position %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == antisenseelement) %>%
    relocate(cpg_rel_start) %$% cpg_rel_start %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)

dfr <- df %>% dplyr::rename(gene_id = Sequence, cpg_rel_start = Sequence_Specific_Position)
merged <- left_join(dfr, mdf, by = c("gene_id", "cpg_rel_start"))
cpg_order <- merged %$% Alignment_Position %>%
    unique() %>%
    sort() %>%
    paste0("cpg", .)
merged <- merged %>% mutate(cpgID = factor(cpgID, levels = cpg_order))
library(tidyHeatmap)


p <- merged %>%
    filter(condition == "AD1") %>%
    complete(gene_id, cpgID) %>%
    group_by(gene_id) %>%
    mutate(
        count_NA = sum(is.na(pctM)), # Count of NA values
        count_nonNA = sum(!is.na(pctM)), # Count of non-NA values
        .groups = "drop" # Drop the grouping
    ) %>%
    ungroup() %>%
    filter(count_nonNA > 50) %>%
    group_by(cpgID) %>%
    mutate(
        count_CpG_NA = sum(is.na(pctM)), # Count of NA values
        count_CpG_nonNA = sum(!is.na(pctM)), # Count of non-NA values
        .groups = "drop" # Drop the grouping
    ) %>%
    ungroup() %>%
    filter(count_CpG_nonNA > 50) %>%
    heatmap(gene_id, cpgID, pctM, cluster_rows = TRUE, cluster_columns = FALSE)

p <- merged %>%
    filter(sample == "AD1") %>%
    group_by(gene_id) %>%
    mutate(cpgs_detected_per_element = n()) %>%
    ungroup() %>%
    filter(cpgs_detected_per_element > 50) %>%
    heatmap(gene_id, cpgID, pctM, cluster_rows = TRUE, cluster_columns = FALSE)

dir.create(outputdir, recursive = TRUE)
mysaveandstore(sprintf("%s/pro5.pdf", outputdir), w = 6, h = 6)

hms <- list()
for (sample in conf$samples) {
    p <- merged %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        heatmap(gene_id, cpgID, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
    hms[[sample]] <- p
    dir.create(outputdir, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir, sample, subfam), w = 6, h = 6)
}

elements_of_interest <- c("L1HS_2p13.2_1", "L1HS_2q21.1_2")
rmann %>%
    filter(gene_id == elements_of_interest[2]) %>%
    print(width = Inf)

p <- merged %>%
    filter(condition == "SEN") %>%
    group_by(cpgID) %>%
    mutate(nCPG = n()) %>%
    ungroup() %>%
    filter(nCPG > 50) %>%
    heatmap(gene_id, cpgID, pctM, cluster_rows = TRUE, cluster_columns = FALSE)

outputdir <- "ldna/results/plots/l1_alignment_meth"
dir.create(outputdir, recursive = TRUE)
mysaveandstore(sprintf("%s/sen.pdf", outputdir), w = 6, h = 6)
