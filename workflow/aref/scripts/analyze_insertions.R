source("~/data/common/myDefaults.r")
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

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")[["aref"]]
)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = "aref/tldr/tldr.table.txt",
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/annotations/repeatmasker_annotation.csv",
            ref = "aref/ref_pre_ins_filtering.fa"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/RefAnalysis/tldr_plots/tldr_plots.rds",
            contigs_to_keep = "aref/contigs_to_keep.txt",
            filtered_tldr = "aref/tldr/tldr.table.kept_in_updated_ref.txt",
            updated_reference
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
plots <- list()


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

df <- read_delim(inputs$tldroutput) %>%
    mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS")


ss <- DNAStringSet(df %>% dplyr::arrange(faName) %$% Consensus)
names(ss) <- df %>% dplyr::arrange(faName) %$% faName
writeXStringSet(ss, paste0(dirname(inputs$ref), "/nonrefcontigs.fa"), append = FALSE, format = "fasta")



genome_lengths <- fasta.seqlengths(inputs$ref)
genome_lengths_df <- data.frame(seqnames = names(genome_lengths), contig_length = genome_lengths) %>% tibble()
refcontigs <- genome_lengths_df %>%
    filter(!grepl("nonref", seqnames)) %>%
    pull(seqnames)
nonrefcontigs <- rmann %>% filter(grepl("nonref", seqnames))
# this will ensure that 1) only elements which are actually picked up by repeat masker (as opposed to 3p transductions with a tiny bit of TE sequence, are retainined in the reference
# and that 2) I can make a gtf with only those elements which inserted, and not elements which are merely contained in the nonref contig (but which are also present in the reference).
nonrefelementspass <- nonrefcontigs %>%
    mutate(seqnames_element_type = gsub("_chr.*", "", gsub("nonrefins_", "", seqnames))) %>%
    mutate(seqnames_element_type = gsub("L1Ta", "L1HS", seqnames_element_type)) %>%
    mutate(seqnames_element_type = gsub("L1preTa", "L1HS", seqnames_element_type)) %>%
    mutate(family_element_type = gsub("^.*/", "", family)) %>%
    filter(seqnames_element_type == family_element_type)

sum(nonrefelementspass %$% seqnames %>% table() > 1)



contigs_to_keep <- c(refcontigs, nonrefelementspass$seqnames %>% unique())
write_delim(as.data.frame(contigs_to_keep), outputs$contigs_to_keep, delim = "\n", col_names = FALSE)

df_filtered <- df %>%
    filter(faName %in% contigs_to_keep)
write_delim(df_filtered, outputs$filtered_tldr, delim = "\t")


dfgrs <- GRanges(seqnames = df$Chrom, ranges = IRanges(start = df$Start, end = df$End), strand = df$Strand, UUID = df$UUID)

df %$% Subfamily %>% table()
df %$% Filter %>% table()
p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysave(sprintf("%s/insertions.png", outputdir), 5, 5)
plots[["insertions"]] <- p


p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_wrap(~Family, scales = "free") +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysave(sprintf("%s/insertions_subfamily.png", outputdir), 5, 5)
plots[["insertions_subfamily"]] <- p

p <- df %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mythemeconditions
mysave(sprintf("%s/l1hs_length.png", outputdir), 5, 5)
plots[["l1hs_length"]] <- p


# only insertions that are in the updated reference
df_filtered %$% Subfamily %>% table()
df_filtered %$% Filter %>% table()
p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Family)) +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysave(sprintf("%s/insertions_in_updated_ref.png", outputdir), 5, 5)
plots[["insertions_in_updated_ref"]] <- p


p <- df_filtered %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_wrap(~Family, scales = "free") +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    scale_fill_manual(values = (my_palette)) +
    ggtitle("Non-reference RTE Insertions")
mysave(sprintf("%s/insertions_subfamily_in_updated_ref.png", outputdir), 6, 4)
plots[["insertions_subfamily_in_updated_ref"]] <- p

p <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mythemeconditions
mysave(sprintf("%s/l1hs_length_in_updated_ref.png", outputdir), 5, 5)
plots[["l1hs_length"]] <- p

###############

trsd <- df %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd %>% head(n = 4) %$% Transduction_3p

trsd %$% TransductionLen %>% summary()

save(plots, file = outputs$plots)
