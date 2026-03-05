module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

library("readr")
library("ggplot2")
library("tibble")
library("RColorBrewer")
library("cowplot")

library("tidyr")
library(stringr)
library("dplyr")
library(ComplexHeatmap)
library(patchwork)
library(rtracklayer)
library(GenomicFeatures)

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "annotation_genes" = conf$annotation_genes,
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("outputs", list(
            tpm = "srna/outs/agg/tpm/tpmdf.tsv"
        ), env = globalenv())
        assign("inputs", list(
            counts = sprintf("%s/outs/agg/featurecounts_genes/counts.txt", conf$prefix),
            rte_counts = sprintf("%s/outs/%s/telescope/telescope-run_stats.tsv", conf$prefix, conf$samples)
        ), env = globalenv())
    }
)

# https://support.bioconductor.org/p/91218/
tpm <- function(counts, len) {
    x <- counts / len
    return(t(t(x) * 1e6 / colSums(x)))
}


counttype <- params[["counttype"]]

print(counttype)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
coldata <- read.csv(params[["sample_table"]])
samples <- conf$samples
coldata <- coldata[match(conf$samples, coldata$sample_name), ]

# repeats
rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE
)
repeat_lengths <- rmann %>% dplyr::select(gene_id, length)
repeat_lengths %>% filter(length == 0)

if (counttype == "telescope_multi") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- gsub("run_stats", "TE_counts", grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE))
        tdf <- read_delim(path, comment = "#", col_names = TRUE) %>%
            dplyr::rename(gene_id = transcript)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

if (counttype == "telescope_unique") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE)
        tdf <- read_delim(path, comment = "#", col_names = FALSE) %>%
            dplyr::select(X1, X6) %>%
            dplyr::rename(gene_id = X1)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]

refseq <- import(params$annotation_genes)
refseqdf <- refseq %>%
    as.data.frame() %>%
    tibble()
refseqexons <- refseq[refseq$type == "exon"]
refseqgrl <- split(refseqexons, refseqexons$transcript_id)
gene_length_vec <- width(refseqgrl) %>% sum()
gene_lengths1 <- tibble(transcript_id = names(gene_length_vec), length = unlist(gene_length_vec))
gene_lengths2 <- gene_lengths1 %>% left_join(refseqdf %>% dplyr::select(gene_id, transcript_id) %>% distinct())
gene_lengths <- gene_lengths2 %>%
    group_by(gene_id) %>%
    summarize(length = median(length))

gene_cts <- read.delim(inputs$counts)
str(gene_cts)
length(gene_cts$gene_id)
cnames <- colnames(gene_cts)

snames <- c()
cvalues <- c()
for (sample in conf$samples) {
    colname_for_sample <- grep(paste0("outs.", sample, ".star"), cnames, value = TRUE)
    snames <- c(snames, sample)
    cvalues <- c(cvalues, colname_for_sample)
}
gene_cts <- gene_cts %>%
    dplyr::rename(!!!setNames(cvalues, snames)) %>%
    dplyr::rename(gene_id = Geneid)
gene_cts <- gene_cts[c("gene_id", sample_table$sample_name)]



# merge genes and repeats
cts <- rbind(gene_cts, as.data.frame(bounddf1 %>% replace(is.na(.), 0)))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

lengths <- rbind(gene_lengths, repeat_lengths)
len_df <- tibble(gene_id = rownames(cts)) %>% left_join(lengths)

filter_out <- len_df %>% filter(if_any(everything(), is.na) | length == 0) %$% gene_id


cts1 <- cts[!(rownames(cts) %in% filter_out), ]
len_df1 <- len_df %>% filter(!(gene_id %in% filter_out))
tpms <- tpm(cts1, len_df1 %$% length)
colSums(tpms)
tpmdf <- tpms %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    tibble() %>%
    relocate(gene_id)
tpmdf %>% filter(grepl("_NI", gene_id))
outpath <- outputs$tpm
dir.create(dirname(outpath), recursive = TRUE)
write_delim(tpmdf, outpath)
