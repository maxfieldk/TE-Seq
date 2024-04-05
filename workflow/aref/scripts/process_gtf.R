source("workflow/scripts/defaults.R")
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(dtplyr)
library(Biostrings)
library(Rsamtools)

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
            gtf = "aref/repeatmasker/repeatmasker_raw.gtf",
            contigs_to_keep = "aref/contigs_to_keep.txt",
            ref_cytobands = "aref/annotations/cytobands.bed",
            tldroutput = "aref/tldr/tldr.table.txt",
            ref = "aref/ref_pre_ins_filtering.fa"
        ), env = globalenv())
        assign("outputs", list(
            contigs_to_keep = "aref/contigs_to_keep.txt",
            filtered_tldr = "aref/tldr/tldr.table.kept_in_updated_ref.txt",
            repmask_gff2 = "aref/annotations/repeatmasker.gff2",
            repmask_gff3 = "aref/annotations/repeatmasker.gff3",
            r_annotation = "aref/annotations/repeatmasker.gtf.rformatted.csv",
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_annotation_families = "aref/annotations/families_annotation.csv"
        ), env = globalenv())
    }
)

# first get GTF into fragment accounted for format
rmgr <- import(inputs$gtf)
rm <- rmgr %>%
    as.data.frame() %>%
    tibble()

# filter out mitochondrial repeats
rm <- rm %>% filter(seqnames != "chrM")

rm %$% Target %>% head()
rm1 <- rm %>% separate(Target, into = c("family", "element_start", "element_end"), sep = " ")
rm2 <- rm1 %>%
    dplyr::select(-transcript_id) %>%
    mutate(length = as.numeric(element_end) - as.numeric(element_start))
rm3 <- rm2 %>%
    mutate(pctdiv = as.numeric(pctdiv)) %>%
    mutate(pctdel = as.numeric(pctdel)) %>%
    mutate(element_start = as.numeric(element_start)) %>%
    mutate(element_end = as.numeric(element_end))



rmfragments <- rm3 %>%
    group_by(gene_id) %>%
    summarise(seqnames = dplyr::first(seqnames), source = dplyr::first(source), type = dplyr::first(type), start = min(start), end = max(end), strand = dplyr::first(strand), phase = dplyr::first(phase), family = dplyr::first(family), element_start = min(element_start), element_end = max(element_end), pctdiv = sum(pctdiv * length) / sum(length), length = sum(length), num_fragments = n())
rmfragments <- rmfragments %>% mutate(refstatus = ifelse(str_detect(seqnames, "nonref"), "NonRef", "Ref"))

# for annotation purposes, I will have to have the location of nonreference inserts be their insertion site
if (any(str_detect(rmfragments$refstatus, "NonRef"))) {
    rmfragments_ref <- rmfragments %>% filter(refstatus == "Ref")
    rmfragments_nonref <- rmfragments %>% filter(refstatus == "NonRef")
    seqnamesNonref <- rmfragments_nonref %$% seqnames
    insert_seqnames <- c()
    insert_start <- c()
    insert_end <- c()
    for (seqnames in seqnamesNonref) {
        split <- str_split(seqnames, "_") %>% unlist()
        splitlen <- length(split)
        insert_seqnames <- c(insert_seqnames, split[splitlen - 2])
        insert_start <- c(insert_start, split[splitlen - 1])
        insert_end <- c(insert_end, split[splitlen])
    }
    rmfragments_nonref$insert_seqnames <- insert_seqnames
    rmfragments_nonref$insert_start <- insert_start
    rmfragments_nonref$insert_end <- insert_end

    rmfragments_nonrefgr <- GRanges(rmfragments_nonref %>% dplyr::select(-seqnames, -start, -end) %>% dplyr::relocate(gene_id, insert_seqnames, source, type, insert_start, insert_end, strand) %>% dplyr::rename(seqnames = insert_seqnames, start = insert_start, end = insert_end))
    rmfragments_refgr <- GRanges(rmfragments_ref)
    rmfragmentsgr_properinsertloc <- c(rmfragments_refgr, rmfragments_nonrefgr)
    gr_for_overlap_analysis <- rmfragmentsgr_properinsertloc
} else {
    gr_for_overlap_analysis <- GRanges(rmfragments)
}


cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
cytobands <- cytobandsdf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()

overlaps <- findOverlaps(gr_for_overlap_analysis, cytobands, select = "first")
gr_for_overlap_analysis$cytobands <- mcols(cytobands[overlaps])$X4

new_id_mapping <- gr_for_overlap_analysis %>%
    as.data.frame() %>%
    tibble() %>%
    mutate(element_basename = str_split(family, "/") %>% map_chr(tail, 1)) %>%
    mutate(cytobands_good = paste0(str_replace(seqnames, "chr", ""), cytobands)) %>%
    group_by(element_basename, cytobands_good) %>%
    mutate(new_id = paste0(element_basename, "_", cytobands_good, "_", row_number()), n = n()) %>%
    ungroup() %>%
    relocate(gene_id, new_id)


new_id_mapping <- new_id_mapping %>%
    dplyr::select(gene_id, new_id) %>%
    dplyr::rename(GiesmaID = new_id)



rmgr <- GRanges(rm %>% full_join(new_id_mapping) %>% relocate(GiesmaID, .after = gene_id) %>% dplyr::select(-gene_id) %>% dplyr::rename(gene_id = GiesmaID))
# write_csv(rm3 %>% full_join(new_id_mapping) %>% relocate(GiesmaID, .after = gene_id) %>% dplyr::select(-gene_id) %>% dplyr::rename(gene_id = GiesmaID), outputs$r_annotation)

rmfragments <- rmfragments %>%
    full_join(new_id_mapping) %>%
    relocate(GiesmaID, .after = gene_id) %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_id = GiesmaID)



write_csv(rmfragments, outputs$r_annotation_fragmentsjoined)

rtracklayer::export(rmgr[seqnames(rmgr) %in% contigs_to_keep], con = outputs$repmask_gff2, format = "GFF2")
rtracklayer::export(rmgr[seqnames(rmgr) %in% contigs_to_keep], con = outputs$repmask_gff3, format = "GFF3")
