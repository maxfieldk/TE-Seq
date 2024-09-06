module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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
# library(ggsankey)
library(circlize)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            filtered_tldr = "aref/A.REF_tldr/A.REF.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/default/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/default/A.REF.fa",
            blast_njs = "aref/default/A.REF.njs"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/default/A.REF_Analysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/default/A.REF_Analysis/transduction_df.csv"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "A.REF"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)
rmann %>% filter(refstatus == "NonRef")
df_filtered <- read_delim(inputs$filtered_tldr)
l1ta <- df_filtered %>%
    filter(Family %in% c("L1")) %>% 
    filter(Filter == "PASS")
l1ta_sankey <- l1ta %$% UUID



#TRANSDUCTION BASED MAPPING
tryCatch({

trsd <- l1ta %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 20)
trsd_sankey <- trsd %$% UUID
trsd %>% head(n = 4) %$% Transduction_3p


trsd_ss <- DNAStringSet(trsd$Transduction_3p)
names(trsd_ss) <- trsd$UUID


at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)

trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]
trsd_ss_for_blast_sankey <-names(trsd_ss_for_blast)

bl <- blast(db = sub('\\.[^.]*$', '', inputs$blast_njs))
bres <- tibble(predict(bl, trsd_ss_for_blast)) 

bres1 <- bres %>% left_join(trsd, by = c("qseqid" = "UUID"))

bres %$% qseqid %>%
    unique() %>%
    length()

bres_hits <- bres1 %>%
    filter(!grepl("NI", sseqid)) %>%
    group_by(qseqid) %>%
    filter(bitscore == max(bitscore)) %>%
    filter(pident == max(pident)) %>%
    filter(gapopen == min(gapopen)) %>%
    filter(length == max(length)) %>% filter(nrow(.) == 1) %>% ungroup()

join_df <- rmann %>% filter(refstatus == "NonRef") %>% dplyr::select(gene_id, old_id) %>% mutate(Chrom = gsub(".*chr", "chr", old_id) %>% gsub("_.*", "", .)) %>% mutate(Start = gsub(".*chr", "chr", old_id) %>% str_split_i("_", 2) %>% as.numeric()) %>% dplyr::select(-old_id)
bres_hits1 <- bres_hits %>% left_join(join_df)
relevant_subfamilies <- bres_hits1 %$% Subfamily %>% unique()
l1grs <- GRanges(rmann %>% filter(grepl(paste(relevant_subfamilies, sep = "|", colapse = TRUE), rte_subfamily)))

bres_hits_grs_prep <- bres_hits1 %>%
    dplyr::select(gene_id, sseqid, sstart, send) %>%
    group_by(gene_id) %>%
    mutate(seqnames = sseqid) %>%
    mutate(start = min(sstart, send), end = max(sstart, send)) %>%
    dplyr::select(seqnames, start, end, gene_id)
bresgrs <- GRanges(bres_hits_grs_prep)
# extend by 500 bp on either side
bresgrs <- resize(bresgrs, width = width(bresgrs) + 20000, fix = "center")

trsd_adjacent_l1 <- mergeByOverlaps(bresgrs, l1grs) %>% as.data.frame() %>% tibble() %>% dplyr::rename(from_gene_id = l1grs.gene_id, to_gene_id = bresgrs.gene_id) %>% dplyr::select(from_gene_id, to_gene_id)

from_df <- trsd_adjacent_l1 %>% dplyr::select(gene_id = from_gene_id) %>% left_join(rmann) %>% dplyr::select(gene_id, seqnames, start, end)
to_df <- trsd_adjacent_l1 %>%
    dplyr::select(gene_id = to_gene_id) %>%
    left_join(rmann ) %>% dplyr::select(gene_id, old_id) %>%
        mutate(seqnames = gsub(".*chr", "chr", old_id) %>%
            gsub("_.*", "", .)) %>%
        mutate(start = gsub(".*chr", "chr", old_id) %>%
            str_split_i("_", 2) %>% as.numeric()) %>% 
        mutate(end = gsub(".*chr", "chr", old_id) %>% 
            str_split_i("_", 3) %>% as.numeric()) %>%
        dplyr::select(gene_id, seqnames, start, end)

from_elements <- from_df %>% pull(gene_id) %>% unique()
from_elements_colors <- tibble(gene_id = from_elements, color = adjust_transparency(as.character(paletteer_c("grDevices::Roma", direction = 1, n = length(from_elements))), alpha = 1))
from_df <- left_join(from_df, from_elements_colors, by = "gene_id")
cytobands <- read.cytoband(
    cytoband = ,
    species = NULL,
    sort.chr = TRUE)

{
pdf(sprintf("%s/transduction_mapping.pdf", outputdir), width = 5, height = 5)

circos.initializeWithIdeogram(
    cytoband = conf$ref_cytobands,
    sort.chr = TRUE,
    major.by = NULL,
    plotType = c("ideogram", "labels"),
    track.height = NULL,
    ideogram.height = convert_height(5, "mm"))

circos.genomicLink(from_df %>% relocate(-gene_id), to_df%>% relocate(-gene_id),col = from_df$color, lwd = 2, direction = 1)
title(sprintf("%s Source Element Mapping", paste(relevant_subfamilies, sep = "|", colapse = TRUE)))
dev.off()
}

number_of_offspring <- from_df %>% group_by(gene_id) %>% summarize(n = n())
p <- number_of_offspring  %>% ggplot(aes(x = fct_reorder(gene_id, n), y = n)) + geom_col() + mtopen + coord_flip() +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    labs(x = "Source Element", y = "Number of Offspring", caption = "Transduction-based Source Element Mapping") +
    anchorbar
mysaveandstore(sprintf("%s/Transduction_based_offspring_per_source_element.pdf", outputdir), 3, 4)
}, error = function(e) {
    print("error in transduction mapping. There may be no non-reference elements")
    print(e)
    
})

#PHYLOGENY BASED MAPPING OF INTACT ELEMENTS
for (alignment in list.files(sprintf("aref/default/%s_Analysis", params$sample_or_ref), pattern = "_intact.aln.fa", full.names = TRUE)) {
    tryCatch({
    element_name <- gsub("_intact.aln.fa", "", alignment) %>% gsub(".*/", "", .)

    aln <- read.alignment(alignment, format = "fasta")
    d <- dist.alignment(aln, "identity", gap = FALSE)
    d_gapped <- dist.alignment(aln, "identity", gap = TRUE)

    # this is the pct difference in identity between sequences
    dsqaured <- 100 * d**2
    d_gappedsquared <- 100 * d_gapped**2

    str(dsqaured)
    tree_raw <- as.treedata(fastme.bal(dsqaured))
    treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
    # to df for further modifications
    ttibble_raw <- as_tibble(tree_raw)
    ttibblegapped_raw <- as_tibble(treegapped_raw)

    
    rmannl1 <- rmann %>% filter(grepl(element_name, rte_subfamily))

    non_ref_tree <- as.phylo(ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")))

    non_ref_nodes <- ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% filter(refstatus == "NonRef") %>% pull(node)
    ref_nodes <- ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% filter(refstatus == "Ref") %>% pull(node)
    treedf <- ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% filter(!is.na(refstatus))
    nearest_tip_per_tip <- castor::find_nearest_tips(as.phylo(ttibble_raw), target_tips = ref_nodes)$nearest_tip_per_tip
    treedf$nearest_ref_tip <- nearest_tip_per_tip

    treedf_nonref <- treedf %>% filter(refstatus == "NonRef") %>% select(node, nearest_ref_tip) %>% dplyr::rename(from = nearest_ref_tip, to = node)
    #initialize new data frame
    from_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
    to_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
    for (index in 1:nrow(treedf_nonref)) {
            from <- treedf_nonref[index,]$from
            to <- treedf_nonref[index,]$to
            from_gene_id <- treedf %>% filter(node == from) %>% pull(label)
            to_gene_id <- treedf %>% filter(node == to) %>% pull(label)
            from_chr <- rmannl1 %>% filter(gene_id == from_gene_id) %>% pull(seqnames)
            to_chr <- rmannl1 %>% filter(gene_id == to_gene_id) %>% pull(seqnames) %>% str_extract("chr[0-9XY]+")
            from_start <- rmannl1 %>% filter(gene_id == from_gene_id) %>% pull(start)
            to_start <- rmannl1 %>% filter(gene_id == to_gene_id) %>% pull(start)
            from_end <- rmannl1 %>% filter(gene_id == from_gene_id) %>% pull(end)
            to_end <- rmannl1 %>% filter(gene_id == to_gene_id) %>% pull(end)
            from_df_row <- tibble(gene_id = from_gene_id, chr = from_chr, start = from_start, end = from_end)
            to_df_row <- tibble(gene_id = to_gene_id, chr = to_chr, start = to_start, end = to_end)
            from_df <- bind_rows(from_df, from_df_row)
            to_df <- bind_rows(to_df, to_df_row)
    }

    from_elements <- from_df %>% pull(gene_id) %>% unique()
    from_elements_colors <- tibble(gene_id = from_elements, color = adjust_transparency(as.character(paletteer_c("grDevices::Roma", direction = 1, n = length(from_elements))), alpha = 1))
    from_df <- left_join(from_df, from_elements_colors, by = "gene_id")
    cytobands <- read.cytoband(
        cytoband = ,
        species = NULL,
        sort.chr = TRUE)

    {
    pdf(sprintf("%s/%s_phylogeny_mapping.pdf", outputdir, element_name), width = 5, height = 5)

    circos.initializeWithIdeogram(
        cytoband = conf$ref_cytobands,
        sort.chr = TRUE,
        major.by = NULL,
        plotType = c("ideogram", "labels"),
        track.height = NULL,
        ideogram.height = convert_height(5, "mm"))

    circos.genomicLink(from_df, to_df, col = from_df$color, lwd = 2, direction = 1)
    title("element_name Source Element Mapping")
    dev.off()
    }



    number_of_offspring <- from_df %>% group_by(gene_id) %>% summarize(n = n())
    p <- number_of_offspring  %>% ggplot(aes(x = fct_reorder(gene_id, n), y = n)) + geom_col() + mtopen + coord_flip() +
        scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
        labs(x = "Source Element", y = "Number of Offspring", caption = "Phylogeny-based Source Element Mapping") +
        anchorbar
    mysaveandstore(sprintf("%s/%s_offspring_per_source_element.pdf", outputdir, element_name), 4, 6) 
}, error = function(e) {
    print("error in transduction mapping. There may be no non-reference elements")
    print(e)
    
})
}

x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)
