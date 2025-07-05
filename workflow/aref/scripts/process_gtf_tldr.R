module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
set.seed(123)

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(Biostrings)
library(Rsamtools)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            gtf = "aref/default/A.REF_repeatmasker/A.REF_repeatmasker_raw.gtf",
            ref_cytobands = "aref/default/A.REF_annotations/cytobands.bed",
            tldroutput = "aref/A.REF_tldr/A.REF.table.txt",
            ref = "aref/default/A.REF-pre-ins-filtering.fa"
        ), env = globalenv())
        assign("params", list(
            tldr_switch = "process_gtf",
            sample_or_ref = "A.REF"
        ), env = globalenv())
        assign("outputs", list(
            contigs_to_keep = "aref/default/contigs_to_keep.txt",
            filtered_tldr = "aref/default/A.REF.table.kept_in_updated_ref.txt",
            repmask_gff2 = "aref/default/A.REF_annotations/A.REF_repeatmasker.gff2",
            repmask_gff3 = "aref/default/A.REF_annotations/A.REF_repeatmasker.gff3",
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv"
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
rm1 <- rm %>% tidyr::separate(Target, into = c("family", "element_start", "element_end", "element_bp_remaining"), sep = " ")
rm2 <- rm1 %>%
    dplyr::select(-transcript_id) %>%
    mutate(element_bp_remaining = as.numeric(gsub("\\(|\\)", "", element_bp_remaining))) %>%
    mutate(element_start = as.numeric(element_start)) %>%
    mutate(element_end = as.numeric(element_end)) %>%
    mutate(length = element_end - element_start) %>%
    mutate(consensus_length = element_end + element_bp_remaining) %>%
    mutate(pctconsensuscovered = 100 * (length / consensus_length))
rm3 <- rm2 %>%
    mutate(pctdiv = as.numeric(pctdiv)) %>%
    mutate(pctdel = as.numeric(pctdel)) %>%
    mutate(element_start = as.numeric(element_start)) %>%
    mutate(element_end = as.numeric(element_end))


rmfragments <- rm3 %>%
    group_by(gene_id) %>%
    summarise(seqnames = dplyr::first(seqnames), source = dplyr::first(source), type = dplyr::first(type), start = min(start), end = max(end), strand = dplyr::first(strand), phase = dplyr::first(phase), family = dplyr::first(family), element_start = min(element_start), element_end = max(element_end), pctdiv = sum(pctdiv * length) / sum(length), length = sum(length), pctconsensuscovered = sum(pctconsensuscovered), num_fragments = n()) %>%
    mutate(pctconsensustruncated = 100 - pctconsensuscovered)
rmfragments <- rmfragments %>% mutate(refstatus = ifelse(str_detect(seqnames, "^NI_"), "NonRef", "Ref"))

# for annotation purposes, I will have to have the location of nonreference inserts be their insertion site
if (any(str_detect(rmfragments$refstatus, "NonRef"))) {
    rmfragments_ref <- rmfragments %>% filter(refstatus == "Ref")
    rmfragments_nonref <- rmfragments %>% filter(refstatus == "NonRef")
    seqnamesNonref <- rmfragments_nonref %$% seqnames
    insert_seqnames <- c()
    insert_start <- c()
    insert_end <- c()
    for (seqname in seqnamesNonref) {
        seqname_no_rte <- seqname %>% gsub(".*_chr", "chr", .)
        split <- str_split(seqname_no_rte, "_") %>% unlist()
        splitlen <- length(split)
        insert_start <- c(insert_start, split[splitlen - 1])
        insert_end <- c(insert_end, split[splitlen])
        insert_seqnames <- c(insert_seqnames, gsub(paste0("_", split[splitlen - 1], "_", split[splitlen]), "", seqname_no_rte))
    }
    rmfragments_nonref$insert_seqnames <- insert_seqnames
    rmfragments_nonref$insert_start <- insert_start
    rmfragments_nonref$insert_end <- insert_end
    rmfragments_nonref %$% seqnames

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

gr_for_overlap_with_cytobands <- subsetByOverlaps(gr_for_overlap_analysis, cytobands)
gr_for_overlap_analysis_lacking_cytobands <- subsetByOverlaps(gr_for_overlap_analysis, cytobands, invert = TRUE)

overlaps <- findOverlaps(gr_for_overlap_with_cytobands, cytobands, select = "first")
gr_for_overlap_with_cytobands$cytobands <- mcols(cytobands[overlaps])$X4
if (length(gr_for_overlap_analysis_lacking_cytobands) != 0) {
    gr_for_overlap_analysis_lacking_cytobands$cytobands <- "NoCytobandInfo"
    gr_for_overlap_analysis <<- c(gr_for_overlap_with_cytobands, gr_for_overlap_analysis_lacking_cytobands)
} else {
    gr_for_overlap_analysis <<- gr_for_overlap_with_cytobands
}
new_id_mapping_ref <- gr_for_overlap_analysis %>%
    as.data.frame() %>%
    tibble() %>%
    filter(refstatus == "Ref") %>%
    mutate(element_basename = str_split(family, "/") %>% map_chr(tail, 1)) %>%
    mutate(cytobands_good = paste0(str_replace(seqnames, "chr", ""), cytobands)) %>%
    group_by(element_basename, cytobands_good) %>%
    mutate(new_id = paste0(element_basename, "_", cytobands_good, "_", row_number()), n = n()) %>%
    ungroup() %>%
    relocate(gene_id, new_id)

new_id_mapping_nonref <- gr_for_overlap_analysis %>%
    as.data.frame() %>%
    tibble() %>%
    filter(refstatus == "NonRef") %>%
    mutate(element_basename = str_split(family, "/") %>% map_chr(tail, 1)) %>%
    mutate(cytobands_good = paste0(str_replace(seqnames, "chr", ""), cytobands)) %>%
    group_by(element_basename, cytobands_good) %>%
    mutate(new_id = paste0(element_basename, "_", "NI_", params$sample_or_ref, "_", cytobands_good, "_", row_number()), n = n()) %>%
    ungroup() %>%
    relocate(gene_id, new_id)

new_id_mapping <- bind_rows(new_id_mapping_ref, new_id_mapping_nonref) %>%
    dplyr::select(gene_id, new_id) %>%
    dplyr::rename(GiesmaID = new_id)

rmgr <- GRanges(rm %>% full_join(new_id_mapping) %>%
    relocate(GiesmaID, .after = gene_id) %>%
    dplyr::rename(old_id = gene_id) %>%
    dplyr::rename(gene_id = GiesmaID) %>%
    dplyr::relocate(-old_id))

rmfragments <- rmfragments %>%
    full_join(new_id_mapping) %>%
    relocate(GiesmaID, .after = gene_id) %>%
    dplyr::rename(old_id = gene_id) %>%
    dplyr::rename(gene_id = GiesmaID)

if (params$tldr_switch == "process_gtf_tldr") {
    # now determine which nonref contigs should be kept in the reference
    df <- read_delim(inputs$tldroutput) %>%
        mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
        filter(!is.na(Family)) %>%
        filter(!is.na(StartTE)) %>%
        filter(Filter == "PASS")

    er <- df %$% EmptyReads %>% str_extract_all("\\|[0-9]+")
    emptyreadsnum <- c()
    for (i in 1:length(er)) {
        if (length(er[[i]]) == 1) {
            val <- er[[i]] %>%
                gsub("\\|", "", .) %>%
                as.numeric()
        } else {
            val <- er[[i]] %>%
                gsub("\\|", "", .) %>%
                as.numeric() %>%
                sum()
        }
        emptyreadsnum <- c(emptyreadsnum, val)
    }

    df$emptyreadsnum <- emptyreadsnum %>% replace_na(0)
    df <- df %>%
        mutate(fraction_reads_count = UsedReads / (UsedReads + emptyreadsnum))
    df <- df[!(df %$% faName %>% duplicated()), ]

    ss <- DNAStringSet(df %>% dplyr::arrange(faName) %$% Consensus)
    names(ss) <- df %>% dplyr::arrange(faName) %$% faName
    writeXStringSet(ss, paste0(dirname(inputs$ref), "/nonrefcontigs.fa"), append = FALSE, format = "fasta")



    genome_lengths <- fasta.seqlengths(inputs$ref)
    genome_lengths_df <- data.frame(seqnames = names(genome_lengths), contig_length = genome_lengths) %>% tibble()
    refcontigs <- genome_lengths_df %>%
        filter(!grepl("^NI_", seqnames)) %>%
        pull(seqnames)
    nonrefcontigs <- rmfragments %>% filter(grepl("^NI_", seqnames))
    # this will ensure that 1) only elements which are actually picked up by repeat masker (as opposed to 3p transductions with a tiny bit of TE sequence, are retainined in the reference
    # and that 2) I can make a gtf with only those elements which inserted, and not elements which are merely contained in the nonref contig (but which are also present in the reference).
    nonrefelementspass <- nonrefcontigs %>%
        mutate(seqnames_element_type = gsub("_chr.*", "", gsub("NI_", "", seqnames))) %>%
        mutate(seqnames_element_type = gsub("L1Ta", "L1HS", seqnames_element_type)) %>%
        mutate(seqnames_element_type = gsub("L1preTa", "L1HS", seqnames_element_type)) %>%
        mutate(family_element_type = gsub("^.*/", "", family)) %>%
        filter(seqnames_element_type == family_element_type)

    sum(nonrefelementspass %$% seqnames %>% table() > 1)
    contigs_to_keep <- c(refcontigs, nonrefelementspass$seqnames %>% unique() %>% as.character())
    write_delim(as.data.frame(contigs_to_keep), outputs$contigs_to_keep, delim = "\n", col_names = FALSE)

    df_filtered <- df %>%
        filter(faName %in% contigs_to_keep)
    write_delim(df_filtered, outputs$filtered_tldr, delim = "\t")


    # now determine which of the nonref inserts are the bona fide nonref inserts for annotation in the gtf
    rmref <- rmfragments %>% filter(!grepl("^NI_", seqnames))
    rmnonref <- rmfragments %>%
        filter(grepl("^NI_", seqnames)) %>%
        filter(seqnames %in% contigs_to_keep)
    rmnonrefkeep <- nonrefelementspass
    # determine which of these is the bona fide nonref ins, and which are just ref ins which are on a nonref contig
    a <- rmnonrefkeep %>%
        mutate(seqname_ins_type = gsub("_chr.*", "", gsub("NI_", "", seqnames))) %>%
        mutate(seqname_ins_type = gsub("L1Ta", "L1HS", gsub("L1preTa", "L1HS", seqname_ins_type))) %>%
        mutate(ins_type = gsub("^.*/", "", gsub("/NI_.*", "", family))) %>%
        filter(ins_type == seqname_ins_type)

    # bonafide element should be the center of the contig, so select the central element
    a <- a %>%
        left_join(genome_lengths_df) %>%
        mutate(center = contig_length / 2) %>%
        mutate(dist = abs(center - (end - start) / 2))

    rmnonrefkeep_central_element <- a %>%
        left_join(df %>% dplyr::rename(seqnames = faName) %>% dplyr::select(seqnames, Strand)) %>%
        group_by(seqnames) %>%
        filter(dist == min(dist)) %>%
        filter(strand == Strand) %>%
        ungroup() %>%
        dplyr::select(-center, -dist, -Strand, -seqname_ins_type, -ins_type, -contig_length, -seqnames_element_type, -family_element_type) %>%
        left_join(df_filtered %>% dplyr::select(faName, UUID), by = c("seqnames" = "faName")) %>%
        dplyr::rename(nonref_UUID = UUID)


    rmnonref_noncentral_elements <- rmnonref %>%
        anti_join(rmnonrefkeep_central_element) %>%
        mutate(refstatus = "NonCentral")


    rmref$nonref_UUID <- "NotApplicable"
    rmnonref_noncentral_elements$nonref_UUID <- "NotApplicable"
    colnames(rmref)
    colnames(rmnonrefkeep_central_element)
    colnames(rmnonref_noncentral_elements)
    rmfragments <- rbind(rmref, rmnonrefkeep_central_element, rmnonref_noncentral_elements)

    # adding L1 antisense promoter annotations
    if (conf$species %in% c("human")) {
        # adding ASP annotations - first to rmfragments
        fulllength_trnc_length_threshold <- conf$fulllength_trnc_length_threshold
        human_l1_antisense_subfamilies <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")
        ASP <- rmfragments %>%
            filter(grepl(paste(human_l1_antisense_subfamilies, collapse = "|"), gene_id, perl = TRUE)) %>%
            filter(pctconsensuscovered >= fulllength_trnc_length_threshold) %>%
            filter(element_start < 300) %>%
            mutate(gene_id = paste0(gene_id, "__AS")) %>%
            mutate(family = family %>%
                str_split("/") %>%
                map(~ paste0(.x, "__AS")) %>%
                map_chr(~ paste(.x, collapse = "/")))
        ASPpos <- ASP %>%
            filter(strand == "+") %>%
            GRanges() %>%
            promoters(upstream = 0, downstream = 500)
        ASPneg <- ASP %>%
            filter(strand == "-") %>%
            GRanges() %>%
            promoters(upstream = 0, downstream = 500)
        strand(ASPpos) <- "-"
        strand(ASPneg) <- "+"
        ASPfinal <- c(ASPpos, ASPneg) %>%
            as.data.frame() %>%
            tibble() %>%
            dplyr::select(-width)
        rmfragments <- bind_rows(rmfragments, ASPfinal)

        # now for the fragmented
        rmgrdf <- rmgr %>%
            as.data.frame() %>%
            tibble()
        ids <- rmfragments %>% filter(grepl("__AS", gene_id, perl = TRUE)) %$% old_id
        ASPgr <- rmgrdf %>%
            filter(old_id %in% ids) %>%
            mutate(Target2 = Target) %>%
            tidyr::separate(Target2, into = c("family", "element_start", "element_end", "element_bp_remaining"), sep = " ") %>%
            dplyr::select(-c("family", "element_end", "element_bp_remaining")) %>%
            group_by(old_id) %>%
            arrange(element_start) %>%
            filter(row_number() == 1) %>%
            ungroup() %>%
            dplyr::select(-element_start) %>%
            mutate(gene_id = paste0(gene_id, "__AS"))

        ASPposgr <- ASPgr %>%
            filter(strand == "+") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        ASPneggr <- ASPgr %>%
            filter(strand == "-") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        strand(ASPposgr) <- "-"
        strand(ASPneggr) <- "+"
        ASPfinalgr <- c(ASPposgr, ASPneggr) %>%
            as.data.frame() %>%
            tibble() %>%
            dplyr::select(-width)
        rmgr <- c(rmgr, GRanges(ASPfinalgr))
    }

    write_csv(rmfragments %>% dplyr::relocate(-old_id), outputs$r_annotation_fragmentsjoined)
    contigs_to_keep <- rmfragments %$% seqnames %>% unique()
    rtracklayer::export(rmgr[seqnames(rmgr) %in% contigs_to_keep], con = outputs$repmask_gff2, format = "GFF2")
    rtracklayer::export(rmgr[seqnames(rmgr) %in% contigs_to_keep], con = outputs$repmask_gff3, format = "GFF3")
} else {
    # adding L1 antisense promoter annotations
    if (conf$species %in% c("human")) {
        # adding ASP annotations - first to rmfragments
        fulllength_trnc_length_threshold <- conf$fulllength_trnc_length_threshold
        human_l1_antisense_subfamilies <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")
        ASP <- rmfragments %>%
            filter(grepl(paste(human_l1_antisense_subfamilies, collapse = "|"), gene_id, perl = TRUE)) %>%
            filter(pctconsensuscovered >= fulllength_trnc_length_threshold) %>%
            filter(element_start < 300) %>%
            mutate(gene_id = paste0(gene_id, "__AS")) %>%
            mutate(family = family %>%
                str_split("/") %>%
                map(~ paste0(.x, "__AS")) %>%
                map_chr(~ paste(.x, collapse = "/")))

        ASPpos <- ASP %>%
            filter(strand == "+") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        ASPneg <- ASP %>%
            filter(strand == "-") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        strand(ASPpos) <- "-"
        strand(ASPneg) <- "+"
        ASPfinal <- c(ASPpos, ASPneg) %>%
            as.data.frame() %>%
            tibble() %>%
            dplyr::select(-width)
        rmfragments <- bind_rows(rmfragments, ASPfinal)

        # now for the fragmented
        rmgrdf <- rmgr %>%
            as.data.frame() %>%
            tibble()
        ids <- rmfragments %>% filter(grepl("__AS", gene_id, perl = TRUE)) %$% old_id
        ASPgr <- rmgrdf %>%
            filter(old_id %in% ids) %>%
            mutate(Target2 = Target) %>%
            tidyr::separate(Target2, into = c("family", "element_start", "element_end", "element_bp_remaining"), sep = " ") %>%
            dplyr::select(-c("family", "element_end", "element_bp_remaining")) %>%
            group_by(old_id) %>%
            arrange(element_start) %>%
            filter(row_number() == 1) %>%
            ungroup() %>%
            dplyr::select(-element_start) %>%
            mutate(gene_id = paste0(gene_id, "__AS"))

        ASPposgr <- ASPgr %>%
            filter(strand == "+") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        ASPneggr <- ASPgr %>%
            filter(strand == "-") %>%
            GRanges() %>%
            promoters(upstream = 30, downstream = 500)
        strand(ASPposgr) <- "-"
        strand(ASPneggr) <- "+"
        ASPfinalgr <- c(ASPposgr, ASPneggr) %>%
            as.data.frame() %>%
            tibble() %>%
            dplyr::select(-width)
        rmgr <- c(rmgr, GRanges(ASPfinalgr))
    }

    write_csv(rmfragments %>% dplyr::relocate(-old_id), outputs$r_annotation_fragmentsjoined)
    rtracklayer::export(rmgr, con = outputs$repmask_gff2, format = "GFF2")
    rtracklayer::export(rmgr, con = outputs$repmask_gff3, format = "GFF3")
}
