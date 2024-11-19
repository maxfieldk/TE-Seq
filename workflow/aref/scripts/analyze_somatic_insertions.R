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
                rep(sprintf("aref/default/%s.table.kept_in_updated_ref.txt", "A.REF"), length(sample_table$sample_name))
            },
            bam = sprintf("aref/intermediates/%s/alignments/5khz/%s.hac.5mCG_5hmCG.sorted.bam", sample_table$sample_name, sample_table$sample_name),
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

# load genome annotation data
rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)


segdupsdf <- read_delim("/users/mkelsey/data/Nanopore/alz/segdups.bed", col_names = TRUE, delim = "\t")
segdups <- segdupsdf %>%
    dplyr::rename(seqnames = `#chrom`, start = chromStart, end = chromEnd) %>%
    GRanges()
cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
cytobands <- cytobandsdf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
centromere <- cytobandsdf %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
telomeredf <- read_delim(conf$ref_telomere, col_names = FALSE, delim = "\t")
telomere <- telomeredf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
mappability <- import("/users/mkelsey/data/Nanopore/alz/k50.Unique.Mappability.bb")

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

# load tldr germinline tldr insertions that wind up in reference
dfs_filtered <- list()
if (conf$update_ref_with_tldr$per_sample == "yes") {
    for (sample in sample_table$sample_name) {
        df <- read.table(grep(sprintf("%s.table", sample), inputs$filtered_tldr, value = TRUE), header = TRUE)
        df$sample_name <- sample
        df <- df %>% left_join(sample_table)
        dfs_filtered[[sample]] <- df
        germline <- do.call(rbind, dfs_filtered) %>%
            tibble() %>%
            left_join(sample_sequencing_data)
    }
} else {
    germline <- read.table(inputs$filtered_tldr, header = TRUE) %>%
        tibble()
}
germline_insert_df <- GRanges(germline) %>%
    as.data.frame() %>%
    tibble()




# load all TLDR insertions
library(BSgenome)
fa <- Rsamtools::FaFile(conf$ref)
dflist <- list()
if (conf$update_ref_with_tldr$per_sample == "yes") {
    for (sample in sample_table$sample_name) {
        df <- read.table(grep(sprintf("%s_tldr", sample), inputs$tldroutput, value = TRUE), header = TRUE) %>%
            mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
            tibble()

        df$sample_name <- sample
        df <- df %>%
            filter(grepl(sample, SampleReads)) %>%
            filter(!grepl(paste(setdiff(sample_table$sample_name, sample), collapse = "|"), SampleReads))
        dflist[[sample]] <- df
    }
    dfall <- do.call(bind_rows, dflist) %>%
        left_join(sample_sequencing_data) %>%
        left_join(sample_table)
    somatic <- dfall %>%
        separate_wider_delim(EmptyReads, delim = "|", names = c("bamname", "emptyreadsnum")) %>%
        mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>%
        filter(fraction_reads_count < 0.1) %>%
        filter(MedianMapQ >= 60) %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(!is.na(TSD))

    somatic_all <- dfall %>%
        separate_wider_delim(EmptyReads, delim = "|", names = c("bamname", "emptyreadsnum")) %>%
        mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>%
        filter(fraction_reads_count < 0.1) %>%
        filter(MedianMapQ >= 60) %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1)
} else {
    dfall <- read.table("aref/A.REF_tldr/A.REF.table.txt", header = TRUE) %>%
        mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
        tibble()

    somatic1 <- dfall %>%
        filter(!is.na(EmptyReads)) %>%
        rowwise() %>%
        mutate(emptyreadsnum = sum(as.numeric(gsub("\\|", "", unlist(str_extract_all(EmptyReads, pattern = "\\|[0-9]+")))))) %>%
        ungroup() %>%
        mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>%
        filter(fraction_reads_count < 0.1) %>%
        filter(MedianMapQ >= 60) %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(!is.na(TSD)) %>%
        mutate(sample_name = str_extract(SampleReads, paste(conf$samples, collapse = "|")))

    somatic_all1 <- dfall %>%
        filter(!is.na(EmptyReads)) %>%
        rowwise() %>%
        mutate(emptyreadsnum = sum(as.numeric(gsub("\\|", "", unlist(str_extract_all(EmptyReads, pattern = "\\|[0-9]+")))))) %>%
        ungroup() %>%
        mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>%
        filter(fraction_reads_count < 0.1) %>%
        filter(MedianMapQ >= 60) %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        mutate(sample_name = str_extract(SampleReads, paste(conf$samples, collapse = "|")))
}

somatic %$% sample_name
somatic_all %$% sample_name

somatic_all_not_read_filtered <- GRanges(somatic_all1) %>%
    subsetByOverlaps(centromere, invert = TRUE) %>%
    subsetByOverlaps(telomere, invert = TRUE) %>%
    subsetByOverlaps(segdups, invert = TRUE) %>%
    subsetByOverlaps(mappability) %>%
    as.data.frame() %>%
    tibble()
somatic_not_read_filtered <- GRanges(somatic1) %>%
    subsetByOverlaps(centromere, invert = TRUE) %>%
    subsetByOverlaps(telomere, invert = TRUE) %>%
    subsetByOverlaps(segdups, invert = TRUE) %>%
    subsetByOverlaps(mappability) %>%
    as.data.frame() %>%
    tibble()

somatic_repregion_not_read_filtered <- GRanges(somatic1) %>%
    subsetByOverlaps(c(centromere, telomere, segdups)) %>%
    subsetByOverlaps(mappability) %>%
    as.data.frame() %>%
    tibble()


filter_by_read_metadata <- function(insertdf) {
    flag <- scanBamFlag(
        isUnmappedQuery = NA,
        isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
        isDuplicate = NA, isSupplementaryAlignment = NA
    )

    insert_mean_mapqs <- list()
    insert_supplementary_status <- list()
    insert_indel_and_clipping_status <- list()
    for (sample in conf$samples) {
        print(sample)
        bampath <- grep(sample, inputs$bam, value = TRUE)
        insert_ids <- insertdf %>% filter(sample_name == sample) %$% UUID
        print(length(insert_ids))
        for (insert_id in insert_ids) {
            insert_row <- insertdf %>%
                filter(sample_name == sample) %>%
                filter(UUID == insert_id)
            whichgr <- GRanges(insertdf %>% filter(sample_name == sample) %>% filter(UUID == insert_id))
            aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
            alndf <- as.data.frame(aln1) %>%
                tibble()
            mean_mapq <- alndf$mapq %>% mean()
            insert_mean_mapqs[[insert_id]] <- mean_mapq

            # insert fails if it was captured in a read which has supplementary alignments
            if (conf$update_ref_with_tldr$per_sample == "yes") {
                tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
            } else {
                tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", "A.REF", "A.REF", insert_id))
            }
            read_name <- tldr_df %>% filter(Useable == TRUE & IsSpanRead == TRUE) %$% ReadName
            read_of_interest <- alndf %>% filter(qname == read_name)
            insert_supplementary_status[[insert_id]] <- ifelse(is.na(read_of_interest$SA), "PASS", "FAIL")

            # insert fails if it was captured in a read which has extensive indels or clipping
            inserts <- str_extract_all(read_of_interest$cigar, "[0-9]+I") %>%
                unlist() %>%
                gsub("I", "", .) %>%
                as.numeric()
            deletions <- str_extract_all(read_of_interest$cigar, "[0-9]+D") %>%
                unlist() %>%
                gsub("D", "", .) %>%
                as.numeric()
            clips <- str_extract_all(read_of_interest$cigar, "([0-9]+S)|([0-9]+H)") %>%
                unlist() %>%
                gsub("S|H", "", .) %>%
                as.numeric()
            insert_indel_and_clipping_status[[insert_id]] <- ifelse((sum(inserts > 50) > 1) | (sum(deletions > 50) > 0) | (sum(clips > 50) > 0), "FAIL", "PASS")
        }
    }
    mean_mapq_filter <- tibble(UUID = names(insert_mean_mapqs), insert_mean_mapqs = unname(insert_mean_mapqs) %>% unlist())
    supplementary_alignment_filter <- tibble(UUID = names(insert_supplementary_status), insert_supplementary_status = unname(insert_supplementary_status) %>% map(~ .x[1]) %>% unlist())
    indel_and_clipping_filter <- tibble(UUID = names(insert_indel_and_clipping_status), insert_indel_and_clipping_status = unname(insert_indel_and_clipping_status) %>% unlist())

    insertdf <- insertdf %>%
        left_join(mean_mapq_filter) %>%
        left_join(supplementary_alignment_filter) %>%
        left_join(indel_and_clipping_filter)
    insertdf <- insertdf %>%
        filter(insert_mean_mapqs > 55) %>%
        filter(insert_supplementary_status == "PASS") %>%
        filter(insert_indel_and_clipping_status == "PASS")
    return(insertdf)
}

filter_by_teend <- function(insertdf) {
    insertdf <- insertdf %>%
        mutate(endte_manual_filter = case_when(
            Family == "ALU" ~ 250,
            Family == "L1" ~ 5800,
        )) %>%
        filter(EndTE >= endte_manual_filter | is.na(endte_manual_filter)) %>%
        return(insertdf)
}

somatic <- filter_by_read_metadata(somatic_not_read_filtered)
somatic %>%
    write_delim(sprintf("%s/somatic_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")
somatic <- read_delim(sprintf("%s/somatic_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")

somatic_filt <- filter_by_teend(somatic)
somatic_filt %>%
    write_delim(sprintf("%s/somatic_filt_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")
somatic_filt <- read_delim(sprintf("%s/somatic_filt_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")

somatic_repregion <- filter_by_read_metadata(somatic_repregion_not_read_filtered)
somatic_repregion %>%
    write_delim(sprintf("%s/somatic_repregion_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")
somatic_repregion <- read_delim(sprintf("%s/somatic_repregion_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")

somatic_repregion_filt <- filter_by_teend(somatic_repregion)
somatic_repregion_filt %>%
    write_delim(sprintf("%s/somatic_repregion_filt_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")
somatic_repregion_filt <- read_delim(sprintf("%s/somatic_repregion_filt_inserts_read_level_filters_passed.tsv", outputdir), delim = "\t")


# ADDTIONAL FILTERS
# I think these are overly restrictive - not applying these
germline_insert_characteristics <- tibble(germline_tsd_95 = numeric(), germline_trsd3_95 = numeric(), germline_trsd5_95 = numeric(), germline_endte_05 = numeric(), Subfamily = character())
for (element_type in somatic %$% Subfamily %>% unique()) {
    germlinetemp <- germline %>% filter(Subfamily == element_type)
    germline_tsd_95 <- germlinetemp %$% TSD %>%
        nchar() %>%
        replace_na(0) %>%
        quantile(0.95, na.rm = TRUE)
    germline_trsd3_95 <- germlinetemp %$% Transduction_3p %>%
        nchar() %>%
        replace_na(0) %>%
        quantile(0.95, na.rm = TRUE)
    germline_trsd5_95 <- germlinetemp %$% Transduction_5p %>%
        nchar() %>%
        replace_na(0) %>%
        quantile(0.95, na.rm = TRUE)
    germline_endte_05 <- germlinetemp %$% EndTE %>%
        quantile(0.05, na.rm = TRUE) %>%
        unname()
    germline_insert_characteristics <- germline_insert_characteristics %>% add_row(germline_tsd_95 = germline_tsd_95, germline_trsd3_95 = germline_trsd3_95, germline_trsd5_95 = germline_trsd5_95, germline_endte_05 = germline_endte_05, Subfamily = element_type)
}

#### transduction mapping
# overall it seems like the transductions are just more duplicated adjacent genomic sequence - perhaps it had too many consequtive mismatches to be called in the TSD.
get_promising_transduction <- function(insertdf) {
    promising_transductions <- list()
    for (sample in conf$samples) {
        sf_with_trsd <- insertdf %>%
            filter(sample_name == sample) %>%
            filter(nchar(Transduction_3p) > 20)
        if (nrow(sf_with_trsd) == 0) {
            next
        }

        trsd_ss <- DNAStringSet(sf_with_trsd$Transduction_3p)
        names(trsd_ss) <- sf_with_trsd$UUID


        at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)

        trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]


        bl <- blast(db = sub("\\.[^.]*$", "", grep(sample, inputs$blast_njs, value = TRUE)))
        bres <- tibble(predict(bl, trsd_ss_for_blast))
        if (nrow(bres) == 0) {
            next
        }
        bres1 <- bres %>% left_join(sf_with_trsd, by = c("qseqid" = "UUID"))

        bres_hits <- bres1 %>%
            group_by(qseqid) %>%
            filter(bitscore == max(bitscore)) %>%
            filter(pident == max(pident)) %>%
            filter(gapopen == min(gapopen)) %>%
            filter(length == max(length)) %>%
            ungroup()

        bres_hits_other_loc <- bres_hits %>% filter(!((sseqid == seqnames) & (abs(sstart - start) < 2000)))

        relevant_subfamilies <- bres_hits_other_loc %$% Subfamily %>% unique()
        relevant_subfamilies_grs <- GRanges(rmann %>% filter(grepl(paste(relevant_subfamilies, sep = "|", colapse = TRUE), rte_subfamily)))

        bresgrs <- GRanges(bres_hits_other_loc)
        # extend by 500 bp on either side
        bresgrs <- resize(bresgrs, width = width(bresgrs) + 500, fix = "center")

        trsd_adjacent_rte <- subsetByOverlaps(bresgrs, relevant_subfamilies_grs) %>%
            as.data.frame() %>%
            tibble()
        promising_transductions[[sample]] <- trsd_adjacent_rte
        if (length(trsd_adjacent_rte) > 0) {
            print("promising_ids")
            trsd_adjacent_rte %$% qseqid
            print("end_promising_hits")
        }
    }
    return(promising_transductions)
}

somatic_repregion_filt_promosing_transductions <- get_promising_transduction(somatic_repregion_filt)

######

# somatic
# nice inserts
insert_id <- "8a2116c5-6b9a-44ee-ba15-b3940511a048"
insert_id <- "49a6f96f-d10a-413d-8627-526aefbc1904"
# odd insert
insert_id <- "1d4c8f79-de32-492b-9ac1-d6b6e5d192db"
# first only elements which pass all filters



analyze_insert <- function(group_frame, insert_id, insert_outputdir) {
    tryCatch(
        {
            insert_row <- group_frame %>% filter(UUID == insert_id)
            insert_length <- insert_row$LengthIns
            tsd_length <- nchar(insert_row$TSD)
            trsd5_length <- nchar(insert_row$Transduction_5p)
            trsd3_length <- nchar(insert_row$Transduction_3p)
            sample <- insert_row$sample_name

            # making a folder with info to assess this insertion
            dir.create(insert_outputdir, recursive = TRUE)
            insert_row %>% write_delim(sprintf("%s/insert_info.tsv", insert_outputdir), delim = "\t")

            # # Creating dot plot of read vs insert site
            bampath <- grep(sample, inputs$bam, value = TRUE)
            whichgr <- GRanges(group_frame %>% filter(sample_name == sample) %>% filter(UUID == insert_id))
            # flag <- scanBamFlag(
            #     isUnmappedQuery = NA,
            #     isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
            #     isDuplicate = NA, isSupplementaryAlignment = NA
            # )
            # aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
            # alndf <- as.data.frame(aln1) %>%
            #     tibble()

            # if (conf$update_ref_with_tldr$per_sample == "yes") {
            #     tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
            # } else {
            #     tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", "A.REF", "A.REF", insert_id))
            # }
            # read_name <- tldr_df %>% filter(Useable == TRUE & IsSpanRead == TRUE) %$% ReadName
            # read_of_interest <- alndf %>% filter(qname == read_name)
            # read_length <- read_of_interest$seq %>% nchar()
            # cigar <- read_of_interest$cigar

            # inserts <- str_extract_all(read_of_interest$cigar, "[0-9]+I") %>%
            #     unlist() %>%
            #     gsub("I", "", .) %>%
            #     as.numeric()
            # insert_start_in_read <- str_split(cigar, paste0(inserts[inserts > 50], "I"))[[1]][1] %>%
            #     str_split(., pattern = "[a-zA-Z]") %>%
            #     unlist() %>%
            #     as.numeric() %>%
            #     sum(na.rm = TRUE)

            # read_seq_sense <- DNAStringSet(read_of_interest$seq)
            # tryCatch(
            #     {
            #         read_seq_sense_1001_window <<- subseq(read_seq_sense, start = insert_start_in_read - 500, end = insert_start_in_read + 500)
            #     },
            #     error = function(e) {
            #         read_seq_sense_1001_window <<- subseq(read_seq_sense, start = insert_start_in_read - 200, end = insert_start_in_read + 400)
            #     }
            # )

            # read_seq_split <- unlist(strsplit(unname(as.character(read_seq_sense_1001_window)), split = ""))
            consensus_small_split <- toupper(unlist(strsplit(unname(as.character(insert_row$Consensus)), split = "")))
            if (insert_row$strand == "+") {
                insert_site_sense <- getSeq(fa, whichgr + round((length(consensus_small_split) - width(whichgr)) / 2))
            } else {
                insert_site_sense <- getSeq(fa, whichgr + round((length(consensus_small_split) - width(whichgr)) / 2)) %>% reverseComplement()
            }

            insert_site_split <- unlist(strsplit(unname(as.character(insert_site_sense)), split = ""))
            # pdf(sprintf("%s/dotplot.pdf", insert_outputdir))
            # dotPlot(read_seq_split, insert_site_split, wsize = 10, nmatch = 9)
            # dev.off()
            pdf(sprintf("%s/dotplot_consensussmall_vs_consensussmall.pdf", insert_outputdir))
            dotPlot(consensus_small_split, consensus_small_split, wsize = 10, nmatch = 9)
            dev.off()
            pdf(sprintf("%s/dotplot_insertsite_vs_insertsite.pdf", insert_outputdir))
            dotPlot(insert_site_split, insert_site_split, wsize = 10, nmatch = 9)
            dev.off()
            pdf(sprintf("%s/dotplot_consensussmall_vs_insertsite.pdf", insert_outputdir))
            dotPlot(consensus_small_split, insert_site_split, wsize = 10, nmatch = 9)
            dev.off()
            writeXStringSet(insert_site_sense, sprintf("%s/insert_site.fa", insert_outputdir))
            # writeXStringSet(read_seq_sense_1001_window, sprintf("%s/supporting_read_1000_window.fa", insert_outputdir))
            # writeXStringSet(read_seq_sense, sprintf("%s/supporting_read.fa", insert_outputdir))
            # Creating an annotated fasta file of insert
            tryCatch(
                {
                    # insert_cons_df <- tibble(Consensus = read_lines(sprintf("aref/%s_tldr/%s/%s.cons.ref.fa", sample, sample, insert_row$UUID))[[2]])
                    insert_cons_df <- tibble(Consensus = insert_row$Consensus)

                    insert_cons_df <- insert_cons_df %>%
                        mutate(Consensus_lower = str_extract_all(Consensus, pattern = "[:lower:]+")) %>%
                        rowwise() %>%
                        mutate(insLenMatch = which(as.vector((nchar(Consensus_lower))) == insert_row$LengthIns)) %>%
                        mutate(ins_consensus_noflank = Consensus_lower[insLenMatch]) %>%
                        mutate(ins_consensus_30flank = str_sub(Consensus, str_locate(Consensus, ins_consensus_noflank)[1] - 30, str_locate(Consensus, ins_consensus_noflank)[2] + 30)) %>%
                        mutate(ins_start = str_locate(Consensus, ins_consensus_noflank)[1]) %>%
                        mutate(ins_end = str_locate(Consensus, ins_consensus_noflank)[2])

                    insert_cons_ss <- DNAStringSet(insert_cons_df$Consensus)
                    names(insert_cons_ss) <- insert_row$UUID
                    writeXStringSet(insert_cons_ss, sprintf("%s/consensus_%s.fa", insert_outputdir, insert_row$UUID))

                    insert_fa_path <- sprintf("%s/all_inserts.fa", dirname(insert_outputdir))
                    if (file.exists(insert_fa_path)) {
                        if (!any(grepl(names(insert_cons_ss), insert_fa_path))) {
                            writeXStringSet(insert_cons_ss, insert_fa_path, append = TRUE)
                        }
                    } else {
                        writeXStringSet(insert_cons_ss, insert_fa_path)
                    }

                    insert_grs <- GRanges(seqnames = insert_row$UUID, ranges = IRanges(start = insert_cons_df$ins_start, insert_cons_df$ins_end), name = "insert")
                    # TSD ANNOT
                    # allowing no mismatches
                    nomismatch <- vmatchPattern(insert_row$TSD, insert_cons_ss, max.mismatch = 0)[[1]]
                    if (length(nomismatch) > 0) {
                        tsd_grs_no_mismatch <- GRanges(seqnames = insert_row$UUID, nomismatch, name = "tsd_no_mismatch")
                    } else {
                        tsd_grs_no_mismatch <- GRanges()
                    }
                    # allowing some mismatches for nanopore error
                    max_mismatch <- 0
                    twomatchingregions <- FALSE
                    while (twomatchingregions == FALSE) {
                        maybe_mismatch <- vmatchPattern(insert_row$TSD, insert_cons_ss, max.mismatch = max_mismatch)[[1]]
                        if (length(maybe_mismatch) > 1) {
                            tsd_grs_maybe_mismatch <<- GRanges(seqnames = insert_row$UUID, maybe_mismatch, name = sprintf("tsd_up_to_%s_mismatch", max_mismatch))
                            twomatchingregions <- TRUE
                        } else {
                            max_mismatch <- max_mismatch + 1
                        }
                        if (max_mismatch > 30) {
                            twomatchingregions <- TRUE
                        }
                    }
                    tsd_grs <- c(tsd_grs_no_mismatch, tsd_grs_maybe_mismatch)

                    grs <- c(insert_grs, tsd_grs)
                    # trd annot
                    if (!is.na(insert_row$Transduction_5p)) {
                        if (length(insert_row$Transduction_5p) > 7) {
                            trd5_grs <- GRanges(seqnames = insert_row$UUID, vmatchPattern(insert_row$Transduction_5p, insert_cons_ss, max.mismatch = 0)[[1]], name = "trd5")
                            grs <- c(grs, trd5_grs)
                        }
                    }
                    if (!is.na(insert_row$Transduction_3p)) {
                        if (length(insert_row$Transduction_3p) > 7) {
                            trd3_grs <- GRanges(seqnames = insert_row$UUID, vmatchPattern(insert_row$Transduction_3p, insert_cons_ss, max.mismatch = 0)[[1]], name = "trd3")
                            grs <- c(grs, trd3_grs)
                        }
                    }
                    rtracklayer::export(grs, sprintf("%s/insert_structure_annot.gtf", insert_outputdir), format = "gtf")
                    system(
                        sprintf(
                            "awk '!/#/ {print}' %s > %s",
                            sprintf("%s/insert_structure_annot.gtf", insert_outputdir), sprintf("%s/insert_structure_annot.good.gtf", insert_outputdir)
                        )
                    )
                },
                error = function(e) {
                    print(sprintf("%s failed", insert_row$UUID))
                    print(e)
                }
            )
        },
        error = function(e) {
            print(insert_id)
            print(e)
        }
    )
    return("done")
}

# repeatmask
callRM <- function(conda_base_path, rm_outputdir, species, insertsfa) {
    system(
        sprintf(
            "source %s/etc/profile.d/conda.sh && conda activate rmtest && RepeatMasker -species %s -pa 10 -gff %s -dir %s",
            conda_base_path, species, insertsfa, rm_outputdir
        )
    )
    system(
        sprintf(
            "awk '!/#/ {print}' %s/all_inserts.fa.out.gff > %s/all_inserts.fa.out.gff.temp; mv %s/all_inserts.fa.out.gff.temp %s/all_inserts.fa.out.gff; rm %s/all_inserts.fa.out.gff.temp",
            rm_outputdir, rm_outputdir, rm_outputdir, rm_outputdir, rm_outputdir
        )
    )
    system(
        sprintf(
            "awk '{IFS=OFS=\"\\t\"} NR>3 {if ($9 == \"C\") $9 = \"-\"; print $5,$6,$7,$10,$1,$9}' %s/all_inserts.fa.out > %s/all_inserts.fa.out.bed",
            rm_outputdir, rm_outputdir
        )
    )
}

extractRMperID <- function(group_frame, insert_id, insert_outputdir) {
    print(insert_id)
    insert_row <- group_frame %>% filter(UUID == insert_id)
    sample <- insert_row$sample_name
    system(
        sprintf(
            "awk '/%s/ {print}' %s/all_inserts.fa.out.bed > %s/rm_%s.bed", insert_id, dirname(insert_outputdir), insert_outputdir, insert_id
        )
    )
    tryCatch(
        {
            teranges <- import(sprintf("%s/rm_%s.bed", insert_outputdir, insert_id))
            export(teranges, sprintf("%s/rm_%s.temp.gtf", insert_outputdir, insert_id))

            system(
                sprintf(
                    "awk '!/#/ {print}' %s/rm_%s.temp.gtf > %s/rm_%s.good.gtf ; rm %s/rm_%s.temp.gtf", insert_outputdir, insert_id, insert_outputdir, insert_id, insert_outputdir, insert_id
                )
            )
        },
        error = function(e) {
            print(e)
        }
    )
}

analyze_inserts <- function(group_frame, element_group_name) {
    print(element_group_name)
    for (insert_id in group_frame$UUID) {
        print(insert_id)
        insert_outputdir <- sprintf("%s/unique_insert_data/%s/%s", outputdir, element_group_name, insert_id)
        analyze_insert(group_frame, insert_id, insert_outputdir)
    }

    conda_base_path <- system("conda info --base", intern = TRUE)
    rm_outputdir <- dirname(insert_outputdir)
    species <- conf$species
    insertsfa <- sprintf("%s/all_inserts.fa", rm_outputdir)
    callRM(conda_base_path, rm_outputdir, species, insertsfa)

    for (insert_id in group_frame$UUID) {
        insert_outputdir <- sprintf("%s/unique_insert_data/%s/%s", outputdir, element_group_name, insert_id)
        extractRMperID(group_frame, insert_id, insert_outputdir)
    }
}

tsd_pass_elements <- somatic_filt %>% filter(nchar(TSD) < 20)
tsd_nonpass_elements <- somatic_filt %>% filter(nchar(TSD) >= 20)

insert_frames <- list(
    "tsd_pass" = tsd_pass_elements,
    "inrepregion" = somatic_repregion_filt
)

imap(insert_frames, ~ analyze_inserts(.x, .y))

insert_frames <- list(
    "germline" = germline_insert_df[1:20, ] %>% filter(Family == "ALU")
)

imap(insert_frames, ~ analyze_inserts(.x, .y))
### learned features from elements that passed manual curation
curation_df <- tibble(UUID = c(
    "996b0eb6-2629-4e8f-94b3-ea1e9befe22d",
    "8a2116c5-6b9a-44ee-ba15-b3940511a048",
    "8ca0fd30-d707-4292-af43-cccc1c264aa4",
    "712ac624-a5d1-44c3-a49e-97c9149cf041",
    "99739fc6-95f2-4e2f-8855-fc3628579e7c",
    "ce5013f1-5000-45ab-9091-bb89b1cad6ce",
    "d5413498-68a1-447e-aae9-10fa564c5197",
    "edba26cc-ae66-49d7-b94c-a7a47dd97aa0"
), pass_curation = TRUE)

tsd_pass_elements %>%
    left_join(curation_df) %>%
    filter(Family != "SVA") %>%
    filter(pass_curation == TRUE) %>%
    dplyr::select(width, strand, StartTE, EndTE, Inversion, UnmapCover, TEMatch, Remappable, Filter)
tsd_pass_elements %>%
    left_join(curation_df) %>%
    filter(Family != "SVA") %>%
    filter(is.na(pass_curation)) %>%
    dplyr::select(width, strand, StartTE, EndTE, Inversion, UnmapCover, TEMatch, Remappable, Filter)

# lessons are that an Alu should completely span StartTE - EndTE, and that TEmatch should be over 90
# also all "SVA" observed were just stretches of poly A
tsd_nonpass_elements %>%
    filter(Family == "ALU") %>%
    mutate(startcheck = ifelse(StartTE > 5, FALSE, TRUE)) %$% startcheck %>%
    table()
tsd_nonpass_elements %>%
    filter(Family == "ALU") %>%
    mutate(matchcheck = ifelse(TEMatch > 90, TRUE, FALSE)) %$% matchcheck %>%
    table()

tsd_nonpass_start_match_pass <- tsd_nonpass_elements %>%
    filter(Family == "ALU") %>%
    filter(LengthIns > 290) %>%
    mutate(startcheck = ifelse(StartTE < 5, TRUE, FALSE)) %>%
    mutate(matchcheck = ifelse(TEMatch > 90, TRUE, FALSE)) %>%
    mutate(both = matchcheck & startcheck) %>%
    filter(both == TRUE)

analyze_inserts(tsd_nonpass_start_match_pass, "tsd_nonpass")



# all elements that passed scrutiny:
curation_df_complete <- tibble(UUID = c(
    "996b0eb6-2629-4e8f-94b3-ea1e9befe22d",
    "8a2116c5-6b9a-44ee-ba15-b3940511a048",
    "8ca0fd30-d707-4292-af43-cccc1c264aa4",
    "712ac624-a5d1-44c3-a49e-97c9149cf041",
    "99739fc6-95f2-4e2f-8855-fc3628579e7c",
    "ce5013f1-5000-45ab-9091-bb89b1cad6ce",
    "d5413498-68a1-447e-aae9-10fa564c5197",
    "edba26cc-ae66-49d7-b94c-a7a47dd97aa0",
    "73787540-8575-4ada-9e04-4c94c46a06ee"
), pass_curation = TRUE)

pass <- somatic_all %>%
    left_join(curation_df_complete) %>%
    filter(pass_curation == TRUE) %>%
    mutate(sample_name = ordered(sample_name, levels = conf$samples)) %>%
    mutate(condition = factor(condition, levels = conf$levels))

p <- pass %>%
    ggplot(aes(x = sample_name, fill = condition)) +
    geom_bar() +
    facet_wrap(~Family) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_family.pdf", outputdir), 10, 5)

p <- pass %>%
    ggplot(aes(x = sample_name, fill = Subfamily)) +
    geom_bar() +
    labs(x = "", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_limited.pdf", outputdir), 5, 5)

p <- pass %>%
    group_by(sample_name, Subfamily) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, Subfamily, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    ggplot(aes(x = sample_name, fill = Subfamily)) +
    geom_col(aes(y = n)) +
    labs(x = "", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar.pdf", outputdir), 5, 5)

p <- pass %>%
    group_by(sample_name, Subfamily) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, Subfamily, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    ggplot(aes(x = sample_name)) +
    geom_point(aes(x = age, y = n, color = sex), size = 3) +
    labs(x = "Age", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_palette +
    scale_y_continuous(limits = c(-0.5, NA)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_age.pdf", outputdir), 5, 5)

p <- pass %>%
    group_by(sample_name, Subfamily) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, Subfamily, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    ggplot(aes(x = sample_name, fill = Subfamily)) +
    geom_point(aes(x = braak, y = n), size = 3) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_palette +
    scale_y_continuous(limits = c(-0.5, NA)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_braak.pdf", outputdir), 5, 5)

p <- pass %>%
    group_by(sample_name, Subfamily) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, Subfamily, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    ggplot(aes(x = sample_name, fill = Subfamily)) +
    geom_point(aes(x = braak, y = n), size = 3) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
    mtclosed +
    anchorbar +
    scale_palette +
    scale_y_continuous(limits = c(-1, NA)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_braak.pdf", outputdir), 5, 5)

p <- pass %>%
    group_by(sample_name) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    left_join(sample_sequencing_data) %>%
    ggplot(aes(x = sample_name)) +
    geom_point(aes(x = bases_number, y = n, color = sex), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Bases Sequenced", title = "RTE Somatic Insertions") +
    mtopen +
    anchorbar +
    scale_palette +
    scale_y_continuous(limits = c(-0.5, NA)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/curated_elements/bar_bases.pdf", outputdir), 5, 5)
########### OLD CODE

# germline
for (insert_id in germline$UUID) {
    insert_row <- germline %>% filter(UUID == insert_id)
    bampath <- inputs$bam[[1]]
    whichgr <- GRanges(germline %>% filter(UUID == insert_id))
    flag <- scanBamFlag(
        isUnmappedQuery = NA,
        isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
        isDuplicate = NA, isSupplementaryAlignment = NA
    )
    aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
    alndf <- as.data.frame(aln1) %>%
        tibble()
    # insert fails if it was captured in a read which has supplementary alignments
    if (conf$update_ref_with_tldr$per_sample == "yes") {
        tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
    } else {
        tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", "A.REF", "A.REF", insert_id))
    }
    read_name <- tldr_df %>%
        filter(Useable == TRUE & IsSpanRead == TRUE) %$% ReadName %>%
        pluck(1)
    read_of_interest <- alndf %>% filter(qname == read_name)
    read_length <- read_of_interest$seq %>% nchar()
    cigar <- read_of_interest$cigar

    inserts <- str_extract_all(read_of_interest$cigar, "[0-9]+I") %>%
        unlist() %>%
        gsub("I", "", .) %>%
        as.numeric()
    insert_start_in_read <- str_split(cigar, paste0(inserts[inserts > 50], "I"))[[1]][1] %>%
        str_split(., pattern = "[a-zA-Z]") %>%
        unlist() %>%
        as.numeric() %>%
        sum(na.rm = TRUE)

    read_seq_sense <- DNAStringSet(read_of_interest$seq)
    read_seq_sense <- subseq(read_seq_sense, start = insert_start_in_read - 600, end = insert_start_in_read + 800)
    tryCatch(
        {
            if (insert_row$strand == "+") {
                insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2))
            } else {
                insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2)) %>% reverseComplement()
            }
        },
        error = function(e) {
            if (insert_row$Strand == "+") {
                insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2))
            } else {
                insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2)) %>% reverseComplement()
            }
        }
    )

    path <- "temp_read.fa"
    writeXStringSet(read_seq_sense, path)
    makeblastdb(path)
    bl <- blast(db = path)
    bres <- tibble(predict(bl, insert_site_sense))
    insert_site_split <- unlist(strsplit(unname(as.character(insert_site_sense)), split = ""))
    read_seq_split <- unlist(strsplit(unname(as.character(read_seq_sense)), split = ""))

    plotpath <- sprintf("%s/germline_dotplots/%s_insertlen_%s_%s.pdf", outputdir, insert_row$Subfamily, insert_row$LengthIns, gsub("-", "_", insert_id))
    dir.create(dirname(plotpath), recursive = TRUE)
    pdf(plotpath)
    dotPlot(read_seq_split, insert_site_split, wsize = 10, nmatch = 9)
    dev.off()
} # TODO STILL HAVE AN ISSUE WITH THESE DOT PLOTS.

tldr_te_ref_path <- conf$update_ref_with_tldr$tldr_te_ref[[conf$species]]

my_tsd_check <- tibble(sample_name = character(), UUID = character(), my_tsd = character(), my_tsd_left_start = numeric(), my_tsd_right_start = numeric(), my_te_start = numeric(), my_te_end = numeric())
for (sample in conf$samples) {
    print(sample)
    bampath <- grep(sample, inputs$bam, value = TRUE)
    insert_ids <- somatic %>% filter(sample_name == sample) %$% UUID
    print(length(insert_ids))
    for (insert_id in insert_ids) {
        insert_row <- somatic %>%
            filter(sample_name == sample) %>%
            filter(UUID == insert_id)
        insert_loc <- GRanges(insert_row)
        css <- DNAStringSet(insert_row$Consensus)

        tryCatch(
            {
                bl <- blast(db = tldr_te_ref_path)
                bres <- tibble(predict(bl, css))
                hit <- bres %>%
                    filter(grepl(insert_row$Subfamily, sseqid)) %>%
                    mutate(length_dif_from_insert = abs(length - insert_row$LengthIns)) %>%
                    arrange(length_dif_from_insert) %>%
                    head(1)

                # I scan the 220 bp region flanking the predicted insert for TSD. Will fail if there is a large five or three prime transduction
                five_start <- hit$qstart - 200
                five_end <- hit$qstart + 20
                three_start <- hit$qend - 20
                three_end <- hit$qend + 200

                five_string <- str_sub(as.character(css), start = five_start, end = five_end)
                three_string <- str_sub(as.character(css), start = three_start, end = three_end)
                if (is.null(five_string) | is.null(three_string)) {
                    print("empty strings")
                } else {
                    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1000, baseOnly = TRUE)
                    paln <- pairwiseAlignment(DNAStringSet(five_string), DNAStringSet(three_string), type = "local", substitutionMatrix = mat, gapOpening = 1000)
                    tsd_loc <- str_locate_all(as.character(css), as.character(paln@subject))[[1]]
                    row <- tibble(sample_name = sample, UUID = insert_id, my_tsd = as.character(paln@pattern), my_tsd_left_start = as.data.frame(tsd_loc)$start[1], my_tsd_right_start = as.data.frame(tsd_loc)$start[2], my_te_start = hit$qstart, my_te_end = hit$qend)
                    my_tsd_check <- bind_rows(my_tsd_check, row)
                }
                # write(paste0(">five\n", five_string), "five.fa")
                # makeblastdb("five.fa")
                # bl <- blast(db = "five.fa")
                # bres <- tibble(predict(bl, DNAStringSet(three_string)))
            },
            error = function(e) {
                "issue"
            }
        )
    }
}






for (element_type in somatic_filt %$% Subfamily %>% unique()) {
    dftemp <- somatic %>% filter(Subfamily == element_type)
    p <- dftemp %>%
        gghistogram(x = "LengthIns", fill = "blue") +
        mtopen +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/filter_set1/length_distribution.pdf", outputdir, element_type))
    dftemp %$% LengthIns %>% quantile()
    p <- dftemp %>%
        mutate(tsdlen = nchar(TSD)) %>%
        gghistogram(x = "tsdlen", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set1/tsd_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd5 = nchar(Transduction_5p)) %>%
        gghistogram(x = "trsd5", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set1/trsd5_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd3 = nchar(Transduction_3p)) %>%
        gghistogram(x = "trsd3", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set1/trsd3_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        arrange(-EndTE) %>%
        mutate(nrow = row_number()) %>%
        ggplot() +
        geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
        mtclosed +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/filter_set1/te_body_distribution.pdf", outputdir, element_type))
}

for (element_type in somatic_filt %$% Subfamily %>% unique()) {
    dftemp <- somatic_filt %>% filter(Subfamily == element_type)
    p <- dftemp %>%
        gghistogram(x = "LengthIns", fill = "blue") +
        mtopen +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/filter_set2/length_distribution.pdf", outputdir, element_type))
    dftemp %$% LengthIns %>% quantile()
    p <- dftemp %>%
        mutate(tsdlen = nchar(TSD)) %>%
        gghistogram(x = "tsdlen", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set2/tsd_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd5 = nchar(Transduction_5p)) %>%
        gghistogram(x = "trsd5", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set2/trsd5_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd3 = nchar(Transduction_3p)) %>%
        gghistogram(x = "trsd3", fill = "blue") +
        labs(title = element_type) + mtopen
    mysaveandstore(sprintf("%s/%s/filter_set2/trsd3_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        arrange(-EndTE) %>%
        mutate(nrow = row_number()) %>%
        ggplot() +
        geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
        mtclosed +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/filter_set2/te_body_distribution.pdf", outputdir, element_type))
}


for (sample in unique(somatic_filt$sample_name)) {
    tempoutputdir <- sprintf("aref/%s_Analysis/tldr_plots/somatic", sample)
    dfallsample <- dfall %>% filter(sample_name == sample)
    somaticsample <- somatic_filt %>% filter(sample_name == sample)

    somaticsample %>%
        write_delim(sprintf("%s/somatic_pass.tsv", tempoutputdir), delim = "\t")

    p <- somaticsample %>%
        ggplot(aes(x = Family)) +
        geom_bar() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/somatic_bar.pdf", tempoutputdir), 5, 4)


    p <- somaticsample %>%
        ggplot(aes(x = LengthIns)) +
        geom_histogram() +
        labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions") +
        facet_wrap(~Family) +
        mtclosed +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/somatic_pass_insertion_length.pdf", tempoutputdir), 5, 3)
}


p <- somatic_filt %>%
    ggplot(aes(x = sample_name, fill = condition)) +
    geom_bar() +
    facet_wrap(~Family) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/somatic_bar.pdf", outputdir), 10, 5)

p <- somatic_filt %>%
    group_by(sample_name, Subfamily, condition) %>%
    summarise(nins = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Subfamily, values_from = nins) %>%
    replace(is.na(.), 0) %>%
    ggplot(aes(x = L1HS, y = AluY)) +
    geom_point(aes(color = condition)) +
    scale_conditions +
    stat_cor() +
    stat_smooth(method = "lm", formula = y ~ x, geom = "smooth") +
    labs(x = "L1HS", y = "AluY", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/somatic_pass_l1hs_aluy_corr_bar.pdf", outputdir), 4, 4)

somatic_filt %>%
    filter(Subfamily %in% c("AluY")) %>%
    filter(sample_name == "AD1") %>%
    filter(UUID == "a72cba1d-f8c4-4dd4-8054-ee7c97a6c704") %>%
    print(width = Inf)

p <- somatic_filt %>%
    filter(Subfamily %in% c("L1HS", "AluY")) %>%
    ggplot(aes(x = TEMatch, fill = condition)) +
    geom_histogram() +
    facet_grid2(rows = vars(sample_name), cols = vars(Subfamily), scale = "free_x", axes = "all", remove_labels = "y") +
    labs(x = "TE Match", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/somatic_histogram.pdf", outputdir), 5, 20)

p <- somatic_filt %>%
    filter(Subfamily %in% c("AluY")) %>%
    ggplot(aes(x = sample_name, fill = condition)) +
    geom_bar() +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/somatic_bar_aluy.pdf", outputdir), 8, 5)



p <- somatic_filt %>%
    ggplot(aes(x = LengthIns)) +
    geom_histogram() +
    labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    facet_grid(rows = vars(sample_name), cols = vars(Family), scales = "free") +
    mtclosed +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/somatic_insertion_length.pdf", outputdir), 10, 10)

somatic_filt %>%
    filter(LengthIns > 5000) %>%
    filter(Subfamily == "L1HS") %>%
    print(width = Inf)

tryCatch(
    {
        metadata_vars <- colnames(sample_table)[!(colnames(sample_table) %in% c("sample_name", "condition"))]
        sequencing_metadata_vars <- c("N50", "reads_number", "bases_number")


        grouping_vars <- c("sample_name", "condition", sequencing_metadata_vars, metadata_vars)
        germline_n <- germline %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(Insert_Type = "Germline")
        somatic_n <- somatic_filt %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(Insert_Type = "Somatic")
        insert_n <- bind_rows(germline_n, somatic_n)
        p <- insert_n %>% ggplot(aes(x = N50, y = nins, color = sample_name)) +
            geom_point(size = 2) +
            facet_wrap(~Insert_Type) +
            labs(x = "N50", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/n50_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- insert_n %>% ggplot(aes(x = reads_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~Insert_Type) +
            labs(x = "Total Reads Count", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/reads_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- insert_n %>% ggplot(aes(x = bases_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~Insert_Type) +
            labs(x = "Total Bases Count", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/bases_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        predictors <- c(sequencing_metadata_vars, metadata_vars, "condition")
        model_formula <- paste("nins ~ ", paste(predictors, collapse = "+"), sep = "")

        model <- lm(as.formula(model_formula), data = insert_n)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_allins2.txt", outputdir), delim = "\t")

        model <- lm(as.formula(model_formula), data = somatic_n)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_somatic.txt", outputdir), delim = "\t")

        predictors <- c(sequencing_metadata_vars, "condition")
        model_formula <- paste("nins ~ ", paste(predictors, collapse = "+"), sep = "")
        model <- lm(as.formula(model_formula), data = somatic_n)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_somatic_model_limited.txt", outputdir), delim = "\t")


        germline_n_bysubfam <- germline %>%
            group_by(across(c(grouping_vars, Subfamily))) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(Insert_Type = "Germline")
        somatic_n_bysubfam <- somatic_filt %>%
            group_by(across(c(grouping_vars, Subfamily))) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(Insert_Type = "Somatic")
        insert_n_bysubfam <- bind_rows(germline_n_bysubfam, somatic_n_bysubfam)
        for (insertion_type in somatic_n_bysubfam %$% Subfamily %>% unique()) {
            insert_n_bysubfam <- somatic_n_bysubfam %>% filter(Subfamily == insertion_type)
            for (mvar in metadata_vars) {
                tryCatch(
                    {
                        if (sample_table[[mvar]] %>% is.numeric()) {
                            p <- insert_n_bysubfam %>% ggscatter(x = mvar, y = "nins", color = "condition", size = 3) +
                                scale_conditions + stat_cor() +
                                stat_smooth(method = "lm", formula = y ~ x, geom = "smooth")
                        } else {
                            p <- insert_n_bysubfam %>% ggplot(aes(x = sample_name, y = nins, fill = condition)) +
                                facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") +
                                geom_col() +
                                labs(x = "", y = "Number of Insertions", , title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") +
                                mtclosed +
                                anchorbar +
                                scale_conditions +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        }
                        mysaveandstore(sprintf("%s/%s_by_%s_somatic_pass.pdf", outputdir, insertion_type, mvar), 5, 4)
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        }
    },
    error = function(e) {
        message("Likely single condition")
    }
)

##########

x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)


somatic_filt %>%
    filter(Subfamily == "L1HS") %>%
    print(width = Inf)
somatic_filt_out
# ctrl6
