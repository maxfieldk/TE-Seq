module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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


tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = sprintf("aref/%s_tldr/%s.table.txt", conf$samples, conf$samples),
            json = sprintf("aref/qc/%s/%spycoQC.json", conf$samples, conf$samples),
            filtered_tldr = sprintf("aref/default/%s.table.kept_in_updated_ref.txt", sample_table$sample_name),
            bam = sprintf("aref/intermediates/%s/alignments/5khz/%s.hac.5mCG_5hmCG.sorted.bam", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/results/somatic_insertions/analyze_nongermline_insertions.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


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

# only tldr germinline tldr insertions that wind up in reference
dfs_filtered <- list()
for (sample in sample_table$sample_name) {
    df <- read.table(grep(sprintf("%s.table", sample), inputs$filtered_tldr, value = TRUE), header = TRUE)
    df$sample_name <- sample
    df <- df %>% left_join(sample_table)
    dfs_filtered[[sample]] <- df
}

germline <- do.call(rbind, dfs_filtered) %>%
    tibble() %>%
    left_join(sample_sequencing_data)


# all TLDR insertions
library(BSgenome)
fa <- Rsamtools::FaFile(conf$ref)
dflist <- list()
for (sample in sample_table$sample_name) {
    if (conf$update_ref_with_tldr$per_sample == "yes") {
        df <- read.table(grep(sprintf("%s_tldr", sample), inputs$tldroutput, value = TRUE), header = TRUE) %>%
            mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
            tibble()
    } else {
        df <- read.table("aref/A.REF_tldr/A.REF.table.txt", header = TRUE) %>%
            mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
            tibble()
    }
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


somatic <- GRanges(somatic) %>%
    subsetByOverlaps(centromere, invert = TRUE) %>%
    subsetByOverlaps(telomere, invert = TRUE) %>%
    subsetByOverlaps(segdups, invert = TRUE) %>%
    subsetByOverlaps(mappability) %>%
    as.data.frame() %>%
    tibble()

library(GenomicAlignments)
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
    insert_ids <- somatic %>% filter(sample_name == sample) %$% UUID
    print(length(insert_ids))
    for (insert_id in insert_ids) {
        insert_row <- somatic %>%
            filter(sample_name == sample) %>%
            filter(UUID == insert_id)
        whichgr <- GRanges(somatic %>% filter(sample_name == sample) %>% filter(UUID == insert_id))
        aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
        alndf <- as.data.frame(aln1) %>%
            tibble()
        mean_mapq <- alndf$mapq %>% mean()
        insert_mean_mapqs[[insert_id]] <- mean_mapq

        # insert fails if it was captured in a read which has supplementary alignments
        tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
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

somatic <- somatic %>%
    left_join(mean_mapq_filter) %>%
    left_join(supplementary_alignment_filter) %>%
    left_join(indel_and_clipping_filter)
somatic <- somatic %>%
    filter(insert_mean_mapqs > 55) %>%
    filter(insert_supplementary_status == "PASS") %>%
    filter(insert_indel_and_clipping_status == "PASS")
somatic %>%
    write_delim(sprintf("%s/somatic_inserts.tsv", outputdir), delim = "\t")



# ADDTIONAL FILTERS
# I think these are overly restrictive given how tsd calling on a single lossy read may occur - doesn't have the benefit of deriving from cluster
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

somatic_filt <- somatic %>%
    left_join(germline_insert_characteristics) %>%
    mutate(endte_manual_filter = case_when(
        Family == "ALU" ~ 250,
        Family == "L1" ~ 5800,
    )) %>%
    filter(EndTE >= endte_manual_filter | is.na(endte_manual_filter)) %>%
    filter(EndTE >= germline_endte_05 - 50)
somatic_filt %>%
    write_delim(sprintf("%s/somatic_inserts_filt.tsv", outputdir), delim = "\t")
# somatic_filt <- read_delim(sprintf("%s/somatic_inserts_filt.tsv", outputdir))
# somatic_filt %>% filter(Family == "L1")
# filter(nchar(TSD) <= germline_tsd_95) %>%
# filter(nchar(Transduction_5p) <= germline_trsd5_95)

# somatic_filt %>% filter(sample_name == "CTRL2")
# filter(nchar(Transduction_3p) <= germline_trsd3_95)
# filter(EndTE - StartTE < LengthIns + 30) %>% #what if you just have a couple of basepairs that spuriously align to the front of TE

for (insert_id in somatic_filt$UUID) {
    insert_row <- somatic_filt %>% filter(UUID == insert_id)
    sample <- insert_row$sample_name
    bampath <- grep(sample, inputs$bam, value = TRUE)
    whichgr <- GRanges(somatic %>% filter(sample_name == sample) %>% filter(UUID == insert_id))
    aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
    alndf <- as.data.frame(aln1) %>%
        tibble()
    # insert fails if it was captured in a read which has supplementary alignments
    tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
    read_name <- tldr_df %>% filter(Useable == TRUE & IsSpanRead == TRUE) %$% ReadName
    read_of_interest <- alndf %>% filter(qname == read_name)
    read_length <- read_of_interest$seq %>% nchar()


    read_seq_sense <- DNAStringSet(read_of_interest$seq)
    if (insert_row$strand == "+") {
        insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(insert_loc)) / 2))
    } else {
        insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(insert_loc)) / 2)) %>% reverseComplement()
    }

    path <- "temp_read.fa"
    writeXStringSet(read_seq_sense, path)
    makeblastdb(path)
    bl <- blast(db = path)
    bres <- tibble(predict(bl, insert_site_sense))
    bres
    library(seqinr)
    insert_site_split <- unlist(strsplit(unname(as.character(insert_site_sense)), split = ""))
    read_seq_split <- unlist(strsplit(unname(as.character(read_seq_sense)), split = ""))


    pdf(sprintf("%s/%s_%s.pdf", outputdir, sample, insert_id))
    dotPlot(read_seq_split, insert_site_split, wsize = 10, nmatch = 9)
    dev.off()
}

library(rBLAST)
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
