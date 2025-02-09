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
            bam = sprintf("aref/intermediates/%s/alignments/%s/%s.%s.%s.sorted.bam", sample_table$sample_name, conf$rate, sample_table$sample_name, conf$type, conf$modification_string),
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

## FUNCTIONS
{
    annotate_mappability <- function(insert_df) {
        df1 <- insert_df %>%
            GRanges() %>%
            subsetByOverlaps(centromere, invert = TRUE) %>%
            subsetByOverlaps(telomere, invert = TRUE) %>%
            subsetByOverlaps(segdups, invert = TRUE) %>%
            as.data.frame() %>%
            tibble()
        df1$inrepregion <- FALSE
        df2 <- insert_df %>%
            GRanges() %>%
            subsetByOverlaps(c(centromere, telomere, segdups)) %>%
            as.data.frame() %>%
            tibble()
        df2$inrepregion <- TRUE
        df <- full_join(df1, df2)
        df1 <- df %>%
            GRanges() %>%
            subsetByOverlaps(mappability) %>%
            as.data.frame() %>%
            tibble()
        df1$k50_mappable <- TRUE
        df2 <- df %>%
            GRanges() %>%
            subsetByOverlaps(mappability, invert = TRUE) %>%
            as.data.frame() %>%
            tibble()
        df2$k50_mappable <- FALSE
        df <- full_join(df1, df2)
        return(df)
    }

    annotate_read_metadata <- function(insertdf) {
        flag <- scanBamFlag(
            isUnmappedQuery = NA,
            isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
            isDuplicate = NA, isSupplementaryAlignment = NA
        )

        insert_mean_mapqs <- list()
        insert_supplementary_status <- list()
        insert_indel_and_clipping_and_mapq_status <- list()
        insert_indel_and_clipping_status <- list()
        insert_status <- list()
        deletion_status <- list()
        clip_status <- list()
        insert_fail_reason_detail <- list()
        deletion_fail_reason_detail <- list()
        clip_fail_reason_detail <- list()

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
                read_of_interest <- alndf %>% filter(qname %in% read_name)
                insert_supplementary_status[[insert_id]] <- ifelse(all(is.na(read_of_interest$SA)), TRUE, FALSE)

                # insert fails if it was captured in a read which has extensive indels or clipping
                indelclipmapq <- list()
                indelclip <- list()
                insert <- list()
                insert_fail_reason <- list()
                deletion <- list()
                deletion_fail_reason <- list()
                clip <- list()
                clip_fail_reason <- list()
                for (read_index in 1:nrow(read_of_interest)) {
                    read <- read_of_interest[read_index, ]

                    inserts <- str_extract_all(read$cigar, "[0-9]+I") %>%
                        unlist() %>%
                        gsub("I", "", .) %>%
                        as.numeric()
                    deletions <- str_extract_all(read$cigar, "[0-9]+D") %>%
                        unlist() %>%
                        gsub("D", "", .) %>%
                        as.numeric()
                    clips <- str_extract_all(read$cigar, "([0-9]+S)|([0-9]+H)") %>%
                        unlist() %>%
                        gsub("S|H", "", .) %>%
                        as.numeric()
                    mapq <- read$mapq
                    indelclipmapq[[read_index]] <- ifelse((mapq < 60) | (sum(inserts > 50) > 1) | (sum(deletions > 50) > 0) | (sum(clips > 50) > 0), FALSE, TRUE)
                    indelclip[[read_index]] <- ifelse((sum(inserts > 50) > 1) | (sum(deletions > 50) > 0) | (sum(clips > 50) > 0), FALSE, TRUE)
                    insert[[read_index]] <- ifelse((sum(inserts > 50) > 1), FALSE, TRUE)
                    insert_fail_reason[[read_index]] <- ifelse((sum(inserts > 50) > 1), paste(inserts[inserts > 50], collapse = " "), "pass")
                    deletion[[read_index]] <- ifelse((sum(deletions > 50) > 0), FALSE, TRUE)
                    deletion_fail_reason[[read_index]] <- ifelse((sum(deletions > 50) > 0), paste(deletions[deletions > 50], collapse = " "), "pass")
                    clip[[read_index]] <- ifelse((sum(clips > 50) > 0), FALSE, TRUE)
                    clip_fail_reason[[read_index]] <- ifelse((sum(clips > 50) > 0), paste(clips[clips > 50], collapse = " "), "pass")
                }
                insert_indel_and_clipping_and_mapq_status[[insert_id]] <- any(unlist(indelclipmapq))
                insert_indel_and_clipping_status[[insert_id]] <- any(unlist(indelclip))
                insert_status[[insert_id]] <- any(unlist(insert))
                deletion_status[[insert_id]] <- any(unlist(deletion))
                clip_status[[insert_id]] <- any(unlist(clip))
                insert_fail_reason_detail[[insert_id]] <- paste(names(insert_fail_reason), insert_fail_reason, collapse = ";", sep = ":")
                deletion_fail_reason_detail[[insert_id]] <- paste(names(deletion_fail_reason), deletion_fail_reason, collapse = ";", sep = ":")
                clip_fail_reason_detail[[insert_id]] <- paste(names(clip_fail_reason), clip_fail_reason, collapse = ";", sep = ":")
            }
        }
        mean_mapq_filter <- tibble(UUID = names(insert_mean_mapqs), insert_mean_mapqs = unname(insert_mean_mapqs) %>% unlist())
        supplementary_alignment_filter <- tibble(UUID = names(insert_supplementary_status), insert_supplementary_status = unname(insert_supplementary_status) %>% map(~ .x[1]) %>% unlist())
        indel_and_clipping_and_mapq_filter <- tibble(UUID = names(insert_indel_and_clipping_and_mapq_status), insert_indel_and_clipping_and_mapq_status = unname(insert_indel_and_clipping_and_mapq_status) %>% unlist())
        indel_and_clipping_filter <- tibble(UUID = names(insert_indel_and_clipping_status), insert_indel_and_clipping_status = unname(insert_indel_and_clipping_status) %>% unlist())
        insert_status_filter <- tibble(UUID = names(insert_status), insert_status = unname(insert_status) %>% unlist())
        deletion_status_filter <- tibble(UUID = names(deletion_status), deletion_status = unname(deletion_status) %>% unlist())
        clip_status_filter <- tibble(UUID = names(clip_status), clip_status = unname(clip_status) %>% unlist())
        insert_fail_reason_detail_df <- tibble(UUID = names(insert_fail_reason_detail), insert_fail_reason_detail = unname(insert_fail_reason_detail) %>% unlist())
        deletion_fail_reason_detail_df <- tibble(UUID = names(deletion_fail_reason_detail), deletion_fail_reason_detail = unname(deletion_fail_reason_detail) %>% unlist())
        clip_fail_reason_detail_df <- tibble(UUID = names(clip_fail_reason_detail), clip_fail_reason_detail = unname(clip_fail_reason_detail) %>% unlist())

        insertdf <- insertdf %>%
            left_join(mean_mapq_filter) %>%
            left_join(supplementary_alignment_filter) %>%
            left_join(indel_and_clipping_and_mapq_filter) %>%
            left_join(indel_and_clipping_filter) %>%
            left_join(insert_status_filter) %>%
            left_join(deletion_status_filter) %>%
            left_join(clip_status_filter) %>%
            left_join(insert_fail_reason_detail_df) %>%
            left_join(deletion_fail_reason_detail_df) %>%
            left_join(clip_fail_reason_detail_df)
        return(insertdf)
    }

    annotate_teend <- function(insertdf) {
        insertdf <- insertdf %>%
            mutate(endte_manual_filter = case_when(
                Family == "ALU" ~ 250,
                Family == "L1" ~ 5800,
            ))
        return(insertdf)
    }

    analyze_insert <- function(group_frame, insert_id, insert_outputdir) {
        tryCatch(
            {
                insert_row <- group_frame %>% filter(UUID == insert_id)
                insert_length <- insert_row$LengthIns
                tsd_length <- nchar(insert_row$TSD)
                trsd5_length <- nchar(insert_row$Transduction_5p)
                trsd3_length <- nchar(insert_row$Transduction_3p)
                sample <- insert_row$sample_name

                dir.create(insert_outputdir, recursive = TRUE)
                insert_row %>% write_delim(sprintf("%s/insert_info.tsv", insert_outputdir), delim = "\t")

                # # Creating dot plot of read vs insert site
                bampath <- grep(sample, inputs$bam, value = TRUE)
                whichgr <- GRanges(group_frame %>% filter(sample_name == sample) %>% filter(UUID == insert_id))
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
                if (length(consensus_small_split) < 2000) {
                    # pdf(sprintf("%s/dotplot_consensussmall_vs_consensussmall.pdf", insert_outputdir))
                    # dotPlot(consensus_small_split, consensus_small_split, wsize = 10, nmatch = 9)
                    # dev.off()
                    # pdf(sprintf("%s/dotplot_insertsite_vs_insertsite.pdf", insert_outputdir))
                    # dotPlot(insert_site_split, insert_site_split, wsize = 10, nmatch = 9)
                    # dev.off()
                    pdf(sprintf("%s/dotplot_consensussmall_vs_insertsite.pdf", insert_outputdir))
                    dotPlot(consensus_small_split, insert_site_split, wsize = 10, nmatch = 9)
                    dev.off()
                }
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
                            tsd_grs <- GRanges(seqnames = insert_row$UUID, nomismatch, name = "tsd_no_mismatch")
                        } else {
                            # allowing some mismatches for nanopore error
                            max_mismatch <- 1
                            twomatchingregions <- FALSE
                            while (twomatchingregions == FALSE) {
                                maybe_mismatch <- vmatchPattern(insert_row$TSD, insert_cons_ss, max.mismatch = max_mismatch)[[1]]
                                if (length(maybe_mismatch) > 1) {
                                    tsd_grs <<- GRanges(seqnames = insert_row$UUID, maybe_mismatch, name = sprintf("tsd_up_to_%s_mismatch", max_mismatch))
                                    twomatchingregions <- TRUE
                                } else {
                                    max_mismatch <- max_mismatch + 1
                                }
                                if (max_mismatch > 10) {
                                    twomatchingregions <- TRUE
                                    tsd_grs <- GRanges()
                                }
                            }
                        }

                        grs <- c(insert_grs, tsd_grs)
                        # trd annot
                        if (!is.na(insert_row$Transduction_5p)) {
                            if (nchar(insert_row$Transduction_5p) > 7) {
                                trd5_grs <- GRanges(seqnames = insert_row$UUID, vmatchPattern(insert_row$Transduction_5p, insert_cons_ss, max.mismatch = 0)[[1]], name = "trd5")
                                grs <- c(grs, trd5_grs)
                            }
                        }
                        if (!is.na(insert_row$Transduction_3p)) {
                            if (nchar(insert_row$Transduction_3p) > 7) {
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

        insert_suspicion_level <- list()
        for (insert_id in group_frame$UUID) {
            insert_outputdir <- sprintf("%s/unique_insert_data/%s/%s", outputdir, element_group_name, insert_id)

            insert_row <- read_tsv(sprintf("%s/insert_info.tsv", insert_outputdir))
            insert_annot <- import(sprintf("%s/insert_structure_annot.gtf", insert_outputdir))
            rm_annot <- import(sprintf("%s/rm_%s.good.gtf", insert_outputdir, insert_id))
            rm_annot
            rm_intersect_insert <- pintersect(rm_annot, insert_annot[insert_annot$name == "insert"])
            likely_insert <- rm_intersect_insert[width(rm_intersect_insert) == max(width(rm_intersect_insert))]
            extended_insert <- resize(insert_annot[insert_annot$name == "insert"],
                width = width(insert_annot[insert_annot$name == "insert"]) + 600, fix = "center"
            )
            rms_near_insert <- subsetByOverlaps(rm_annot, extended_insert)
            if (grepl("Alu", insert_row$Family, ignore.case = TRUE)) {
                if (any(table(grep(likely_insert$name, rms_near_insert$name, ignore.case = TRUE, value = TRUE)) > 1)) {
                    insert_suspicion_level[[insert_id]] <- "suspicious"
                } else {
                    insert_suspicion_level[[insert_id]] <- "ok"
                }
            } else if (grepl("L1", insert_row$Family, ignore.case = TRUE)) {
                if (length(grep(likely_insert$name, rms_near_insert$name, ignore.case = TRUE, value = TRUE)) > 1) {
                    insert_suspicion_level[[insert_id]] <- "suspicious"
                } else {
                    insert_suspicion_level[[insert_id]] <- "ok"
                }
            } else {
                insert_suspicion_level[[insert_id]] <- "suspicious"
            }
        }
        susdf <- tibble(UUID = names(insert_suspicion_level), insert_te_adjacent_to_same_subfam_te = unlist(insert_suspicion_level))
        susdf <- susdf %>% left_join(group_frame %>% dplyr::select(UUID, TEMatch))
        write_csv(susdf, sprintf("%s/suspicion_level.csv", rm_outputdir))
    }

    get_promising_transduction <- function(sf_with_trsd, sample_or_aref) {
        print(sample_or_aref)
        promising_transductions <- list()
        if (nrow(sf_with_trsd) == 0) {
            return(tibble())
        }

        trsd_ss <- DNAStringSet(sf_with_trsd$Transduction_3p)
        names(trsd_ss) <- sf_with_trsd$UUID


        at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)

        trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]


        bl <- blast(db = sub("\\.[^.]*$", "", grep(sample_or_aref, inputs$blast_njs, value = TRUE)))
        bres <- tibble(predict(bl, trsd_ss_for_blast))
        if (nrow(bres) == 0) {
            return(tibble())
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

        hitlist <- list()
        for (qseqid in bres_hits_other_loc$qseqid) {
            hitgrs <- bres_hits_other_loc[bres_hits_other_loc$qseqid == qseqid, ]

            relevant_subfamilies <- hitgrs %$% Subfamily %>% unique()
            relevant_subfamilies_grs <- GRanges(rmann %>% filter(grepl(paste(relevant_subfamilies, sep = "|", colapse = TRUE), rte_subfamily)))

            bresgrs <- GRanges(
                seqnames = hitgrs$sseqid,
                ranges = IRanges(
                    start = min(hitgrs$sstart, hitgrs$send),
                    end = max(hitgrs$sstart, hitgrs$send)
                )
            )
            mcols(bresgrs)$UUID <- hitgrs$qseqid
            hits <- as.data.frame(distanceToNearest(bresgrs, relevant_subfamilies_grs)) %>% tibble()
            trsd_adjacent_rte <- tibble(
                UUID = bresgrs[hits$queryHits, ]$UUID,
                gene_id = relevant_subfamilies_grs[hits$subjectHits, ]$gene_id,
                distance = hits$distance
            ) %>% mutate(adjacent = ifelse(distance < 100, TRUE, FALSE))

            hitlist[[qseqid]] <- trsd_adjacent_rte
        }
        if (length(hitlist) != 0) {
            transductions <- purrr::reduce(hitlist, bind_rows) %>%
                full_join(tibble(UUID = bres_hits_other_loc$qseqid)) %>%
                left_join(rmann) %>%
                relocate(UUID)
        } else {
            transductions <- tibble()
        }
    }


    plot_group_characteristics <- function(insert_df, group_name) {
        p <- insert_df %>%
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
        mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/bar.pdf", outputdir, group_name, "all"))

        for (element_type in insert_df %$% Subfamily %>% unique()) {
            dftemp <- insert_df %>% filter(Subfamily == element_type)
            p <- dftemp %>%
                gghistogram(x = "LengthIns", fill = "blue") +
                mtopen +
                labs(title = element_type)
            mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/length_distribution.pdf", outputdir, group_name, element_type))

            p <- dftemp %>%
                mutate(tsdlen = nchar(TSD)) %>%
                gghistogram(x = "tsdlen", fill = "blue") +
                labs(title = element_type) + mtopen
            mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/tsd_length_distribution.pdf", outputdir, group_name, element_type))

            p <- dftemp %>%
                mutate(trsd5 = nchar(Transduction_5p)) %>%
                gghistogram(x = "trsd5", fill = "blue") +
                labs(title = element_type) + mtopen
            mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/trsd5_length_distribution.pdf", outputdir, group_name, element_type))

            p <- dftemp %>%
                mutate(trsd3 = nchar(Transduction_3p)) %>%
                gghistogram(x = "trsd3", fill = "blue") +
                labs(title = element_type) + mtopen
            mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/trsd3_length_distribution.pdf", outputdir, group_name, element_type))

            p <- dftemp %>%
                arrange(-EndTE) %>%
                mutate(nrow = row_number()) %>%
                ggplot() +
                geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
                mtclosed +
                labs(title = element_type, x = "Consensus Position")
            mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/te_body_distribution.pdf", outputdir, group_name, element_type))

            if (element_type == "AluY") {
                p <- dftemp %>%
                    arrange(-EndTE) %>%
                    mutate(nrow = row_number()) %>%
                    ggplot() +
                    geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
                    geom_vline(xintercept = 250) +
                    mtclosed +
                    labs(title = element_type, x = "Consensus Position")
                mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/te_body_distribution_withfilterguide.pdf", outputdir, group_name, element_type))
            } else if (element_type == "L1HS") {
                p <- dftemp %>%
                    arrange(-EndTE) %>%
                    mutate(nrow = row_number()) %>%
                    ggplot() +
                    geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
                    geom_vline(xintercept = 5800) +
                    mtclosed +
                    labs(title = element_type, x = "Consensus Position")
                mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/%s/te_body_distribution_withfilterguide.pdf", outputdir, group_name, element_type))
            }
        }
    }
}


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



sample_table %>%
    left_join(sample_sequencing_data) %>%
    dplyr::select(-condition, -nanopore_rawdata_dir) %>%
    mutate(coverage = bases_number / (3.1 * 10**9)) %>%
    write_delim(sprintf("aref/results/sample_chars.txt"))
p <- sample_table %>%
    left_join(sample_sequencing_data) %>%
    dplyr::select(-condition, -nanopore_rawdata_dir) %>%
    mutate(coverage = signif(bases_number / (3.1 * 10**9), 3)) %>%
    mutate(across(
        c(reads_number, N50, bases_number),
        ~ formatC(as.numeric(.x), format = "e", digits = 2)
    )) %>%
    ggtexttable(theme = ttheme("minimal"))
mysaveandstore(pl = p, sprintf("aref/results/sample_chars.pdf"), 8, 5)
p <- sample_table %>%
    left_join(sample_sequencing_data) %>%
    dplyr::select(-nanopore_rawdata_dir) %>%
    group_by(condition) %>%
    summarise(mean_bases = mean(bases_number), mean_N50 = mean(N50), mean_reads = mean(reads_number)) %>%
    mutate(mean_coverage = signif(mean_bases / (3.1 * 10**9), 3)) %>%
    mutate(across(
        c(mean_bases, mean_reads, mean_N50),
        ~ formatC(as.numeric(.x), format = "e", digits = 2)
    )) %>%
    ggtexttable(theme = ttheme("minimal"))
mysaveandstore(pl = p, sprintf("aref/results/sample_chars_grouped_means.pdf"), 6, 3)
p <- sample_table %>%
    left_join(sample_sequencing_data) %>%
    dplyr::select(-nanopore_rawdata_dir) %>%
    summarise(mean_bases = mean(bases_number), mean_N50 = mean(N50), mean_reads = mean(reads_number)) %>%
    mutate(mean_coverage = signif(mean_bases / (3.1 * 10**9), 3)) %>%
    mutate(across(
        c(mean_bases, mean_reads, mean_N50),
        ~ formatC(as.numeric(.x), format = "e", digits = 2)
    )) %>%
    ggtexttable(theme = ttheme("minimal"))
mysaveandstore(pl = p, sprintf("aref/results/sample_chars_means.pdf"), 6, 3)
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
} else {
    dfall <- read.table("aref/A.REF_tldr/A.REF.table.txt", header = TRUE) %>%
        mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
        tibble()
}
somatic_naught <- dfall %>%
    filter(!is.na(EmptyReads)) %>%
    rowwise() %>%
    mutate(emptyreadsnum = sum(as.numeric(gsub("\\|", "", unlist(str_extract_all(EmptyReads, pattern = "\\|[0-9]+")))))) %>%
    ungroup() %>%
    filter(MedianMapQ >= 60) %>%
    filter(Filter == "PASS") %>%
    mutate(sample_name = gsub("\\..*", "", str_extract(SampleReads, paste(conf$samples, collapse = "|")))) %>%
    filter(SpanReads > 0) %>%
    filter(SpanReads < 6) %>%
    mutate(TSD_OK = ifelse(is.na(TSD), FALSE, ifelse(nchar(TSD) < 21, TRUE, FALSE))) %>%
    # binomial probability that we observe this few reads were it a heterozygous insert
    mutate(prob_fp = ifelse(UsedReads > as.numeric(emptyreadsnum), 1, pbinom(q = UsedReads, size = UsedReads + as.numeric(emptyreadsnum), prob = 0.5))) %>%
    filter(prob_fp < 0.001) %>%
    filter(is.na(NonRef)) %>%
    relocate(prob_fp, UsedReads, emptyreadsnum)

dfall %>%
    filter(Filter == "PASS") %>%
    filter(!is.na(EmptyReads)) %>%
    rowwise() %>%
    mutate(emptyreadsnum = sum(as.numeric(gsub("\\|", "", unlist(str_extract_all(EmptyReads, pattern = "\\|[0-9]+")))))) %>%
    ungroup() %>%
    mutate(sample_name = gsub("\\..*", "", str_extract(SampleReads, paste(conf$samples, collapse = "|")))) %>%
    mutate(prob_fp = ifelse(UsedReads > as.numeric(emptyreadsnum), 1, dbinom(x = UsedReads, size = UsedReads + as.numeric(emptyreadsnum), prob = 0.5))) %>%
    filter(prob_fp < 0.001) %>%
    relocate(prob_fp, UsedReads, emptyreadsnum)

somatic_alpha <- dfall %>%
    filter(!is.na(EmptyReads)) %>%
    rowwise() %>%
    mutate(emptyreadsnum = sum(as.numeric(gsub("\\|", "", unlist(str_extract_all(EmptyReads, pattern = "\\|[0-9]+")))))) %>%
    ungroup() %>%
    filter(MedianMapQ >= 60) %>%
    filter(Filter == "PASS") %>%
    filter(!is.na(TSD)) %>%
    mutate(sample_name = gsub("\\..*", "", str_extract(SampleReads, paste(conf$samples, collapse = "|")))) %>%
    filter(SpanReads > 0) %>%
    filter(SpanReads < 6) %>%
    mutate(TSD_OK = ifelse(nchar(TSD) < 21, TRUE, FALSE)) %>%
    # binomial probability that we observe this few reads were it a heterozygous insert
    mutate(prob_fp = ifelse(UsedReads > as.numeric(emptyreadsnum), 1, dbinom(x = UsedReads, size = UsedReads + as.numeric(emptyreadsnum), prob = 0.5))) %>%
    filter(prob_fp < 0.001) %>%
    relocate(prob_fp, UsedReads, emptyreadsnum)

somatic_alpha %>%
    group_by(UsedReads, SpanReads, TSD_OK) %>%
    summarize(n = n())

somatic_alpha_annotated <- somatic_alpha %>%
    annotate_mappability() %>%
    annotate_read_metadata() %>%
    annotate_teend()

somatic_alpha_annotated %$% insert_fail_reason_detail
somatic_alpha_annotated %$% k50_mappable %>% table()
somatic_alpha_annotated %$% insert_mean_mapqs %>% quantile()
somatic_alpha_annotated %$% insert_supplementary_status %>% table()
somatic_alpha_annotated %$% insert_indel_and_clipping_status %>% table()
somatic_alpha_annotated %$% insert_status %>% table()
somatic_alpha_annotated %$% deletion_status %>% table()
somatic_alpha_annotated %$% clip_status %>% table()
somatic_alpha_annotated %$% inrepregion %>% table()


somatic_alpha_annotated %$% inrepregion %>% table()



somatic_mappable %$% inrepregion %>% table()
somatic_mappable_less_stringent %$% inrepregion %>% table()

somatic_readfiltered <- somatic_alpha_annotated %>%
    filter(insert_supplementary_status == TRUE) %>%
    filter(insert_indel_and_clipping_and_mapq_status == TRUE)

somatic_readfiltered_teendfiltered <- somatic_readfiltered %>%
    filter(EndTE >= endte_manual_filter | is.na(endte_manual_filter)) %>%
    filter(TEMatch >= 90)

somatic_readfiltered_l1hs_extended <- somatic_readfiltered %>%
    filter(EndTE >= 5600 & EndTE < 5800) %>%
    filter(TEMatch >= 90)

somatic_mappable_readfiltered_teendfiltered <- somatic_readfiltered_teendfiltered %>%
    filter(is.na(NonRef)) %>%
    filter(insert_mean_mapqs > 55) %>%
    filter(k50_mappable == TRUE) %>%
    mutate(mappability_stringency = "high")
somatic_mappable_less_stringentreadfiltered_teendfiltered <- somatic_readfiltered_teendfiltered %>%
    filter(is.na(NonRef)) %>%
    filter((insert_mean_mapqs > 40 & insert_mean_mapqs <= 55) | (k50_mappable == FALSE & insert_mean_mapqs > 40)) %>%
    mutate(mappability_stringency = "medium")


somatic_mappable_readfiltered_teendfiltered_l1hs_extended <- somatic_readfiltered_l1hs_extended %>%
    filter(is.na(NonRef)) %>%
    filter(insert_mean_mapqs > 55) %>%
    filter(k50_mappable == TRUE) %>%
    mutate(mappability_stringency = "high")
somatic_mappable_less_stringentreadfiltered_teendfiltered_l1hs_extended <- somatic_readfiltered_l1hs_extended %>%
    filter(is.na(NonRef)) %>%
    filter((insert_mean_mapqs > 40 & insert_mean_mapqs <= 55) | (k50_mappable == FALSE & insert_mean_mapqs > 40)) %>%
    mutate(mappability_stringency = "medium")

sdf <- rbind(somatic_mappable_readfiltered_teendfiltered, somatic_mappable_less_stringentreadfiltered_teendfiltered)
sdf_l1hs_extended <- rbind(somatic_mappable_readfiltered_teendfiltered_l1hs_extended, somatic_mappable_less_stringentreadfiltered_teendfiltered_l1hs_extended)


sdf %$% inrepregion %>% table()
sdf_l1hs_extended %$% inrepregion %>% table()

plot_group_characteristics(somatic_naught, "no_filter")
plot_group_characteristics(somatic_alpha, "has_tsd")
plot_group_characteristics(somatic_readfiltered, "has_tsd_is_readfiltered")
plot_group_characteristics(somatic_readfiltered_teendfiltered, "has_tsd_is_readfiltered_teendfiltered")
plot_group_characteristics(somatic_mappable_readfiltered_teendfiltered, "has_tsd_is_mappable_readfiltered_teendfiltered")
plot_group_characteristics(somatic_mappable_less_stringentreadfiltered_teendfiltered, "has_tsd_is_mappablelowerstringency_readfiltered_teendfiltered")
plot_group_characteristics(somatic_mappable_readfiltered_teendfiltered_l1hs_extended, "has_tsd_is_mappable_readfiltered_teendfiltered_l1hs_extended")
plot_group_characteristics(somatic_mappable_less_stringentreadfiltered_teendfiltered_l1hs_extended, "has_tsd_is_mappablelowerstringency_readfiltered_teendfiltered_l1hs_extended")



insert_frames <- sdf %>%
    mutate(tsd_filter = ifelse(TSD_OK, "tsd_pass", "tsd_fail")) %>%
    mutate(sup_read_combination = paste0(mappability_stringency, "_", tsd_filter, "_", SpanReads, ".", UsedReads)) %>%
    split(.$sup_read_combination)
names(insert_frames)

imap(insert_frames, ~ analyze_inserts(.x, .y))

insert_frames_l1hs_extended <- sdf_l1hs_extended %>%
    mutate(tsd_filter = ifelse(TSD_OK, "tsd_pass", "tsd_fail")) %>%
    mutate(sup_read_combination = paste0("l1hs_extended", "_", mappability_stringency, "_", tsd_filter, "_", SpanReads, ".", UsedReads)) %>%
    split(.$sup_read_combination)
names(insert_frames_l1hs_extended)

imap(insert_frames_l1hs_extended, ~ analyze_inserts(.x, .y))

# insert_frames <- list(
#     "germline" = germline_insert_df[1:20, ] %>% filter(Family == "ALU")
# )

# imap(insert_frames, ~ analyze_inserts(.x, .y))


#### transduction mapping
if (conf$update_ref_with_tldr$per_sample == "yes") {
    transductions_list <- list()
    for (sample in conf$samples) {
        sf_with_trsd <- sdf %>%
            filter(sample_name == sample) %>%
            filter(nchar(Transduction_3p) > 10)
        transductions_list[[sample]] <- get_promising_transduction(sf_with_trsd, sample) %>%
            mutate(sample_name = sample)
    }
    transductions_df <- purrr::reduce(transductions_list, bind_rows) %>%
        group_by(UUID) %>%
        filter(row_number() == 1)
} else {
    sf_with_trsd <- sdf %>%
        filter(nchar(Transduction_3p) > 10)
    transductions_df <- get_promising_transduction(sf_with_trsd, "A.REF") %>%
        mutate(sample_name = "A.REF") %>%
        group_by(UUID) %>%
        filter(row_number() == 1)
}
write_csv(transductions_df %>% filter(adjacent == TRUE), sprintf("%s/promising_transductions.csv", outputdir))


# NOW WITH ALL THIS DATA PAUSE - MANUALLY INSPECT EACH PUTATIVE INSERT, AND WRITE UUIDS OF PASSING ELEMENTS INTO THE CURATED_ELEMENTS_PATH FILE BELOW
# all elements that passed scrutiny:
curated_elements_path <- sprintf("%s/insert_characteristics/curated_inserts.txt", outputdir)
if (file.exists(curated_elements_path)) {
    ce <- read_csv(curated_elements_path, col_names = FALSE)
    curation_df_complete <- tibble(UUID = ce$X1, pass_curation = TRUE)

    for (UUID in curation_df_complete$UUID) {
        path <- system(sprintf("find aref/results/somatic_insertions/unique_insert_data -name %s", UUID), intern = TRUE)
        system("mkdir -p aref/results/somatic_insertions/insert_characteristics/curated_elements/insert_data")
        system(sprintf("cp -r %s aref/results/somatic_insertions/insert_characteristics/curated_elements/insert_data", path))
        file.exists(path)
    }


    if (is.null(dfall$condition)) {
        dfall$condition <- conf$levels
    }
    pass <- somatic_alpha_annotated %>%
        left_join(curation_df_complete) %>%
        filter(pass_curation == TRUE) %>%
        mutate(sample_name = factor(sample_name, levels = conf$samples)) %>%
        mutate(condition = factor(condition, levels = conf$levels))

    pass %>% arrange(UUID)

    group_name <- "curated_elements"
    plot_group_characteristics(pass, group_name)
    library(ggrepel)

    pf <- pass %>%
        group_by(sample_name, condition, apoe, sex, age, Subfamily) %>%
        summarise(n = n()) %>%
        full_join(sample_table) %>%
        mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n)) %>%
        ungroup()
    p <- pass %>%
        group_by(sample_name, condition, apoe, sex, age) %>%
        summarise(n = n()) %>%
        full_join(sample_table) %>%
        mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n)) %>%
        ungroup() %>%
        mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
        ggplot(aes(y = sample_name, x = n, color = condition, shape = sex)) +
        geom_segment(aes(x = 0, xend = n, yend = sample_name), color = "green") +
        geom_point(size = 3) +
        scale_conditions +
        # geom_vline(xintercept = pf %>% filter(condition == "AD") %$% n %>% mean(), color = "blue", linetype = "dashed") +
        # geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% n %>% mean(), color = "grey", linetype = "dashed") +
        geom_text_repel(aes(label = apoe)) +
        # new_scale_fill() +
        # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
        scale_conditions +
        mtopen
    library(broom)
    stats <- summary(lm(n ~ condition + sex + age, pf)) %>% tidy()
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/point_char.pdf", outputdir, group_name), 5, 4, sf = stats)

    p <- pass %>%
        group_by(sample_name, condition, apoe, sex, age) %>%
        summarise(n = n()) %>%
        full_join(sample_table) %>%
        mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n)) %>%
        ungroup() %>%
        mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
        ggplot(aes(x = condition, y = n)) +
        geom_boxplot(aes(color = condition)) +
        scale_conditions +
        geom_point(size = 3) +
        mtopen
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/point_boxplot_char.pdf", outputdir, group_name), 5, 4, sf = stats)


    pf <- pass %>%
        group_by(sample_name, condition, apoe, sex, age, Subfamily) %>%
        summarise(n = n()) %>%
        full_join(sample_table %>% crossing(Subfamily = c("L1HS", "AluY"))) %>%
        mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n)) %>%
        ungroup()
    p <- pf %>%
        mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
        ggplot(aes(y = sample_name, x = n, color = condition, shape = sex)) +
        geom_segment(aes(x = 0, xend = n, yend = sample_name), color = "green") +
        geom_point(size = 3) +
        scale_conditions +
        # geom_vline(xintercept = pf %>% filter(condition == "AD") %$% n %>% mean(), color = "blue", linetype = "dashed") +
        # geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% n %>% mean(), color = "grey", linetype = "dashed") +
        geom_text_repel(aes(label = apoe)) +
        facet_wrap(~Subfamily) +
        # new_scale_fill() +
        # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
        scale_conditions +
        mtclosed
    library(broom)
    statsAlu <- summary(lm(n ~ condition + sex + age, pf %>% filter(Subfamily == "AluY"))) %>%
        tidy() %>%
        mutate(Subfamily = "AluY")
    statsL1 <- summary(lm(n ~ condition + sex + age, pf %>% filter(Subfamily == "L1HS"))) %>%
        tidy() %>%
        mutate(Subfamily = "L1HS")
    stats <- full_join(statsAlu, statsL1)
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/point_char_faceted.pdf", outputdir, group_name), 5, 4, sf = stats)

    p <- pass %>%
        ggplot(aes(x = sample_name, fill = condition)) +
        geom_bar() +
        facet_wrap(~Family) +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtclosed +
        anchorbar +
        scale_conditions +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar_family.pdf", outputdir, group_name))

    p <- pass %>%
        ggplot(aes(x = sample_name, fill = Subfamily)) +
        geom_bar() +
        labs(x = "", title = "RTE Somatic Insertions") +
        mtclosed +
        anchorbar +
        scale_palette +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar_limited.pdf", outputdir, group_name))

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
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar1.pdf", outputdir, group_name))

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
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar_age.pdf", outputdir, group_name))

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
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar_braak.pdf", outputdir, group_name))

    p <- pass %>%
        group_by(sample_name) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        complete(sample_name, fill = list(n = 0)) %>%
        left_join(sample_table) %>%
        left_join(sample_sequencing_data) %>%
        ggplot(aes(x = sample_name)) +
        geom_point(aes(x = bases_number, y = n), size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "Bases Sequenced", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette +
        scale_y_continuous(limits = c(-0.5, NA)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/bar_bases.pdf", outputdir, group_name))



    pass %>%
        group_by(Subfamily) %>%
        summary()
    library(dplyr)
    library(tidyr)
    data <- pass %>%
        filter(Subfamily == "AluY") %>%
        dplyr::select(LengthIns, StartTE, EndTE)
    p <- as.data.frame(sapply(data, summary))
    df <- as.data.frame(lapply(p, function(x) if (is.numeric(x)) round(x, 0) else x))

    cut_site_df <- read_csv("aref/results/somatic_insertions/insert_characteristics/curated_elements/curated_element_cutsite_region.csv")
    p <- pass %>%
        left_join(cut_site_df) %>%
        mutate(tsd_length = nchar(TSD)) %>%
        dplyr::select(Subfamily, sample_name, sample_name, seqnames, LengthIns, start, StartTE, EndTE, tsd_length, UsedReads, emptyreadsnum, `cut_site_seq_-3_+7`) %>%
        dplyr::rename(`Used Reads` = UsedReads) %>%
        dplyr::rename(`Empty Reads` = emptyreadsnum) %>%
        arrange(sample_name, LengthIns) %>%
        ggtexttable(theme = ttheme("minimal"))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/insert_chars.pdf", outputdir, group_name), 12, 5)

    p <- pass %>%
        left_join(cut_site_df) %>%
        mutate(tsd_length = nchar(TSD)) %>%
        dplyr::select(Subfamily, sample_name, seqnames, start, LengthIns, StartTE, EndTE, tsd_length, UsedReads, emptyreadsnum, `cut_site_seq_-3_+7`, TEMatch, insert_mean_mapqs, MedianMapQ) %>%
        dplyr::rename(`Used Reads` = UsedReads) %>%
        dplyr::rename(`Empty Reads` = emptyreadsnum) %>%
        arrange(sample_name, LengthIns) %>%
        ggtexttable(theme = ttheme("minimal"))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/insert_chars_detailed.pdf", outputdir, group_name), 15, 5)

    p <- pass %>%
        left_join(cut_site_df) %>%
        mutate(tsd_length = nchar(TSD)) %>%
        dplyr::select(UUID, Subfamily, sample_name, LengthIns, StartTE, EndTE, tsd_length, UsedReads, emptyreadsnum, `cut_site_seq_-3_+7`, TEMatch, insert_mean_mapqs, MedianMapQ) %>%
        dplyr::rename(`Used Reads` = UsedReads) %>%
        dplyr::rename(`Empty Reads` = emptyreadsnum) %>%
        arrange(sample_name, LengthIns) %>%
        ggtexttable(theme = ttheme("minimal"))
    mysaveandstore(pl = p, sprintf("%s/insert_characteristics/%s/insert_chars_detailed1.pdf", outputdir, group_name), 18, 5)
    p1 <- pass %>%
        left_join(cut_site_df) %>%
        mutate(tsd_length = nchar(TSD)) %>%
        dplyr::select(Subfamily, LengthIns, StartTE, EndTE, tsd_length) %>%
        arrange(StartTE) %>%
        ggtexttable(theme = ttheme("minimal"))
    mysaveandstore(pl = p1, sprintf("%s/insert_characteristics/%s/insert_chars_1.pdf", outputdir, group_name), 6, 4)
    p2 <- pass %>%
        left_join(cut_site_df) %>%
        mutate(tsd_length = nchar(TSD)) %>%
        dplyr::select(StartTE, UsedReads, emptyreadsnum, `cut_site_seq_-3_+7`) %>%
        dplyr::rename(`Used Reads` = UsedReads) %>%
        dplyr::rename(`Empty Reads` = emptyreadsnum) %>%
        arrange(StartTE) %>%
        dplyr::select(-StartTE) %>%
        ggtexttable(theme = ttheme("minimal"))
    mysaveandstore(pl = p2, sprintf("%s/insert_characteristics/%s/insert_chars_2.pdf", outputdir, group_name), 6, 4)
}
tdf <- pass %>%
    group_by(sample_name) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    complete(sample_name, fill = list(n = 0)) %>%
    left_join(sample_table) %>%
    left_join(sample_sequencing_data)

fa
tdf$bases_number / (3 * 10**9)
haploid_genome_length <- seqinfo(fa) %>%
    data.frame() %$% seqlengths %>%
    sum()
tdf %>% mutate(fraction_genomes_with_a_somatic_insert = n / (bases_number / (2 * haploid_genome_length)))

########### OLD CODE

# # germline
# for (insert_id in germline$UUID) {
#     insert_row <- germline %>% filter(UUID == insert_id)
#     bampath <- inputs$bam[[1]]
#     whichgr <- GRanges(germline %>% filter(UUID == insert_id))
#     flag <- scanBamFlag(
#         isUnmappedQuery = NA,
#         isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
#         isDuplicate = NA, isSupplementaryAlignment = NA
#     )
#     aln1 <- readGAlignments(bampath, param = ScanBamParam(which = whichgr, what = scanBamWhat(), tag = c("SA"), flag = flag))
#     alndf <- as.data.frame(aln1) %>%
#         tibble()
#     # insert fails if it was captured in a read which has supplementary alignments
#     if (conf$update_ref_with_tldr$per_sample == "yes") {
#         tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", sample, sample, insert_id))
#     } else {
#         tldr_df <- read_delim(sprintf("aref/%s_tldr/%s/%s.detail.out", "A.REF", "A.REF", insert_id))
#     }
#     read_name <- tldr_df %>%
#         filter(Useable == TRUE & IsSpanRead == TRUE) %$% ReadName %>%
#         pluck(1)
#     read_of_interest <- alndf %>% filter(qname == read_name)
#     read_length <- read_of_interest$seq %>% nchar()
#     cigar <- read_of_interest$cigar

#     inserts <- str_extract_all(read_of_interest$cigar, "[0-9]+I") %>%
#         unlist() %>%
#         gsub("I", "", .) %>%
#         as.numeric()
#     insert_start_in_read <- str_split(cigar, paste0(inserts[inserts > 50], "I"))[[1]][1] %>%
#         str_split(., pattern = "[a-zA-Z]") %>%
#         unlist() %>%
#         as.numeric() %>%
#         sum(na.rm = TRUE)

#     read_seq_sense <- DNAStringSet(read_of_interest$seq)
#     read_seq_sense <- subseq(read_seq_sense, start = insert_start_in_read - 600, end = insert_start_in_read + 800)
#     tryCatch(
#         {
#             if (insert_row$strand == "+") {
#                 insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2))
#             } else {
#                 insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2)) %>% reverseComplement()
#             }
#         },
#         error = function(e) {
#             if (insert_row$Strand == "+") {
#                 insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2))
#             } else {
#                 insert_site_sense <- getSeq(fa, whichgr + round((width(read_seq_sense) - width(whichgr)) / 2)) %>% reverseComplement()
#             }
#         }
#     )

#     path <- "temp_read.fa"
#     writeXStringSet(read_seq_sense, path)
#     makeblastdb(path)
#     bl <- blast(db = path)
#     bres <- tibble(predict(bl, insert_site_sense))
#     insert_site_split <- unlist(strsplit(unname(as.character(insert_site_sense)), split = ""))
#     read_seq_split <- unlist(strsplit(unname(as.character(read_seq_sense)), split = ""))

#     plotpath <- sprintf("%s/germline_dotplots/%s_insertlen_%s_%s.pdf", outputdir, insert_row$Subfamily, insert_row$LengthIns, gsub("-", "_", insert_id))
#     dir.create(dirname(plotpath), recursive = TRUE)
#     pdf(plotpath)
#     dotPlot(read_seq_split, insert_site_split, wsize = 10, nmatch = 9)
#     dev.off()
# } # TODO STILL HAVE AN ISSUE WITH THESE DOT PLOTS.

# tldr_te_ref_path <- conf$update_ref_with_tldr$tldr_te_ref[[conf$species]]

# my_tsd_check <- tibble(sample_name = character(), UUID = character(), my_tsd = character(), my_tsd_left_start = numeric(), my_tsd_right_start = numeric(), my_te_start = numeric(), my_te_end = numeric())
# for (sample in conf$samples) {
#     print(sample)
#     bampath <- grep(sample, inputs$bam, value = TRUE)
#     insert_ids <- somatic %>% filter(sample_name == sample) %$% UUID
#     print(length(insert_ids))
#     for (insert_id in insert_ids) {
#         insert_row <- somatic %>%
#             filter(sample_name == sample) %>%
#             filter(UUID == insert_id)
#         insert_loc <- GRanges(insert_row)
#         css <- DNAStringSet(insert_row$Consensus)

#         tryCatch(
#             {
#                 bl <- blast(db = tldr_te_ref_path)
#                 bres <- tibble(predict(bl, css))
#                 hit <- bres %>%
#                     filter(grepl(insert_row$Subfamily, sseqid)) %>%
#                     mutate(length_dif_from_insert = abs(length - insert_row$LengthIns)) %>%
#                     arrange(length_dif_from_insert) %>%
#                     head(1)

#                 # I scan the 220 bp region flanking the predicted insert for TSD. Will fail if there is a large five or three prime transduction
#                 five_start <- hit$qstart - 200
#                 five_end <- hit$qstart + 20
#                 three_start <- hit$qend - 20
#                 three_end <- hit$qend + 200

#                 five_string <- str_sub(as.character(css), start = five_start, end = five_end)
#                 three_string <- str_sub(as.character(css), start = three_start, end = three_end)
#                 if (is.null(five_string) | is.null(three_string)) {
#                     print("empty strings")
#                 } else {
#                     mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1000, baseOnly = TRUE)
#                     paln <- pairwiseAlignment(DNAStringSet(five_string), DNAStringSet(three_string), type = "local", substitutionMatrix = mat, gapOpening = 1000)
#                     tsd_loc <- str_locate_all(as.character(css), as.character(paln@subject))[[1]]
#                     row <- tibble(sample_name = sample, UUID = insert_id, my_tsd = as.character(paln@pattern), my_tsd_left_start = as.data.frame(tsd_loc)$start[1], my_tsd_right_start = as.data.frame(tsd_loc)$start[2], my_te_start = hit$qstart, my_te_end = hit$qend)
#                     my_tsd_check <- bind_rows(my_tsd_check, row)
#                 }
#                 # write(paste0(">five\n", five_string), "five.fa")
#                 # makeblastdb("five.fa")
#                 # bl <- blast(db = "five.fa")
#                 # bres <- tibble(predict(bl, DNAStringSet(three_string)))
#             },
#             error = function(e) {
#                 "issue"
#             }
#         )
#     }
# }


### KEEP WORKING HERE - ALSO analyze the tsd_pass_othernumber


p <- pass %>%
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

tryCatch(
    {
        metadata_vars <- colnames(sample_table)[!(colnames(sample_table) %in% c("sample_name", "condition"))]
        sequencing_metadata_vars <- c("N50", "reads_number", "bases_number")


        grouping_vars <- c("sample_name", "condition", sequencing_metadata_vars, metadata_vars)
        germline_n <- germline %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            full_join(sample_table %>% full_join(sample_sequencing_data)) %>%
            replace_na(list(nins = 0)) %>%
            mutate(Insert_Type = "Germline")
        somatic_n <- pass %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            full_join(sample_table %>% full_join(sample_sequencing_data)) %>%
            replace_na(list(nins = 0)) %>%
            mutate(Insert_Type = "Somatic")
        insert_n <- bind_rows(germline_n, somatic_n)
        p <- insert_n %>% ggplot(aes(x = N50, y = nins, color = sample_name)) +
            geom_point(size = 2) +
            facet_wrap(~Insert_Type, scales = "free_y") +
            labs(x = "N50", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/n50_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- insert_n %>% ggplot(aes(x = reads_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~Insert_Type, scales = "free_y") +
            labs(x = "Total Reads Count", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/reads_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- insert_n %>% ggplot(aes(x = bases_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~Insert_Type, scales = "free_y") +
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
            complete(sample_name, Subfamily, fill = list(nins = 0)) %>%
            dplyr::select(sample_name, Subfamily, nins) %>%
            full_join(sample_table %>% full_join(sample_sequencing_data)) %>%
            mutate(Insert_Type = "Germline")
        somatic_n_bysubfam <- pass %>%
            group_by(across(c(grouping_vars, Subfamily))) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            complete(sample_name, Subfamily, fill = list(nins = 0)) %>%
            dplyr::select(sample_name, Subfamily, nins) %>%
            full_join(sample_table %>% full_join(sample_sequencing_data)) %>%
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


pass %>%
    filter(Subfamily == "L1HS") %>%
    print(width = Inf)
pass_out
# ctrl6




# # ADDTIONAL FILTERS
# # I think these are overly restrictive - not applying these
# germline_insert_characteristics <- tibble(germline_tsd_95 = numeric(), germline_trsd3_95 = numeric(), germline_trsd5_95 = numeric(), germline_endte_05 = numeric(), Subfamily = character())
# for (element_type in somatic %$% Subfamily %>% unique()) {
#     germlinetemp <- germline %>% filter(Subfamily == element_type)
#     germline_tsd_95 <- germlinetemp %$% TSD %>%
#         nchar() %>%
#         replace_na(0) %>%
#         quantile(0.95, na.rm = TRUE)
#     germline_trsd3_95 <- germlinetemp %$% Transduction_3p %>%
#         nchar() %>%
#         replace_na(0) %>%
#         quantile(0.95, na.rm = TRUE)
#     germline_trsd5_95 <- germlinetemp %$% Transduction_5p %>%
#         nchar() %>%
#         replace_na(0) %>%
#         quantile(0.95, na.rm = TRUE)
#     germline_endte_05 <- germlinetemp %$% EndTE %>%
#         quantile(0.05, na.rm = TRUE) %>%
#         unname()
#     germline_insert_characteristics <- germline_insert_characteristics %>% add_row(germline_tsd_95 = germline_tsd_95, germline_trsd3_95 = germline_trsd3_95, germline_trsd5_95 = germline_trsd5_95, germline_endte_05 = germline_endte_05, Subfamily = element_type)
# }



rmann_nr_list <- list()
for (sample in sample_table$sample_name) {
    rmann_nr_temp <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
    rmann_nr_temp$sample_name <- sample
    rmann_nr_list[[sample]] <- rmann_nr_temp
}

rmann_nr <- do.call(rbind, rmann_nr_list) %>%
    tibble() %>%
    mutate(gene_id = paste0(sample_name, "___", gene_id)) %>%
    mutate(seqnames = paste0(sample_name, "___", seqnames))


germline_insert_grs <- rmann_nr %$%
    seqnames %>%
    str_extract(., "chr[1-9XY].*") %>%
    str_split("_") %>%
    map(~ GRanges(tibble(
        seqnames = .x[1],
        start = .x[2],
        stop = .x[3]
    ))) %>%
    purrr::reduce(., c)

pass_grs <- pass %>%
    GRanges() %>%
    flank(width = 100, both = TRUE)

pass_grs %>% subsetByOverlaps(germline_insert_grs)



total_cov <- sample_table %>%
    left_join(sample_sequencing_data) %>%
    dplyr::select(-condition, -nanopore_rawdata_dir) %>%
    mutate(coverage = bases_number / (3.1 * 10**9)) %$% coverage %>%
    sum()

frost_total_cov <- 7.4 * 18

pass_grs %>% subsetByOverlaps(segdups)
