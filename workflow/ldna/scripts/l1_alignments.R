modcode <- "m"
module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(Rsamtools)
library(cowplot)
library(patchwork)
library(forcats)
library(parallel)

tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/extended/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/extended/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/extended/A.REF.fa"
        ), env = globalenv())
        assign("params", list(
            l13 = conf$l13fasta,
            mod_code = modcode
        ), env = globalenv())
        assign("outputs", list(), env = globalenv())
    }
)

fa <- Rsamtools::FaFile(inputs$ref)

if (confALL$aref$update_ref_with_tldr$response == "yes") {
    if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
        rmannShared <- read_csv(conf$rmann)
        rmannSamples <- list()
        for (sample in sample_table$sample_name) {
            df <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
            df$sample_name <- sample
            df <- df %>% left_join(sample_table)
            rmannSamples[[sample]] <- df
        }
        rmannnonref <- do.call(rbind, rmannSamples) %>% tibble()
        rmann <- bind_rows(rmannShared, rmannnonref) %>%
            filter(refstatus != "NonCentral")
    } else if (confALL$aref$update_ref_with_tldr$per_sample == "no") {
        rmann <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann.csv", "A.REF", "A.REF")) %>% filter(refstatus != "NonCentral")
    }
} else {
    rmann <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann.csv", "A.REF", "A.REF")) %>% filter(refstatus != "NonCentral")
}

outputdir <- sprintf("ldna/results/%s/plots/l1_alignment", modcode)
dir.create(sprintf("%s/alignments", outputdir), recursive = TRUE)
dir.create(sprintf("%s/consensus", outputdir), recursive = TRUE)

# Helper: extract sequences for elements, handling per-sample nonref fasta files
get_element_sequences <- function(element_df, fa, confALL) {
    ref_df <- element_df %>% filter(refstatus == "Ref")
    ref_ss <- DNAStringSet()
    if (nrow(ref_df) > 0) {
        ref_gr <- GRanges(ref_df)
        ref_ss <- getSeq(fa, ref_gr)
        names(ref_ss) <- mcols(ref_gr)$gene_id
    }

    nonref_ss <- DNAStringSet()
    if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
        nonref_list <- list()
        for (sample in sample_table$sample_name) {
            sample_df <- element_df %>% filter(sample_name == sample)
            if (nrow(sample_df) > 0) {
                sample_gr <- GRanges(sample_df)
                sample_fa <- Rsamtools::FaFile(sprintf("aref/extended/%s.fa", sample))
                sample_ss <- getSeq(sample_fa, sample_gr)
                names(sample_ss) <- paste0(sample, "___", mcols(sample_gr)$gene_id)
                nonref_list[[sample]] <- sample_ss
            }
        }
        if (length(nonref_list) > 0) {
            nonref_ss <- purrr::reduce(nonref_list, c)
        }
    } else {
        nonref_df <- element_df %>% filter(refstatus == "NonRef")
        if (nrow(nonref_df) > 0) {
            nonref_gr <- GRanges(nonref_df)
            nonref_ss <- getSeq(fa, nonref_gr)
            names(nonref_ss) <- mcols(nonref_gr)$gene_id
        }
    }

    return(c(ref_ss, nonref_ss))
}

# Map pairwise alignment result to consensus positions for a single element
map_to_consensus <- function(aln_result, element_name) {
    aln_subject <- as.character(aligned(subject(aln_result)))
    aln_pattern <- as.character(aligned(pattern(aln_result)))
    subject_chars <- str_split_1(aln_subject, "")
    pattern_chars <- str_split_1(aln_pattern, "")

    consensus_pos <- integer(length(subject_chars))
    sequence_pos <- integer(length(subject_chars))
    cons_idx <- start(subject(aln_result)) - 1L
    seq_idx <- start(pattern(aln_result)) - 1L

    for (i in seq_along(subject_chars)) {
        has_cons <- subject_chars[i] != "-"
        has_seq <- pattern_chars[i] != "-"
        if (has_cons) cons_idx <- cons_idx + 1L
        if (has_seq) seq_idx <- seq_idx + 1L
        consensus_pos[i] <- if (has_cons) cons_idx else NA_integer_
        sequence_pos[i] <- if (has_seq) seq_idx else NA_integer_
    }

    tibble(
        gene_id = element_name,
        consensus_pos = consensus_pos,
        sequence_pos = sequence_pos,
        consensus_seq = subject_chars,
        seq = pattern_chars
    ) %>% filter(!(is.na(sequence_pos) & is.na(consensus_pos)))
}

###############################################################################
# STEP 1: Derive consensus sequences from FL elements via MSA (per subfamily)
###############################################################################
subfams <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")
consensus_seqs <- list()

# Set to TRUE to skip MSA and read pre-computed consensus sequences
read_precomputed_consensus <- TRUE

if (read_precomputed_consensus) {
    message("=== Reading pre-computed consensus sequences ===")
    for (subfam in subfams) {
        consensus_path <- sprintf("%s/consensus/%s_consensus.fa", outputdir, subfam)
        consensus_seqs[[subfam]] <- readDNAStringSet(consensus_path)
        message(sprintf("  %s consensus length: %d bp", subfam, width(consensus_seqs[[subfam]])))
    }
} else {
    for (subfam in subfams) {
        message(sprintf("=== %s: Building consensus from FL elements ===", subfam))
        grs_fl <- rmann %>%
            filter(rte_length_req == "FL") %>%
            filter(rte_subfamily == subfam) %>%
            group_by(rte_subfamily) %>%
            mutate(length99 = quantile(length, 0.99)) %>%
            ungroup() %>%
            filter(length < length99 * 1.05)

        fl_ss <- get_element_sequences(grs_fl, fa, confALL)
        fl_fa_path <- sprintf("%s/alignments/%s_fl.fa", outputdir, subfam)
        writeXStringSet(fl_ss, fl_fa_path)

        fl_aln_path <- sprintf("%s/alignments/%s_fl.aln.fa", outputdir, subfam)
        system(sprintf("mafft --auto %s > %s", fl_fa_path, fl_aln_path))

        aln <- readDNAMultipleAlignment(fl_aln_path)
        consensusMat <- consensusMatrix(aln)
        base_names <- rownames(consensusMat)
        max_indices <- apply(consensusMat, 2, which.max)
        consensus_vec <- base_names[max_indices]
        consensus_str <- paste(consensus_vec, collapse = "")

        consensus_nogaps <- gsub("-", "", consensus_str)
        consensus_ss <- DNAStringSet(consensus_nogaps)
        names(consensus_ss) <- sprintf("%s_consensus", subfam)
        consensus_path <- sprintf("%s/consensus/%s_consensus.fa", outputdir, subfam)
        writeXStringSet(consensus_ss, consensus_path)

        consensus_seqs[[subfam]] <- consensus_ss
        message(sprintf("  %s consensus length: %d bp", subfam, nchar(consensus_nogaps)))
    }
}

###############################################################################
# STEP 2: Align ALL elements to their subfamily consensus (pairwise)
###############################################################################
n_cores <- 40
all_consensus_maps <- list()
alignment_stats_list <- list()

align_one_element <- function(i, all_ss, consensus_seq) {
    element_ss <- all_ss[i]
    aln_result <- pairwiseAlignment(
        element_ss[[1]], consensus_seq,
        type = "local"
    )
    mapping <- map_to_consensus(aln_result, names(all_ss)[i])
    stats <- tibble(
        gene_id = names(all_ss)[i],
        element_length = width(element_ss),
        aln_score = score(aln_result),
        aln_nchar = nchar(aln_result),
        pct_identity = pid(aln_result),
        consensus_start = start(subject(aln_result)),
        consensus_end = end(subject(aln_result)),
        consensus_length = length(consensus_seq)
    )
    list(mapping = mapping, stats = stats)
}

for (subfam in subfams) {
    message(sprintf("=== %s: Pairwise alignment of all elements to consensus ===", subfam))
    consensus_ss <- consensus_seqs[[subfam]]

    grs_all <- rmann %>%
        filter(rte_subfamily == subfam)

    all_ss <- get_element_sequences(grs_all, fa, confALL)
    message(sprintf("  %d elements to align using %d cores", length(all_ss), n_cores))

    results <- mclapply(seq_along(all_ss), align_one_element,
        all_ss = all_ss,
        consensus_seq = consensus_ss[[1]],
        mc.cores = n_cores
    )

    consensus_index_long <- bind_rows(lapply(results, `[[`, "mapping"))
    stats_df <- bind_rows(lapply(results, `[[`, "stats")) %>%
        mutate(rte_subfamily = subfam)

    all_consensus_maps[[subfam]] <- consensus_index_long
    alignment_stats_list[[subfam]] <- stats_df

    write_csv(consensus_index_long, sprintf("%s/%s_all_mapping_to_consensus_table.csv", outputdir, subfam))
    message(sprintf("  wrote %s mapping table (%d rows)", subfam, nrow(consensus_index_long)))
}

alignment_stats <- bind_rows(alignment_stats_list)
alignment_stats <- alignment_stats %>%
    left_join(rmann %>% select(gene_id, rte_length_req) %>% distinct(), by = "gene_id")
write_csv(alignment_stats, sprintf("%s/alignment_stats.csv", outputdir))

###############################################################################
# STEP 3: Build harmonized positions across subfamilies
#   MSA of all subfamily consensus sequences. The alignment column index
#   becomes harmonized_pos — a unified coordinate system that captures
#   positions present in ANY subfamily, not just L1HS.
###############################################################################
message("=== Building cross-subfamily harmonized positions ===")

# Combine and align all subfamily consensus sequences
all_consensus_unaligned <- DNAStringSet(lapply(consensus_seqs[subfams], `[[`, 1))
names(all_consensus_unaligned) <- subfams
consensus_input_path <- sprintf("%s/consensus/all_subfamilies_unaligned.fa", outputdir)
consensus_aligned_path <- sprintf("%s/consensus/all_subfamilies_aligned.fa", outputdir)
writeXStringSet(all_consensus_unaligned, consensus_input_path)

system(sprintf("mafft --auto %s > %s", consensus_input_path, consensus_aligned_path))

aligned_consensuses <- readDNAStringSet(consensus_aligned_path)

# Build harmonization map from the MSA: each alignment column = one harmonized_pos
harmonization_maps <- list()
for (subfam in subfams) {
    gapped_seq <- as.character(aligned_consensuses[subfam])
    chars <- str_split_1(gapped_seq, "")
    cons_idx <- 0L
    rows <- vector("list", length(chars))
    for (i in seq_along(chars)) {
        has_base <- chars[i] != "-"
        if (has_base) cons_idx <- cons_idx + 1L
        if (has_base) {
            rows[[i]] <- tibble(
                rte_subfamily = subfam,
                consensus_pos = cons_idx,
                harmonized_pos = i
            )
        }
    }
    harmonization_maps[[subfam]] <- bind_rows(rows)
}

harmonization_df <- bind_rows(harmonization_maps)
write_csv(harmonization_df, sprintf("%s/consensus_harmonization_map.csv", outputdir))
message(sprintf("  Harmonized coordinate space: %d positions", max(harmonization_df$harmonized_pos)))

###############################################################################
# STEP 4: Join harmonized positions into the per-element mapping tables
###############################################################################
message("=== Adding harmonized_pos to element mapping tables ===")

all_results <- list()
for (subfam in subfams) {
    subfam_map <- all_consensus_maps[[subfam]] %>%
        mutate(rte_subfamily = subfam) %>%
        left_join(harmonization_df, by = c("rte_subfamily", "consensus_pos"))

    all_results[[subfam]] <- subfam_map
    write_csv(subfam_map, sprintf("%s/%s_all_mapping_to_consensus_harmonized.csv", outputdir, subfam))
    message(sprintf("  %s: %d rows with harmonized positions", subfam, nrow(subfam_map)))
}

# Combined output across all subfamilies
all_combined <- bind_rows(all_results)
write_csv(all_combined, sprintf("%s/all_subfamilies_mapping_to_consensus_harmonized.csv", outputdir))
message(sprintf("=== Done. Combined table: %d rows across %d subfamilies ===", nrow(all_combined), length(subfams)))


###############################################################################
# STEP 5: Visualizations and QC
###############################################################################
plotdir <- sprintf("%s/qc_plots", outputdir)
dir.create(plotdir, recursive = TRUE)

# --- 5A: Consensus length comparison across subfamilies ---
consensus_lengths <- tibble(
    rte_subfamily = factor(subfams, levels = subfams),
    consensus_length = sapply(consensus_seqs, width)
)
print(consensus_lengths)

p <- ggplot(consensus_lengths, aes(x = rte_subfamily, y = consensus_length)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = consensus_length), vjust = -0.3, size = 3) +
    labs(title = "Consensus sequence length by subfamily", x = "Subfamily", y = "Length (bp)") +
    mtclosed
mysaveandstore(sprintf("%s/consensus_lengths.pdf", plotdir), w = 5, h = 4)

# --- 5B: Alignment quality - percent identity distribution by FL/Trnc ---
p <- ggplot(alignment_stats, aes(x = pct_identity, fill = rte_length_req)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    labs(
        title = "Pairwise alignment percent identity to consensus",
        x = "Percent identity", y = "Count", fill = "Length status"
    ) +
    mtclosed
mysaveandstore(sprintf("%s/pct_identity_distribution.pdf", plotdir), w = 10, h = 7)

# --- 5C: Consensus coverage - what fraction of the consensus does each element span? ---
alignment_stats <- alignment_stats %>%
    mutate(
        consensus_span = consensus_end - consensus_start + 1,
        pct_consensus_covered = consensus_span / consensus_length * 100
    )

p <- ggplot(alignment_stats, aes(x = pct_consensus_covered, fill = rte_length_req)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    labs(
        title = "Consensus coverage per element",
        x = "% of consensus covered by alignment", y = "Count", fill = "Length status"
    ) +
    mtclosed
mysaveandstore(sprintf("%s/consensus_coverage_distribution.pdf", plotdir), w = 10, h = 7)

# --- 5D: Where do alignments land on the consensus? (start/end positions) ---
p <- ggplot(alignment_stats, aes(x = consensus_start, y = consensus_end, color = rte_length_req)) +
    geom_point(alpha = 0.3, size = 0.5) +
    facet_wrap(~rte_subfamily, scales = "free") +
    labs(
        title = "Alignment start vs end on consensus",
        x = "Consensus start", y = "Consensus end", color = "Length status"
    ) +
    mtclosed
mysaveandstore(sprintf("%s/alignment_start_end.pdf", plotdir), w = 10, h = 7)

# --- 5E: Element count summary ---
element_counts <- alignment_stats %>%
    group_by(rte_subfamily, rte_length_req) %>%
    summarise(n = n(), .groups = "drop")
print(element_counts)

p <- ggplot(element_counts, aes(x = rte_subfamily, y = n, fill = rte_length_req)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5) +
    labs(title = "Element counts by subfamily and length status", x = "Subfamily", y = "Count", fill = "Length status") +
    mtclosed
mysaveandstore(sprintf("%s/element_counts.pdf", plotdir), w = 6, h = 4)

# --- 5F: Mismatch rate along consensus (sliding window) ---
# Computed per-subfamily using the mapping tables
window_size <- 50
mismatch_profiles <- list()
for (subfam in subfams) {
    subfam_map <- all_consensus_maps[[subfam]] %>%
        filter(!is.na(consensus_pos), !is.na(sequence_pos)) %>%
        mutate(is_mismatch = seq != consensus_seq)

    if (nrow(subfam_map) == 0) next

    max_cons_pos <- max(subfam_map$consensus_pos, na.rm = TRUE)
    window_starts <- seq(1, max_cons_pos - window_size + 1, by = window_size %/% 2)

    window_stats <- lapply(window_starts, function(ws) {
        window_data <- subfam_map %>%
            filter(consensus_pos >= ws, consensus_pos < ws + window_size)
        tibble(
            rte_subfamily = subfam,
            window_mid = ws + window_size / 2,
            mismatch_rate = mean(window_data$is_mismatch, na.rm = TRUE),
            n_bases = nrow(window_data)
        )
    }) %>% bind_rows()
    mismatch_profiles[[subfam]] <- window_stats
}
mismatch_profile_df <- bind_rows(mismatch_profiles)

p <- ggplot(mismatch_profile_df, aes(x = window_mid, y = mismatch_rate, color = rte_subfamily)) +
    geom_line(alpha = 0.8) +
    labs(
        title = sprintf("Mismatch rate along consensus (window = %d bp)", window_size),
        x = "Consensus position", y = "Mismatch rate", color = "Subfamily"
    ) +
    mtclosedgridh
mysaveandstore(sprintf("%s/mismatch_rate_along_consensus.pdf", plotdir), w = 8, h = 4)

p <- ggplot(mismatch_profile_df, aes(x = window_mid, y = mismatch_rate)) +
    geom_line(alpha = 0.8, color = "steelblue") +
    facet_wrap(~rte_subfamily, ncol = 2) +
    labs(
        title = sprintf("Mismatch rate along consensus (window = %d bp)", window_size),
        x = "Consensus position", y = "Mismatch rate"
    ) +
    mtclosedgridh
mysaveandstore(sprintf("%s/mismatch_rate_along_consensus_faceted.pdf", plotdir), w = 8, h = 10)

# --- 5G: Harmonization QC - consensus-to-consensus alignment dot plots ---
# Show where each subfamily's consensus positions land in the harmonized space
p <- ggplot(harmonization_df, aes(x = harmonized_pos, y = consensus_pos, color = rte_subfamily)) +
    geom_point(size = 0.2, alpha = 0.5) +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    labs(
        title = "Consensus harmonization: subfamily consensus pos vs harmonized pos (MSA column)",
        x = "Harmonized position (MSA alignment column)", y = "Subfamily consensus position"
    ) +
    guides(color = "none") +
    mtclosed
mysaveandstore(sprintf("%s/harmonization_dotplot.pdf", plotdir), w = 10, h = 7)

# --- 5H: Harmonized position occupancy - which subfamilies have a base at each position ---
occupancy_df <- harmonization_df %>%
    select(rte_subfamily, harmonized_pos) %>%
    distinct() %>%
    mutate(present = 1) %>%
    complete(rte_subfamily = subfams, harmonized_pos = seq_len(max(harmonization_df$harmonized_pos)), fill = list(present = 0))

occupancy_summary <- occupancy_df %>%
    group_by(harmonized_pos) %>%
    summarise(n_subfams = sum(present), .groups = "drop")

p <- ggplot(occupancy_summary, aes(x = harmonized_pos, y = n_subfams)) +
    geom_line(alpha = 0.5, linewidth = 0.3) +
    labs(
        title = "Number of subfamilies with a base at each harmonized position",
        x = "Harmonized position", y = "Subfamilies present"
    ) +
    scale_y_continuous(breaks = 1:6) +
    mtclosedgridh
mysaveandstore(sprintf("%s/harmonized_position_occupancy.pdf", plotdir), w = 8, h = 3)

# Per-subfamily: how many consensus positions map to the harmonized space (should be 100% by construction)
harm_summary <- harmonization_df %>%
    group_by(rte_subfamily) %>%
    summarise(
        consensus_positions = max(consensus_pos),
        harmonized_positions_used = n_distinct(harmonized_pos),
        total_harmonized_space = max(harmonization_df$harmonized_pos),
        .groups = "drop"
    )
print(harm_summary)

# --- 5I: Spot check - verify FL elements cover most of the consensus ---
fl_coverage_check <- alignment_stats %>%
    filter(rte_length_req == "FL") %>%
    group_by(rte_subfamily) %>%
    summarise(
        n_fl = n(),
        median_pct_covered = median(pct_consensus_covered),
        mean_pct_identity = mean(pct_identity),
        .groups = "drop"
    )
message("=== FL element alignment QC ===")
print(fl_coverage_check)

# --- 5J: Sanity check - compare element sequence position range vs element length ---
seq_pos_range_check <- all_combined %>%
    filter(!is.na(sequence_pos)) %>%
    group_by(gene_id, rte_subfamily) %>%
    summarise(
        min_seq_pos = min(sequence_pos),
        max_seq_pos = max(sequence_pos),
        n_mapped_bases = sum(!is.na(sequence_pos)),
        .groups = "drop"
    ) %>%
    left_join(alignment_stats %>% select(gene_id, element_length, rte_length_req), by = "gene_id") %>%
    mutate(pct_element_mapped = n_mapped_bases / element_length * 100)

p <- ggplot(seq_pos_range_check, aes(x = element_length, y = n_mapped_bases, color = rte_length_req)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~rte_subfamily, scales = "free") +
    labs(
        title = "Element length vs mapped bases (should track diagonal)",
        x = "Element length (bp)", y = "Mapped bases in alignment", color = "Length status"
    ) +
    mtclosed
mysaveandstore(sprintf("%s/element_length_vs_mapped_bases.pdf", plotdir), w = 10, h = 7)

message("=== All QC plots saved to: ", plotdir, " ===")
