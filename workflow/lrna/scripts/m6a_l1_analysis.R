# ============================================================
# m6A Analysis on Young LINE-1 Elements
# Comprehensive analysis: IVT filtering, reliable site ID,
# mixed-effects modeling, and evolutionary analysis
# ============================================================

# //ANCHOR - PREAMBLE
module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)
library(ggrepel)
library(tidyHeatmap)
library(glmmTMB)
library(broom.mixed)
library(furrr)
library(Rsamtools)
library(conflicted)
conflicts_prefer(
    dplyr::filter,
    dplyr::select,
    dplyr::first,
    dplyr::setdiff,
    dplyr::intersect,
    dplyr::union
)

current_mod_name <- "m6A"
samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]
sample_table_lrna <- read_csv("conf/sample_table_lrna.csv")

# Genome setup
genome_lengths <- fasta.seqlengths(conf$reference)
chromosomesAll <- names(genome_lengths)
nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, "chrX", "chrY", nonrefchromosomes)
MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS

# Mod code mapping
MOD_CODE_MAP <- tibble::tibble(
    mod_code = c("a", "m", "17596", "17802", "19227", "19228", "69426"),
    mod_name = c("m6A", "m5C", "Am", "pseU", "Um", "Cm", "inosine")
)
MOD_CANONICAL_BASE <- c(m6A = "A", m5C = "C", Am = "A", pseU = "T", Um = "T", Cm = "C", inosine = "A")

# Young L1 subfamilies of interest
YOUNG_L1_SUBFAMILIES <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")
SUBFAM_ORDER <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")

# Load repeat annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
rm(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
rmannextended <- get_repeat_annotations(
    default_or_extended = "extended",
    keep_non_central = FALSE
)
# Load pre-computed intermediates
yl1annm <- read_delim("lrna/Rintermediates/yl1annm.fa", col_names = TRUE)

# rtedf <- read_delim(sprintf("lrna/Rintermediates/%s/rtedf.tsv", current_mod_name), col_names = TRUE)
yl1annm_alllength <- read_delim("lrna/Rintermediates/yl1annm_alllength.tsv")
# yl1annm_alllength <- yl1annm_alllength %>% dplyr::rename(rte_subfamily = rte_subfamily.x)
yl1annm %>% filter(rte_length_req == "FL")

reads_path <- sprintf("lrna/Rintermediates/%s/reads_context_all.tsv", current_mod_name)
reads <- read_delim(reads_path, col_names = TRUE)
reads %>% pw()

annot
readslong <- reads %>%
    filter(read_length > 900) %>%
    filter(sample != "N_ORF1_4")

readsvlong <- readslong %>% filter(read_length > 5000)
readsvlong_mapped <- readsvlong %>%
    left_join(yl1annm_filtered %>% dplyr::select(start = genome_pos, harmonized_pos, gene_id)) %>%
    left_join(annot) %>%
    filter(orf1_intactness == "ORF1Intact") %>%
    filter(harmonized_pos %in% reliable_sites_harmonized_pos) %>%
    mutate(is_mod = ifelse(mod_qual > 0.5, 1, 0)) %>%
    mutate(bcondition = case_when(
        grepl("ORF1", condition) ~ "ORF1",
        grepl("TOT", condition) ~ "TOT",
        TRUE ~ condition
    ))

# Per-read mean m6A
readsvlong_mapped %>%
    group_by(read_id, condition, sample) %>%
    summarise(mm = mean(is_mod)) %>%
    group_by(condition, sample) %>%
    summarise(mm = mean(mm))

outdir_vlong <- str_glue("{outdir_base}/vlong_per_site")
dir.create(outdir_vlong, recursive = TRUE, showWarnings = FALSE)

# --- 1. Per-site dotplot: %m6A for ORF1 vs TOT ---
vlong_site_stats <- readsvlong_mapped %>%
    group_by(harmonized_pos, bcondition) %>%
    summarise(pctM = mean(is_mod) * 100, n_reads = n(), .groups = "drop")

vlong_site_wide <- vlong_site_stats %>%
    pivot_wider(names_from = bcondition, values_from = c(pctM, n_reads), names_sep = "_") %>%
    filter(!is.na(pctM_ORF1) & !is.na(pctM_TOT))

p <- ggplot(vlong_site_wide, aes(
    x = pctM_TOT, y = pctM_ORF1,
    size = log10(n_reads_ORF1 + n_reads_TOT)
)) +
    geom_point(alpha = 0.7, color = "#4393C3") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = harmonized_pos), size = 2.5, max.overlaps = 20) +
    scale_size_continuous(name = "log10(total cov)", range = c(1, 6)) +
    labs(
        x = "% m6A (TOT)", y = "% m6A (ORF1)",
        title = "L1HS ORF1-intact (>5kb reads): Per-site m6A"
    ) +
    coord_equal() +
    mtclosed
mss(
    pl = p, fn = str_glue("{outdir_vlong}/per_site_m6a_orf1_vs_tot1.pdf"),
    fw = 2, fh = 2, plus_void = TRUE
)

# Bar version
vlong_site_plot <- vlong_site_stats %>%
    mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))))

p <- ggplot(vlong_site_plot, aes(x = harmonized_pos, y = pctM, fill = bcondition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = scales::comma(n_reads)),
        position = position_dodge(0.8),
        vjust = -0.5, size = 2
    ) +
    labs(
        x = "Harmonized position", y = "% m6A", fill = "Condition",
        title = "L1HS ORF1-intact (>5kb reads): Per-site m6A"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    mtclosed +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
mss(
    pl = p, fn = str_glue("{outdir_vlong}/per_site_m6a_barplot111.pdf"),
    fw = 4, fh = 1, plus_void = TRUE
)

# --- 2. Intra-read pairwise correlation heatmap ---
# Pivot to wide: one column per site, one row per read
vlong_wide <- readsvlong_mapped %>%
    dplyr::select(read_id, bcondition, harmonized_pos, is_mod) %>%
    distinct(read_id, harmonized_pos, .keep_all = TRUE) %>%
    pivot_wider(names_from = harmonized_pos, values_from = is_mod, names_prefix = "site_")

vlong_site_cols <- grep("^site_", colnames(vlong_wide), value = TRUE)
vlong_site_mat <- vlong_wide %>%
    dplyr::select(all_of(vlong_site_cols)) %>%
    as.matrix()

cor_mat_vl <- cor(vlong_site_mat, use = "pairwise.complete.obs")
colnames(cor_mat_vl) <- gsub("^site_", "", colnames(cor_mat_vl))
rownames(cor_mat_vl) <- gsub("^site_", "", rownames(cor_mat_vl))
pos_order_vl <- order(as.numeric(colnames(cor_mat_vl)))
cor_mat_vl <- cor_mat_vl[pos_order_vl, pos_order_vl]

p <- Heatmap(
    cor_mat_vl,
    name = "Pearson r",
    col = circlize::colorRamp2(c(-0.3, 0, 0.3, 0.6, 1), c("#2166AC", "#F7F7F7", "#F7F7F7", "#FDDBC7", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_title = "L1HS ORF1-intact (>5kb): Intra-read m6A correlation (all)",
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7)
)
pdf(str_glue("{outdir_vlong}/intraread_correlation_heatmap_all1.pdf"), width = 6, height = 5)
draw(p)
dev.off()

# Per-condition correlation heatmaps
for (cond in c("ORF1", "TOT")) {
    cond_mat <- vlong_wide %>%
        filter(bcondition == cond) %>%
        dplyr::select(all_of(vlong_site_cols)) %>%
        as.matrix()
    if (nrow(cond_mat) < 10) next

    cor_cond <- cor(cond_mat, use = "pairwise.complete.obs")
    colnames(cor_cond) <- gsub("^site_", "", colnames(cor_cond))
    rownames(cor_cond) <- gsub("^site_", "", rownames(cor_cond))
    cor_cond <- cor_cond[pos_order_vl, pos_order_vl]

    p <- Heatmap(
        cor_cond,
        name = "Pearson r",
        col = circlize::colorRamp2(c(-0.3, 0, 0.3, 0.6, 1), c("#2166AC", "#F7F7F7", "#F7F7F7", "#FDDBC7", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_title = str_glue("L1HS ORF1-intact (>5kb): Intra-read m6A correlation ({cond})"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7)
    )
    pdf(str_glue("{outdir_vlong}/intraread_correlation_heatmap_{tolower(cond)}.pdf"), width = 6, height = 5)
    draw(p)
    dev.off()
}

# Difference heatmap: ORF1 - TOT
orf1_mat_vl <- vlong_wide %>%
    filter(bcondition == "ORF1") %>%
    dplyr::select(all_of(vlong_site_cols)) %>%
    as.matrix()
tot_mat_vl <- vlong_wide %>%
    filter(bcondition == "TOT") %>%
    dplyr::select(all_of(vlong_site_cols)) %>%
    as.matrix()
if (nrow(orf1_mat_vl) >= 10 && nrow(tot_mat_vl) >= 10) {
    cor_diff_vl <- cor(orf1_mat_vl, use = "pairwise.complete.obs") - cor(tot_mat_vl, use = "pairwise.complete.obs")
    colnames(cor_diff_vl) <- gsub("^site_", "", colnames(cor_diff_vl))
    rownames(cor_diff_vl) <- gsub("^site_", "", rownames(cor_diff_vl))
    cor_diff_vl <- cor_diff_vl[pos_order_vl, pos_order_vl]

    p <- Heatmap(
        cor_diff_vl,
        name = "Δr (ORF1-TOT)",
        col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "#F7F7F7", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_title = "L1HS ORF1-intact (>5kb): Correlation diff (ORF1 - TOT)",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (i != j) grid.text(sprintf("%.2f", cor_diff_vl[i, j]), x, y, gp = gpar(fontsize = 5))
        }
    )
    pdf(str_glue("{outdir_vlong}/intraread_correlation_diff_orf1_minus_tot.pdf"), width = 8, height = 7)
    draw(p)
    dev.off()
}

# Load consensus sequences
l1hs_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1HS_fl_consensus.fa")
l1pa2_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1PA2_fl_consensus.fa")
l1pa3_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1PA3_fl_consensus.fa")
l1pa4_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1PA4_fl_consensus.fa")
l1pa5_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1PA5_fl_consensus.fa")
l1pa6_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment_meth/alignments/L1PA6_fl_consensus.fa")
consensus_seqs <- list(L1HS = l1hs_cs, L1PA2 = l1pa2_cs, L1PA3 = l1pa3_cs, L1PA4 = l1pa4_cs, L1PA5 = l1pa5_cs, L1PA6 = l1pa6_cs)

# Load aligned consensus
aligned <- readDNAStringSet("lrna/Rintermediates/l1_subfamily_consensus_aligned.fa")

# Load mapping tables
consensus_index_long <- read_csv(sprintf("ldna/results/m/plots/l1_alignment/all_subfamilies_mapping_to_consensus_harmonized.csv"))

# Build pos_to_aln_df: consensus_pos -> harmonized_pos for each subfamily
subfam_pos_to_aln <- list()
for (sf in names(aligned)) {
    chars <- strsplit(as.character(aligned[[sf]]), "")[[1]]
    ungapped <- 0L
    mapping <- integer(0)
    for (col in seq_along(chars)) {
        if (chars[col] != "-") {
            ungapped <- ungapped + 1L
            mapping[ungapped] <- col
        }
    }
    subfam_pos_to_aln[[sf]] <- tibble(
        subfamily_consensus = sf,
        consensus_pos = seq_along(mapping),
        harmonized_pos = mapping
    )
}
pos_to_aln_df <- bind_rows(subfam_pos_to_aln)

# Load expression / enrichment data



# contrasts_lr <- c("condition_N_ORF1_vs_N_TOT")
# elements_for_enrichment <- rmann %>% filter(rte_subfamily %in% YOUNG_L1_SUBFAMILIES, rte_length_req == "FL") %$% gene_id
# elements_for_enrichment_all <- rmann %>% filter(rte_subfamily %in% YOUNG_L1_SUBFAMILIES) %$% gene_id
# enrichment_lr <- tibble(gene_id = elements_for_enrichment)
# for (contrast in contrasts_lr) {
#     filepath <- str_glue("lrna/results/agg/deseq/telescope_multi/{contrast}/results_rtes.csv")
#     df <- read_csv(filepath, show_col_types = FALSE)
#     colnames(df)[1] <- "gene_id"
#     df <- df %>%
#         filter(gene_id %in% elements_for_enrichment) %>%
#         dplyr::select(gene_id, baseMean, log2FoldChange, lfcSE, padj) %>%
#         dplyr::rename(
#             !!str_glue("baseMean_{contrast}") := baseMean,
#             !!str_glue("l2fc_{contrast}") := log2FoldChange,
#             !!str_glue("lfcSE_{contrast}") := lfcSE,
#             !!str_glue("padj_{contrast}") := padj
#         )
#     enrichment_lr <- enrichment_lr %>% left_join(df, by = "gene_id")
# }
# for (contrast in contrasts_lr) {
#     enrichment_lr <- enrichment_lr %>%
#         mutate(!!str_glue("l2fc_{contrast}") := replace_na(!!sym(str_glue("l2fc_{contrast}")), 0)) %>%
#         mutate(!!str_glue("padj_{contrast}") := replace_na(!!sym(str_glue("padj_{contrast}")), 1))
# }

# normed_lr_raw <- read_csv("lrna/results/agg/deseq/telescope_multi/counttablesizenormed.csv", show_col_types = FALSE)
# colnames(normed_lr_raw)[1] <- "gene_id"
# normed_lr_long_all <- normed_lr_raw %>%
#     filter(gene_id %in% elements_for_enrichment_all) %>%
#     pivot_longer(-gene_id, names_to = "sample_name", values_to = "normed_count") %>%
#     left_join(sample_table_lrna %>% dplyr::select(sample_name, condition), by = "sample_name")
# normed_lr_long <- normed_lr_long_all %>% filter(gene_id %in% elements_for_enrichment)
# normed_lr_mean <- normed_lr_long %>%
#     group_by(gene_id, condition) %>%
#     summarise(mean_normed = mean(normed_count, na.rm = TRUE), .groups = "drop") %>%
#     pivot_wider(names_from = condition, values_from = mean_normed, names_prefix = "mean_lr_")
# enrichment_lr <- enrichment_lr %>%
#     left_join(normed_lr_mean, by = "gene_id") %>%
#     mutate(across(starts_with("mean_lr_"), ~ replace_na(., 0)))

# # Build enrichment for all-length elements (FL + trnc)
# enrichment_lr_all <- tibble(gene_id = elements_for_enrichment_all)
# for (contrast in contrasts_lr) {
#     filepath <- str_glue("lrna/results/agg/deseq/telescope_multi/{contrast}/results_rtes.csv")
#     df <- read_csv(filepath, show_col_types = FALSE)
#     colnames(df)[1] <- "gene_id"
#     df <- df %>%
#         filter(gene_id %in% elements_for_enrichment_all) %>%
#         dplyr::select(gene_id, baseMean, log2FoldChange, lfcSE, padj) %>%
#         dplyr::rename(
#             !!str_glue("baseMean_{contrast}") := baseMean,
#             !!str_glue("l2fc_{contrast}") := log2FoldChange,
#             !!str_glue("lfcSE_{contrast}") := lfcSE,
#             !!str_glue("padj_{contrast}") := padj
#         )
#     enrichment_lr_all <- enrichment_lr_all %>% left_join(df, by = "gene_id")
# }
# for (contrast in contrasts_lr) {
#     enrichment_lr_all <- enrichment_lr_all %>%
#         mutate(!!str_glue("l2fc_{contrast}") := replace_na(!!sym(str_glue("l2fc_{contrast}")), 0)) %>%
#         mutate(!!str_glue("padj_{contrast}") := replace_na(!!sym(str_glue("padj_{contrast}")), 1))
# }
# normed_lr_mean_all <- normed_lr_long_all %>%
#     group_by(gene_id, condition) %>%
#     summarise(mean_normed = mean(normed_count, na.rm = TRUE), .groups = "drop") %>%
#     pivot_wider(names_from = condition, values_from = mean_normed, names_prefix = "mean_lr_")
# enrichment_lr_all <- enrichment_lr_all %>%
#     left_join(normed_lr_mean_all, by = "gene_id") %>%
#     mutate(across(starts_with("mean_lr_"), ~ replace_na(., 0)))

count_data <- readRDS("custom/counts/count_data.rds")
method_data <- count_data$method_data
quant_all <- count_data$quant_all
quant_all_raw <- count_data$quant_all_raw
annot <- count_data$annot

annot %$% orf1_intactness %>% table()
# Create output directories
outdir_base <- str_glue("lrna/results/plots/m6a_l1_analysis")
outdir_ivt <- str_glue("{outdir_base}/ivt_filtering")
outdir_sites <- str_glue("{outdir_base}/reliable_sites")
outdir_glmm <- str_glue("{outdir_base}/glmm")
outdir_evo <- str_glue("{outdir_base}/evolution")
intdir <- "lrna/Rintermediates/m6a_l1_analysis"
for (d in c(outdir_ivt, outdir_sites, outdir_glmm, outdir_evo, intdir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# L1 gene structure coordinates (approximate, for annotation on plots)
# L1HS consensus is ~6kb: 5'UTR (1-910), ORF1 (911-2811), inter-ORF (2812-2981), ORF2 (2982-5817), 3'UTR (5818-6064)
L1_REGIONS <- tibble(
    region = c("5'UTR", "ORF1", "Inter-ORF", "ORF2", "3'UTR"),
    start = c(1, 907, 1924, 1987, 5815),
    end = c(906, 1923, 1986, 5814, 6032),
    color = c("#E69F00", "#56B4E9", "grey70", "#009E73", "#F0E442")
)

# Helper: add L1 gene structure annotation to a ggplot with x = position
add_l1_structure <- function(p, y_pos = -5, height = 3, alpha = 0.3) {
    p + annotate("rect",
        xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
        ymin = y_pos - height / 2, ymax = y_pos + height / 2,
        fill = L1_REGIONS$color, alpha = alpha
    ) +
        annotate("text",
            x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
            y = y_pos, label = L1_REGIONS$region, size = 2
        )
}


# ============================================================
# //ANCHOR - SECTION 1: IVT FALSE POSITIVE FILTERING
# ============================================================
cat("\n=== Section 1: IVT False Positive Filtering ===\n")

# Parameters
FP_THRESHOLD <- 5 # % modification to call a position a false positive
FP_MIN_COV <- 1000 # minimum coverage in IVT to trust the call

# 1a. Load IVT bedMethyl
ivt_bed_path <- "/users/mkelsey/data/n2102ep/RTE/lrna/intermediates/N_ORF1_4/methylation/bymod/N_ORF1_4_construct_m6A_bedMethyl.bed"
# bedMethyl column 10 is space-separated (Nvalid pctMod Nmod ...); read_table splits on all whitespace
ivt_raw <- read_table(ivt_bed_path, col_names = FALSE, show_col_types = FALSE)
# After whitespace splitting: X1=chrom, X2=start, X3=end, X4=name, X5=score, X6=strand,
# X7=thickStart, X8=thickEnd, X9=itemRgb, X10=cov(Nvalid), X11=pctM, ...
ivt_raw <- ivt_raw %>%
    dplyr::rename(
        chrom = X1, start = X2, end = X3, name = X4, score = X5, strand = X6,
        thickStart = X7, thickEnd = X8, itemRgb = X9, cov = X10, pctM = X11
    ) %>%
    mutate(pctM = as.double(pctM))
cat(sprintf("IVT bedMethyl: %d positions on %s\n", nrow(ivt_raw), paste(unique(ivt_raw$chrom), collapse = ", ")))

# 1b. Load construct FASTA and align to L1HS consensus
construct_fa <- readDNAStringSet("/users/mkelsey/data/n2102ep/RTE/conf/construct.fa")
construct_seq <- construct_fa[[1]]
# cat(sprintf("Construct sequence length: %d bp\n", BiocGenerics::width(construct_seq)))

# Pairwise align construct to L1HS consensus to find L1 insert and build position mapping
aln_result <- pairwiseAlignment(construct_seq, l1hs_cs[[1]], type = "local")
construct_aln_start <- start(Biostrings::pattern(aln_result))
construct_aln_end <- end(Biostrings::pattern(aln_result))
consensus_aln_start <- start(Biostrings::subject(aln_result))
cat(sprintf(
    "L1 insert in construct: positions %d-%d (aligned to L1HS consensus starting at %d)\n",
    construct_aln_start, construct_aln_end, consensus_aln_start
))

# Build position mapping: construct position -> L1HS consensus position
# Walk through the aligned sequences, tracking ungapped positions on both sides
pattern_chars <- strsplit(as.character(Biostrings::pattern(aln_result)), "")[[1]]
subject_chars <- strsplit(as.character(Biostrings::subject(aln_result)), "")[[1]]

construct_pos_vec <- integer(0)
consensus_pos_vec <- integer(0)
c_pos <- construct_aln_start - 1L
s_pos <- consensus_aln_start - 1L

for (i in seq_along(pattern_chars)) {
    p_gap <- pattern_chars[i] == "-"
    s_gap <- subject_chars[i] == "-"
    if (!p_gap) c_pos <- c_pos + 1L
    if (!s_gap) s_pos <- s_pos + 1L
    if (!p_gap && !s_gap) {
        construct_pos_vec <- c(construct_pos_vec, c_pos)
        consensus_pos_vec <- c(consensus_pos_vec, s_pos)
    }
}
ivt_to_consensus <- tibble(construct_pos = construct_pos_vec, consensus_pos = consensus_pos_vec)
cat(sprintf("Mapped %d construct positions to L1HS consensus\n", nrow(ivt_to_consensus)))

# 1c. Map IVT modification to consensus space
# IVT bedmethyl uses 0-based coordinates; construct_pos is 1-based
ivt_profile <- ivt_raw %>%
    mutate(construct_pos = start + 1L) %>%
    left_join(ivt_to_consensus, by = "construct_pos") %>%
    filter(!is.na(consensus_pos))

cat(sprintf("IVT positions mapped to L1HS consensus: %d\n", nrow(ivt_profile)))

# 1d. Identify false positive positions
ivt_fp <- ivt_profile %>%
    filter(cov >= FP_MIN_COV, pctM > FP_THRESHOLD)
fp_consensus_positions <- unique(ivt_fp$consensus_pos)
cat(sprintf("False positive positions (pctM > %g%%, cov >= %d): %d\n", FP_THRESHOLD, FP_MIN_COV, length(fp_consensus_positions)))

# 1e. Map FP positions to harmonized space
l1hs_mapping <- pos_to_aln_df %>% filter(subfamily_consensus == "L1HS")
fp_harmonized_positions <- l1hs_mapping %>%
    filter(consensus_pos %in% fp_consensus_positions) %>%
    pull(harmonized_pos) %>%
    unique()
cat(sprintf("FP harmonized positions: %d\n", length(fp_harmonized_positions)))

# Save FP positions
write_tsv(tibble(
    consensus_pos = fp_consensus_positions,
    harmonized_pos = l1hs_mapping$harmonized_pos[match(fp_consensus_positions, l1hs_mapping$consensus_pos)]
), str_glue("{intdir}/fp_positions.tsv"))

# 1f. Filter yl1annm (FL only) and yl1annm_alllength (FL + trnc)
n_before <- nrow(yl1annm)
yl1annm_filtered <- yl1annm %>% filter(!(harmonized_pos %in% fp_harmonized_positions))
n_after <- nrow(yl1annm_filtered)
cat(sprintf(
    "yl1annm (FL): %d rows before -> %d rows after FP filtering (%d removed, %.1f%%)\n",
    n_before, n_after, n_before - n_after, 100 * (n_before - n_after) / n_before
))

n_before_all <- nrow(yl1annm_alllength)
yl1annm_alllength_filtered <- yl1annm_alllength %>% filter(!(harmonized_pos %in% fp_harmonized_positions))
n_after_all <- nrow(yl1annm_alllength_filtered)
cat(sprintf(
    "yl1annm_alllength (FL+trnc): %d rows before -> %d rows after FP filtering (%d removed, %.1f%%)\n",
    n_before_all, n_after_all, n_before_all - n_after_all, 100 * (n_before_all - n_after_all) / n_before_all
))

# Update position sets
posmeth_filtered <- setdiff(
    yl1annm_filtered %>%
        group_by(harmonized_pos, rte_subfamily) %>%
        summarise(sc = sum(cov), n = n(), mm = mean(pctM), .groups = "drop") %>%
        filter(mm >= 10) %$% harmonized_pos %>% unique(),
    fp_harmonized_positions
)
poskeep_filtered <- setdiff(
    yl1annm_filtered %>%
        group_by(harmonized_pos, rte_subfamily) %>%
        summarise(sc = sum(cov), n = n(), .groups = "drop") %>%
        filter(n >= quantile(n, 0.25), sc >= quantile(sc, 0.25)) %$% harmonized_pos %>% unique(),
    fp_harmonized_positions
)

# Save IVT profile
write_tsv(ivt_profile, str_glue("{intdir}/ivt_profile.tsv"))

# --- IVT PLOTS ---

# Plot 1: IVT modification profile along element
p <- ivt_profile %>%
    filter(cov >= FP_MIN_COV) %>%
    mutate(is_fp = consensus_pos %in% fp_consensus_positions) %>%
    ggplot(aes(x = consensus_pos, y = pctM)) +
    geom_point(aes(color = is_fp), size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = FP_THRESHOLD, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red"), labels = c("Pass", "False Positive")) +
    annotate("rect",
        xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
        ymin = -3, ymax = -1, fill = L1_REGIONS$color, alpha = 0.5
    ) +
    annotate("text",
        x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
        y = -2, label = L1_REGIONS$region, size = 2.5
    ) +
    labs(
        x = "L1HS consensus position", y = "% m6A (IVT)", color = NULL
    ) +
    coord_cartesian(ylim = c(-4, max(ivt_profile$pctM[ivt_profile$cov >= 5], na.rm = TRUE) * 1.1)) +
    theme_minimal(base_size = 12)
mss(str_glue("{outdir_ivt}/ivt_modification_profile11.pdf"), fw = 4, fh = 2, plus_void = TRUE)

# Plot 2: IVT coverage along element
p <- ivt_profile %>%
    ggplot(aes(x = consensus_pos, y = cov)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, width = 1) +
    geom_hline(yintercept = FP_MIN_COV, linetype = "dashed", color = "red") +
    labs(x = "L1HS consensus position", y = "Coverage", title = "IVT coverage along L1HS construct") +
    theme_minimal(base_size = 12)
mss(str_glue("{outdir_ivt}/ivt_coverage_profile.pdf"), fw = 4, fh = 2)

# Plot 3: Histogram of IVT pctM
p <- ivt_profile %>%
    filter(cov >= FP_MIN_COV) %>%
    ggplot(aes(x = pctM)) +
    geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white") +
    geom_vline(xintercept = FP_THRESHOLD, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = FP_THRESHOLD + 1, y = Inf, vjust = 2, label = str_glue("FP threshold = {FP_THRESHOLD}%"), color = "red", size = 3.5) +
    labs(x = "% m6A in IVT", y = "# positions", title = "Distribution of apparent m6A in IVT (no true modification expected)") +
    theme_minimal(base_size = 12)
mss(str_glue("{outdir_ivt}/ivt_pctM_histogram111.pdf"), fw = 1.5, fh = 1.5, plus_void = TRUE)

# Plot 4: Before/after filtering comparison
orf1_mean_before <- yl1annm %>%
    filter(condition == "N_ORF1", rte_subfamily == "L1HS") %>%
    group_by(harmonized_pos) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    mutate(filter_status = "Before IVT filter")

orf1_mean_after <- yl1annm_filtered %>%
    filter(condition == "N_ORF1", rte_subfamily == "L1HS") %>%
    group_by(harmonized_pos) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    mutate(filter_status = "After IVT filter")

# Map harmonized_pos back to approximate consensus_pos for x-axis
hp_to_cp <- l1hs_mapping %>% dplyr::select(harmonized_pos, consensus_pos)

before_after <- bind_rows(orf1_mean_before, orf1_mean_after) %>%
    left_join(hp_to_cp, by = "harmonized_pos") %>%
    mutate(is_removed = harmonized_pos %in% fp_harmonized_positions)

p <- before_after %>%
    ggplot(aes(x = consensus_pos, y = mm)) +
    geom_point(aes(color = is_removed), size = 1, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey30", "TRUE" = "red"), labels = c("Retained", "Removed by IVT filter")) +
    facet_wrap(~filter_status, ncol = 1) +
    labs(
        x = "L1HS consensus position", y = "Mean % m6A (ORF1 RIP)", color = NULL,
        title = "ORF1 RIP L1HS modification: before vs after IVT filtering"
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_ivt}/before_after_filtering.pdf"), w = 12, h = 8)

# Plot 5: IVT pctM vs ORF1-RIP pctM scatter
orf1_at_ivt_positions <- yl1annm %>%
    filter(condition == "N_ORF1", rte_subfamily == "L1HS") %>%
    group_by(harmonized_pos) %>%
    summarise(orf1_mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    left_join(hp_to_cp, by = "harmonized_pos")

ivt_for_scatter <- ivt_profile %>%
    filter(cov >= FP_MIN_COV) %>%
    dplyr::select(consensus_pos, ivt_pctM = pctM)

scatter_df <- orf1_at_ivt_positions %>%
    inner_join(ivt_for_scatter, by = "consensus_pos")

p <- scatter_df %>%
    mutate(is_fp = consensus_pos %in% fp_consensus_positions) %>%
    ggplot(aes(x = ivt_pctM, y = orf1_mm)) +
    geom_point(aes(color = is_fp), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
    scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red"), labels = c("Retained", "False Positive")) +
    geom_text_repel(
        data = . %>% filter(is_fp | orf1_mm > 30),
        aes(label = consensus_pos), size = 2.5, max.overlaps = 20
    ) +
    labs(
        x = "% m6A in IVT (no real modification)", y = "% m6A in ORF1 RIP",
        color = NULL
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_ivt}/ivt_vs_orf1_scatter1.pdf"), w = 6, h = 5)

# Plot 6: Summary stats
cat(sprintf("\n--- IVT Filtering Summary ---\n"))
cat(sprintf("Total IVT positions mapped to L1HS: %d\n", nrow(ivt_profile)))
cat(sprintf("Positions with adequate coverage (>=%d): %d\n", FP_MIN_COV, sum(ivt_profile$cov >= FP_MIN_COV)))
cat(sprintf("False positive positions (>%g%% mod): %d\n", FP_THRESHOLD, length(fp_consensus_positions)))
cat(sprintf("FP harmonized positions removed: %d\n", length(fp_harmonized_positions)))
cat(sprintf(
    "yl1annm rows removed: %d / %d (%.1f%%)\n",
    n_before - n_after, n_before, 100 * (n_before - n_after) / n_before
))


# ============================================================
# //ANCHOR - SECTIONS 2+3 WRAPPER FUNCTION
# ============================================================

# Define fit_site_model at global scope so future_map workers don't serialize
# the entire run_sections_2_and_3 environment
GLMM_TIMEOUT_SECONDS <- 120 # max seconds per model fit

fit_site_model <- function(target_site, data, site_cols, include_other_sites = TRUE) {
    model_data <- data %>%
        filter(!is.na(.data[[target_site]]))

    other_sites_keep <- character(0)
    if (include_other_sites) {
        other_sites <- setdiff(site_cols, target_site)
        other_sites_keep <- other_sites[colMeans(!is.na(model_data[other_sites])) > 0.5]
        model_data <- model_data %>% drop_na(all_of(other_sites_keep))
    }

    if (nrow(model_data) < 50) {
        return(NULL)
    }

    fixed_terms <- c(other_sites_keep, "l2fc_scaled", "baseMean_scaled", "intactness_req", "genic_loc")
    for (term in fixed_terms) {
        if (is.numeric(model_data[[term]])) {
            if (var(model_data[[term]], na.rm = TRUE) == 0) fixed_terms <- setdiff(fixed_terms, term)
        } else {
            if (n_distinct(model_data[[term]], na.rm = TRUE) < 2) fixed_terms <- setdiff(fixed_terms, term)
        }
    }

    formula_str <- paste0(
        target_site, " ~ ",
        paste(fixed_terms, collapse = " + "),
        " + (1|gene_id) + (1|sample)"
    )

    # Wrap in timeout to avoid hanging on convergence issues
    fit_with_timeout <- function(formula_str, model_data) {
        setTimeLimit(elapsed = GLMM_TIMEOUT_SECONDS, transient = TRUE)
        on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
        glmmTMB(as.formula(formula_str),
            data = model_data,
            family = binomial(link = "logit")
        )
    }

    model <- tryCatch(
        fit_with_timeout(formula_str, model_data),
        error = function(e) {
            cat(sprintf("    Primary model failed for %s: %s\n", target_site, conditionMessage(e)))
            formula_str2 <- paste0(
                target_site, " ~ ",
                paste(fixed_terms, collapse = " + "),
                " + sample + (1|gene_id)"
            )
            tryCatch(
                fit_with_timeout(formula_str2, model_data),
                error = function(e2) {
                    cat(sprintf("    Fallback model also failed for %s: %s\n", target_site, conditionMessage(e2)))
                    NULL
                }
            )
        }
    )
    return(model)
}

# Wrap Sections 2 and 3 in a function to run for different subfamily sets and length filters
run_sections_2_and_3 <- function(subfamilies, label, suffix, length_filter = "FL") {
    # label: for plot titles (e.g. "Young L1" or "L1HS")
    # suffix: for file names (e.g. "" or "_L1HS" or "_alllength")
    # length_filter: "FL" for full-length only, "all" for FL + truncated
    # mysaveandstore's default pl=p resolves p in the global env (where it was defined),
    # not our local scope. Shadow it with a version that captures p from this scope.
    .global_mysaveandstore <- get("mysaveandstore", envir = globalenv())
    mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, store = store_var, raster = FALSE, ...) {
        .global_mysaveandstore(fn = fn, w = w, h = h, res = res, pl = p, store = store, raster = raster, ...)
    }

    # Select the appropriate data source based on length_filter
    if (length_filter == "all") {
        meth_data <- yl1annm_alllength_filtered
        enrich_data <- enrichment_lr_all
        length_label <- " (FL+trnc)"
    } else {
        meth_data <- yl1annm_filtered
        enrich_data <- enrichment_lr
        length_label <- " (FL)"
    }
    full_label <- paste0(label, length_label)

    # Shadow global output dirs so all str_glue references within this function write to suffixed paths
    outdir_sites <- str_glue("{outdir_sites}{suffix}")
    outdir_glmm <- str_glue("{outdir_glmm}{suffix}")
    intdir <- str_glue("{intdir}{suffix}")
    for (d in c(outdir_sites, outdir_glmm, intdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    # ============================================================
    # //ANCHOR - SECTION 2: IDENTIFY RELIABLY MODIFIED SITES
    # ============================================================
    cat(sprintf("\n=== Section 2: Identify Reliably Modified Sites [%s] ===\n", full_label))

    # Parameters
    MIN_MEAN_PCTM <- 20
    MIN_REP_PCTM <- 2
    MIN_FRAC_ELEMENTS <- 0.1
    MIN_N_ELEMENTS <- 10
    MIN_TOTAL_COV <- 50

    # 2a. Subset to requested subfamilies and length filter (pool all conditions)
    filtered_meth <- meth_data %>%
        filter(rte_subfamily %in% subfamilies)
    if (length_filter == "FL") {
        filtered_meth <- filtered_meth %>% filter(rte_length_req == "FL")
    }
    conditions_present <- unique(filtered_meth$condition)
    cat(sprintf("Conditions pooled: %s\n", paste(conditions_present, collapse = ", ")))

    # 2b. Per-position, per-sample summary (pooling all conditions for max power)
    pos_sample_summary <- filtered_meth %>%
        filter(cov >= 5) %>%
        group_by(harmonized_pos, sample) %>%
        summarise(
            condition = first(condition),
            n_elements = n_distinct(gene_id),
            mean_pctM = mean(pctM, na.rm = TRUE),
            median_pctM = median(pctM, na.rm = TRUE),
            frac_mod = mean(pctM > 10, na.rm = TRUE),
            total_cov = sum(cov, na.rm = TRUE),
            weighted_pctM = sum(pctM * cov) / sum(cov),
            .groups = "drop"
        )
    write_tsv(pos_sample_summary, str_glue("{intdir}/pos_sample_summary.tsv"))

    # 2c. Cross-sample summary (all samples pooled regardless of condition)
    pos_summary <- pos_sample_summary %>%
        group_by(harmonized_pos) %>%
        summarise(
            n_samples = n(),
            mean_pctM = mean(weighted_pctM, na.rm = TRUE),
            sd_pctM = sd(weighted_pctM, na.rm = TRUE),
            cv_pctM = ifelse(mean(weighted_pctM) > 0, sd(weighted_pctM) / mean(weighted_pctM), NA_real_),
            min_sample_pctM = min(weighted_pctM, na.rm = TRUE),
            max_sample_pctM = max(weighted_pctM, na.rm = TRUE),
            mean_frac_mod = mean(frac_mod, na.rm = TRUE),
            mean_n_elements = mean(n_elements, na.rm = TRUE),
            mean_total_cov = mean(total_cov, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        left_join(hp_to_cp, by = "harmonized_pos")

    # 2d. Multi-criteria filtering
    reliable_sites <- pos_summary %>%
        filter(
            mean_pctM >= MIN_MEAN_PCTM,
            min_sample_pctM >= MIN_REP_PCTM,
            mean_frac_mod >= MIN_FRAC_ELEMENTS,
            mean_n_elements >= MIN_N_ELEMENTS,
            mean_total_cov >= MIN_TOTAL_COV
        )
    reliable_sites_harmonized_pos <- reliable_sites %$% harmonized_pos

    reliable_sites2 <- pos_summary %>%
        filter(
            mean_pctM <= MIN_MEAN_PCTM,
            mean_pctM >= 10,
            mean_n_elements >= MIN_N_ELEMENTS,
            mean_total_cov >= 100
        )

    cat(sprintf("Total positions assessed: %d\n", nrow(pos_summary)))
    cat(sprintf("Reliable sites (pooled): %d\n", nrow(reliable_sites)))

    # 2e. Classify sites as reliable or not
    site_classification <- pos_summary %>%
        mutate(
            reliable = mean_pctM >= MIN_MEAN_PCTM & min_sample_pctM >= MIN_REP_PCTM &
                mean_frac_mod >= MIN_FRAC_ELEMENTS & mean_n_elements >= MIN_N_ELEMENTS &
                mean_total_cov >= MIN_TOTAL_COV
        )

    write_tsv(pos_summary, str_glue("{intdir}/pos_summary.tsv"))
    write_tsv(reliable_sites, str_glue("{intdir}/reliable_sites.tsv"))
    write_tsv(site_classification, str_glue("{intdir}/site_classification.tsv"))

    cat(sprintf(
        "Reliable: %d, Not reliable: %d\n",
        sum(site_classification$reliable), sum(!site_classification$reliable)
    ))

    # --- RELIABLE SITES PLOTS ---

    # Plot 1: Modification profile with reliable sites highlighted
    p <- site_classification %>%
        filter(!is.na(consensus_pos)) %>%
        ggplot(aes(x = consensus_pos, y = mean_pctM)) +
        geom_point(aes(color = reliable, size = mean_n_elements), alpha = 0.7) +
        scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey70")) +
        scale_size_continuous(range = c(0.5, 4)) +
        geom_text_repel(
            data = . %>% filter(reliable),
            aes(label = harmonized_pos), size = 2.5, max.overlaps = 25
        ) +
        annotate("rect",
            xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
            ymin = -4, ymax = -2, fill = L1_REGIONS$color, alpha = 0.5
        ) +
        annotate("text",
            x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
            y = -3, label = L1_REGIONS$region, size = 2.5
        ) +
        labs(
            x = "L1HS consensus position", y = "Mean % m6A (all conditions pooled)",
            color = "Reliable", size = "# elements",
            title = sprintf("m6A modification profile: reliable sites [%s]", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_sites}/modification_profile_reliable_sites.pdf"), w = 14, h = 6)

    # Plot 2: Per-condition profiles at reliable sites
    pos_by_condition <- filtered_meth %>%
        filter(cov >= 5) %>%
        group_by(condition, harmonized_pos) %>%
        summarise(
            mean_pctM = weighted.mean(pctM, cov, na.rm = TRUE),
            n_elements = n_distinct(gene_id),
            .groups = "drop"
        ) %>%
        left_join(hp_to_cp, by = "harmonized_pos")

    p <- pos_by_condition %>%
        filter(!is.na(consensus_pos)) %>%
        ggplot(aes(x = consensus_pos, y = mean_pctM, color = condition)) +
        geom_point(alpha = 0.4, size = 0.8) +
        geom_smooth(method = "loess", se = FALSE, span = 0.1) +
        geom_point(
            data = . %>% filter(harmonized_pos %in% reliable_sites$harmonized_pos),
            size = 2, shape = 1
        ) +
        annotate("rect",
            xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
            ymin = -4, ymax = -2, fill = L1_REGIONS$color, alpha = 0.5
        ) +
        annotate("text",
            x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
            y = -3, label = L1_REGIONS$region, size = 2.5
        ) +
        labs(
            x = "L1HS consensus position", y = "Mean % m6A",
            color = "Condition",
            title = sprintf("m6A profile by condition [%s]", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_sites}/modification_profile_by_condition.pdf"), w = 14, h = 6)

    # Plot 3: Replicate reproducibility scatter (all samples)
    all_samples <- unique(pos_sample_summary$sample)
    rep_pairs <- combn(all_samples, 2, simplify = FALSE)
    rep_scatter_list <- list()
    for (pair in rep_pairs) {
        pair_df <- pos_sample_summary %>%
            filter(sample %in% pair) %>%
            dplyr::select(harmonized_pos, sample, weighted_pctM) %>%
            pivot_wider(names_from = sample, values_from = weighted_pctM)
        colnames(pair_df)[2:3] <- c("rep1", "rep2")
        r_val <- cor(pair_df$rep1, pair_df$rep2, use = "pairwise.complete.obs")
        pair_conds <- pos_sample_summary %>%
            filter(sample %in% pair) %>%
            distinct(sample, condition) %>%
            pull(condition)
        pair_label <- sprintf("%s (%s) vs %s (%s)", pair[1], pair_conds[1], pair[2], pair_conds[2])
        pair_df <- pair_df %>%
            left_join(site_classification %>% dplyr::select(harmonized_pos, reliable), by = "harmonized_pos")
        rep_scatter_list[[pair_label]] <- pair_df %>%
            ggplot(aes(x = rep1, y = rep2, color = reliable)) +
            geom_point(alpha = 0.5, size = 1) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
            scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey70")) +
            annotate("text", x = 5, y = 90, label = sprintf("r = %.3f", r_val), size = 4) +
            labs(x = pair[1], y = pair[2], title = pair_label) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none")
    }
    p <- plot_grid(plotlist = rep_scatter_list, ncol = min(length(rep_scatter_list), 3))
    mysaveandstore(str_glue("{outdir_sites}/replicate_reproducibility_scatter.pdf"),
        w = 6 * min(length(rep_scatter_list), 3), h = 5 * ceiling(length(rep_scatter_list) / 3)
    )

    # Plot 4: Volcano-style: mean modification vs consistency (inverted CV)
    p <- site_classification %>%
        filter(!is.na(cv_pctM), cv_pctM > 0) %>%
        ggplot(aes(x = mean_pctM, y = -log10(cv_pctM + 0.01), color = reliable)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey70")) +
        geom_text_repel(
            data = . %>% filter(reliable),
            aes(label = harmonized_pos), size = 2.5, max.overlaps = 20
        ) +
        labs(
            x = "Mean % m6A (pooled)", y = "-log10(CV)", color = "Reliable",
            title = "Site reliability: modification level vs cross-sample consistency"
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_sites}/volcano_style_reliability.pdf"), w = 8, h = 6)

    # Plot 5: Element-level heatmap at reliable sites
    if (nrow(reliable_sites) > 0) {
        element_hm <- filtered_meth %>%
            filter(harmonized_pos %in% reliable_sites$harmonized_pos, cov >= 5) %>%
            group_by(gene_id, harmonized_pos, rte_subfamily) %>%
            summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            group_by(gene_id) %>%
            filter(n() >= 0.5 * nrow(reliable_sites)) %>%
            ungroup()

        if (nrow(element_hm) > 0 && n_distinct(element_hm$gene_id) >= 5) {
            hm_mat <- element_hm %>%
                pivot_wider(id_cols = gene_id, names_from = harmonized_pos, values_from = mean_pctM, values_fn = mean) %>%
                column_to_rownames("gene_id") %>%
                as.matrix()
            hm_mat <- hm_mat[, order(as.numeric(colnames(hm_mat))), drop = FALSE]
            hm_mat <- hm_mat[rowSums(!is.na(hm_mat)) > 0, , drop = FALSE]

            row_ann_df <- rmann %>%
                filter(gene_id %in% rownames(hm_mat)) %>%
                dplyr::select(gene_id, intactness_req, genic_loc, rte_subfamily, rte_length_req) %>%
                distinct() %>%
                column_to_rownames("gene_id")
            row_ann_df <- row_ann_df[rownames(hm_mat), , drop = FALSE]
            row_ha <- rowAnnotation(
                Length = row_ann_df$rte_length_req,
                Intactness = row_ann_df$intactness_req,
                Location = row_ann_df$genic_loc,
                Subfamily = row_ann_df$rte_subfamily,
                col = list(Length = c("FL" = "#D62728", "Trnc" = "#1F77B4"))
            )

            # Build column annotation: map L1_REGIONS (consensus_pos) to harmonized_pos
            make_col_region_annotation <- function(hm_cols) {
                hp_vals <- as.numeric(hm_cols)
                col_regions <- sapply(hp_vals, function(hp) {
                    cp <- hp_to_cp$consensus_pos[match(hp, hp_to_cp$harmonized_pos)]
                    if (is.na(cp)) {
                        return("Unknown")
                    }
                    region_hit <- L1_REGIONS %>% filter(start <= cp, end >= cp)
                    if (nrow(region_hit) == 0) {
                        return("Unknown")
                    }
                    region_hit$region[1]
                })
                region_colors <- setNames(L1_REGIONS$color, L1_REGIONS$region)
                region_colors["Unknown"] <- "white"
                HeatmapAnnotation(
                    Region = col_regions,
                    col = list(Region = region_colors),
                    annotation_name_side = "left"
                )
            }
            col_ha <- make_col_region_annotation(colnames(hm_mat))

            p <- Heatmap(hm_mat,
                name = "% m6A",
                col = circlize::colorRamp2(c(0, 50, 100), c("blue", "white", "red")),
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = FALSE, na_col = "grey90",
                use_raster = TRUE,
                left_annotation = row_ha,
                top_annotation = col_ha,
                column_title = sprintf("Reliable m6A sites on %s elements (pooled)", full_label)
            )
            mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites.pdf"),
                w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
            )

            p <- Heatmap(hm_mat,
                name = "% m6A",
                col = circlize::colorRamp2(
                    seq(0, 100, length.out = 9),
                    viridisLite::viridis(9)
                ),
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = FALSE, na_col = "grey90",
                use_raster = TRUE,
                left_annotation = row_ha,
                top_annotation = col_ha,
                column_title = sprintf("Reliable m6A sites on %s elements (pooled) [viridis]", full_label)
            )
            mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_viridis.pdf"),
                w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
            )

            p <- Heatmap(hm_mat,
                name = "% m6A",
                col = circlize::colorRamp2(
                    seq(0, 100, length.out = 9),
                    viridisLite::plasma(9)
                ),
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = FALSE, na_col = "grey90",
                use_raster = TRUE,
                left_annotation = row_ha,
                top_annotation = col_ha,
                column_title = sprintf("Reliable m6A sites on %s elements (pooled) [plasma]", full_label)
            )
            mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_plasma.pdf"),
                w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
            )

            # Column-scaled (Z-score) version
            col_scale_mat <- function(mat) {
                scaled <- apply(mat, 2, function(x) {
                    s <- sd(x, na.rm = TRUE)
                    if (is.na(s) || s == 0) {
                        return(rep(0, length(x)))
                    }
                    (x - mean(x, na.rm = TRUE)) / s
                })
                rownames(scaled) <- rownames(mat)
                scaled
            }

            hm_mat_z <- col_scale_mat(hm_mat)
            z_max <- min(max(abs(hm_mat_z), na.rm = TRUE), 3)

            p <- Heatmap(hm_mat_z,
                name = "Z-score",
                col = circlize::colorRamp2(c(-z_max, 0, z_max), c("blue", "white", "red")),
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = FALSE, na_col = "grey90",
                use_raster = TRUE,
                left_annotation = row_ha,
                top_annotation = col_ha,
                column_title = sprintf("Reliable m6A sites (column Z-score) [%s]", full_label)
            )
            mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_zscore.pdf"),
                w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
            )

            p <- Heatmap(hm_mat_z,
                name = "Z-score",
                col = circlize::colorRamp2(
                    seq(-z_max, z_max, length.out = 9),
                    viridisLite::viridis(9)
                ),
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = FALSE, na_col = "grey90",
                use_raster = TRUE,
                left_annotation = row_ha,
                top_annotation = col_ha,
                column_title = sprintf("Reliable m6A sites (column Z-score) [%s, viridis]", full_label)
            )
            mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_zscore_viridis.pdf"),
                w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
            )

            # Heatmap with strict 90% site coverage filter
            element_hm_strict <- filtered_meth %>%
                filter(harmonized_pos %in% reliable_sites$harmonized_pos, cov >= 5) %>%
                group_by(gene_id, harmonized_pos, rte_subfamily) %>%
                summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
                group_by(gene_id) %>%
                filter(n() >= 0.9 * nrow(reliable_sites)) %>%
                ungroup()

            if (nrow(element_hm_strict) > 0 && n_distinct(element_hm_strict$gene_id) >= 5) {
                hm_mat_strict <- element_hm_strict %>%
                    pivot_wider(id_cols = gene_id, names_from = harmonized_pos, values_from = mean_pctM, values_fn = mean) %>%
                    column_to_rownames("gene_id") %>%
                    as.matrix()
                hm_mat_strict <- hm_mat_strict[, order(as.numeric(colnames(hm_mat_strict))), drop = FALSE]
                hm_mat_strict <- hm_mat_strict[rowSums(!is.na(hm_mat_strict)) > 0, , drop = FALSE]

                row_ann_df_strict <- rmann %>%
                    filter(gene_id %in% rownames(hm_mat_strict)) %>%
                    dplyr::select(gene_id, intactness_req, genic_loc, rte_subfamily, rte_length_req) %>%
                    distinct() %>%
                    column_to_rownames("gene_id")
                row_ann_df_strict <- row_ann_df_strict[rownames(hm_mat_strict), , drop = FALSE]
                row_ha_strict <- rowAnnotation(
                    Length = row_ann_df_strict$rte_length_req,
                    Intactness = row_ann_df_strict$intactness_req,
                    Location = row_ann_df_strict$genic_loc,
                    Subfamily = row_ann_df_strict$rte_subfamily,
                    col = list(Length = c("FL" = "#D62728", "Trnc" = "#1F77B4"))
                )

                col_ha_strict <- make_col_region_annotation(colnames(hm_mat_strict))

                p <- Heatmap(hm_mat_strict,
                    name = "% m6A",
                    col = circlize::colorRamp2(
                        seq(0, 100, length.out = 9),
                        viridisLite::viridis(9)
                    ),
                    cluster_rows = TRUE, cluster_columns = FALSE,
                    show_row_names = FALSE, na_col = "grey90",
                    use_raster = TRUE,
                    left_annotation = row_ha_strict,
                    top_annotation = col_ha_strict,
                    column_title = sprintf("Reliable m6A sites on %s elements (>=90%% site coverage) [viridis]", full_label)
                )
                mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_strict90_viridis.pdf"),
                    w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
                )
                cat(sprintf("Strict 90%% heatmap: %d elements\n", nrow(hm_mat_strict)))

                # Column-scaled (Z-score) version of strict 90%
                hm_mat_strict_z <- col_scale_mat(hm_mat_strict)
                z_max_strict <- min(max(abs(hm_mat_strict_z), na.rm = TRUE), 3)

                p <- Heatmap(hm_mat_strict_z,
                    name = "Z-score",
                    col = circlize::colorRamp2(c(-z_max_strict, 0, z_max_strict), c("blue", "white", "red")),
                    cluster_rows = TRUE, cluster_columns = FALSE,
                    show_row_names = FALSE, na_col = "grey90",
                    use_raster = TRUE,
                    left_annotation = row_ha_strict,
                    top_annotation = col_ha_strict,
                    column_title = sprintf("Reliable m6A sites (>=90%% coverage, column Z-score) [%s]", full_label)
                )
                mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_strict90_zscore.pdf"),
                    w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
                )

                p <- Heatmap(hm_mat_strict_z,
                    name = "Z-score",
                    col = circlize::colorRamp2(
                        seq(-z_max_strict, z_max_strict, length.out = 9),
                        viridisLite::viridis(9)
                    ),
                    cluster_rows = TRUE, cluster_columns = FALSE,
                    show_row_names = FALSE, na_col = "grey90",
                    use_raster = TRUE,
                    left_annotation = row_ha_strict,
                    top_annotation = col_ha_strict,
                    column_title = sprintf("Reliable m6A sites (>=90%% coverage, column Z-score) [%s, viridis]", full_label)
                )
                mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_strict90_zscore_viridis.pdf"),
                    w = max(8, nrow(reliable_sites) * 0.3 + 4), h = 10
                )
            }

            # Heatmap restricted to harmonized_pos < 1600
            reliable_sub <- reliable_sites %>% filter(harmonized_pos < 1600)
            if (nrow(reliable_sub) > 0) {
                element_hm_sub <- filtered_meth %>%
                    filter(harmonized_pos %in% reliable_sub$harmonized_pos, cov >= 5) %>%
                    group_by(gene_id, harmonized_pos, rte_subfamily) %>%
                    summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
                    group_by(gene_id) %>%
                    filter(n() >= 0.5 * nrow(reliable_sub)) %>%
                    ungroup()

                if (nrow(element_hm_sub) > 0 && n_distinct(element_hm_sub$gene_id) >= 5) {
                    hm_mat_sub <- element_hm_sub %>%
                        pivot_wider(id_cols = gene_id, names_from = harmonized_pos, values_from = mean_pctM, values_fn = mean) %>%
                        column_to_rownames("gene_id") %>%
                        as.matrix()
                    hm_mat_sub <- hm_mat_sub[, order(as.numeric(colnames(hm_mat_sub))), drop = FALSE]
                    hm_mat_sub <- hm_mat_sub[rowSums(!is.na(hm_mat_sub)) > 0, , drop = FALSE]

                    row_ann_df_sub <- rmann %>%
                        filter(gene_id %in% rownames(hm_mat_sub)) %>%
                        dplyr::select(gene_id, intactness_req, genic_loc, rte_subfamily, rte_length_req) %>%
                        distinct() %>%
                        column_to_rownames("gene_id")
                    row_ann_df_sub <- row_ann_df_sub[rownames(hm_mat_sub), , drop = FALSE]
                    row_ha_sub <- rowAnnotation(
                        Length = row_ann_df_sub$rte_length_req,
                        Intactness = row_ann_df_sub$intactness_req,
                        Location = row_ann_df_sub$genic_loc,
                        Subfamily = row_ann_df_sub$rte_subfamily,
                        col = list(Length = c("FL" = "#D62728", "Trnc" = "#1F77B4"))
                    )

                    col_ha_sub <- make_col_region_annotation(colnames(hm_mat_sub))

                    p <- Heatmap(hm_mat_sub,
                        name = "% m6A",
                        col = circlize::colorRamp2(
                            seq(0, 100, length.out = 9),
                            viridisLite::viridis(9)
                        ),
                        cluster_rows = TRUE, cluster_columns = FALSE,
                        show_row_names = FALSE, na_col = "grey90",
                        use_raster = TRUE,
                        left_annotation = row_ha_sub,
                        top_annotation = col_ha_sub,
                        column_title = sprintf("Reliable m6A sites < 1600 on %s elements [viridis]", full_label)
                    )
                    mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_lt1600_viridis.pdf"),
                        w = max(8, nrow(reliable_sub) * 0.3 + 4), h = 10
                    )

                    # Column-scaled (Z-score) version of lt1600
                    hm_mat_sub_z <- col_scale_mat(hm_mat_sub)
                    z_max_sub <- min(max(abs(hm_mat_sub_z), na.rm = TRUE), 3)

                    p <- Heatmap(hm_mat_sub_z,
                        name = "Z-score",
                        col = circlize::colorRamp2(c(-z_max_sub, 0, z_max_sub), c("blue", "white", "red")),
                        cluster_rows = TRUE, cluster_columns = FALSE,
                        show_row_names = FALSE, na_col = "grey90",
                        use_raster = TRUE,
                        left_annotation = row_ha_sub,
                        top_annotation = col_ha_sub,
                        column_title = sprintf("Reliable m6A sites < 1600 (column Z-score) [%s]", full_label)
                    )
                    mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_lt1600_zscore.pdf"),
                        w = max(8, nrow(reliable_sub) * 0.3 + 4), h = 10
                    )

                    p <- Heatmap(hm_mat_sub_z,
                        name = "Z-score",
                        col = circlize::colorRamp2(
                            seq(-z_max_sub, z_max_sub, length.out = 9),
                            viridisLite::viridis(9)
                        ),
                        cluster_rows = TRUE, cluster_columns = FALSE,
                        show_row_names = FALSE, na_col = "grey90",
                        use_raster = TRUE,
                        left_annotation = row_ha_sub,
                        top_annotation = col_ha_sub,
                        column_title = sprintf("Reliable m6A sites < 1600 (column Z-score) [%s, viridis]", full_label)
                    )
                    mysaveandstore(str_glue("{outdir_sites}/element_heatmap_reliable_sites_lt1600_zscore_viridis.pdf"),
                        w = max(8, nrow(reliable_sub) * 0.3 + 4), h = 10
                    )
                }
            }
        }
    }

    # Plot 6: Filtering funnel
    funnel_df <- tibble(
        step = c(
            "All positions", sprintf("Mean pctM >= %g%%", MIN_MEAN_PCTM), sprintf("Min sample pctM >= %g%%", MIN_REP_PCTM),
            sprintf("Frac elements >= %g%%", MIN_FRAC_ELEMENTS * 100), sprintf("N elements >= %d", MIN_N_ELEMENTS), sprintf("Total cov >= %d", MIN_TOTAL_COV)
        ),
        n = c(
            nrow(pos_summary),
            sum(pos_summary$mean_pctM >= MIN_MEAN_PCTM),
            sum(pos_summary$mean_pctM >= MIN_MEAN_PCTM & pos_summary$min_sample_pctM >= MIN_REP_PCTM),
            sum(pos_summary$mean_pctM >= MIN_MEAN_PCTM & pos_summary$min_sample_pctM >= MIN_REP_PCTM & pos_summary$mean_frac_mod >= MIN_FRAC_ELEMENTS),
            sum(pos_summary$mean_pctM >= MIN_MEAN_PCTM & pos_summary$min_sample_pctM >= MIN_REP_PCTM & pos_summary$mean_frac_mod >= MIN_FRAC_ELEMENTS & pos_summary$mean_n_elements >= MIN_N_ELEMENTS),
            nrow(reliable_sites)
        )
    ) %>% mutate(step = fct_inorder(step))

    p <- funnel_df %>%
        ggplot(aes(x = step, y = n)) +
        geom_col(fill = "steelblue", alpha = 0.8) +
        geom_text(aes(label = n), vjust = -0.5, size = 4) +
        labs(x = NULL, y = "# positions", title = sprintf("Filtering funnel for reliable m6A sites [%s]", full_label)) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    mysaveandstore(str_glue("{outdir_sites}/filtering_funnel.pdf"), w = 8, h = 5)

    # Plot 7: Density of pctM at reliable vs non-reliable sites
    p <- filtered_meth %>%
        filter(cov >= 5) %>%
        mutate(site_type = ifelse(harmonized_pos %in% reliable_sites$harmonized_pos, "Reliable", "Non-reliable")) %>%
        ggplot(aes(x = pctM, fill = site_type)) +
        geom_density(alpha = 0.4) +
        scale_fill_manual(values = c("Reliable" = "#D62728", "Non-reliable" = "grey70")) +
        facet_wrap(~condition) +
        labs(
            x = "% m6A per element-position", y = "Density", fill = NULL,
            title = sprintf("Modification distribution at reliable vs non-reliable sites [%s]", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_sites}/density_reliable_vs_nonreliable.pdf"),
        w = 5 * length(conditions_present), h = 5
    )

    # Plot 8: Per-site dot plot for reliable sites, colored by condition
    reliable_by_cond <- pos_by_condition %>%
        filter(harmonized_pos %in% reliable_sites$harmonized_pos)
    if (nrow(reliable_by_cond) > 0) {
        p <- reliable_by_cond %>%
            filter(!is.na(consensus_pos)) %>%
            mutate(harmonized_pos = fct_reorder(factor(harmonized_pos), mean_pctM, .fun = mean)) %>%
            ggplot(aes(x = mean_pctM, y = harmonized_pos, color = condition, size = n_elements)) +
            geom_point(alpha = 0.8) +
            scale_size_continuous(range = c(1, 5)) +
            labs(
                x = "Mean % m6A", y = "Harmonized position",
                size = "# elements", color = "Condition",
                title = sprintf("Reliable m6A sites by condition [%s]", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/site_dotplot_reliable.pdf"),
            w = 9, h = max(4, nrow(reliable_sites) * 0.25)
        )
    }

    # --- Plot 9: IVT vs Sample signal comparison for FP positions ---
    # Some FP positions may have genuine signal far exceeding the IVT artifact
    # Use pre-FP-filtered data since both meth_data and filtered_meth already excluded FP positions
    unfiltered_source <- if (length_filter == "all") yl1annm_alllength else yl1annm
    unfiltered_meth <- unfiltered_source %>% filter(rte_subfamily %in% subfamilies)
    if (length_filter == "FL") unfiltered_meth <- unfiltered_meth %>% filter(rte_length_req == "FL")

    sample_at_fp <- unfiltered_meth %>%
        filter(cov >= 5, harmonized_pos %in% fp_harmonized_positions) %>%
        group_by(harmonized_pos) %>%
        summarise(
            sample_mean_pctM = weighted.mean(pctM, cov, na.rm = TRUE),
            sample_median_pctM = median(pctM, na.rm = TRUE),
            n_elements = n_distinct(gene_id),
            n_samples = n_distinct(sample),
            total_cov = sum(cov),
            .groups = "drop"
        ) %>%
        left_join(hp_to_cp, by = "harmonized_pos")

    ivt_at_fp <- ivt_profile %>%
        filter(cov >= FP_MIN_COV) %>%
        dplyr::select(consensus_pos, ivt_pctM = pctM)

    fp_comparison <- sample_at_fp %>%
        filter(!is.na(consensus_pos)) %>%
        left_join(ivt_at_fp, by = "consensus_pos") %>%
        filter(!is.na(ivt_pctM)) %>%
        mutate(
            ratio_sample_ivt = sample_mean_pctM / pmax(ivt_pctM, 0.1),
            diff_sample_ivt = sample_mean_pctM - ivt_pctM,
            potentially_real = ratio_sample_ivt > 3 & sample_mean_pctM > 15
        )

    if (nrow(fp_comparison) > 0) {
        p <- fp_comparison %>%
            ggplot(aes(x = ivt_pctM, y = sample_mean_pctM)) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
            geom_abline(slope = 3, intercept = 0, linetype = "dotted", color = "blue", alpha = 0.4) +
            geom_point(aes(color = potentially_real, size = n_elements), alpha = 0.7) +
            scale_color_manual(
                values = c("TRUE" = "#E64B35", "FALSE" = "grey50"),
                labels = c("TRUE" = "Sample >> IVT (ratio > 3x)", "FALSE" = "Likely artifact")
            ) +
            scale_size_continuous(range = c(1.5, 5)) +
            geom_text_repel(
                data = . %>% filter(potentially_real | sample_mean_pctM > 20),
                aes(label = sprintf("hp%d (%.0f%%)", harmonized_pos, sample_mean_pctM)),
                size = 2.5, max.overlaps = 20
            ) +
            labs(
                x = "% m6A in IVT (artifact level)",
                y = sprintf("Mean %% m6A in samples (pooled, %s)", full_label),
                color = NULL, size = "# elements",
                title = "FP positions: sample signal vs IVT artifact",
                subtitle = "Blue dotted line = 3x IVT level. Points above may have real signal despite artifact."
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/fp_sample_vs_ivt_signal.pdf"), w = 9, h = 7)

        # Also save a table of potentially rescuable FP positions
        if (any(fp_comparison$potentially_real)) {
            write_tsv(
                fp_comparison %>% filter(potentially_real) %>% arrange(desc(ratio_sample_ivt)),
                str_glue("{intdir}/fp_potentially_rescuable.tsv")
            )
            cat(sprintf(
                "Potentially rescuable FP positions (sample >> IVT): %d\n",
                sum(fp_comparison$potentially_real)
            ))
        }
    }

    # --- Plot 10: Modification profile colored by DRACH context ---
    # DRACH motif: D=[AGT], R=[AG], A, C, H=[ACT]
    # Scan L1HS consensus for DRACH, map to harmonized positions
    l1hs_seq_str <- as.character(l1hs_cs[[1]])
    drach_pattern <- DNAString("DRACH")
    drach_hits <- matchPattern(drach_pattern, l1hs_cs[[1]], fixed = FALSE)
    # The A in DRACH is at position 3 of the 5-mer
    drach_a_positions <- start(drach_hits) + 2L
    # Map to harmonized positions
    drach_harmonized <- l1hs_mapping %>%
        filter(consensus_pos %in% drach_a_positions) %>%
        pull(harmonized_pos) %>%
        unique()
    cat(sprintf(
        "DRACH motif A-positions in L1HS consensus: %d, mapped to %d harmonized positions\n",
        length(drach_a_positions), length(drach_harmonized)
    ))

    site_classification_drach <- site_classification %>%
        mutate(is_drach = harmonized_pos %in% drach_harmonized)

    p <- site_classification_drach %>%
        filter(!is.na(consensus_pos)) %>%
        ggplot(aes(x = consensus_pos, y = mean_pctM)) +
        geom_point(aes(color = interaction(reliable, is_drach), size = mean_n_elements), alpha = 0.7) +
        scale_color_manual(
            values = c(
                "TRUE.TRUE" = "#E64B35", "TRUE.FALSE" = "#4DBBD5",
                "FALSE.TRUE" = "#F39B7F", "FALSE.FALSE" = "grey80"
            ),
            labels = c(
                "TRUE.TRUE" = "Reliable + DRACH", "TRUE.FALSE" = "Reliable + non-DRACH",
                "FALSE.TRUE" = "Non-reliable + DRACH", "FALSE.FALSE" = "Non-reliable + non-DRACH"
            )
        ) +
        scale_size_continuous(range = c(0.5, 4)) +
        geom_text_repel(
            data = . %>% filter(reliable),
            aes(label = sprintf("%d", harmonized_pos)),
            size = 2.5, max.overlaps = 25
        ) +
        annotate("rect",
            xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
            ymin = -4, ymax = -2, fill = L1_REGIONS$color, alpha = 0.5
        ) +
        annotate("text",
            x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
            y = -3, label = L1_REGIONS$region, size = 2.5
        ) +
        labs(
            x = "L1HS consensus position", y = "Mean % m6A (pooled)",
            color = NULL, size = "# elements",
            title = sprintf("m6A profile: DRACH context [%s]", full_label),
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_sites}/modification_profile_drach.pdf"), w = 14, h = 6)

    # Summary of DRACH vs non-DRACH among reliable sites
    drach_summary <- site_classification_drach %>%
        filter(reliable) %>%
        count(is_drach) %>%
        mutate(pct = 100 * n / sum(n))
    cat(sprintf(
        "Reliable sites in DRACH context: %d/%d (%.1f%%)\n",
        sum(drach_summary$n[drach_summary$is_drach]),
        sum(drach_summary$n),
        ifelse(any(drach_summary$is_drach), drach_summary$pct[drach_summary$is_drach], 0)
    ))

    # --- Plot 11: Motif enrichment - modified vs unmodified DRACH sites ---
    # Compare sequence context around modified vs unmodified DRACH A-positions
    # Extract windows around each DRACH A position
    WINDOW_SIZE <- 10 # bases on each side of the A
    l1hs_len <- nchar(l1hs_seq_str)

    all_drach_df <- tibble(
        consensus_pos = drach_a_positions
    ) %>%
        left_join(l1hs_mapping %>% dplyr::select(consensus_pos, harmonized_pos), by = "consensus_pos") %>%
        filter(!is.na(harmonized_pos)) %>%
        mutate(
            is_modified = harmonized_pos %in% reliable_sites$harmonized_pos,
            window_start = pmax(1L, consensus_pos - WINDOW_SIZE),
            window_end = pmin(l1hs_len, consensus_pos + WINDOW_SIZE),
            window_seq = substring(l1hs_seq_str, window_start, window_end),
            # Centered window (may be shorter at edges)
            center_offset = consensus_pos - window_start + 1L
        )

    # Filter to windows with full length
    full_window <- all_drach_df %>%
        filter(nchar(window_seq) == 2 * WINDOW_SIZE + 1)

    if (nrow(full_window) > 0 && any(full_window$is_modified) && any(!full_window$is_modified)) {
        # Build position frequency matrices for modified vs unmodified
        mod_seqs <- DNAStringSet(full_window$window_seq[full_window$is_modified])
        unmod_seqs <- DNAStringSet(full_window$window_seq[!full_window$is_modified])

        mod_pfm <- consensusMatrix(mod_seqs, as.prob = TRUE)[1:4, ]
        unmod_pfm <- consensusMatrix(unmod_seqs, as.prob = TRUE)[1:4, ]

        # Convert to long format for plotting
        pfm_to_long <- function(pfm, group_label) {
            as_tibble(t(pfm)) %>%
                mutate(position = seq_len(ncol(pfm)) - (WINDOW_SIZE + 1L)) %>%
                pivot_longer(-position, names_to = "base", values_to = "freq") %>%
                mutate(group = group_label)
        }

        pfm_long <- bind_rows(
            pfm_to_long(mod_pfm, sprintf("Modified DRACH (n=%d)", length(mod_seqs))),
            pfm_to_long(unmod_pfm, sprintf("Unmodified DRACH (n=%d)", length(unmod_seqs)))
        )

        p <- pfm_long %>%
            ggplot(aes(x = position, y = freq, fill = base)) +
            geom_col(position = "stack") +
            facet_wrap(~group, ncol = 1) +
            scale_fill_manual(values = c(A = "#E64B35", C = "#4DBBD5", G = "#00A087", T = "#F39B7F")) +
            geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
            labs(
                x = "Position relative to modified A", y = "Base frequency",
                fill = "Base",
                title = sprintf("Sequence context: modified vs unmodified DRACH sites [%s]", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/drach_motif_context_comparison.pdf"), w = 10, h = 7)

        # Per-position enrichment: Fisher's test at each position for each base
        enrichment_results <- list()
        for (pos_idx in seq_len(2 * WINDOW_SIZE + 1)) {
            rel_pos <- pos_idx - (WINDOW_SIZE + 1L)
            for (base in c("A", "C", "G", "T")) {
                n_mod_base <- mod_pfm[base, pos_idx] * length(mod_seqs)
                n_mod_other <- length(mod_seqs) - n_mod_base
                n_unmod_base <- unmod_pfm[base, pos_idx] * length(unmod_seqs)
                n_unmod_other <- length(unmod_seqs) - n_unmod_base
                ft <- fisher.test(matrix(c(
                    round(n_mod_base), round(n_mod_other),
                    round(n_unmod_base), round(n_unmod_other)
                ), nrow = 2))
                enrichment_results[[length(enrichment_results) + 1]] <- tibble(
                    position = rel_pos,
                    base = base,
                    odds_ratio = ft$estimate,
                    p_value = ft$p.value,
                    mod_freq = mod_pfm[base, pos_idx],
                    unmod_freq = unmod_pfm[base, pos_idx]
                )
            }
        }
        enrichment_df <- bind_rows(enrichment_results) %>%
            mutate(p_adj = p.adjust(p_value, method = "BH"))

        write_tsv(enrichment_df, str_glue("{intdir}/drach_motif_enrichment.tsv"))

        # Plot enrichment/depletion
        p <- enrichment_df %>%
            mutate(log2_or = log2(pmax(odds_ratio, 0.01))) %>%
            ggplot(aes(x = position, y = log2_or, fill = base)) +
            geom_col(position = "dodge") +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
            scale_fill_manual(values = c(A = "#E64B35", C = "#4DBBD5", G = "#00A087", T = "#F39B7F")) +
            geom_point(
                data = . %>% filter(p_adj < 0.05),
                aes(x = position, y = log2_or),
                shape = 8, size = 1.5, color = "black", show.legend = FALSE,
                position = position_dodge(width = 0.9)
            ) +
            labs(
                x = "Position relative to modified A", y = "log2(Odds Ratio)",
                fill = "Base",
                title = sprintf("Base enrichment: modified vs unmodified DRACH [%s]", full_label),
                subtitle = "Stars = BH-adjusted p < 0.05"
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/drach_motif_enrichment.pdf"), w = 10, h = 5)

        # Summary statistics
        sig_positions <- enrichment_df %>% filter(p_adj < 0.05)
        cat(sprintf(
            "Motif enrichment: %d significant base-position combinations (padj < 0.05)\n",
            nrow(sig_positions)
        ))
        if (nrow(sig_positions) > 0) {
            cat("Top enriched/depleted:\n")
            top_sig <- sig_positions %>%
                arrange(p_adj) %>%
                head(10)
            for (i in seq_len(nrow(top_sig))) {
                cat(sprintf(
                    "  Position %+d, base %s: OR=%.2f, padj=%.3g\n",
                    top_sig$position[i], top_sig$base[i], top_sig$odds_ratio[i], top_sig$p_adj[i]
                ))
            }
        }

        write_tsv(all_drach_df, str_glue("{intdir}/drach_positions_annotated.tsv"))
    } else {
        cat("Insufficient DRACH sites for motif enrichment (need both modified and unmodified)\n")
    }

    # --- FACETED BY LENGTH PLOTS (all-length mode only) ---
    if (length_filter == "all") {
        by_length <- filtered_meth %>%
            filter(rte_subfamily %in% subfamilies, cov >= 5) %>%
            group_by(condition, harmonized_pos, rte_length_req) %>%
            summarise(
                mean_pctM = mean(pctM, na.rm = TRUE),
                n_elements = n_distinct(gene_id),
                .groups = "drop"
            ) %>%
            left_join(hp_to_cp, by = "harmonized_pos")

        p <- by_length %>%
            filter(!is.na(consensus_pos)) %>%
            ggplot(aes(x = consensus_pos, y = mean_pctM, color = rte_length_req)) +
            geom_point(aes(size = n_elements), alpha = 0.5) +
            geom_smooth(method = "loess", se = FALSE, span = 0.1) +
            scale_color_manual(values = c("FL" = "#D62728", "Trnc" = "#1F77B4")) +
            scale_size_continuous(range = c(0.5, 3)) +
            geom_rect(
                data = L1_REGIONS, inherit.aes = FALSE,
                aes(xmin = start, xmax = end, ymin = -4, ymax = -2, fill = region), alpha = 0.5
            ) +
            scale_fill_manual(values = setNames(L1_REGIONS$color, L1_REGIONS$region), guide = "none") +
            geom_text(
                data = L1_REGIONS, inherit.aes = FALSE,
                aes(x = (start + end) / 2, y = -3, label = region), size = 2.5
            ) +
            facet_wrap(~condition, ncol = 1) +
            labs(
                x = "L1HS consensus position", y = "Mean % m6A",
                color = "Length", size = "# elements",
                title = sprintf("m6A profile: FL vs truncated (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/modification_profile_by_length.pdf"),
            w = 14, h = 5 * length(conditions_present)
        )

        # Faceted by condition x length
        p <- by_length %>%
            filter(!is.na(consensus_pos)) %>%
            ggplot(aes(x = consensus_pos, y = mean_pctM)) +
            geom_point(aes(size = n_elements), alpha = 0.5, color = "grey40") +
            geom_smooth(method = "loess", se = FALSE, span = 0.1, color = "steelblue") +
            geom_point(
                data = . %>% filter(harmonized_pos %in% reliable_sites$harmonized_pos),
                color = "#D62728", size = 2
            ) +
            scale_size_continuous(range = c(0.5, 3)) +
            facet_grid(condition ~ rte_length_req) +
            geom_rect(
                data = L1_REGIONS, inherit.aes = FALSE,
                aes(xmin = start, xmax = end, ymin = -4, ymax = -2, fill = region), alpha = 0.5
            ) +
            scale_fill_manual(values = setNames(L1_REGIONS$color, L1_REGIONS$region), guide = "none") +
            geom_text(
                data = L1_REGIONS, inherit.aes = FALSE,
                aes(x = (start + end) / 2, y = -3, label = region), size = 2.5
            ) +
            labs(
                x = "L1HS consensus position", y = "Mean % m6A",
                size = "# elements",
                title = sprintf("m6A profile by condition x length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/modification_profile_faceted_condition_x_length.pdf"),
            w = 14, h = 4 * length(conditions_present)
        )

        # Per-reliable-site: FL vs trnc comparison, per condition
        site_by_length <- filtered_meth %>%
            filter(rte_subfamily %in% subfamilies) %>%
            filter(harmonized_pos %in% reliable_sites$harmonized_pos) %>%
            filter(cov >= 5) %>%
            group_by(condition, harmonized_pos, rte_length_req) %>%
            summarise(
                mean_pctM = mean(pctM, na.rm = TRUE),
                n_elements = n_distinct(gene_id),
                .groups = "drop"
            )

        if (nrow(site_by_length) > 0) {
            p <- site_by_length %>%
                mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
                ggplot(aes(x = harmonized_pos, y = mean_pctM, fill = rte_length_req)) +
                geom_col(position = "dodge") +
                geom_text(aes(label = n_elements), position = position_dodge(width = 0.9), vjust = -0.3, size = 2) +
                scale_fill_manual(values = c("FL" = "#D62728", "Trnc" = "#1F77B4")) +
                facet_wrap(~condition, ncol = 1) +
                labs(
                    x = "Reliable site (harmonized position)", y = "Mean % m6A",
                    fill = "Length",
                    title = sprintf("FL vs truncated at reliable sites (%s)", full_label)
                ) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            mysaveandstore(str_glue("{outdir_sites}/reliable_sites_fl_vs_trnc.pdf"),
                w = max(7, nrow(reliable_sites) * 0.8 + 2), h = 4 * length(conditions_present)
            )
        }

        # Element count comparison
        element_counts <- filtered_meth %>%
            filter(rte_subfamily %in% subfamilies, cov >= 5) %>%
            distinct(gene_id, rte_length_req, rte_subfamily, condition) %>%
            count(condition, rte_length_req, rte_subfamily)

        p <- element_counts %>%
            ggplot(aes(x = rte_subfamily, y = n, fill = rte_length_req)) +
            geom_col(position = "dodge") +
            geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5) +
            scale_fill_manual(values = c("FL" = "#D62728", "Trnc" = "#1F77B4")) +
            facet_wrap(~condition) +
            labs(
                x = "Subfamily", y = "# elements with data",
                fill = "Length",
                title = sprintf("Elements by subfamily and length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_sites}/element_counts_by_length.pdf"),
            w = 5 * length(conditions_present), h = 5
        )
    }


    # ============================================================
    # //ANCHOR - SECTION 3: MIXED EFFECTS MODELING (glmmTMB)
    # ============================================================
    cat(sprintf("\n=== Section 3: Mixed Effects Modeling [%s] ===\n", full_label))

    # 3a. Map read-level data to harmonized positions (all conditions pooled)
    reads_filtered <- reads %>%
        filter(!is.na(gene_id))

    rmann_sub <- rmann %>% filter(rte_subfamily %in% subfamilies)
    if (length_filter == "FL") {
        rmann_sub <- rmann_sub %>% filter(rte_length_req == "FL")
    }
    reads_with_element <- reads_filtered %>%
        left_join(
            rmann_sub %>%
                dplyr::select(gene_id, rte_start = start, rte_end = end, rte_strand = strand),
            by = "gene_id"
        ) %>%
        filter(!is.na(rte_start)) %>%
        mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end)))

    reads_mapped <- reads_with_element %>%
        left_join(
            consensus_index_long %>%
                dplyr::select(gene_id, sequence_pos, consensus_pos) %>%
                distinct(gene_id, sequence_pos, .keep_all = TRUE),
            by = c("gene_id", "sequence_pos")
        ) %>%
        filter(!is.na(consensus_pos)) %>%
        left_join(
            l1hs_mapping %>% dplyr::select(consensus_pos, harmonized_pos),
            by = "consensus_pos"
        ) %>%
        filter(!is.na(harmonized_pos)) %>%
        filter(harmonized_pos %in% reliable_sites$harmonized_pos) %>%
        filter(!(harmonized_pos %in% fp_harmonized_positions))

    # 3b. Binarize modification
    reads_mapped <- reads_mapped %>%
        mutate(is_mod = ifelse(mod_qual > 0.5, 1, 0))

    cat(sprintf("Read-level observations at reliable sites: %d\n", nrow(reads_mapped)))
    cat(sprintf(
        "Unique reads: %d, Unique elements: %d, Conditions: %s\n",
        n_distinct(reads_mapped$read_id), n_distinct(reads_mapped$gene_id),
        paste(unique(reads_mapped$condition), collapse = ", ")
    ))

    # 3c. Pivot to wide format
    reads_wide <- reads_mapped %>%
        dplyr::select(read_id, gene_id, sample, condition, harmonized_pos, is_mod) %>%
        distinct(read_id, harmonized_pos, .keep_all = TRUE) %>%
        pivot_wider(
            names_from = harmonized_pos,
            values_from = is_mod,
            names_prefix = "site_"
        )

    # Keep reads spanning at least 50% of reliable sites
    site_cols <- grep("^site_", colnames(reads_wide), value = TRUE)
    n_reliable <- length(site_cols)
    reads_wide <- reads_wide %>%
        mutate(n_sites_covered = rowSums(!is.na(across(all_of(site_cols))))) %>%
        filter(n_sites_covered >= 0.5 * n_reliable)

    cat(sprintf("Reads spanning >= 50%% of reliable sites: %d\n", nrow(reads_wide)))

    # 3d. Join element-level annotations
    reads_annotated <- reads_wide %>%
        left_join(
            rmann %>%
                dplyr::select(gene_id, intactness_req, genic_loc, loc_superlowres_integrative_stranded, rte_subfamily, rte_length_req) %>%
                distinct(),
            by = "gene_id"
        ) %>%
        left_join(
            enrich_data %>%
                dplyr::select(
                    gene_id, starts_with("l2fc_"), starts_with("baseMean_"),
                    starts_with("mean_lr_")
                ),
            by = "gene_id"
        )

    # Scale continuous predictors (use first available contrast)
    l2fc_col <- grep("^l2fc_", colnames(reads_annotated), value = TRUE)[1]
    basemean_col <- grep("^baseMean_", colnames(reads_annotated), value = TRUE)[1]
    if (!is.na(l2fc_col)) {
        reads_annotated <- reads_annotated %>%
            mutate(l2fc_scaled = as.numeric(scale(.data[[l2fc_col]])))
    }
    if (!is.na(basemean_col)) {
        reads_annotated <- reads_annotated %>%
            mutate(baseMean_scaled = as.numeric(scale(log1p(.data[[basemean_col]]))))
    }

    # Cap reads per element to prevent dominant elements
    MAX_READS_PER_ELEMENT <- 100
    reads_annotated_sub <- reads_annotated %>%
        group_by(gene_id) %>%
        dplyr::slice_sample(prop = 1) %>%
        dplyr::slice_head(n = MAX_READS_PER_ELEMENT) %>%
        ungroup()

    cat(sprintf("Reads after subsampling (max %d per element): %d\n", MAX_READS_PER_ELEMENT, nrow(reads_annotated_sub)))

    write_tsv(reads_annotated_sub %>% dplyr::select(-n_sites_covered), str_glue("{intdir}/reads_wide_annotated.tsv"))

    # 3e. Return modeling inputs
    modeling_inputs <- list(
        reads_annotated_sub = reads_annotated_sub,
        site_cols = site_cols,
        reliable_sites = reliable_sites
    )

    # Return Section 2 results + modeling inputs
    list(
        reliable_sites = reliable_sites,
        site_classification = site_classification,
        modeling_inputs = modeling_inputs
    )
} # end run_sections_2_and_3

# --- Run Sections 2+3 data prep for all subfamily/length combinations ---
# FL-only (original)
results_all_yl1 <- run_sections_2_and_3(subfamilies = YOUNG_L1_SUBFAMILIES, label = "Young L1", suffix = "", length_filter = "FL")
results_l1hs <- run_sections_2_and_3(subfamilies = "L1HS", label = "L1HS", suffix = "_L1HS", length_filter = "FL")
# All-length (FL + trnc)
results_all_yl1_alllength <- run_sections_2_and_3(YOUNG_L1_SUBFAMILIES, label = "Young L1", suffix = "_alllength", length_filter = "all")
results_l1hs_alllength <- run_sections_2_and_3(subfamilies = "L1HS", label = "L1HS", suffix = "_L1HS_alllength", length_filter = "all")

# ============================================================
# //ANCHOR - PER-SITE m6A DOTPLOT + INTRA-READ CORRELATION HEATMAP
# ============================================================
# Uses reads_wide_annotated from Section 2/3 (FL L1HS)

outdir_m6a_sites <- str_glue("{outdir_base}/per_site_plots")
dir.create(outdir_m6a_sites, recursive = TRUE, showWarnings = FALSE)

# Load read-level data at reliable sites
reads_ann <- read_tsv(str_glue("{intdir}_L1HS/reads_wide_annotated.tsv"))
reliable_sites_l1hs <- results_l1hs$reliable_sites
site_cols <- grep("^site_", colnames(reads_ann), value = TRUE)

# Derive condition label (ORF1 vs TOT)
reads_ann <- reads_ann %>%
    mutate(bcondition = case_when(
        grepl("ORF1", condition) ~ "ORF1",
        grepl("TOT", condition) ~ "TOT",
        TRUE ~ condition
    ))

# --- 1. Per-site dotplot: %m6A for ORF1 vs TOT ---
# Pivot back to long for per-site stats
reads_long <- reads_ann %>%
    pivot_longer(
        cols = all_of(site_cols),
        names_to = "site", values_to = "is_mod",
        names_prefix = "site_"
    ) %>%
    filter(!is.na(is_mod)) %>%
    mutate(harmonized_pos = as.numeric(site))

site_stats <- reads_long %>%
    group_by(harmonized_pos, bcondition) %>%
    summarise(
        pctM = mean(is_mod) * 100,
        n_reads = n(),
        .groups = "drop"
    )

site_stats_wide <- site_stats %>%
    pivot_wider(
        names_from = bcondition,
        values_from = c(pctM, n_reads),
        names_sep = "_"
    ) %>%
    filter(!is.na(pctM_ORF1) & !is.na(pctM_TOT))

# Join with consensus_pos for x-axis labeling
if (exists("hp_to_cp")) {
    site_stats_wide <- site_stats_wide %>%
        left_join(hp_to_cp, by = "harmonized_pos")
}

p <- ggplot(site_stats_wide, aes(
    x = pctM_TOT, y = pctM_ORF1,
    size = log10(n_reads_ORF1 + n_reads_TOT)
)) +
    geom_point(alpha = 0.7, color = "#4393C3") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = harmonized_pos), size = 2.5, max.overlaps = 20) +
    scale_size_continuous(name = "log10(total cov)", range = c(1, 6)) +
    labs(
        x = "% m6A (TOT)", y = "% m6A (ORF1)",
        title = "L1HS FL: Per-site m6A — ORF1 vs TOT"
    ) +
    coord_equal() +
    mtclosed
mss(
    pl = p, fn = str_glue("{outdir_m6a_sites}/per_site_m6a_orf1_vs_tot.pdf"),
    fw = 2, fh = 2, plus_void = TRUE
)

# Bar version: each site as a bar, ORF1 and TOT side by side
site_stats_plot <- site_stats %>%
    mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))))

p <- ggplot(site_stats_plot, aes(x = harmonized_pos, y = pctM, fill = bcondition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.0f", pctM)),
        position = position_dodge(0.8),
        vjust = -0.5, size = 2
    ) +
    labs(
        x = "Harmonized position", y = "% m6A", fill = "Condition",
        title = "L1HS FL: Per-site m6A levels"
    ) +
    mtclosed
mss(
    pl = p, fn = str_glue("{outdir_m6a_sites}/per_site_m6a_barplot.pdf"),
    fw = 3, fh = 1.5, plus_void = TRUE
)

# --- 2. Intra-read pairwise correlation heatmap ---
# For each pair of reliable sites, compute correlation of is_mod across reads
site_mat <- reads_ann %>%
    dplyr::select(all_of(site_cols)) %>%
    as.matrix()

# Pairwise phi (Pearson on binary) correlation
cor_mat <- cor(site_mat, use = "pairwise.complete.obs")
colnames(cor_mat) <- gsub("^site_", "", colnames(cor_mat))
rownames(cor_mat) <- gsub("^site_", "", rownames(cor_mat))

# Order by position
pos_order <- order(as.numeric(colnames(cor_mat)))
cor_mat <- cor_mat[pos_order, pos_order]

p <- Heatmap(
    cor_mat,
    name = "Pearson r",
    col = circlize::colorRamp2(c(-0.3, 0, 0.3, 0.6, 1), c("#2166AC", "#F7F7F7", "#F7F7F7", "#FDDBC7", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_title = "L1HS FL: Intra-read m6A correlation (all conditions)",
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (i != j) grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 5))
    }
)
pdf(str_glue("{outdir_m6a_sites}/intraread_correlation_heatmap_all.pdf"), width = 8, height = 7)
draw(p)
dev.off()

# Per-condition correlation heatmaps
for (cond in c("ORF1", "TOT")) {
    cond_reads <- reads_ann %>%
        filter(bcondition == cond) %>%
        dplyr::select(all_of(site_cols)) %>%
        as.matrix()

    if (nrow(cond_reads) < 10) next

    cor_cond <- cor(cond_reads, use = "pairwise.complete.obs")
    colnames(cor_cond) <- gsub("^site_", "", colnames(cor_cond))
    rownames(cor_cond) <- gsub("^site_", "", rownames(cor_cond))
    cor_cond <- cor_cond[pos_order, pos_order]

    p <- Heatmap(
        cor_cond,
        name = "Pearson r",
        col = circlize::colorRamp2(c(-0.3, 0, 0.3, 0.6, 1), c("#2166AC", "#F7F7F7", "#F7F7F7", "#FDDBC7", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_title = str_glue("L1HS FL: Intra-read m6A correlation ({cond})"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (i != j) grid.text(sprintf("%.2f", cor_cond[i, j]), x, y, gp = gpar(fontsize = 5))
        }
    )
    pdf(str_glue("{outdir_m6a_sites}/intraread_correlation_heatmap_{tolower(cond)}.pdf"), width = 8, height = 7)
    draw(p)
    dev.off()
}

# Difference heatmap: ORF1 - TOT correlation
orf1_reads <- reads_ann %>%
    filter(bcondition == "ORF1") %>%
    dplyr::select(all_of(site_cols)) %>%
    as.matrix()
tot_reads <- reads_ann %>%
    filter(bcondition == "TOT") %>%
    dplyr::select(all_of(site_cols)) %>%
    as.matrix()
if (nrow(orf1_reads) >= 10 && nrow(tot_reads) >= 10) {
    cor_orf1 <- cor(orf1_reads, use = "pairwise.complete.obs")
    cor_tot <- cor(tot_reads, use = "pairwise.complete.obs")
    cor_diff <- cor_orf1 - cor_tot
    colnames(cor_diff) <- gsub("^site_", "", colnames(cor_diff))
    rownames(cor_diff) <- gsub("^site_", "", rownames(cor_diff))
    cor_diff <- cor_diff[pos_order, pos_order]

    p <- Heatmap(
        cor_diff,
        name = "Δr (ORF1-TOT)",
        col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "#F7F7F7", "#B2182B")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_title = "L1HS FL: Correlation difference (ORF1 - TOT)",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (i != j) grid.text(sprintf("%.2f", cor_diff[i, j]), x, y, gp = gpar(fontsize = 5))
        }
    )
    pdf(str_glue("{outdir_m6a_sites}/intraread_correlation_diff_orf1_minus_tot.pdf"), width = 8, height = 7)
    draw(p)
    dev.off()
}


# ============================================================
# //ANCHOR - SECTION 3: GLMM FITTING + PLOTS (at global scope for parallelism)
# ============================================================

run_glmm <- function(modeling_inputs, suffix = "") {
    reads_annotated_sub <- modeling_inputs$reads_annotated_sub
    site_cols <- modeling_inputs$site_cols
    reliable_sites <- modeling_inputs$reliable_sites
    outdir_glmm_local <- str_glue("{outdir_glmm}{suffix}")
    intdir_local <- str_glue("{intdir}{suffix}")
    for (d in c(outdir_glmm_local, intdir_local)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    # Shadow mysaveandstore to capture function-local p instead of global p
    .global_mysaveandstore <- get("mysaveandstore", envir = globalenv())
    mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, store = store_var, raster = FALSE, ...) {
        .global_mysaveandstore(fn = fn, w = w, h = h, res = res, pl = p, store = store, raster = raster, ...)
    }

    cat(sprintf("\n=== Fitting GLMMs%s ===\n", ifelse(suffix == "", "", sprintf(" [%s]", suffix))))

    # Fit full models (parallel)
    plan(multisession, workers = min(length(site_cols), 16))

    model_results <- site_cols %>%
        purrr::set_names() %>%
        future_map(~ {
            fit_site_model(.x, reads_annotated_sub, site_cols)
        }, .progress = TRUE)

    plan(sequential)

    n_success <- sum(!sapply(model_results, is.null))
    cat(sprintf("Models fit: %d / %d\n", n_success, length(site_cols)))
    saveRDS(model_results, str_glue("{intdir_local}/glmm_models.rds"))

    # Fit context-only models (parallel)
    plan(multisession, workers = min(length(site_cols), 8))

    model_results_nosite <- site_cols %>%
        purrr::set_names() %>%
        future_map(~ {
            fit_site_model(.x, reads_annotated_sub, site_cols, include_other_sites = FALSE)
        }, .progress = TRUE)

    plan(sequential)

    n_success_nosite <- sum(!sapply(model_results_nosite, is.null))
    cat(sprintf("Context-only models fit: %d / %d\n", n_success_nosite, length(site_cols)))
    saveRDS(model_results_nosite, str_glue("{intdir_local}/glmm_models_nosite.rds"))

    # Extract coefficients (full models)
    coef_table <- model_results %>%
        keep(~ !is.null(.x)) %>%
        imap_dfr(~ {
            tryCatch(
                broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE) %>%
                    mutate(target_site = .y),
                error = function(e) tibble()
            )
        })

    ranef_table <- model_results %>%
        keep(~ !is.null(.x)) %>%
        imap_dfr(~ {
            tryCatch(
                broom.mixed::tidy(.x, effects = "ran_pars") %>%
                    mutate(target_site = .y),
                error = function(e) tibble()
            )
        })

    if (nrow(coef_table) > 0) {
        coef_table <- coef_table %>%
            group_by(term) %>%
            mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
            ungroup()
    }

    write_tsv(coef_table, str_glue("{intdir_local}/glmm_coef_table.tsv"))
    write_tsv(ranef_table, str_glue("{intdir_local}/glmm_ranef_table.tsv"))

    cat(sprintf(
        "Fixed effects: %d, Significant (padj < 0.05): %d\n",
        nrow(coef_table), sum(coef_table$p_adj < 0.05, na.rm = TRUE)
    ))

    # Extract coefficients (context-only models)
    coef_table_nosite <- model_results_nosite %>%
        keep(~ !is.null(.x)) %>%
        imap_dfr(~ {
            tryCatch(
                broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE) %>%
                    mutate(target_site = .y),
                error = function(e) tibble()
            )
        })

    ranef_table_nosite <- model_results_nosite %>%
        keep(~ !is.null(.x)) %>%
        imap_dfr(~ {
            tryCatch(
                broom.mixed::tidy(.x, effects = "ran_pars") %>%
                    mutate(target_site = .y),
                error = function(e) tibble()
            )
        })

    if (nrow(coef_table_nosite) > 0) {
        coef_table_nosite <- coef_table_nosite %>%
            group_by(term) %>%
            mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
            ungroup()
    }

    write_tsv(coef_table_nosite, str_glue("{intdir_local}/glmm_coef_table_nosite.tsv"))
    write_tsv(ranef_table_nosite, str_glue("{intdir_local}/glmm_ranef_table_nosite.tsv"))

    # --- GLMM PLOTS ---
    # Use outdir_glmm_local and intdir_local for all paths
    outdir_glmm <- outdir_glmm_local
    intdir <- intdir_local

    if (nrow(coef_table) > 0 && n_success > 0) {
        coef_table <- coef_table %>%
            mutate(
                predictor_category = case_when(
                    grepl("^site_", term) ~ "Co-modification",
                    term %in% c("l2fc_scaled", "baseMean_scaled", "orf1_expression_scaled") ~ "Expression",
                    term %in% c("intactness_req", "genic_loc", "loc_superlowres_integrative_stranded") ~ "Genomic context",
                    term == "(Intercept)" ~ "Intercept",
                    grepl("^sample", term) ~ "Sample",
                    TRUE ~ "Other"
                )
            )

        # Plot 1: Coefficient heatmap
        coef_for_hm <- coef_table %>%
            filter(term != "(Intercept)", !grepl("^sample", term)) %>%
            dplyr::select(target_site, term, estimate) %>%
            pivot_wider(names_from = target_site, values_from = estimate) %>%
            column_to_rownames("term") %>%
            as.matrix()

        sig_for_hm <- coef_table %>%
            filter(term != "(Intercept)", !grepl("^sample", term)) %>%
            dplyr::select(target_site, term, p_adj) %>%
            pivot_wider(names_from = target_site, values_from = p_adj) %>%
            column_to_rownames("term") %>%
            as.matrix()

        star_mat <- matrix("", nrow = nrow(sig_for_hm), ncol = ncol(sig_for_hm))
        star_mat[sig_for_hm < 0.05] <- "*"
        star_mat[sig_for_hm < 0.01] <- "**"
        star_mat[sig_for_hm < 0.001] <- "***"

        if (nrow(coef_for_hm) > 1 && ncol(coef_for_hm) > 1) {
            row_cats <- coef_table %>%
                filter(term != "(Intercept)", !grepl("^sample", term)) %>%
                dplyr::select(term, predictor_category) %>%
                distinct() %>%
                column_to_rownames("term")
            row_cats <- row_cats[rownames(coef_for_hm), , drop = FALSE]

            col_ann <- reliable_sites %>%
                mutate(site_name = paste0("site_", harmonized_pos)) %>%
                filter(site_name %in% colnames(coef_for_hm)) %>%
                dplyr::select(site_name, mean_pctM) %>%
                column_to_rownames("site_name")
            col_ann <- col_ann[colnames(coef_for_hm), , drop = FALSE]

            max_abs <- min(max(abs(coef_for_hm), na.rm = TRUE), 5)

            p <- Heatmap(coef_for_hm,
                name = "log-odds",
                col = circlize::colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red")),
                na_col = "grey90",
                cluster_rows = TRUE, cluster_columns = TRUE,
                left_annotation = rowAnnotation(Category = row_cats$predictor_category),
                top_annotation = HeatmapAnnotation(
                    `Mean m6A` = col_ann$mean_pctM,
                    col = list(`Mean m6A` = circlize::colorRamp2(c(0, 50, 100), c("white", "orange", "red")))
                ),
                cell_fun = function(j, i, x, y, width, height, fill) {
                    if (!is.na(star_mat[i, j]) && star_mat[i, j] != "") {
                        grid::grid.text(star_mat[i, j], x, y, gp = grid::gpar(fontsize = 8))
                    }
                },
                column_title = "GLMM coefficients: predictors of m6A at each site"
            )
            mysaveandstore(str_glue("{outdir_glmm}/coefficient_heatmap.pdf"),
                w = max(8, ncol(coef_for_hm) * 0.8 + 4), h = max(6, nrow(coef_for_hm) * 0.4 + 2)
            )
        }

        # Plot 2: Forest plots per target site
        for (site in unique(coef_table$target_site)) {
            site_coefs <- coef_table %>%
                filter(target_site == site, term != "(Intercept)", !grepl("^sample", term))
            if (nrow(site_coefs) == 0) next

            p <- site_coefs %>%
                mutate(term = fct_reorder(term, estimate), significant = p_adj < 0.05) %>%
                ggplot(aes(x = estimate, y = term, color = significant)) +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
                geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
                scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey50")) +
                labs(
                    x = "Log-odds (estimate)", y = NULL,
                    title = sprintf("GLMM: predictors of m6A at %s", site), color = "padj < 0.05"
                ) +
                theme_minimal(base_size = 12)
            mysaveandstore(str_glue("{outdir_glmm}/forest_{site}.pdf"), w = 8, h = max(4, nrow(site_coefs) * 0.35 + 1))
        }

        # Plot 3: Co-modification network
        comod_edges <- coef_table %>%
            filter(grepl("^site_", term), p_adj < 0.05) %>%
            mutate(from = gsub("site_", "", term), to = gsub("site_", "", target_site)) %>%
            dplyr::select(from, to, estimate, p_adj) %>%
            filter(from != to)

        if (nrow(comod_edges) > 0) {
            nodes <- unique(c(comod_edges$from, comod_edges$to))
            node_pos <- tibble(
                node = nodes,
                angle = seq(0, 2 * pi, length.out = length(nodes) + 1)[1:length(nodes)],
                x = cos(angle), y = sin(angle)
            )
            edge_df <- comod_edges %>%
                left_join(node_pos %>% dplyr::rename(from = node, x_from = x, y_from = y, angle_from = angle), by = "from") %>%
                left_join(node_pos %>% dplyr::rename(to = node, x_to = x, y_to = y, angle_to = angle), by = "to")

            p <- ggplot() +
                geom_segment(data = edge_df, aes(
                    x = x_from, y = y_from, xend = x_to, yend = y_to,
                    color = estimate, linewidth = abs(estimate)
                ), alpha = 0.6) +
                geom_point(data = node_pos, aes(x = x, y = y), size = 5, color = "black") +
                geom_text(data = node_pos, aes(x = x * 1.15, y = y * 1.15, label = node), size = 3) +
                scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0, name = "Effect\n(log-odds)") +
                scale_linewidth_continuous(range = c(0.5, 3), guide = "none") +
                coord_fixed() +
                labs(title = "Co-modification network (significant associations)") +
                theme_void() +
                theme(legend.position = "right")
            mysaveandstore(str_glue("{outdir_glmm}/comodification_network.pdf"), w = 8, h = 8)
        }

        # Plot 4: Random effect variance
        if (nrow(ranef_table) > 0) {
            p <- ranef_table %>%
                filter(effect == "ran_pars") %>%
                ggplot(aes(x = target_site, y = estimate, fill = group)) +
                geom_col(position = "dodge") +
                labs(
                    x = "Target site", y = "Variance (SD)", fill = "Random effect",
                    title = "Random effect variance components across models"
                ) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            mysaveandstore(str_glue("{outdir_glmm}/random_effect_variance.pdf"), w = max(6, n_success * 0.8 + 2), h = 5)
        }

        # Plot 5: Model fit summary
        model_summary <- model_results %>%
            keep(~ !is.null(.x)) %>%
            imap_dfr(~ tibble(target_site = .y, AIC = AIC(.x), n_obs = nobs(.x), converged = .x$sdr$pdHess))
        write_tsv(model_summary, str_glue("{intdir}/glmm_model_summary.tsv"))

        p <- model_summary %>%
            mutate(target_site = fct_reorder(target_site, AIC)) %>%
            ggplot(aes(x = target_site, y = AIC, fill = converged)) +
            geom_col() +
            scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "red")) +
            labs(x = NULL, y = "AIC", fill = "Converged", title = "Model fit (AIC) for each target site") +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outdir_glmm}/model_fit_summary.pdf"), w = max(6, n_success * 0.8 + 2), h = 5)

        # Plot 6: Predicted probability curves
        top_predictors <- coef_table %>%
            filter(term != "(Intercept)", !grepl("^sample", term), p_adj < 0.05) %>%
            group_by(term) %>%
            summarise(mean_abs_effect = mean(abs(estimate)), .groups = "drop") %>%
            arrange(desc(mean_abs_effect)) %>%
            head(5)

        if (nrow(top_predictors) > 0) {
            tryCatch(
                {
                    library(ggeffects)
                    for (pred in top_predictors$term) {
                        target <- coef_table %>%
                            filter(term == pred, p_adj < 0.05) %>%
                            arrange(p_adj) %>%
                            head(1) %>%
                            pull(target_site)
                        if (length(target) == 0) next
                        mod <- model_results[[target]]
                        if (is.null(mod)) next
                        pred_effects <- tryCatch(ggpredict(mod, terms = pred), error = function(e) NULL)
                        if (is.null(pred_effects)) next
                        p <- plot(pred_effects) +
                            labs(title = sprintf("Predicted P(m6A) at %s by %s", target, pred)) + theme_minimal(base_size = 12)
                        mysaveandstore(str_glue("{outdir_glmm}/predicted_prob_{target}_{pred}.pdf"), w = 7, h = 5)
                    }
                },
                error = function(e) cat(sprintf("ggeffects plotting failed: %s\n", e$message))
            )
        }
    }

    # --- CONTEXT-ONLY GLMM PLOTS ---
    if (nrow(coef_table_nosite) > 0 && n_success_nosite > 0) {
        coef_table_nosite <- coef_table_nosite %>%
            mutate(predictor_category = case_when(
                term %in% c("l2fc_scaled", "baseMean_scaled", "orf1_expression_scaled") ~ "Expression",
                term %in% c("intactness_req", "genic_loc", "loc_superlowres_integrative_stranded") ~ "Genomic context",
                term == "(Intercept)" ~ "Intercept",
                grepl("^sample", term) ~ "Sample",
                TRUE ~ "Other"
            ))

        coef_nosite_hm <- coef_table_nosite %>%
            filter(term != "(Intercept)", !grepl("^sample", term)) %>%
            dplyr::select(target_site, term, estimate) %>%
            pivot_wider(names_from = target_site, values_from = estimate) %>%
            column_to_rownames("term") %>%
            as.matrix()

        sig_nosite_hm <- coef_table_nosite %>%
            filter(term != "(Intercept)", !grepl("^sample", term)) %>%
            dplyr::select(target_site, term, p_adj) %>%
            pivot_wider(names_from = target_site, values_from = p_adj) %>%
            column_to_rownames("term") %>%
            as.matrix()

        star_mat_nosite <- matrix("", nrow = nrow(sig_nosite_hm), ncol = ncol(sig_nosite_hm))
        star_mat_nosite[sig_nosite_hm < 0.05] <- "*"
        star_mat_nosite[sig_nosite_hm < 0.01] <- "**"
        star_mat_nosite[sig_nosite_hm < 0.001] <- "***"

        if (nrow(coef_nosite_hm) > 1 && ncol(coef_nosite_hm) > 1) {
            row_cats_ns <- coef_table_nosite %>%
                filter(term != "(Intercept)", !grepl("^sample", term)) %>%
                dplyr::select(term, predictor_category) %>%
                distinct() %>%
                column_to_rownames("term")
            row_cats_ns <- row_cats_ns[rownames(coef_nosite_hm), , drop = FALSE]
            max_abs_ns <- min(max(abs(coef_nosite_hm), na.rm = TRUE), 5)

            p <- Heatmap(coef_nosite_hm,
                name = "log-odds",
                col = circlize::colorRamp2(c(-max_abs_ns, 0, max_abs_ns), c("blue", "white", "red")),
                na_col = "grey90", cluster_rows = TRUE, cluster_columns = TRUE,
                left_annotation = rowAnnotation(Category = row_cats_ns$predictor_category),
                cell_fun = function(j, i, x, y, width, height, fill) {
                    if (!is.na(star_mat_nosite[i, j]) && star_mat_nosite[i, j] != "") {
                        grid::grid.text(star_mat_nosite[i, j], x, y, gp = grid::gpar(fontsize = 8))
                    }
                },
                column_title = "Context-only GLMM coefficients (no co-modification predictors)"
            )
            mysaveandstore(str_glue("{outdir_glmm}/coefficient_heatmap_nosite.pdf"),
                w = max(8, ncol(coef_nosite_hm) * 0.8 + 4), h = max(6, nrow(coef_nosite_hm) * 0.4 + 2)
            )
        }

        # Forest plots (context-only)
        for (site in unique(coef_table_nosite$target_site)) {
            site_coefs <- coef_table_nosite %>%
                filter(target_site == site, term != "(Intercept)", !grepl("^sample", term))
            if (nrow(site_coefs) == 0) next
            p <- site_coefs %>%
                mutate(term = fct_reorder(term, estimate), significant = p_adj < 0.05) %>%
                ggplot(aes(x = estimate, y = term, color = significant)) +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
                geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
                scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey50")) +
                labs(
                    x = "Log-odds (estimate)", y = NULL,
                    title = sprintf("Context-only GLMM: predictors of m6A at %s", site), color = "padj < 0.05"
                ) +
                theme_minimal(base_size = 12)
            mysaveandstore(str_glue("{outdir_glmm}/forest_nosite_{site}.pdf"), w = 8, h = max(4, nrow(site_coefs) * 0.35 + 1))
        }

        # AIC comparison
        model_summary_nosite <- model_results_nosite %>%
            keep(~ !is.null(.x)) %>%
            imap_dfr(~ tibble(target_site = .y, AIC = AIC(.x), n_obs = nobs(.x), converged = .x$sdr$pdHess))
        write_tsv(model_summary_nosite, str_glue("{intdir}/glmm_model_summary_nosite.tsv"))

        if (exists("model_summary") && nrow(model_summary) > 0) {
            aic_compare <- model_summary %>%
                dplyr::select(target_site, AIC_full = AIC) %>%
                inner_join(model_summary_nosite %>% dplyr::select(target_site, AIC_nosite = AIC), by = "target_site") %>%
                mutate(delta_AIC = AIC_nosite - AIC_full)
            write_tsv(aic_compare, str_glue("{intdir}/glmm_aic_comparison.tsv"))

            p <- aic_compare %>%
                mutate(target_site = fct_reorder(target_site, delta_AIC)) %>%
                ggplot(aes(x = target_site, y = delta_AIC)) +
                geom_col(aes(fill = delta_AIC > 0)) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                scale_fill_manual(
                    values = c("TRUE" = "#D62728", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "Co-mod improves fit", "FALSE" = "Context-only better")
                ) +
                labs(
                    x = NULL, y = expression(Delta * "AIC (context-only - full)"), fill = NULL,
                    title = "Model comparison: does co-modification information improve fit?"
                ) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            mysaveandstore(str_glue("{outdir_glmm}/aic_full_vs_nosite.pdf"), w = max(6, nrow(aic_compare) * 0.8 + 2), h = 5)
        }

        # Per-predictor context plots
        context_predictors <- c("l2fc_scaled", "baseMean_scaled", "intactness_req", "genic_loc")
        context_predictors <- context_predictors[context_predictors %in% coef_table_nosite$term]

        coef_nosite_ctx <- coef_table_nosite %>%
            filter(term %in% context_predictors) %>%
            mutate(
                pos = as.integer(gsub("site_", "", target_site)),
                target_site = fct_reorder(target_site, pos)
            )

        plot_ctx_predictor_panels <- function(ctx_data, suffix = "", title_extra = "") {
            for (pred in context_predictors) {
                pred_data <- ctx_data %>% filter(term == pred)
                if (nrow(pred_data) == 0) next
                pred_label <- case_when(
                    pred == "l2fc_scaled" ~ "ORF1 enrichment (log2FC)",
                    pred == "baseMean_scaled" ~ "Expression level (baseMean)",
                    pred == "intactness_req" ~ "Element intactness",
                    pred == "genic_loc" ~ "Genic location", TRUE ~ pred
                )
                pred_title <- paste0(pred_label, title_extra)

                p <- pred_data %>%
                    mutate(sig_label = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "")) %>%
                    ggplot(aes(x = target_site, y = estimate)) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                    geom_segment(aes(xend = target_site, y = 0, yend = estimate, color = estimate > 0), linewidth = 0.8) +
                    geom_point(aes(fill = estimate > 0), shape = 21, size = 3, color = "black") +
                    geom_text(aes(label = sig_label), vjust = -1, size = 4) +
                    scale_fill_manual(values = c("TRUE" = "#D62728", "FALSE" = "steelblue"), guide = "none") +
                    scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "steelblue"), guide = "none") +
                    labs(
                        x = "Site (linear position)", y = "Coefficient (log-odds)",
                        title = sprintf("%s: effect on m6A across sites", pred_title),
                        subtitle = "Stars = BH-adjusted p-value (* <0.05, ** <0.01, *** <0.001)"
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outdir_glmm}/ctx_lollipop_{pred}{suffix}.pdf"),
                    w = max(7, nrow(pred_data) * 0.6 + 2), h = 5
                )

                p <- pred_data %>%
                    mutate(padj_sig = p_adj < 0.05, neg_log10_p = -log10(p.value)) %>%
                    ggplot(aes(x = target_site, y = estimate)) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                    geom_pointrange(aes(ymin = conf.low, ymax = conf.high, color = neg_log10_p, shape = padj_sig), size = 0.8) +
                    scale_color_viridis_c(option = "inferno", name = expression(-log[10] * "(p)")) +
                    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "padj < 0.05") +
                    labs(
                        x = "Site (linear position)", y = "Coefficient (log-odds)",
                        title = sprintf("%s: effect on m6A across sites", pred_title)
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outdir_glmm}/ctx_pointrange_{pred}{suffix}.pdf"),
                    w = max(7, nrow(pred_data) * 0.6 + 2), h = 5
                )

                p <- pred_data %>%
                    mutate(padj_bin = cut(p_adj,
                        breaks = c(0, 0.001, 0.01, 0.05, 1),
                        labels = c("<0.001", "<0.01", "<0.05", "n.s."), include.lowest = TRUE
                    )) %>%
                    ggplot(aes(x = target_site, y = estimate, fill = padj_bin)) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
                    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3) +
                    scale_fill_manual(values = c("<0.001" = "#D62728", "<0.01" = "#FF7F0E", "<0.05" = "#FFBB78", "n.s." = "grey80"), name = "Adjusted p") +
                    labs(
                        x = "Site (linear position)", y = "Coefficient (log-odds)",
                        title = sprintf("%s: effect on m6A across sites", pred_title)
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outdir_glmm}/ctx_barplot_{pred}{suffix}.pdf"),
                    w = max(7, nrow(pred_data) * 0.6 + 2), h = 5
                )

                p <- pred_data %>%
                    ggplot(aes(x = target_site, y = 1)) +
                    geom_tile(aes(fill = estimate), color = "white", linewidth = 0.5) +
                    geom_text(aes(label = sprintf("%.2f", estimate)), size = 2.5) +
                    geom_point(
                        data = pred_data %>% filter(p_adj < 0.05),
                        aes(x = target_site, y = 1), shape = 8, size = 2, color = "black", nudge_y = 0.35
                    ) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "#D62728", midpoint = 0, name = "log-odds") +
                    labs(
                        x = "Site (linear position)", y = NULL,
                        title = sprintf("%s: effect on m6A%s", pred_label, title_extra), caption = "Star = padj < 0.05"
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(
                        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1)
                    )
                mysaveandstore(str_glue("{outdir_glmm}/ctx_tile_{pred}{suffix}.pdf"),
                    w = max(7, nrow(pred_data) * 0.6 + 2), h = 3
                )
            }

            if (nrow(ctx_data) > 0) {
                p <- ctx_data %>%
                    mutate(
                        pred_label = case_when(
                            term == "l2fc_scaled" ~ "ORF1 enrichment",
                            term == "baseMean_scaled" ~ "Expression", term == "intactness_req" ~ "Intactness",
                            term == "genic_loc" ~ "Genic location", TRUE ~ term
                        ),
                        sig_star = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "")
                    ) %>%
                    ggplot(aes(x = target_site, y = pred_label)) +
                    geom_tile(aes(fill = estimate), color = "white", linewidth = 0.5) +
                    geom_text(aes(label = sig_star), size = 3.5) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "#D62728", midpoint = 0, name = "log-odds") +
                    labs(
                        x = "Site (linear position)", y = NULL,
                        title = paste0("Context predictors of m6A across L1 sites", title_extra)
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                mysaveandstore(str_glue("{outdir_glmm}/ctx_combined_tile{suffix}.pdf"),
                    w = max(7, n_distinct(ctx_data$target_site) * 0.6 + 3), h = 4
                )
            }
        }

        plot_ctx_predictor_panels(coef_nosite_ctx, suffix = "", title_extra = "")
        coef_nosite_ctx_no5989 <- coef_nosite_ctx %>%
            filter(target_site != "site_5989") %>%
            mutate(target_site = fct_drop(target_site))
        plot_ctx_predictor_panels(coef_nosite_ctx_no5989, suffix = "_no5989", title_extra = " (excl. site 5989)")
    }

    list(
        coef_table = coef_table, coef_table_nosite = coef_table_nosite,
        model_results = model_results, model_results_nosite = model_results_nosite
    )
}

# --- Run GLMMs at global scope (parallelization works cleanly here) ---
# FL-only
glmm_all_yl1 <- run_glmm(results_all_yl1$modeling_inputs, suffix = "")
glmm_l1hs <- run_glmm(results_l1hs$modeling_inputs, suffix = "_L1HS")
# All-length
glmm_all_yl1_alllength <- run_glmm(results_all_yl1_alllength$modeling_inputs, suffix = "_alllength")
glmm_l1hs_alllength <- run_glmm(results_l1hs_alllength$modeling_inputs, suffix = "_L1HS_alllength")

# Use young L1 FL results as defaults for Section 4
reliable_sites <- results_all_yl1$reliable_sites
site_classification <- results_all_yl1$site_classification


# ============================================================
# //ANCHOR - SECTION 4: EVOLUTIONARY ANALYSIS
# ============================================================
cat("\n=== Section 4: Evolutionary Analysis ===\n")

DRACH_REGEX <- "^[AGT][AG]AC[ACT]$"
DRACH_REGEX_MATCH <- "[AGT][AG]AC[ACT]"

# 4a. Map reliable sites across subfamilies
reliable_harmonized <- reliable_sites$harmonized_pos

site_across_subfamilies <- tibble(harmonized_pos = reliable_harmonized) %>%
    left_join(pos_to_aln_df, by = "harmonized_pos", relationship = "many-to-many") %>%
    filter(subfamily_consensus %in% SUBFAM_ORDER) %>%
    dplyr::rename(subfamily = subfamily_consensus)

# 4b. Extract 5-mer context from each subfamily consensus
site_across_subfamilies <- site_across_subfamilies %>%
    rowwise() %>%
    mutate(
        cs_str = as.character(consensus_seqs[[subfamily]][[1]]),
        cs_len = nchar(cs_str),
        context_5mer = ifelse(
            consensus_pos >= 3 & consensus_pos <= cs_len - 2,
            substr(cs_str, consensus_pos - 2, consensus_pos + 2),
            NA_character_
        ),
        central_base = substr(cs_str, consensus_pos, consensus_pos),
        is_A = central_base == "A",
        is_DRACH = !is.na(context_5mer) & is_A & grepl(DRACH_REGEX, context_5mer)
    ) %>%
    ungroup() %>%
    dplyr::select(-cs_str, -cs_len)

write_tsv(site_across_subfamilies, str_glue("{intdir}/site_context_across_subfamilies.tsv"))

# 4c. Categorize evolutionary trajectory
site_evolution <- site_across_subfamilies %>%
    mutate(subfamily = factor(subfamily, levels = SUBFAM_ORDER)) %>%
    arrange(harmonized_pos, subfamily) %>%
    group_by(harmonized_pos) %>%
    summarise(
        drach_pattern = paste(ifelse(is_DRACH, "Y", "N"), collapse = "->"),
        n_subfams_drach = sum(is_DRACH, na.rm = TRUE),
        n_subfams_total = n(),
        youngest_drach = ifelse(any(is_DRACH), as.character(dplyr::first(subfamily[is_DRACH])), NA_character_),
        oldest_drach = ifelse(any(is_DRACH), as.character(dplyr::last(subfamily[is_DRACH])), NA_character_),
        l1hs_drach = is_DRACH[subfamily == "L1HS"],
        .groups = "drop"
    ) %>%
    mutate(
        category = case_when(
            n_subfams_drach == n_subfams_total ~ "Ancestrally conserved",
            l1hs_drach & n_subfams_drach < n_subfams_total ~ "Recently gained",
            !l1hs_drach & n_subfams_drach > 0 ~ "Lost in young",
            n_subfams_drach == 0 ~ "No DRACH",
            TRUE ~ "Complex"
        )
    )

write_tsv(site_evolution, str_glue("{intdir}/site_evolution_categories.tsv"))

cat("Evolutionary categories:\n")
print(table(site_evolution$category))

# 4d. Within-subfamily DRACH polymorphism
# Scan each FL element's actual genomic sequence for DRACH motifs using IUPAC-aware
# pattern matching, then map back to harmonized positions.
# This approach (from the validated bedmethylanalysis_process_new.R) is more robust
# than extracting 5-mers at specific positions.

genome_fa_file <- FaFile(conf$reference)

# Get FL elements for each subfamily
fl_elements <- rmann %>%
    filter(rte_subfamily %in% SUBFAM_ORDER, rte_length_req == "FL") %>%
    dplyr::select(gene_id, seqnames, start, end, strand, rte_subfamily) %>%
    distinct()

# Subfamily-specific consensus-to-harmonized mapping
subfam_mappings <- lapply(SUBFAM_ORDER, function(sf) {
    pos_to_aln_df %>%
        filter(subfamily_consensus == sf) %>%
        dplyr::select(consensus_pos, harmonized_pos)
})
names(subfam_mappings) <- SUBFAM_ORDER

# Process each subfamily
drach_all_subfamilies <- list()

for (sf in SUBFAM_ORDER) {
    sf_elements <- fl_elements %>% filter(rte_subfamily == sf)
    if (nrow(sf_elements) == 0) next
    cat(sprintf("Scanning %d FL %s elements for DRACH...\n", nrow(sf_elements), sf))

    # Extract full element sequences in one batch
    el_gr <- GRanges(
        seqnames = sf_elements$seqnames,
        ranges = IRanges(start = sf_elements$start, end = sf_elements$end)
    )
    el_seqs <- getSeq(genome_fa_file, el_gr)
    # Reverse complement minus-strand elements so sequence is in consensus orientation
    minus_idx <- as.character(sf_elements$strand) == "-"
    el_seqs[minus_idx] <- reverseComplement(el_seqs[minus_idx])
    names(el_seqs) <- sf_elements$gene_id

    # Scan each element for DRACH motifs using IUPAC-aware matching
    sf_mapping <- subfam_mappings[[sf]]

    drach_per_element <- list()
    for (i in seq_along(el_seqs)) {
        gid <- names(el_seqs)[i]
        seq_i <- el_seqs[[i]]
        dm <- matchPattern("DRACH", seq_i, fixed = FALSE)
        if (length(dm) == 0) next
        a_positions <- start(dm) + 2L # A is 3rd base in DRACH
        motifs <- as.character(dm)

        drach_per_element <- c(drach_per_element, list(tibble(
            gene_id = gid,
            element_pos = a_positions,
            motif = motifs,
            is_drach = TRUE
        )))
    }
    if (length(drach_per_element) == 0) next
    drach_element_df <- bind_rows(drach_per_element)

    # Map element positions → consensus positions via consensus_index_long
    drach_element_mapped <- drach_element_df %>%
        left_join(
            consensus_index_long %>%
                filter(gene_id %in% sf_elements$gene_id) %>%
                dplyr::select(gene_id, sequence_pos, consensus_pos),
            by = c("gene_id", "element_pos" = "sequence_pos")
        ) %>%
        filter(!is.na(consensus_pos)) %>%
        left_join(sf_mapping, by = "consensus_pos") %>%
        filter(!is.na(harmonized_pos)) %>%
        mutate(rte_subfamily = sf)

    # Build complete element × harmonized_pos matrix (DRACH or not at each position)
    # All FL elements × all harmonized positions they have data for
    sf_element_hpos <- consensus_index_long %>%
        filter(gene_id %in% sf_elements$gene_id) %>%
        dplyr::select(gene_id, sequence_pos, consensus_pos) %>%
        left_join(sf_mapping, by = "consensus_pos") %>%
        filter(!is.na(harmonized_pos)) %>%
        dplyr::select(gene_id, harmonized_pos) %>%
        distinct()

    drach_poly_df_sf <- sf_element_hpos %>%
        left_join(
            drach_element_mapped %>%
                dplyr::select(gene_id, harmonized_pos, motif) %>%
                mutate(is_drach = TRUE) %>%
                distinct(gene_id, harmonized_pos, .keep_all = TRUE),
            by = c("gene_id", "harmonized_pos")
        ) %>%
        mutate(
            is_drach = replace_na(is_drach, FALSE),
            motif = replace_na(motif, ""),
            rte_subfamily = sf
        )

    drach_all_subfamilies <- c(drach_all_subfamilies, list(drach_poly_df_sf))
    cat(sprintf(
        "  %s: %d DRACH hits in %d elements, %d element×position pairs\n",
        sf, nrow(drach_element_df), n_distinct(drach_element_df$gene_id), nrow(drach_poly_df_sf)
    ))
}

drach_poly_df <- bind_rows(drach_all_subfamilies)

# Summarize polymorphism per harmonized position × subfamily
polymorphism_summary <- drach_poly_df %>%
    group_by(harmonized_pos, rte_subfamily) %>%
    summarise(
        n_elements = n_distinct(gene_id),
        n_drach = sum(is_drach),
        n_nondrach = sum(!is_drach),
        frac_drach = n_drach / n_elements,
        pct_drach = 100 * n_drach / n_elements,
        is_polymorphic = pmin(pct_drach, 100 - pct_drach) >= 10,
        .groups = "drop"
    )

write_tsv(polymorphism_summary, str_glue("{intdir}/element_drach_polymorphism.tsv"))

# Also save the full element × position DRACH status
write_tsv(drach_poly_df, str_glue("{intdir}/element_drach_status_full.tsv"))

# Summary of polymorphic positions
poly_summary_overall <- polymorphism_summary %>%
    group_by(rte_subfamily) %>%
    summarise(
        n_positions = n(),
        n_polymorphic = sum(is_polymorphic),
        n_conserved_drach = sum(pct_drach >= 90),
        n_conserved_nondrach = sum(pct_drach <= 10),
        .groups = "drop"
    )
cat("DRACH polymorphism summary by subfamily:\n")
print(poly_summary_overall)

# Polymorphic positions at reliable sites specifically
polymorphism_at_reliable <- polymorphism_summary %>%
    filter(harmonized_pos %in% reliable_harmonized)
cat(sprintf(
    "\nPolymorphic DRACH at reliable sites: %d / %d position×subfamily combinations\n",
    sum(polymorphism_at_reliable$is_polymorphic), nrow(polymorphism_at_reliable)
))

# 4e. Link polymorphism to modification
drach_mod_link <- drach_poly_df %>%
    dplyr::select(gene_id, harmonized_pos, is_drach) %>%
    distinct() %>%
    inner_join(
        yl1annm_filtered %>%
            filter(condition == "N_ORF1") %>%
            dplyr::select(gene_id, harmonized_pos, pctM, sample),
        by = c("gene_id", "harmonized_pos")
    )

# --- EVOLUTIONARY PLOTS ---

# Plot 0a: Global DRACH polymorphism histogram (all positions, L1HS)
l1hs_poly <- polymorphism_summary %>% filter(rte_subfamily == "L1HS")
p <- l1hs_poly %>%
    ggplot(aes(x = pct_drach)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    geom_vline(xintercept = c(10, 90), linetype = "dashed", color = "red") +
    labs(
        x = "% elements with DRACH motif", y = "# positions",
        title = sprintf(
            "L1HS DRACH polymorphism (%d positions, %d polymorphic)",
            nrow(l1hs_poly), sum(l1hs_poly$is_polymorphic)
        )
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_evo}/drach_polymorphism_histogram.pdf"), w = 7, h = 5)

# Plot 0b: Polymorphism along the element (all positions)
p <- l1hs_poly %>%
    left_join(pos_to_aln_df %>% filter(subfamily_consensus == "L1HS") %>%
        dplyr::select(harmonized_pos, consensus_pos), by = "harmonized_pos") %>%
    mutate(status = case_when(
        is_polymorphic ~ "Polymorphic",
        pct_drach >= 90 ~ "Conserved DRACH",
        TRUE ~ "Conserved non-DRACH"
    )) %>%
    filter(!is.na(consensus_pos)) %>%
    ggplot(aes(x = consensus_pos, y = pct_drach, color = status)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = c(10, 90), linetype = "dashed", alpha = 0.5) +
    geom_point(
        data = . %>% filter(harmonized_pos %in% reliable_harmonized),
        shape = 1, size = 4, color = "black", stroke = 1
    ) +
    annotate("rect",
        xmin = L1_REGIONS$start, xmax = L1_REGIONS$end,
        ymin = -8, ymax = -4, fill = L1_REGIONS$color, alpha = 0.5
    ) +
    annotate("text",
        x = (L1_REGIONS$start + L1_REGIONS$end) / 2,
        y = -6, label = L1_REGIONS$region, size = 2.5
    ) +
    labs(
        x = "L1HS consensus position", y = "% elements with DRACH",
        color = NULL, title = "DRACH polymorphism along L1HS (circles = reliable m6A sites)"
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_evo}/drach_polymorphism_along_element.pdf"), w = 14, h = 6)

# Plot 1: DRACH fraction evolution across subfamilies (reliable sites only)
reliable_poly <- polymorphism_summary %>%
    filter(harmonized_pos %in% reliable_harmonized)

p <- reliable_poly %>%
    left_join(site_evolution %>% dplyr::select(harmonized_pos, category), by = "harmonized_pos") %>%
    mutate(rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)) %>%
    ggplot(aes(x = rte_subfamily, y = frac_drach, group = factor(harmonized_pos), color = category)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        x = "Subfamily (young -> old)", y = "Fraction of elements with DRACH",
        color = "Evolutionary category",
        title = "DRACH motif evolution across L1 subfamilies at reliable m6A sites"
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_evo}/drach_fraction_evolution.pdf"), w = 10, h = 6)

# Same but faceted by site for clarity
n_reliable_sites <- n_distinct(reliable_poly$harmonized_pos)
p <- reliable_poly %>%
    left_join(site_evolution %>% dplyr::select(harmonized_pos, category), by = "harmonized_pos") %>%
    mutate(rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)) %>%
    ggplot(aes(x = rte_subfamily, y = frac_drach, group = 1)) +
    geom_line(color = "steelblue") +
    geom_point(aes(color = category), size = 2) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~harmonized_pos, scales = "free_y") +
    labs(x = "Subfamily", y = "Fraction DRACH", color = "Category") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
mysaveandstore(str_glue("{outdir_evo}/drach_fraction_evolution_faceted.pdf"),
    w = max(8, ceiling(sqrt(n_reliable_sites)) * 3),
    h = max(6, ceiling(n_reliable_sites / ceiling(sqrt(n_reliable_sites))) * 3)
)

# Plot 2: Alignment window tiles (reusing existing pattern)
flank <- 2
subfams <- names(aligned)
window_data_evo <- list()

for (hpos in reliable_harmonized) {
    col_start <- max(1, hpos - flank)
    col_end <- min(width(aligned)[1], hpos + flank)
    for (sf in subfams) {
        aln_seq <- as.character(subseq(aligned[[sf]], start = col_start, end = col_end))
        bases <- strsplit(aln_seq, "")[[1]]
        relative_pos <- seq(col_start - hpos, col_end - hpos)
        window_data_evo <- c(window_data_evo, list(tibble(
            harmonized_pos = hpos,
            subfamily = sf,
            rel_pos = relative_pos,
            base = bases
        )))
    }
}
window_df_evo <- bind_rows(window_data_evo)

drach_check_evo <- window_df_evo %>%
    arrange(harmonized_pos, subfamily, rel_pos) %>%
    group_by(harmonized_pos, subfamily) %>%
    summarise(motif = paste0(base, collapse = ""), .groups = "drop") %>%
    mutate(is_drach = grepl(DRACH_REGEX, motif))

base_colors <- c("A" = "#4DAF4A", "T" = "#E41A1C", "C" = "#377EB8", "G" = "#FF7F00", "-" = "grey90")

# Get mean modification at each site x subfamily
mod_by_site_subfam <- yl1annm_filtered %>%
    filter(condition == "N_ORF1") %>%
    filter(harmonized_pos %in% reliable_harmonized) %>%
    group_by(harmonized_pos, rte_subfamily) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop")

plot_df_evo <- window_df_evo %>%
    mutate(
        subfamily = factor(subfamily, levels = rev(subfams)),
        harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))
    ) %>%
    left_join(
        drach_check_evo %>% mutate(
            subfamily = factor(subfamily, levels = rev(subfams)),
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))
        ),
        by = c("harmonized_pos", "subfamily")
    )

star_df_evo <- plot_df_evo %>%
    filter(!is_drach) %>%
    group_by(harmonized_pos, subfamily) %>%
    dplyr::slice(1)

meth_annot_evo <- mod_by_site_subfam %>%
    mutate(
        subfamily = factor(rte_subfamily, levels = rev(subfams)),
        harmonized_pos = factor(harmonized_pos, levels = sort(unique(as.numeric(as.character(plot_df_evo$harmonized_pos))))),
        mm_label = sprintf("%.0f%%", mm)
    ) %>%
    filter(!is.na(harmonized_pos))

p <- plot_df_evo %>%
    ggplot(aes(x = rel_pos, y = subfamily)) +
    geom_tile(aes(fill = base), color = "white", linewidth = 0.3) +
    geom_text(aes(label = base), size = 3) +
    geom_text(
        data = star_df_evo, aes(x = flank + 0.8, y = subfamily, label = "*"),
        size = 9, color = "black", inherit.aes = FALSE
    ) +
    geom_text(
        data = meth_annot_evo, aes(x = flank + 1.8, y = subfamily, label = mm_label, color = mm),
        size = 5, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_color_gradient2(low = "blue", mid = "grey40", high = "red", midpoint = 25, guide = "none") +
    scale_fill_manual(values = base_colors) +
    facet_wrap(~harmonized_pos, ncol = 4, labeller = labeller(harmonized_pos = function(x) paste0("pos ", x))) +
    scale_x_continuous(breaks = seq(-flank, flank)) +
    coord_cartesian(xlim = c(-flank - 0.5, flank + 2.5), clip = "off") +
    labs(
        x = "Position relative to site", y = NULL,
        title = "Sequence context at reliable m6A sites across L1 subfamilies"
    ) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(size = 8))
mysaveandstore(str_glue("{outdir_evo}/alignment_windows_reliable_sites.pdf"),
    w = 4 * min(length(reliable_harmonized), 4) + 2,
    h = 2 + length(subfams) * 0.4 * ceiling(length(reliable_harmonized) / 4)
)

# Plot 3: Evolutionary category summary
p <- site_evolution %>%
    count(category) %>%
    mutate(category = fct_reorder(category, n)) %>%
    ggplot(aes(x = category, y = n, fill = category)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    labs(
        x = NULL, y = "# reliable m6A sites",
        title = "Evolutionary classification of DRACH context at reliable m6A sites"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
mysaveandstore(str_glue("{outdir_evo}/evolutionary_category_summary.pdf"), w = 7, h = 5)

# Plot 4: Within-subfamily polymorphism violin (DRACH vs non-DRACH elements -> modification)
if (exists("drach_mod_link") && nrow(drach_mod_link) > 0) {
    # Identify positions with polymorphism (10-90% DRACH)
    poly_positions <- polymorphism_summary %>%
        filter(frac_drach >= 0.1, frac_drach <= 0.9, n_elements >= 10) %>%
        pull(harmonized_pos) %>%
        unique()

    if (length(poly_positions) > 0) {
        p <- drach_mod_link %>%
            filter(harmonized_pos %in% poly_positions) %>%
            mutate(drach_status = ifelse(is_drach, "DRACH", "non-DRACH")) %>%
            group_by(gene_id, harmonized_pos, drach_status) %>%
            summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            ggplot(aes(x = drach_status, y = mean_pctM, fill = drach_status)) +
            geom_violin(alpha = 0.3, scale = "width") +
            geom_boxplot(width = 0.2, outlier.shape = NA) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "Mean % m6A", fill = NULL,
                title = "m6A modification at polymorphic DRACH sites"
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none")
        mysaveandstore(str_glue("{outdir_evo}/drach_polymorphism_vs_modification.pdf"),
            w = 3 * min(length(poly_positions), 4) + 1,
            h = 3 * ceiling(length(poly_positions) / 4) + 1
        )
    }
}

# Plot 5: Modification vs DRACH conservation
mod_at_sites <- yl1annm_filtered %>%
    filter(condition == "N_ORF1", harmonized_pos %in% reliable_harmonized) %>%
    group_by(harmonized_pos) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop")

p <- site_evolution %>%
    left_join(mod_at_sites, by = "harmonized_pos") %>%
    filter(!is.na(mm)) %>%
    ggplot(aes(x = n_subfams_drach, y = mm, color = category)) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_text_repel(aes(label = harmonized_pos), size = 2.5, max.overlaps = 15) +
    labs(
        x = "# subfamilies with DRACH context", y = "Mean % m6A (ORF1 RIP)",
        color = "Category",
        title = "Does DRACH conservation predict modification level?"
    ) +
    theme_minimal(base_size = 12)
mysaveandstore(str_glue("{outdir_evo}/modification_vs_drach_conservation.pdf"), w = 8, h = 6)

# Plot 6: Consensus evolution heatmap (5-mer across subfamilies)
consensus_5mer_hm <- site_across_subfamilies %>%
    filter(!is.na(context_5mer)) %>%
    mutate(
        subfamily = factor(subfamily, levels = SUBFAM_ORDER),
        label = ifelse(is_DRACH, context_5mer, paste0(context_5mer, "*"))
    ) %>%
    dplyr::select(harmonized_pos, subfamily, context_5mer, is_DRACH)

p <- consensus_5mer_hm %>%
    mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
    ggplot(aes(x = harmonized_pos, y = subfamily, fill = is_DRACH)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = context_5mer), size = 2.5) +
    scale_fill_manual(values = c("TRUE" = "#90EE90", "FALSE" = "#FFB6C1"), labels = c("non-DRACH", "DRACH")) +
    labs(
        x = "Harmonized position", y = NULL, fill = "Motif context",
        title = "DRACH motif evolution at reliable m6A sites"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
mysaveandstore(str_glue("{outdir_evo}/consensus_5mer_evolution_heatmap.pdf"),
    w = max(8, n_distinct(consensus_5mer_hm$harmonized_pos) * 0.8 + 3), h = 5
)

# Plot 7: Motif variant frequency at polymorphic sites
if (nrow(drach_poly_df) > 0 && length(poly_positions) > 0) {
    # Get mean modification per motif variant using drach_poly_df
    motif_mod <- drach_poly_df %>%
        filter(harmonized_pos %in% poly_positions, motif != "") %>%
        left_join(
            yl1annm_filtered %>%
                filter(condition == "N_ORF1") %>%
                group_by(gene_id, harmonized_pos) %>%
                summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop"),
            by = c("gene_id", "harmonized_pos")
        ) %>%
        group_by(harmonized_pos, motif, is_drach) %>%
        summarise(mm = mean(mean_pctM, na.rm = TRUE), n = n(), .groups = "drop") %>%
        filter(!is.na(mm)) %>%
        group_by(harmonized_pos) %>%
        mutate(pct = 100 * n / sum(n)) %>%
        arrange(harmonized_pos, -n) %>%
        slice_head(n = 5) %>%
        ungroup()

    if (nrow(motif_mod) > 0) {
        p <- motif_mod %>%
            mutate(
                drach_label = ifelse(is_drach, "DRACH", "*"),
                motif_label = sprintf("%s (%d%%, %.0f%%mod)", motif, round(pct), mm),
                harmonized_pos = factor(harmonized_pos)
            ) %>%
            ggplot(aes(x = 1, y = reorder(motif_label, mm), fill = mm)) +
            geom_tile(color = "white") +
            geom_text(aes(label = drach_label), size = 3) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 25, name = "% mod") +
            facet_wrap(~harmonized_pos, scales = "free_y", ncol = 3) +
            labs(x = NULL, y = NULL, title = "Motif variants at polymorphic DRACH m6A sites") +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        mysaveandstore(str_glue("{outdir_evo}/motif_variant_frequency.pdf"),
            w = 4 * min(length(poly_positions), 3) + 1, h = 3 * ceiling(length(poly_positions) / 3) + 1
        )
    }
}

# Plot 8: DRACH gain/loss timeline (stylized cladogram)
# Simple horizontal layout: subfamilies as branches, colored nodes for DRACH status
timeline_data <- site_across_subfamilies %>%
    dplyr::select(harmonized_pos, subfamily, is_DRACH) %>%
    mutate(
        subfam_x = match(subfamily, SUBFAM_ORDER),
        harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))
    )

p <- timeline_data %>%
    ggplot(aes(x = subfam_x, y = harmonized_pos, fill = is_DRACH)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is_DRACH, "D", "")), size = 3) +
    scale_x_continuous(breaks = seq_along(SUBFAM_ORDER), labels = SUBFAM_ORDER) +
    scale_fill_manual(
        values = c("TRUE" = "#2ECC71", "FALSE" = "#E74C3C"),
        labels = c("non-DRACH", "DRACH")
    ) +
    labs(
        x = "Subfamily (young -> old)", y = "Harmonized position", fill = "DRACH context",
        title = "DRACH context evolution across L1 phylogeny at reliable m6A sites"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(str_glue("{outdir_evo}/drach_evolution_timeline.pdf"),
    w = 8, h = max(4, n_distinct(timeline_data$harmonized_pos) * 0.3 + 2)
)


# ============================================================
# //ANCHOR - SECTION 4b: EVOLUTIONARY ANALYSIS (ALL-LENGTH)
# ============================================================
cat("\n=== Section 4b: Evolutionary Analysis (FL + trnc) ===\n")
outdir_evo_al <- str_glue("{outdir_evo}_alllength")
intdir_evo_al <- str_glue("{intdir}_alllength")
for (d in c(outdir_evo_al, intdir_evo_al)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Scan ALL elements (FL + trnc) for DRACH motifs using the same approach as Section 4
all_elements <- rmann %>%
    filter(rte_subfamily %in% SUBFAM_ORDER) %>%
    dplyr::select(gene_id, seqnames, start, end, strand, rte_subfamily, rte_length_req) %>%
    distinct()

drach_all_subfamilies_al <- list()
for (sf in SUBFAM_ORDER) {
    sf_elements <- all_elements %>% filter(rte_subfamily == sf)
    if (nrow(sf_elements) == 0) next
    cat(sprintf("Scanning %d %s elements (FL+trnc) for DRACH...\n", nrow(sf_elements), sf))

    el_gr <- GRanges(
        seqnames = sf_elements$seqnames,
        ranges = IRanges(start = sf_elements$start, end = sf_elements$end)
    )
    el_seqs <- getSeq(genome_fa_file, el_gr)
    minus_idx <- as.character(sf_elements$strand) == "-"
    el_seqs[minus_idx] <- reverseComplement(el_seqs[minus_idx])
    names(el_seqs) <- sf_elements$gene_id

    sf_mapping <- subfam_mappings[[sf]]

    drach_per_element <- list()
    for (i in seq_along(el_seqs)) {
        gid <- names(el_seqs)[i]
        seq_i <- el_seqs[[i]]
        dm <- matchPattern("DRACH", seq_i, fixed = FALSE)
        if (length(dm) == 0) next
        a_positions <- start(dm) + 2L
        motifs <- as.character(dm)
        drach_per_element <- c(drach_per_element, list(tibble(
            gene_id = gid, element_pos = a_positions, motif = motifs, is_drach = TRUE
        )))
    }
    if (length(drach_per_element) == 0) next
    drach_element_df <- bind_rows(drach_per_element)

    drach_element_mapped <- drach_element_df %>%
        left_join(
            consensus_index_long %>%
                filter(gene_id %in% sf_elements$gene_id) %>%
                dplyr::select(gene_id, sequence_pos, consensus_pos),
            by = c("gene_id", "element_pos" = "sequence_pos")
        ) %>%
        filter(!is.na(consensus_pos)) %>%
        left_join(sf_mapping, by = "consensus_pos") %>%
        filter(!is.na(harmonized_pos)) %>%
        mutate(rte_subfamily = sf)

    sf_element_hpos <- consensus_index_long %>%
        filter(gene_id %in% sf_elements$gene_id) %>%
        dplyr::select(gene_id, sequence_pos, consensus_pos) %>%
        left_join(sf_mapping, by = "consensus_pos") %>%
        filter(!is.na(harmonized_pos)) %>%
        dplyr::select(gene_id, harmonized_pos) %>%
        distinct()

    drach_poly_df_sf <- sf_element_hpos %>%
        left_join(
            drach_element_mapped %>%
                dplyr::select(gene_id, harmonized_pos, motif) %>%
                mutate(is_drach = TRUE) %>%
                distinct(gene_id, harmonized_pos, .keep_all = TRUE),
            by = c("gene_id", "harmonized_pos")
        ) %>%
        mutate(
            is_drach = replace_na(is_drach, FALSE),
            motif = replace_na(motif, ""),
            rte_subfamily = sf
        )

    drach_all_subfamilies_al <- c(drach_all_subfamilies_al, list(drach_poly_df_sf))
}

drach_poly_df_al <- bind_rows(drach_all_subfamilies_al) %>%
    left_join(all_elements %>% dplyr::select(gene_id, rte_length_req), by = "gene_id")

# Polymorphism summary by length
polymorphism_summary_all <- drach_poly_df_al %>%
    group_by(harmonized_pos, rte_subfamily, rte_length_req) %>%
    summarise(
        n_elements = n_distinct(gene_id),
        n_drach = sum(is_drach),
        frac_drach = n_drach / n_elements,
        pct_drach = 100 * n_drach / n_elements,
        .groups = "drop"
    )
write_tsv(polymorphism_summary_all, str_glue("{intdir_evo_al}/element_drach_polymorphism_alllength.tsv"))

# Combined (not split by length)
polymorphism_summary_all_combined <- drach_poly_df_al %>%
    group_by(harmonized_pos, rte_subfamily) %>%
    summarise(
        n_elements = n_distinct(gene_id),
        n_drach = sum(is_drach),
        frac_drach = n_drach / n_elements,
        .groups = "drop"
    )

# Link DRACH to modification using all-length data
drach_mod_link_all <- drach_poly_df_al %>%
    dplyr::select(gene_id, harmonized_pos, is_drach, rte_length_req) %>%
    distinct() %>%
    inner_join(
        yl1annm_alllength_filtered %>%
            filter(condition == "N_ORF1") %>%
            dplyr::select(gene_id, harmonized_pos, pctM, sample),
        by = c("gene_id", "harmonized_pos")
    )

if (nrow(drach_poly_df_al) > 0) {
    # --- ALL-LENGTH EVOLUTIONARY PLOTS ---

    # Plot 1: DRACH fraction evolution, faceted by length
    p <- polymorphism_summary_all %>%
        left_join(site_evolution %>% dplyr::select(harmonized_pos, category), by = "harmonized_pos") %>%
        mutate(rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)) %>%
        ggplot(aes(x = rte_subfamily, y = frac_drach, group = factor(harmonized_pos), color = category)) +
        geom_line(alpha = 0.6) +
        geom_point(size = 2) +
        scale_y_continuous(labels = scales::percent) +
        facet_wrap(~rte_length_req) +
        labs(
            x = "Subfamily (young -> old)", y = "Fraction of elements with DRACH",
            color = "Evolutionary category",
            title = "DRACH motif evolution: FL vs truncated elements"
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_evo_al}/drach_fraction_evolution_by_length.pdf"), w = 14, h = 6)

    # Plot 2: FL vs trnc DRACH fraction comparison
    drach_fl_vs_trnc <- polymorphism_summary_all %>%
        dplyr::select(harmonized_pos, rte_subfamily, rte_length_req, frac_drach) %>%
        pivot_wider(names_from = rte_length_req, values_from = frac_drach, names_prefix = "frac_")

    if ("frac_FL" %in% names(drach_fl_vs_trnc) && "frac_trnc" %in% names(drach_fl_vs_trnc)) {
        p <- drach_fl_vs_trnc %>%
            filter(!is.na(frac_FL), !is.na(frac_trnc)) %>%
            ggplot(aes(x = frac_FL, y = frac_trnc, color = rte_subfamily)) +
            geom_point(size = 2, alpha = 0.7) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            scale_x_continuous(labels = scales::percent) +
            scale_y_continuous(labels = scales::percent) +
            labs(
                x = "DRACH fraction (FL elements)", y = "DRACH fraction (trnc elements)",
                color = "Subfamily",
                title = "DRACH polymorphism: FL vs truncated elements"
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_evo_al}/drach_fraction_fl_vs_trnc.pdf"), w = 8, h = 7)
    }

    # Plot 3: DRACH vs modification, faceted by length
    poly_positions_all <- polymorphism_summary_all %>%
        filter(frac_drach >= 0.1, frac_drach <= 0.9, n_elements >= 10) %>%
        pull(harmonized_pos) %>%
        unique()

    if (length(poly_positions_all) > 0 && nrow(drach_mod_link_all) > 0) {
        p <- drach_mod_link_all %>%
            filter(harmonized_pos %in% poly_positions_all) %>%
            mutate(drach_status = ifelse(is_drach, "DRACH", "non-DRACH")) %>%
            group_by(gene_id, harmonized_pos, drach_status, rte_length_req) %>%
            summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            ggplot(aes(x = drach_status, y = mean_pctM, fill = drach_status)) +
            geom_violin(alpha = 0.3, scale = "width") +
            geom_boxplot(width = 0.2, outlier.shape = NA) +
            facet_grid(rte_length_req ~ harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "Mean % m6A", fill = NULL,
                title = "m6A at polymorphic DRACH sites: FL vs truncated"
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none")
        mysaveandstore(str_glue("{outdir_evo_al}/drach_polymorphism_vs_modification_by_length.pdf"),
            w = 3 * min(length(poly_positions_all), 4) + 1,
            h = 6 + ceiling(length(poly_positions_all) / 4)
        )
    }

    # Plot 4: Modification at reliable sites, FL vs trnc comparison
    mod_by_site_length <- yl1annm_alllength_filtered %>%
        filter(condition == "N_ORF1", harmonized_pos %in% reliable_harmonized) %>%
        group_by(harmonized_pos, rte_length_req) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), n = n_distinct(gene_id), .groups = "drop")

    p <- mod_by_site_length %>%
        mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
        ggplot(aes(x = harmonized_pos, y = mm, fill = rte_length_req)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5) +
        scale_fill_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
        labs(
            x = "Reliable site (harmonized position)", y = "Mean % m6A (ORF1 RIP)",
            fill = "Length",
            title = "Modification at reliable m6A sites: FL vs truncated elements"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outdir_evo_al}/modification_reliable_sites_fl_vs_trnc.pdf"),
        w = max(7, length(reliable_harmonized) * 0.8 + 2), h = 5
    )

    # Plot 5: Modification vs DRACH conservation, all-length combined
    mod_at_sites_all <- yl1annm_alllength_filtered %>%
        filter(condition == "N_ORF1", harmonized_pos %in% reliable_harmonized) %>%
        group_by(harmonized_pos, rte_length_req) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop")

    p <- site_evolution %>%
        left_join(mod_at_sites_all, by = "harmonized_pos") %>%
        filter(!is.na(mm)) %>%
        ggplot(aes(x = n_subfams_drach, y = mm, color = category, shape = rte_length_req)) +
        geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        geom_text_repel(aes(label = harmonized_pos), size = 2.5, max.overlaps = 15) +
        labs(
            x = "# subfamilies with DRACH context", y = "Mean % m6A (ORF1 RIP)",
            color = "Category", shape = "Length",
            title = "DRACH conservation vs modification (FL + trnc)"
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_evo_al}/modification_vs_drach_conservation_alllength.pdf"), w = 9, h = 6)

    # Plot 6: Element count by length and subfamily at reliable sites
    element_count_at_sites <- drach_poly_df_al %>%
        filter(harmonized_pos %in% reliable_harmonized) %>%
        distinct(gene_id, harmonized_pos, rte_subfamily, rte_length_req) %>%
        count(harmonized_pos, rte_length_req)

    p <- element_count_at_sites %>%
        mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
        ggplot(aes(x = harmonized_pos, y = n, fill = rte_length_req)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
        labs(
            x = "Reliable site", y = "# elements with data",
            fill = "Length",
            title = "Element coverage at reliable sites: FL vs truncated"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outdir_evo_al}/element_coverage_by_length.pdf"),
        w = max(7, length(reliable_harmonized) * 0.8 + 2), h = 5
    )
}

cat("\n=== Section 4b complete ===\n")


# ============================================================
# //ANCHOR - SECTION 5: ELEMENT-LEVEL METHYLATION ANALYSIS
# ============================================================

run_element_analysis <- function(meth_data, enrich_data, subfamilies, reliable_sites,
                                 label, suffix, length_filter = "FL") {
    .global_mysaveandstore <- get("mysaveandstore", envir = globalenv())
    mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, store = store_var, raster = FALSE, ...) {
        .global_mysaveandstore(fn = fn, w = w, h = h, res = res, pl = p, store = store, raster = raster, ...)
    }

    outdir_elem <- str_glue("{outdir_base}/element_level{suffix}")
    intdir_elem <- str_glue("{intdir}{suffix}/element_level")
    for (d in c(outdir_elem, intdir_elem)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    length_label <- if (length_filter == "all") " (FL+trnc)" else " (FL)"
    full_label <- paste0(label, length_label)
    cat(sprintf("\n=== Section 5: Element-Level Methylation [%s] ===\n", full_label))

    reliable_hpos <- reliable_sites$harmonized_pos

    # 5a. Compute per-element per-site methylation at reliable sites
    elem_site <- meth_data %>%
        filter(condition == "N_ORF1") %>%
        filter(rte_subfamily %in% subfamilies) %>%
        filter(harmonized_pos %in% reliable_hpos) %>%
        filter(cov >= 5)
    if (length_filter == "FL") {
        elem_site <- elem_site %>% filter(rte_length_req == "FL")
    }

    elem_site_summary <- elem_site %>%
        group_by(gene_id, harmonized_pos, rte_subfamily, rte_length_req) %>%
        summarise(
            mean_pctM = mean(pctM, na.rm = TRUE),
            mean_cov = mean(cov, na.rm = TRUE),
            n_samples = n_distinct(sample),
            .groups = "drop"
        )

    # Per-subfamily per-site coverage summary (for annotation on plots)
    site_cov_summary <- elem_site_summary %>%
        group_by(rte_subfamily, harmonized_pos) %>%
        summarise(
            n_elements = n_distinct(gene_id),
            total_cov = sum(mean_cov, na.rm = TRUE),
            mean_cov_across_elements = mean(mean_cov, na.rm = TRUE),
            .groups = "drop"
        )

    # 5b. Per-element aggregate across all reliable sites
    elem_agg <- elem_site_summary %>%
        group_by(gene_id, rte_subfamily, rte_length_req) %>%
        summarise(
            n_sites_covered = n_distinct(harmonized_pos),
            mean_m6a = mean(mean_pctM, na.rm = TRUE),
            median_m6a = median(mean_pctM, na.rm = TRUE),
            max_m6a = max(mean_pctM, na.rm = TRUE),
            frac_sites_mod = mean(mean_pctM > 10, na.rm = TRUE),
            total_cov = sum(mean_cov, na.rm = TRUE),
            .groups = "drop"
        )

    # Join element annotations
    elem_agg <- elem_agg %>%
        left_join(
            rmann %>%
                dplyr::select(
                    gene_id, intactness_req, genic_loc, loc_superlowres_integrative_stranded,
                    pctdiv, length, pctconsensuscovered, rte_length_req
                ) %>%
                distinct(),
            by = c("gene_id", "rte_length_req")
        ) %>%
        left_join(
            enrich_data %>%
                dplyr::select(
                    gene_id, l2fc_condition_N_ORF1_vs_N_TOT,
                    baseMean_condition_N_ORF1_vs_N_TOT,
                    starts_with("mean_lr_")
                ),
            by = "gene_id"
        )

    write_tsv(elem_site_summary, str_glue("{intdir_elem}/element_site_methylation.tsv"))
    write_tsv(elem_agg, str_glue("{intdir_elem}/element_aggregate_methylation.tsv"))

    cat(sprintf("Elements with data: %d\n", nrow(elem_agg)))
    cat(sprintf(
        "Sites covered per element: median=%g, range=%d-%d\n",
        median(elem_agg$n_sites_covered), min(elem_agg$n_sites_covered), max(elem_agg$n_sites_covered)
    ))

    # --- PLOTS ---

    # Plot 1: Strip/violin of element-level pctM per reliable site
    # n labels per site (total across subfamilies)
    n_per_site <- site_cov_summary %>%
        group_by(harmonized_pos) %>%
        summarise(n_total = sum(n_elements), cov_total = round(sum(total_cov)), .groups = "drop") %>%
        mutate(
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
            n_label = paste0("n=", n_total, "\ncov=", cov_total)
        )
    p <- elem_site_summary %>%
        mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
        ggplot(aes(x = harmonized_pos, y = mean_pctM)) +
        geom_violin(fill = "grey85", alpha = 0.5, scale = "width") +
        geom_jitter(aes(color = rte_subfamily), width = 0.2, size = 0.6, alpha = 0.4) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "red", linewidth = 0.4) +
        geom_text(
            data = n_per_site, aes(x = harmonized_pos, label = n_label),
            y = -Inf, vjust = -0.3, size = 2.5, inherit.aes = FALSE
        ) +
        labs(
            x = "Reliable site (harmonized position)", y = "Mean % m6A per element",
            color = "Subfamily",
            title = sprintf("Element-level m6A at reliable sites (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outdir_elem}/element_m6a_per_site_violin.pdf"),
        w = max(8, length(reliable_hpos) * 0.8 + 2), h = 6
    )

    # Plot 2: Same but faceted by length (if all-length)
    if (length_filter == "all") {
        p <- elem_site_summary %>%
            mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
            ggplot(aes(x = harmonized_pos, y = mean_pctM)) +
            geom_violin(fill = "grey85", alpha = 0.5, scale = "width") +
            geom_jitter(aes(color = rte_subfamily), width = 0.2, size = 0.6, alpha = 0.4) +
            stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "red", linewidth = 0.4) +
            facet_wrap(~rte_length_req, ncol = 1) +
            labs(
                x = "Reliable site", y = "Mean % m6A per element",
                color = "Subfamily",
                title = sprintf("Element-level m6A by length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/element_m6a_per_site_violin_by_length.pdf"),
            w = max(8, length(reliable_hpos) * 0.8 + 2), h = 10
        )
    }

    # Plot 2b: Per-site element count by subfamily
    p <- site_cov_summary %>%
        mutate(
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
            rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)
        ) %>%
        ggplot(aes(x = rte_subfamily, y = n_elements, fill = rte_subfamily)) +
        geom_col(alpha = 0.7) +
        geom_text(aes(label = paste0(n_elements, "\n(", round(total_cov), "x)")), vjust = -0.1, size = 2) +
        facet_wrap(~harmonized_pos, scales = "free_y") +
        labs(
            x = NULL, y = "# elements with data",
            subtitle = "Labels: # elements (total coverage)",
            title = sprintf("Per-site coverage by subfamily (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none"
        )
    mysaveandstore(str_glue("{outdir_elem}/per_site_coverage_by_subfamily.pdf"),
        w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
        h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
    )

    # Plot 2c: Total read coverage by subfamily per site
    p <- site_cov_summary %>%
        mutate(
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
            rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)
        ) %>%
        ggplot(aes(x = rte_subfamily, y = total_cov, fill = rte_subfamily)) +
        geom_col(alpha = 0.7) +
        geom_text(aes(label = round(total_cov)), vjust = -0.3, size = 2.5) +
        facet_wrap(~harmonized_pos, scales = "free_y") +
        labs(
            x = NULL, y = "Total read coverage",
            title = sprintf("Per-site total coverage by subfamily (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none"
        )
    mysaveandstore(str_glue("{outdir_elem}/per_site_total_coverage_by_subfamily.pdf"),
        w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
        h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
    )

    # Plot 2d: Coverage summary heatmap - elements and total coverage per site per subfamily
    p <- site_cov_summary %>%
        mutate(
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
            rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER),
            cell_label = paste0(n_elements, " el\n", round(total_cov), "x")
        ) %>%
        ggplot(aes(x = harmonized_pos, y = rte_subfamily, fill = total_cov)) +
        geom_tile(color = "white") +
        geom_text(aes(label = cell_label), size = 2.5, lineheight = 0.85) +
        scale_fill_gradient(low = "white", high = "steelblue", name = "Total cov") +
        labs(
            x = "Reliable site (harmonized position)", y = "Subfamily",
            title = sprintf("Coverage per site per subfamily (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(str_glue("{outdir_elem}/site_subfamily_coverage_heatmap.pdf"),
        w = max(7, length(reliable_hpos) * 0.6 + 3), h = max(4, length(subfamilies) * 0.5 + 2)
    )

    # Plot 3: Heatmap of element x site methylation (top elements by coverage)
    hm_data <- elem_site_summary %>%
        group_by(gene_id) %>%
        filter(n_distinct(harmonized_pos) >= max(2, length(reliable_hpos) * 0.5)) %>%
        ungroup()

    if (n_distinct(hm_data$gene_id) >= 5) {
        hm_wide <- hm_data %>%
            dplyr::select(gene_id, harmonized_pos, mean_pctM) %>%
            pivot_wider(names_from = harmonized_pos, values_from = mean_pctM, names_prefix = "pos_")

        hm_mat <- as.matrix(hm_wide %>% dplyr::select(-gene_id))
        rownames(hm_mat) <- hm_wide$gene_id

        # Limit to manageable size
        if (nrow(hm_mat) > 200) {
            row_vars <- apply(hm_mat, 1, var, na.rm = TRUE)
            keep_rows <- names(sort(row_vars, decreasing = TRUE))[1:200]
            hm_mat <- hm_mat[keep_rows, , drop = FALSE]
        }

        # Row annotations
        hm_ann_df <- elem_agg %>%
            filter(gene_id %in% rownames(hm_mat)) %>%
            dplyr::select(gene_id, rte_subfamily, intactness_req, genic_loc, rte_length_req) %>%
            distinct() %>%
            column_to_rownames("gene_id")
        hm_ann_df <- hm_ann_df[rownames(hm_mat), , drop = FALSE]

        ha <- ComplexHeatmap::rowAnnotation(
            Subfamily = hm_ann_df$rte_subfamily,
            Intactness = hm_ann_df$intactness_req,
            Location = hm_ann_df$genic_loc,
            Length = hm_ann_df$rte_length_req
        )

        p <- ComplexHeatmap::Heatmap(
            hm_mat,
            name = "% m6A",
            col = circlize::colorRamp2(c(0, 25, 50, 100), c("white", "#FFEDA0", "#FEB24C", "#F03B20")),
            cluster_columns = FALSE,
            show_row_names = FALSE,
            left_annotation = ha,
            column_title = sprintf("Element x site m6A (%s)", full_label),
            na_col = "grey95"
        )
        pdf(str_glue("{outdir_elem}/element_site_heatmap.pdf"),
            width = max(8, length(reliable_hpos) * 0.6 + 4),
            height = max(8, min(nrow(hm_mat) * 0.04 + 3, 20))
        )
        draw(p)
        dev.off()
    }

    # Plot 3b: Pairwise correlation heatmap of methylation across reliable sites
    cor_wide <- elem_site_summary %>%
        dplyr::select(gene_id, harmonized_pos, mean_pctM) %>%
        pivot_wider(names_from = harmonized_pos, values_from = mean_pctM, names_prefix = "pos_")
    cor_mat_vals <- cor_wide %>%
        dplyr::select(-gene_id) %>%
        as.matrix()
    if (ncol(cor_mat_vals) >= 2) {
        cor_mat <- cor(cor_mat_vals, use = "pairwise.complete.obs", method = "spearman")
        # n observations per pair
        n_mat <- crossprod(!is.na(cor_mat_vals))

        col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
        # Cell labels: r value + n
        cell_labels <- matrix(
            sprintf("%.2f\n(%d)", as.vector(cor_mat), as.vector(n_mat)),
            nrow = nrow(cor_mat)
        )
        p <- ComplexHeatmap::Heatmap(
            cor_mat,
            name = "Spearman r",
            col = col_fun,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid::grid.text(cell_labels[i, j], x, y, gp = grid::gpar(fontsize = 7))
            },
            cluster_rows = nrow(cor_mat) > 2,
            cluster_columns = ncol(cor_mat) > 2,
            column_title = sprintf("Site-site methylation correlation (%s)", full_label),
            row_names_gp = grid::gpar(fontsize = 8),
            column_names_gp = grid::gpar(fontsize = 8)
        )
        pdf(str_glue("{outdir_elem}/site_methylation_correlation_heatmap.pdf"),
            width = max(6, ncol(cor_mat) * 0.8 + 3),
            height = max(6, nrow(cor_mat) * 0.8 + 3)
        )
        draw(p)
        dev.off()
    }

    # Plot 4: Distribution of aggregate element methylation
    p <- elem_agg %>%
        ggplot(aes(x = mean_m6a)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
        geom_vline(
            xintercept = median(elem_agg$mean_m6a, na.rm = TRUE),
            linetype = "dashed", color = "red"
        ) +
        labs(
            x = "Mean % m6A across reliable sites", y = "# elements",
            title = sprintf("Distribution of aggregate element methylation (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_distribution.pdf"), w = 7, h = 5)

    # Plot 5: Aggregate methylation by subfamily
    p <- elem_agg %>%
        mutate(rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)) %>%
        ggplot(aes(x = rte_subfamily, y = mean_m6a, fill = rte_subfamily)) +
        geom_violin(alpha = 0.4, scale = "width") +
        geom_boxplot(width = 0.15, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 0.5, alpha = 0.2) +
        labs(
            x = "Subfamily", y = "Mean % m6A across reliable sites",
            title = sprintf("Aggregate element m6A by subfamily (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")
    mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_by_subfamily.pdf"), w = 7, h = 5)

    # Plot 6: Aggregate methylation by intactness
    if (n_distinct(na.omit(elem_agg$intactness_req)) > 1) {
        p <- elem_agg %>%
            filter(!is.na(intactness_req)) %>%
            ggplot(aes(x = intactness_req, y = mean_m6a, fill = intactness_req)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_jitter(width = 0.15, size = 0.5, alpha = 0.2) +
            labs(
                x = "Intactness", y = "Mean % m6A across reliable sites",
                title = sprintf("Aggregate m6A by intactness (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_by_intactness.pdf"), w = 7, h = 5)
    }

    # Plot 7: Aggregate methylation by genic location
    if (n_distinct(na.omit(elem_agg$genic_loc)) > 1) {
        p <- elem_agg %>%
            filter(!is.na(genic_loc)) %>%
            ggplot(aes(x = genic_loc, y = mean_m6a, fill = genic_loc)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_jitter(width = 0.15, size = 0.5, alpha = 0.2) +
            labs(
                x = "Genic location", y = "Mean % m6A across reliable sites",
                title = sprintf("Aggregate m6A by genic location (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_by_genic_loc.pdf"), w = 7, h = 5)
    }

    # Plot 8: Aggregate methylation by loc_superlowres_integrative_stranded
    if (n_distinct(na.omit(elem_agg$loc_superlowres_integrative_stranded)) > 1) {
        p <- elem_agg %>%
            filter(!is.na(loc_superlowres_integrative_stranded)) %>%
            ggplot(aes(x = loc_superlowres_integrative_stranded, y = mean_m6a, fill = loc_superlowres_integrative_stranded)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_jitter(width = 0.15, size = 0.5, alpha = 0.2) +
            labs(
                x = "Integrative location", y = "Mean % m6A across reliable sites",
                title = sprintf("Aggregate m6A by integrative location (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_by_loc_superlowres_integrative_stranded.pdf"), w = 8, h = 5)
    }

    # Plot 8b: Per-site violin by subfamily (when multiple subfamilies)
    if (length(subfamilies) > 1) {
        # Prepare n labels for each subfamily x site
        n_labels_8b <- site_cov_summary %>%
            mutate(
                harmonized_pos = factor(harmonized_pos, levels = sort(unique(site_cov_summary$harmonized_pos))),
                rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER),
                n_label = paste0("n=", n_elements, "\n", round(total_cov), "x")
            )
        p <- elem_site_summary %>%
            mutate(
                harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
                rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)
            ) %>%
            ggplot(aes(x = rte_subfamily, y = mean_pctM, fill = rte_subfamily)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_text(
                data = n_labels_8b, aes(x = rte_subfamily, label = n_label),
                y = -Inf, vjust = -0.3, size = 2, inherit.aes = FALSE
            ) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "% m6A per element", fill = "Subfamily",
                title = sprintf("Per-site m6A by subfamily (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
        mysaveandstore(str_glue("{outdir_elem}/per_site_m6a_by_subfamily.pdf"),
            w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
            h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
        )
    }

    # Join element annotations for per-site plots and differential methylation
    elem_site_ann <- elem_site_summary %>%
        left_join(
            rmann %>%
                dplyr::select(gene_id, intactness_req, genic_loc, loc_superlowres_integrative_stranded, rte_length_req) %>%
                distinct(),
            by = c("gene_id", "rte_length_req")
        )

    # Helper to create n-label data for faceted per-site plots
    make_n_labels <- function(data, group_var) {
        data %>%
            filter(!is.na(.data[[group_var]])) %>%
            group_by(harmonized_pos, .data[[group_var]]) %>%
            summarise(n = n(), cov = round(sum(mean_cov, na.rm = TRUE)), .groups = "drop") %>%
            mutate(
                harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
                n_label = paste0("n=", n, "\n", cov, "x")
            )
    }

    # Plot 8c: Per-site violin by intactness
    if (n_distinct(na.omit(elem_site_ann$intactness_req)) > 1) {
        n_labels_8c <- make_n_labels(elem_site_ann, "intactness_req")
        p <- elem_site_ann %>%
            filter(!is.na(intactness_req)) %>%
            mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
            ggplot(aes(x = intactness_req, y = mean_pctM, fill = intactness_req)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_text(
                data = n_labels_8c, aes(x = .data[["intactness_req"]], label = n_label),
                y = -Inf, vjust = -0.3, size = 2, inherit.aes = FALSE
            ) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "% m6A per element", fill = "Intactness",
                title = sprintf("Per-site m6A by intactness (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
        mysaveandstore(str_glue("{outdir_elem}/per_site_m6a_by_intactness.pdf"),
            w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
            h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
        )
    }

    # Plot 8d: Per-site violin by location
    if (n_distinct(na.omit(elem_site_ann$loc_superlowres_integrative_stranded)) > 1) {
        n_labels_8d <- make_n_labels(elem_site_ann, "loc_superlowres_integrative_stranded")
        p <- elem_site_ann %>%
            filter(!is.na(loc_superlowres_integrative_stranded)) %>%
            mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
            ggplot(aes(x = loc_superlowres_integrative_stranded, y = mean_pctM, fill = loc_superlowres_integrative_stranded)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_text(
                data = n_labels_8d, aes(x = .data[["loc_superlowres_integrative_stranded"]], label = n_label),
                y = -Inf, vjust = -0.3, size = 1.8, inherit.aes = FALSE
            ) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "% m6A per element", fill = "Location",
                title = sprintf("Per-site m6A by genomic location (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                legend.position = "none"
            )
        mysaveandstore(str_glue("{outdir_elem}/per_site_m6a_by_location.pdf"),
            w = max(10, ceiling(sqrt(length(reliable_hpos))) * 4 + 2),
            h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
        )
    }

    # Plot 8e: Per-site violin by genic location
    if (n_distinct(na.omit(elem_site_ann$genic_loc)) > 1) {
        n_labels_8e <- make_n_labels(elem_site_ann, "genic_loc")
        p <- elem_site_ann %>%
            filter(!is.na(genic_loc)) %>%
            mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
            ggplot(aes(x = genic_loc, y = mean_pctM, fill = genic_loc)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_text(
                data = n_labels_8e, aes(x = .data[["genic_loc"]], label = n_label),
                y = -Inf, vjust = -0.3, size = 2, inherit.aes = FALSE
            ) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = NULL, y = "% m6A per element", fill = "Genic loc",
                title = sprintf("Per-site m6A by genic location (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
        mysaveandstore(str_glue("{outdir_elem}/per_site_m6a_by_genic_loc.pdf"),
            w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
            h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
        )
    }

    # Plot 9: Scatter - aggregate m6A vs ORF1 enrichment (l2fc)
    p <- elem_agg %>%
        filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT)) %>%
        ggplot(aes(x = l2fc_condition_N_ORF1_vs_N_TOT, y = mean_m6a, color = rte_subfamily)) +
        geom_point(size = 1.5, alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        labs(
            x = "ORF1 RIP enrichment (log2FC vs Total)", y = "Mean % m6A across reliable sites",
            color = "Subfamily",
            title = sprintf("Element m6A vs ORF1 enrichment (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/m6a_vs_orf1_enrichment.pdf"), w = 8, h = 6)

    # Plot 10: Scatter - aggregate m6A vs expression (baseMean)
    p <- elem_agg %>%
        filter(
            !is.na(baseMean_condition_N_ORF1_vs_N_TOT),
            baseMean_condition_N_ORF1_vs_N_TOT > 0
        ) %>%
        ggplot(aes(x = log10(baseMean_condition_N_ORF1_vs_N_TOT), y = mean_m6a, color = rte_subfamily)) +
        geom_point(size = 1.5, alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        labs(
            x = "Expression (log10 baseMean)", y = "Mean % m6A across reliable sites",
            color = "Subfamily",
            title = sprintf("Element m6A vs expression (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/m6a_vs_expression.pdf"), w = 8, h = 6)

    # Plot 11: Scatter - aggregate m6A vs divergence
    p <- elem_agg %>%
        filter(!is.na(pctdiv)) %>%
        ggplot(aes(x = pctdiv, y = mean_m6a, color = rte_subfamily)) +
        geom_point(size = 1.5, alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        labs(
            x = "% divergence from consensus", y = "Mean % m6A across reliable sites",
            color = "Subfamily",
            title = sprintf("Element m6A vs sequence divergence (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/m6a_vs_divergence.pdf"), w = 8, h = 6)

    # Plot 12: Scatter - aggregate m6A vs element length
    p <- elem_agg %>%
        filter(!is.na(length)) %>%
        ggplot(aes(x = length, y = mean_m6a, color = rte_subfamily)) +
        geom_point(size = 1.5, alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        labs(
            x = "Element length (bp)", y = "Mean % m6A across reliable sites",
            color = "Subfamily",
            title = sprintf("Element m6A vs element length (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/m6a_vs_length.pdf"), w = 8, h = 6)

    # Plot 12b: Scatter - aggregate m6A vs Total RNA expression
    if ("mean_lr_N_TOT" %in% colnames(elem_agg)) {
        p <- elem_agg %>%
            filter(!is.na(mean_lr_N_TOT), mean_lr_N_TOT > 0) %>%
            ggplot(aes(x = log10(mean_lr_N_TOT), y = mean_m6a, color = rte_subfamily)) +
            geom_point(size = 1.5, alpha = 0.5) +
            geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
            labs(
                x = "Total RNA expression (log10 normed count)", y = "Mean % m6A across reliable sites",
                color = "Subfamily",
                title = sprintf("Element m6A vs Total RNA expression (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/m6a_vs_total_rna_expression.pdf"), w = 8, h = 6)
    }

    # Plot 12c: Scatter - aggregate m6A vs ORF1 RIP expression
    if ("mean_lr_N_ORF1" %in% colnames(elem_agg)) {
        p <- elem_agg %>%
            filter(!is.na(mean_lr_N_ORF1), mean_lr_N_ORF1 > 0) %>%
            ggplot(aes(x = log10(mean_lr_N_ORF1), y = mean_m6a, color = rte_subfamily)) +
            geom_point(size = 1.5, alpha = 0.5) +
            geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
            labs(
                x = "ORF1 RIP expression (log10 normed count)", y = "Mean % m6A across reliable sites",
                color = "Subfamily",
                title = sprintf("Element m6A vs ORF1 RIP expression (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/m6a_vs_orf1_rip_expression.pdf"), w = 8, h = 6)
    }

    # Plot 12d: Total vs ORF1 expression colored by m6A level
    if ("mean_lr_N_TOT" %in% colnames(elem_agg) && "mean_lr_N_ORF1" %in% colnames(elem_agg)) {
        p <- elem_agg %>%
            filter(
                !is.na(mean_lr_N_TOT), !is.na(mean_lr_N_ORF1),
                mean_lr_N_TOT > 0, mean_lr_N_ORF1 > 0
            ) %>%
            ggplot(aes(x = log10(mean_lr_N_TOT), y = log10(mean_lr_N_ORF1), color = mean_m6a)) +
            geom_point(size = 1.5, alpha = 0.6) +
            scale_color_gradient2(
                low = "steelblue", mid = "grey80", high = "#D62728",
                midpoint = median(elem_agg$mean_m6a, na.rm = TRUE), name = "Mean % m6A"
            ) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                x = "Total RNA (log10 normed count)", y = "ORF1 RIP (log10 normed count)",
                title = sprintf("Expression landscape colored by m6A (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/total_vs_orf1_expression_by_m6a.pdf"), w = 8, h = 7)
    }

    # --- FL vs trnc comparisons (all-length mode only) ---
    if (length_filter == "all" && n_distinct(na.omit(elem_agg$rte_length_req)) > 1) {
        # Aggregate m6A by FL vs trnc
        p <- elem_agg %>%
            ggplot(aes(x = rte_length_req, y = mean_m6a, fill = rte_length_req)) +
            geom_violin(alpha = 0.4, scale = "width") +
            geom_boxplot(width = 0.15, outlier.shape = NA) +
            geom_jitter(width = 0.15, size = 0.5, alpha = 0.2) +
            scale_fill_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
            labs(
                x = NULL, y = "Mean % m6A across reliable sites",
                title = sprintf("Aggregate m6A: FL vs truncated (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none")
        mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_fl_vs_trnc.pdf"), w = 5, h = 5)

        # By subfamily, faceted by length
        p <- elem_agg %>%
            mutate(rte_subfamily = factor(rte_subfamily, levels = SUBFAM_ORDER)) %>%
            ggplot(aes(x = rte_subfamily, y = mean_m6a, fill = rte_length_req)) +
            geom_violin(alpha = 0.4, scale = "width", position = position_dodge(width = 0.8)) +
            geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
            labs(
                x = "Subfamily", y = "Mean % m6A across reliable sites",
                fill = "Length",
                title = sprintf("Aggregate m6A by subfamily and length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/aggregate_m6a_by_subfamily_and_length.pdf"), w = 9, h = 5)

        # Key scatter plots faceted by length
        # m6A vs enrichment by length
        p <- elem_agg %>%
            filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT)) %>%
            ggplot(aes(x = l2fc_condition_N_ORF1_vs_N_TOT, y = mean_m6a, color = rte_length_req)) +
            geom_point(size = 1, alpha = 0.4) +
            geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
            scale_color_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
            labs(
                x = "ORF1 RIP enrichment (log2FC vs Total)", y = "Mean % m6A",
                color = "Length",
                title = sprintf("m6A vs enrichment by length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/m6a_vs_enrichment_by_length.pdf"), w = 8, h = 6)

        # m6A vs divergence by length
        p <- elem_agg %>%
            filter(!is.na(pctdiv)) %>%
            ggplot(aes(x = pctdiv, y = mean_m6a, color = rte_length_req)) +
            geom_point(size = 1, alpha = 0.4) +
            geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
            scale_color_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
            labs(
                x = "% divergence from consensus", y = "Mean % m6A",
                color = "Length",
                title = sprintf("m6A vs divergence by length (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/m6a_vs_divergence_by_length.pdf"), w = 8, h = 6)

        # m6A vs Total RNA expression by length
        if ("mean_lr_N_TOT" %in% colnames(elem_agg)) {
            p <- elem_agg %>%
                filter(!is.na(mean_lr_N_TOT), mean_lr_N_TOT > 0) %>%
                ggplot(aes(x = log10(mean_lr_N_TOT), y = mean_m6a, color = rte_length_req)) +
                geom_point(size = 1, alpha = 0.4) +
                geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
                scale_color_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
                labs(
                    x = "Total RNA expression (log10)", y = "Mean % m6A",
                    color = "Length",
                    title = sprintf("m6A vs Total RNA by length (%s)", full_label)
                ) +
                theme_minimal(base_size = 12)
            mysaveandstore(str_glue("{outdir_elem}/m6a_vs_total_rna_by_length.pdf"), w = 8, h = 6)
        }

        # m6A vs ORF1 RIP expression by length
        if ("mean_lr_N_ORF1" %in% colnames(elem_agg)) {
            p <- elem_agg %>%
                filter(!is.na(mean_lr_N_ORF1), mean_lr_N_ORF1 > 0) %>%
                ggplot(aes(x = log10(mean_lr_N_ORF1), y = mean_m6a, color = rte_length_req)) +
                geom_point(size = 1, alpha = 0.4) +
                geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
                scale_color_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
                labs(
                    x = "ORF1 RIP expression (log10)", y = "Mean % m6A",
                    color = "Length",
                    title = sprintf("m6A vs ORF1 RIP by length (%s)", full_label)
                ) +
                theme_minimal(base_size = 12)
            mysaveandstore(str_glue("{outdir_elem}/m6a_vs_orf1_rip_by_length.pdf"), w = 8, h = 6)
        }

        # Multi-feature panel faceted by length
        if (nrow(elem_long) > 0) {
            p <- elem_long %>%
                left_join(elem_agg %>% dplyr::select(gene_id, rte_length_req), by = "gene_id") %>%
                ggplot(aes(x = value, y = mean_m6a, color = rte_length_req)) +
                geom_point(size = 0.6, alpha = 0.3) +
                geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.5) +
                scale_color_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
                facet_wrap(~feature, scales = "free_x", ncol = 3) +
                labs(
                    x = NULL, y = "Mean % m6A across reliable sites",
                    color = "Length",
                    title = sprintf("Element m6A vs features by length (%s)", full_label)
                ) +
                theme_minimal(base_size = 12) +
                theme(strip.text = element_text(size = 9))
            mysaveandstore(str_glue("{outdir_elem}/m6a_vs_features_panel_by_length.pdf"), w = 14, h = 10)
        }
    }

    # Plot 13: Correlation matrix of per-site methylation across elements
    if (n_distinct(hm_data$gene_id) >= 10 && length(reliable_hpos) >= 2) {
        cor_mat <- hm_data %>%
            dplyr::select(gene_id, harmonized_pos, mean_pctM) %>%
            pivot_wider(names_from = harmonized_pos, values_from = mean_pctM) %>%
            dplyr::select(-gene_id) %>%
            cor(use = "pairwise.complete.obs")

        if (all(!is.na(cor_mat))) {
            col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("steelblue", "white", "#D62728"))
            p <- ComplexHeatmap::Heatmap(
                cor_mat,
                name = "Pearson r",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_title = sprintf("Site-site m6A correlation across elements (%s)", full_label),
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = grid::gpar(fontsize = 8))
                }
            )
            pdf(str_glue("{outdir_elem}/site_correlation_matrix.pdf"),
                width = max(6, length(reliable_hpos) * 0.8 + 3),
                height = max(6, length(reliable_hpos) * 0.8 + 3)
            )
            draw(p)
            dev.off()
        }
    }

    # Plot 14: Multi-feature scatter panel (cowplot grid)
    # Faceted scatter of mean_m6a vs continuous features
    elem_long <- elem_agg %>%
        dplyr::select(
            gene_id, mean_m6a, rte_subfamily,
            l2fc_condition_N_ORF1_vs_N_TOT, baseMean_condition_N_ORF1_vs_N_TOT,
            pctdiv, length, pctconsensuscovered,
            any_of(c("mean_lr_N_TOT", "mean_lr_N_ORF1"))
        ) %>%
        mutate(
            log10_baseMean = log10(pmax(baseMean_condition_N_ORF1_vs_N_TOT, 0.1)),
            log10_total_rna = if ("mean_lr_N_TOT" %in% names(.)) log10(pmax(mean_lr_N_TOT, 0.1)) else NA_real_,
            log10_orf1_rip = if ("mean_lr_N_ORF1" %in% names(.)) log10(pmax(mean_lr_N_ORF1, 0.1)) else NA_real_
        ) %>%
        dplyr::select(-baseMean_condition_N_ORF1_vs_N_TOT, -any_of(c("mean_lr_N_TOT", "mean_lr_N_ORF1"))) %>%
        pivot_longer(
            cols = c(
                l2fc_condition_N_ORF1_vs_N_TOT, log10_baseMean, pctdiv, length, pctconsensuscovered,
                log10_total_rna, log10_orf1_rip
            ),
            names_to = "feature", values_to = "value"
        ) %>%
        mutate(feature = case_when(
            feature == "l2fc_condition_N_ORF1_vs_N_TOT" ~ "ORF1 enrichment (l2fc)",
            feature == "log10_baseMean" ~ "DESeq2 baseMean (log10)",
            feature == "log10_total_rna" ~ "Total RNA (log10)",
            feature == "log10_orf1_rip" ~ "ORF1 RIP (log10)",
            feature == "pctdiv" ~ "% divergence",
            feature == "length" ~ "Element length",
            feature == "pctconsensuscovered" ~ "% consensus covered",
            TRUE ~ feature
        )) %>%
        filter(!is.na(value))

    p <- elem_long %>%
        ggplot(aes(x = value, y = mean_m6a)) +
        geom_point(aes(color = rte_subfamily), size = 0.8, alpha = 0.3) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
        facet_wrap(~feature, scales = "free_x", ncol = 3) +
        labs(
            x = NULL, y = "Mean % m6A across reliable sites",
            color = "Subfamily",
            title = sprintf("Element m6A vs features (%s)", full_label)
        ) +
        theme_minimal(base_size = 12) +
        theme(strip.text = element_text(size = 9))
    mysaveandstore(str_glue("{outdir_elem}/m6a_vs_features_panel.pdf"), w = 14, h = 10)

    # Plot 15: Per-site m6A vs ORF1 enrichment (faceted by site)
    elem_site_enrich <- elem_site_summary %>%
        left_join(
            enrich_data %>%
                dplyr::select(gene_id, l2fc_condition_N_ORF1_vs_N_TOT),
            by = "gene_id"
        ) %>%
        filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT))

    if (nrow(elem_site_enrich) > 0) {
        p <- elem_site_enrich %>%
            mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
            ggplot(aes(x = l2fc_condition_N_ORF1_vs_N_TOT, y = mean_pctM)) +
            geom_point(aes(color = rte_subfamily), size = 0.8, alpha = 0.4) +
            geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(
                x = "ORF1 enrichment (log2FC)", y = "% m6A",
                color = "Subfamily",
                title = sprintf("Per-site m6A vs ORF1 enrichment (%s)", full_label)
            ) +
            theme_minimal(base_size = 12)
        mysaveandstore(str_glue("{outdir_elem}/per_site_m6a_vs_enrichment.pdf"),
            w = max(8, ceiling(sqrt(length(reliable_hpos))) * 3 + 2),
            h = max(6, ceiling(length(reliable_hpos) / ceiling(sqrt(length(reliable_hpos)))) * 3 + 1)
        )
    }

    # Plot 16: Number of sites covered per element (histogram)
    p <- elem_agg %>%
        ggplot(aes(x = n_sites_covered)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
        labs(
            x = "# reliable sites covered", y = "# elements",
            title = sprintf("Site coverage per element (%s)", full_label)
        ) +
        theme_minimal(base_size = 12)
    mysaveandstore(str_glue("{outdir_elem}/sites_covered_per_element.pdf"), w = 7, h = 5)

    # ============================================================
    # 5c. DIFFERENTIAL METHYLATION PER SITE BY ELEMENT FEATURES
    # ============================================================
    cat(sprintf("--- Differential methylation testing [%s] ---\n", full_label))

    # Helper: per-site Wilcoxon test for a binary/categorical grouping variable
    run_site_diffmeth <- function(data, group_var, group_label) {
        data <- data %>% filter(!is.na(!!sym(group_var)))
        groups <- unique(data[[group_var]])
        if (length(groups) < 2) {
            return(NULL)
        }

        results <- map_dfr(sort(unique(data$harmonized_pos)), function(hp) {
            site_data <- data %>% filter(harmonized_pos == hp)
            site_groups <- unique(site_data[[group_var]])
            if (length(site_groups) < 2 || min(table(site_data[[group_var]])) < 3) {
                return(tibble(
                    harmonized_pos = hp, group_var = group_label,
                    p_value = NA_real_, n_total = nrow(site_data)
                ))
            }
            # For 2 groups: Wilcoxon, for >2 groups: Kruskal-Wallis
            if (length(site_groups) == 2) {
                tt <- wilcox.test(mean_pctM ~ get(group_var), data = site_data)
            } else {
                tt <- kruskal.test(mean_pctM ~ get(group_var), data = site_data)
            }
            # Per-group summaries
            grp_summary <- site_data %>%
                group_by(!!sym(group_var)) %>%
                summarise(
                    mean_m6a = mean(mean_pctM, na.rm = TRUE),
                    median_m6a = median(mean_pctM, na.rm = TRUE),
                    n = n(), .groups = "drop"
                )
            tibble(
                harmonized_pos = hp,
                group_var = group_label,
                p_value = tt$p.value,
                n_total = nrow(site_data),
                group_summary = list(grp_summary)
            )
        })
        if (nrow(results) > 0) {
            results <- results %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
        }
        results
    }

    # Test 1: Intact vs Not Intact
    diffmeth_intactness <- run_site_diffmeth(elem_site_ann, "intactness_req", "Intactness")

    # Test 2: Genic vs Intergenic
    diffmeth_genic <- run_site_diffmeth(elem_site_ann, "genic_loc", "Genic location")

    # Test 3: FL vs trnc (only meaningful in all-length mode)
    diffmeth_length <- NULL
    if (length_filter == "all") {
        diffmeth_length <- run_site_diffmeth(elem_site_ann, "rte_length_req", "Length (FL vs trnc)")
    }

    # Combine all differential methylation results
    diffmeth_all <- bind_rows(diffmeth_intactness, diffmeth_genic, diffmeth_length) %>%
        filter(!is.na(p_value))

    write_tsv(
        diffmeth_all %>% dplyr::select(-group_summary),
        str_glue("{intdir_elem}/differential_methylation_tests.tsv")
    )

    cat(sprintf(
        "Differential tests: %d total, %d significant (padj < 0.05)\n",
        nrow(diffmeth_all), sum(diffmeth_all$p_adj < 0.05, na.rm = TRUE)
    ))

    # --- DIFFERENTIAL METHYLATION PLOTS ---

    # Plot DM1: Summary tile - significance across sites and features
    if (nrow(diffmeth_all) > 0) {
        p <- diffmeth_all %>%
            mutate(
                harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
                sig_label = case_when(
                    p_adj < 0.001 ~ "***",
                    p_adj < 0.01 ~ "**",
                    p_adj < 0.05 ~ "*",
                    TRUE ~ ""
                ),
                neg_log_padj = -log10(pmax(p_adj, 1e-20))
            ) %>%
            ggplot(aes(x = harmonized_pos, y = group_var, fill = neg_log_padj)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = sig_label), size = 4) +
            scale_fill_gradient(low = "grey95", high = "#D62728", name = "-log10(padj)") +
            labs(
                x = "Site (harmonized position)", y = NULL,
                title = sprintf("Differential m6A by element features (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/diffmeth_summary_tile.pdf"),
            w = max(7, length(reliable_hpos) * 0.6 + 3), h = max(3, n_distinct(diffmeth_all$group_var) * 0.8 + 2)
        )
    }

    # Plot DM2: Per-site violin by intactness
    if (!is.null(diffmeth_intactness) && any(diffmeth_intactness$p_adj < 0.1, na.rm = TRUE)) {
        sig_sites_int <- diffmeth_intactness %>%
            filter(p_adj < 0.1) %>%
            pull(harmonized_pos)
        if (length(sig_sites_int) > 0) {
            p <- elem_site_ann %>%
                filter(harmonized_pos %in% sig_sites_int, !is.na(intactness_req)) %>%
                mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
                ggplot(aes(x = intactness_req, y = mean_pctM, fill = intactness_req)) +
                geom_violin(alpha = 0.4, scale = "width") +
                geom_boxplot(width = 0.15, outlier.shape = NA) +
                geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
                facet_wrap(~harmonized_pos, scales = "free_y") +
                labs(
                    x = NULL, y = "% m6A", fill = "Intactness",
                    title = sprintf("Differential m6A by intactness (padj < 0.1) (%s)", full_label)
                ) +
                theme_minimal(base_size = 12) +
                theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
            mysaveandstore(str_glue("{outdir_elem}/diffmeth_intactness_violin.pdf"),
                w = max(6, ceiling(sqrt(length(sig_sites_int))) * 3 + 1),
                h = max(5, ceiling(length(sig_sites_int) / ceiling(sqrt(length(sig_sites_int)))) * 3 + 1)
            )
        }
    }

    # Plot DM3: Per-site violin by genic location
    if (!is.null(diffmeth_genic) && any(diffmeth_genic$p_adj < 0.1, na.rm = TRUE)) {
        sig_sites_gen <- diffmeth_genic %>%
            filter(p_adj < 0.1) %>%
            pull(harmonized_pos)
        if (length(sig_sites_gen) > 0) {
            p <- elem_site_ann %>%
                filter(harmonized_pos %in% sig_sites_gen, !is.na(genic_loc)) %>%
                mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
                ggplot(aes(x = genic_loc, y = mean_pctM, fill = genic_loc)) +
                geom_violin(alpha = 0.4, scale = "width") +
                geom_boxplot(width = 0.15, outlier.shape = NA) +
                geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
                facet_wrap(~harmonized_pos, scales = "free_y") +
                labs(
                    x = NULL, y = "% m6A", fill = "Genic location",
                    title = sprintf("Differential m6A by genic location (padj < 0.1) (%s)", full_label)
                ) +
                theme_minimal(base_size = 12) +
                theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
            mysaveandstore(str_glue("{outdir_elem}/diffmeth_genic_loc_violin.pdf"),
                w = max(6, ceiling(sqrt(length(sig_sites_gen))) * 3 + 1),
                h = max(5, ceiling(length(sig_sites_gen) / ceiling(sqrt(length(sig_sites_gen)))) * 3 + 1)
            )
        }
    }

    # Plot DM4: Per-site violin by FL vs trnc (all-length mode only)
    if (!is.null(diffmeth_length) && any(diffmeth_length$p_adj < 0.1, na.rm = TRUE)) {
        sig_sites_len <- diffmeth_length %>%
            filter(p_adj < 0.1) %>%
            pull(harmonized_pos)
        if (length(sig_sites_len) > 0) {
            p <- elem_site_ann %>%
                filter(harmonized_pos %in% sig_sites_len) %>%
                mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))) %>%
                ggplot(aes(x = rte_length_req, y = mean_pctM, fill = rte_length_req)) +
                geom_violin(alpha = 0.4, scale = "width") +
                geom_boxplot(width = 0.15, outlier.shape = NA) +
                geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
                scale_fill_manual(values = c("FL" = "#D62728", "trnc" = "#1F77B4")) +
                facet_wrap(~harmonized_pos, scales = "free_y") +
                labs(
                    x = NULL, y = "% m6A", fill = "Length",
                    title = sprintf("Differential m6A: FL vs truncated (padj < 0.1) (%s)", full_label)
                ) +
                theme_minimal(base_size = 12) +
                theme(legend.position = "bottom")
            mysaveandstore(str_glue("{outdir_elem}/diffmeth_fl_vs_trnc_violin.pdf"),
                w = max(6, ceiling(sqrt(length(sig_sites_len))) * 3 + 1),
                h = max(5, ceiling(length(sig_sites_len) / ceiling(sqrt(length(sig_sites_len)))) * 3 + 1)
            )
        }
    }

    # Plot DM5: All sites, all features - dot plot with effect direction
    # Show mean difference (or largest group mean - smallest) plus significance
    effect_summary <- diffmeth_all %>%
        rowwise() %>%
        mutate(
            max_grp_mean = max(group_summary$mean_m6a),
            min_grp_mean = min(group_summary$mean_m6a),
            effect_size = max_grp_mean - min_grp_mean
        ) %>%
        ungroup()

    if (nrow(effect_summary) > 0) {
        p <- effect_summary %>%
            mutate(
                harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))),
                significant = p_adj < 0.05
            ) %>%
            ggplot(aes(x = harmonized_pos, y = effect_size, color = significant, shape = group_var)) +
            geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.5)) +
            scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "grey60")) +
            labs(
                x = "Site (harmonized position)",
                y = "Effect size (max - min group mean % m6A)",
                color = "Significant (padj<0.05)", shape = "Feature",
                title = sprintf("Differential m6A effect sizes (%s)", full_label)
            ) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(str_glue("{outdir_elem}/diffmeth_effect_sizes.pdf"),
            w = max(8, length(reliable_hpos) * 0.6 + 3), h = 6
        )
    }

    cat(sprintf("Section 5 complete [%s]: %d elements\n", full_label, nrow(elem_agg)))

    list(elem_agg = elem_agg, elem_site_summary = elem_site_summary, diffmeth = diffmeth_all)
}

# --- Run Section 5 for all combinations ---
# FL-only
elem_all_yl1 <- run_element_analysis(
    yl1annm_filtered, enrichment_lr, YOUNG_L1_SUBFAMILIES,
    results_all_yl1$reliable_sites, "Young L1",
    suffix = "", length_filter = "FL"
)
elem_l1hs <- run_element_analysis(
    yl1annm_filtered, enrichment_lr, "L1HS",
    results_l1hs$reliable_sites, "L1HS",
    suffix = "_L1HS", length_filter = "FL"
)
# All-length
elem_all_yl1_alllength <- run_element_analysis(
    yl1annm_alllength_filtered, enrichment_lr_all, YOUNG_L1_SUBFAMILIES,
    results_all_yl1_alllength$reliable_sites, "Young L1",
    suffix = "_alllength", length_filter = "all"
)
elem_l1hs_alllength <- run_element_analysis(
    yl1annm_alllength_filtered, enrichment_lr_all, "L1HS",
    results_l1hs_alllength$reliable_sites, "L1HS",
    suffix = "_L1HS_alllength", length_filter = "all"
)


# ============================================================
# //ANCHOR - FINAL SUMMARY
# ============================================================
cat("\n=== Analysis Complete ===\n")
cat(sprintf("IVT false positive positions removed: %d\n", length(fp_harmonized_positions)))
cat(sprintf("Reliable m6A sites identified (FL): %d\n", nrow(reliable_sites)))
cat(sprintf("Reliable m6A sites (L1HS FL): %d\n", nrow(results_l1hs$reliable_sites)))
cat(sprintf("Reliable m6A sites (all-length): %d\n", nrow(results_all_yl1_alllength$reliable_sites)))
cat(sprintf("Reliable m6A sites (L1HS all-length): %d\n", nrow(results_l1hs_alllength$reliable_sites)))
cat(sprintf("Evolutionary categories:\n"))
print(table(site_evolution$category))
cat(sprintf("\nOutputs saved to:\n  Intermediates: %s\n  Plots: %s\n", intdir, outdir_base))
