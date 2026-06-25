############## Region-specific modification analysis
# Uses region-restricted modkit pileup outputs (higher sensitivity for specific RTE subfamilies)
cat(sprintf("\n=== Region-specific analysis for %s ===\n", current_mod_name))

regions <- if (!is.null(params$regions)) params$regions else conf$rte_subfamily_read_level_analysis
dir.create(sprintf("lrna/Rintermediates/%s/regions", current_mod_name), recursive = TRUE, showWarnings = FALSE)

for (region in regions) {
    cat(sprintf("  Processing region: %s\n", region))

    region_grs_list <- list()
    for (sample_name in samples) {
        filepath <- grep(region,
            grep(sprintf("/%s/", sample_name),
                inputs$region_bedmethylpaths,
                value = TRUE
            ),
            value = TRUE
        )
        if (length(filepath) == 0) next
        df <- read_table(filepath, col_names = FALSE)
        if (nrow(df) == 0) next
        gr <- GRanges(
            seqnames = df$X1,
            ranges = IRanges(start = df$X2, end = df$X2),
            cov = df$X10,
            pctM = as.double(df$X11)
        )
        gr$sample <- sample_name
        gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
        region_grs_list[[sample_name]] <- gr
        rm(df, gr)
    }
    if (length(region_grs_list) == 0) next

    region_grs <- Reduce(c, region_grs_list)
    rm(region_grs_list)

    # Filter by coverage
    region_grs <- region_grs[region_grs$cov > MINIMUMCOVERAGE]

    # Overlap with the specific RTE annotation bed
    eoi <- import(paste0("aref/default/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
    strand(eoi) <- "*"

    mbo <- mergeByOverlaps(region_grs, eoi)
    region_methdf <- mbo$region_grs %>%
        as.data.frame() %>%
        tibble()
    region_methdf$gene_id <- mbo$eoi %>%
        as.data.frame() %>%
        tibble() %$% name
    write_delim(region_methdf, sprintf("lrna/Rintermediates/%s/regions/%s_methdf.tsv", current_mod_name, region), col_names = TRUE)

    # Per-element summary
    region_perelementdf <- region_methdf %>%
        group_by(gene_id, sample, condition) %>%
        summarize(mean_meth = mean(pctM), n_sites = n(), .groups = "drop") %>%
        left_join(joinframe)
    region_perelementdf <- region_perelementdf %>% filter(!is.na(rte_length_req))
    write_delim(region_perelementdf, sprintf("lrna/Rintermediates/%s/regions/%s_perelementdf.tsv", current_mod_name, region), col_names = TRUE)

    rm(region_grs, region_methdf, region_perelementdf)
    gc()
}

cat(sprintf("\n=== bedmethylanalysis_process_new.R complete for %s ===\n", current_mod_name))




outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

# rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)

cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "A") %>%
    start() %>%
    unlist() %>%
    as.numeric()
consensus_ss[[1]][5762:5763]


cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

methdf <- region_methdf %>% left_join(joinframe)
mdf <- methdf %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end)))

senseelement <- mdf %>%
    filter(rte_strand == "+") %$% gene_id %>%
    pluck(1)

antisenseelement <- mdf %>%
    filter(rte_strand == "-") %$% gene_id %>%
    pluck(1)

cpgmapping_check <- cg_positions_df %>%
    filter(gene_id == senseelement) %$% sequence_pos %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == senseelement) %>%
    relocate(sequence_pos) %$% sequence_pos %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)
cpgmapping_check <- cg_positions_df %>%
    filter(gene_id == antisenseelement) %$% sequence_pos %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == antisenseelement) %>%
    relocate(sequence_pos) %$% sequence_pos %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)

merged <- left_join(cg_positions_df, mdf, by = c("gene_id", "sequence_pos"))
merged %>% filter(!is.na(seqnames))


library(ggrepel)

# detect DRACH motifs in consensus
# DRACH = [D=A/G/T][R=A/G][A][C][H=A/C/T], m6A target is the A at position 3
drach_matches <- vmatchPattern("DRACH", consensus_ss, fixed = FALSE)
# Get the A position (3rd base in each DRACH match)
drach_a_positions <- start(drach_matches[[1]]) + 2
cat(sprintf("Found %d DRACH motifs in L1HS consensus\n", length(drach_a_positions)))


p <- merged %>%
    filter(!is.na(seqnames)) %>%
    group_by(consensus_pos, condition) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    mutate(
        label = ifelse(mm > 50, as.character(consensus_pos), NA_character_),
        is_drach = consensus_pos %in% drach_a_positions
    ) %>%
    ggplot(aes(x = consensus_pos, y = mm)) +
    geom_point(data = . %>% filter(is_drach), color = "black", size = 2.5) +
    geom_point(aes(color = condition)) +
    facet_wrap(~condition) +
    geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/mod_along_element11111.pdf", w = 8, h = 4)

p <- merged %>%
    filter(!is.na(seqnames)) %>%
    mutate(motif = ifelse(consensus_pos %in% drach_a_positions, "DRACH", "non-DRACH")) %>%
    group_by(motif, condition) %>%
    summarise(
        mean_meth = mean(pctM, na.rm = TRUE),
        se_meth = sd(pctM, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = motif, y = mean_meth, fill = motif)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_meth - se_meth, ymax = mean_meth + se_meth), width = 0.2) +
    facet_wrap(~condition) +
    labs(x = NULL, y = "Mean % modified") +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/drach_vs_nondrach_barplot.pdf", w = 6, h = 4)


p <- merged %>%
    filter(!is.na(seqnames)) %>%
    filter(consensus_pos %in% drach_a_positions) %>%
    group_by(consensus_pos, condition) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    mutate(label = ifelse(mm > 50, as.character(consensus_pos), NA_character_)) %>%
    ggplot(aes(x = consensus_pos, y = mm)) +
    geom_point(aes(color = condition)) +
    facet_wrap(~condition) +
    geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/mod_along_element_drach_only1111.pdf", w = 8, h = 4)


top_pos <- merged %>%
    filter(!is.na(seqnames)) %>%
    filter(consensus_pos %in% drach_a_positions) %>%
    group_by(consensus_pos, condition) %>%
    summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    filter(mm > 30) %$% consensus_pos %>%
    unique()


library(tidyHeatmap)

# Heatmap of modification across all L1HS elements
heatmap_df <- merged %>%
    filter(!is.na(seqnames)) %>%
    group_by(gene_id, consensus_pos, condition) %>%
    summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>%
    mutate(consensus_pos = factor(consensus_pos, levels = sort(unique(consensus_pos))))

# Filter to elements with >= 50% position coverage to allow clustering
n_positions <- n_distinct(heatmap_df$consensus_pos)
well_covered_elements <- heatmap_df %>%
    group_by(gene_id, condition) %>%
    summarise(n_pos = n(), .groups = "drop") %>%
    filter(n_pos >= n_positions * 0.5) %>%
    pull(gene_id) %>%
    unique()
heatmap_df_filt <- heatmap_df %>% filter(gene_id %in% well_covered_elements)

p <- heatmap_df_filt %>%
    filter(condition == "N_ORF1") %>%
    heatmap(
        .row = gene_id,
        .column = consensus_pos,
        .value = mean_pctM,
        scale = "none",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE
    )
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/element_heatmap1.pdf", w = 10, h = 8)

# Heatmap filtered to DRACH motif positions only
heatmap_df_drach <- heatmap_df %>%
    filter(as.numeric(consensus_pos) %in% drach_a_positions)

n_drach_positions <- n_distinct(heatmap_df_drach$consensus_pos)
well_covered_drach <- heatmap_df_drach %>%
    group_by(gene_id, condition) %>%
    summarise(n_pos = n(), .groups = "drop") %>%
    filter(n_pos >= n_drach_positions * 0.5) %>%
    pull(gene_id) %>%
    unique()
heatmap_df_drach_filt <- heatmap_df_drach %>% filter(gene_id %in% well_covered_drach)

p <- heatmap_df_drach_filt %>%
    filter(condition == "N_ORF1") %>%
    heatmap(
        .row = gene_id,
        .column = consensus_pos,
        .value = mean_pctM,
        scale = "none",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        palette_value = circlize::colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
    )
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/element_heatmap_drach111.pdf", w = 14, h = 8)




readslist <- list()
for (region in conf$rte_subfamily_read_level_analysis) {
    for (sample_name in samples) {
        filepath <- grep(region,
            grep(sprintf("/%s/", sample_name),
                inputs$read_mods,
                value = TRUE
            ),
            value = TRUE
        )
        if (length(filepath) == 0) next
        if (file.size(filepath) == 0) next
        df <- read_delim(filepath)
        if (nrow(df) == 0) next
        df$region <- region
        df$sample <- sample_name
        df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grsx <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0("aref/default/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
        strand(eoi) <- "*"
        mbo <- mergeByOverlaps(grsx, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grsx"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grsx.", "", colnames(df1))
        readslist <- c(readslist, list(df1))
    }
}
if (length(readslist) > 0) {
    reads <- Reduce(rbind, readslist)
    write_delim(reads, sprintf("lrna/Rintermediates/%s/reads_context_all.tsv", current_mod_name), col_names = TRUE)
}


reads



readssmall <- reads %>% head(n = 100000)



# cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

mdfreadsfull <- reads %>%
    left_join(joinframe) %>%
    filter(strand == rte_strand) %>%
    mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end)))

mdfreads <- mdfreadsfull %>% dplyr::select(gene_id, sample, condition, sequence_pos, read_id, mod_qual)

mergedreads <- left_join(cg_positions_df, mdfreads, by = c("gene_id", "sequence_pos")) %>%
    filter(!is.na(seqnames))
mergedreads %>% filter(!is.na(seqnames))


p <- mergedreads %>%
    filter(!is.na(seqnames)) %>%
    dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0)) %>%
    group_by(consensus_pos, condition) %>%
    summarise(mm = mean(mod_indicator, na.rm = TRUE), .groups = "drop") %>%
    mutate(
        is_drach = consensus_pos %in% drach_a_positions
    ) %>%
    ggplot(aes(x = consensus_pos, y = mm)) +
    geom_point(data = . %>% filter(is_drach), color = "black", size = 2.5) +
    geom_point(aes(color = condition)) +
    facet_wrap(~condition) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/mod_along_element111111.pdf", w = 8, h = 4)


mergedreads %$% read_id %>% unique()


mergedreads %>% filter(consensus_pos %in% top_pos)

# Per-read heatmap at top DRACH positions
read_heatmap_df <- mergedreads %>%
    filter(!is.na(seqnames)) %>%
    filter(consensus_pos %in% top_pos) %>%
    mutate(
        mod_indicator = ifelse(mod_qual > 0.5, 1, 0),
        consensus_pos = factor(consensus_pos, levels = sort(unique(consensus_pos)))
    ) %>%
    dplyr::select(read_id, consensus_pos, mod_indicator, sample, condition) %>%
    distinct()

# Keep reads that cover at least 50% of top positions
n_top <- n_distinct(read_heatmap_df$consensus_pos)
good_reads <- read_heatmap_df %>%
    group_by(read_id) %>%
    summarise(n_pos = n(), .groups = "drop") %>%
    filter(n_pos >= n_top * 0.5) %>%
    pull(read_id)
read_heatmap_df_filt <- read_heatmap_df %>% filter(read_id %in% good_reads)

# p <- read_heatmap_df_filt %>%
#     heatmap(
#         .row = read_id,
#         .column = consensus_pos,
#         .value = mod_indicator,
#         scale = "none",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = FALSE,
#         palette_value = circlize::colorRamp2(c(0, 1), c("white", "red"))
#     )
# mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/per_read_heatmap_top_drach.pdf", w = 10, h = 10)

# # Scalable alternatives for within-read modification correlation

# 1. Pairwise position correlation matrix
read_mod_wide <- read_heatmap_df_filt %>%
    pivot_wider(id_cols = c(read_id, condition), names_from = consensus_pos, values_from = mod_indicator)
mod_mat <- read_mod_wide %>%
    dplyr::select(-read_id, -condition) %>%
    as.matrix()
cor_mat <- cor(mod_mat, use = "pairwise.complete.obs")

p <- cor_mat %>%
    as.data.frame() %>%
    rownames_to_column("pos1") %>%
    pivot_longer(-pos1, names_to = "pos2", values_to = "correlation") %>%
    mutate(
        pos1 = factor(pos1, levels = sort(as.numeric(unique(pos1)))),
        pos2 = factor(pos2, levels = sort(as.numeric(unique(pos2))))
    ) %>%
    ggplot(aes(x = pos1, y = pos2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    labs(x = "Consensus position", y = "Consensus position", title = "Pairwise modification correlation across reads") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/position_correlation_matrix1.pdf", w = 8, h = 7)

# 2. Per-read modification count distribution
p <- read_heatmap_df_filt %>%
    group_by(read_id, condition) %>%
    summarise(n_modified = sum(mod_indicator), n_total = n(), .groups = "drop") %>%
    mutate(frac_modified = n_modified / n_total) %>%
    ggplot(aes(x = frac_modified, fill = condition)) +
    geom_histogram(bins = 20, position = "dodge") +
    labs(x = "Fraction of top DRACH positions modified per read", y = "Number of reads") +
    facet_wrap(~condition) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/per_read_frac_modified_hist1.pdf", w = 8, h = 4)

# 3. Co-modification odds ratios for all position pairs
positions <- levels(read_heatmap_df_filt$consensus_pos)
comod_results <- list()
for (i in seq_along(positions)) {
    for (j in seq_along(positions)) {
        if (j <= i) next
        p1 <- positions[i]
        p2 <- positions[j]
        both <- read_mod_wide %>% filter(!is.na(!!sym(p1)) & !is.na(!!sym(p2)))
        a <- sum(both[[p1]] == 1 & both[[p2]] == 1)
        b <- sum(both[[p1]] == 1 & both[[p2]] == 0)
        cc <- sum(both[[p1]] == 0 & both[[p2]] == 1)
        d <- sum(both[[p1]] == 0 & both[[p2]] == 0)
        or_val <- ifelse((b * cc) > 0, (a * d) / (b * cc), NA)
        comod_results <- c(comod_results, list(tibble(pos1 = p1, pos2 = p2, odds_ratio = or_val, n = nrow(both))))
    }
}
comod_df <- bind_rows(comod_results)

p <- comod_df %>%
    mutate(log2_OR = log2(odds_ratio)) %>%
    mutate(
        pos1 = factor(pos1, levels = sort(as.numeric(unique(c(pos1, pos2))))),
        pos2 = factor(pos2, levels = sort(as.numeric(unique(c(pos1, pos2)))))
    ) %>%
    ggplot(aes(x = pos1, y = pos2, fill = log2_OR)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2(OR)") +
    labs(x = "Consensus position", y = "Consensus position", title = "Co-modification odds ratios between DRACH positions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/comodification_odds_ratios1.pdf", w = 8, h = 7)

# 4. PCA of per-read modification patterns
pca_mat <- mod_mat
pca_mat[is.na(pca_mat)] <- 0
# Remove zero-variance columns
pca_mat <- pca_mat[, apply(pca_mat, 2, var) > 0, drop = FALSE]
pca_res <- prcomp(pca_mat, center = TRUE, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:min(3, ncol(pca_res$x))]) %>%
    mutate(
        read_id = read_mod_wide$read_id,
        condition = read_mod_wide$condition
    )
var_explained <- summary(pca_res)$importance[2, 1:2] * 100

p <- pca_df %>%
    ggplot(aes(x = PC1, y = PC2, color = condition)) +
    geom_point(alpha = 0.3, size = 0.5) +
    labs(
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC2 (%.1f%%)", var_explained[2]),
        title = "PCA of per-read modification at top DRACH positions"
    ) +
    mtopen
mysaveandstore("lrna/results/plots/modification_analysis/rtes/l1hs/reads/per_read_pca1.pdf", w = 7, h = 5)
