if (interactive()) {
    module_name <- "srna"
} else {
    module_name <- snakemake@params$module_name
}
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)
set.seed(123)

library(magrittr)
library(stringr)
library(cowplot)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(dplyr)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
library(vcfR)

# analysis parameters
{
    tryCatch(
        {
            params <- snakemake@params
            inputs <- snakemake@input
            outputs <- snakemake@output
        },
        error = function(e) {
            assign("params", list(
                "outputdir" = "srna/results/agg/loi_analysis",
                "imprinted_genes" = "resources/imprinted_genes.csv",
                "annotation_genes" = conf$annotation_genes
            ), env = globalenv())
            assign("inputs", list(
                ase_counts = paste0("srna/outs/", conf$samples, "/ase/ase_counts.table"),
                phased_vcf = paste0(
                    "aref/intermediates/A.REF/clair3/", confALL$aref$rate, "/",
                    confALL$aref$type, ".", confALL$aref$modification_string, "/phased_hets.filtered.vcf.gz"
                )
            ), env = globalenv())
            assign("outputs", list(
                plots = "srna/results/agg/loi_analysis/loi_plots.RData"
            ), env = globalenv())
        }
    )
}

outputdir <- params[["outputdir"]]
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

# load gene annotation GTF for SNP-to-gene mapping (exons only to exclude intronic signal)
refseq_gr <- import(confALL$ldna$refseq_unaltered)
exons_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "exon", ]
exons_gr <- exons_gr[mcols(exons_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(exons_gr)$gene_id <- mcols(exons_gr)$gene
# collapse overlapping exons per gene to avoid double-counting
exons_gr <- reduce(split(exons_gr, exons_gr$gene_id))
exons_gr <- unlist(exons_gr)
mcols(exons_gr) <- DataFrame(gene_id = names(exons_gr))
names(exons_gr) <- NULL


# load imprinted gene list
imprinted <- read_csv("/users/mkelsey/data/LF1_newnanoporeset/TE-Seq/custom/resources/akbarisup7_imprintedgenes.csv")
# expect at minimum a column called "gene" with gene symbols
imprinted_genes <- unique(imprinted$gene_id)

# load ASE count tables
ase_list <- list()
for (i in seq_along(inputs$ase_counts)) {
    fn <- inputs$ase_counts[[i]]
    sample_name <- basename(dirname(dirname(fn)))
    df <- read_delim(fn, delim = "\t", show_col_types = FALSE)
    df$sample_name <- sample_name
    ase_list[[sample_name]] <- df
}
ase <- bind_rows(ase_list)

# GATK ASEReadCounter columns: contig, position, variantID, refAllele, altAllele,
# refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs
ase <- ase %>%
    filter(totalCount >= 10) %>%
    mutate(
        ref_ratio = refCount / totalCount,
        alt_ratio = altCount / totalCount
    )

# load phased VCF to get GT orientation and phase sets
# we need to know 0|1 vs 1|0 to orient ref/alt to haplotype 1/2
phased_vcf_path <- inputs$phased_vcf
if (is.list(phased_vcf_path)) phased_vcf_path <- phased_vcf_path[[1]]
vcf <- vcfR::read.vcfR(phased_vcf_path, verbose = FALSE)
gt_mat <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
ps_mat <- vcfR::extract.gt(vcf, element = "PS", as.numeric = FALSE)
phase_df <- tibble(
    contig = vcfR::getCHROM(vcf),
    position = as.integer(vcfR::getPOS(vcf)),
    GT = gt_mat[, 1],
    PS = ps_mat[, 1]
) %>%
    mutate(
        phased = grepl("\\|", GT),
        ref_is_hap1 = GT == "0|1"
    )
rm(vcf, gt_mat, ps_mat)

# map SNPs to exons to get gene assignment (excludes intronic SNPs)
snp_gr <- GRanges(seqnames = ase$contig, ranges = IRanges(start = ase$position, end = ase$position))
hits <- findOverlaps(snp_gr, exons_gr)
ase$gene_id <- NA_character_
ase$gene_id[queryHits(hits)] <- exons_gr$gene_id[subjectHits(hits)]

# remove SNPs not in any exon
ase <- ase %>% filter(!is.na(gene_id))

# add condition from sample table
ase <- ase %>%
    left_join(sample_table %>% select(sample_name, condition))

# join phase info and orient counts to haplotypes
ase <- ase %>%
    left_join(phase_df %>% select(contig, position, GT, PS, phased, ref_is_hap1),
        by = c("contig", "position")
    )
# for phased sites: assign counts to hap1/hap2 based on GT orientation
# for unphased sites (/ genotype): fall back to ref/alt (single-SNP genes still work)
ase <- ase %>%
    mutate(
        phased = ifelse(is.na(phased), FALSE, phased),
        ref_is_hap1 = ifelse(is.na(ref_is_hap1), TRUE, ref_is_hap1),
        PS = ifelse(is.na(PS) | !phased, paste0(contig, ":", position), PS),
        hap1Count = ifelse(ref_is_hap1, refCount, altCount),
        hap2Count = ifelse(ref_is_hap1, altCount, refCount)
    )


ase %>% filter(gene_id == "H19")
# aggregate allelic counts per gene per sample per phase block
# within a phase block, hap1/hap2 assignments are consistent
gene_ase <- ase %>%
    group_by(sample_name, condition, gene_id, PS) %>%
    summarise(
        hap1Count = sum(hap1Count),
        hap2Count = sum(hap2Count),
        totalCount = sum(totalCount),
        n_snps = n(),
        .groups = "drop"
    ) %>%
    # for genes spanning multiple phase blocks, keep the block with most SNPs
    group_by(sample_name, condition, gene_id) %>%
    slice_max(n_snps, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(hap1_ratio = hap1Count / totalCount) %>%
    filter(!(condition %in% c("AP", "AE", "AL")))

gene_ase %>% filter(gene_id == "H19")
# determine which haplotype is the minor one ONCE per gene, using pooled counts across all samples
# this ensures consistent orientation when comparing conditions
gene_hap_orientation <- gene_ase %>%
    group_by(gene_id) %>%
    summarise(global_hap1_ratio = sum(hap1Count) / sum(totalCount), .groups = "drop") %>%
    mutate(hap1_is_minor = global_hap1_ratio < 0.5)
gene_ase <- gene_ase %>%
    left_join(gene_hap_orientation, by = "gene_id") %>%
    mutate(
        minor_count = ifelse(hap1_is_minor, hap1Count, hap2Count),
        major_count = ifelse(hap1_is_minor, hap2Count, hap1Count),
        minor_allele_ratio = minor_count / totalCount
    )

# flag imprinted genes and chrX genes
chrX_genes <- unique(exons_gr$gene_id[as.character(seqnames(exons_gr)) == "chrX"])
gene_ase <- gene_ase %>%
    mutate(
        imprinted = gene_id %in% imprinted_genes,
        chrX = gene_id %in% chrX_genes
    )

# binomial test for biallelic expression per gene per sample
# under monoallelic imprinting, expect ratio ~0 or ~1
# under LOI, expect ratio ~0.5
gene_ase <- gene_ase %>%
    rowwise() %>%
    mutate(
        binom_p = binom.test(hap1Count, totalCount, p = 0.5)$p.value
    ) %>%
    ungroup()

# classify allelic expression
gene_ase <- gene_ase %>%
    mutate(allelic_class = case_when(
        minor_allele_ratio >= 0.35 ~ "biallelic",
        minor_allele_ratio >= 0.1 ~ "skewed",
        TRUE ~ "monoallelic"
    ))

# summary per gene across conditions
# pool counts across samples so high-expression SNPs (exonic) dominate over low-expression (intronic)
gene_condition_summary <- gene_ase %>%
    group_by(gene_id, condition, imprinted) %>%
    summarise(
        pooled_minor = sum(minor_count),
        pooled_major = sum(major_count),
        pooled_total = sum(totalCount),
        pooled_minor_ratio = pooled_minor / pooled_total,
        mean_minor_ratio = mean(minor_allele_ratio),
        sd_minor_ratio = sd(minor_allele_ratio),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    filter(n_samples == (sample_table %>% group_by(condition) %>% mutate(n = n()) %$% n %>% max())) %>%
    filter(pooled_total > 200) %>%
    group_by(gene_id) %>%
    mutate(n = n()) %>%
    filter(n == length(unique(sample_table$condition)))




# test for LOI: compare minor allele ratio between conditions
# for imprinted genes, test whether SEN has higher minor allele ratio than PRO
loi_test <- gene_ase %>%
    filter(imprinted) %>%
    group_by(gene_id) %>%
    filter(n_distinct(condition) == 2) %>%
    summarise(
        mean_ratio_PRO = mean(minor_allele_ratio[condition == conf$levels[1]]),
        mean_ratio_SEN = mean(minor_allele_ratio[condition == conf$levels[2]]),
        delta_ratio = mean_ratio_SEN - mean_ratio_PRO,
        wilcox_p = tryCatch(
            wilcox.test(
                minor_allele_ratio[condition == conf$levels[2]],
                minor_allele_ratio[condition == conf$levels[1]]
            )$p.value,
            error = function(e) NA_real_
        ),
        n_PRO = sum(condition == conf$levels[1]),
        n_SEN = sum(condition == conf$levels[2]),
        .groups = "drop"
    ) %>%
    mutate(padj = p.adjust(wilcox_p, method = "BH"))

# save tables
write_tsv(gene_ase, file.path(outputdir, "gene_ase_all.tsv"))
write_tsv(gene_ase %>% filter(imprinted), file.path(outputdir, "gene_ase_imprinted.tsv"))
write_tsv(loi_test, file.path(outputdir, "loi_test_results.tsv"))
write_tsv(gene_condition_summary, file.path(outputdir, "gene_condition_summary.tsv"))

######################
# PLOTS
######################

# restrict all plots to genes that passed the gene_condition_summary filters
passing_genes <- unique(gene_condition_summary$gene_id)
gene_ase_plot <- gene_ase %>% filter(gene_id %in% passing_genes)

# 1. Allelic ratio distribution for imprinted genes by condition
imprinted_ase <- gene_ase_plot %>% filter(imprinted)

if (nrow(imprinted_ase) > 0) {
    p <- ggplot(imprinted_ase, aes(x = condition, y = minor_allele_ratio, fill = condition)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
        facet_wrap(~gene_id, scales = "free_y") +
        geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = 0.1, linetype = "dotted", alpha = 0.3) +
        labs(
            title = "Allelic Ratios at Imprinted Genes",
            subtitle = "Minor allele ratio (0 = monoallelic, 0.5 = biallelic)",
            y = "Minor Allele Ratio",
            x = NULL
        ) +
        mtclosed
    mysaveandstore(file.path(outputdir, "imprinted_allelic_ratios_boxplot.pdf"),
        w = max(6, 3 * ceiling(sqrt(length(unique(imprinted_ase$gene_id))))),
        h = max(5, 3 * ceiling(length(unique(imprinted_ase$gene_id)) / ceiling(sqrt(length(unique(imprinted_ase$gene_id)))))),
        pl = p
    )

    # 2. Heatmap of minor allele ratios at imprinted genes
    ratio_mat <- imprinted_ase %>%
        dplyr::select(sample_name, gene_id, minor_allele_ratio) %>%
        pivot_wider(names_from = sample_name, values_from = minor_allele_ratio) %>%
        tibble::column_to_rownames("gene_id") %>%
        as.matrix()
    ratio_mat <- ratio_mat[complete.cases(ratio_mat), , drop = FALSE]

    if (nrow(ratio_mat) > 1 && ncol(ratio_mat) > 1) {
        sample_cond <- sample_table %>%
            filter(sample_name %in% colnames(ratio_mat)) %>%
            arrange(match(sample_name, colnames(ratio_mat)))
        ha <- ComplexHeatmap::HeatmapAnnotation(
            Condition = sample_cond$condition
        )
        p <- ComplexHeatmap::Heatmap(
            ratio_mat,
            name = "Minor\nAllele\nRatio",
            top_annotation = ha,
            cluster_rows = nrow(ratio_mat) > 2,
            cluster_columns = FALSE,
            column_split = sample_cond$condition,
            col = circlize::colorRamp2(c(0, 0.25, 0.5), c("navy", "white", "firebrick")),
            row_names_gp = grid::gpar(fontsize = 10),
            column_names_gp = grid::gpar(fontsize = 10),
            column_title = "Imprinted Gene Allelic Ratios"
        )
        mysaveandstore(file.path(outputdir, "imprinted_allelic_ratios_heatmap.pdf"),
            w = max(6, ncol(ratio_mat) * 0.8 + 3),
            h = max(4, nrow(ratio_mat) * 0.4 + 3),
            pl = p
        )
    }

    # 3. LOI volcano-style plot: delta ratio vs significance
    loi_test_plot <- loi_test %>% filter(gene_id %in% passing_genes)
    if (nrow(loi_test_plot) > 0 && any(!is.na(loi_test_plot$wilcox_p))) {
        p <- loi_test_plot %>%
            mutate(neg_log10_p = -log10(wilcox_p)) %>%
            ggplot(aes(x = delta_ratio, y = neg_log10_p)) +
            geom_point(aes(color = padj < 0.1), size = 3) +
            geom_text_repel(aes(label = gene_id), max.overlaps = 20, size = 3) +
            geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
            geom_hline(yintercept = -log10(0.05), linetype = "dotted", alpha = 0.3) +
            labs(
                title = "Loss of Imprinting Test",
                subtitle = sprintf("%s vs %s", conf$levels[2], conf$levels[1]),
                x = sprintf("Delta Minor Allele Ratio (%s - %s)", conf$levels[2], conf$levels[1]),
                y = "-log10(p-value)",
                color = "padj < 0.1"
            ) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
            mtclosed
        mysaveandstore(file.path(outputdir, "loi_volcano.pdf"), w = 7, h = 6, pl = p)
    }

    # 4. Per-gene barplot of allelic counts for top imprinted genes
    top_genes <- imprinted_ase %>%
        group_by(gene_id) %>%
        summarise(mean_total = mean(totalCount), .groups = "drop") %>%
        arrange(desc(mean_total)) %>%
        head(1000) %>%
        pull(gene_id)

    if (length(top_genes) > 0) {
        plot_df <- imprinted_ase %>%
            filter(gene_id %in% top_genes) %>%
            dplyr::select(sample_name, condition, gene_id, hap1Count, hap2Count) %>%
            pivot_longer(cols = c(hap1Count, hap2Count), names_to = "haplotype", values_to = "count") %>%
            mutate(haplotype = gsub("Count", "", haplotype))

        p <- ggplot(plot_df, aes(x = sample_name, y = count, fill = haplotype)) +
            geom_bar(stat = "identity", position = "stack") +
            facet_wrap(~gene_id, scales = "free_y") +
            labs(
                title = "Haplotype-Specific Read Counts at Imprinted Genes",
                y = "Read Count",
                x = NULL,
                fill = "Haplotype"
            ) +
            mtclosed +
            scale_fill_manual(values = c("hap1" = "#4393C3", "hap2" = "#D6604D")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(file.path(outputdir, "imprinted_allelic_counts_barplot.pdf"),
            w = max(8, length(top_genes) * 2), h = max(6, ceiling(length(top_genes) / 4) * 3), pl = p
        )
    }
}

# 5. Genome-wide allelic ratio scatter: PRO vs SEN mean minor allele ratio
gw_summary <- gene_ase_plot %>%
    group_by(gene_id, condition, imprinted) %>%
    summarise(mean_ratio = mean(minor_allele_ratio), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_ratio)

level_pairs <- combn(conf$levels, 2, simplify = FALSE)
for (pair in level_pairs) {
    if (all(pair %in% colnames(gw_summary))) {
        p <- ggplot(gw_summary, aes(x = !!sym(pair[1]), y = !!sym(pair[2]))) +
            geom_abline(intercept = 0, slope = 1, alpha = 0.4) +
            geom_point(aes(color = imprinted), size = 1.5, alpha = 0.6) +
            geom_text_repel(
                data = gw_summary %>% filter(imprinted),
                aes(label = gene_id), size = 3, max.overlaps = 20
            ) +
            labs(
                title = "Genome-wide Allelic Ratios",
                x = sprintf("%s Mean Minor Allele Ratio", pair[1]),
                y = sprintf("%s Mean Minor Allele Ratio", pair[2]),
                color = "Imprinted"
            ) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
            mtclosed +
            theme(aspect.ratio = 1)
        mysaveandstore(file.path(outputdir, sprintf("genomewide_allelic_ratio_scatter_%s_vs_%s.pdf", pair[1], pair[2])), w = 7, h = 7, pl = p)
    }
}

######################
# XCI EROSION ANALYSIS
######################
xci_outputdir <- file.path(outputdir, "xci")
dir.create(xci_outputdir, showWarnings = FALSE, recursive = TRUE)

# condition summary for chrX genes (same filters as LOI)
xci_condition_summary <- gene_ase %>%
    filter(chrX) %>%
    group_by(gene_id, condition) %>%
    summarise(
        pooled_minor = sum(minor_count),
        pooled_major = sum(major_count),
        pooled_total = sum(totalCount),
        pooled_minor_ratio = pooled_minor / pooled_total,
        mean_minor_ratio = mean(minor_allele_ratio),
        sd_minor_ratio = sd(minor_allele_ratio),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    filter(n_samples == (sample_table %>% group_by(condition) %>% mutate(n = n()) %$% n %>% max())) %>%
    filter(pooled_total > 200) %>%
    group_by(gene_id) %>%
    mutate(n = n()) %>%
    filter(n == length(unique(sample_table$condition)))


xci_passing_genes <- unique(xci_condition_summary$gene_id)
write_tsv(xci_condition_summary, file.path(xci_outputdir, "xci_condition_summary.tsv"))

# test for XCI erosion: compare minor allele ratio between conditions
xci_test <- gene_ase %>%
    filter(chrX, gene_id %in% xci_passing_genes) %>%
    group_by(gene_id) %>%
    filter(n_distinct(condition) == 2) %>%
    summarise(
        mean_ratio_PRO = mean(minor_allele_ratio[condition == conf$levels[1]]),
        mean_ratio_SEN = mean(minor_allele_ratio[condition == conf$levels[2]]),
        delta_ratio = mean_ratio_SEN - mean_ratio_PRO,
        wilcox_p = tryCatch(
            wilcox.test(
                minor_allele_ratio[condition == conf$levels[2]],
                minor_allele_ratio[condition == conf$levels[1]]
            )$p.value,
            error = function(e) NA_real_
        ),
        n_PRO = sum(condition == conf$levels[1]),
        n_SEN = sum(condition == conf$levels[2]),
        .groups = "drop"
    ) %>%
    mutate(padj = p.adjust(wilcox_p, method = "BH"))
write_tsv(xci_test, file.path(xci_outputdir, "xci_test_results.tsv"))

# XCI plots
xci_ase <- gene_ase %>% filter(chrX, gene_id %in% xci_passing_genes)

if (nrow(xci_ase) > 0) {
    # 1. Distribution of minor allele ratios for all chrX genes by condition
    p <- ggplot(xci_ase, aes(x = minor_allele_ratio, fill = condition)) +
        geom_density(alpha = 0.5) +
        geom_vline(xintercept = 0.1, linetype = "dotted", alpha = 0.3) +
        geom_vline(xintercept = 0.35, linetype = "dashed", alpha = 0.3) +
        labs(
            title = "X Chromosome Allelic Ratio Distribution",
            subtitle = "Shift toward 0.5 in SEN suggests XCI erosion",
            x = "Minor Allele Ratio",
            y = "Density"
        ) +
        mtclosed
    mysaveandstore(file.path(xci_outputdir, "xci_ratio_density.pdf"), w = 7, h = 5, pl = p)

    # 2. Paired boxplot: per-gene minor allele ratio PRO vs SEN
    p <- ggplot(xci_ase, aes(x = condition, y = minor_allele_ratio, fill = condition)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
        labs(
            title = "X Chromosome Allelic Ratios by Condition",
            y = "Minor Allele Ratio",
            x = NULL
        ) +
        mtclosed
    mysaveandstore(file.path(xci_outputdir, "xci_ratio_boxplot.pdf"), w = 5, h = 5, pl = p)

    # 3. Heatmap of chrX minor allele ratios
    xci_ratio_mat <- xci_ase %>%
        dplyr::select(sample_name, gene_id, minor_allele_ratio) %>%
        pivot_wider(names_from = sample_name, values_from = minor_allele_ratio) %>%
        tibble::column_to_rownames("gene_id") %>%
        as.matrix()
    xci_ratio_mat <- xci_ratio_mat[complete.cases(xci_ratio_mat), , drop = FALSE]

    if (nrow(xci_ratio_mat) > 1 && ncol(xci_ratio_mat) > 1) {
        sample_cond <- sample_table %>%
            filter(sample_name %in% colnames(xci_ratio_mat)) %>%
            arrange(match(sample_name, colnames(xci_ratio_mat)))
        ha <- ComplexHeatmap::HeatmapAnnotation(
            Condition = sample_cond$condition
        )
        p <- ComplexHeatmap::Heatmap(
            xci_ratio_mat,
            name = "Minor\nAllele\nRatio",
            top_annotation = ha,
            cluster_rows = nrow(xci_ratio_mat) > 2,
            cluster_columns = FALSE,
            column_split = sample_cond$condition,
            col = circlize::colorRamp2(c(0, 0.25, 0.5), c("navy", "white", "firebrick")),
            row_names_gp = grid::gpar(fontsize = if (nrow(xci_ratio_mat) > 50) 4 else 8),
            column_names_gp = grid::gpar(fontsize = 10),
            column_title = "X Chromosome Allelic Ratios",
            show_row_names = nrow(xci_ratio_mat) <= 100
        )
        mysaveandstore(file.path(xci_outputdir, "xci_allelic_ratios_heatmap.pdf"),
            w = max(6, ncol(xci_ratio_mat) * 0.8 + 3),
            h = max(4, nrow(xci_ratio_mat) * 0.25 + 3),
            pl = p
        )
    }

    # 4. Volcano plot: XCI erosion
    if (nrow(xci_test) > 0 && any(!is.na(xci_test$wilcox_p))) {
        p <- xci_test %>%
            mutate(neg_log10_p = -log10(wilcox_p)) %>%
            ggplot(aes(x = delta_ratio, y = neg_log10_p)) +
            geom_point(aes(color = padj < 0.1), size = 2) +
            geom_text_repel(
                data = . %>% filter(padj < 0.1 | abs(delta_ratio) > 0.1),
                aes(label = gene_id), max.overlaps = 30, size = 3
            ) +
            geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
            geom_hline(yintercept = -log10(0.05), linetype = "dotted", alpha = 0.3) +
            labs(
                title = "XCI Erosion Test (X Chromosome Genes)",
                subtitle = sprintf("%s vs %s", conf$levels[2], conf$levels[1]),
                x = sprintf("Delta Minor Allele Ratio (%s - %s)", conf$levels[2], conf$levels[1]),
                y = "-log10(p-value)",
                color = "padj < 0.1"
            ) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
            mtclosed
        mysaveandstore(file.path(xci_outputdir, "xci_volcano.pdf"), w = 7, h = 6, pl = p)
    }

    # 5. Scatter: PRO vs SEN minor allele ratio for chrX genes
    xci_scatter <- xci_condition_summary %>%
        dplyr::select(gene_id, condition, pooled_minor_ratio) %>%
        pivot_wider(names_from = condition, values_from = pooled_minor_ratio)

    if (all(conf$levels %in% colnames(xci_scatter))) {
        p <- ggplot(xci_scatter, aes(x = !!sym(conf$levels[1]), y = !!sym(conf$levels[2]))) +
            geom_abline(intercept = 0, slope = 1, alpha = 0.4) +
            geom_point(size = 2, alpha = 0.6) +
            geom_text_repel(aes(label = gene_id), size = 2.5, max.overlaps = 20) +
            labs(
                title = "X Chromosome: Pooled Allelic Ratios",
                x = sprintf("%s Pooled Minor Allele Ratio", conf$levels[1]),
                y = sprintf("%s Pooled Minor Allele Ratio", conf$levels[2])
            ) +
            mtclosed +
            theme(aspect.ratio = 1) +
            coord_fixed(xlim = c(0, 0.5), ylim = c(0, 0.5))
        mysaveandstore(file.path(xci_outputdir, "xci_scatter_pro_vs_sen.pdf"), w = 7, h = 7, pl = p)
    }

    # 6. Per-gene barplot of haplotype counts for top XCI erosion candidates
    top_xci <- xci_test %>%
        filter(!is.na(delta_ratio)) %>%
        arrange(desc(delta_ratio)) %>%
        head(1000) %>%
        pull(gene_id)

    if (length(top_xci) > 0) {
        plot_df <- xci_ase %>%
            filter(gene_id %in% top_xci) %>%
            dplyr::select(sample_name, condition, gene_id, minor_count, major_count) %>%
            pivot_longer(cols = c(minor_count, major_count), names_to = "allele", values_to = "count") %>%
            mutate(allele = ifelse(allele == "minor_count", "Minor", "Major"))

        p <- ggplot(plot_df, aes(x = sample_name, y = count, fill = allele)) +
            geom_bar(stat = "identity", position = "stack") +
            facet_wrap(~gene_id, scales = "free_y") +
            labs(
                title = "Top XCI Erosion Candidates: Allele-Specific Counts",
                y = "Read Count",
                x = NULL,
                fill = "Allele"
            ) +
            mtclosed +
            scale_fill_manual(values = c("Minor" = "#D6604D", "Major" = "#4393C3")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(file.path(xci_outputdir, "xci_top_candidates_barplot2.pdf"),
            w = max(8, length(top_xci) * 2), h = max(6, ceiling(length(top_xci) / 4) * 3), pl = p
        )
    }

    # 7. Genes in the same phase set as XIST — hap1/hap2 coloring to check inversion
    # XIST is expressed from Xi; other X genes from Xa — expect opposite haplotype dominance
    xist_ps <- gene_ase %>%
        filter(gene_id == "XIST", chrX) %>%
        pull(PS) %>%
        unique()

    if (length(xist_ps) > 0) {
        xist_block_genes <- gene_ase %>%
            filter(chrX, PS %in% xist_ps, gene_id %in% xci_passing_genes) %>%
            pull(gene_id) %>%
            unique()

        if (length(xist_block_genes) > 1) {
            plot_df <- gene_ase %>%
                filter(gene_id %in% xist_block_genes) %>%
                dplyr::select(sample_name, condition, gene_id, hap1Count, hap2Count) %>%
                pivot_longer(cols = c(hap1Count, hap2Count), names_to = "haplotype", values_to = "count") %>%
                mutate(
                    haplotype = gsub("Count", "", haplotype),
                    gene_id = factor(gene_id, levels = c("XIST", setdiff(xist_block_genes, "XIST")))
                )

            p <- ggplot(plot_df, aes(x = sample_name, y = count, fill = haplotype)) +
                geom_bar(stat = "identity", position = "stack") +
                facet_wrap(~gene_id, scales = "free_y") +
                labs(
                    title = "X Genes in XIST Phase Block: Haplotype Counts",
                    subtitle = "XIST (from Xi) should show opposite haplotype dominance vs other genes (from Xa)",
                    y = "Read Count",
                    x = NULL,
                    fill = "Haplotype"
                ) +
                mtclosed +
                scale_fill_manual(values = c("hap1" = "#4393C3", "hap2" = "#D6604D")) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            mysaveandstore(file.path(xci_outputdir, "xist_phaseblock_haplotype_barplot1.pdf"),
                w = max(8, length(xist_block_genes) * 2),
                h = max(6, ceiling(length(xist_block_genes) / 4) * 3),
                pl = p
            )
        }
    }
}

# save all plots
if (exists("mysaveandstoreplots")) {
    save(mysaveandstoreplots, file = outputs$plots)
}
