module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>% mutate(sample = sample_name)
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(zoo)
library(pryr)
library(circlize)
library(rGREAT)
library(reactome.db)
library(msigdb)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)
# library(ReMapEnrich)
library(msigdbr)
library(Biostrings)
library(ggpubr)
library(PCAtools)
library(patchwork)
library(ggh4x)
library(ggnewscale)
library(betareg)
library(scales)
library(ggnewscale)

samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]
conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name
enough_samples_per_condition_for_stats <- ifelse(length(condition1samples) < 3 | length(condition2samples) < 3, FALSE, TRUE)

adjustment_set <- c(conf$linear_model_adjustment_set)
lm_right_hand_side <- ifelse(is.null(adjustment_set), "condition", paste0(c("condition", adjustment_set), collapse = " + "))
adjustment_set_categorical <- if (is.null(adjustment_set)) {
    NULL
} else {
    sample_table %>%
        dplyr::select(all_of(adjustment_set)) %>%
        map(class) %>%
        unlist() %>%
        grep(pattern = "character|factor", value = TRUE) %>%
        names()
}

asc <- ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1])


{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")
    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
    if ("chrY" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, grep("_chrX_|_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes)
        } else {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, grep("_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX)
        }
    } else if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrY, grep("_chrX_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrY)
    } else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
    }
}
#################### functions and themes
named_group_split <- function(.tbl, ...) {
    grouped <- group_by(.tbl, ...)
    names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

    grouped %>%
        group_split() %>%
        rlang::set_names(names)
}
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", sample_table$sample_name, sample_table$sample_name),
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/m/tables/dmrs.tsv",
            dmls = "ldna/results/m/tables/dmls.tsv"
        ), env = globalenv())
        assign("params", list(
            dmrs = "ldna/results/m/tables/dmrs.tsv",
            dmls = "ldna/results/m/tables/dmls.tsv",
            mod_code = "m"
        ), env = globalenv())
        assign("outputs", list(
            promoters_bed = "ldna/Rintermediates/m/promoters_t05.bed",
            dmrpromoterhyper_bed = "ldna/Rintermediates/m/promoters_dmhyperregions_t05.bed",
            dmrpromoterhypo_bed = "ldna/Rintermediates/m/promoters_dmhyporegions_t05.bed"
        ), env = globalenv())
    }
)


merge_with_grs <- function(grs, rte_frame) {
    mbo <- mergeByOverlaps(grs, rte_frame)
    methdf <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    rte_only_frame <- mbo$rte_frame %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
    rtedf_promoters <- bind_cols(methdf, rte_only_frame)
    return(rtedf_promoters)
}

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)


flRTEpromoter <- read_delim(sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)
RMdf <- read_delim(sprintf("ldna/Rintermediates/%s/RMdf.tsv", params$mod_code), col_names = TRUE)

rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)
rtedf$sample <- factor(rtedf$sample, levels = conf$samples)
rtedf$condition <- factor(rtedf$condition, levels = conf$levels)

rtedf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf_promoters.tsv", params$mod_code), col_names = TRUE)
rtedf_promoters$sample <- factor(rtedf_promoters$sample, levels = conf$samples)
rtedf_promoters$condition <- factor(rtedf_promoters$condition, levels = conf$levels)

l1hs_intrautr <- read_delim(sprintf("ldna/Rintermediates/%s/l1hs_intrautr.tsv", params$mod_code), col_names = TRUE)
l1hs_intrautr$sample <- factor(l1hs_intrautr$sample, levels = conf$samples)
l1hs_intrautr$condition <- factor(l1hs_intrautr$condition, levels = conf$levels)


perelementdf <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf.tsv", params$mod_code), col_names = TRUE)
perelementdf$sample <- factor(perelementdf$sample, levels = conf$samples)
perelementdf$condition <- factor(perelementdf$condition, levels = conf$levels)

perelementdf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf_promoters.tsv", params$mod_code), col_names = TRUE)
perelementdf_promoters$sample <- factor(perelementdf_promoters$sample, levels = conf$samples)
perelementdf_promoters$condition <- factor(perelementdf_promoters$condition, levels = conf$levels)

perl1hs_5utr_region <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE) %>% mutate(region = ordered(region, levels = c("328", "500", "909", "ASP")))
perl1hs_5utr_region$sample <- factor(perl1hs_5utr_region$sample, levels = conf$samples)
perl1hs_5utr_region$condition <- factor(perl1hs_5utr_region$condition, levels = conf$levels)

cpg_islands <- rtracklayer::import(conf$cpg_islands)
cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)

readscg <- read_delim(sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE) %>%
    mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels))

refseq_gr <- import(conf$refseq_unaltered)
genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
mcols(genes_gr) %>% colnames()
mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
write_delim(tibble(as.data.frame(promoters)) %>% mutate(score = 1000) %>% dplyr::select(seqnames, start, end, gene_id, score, strand), sprintf("ldna/Rintermediates/%s/promoters.bed", params$mod_code), col_names = FALSE, delim = "\t")


if ((conf$single_condition == "no")) {
    dmrs <- read_delim(params$dmrs, delim = "\t", col_names = TRUE) %>% filter(dmr_type %in% c("t01", "t05"))
    dmls <- read_delim(params$dmls, delim = "\t", col_names = TRUE) %>% filter(fdrs <= 0.2)

    dmrsgr <- GRanges(dmrs)
    dmlsgr <- GRanges(
        seqnames = dmls$chr,
        ranges = IRanges(start = dmls$pos, end = dmls$pos),
        stat = dmls$stat,
        pval = dmls$pvals,
        fdr = dmls$fdrs,
        direction = dmls$direction
    )
    dmrsannot <- dmrs %>%
        mutate(direction_threshold = paste(direction, gsub("t", "", dmr_type), sep = "_")) %>%
        GRanges()
    dmrsgr_split <- split(dmrsannot, dmrsannot$direction_threshold)
}

outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)
cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()
cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)


##################################### DML / DMR analysis
if ((conf$single_condition == "no")) {
    dmrtypes <- dmrs$dmr_type %>% unique()
    dmr_grs_cpg_islands <- dmrsgr %>% subsetByOverlaps(cpg_islands)
    dmr_grs_cpg_islands$islandStatus <- "island"
    dmr_grs_cpgi_shores <- dmrsgr %>% subsetByOverlaps(cpgi_shores)
    dmr_grs_cpgi_shores_filtered <- dmr_grs_cpgi_shores %>% subsetByOverlaps(dmr_grs_cpg_islands, invert = TRUE)
    dmr_grs_cpgi_shores_filtered$islandStatus <- "shore"
    dmr_grs_cpgi_shelves <- dmrsgr %>% subsetByOverlaps(cpgi_shelves)
    dmr_grs_cpgi_shelves_filtered <- dmr_grs_cpgi_shelves %>% subsetByOverlaps(dmr_grs_cpgi_shores, invert = TRUE)
    dmr_grs_cpgi_shelves_filtered$islandStatus <- "shelf"
    dmr_grs_cpg_opensea <- dmrsgr %>% subsetByOverlaps(cpgi_features, invert = TRUE)
    dmr_grs_cpg_opensea$islandStatus <- "opensea"
    dmrsgrislandStatusdf <- c(dmr_grs_cpg_islands, dmr_grs_cpgi_shores_filtered, dmr_grs_cpgi_shelves_filtered, dmr_grs_cpg_opensea) %>%
        as.data.frame() %>%
        tibble()

    for (dmrtype in dmrtypes) {
        dmrs_temp <- dmrsgrislandStatusdf %>% filter(dmr_type == dmrtype)
        dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
        p <- dmrs_temp %>%
            ggplot() +
            geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
            labs(x = "") +
            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
            ggtitle("DMR Counts") +
            mtopen +
            scale_methylation
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_number.pdf", params$mod_code, dmrtype), 4, 4)

        # what is their average length
        p <- ggplot(data = dmrs_temp) +
            geom_histogram(aes(length), fill = mycolor, color = "black") +
            ggtitle("DMR Lengths") +
            labs(x = "length (bp)") +
            xlim(0, 3000) +
            mtopen +
            anchorbar
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_length.pdf", params$mod_code, dmrtype), w = 4, h = 4)

        p <- ggplot(data = dmrs_temp) +
            geom_histogram(aes(length, fill = direction), alpha = 0.7, color = "black") +
            ggtitle("DMR Lengths") +
            labs(x = "length (bp)") +
            xlim(0, 3000) +
            mtopen +
            scale_methylation +
            anchorbar
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_length_stratified.pdf", params$mod_code, dmrtype), w = 4, h = 4)

        ## where are they?

        p <- dmrsgrislandStatusdf %>%
            group_by(islandStatus, direction) %>%
            summarize(n = n()) %>%
            ggplot() +
            geom_col(aes(x = islandStatus, y = n, fill = direction), position = "dodge", color = "black") +
            mtopen +
            scale_methylation +
            anchorbar +
            labs(x = "", y = "Count") +
            ggtitle("DMR") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_count_islandstatus.pdf", params$mod_code, dmrtype), 5, 4)


        dmrlocdf <- dmrs_temp %>%
            group_by(seqnames, direction) %>%
            summarize(n = n())
        dmrlocdf$seqnames <- factor(dmrlocdf$seqnames, levels = chromosomes)
        p <- ggplot(data = dmrlocdf) +
            geom_col(aes(y = seqnames, x = n, fill = direction), position = "dodge", color = "black") +
            theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
            labs(x = "count", y = "") +
            ggtitle("DMR Location") +
            mtopen +
            scale_methylation
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_count.pdf", params$mod_code, dmrtype), 5, 5)

        # p <- dmrs %>%
        #     ggplot(aes(x = meanMethy_c1, y = meanMethy_c2)) +
        #     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
        #     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
        #     scale_x_continuous(expand = c(0, 0)) +
        #     scale_y_continuous(expand = c(0, 0)) +
        #     scale_fill_distiller(palette = "Spectral", direction = 1) +
        #     xlab(sprintf("CpG Methylation %s", condition1)) +
        #     ylab(sprintf("CpG Methylation %s", condition2)) +
        #     ggtitle("DMR Density") +
        #     mtclosed
        # mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrdensity.pdf", params$mod_code), 5, 5)
    }

    p <- dmrsgrislandStatusdf %>%
        filter(dmr_type != "t05CG10") %>%
        filter(dmr_type != "t001") %>%
        group_by(direction, dmr_type) %>%
        summarize(n = n()) %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        ggplot() +
        geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
        geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
        mtopen +
        scale_methylation_thresholds +
        annotation_logticks(sides = "l") +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", math_format(10^.x))
        ) +
        labs(x = "", y = "Count") +
        ggtitle("DMR") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_count.pdf", params$mod_code, "both"), 5, 4)

    p <- dmrsgrislandStatusdf %>%
        filter(dmr_type != "t05CG10", dmr_type != "t001") %>%
        group_by(direction, dmr_type) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        ggplot() +
        geom_col(
            data = . %>% filter(dmr_type == "t05"),
            aes(x = direction, y = n, group = direction, fill = direction_threshold),
            position = position_dodge(), color = "black"
        ) +
        geom_col(
            data = . %>% filter(dmr_type == "t01"),
            aes(x = direction, y = n, group = direction, fill = direction_threshold),
            position = position_dodge(), color = "black"
        ) +
        # Add text labels on top of bars
        geom_text(
            data = . %>% filter(dmr_type %in% c("t01", "t05")),
            aes(x = direction, y = n, label = n, group = direction),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3
        ) +
        mtopen +
        scale_methylation_thresholds +
        annotation_logticks(sides = "l") +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", math_format(10^.x))
        ) +
        labs(x = "", y = "Count") +
        ggtitle("DMR")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_count_numannot.pdf", params$mod_code, "both"), 3.75, 4)

    p <- dmrsgrislandStatusdf %>%
        filter(dmr_type != "t05CG10") %>%
        filter(dmr_type != "t001") %>%
        group_by(islandStatus, direction, dmr_type) %>%
        summarize(n = n()) %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        ggplot() +
        geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = islandStatus, y = n, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
        geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = islandStatus, y = n, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
        mtopen +
        scale_methylation_thresholds +
        annotation_logticks(sides = "l") +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", math_format(10^.x))
        ) +
        labs(x = "", y = "Count") +
        ggtitle("DMR") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_count_islandstatus.pdf", params$mod_code, "both"), 5, 4)

    p <- dmls %>%
        ggplot() +
        geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
        ggtitle("DML Counts") +
        labs(x = "", y = "Count") +
        anchorbar +
        mtclosed +
        scale_methylation
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dml_count.pdf", params$mod_code), 4, 4)
    # dmls
    # p <- ggplot() +
    #     geom_density(data = dmls, aes(x = diff_c2_minus_c1), fill = mycolor) +
    #     geom_vline(xintercept = 0, linetype = "dashed") +
    #     labs(x = sprintf("DML Methylation (%s-%s)", condition2, condition1), y = "Density") +
    #     ggtitle("DML Methylation Density") +
    #     annotate("label", x = -Inf, y = Inf, label = "Hypo", hjust = 0, vjust = 1) +
    #     annotate("label", x = Inf, y = Inf, label = "Hyper", hjust = 1, vjust = 1) +
    #     mtopen +
    #     scale_contrasts +
    #     anchorbar

    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dml_delta.pdf", params$mod_code), 4, 4)

    dmllocdf <- dmls %>%
        group_by(chr, direction) %>%
        summarize(n = n())
    dmllocdf$chr <- factor(dmllocdf$chr, levels = chromosomes)

    # p <- dmls %>%
    #     filter(chr %in% chromosomes) %>%
    #     ggplot(aes(x = mu_c1, y = mu_c2)) +
    #     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    #     scale_x_continuous(expand = c(0, 0)) +
    #     scale_y_continuous(expand = c(0, 0)) +
    #     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    #     scale_fill_distiller(palette = "Spectral", direction = 1) +
    #     xlab(sprintf("CpG Methylation %s", condition1)) +
    #     ylab(sprintf("CpG Methylation %s", condition2)) +
    #     ggtitle("DML Density") +
    #     mtclosed
    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmldensity.pdf", params$mod_code), w = 5, h = 5)



    dml_grs_cpg_islands <- dmlsgr %>% subsetByOverlaps(cpg_islands)
    dml_grs_cpg_islands$islandStatus <- "island"
    dml_grs_cpgi_shores <- dmlsgr %>% subsetByOverlaps(cpgi_shores)
    dml_grs_cpgi_shores_filtered <- dml_grs_cpgi_shores %>% subsetByOverlaps(dml_grs_cpg_islands, invert = TRUE)
    dml_grs_cpgi_shores_filtered$islandStatus <- "shore"
    dml_grs_cpgi_shelves <- dmlsgr %>% subsetByOverlaps(cpgi_shelves)
    dml_grs_cpgi_shelves_filtered <- dml_grs_cpgi_shelves %>% subsetByOverlaps(dml_grs_cpgi_shores, invert = TRUE)
    dml_grs_cpgi_shelves_filtered$islandStatus <- "shelf"
    dml_grs_cpg_opensea <- dmlsgr %>% subsetByOverlaps(cpgi_features, invert = TRUE)
    dml_grs_cpg_opensea$islandStatus <- "opensea"
    dmlsgrislandStatusdf <- c(dml_grs_cpg_islands, dml_grs_cpgi_shores_filtered, dml_grs_cpgi_shelves_filtered, dml_grs_cpg_opensea) %>%
        as.data.frame() %>%
        tibble()

    p <- dmlsgrislandStatusdf %>%
        group_by(islandStatus, direction) %>%
        summarize(n = n()) %>%
        ggplot() +
        geom_col(aes(x = islandStatus, y = n, fill = direction), position = "dodge", color = "black") +
        mtopen +
        scale_methylation +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        labs(x = "", y = "Count") +
        ggtitle("DML Island Status")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dml_count_islandstatus.pdf", params$mod_code), 4, 4)
}



if ((conf$single_condition == "no")) {
    dir.create(sprintf("ldna/results/%s/plots/figs", params$mod_code))
    dmrtypes <- dmrs$dmr_type %>% unique()


    flRTEpromoterlong <- pivot_longer(data = flRTEpromoter %>% dplyr::select(-t001, -t05CG10), cols = dmrtypes, names_to = "dmr_type", values_to = "direction")

    pff <- flRTEpromoterlong %>%
        group_by(rte_subfamily, dmr_type) %>%
        mutate(group_n = n()) %>%
        ungroup() %>%
        filter(!is.na(direction)) %>%
        group_by(rte_subfamily, dmr_type, direction) %>%
        summarise(n = n(), group_n = dplyr::first(group_n)) %>%
        mutate(frac_dm = n / group_n) %>%
        filter(!is.na(rte_subfamily)) %>%
        filter(direction != "discordant") %>%
        ungroup() %>%
        complete(rte_subfamily, dmr_type, direction, fill = list(n = 0, group_n = 0, frac_dm = 0)) %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type)))

    numdf <- flRTEpromoterlong %>%
        group_by(rte_subfamily, dmr_type) %>%
        summarise(group_n_accurate = n())
    flRTEpromoterlong %>%
        filter(rte_subfamily == "L1HS") %>%
        pw()
    p <- pff %>%
        left_join(numdf) %>%
        mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
        ggplot() +
        geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        labs(x = "", y = "Fraction Differentially Methylated") +
        ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
        mtclosed +
        scale_methylation_thresholds
    mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter.pdf", params$mod_code, "all", "all"), 12, 4)

    p <- pff %>%
        left_join(numdf) %>%
        mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
        filter(rte_subfamily != "Other") %>%
        ggplot() +
        geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        labs(x = "", y = "Fraction Differentially Methylated") +
        ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
        mtclosed +
        scale_methylation_thresholds
    mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter_noOther.pdf", params$mod_code, "all", "all"), 12, 4)

    p <- pff %>%
        left_join(numdf) %>%
        mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
        filter(grepl("^L1", rte_subfamily)) %>%
        ggplot() +
        geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
        labs(x = "", y = "Fraction Differentially Methylated") +
        ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
        mtclosed +
        scale_methylation_thresholds
    mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter_l1s.pdf", params$mod_code, "all", "all"), 10, 4)


    for (dmrtype in dmrtypes) {
        pff <- flRTEpromoter %>%
            group_by(rte_subfamily, genic_loc) %>%
            mutate(group_n = n()) %>%
            group_by(rte_subfamily, !!sym(dmrtype), genic_loc) %>%
            summarise(n = n(), group_n = dplyr::first(group_n)) %>%
            mutate(frac_dm = n / group_n) %>%
            filter(!is.na(rte_subfamily)) %>%
            filter(!is.na(!!sym(dmrtype))) %>%
            filter(!!sym(dmrtype) != "discordant") %>%
            ungroup() %>%
            complete(rte_subfamily, !!sym(dmrtype), genic_loc, fill = list(n = 0, group_n = 0, frac_dm = 0))
        numdf <- flRTEpromoter %>%
            group_by(rte_subfamily, genic_loc) %>%
            summarise(group_n_accurate = n())
        p <- pff %>%
            left_join(numdf) %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            ggplot() +
            geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
            facet_wrap(~genic_loc, scales = "free_x", nrow = 2) +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
            mtclosed +
            scale_methylation
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter_regionstrat.pdf", params$mod_code, dmrtype, "all"), 12, 7)

        pff <- flRTEpromoter %>%
            group_by(rte_subfamily) %>%
            mutate(group_n = n()) %>%
            group_by(rte_subfamily, !!sym(dmrtype)) %>%
            summarise(n = n(), group_n = dplyr::first(group_n)) %>%
            mutate(frac_dm = n / group_n) %>%
            filter(!is.na(rte_subfamily)) %>%
            filter(!is.na(!!sym(dmrtype))) %>%
            filter(!!sym(dmrtype) != "discordant") %>%
            ungroup() %>%
            complete(rte_subfamily, !!sym(dmrtype), fill = list(n = 0, group_n = 0, frac_dm = 0))
        numdf <- flRTEpromoter %>%
            group_by(rte_subfamily) %>%
            summarise(group_n_accurate = n())
        p <- pff %>%
            left_join(numdf) %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            ggplot() +
            geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
            mtclosed +
            scale_methylation
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter.pdf", params$mod_code, dmrtype, "all"), 12, 4)

        p <- pff %>%
            left_join(numdf) %>%
            filter(rte_subfamily != "Other") %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            ggplot() +
            geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
            mtclosed +
            scale_methylation
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter1.pdf", params$mod_code, dmrtype, "all"), 10, 4)

        p <- pff %>%
            left_join(numdf) %>%
            filter(grepl("^L1", rte_subfamily)) %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            ggplot() +
            geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
            mtclosed +
            scale_methylation
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/dm_%s_fl%s_promoter2.pdf", params$mod_code, dmrtype, "all"), 8, 4)
    }
}


### CUSTOM
p <- perelementdf_promoters %>%
    filter(sample == conf$samples[[1]]) %>%
    filter(grepl("^L1", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = rte_subfamily, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1 CpG Methylation") +
    mtopen +
    scale_conditions
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s_sample1only.pdf", params$mod_code), raster = TRUE, 12, 4)
#####

p <- perelementdf_promoters %>%
    group_by(rte_subfamily, sample) %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_subfamily, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_subfamily, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        group_by(sample, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters.pdf", params$mod_code), 12, 5, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters.pdf", params$mod_code), raster = TRUE, 12, 5)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters.pdf", params$mod_code), raster = TRUE, 12, 5)
}

p <- perelementdf_promoters %>%
    filter(!grepl("HERV", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_subfamily, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_subfamily, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        group_by(sample, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_1.pdf", params$mod_code), 14, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_1.pdf", params$mod_code), raster = TRUE, 14, 6)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_1.pdf", params$mod_code), raster = TRUE, 14, 6)
}
### CUSTOM
p <- perelementdf_promoters %>%
    filter(condition == condition1) %>%
    filter(grepl("^L1", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = rte_subfamily, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1 CpG Methylation") +
    mtopen +
    scale_conditions
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s_condition1only.pdf", params$mod_code), raster = TRUE, 12, 4)
#####

p <- perelementdf_promoters %>%
    filter(grepl("^L1", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = rte_subfamily, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1 CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1", rte_subfamily)) %>%
        group_by(sample, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s.pdf", params$mod_code), 14, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s.pdf", params$mod_code), raster = TRUE, 12, 4)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = sample, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    facet_wrap(~loc_lowres_integrative_stranded) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(sample, condition, rte_subfamily, loc_lowres_integrative_stranded) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS_loc.pdf", params$mod_code), 12, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS_loc.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS_loc.pdf", params$mod_code), raster = TRUE, 12, 4)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = loc_lowres_integrative_stranded, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(sample, condition, rte_subfamily, loc_lowres_integrative_stranded) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ loc_lowres_integrative_stranded, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc.pdf", params$mod_code), 6, 4, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoterse_by_condition_L1HS_loc.pdf", params$mod_code), raster = TRUE, 6, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc.pdf", params$mod_code), raster = TRUE, 6, 6)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    mutate(loc_lowres_integrative_stranded = case_when(
        loc_lowres_integrative_stranded == "Gene_Adj_Antisense" | loc_lowres_integrative_stranded == "Gene_Adj_Sense" ~ "Intergenic",
        TRUE ~ loc_lowres_integrative_stranded
    )) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = loc_lowres_integrative_stranded, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        mutate(loc_lowres_integrative_stranded = case_when(
            loc_lowres_integrative_stranded == "Gene_Adj_Antisense" | loc_lowres_integrative_stranded == "Gene_Adj_Sense" ~ "Intergenic",
            TRUE ~ loc_lowres_integrative_stranded
        )) %>%
        group_by(sample, condition, rte_subfamily, loc_lowres_integrative_stranded) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ loc_lowres_integrative_stranded, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc_lowres.pdf", params$mod_code), 5, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoterse_by_condition_L1HS_loc_lowres.pdf", params$mod_code), raster = TRUE, 5, 6)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc_lowres.pdf", params$mod_code), raster = TRUE, 5, 6)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    ggplot(aes(x = sample, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(sample, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS.pdf", params$mod_code), 10, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS.pdf", params$mod_code), raster = TRUE, 12, 4)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    filter(intactness_req == "Intact") %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("Intact L1HS CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        filter(intactness_req == "Intact") %>%
        group_by(sample, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        ungroup() %>%
        compare_means(mean_meth ~ condition, data = ., method = "t.test", group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_intactL1s.pdf", params$mod_code), 10, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
}


p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    left_join(sample_table) %>%
    filter(intactness_req == "Intact") %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot(aes(x = braak, y = mean_meth, color = sample)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    stats <- perelementdf_promoters %>%
        compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), 10, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    left_join(sample_table) %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = mean_meth, color = braak)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    scale_palette_alt +
    mtopen
if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(sample) %>%
        summarize(median = median(mean_meth)) %>%
        left_join(sample_table) %>%
        ungroup() %>%
        group_by(condition) %>%
        summarize(mean_of_median = mean(median))
    stats <- perelementdf_promoters %>%
        compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_L1s.pdf", params$mod_code), 10, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_L1s.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_L1s.pdf", params$mod_code), raster = TRUE, 12, 4)
}

if ((conf$single_condition == "no")) {
    pfl1 <- perelementdf_promoters %>%
        filter(grepl("^L1", rte_subfamily)) %>%
        dplyr::select(-t05CG10, -t001)
    p <- pfl1 %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        group_by(rte_subfamily) %>%
        mutate(rank_change = row_number()) %>%
        mutate(top_change = ifelse(rank_change <= 10, "Top", "NotTop")) %>%
        arrange(abs_dif) %>%
        ungroup() %>%
        ggpaired(cond1 = condition1, cond2 = condition2, line.color = "top_change", alpha = "top_change", facet.by = "rte_subfamily") +
        scale_alpha_manual(values = c(1, 0.5)) +
        scale_color_manual(values = c("Top" = "red", "NotTop" = "grey")) +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_L1s.pdf", params$mod_code), 14, 6, raster = TRUE)

    pfl1hs <- pfl1 %>%
        filter(rte_subfamily == "L1HS")
    pfl1hs %>% arrange(mean_meth)
    l1hs_paired_dif_frame <- pfl1 %>%
        filter(rte_subfamily == "L1HS") %>%
        pivot_longer(cols = dmrtypes, names_to = "dmr_type", values_to = "direction") %>%
        mutate(direction_threshold = ifelse(is.na(direction), "NS", paste0(direction, "_", gsub("t", "", dmr_type)))) %>%
        filter(!(dmr_type == "t01" & is.na(direction))) %>%
        # filter(!(dmr_type == "t05" & is.na(direction))) %>%
        group_by(gene_id, rte_subfamily, condition, dmr_type, direction_threshold) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        mutate(
            x_1 = factor(condition1, levels = conf$levels),
            x_2 = factor(condition2, levels = conf$levels),
            y_1 = !!sym(condition1),
            y_2 = !!sym(condition2)
        )

    p <- ggplot() +
        geom_boxplot(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold != "Hypo_01") %>% filter(direction_threshold != "Hyper_01"),
            aes(x = x_1, y = y_1)
        ) +
        geom_boxplot(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold != "Hypo_01") %>% filter(direction_threshold != "Hyper_01"),
            aes(x = x_2, y = y_2)
        ) +
        geom_segment(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold == "NS"),
            aes(x = x_1, y = y_1, xend = x_2, yend = y_2, color = direction_threshold),
            alpha = 0.35
        ) +
        geom_segment(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold == "Hyper_05"),
            aes(x = x_1, y = y_1, xend = x_2, yend = y_2, color = direction_threshold)
        ) +
        geom_segment(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold == "Hypo_05"),
            aes(x = x_1, y = y_1, xend = x_2, yend = y_2, color = direction_threshold)
        ) +
        geom_segment(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold == "Hyper_01"),
            aes(x = x_1, y = y_1, xend = x_2, yend = y_2, color = direction_threshold)
        ) +
        geom_segment(
            data = l1hs_paired_dif_frame %>% filter(direction_threshold == "Hypo_01"),
            aes(x = x_1, y = y_1, xend = x_2, yend = y_2, color = direction_threshold)
        ) +
        labs(y = "L1HS 5UTR Methylation", x = "Condition") +
        scale_color_manual(values = c("Hypo_01" = "#ffa200", "Hyper_01" = "49fdfa", "Hypo_05" = "red", "Hyper_05" = "blue", "NS" = "grey")) +
        mtclosedgridh
    mysaveandstore(pl = p, fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_l1hs_%s.pdf", params$mod_code, "all"), 5, 4)

    for (dmrtype in dmrs$dmr_type %>% unique()) {
        p <- pfl1 %>%
            mutate(!!sym(dmrtype) := ifelse(is.na(!!sym(dmrtype)), "NS", !!sym(dmrtype))) %>%
            group_by(gene_id, rte_subfamily, condition, !!sym(dmrtype)) %>%
            summarize(mean_meth = mean(mean_meth)) %>%
            pivot_wider(names_from = condition, values_from = mean_meth) %>%
            ungroup() %>%
            mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
            mutate(abs_dif = abs(dif)) %>%
            arrange(-abs_dif) %>%
            group_by(rte_subfamily) %>%
            arrange(abs_dif) %>%
            ungroup() %>%
            ggpaired(cond1 = condition1, cond2 = condition2, line.color = dmrtype, alpha = dmrtype, facet.by = "rte_subfamily") +
            scale_color_manual(values = c("Hypo" = "red", "Hyper" = "blue", "NS" = "grey")) +
            scale_alpha_manual(values = c(1, 0.1)) +
            mtclosedgridh
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_L1s2_%s.pdf", params$mod_code, dmrtype), 14, 6, raster = TRUE)

        p <- pfl1 %>%
            filter(rte_subfamily == "L1HS") %>%
            mutate(!!sym(dmrtype) := ifelse(is.na(!!sym(dmrtype)), "NS", !!sym(dmrtype))) %>%
            group_by(gene_id, rte_subfamily, condition, !!sym(dmrtype)) %>%
            summarize(mean_meth = mean(mean_meth)) %>%
            pivot_wider(names_from = condition, values_from = mean_meth) %>%
            ungroup() %>%
            mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
            mutate(abs_dif = abs(dif)) %>%
            arrange(-abs_dif) %>%
            group_by(rte_subfamily) %>%
            arrange(abs_dif) %>%
            ungroup() %>%
            ggpaired(cond1 = condition1, cond2 = condition2, line.color = dmrtype, alpha = 0.85, ylab = "L1HS 5UTR Methylation") +
            scale_color_manual(values = c("Hypo" = "red", "Hyper" = "blue", "NS" = "grey")) +
            mtclosedgridh
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_l1hs_%s.pdf", params$mod_code, dmrtype), 4, 4, raster = FALSE)
    }

    top_l1hs_movers <- pfl1 %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        group_by(rte_subfamily) %>%
        mutate(rank_change = row_number()) %>%
        mutate(top_change = ifelse(rank_change <= 10, "Top", "NotTop")) %>%
        ungroup() %>%
        filter(rte_subfamily == "L1HS") %$% gene_id %>%
        head(n = 10)

    top_l1hs_movers_intact <- pfl1 %>%
        filter(intactness_req == "Intact") %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        group_by(rte_subfamily) %>%
        mutate(rank_change = row_number()) %>%
        mutate(top_change = ifelse(rank_change <= 10, "Top", "NotTop")) %>%
        ungroup() %>%
        filter(rte_subfamily == "L1HS") %$% gene_id %>%
        head(n = 10)


    pf <- perelementdf_promoters %>%
        filter(rte_subfamily == "L1HS")
    p <- pf %>%
        ggplot() +
        geom_quasirandom(aes(x = intactness_req, y = mean_meth, color = condition), dodge.width = 0.75) +
        geom_boxplot(aes(x = intactness_req, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
        xlab("") +
        ylab("Average CpG Methylation Per Element") +
        ggtitle("RTE CpG Methylation") +
        geom_pwc(
            data = pf %>% group_by(sample, condition, intactness_req) %>% summarize(mean_meth = mean(mean_meth)), aes(x = intactness_req, y = mean_meth, group = condition), tip.length = 0,
            method = "t.test", label = "{p.adj.format}",
            p.adjust.method = "fdr", p.adjust.by = "panel",
            hide.ns = FALSE
        ) +
        mtopen +
        scale_conditions
    p <- pf %>%
        ggplot(aes(x = intactness_req, y = mean_meth, color = condition)) +
        geom_quasirandom(dodge.width = 0.75) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        xlab("") +
        ylab("Average CpG Methylation Per Element") +
        ggtitle("RTE CpG Methylation") +
        geom_pwc(aes(group = condition),
            tip.length = 0,
            method = "t.test", label = "{p.adj.format}",
            p.adjust.method = "fdr", p.adjust.by = "panel",
            hide.ns = FALSE
        ) +
        mtopen +
        scale_conditions
    tryCatch(
        {
            stats <- pf %>%
                compare_means(mean_meth ~ condition, group.by = "intactness_req", data = ., method = "t.test", p.adjust.method = "fdr")
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_promoters.pdf", params$mod_code), 5, 4, sf = stats)
        },
        error = function(e) {
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_promoters.pdf", params$mod_code), 5, 4)
        }
    )
}


pf <- perl1hs_5utr_region

p <- pf %>%
    group_by(gene_id, region) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ggplot(aes(x = mean_meth)) +
    geom_histogram(color = "black") +
    facet_wrap(~region, nrow = 2) +
    mtclosedgrid
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_histogram_5utr_regions.pdf", params$mod_code), 5, 5)


# Define function to compute ECDF percentages
compute_ecdf_df <- function(df, region, group_vars, breakpoints, by_sample) {
    df |>
        filter(region == region) |>
        group_by(across(all_of(group_vars))) |>
        summarise(mean_meth = mean(mean_meth), .groups = "drop") |>
        pull(mean_meth) |>
        (\(x) {
            pct_below <- ecdf(x)(breakpoints) * 100
            tibble(threshold = breakpoints, percent_below = pct_below)
        })() |>
        mutate(region = region, by_sample = by_sample)
}

# Define regions and breakpoints
regions <- c("909", "500", "328")
breakpoints <- seq(0, 100, 5)

# Compute ECDF for both grouping methods
ecdf_gene <- map_dfr(regions, ~ compute_ecdf_df(pf, .x, "gene_id", breakpoints, FALSE))
ecdf_sample <- map_dfr(regions, ~ compute_ecdf_df(pf, .x, c("gene_id", "sample"), breakpoints, TRUE))

# Combine both results
ecdf_df <- bind_rows(ecdf_gene, ecdf_sample)
ecdf_df %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_meth_ecdf.csv", params$mod_code))


# quantiles
# Define function to compute quantiles
compute_quantile_df <- function(df, region, group_vars, probs, by_sample) {
    df |>
        filter(region == region) |>
        group_by(across(all_of(group_vars))) |>
        summarise(mean_meth = mean(mean_meth), .groups = "drop") |>
        pull(mean_meth) |>
        (\(x) {
            quant_values <- quantile(x, probs)
            tibble(quantile = probs, value = quant_values)
        })() |>
        mutate(region = region, by_sample = by_sample)
}

# Define regions and quantile probabilities
regions <- c("909", "500", "328")
quantile_probs <- seq(0, 1, 0.05) # 5% increments

# Compute quantiles for both grouping methods
quantile_gene <- map_dfr(regions, ~ compute_quantile_df(pf, .x, "gene_id", quantile_probs, FALSE))
quantile_sample <- map_dfr(regions, ~ compute_quantile_df(pf, .x, c("gene_id", "sample"), quantile_probs, TRUE))

# Combine both results
quantile_df <- bind_rows(quantile_gene, quantile_sample)
quantile_df %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_meth_quantiles.csv", params$mod_code))

p <- pf %>%
    filter(region == 909) %>%
    group_by(sample, condition) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggboxplot(x = "condition", y = "mean_meth", fill = "condition", add = c("mean_se", "dotplot")) +
    scale_conditions +
    scale_conditions +
    mtopen
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/mean_meth_point_fll1hs_box.pdf", params$mod_code), 4.5, 3.75, sf = stats)


p <- pf %>%
    ggplot(aes(x = region, y = mean_meth, color = condition)) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
tryCatch(
    {
        mdf <- l1hs_intrautr %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

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

        merged <- left_join(mdf, cg_positions_df, by = c("gene_id", "sequence_pos"))
        cpg_order <- merged %$% consensus_pos %>%
            unique() %>%
            sort()
        merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))

        dat <- merged %>%
            filter(!is.na(pctM))
        dat1 <- merged %>%
            filter(!is.na(consensus_pos))
        merged %>%
            filter(is.na(pctM)) %>%
            pw()
        # write_csv(data_model, "meth_for_sarah.csv")

        data_model <- dat %>%
            mutate(consensus_pos = as.character(consensus_pos)) %>%
            filter(region == "909") %>%
            dplyr::rename(sample_name = sample) %>%
            left_join(sample_table) %>%
            mutate(
                total_sites = cov,
                methylated_sites = round(cov * pctM / 100)
            ) %>%
            mutate(age_z = as.numeric(scale(age))) %>%
            mutate(condition = factor(condition, levels = conf$levels)) %>%
            dplyr::select(sample_name, condition, sex, age_z, seqnames, start, gene_id, consensus_pos, total_sites, methylated_sites)

        library(glmmTMB)
        library(broom.mixed)

        global_model <- glmmTMB(
            cbind(methylated_sites, total_sites - methylated_sites) ~
                condition + sex + age_z + (1 | sample_name) + (1 | gene_id) + (1 | consensus_pos),
            data = data_model,
            family = binomial()
        )


        broom::tidy(global_model) %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909.csv", params$mod_code))
        # now data for bayes
        data_model %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909_data.csv", params$mod_code))


        global_model_with_interaction <- glmmTMB(
            cbind(methylated_sites, total_sites - methylated_sites) ~
                condition * sex + age_z + (1 | sample_name) + (1 | gene_id) + (1 | consensus_pos),
            data = data_model,
            family = binomial()
        )
        broom::tidy(global_model_with_interaction) %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_global_hierarchical_model_909_withsexinteraction.csv", params$mod_code))

        gene_model_results <- data_model %>%
            group_by(gene_id) %>%
            group_split() %>%
            map_df(~ {
                df <- .x
                if (n_distinct(df$condition) < 2) {
                    return(NULL)
                } # skip genes with only one condition

                model <- tryCatch(
                    {
                        glmmTMB(
                            cbind(methylated_sites, total_sites - methylated_sites) ~
                                condition + sex + age_z + (1 | sample_name) + (1 | consensus_pos),
                            data = df,
                            family = binomial()
                        )
                    },
                    error = function(e) {
                        return(NULL)
                    }
                )

                if (!is.null(model)) {
                    est <- summary(model)$coefficients$cond
                    data.frame(
                        gene_id = unique(df$gene_id),
                        condition_effect = est["conditionAD", "Estimate"],
                        std_error = est["conditionAD", "Std. Error"],
                        z_value = est["conditionAD", "z value"],
                        p_value = est["conditionAD", "Pr(>|z|)"]
                    )
                } else {
                    NULL
                }
            }) %>%
            mutate(p_adj = p.adjust(p_value, method = "fdr")) # adjust for multiple testing

        gene_model_results %>% write_mycsv(sprintf("ldna/results/%s/tables/rte/fl_l1hs_gene_level_hierarchical_model_909.csv", params$mod_code))

        gene_model_results %>%
            mutate(sig = ifelse(p_adj <= 0.05, TRUE, FALSE)) %>%
            mutate(sig_dir = case_when(
                condition_effect > 0 & sig == TRUE ~ "sig_up",
                condition_effect <= 0 & sig == TRUE ~ "sig_down",
                condition_effect > 0 & sig == FALSE ~ "notsig_up",
                condition_effect <= 0 & sig == FALSE ~ "notsig_down",
            )) %$% sig_dir %>%
            table()


        ####


        res909 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 909) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
            broom::tidy() %>%
            mutate(region = "909") %>%
            mutate(model_type = "no_interaction")

        res500 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 500) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
            broom::tidy() %>%
            mutate(region = "500") %>%
            mutate(model_type = "no_interaction")
        res328 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 328) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
            broom::tidy() %>%
            mutate(region = "328") %>%
            mutate(model_type = "no_interaction")
        resASP <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == "ASP") %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .) %>%
            broom::tidy() %>%
            mutate(region = "ASP") %>%
            mutate(model_type = "no_interaction")

        res909interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 909) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
            broom::tidy() %>%
            mutate(region = "909") %>%
            mutate(model_type = "interaction")
        res500interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 500) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
            broom::tidy() %>%
            mutate(region = "500") %>%
            mutate(model_type = "interaction")
        res328interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == 328) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
            broom::tidy() %>%
            mutate(region = "328") %>%
            mutate(model_type = "interaction")
        resASPinteraction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth) / 100) %>%
            filter(region == "ASP") %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            betareg(formula(sprintf("%s ~ %s", "pctM", ifelse(is.null(conf$linear_model_adjustment_interactions), lm_right_hand_side, paste0(lm_right_hand_side, " + ", paste0(conf$linear_model_adjustment_interactions, collapse = " + "))))), data = .) %>%
            broom::tidy() %>%
            mutate(region = "ASP") %>%
            mutate(model_type = "interaction")

        stats <- bind_rows(res909, res500, res328, resASP, res909interaction, res500interaction, res328interaction, resASPinteraction)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4, sf = stats)
    },
    error = function(e) {
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4)
    }
)


p <- pf %>%
    left_join(rmann) %>%
    filter(region == "909") %>%
    ggplot(aes(x = loc_lowres_integrative_stranded, y = mean_meth, color = condition)) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
tryCatch(
    {
        res <- pf %>%
            left_join(rmann) %>% # Join with rmann
            filter(region == "909") %>% # Filter for region 909
            group_by(sample, loc_lowres_integrative_stranded) %>%
            summarise(pctM = mean(mean_meth), .groups = "drop") %>% # Summarize mean methylation
            left_join(sample_table %>% mutate(sample = sample_name)) %>% # Join sample table
            nest(data = -loc_lowres_integrative_stranded) %>% # Nest data for each locus
            mutate(
                model_results = map(
                    data,
                    ~ lm(formula(sprintf("%s ~ %s", "pctM", lm_right_hand_side)), data = .x) %>%
                        summary() %>%
                        broom::tidy()
                )
            ) %>%
            dplyr::select(-data) %>% # Drop nested data
            unnest(model_results)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_gene_loc.pdf", params$mod_code), 5, 4, sf = res)
    },
    error = function(e) {
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_gene_loc.pdf", params$mod_code), 5, 4)
    }
)


########## PCA
library(PCAtools)
pcaframe <- perl1hs_5utr_region %>%
    filter(region == "909") %>%
    dplyr::select(sample, mean_meth, gene_id) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    pivot_wider(names_from = gene_id, values_from = mean_meth) %>%
    arrange(sample) %>%
    column_to_rownames(var = "sample") %>%
    as.matrix() %>%
    t()
pcaObj <- pca(pcaframe, center = TRUE, scale = FALSE, metadata = sample_table %>% column_to_rownames(var = "sample_name"))

p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/screeplot.pdf", params$mod_code), 5, 4)


# p <- plotloadings(pcaObj,
#     components = getComponents(pcaObj, seq_len(3)),
#     rangeRetain = 0.045, labSize = 2
# ) + mtopen
# mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/loadings.pdf", params$mod_code), 5, 4)

if (is.null(adjustment_set_categorical)) {
    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", colby = "condition",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
} else {
    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1]), colby = "condition",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)
}





pf <- pcaframe %>%
    colMeans() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::rename(sample_name = "rowname", mean_meth = ".") %>%
    tibble() %>%
    left_join(sample_table) %>%
    mutate(sample_name = fct_reorder(sample_name, mean_meth))
p <- pf %>%
    ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
    geom_col() +
    scale_conditions +
    new_scale_fill() +
    geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
    scale_palette +
    mtopen
# stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar1.pdf", params$mod_code), 5, 4)

tryCatch(
    {
        p <- pf %>%
            mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
            ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
            geom_col() +
            scale_conditions +
            new_scale_fill() +
            geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
            scale_palette +
            mtopen
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage.pdf", params$mod_code), 5, 4)

        p <- pf %>%
            mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
            ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
            geom_col() +
            scale_conditions +
            new_scale_fill() +
            geom_tile(aes(x = -1, fill = !!sym(asc)), width = 2) + # Add metadata strip
            scale_palette +
            mtopen
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage_ordered.pdf", params$mod_code), 5, 4)


        p <- pf %>%
            mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
            ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +
            # new_scale_fill() +
            # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
            scale_conditions +
            mtopen
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% broom::tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_ordered.pdf", params$mod_code), 5, 4)

        p <- pf %>%
            mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
            ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +
            # new_scale_fill() +
            # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
            scale_conditions +
            mtopen
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_orderedmeth.pdf", params$mod_code), 5, 4)
    },
    error = function(e) {

    }
)
# split_frame <- pf %>% group_by(sample, condition, rte_subfamily) %>%
#     summarise(mean_meth = mean(mean_meth)) %>%
#     ungroup() %>%
#     pivot_wider(names_from = sample, values_from = mean_meth) %>%
#     dplyr::select(-condition) %>% split(.$rte_subfamily)
# results_frame <- data.frame()
# for (frame in split_frame) {
#     condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
#     condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
#     sample_means <- frame %>% dplyr::select(-rte_subfamily) %>% apply(2, mean, na.rm = TRUE)
#     sample_means <- sample_means[c(condition1samples, condition2samples)]
#     condition_vec <- c(condition1samples, condition2samples)
#     condition <- sample_table[match(condition_vec, sample_table$sample_name), ] %$% condition
#     obs_stat <- mean(sample_means[condition == condition2]) - mean(sample_means[condition == condition1])
#     set.seed(123)  # For reproducibility
#     n_permutations <- 10000
#     perm_stats <- replicate(n_permutations, {
#         permuted_condition <- sample(condition)
#         perm_stat <- mean(sample_means[permuted_condition == condition2]) - mean(sample_means[permuted_condition == condition1])
#         perm_stat
#     })
#     p_value <- mean(abs(perm_stats) >= abs(obs_stat))
#     cat("Observed statistic:", obs_stat, "\n")
#     cat("P-value:", p_value, "\n")
#     #store the values
#     results_frame <- rbind(results_frame, data.frame(rte_subfamily = frame$rte_subfamily[1], obs_stat = obs_stat, p_value = p_value))
# }
# sample_means <- sample_means[c(condition1samples, condition2samples)]
# sample_medians <- aa %>%
#     dplyr::select(read_id, sample, fraction_meth) %>%
#     pivot_wider(names_from = sample, values_from = fraction_meth) %>%
#     dplyr::select(-read_id) %>%
#     apply(2, median, na.rm = TRUE)
# sample_medians <- sample_medians[c(condition1samples, condition2samples)]
# condition_vec <- c(condition1samples, condition2samples)
# condition <- sample_table[match(condition_vec, sample_table$sample_name), ] %$% condition
# obs_stat <- mean(sample_medians[condition == condition2]) - mean(sample_medians[condition == condition1])

# # Permutation test
# set.seed(123) # For reproducibility
# n_permutations <- 10000
# perm_stats <- replicate(n_permutations, {
#     permuted_condition <- sample(condition)
#     perm_stat <- mean(sample_medians[permuted_condition == "AD"]) - mean(sample_medians[permuted_condition == "CTRL"])
#     perm_stat
# })

# # Calculate p-value
# p_value <- mean(abs(perm_stats) >= abs(obs_stat))

# # Display result
# cat("Observed statistic:", obs_stat, "\n")
# cat("P-value:", p_value, "\n")

#################

l1hsintactmethgr <- rtedf %>%
    filter(intactness_req == "Intact")
l1hsintactmethgr <- l1hsintactmethgr %>%
    mutate(rel_start = start - rte_start) %>%
    mutate(rel_end = end - rte_start)
write_delim(l1hsintactmethgr, sprintf("ldna/Rintermediates/%s/l1hsintactdf.tsv", params$mod_code), col_names = TRUE)

library(zoo)
pf_pos <- l1hsintactmethgr %>%
    filter(rte_strand == "+") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()
pf_neg <- l1hsintactmethgr %>%
    filter(rte_strand == "-") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()

p <- pf_pos %>% ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtclosed +
    scale_conditions

mysaveandstore(sprintf("ldna/results/%s/plots/rte/l1intact_Lines_pos_strand.pdf", params$mod_code), 12, 30)


p <- pf_pos %>%
    filter(rel_start < 910) %>%
    ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    xlim(c(1, 910)) +
    ylim(c(0, 100)) +
    mtclosed +
    scale_conditions
mysaveandstore(sprintf("ldna/results/%s/plots/rte/l1intact_Lines_pos_strand_promoter.pdf", params$mod_code), 12, 30)

p <- pf_neg %>% ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtclosed +
    scale_conditions
mysaveandstore(sprintf("ldna/results/%s/plots/rte/l1intact_Lines_neg_strand.pdf", params$mod_code), 12, 30)

p <- pf_neg %>%
    filter(rel_start < 910) %>%
    ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    xlim(c(1, 910)) +
    ylim(c(0, 100)) +
    mtclosed +
    scale_conditions
mysaveandstore(sprintf("ldna/results/%s/plots/rte/l1intact_Lines_neg_strand_promoter.pdf", params$mod_code), 12, 30)


if ((conf$single_condition == "no")) {
    element_anatomy <- read_delim("aref/default/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")

    dm_intact_l1hs_elements <- flRTEpromoter %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(intactness_req == "Intact") %>%
        filter(t05 == "Hypo")
    dm_fl_l1hs_elements <- flRTEpromoter %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(t05 == "Hypo")
    topmovers_l1hs_elements <- flRTEpromoter %>%
        filter(gene_id %in% top_l1hs_movers_intact)
    allfl_l1hs_elements <- flRTEpromoter %>%
        filter(rte_subfamily == "L1HS")
    element_sets_of_interst <- list("dm_fl_l1hs" = dm_fl_l1hs_elements, "dm_intact_l1hs" = dm_intact_l1hs_elements, "Top_Movers" = topmovers_l1hs_elements, "all_elements" = allfl_l1hs_elements)



    l1hsflmethgr <- rtedf %>%
        filter(rte_subfamily == "L1HS")
    l1hsflmethgr <- l1hsflmethgr %>%
        mutate(rel_start = start - rte_start) %>%
        mutate(rel_end = end - rte_start)
    write_delim(l1hsflmethgr, sprintf("ldna/Rintermediates/%s/l1hsfldf.tsv", params$mod_code), col_names = TRUE)
    l1hsflmethgr <- read_delim(sprintf("ldna/Rintermediates/%s/l1hsfldf.tsv", params$mod_code), col_names = TRUE)

    library(zoo)
    pf_pos <- l1hsflmethgr %>%
        filter(rte_strand == "+") %>%
        as.data.frame() %>%
        tibble() %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        group_by(gene_id, condition) %>%
        mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
        filter(!is.na(rM)) %>%
        ungroup()
    pf_neg <- l1hsflmethgr %>%
        filter(rte_strand == "-") %>%
        as.data.frame() %>%
        tibble() %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        group_by(gene_id, condition) %>%
        mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
        filter(!is.na(rM)) %>%
        ungroup()

    for (element_type in names(element_sets_of_interst)) {
        df <- element_sets_of_interst[[element_type]]
        dir.create(sprintf("ldna/Rintermediates/%s/l1hs/", params$mod_code), recursive = TRUE)
        write_delim(df %>% dplyr::select(gene_id), sprintf("ldna/Rintermediates/%s/l1hs/%s_gene_id.tsv", params$mod_code, element_type), col_names = FALSE)
        write_delim(df %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/%s/l1hs/%s_promoters.bed", params$mod_code, element_type), col_names = FALSE, delim = "\t")
        write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ] %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/%s/l1hs/%s_full_elements.bed", params$mod_code, element_type), col_names = FALSE, delim = "\t")
        write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ], sprintf("ldna/Rintermediates/%s/l1hs/%s_full_elements.tsv", params$mod_code, element_type), col_names = TRUE, delim = "\t")

        for (element in df$gene_id) {
            y_lim_lower <- 50
            y_lim_upper <- 100
            y_valmin <- y_lim_lower
            y_valmax <- y_lim_lower + ((y_lim_upper - y_lim_lower) / 10)

            if (rmann %>% filter(gene_id == element) %$% strand == "+") {
                modifier <- rmann %>% filter(gene_id == element) %$% start
                color_intervals <- element_anatomy %>%
                    filter(!(feature %in% c("EN", "RT"))) %>%
                    filter(gene_id == element) %>%
                    mutate(across(where(is.numeric), ~ . + modifier))
                pf <- pf_pos
                start_vec <- pf %>%
                    filter(gene_id == element) %$% start
                offset <- min(start_vec)
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    mutate(start = start - offset) %>%
                    ggplot() +
                    geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
                    geom_line(data = . %>%
                        group_by(start, condition) %>%
                        summarise(rM = mean(rM)), aes(x = start, y = rM, color = condition)) +
                    scale_samples_unique +
                    labs(y = "Methylation Rolling Mean") +
                    mtclosed +
                    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
            } else {
                modifier <- rmann %>% filter(gene_id == element) %$% end
                color_intervals <- element_anatomy %>%
                    filter(!(feature %in% c("EN", "RT"))) %>%
                    filter(gene_id == element) %>%
                    mutate(across(where(is.numeric), ~ modifier - .))
                pf <- pf_neg
                start_vec <- pf %>%
                    filter(gene_id == element) %$% start
                offset <- min(start_vec * -1) * -1
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    mutate(start = (-start) + offset) %>%
                    ggplot() +
                    geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
                    geom_line(data = . %>%
                        group_by(start, condition) %>%
                        summarise(rM = mean(rM)), aes(x = start, y = rM, color = condition)) +
                    geom_vline(xintercept = c(0, 909)) +
                    scale_samples_unique +
                    labs(y = "Methylation Rolling Mean") +
                    mtclosed +
                    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
            }

            p <- p1 + plot_layout(heights = c(1))

            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/%s_methylation.pdf", params$mod_code, element_type, element), 5, 5)
            mysaveandstore(pl = p1 + ggtitle(element) + mtclosed, sprintf("ldna/results/%s/plots/rte/%s/%s_methylation_nc.pdf", params$mod_code, element_type, element), 5, 4)
        }

        for (element in df$gene_id) {
            y_lim_lower <- 50
            y_lim_upper <- 100
            y_valmin <- y_lim_lower
            y_valmax <- y_lim_lower + ((y_lim_upper - y_lim_lower) / 10)

            if (rmann %>% filter(gene_id == element) %$% strand == "+") {
                modifier <- rmann %>% filter(gene_id == element) %$% start
                color_intervals <- element_anatomy %>%
                    filter(!(feature %in% c("EN", "RT"))) %>%
                    filter(gene_id == element) %>%
                    mutate(across(where(is.numeric), ~ . + modifier))
                pf <- pf_pos
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    group_by(seqnames, start, end, condition) %>%
                    summarise(rM = mean(rM)) %>%
                    ggplot() +
                    geom_point(aes(x = start, y = rM, fill = condition, color = condition)) +
                    scale_samples_unique +
                    labs(y = "Methylation Rolling Mean") +
                    mtclosed +
                    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
                p2 <- color_intervals %>%
                    ggplot() +
                    geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
                    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
                    geom_text(aes(x = -200 + ((start + end) / 2), y = 1.5, label = feature)) +
                    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
                    ggtitle(element) +
                    scale_fill_paletteer_d("dutchmasters::milkmaid") +
                    theme_map() +
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
                    scale_y_continuous(expand = c(0, 0.4)) +
                    theme(legend.position = "none")
            } else {
                modifier <- rmann %>% filter(gene_id == element) %$% end
                color_intervals <- element_anatomy %>%
                    filter(!(feature %in% c("EN", "RT"))) %>%
                    filter(gene_id == element) %>%
                    mutate(across(where(is.numeric), ~ modifier - .))
                pf <- pf_neg
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    group_by(seqnames, start, end, condition) %>%
                    summarise(rM = mean(rM)) %>%
                    ggplot() +
                    geom_point(aes(x = start, y = rM, fill = condition, color = condition)) +
                    scale_conditions +
                    labs(y = "Methylation Rolling Mean") +
                    scale_x_reverse() +
                    mtclosed +
                    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
                p2 <- color_intervals %>%
                    ggplot() +
                    geom_rect(aes(xmin = -start, xmax = -end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
                    geom_rect(aes(xmin = -start, xmax = -end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
                    geom_text(aes(x = -200 + ((-start + -end) / 2), y = 1.5, label = feature)) +
                    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
                    ggtitle(element) +
                    scale_fill_paletteer_d("dutchmasters::milkmaid") +
                    theme_map() +
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
                    scale_y_continuous(expand = c(0, 0.4)) +
                    theme(legend.position = "none")
            }

            p <- p2 / p1 + plot_layout(heights = c(0.2, 1))

            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/%s_methylation_conditionaveraged.pdf", params$mod_code, element_type, element), 5, 5)
        }
    }
}


l1hsintactmethdf <- l1hsintactmethgr %>%
    as.data.frame() %>%
    tibble()

for (gene_id in l1hsintactmethdf %$% gene_id %>% unique()) {
    tryCatch(
        {
            pf <- l1hsintactmethdf %>%
                filter(gene_id == !!gene_id) %>%
                filter(cov > MINIMUMCOVERAGE) %>%
                group_by(sample) %>%
                mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
                filter(!is.na(rM)) %>%
                ungroup()
            p <- pf %>% ggplot() +
                geom_line(aes(x = start, y = rM, color = condition)) +
                scale_x_continuous(breaks = scales::breaks_pretty(3)) +
                facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
                ylim(c(0, 100)) +
                mtopen +
                scale_conditions
            dir.create(sprintf("ldna/results/%s/plots/rte/l1hsintact", params$mod_code))
            png(paste0(sprintf("ldna/results/%s/plots/rte/l1hsintact/", params$mod_code), gene_id, ".png"), 8, 3, units = "in", res = 300)
            print(p)
            dev.off()
        },
        error = function(e) {

        }
    )
}



# heatmap 5UTR
heatmapprep <- l1hsintactmethdf %>%
    filter(case_when(
        rte_strand == "+" ~ (start > rte_start) & (start < rte_start + 909),
        rte_strand == "-" ~ (start > rte_end - 909) & (start < rte_end)
    )) %>%
    group_by(gene_id, sample) %>%
    summarise(mean = mean(pctM)) %>%
    pivot_wider(names_from = sample, values_from = mean) %>%
    arrange(gene_id)
m <- as.matrix(heatmapprep %>% ungroup() %>% dplyr::select(-gene_id))
rownames(m) <- heatmapprep %$% gene_id
remove_rows <- rep(FALSE, nrow(m))
for (i in 1:nrow(m)) {
    print(m[i, ])
    if (any(is.na(m[i, ]))) {
        remove_rows[i] <- TRUE
    }
}
m <- m[!remove_rows, ]

rownames(m)
library(circlize)
col_fun <- colorRamp2(c(50, 75, 100), c("red", "white", "blue"))
col_fun(seq(50, 100, by = 12.5))

if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    l1hsintactdf <- l1hsintactmethdf %>%
        group_by(gene_id) %>%
        summarise(concordance = dplyr::first(concordance), genic_loc = dplyr::first(genic_loc))
    pvals <- l1hsintactdf %>%
        arrange(gene_id) %>%
        filter(gene_id %in% rownames(m)) %$% concordance
    length(pvals)
    is_sig <- !is.na(pvals)
    pch <- rep("*", length(pvals))
    pch[!is_sig] <- NA
    genic_locs <- l1hsintactdf %>%
        arrange(gene_id) %>%
        filter(gene_id %in% rownames(m)) %$% genic_loc
    row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), genic_loc = genic_locs, col = list(genic_loc = c("Genic" = "brown", "Intergenic" = "tan")))
    conditions <- c(sample_table %>% filter(condition == condition1) %$% condition, sample_table %>% filter(condition == condition2) %$% condition)

    conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
    topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))


    heatmapL1UTR <<- m %>%
        Heatmap(
            name = "CpG Methylation",
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_rot = 45,
            col = col_fun,
            split = pvals,
            top_annotation = topAnn,
            right_annotation = row_ha,
            row_title = "Intact L1HS"
        )

    col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
    col_fun(seq(50, 100, by = 12.5))

    heatmapL1UTR2 <- m %>%
        Heatmap(
            name = "CpG Methylation",
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_rot = 45,
            split = pvals,
            col = col_fun,
            top_annotation = topAnn,
            right_annotation = row_ha,
            row_title = "Intact L1HS"
        )
} else {
    heatmapL1UTR <<- m %>%
        Heatmap(
            name = "CpG Methylation",
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_rot = 45,
            col = col_fun,
            row_title = "Intact L1HS"
        )

    col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
    col_fun(seq(50, 100, by = 12.5))

    heatmapL1UTR2 <- m %>%
        Heatmap(
            name = "CpG Methylation",
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_rot = 45,
            col = col_fun,
            row_title = "Intact L1HS"
        )
}

p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR, heatmap_legend_side = "right", annotation_legend_side = "right")))
mysaveandstore(sprintf("ldna/results/%s/plots/l1intactheatmap_5utr.pdf", params$mod_code), 7, 14)


p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR2, heatmap_legend_side = "right", annotation_legend_side = "right")))
mysaveandstore(sprintf("ldna/results/%s/plots/l1intactheatmap_5utr_fullrange.pdf", params$mod_code), 7, 14)
# plgrob <- grid.grabExpr(ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right"))
# plots[["l1intactheatmap_5utr"]] <- plgrob

outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)
cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()

read_analysis2 <- function(
    readscg,
    cg_indices,
    mod_code_var = "m",
    regions_of_interest = list(c(0, 328), c(0, 500), c(0, 909), c(400, 600)),
    required_fraction_of_total_cg = 0.75,
    meth_thresholds = c(0.1, 0.25, 0.5),
    context = "CpG") {
    region <- "L1HS_FL"
    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, region, required_fraction_of_total_cg)
    dir.create(outputdirtables, recursive = TRUE)

    readsdf1 <- readscg %>%
        left_join(rmann %>%
            dplyr::select(gene_id, start, end, strand, rte_length_req, intactness_req) %>%
            dplyr::rename(element_strand = strand, element_start = start, element_end = end)) %>%
        filter(rte_length_req == "FL")


    by_cpg_l <- list()
    by_read_l <- list()
    by_sample_l <- list()
    by_gene_id_l <- list()

    for (region_of_interest in regions_of_interest) {
        roistart <- region_of_interest[1]
        roiend <- region_of_interest[2]
        roistring <- paste0(roistart, "to", roiend)

        numCGneeded <- ceiling(length(cg_indices[(cg_indices <= roiend) & (cg_indices >= roistart)]) * required_fraction_of_total_cg)

        utr1 <- readsdf1 %>%
            filter(mod_code == mod_code_var) %>%
            filter(case_when(
                element_strand == "+" ~ (start > element_start + roistart) & (start < element_start + roiend),
                element_strand == "-" ~ (start > element_end - roiend) & (start < element_end - roistart)
            )) %>%
            dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

        by_cpg_temp <- utr1 %>%
            group_by(gene_id, read_id, condition, sample) %>%
            mutate(num_cpgs_in_read = n()) %>%
            mutate(fraction_meth = mean(mod_indicator)) %>%
            relocate(gene_id) %>%
            ungroup()

        by_read_temp <- by_cpg_temp %>%
            filter(num_cpgs_in_read >= numCGneeded) %>%
            group_by(read_id, gene_id, sample, condition, region) %>%
            summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read), strand = dplyr::first(element_strand), numCGneeded = dplyr::first(numCGneeded)) %>%
            ungroup()

        by_cpg_temp$subset <- as.character(roistring)
        by_read_temp$subset <- as.character(roistring)

        by_cpg_l[[as.character(roistring)]] <- by_cpg_temp
        by_read_l[[as.character(roistring)]] <- by_read_temp

        for (meth_threshold in meth_thresholds) {
            subset_threshold <- paste0(roistring, "_", meth_threshold)
            by_sample_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth > meth_threshold, 0, 1)) %>%
                group_by(sample, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                group_by(condition, region) %>%
                mutate(meanProp = mean(propUnmeth))
            by_sample_temp$meth_threshold <- meth_threshold
            by_sample_temp$subset <- as.character(roistring)
            by_sample_temp$subset_threshold <- subset_threshold

            by_gene_id_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth >= meth_threshold, 0, 1)) %>%
                group_by(sample, gene_id, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                ungroup() %>%
                group_by(gene_id) %>%
                mutate(group_size = n()) %>%
                ungroup()
            by_gene_id_temp$meth_threshold <- meth_threshold
            by_gene_id_temp$subset <- as.character(roistring)
            by_gene_id_temp$subset_threshold <- subset_threshold

            by_sample_l[[subset_threshold]] <- by_sample_temp
            by_gene_id_l[[subset_threshold]] <- by_gene_id_temp
        }
    }

    by_cpg <- purrr::reduce(by_cpg_l, bind_rows)
    by_read <- purrr::reduce(by_read_l, bind_rows)
    by_gene_id <- purrr::reduce(by_gene_id_l, bind_rows)
    by_sample <- purrr::reduce(by_sample_l, bind_rows)

    dir.create(dirname(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)), recursive = TRUE)
    by_cpg %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region))
    by_read %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region))
    by_gene_id %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region))
    by_sample %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region))

    # by_cpg <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_read <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_gene_id <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_sample <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region)) %>%
    #         mutate(condition = factor(condition, levels = conf$levels)) %>%
    #         mutate(sample = factor(sample, levels = conf$samples))

    get_read_ecdf <- function(df, subset_val, breakpoints, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by specified variable
            summarise(
                percent_below = list(ecdf(fraction_meth)(breakpoints) * 100), # Compute ECDF
                .groups = "drop"
            ) %>%
            unnest_longer(percent_below) %>%
            mutate(
                threshold = rep(breakpoints, times = n() / length(breakpoints)), # Expand breakpoints
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), threshold, percent_below) # Keep relevant columns
    }


    get_read_quantiles <- function(df, subset_val, probs, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by sample_name or another variable
            summarise(
                quantiles = list(quantile(fraction_meth, probs = probs)), # Store as list
                .groups = "drop"
            ) %>%
            unnest_longer(quantiles) %>%
            mutate(
                quantile = rep(probs, times = n() / length(probs)), # Expand probs for each group
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), quantile, mean_meth = quantiles) # Keep relevant columns
    }
    subsets <- c("0to909", "0to500", "0to328")
    breakpoints <- seq(0, 1, 0.05)

    ecdf_reads <- map_dfr(subsets, ~ get_read_ecdf(by_cpg, .x, breakpoints, "sample"))
    quantile_reads <- map_dfr(subsets, ~ get_read_quantiles(by_cpg, .x, breakpoints, "sample"))


    ecdf_reads %>% write_mycsv(sprintf("%s/read_ecdf.csv", outputdirtables))
    quantile_reads %>% write_mycsv(sprintf("%s/read_quantiles.csv", outputdirtables))

    ecdf_reads_acrosssamplemean <- ecdf_reads %>%
        group_by(subset, threshold) %>%
        summarise(percent_below = mean(percent_below))
    quantile_reads_acrosssamplemean <- quantile_reads %>%
        group_by(subset, quantile) %>%
        summarise(mean_meth = mean(mean_meth))

    ecdf_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_ecdf_acrosssamplemean.csv", outputdirtables))
    quantile_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_quantiles_acrosssamplemean.csv", outputdirtables))


    # Example ECDF plot function
    plot_ecdf <- function(ecdf_data) {
        ggplot(ecdf_data, aes(x = threshold, y = percent_below, color = as.factor(subset))) +
            geom_line() + # ECDF curve
            geom_point() + # Add points for clarity
            labs(
                title = "Empirical CDF of Fraction Methylation",
                x = "Fraction Methylation Threshold",
                y = "Percent Below Threshold",
                color = "Subset"
            ) +
            theme_minimal()
    }


    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)


    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        # stats <- by_sample %>%
        #     filter(subset != "400to600") %>%
        #     ungroup() %>%
        #     left_join(sample_table) %>%
        #     group_split(subset_threshold) %>%
        #     set_names(unique(by_sample %>% filter(subset != "400to600") %$% subset_threshold)) %>%
        #     map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
        #     imap_dfr(~ .x %>% mutate(subset = .y))
        stats_list <- list()
        i <- 1
        for (subset in unique(by_cpg$subset)) {
            for (threshold in unique(by_sample$meth_threshold)) {
                by_cpg_tmp <- by_cpg %>%
                    filter(subset == !!subset) %>%
                    mutate(unmeth = ifelse(fraction_meth > threshold, 0, 1)) %>%
                    dplyr::rename(sample_name = sample) %>%
                    group_by(sample_name, condition, gene_id) %>%
                    summarise(unmeth = sum(unmeth), total = n()) %>%
                    ungroup() %>%
                    left_join(sample_table) %>%
                    mutate(age_z = as.numeric(scale(age)))
                model_tmp <- glmmTMB(
                    cbind(unmeth, total - unmeth) ~
                        condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
                    data = by_cpg_tmp,
                    family = binomial()
                )
                res_tmp <- broom::tidy(model_tmp) %>%
                    mutate(threshold = !!threshold) %>%
                    mutate(subset = !!subset)
                stats_list[[i]] <- res_tmp
                i <- i + 1
            }
        }
        stats <- purrr::reduce(stats_list, bind_rows)
        stats_padj <- stats %>%
            filter(term == "conditionAD") %>%
            filter(subset != "400to600") %>%
            mutate(padj = p.adjust(p.value, method = "fdr"))
        statswpadj <- stats %>% left_join(stats_padj)
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = statswpadj)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }
    p <- by_sample %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        stats <- by_sample %>%
            ungroup() %>%
            left_join(sample_table) %>%
            group_split(subset_threshold) %>%
            set_names(unique(by_sample$subset_threshold)) %>%
            map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
            imap_dfr(~ .x %>% mutate(subset = .y))
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = stats)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }

    p <- by_read %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)


    # Get group-wise mean fraction_meth
    # beta_overlays <- by_read %>%
    # filter(subset != "400to600") %>%
    # group_by(condition, subset) %>%
    # summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop")

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == "AD"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == "CTRL"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))
    p <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == "AD"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == "CTRL"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay_L1HS_4q28.3_9.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        filter(subset == "328") %$% condition %>%
        table()


    p <- by_gene_id %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        filter(meth_threshold == "0.5") %>%
        filter(subset == "0to909") %>%
        group_by(gene_id) %>%
        mutate(samples_detected_per_element = n()) %>%
        ungroup() %>%
        filter(samples_detected_per_element > 11) %>%
        distinct() %>%
        # group_by(direction) %>%
        tidyHeatmap::heatmap(gene_id, sample, propUnmeth,
            cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE,
            show_row_dend = FALSE,
            palette_value = circlize::colorRamp2(
                c(0, 0.00001, seq(0.1, 1, length.out = 3)),
                c("black", rev(RColorBrewer::brewer.pal(4, "Oranges")))
            )
        )
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_gene_id_consistent_across_samples.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)

    #
    df_filtered <- by_gene_id %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        filter(meth_threshold == "0.5", subset == "0to909") %>%
        group_by(gene_id) %>%
        mutate(samples_detected_per_element = n()) %>%
        ungroup() %>%
        filter(samples_detected_per_element > 11) %>%
        distinct()
    df_sorted <- df_filtered %>%
        group_by(sample) %>%
        arrange(propUnmeth) %>%
        mutate(gene_rank = row_number()) %>%
        ungroup()
    p <- ggplot(df_sorted, aes(x = sample, y = gene_rank, fill = propUnmeth)) +
        geom_tile() +
        scale_fill_gradientn(
            colours = c("lightblue", rev(RColorBrewer::brewer.pal(4, "Oranges"))),
            values = rescale(c(0, 0.00001, seq(0.1, 1, length.out = 3))),
            name = "propUnmeth"
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) + # No padding; bottom = high rank
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(x = "sample", y = "genes (sorted by propUnmeth per sample)") +
        mtclosed +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p, raster = TRUE)

    #

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    named_group_split <- function(.tbl, ...) {
        grouped <- group_by(.tbl, ...)
        names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

        grouped %>%
            group_split() %>%
            rlang::set_names(names)
    }

    groups_needed <- length(condition1samples) + 1
    groups_needed <- 10
    n_under_consideration <- by_gene_id %>%
        ungroup() %>%
        mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
        left_join(sample_table) %>%
        filter(group_size >= groups_needed) %>%
        group_by(subset_threshold) %>%
        dplyr::select(subset_threshold, gene_id) %>%
        distinct() %>%
        summarise(n_under_consideration = n())
    group_size_df <- by_gene_id %>% dplyr::select(gene_id, group_size, subset, meth_threshold)

    dat_tmp <- meth_thresholds %>%
        map(~ by_cpg %>%
            mutate(unmeth = ifelse(fraction_meth > .x, 0, 1), meth_threshold = .x) %>%
            filter(subset != "400to600") %>%
            group_by(sample, gene_id, subset) %>%
            mutate(total_read = n()) %>%
            ungroup() %>%
            group_by(sample, gene_id, subset, meth_threshold) %>%
            summarise(unmeth = sum(unmeth), total_read = dplyr::first(total_read)) %>%
            ungroup()) %>%
        list_rbind() %>%
        dplyr::rename(sample_name = sample) %>%
        mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
        left_join(group_size_df %>% distinct()) %>%
        left_join(sample_table) %>%
        filter(group_size >= groups_needed) %>%
        mutate(age_z = as.numeric(scale(age)))




    # Create an empty list to hold tidy model results
    stats_list <- list()

    # Split data by group
    dat_groups <- named_group_split(dat_tmp, split_var)
    library(broom.mixed)
    # Loop through each group
    for (i in seq_along(dat_groups)) {
        print(i)
        group_data <- dat_groups[[i]]
        group_name <- names(dat_groups)[i]

        print(group_name)

        # Check if unmeth values are constant (e.g. all 0s or all the same)
        if (length(unique(group_data$unmeth)) == 1 | sum(group_data$unmeth > 0) < 2) {
            print("no info")
            tidy_model <- tibble(
                term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                subset = group_name,
                note = "No variation in unmeth"
            )
        } else {
            model <- glmmTMB(
                cbind(unmeth, total_read - unmeth) ~
                    condition + sex + age_z + (1 | sample_name),
                data = group_data,
                family = binomial()
            )
            model <- glmmTMB(
                cbind(unmeth, total_read - unmeth) ~
                    condition + sex + age_z + (1 | sample_name),
                data = group_data,
                family = binomial()
            )

            tidy_model <- broom::tidy(model) %>%
                mutate(subset = group_name)
        }
        stats_list[[i]] <- tidy_model
    }

    # Combine the results into a single data frame
    stats <- purrr::reduce(stats_list, bind_rows) %>%
        tidyr::separate(subset, into = c("gene_id", "subset", "meth_threshold"), sep = "/", convert = TRUE)


    write_csv(stats, sprintf("ldna/results/%s/tables/reads_new/%s_%s/by_gene_new.csv", params$mod_code, region, required_fraction_of_total_cg))

    gene_condition_stats <- stats %>%
        filter(term == "conditionAD") %>%
        mutate(p.value = case_when(
            is.nan(p.value) ~ 1,
            TRUE ~ p.value
        )) %>%
        mutate(statistic = case_when(
            is.nan(statistic) ~ 0,
            TRUE ~ statistic
        ))

    tryCatch(
        {
            p <- stats %>%
                mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                filter(p.value <= 0.05) %>%
                filter(grepl("condition", term)) %>%
                mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                count(meth_threshold, subset, dir_stat) %>%
                complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                geom_col(position = "dodge", color = "black") +
                facet_wrap(vars(subset), nrow = 1) +
                ggtitle(sprintf("Significant Loci")) +
                labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                mtclosed +
                anchorbar +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 5.5, 4, pl = p)

            p <- stats %>%
                filter(subset != "400to600") %>%
                mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                filter(p.value <= 0.05) %>%
                filter(grepl("condition", term)) %>%
                mutate(meth_threshold = factor(as.character(meth_threshold), levels = c("0.1", "0.25", "0.5"))) %>%
                mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                count(meth_threshold, subset, dir_stat) %>%
                complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                geom_col(position = "dodge", color = "black") +
                facet_wrap(vars(subset), nrow = 1) +
                ggtitle(sprintf("Significant Loci")) +
                labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                mtclosed +
                anchorbar +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 4, pl = p)
        },
        error = function(e) {
            print("no sig de")
        }
    )

    dispersion_models <- list()
    for (subsetofinterest in by_cpg$subset %>%
        unique() %>%
        grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
        dispersion_model <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = by_cpg %>% filter(subset == subsetofinterest) %>%
                mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
        )
        dispersion_models[[subsetofinterest]] <- dispersion_model
    }

    res_text <- capture.output(map(dispersion_models, summary))
    writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))
    summary_disp <- summary(mod)$dispersion

    # # Extract estimates and standard errors
    # log_phi_control <- 1.95678
    # se_log_phi_control <- 0.01539

    # log_phi_diff_AD <- -0.32685
    # se_log_phi_diff_AD <- 0.01990

    # # Compute log(φ) for AD and its SE
    # log_phi_AD <- log_phi_control + log_phi_diff_AD
    # se_log_phi_AD <- sqrt(se_log_phi_control^2 + se_log_phi_diff_AD^2)

    # # 95% CI on log-scale
    # z <- 1.96
    # ci_log_phi_control <- log_phi_control + c(-1, 1) * z * se_log_phi_control
    # ci_log_phi_AD <- log_phi_AD + c(-1, 1) * z * se_log_phi_AD

    # # Exponentiate to get CIs for φ
    # phi_control <- exp(log_phi_control)
    # phi_AD <- exp(log_phi_AD)

    # ci_phi_control <- exp(ci_log_phi_control)
    # ci_phi_AD <- exp(ci_log_phi_AD)


    # dispersion_models_withreadid <- list()
    # for (subsetofinterest in by_cpg$subset %>%
    #     unique() %>%
    #     grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
    #     dispersion_model <- glmmTMB(
    #         cbind(meth, unmeth) ~ 1 + (1 | read_id),
    #         dispformula = ~condition,
    #         family = glmmTMB::betabinomial(),
    #         data = by_cpg %>% filter(subset == subsetofinterest) %>%
    #             mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
    #     )
    #     dispersion_models[[subsetofinterest]] <- dispersion_model
    # }

    # res_text <- capture.output(map(dispersion_models, summary))
    # writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))


    # dat <- by_cpg %>%
    #     filter(subset == "0to909") %>%
    #     mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
    # mod0 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
    #     data = dat, family = binomial()
    # )

    # mod05 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = dat, family = glmmTMB::betabinomial()
    # )
    # mod05disp <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     dispformula = ~1,
    #     data = dat, family = glmmTMB::betabinomial()
    # )
    # datmorebinom <- dat %>%
    #     filter(fraction_meth < 0.95) %>%
    #     filter(fraction_meth > 0.35)
    # mod05morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = datmorebinom, family = glmmTMB::betabinomial()
    # )
    # mod05binom_morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = datmorebinom, family = binomial()
    # )


    # dispersion_model_morebinomdat <- glmmTMB(
    #     cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
    #     dispformula = ~condition,
    #     family = glmmTMB::betabinomial(),
    #     data = datmorebinom
    # )

    # anova(mod0, mod05)
    # anova(mod05binom_morebinom, mod05morebinom)
}


read_analysis2(readscg, cg_indices)


##################
grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
grsdf %$% sample %>% unique()
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
grs <- GRanges(grsdf)

grs_cpg_islands <- grs %>% subsetByOverlaps(cpg_islands)
grs_cpg_islands$islandStatus <- "island"
grs_cpgi_shelves <- grs %>% subsetByOverlaps(cpgi_shelves)
grs_cpgi_shelves$islandStatus <- "shelf"
grs_cpgi_shores <- grs %>% subsetByOverlaps(cpgi_shores)
grs_cpgi_shores$islandStatus <- "shore"
grs_cpg_opensea <- grs %>% subsetByOverlaps(cpgi_features, invert = TRUE)
grs_cpg_opensea$islandStatus <- "opensea"
# SETTING UP SOME SUBSETS FOR EXPLORATION
set.seed(75)
possample <- sample(grsdf$pos, size = 1000, replace = FALSE)
grsdfs <- grsdf %>% filter(pos %in% possample)
grss <- GRanges(grsdfs)


grsdf %>%
    ungroup() %>%
    group_by(pos, condition) %>%
    summarise(variance = var(pctM)) %>%
    group_by(condition) %>%
    summarise(meanvar = mean(variance))



grsdfs %>%
    ungroup() %>%
    dplyr::select(pos, pctM, sample) %>%
    pivot_wider(id_cols = pos, names_from = sample, values_from = pctM) %>%
    na.omit()
############
# GLOBAL
dir.create(sprintf("ldna/results/%s/plots/genomewide", params$mod_code), showWarnings = FALSE)

grsdfsummary <- grsdf %>%
    mutate(methylated_sites = round(cov * pctM / 100)) %>%
    group_by(sample) %>%
    summarise(total_cov = sum(cov), total_meth = sum(methylated_sites))

pf <- grsdfsummary %>%
    mutate(mean_meth = total_meth / total_cov) %>%
    left_join(sample_table)
pf %>%
    group_by(condition) %>%
    summarise(mean_meth = mean(mean_meth))

p <- pf %>%
    ggplot(aes(y = sample, x = mean_meth, color = condition, shape = !!sym(asc))) +
    geom_point(size = 3) +
    scale_conditions +
    geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
    geom_text_repel(aes(label = apoe)) +
    # new_scale_fill() +
    # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_conditions +
    mtopen
library(broom)
if (enough_samples_per_condition_for_stats) {
    # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
    genomewide_model <- glmmTMB(
        cbind(total_meth, total_cov - total_meth) ~
            condition + sex + age_z + (1 | sample),
        data = grsdfsummary %>% left_join(sample_table %>%
            mutate(age_z = as.numeric(scale(age)))),
        family = binomial()
    )
    stats <- broom::tidy(genomewide_model)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4)
}

tryCatch(
    {
        p <- pf %>%
            mutate(sample = fct_reorder(paste0(sample, "_", age), age)) %>%
            ggplot(aes(y = sample, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +
            # new_scale_fill() +
            # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
            scale_conditions +
            mtopen
        library(broom)
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_withage_ordered_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
    },
    error = function(e) {

    }
)

# t.test(pctM ~ condition, data = grsdf, var.equal = TRUE)
# t.test(pctM ~ condition, data = grsdf %>% filter(seqnames %in% chromosomesNoX), var.equal = TRUE)

p <- grsdfs %>% ggplot() +
    geom_boxplot(aes(x = islandStatus, y = pctM, fill = condition), outliers = FALSE) +
    mtopen +
    scale_conditions
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/cpgislandstatusbox.pdf", params$mod_code), w = 5, h = 4, res = 300, pl = p)


p <- grsdfs %>%
    group_by(islandStatus, condition) %>%
    summarize(pctM = mean(pctM)) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = pctM, fill = condition), position = "dodge", color = "black") +
    mtopen +
    scale_conditions +
    anchorbar
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/cpgislandstatusbar_1000.pdf", params$mod_code), w = 4, h = 4, res = 300, pl = p)



if ((conf$single_condition == "no") & enough_samples_per_condition_for_stats) {
    ######### GENES
    {
        directions <- c("Hypo", "Hyper", "Dif")
        mydir <- sprintf("ldna/results/%s/plots/great", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great", params$mod_code)
        dir.create(mydir, recursive = TRUE, showWarnings = FALSE)
        dir.create(mydirtables, recursive = TRUE, showWarnings = FALSE)
        # convert to entrez
        # symbols <- mcols(gr)$Name
        # library(org.Hs.eg.db)
        # columns(org.Hs.eg.db)
        # anno.result <- mapIds(org.Hs.eg.db,
        #     keys = symbols,
        #     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
        # )
        # mask <- !is.na(mcols(genes))
        # genes_notNA <- genes[mask]
        et <- extendTSS(genes_gr %>% resize(width = 1), genome_lengths, gene_id_type = "SYMBOL")
        et_noextension <- extendTSS(genes_gr %>% resize(width = 1), genome_lengths, gene_id_type = "SYMBOL", extension = 0)
    }



    {
        library(clusterProfiler)
        library(msigdbr)
        gs <- msigdbr("human")
        # gs %>% filter(grepl("GOBP_HISTONE_H3_K27_TRI", gs_name))
        get_gs_enrichments <- function(gs, gs_ontology_level, outputdir, regionsgrs, etparam, et_mode_string, directions = c("Hypo", "Hyper", "Dif"), background = NULL) {
            outputdirplots <- file.path(outputdir, gs_ontology_level, et_mode_string)
            outputdirtables <- file.path(gsub("m/plots/great", "m/tables/great", outputdir), gs_ontology_level, et_mode_string)
            print(outputdirplots)
            tablesMsigdb <- list()
            genecollections <- gs %>%
                pluck(gs_ontology_level) %>%
                unique()
            for (collection in genecollections) {
                print(collection)
                tryCatch(
                    {
                        dir.create(paste(outputdirplots, collection, sep = "/"), recursive = TRUE)
                        dir.create(paste(outputdirtables, collection, sep = "/"), recursive = TRUE)
                        genesets <- gs %>%
                            filter(!!sym(gs_ontology_level) == collection) %>%
                            dplyr::select(gs_name, gene_symbol) %>%
                            dplyr::rename(term = gs_name, gene = gene_symbol) %>%
                            as.data.frame()
                        for (direction in directions) {
                            if (direction == "Hypo") {
                                regionstemp <<- regionsgrs[grepl("Hypo", regionsgrs$direction)]
                            }
                            if (direction == "Hyper") {
                                regionstemp <<- regionsgrs[grepl("Hyper", regionsgrs$direction)]
                            }
                            if (direction == "Dif") {
                                regionstemp <<- regionsgrs
                            }

                            if (is.null(background)) {
                                res <- great(regionstemp, gene_sets = genesets, extended_tss = etparam, background = CHROMOSOMESINCLUDEDINANALYSIS_REF)
                            } else {
                                res <- great(regionstemp, gene_sets = genesets, extended_tss = etparam, background = background)
                            }
                            tb <- getEnrichmentTable(res)
                            if (nrow(tb) != 0) {
                                tb <- tb %>% dplyr::arrange(p_adjust)
                                tablesMsigdb[[collection]][[direction]] <- tb
                                write_delim(tb, paste(outputdirtables, collection, paste0(direction, "great_enrichment.tsv"), sep = "/"))

                                png(paste(outputdirplots, collection, paste0(direction, "volcano.png"), sep = "/"), height = 5, width = 5, res = 300, units = "in")
                                plotVolcano(res)
                                dev.off()

                                png(paste(outputdirplots, collection, paste0(direction, "associations.png"), sep = "/"), height = 5, width = 10, res = 300, units = "in")
                                plotRegionGeneAssociations(res)
                                dev.off()
                            }
                        }
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
            # save(file = sprintf("ldna/Rintermediates/%s/tablesMsigdb_%s_%s.rds", params$mod_code, gs_ontology_level, et_mode_string), tablesMsigdb)
            save(file = sprintf("%s/tablesMsigdb.rds", outputdirtables), tablesMsigdb)
            return(tablesMsigdb)
        }

        make_enrich_plots <- function(tempdirectory) {
            tsv_files <- list.files(
                path = tempdirectory,
                pattern = "\\.tsv$",
                recursive = TRUE,
                full.names = TRUE
            )
            for (tempfile in tsv_files) {
                tb <- read_delim(tempfile)
                outputdirplots <- gsub("tables", "plots", dirname(dirname(tempfile)))
                collection <- basename(dirname(tempfile))
                direction <- gsub("great_enrichment.tsv", "", basename(tempfile))

                tbnames <- tb %>%
                    tibble() %>%
                    mutate(id_nchar = nchar(id)) %>%
                    mutate(id = case_when(
                        id_nchar < 40 ~ paste0(strrep("-", pmax(0, 40 - id_nchar)), id),
                        TRUE ~ id
                    )) %>%
                    mutate(mean_padj = (p_adjust + p_adjust_hyper) / 2)

                binom_df <- tbnames %>%
                    dplyr::select(id, fold_enrichment, p_adjust, mean_padj) %>%
                    mutate(type = "Binom")
                hyper_df <- tbnames %>%
                    dplyr::select(id, fold_enrichment = fold_enrichment_hyper, p_adjust = p_adjust_hyper, mean_padj) %>%
                    mutate(type = "Hyper")
                combined_df <- bind_rows(
                    binom_df,
                    hyper_df %>% mutate(fold_enrichment)
                )

                terms_to_plot <- tbnames %>%
                    arrange(mean_padj) %>%
                    head(n = 5) %$% id
                p <- combined_df %>%
                    filter(id %in% terms_to_plot) %>%
                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                    mutate(id = fct_reorder(id, -mean_padj)) %>%
                    ggplot(aes(x = id)) +
                    geom_col(data = . %>% filter(type == "Binom"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position_nudge(x = 0.45 / 2), width = 0.45) +
                    coord_flip() +
                    scale_fill_distiller(
                        name = "BinomP",
                        palette = "Blues",
                        direction = -1,
                        limits = c(0, 1),
                        oob = scales::squish
                    ) +
                    new_scale_fill() +
                    geom_col(data = . %>% filter(type == "Hyper"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position_nudge(x = -0.45 / 2), width = 0.45) +
                    scale_fill_distiller(
                        name = "HyperP",
                        palette = "Greens",
                        direction = -1,
                        limits = c(0, 1),
                        oob = scales::squish
                    ) +
                    geom_vline(xintercept = 0, color = "black") +
                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                    mtclosedgridv +
                    anchorbar +
                    geom_hline(yintercept = 0, color = "black") +
                    theme(axis.text.y = element_text(family = "mono"))
                mysaveandstore(pl = p, fn = paste(outputdirplots, collection, paste0(direction, "lollipop5.pdf"), sep = "/"), 6, 4)

                terms_to_plot <- tbnames %>%
                    arrange(mean_padj) %>%
                    head(n = 10) %$% id
                p <- combined_df %>%
                    filter(id %in% terms_to_plot) %>%
                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                    mutate(id = fct_reorder(id, -mean_padj)) %>%
                    ggplot(aes(x = id)) +
                    geom_col(data = . %>% filter(type == "Binom"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position_nudge(x = 0.45 / 2), width = 0.45) +
                    coord_flip() +
                    scale_fill_distiller(
                        name = "BinomP",
                        palette = "Blues",
                        direction = -1,
                        limits = c(0, 1),
                        oob = scales::squish
                    ) +
                    new_scale_fill() +
                    geom_col(data = . %>% filter(type == "Hyper"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position_nudge(x = -0.45 / 2), width = 0.45) +
                    scale_fill_distiller(
                        name = "HyperP",
                        palette = "Greens",
                        direction = -1,
                        limits = c(0, 1),
                        oob = scales::squish
                    ) +
                    geom_vline(xintercept = 0, color = "black") +
                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                    mtclosedgridv +
                    anchorbar +
                    geom_hline(yintercept = 0, color = "black") +
                    theme(axis.text.y = element_text(family = "mono"))

                mysaveandstore(pl = p, fn = paste(outputdirplots, collection, paste0(direction, "lollipop.pdf"), sep = "/"), 6.5, 6)
            }
        }

        #####
        # promoters but with background
        mydir <- sprintf("ldna/results/%s/plots/great_promoters", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great_promoters", params$mod_code)
        for (dmrtype in dmrs$dmr_type %>% unique()) {
            regions1 <- mergeByOverlaps(promoters, dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype])
            regions2 <- regions1$promoters
            mcols(regions2)$direction <- as.data.frame(regions1)$direction
            regions <- as.data.frame(regions2) %>%
                tibble() %>%
                distinct() %>%
                GRanges()


            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_cat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et_noextension,
                        et_mode_string = "et_noextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = promoters
                    )
                },
                error = function(e) {}
            )
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et_noextension,
                        et_mode_string = "et_noextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = promoters
                    )
                },
                error = function(e) {}
            )
        }

        ##### cpg island background
        mydir <- sprintf("ldna/results/%s/plots/great_cpgislands", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great_cpgislands", params$mod_code)
        for (dmrtype in dmrs$dmr_type %>% unique()) {
            regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_cat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = cpg_islands
                    )
                },
                error = function(e) {}
            )
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = cpg_islands
                    )
                },
                error = function(e) {}
            )
        }

        ####

        ##### promoters and enhancers background
        chromHMMgr <- import(conf$chromHMM)
        chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]
        prom_no_mcols <- promoters
        mcols(prom_no_mcols) <- NULL
        enh_no_mcols <- chromHMM_enhancers_grs
        mcols(enh_no_mcols) <- NULL
        background <- c(enh_no_mcols, prom_no_mcols)
        mydir <- sprintf("ldna/results/%s/plots/great_prom_enh", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great_prom_enh", params$mod_code)

        for (dmrtype in dmrs$dmr_type %>% unique()) {
            regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_cat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
        }

        make_enrich_plots("ldna/results/m/tables/great_prom_enh")



        ##### promoters and enhancers island background
        chromHMMgr <- import(conf$chromHMM)
        chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]
        prom_no_mcols <- promoters
        mcols(prom_no_mcols) <- NULL
        enh_no_mcols <- chromHMM_enhancers_grs
        mcols(enh_no_mcols) <- NULL
        prom_enh <- c(enh_no_mcols, prom_no_mcols)

        background <- cpg_islands %>% subsetByOverlaps(prom_enh)
        mydir <- sprintf("ldna/results/%s/plots/great_prom_enh_intersect_cpgI", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great_prom_enh_intersect_cpgI", params$mod_code)
        for (dmrtype in dmrs$dmr_type %>% unique()) {
            regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_cat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
        }
    }
    ####


    {
        grsdf_island_status <- grsdf %>%
            group_by(sample, islandStatus) %>%
            mutate(methylated_sites = round(cov * pctM / 100)) %>%
            summarise(methylated_sites = sum(methylated_sites), cov = sum(cov)) %>%
            mutate(pctM = methylated_sites / cov)

        dat_tmp <- grsdf_island_status %>%
            ungroup() %>%
            dplyr::rename(sample_name = sample) %>%
            left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))
        genomewide_model_islandstatus <- glmmTMB(
            cbind(methylated_sites, cov - methylated_sites) ~
                condition * islandStatus + sex + age_z + (1 | sample_name),
            data = dat_tmp,
            family = binomial()
        )
        stats <- broom::tidy(genomewide_model_islandstatus)

        p <- grsdf_island_status %>%
            mutate(pctM = methylated_sites / cov) %>%
            left_join(sample_table) %>%
            ggplot() +
            geom_point(aes(x = islandStatus, y = pctM, color = condition), position = position_dodge(width = 0.75)) +
            mtopen +
            scale_conditions +
            anchorbar
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_islandstatus_genomewide.pdf", params$mod_code), 5, 4, sf = stats, pl = p)
    }

    {
        island_promoters <- merge_with_grs(grs[mcols(grs)$islandStatus == "island"], promoters)

        ip <- island_promoters %>%
            group_by(gene_id, sample, condition) %>%
            summarise(mean_meth = mean(pctM))


        ip_mean <- island_promoters %>%
            group_by(sample, condition) %>%
            summarise(mean_meth = mean(pctM)) %>%
            left_join(sample_table %>% mutate(sample = sample_name, age_z = as.numeric(scale(age))))


        #### multilevel STATS
        dat_tmp <- island_promoters %>%
            mutate(methylated_sites = round(cov * pctM / 100), total_sites = cov) %>%
            dplyr::rename(sample_name = sample) %>%
            group_by(gene_id, sample_name) %>%
            summarise(methylated_sites = sum(methylated_sites), total_sites = sum(total_sites)) %>%
            left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))


        p <- ip_mean %>%
            mutate(condition = factor(condition, levels = conf$levels)) %>%
            ggviolin(x = "condition", y = "mean_meth", fill = "condition", add = c("mean_se", "dotplot")) +
            scale_conditions +
            scale_conditions +
            mtopen
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_violin.pdf", params$mod_code), 4.5, 3.75)

        p <- ip_mean %>%
            ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = ip_mean %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = ip_mean %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +

            # new_scale_fill() +
            # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
            scale_conditions +
            mtopen

        library(broom)
        if (enough_samples_per_condition_for_stats) {
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()

            gene_model <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
                data = dat_tmp,
                family = binomial()
            )
            stats <- broom::tidy(gene_model)
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
        } else {
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide.pdf", params$mod_code), 5, 3.7)
        }

        island_promotersnoX <- island_promoters %>% filter(seqnames != "chrX")

        ipnoX <- island_promotersnoX %>%
            group_by(gene_id, sample, condition) %>%
            summarise(mean_meth = mean(pctM))


        ip_meannoX <- island_promotersnoX %>%
            group_by(sample, condition) %>%
            summarise(mean_meth = mean(pctM)) %>%
            left_join(sample_table %>% mutate(sample = sample_name, age_z = as.numeric(scale(age))))


        #### multilevel STATS THIS UP
        dat_tmp <- island_promotersnoX %>%
            mutate(methylated_sites = round(cov * pctM / 100), total_sites = cov) %>%
            dplyr::rename(sample_name = sample) %>%
            group_by(gene_id, sample_name) %>%
            summarise(methylated_sites = sum(methylated_sites), total_sites = sum(total_sites)) %>%
            left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))


        p <- ip_meannoX %>%
            ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = ip_mean %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = ip_mean %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +
            scale_conditions +
            mtopen
        library(broom)
        if (enough_samples_per_condition_for_stats) {
            # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()

            gene_model <- glmmTMB(
                cbind(methylated_sites, total_sites - methylated_sites) ~
                    condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
                data = dat_tmp,
                family = binomial()
            )
            stats <- broom::tidy(gene_model)
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_noX.pdf", params$mod_code), 5, 4, sf = stats)
        } else {
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_noX.pdf", params$mod_code), 5, 3.7)
        }

        #####

        results <- ip %>%
            left_join(sample_table) %>%
            group_by(gene_id) %>%
            summarise(
                model = list(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), data = cur_data())),
                .groups = "drop"
            ) %>%
            mutate(
                tidied = map(model, broom::tidy)
            ) %>%
            unnest(tidied) %>%
            filter(term == "conditionCTRL") %>%
            dplyr::select(gene_id, estimate, p.value) %>%
            mutate(
                padj = p.adjust(p.value, method = "fdr") # Adjust p-values for multiple testing
            )

        res <- results %>%
            mutate(signed_log10p = sign(estimate) * abs(log10(p.value))) %>%
            arrange(-signed_log10p)
        ordered_by_stat <- setNames(res[["signed_log10p"]], res$gene_id) %>% na.omit()
        # GSEA untargeted
        tryCatch(
            {
                gene_sets <- msigdbr(species = confALL$aref$species)
            },
            error = function(e) {
                gene_sets <<- msigdbr(species = "human")
            }
        )
        library(clusterProfiler)
        rm(gse_df)

        for (category in gene_sets %$% gs_cat %>% unique()) {
            cat(category, "\n")
            tryCatch({
                collection <- category
                msigdbr_df <- gene_sets %>% filter(gs_cat == category)
                msigdbr_t2g <- msigdbr_df %>%
                    dplyr::distinct(gs_name, gene_symbol) %>%
                    as.data.frame()
                gse <- GSEA(ordered_by_stat, TERM2GENE = msigdbr_t2g, maxGSSize = 100000, minGSSize = 1)
                df <- gse@result %>% tibble()
                df$collection <- collection
                df$contrast <- str_glue("condition_{condition2}_vs_{condition1}")
                # gse_results[[contrast]][[collection]] <- as.data.frame(df) %>% tibble()
                if (!exists("gse_df")) {
                    gse_df <<- df
                } else {
                    gse_df <<- rbind(gse_df, df)
                }
                genesettheme <- theme_gray() + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

                tryCatch(
                    {
                        for (num in c(5, 10, 15)) {
                            dftemp <- arrange(df, -abs(NES)) %>%
                                group_by(sign(NES)) %>%
                                slice(1:num) %>%
                                mutate(log10padj = -log10(p.adjust))
                            dftemp <- dftemp %>% mutate(Description = str_wrap(as.character(Description) %>% gsub("_", " ", .), width = 40))
                            dftemp <- dftemp %>% mutate(`Gene Ratio` = 0.01 * as.numeric(gsub("%", "", gsub(",.*", "", gsub("tags=", "", leading_edge)))))
                            p <- ggplot(dftemp, aes(NES, fct_reorder(Description, NES), fill = log10padj)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient(high = "red", low = "white", limits = c(0, 5), oob = scales::squish) +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)

                            p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust) * `sign(NES)`)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                guides(fill = guide_legend(title = "Signed \n-log10(p.adjust)")) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )
            })
        }

        # ipcondition <- ip %>%
        #     group_by(gene_id, condition) %>%
        #     mutate(mean_meth = mean(mean_meth)) %>%
        #     group_by(gene_id, condition) %>%
        #     mutate(nrow = row_number()) %>%
        #     filter(nrow == 1) %>%
        #     dplyr::select(-nrow, -sample) %>%
        #     pivot_wider(id_cols = gene_id, names_from = condition, values_from = mean_meth) %>%
        #     mutate(dif = !!sym(condition2) - !!sym(condition1))


        # ipcondition %>%
        #     arrange(dif) %>%
        #     head(n = 30) %>%
        #     print(n = 30)
        # ipcondition %>%
        #     arrange(-dif) %>%
        #     head(n = 30) %>%
        #     print(n = 30)


        # ipwide <- ip %>%
        #     dplyr::select(-condition) %>%
        #     pivot_wider(names_from = sample, values_from = mean_meth)

        # ipwide %>% filter(gene_id == "LOC112268283")

        # library(dplyr)
        # library(broom)

        # results <- ip %>%
        #     group_by(gene_id) %>%
        #     summarise(
        #         model = list(lm(mean_meth ~ condition, data = cur_data())),
        #         .groups = "drop"
        #     ) %>%
        #     mutate(
        #         tidied = map(model, broom::tidy)
        #     ) %>%
        #     unnest(tidied) %>%
        #     filter(term == "conditionCTRL") %>%
        #     dplyr::select(gene_id, estimate, p.value) %>%
        #     mutate(
        #         log2FC = log2(exp(estimate)), # Convert from natural log to log2 fold-change
        #         padj = p.adjust(p.value, method = "BH") # Adjust p-values
        #     )

        # results %>% filter(p.value < 0.05)

        # ipwide %>% filter(gene_id == "AGBL4-AS1")

        # # View results
        # head(results)
        ## END NEW CODE



        # threshold_dfs <- list()
        # for (dmrtype in dmrs$dmr_type %>% unique()) {
        #     dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
        #     hyporegions <- dmrsgr_temp[grepl("Hypo", dmrsgr_temp$direction)]
        #     hyperregions <- dmrsgr_temp[grepl("Hyper", dmrsgr_temp$direction)]

        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, dmrsgr_temp, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyporegions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hypo_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyperregions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hyper_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")

        #     hypo <- mcols(subsetByOverlaps(promoters, hyporegions))$gene_id
        #     hyper <- mcols(subsetByOverlaps(promoters, hyperregions))$gene_id
        #     disc <- intersect(hypo, hyper)
        #     hypo_nd <- setdiff(hypo, disc)
        #     hyper_nd <- setdiff(hyper, disc)



        #     threshold_df <- bind_rows(
        #         tibble(gene_id = hypo_nd, !!sym(dmrtype) := "Hypo"),
        #         tibble(gene_id = hyper_nd, !!sym(dmrtype) := "Hyper"),
        #         tibble(gene_id = disc, !!sym(dmrtype) := "Discordant")
        #     )
        #     threshold_dfs[[dmrtype]] <- threshold_df
        # }
        # thresholddf <- purrr::reduce(threshold_dfs, full_join)
        # promoters_df <- full_join(tibble(as.data.frame(promoters)), thresholddf) %>%
        #     filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
        #     pivot_longer(cols = c("t05", "t01"), names_to = "dmr_type", values_to = "direction") %>%
        #     mutate(direction = factor(direction, levels = c("Discordant", "Hyper", "Hypo")))

        # total_possible <- nrow(mcols(promoters))

        # p <- promoters_df %>%
        #     group_by(dmr_type, direction) %>%
        #     summarise(n = n()) %>%
        #     ungroup() %>%
        #     tidyr::complete(dmr_type, direction, fill = list(n = 0)) %>%
        #     filter(!is.na(direction)) %>%
        #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        #     ggplot() +
        #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
        #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
        #     ggtitle("Promoter Methylation") +
        #     labs(x = "Direction", y = str_glue("Count (out of {total_possible})")) +
        #     theme(legend.position = "none") +
        #     annotation_logticks(sides = "l") +
        #     scale_y_log10(
        #         breaks = scales::trans_breaks("log10", function(x) 10^x),
        #         labels = scales::trans_format("log10", math_format(10^.x))
        #     ) +
        #     mtopen +
        #     scale_methylation_thresholds
        # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/genes_concordance.pdf", params$mod_code), 5, 4)

        # threshold_dfs <- list()
        # for (dmrtype in dmrs$dmr_type %>% unique()) {
        #     dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
        #     hyporegions <- dmrsgr_temp[grepl("Hypo", dmrsgr_temp$direction)]
        #     hyperregions <- dmrsgr_temp[grepl("Hyper", dmrsgr_temp$direction)]
        #     chromHMM_enhancers_grs_with_id <- chromHMM_enhancers_grs
        #     mcols(chromHMM_enhancers_grs_with_id)$ID <- paste0("ID", 1:nrow(mcols(chromHMM_enhancers_grs)))
        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, dmrsgr_temp, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, hyporegions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_hypo_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
        #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, hyperregions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_hyper_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")

        #     hypo <- mcols(subsetByOverlaps(chromHMM_enhancers_grs_with_id, hyporegions))$ID
        #     hyper <- mcols(subsetByOverlaps(chromHMM_enhancers_grs_with_id, hyperregions))$ID
        #     disc <- intersect(hypo, hyper)
        #     hypo_nd <- setdiff(hypo, disc)
        #     hyper_nd <- setdiff(hyper, disc)

        #     threshold_df <- bind_rows(
        #         tibble(ID = hypo_nd, !!sym(dmrtype) := "Hypo"),
        #         tibble(ID = hyper_nd, !!sym(dmrtype) := "Hyper"),
        #         tibble(ID = disc, !!sym(dmrtype) := "Discordant")
        #     )
        #     threshold_dfs[[dmrtype]] <- threshold_df
        # }
        # thresholddf <- purrr::reduce(threshold_dfs, full_join)
        # enh_df <- full_join(tibble(as.data.frame(chromHMM_enhancers_grs_with_id)), thresholddf) %>%
        #     filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
        #     pivot_longer(cols = c("t05", "t01"), names_to = "dmr_type", values_to = "direction") %>%
        #     mutate(direction = factor(direction, levels = c("Discordant", "Hyper", "Hypo")))

        # total_possible <- nrow(mcols(chromHMM_enhancers_grs))

        # p <- enh_df %>%
        #     group_by(dmr_type, direction) %>%
        #     summarise(n = n()) %>%
        #     ungroup() %>%
        #     tidyr::complete(dmr_type, direction, fill = list(n = 0)) %>%
        #     filter(!is.na(direction)) %>%
        #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        #     ggplot() +
        #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
        #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
        #     ggtitle("Enhancer Methylation") +
        #     labs(x = "Direction", y = str_glue("Count (out of {total_possible})")) +
        #     theme(legend.position = "none") +
        #     annotation_logticks(sides = "l") +
        #     scale_y_log10(
        #         breaks = scales::trans_breaks("log10", function(x) 10^x),
        #         labels = scales::trans_format("log10", math_format(10^.x))
        #     ) +
        #     mtopen +
        #     scale_methylation_thresholds
        # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/enhancer_concordance.pdf", params$mod_code), 5, 4)
        # #####
        # enh_df %>% bind_rows()
        # prom_enh_df <- promoters_df %>%
        #     mutate(biotype = "Promoter") %>%
        #     filter(!is.na(direction)) %>%
        #     dplyr::select(direction, biotype, dmr_type, ID) %>%
        #     bind_rows(enh_df %>% mutate(biotype = "Enhancer") %>% filter(!is.na(direction)) %>% dplyr::select(direction, dmr_type, biotype, ID))
        # p <- prom_enh_df %>%
        #     group_by(dmr_type, direction, biotype) %>%
        #     summarise(n = n()) %>%
        #     ungroup() %>%
        #     tidyr::complete(dmr_type, direction, biotype, fill = list(n = 0)) %>%
        #     filter(direction != "Discordant") %>%
        #     filter(!is.na(direction)) %>%
        #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        #     ggplot() +
        #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = biotype, y = n, fill = direction_threshold), color = "black", position = "dodge2") +
        #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = biotype, y = n, fill = direction_threshold), color = "black", position = "dodge2") +
        #     ggtitle("Enhancer Methylation") +
        #     labs(x = "Direction", y = str_glue("Count")) +
        #     theme(legend.position = "none") +
        #     annotation_logticks(sides = "l") +
        #     scale_y_log10(
        #         breaks = scales::trans_breaks("log10", function(x) 10^x),
        #         labels = scales::trans_format("log10", math_format(10^.x))
        #     ) +
        #     mtopen +
        #     scale_methylation_thresholds
        # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/prom_enhancer_concordance.pdf", params$mod_code), 5, 4)

        # # dmrsgr_enh <- dmrsgr %>% subsetByOverlaps(chromHMM_enhancers_grs)
        # # dmrsgr_prom <- dmrsgr %>% subsetByOverlaps(promoters)
        # # length_overlap <- dmrsgr_enh %>% subsetByOverlaps(dmrsgr_prom) %>% mcols() %>% nrow()
        # # perc_prom_in_enh <- length_overlap/nrow(mcols(dmrsgr_prom))
        # # perc_enh_in_prom <- length_overlap/nrow(mcols(dmrsgr_enh))
        # enh_in_dmrs <- chromHMM_enhancers_grs %>% subsetByOverlaps(dmrsgr)
        # prom_in_dmrs <- promoters %>% subsetByOverlaps(dmrsgr)
        # enh_intersect <- enh_in_dmrs %>%
        #     subsetByOverlaps(prom_in_dmrs) %>%
        #     mcols() %>%
        #     nrow()
        # perc_prom_in_enh <- enh_intersect / nrow(mcols(dmrsgr_prom))
        # prom_intersect <- prom_in_dmrs %>%
        #     subsetByOverlaps(enh_in_dmrs) %>%
        #     mcols() %>%
        #     nrow()
        # perc_enh_in_prom <- prom_intersect / nrow(mcols(dmrsgr_enh))

        #####


        # for (dmrtype in dmrs$dmr_type %>% unique()) {
        #     library(clusterProfiler)
        #     hypo_genes <- promoters_df %>% filter(direction == "Hypo") %$% gene_id
        #     hyper_genes <- promoters_df %>% filter(direction == "Hyper") %$% gene_id
        #     background <- mcols(promoters)$gene_id

        #     gs <- msigdbr("human")
        #     tablesORA_subcollection <- list()
        #     results_ORA_hypo <- list()
        #     results_ORA_hyper <- list()
        #     genesubcollections <- gs$gs_subcat %>% unique()
        #     for (collection in genesubcollections) {
        #         term2gene <- gs %>%
        #             filter(gs_subcat == collection) %>%
        #             dplyr::rename(term = gs_name, gene = gene_symbol) %>%
        #             select(term, gene)
        #         res_hypo <- enricher(hypo_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

        #         results_ORA_hypo[[collection]] <- res_hypo %>% as.data.frame()

        #         res_hyper <- enricher(hyper_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

        #         results_ORA_hyper[[collection]] <- res_hyper %>% as.data.frame()
        #     }



        ## DMR intersection
        {
            get_region_enrichment <- function(dmrstemp, rangestemp, genome_size = 3.1e9) {
                n_overlap <- dmrstemp %>%
                    subsetByOverlaps(rangestemp) %>%
                    as.data.frame() %>%
                    nrow()
                n_total <- dmrstemp %>%
                    as.data.frame() %>%
                    nrow()
                prop_genome <- (1 / genome_size) * width(rangestemp) %>% sum()
                fold_enrichment <- (n_overlap / n_total) / prop_genome
                p_value <- binom.test(n_overlap, n_total, p = prop_genome)$p.value
                return(c("fold_enrichment" = fold_enrichment, "p_value" = p_value, "n_overlap" = n_overlap, "n_total" = n_total, "prop_genome" = prop_genome))
            }

            # i.e hypodmrs05 and cpgislands
            # note subject must have name attribute
            get_region_enrichment2 <- function(querygrs, subjectgrs, genome_size = 3.1e9) {
                hits <- findOverlaps(querygrs, subjectgrs)

                query_hits <- queryHits(hits)
                subject_hits <- subjectHits(hits)

                overlaps <- pintersect(querygrs[query_hits], subjectgrs[subject_hits])
                overlap_widths <- width(overlaps)

                hitsdf <- data.frame(
                    query = query_hits,
                    subject = subject_hits,
                    overlap_bp = overlap_widths
                )

                max_overlap_df <- hitsdf %>%
                    group_by(query) %>%
                    slice_max(overlap_bp, with_ties = TRUE) %>%
                    ungroup()

                subjectbesthit_grs <- subjectgrs[max_overlap_df$subject]

                n_overlap_table <- subjectbesthit_grs$name %>% table()

                n_overlap_df <- tibble(subject_region = names(n_overlap_table), n_overlap = as.integer(n_overlap_table))
                if (nrow(n_overlap_df) == 0) {
                    n_overlap_df <- tibble(
                        subject_region = mcols(subjectgrs)$name %>% unique(),
                        n_overlap = 0
                    ) %>%
                        mutate(n_total = querygrs %>%
                            as.data.frame() %>%
                            nrow())
                } else {
                    n_overlap_df <- n_overlap_df %>%
                        mutate(subject_region = factor(subject_region, levels = mcols(subjectgrs)$name %>% unique())) %>%
                        complete(subject_region, fill = list(n_overlap = 0)) %>%
                        mutate(n_total = querygrs %>%
                            as.data.frame() %>%
                            nrow())
                }

                subject_split <- split(subjectgrs, subjectgrs$name)
                subjectsize_list <- map(subject_split, ~ sum(width(.x)))
                subjectsize_df <- tibble(subject_region = names(subjectsize_list), size = unlist(subjectsize_list))



                resdf <- n_overlap_df %>%
                    left_join(subjectsize_df) %>%
                    mutate(genome_size = genome_size) %>%
                    mutate(prop_genome = size / genome_size) %>%
                    mutate(fold_enrichment = (n_overlap / n_total) / prop_genome) %>%
                    rowwise() %>%
                    mutate(p_value = binom.test(n_overlap, n_total, p = prop_genome)$p.value) %>%
                    ungroup() %>%
                    dplyr::rename(name = subject_region)
                return(resdf)
            }
            # load annotations of interest
            ccresdf <- read_delim(conf$ccres, col_names = FALSE)
            ccresgr <- GRanges(
                seqnames = ccresdf$X1,
                ranges = IRanges(start = ccresdf$X2, end = ccresdf$X3),
                name = ccresdf$X10,
                UID = ccresdf$X4,
            )

            chromHMMgr <- import(conf$chromHMM)
            width(chromHMMgr) %>% sum()

            prom_no_mcols <- promoters
            mcols(prom_no_mcols) <- NULL
            mcols(prom_no_mcols)$name <- "prom"
            chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]

            enh_no_mcols <- chromHMM_enhancers_grs
            mcols(enh_no_mcols) <- NULL
            mcols(enh_no_mcols)$name <- "enh"
            prom_enh <- c(enh_no_mcols, prom_no_mcols)


            fa <- Rsamtools::FaFile(confALL$aref$ref)
            genome_size_chrfiltered <- seqinfo(fa) %>%
                data.frame() %>%
                rownames_to_column("seqnames") %>%
                filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %$%
                seqlengths %>%
                sum()

            cpgi_features_temp <- cpgi_features
            mcols(cpgi_features_temp)$name <- mcols(cpgi_features_temp)$name %>% gsub("CpG.*", "Island", .)

            genome_ranges <- seqinfo(fa) %>%
                data.frame() %>%
                rownames_to_column("seqnames") %>%
                filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
                tibble() %>%
                mutate(start = 1, end = seqlengths) %>%
                dplyr::select(seqnames, start, end) %>%
                GRanges()

            open_sea_regions <- GenomicRanges::setdiff(genome_ranges, reduce(cpgi_features))
            mcols(open_sea_regions)$name <- "OpenSea"
            cpgi_features_forenrichment <- c(cpgi_features_temp, open_sea_regions)

            flrtepromoter_for_enrichment <- flRTEpromoter %>%
                mutate(name = rte_subfamily) %>%
                filter(name != "Other") %>%
                GRanges()

            annotations_of_interest <- list(cpgi = cpgi_features_forenrichment, flrteprom = flrtepromoter_for_enrichment, chromHMM = chromHMMgr, cCREs = ccresgr, prom_enh = prom_enh)
            # note I should be testing for enrichment of these sets too
            levelslist <- list()
            levelslist[["chromHMM"]] <- c(
                "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "TssBiv",
                "Tx", "TxWk",
                "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "EnhBiv", "Het", "ZNF/Rpts", "ReprPC", "ReprPCWk", "Quies"
            )
            levelslist[["rtesubfamily"]] <- c(
                "AluY", "HERVK_LTR", "HERVL_LTR",
                "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"
            )
            levelslist[["ccres"]] <- c("PLS,CTCF-bound", "PLS", "pELS,CTCF-bound", "pELS", "dELS,CTCF-bound", "dELS", "DNase-H3K4me3,CTCF-bound", "CTCF-only,CTCF-bound", "DNase-H3K4me3", "DNase-only")

            for (annotation_of_interest in names(annotations_of_interest)) {
                annotdf <- as.data.frame(annotations_of_interest[[annotation_of_interest]]) %>%
                    tibble() %>%
                    filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
                annot <- GRanges(annotdf)
                annot_split <- split(annot, annot$name)

                enrichdf <- map(
                    dmrsgr_split,
                    ~ map(annot_split, function(region) get_region_enrichment(.x, region, genome_size = genome_size_chrfiltered)) %>%
                        as.data.frame() %>%
                        as.matrix() %>%
                        t() %>%
                        as.data.frame() %>%
                        rownames_to_column("name")
                ) %>%
                    bind_rows(.id = "type") %>%
                    tibble() %>%
                    mutate(direction_threshold = sub("t", "", type)) %>%
                    {
                        if (annotation_of_interest == "chromHMM") {
                            mutate(., name = gsub("ZNF.Rpts", "ZNF/Rpts", name)) %>%
                                mutate(name = factor(name, levels = levelslist[["chromHMM"]])) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["chromHMM"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "flrteprom") {
                            mutate(., name = factor(name, levels = rev(levelslist[["rtesubfamily"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["rtesubfamily"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "cCREs") {
                            mutate(., name = factor(name, levels = rev(levelslist[["ccres"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["ccres"]], width = 20, side = "left", pad = "_")))
                        } else {
                            mutate(., name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_")))
                        }
                    } %>%
                    mutate(enrichment_direction = "dmrs_enriched_for_a_region")

                enrichdfrev <- map(
                    annot_split,
                    ~ map(dmrsgr_split, function(region) get_region_enrichment(.x, region, genome_size = genome_size_chrfiltered)) %>%
                        as.data.frame() %>%
                        as.matrix() %>%
                        t() %>%
                        as.data.frame() %>%
                        rownames_to_column("type")
                ) %>%
                    bind_rows(.id = "name") %>%
                    tibble() %>%
                    mutate(direction_threshold = sub("t", "", type)) %>%
                    {
                        if (annotation_of_interest == "chromHMM") {
                            mutate(., name = gsub("ZNF.Rpts", "ZNF/Rpts", name)) %>%
                                mutate(name = factor(name, levels = levelslist[["chromHMM"]])) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["chromHMM"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "flrteprom") {
                            mutate(., name = factor(name, levels = rev(levelslist[["rtesubfamily"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["rtesubfamily"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "cCREs") {
                            mutate(., name = factor(name, levels = rev(levelslist[["ccres"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["ccres"]], width = 20, side = "left", pad = "_")))
                        } else {
                            mutate(., name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_")))
                        }
                    } %>%
                    mutate(enrichment_direction = "region_enriched_for_a_dmr")

                bidirec_enrich_df <- bind_rows(enrichdf, enrichdfrev)

                enrichdf_bestoverlap <- map(
                    dmrsgr_split,
                    ~ get_region_enrichment2(.x, annot, genome_size = genome_size_chrfiltered)
                ) %>%
                    bind_rows(.id = "type") %>%
                    tibble() %>%
                    mutate(direction_threshold = sub("t", "", type)) %>%
                    {
                        if (annotation_of_interest == "chromHMM") {
                            mutate(., name = gsub("ZNF.Rpts", "ZNF/Rpts", name)) %>%
                                mutate(name = factor(name, levels = levelslist[["chromHMM"]])) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["chromHMM"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "flrteprom") {
                            mutate(., name = factor(name, levels = rev(levelslist[["rtesubfamily"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["rtesubfamily"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "cCREs") {
                            mutate(., name = factor(name, levels = rev(levelslist[["ccres"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["ccres"]], width = 20, side = "left", pad = "_")))
                        } else {
                            mutate(., name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_")))
                        }
                    }


                max_abs <- log(max(abs(bidirec_enrich_df %>% filter(fold_enrichment != 0) %$% fold_enrichment))) * 1.1
                p <- bidirec_enrich_df %>%
                    mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                    ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                    geom_col(position = "dodge", color = "black") +
                    geom_vline(xintercept = 0, color = "green") +
                    coord_cartesian(xlim = c(-max_abs, max_abs)) +
                    facet_grid(rows = vars(direction), cols = vars(enrichment_direction)) +
                    scale_methylation_thresholds +
                    ggtitle("Enrichments") +
                    mtclosedgridv
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet.pdf", params$mod_code, annotation_of_interest, "bidirectional_all", annotation_of_interest), 8, 6)

                p <- bidirec_enrich_df %>%
                    filter(grepl("05", direction_threshold)) %>%
                    mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                    ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                    geom_col(position = "dodge", color = "black") +
                    geom_vline(xintercept = 0, color = "green") +
                    coord_cartesian(xlim = c(-max_abs, max_abs)) +
                    facet_grid(rows = vars(direction), cols = vars(enrichment_direction)) +
                    scale_methylation_thresholds +
                    ggtitle("Enrichments") +
                    mtclosedgridv
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet_05.pdf", params$mod_code, annotation_of_interest, "bidirectional_all", annotation_of_interest), 8, 6)

                p <- bidirec_enrich_df %>%
                    filter(grepl("05", direction_threshold)) %>%
                    mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                    ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                    geom_col(position = "dodge", color = "black") +
                    geom_vline(xintercept = 0, color = "green") +
                    coord_cartesian(xlim = c(-max_abs, max_abs)) +
                    facet_grid(cols = vars(enrichment_direction)) +
                    scale_methylation_thresholds +
                    ggtitle("Enrichments") +
                    mtclosedgridv
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet_05_together.pdf", params$mod_code, annotation_of_interest, "bidirectional_all", annotation_of_interest), 7, 4)


                enrichdflist <- list("alloverlaps" = enrichdf, "alloverlapsrev" = enrichdfrev, "bestoverlaps" = enrichdf_bestoverlap)
                for (enrichtype in names(enrichdflist)) {
                    enrichdftemp <- enrichdflist[[enrichtype]]

                    p <- enrichdftemp %>%
                        ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                        geom_col(position = "dodge", color = "black") +
                        scale_methylation_thresholds +
                        mtclosedgridv
                    mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s.pdf", params$mod_code, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

                    p <- enrichdftemp %>%
                        filter(direction_threshold %in% c("Hyper_05", "Hypo_05")) %>%
                        ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                        geom_col(position = "dodge", color = "black") +
                        scale_methylation_thresholds +
                        mtclosedgridv
                    mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_05.pdf", params$mod_code, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

                    max_abs <- log(max(abs(enrichdftemp %>% filter(fold_enrichment != 0) %$% fold_enrichment))) * 1.1
                    p <- enrichdftemp %>%
                        mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                        ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                        geom_col(position = "dodge", color = "black") +
                        geom_vline(xintercept = 0, color = "green") +
                        coord_cartesian(xlim = c(-max_abs, max_abs)) +
                        facet_wrap(~direction) +
                        scale_methylation_thresholds +
                        ggtitle("Enrichments") +
                        mtclosedgridv
                    mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet.pdf", params$mod_code, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

                    p <- enrichdftemp %>%
                        mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                        ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                        geom_col(position = "dodge", color = "black") +
                        geom_vline(xintercept = 0, color = "green") +
                        facet_wrap(~direction) +
                        scale_methylation_thresholds +
                        ggtitle("Enrichments") +
                        mtclosedgridv
                    mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet_noaxislimmirror.pdf", params$mod_code, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

                    p <- enrichdftemp %>%
                        mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                        ggplot(aes(y = name_padded, x = log(fold_enrichment), color = direction_threshold)) +
                        geom_point() +
                        geom_vline(xintercept = 0, linetype = "dashed") +
                        coord_cartesian(xlim = c(-max_abs, max_abs)) +
                        facet_wrap(~direction) +
                        scale_methylation_thresholds +
                        ggtitle("Enrichments") +
                        mtclosedgridv
                    mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmr_enrichments_in_%s_facet_dot.pdf", params$mod_code, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)
                }
                # counts

                mbo <- mergeByOverlaps(annot, dmrsgr)
                mbodf <- tibble(as.data.frame(mbo)) %>%
                    {
                        if (annotation_of_interest == "chromHMM") {
                            mutate(., name = gsub("ZNF.Rpts", "ZNF/Rpts", name)) %>%
                                mutate(name = factor(name, levels = levelslist[["chromHMM"]])) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["chromHMM"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "flrteprom") {
                            mutate(., name = factor(name, levels = rev(levelslist[["rtesubfamily"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["rtesubfamily"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "cCREs") {
                            mutate(., name = factor(name, levels = rev(levelslist[["ccres"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["ccres"]], width = 20, side = "left", pad = "_")))
                        } else {
                            mutate(., name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_")))
                        }
                    }

                p <- mbodf %>%
                    group_by(name_padded, dmr_type, direction) %>%
                    summarise(n = n()) %>%
                    ggplot() +
                    geom_col(aes(x = name_padded, y = n, fill = direction), position = "dodge", color = "black") +
                    labs(x = "") +
                    coord_flip() +
                    ggtitle(str_glue("{annotation_of_interest} Methylation")) +
                    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                    mtopen +
                    scale_methylation
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/dmrs_in_%s.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                library(scales)
                p <- mbodf %>%
                    filter(dmr_type != "t05CG10") %>%
                    filter(dmr_type != "t001") %>%
                    group_by(name_padded, direction, dmr_type) %>%
                    summarize(n = n()) %>%
                    ungroup() %>% # Summarize with groups dropped for completeness
                    tidyr::complete(name_padded, direction, dmr_type, fill = list(n = 0)) %>%
                    mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
                    ggplot() +
                    geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = name_padded, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                    geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = name_padded, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                    labs(x = "") +
                    coord_flip() +
                    ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                    annotation_logticks(sides = "b") +
                    scale_y_log10(
                        breaks = scales::trans_breaks("log10", function(x) 10^x),
                        labels = scales::trans_format("log10", math_format(10^.x))
                    ) +
                    mtopen +
                    scale_methylation_thresholds
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/dmrs_in_%s_log.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                p <- mbodf %>%
                    filter(dmr_type != "t05CG10") %>%
                    filter(dmr_type != "t001") %>%
                    group_by(name_padded, direction, dmr_type) %>%
                    summarize(n = n()) %>%
                    ungroup() %>% # Summarize with groups dropped for completeness
                    tidyr::complete(name_padded, direction, dmr_type, fill = list(n = 0)) %>%
                    mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
                    ggplot() +
                    geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = name_padded, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                    geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = name_padded, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                    labs(x = "") +
                    geom_text(
                        data = . %>% filter(dmr_type %in% c("t01", "t05")),
                        aes(x = name_padded, y = n * 1.1, label = n, group = direction_threshold), # use direction_threshold for dodging
                        position = position_dodge2(width = 0.9, preserve = "single"),
                        vjust = 0, # 0 = top-aligned; more predictable than vjust=1 on log scale
                        size = 3
                    ) +
                    coord_flip() +
                    ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                    annotation_logticks(sides = "b") +
                    scale_y_log10(
                        breaks = scales::trans_breaks("log10", function(x) 10^x),
                        labels = scales::trans_format("log10", math_format(10^.x))
                    ) +
                    mtopen +
                    scale_methylation_thresholds
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/dmrs_in_%s_log_numannot.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                total <- annotdf %>%
                    group_by(name) %>%
                    summarize(n = n())
                totaldm <- mbodf %>%
                    group_by(name, name_padded, direction, dmr_type) %>%
                    summarize(n = n()) %>%
                    ungroup() %>% # Summarize with groups dropped for completeness
                    tidyr::complete(name_padded, direction, dmr_type, fill = list(n = 0)) %>%
                    mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type)))

                pctdm <- left_join(totaldm, total, by = c("name")) %>%
                    mutate(pct = 100 * n.x / n.y) %>%
                    mutate(myaxis = paste0(name_padded, "\n", "n=", n.y)) %>%
                    drop_na() %>%
                    {
                        if (annotation_of_interest == "chromHMM") {
                            mutate(., name = gsub("ZNF.Rpts", "ZNF/Rpts", name)) %>%
                                mutate(name = factor(name, levels = levelslist[["chromHMM"]])) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["chromHMM"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "flrteprom") {
                            mutate(., name = factor(name, levels = rev(levelslist[["rtesubfamily"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["rtesubfamily"]], width = 20, side = "left", pad = "_")))
                        } else if (annotation_of_interest == "cCREs") {
                            mutate(., name = factor(name, levels = rev(levelslist[["ccres"]]))) %>%
                                mutate(name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_"), levels = str_pad(levelslist[["ccres"]], width = 20, side = "left", pad = "_")))
                        } else {
                            mutate(., name_padded = factor(str_pad(name, width = 20, side = "left", pad = "_")))
                        }
                    }

                p <- pctdm %>% ggplot() +
                    geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = myaxis, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                    geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = myaxis, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                    labs(x = "", y = "Pct Differentially Methylated") +
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                    ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                    coord_flip() +
                    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                    mtopen +
                    scale_methylation_thresholds
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/dmrs_in_%s_pct.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                p <- pctdm %>% ggplot() +
                    geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = name_padded, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                    geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = name_padded, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                    labs(x = "", y = "Pct Differentially Methylated") +
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                    ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                    coord_flip() +
                    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                    mtopen +
                    scale_methylation_thresholds
                mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/dmrs_in_%s_pct_orderedaxis.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)
            }
        }
    }



    # dmrs heatmap
    dmrsuuid <- dmrs %>%
        mutate(dmrid = paste0(direction, row_number()))
    dmrs05 <- dmrsuuid %>%
        filter(dmr_type == "t05") %>%
        GRanges()
    dmrs01 <- dmrsuuid %>%
        filter(dmr_type == "t01") %>%
        GRanges()
    dmrs05setdiff <- dmrs05 %>% subsetByOverlaps(dmrs01, invert = TRUE)

    dmrsgrsuuid <- c(dmrs01, dmrs05setdiff)
    mbo <- mergeByOverlaps(grs, dmrsgrsuuid)
    dmrs_meth_df <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    dmrs_meth_df$dmrid <- mbo$dmrsgrsuuid$dmrid
    dmrs_meth_df$areaStat <- mbo$dmrsgrsuuid$areaStat
    dmrs_meth_df$direction <- mbo$dmrsgrsuuid$direction
    dmrs_meth_df$dmr_type <- mbo$dmrsgrsuuid$dmr_type
    dmrs_meth_df$length <- mbo$dmrsgrsuuid$length
    dmrs_meth_df$nCG <- mbo$dmrsgrsuuid$nCG

    # dmr av meth stats
    ac <- dmrs_meth_df %>%
        group_by(pos, direction, dmr_type, condition) %>%
        summarise(mean_meth = mean(pctM))
    ad <- ac %>%
        filter(dmr_type %in% c("t01", "t05")) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        mutate(dif = !!sym(condition2) - !!sym(condition1)) %>%
        group_by(dmr_type, direction) %>%
        summarize(mean_dif = mean(dif))
    dir.create(sprintf("ldna/results/%s/tables/dmrs", params$mod_code), recursive = TRUE)
    ad %>% write_delim(sprintf("ldna/results/%s/tables/dmrs/dmrs_meth_dif", params$mod_code))

    # topdmrshypo <- dmrsgrsuuid %>%
    #     as.data.frame() %>%
    #     tibble() %>%
    #     filter(direction == "Hypo") %>%
    #     arrange(areaStat) %$% dmrid %>%
    #     head(n = 1000)
    # topdmrshyper <- dmrsgrsuuid %>%
    #     as.data.frame() %>%
    #     tibble() %>%
    #     filter(direction == "Hyper") %>%
    #     arrange(-areaStat) %$% dmrid %>%
    #     head(n = 1000)

    pf1 <- dmrs_meth_df %>%
        group_by(sample, dmrid, direction, dmr_type, areaStat) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup()

    topdmrshypo <- pf1 %>%
        as.data.frame() %>%
        tibble() %>%
        filter(direction == "Hypo") %>%
        arrange(areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = sample_table$sample_name, names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    topdmrshyper <- pf1 %>%
        as.data.frame() %>%
        tibble() %>%
        filter(direction == "Hyper") %>%
        arrange(-areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = sample_table$sample_name, names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    pf <- bind_rows(topdmrshyper, topdmrshypo)

    library(tidyHeatmap)
    p <- pf %>%
        mutate(sample = factor(sample, levels = conf$samples)) %>%
        group_by(direction) %>%
        heatmap(dmrid, sample, mean_meth,
            cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_row_dend = FALSE # palette_value = circlize::colorRamp2(
            # seq(0, 100, length.out = 11),
            # RColorBrewer::brewer.pal(11, "RdBu")
            # )
        ) %>%
        add_tile(sex)

    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap.pdf", params$mod_code), w = 5, h = 10, res = 300, pl = p, raster = FALSE)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap.pdf", params$mod_code), w = 5, h = 10, res = 300, pl = p, raster = TRUE)

    p <- pf %>%
        mutate(sample = factor(sample, levels = conf$samples)) %>%
        group_by(direction) %>%
        heatmap(dmrid, sample, mean_meth,
            cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_row_dend = FALSE, column_km = 2 # palette_value = circlize::colorRamp2(
            # seq(0, 100, length.out = 11),
            # RColorBrewer::brewer.pal(11, "RdBu")
            # )
        ) %>%
        add_tile(sex)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap_colclust.pdf", params$mod_code), w = 4.5, h = 9, res = 300, pl = p, raster = FALSE)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap_colclust.pdf", params$mod_code), w = 4.5, h = 9, res = 300, pl = p, raster = TRUE)


    pf <- dmrs_meth_df %>%
        group_by(sample, dmrid, condition, direction, dmr_type, nCG, length) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup()

    p <- pf %>%
        filter(condition == condition1) %>%
        group_by(direction, dmrid, dmr_type) %>%
        summarise(mean_meth = mean(mean_meth)) %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        mutate(direction_threshold = factor(direction_threshold, levels = c("Hypo_05", "Hypo_01", "Hyper_05", "Hyper_01"))) %>%
        ggplot(aes(x = mean_meth, fill = direction_threshold)) +
        xlab(sprintf("Average DMR Meth (%s)", condition1)) +
        geom_histogram(color = "black") +
        facet_wrap(~direction, scales = "free_y") +
        scale_methylation_thresholds +
        mtclosed
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_histogram.pdf", params$mod_code), w = 6, h = 4, res = 300, pl = p)


    p <- dmrsgrsuuid %>%
        as.data.frame() %>%
        tibble() %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        mutate(direction_threshold = factor(direction_threshold, levels = c("Hypo_05", "Hypo_01", "Hyper_05", "Hyper_01"))) %>%
        ggplot(aes(x = nCG, fill = direction_threshold)) +
        xlab(sprintf("Average DMR Meth (%s)", condition1)) +
        geom_histogram(color = "black", bins = 10) +
        facet_wrap(~direction, scales = "free_y") +
        scale_methylation_thresholds +
        lims(x = c(0, 15)) +
        mtclosed
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_histogram_nCG.pdf", params$mod_code), w = 6.5, h = 4, res = 300, pl = p)
    p <- dmrsgrsuuid %>%
        as.data.frame() %>%
        tibble() %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        mutate(direction_threshold = factor(direction_threshold, levels = c("Hypo_05", "Hypo_01", "Hyper_05", "Hyper_01"))) %>%
        ggplot(aes(x = length, fill = direction_threshold)) +
        xlab("DMR Length (bp)") +
        geom_histogram(color = "black", bins = 20) +
        facet_wrap(~direction, scales = "free_y") +
        scale_methylation_thresholds +
        lims(x = c(0, 1500)) +
        mtclosed
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_histogram_length.pdf", params$mod_code), w = 6, h = 4, res = 300, pl = p)

    dmrsgrsuuid %>%
        as.data.frame() %>%
        tibble() %$% nCG %>%
        mean()

    dmrsgrsuuid %>%
        as.data.frame() %>%
        tibble() %$% length %>%
        mean()


    pcaframe <- dmrs_meth_df %>%
        group_by(sample, dmrid, condition, direction) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup() %>%
        dplyr::select(sample, mean_meth, dmrid) %>%
        pivot_wider(names_from = dmrid, values_from = mean_meth) %>%
        mutate(sample = factor(sample, levels = conf$samples)) %>%
        arrange(sample) %>%
        column_to_rownames(var = "sample") %>%
        as.matrix() %>%
        t()
    library(PCAtools)

    pcaObj <- pca(pcaframe, center = TRUE, scale = FALSE, metadata = sample_table %>% column_to_rownames(var = "sample_name"))

    p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_scree.pdf", params$mod_code), w = 5, h = 10, res = 300, pl = p)


    p <- plotloadings(pcaObj,
        components = getComponents(pcaObj, seq_len(3)),
        rangeRetain = 0.045, labSize = 2
    ) + mtopen
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_loadings.pdf", params$mod_code), w = 5, h = 10, res = 300, pl = p)

    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = "sex", colby = "condition",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmr_biplot.pdf", params$mod_code), w = 5, h = 5, res = 300, pl = p)

    # p <- dmrs %>%
    #     ggplot(aes(x = meanMethy_c1, y = meanMethy_c2)) +
    #     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    #     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    #     scale_x_continuous(expand = c(0, 0)) +
    #     scale_y_continuous(expand = c(0, 0)) +
    #     scale_fill_distiller(palette = "Spectral", direction = 1) +
    #     xlab(sprintf("CpG Methylation %s", condition1)) +
    #     ylab(sprintf("CpG Methylation %s", condition2)) +
    #     ggtitle("DMR Density") +
    #     mtclosed
    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrdensity.pdf", params$mod_code), 5, 5)
}



## CENTROMERE
# ANNOTATIONS
##########################################

giesma <- read_delim(conf$cytobands, col_names = FALSE)
cent_gr <- giesma %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
cenRegion <- reduce(cent_gr)
lengths(cenRegion) %>% sum()

cenRegionMeth <- subsetByOverlaps(grs, cenRegion)
cenRegionMethdf <- tibble(as.data.frame(cenRegionMeth))
cenRegionMethdf %$% seqnames %>% unique()
colnames(grsdf)
colnames(mcols(grs))
# group Av
p <- cenRegionMethdf %>%
    filter(cov > 4) %>%
    ggplot() +
    geom_boxplot(aes(x = seqnames, y = pctM, fill = condition), outlier.shape = NA) +
    theme_cowplot() +
    scale_conditions +
    theme(aspect.ratio = 0.33) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") +
    ylab("CpG Fraction Methylated") +
    ggtitle("", ) +
    theme(plot.title = element_text(hjust = 0.5))
mysaveandstore(sprintf("ldna/results/%s/plots/centromere/cenRegion_boxplot.pdf", params$mod_code), 12, 5)

# group rolling av
pf <- cenRegionMethdf %>%
    filter(cov > 4) %>%
    group_by(seqnames, condition) %>%
    mutate(rM = rollmean(pctM, 100, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()
p <- pf %>% ggplot() +
    geom_line(aes(x = start, y = rM, color = condition), alpha = 0.5) +
    scale_conditions +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~seqnames, scales = "free_x")
mysaveandstore(sprintf("ldna/results/%s/plots/centromere/cenRegion_Lines.pdf", params$mod_code), 12, 12)

###########
censat <- import.bed("resources/genomes/hs1/annotations/censat.bed")
censatdf <- censat %>%
    as.data.frame() %>%
    tibble() %>%
    mutate(types = gsub("_.*", "", name))

censat <- GRanges(censatdf)
censatMeth <- subsetByOverlaps(grs, censat)
censatMethOL <- findOverlaps(grs, censat)
typesOL <- censatdf$types[censatMethOL@to]
censatMeth$types <- typesOL
censatMethdf <- tibble(as.data.frame(censatMeth))
censatMethdf %>%
    group_by(types, condition) %>%
    filter(cov > 4) %>%
    summarise(avMeth = mean(pctM))


sample_size <- censatMethdf %>%
    group_by(types) %>%
    summarize(num = n())

# Plot
p <- censatMethdf %>%
    left_join(sample_size) %>%
    mutate(myaxis = paste0(types, "\n", "n=", num)) %>%
    filter(cov > 4) %>%
    ggplot() +
    geom_boxplot(aes(x = myaxis, y = pctM, fill = condition), outlier.shape = NA) +
    theme_cowplot() +
    scale_conditions +
    theme(aspect.ratio = 0.33) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") +
    ylab("CpG Fraction Methylated") +
    ggtitle("", ) +
    theme(plot.title = element_text(hjust = 0.5))
mysaveandstore(sprintf("ldna/results/%s/plots/centromere/censatTypes_boxplot.pdf", params$mod_code), 12, 6)



#######

# biotype_levels <- c("repeat", "Other", "lncRNA", "miRNA", "snRNA", "tRNA", "rRNA", "protein_coding")


# mutate(Class = factor(Class, levels = c("Gene", "LowComp", "SimpleRep", "Retroposon", "SAT", "DNA", "LTR", "LINE", "SINE", "Other"))) %>%




##############
# METHYLATION CLUSTERING

outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

# rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)

cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()
consensus_ss[[1]][5762:5763]


cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

methdf <- rtedf %>% filter(rte_subfamily == subfam)
mdf <- methdf %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

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
cpg_order <- merged %$% consensus_pos %>%
    unique() %>%
    sort()
merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))
library(tidyHeatmap)

dat <- merged %>%
    filter(!is.na(pctM))
dat %$% gene_id %>%
    unique() %>%
    length()
# write_csv(dat, "meth_for_bayes.csv")


# methylation heatmaps
hms <- list()
for (sample in conf$samples) {
    p <- merged %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE) %>%
        add_tile(genic_loc) %>%
        add_tile(intactness_req) %>%
        add_tile(loc_highres_integrative_stranded) %>%
        as_ComplexHeatmap()
    hms[[sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 6)
}
# Generate the expression as a string and parse it
p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
rm(p)

# methylation heatmaps
hms <- list()
for (sample in conf$samples) {
    p <- merged %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE) %>%
        as_ComplexHeatmap()
    hms[[sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s_noannot.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 6)
}
# Generate the expression as a string and parse it
p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
mysaveandstore(sprintf("%s/%s_methylation_%s_noannot.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
rm(p)


# nonref analysis

grsdf_nr <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf_nonref.tsv", params$mod_code))
grs_nr <- GRanges(grsdf_nr)
rmann_nr_list <- list()
merged_nr_list <- list()
for (sample in sample_table$sample_name) {
    rmann_nr_temp <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
    rmann_nr_temp$sample_name <- sample
    rmann_nr_list[[sample]] <- rmann_nr_temp
    grs_nr_temp <- grs_nr[mcols(grs_nr)$sample == sample]
    merged_temp <- merge_with_grs(grs_nr_temp, GRanges(rmann_nr_temp))
    merged_nr_list[[sample]] <- merged_temp
}

rmann_nr <- do.call(rbind, rmann_nr_list) %>%
    tibble() %>%
    mutate(gene_id = paste0(sample_name, "___", gene_id)) %>%
    mutate(seqnames = paste0(sample_name, "___", seqnames))

merged_nr <- do.call(rbind, merged_nr_list) %>%
    tibble() %>%
    mutate(gene_id = paste0(sample_name, "___", gene_id)) %>%
    mutate(seqnames = paste0(sample_name, "___", seqnames))

l1hs_nr <- GRanges(merged_nr %>% filter(rte_subfamily == "L1HS")) %>%
    as.data.frame() %>%
    tibble()


outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

# rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)

cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()
cg_indices %>% length()
cg_indices[(909 >= cg_indices)] %>% length()
cg_indices[(500 > cg_indices)] %>% length()
cg_indices[(328 > cg_indices)] %>% length()

p <- tibble(cg_site = cg_indices) %>%
    ggplot(aes(x = cg_site)) +
    geom_rect(aes(xmin = 0, xmax = 328, ymin = 0, ymax = 0.25), fill = "yellow") +
    geom_rect(aes(xmin = 0, xmax = 500, ymin = 0.25, ymax = 0.5), fill = "orange") +
    geom_rect(aes(xmin = 0, xmax = 909, ymin = 0.5, ymax = 0.75), fill = "red") +
    geom_rect(aes(xmin = 0, xmax = 6031, ymin = 0.75, ymax = 1), fill = "blue") +
    geom_segment(aes(y = 0, yend = 1)) +
    scale_x_continuous(breaks = seq(0, 6000, by = 500)) +
    mtclosed +
    theme(
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.y = element_blank(), # Remove y-axis text
        axis.ticks.y = element_blank() # Remove y-axis ticks
    )
mysaveandstore(sprintf("%s/l1_consensus_cpg.pdf", outputdir_meth_clustering), w = 12, h = 2)

p <- tibble(cg_site = cg_indices[(909 >= cg_indices)]) %>%
    ggplot(aes(x = cg_site)) +
    geom_rect(aes(xmin = 0, xmax = 328, ymin = 0, ymax = 0.25), fill = "yellow") +
    geom_rect(aes(xmin = 0, xmax = 500, ymin = 0.25, ymax = 0.5), fill = "orange") +
    geom_rect(aes(xmin = 0, xmax = 909, ymin = 0.5, ymax = 0.75), fill = "red") +
    geom_rect(aes(xmin = 0, xmax = 909, ymin = 0.75, ymax = 1), fill = "blue") +
    geom_segment(aes(y = 0, yend = 1)) +
    scale_x_continuous(breaks = seq(0, 900, by = 150)) +
    mtclosed +
    theme(
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.y = element_blank(), # Remove y-axis text
        axis.ticks.y = element_blank() # Remove y-axis ticks
    )
mysaveandstore(sprintf("%s/l1_5utrconsensus_cpg.pdf", outputdir_meth_clustering), w = 12, h = 2)


consensus_ss[[1]][5762:5763]
cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

methdf <- l1hs_nr %>% filter(rte_subfamily == subfam)
mdf <- methdf %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

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

merged <- left_join(cg_positions_df, mdf, by = c("gene_id", "sequence_pos")) %>% filter()
cpg_order <- merged %$% consensus_pos %>%
    unique() %>%
    sort()
merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))
library(tidyHeatmap)

dat <- merged %>%
    filter(!is.na(pctM))
dat %$% gene_id %>%
    unique() %>%
    length()
# write_csv(dat, "meth_for_bayes.csv")


# methylation heatmaps
hms <- list()
for (sample in conf$samples) {
    p <- merged %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        distinct() %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE) %>%
        as_ComplexHeatmap()
    hms[[sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/nonref_%s_methylation_%s.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 3)
}
merged %$% gene_id
p <- merged %>%
    filter(!is.na(sample)) %>%
    group_by(gene_id) %>%
    mutate(cpgs_detected_per_element = n()) %>%
    ungroup() %>%
    filter(cpgs_detected_per_element > 50) %>%
    distinct() %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    group_by(sample_name) %>%
    heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE) %>%
    as_ComplexHeatmap()
dir.create(outputdir_meth_clustering, recursive = TRUE)
mysaveandstore(sprintf("%s/nonref_%s_methylation_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 6, h = 8)


filter_by_consensus_pos <- function(fl_grs, pos_mapping, include_up_to_pos) {
    pos_genes <- pos_mapping %$% gene_id %>% unique()
    grs_genes <- mcols(fl_grs)$gene_id %>% unique()
    genes_to_map <- intersect(pos_genes, grs_genes)
    filter_pos_list <- list()
    i <- 0
    for (element in genes_to_map) {
        i <- i + 1
        print(i)
        print(element)
        dfs <- pos_mapping %>% filter(gene_id == element)
        seqval <- dfs %>%
            filter(consensus_pos == include_up_to_pos) %$% sequence_pos %>%
            pluck(1)
        if (!is.na(seqval)) {
            filter_pos <- seqval
        } else {
            start_pos <- include_up_to_pos - 1
            match <- FALSE
            while (match == FALSE) {
                seqval <- dfs %>% filter(consensus_pos == start_pos) %$% sequence_pos
                if (length(seqval) != 0) {
                    if (!is.na(seqval)) {
                        filter_pos <- seqval
                        match <- TRUE
                    } else {
                        start_pos <- start_pos - 1
                    }
                } else {
                    filter_pos <- "NoRegionHomology"
                    match <- TRUE
                }
            }
        }
        filter_pos_list[[element]] <- filter_pos
    }
    print("loop done")
    mapping <- tibble(gene_id = names(filter_pos_list), filter_pos = unlist(filter_pos_list))
    gene_ids_with_homology <- mapping %>% filter(filter_pos != "NoRegionHomology") %$% gene_id
    fl_grs_with_homology <- fl_grs[mcols(fl_grs)$gene_id %in% gene_ids_with_homology]
    grlist <- map(seq_along(fl_grs_with_homology), function(x) fl_grs_with_homology[x])
    grlistresized <- map(grlist, function(x) {
        if (mcols(x)$gene_id %in% mapping$gene_id) {
            return(resize(x, width = (mapping %>% filter(gene_id == mcols(x)$gene_id) %$% filter_pos)))
        }
    })
    number_omitted <- length(grlist) - nrow(mapping)
    l1hs_resized <- purrr::reduce(grlistresized, c)
    return(l1hs_resized)
}

flL1HS5UTR <- filter_by_consensus_pos(fl_grs = rmann_nr %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %>% GRanges(), pos_mapping = consensus_index_long, include_up_to_pos = 909)

nr_5utr <- GRanges(merged_nr %>% filter(rte_subfamily == "L1HS")) %>% subsetByOverlaps(flL1HS5UTR)

pf <- nr_5utr %>%
    as.data.frame() %>%
    tibble() %>%
    group_by(gene_id, sample, condition) %>%
    summarise(mean_meth = mean(pctM), nCpG = n(), mean_cov = mean(cov)) %>%
    ungroup()
p <- pf %>%
    filter(nCpG > 25, mean_cov > 6) %>%
    mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
    ggplot(aes(x = sample, y = mean_meth, color = sample)) +
    geom_beeswarm() +
    labs(title = "NonRef L1HS 5UTR") +
    scale_samples_unique +
    mtopen +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/nonref/l1hs_meth.pdf", params$mod_code), w = 5, h = 4, res = 300, pl = p)


# srna_seq_table <- read_csv("conf/srna_seq_chars.csv")
# p <- srna_seq_table %>%
#     mutate(across(
#         c(`# Reads`, `Yield (Mbases)`),
#         ~ formatC(as.numeric(.x), format = "e", digits = 2)
#     )) %>%
#     ggtexttable(theme = ttheme("minimal"))
# mysaveandstore(pl = p, sprintf("aref/results/sample_srna_chars.pdf"), 8, 5)


########


## READ DISPERSION ANALYSIS
{
    flyngl1 <- rmann %>%
        filter(rte_length_req == "FL") %>%
        filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2")

    flyngl1grs <- GRanges(flyngl1)
    flyngl1grspromoter <- promoters(flyngl1grs, upstream = 0, downstream = 909)


    l1cpggrs <- mergeByOverlaps(cpg_islands, flyngl1grspromoter)
    l1cpggrs1 <- l1cpggrs$cpg_islands
    mcols(l1cpggrs1)$gene_id <- mcols(l1cpggrs$flyngl1grspromoter)$gene_id

    promoterscpggrs <- mergeByOverlaps(cpg_islands, promoters)
    promoterscpggrs1 <- promoterscpggrs$cpg_islands
    mcols(promoterscpggrs1)$gene_id <- mcols(promoterscpggrs$promoters)$gene

    flanked_promoters <- flank(promoterscpggrs1, width = 10000, both = TRUE)

    merged <- mergeByOverlaps(flanked_promoters, l1cpggrs1)
    merged[!(mcols(merged$flanked_promoters)$name == mcols(merged$l1cpggrs1)$name), ]


    rmann %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        pw()


    reads_cgI_chr1 <- read_delim("ldna/intermediates/merged/merged_cpgI_reads1.tsv")
    reads_cgI_chr1_grs <- GRanges(reads_cgI_chr1 %>% dplyr::rename(seqnames = chrom, start = ref_position) %>% mutate(end = start))

    mbo <- mergeByOverlaps(reads_cgI_chr1_grs, promoters)
    mboreadsdf <- mbo$reads_cgI_chr1_grs %>%
        as.data.frame() %>%
        tibble()
    mbopromoters <- mbo$promoters %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width)

    readscg_genepromoters <- bind_cols(mboreadsdf, mbopromoters)
    rm(mbopromoters)
    rm(readsdf)
    rm(mbo)

    chr1rtes_grs <- rmann %>%
        filter(seqnames == "chr1") %>%
        filter(rte_subfamily != "Other") %>%
        GRanges()
    mbo <- mergeByOverlaps(reads_cgI_chr1_grs, chr1rtes_grs)
    mboreadsdf <- mbo$reads_cgI_chr1_grs %>%
        as.data.frame() %>%
        tibble()
    mbortes <- mbo$chr1rtes_grs %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width)

    readscg_rtes <- bind_cols(mboreadsdf, mbortes)
    rm(mboreadsdf)
    rm(mbortes)
    rm(mbo)


    by_cpg_genepromoters <- readscg_genepromoters %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    numCGneeded_genes <- by_cpg_genepromoters %>%
        group_by(gene_id, read_id) %>%
        summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        group_by(gene_id) %>%
        summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    by_cpg_genepromoters_ds <- by_cpg_genepromoters %>%
        filter(num_cpgs_in_read >= 25) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = 25, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()

    by_read_genepromoters_ds <- by_cpg_genepromoters_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id) %>%
        mutate(mean_fm = mean(fraction_meth))

    # by_read_genepromoters_ds %>%
    #     group_by(gene_id) %>%
    #     summarise(mean_fm = mean(fraction_meth)) %$% mean_fm %>%
    #     quantile(., probs = seq(0, 1, 0.05))
    # by_read_genepromoters_ds %$% gene_id %>% unique()


    # p <- by_read_genepromoters_ds %>%
    #     group_by(gene_id) %>%
    #     summarise(mean_fm = mean(fraction_meth)) %>%
    #     ggplot() +
    #     geom_histogram(aes(x = mean_fm))
    # mysaveandstore()

    # by_read_genepromoters_ds %>%
    #     filter(mean_fm > 0.7) %>%
    #     filter(mean_fm < 0.98) %$% gene_id %>%
    #     unique()


    dispersion_model_genes_high_meth <- glmmTMB(
        cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
        dispformula = ~condition,
        family = glmmTMB::betabinomial(),
        data = by_read_genepromoters_ds %>% filter(mean_fm > 0.7) %>% filter(mean_fm < 0.98) %>%
            mutate(meth = fraction_meth * 25, unmeth = 25 - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
    )

    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, "genes", "cpg25req")
    dir.create(outputdirtables, recursive = TRUE)
    summary(dispersion_model_genes_high_meth)
    res_text <- capture.output(dispersion_model_genes_high_meth %>% summary())
    writeLines(res_text, sprintf("%s/gene_high_meth_dispersion_model_summary.txt", outputdirtables))




    by_cpg_rtes <- readscg_rtes %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    numCGneeded_rtes <- by_cpg_rtes %>%
        group_by(gene_id, read_id) %>%
        summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        group_by(gene_id) %>%
        summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    by_cpg_rtes_ds <- by_cpg_rtes %>%
        filter(num_cpgs_in_read >= 25) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = 25, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()

    by_read_rtes_ds <- by_cpg_rtes_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id) %>%
        mutate(mean_fm = mean(fraction_meth)) %>%
        left_join(rmann %>% dplyr::select(gene_id, rte_subfamily, rte_length_req))




    rtemodels <- list()
    for (subfam in c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")) {
        dispersion_model_high_meth_fl <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = by_read_rtes_ds %>% filter(mean_fm > 0.7) %>% filter(mean_fm < 0.98) %>%
                filter(rte_subfamily == subfam) %>%
                filter(rte_length_req == "FL") %>%
                mutate(meth = fraction_meth * 25, unmeth = 25 - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
        )
        rtemodels[[subfam]] <- dispersion_model_high_meth_fl
    }

    res_text <- capture.output(map(dispersion_model_genes_high_meth, summary))
    writeLines(res_text, sprintf("%s/gene_high_meth_dispersion_model_summary.txt", outputdirtables))



    ################### reads to extract
    imprinted_genes_df <- read_csv("/users/mkelsey/data/Nanopore/alz/imprinted_genes.csv", comment = "#")
    imprinted_genes_df %$% Aliases %>%
        str_split(., ", ") %>%
        map(., ~ trimws(.x))

    temppromoteranalysis <- ip %>%
        group_by(gene_id) %>%
        summarise(mm = mean(mean_meth)) %>%
        filter(!grepl("-AS", gene_id))
    highmethgenes <- temppromoteranalysis %>% filter(mm > 75) %$% gene_id
    lowmethgenes <- temppromoteranalysis %>% filter(mm < 15) %$% gene_id
    imprinted_genes <- temppromoteranalysis %>% filter(gene_id %in% (imprinted_genes_df %$% gene_id)) %$% gene_id

    sex_chromosome_genes <- mcols(promoters[seqnames(promoters) == "chrX" | seqnames(promoters) == "chrY"])$gene_id

    # non interesting cpgIs
    # filter out any related to genes, enhancers, ccres, L1s
    boring_islands <- cpg_islands %>%
        subsetByOverlaps(promoters, invert = TRUE) %>%
        subsetByOverlaps(refseq_gr[mcols(refseq_gr)$type == "gene"], invert = TRUE) %>%
        subsetByOverlaps(rmann %>% filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")) %>% filter(element_start < 1000) %>% GRanges(), invert = TRUE) %>%
        subsetByOverlaps(ccresgr, invert = TRUE) %>%
        subsetByOverlaps(chromHMMgr[mcols(chromHMMgr)$name == "Quies"], invert = FALSE)

    set.seed(72)
    boring_islands_to_extract_reads_from <- boring_islands[width(boring_islands) < 1000 & !(seqnames(boring_islands) %in% c("chrX", "chrY"))] %>%
        sample(., size = 500)
    mcols(boring_islands_to_extract_reads_from)$score <- NULL
    mcols(boring_islands_to_extract_reads_from)$gene_id <- mcols(boring_islands_to_extract_reads_from)$name
    mcols(boring_islands_to_extract_reads_from)$name <- NULL

    set.seed(73)
    genes_to_extract_reads_from <- c(
        highmethgenes[!(highmethgenes %in% sex_chromosome_genes)] %>% sample(., size = 500),
        lowmethgenes[!(lowmethgenes %in% sex_chromosome_genes)] %>% sample(., size = 500),
        imprinted_genes[!(imprinted_genes %in% sex_chromosome_genes)]
    )

    genes_to_extract_reads_from_grs <- promoters[mcols(promoters)$gene_id %in% genes_to_extract_reads_from]

    set.seed(74)
    rtes_to_extract_reads_from_grs <- rmann %>%
        filter(!(seqnames %in% c("chrX", "chrY"))) %>%
        filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")) %>%
        filter(rte_length_req == "FL") %>%
        filter(element_start < 150) %>%
        GRanges() %>%
        promoters(., upstream = 0, downstream = 909) %>%
        as.data.frame() %>%
        tibble() %>%
        group_by(rte_subfamily) %>%
        slice_sample(n = 500, replace = FALSE) %>%
        GRanges()

    regions_to_extract_reads_from_with_mcols <- c(boring_islands_to_extract_reads_from, genes_to_extract_reads_from_grs, rtes_to_extract_reads_from_grs)
    strand(regions_to_extract_reads_from_with_mcols) <- "*"
    mcols(regions_to_extract_reads_from_with_mcols)$name <- mcols(regions_to_extract_reads_from_with_mcols)$gene_id
    rtracklayer::export.bed(regions_to_extract_reads_from_with_mcols, "ldna/Rintermediates/m/regions_to_extract_reads_from.bed")


    ###


    reads_rois <- read_delim("ldna/intermediates/merged/merged_rois_reads.tsv")
    reads_rois_grs <- GRanges(reads_rois %>% dplyr::rename(seqnames = chrom, start = ref_position) %>% mutate(end = start))


    mbo <- mergeByOverlaps(reads_rois_grs, regions_to_extract_reads_from_with_mcols)
    mboreadsdf <- mbo$reads_rois_grs %>%
        as.data.frame() %>%
        tibble()
    mborois <- mbo$regions_to_extract_reads_from_with_mcols %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width, gene_id)

    readscg_rois <- bind_cols(mboreadsdf, mborois)
    rm(mboreadsdf)
    rm(mborois)
    rm(mbo)


    by_cpg_readscg_rois <- readscg_rois %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    # numCGneeded_rois <- by_cpg_readscg_rois %>%
    #     group_by(gene_id, read_id) %>%
    #     summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
    #     group_by(gene_id) %>%
    #     summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    downsample_to_n <- 20
    set.seed(10)
    by_cpg_rois_ds <- by_cpg_readscg_rois %>%
        filter(num_cpgs_in_read >= downsample_to_n) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = downsample_to_n, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()
    by_cpg_rois_ds %>% write_csv("ldna/Rintermediates/m/by_cpg_rois_ds")

    by_read_rois_ds <- by_cpg_rois_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id, sample_name) %>%
        mutate(mean_fm = mean(fraction_meth)) #can use this to remove high meth genes that are actually low in one sample 

    by_read_imprinted <- by_read_rois_ds %>% filter(gene_id %in% imprinted_genes) %>% 
        mutate(roi = "imprinted_gene") %>%
        filter(mean_fm > .80)
    by_read_highmeth_gene <- by_read_rois_ds %>%
        filter(!(gene_id %in% imprinted_genes)) %>%
        filter(gene_id %in% highmethgenes) %>% 
        mutate(roi = "high_meth_gene") %>%
        filter(mean_fm > .80)
    by_read_lowmeth_gene <- by_read_rois_ds %>%
        filter(!(gene_id %in% imprinted_genes)) %>%
        filter(gene_id %in% lowmethgenes) %>% 
        mutate(roi = "low_meth_gene") %>%
        filter(mean_fm < .20)
    by_read_boringcpgi <- by_read_rois_ds %>% filter(grepl("CpG:", gene_id)) %>% 
        mutate(roi = "isolated_cpgI") %>%
        filter(mean_fm > .80)
    by_read_yngl1s <- by_read_rois_ds %>%
        filter(gene_id %in% mcols(rtes_to_extract_reads_from_grs)$gene_id) %>%
        left_join(rmann) %>% 
        mutate(roi = rte_subfamily) %>%
        filter(mean_fm > .80)

    dfs <- list(
        "imprinted" = by_read_imprinted,
        "highmethgene" = by_read_highmeth_gene,
        "lowmethgene" = by_read_lowmeth_gene,
        "boringcpgi" = by_read_boringcpgi
    )

    dfsl1s <- split(by_read_yngl1s %>% ungroup() %>% dplyr::select(-colnames(rmann)[!colnames(rmann) %in% c("gene_id", "rte_subfamily")]), by_read_yngl1s$rte_subfamily)

    dfsall <- c(dfs, dfsl1s)
    dfsallbound <- bind_rows(dfsall)


cpg_fraction_reads_highly_demeth <- dfsallbound %>% 
mutate(fraction_meth_lt50 = case_when(
    fraction_meth < 0.50 ~ 1,
    TRUE ~ 0,
)) %>%
mutate(fraction_meth_lt25 = case_when(
    fraction_meth < 0.25 ~ 1,
    TRUE ~ 0,
)) %>%
mutate(fraction_meth_lt10 = case_when(
    fraction_meth < 0.1 ~ 1,
    TRUE ~ 0,
)) %>% group_by(roi) %>%
 summarise(nreads = n(), 
    mlt50 = mean(fraction_meth_lt50), 
    mlt25 = mean(fraction_meth_lt25), 
    mtl10 = mean(fraction_meth_lt10))

    p <- dfsallbound %>% ggplot() +
        geom_density(aes(x = fraction_meth)) +
        facet_wrap(~roi) +
        mtclosed
    mysaveandstore(fn = "zzzte3.pdf", w = 8, h = 6)


    p <- dfsallbound %>% ggplot() +
        geom_histogram(
            aes(x = fraction_meth, y = after_stat(count / sum(count))),
            position = "identity"
        ) +
        facet_wrap(~roi) +
        mtclosed
    mysaveandstore(fn = "zzzte4.pdf", w = 8, h = 6)


    getdispmodel <- function(df) {
        dispersion_model_genes_high_meth <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = df %>%
                mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
        )
    }

    models_dfs <- map(dfs, getdispmodel)
    models_l1s <- map(dfsl1s, getdispmodel)

    modelsall <- c(models_dfs, models_l1s)


    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, "allrois_types", sprintf("cpg%sreq", downsample_to_n))
    dir.create(outputdirtables, recursive = TRUE)
    res_text <- capture.output(map(models_l1s, summary))
    writeLines(res_text, sprintf("%s/model_summary.txt", outputdirtables))

    modelsall[[1]] %>% broom::tidy()
    fixef(modelsall[[1]])$disp
    modelsall[[1]]$fit$par["thetaf"] # Overdispersion theta (log scale)
    sigma(modelsall[[1]], type = "dispersion")
    VarCorr(modelsall[[1]])
    summary(modelsall[[1]])
    summary(modelsall[[1]]) %>% as.data.frame()
    summary(modelsall[[1]])$coefficients$disp

    disp_df <- map(modelsall, ~ summary(.x)$coefficients$disp %>%
        as.data.frame() %>%
        rownames_to_column("term")) %>%
        bind_rows(., .id = "roi")
    write_csv(disp_df, sprintf("%s/model_disp_coef.csv", outputdirtables))



    ### mixture model
    library(brms)

    # Define mixture of two binomial components
    mix <- mixture(binomial(), binomial())

    fit <- brm(
        bf(meth | trials(20) ~ 1), # Intercept-only for now
        family = mix,
        data = dfsall[["L1HS"]] %>%
            mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
            left_join(sample_table %>%
                dplyr::select(sample_name, condition)) %>%
            dplyr::select(gene_id, sample_name, meth),
        chains = 4,
        iter = 4000,
        #   control = list(adapt_delta = 0.95),
        cores = 4
    )

    library(flexmix)
    fitflex <- flexmix(
        cbind(meth, 20 - meth) ~ 1,
        data = dfsall[["L1HS"]] %>%
            mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
            left_join(sample_table %>%
                dplyr::select(sample_name, condition)) %>%
            dplyr::select(gene_id, sample_name, meth) %>%
            filter(gene_id == "L1HS_4q28.3_9"),
        k = 3, # number of components
        model = FLXglm(family = "binomial")
    )

    summary(fitflex)

    parameters(fitflex)

    posterior_probs <- posterior(fitflex)

    mixturemodelk2 <- function(df) {
        fit <- flexmix(
            cbind(meth, 20 - meth) ~ 1,
            data = df %>%
                ungroup() %>%
                mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
                left_join(sample_table %>%
                    dplyr::select(sample_name, condition)) %>%
                dplyr::select(gene_id, sample_name, meth),
            k = 2, # number of components
            model = FLXglm(family = "binomial")
        )
        tibble(proportions = prior(fit), mean_meth = unname(plogis(parameters(fit))))
    }


    mix_res <- bind_rows(map(dfsall, mixturemodelk2), .id = "type")

    dattemp <- dfsall[["imprinted"]] %>%
        mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
        left_join(sample_table %>%
            dplyr::select(sample_name, condition)) %>%
        dplyr::select(gene_id, sample_name, meth) %>%
        filter(gene_id == "L1HS_4q28.3_9")


    mixturemodelk2(dfsall[["imprinted"]])
}






read_analysis_alt_regions <- function(
    readscg,
    mod_code_var = "m",
    required_number_cg = 0.75,
    meth_thresholds = c(0.1, 0.25, 0.5),
    context = "CpG"
    region = "L1HS_FL"
) {
    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, region, required_fraction_of_total_cg)
    dir.create(outputdirtables, recursive = TRUE)

    readsdf1 <- readscg %>%
        left_join(rmann %>%
            dplyr::select(gene_id, start, end, strand, rte_length_req, intactness_req) %>%
            dplyr::rename(element_strand = strand, element_start = start, element_end = end)) %>%
        filter(rte_length_req == "FL")


    by_cpg_l <- list()
    by_read_l <- list()
    by_sample_l <- list()
    by_gene_id_l <- list()

    for (region_of_interest in regions_of_interest) {
        roistart <- region_of_interest[1]
        roiend <- region_of_interest[2]
        roistring <- paste0(roistart, "to", roiend)

        numCGneeded <- ceiling(length(cg_indices[(cg_indices <= roiend) & (cg_indices >= roistart)]) * required_fraction_of_total_cg)

        utr1 <- readsdf1 %>%
            filter(mod_code == mod_code_var) %>%
            filter(case_when(
                element_strand == "+" ~ (start > element_start + roistart) & (start < element_start + roiend),
                element_strand == "-" ~ (start > element_end - roiend) & (start < element_end - roistart)
            )) %>%
            dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

        by_cpg_temp <- utr1 %>%
            group_by(gene_id, read_id, condition, sample) %>%
            mutate(num_cpgs_in_read = n()) %>%
            mutate(fraction_meth = mean(mod_indicator)) %>%
            relocate(gene_id) %>%
            ungroup()

        by_read_temp <- by_cpg_temp %>%
            filter(num_cpgs_in_read >= numCGneeded) %>%
            group_by(read_id, gene_id, sample, condition, region) %>%
            summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read), strand = dplyr::first(element_strand), numCGneeded = dplyr::first(numCGneeded)) %>%
            ungroup()

        by_cpg_temp$subset <- as.character(roistring)
        by_read_temp$subset <- as.character(roistring)

        by_cpg_l[[as.character(roistring)]] <- by_cpg_temp
        by_read_l[[as.character(roistring)]] <- by_read_temp

        for (meth_threshold in meth_thresholds) {
            subset_threshold <- paste0(roistring, "_", meth_threshold)
            by_sample_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth > meth_threshold, 0, 1)) %>%
                group_by(sample, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                group_by(condition, region) %>%
                mutate(meanProp = mean(propUnmeth))
            by_sample_temp$meth_threshold <- meth_threshold
            by_sample_temp$subset <- as.character(roistring)
            by_sample_temp$subset_threshold <- subset_threshold

            by_gene_id_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth >= meth_threshold, 0, 1)) %>%
                group_by(sample, gene_id, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                ungroup() %>%
                group_by(gene_id) %>%
                mutate(group_size = n()) %>%
                ungroup()
            by_gene_id_temp$meth_threshold <- meth_threshold
            by_gene_id_temp$subset <- as.character(roistring)
            by_gene_id_temp$subset_threshold <- subset_threshold

            by_sample_l[[subset_threshold]] <- by_sample_temp
            by_gene_id_l[[subset_threshold]] <- by_gene_id_temp
        }
    }

    by_cpg <- purrr::reduce(by_cpg_l, bind_rows)
    by_read <- purrr::reduce(by_read_l, bind_rows)
    by_gene_id <- purrr::reduce(by_gene_id_l, bind_rows)
    by_sample <- purrr::reduce(by_sample_l, bind_rows)

    dir.create(dirname(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)), recursive = TRUE)
    by_cpg %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region))
    by_read %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region))
    by_gene_id %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region))
    by_sample %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region))

    # by_cpg <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_read <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_gene_id <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_sample <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region)) %>%
    #         mutate(condition = factor(condition, levels = conf$levels)) %>%
    #         mutate(sample = factor(sample, levels = conf$samples))

    get_read_ecdf <- function(df, subset_val, breakpoints, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by specified variable
            summarise(
                percent_below = list(ecdf(fraction_meth)(breakpoints) * 100), # Compute ECDF
                .groups = "drop"
            ) %>%
            unnest_longer(percent_below) %>%
            mutate(
                threshold = rep(breakpoints, times = n() / length(breakpoints)), # Expand breakpoints
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), threshold, percent_below) # Keep relevant columns
    }


    get_read_quantiles <- function(df, subset_val, probs, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by sample_name or another variable
            summarise(
                quantiles = list(quantile(fraction_meth, probs = probs)), # Store as list
                .groups = "drop"
            ) %>%
            unnest_longer(quantiles) %>%
            mutate(
                quantile = rep(probs, times = n() / length(probs)), # Expand probs for each group
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), quantile, mean_meth = quantiles) # Keep relevant columns
    }
    subsets <- c("0to909", "0to500", "0to328")
    breakpoints <- seq(0, 1, 0.05)

    ecdf_reads <- map_dfr(subsets, ~ get_read_ecdf(by_cpg, .x, breakpoints, "sample"))
    quantile_reads <- map_dfr(subsets, ~ get_read_quantiles(by_cpg, .x, breakpoints, "sample"))


    ecdf_reads %>% write_mycsv(sprintf("%s/read_ecdf.csv", outputdirtables))
    quantile_reads %>% write_mycsv(sprintf("%s/read_quantiles.csv", outputdirtables))

    ecdf_reads_acrosssamplemean <- ecdf_reads %>%
        group_by(subset, threshold) %>%
        summarise(percent_below = mean(percent_below))
    quantile_reads_acrosssamplemean <- quantile_reads %>%
        group_by(subset, quantile) %>%
        summarise(mean_meth = mean(mean_meth))

    ecdf_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_ecdf_acrosssamplemean.csv", outputdirtables))
    quantile_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_quantiles_acrosssamplemean.csv", outputdirtables))


    # Example ECDF plot function
    plot_ecdf <- function(ecdf_data) {
        ggplot(ecdf_data, aes(x = threshold, y = percent_below, color = as.factor(subset))) +
            geom_line() + # ECDF curve
            geom_point() + # Add points for clarity
            labs(
                title = "Empirical CDF of Fraction Methylation",
                x = "Fraction Methylation Threshold",
                y = "Percent Below Threshold",
                color = "Subset"
            ) +
            theme_minimal()
    }


    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)


    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        # stats <- by_sample %>%
        #     filter(subset != "400to600") %>%
        #     ungroup() %>%
        #     left_join(sample_table) %>%
        #     group_split(subset_threshold) %>%
        #     set_names(unique(by_sample %>% filter(subset != "400to600") %$% subset_threshold)) %>%
        #     map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
        #     imap_dfr(~ .x %>% mutate(subset = .y))
        stats_list <- list()
        i <- 1
        for (subset in unique(by_cpg$subset)) {
            for (threshold in unique(by_sample$meth_threshold)) {
                by_cpg_tmp <- by_cpg %>%
                    filter(subset == !!subset) %>%
                    mutate(unmeth = ifelse(fraction_meth > threshold, 0, 1)) %>%
                    dplyr::rename(sample_name = sample) %>%
                    group_by(sample_name, condition, gene_id) %>%
                    summarise(unmeth = sum(unmeth), total = n()) %>%
                    ungroup() %>%
                    left_join(sample_table) %>%
                    mutate(age_z = as.numeric(scale(age)))
                model_tmp <- glmmTMB(
                    cbind(unmeth, total - unmeth) ~
                        condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
                    data = by_cpg_tmp,
                    family = binomial()
                )
                res_tmp <- broom::tidy(model_tmp) %>%
                    mutate(threshold = !!threshold) %>%
                    mutate(subset = !!subset)
                stats_list[[i]] <- res_tmp
                i <- i + 1
            }
        }
        stats <- purrr::reduce(stats_list, bind_rows)
        stats_padj <- stats %>%
            filter(term == "conditionAD") %>%
            filter(subset != "400to600") %>%
            mutate(padj = p.adjust(p.value, method = "fdr"))
        statswpadj <- stats %>% left_join(stats_padj)
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = statswpadj)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }
    p <- by_sample %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        stats <- by_sample %>%
            ungroup() %>%
            left_join(sample_table) %>%
            group_split(subset_threshold) %>%
            set_names(unique(by_sample$subset_threshold)) %>%
            map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
            imap_dfr(~ .x %>% mutate(subset = .y))
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = stats)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }

    p <- by_read %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)


    # Get group-wise mean fraction_meth
    # beta_overlays <- by_read %>%
    # filter(subset != "400to600") %>%
    # group_by(condition, subset) %>%
    # summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop")

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == "AD"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == "CTRL"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))
    p <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == "AD"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == "CTRL"),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay_L1HS_4q28.3_9.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        filter(subset == "328") %$% condition %>%
        table()


    p <- by_gene_id %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        filter(meth_threshold == "0.5") %>%
        filter(subset == "0to909") %>%
        group_by(gene_id) %>%
        mutate(samples_detected_per_element = n()) %>%
        ungroup() %>%
        filter(samples_detected_per_element > 11) %>%
        distinct() %>%
        # group_by(direction) %>%
        tidyHeatmap::heatmap(gene_id, sample, propUnmeth,
            cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE,
            show_row_dend = FALSE,
            palette_value = circlize::colorRamp2(
                c(0, 0.00001, seq(0.1, 1, length.out = 3)),
                c("black", rev(RColorBrewer::brewer.pal(4, "Oranges")))
            )
        )
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_gene_id_consistent_across_samples.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)

    #
    df_filtered <- by_gene_id %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        filter(meth_threshold == "0.5", subset == "0to909") %>%
        group_by(gene_id) %>%
        mutate(samples_detected_per_element = n()) %>%
        ungroup() %>%
        filter(samples_detected_per_element > 11) %>%
        distinct()
    df_sorted <- df_filtered %>%
        group_by(sample) %>%
        arrange(propUnmeth) %>%
        mutate(gene_rank = row_number()) %>%
        ungroup()
    p <- ggplot(df_sorted, aes(x = sample, y = gene_rank, fill = propUnmeth)) +
        geom_tile() +
        scale_fill_gradientn(
            colours = c("lightblue", rev(RColorBrewer::brewer.pal(4, "Oranges"))),
            values = rescale(c(0, 0.00001, seq(0.1, 1, length.out = 3))),
            name = "propUnmeth"
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) + # No padding; bottom = high rank
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(x = "sample", y = "genes (sorted by propUnmeth per sample)") +
        mtclosed +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p, raster = TRUE)

    #

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmann) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    named_group_split <- function(.tbl, ...) {
        grouped <- group_by(.tbl, ...)
        names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

        grouped %>%
            group_split() %>%
            rlang::set_names(names)
    }

    groups_needed <- length(condition1samples) + 1
    groups_needed <- 10
    n_under_consideration <- by_gene_id %>%
        ungroup() %>%
        mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
        left_join(sample_table) %>%
        filter(group_size >= groups_needed) %>%
        group_by(subset_threshold) %>%
        dplyr::select(subset_threshold, gene_id) %>%
        distinct() %>%
        summarise(n_under_consideration = n())
    group_size_df <- by_gene_id %>% dplyr::select(gene_id, group_size, subset, meth_threshold)

    dat_tmp <- meth_thresholds %>%
        map(~ by_cpg %>%
            mutate(unmeth = ifelse(fraction_meth > .x, 0, 1), meth_threshold = .x) %>%
            filter(subset != "400to600") %>%
            group_by(sample, gene_id, subset) %>%
            mutate(total_read = n()) %>%
            ungroup() %>%
            group_by(sample, gene_id, subset, meth_threshold) %>%
            summarise(unmeth = sum(unmeth), total_read = dplyr::first(total_read)) %>%
            ungroup()) %>%
        list_rbind() %>%
        dplyr::rename(sample_name = sample) %>%
        mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
        left_join(group_size_df %>% distinct()) %>%
        left_join(sample_table) %>%
        filter(group_size >= groups_needed) %>%
        mutate(age_z = as.numeric(scale(age)))




    # Create an empty list to hold tidy model results
    stats_list <- list()

    # Split data by group
    dat_groups <- named_group_split(dat_tmp, split_var)
    library(broom.mixed)
    # Loop through each group
    for (i in seq_along(dat_groups)) {
        print(i)
        group_data <- dat_groups[[i]]
        group_name <- names(dat_groups)[i]

        print(group_name)

        # Check if unmeth values are constant (e.g. all 0s or all the same)
        if (length(unique(group_data$unmeth)) == 1 | sum(group_data$unmeth > 0) < 2) {
            print("no info")
            tidy_model <- tibble(
                term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                subset = group_name,
                note = "No variation in unmeth"
            )
        } else {
            model <- glmmTMB(
                cbind(unmeth, total_read - unmeth) ~
                    condition + sex + age_z + (1 | sample_name),
                data = group_data,
                family = binomial()
            )
            model <- glmmTMB(
                cbind(unmeth, total_read - unmeth) ~
                    condition + sex + age_z + (1 | sample_name),
                data = group_data,
                family = binomial()
            )

            tidy_model <- broom::tidy(model) %>%
                mutate(subset = group_name)
        }
        stats_list[[i]] <- tidy_model
    }

    # Combine the results into a single data frame
    stats <- purrr::reduce(stats_list, bind_rows) %>%
        tidyr::separate(subset, into = c("gene_id", "subset", "meth_threshold"), sep = "/", convert = TRUE)


    write_csv(stats, sprintf("ldna/results/%s/tables/reads_new/%s_%s/by_gene_new.csv", params$mod_code, region, required_fraction_of_total_cg))

    gene_condition_stats <- stats %>%
        filter(term == "conditionAD") %>%
        mutate(p.value = case_when(
            is.nan(p.value) ~ 1,
            TRUE ~ p.value
        )) %>%
        mutate(statistic = case_when(
            is.nan(statistic) ~ 0,
            TRUE ~ statistic
        ))

    tryCatch(
        {
            p <- stats %>%
                mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                filter(p.value <= 0.05) %>%
                filter(grepl("condition", term)) %>%
                mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                count(meth_threshold, subset, dir_stat) %>%
                complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                geom_col(position = "dodge", color = "black") +
                facet_wrap(vars(subset), nrow = 1) +
                ggtitle(sprintf("Significant Loci")) +
                labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                mtclosed +
                anchorbar +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 5.5, 4, pl = p)

            p <- stats %>%
                filter(subset != "400to600") %>%
                mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                filter(p.value <= 0.05) %>%
                filter(grepl("condition", term)) %>%
                mutate(meth_threshold = factor(as.character(meth_threshold), levels = c("0.1", "0.25", "0.5"))) %>%
                mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                count(meth_threshold, subset, dir_stat) %>%
                complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                geom_col(position = "dodge", color = "black") +
                facet_wrap(vars(subset), nrow = 1) +
                ggtitle(sprintf("Significant Loci")) +
                labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                mtclosed +
                anchorbar +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 4, pl = p)
        },
        error = function(e) {
            print("no sig de")
        }
    )

    dispersion_models <- list()
    for (subsetofinterest in by_cpg$subset %>%
        unique() %>%
        grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
        dispersion_model <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = by_cpg %>% filter(subset == subsetofinterest) %>%
                mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
        )
        dispersion_models[[subsetofinterest]] <- dispersion_model
    }

    res_text <- capture.output(map(dispersion_models, summary))
    writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))
    summary_disp <- summary(mod)$dispersion

    # # Extract estimates and standard errors
    # log_phi_control <- 1.95678
    # se_log_phi_control <- 0.01539

    # log_phi_diff_AD <- -0.32685
    # se_log_phi_diff_AD <- 0.01990

    # # Compute log(φ) for AD and its SE
    # log_phi_AD <- log_phi_control + log_phi_diff_AD
    # se_log_phi_AD <- sqrt(se_log_phi_control^2 + se_log_phi_diff_AD^2)

    # # 95% CI on log-scale
    # z <- 1.96
    # ci_log_phi_control <- log_phi_control + c(-1, 1) * z * se_log_phi_control
    # ci_log_phi_AD <- log_phi_AD + c(-1, 1) * z * se_log_phi_AD

    # # Exponentiate to get CIs for φ
    # phi_control <- exp(log_phi_control)
    # phi_AD <- exp(log_phi_AD)

    # ci_phi_control <- exp(ci_log_phi_control)
    # ci_phi_AD <- exp(ci_log_phi_AD)


    # dispersion_models_withreadid <- list()
    # for (subsetofinterest in by_cpg$subset %>%
    #     unique() %>%
    #     grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
    #     dispersion_model <- glmmTMB(
    #         cbind(meth, unmeth) ~ 1 + (1 | read_id),
    #         dispformula = ~condition,
    #         family = glmmTMB::betabinomial(),
    #         data = by_cpg %>% filter(subset == subsetofinterest) %>%
    #             mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
    #     )
    #     dispersion_models[[subsetofinterest]] <- dispersion_model
    # }

    # res_text <- capture.output(map(dispersion_models, summary))
    # writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))


    # dat <- by_cpg %>%
    #     filter(subset == "0to909") %>%
    #     mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
    # mod0 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
    #     data = dat, family = binomial()
    # )

    # mod05 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = dat, family = glmmTMB::betabinomial()
    # )
    # mod05disp <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     dispformula = ~1,
    #     data = dat, family = glmmTMB::betabinomial()
    # )
    # datmorebinom <- dat %>%
    #     filter(fraction_meth < 0.95) %>%
    #     filter(fraction_meth > 0.35)
    # mod05morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = datmorebinom, family = glmmTMB::betabinomial()
    # )
    # mod05binom_morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
    #     data = datmorebinom, family = binomial()
    # )


    # dispersion_model_morebinomdat <- glmmTMB(
    #     cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
    #     dispformula = ~condition,
    #     family = glmmTMB::betabinomial(),
    #     data = datmorebinom
    # )

    # anova(mod0, mod05)
    # anova(mod05binom_morebinom, mod05morebinom)
}
