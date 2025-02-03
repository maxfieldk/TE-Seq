module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
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


samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]
conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

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

perl1hs_5utr_region <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE) %>% mutate(region = ordered(region, levels = c("328", "500", "909")))
perl1hs_5utr_region$sample <- factor(perl1hs_5utr_region$sample, levels = conf$samples)
perl1hs_5utr_region$condition <- factor(perl1hs_5utr_region$condition, levels = conf$levels)

cpg_islands <- rtracklayer::import(conf$cpg_islands)
cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)

readscg <- read_delim(sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE) %>%
    mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels))

if (conf$single_condition == "no") {
    dmrs <- read_delim(params$dmrs, delim = "\t", col_names = TRUE) %>% filter(dmr_type != "t05CG10")
    dmls <- read_delim(params$dmls, delim = "\t", col_names = TRUE) %>% filter(fdrs <= 0.05)

    dmrsgr <- GRanges(dmrs)
    dmlsgr <- GRanges(
        seqnames = dmls$chr,
        ranges = IRanges(start = dmls$pos, end = dmls$pos),
        stat = dmls$stat,
        pval = dmls$pvals,
        fdr = dmls$fdrs,
        direction = dmls$direction
    )
}

##################################### DML / DMR analysis
if (conf$single_condition == "no") {
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
            group_by(chr, direction) %>%
            summarize(n = n())
        dmrlocdf$chr <- factor(dmrlocdf$chr, levels = chromosomes)
        p <- ggplot(data = dmrlocdf) +
            geom_col(aes(y = chr, x = n, fill = direction), position = "dodge", color = "black") +
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
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
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



if (conf$single_condition == "no") {
    dir.create(sprintf("ldna/results/%s/plots/figs", params$mod_code))
    dmrtypes <- dmrs$dmr_type %>% unique()


    flRTEpromoterlong <- pivot_longer(data = flRTEpromoter, cols = dmrtypes, names_to = "dmr_type", values_to = "direction")

    pff <- flRTEpromoterlong %>%
        filter(dmr_type != "t05CG10") %>%
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
if (conf$single_condition == "no") {
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
if (conf$single_condition == "no") {
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
if (conf$single_condition == "no") {
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
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
if (conf$single_condition == "no") {
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
if (conf$single_condition == "no") {
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
    left_join(sample_table %>% dplyr::rename(sample = sample_name)) %>%
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
if (conf$single_condition == "no") {
    stats <- perelementdf_promoters %>%
        compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), 10, 6, sf = stats)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", params$mod_code), raster = TRUE, 12, 4)
}

p <- perelementdf_promoters %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    left_join(sample_table %>% dplyr::rename(sample = sample_name)) %>%
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
if (conf$single_condition == "no") {
    perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(sample) %>%
        summarize(median = median(mean_meth)) %>%
        left_join(sample_table %>% dplyr::rename(sample = sample_name)) %>%
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

if (conf$single_condition == "no") {
    pfl1 <- perelementdf_promoters %>% filter(grepl("^L1", rte_subfamily))
    p <- pfl1 %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = CTRL - AD) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        group_by(rte_subfamily) %>%
        mutate(rank_change = row_number()) %>%
        mutate(top_change = ifelse(rank_change <= 10, "Top", "NotTop")) %>%
        arrange(abs_dif) %>%
        ungroup() %>%
        ggpaired(cond1 = "CTRL", cond2 = "AD", line.color = "top_change", alpha = "top_change", facet.by = "rte_subfamily") +
        scale_alpha_manual(values = c(1, 0.5)) +
        scale_color_manual(values = c("Top" = "red", "NotTop" = "grey")) +
        mtclosedgridh
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_L1s.pdf", params$mod_code), 14, 6, raster = TRUE)

    l1hs_paired_dif_frame <- pfl1 %>%
        filter(rte_subfamily == "L1HS") %>%
        pivot_longer(cols = dmrtypes, names_to = "dmr_type", values_to = "direction") %>%
        filter(dmr_type != "t05CG10") %>%
        mutate(direction_threshold = ifelse(is.na(direction), "NS", paste0(direction, "_", gsub("t", "", dmr_type)))) %>%
        filter(!(dmr_type == "t01" & is.na(direction))) %>%
        group_by(gene_id, rte_subfamily, condition, dmr_type, direction_threshold) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = CTRL - AD) %>%
        mutate(abs_dif = abs(dif)) %>%
        arrange(-abs_dif) %>%
        mutate(
            x_1 = factor("CTRL", levels = c("CTRL", "AD")),
            x_2 = factor("AD", levels = c("CTRL", "AD")),
            y_1 = CTRL,
            y_2 = AD
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
    mysaveandstore(pl = p, fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_L1s_%s.pdf", params$mod_code, "all"), 5, 4)

    for (dmrtype in dmrs$dmr_type %>% unique()) {
        p <- pfl1 %>%
            mutate(!!sym(dmrtype) := ifelse(is.na(!!sym(dmrtype)), "NS", !!sym(dmrtype))) %>%
            group_by(gene_id, rte_subfamily, condition, !!sym(dmrtype)) %>%
            summarize(mean_meth = mean(mean_meth)) %>%
            pivot_wider(names_from = condition, values_from = mean_meth) %>%
            ungroup() %>%
            mutate(dif = CTRL - AD) %>%
            mutate(abs_dif = abs(dif)) %>%
            arrange(-abs_dif) %>%
            group_by(rte_subfamily) %>%
            arrange(abs_dif) %>%
            ungroup() %>%
            ggpaired(cond1 = "CTRL", cond2 = "AD", line.color = dmrtype, alpha = dmrtype, facet.by = "rte_subfamily") +
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
            mutate(dif = CTRL - AD) %>%
            mutate(abs_dif = abs(dif)) %>%
            arrange(-abs_dif) %>%
            group_by(rte_subfamily) %>%
            arrange(abs_dif) %>%
            ungroup() %>%
            ggpaired(cond1 = "CTRL", cond2 = "AD", line.color = dmrtype, alpha = 0.85, ylab = "L1HS 5UTR Methylation") +
            scale_color_manual(values = c("Hypo" = "red", "Hyper" = "blue", "NS" = "grey")) +
            mtclosedgridh
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_paired_promoters_l1hs_%s.pdf", params$mod_code, dmrtype), 4, 4, raster = FALSE)
    }

    top_l1hs_movers <- pfl1 %>%
        group_by(gene_id, rte_subfamily, condition) %>%
        summarize(mean_meth = mean(mean_meth)) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        ungroup() %>%
        mutate(dif = CTRL - AD) %>%
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
        mutate(dif = CTRL - AD) %>%
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
    ggplot(aes(x = region, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("L1HS CpG Methylation") +
    mtopen +
    scale_conditions
tryCatch(
    {
        res909 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 909) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 909) %>%
            mutate(model_type = "no_interaction")
        res500 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 500) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 500) %>%
            mutate(model_type = "no_interaction")
        res328 <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 328) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 328) %>%
            mutate(model_type = "no_interaction")

        res909interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 909) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age + condition * sex, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 909) %>%
            mutate(model_type = "interaction")
        res500interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 500) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age + condition * sex, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 500) %>%
            mutate(model_type = "interaction")
        res328interaction <- pf %>%
            group_by(sample, region) %>%
            summarise(pctM = mean(mean_meth)) %>%
            filter(region == 328) %>%
            left_join(sample_table %>% mutate(sample = sample_name)) %>%
            lm(pctM ~ condition + sex + age + condition * sex, data = .) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(region = 328) %>%
            mutate(model_type = "interaction")
        stats <- bind_rows(res909, res500, res328, res909interaction, res500interaction, res328interaction)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4, sf = stats)
    },
    error = function(e) {
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/l1hs_boxplot_5utr_regions.pdf", params$mod_code), 5, 4)
    }
)



########## PCA
library(PCAtools)
pcaframe <- perl1hs_5utr_region %>%
    filter(region == "909") %>%
    dplyr::select(sample, mean_meth, gene_id) %>%
    pivot_wider(names_from = gene_id, values_from = mean_meth) %>%
    arrange(sample) %>%
    column_to_rownames(var = "sample") %>%
    as.matrix() %>%
    t()
pcaObj <- pca(pcaframe, center = TRUE, scale = FALSE, metadata = sample_table %>% column_to_rownames(var = "sample_name"))

p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/screeplot.pdf", params$mod_code), 5, 4)


p <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 2
) + mtopen
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/loadings.pdf", params$mod_code), 5, 4)

p <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = "sex", colby = "condition",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) + mtopen + scale_conditions
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/biplot.pdf", params$mod_code), 5, 5)


library(ggh4x)
library(ggnewscale)

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
    geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_palette +
    mtopen
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar1.pdf", params$mod_code), 5, 4, sf = stats)

p <- pf %>%
    mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
    ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
    geom_col() +
    scale_conditions +
    new_scale_fill() +
    geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_palette +
    mtopen
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage.pdf", params$mod_code), 5, 4, sf = stats)

p <- pf %>%
    mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
    ggplot(aes(y = sample_name, x = mean_meth, fill = condition)) +
    geom_col() +
    scale_conditions +
    new_scale_fill() +
    geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_palette +
    mtopen
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_bar_withage_ordered.pdf", params$mod_code), 5, 4, sf = stats)


p <- pf %>%
    mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
    ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = sex)) +
    geom_point(size = 3) +
    scale_conditions +
    geom_vline(xintercept = pf %>% filter(condition == "AD") %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
    geom_text_repel(aes(label = apoe)) +
    # new_scale_fill() +
    # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_conditions +
    mtopen
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_ordered.pdf", params$mod_code), 5, 4, sf = stats)

p <- pf %>%
    mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), mean_meth)) %>%
    ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = sex)) +
    geom_point(size = 3) +
    scale_conditions +
    geom_vline(xintercept = pf %>% filter(condition == "AD") %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
    geom_text_repel(aes(label = apoe)) +
    # new_scale_fill() +
    # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_conditions +
    mtopen
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_orderedmeth.pdf", params$mod_code), 5, 4, sf = stats)

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


if (conf$single_condition == "no") {
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


    library(ggnewscale)
    library(patchwork)
    l1hsflmethgr <- rtedf %>%
        filter(rte_subfamily == "L1HS")
    l1hsflmethgr <- l1hsflmethgr %>%
        mutate(rel_start = start - rte_start) %>%
        mutate(rel_end = end - rte_start)
    write_delim(l1hsflmethgr, sprintf("ldna/Rintermediates/%s/l1hsfldf.tsv", params$mod_code), col_names = TRUE)

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
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    ggplot() +
                    geom_point(aes(x = rel_start, y = rM, fill = sample, color = sample)) +
                    scale_samples_unique +
                    labs(y = "Methylation Rolling Mean") +
                    mtclosed +
                    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
                # p2 <- color_intervals %>%
                #     ggplot() +
                #     geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
                #     geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
                #     geom_text(aes(x = -200 + ((start + end) / 2), y = 1.5, label = feature)) +
                #     coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
                #     ggtitle(element) +
                #     scale_fill_paletteer_d("dutchmasters::milkmaid") +
                #     theme_map() +
                #     theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
                #     scale_y_continuous(expand = c(0, 0.4)) +
                #     theme(legend.position = "none")
            } else {
                modifier <- rmann %>% filter(gene_id == element) %$% end
                color_intervals <- element_anatomy %>%
                    filter(!(feature %in% c("EN", "RT"))) %>%
                    filter(gene_id == element) %>%
                    mutate(across(where(is.numeric), ~ modifier - .))
                pf <- pf_neg
                p1 <- pf %>%
                    filter(gene_id == element) %>%
                    ggplot() +
                    geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
                    scale_samples_unique +
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

if (conf$single_condition == "no") {
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


read_analysis <- function(
    readsdf,
    region = "L1HS_intactness_req_ALL",
    mod_code_var = "m",
    region_of_interest_from_start = 909,
    required_fraction_of_total_cg = 0.8,
    fraction_meth_threshold = 0.5,
    context = "CG") {
    readsdf1 <- readsdf %>% left_join(rmann %>% dplyr::select(gene_id, start, end, strand, rte_length_req) %>% dplyr::rename(element_strand = strand, element_start = start, element_end = end))
    readsdf2 <- readsdf1 %>% filter(rte_length_req == "FL")
    utr1 <- readsdf1 %>%
        filter(mod_code == mod_code_var) %>%
        filter(case_when(
            element_strand == "+" ~ (start > element_start) & (start < element_start + region_of_interest_from_start),
            element_strand == "-" ~ (start > element_end - region_of_interest_from_start) & (start < element_end)
        )) %>%
        dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

    utr <- utr1 %>%
        group_by(gene_id, read_id, condition, sample) %>%
        mutate(num_cpgs_in_read = n()) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        relocate(gene_id) %>%
        ungroup()


    outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
    subfam <- "L1HS"
    consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
    consensus_ss <- readDNAStringSet(consensus_path)
    cg_indices <- consensus_ss %>%
        vmatchPattern(pattern = "CG") %>%
        start() %>%
        unlist() %>%
        as.numeric()


    numCGneeded <- ceiling(length(cg_indices[cg_indices <= region_of_interest_from_start]) * required_fraction_of_total_cg)

    # write_delim(utr, sprintf("ldna/Rintermediates/%s/%s_to_%s_considering_reads_%s_fraction_%s_%s_%s_reads.tsv", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), delim = "\t")

    # p <- utr %>%
    #     group_by(gene_id, read_id, condition, region) %>%
    #     summarise(nc = max(read_span), strand = dplyr::first(element_strand)) %>%
    #     ggplot() +
    #     geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
    #     facet_wrap(vars(region)) +
    #     mtclosed +
    #     scale_conditions
    # mysaveandstore(sprintf("ldna/results/%s/plots/reads/read_span_distribution_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 5, 5, pl = p)


    # p <- utr %>%
    #     filter(read_span > 250) %>%
    #     group_by(gene_id, read_id, condition, region) %>%
    #     summarise(nc = max(num_cpgs_in_read), strand = dplyr::first(element_strand)) %>%
    #     ggplot() +
    #     geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
    #     facet_wrap(vars(region)) +
    #     mtclosed +
    #     scale_conditions
    # mysaveandstore(sprintf("ldna/results/%s/plots/reads/read_num_cpg_distribution_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 5, 5, pl = p)



    # p <- utr %>%
    #     filter(read_span > 250) %>%
    #     group_by(read_id, condition, region) %>%
    #     summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
    #     ggplot() +
    #     geom_point(aes(x = read_span, y = fraction_meth, color = condition)) +
    #     facet_wrap(vars(region)) +
    #     mtclosed +
    #     scale_conditions
    # mysaveandstore(sprintf("ldna/results/%s/plots/reads/fraction_meth_distribution_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 5, 5, pl = p)

    # utr %>%
    #     group_by(read_id, sample, condition, region) %>%
    #     summarise(fraction_meth = dplyr::first(fraction_meth), strand = dplyr::first(element_strand)) %>%
    #     ungroup()
    aa <- utr %>%
        filter(num_cpgs_in_read >= numCGneeded) %>%
        group_by(read_id, sample, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), strand = dplyr::first(element_strand)) %>%
        ungroup()
    aa %$% fraction_meth %>%
        quantile() %>%
        tibble() %>%
        tibble() %>%
        write_delim(sprintf("ldna/results/%s/plots/reads/number_of_reads_assessed_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s_read_methylation_quartiles.txt", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context))

    # dir.create("bayes/data", recursive = TRUE)
    # aa %>%
    #     write_delim(sprintf("bayes/data/reads_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.tsv", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context),
    #         col_names = TRUE
    #     )

    p <- aa %>%
        group_by(sample) %>%
        summarise(number_of_reads = n()) %>%
        ggplot(aes(x = sample, y = number_of_reads)) +
        geom_col() +
        coord_flip() +
        mtopen
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/reads/number_of_reads_assessed_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 5, 5, pl = p)

    ab <- aa %>%
        mutate(unmeth = ifelse(fraction_meth > fraction_meth_threshold, 0, 1)) %>%
        group_by(sample, region, condition) %>%
        summarise(propUnmeth = mean(unmeth)) %>%
        group_by(condition, region) %>%
        mutate(meanProp = mean(propUnmeth))

    p <- ab %>%
        ggplot(aes(x = condition)) +
        stat_summary(aes(y = propUnmeth, fill = condition), color = "black", fun = "mean", geom = "bar") +
        geom_point(aes(y = propUnmeth)) +
        facet_wrap(vars(region)) +
        geom_pwc(aes(x = condition, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "", y = sprintf("Reads Fraction < %s methylated", fraction_meth_threshold)) +
        ggtitle(sprintf("1-%s Methylation", region_of_interest_from_start)) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    stats <- ab %>%
        ungroup() %>%
        left_join(sample_table %>% dplyr::rename(sample = sample_name)) %>%
        lm(propUnmeth ~ condition + sex + age, .) %>%
        summary() %>%
        broom::tidy()
    mysaveandstore(sprintf("ldna/results/%s/plots/reads/barplot_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 4, 4, pl = p, sf = stats)

    p <- ab %>%
        ggplot(aes(x = condition)) +
        stat_summary(aes(y = propUnmeth, fill = condition), color = "black", fun = "mean", geom = "bar") +
        geom_point(aes(y = propUnmeth)) +
        geom_pwc(aes(x = condition, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "", y = sprintf("Reads Fraction < %s methylated", fraction_meth_threshold)) +
        ggtitle(sprintf("1-%s Methylation", region_of_interest_from_start)) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads/barplot_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s_nosubtitle.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 4, 4, pl = p)

    p <- aa %>%
        mutate(condition = factor(condition, levels = c("AD", "CTRL"))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        labs(x = "Pct CpG Methylation per Read") +
        ggtitle(sprintf("1-%s Methylation", region_of_interest_from_start)) +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads/fraction_meth_density_distribution_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 4, 4, pl = p)


    # aa <- utr %>%
    #     filter(read_span > region_of_interest_from_start * required_fraction_of_total_cg) %>%
    #     group_by(read_id, sample, condition, region, gene_id) %>%
    #     summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
    #     ungroup()

    # p <- aa %>%
    #     group_by(sample, gene_id) %>%
    #     summarise(number_of_reads = n()) %>%
    #     ggplot(aes(x = sample, y = number_of_reads)) +
    #     geom_col() +
    #     facet_wrap(~gene_id) +
    #     coord_flip() +
    #     mtopen
    # mysaveandstore(fn = sprintf("ldna/results/%s/plots/reads/number_of_reads_assessed_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s_split_by_geneid.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 30, 30, pl = p)

    # ab <- aa %>%
    #     mutate(unmeth = ifelse(fraction_meth > fraction_meth_threshold, 0, 1)) %>%
    #     group_by(sample, region, condition, gene_id) %>%
    #     summarise(propUnmeth = mean(unmeth)) %>%
    #     group_by(condition, region, gene_id) %>%
    #     mutate(meanProp = mean(propUnmeth))

    # p <- ab %>%
    #     ggplot(aes(x = condition)) +
    #     stat_summary(aes(y = propUnmeth, fill = condition), color = "black", fun = "mean", geom = "bar") +
    #     geom_point(aes(y = propUnmeth)) +
    #     facet_wrap(~gene_id) +
    #     geom_pwc(aes(x = condition, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
    #     labs(x = "", y = sprintf("Pct Reads < %s methylated", fraction_meth_threshold)) +
    #     ggtitle(sprintf("1-%s Methylation", region_of_interest_from_start)) +
    #     mtclosedgridh +
    #     scale_conditions +
    #     anchorbar
    # mysaveandstore(sprintf("ldna/results/%s/plots/reads/barplot_%s_to_%s_considering_reads_%s_fraction_%s_%s_%s_split_by_geneid.pdf", params$mod_code, region, region_of_interest_from_start, required_fraction_of_total_cg, fraction_meth_threshold, mod_code_var, context), 30, 30, pl = p)
}

tryCatch(
    {
        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 909,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.5,
            context = "CpG"
        )

        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 909,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.25,
            context = "CpG"
        )

        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 500,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.5,
            context = "CpG"
        )

        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 500,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.25,
            context = "CpG"
        )

        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 328,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.5,
            context = "CpG"
        )

        read_analysis(
            readsdf = readscg,
            region = "L1HS_intactness_req_ALL",
            mod_code_var = "m",
            region_of_interest_from_start = 328,
            required_fraction_of_total_cg = 0.9,
            fraction_meth_threshold = 0.25,
            context = "CpG"
        )

        # read_analysis(reads, "L1HS_intactness_req_ALL", "m", "NoContext")
        # read_analysis(readscg, "L1HS_intactness_req_ALL", "h", "CpG")
        # read_analysis(reads, "L1HS_intactness_req_ALL", "h", "NoContext")
        # read_analysis(reads, "L1HS_intactness_req_ALL", "a", "NoContext")
    },
    error = function(e) {
        print(e)
    }
)


if (conf$single_condition == "no") {
    ######### GENES


    {
        directions <- c("Hypo", "Hyper", "Dif")
        mydir <- sprintf("ldna/results/%s/plots/great", params$mod_code)
        mydirtables <- sprintf("ldna/results/%s/tables/great", params$mod_code)
        dir.create(mydir, recursive = TRUE, showWarnings = FALSE)
        dir.create(mydirtables, recursive = TRUE, showWarnings = FALSE)

        refseq_gr <- import(conf$refseq_unaltered)
        genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
        genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
        genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
        mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
        mcols(genes_gr) %>% colnames()
        mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]

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

    library(clusterProfiler)
    library(scales)
    library(msigdbr)


    {
        gs <- msigdbr("human")
        get_gs_enrichments <- function(gs, gs_ontology_level, outputdir, dmrsgr, et, et_mode_string) {
            mydir <- file.path(outputdir, gs_ontology_level, et_mode_string)
            print(mydir)
            tablesMsigdb <- list()
            genecollections <- gs %>%
                pluck(gs_ontology_level) %>%
                unique()
            for (collection in genecollections) {
                print(collection)
                tryCatch(
                    {
                        dir.create(paste(mydir, collection, sep = "/"), recursive = TRUE)
                        dir.create(paste(mydirtables, collection, sep = "/"), recursive = TRUE)
                        genesets <- gs %>%
                            filter(!!sym(gs_ontology_level) == collection) %>%
                            dplyr::select(gs_name, gene_symbol) %>%
                            dplyr::rename(term = gs_name, gene = gene_symbol) %>%
                            as.data.frame()
                        for (direction in directions) {
                            if (direction == "Hypo") {
                                regions <- dmrsgr[grepl("Hypo", dmrsgr$direction)]
                            } else if (direction == "Hyper") {
                                regions <- dmrsgr[grepl("Hyper", dmrsgr$direction)]
                            } else {
                                regions <- dmrsgr
                            }


                            res <- great(regions, gene_sets = genesets, extended_tss = et, background = CHROMOSOMESINCLUDEDINANALYSIS_REF)
                            tb <- getEnrichmentTable(res)
                            if (nrow(tb) != 0) {
                                tb <- tb %>% dplyr::arrange(p_adjust)
                                tablesMsigdb[[collection]][[direction]] <- tb
                                write_delim(tb, paste(mydirtables, collection, paste0(direction, "great_enrichment.tsv"), sep = "/"))

                                png(paste(mydir, collection, paste0(direction, "volcano.png"), sep = "/"), height = 5, width = 5, res = 300, units = "in")
                                plotVolcano(res)
                                dev.off()

                                png(paste(mydir, collection, paste0(direction, "associations.png"), sep = "/"), height = 5, width = 10, res = 300, units = "in")
                                plotRegionGeneAssociations(res)
                                dev.off()


                                tbnames <- tb %>%
                                    tibble() %>%
                                    mutate(id_nchar = nchar(id)) %>%
                                    mutate(id = case_when(
                                        id_nchar < 40 ~ paste0(strrep("-", pmax(0, 40 - id_nchar)), id),
                                        TRUE ~ id
                                    ))
                                p <- tbnames %>%
                                    head(n = 10) %>%
                                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                                    mutate(id = fct_reorder(id, fold_enrichment)) %>%
                                    ggplot(aes(x = id, y = fold_enrichment)) +
                                    geom_col(aes(fill = p_adjust), color = "black") +
                                    coord_flip() +
                                    scale_color_continuous(trans = "reverse") +
                                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                                    mtclosed +
                                    anchorbar
                                mysaveandstore(pl = p, fn = paste(mydir, collection, paste0(direction, "lollipop.pdf"), sep = "/"), 6.5, 6)
                                p <- tbnames %>%
                                    head(n = 10) %>%
                                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                                    mutate(id = fct_reorder(id, fold_enrichment)) %>%
                                    ggplot(aes(x = id, y = fold_enrichment)) +
                                    geom_col(fill = "skyblue", color = "black") +
                                    coord_flip() +
                                    scale_color_continuous(trans = "reverse") +
                                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                                    mtclosed +
                                    anchorbar
                                mysaveandstore(pl = p, fn = paste(mydir, collection, paste0(direction, "lollipop_no_color.pdf"), sep = "/"), 6, 6)

                                p <- tbnames %>%
                                    head(n = 5) %>%
                                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                                    mutate(id = fct_reorder(id, fold_enrichment)) %>%
                                    ggplot(aes(x = id, y = fold_enrichment)) +
                                    geom_col(aes(fill = p_adjust), color = "black") +
                                    coord_flip() +
                                    scale_color_continuous(trans = "reverse") +
                                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                                    mtclosed +
                                    anchorbar +
                                    theme(axis.text.y = element_text(family = "mono"))

                                mysaveandstore(pl = p, fn = paste(mydir, collection, paste0(direction, "lollipop5.pdf"), sep = "/"), 6.5, 4)
                                p <- tbnames %>%
                                    head(n = 5) %>%
                                    mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                                    mutate(id = fct_reorder(id, fold_enrichment)) %>%
                                    ggplot(aes(x = id, y = fold_enrichment)) +
                                    geom_col(fill = "skyblue", color = "black") +
                                    coord_flip() +
                                    scale_color_continuous(trans = "reverse") +
                                    labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                                    mtclosed +
                                    anchorbar +
                                    theme(axis.text.y = element_text(family = "mono"))
                                mysaveandstore(pl = p, fn = paste(mydir, collection, paste0(direction, "lollipop_no_color5.pdf"), sep = "/"), 6, 4)
                            }
                        }
                    },
                    error = function(e) {}
                )
            }
            save(file = sprintf("ldna/Rintermediates/%s/tablesMsigdb_%s_%s.rds", params$mod_code, gs_ontology_level, et_mode_string), tablesMsigdb)
            return(tablesMsigdb)
        }

        for (dmrtype in dmrs$dmr_type %>% unique()) {
            tryCatch(
                {
                    get_gs_enrichments(gs, "gs_cat", sprintf("%s/%s", mydir, dmrtype), dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype], et_noextension, "et_noextension")
                },
                error = function(e) {}
            )

            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcat",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        dmrsgr = dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype],
                        et = et_noextension,
                        et_mode_string = "et_noextension"
                    )
                },
                error = function(e) {}
            )

            # get_gs_enrichments(gs, "gs_cat", sprintf("%s/%s", mydir, dmrtype), dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype], et, "et")
            # get_gs_enrichments(gs, "gs_subcat", sprintf("%s/%s", mydir, dmrtype), dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype], et, "et")
        }

        {
            promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
            write_delim(tibble(as.data.frame(promoters)) %>% mutate(score = 1000) %>% dplyr::select(seqnames, start, end, gene_id, score, strand), sprintf("ldna/Rintermediates/%s/promoters.bed", params$mod_code), col_names = FALSE, delim = "\t")
            threshold_dfs <- list()
            for (dmrtype in dmrs$dmr_type %>% unique()) {
                dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
                hyporegions <- dmrsgr_temp[grepl("Hypo", dmrsgr_temp$direction)]
                hyperregions <- dmrsgr_temp[grepl("Hyper", dmrsgr_temp$direction)]

                write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, dmrsgr_temp, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
                write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyporegions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hypo_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
                write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyperregions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hyper_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")

                hypo <- mcols(subsetByOverlaps(promoters, hyporegions))$gene_id
                hyper <- mcols(subsetByOverlaps(promoters, hyperregions))$gene_id
                disc <- intersect(hypo, hyper)
                hypo_nd <- setdiff(hypo, disc)
                hyper_nd <- setdiff(hyper, disc)



                threshold_df <- bind_rows(
                    tibble(gene_id = hypo_nd, !!sym(dmrtype) := "Hypo"),
                    tibble(gene_id = hyper_nd, !!sym(dmrtype) := "Hyper"),
                    tibble(gene_id = disc, !!sym(dmrtype) := "Discordant")
                )
                threshold_dfs[[dmrtype]] <- threshold_df
            }
            thresholddf <- purrr::reduce(threshold_dfs, full_join) %>% dplyr::select(-"t05CG10")
            promoters_df <- full_join(tibble(as.data.frame(promoters)), thresholddf) %>%
                filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
                pivot_longer(cols = c("t05", "t01"), names_to = "dmr_type", values_to = "direction")
            p <- promoters_df %>%
                group_by(dmr_type, direction) %>%
                summarise(n = n()) %>%
                ungroup() %>%
                tidyr::complete(dmr_type, direction, fill = list(n = 0)) %>%
                filter(!is.na(direction)) %>%
                mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
                ggplot() +
                geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
                geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
                ggtitle("Promoter Methylation") +
                labs(x = "Direction", y = "count") +
                theme(legend.position = "none") +
                annotation_logticks(sides = "l") +
                scale_y_log10(
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))
                ) +
                mtopen +
                scale_methylation_thresholds
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/genes_concordance.pdf", params$mod_code), 5, 4)

            for (dmrtype in dmrs$dmr_type %>% unique()) {
                library(clusterProfiler)
                hypo_genes <- mergeddf %>% filter(direction == "Hypo") %$% promoters.gene_id
                hyper_genes <- mergeddf %>% filter(direction == "Hyper") %$% promoters.gene_id
                background <- mcols(promoters)$gene_id

                gs <- msigdbr("human")
                tablesORA_subcollection <- list()
                results_ORA_hypo <- list()
                results_ORA_hyper <- list()
                genesubcollections <- gs$gs_subcat %>% unique()
                for (collection in genesubcollections) {
                    term2gene <- gs %>%
                        filter(gs_subcat == collection) %>%
                        dplyr::rename(term = gs_name, gene = gene_symbol) %>%
                        select(term, gene)
                    res_hypo <- enricher(hypo_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

                    results_ORA_hypo[[collection]] <- res_hypo %>% as.data.frame()

                    res_hyper <- enricher(hyper_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

                    results_ORA_hyper[[collection]] <- res_hyper %>% as.data.frame()
                }



                ## DMR intersection
                {
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
                    annotations_of_interest <- list(chromHMM = chromHMMgr, cCREs = ccresgr)
                    # note I should be testing for enrichment of these sets too
                    for (annotation_of_interest in names(annotations_of_interest)) {
                        annot <- annotations_of_interest[[annotation_of_interest]]
                        annotdf <- as.data.frame(annot) %>% tibble()
                        mbo <- mergeByOverlaps(annot, dmrsgr)
                        mbodf <- tibble(as.data.frame(mbo))
                        if (annotation_of_interest == "chromHMM") {
                            mbodf$annot.name <- factor(mbodf$annot.name,
                                levels =
                                    c(
                                        "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "TssBiv",
                                        "Tx", "TxWk",
                                        "EnhG1", "EnhG2",
                                        "EnhA1", "EnhA2",
                                        "EnhWk",
                                        "EnhBiv",
                                        "Het",
                                        "ZNF/Rpts",
                                        "ReprPC",
                                        "ReprPCWk",
                                        "Quies"
                                    )
                            )
                        }
                        if (annotation_of_interest == "cCREs") {
                            mbodf$annot.name <- factor(mbodf$annot.name,
                                levels =
                                    c(
                                        "PLS,CTCF-bound", "PLS",
                                        "pELS,CTCF-bound", "pELS",
                                        "dELS,CTCF-bound",
                                        "dELS",
                                        "DNase-H3K4me3,CTCF-bound",
                                        "CTCF-only,CTCF-bound",
                                        "DNase-H3K4me3",
                                        "DNase-only"
                                    )
                            )
                        }

                        p <- mbodf %>%
                            group_by(annot.name, dmr_type, direction) %>%
                            summarise(n = n()) %>%
                            ggplot() +
                            geom_col(aes(x = annot.name, y = n, fill = direction), position = "dodge", color = "black") +
                            labs(x = "") +
                            coord_flip() +
                            ggtitle("cCRE Methylation") +
                            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                            mtopen +
                            scale_methylation
                        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/%s/dmrs_in_%s.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                        library(scales)
                        p <- mbodf %>%
                            filter(dmr_type != "t05CG10") %>%
                            group_by(annot.name, direction, dmr_type) %>%
                            summarize(n = n()) %>%
                            ungroup() %>% # Summarize with groups dropped for completeness
                            tidyr::complete(annot.name, direction, dmr_type, fill = list(n = 0)) %>%
                            mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
                            ggplot() +
                            geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = annot.name, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                            geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = annot.name, y = n, group = direction, fill = direction_threshold), position = position_dodge(preserve = "single"), color = "black") +
                            labs(x = "") +
                            coord_flip() +
                            ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                            annotation_logticks(sides = "b") +
                            scale_y_log10(
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))
                            ) +
                            mtopen +
                            scale_methylation_thresholds
                        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/%s/dmrs_in_%s_log1.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                        total <- annotdf %>%
                            group_by(name) %>%
                            summarize(n = n())
                        totaldm <- mbodf %>%
                            group_by(annot.name, direction, dmr_type) %>%
                            summarize(n = n()) %>%
                            ungroup() %>% # Summarize with groups dropped for completeness
                            tidyr::complete(annot.name, direction, dmr_type, fill = list(n = 0)) %>%
                            mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type)))

                        pctdm <- left_join(totaldm, total, by = c("annot.name" = "name")) %>%
                            mutate(pct = 100 * n.x / n.y) %>%
                            mutate(myaxis = paste0(annot.name, "\n", "n=", n.y)) %>%
                            drop_na()

                        if (annotation_of_interest == "chromHMM") {
                            pctdm$annot.name <- factor(pctdm$annot.name,
                                levels =
                                    c(
                                        "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "TssBiv",
                                        "Tx", "TxWk",
                                        "EnhG1", "EnhG2",
                                        "EnhA1", "EnhA2",
                                        "EnhWk",
                                        "EnhBiv",
                                        "Het",
                                        "ZNF/Rpts",
                                        "ReprPC",
                                        "ReprPCWk",
                                        "Quies"
                                    )
                            )
                        }
                        if (annotation_of_interest == "cCREs") {
                            pctdm$annot.name <- factor(pctdm$annot.name,
                                levels =
                                    c(
                                        "PLS,CTCF-bound", "PLS",
                                        "pELS,CTCF-bound", "pELS",
                                        "dELS,CTCF-bound",
                                        "dELS",
                                        "DNase-H3K4me3,CTCF-bound",
                                        "CTCF-only,CTCF-bound",
                                        "DNase-H3K4me3",
                                        "DNase-only"
                                    )
                            )
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
                        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/%s/dmrs_in_%s_pct.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 6, 4)

                        p <- pctdm %>% ggplot() +
                            geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = annot.name, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                            geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = annot.name, y = pct, group = direction, fill = direction_threshold), position = position_dodge(), color = "black") +
                            labs(x = "", y = "Pct Differentially Methylated") +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                            ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
                            coord_flip() +
                            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                            mtopen +
                            scale_methylation_thresholds
                        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/%s/dmrs_in_%s_pct_clean.pdf", params$mod_code, annotation_of_interest, annotation_of_interest), 7, 4)
                    }
                }
            }
        }
    }
}



# {#analysis specific  look at genes

# dmls
# dmrsgr %>% subsetByOverlaps(
# genes_gr[mcols(genes_gr)$gene_id == "RMI2"] %>% resize(width = 1000)
# )

# }



# gs %>%
#     filter(gs_cat == "C5") %>%
#     filter(grepl("_ION$", gs_name, perl = TRUE)) %$% gs_name %>%
#     unique()

# metal_terms <- c(
#     "GOBP_DETOXIFICATION_OF_COPPER_ION",
#     "GOBP_STRESS_RESPONSE_TO_METAL_ION",
#     "GOBP_CELLULAR_RESPONSE_TO_ZINC_ION",
#     "GOBP_CELLULAR_RESPONSE_TO_CADMIUM_ION",
#     "GOBP_CELLULAR_RESPONSE_TO_COPPER_ION",
#     "GOBP_ZINC_ION_HOMEOSTASIS"
# )
# gs %>% filter(gs_name %in% metal_terms)
# gs %>% filter(gs_name %in% metal_terms[1])

# gs %>% filter(gs_name %in% metal_terms[2])

# gs %>% filter(gs_name %in% metal_terms[3])
# intersect(
#     gs %>% filter(gs_name %in% metal_terms[1]) %$% gene_symbol,
#     gs %>% filter(gs_name %in% metal_terms[2]) %$% gene_symbol
# )
# intersect(
#     gs %>% filter(gs_name %in% metal_terms[1]) %$% gene_symbol,
#     gs %>% filter(gs_name %in% metal_terms[6]) %$% gene_symbol
# )



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
grsdfs <- grsdf %>%
    group_by(sample, seqnames, islandStatus) %>%
    slice_sample(n = 1000)
grss <- GRanges(grsdfs)


############
# GLOBAL
dir.create(sprintf("ldna/results/%s/plots/genomewide", params$mod_code), showWarnings = FALSE)
pf <- grsdf %>%
    group_by(sample) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(sample_table %>% dplyr::rename(sample = sample_name))

p <- pf %>%
    mutate(sample = fct_reorder(paste0(sample, "_", age), age)) %>%
    ggplot(aes(y = sample, x = mean_meth, color = condition, shape = sex)) +
    geom_point(size = 3) +
    scale_conditions +
    geom_vline(xintercept = pf %>% filter(condition == "AD") %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
    geom_text_repel(aes(label = apoe)) +
    # new_scale_fill() +
    # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_conditions +
    mtopen
library(broom)
stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_withage_ordered_genomewide.pdf", params$mod_code), 5, 4, sf = stats)

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
    pivot_wider(names_from = condition, values_from = mean_meth) %>%
    mutate(dif = AD - CTRL) %>%
    group_by(dmr_type, direction) %>%
    summarize(mean_dif = mean(dif))
dir.create(sprintf("ldna/results/%s/tables/dmrs", params$mod_code), recursive = TRUE)
ad %>% write_delim(sprintf("ldna/results/%s/tables/dmrs/dmrs_meth_dif", params$mod_code))


topdmrshypo <- dmrsgrsuuid %>%
    as.data.frame() %>%
    tibble() %>%
    filter(direction == "Hypo") %>%
    arrange(areaStat) %$% dmrid %>%
    head(n = 500)
topdmrshyper <- dmrsgrsuuid %>%
    as.data.frame() %>%
    tibble() %>%
    filter(direction == "Hyper") %>%
    arrange(-areaStat) %$% dmrid %>%
    head(n = 500)

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
    head(n = 500) %>%
    pivot_longer(cols = sample_table$sample_name, names_to = "sample_name", values_to = "mean_meth") %>%
    left_join(sample_table)

topdmrshyper <- pf1 %>%
    as.data.frame() %>%
    tibble() %>%
    filter(direction == "Hyper") %>%
    arrange(-areaStat) %>%
    pivot_wider(names_from = sample, values_from = mean_meth) %>%
    head(n = 500) %>%
    pivot_longer(cols = sample_table$sample_name, names_to = "sample_name", values_to = "mean_meth") %>%
    left_join(sample_table)

pf <- bind_rows(topdmrshyper, topdmrshypo) %>% dplyr::rename(sample = sample_name)

library(tidyHeatmap)
p <- pf %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    group_by(direction) %>%
    heatmap(dmrid, sample, mean_meth,
        cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_row_dend = FALSE # palette_value = circlize::colorRamp2(
        # seq(0, 100, length.out = 11),
        # RColorBrewer::brewer.pal(11, "RdBu")
        # )
    )

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
    )
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap_colclust.pdf", params$mod_code), w = 4.5, h = 9, res = 300, pl = p, raster = FALSE)
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/dmrheatmap_colclust.pdf", params$mod_code), w = 4.5, h = 9, res = 300, pl = p, raster = TRUE)


pf <- dmrs_meth_df %>%
    group_by(sample, dmrid, condition, direction, dmr_type, nCG, length) %>%
    summarise(mean_meth = mean(pctM)) %>%
    ungroup()

p <- pf %>%
    filter(condition == "CTRL") %>%
    group_by(direction, dmrid, dmr_type) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
    mutate(direction_threshold = factor(direction_threshold, levels = c("Hypo_05", "Hypo_01", "Hyper_05", "Hyper_01"))) %>%
    ggplot(aes(x = mean_meth, fill = direction_threshold)) +
    xlab("Average DMR Meth (CTRL)") +
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
    xlab("Average DMR Meth (CTRL)") +
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
        as_ComplexHeatmap()
    hms[[sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 6)
}
# Generate the expression as a string and parse it
p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
mysaveandstore(sprintf("%s/%s_methylation_%s1.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
rm(p)



# nonref analysis



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
