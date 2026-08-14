module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = sample_name)
samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")

cell_fractions <- read_delim("ldna/results/m/tables/scMD_cell_type_fractions_all.csv")
cell_fractions_scaled <- cell_fractions %>%
    mutate(Neuron = Inh + Exc) %>%
    mutate(across(all_of(c(cell_types_cols, "Neuron")), ~ as.numeric(scale(.)), .names = "{.col}_z"))

sample_table <- sample_table %>% left_join(cell_fractions_scaled)

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
library(glmmTMB)
library(broom.mixed)

conditions <- conf$levels
contrasts <- conf$contrasts

# Helper to parse a contrast string into condition1 (reference) and condition2 (test)
parse_contrast <- function(contrast_string) {
    parts <- str_match(contrast_string, "^condition_(.+)_vs_(.+)$")
    list(condition2 = parts[1, 2], condition1 = parts[1, 3])
}

# For backwards compatibility, set condition1/condition2 from first contrast
first_contrast <- parse_contrast(contrasts[[1]])
condition1 <- first_contrast$condition1
condition2 <- first_contrast$condition2
condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
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
    nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
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
            contrasts = conf$contrasts,
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

# rmannextended <- get_repeat_annotations(
#     default_or_extended = "default",
#     keep_non_central = FALSE
# )

# rmannextended %>% filter(rte_subfamily == "L1HS") %>% filter(refstatus == "Ref") %>% filter(intactness_req == "Intact")
rmannextended <- get_repeat_annotations(
    default_or_extended = "extended",
    keep_non_central = FALSE
)

flRTEpromoter <- read_delim(sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)
RMdf <- read_delim(sprintf("ldna/Rintermediates/%s/RMdf.tsv", params$mod_code), col_names = TRUE)




##################


cpg_islands <- rtracklayer::import(conf$cpg_islands)
cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)


if ((conf$single_condition == "no")) {
    # Read per-contrast DMR/DML data
    dmrs_per_contrast <- list()
    dmls_per_contrast <- list()
    dmrsgr_per_contrast <- list()
    dmlsgr_per_contrast <- list()
    dmrsannot_per_contrast <- list()
    dmrsgr_split_per_contrast <- list()
    for (contrast in contrasts) {
        dmr_path <- sprintf("ldna/results/%s/tables/%s/dmrs.tsv", params$mod_code, contrast)
        dml_path <- sprintf("ldna/results/%s/tables/%s/dmls.tsv", params$mod_code, contrast)
        if (!file.exists(dmr_path) || !file.exists(dml_path)) {
            message(sprintf("Skipping contrast %s: DMR/DML files not found", contrast))
            next
        }
        dmrs_per_contrast[[contrast]] <- read_delim(dmr_path, delim = "\t", col_names = TRUE) %>% filter(dmr_type %in% c("t01", "t05"))
        dmls_per_contrast[[contrast]] <- read_delim(dml_path, delim = "\t", col_names = TRUE) %>% filter(fdrs <= 0.2)

        dmrsgr_per_contrast[[contrast]] <- GRanges(dmrs_per_contrast[[contrast]])
        dmlsgr_per_contrast[[contrast]] <- GRanges(
            seqnames = dmls_per_contrast[[contrast]]$chr,
            ranges = IRanges(start = dmls_per_contrast[[contrast]]$pos, end = dmls_per_contrast[[contrast]]$pos),
            stat = dmls_per_contrast[[contrast]]$stat,
            pval = dmls_per_contrast[[contrast]]$pvals,
            fdr = dmls_per_contrast[[contrast]]$fdrs,
            direction = dmls_per_contrast[[contrast]]$direction
        )
        dmrsannot_per_contrast[[contrast]] <- dmrs_per_contrast[[contrast]] %>%
            mutate(direction_threshold = paste(direction, gsub("t", "", dmr_type), sep = "_")) %>%
            GRanges()
        dmrsgr_split_per_contrast[[contrast]] <- split(dmrsannot_per_contrast[[contrast]], dmrsannot_per_contrast[[contrast]]$direction_threshold)
    }
}



refseq_gr <- import(conf$refseq_unaltered)
genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
mcols(genes_gr) %>% colnames()
mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)


ccresdf <- read_delim(conf$ccres, col_names = FALSE)
ccresgr <- GRanges(
    seqnames = ccresdf$X1,
    ranges = IRanges(start = ccresdf$X2, end = ccresdf$X3),
    name = ccresdf$X10,
    UID = ccresdf$X4,
)

chromHMMgr <- import(conf$chromHMM)
chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]
prom_no_mcols <- promoters
mcols(prom_no_mcols) <- NULL
enh_no_mcols <- chromHMM_enhancers_grs
mcols(enh_no_mcols) <- NULL


###########################


##################################### DML / DMR analysis
if ((conf$single_condition == "no")) {
    for (contrast in contrasts) {
        if (is.null(dmrs_per_contrast[[contrast]])) next
        cp <- parse_contrast(contrast)
        condition1 <- cp$condition1
        condition2 <- cp$condition2
        condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
        condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
        dmrs <- dmrs_per_contrast[[contrast]]
        dmls <- dmls_per_contrast[[contrast]]
        dmrsgr <- dmrsgr_per_contrast[[contrast]]
        dmlsgr <- dmlsgr_per_contrast[[contrast]]
        dmrsgr %>%
            as.data.frame() %>%
            tibble() %>%
            arrange(-areaStat)
        dmrtypes <- dmrs$dmr_type %>% unique()
        dmr_grs_cpg_islands <- dmrsgr %>% subsetByOverlaps(cpg_islands)
        mcols(dmr_grs_cpg_islands)$islandStatus <- rep("island", length(dmr_grs_cpg_islands))
        dmr_grs_cpgi_shores <- dmrsgr %>% subsetByOverlaps(cpgi_shores)
        dmr_grs_cpgi_shores_filtered <- dmr_grs_cpgi_shores %>% subsetByOverlaps(dmr_grs_cpg_islands, invert = TRUE)
        mcols(dmr_grs_cpgi_shores_filtered)$islandStatus <- rep("shore", length(dmr_grs_cpgi_shores_filtered))
        dmr_grs_cpgi_shelves <- dmrsgr %>% subsetByOverlaps(cpgi_shelves)
        dmr_grs_cpgi_shelves_filtered <- dmr_grs_cpgi_shelves %>% subsetByOverlaps(dmr_grs_cpgi_shores, invert = TRUE)
        mcols(dmr_grs_cpgi_shelves_filtered)$islandStatus <- rep("shelf", length(dmr_grs_cpgi_shelves_filtered))
        dmr_grs_cpg_opensea <- dmrsgr %>% subsetByOverlaps(cpgi_features, invert = TRUE)
        mcols(dmr_grs_cpg_opensea)$islandStatus <- rep("opensea", length(dmr_grs_cpg_opensea))
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
                ggtitle(sprintf("DMR Counts (%s)", contrast)) +
                mtopen +
                scale_methylation
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_number.pdf", params$mod_code, contrast, dmrtype), 4, 4)

            p <- ggplot(data = dmrs_temp) +
                geom_histogram(aes(length), fill = mycolor, color = "black") +
                ggtitle(sprintf("DMR Lengths (%s)", contrast)) +
                labs(x = "length (bp)") +
                xlim(0, 3000) +
                mtopen +
                anchorbar
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_length.pdf", params$mod_code, contrast, dmrtype), w = 4, h = 4)

            p <- ggplot(data = dmrs_temp) +
                geom_histogram(aes(length, fill = direction), alpha = 0.7, color = "black") +
                ggtitle(sprintf("DMR Lengths (%s)", contrast)) +
                labs(x = "length (bp)") +
                xlim(0, 3000) +
                mtopen +
                scale_methylation +
                anchorbar
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_length_stratified.pdf", params$mod_code, contrast, dmrtype), w = 4, h = 4)

            p <- dmrsgrislandStatusdf %>%
                group_by(islandStatus, direction) %>%
                summarize(n = n()) %>%
                ggplot() +
                geom_col(aes(x = islandStatus, y = n, fill = direction), position = "dodge", color = "black") +
                mtopen +
                scale_methylation +
                anchorbar +
                labs(x = "", y = "Count") +
                ggtitle(sprintf("DMR (%s)", contrast)) +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_count_islandstatus.pdf", params$mod_code, contrast, dmrtype), 5, 4)


            dmrlocdf <- dmrs_temp %>%
                group_by(seqnames, direction) %>%
                summarize(n = n())
            dmrlocdf$seqnames <- factor(dmrlocdf$seqnames, levels = chromosomes)
            p <- ggplot(data = dmrlocdf) +
                geom_col(aes(y = seqnames, x = n, fill = direction), position = "dodge", color = "black") +
                theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
                labs(x = "count", y = "") +
                ggtitle(sprintf("DMR Location (%s)", contrast)) +
                mtopen +
                scale_methylation
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_count.pdf", params$mod_code, contrast, dmrtype), 5, 5)
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
            ggtitle(sprintf("DMR (%s)", contrast)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_count.pdf", params$mod_code, contrast, "both"), 5, 4)

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
            ggtitle(sprintf("DMR (%s)", contrast))
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_count_numannot.pdf", params$mod_code, contrast, "both"), 3.75, 4)

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
            ggtitle(sprintf("DMR (%s)", contrast)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/%s/dmr_count_islandstatus.pdf", params$mod_code, contrast, "both"), 5, 4)

        p <- dmls %>%
            ggplot() +
            geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
            ggtitle(sprintf("DML Counts (%s)", contrast)) +
            labs(x = "", y = "Count") +
            anchorbar +
            mtclosed +
            scale_methylation
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dml_count.pdf", params$mod_code, contrast), 4, 4)

        dmllocdf <- dmls %>%
            group_by(chr, direction) %>%
            summarize(n = n())
        dmllocdf$chr <- factor(dmllocdf$chr, levels = chromosomes)

        dml_grs_cpg_islands <- dmlsgr %>% subsetByOverlaps(cpg_islands)
        mcols(dml_grs_cpg_islands)$islandStatus <- rep("island", length(dml_grs_cpg_islands))
        dml_grs_cpgi_shores <- dmlsgr %>% subsetByOverlaps(cpgi_shores)
        dml_grs_cpgi_shores_filtered <- dml_grs_cpgi_shores %>% subsetByOverlaps(dml_grs_cpg_islands, invert = TRUE)
        mcols(dml_grs_cpgi_shores_filtered)$islandStatus <- rep("shore", length(dml_grs_cpgi_shores_filtered))
        dml_grs_cpgi_shelves <- dmlsgr %>% subsetByOverlaps(cpgi_shelves)
        dml_grs_cpgi_shelves_filtered <- dml_grs_cpgi_shelves %>% subsetByOverlaps(dml_grs_cpgi_shores, invert = TRUE)
        mcols(dml_grs_cpgi_shelves_filtered)$islandStatus <- rep("shelf", length(dml_grs_cpgi_shelves_filtered))
        dml_grs_cpg_opensea <- dmlsgr %>% subsetByOverlaps(cpgi_features, invert = TRUE)
        mcols(dml_grs_cpg_opensea)$islandStatus <- rep("opensea", length(dml_grs_cpg_opensea))
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
            ggtitle(sprintf("DML Island Status (%s)", contrast))
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dml_count_islandstatus.pdf", params$mod_code, contrast), 4, 4)

        # RTE promoter DM analysis for this contrast
        dir.create(sprintf("ldna/results/%s/plots/figs/%s", params$mod_code, contrast), recursive = TRUE, showWarnings = FALSE)

        # Get the contrast-specific DMR columns from flRTEpromoter
        contrast_dmrtype_cols <- grep(paste0("_", contrast, "$"), colnames(flRTEpromoter), value = TRUE)
        if (length(contrast_dmrtype_cols) == 0) next
        # Create a temp version with simple column names for this contrast
        dmrtypes_simple <- gsub(paste0("_", contrast), "", contrast_dmrtype_cols)
        flRTEpromoter_c <- flRTEpromoter
        for (i in seq_along(contrast_dmrtype_cols)) {
            flRTEpromoter_c <- flRTEpromoter_c %>% dplyr::rename(!!sym(dmrtypes_simple[i]) := !!sym(contrast_dmrtype_cols[i]))
        }
        # Drop other contrast columns
        other_contrast_cols <- grep("^t0[0-9]", colnames(flRTEpromoter_c), value = TRUE)
        other_contrast_cols <- other_contrast_cols[!(other_contrast_cols %in% dmrtypes_simple)]
        if (length(other_contrast_cols) > 0) {
            flRTEpromoter_c <- flRTEpromoter_c %>% dplyr::select(-all_of(other_contrast_cols))
        }

        dmrtypes <- dmrtypes_simple[dmrtypes_simple %in% c("t05", "t01")]
        flRTEpromoterlong <- pivot_longer(data = flRTEpromoter_c %>% dplyr::select(-any_of(c("t001", "t05CG10"))), cols = any_of(dmrtypes), names_to = "dmr_type", values_to = "direction")

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
            pw() %$% direction
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
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter.pdf", params$mod_code, contrast, "all", "all"), 12, 4)

        p <- pff %>%
            left_join(numdf) %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            filter(rte_subfamily != "Other") %>%
            ggplot() +
            geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
            geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
            mtclosed +
            scale_methylation_thresholds
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter_noOther.pdf", params$mod_code, contrast, "all", "all"), 12, 4)

        p <- pff %>%
            left_join(numdf) %>%
            mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
            filter(grepl("^L1", rte_subfamily)) %>%
            ggplot() +
            geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
            geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = ann_axis, y = frac_dm, group = direction, fill = direction_threshold), position = "dodge", color = "black") +
            labs(x = "", y = "Fraction Differentially Methylated") +
            ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
            mtclosed +
            scale_methylation_thresholds
        mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter_l1s.pdf", params$mod_code, contrast, "all", "all"), 10, 4)


        for (dmrtype in dmrtypes) {
            pff <- flRTEpromoter_c %>%
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
            numdf <- flRTEpromoter_c %>%
                group_by(rte_subfamily, genic_loc) %>%
                summarise(group_n_accurate = n())
            p <- pff %>%
                left_join(numdf) %>%
                mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
                ggplot() +
                geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
                facet_wrap(~genic_loc, scales = "free_x", nrow = 2) +
                labs(x = "", y = "Fraction Differentially Methylated") +
                ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
                mtclosed +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter_regionstrat.pdf", params$mod_code, contrast, dmrtype, "all"), 12, 7)

            pff <- flRTEpromoter_c %>%
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
            numdf <- flRTEpromoter_c %>%
                group_by(rte_subfamily) %>%
                summarise(group_n_accurate = n())
            p <- pff %>%
                left_join(numdf) %>%
                mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
                ggplot() +
                geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
                labs(x = "", y = "Fraction Differentially Methylated") +
                ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
                mtclosed +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter.pdf", params$mod_code, contrast, dmrtype, "all"), 12, 4)

            p <- pff %>%
                left_join(numdf) %>%
                filter(rte_subfamily != "Other") %>%
                mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
                ggplot() +
                geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
                labs(x = "", y = "Fraction Differentially Methylated") +
                ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
                mtclosed +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter1.pdf", params$mod_code, contrast, dmrtype, "all"), 10, 4)

            p <- pff %>%
                left_join(numdf) %>%
                filter(grepl("^L1", rte_subfamily)) %>%
                mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
                ggplot() +
                geom_col(aes(x = ann_axis, y = frac_dm, fill = !!sym(dmrtype)), position = "dodge", color = "black") +
                labs(x = "", y = "Fraction Differentially Methylated") +
                ggtitle(sprintf("Full Length %s Promoter Differential Methylation (%s)", "RTE", contrast)) +
                mtclosed +
                scale_methylation
            mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/dm_%s_fl%s_promoter2.pdf", params$mod_code, contrast, dmrtype, "all"), 8, 4)
        }
    } # end contrast loop
} # end single_condition check


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

    for (contrast in contrasts) {
    if (is.null(dmrsgr_split_per_contrast[[contrast]])) next
    dmrsgr_split <- dmrsgr_split_per_contrast[[contrast]]
    dmrsgr <- dmrsgr_per_contrast[[contrast]]
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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet.pdf", params$mod_code, contrast, annotation_of_interest, "bidirectional_all", annotation_of_interest), 8, 6)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet_05.pdf", params$mod_code, contrast, annotation_of_interest, "bidirectional_all", annotation_of_interest), 8, 6)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet_05_together.pdf", params$mod_code, contrast, annotation_of_interest, "bidirectional_all", annotation_of_interest), 7, 4)


        enrichdflist <- list("alloverlaps" = enrichdf, "alloverlapsrev" = enrichdfrev, "bestoverlaps" = enrichdf_bestoverlap)
        for (enrichtype in names(enrichdflist)) {
            enrichdftemp <- enrichdflist[[enrichtype]]

            p <- enrichdftemp %>%
                ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                geom_col(position = "dodge", color = "black") +
                scale_methylation_thresholds +
                mtclosedgridv
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s.pdf", params$mod_code, contrast, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

            p <- enrichdftemp %>%
                filter(direction_threshold %in% c("Hyper_05", "Hypo_05")) %>%
                ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                geom_col(position = "dodge", color = "black") +
                scale_methylation_thresholds +
                mtclosedgridv
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_05.pdf", params$mod_code, contrast, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

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
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet.pdf", params$mod_code, contrast, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

            p <- enrichdftemp %>%
                mutate(direction = gsub("_.*", "", direction_threshold)) %>%
                ggplot(aes(y = name_padded, x = log(fold_enrichment), fill = direction_threshold)) +
                geom_col(position = "dodge", color = "black") +
                geom_vline(xintercept = 0, color = "green") +
                facet_wrap(~direction) +
                scale_methylation_thresholds +
                ggtitle("Enrichments") +
                mtclosedgridv
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet_noaxislimmirror.pdf", params$mod_code, contrast, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)

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
            mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/%s/dmr_enrichments_in_%s_facet_dot.pdf", params$mod_code, contrast, annotation_of_interest, enrichtype, annotation_of_interest), 7, 4)
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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmrs_in_%s.pdf", params$mod_code, contrast, annotation_of_interest, annotation_of_interest), 6, 4)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmrs_in_%s_log.pdf", params$mod_code, contrast, annotation_of_interest, annotation_of_interest), 6, 4)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmrs_in_%s_log_numannot.pdf", params$mod_code, contrast, annotation_of_interest, annotation_of_interest), 6, 4)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmrs_in_%s_pct.pdf", params$mod_code, contrast, annotation_of_interest, annotation_of_interest), 6, 4)

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
        mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/dmr_counts_and_enrichments/%s/%s/dmrs_in_%s_pct_orderedaxis.pdf", params$mod_code, contrast, annotation_of_interest, annotation_of_interest), 6, 4)
    }
    } # end contrast loop for enrichments
}



grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
grsdf %$% sample %>% unique()
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
grs <- GRanges(grsdf)


# FINAL SECTION
library(tidyHeatmap)

for (contrast in contrasts) {
    if (is.null(dmrs_per_contrast[[contrast]])) next
    cp <- parse_contrast(contrast)
    condition1 <- cp$condition1
    condition2 <- cp$condition2
    contrast_samples <- sample_table %>% filter(condition %in% c(condition1, condition2)) %>% pull(sample_name)
    dmrs <- dmrs_per_contrast[[contrast]]
    dmrsgr <- dmrsgr_per_contrast[[contrast]]

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

    # filter to contrast samples
    dmrs_meth_df_contrast <- dmrs_meth_df %>% filter(sample %in% contrast_samples)

    # dmr av meth stats
    ac <- dmrs_meth_df_contrast %>%
        group_by(pos, direction, dmr_type, condition) %>%
        summarise(mean_meth = mean(pctM))
    ad <- ac %>%
        filter(dmr_type %in% c("t01", "t05")) %>%
        pivot_wider(names_from = condition, values_from = mean_meth) %>%
        mutate(dif = !!sym(condition2) - !!sym(condition1)) %>%
        group_by(dmr_type, direction) %>%
        summarize(mean_dif = mean(dif))
    dir.create(sprintf("ldna/results/%s/tables/dmrs/%s", params$mod_code, contrast), recursive = TRUE)
    ad %>% write_delim(sprintf("ldna/results/%s/tables/dmrs/%s/dmrs_meth_dif", params$mod_code, contrast))

    pf1 <- dmrs_meth_df_contrast %>%
        group_by(sample, dmrid, direction, dmr_type, areaStat) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup()

    topdmrshypo <- pf1 %>%
        filter(direction == "Hypo") %>%
        arrange(areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = any_of(contrast_samples), names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    topdmrshyper <- pf1 %>%
        filter(direction == "Hyper") %>%
        arrange(-areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = any_of(contrast_samples), names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    pf <- bind_rows(topdmrshyper, topdmrshypo)

    p <- pf %>%
        mutate(sample = factor(sample, levels = contrast_samples)) %>%
        group_by(direction) %>%
        heatmap(dmrid, sample, mean_meth,
            cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_row_dend = FALSE
        ) %>%
        annotation_tile(condition) %>%
        annotation_tile(sex)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap.pdf", params$mod_code, contrast), w = 5, h = 10, res = 300, pl = p, raster = FALSE)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap.pdf", params$mod_code, contrast), w = 5, h = 10, res = 300, pl = p, raster = TRUE)

    p <- pf %>%
        mutate(sample = factor(sample, levels = contrast_samples)) %>%
        group_by(direction) %>%
        heatmap(dmrid, sample, mean_meth,
            cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_row_dend = FALSE, column_km = 2
        ) %>%
        annotation_tile(condition) %>%
        annotation_tile(sex)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap_colclust.pdf", params$mod_code, contrast), w = 4.5, h = 9, res = 300, pl = p, raster = FALSE)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap_colclust.pdf", params$mod_code, contrast), w = 4.5, h = 9, res = 300, pl = p, raster = TRUE)

    pf <- dmrs_meth_df_contrast %>%
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
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_histogram.pdf", params$mod_code, contrast), w = 6, h = 4, res = 300, pl = p)

    p <- dmrsgrsuuid %>%
        as.data.frame() %>%
        tibble() %>%
        mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
        mutate(direction_threshold = factor(direction_threshold, levels = c("Hypo_05", "Hypo_01", "Hyper_05", "Hyper_01"))) %>%
        ggplot(aes(x = nCG, fill = direction_threshold)) +
        xlab("nCG per DMR") +
        geom_histogram(color = "black", bins = 10) +
        facet_wrap(~direction, scales = "free_y") +
        scale_methylation_thresholds +
        lims(x = c(0, 15)) +
        mtclosed
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_histogram_nCG.pdf", params$mod_code, contrast), w = 6.5, h = 4, res = 300, pl = p)

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
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_histogram_length.pdf", params$mod_code, contrast), w = 6, h = 4, res = 300, pl = p)

    # PCA on contrast samples only
    contrast_sample_table <- sample_table %>% filter(sample_name %in% contrast_samples)
    pcaframe <- dmrs_meth_df_contrast %>%
        group_by(sample, dmrid, condition, direction) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup() %>%
        dplyr::select(sample, mean_meth, dmrid) %>%
        pivot_wider(names_from = dmrid, values_from = mean_meth) %>%
        mutate(sample = factor(sample, levels = contrast_samples)) %>%
        arrange(sample) %>%
        column_to_rownames(var = "sample") %>%
        as.matrix() %>%
        t()

    pcaObj <- pca(pcaframe, center = TRUE, scale = FALSE, metadata = contrast_sample_table %>% column_to_rownames(var = "sample_name"))

    p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_scree.pdf", params$mod_code, contrast), w = 5, h = 10, res = 300, pl = p)

    p <- plotloadings(pcaObj,
        components = getComponents(pcaObj, seq_len(min(3, ncol(pcaframe)))),
        rangeRetain = 0.045, labSize = 2
    ) + mtopen
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_loadings.pdf", params$mod_code, contrast), w = 5, h = 10, res = 300, pl = p)

    p <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0, legendPosition = "right", shape = "sex", colby = "condition",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    ) + mtopen + scale_conditions
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmr_biplot.pdf", params$mod_code, contrast), w = 5, h = 5, res = 300, pl = p)

    # All-condition heatmap: show all samples, split by condition
    pf1_all <- dmrs_meth_df %>%
        group_by(sample, dmrid, direction, dmr_type, areaStat) %>%
        summarise(mean_meth = mean(pctM)) %>%
        ungroup()

    topdmrshypo_all <- pf1_all %>%
        filter(direction == "Hypo") %>%
        arrange(areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = any_of(as.character(sample_table$sample_name)), names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    topdmrshyper_all <- pf1_all %>%
        filter(direction == "Hyper") %>%
        arrange(-areaStat) %>%
        pivot_wider(names_from = sample, values_from = mean_meth) %>%
        head(n = 250) %>%
        pivot_longer(cols = any_of(as.character(sample_table$sample_name)), names_to = "sample_name", values_to = "mean_meth") %>%
        left_join(sample_table)

    pf_all <- bind_rows(topdmrshyper_all, topdmrshypo_all)

    p <- pf_all %>%
        mutate(sample = factor(sample, levels = conf$samples)) %>%
        mutate(condition = factor(condition, levels = conf$levels)) %>%
        group_by(direction) %>%
        heatmap(dmrid, sample, mean_meth,
            cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_row_dend = FALSE
        ) %>%
        annotation_tile(condition) %>%
        annotation_tile(sex)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap_allconditions.pdf", params$mod_code, contrast), w = 8, h = 10, res = 300, pl = p, raster = FALSE)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/%s/dmrheatmap_allconditions.pdf", params$mod_code, contrast), w = 8, h = 10, res = 300, pl = p, raster = TRUE)

} # end contrast loop for FINAL SECTION

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "14_dmrs", params$mod_code))
