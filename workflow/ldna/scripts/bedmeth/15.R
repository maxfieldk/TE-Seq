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





##################


perelementdf <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf.tsv", params$mod_code), col_names = TRUE)
perelementdf$sample <- factor(perelementdf$sample, levels = conf$samples)
perelementdf$condition <- factor(perelementdf$condition, levels = conf$levels)

perelementdf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf_promoters.tsv", params$mod_code), col_names = TRUE)
perelementdf_promoters$sample <- factor(perelementdf_promoters$sample, levels = conf$samples)
perelementdf_promoters$condition <- factor(perelementdf_promoters$condition, levels = conf$levels)


###########################



{
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
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_L1s.pdf", params$mod_code), raster = FALSE, 12, 4)
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
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS_loc.pdf", params$mod_code), raster = FALSE, 12, 4)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS_loc.pdf", params$mod_code), raster = TRUE, 12, 4)
    }

    p <- perelementdf_promoters %>%
        filter(grepl("^L1HS", rte_subfamily)) %>%
        group_by(rte_subfamily) %>%
        ggplot(aes(x = loc_superlowres_integrative_stranded, y = mean_meth, color = condition)) +
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
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc.pdf", params$mod_code), raster = FALSE, 6, 4)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc.pdf", params$mod_code), raster = TRUE, 6, 4)
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
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_condition_L1HS_loc_lowres.pdf", params$mod_code), raster = FALSE, 5, 6)
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
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/repmasker_boxplot_promoters_by_sample_L1HS.pdf", params$mod_code), raster = FALSE, 12, 4)
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
    # rmannextended %>%         filter(rte_subfamily == "L1HS") %>% filter(intactness_req == "Intact") %>% filter(refstatus == "Ref")
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
        for (contrast in contrasts) {
            cp <- parse_contrast(contrast)
            condition1 <- cp$condition1
            condition2 <- cp$condition2
            condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
            condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
            dmrs <- dmrs_per_contrast[[contrast]]

            # Get contrast-specific DMR columns and rename to simple names
            contrast_dmrtype_cols <- grep(paste0("_", contrast, "$"), colnames(perelementdf_promoters), value = TRUE)
            dmrtypes_simple <- gsub(paste0("_", contrast), "", contrast_dmrtype_cols)
            perelementdf_promoters_c <- perelementdf_promoters %>%
                filter(condition %in% c(condition1, condition2))
            for (i in seq_along(contrast_dmrtype_cols)) {
                if (contrast_dmrtype_cols[i] %in% colnames(perelementdf_promoters_c)) {
                    perelementdf_promoters_c <- perelementdf_promoters_c %>% dplyr::rename(!!sym(dmrtypes_simple[i]) := !!sym(contrast_dmrtype_cols[i]))
                }
            }
            other_contrast_cols <- grep("^t0[0-9].*_condition_", colnames(perelementdf_promoters_c), value = TRUE)
            if (length(other_contrast_cols) > 0) {
                perelementdf_promoters_c <- perelementdf_promoters_c %>% dplyr::select(-all_of(other_contrast_cols))
            }

            dmrtypes <- dmrtypes_simple[dmrtypes_simple %in% c("t05", "t01")]

            pfl1 <- perelementdf_promoters_c %>%
                filter(grepl("^L1", rte_subfamily)) %>%
                dplyr::select(-any_of(c("t05CG10", "t001")))
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
            mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/%s/repmasker_paired_promoters_L1s.pdf", params$mod_code, contrast), 14, 6, raster = TRUE)

            pfl1hs <- pfl1 %>%
                filter(rte_subfamily == "L1HS")
            pfl1hs %>% arrange(mean_meth)
            l1hs_paired_dif_frame <- pfl1 %>%
                filter(rte_subfamily == "L1HS") %>%
                pivot_longer(cols = any_of(dmrtypes), names_to = "dmr_type", values_to = "direction") %>%
                mutate(direction_threshold = ifelse(is.na(direction), "NS", paste0(direction, "_", gsub("t", "", dmr_type)))) %>%
                filter(!(dmr_type == "t01" & is.na(direction))) %>%
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
            mysaveandstore(pl = p, fn = sprintf("ldna/results/%s/plots/rte/%s/repmasker_paired_promoters_l1hs_%s.pdf", params$mod_code, contrast, "all"), 5, 4)

            for (dmrtype in dmrs$dmr_type %>% unique()) {
                if (!(dmrtype %in% dmrtypes)) next
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
                mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/%s/repmasker_paired_promoters_L1s2_%s.pdf", params$mod_code, contrast, dmrtype), 14, 6, raster = TRUE)

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
                mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/%s/repmasker_paired_promoters_l1hs_%s.pdf", params$mod_code, contrast, dmrtype), 4, 4, raster = FALSE)
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


            pf <- perelementdf_promoters_c %>%
                filter(rte_subfamily == "L1HS")
            p <- pf %>%
                ggplot() +
                geom_quasirandom(aes(x = intactness_req, y = mean_meth, color = condition), dodge.width = 0.75) +
                geom_boxplot(aes(x = intactness_req, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
                xlab("") +
                ylab("Average CpG Methylation Per Element") +
                ggtitle(sprintf("RTE CpG Methylation (%s)", contrast)) +
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
                ggtitle(sprintf("RTE CpG Methylation (%s)", contrast)) +
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
                    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/%s/l1hs_boxplot_promoters.pdf", params$mod_code, contrast), 5, 4, sf = stats)
                },
                error = function(e) {
                    mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/%s/l1hs_boxplot_promoters.pdf", params$mod_code, contrast), 5, 4)
                }
            )
        } # end contrast loop
    }

    # Reset condition1/condition2 to first contrast for remaining code
    cp <- parse_contrast(contrasts[[1]])
    condition1 <- cp$condition1
    condition2 <- cp$condition2
    condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
    condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
    dmrs <- dmrs_per_contrast[[contrasts[[1]]]]
    dmls <- dmls_per_contrast[[contrasts[[1]]]]
    dmrsgr <- dmrsgr_per_contrast[[contrasts[[1]]]]
    dmlsgr <- dmlsgr_per_contrast[[contrasts[[1]]]]
    dmrsannot <- dmrsannot_per_contrast[[contrasts[[1]]]]
    dmrsgr_split <- dmrsgr_split_per_contrast[[contrasts[[1]]]]
}

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "15", params$mod_code))
