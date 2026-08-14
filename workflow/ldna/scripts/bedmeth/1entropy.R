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



###########################




# global entropy
entropydfs <- list()
for (sample in sample_table$sample_name) {
    enttemp <- read_delim(str_glue("ldna/results/m/tables/entropy/{sample}/regions.bed"))
    enttemp$sample_name <- sample
    entropydfs[[sample]] <- enttemp
}
entdf <- entropydfs %>% purrr::reduce(bind_rows)


entdf %>%
    group_by(region_name) %>%
    summarise(n_represented = n())


pf <- entdf %>%
    group_by(region_name) %>%
    mutate(n_represented = n()) %>%
    ungroup() %>%
    filter(n_represented == 12) %>%
    left_join(sample_table) %>%
    group_by(sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    ungroup()
p <- pf %>% ggplot() +
    geom_point(aes(x = condition, y = me))
mysaveandstore("zztemp.32411111111111111111111111111111111.pdf")

summary(lm(me ~ condition + mean_num_reads + age + sex, pf))


entdf %>%
    group_by(sample_name) %>%
    p() <- entdf %>%
    group_by(sample_name) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    ggplot() +
    geom_point(aes(x = mean_num_reads, y = me)) +
    mtclosed
mysaveandstore("zztemp.3241111111111111111111111111111111.pdf")


# dmr entropy
entropydfs <- list()
for (subset in c("t05", "hyper_t05", "hypo_t05", "t01", "hyper_t01", "hypo_t01")) {
    for (sample in sample_table$sample_name) {
        enttemp <- read_delim(str_glue("ldna/results/m/tables/entropy_{subset}/{sample}/regions.bed"))
        enttemp$sample_name <- sample
        enttemp$type <- subset
        entropydfs[[paste0(sample, subset)]] <- enttemp
    }
}
entdf <- entropydfs %>% purrr::reduce(bind_rows)


entdf %>%
    group_by(region_name) %>%
    summarise(n_represented = n())


sflist <- list()
for (subset in c("t05", "hyper_t05", "hypo_t05", "t01", "hyper_t01", "hypo_t01")) {
    tempdf <- entdf %>% filter(type == subset)

    pf <- tempdf %>%
        group_by(region_name) %>%
        mutate(n_represented = n()) %>%
        ungroup() %>%
        filter(n_represented == 12) %>%
        left_join(sample_table) %>%
        group_by(sample_name, condition, sex, age) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ungroup()
    p <- pf %>%
        ggviolin(x = "condition", y = "me", fill = "condition", add = c("mean_se", "dotplot")) +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("ldna/results/m/plots/entropy_{subset}/mean_entropy.pdf"), 5, 5, pl = p)

    sftemp <- broom::tidy(summary(lm(me ~ condition + age + sex, pf)))
    sftemp$type <- subset
    sflist[[subset]] <- sftemp

    p <- entdf %>%
        group_by(sample_name) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ggplot() +
        geom_point(aes(x = mean_num_reads, y = me)) +
        mtclosed
    mysaveandstore("zztemp.324111111111111111111111111111111111.pdf")

    pfregions <- tempdf %>%
        group_by(region_name) %>%
        mutate(n_represented = n()) %>%
        ungroup() %>%
        filter(n_represented == 12) %>%
        left_join(sample_table) %>%
        group_by(region_name, sample_name, condition, sex, age) %>%
        summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
        ungroup()

    library(ggridges)

    pfmeans <- pf %>%
        group_by(condition) %>%
        summarise(me = mean(me))
    p <- pfregions %>%
        group_by(region_name, sample_name, condition) %>%
        summarise(me = mean(me)) %>%
        ggplot(aes(y = sample_name, x = me, fill = sample_name)) +
        geom_density_ridges2(quantile_lines = TRUE, quantiles = c(0.5)) +
        # geom_vline(data = pfmeans %>% filter(condition == condition2), mapping = aes(xintercept = me)) +
        # geom_vline(data = pfmeans %>% filter(condition == condition1), mapping = aes(xintercept = me)) +
        scale_samples_unique +
        ggtitle(sprintf("DMR Entropy %s - %s", condition2, condition1)) +
        mtclosed
    mysaveandstore(str_glue("ldna/results/m/plots/entropy_{subset}/entropy_dif_density_by_sample_unique111.pdf"), 5, 5, pl = p)


    p <- pfregions %>%
        group_by(region_name, condition) %>%
        summarise(me = mean(me)) %>%
        pivot_wider(names_from = condition, values_from = me) %>%
        mutate(dif = !!sym(condition2) - !!sym(condition1)) %>%
        ggplot() +
        geom_histogram(aes(x = dif)) +
        geom_vline(xintercept = 0) +
        ggtitle(sprintf("DMR Entropy %s - %s", condition2, condition1)) +
        mtclosed
    mysaveandstore(str_glue("ldna/results/m/plots/entropy_{subset}/entropy_dif_hist.pdf"), 5, 3.5, pl = p)
}

tempdf <- entdf %>%
    filter(type != "t05") %>%
    filter(type != "t01")
sf <- purrr::reduce(sflist, bind_rows)
pf <- tempdf %>%
    group_by(type, region_name) %>%
    mutate(n_represented = n()) %>%
    ungroup() %>%
    filter(n_represented == 12) %>%
    left_join(sample_table) %>%
    group_by(type, sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    mutate(stat_thresh = gsub(".*_", "", type)) %>%
    mutate(direction = gsub("_.*", "", type)) %>%
    ungroup()
p <- pf %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", facet.by = "direction", add = c("dotplot")) +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/m/plots/entropy_both/entropy_by_direction.pdf"), 6, 4, pl = p, sf = sf)
p <- pf %>%
    filter(direction == "hyper") %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", add = c("dotplot", "mean_se")) +
    labs(y = "Mean Entropy", x = "Condition") +
    ggtitle("Hypo") +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/m/plots/entropy_both/entropy_hyper.pdf"), 3.5, 4, pl = p, sf = sf)
p <- pf %>%
    filter(direction == "hypo") %>%
    ggviolin(x = "stat_thresh", y = "me", fill = "condition", add = c("dotplot", "mean_se")) +
    labs(y = "Mean Entropy", x = "Condition") +
    ggtitle("Hypo") +
    scale_conditions +
    mtopen
mysaveandstore(str_glue("ldna/results/m/plots/entropy_both/entropy_hypo.pdf"), 3.5, 4, pl = p, sf = sf)

pfregions <- tempdf %>%
    group_by(type, region_name) %>%
    mutate(n_represented = n()) %>%
    ungroup() %>%
    filter(n_represented == 12) %>%
    left_join(sample_table) %>%
    group_by(type, region_name, sample_name, condition, sex, age) %>%
    summarise(me = mean(mean_entropy), med_ent = mean(median_entropy), mean_num_reads = mean(mean_num_reads)) %>%
    ungroup()

p <- pfregions %>%
    group_by(type, region_name, condition) %>%
    summarise(me = mean(me)) %>%
    pivot_wider(names_from = condition, values_from = me) %>%
    mutate(dif = !!sym(condition2) - !!sym(condition1)) %>%
    ggplot() +
    geom_histogram(aes(x = dif)) +
    geom_vline(xintercept = 0) +
    facet_wrap(~type, scales = "free_y") +
    ggtitle(sprintf("DMR Entropy %s - %s", condition2, condition1)) +
    mtclosed
mysaveandstore(str_glue("ldna/results/m/plots/entropy_t05/entropy_dif_hist1.pdf"), 5, 5, pl = p)

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "1entropy", params$mod_code))
