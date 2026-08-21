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



##############
# METHYLATION CLUSTERING

outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

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
merged <- merged %>%
    mutate(consensus_pos = factor(consensus_pos, levels = cpg_order)) %>%
    filter(!is.na(seqnames))
library(tidyHeatmap)

dat <- merged %>%
    filter(!is.na(pctM))
dat %$% gene_id %>%
    unique() %>%
    length()
# write_csv(dat, "meth_for_bayes.csv")


# methylation heatmaps
hms <- list()
for (cur_sample in conf$samples) {
    hm_data <- merged %>%
        filter(sample == cur_sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50)
    if (nrow(hm_data) == 0) next
    p <- hm_data %>%
        heatmap(gene_id, consensus_pos, pctM,
            cluster_rows = TRUE, cluster_columns = FALSE,
            clustering_distance_rows = function(m) {
                d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                d[is.na(d)] <- 2
                d
            }
        ) %>%
        annotation_tile(intactness_req) %>%
        annotation_tile(loc_superlowres_integrative_stranded) %>%
        as_ComplexHeatmap()
    hms[[cur_sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 6)
}
# Generate the expression as a string and parse it
p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
rm(p)

hms <- list()
for (cur_sample in conf$samples) {
    hm_data <- merged %>%
        filter(sample == cur_sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50)
    if (nrow(hm_data) == 0) next
    p <- hm_data %>%
        heatmap(gene_id, consensus_pos, pctM,
            cluster_rows = TRUE, cluster_columns = FALSE,
            clustering_distance_rows = function(m) {
                d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                d[is.na(d)] <- 2
                d
            }
        ) %>%
        annotation_tile(genic_loc) %>%
        annotation_tile(intactness_req) %>%
        annotation_tile(loc_superlowres_integrative_stranded) %>%
        annotation_tile(refstatus) %>%
        as_ComplexHeatmap()
    hms[[cur_sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s_refstatus.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 8, h = 10)
}
# Generate the expression as a string and parse it
p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
mysaveandstore(sprintf("%s/%s_methylation_%s1_refstatus.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
rm(p)

# //ANCHOR - Meth clustering heatmaps
{
    # get clusters
    mat <- merged %>%
        group_by(gene_id, consensus_pos) %>%
        summarise(pctM = mean(pctM)) %>% # <-- collapse to 1 value!
        pivot_wider(names_from = consensus_pos, values_from = pctM) %>%
        column_to_rownames("gene_id") %>%
        as.matrix()
    dim(mat)
    mat <- mat[, colSums(!is.na(mat)) > 0.5 * nrow(mat)] # keep positions observed in >90% elements
    dim(mat)
    mat <- mat[rowSums(is.na(mat)) < (0.75 * ncol(mat)), ]
    dim(mat)
    mat <- mat %>%
        apply(1, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)) %>%
        t()

    library(mclust)
    mc <- Mclust(mat, G = 5)
    clusters <- mc$classification
    max(clusters)
    clustersdf <- tibble(gene_id = names(clusters), clusterNum = clusters)
    clustermap <- clustersdf %>%
        group_by(clusterNum) %>%
        summarise(cluster_n = n()) %>%
        arrange(-cluster_n) %>%
        mutate(rownum = row_number()) %>%
        mutate(cluster = LETTERS[rownum]) %>%
        dplyr::select(-rownum)
    clustersdf <- clustersdf %>% left_join(clustermap)
    merged <- merged %>% left_join(clustersdf)
    merged <- merged %>% mutate(pos_num = as.numeric(as.character(consensus_pos)))

    # methylation heatmaps
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- merged %>%
            filter(sample == cur_sample) %>%
            group_by(gene_id) %>%
            mutate(cpgs_detected_per_element = n()) %>%
            ungroup() %>%
            filter(cpgs_detected_per_element > 50)
        if (nrow(hm_data) == 0) next
        p <- hm_data %>%
            heatmap(gene_id, consensus_pos, pctM,
                cluster_rows = TRUE, cluster_columns = FALSE,
                clustering_distance_rows = function(m) {
                    d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                    d[is.na(d)] <- 2
                    d
                }
            ) %>%
            annotation_tile(intactness_req) %>%
            annotation_tile(loc_superlowres_integrative_stranded) %>%
            as_ComplexHeatmap()
        hms[[cur_sample]] <- p
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
    }
    # Generate the expression as a string and parse it
    p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
    rm(p)

    # methylation heatmaps
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- merged %>%
            filter(sample == cur_sample) %>%
            group_by(gene_id) %>%
            mutate(cpgs_detected_per_element = n()) %>%
            ungroup() %>%
            filter(cpgs_detected_per_element > 50)
        if (nrow(hm_data) == 0) next
        p <- hm_data %>%
            heatmap(gene_id, consensus_pos, pctM,
                cluster_rows = TRUE, cluster_columns = FALSE,
                clustering_distance_rows = function(m) {
                    d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                    d[is.na(d)] <- 2
                    d
                }
            ) %>%
            annotation_tile(intactness_req) %>%
            annotation_group(loc_superlowres_integrative_stranded) %>%
            as_ComplexHeatmap()
        hms[[cur_sample]] <- p
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(sprintf("%s/%s_methylation_%s_split.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
    }
    # Generate the expression as a string and parse it
    p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s_split.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
    rm(p)

    # methylation heatmaps
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- merged %>%
            filter(sample == cur_sample) %>%
            group_by(gene_id) %>%
            mutate(cpgs_detected_per_element = n()) %>%
            ungroup() %>%
            filter(cpgs_detected_per_element > 50)
        if (nrow(hm_data) == 0) next
        p <- hm_data %>%
            heatmap(gene_id, consensus_pos, pctM,
                cluster_rows = TRUE, cluster_columns = FALSE,
                clustering_distance_rows = function(m) {
                    d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                    d[is.na(d)] <- 2
                    d
                }
            ) %>%
            annotation_tile(intactness_req) %>%
            annotation_group(cluster) %>%
            as_ComplexHeatmap()
        hms[[cur_sample]] <- p
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(sprintf("%s/%s_methylation_%s_split_bycluster.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
    }
    # Generate the expression as a string and parse it
    p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s_split_bycluster.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
    rm(p)

    library(viridis)

    # not cluster split
    pf <- merged %>%
        group_by(condition, consensus_pos) %>%
        summarise(pctM = mean(pctM)) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos)))
    pf909 <- merged %>%
        group_by(condition, consensus_pos) %>%
        summarise(pctM = mean(pctM)) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
        filter(pos_num < 1000)
    pabs <- ggplot(pf, aes(x = pos_num, y = pctM, color = condition)) +
        geom_ribbon(
            data = pf %>% group_by(pos_num) %>%
                summarise(ymin = min(pctM), ymax = max(pctM)),
            aes(x = pos_num, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            alpha = 0.3,
            fill = "grey70"
        ) +
        geom_line(
            aes(color = condition),
            linewidth = 0.75
        ) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("%s/signal_wpoint_%s.pdf", outputdir_meth_clustering, subfam), pl = pabs, w = 5, h = 2)
    pabs909 <- ggplot(pf909, aes(x = pos_num, y = pctM, color = condition)) +
        geom_ribbon(
            data = pf909 %>% group_by(pos_num) %>%
                summarise(ymin = min(pctM), ymax = max(pctM)),
            aes(x = pos_num, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            alpha = 0.3,
            fill = "grey70"
        ) +
        geom_line(
            aes(color = condition),
            linewidth = 0.75
        ) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("%s/signal_909wpoint_%s.pdf", outputdir_meth_clustering, subfam), pl = pabs909, w = 5, h = 2)

    pf_wide <- merged %>%
        group_by(condition, consensus_pos) %>%
        summarise(pctM = mean(pctM)) %>%
        pivot_wider(names_from = condition, values_from = pctM) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos)))
    for (contrast in contrasts) {
        cp <- parse_contrast(contrast)
        cond1 <- cp$condition1
        cond2 <- cp$condition2
        pfdif <- pf_wide %>% mutate(dif = !!sym(cond1) - !!sym(cond2))
        pfdif909 <- pfdif %>% filter(pos_num < 1000)

        pdif <- pfdif %>%
            ggplot(aes(x = pos_num, y = dif)) +
            geom_line(color = "grey") +
            geom_point(
                aes(fill = dif),
                shape = 21, color = "black", size = 2, stroke = 0.4
            ) +
            scale_fill_viridis(option = "A") +
            mtclosed
        mysaveandstore(sprintf("%s/%s/signaldif_wpoint_%s.pdf", outputdir_meth_clustering, contrast, subfam), pl = pdif, w = 5, h = 2)
        patch <- wrap_plots(list(pabs, pdif), ncol = 1, guides = "collect", axes = "collect")
        mysaveandstore(sprintf("%s/%s/patch_full_%s.pdf", outputdir_meth_clustering, contrast, subfam), pl = patch, w = 5, h = 4)

        pdif909 <- pfdif909 %>%
            ggplot(aes(x = pos_num, y = dif)) +
            geom_line(color = "grey") +
            geom_point(
                aes(fill = dif),
                shape = 21, color = "black", size = 2, stroke = 0.4
            ) +
            scale_fill_viridis(option = "A") +
            mtclosed
        mysaveandstore(sprintf("%s/%s/signaldif_909wpoint_%s.pdf", outputdir_meth_clustering, contrast, subfam), pl = pdif909, w = 5, h = 2)
        patch <- wrap_plots(list(pabs909, pdif909), ncol = 1, guides = "collect", axes = "collect")
        mysaveandstore(sprintf("%s/%s/patch_909_%s.pdf", outputdir_meth_clustering, contrast, subfam), pl = patch, w = 5, h = 4)
    }


    # with clustering
    pf <- merged %>%
        group_by(condition, consensus_pos, cluster) %>%
        summarise(pctM = mean(pctM)) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos)))
    pf909 <- merged %>%
        group_by(condition, consensus_pos, cluster) %>%
        summarise(pctM = mean(pctM)) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
        filter(pos_num < 1000)
    pabs_clust <- ggplot(pf, aes(x = pos_num, y = pctM, color = condition)) +
        geom_ribbon(
            data = pf %>% group_by(pos_num, cluster) %>%
                summarise(ymin = min(pctM), ymax = max(pctM)),
            aes(x = pos_num, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            alpha = 0.3,
            fill = "grey70"
        ) +
        geom_line(
            aes(color = condition),
            linewidth = 0.75
        ) +
        facet_wrap(~cluster, nrow = 1) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("%s/signal_wpoint_%s_clustered.pdf", outputdir_meth_clustering, subfam), pl = pabs_clust, w = 12, h = 2)
    pabs909_clust <- ggplot(pf909, aes(x = pos_num, y = pctM, color = condition)) +
        geom_ribbon(
            data = pf909 %>% group_by(pos_num, cluster) %>%
                summarise(ymin = min(pctM), ymax = max(pctM)),
            aes(x = pos_num, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            alpha = 0.2,
            fill = "grey50"
        ) +
        geom_line(
            aes(color = condition),
            linewidth = 0.75
        ) +
        facet_wrap(~cluster, nrow = 1) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("%s/signal_909wpoint_%s_clustered.pdf", outputdir_meth_clustering, subfam), pl = pabs909_clust, w = 20, h = 2)

    pf_wide_clust <- merged %>%
        group_by(condition, consensus_pos, cluster) %>%
        summarise(pctM = mean(pctM)) %>%
        pivot_wider(names_from = condition, values_from = pctM) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos)))
    for (contrast in contrasts) {
        cp <- parse_contrast(contrast)
        cond1 <- cp$condition1
        cond2 <- cp$condition2
        pfdif <- pf_wide_clust %>% mutate(dif = !!sym(cond1) - !!sym(cond2))
        pfdif909 <- pfdif %>% filter(pos_num < 1000)

        pdif <- pfdif %>%
            ggplot(aes(x = pos_num, y = dif)) +
            geom_hline(yintercept = 0, color = "darkgrey") +
            geom_line(color = "grey") +
            geom_point(
                aes(fill = dif),
                shape = 21, color = "black", size = 2, stroke = 0.4
            ) +
            facet_wrap(~cluster, nrow = 1) +
            scale_fill_viridis(option = "A") +
            mtclosed
        mysaveandstore(sprintf("%s/%s/signaldif_wpoint_%s_clustered.pdf", outputdir_meth_clustering, contrast, subfam), pl = pdif, w = 12, h = 2)
        patch <- wrap_plots(list(pabs_clust, pdif), ncol = 1, guides = "collect", axes = "collect")
        mysaveandstore(sprintf("%s/%s/patch_full_%s_clustered.pdf", outputdir_meth_clustering, contrast, subfam), pl = patch, w = 12, h = 4)

        pdif909 <- pfdif909 %>%
            ggplot(aes(x = pos_num, y = dif)) +
            geom_hline(yintercept = 0, color = "darkgrey") +
            geom_line(color = "grey") +
            geom_point(
                aes(fill = dif),
                shape = 21, color = "black", size = 2, stroke = 0.4
            ) +
            facet_wrap(~cluster, nrow = 1) +
            scale_fill_viridis(option = "A") +
            mtclosed
        mysaveandstore(sprintf("%s/%s/signaldif_909wpoint_%s_clustered.pdf", outputdir_meth_clustering, contrast, subfam), pl = pdif909, w = 20, h = 2)
        patch <- wrap_plots(list(pabs909_clust, pdif909), ncol = 1, guides = "collect", axes = "collect")
        mysaveandstore(sprintf("%s/%s/patch_909_%s_clustered.pdf", outputdir_meth_clustering, contrast, subfam), pl = patch, w = 12, h = 4)
    }



    pfvar <- merged %>%
        group_by(condition, consensus_pos) %>%
        summarise(mv = var(pctM)) %>%
        group_by(consensus_pos) %>%
        summarise(mv = mean(mv)) %>%
        arrange(-mv) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
        mutate(mvZ = as.numeric(scale(mv)))
    p1 <- pfvar %>%
        ggplot(aes(x = pos_num, y = mvZ)) +
        geom_hline(yintercept = 0, color = "darkgrey") +
        geom_line(color = "grey") +
        geom_point(
            aes(fill = mvZ),
            shape = 21, # allows fill + color
            color = "black", # outline
            size = 2,
            stroke = 0.4 # outline thickness
        ) +
        scale_fill_viridis(option = "D") +
        mtclosed
    mysaveandstore(sprintf("%s/cpg_within_cond_var_%s.pdf", outputdir_meth_clustering, subfam), pl = p1, w = 5, h = 2)
    pfvar <- merged %>%
        group_by(condition, consensus_pos) %>%
        summarise(mean_pctM = mean(pctM)) %>%
        group_by(consensus_pos) %>%
        summarise(mv = var(mean_pctM)) %>%
        arrange(-mv) %>%
        mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
        mutate(mvZ = as.numeric(scale(mv)))
    p2 <- pfvar %>%
        ggplot(aes(x = pos_num, y = mvZ)) +
        geom_hline(yintercept = 0, color = "darkgrey") +
        geom_line(color = "grey") +
        geom_point(
            aes(fill = mvZ),
            shape = 21, # allows fill + color
            color = "black", # outline
            size = 2,
            stroke = 0.4 # outline thickness
        ) +
        scale_fill_viridis(option = "A") +
        mtclosed
    mysaveandstore(sprintf("%s/cpg_between_cond_var_%s.pdf", outputdir_meth_clustering, subfam), pl = p2, w = 5, h = 2)

    patch <- wrap_plots(list(p1, p2), ncol = 1, guides = "collect", axes = "collect")
    mysaveandstore(sprintf("%s/cpg_var_%s.pdf", outputdir_meth_clustering, subfam), pl = patch, w = 5, h = 4)



    cdf <- merged %>%
        mutate(region = case_when(
            pos_num <= 329 ~ "HM",
            (pos_num > 329) & (pos_num <= 600) ~ "ASP",
            (pos_num <= 910) & (pos_num > 600) ~ "PostASP",
            pos_num > 910 ~ "Body",
        )) %>%
        group_by(sample, condition, gene_id, region) %>%
        summarise(pctM = mean(pctM)) %>%
        ungroup()

    cdf_wide <- cdf %>%
        pivot_wider(names_from = region, values_from = pctM)

    cor(cdf_wide %>% dplyr::select(HM, ASP, PostASP, Body), use = "pairwise.complete.obs")
    for (cond in unique(cdf_wide$condition)) {
        cat(sprintf("\nCorrelation matrix for %s:\n", cond))
        cormat <- cor(cdf_wide %>% filter(condition == cond) %>% dplyr::select(HM, ASP, PostASP, Body), use = "pairwise.complete.obs")
        print(cormat)
    }

    cor_pair <- function(df, x, y) {
        ct <- cor.test(df[[x]], df[[y]], use = "pairwise.complete.obs", method = "spearman")
        tibble(
            region1 = x,
            region2 = y,
            cor = ct$estimate,
            p.value = ct$p.value
        )
    }
    region_cols <- c("HM", "ASP", "Body")
    region_pairs <- expand.grid(region1 = region_cols, region2 = region_cols, stringsAsFactors = FALSE)

    cor_df <- cdf_wide %>%
        group_by(condition) %>%
        group_modify(~ bind_rows(
            lapply(1:nrow(region_pairs), function(i) {
                cor_pair(.x, region_pairs$region1[i], region_pairs$region2[i])
            })
        )) %>%
        ungroup()


    p <- cor_df %>% ggplot(aes(x = region1, y = region2, fill = cor)) +
        geom_tile() +
        geom_text(aes(label = round(cor, 2)), color = "black", size = 4) + # overlay R²
        facet_wrap(~condition) +
        scale_fill_gradientn(
            colours = RColorBrewer::brewer.pal(4, "Oranges")
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) + # No padding; bottom = high rank
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(x = "sample", y = "genes (sorted by prop_in_bin per sample)") +
        mtclosed +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )
    mysaveandstore(sprintf("%s/region_cor_%s.pdf", outputdir_meth_clustering, subfam), w = 5, h = 4)
    cor_df_tri <- cor_df %>%
        filter(as.numeric(factor(region1)) >= as.numeric(factor(region2)))
    p <- cor_df_tri %>% ggplot(aes(x = region1, y = region2, fill = cor)) +
        geom_tile() +
        geom_text(aes(label = round(cor, 2)), color = "black", size = 4) + # overlay R²
        facet_wrap(~condition) +
        scale_fill_gradientn(
            colours = RColorBrewer::brewer.pal(4, "Oranges")
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) + # No padding; bottom = high rank
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(x = "sample", y = "genes (sorted by prop_in_bin per sample)") +
        mtclosed +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )
    mysaveandstore(sprintf("%s/region_cor_triangle_%s.pdf", outputdir_meth_clustering, subfam), w = 5, h = 4)

    # Correlation heatmap across CpG consensus positions (per condition)
    # Each CpG is one column/row, ordered by numeric position but spaced evenly
    cpg_levels <- cg_positions_df %>%
        arrange(consensus_pos) %>%
        pull(consensus_pos) %>%
        unique()
    pos_wide <- merged %>%
        filter(consensus_pos %in% cpg_levels) %>%
        group_by(gene_id, condition, consensus_pos) %>%
        summarise(pctM = mean(pctM)) %>%
        pivot_wider(names_from = consensus_pos, values_from = pctM, names_sort = TRUE)

    for (cond in unique(pos_wide$condition)) {
        mat <- pos_wide %>%
            filter(condition == cond) %>%
            ungroup() %>%
            dplyr::select(-gene_id, -condition) %>%
            as.matrix()
        # Reorder columns by numeric position
        col_order <- order(as.numeric(colnames(mat)))
        mat <- mat[, col_order]
        pos_cor <- cor(mat, use = "pairwise.complete.obs")
        cpg_factor <- factor(colnames(pos_cor), levels = colnames(pos_cor))

        p <- as.data.frame(as.table(pos_cor)) %>%
            mutate(
                Var1 = factor(Var1, levels = levels(cpg_factor)),
                Var2 = factor(Var2, levels = levels(cpg_factor))
            ) %>%
            ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
            geom_raster() +
            scale_fill_gradient2(
                low = "blue", mid = "white", high = "red", midpoint = 0,
                name = "r", limits = c(-1, 1)
            ) +
            labs(
                x = "CpG", y = "CpG",
                title = sprintf("%s CpG-CpG Correlation (%s)", subfam, cond)
            ) +
            coord_fixed() +
            mtclosed +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
                axis.text.y = element_text(size = 4)
            )
        mysaveandstore(sprintf("%s/position_cor_heatmap_%s_%s.pdf", outputdir_meth_clustering, cond, subfam), pl = p, w = 10, h = 9)
    }
}








# nonref analysis
grs_nr_df <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/ldna/Rintermediates/m/grsdf_nonref.tsv")
grs_nr <- GRanges(grs_nr_df)

rmannextended %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(refstatus == "Ref") %>%
    filter(intactness_req == "Intact") %>%
    pl()
rmannextended_nr_list <- list()
merged_nr_list <- list()
for (sample in sample_table$sample_name) {
    rmannextended_nr_temp <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
    rmannextended_nr_temp$sample_name <- sample
    rmannextended_nr_list[[sample]] <- rmannextended_nr_temp
    grs_nr_temp <- grs_nr[mcols(grs_nr)$sample == sample]
    merged_temp <- merge_with_grs(grs_nr_temp, GRanges(rmannextended_nr_temp))
    merged_nr_list[[sample]] <- merged_temp
}

rmannextended_nr <- bind_rows(rmannextended_nr_list) %>%
    tibble() %>%
    mutate(gene_id = paste0(sample_name, "___", gene_id)) %>%
    mutate(seqnames = paste0(sample_name, "___", seqnames))

merged_nr <- bind_rows(merged_nr_list) %>%
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
for (cur_sample in conf$samples) {
    hm_data <- merged %>%
        filter(sample == cur_sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        distinct()
    if (nrow(hm_data) == 0) next
    p <- hm_data %>%
        heatmap(gene_id, consensus_pos, pctM,
            cluster_rows = TRUE, cluster_columns = FALSE,
            clustering_distance_rows = function(m) {
                d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
                d[is.na(d)] <- 2
                d
            }
        ) %>%
        as_ComplexHeatmap()
    hms[[cur_sample]] <- p
    dir.create(outputdir_meth_clustering, recursive = TRUE)
    mysaveandstore(sprintf("%s/nonref_%s_methylation_%s.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 3)
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
    heatmap(gene_id, consensus_pos, pctM,
        cluster_rows = TRUE, cluster_columns = FALSE,
        clustering_distance_rows = function(m) {
            d <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
            d[is.na(d)] <- 2
            d
        }
    ) %>%
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
                seqval <- dfs %>% filter(consensus_pos == start_pos) %$% sequence_pos %>% pluck(1)
                if (!is.null(seqval)) {
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
    mapping <- tibble(gene_id = names(filter_pos_list), filter_pos = unlist(filter_pos_list)) %>%
        filter(filter_pos != "NoRegionHomology") %>%
        mutate(filter_pos = as.numeric(filter_pos))
    gene_ids_with_homology <- mapping %$% gene_id
    fl_grs_with_homology <- fl_grs[mcols(fl_grs)$gene_id %in% gene_ids_with_homology]
    grlist <- map(seq_along(fl_grs_with_homology), function(x) fl_grs_with_homology[x])
    grlistresized <- map(grlist, function(x) {
        if (mcols(x)$gene_id %in% mapping$gene_id) {
            new_width <- mapping %>% filter(gene_id == mcols(x)$gene_id) %$% filter_pos %>% pluck(1)
            if (!is.null(new_width) && !is.na(new_width) && new_width > 0) {
                return(resize(x, width = new_width))
            }
        }
        return(NULL)
    })
    grlistresized <- purrr::compact(grlistresized)
    number_omitted <- length(grlist) - length(grlistresized)
    if (length(grlistresized) == 0) {
        print("No resized GRanges produced")
        return(GRanges())
    }
    l1hs_resized <- purrr::reduce(grlistresized, c)
    return(l1hs_resized)
}

flL1HS5UTR <- filter_by_consensus_pos(fl_grs = rmannextended_nr %>% filter(rte_subfamily == "L1HS") %>% filter(rte_length_req == "FL") %>% GRanges(), pos_mapping = consensus_index_long, include_up_to_pos = 909)

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

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "17", params$mod_code))
