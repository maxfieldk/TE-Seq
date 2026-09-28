module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")

cell_fractions <- read_delim("ldna/results/m/tables/scMD_cell_type_fractions_all.csv")
cell_fractions_scaled <- cell_fractions %>%
    mutate(Neuron = Inh + Exc) %>%
    mutate(across(all_of(c(cell_types_cols, "Neuron")), ~ as.numeric(scale(.)), .names = "{.col}_z"))

sample_table <- sample_table %>% left_join(cell_fractions_scaled)
sample_table <- sample_table %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = sample_name)
samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

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
library(emmeans)

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
rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)
yngl1unfiltdf <- read_delim(sprintf("ldna/Rintermediates/%s/yl1sunfilt", params$mod_code), col_names = TRUE)


library(tidyHeatmap)
library(emmeans)

consensus_dir <- "ldna/results/m/plots/l1_alignment_meth"
lm_categorical_terms <- c("condition", adjustment_set_categorical)
sample_order <- sample_table %>%
    arrange(condition) %$%
    sample_name

# the per sample heatmaps get composed side by side, which ComplexHeatmap only allows when every
# heatmap has the same rows. an element covered in one sample and not another would otherwise drop
# out of that sample's heatmap, so it is given a single NA cell there instead and renders blank.
pad_heatmap_elements <- function(df, min_cpgs = 50) {
    annotation_cols <- c(
        "intactness_req", "loc_superlowres_integrative_stranded",
        "genic_loc", "refstatus", "cluster", "meth_group"
    )
    kept <- df %>%
        group_by(sample, gene_id) %>%
        filter(n() > min_cpgs) %>%
        ungroup()
    pad <- expand_grid(sample = unique(kept$sample), gene_id = unique(kept$gene_id)) %>%
        anti_join(distinct(kept, sample, gene_id), by = c("sample", "gene_id")) %>%
        mutate(consensus_pos = kept$consensus_pos[1], pctM = NA_real_)
    # these columns describe the element, not the measurement, so they are stripped off and
    # rejoined once after padding. carrying them onto the pad rows separately let a gene end up
    # labelled on its observed rows and unlabelled on its padded ones.
    annotations <- kept %>%
        dplyr::select(gene_id, any_of(annotation_cols)) %>%
        distinct(gene_id, .keep_all = TRUE)
    detection <- kept %>%
        group_by(gene_id) %>%
        summarise(element_detection = ifelse(n_distinct(sample) == n_distinct(kept$sample),
            "All samples", "Absent from some"
        ), .groups = "drop")
    bind_rows(kept %>% dplyr::select(-any_of(annotation_cols)), pad) %>%
        left_join(annotations, by = "gene_id") %>%
        left_join(detection, by = "gene_id")
}
# unfilt keeps non-ref elements that only a handful of samples carry, so every subfamily is
# analysed twice: over all elements, and over only those with data in every sample
analysis_grid <- expand_grid(
    subfam = c("L1HS", "L1PA2", "L1PA3"),
    filt_type = c("filt", "unfilt"),
    element_set = c("all_elements", "complete_elements")
)

for (grid_row in seq_len(nrow(analysis_grid))) {
    subfam <- analysis_grid$subfam[grid_row]
    filt_type <- analysis_grid$filt_type[grid_row]
    element_set <- analysis_grid$element_set[grid_row]
    outputdir_meth_clustering <- sprintf("%s/%s/%s", consensus_dir, filt_type, element_set)
    print(sprintf("%s %s %s", subfam, filt_type, element_set))

    consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", consensus_dir, subfam))


    consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", consensus_dir, subfam)
    consensus_ss <- readDNAStringSet(consensus_path)

    cg_indices <- consensus_ss %>%
        vmatchPattern(pattern = "CG") %>%
        start() %>%
        unlist() %>%
        as.numeric()


    cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

    methdf <- (if (filt_type == "filt") rtedf else yngl1unfiltdf) %>%
        filter(rte_subfamily == subfam) %>%
        filter(refstatus == "Ref", rte_length_req == "FL")

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
        mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
        # condition arrives as character off the tsv, which orders it alphabetically in every plot
        mutate(condition = factor(condition, levels = conf$levels)) %>%
        filter(!is.na(seqnames))
    merged %$% gene_id %>% unique()


    if (element_set == "complete_elements") {
        merged <- merged %>%
            group_by(gene_id) %>%
            filter(n_distinct(sample) == length(conf$samples)) %>%
            ungroup()
    }
    if (nrow(merged) == 0) next

    hm_all <- pad_heatmap_elements(merged)


    for (subset in c(150, 328, 500, 909)) {
        gene_meth_summary <- merged %>%
            filter(as.numeric(as.character(consensus_pos)) <= subset) %>%
            group_by(gene_id, sample) %>%
            summarise(mean_meth = mean(pctM), .groups = "drop") %>%
            left_join(sample_table)

        gene_lm_results <- gene_meth_summary %>%
            group_by(gene_id) %>%
            filter(n_distinct(condition) > 1) %>%
            group_modify(~ {
                gene_df <- .x
                # a sparsely covered element can leave a categorical predictor with a single
                # value, which lm cannot build contrasts from
                if (any(map_dbl(lm_categorical_terms, ~ n_distinct(gene_df[[.x]], na.rm = TRUE)) < 2)) {
                    return(tibble())
                }
                fit <- lm(reformulate(lm_right_hand_side, "mean_meth"), gene_df)
                # too few samples to support the adjusted model; emmeans errors on a rank
                # deficient fit rather than returning NA
                if (fit$df.residual < 1 || anyNA(coef(fit))) {
                    return(tibble())
                }
                emmeans(fit, ~condition) %>%
                    contrast(method = "pairwise") %>%
                    as.data.frame()
            }) %>%
            ungroup() %>%
            dplyr::rename(term = contrast) %>%
            group_by(term) %>%
            mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
            ungroup() %>%
            arrange(term, p.adjust)

        gene_meth_wide <- gene_meth_summary %>%
            dplyr::select(gene_id, sample, mean_meth) %>%
            pivot_wider(names_from = sample, values_from = mean_meth)

        gene_lm_results <- gene_lm_results %>%
            left_join(gene_meth_wide, by = "gene_id")

        lm_table_dir <- str_glue("ldna/results/{params$mod_code}/tables/l1_alignment_meth/{filt_type}/{element_set}")
        dir.create(lm_table_dir, recursive = TRUE, showWarnings = FALSE)
        write_delim(gene_lm_results, str_glue("{lm_table_dir}/gene_meth_promoter{subset}_lm_{subfam}.tsv"), delim = "\t")

        gene_lm_results %>% arrange(p.value)

        top_genes <- gene_lm_results %>%
            arrange(p.value) %>%
            distinct(gene_id) %>%
            head(n = 15) %$%
            gene_id

        p <- gene_meth_summary %>%
            filter(gene_id %in% top_genes) %>%
            mutate(sample = factor(sample, levels = sample_order)) %>%
            group_by(condition) %>%
            heatmap(gene_id, sample, mean_meth,
                cluster_rows = TRUE, cluster_columns = FALSE,
                palette_value = circlize::colorRamp2(seq(0, 100, length.out = 100), viridis::viridis(100, option = "plasma"))
            ) %>%
            annotation_tile(condition)
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(fn = str_glue("{outputdir_meth_clustering}/top15_gene_meth_heatmap1_promoter{subset}_{subfam}.pdf"), w = 7, h = 6, res = 300, pl = p)

        p <- gene_meth_summary %>%
            filter(gene_id %in% top_genes) %>%
            group_by(gene_id) %>%
            mutate(mean_meth_z = as.numeric(scale(mean_meth))) %>%
            ungroup() %>%
            mutate(sample = factor(sample, levels = sample_order)) %>%
            group_by(condition) %>%
            heatmap(gene_id, sample, mean_meth_z,
                cluster_rows = TRUE, cluster_columns = FALSE,
                palette_value = circlize::colorRamp2(c(-2, 0, 2), c("red", "white", "blue"))
            ) %>%
            annotation_tile(condition)
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(fn = str_glue("{outputdir_meth_clustering}/top15_gene_meth_heatmap1_promoter{subset}_zscore_{subfam}.pdf"), w = 7, h = 6, res = 300, pl = p)

        contrast_levels <- map_chr(conf$contrasts, function(cs) {
            cp <- parse_contrast(cs)
            paste(cp$condition1, "-", cp$condition2)
        })

        gene_order <- gene_lm_results %>%
            filter(term == contrast_levels[1]) %>%
            arrange(estimate) %$%
            gene_id

        p <- gene_lm_results %>%
            mutate(
                term = factor(term, levels = contrast_levels),
                gene_id = factor(gene_id, levels = gene_order)
            ) %>%
            ggplot(aes(x = estimate, y = gene_id)) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_pointrange(aes(xmin = estimate - 1.96 * SE, xmax = estimate + 1.96 * SE, color = p.adjust < 0.05), fatten = 0.5) +
            scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), name = "p.adjust < 0.05") +
            facet_wrap(~term, scales = "free_x") +
            labs(x = "Methylation difference (95% CI)", y = "") +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            mtclosed
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(fn = str_glue("{outputdir_meth_clustering}/gene_meth_lm_coefficients_promoter{subset}_{subfam}1.pdf"), w = 9, h = 9, res = 300, pl = p)
    }



    # methylation heatmaps
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- hm_all %>% filter(sample == cur_sample)
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
    p <- base::eval(base::parse(text = paste0("hms[['", names(hms), "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
    rm(p)

    # same heatmap, rows split by whether the element is detected in every sample
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- hm_all %>% filter(sample == cur_sample)
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
            annotation_group(element_detection) %>%
            as_ComplexHeatmap()
        hms[[cur_sample]] <- p
        dir.create(outputdir_meth_clustering, recursive = TRUE)
        mysaveandstore(sprintf("%s/%s_methylation_%s_split_bydetection.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
    }
    p <- base::eval(base::parse(text = paste0("hms[['", names(hms), "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s_split_bydetection.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
    rm(p)

    # methylation heatmaps
    hms <- list()
    for (cur_sample in conf$samples) {
        hm_data <- hm_all %>% filter(sample == cur_sample)
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
        mysaveandstore(sprintf("%s/%s_methylation_%s_split_by_loc.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
    }
    # Generate the expression as a string and parse it
    p <- base::eval(base::parse(text = paste0("hms[['", names(hms), "']]", collapse = " + ")))
    mysaveandstore(sprintf("%s/%s_methylation_%s_split_by_loc.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
    rm(p)



    # //ANCHOR - Meth clustering heatmaps
    {
        # get clusters
        # one row per element, averaged over samples at each consensus CpG, so an element missing
        # from some samples still contributes whatever CpGs it does have
        mat <- merged %>%
            group_by(gene_id, consensus_pos) %>%
            summarise(pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            pivot_wider(names_from = consensus_pos, values_from = pctM) %>%
            column_to_rownames("gene_id") %>%
            as.matrix()
        dim(mat)
        mat <- mat[, colSums(!is.na(mat)) > 0.5 * nrow(mat), drop = FALSE] # keep positions observed in >50% elements
        dim(mat)
        mat <- mat[rowSums(is.na(mat)) < (0.75 * ncol(mat)), , drop = FALSE]
        dim(mat)
        # an unobserved CpG takes that position's average over elements rather than the element's
        # own mean, which would flatten a sparsely covered element into a near constant profile
        # and leave it clustering on overall level instead of shape
        position_means <- colMeans(mat, na.rm = TRUE)
        mat <- t(apply(mat, 1, function(x) ifelse(is.na(x), position_means, x)))

        library(mclust)
        # Mclust cannot estimate covariances over ~90 CpG dimensions from a few hundred elements,
        # so it is fitted on the leading principal components instead
        pca <- prcomp(mat, center = TRUE, scale. = FALSE)
        pc_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
        # keep the PCs holding the first 90% of the variance rather than a fixed count, so trailing
        # noise directions are not handed to Mclust as if they were signal. left unscaled, since
        # every CpG is on the same percent methylation scale and the spread across PCs is real.
        # capped, because Mclust cannot estimate covariances once the dimension count climbs.
        n_pcs <- min(which(cumsum(pc_var) >= 90)[1], 10)
        pcs <- pca$x[, seq_len(n_pcs), drop = FALSE]
        # both the cluster count and the covariance model are chosen by BIC rather than fixed
        mc <- Mclust(pcs, G = 1:9)
        print(sprintf(
            "%s: %s PCs (%.0f%% of variance), BIC picked %s clusters, model %s",
            subfam, n_pcs, cumsum(pc_var)[n_pcs], mc$G, mc$modelName
        ))
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

        # BIC against cluster count, per covariance model, so the choice above can be judged
        p <- mc$BIC %>%
            unclass() %>%
            as.data.frame() %>%
            rownames_to_column("G") %>%
            pivot_longer(-G, names_to = "model", values_to = "BIC") %>%
            filter(!is.na(BIC)) %>%
            mutate(G = as.integer(G)) %>%
            ggplot(aes(x = G, y = BIC, color = model)) +
            geom_vline(xintercept = mc$G, linetype = "dashed", color = "grey40") +
            geom_line() +
            geom_point(size = 1) +
            labs(
                x = "Clusters", y = "BIC",
                title = sprintf(
                    "%s: %s PCs, BIC picked %s clusters (%s)",
                    subfam, n_pcs, mc$G, mc$modelName
                )
            ) +
            mtclosed
        mysaveandstore(sprintf("%s/cluster_bic_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 7, h = 4)

        pcdf <- pca$x[, 1:2] %>%
            as.data.frame() %>%
            rownames_to_column("gene_id") %>%
            left_join(clustersdf, by = "gene_id")
        # arrows for the consensus CpGs loading most strongly on PC1/PC2, scaled onto the scores
        loadings <- pca$rotation[, 1:2] %>%
            as.data.frame() %>%
            rownames_to_column("consensus_pos") %>%
            mutate(load = sqrt(PC1^2 + PC2^2)) %>%
            slice_max(load, n = 10) %>%
            mutate(across(c(PC1, PC2), ~ .x * max(abs(pcdf$PC1)) / max(load)))
        p <- pcdf %>%
            ggplot(aes(x = PC1, y = PC2)) +
            geom_point(aes(color = cluster), size = 1.5, alpha = 0.8) +
            geom_segment(
                data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                arrow = arrow(length = unit(0.15, "cm")), color = "grey30"
            ) +
            ggrepel::geom_text_repel(data = loadings, aes(label = consensus_pos), size = 2.5, color = "grey20") +
            labs(
                x = sprintf("PC1 (%.1f%%)", pc_var[1]), y = sprintf("PC2 (%.1f%%)", pc_var[2]),
                title = sprintf("%s element methylation profiles", subfam)
            ) +
            mtclosed
        mysaveandstore(sprintf("%s/cluster_pca_biplot_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 6, h = 5)

        library(uwot)
        set.seed(123)
        p <- umap(pcs, n_neighbors = min(15, nrow(pcs) - 1), min_dist = 0.1) %>%
            as.data.frame() %>%
            setNames(c("UMAP1", "UMAP2")) %>%
            mutate(gene_id = rownames(pcs)) %>%
            left_join(clustersdf, by = "gene_id") %>%
            ggplot(aes(x = UMAP1, y = UMAP2, color = cluster)) +
            geom_point(size = 1.5, alpha = 0.8) +
            labs(title = sprintf("%s element methylation profiles", subfam)) +
            mtclosed
        mysaveandstore(sprintf("%s/cluster_umap_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 6, h = 5)
        # joined into the heatmap frame rather than into merged, so merged stays written once per
        # grid iteration and rerunning any part of this block cannot produce cluster.x/cluster.y
        hm_all <- pad_heatmap_elements(merged %>% left_join(clustersdf, by = "gene_id"))

        # methylation heatmaps
        hms <- list()
        for (cur_sample in conf$samples) {
            hm_data <- hm_all %>% filter(sample == cur_sample)
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
        p <- base::eval(base::parse(text = paste0("hms[['", names(hms), "']]", collapse = " + ")))
        mysaveandstore(sprintf("%s/%s_methylation_%s_split_bycluster.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
        rm(p)

        # //ANCHOR - methylation over the three consensus regions
        region_labels <- c(r1 = "1-329", r2 = "330-720", r3 = "721+")
        add_region <- function(df) {
            df %>% mutate(region = factor(
                case_when(
                    pos_num <= 329 ~ "r1",
                    pos_num <= 720 ~ "r2",
                    TRUE ~ "r3"
                ),
                levels = names(region_labels)
            ))
        }
        region_meth_per_sample <- merged %>%
            add_region() %>%
            group_by(sample, condition, gene_id, region) %>%
            summarise(region_meth = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            mutate(sample = factor(sample, levels = sample_order))
        region_meth_per_element <- merged %>%
            add_region() %>%
            group_by(gene_id, region) %>%
            summarise(region_meth = mean(pctM, na.rm = TRUE), .groups = "drop")
        region_facet <- facet_wrap(~region, nrow = 1, labeller = labeller(region = region_labels))

        # where the mass sits in each region, and whether 50% is a sensible place to cut
        p <- region_meth_per_sample %>%
            ggplot(aes(x = region_meth, color = condition)) +
            geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
            geom_density() +
            region_facet +
            scale_conditions +
            labs(x = "Element methylation over region (%)", y = "Density") +
            mtclosed
        mysaveandstore(sprintf("%s/region_meth_density_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 3)

        # the same as a cumulative distribution, which reads off what fraction any cut would take
        p <- region_meth_per_sample %>%
            ggplot(aes(x = region_meth, color = condition)) +
            geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
            stat_ecdf() +
            region_facet +
            scale_conditions +
            labs(x = "Element methylation over region (%)", y = "Cumulative fraction of elements") +
            mtclosed
        mysaveandstore(sprintf("%s/region_meth_ecdf_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 3)

        # per sample, to catch a sample sitting away from the rest of its condition
        p <- region_meth_per_sample %>%
            ggplot(aes(x = sample, y = region_meth)) +
            geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
            geom_violin(aes(fill = condition), color = "black", linewidth = 0.2) +
            geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
            facet_wrap(~region, nrow = 3, labeller = labeller(region = region_labels)) +
            scale_conditions +
            labs(x = "", y = "Element methylation over region (%)") +
            mtclosed +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        mysaveandstore(sprintf("%s/region_meth_per_sample_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 8)

        # one value per element per condition, so the condition comparison is not inflated by
        # counting every sample separately
        p <- region_meth_per_sample %>%
            group_by(condition, gene_id, region) %>%
            summarise(region_meth = mean(region_meth), .groups = "drop") %>%
            ggplot(aes(x = condition, y = region_meth)) +
            geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
            geom_violin(aes(fill = condition), color = "black", show.legend = FALSE) +
            geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2) +
            region_facet +
            scale_conditions +
            labs(x = "", y = "Element methylation over region (%)") +
            mtclosed
        mysaveandstore(sprintf("%s/region_meth_per_condition_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 8, h = 3.5)

        # how the regions relate within an element, which is what decides whether all eight
        # groups actually get populated or only the diagonal ones do
        region_wide <- region_meth_per_element %>%
            pivot_wider(names_from = region, values_from = region_meth, names_expand = TRUE)
        p <- bind_rows(
            region_wide %>% transmute(x = r1, y = r2, pair = "1-329 vs 330-720"),
            region_wide %>% transmute(x = r1, y = r3, pair = "1-329 vs 721+"),
            region_wide %>% transmute(x = r2, y = r3, pair = "330-720 vs 721+")
        ) %>%
            ggplot(aes(x = x, y = y)) +
            geom_hline(yintercept = 50, linetype = "dashed", color = "firebrick") +
            geom_vline(xintercept = 50, linetype = "dashed", color = "firebrick") +
            geom_point(size = 0.6, alpha = 0.4) +
            facet_wrap(~pair, nrow = 1) +
            coord_fixed() +
            labs(x = "Methylation (%)", y = "Methylation (%)") +
            mtclosed
        mysaveandstore(sprintf("%s/region_meth_pairs_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 4)

        # the split a 50% cut would actually produce, one point per sample
        p <- region_meth_per_sample %>%
            group_by(sample, condition, region) %>%
            summarise(frac_below_50 = mean(region_meth <= 70), .groups = "drop") %>%
            ggplot(aes(x = condition, y = frac_below_50)) +
            geom_violin(aes(fill = condition), color = "black", show.legend = FALSE) +
            geom_quasirandom(width = 0.15, size = 1) +
            stat_summary(fun.data = mean_se, geom = "pointrange", color = "black", size = 0.3) +
            region_facet +
            scale_conditions +
            labs(x = "", y = "Fraction of elements at or below 50%") +
            mtclosed
        mysaveandstore(sprintf("%s/region_meth_frac_below70_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 8, h = 3.5)

        # each element is called methylated (M) or not (U) over the three regions, and the three
        # calls concatenated give eight groups, ordered UUU through MMM
        meth_group_levels <- expand_grid(r1 = c("U", "M"), r2 = c("U", "M"), r3 = c("U", "M")) %$%
            paste0(r1, r2, r3)
        assign_meth_group <- function(region_df) {
            region_df %>%
                mutate(state = ifelse(region_meth > 70, "M", "U")) %>%
                dplyr::select(-region_meth) %>%
                pivot_wider(names_from = region, values_from = state, names_expand = TRUE) %>%
                # nothing covered in one of the three regions leaves the element uncallable
                mutate(meth_group = factor(paste0(r1, r2, r3), levels = meth_group_levels)) %>%
                filter(!is.na(meth_group))
        }

        # averaged over samples, so one group per element, used here to split the heatmap rows
        element_meth_groups <- assign_meth_group(region_meth_per_element)
        hm_all <- pad_heatmap_elements(
            merged %>% left_join(element_meth_groups %>% dplyr::select(gene_id, meth_group), by = "gene_id")
        )

        hms <- list()
        for (cur_sample in conf$samples) {
            # elements with no CpG covered in one of the three regions cannot be called, but they
            # stay in as their own group so the row set matches the other heatmaps and composes
            hm_data <- hm_all %>%
                filter(sample == cur_sample) %>%
                mutate(meth_group = replace_na(as.character(meth_group), "uncalled"))
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
                annotation_group(meth_group) %>%
                as_ComplexHeatmap()
            hms[[cur_sample]] <- p
            mysaveandstore(sprintf("%s/%s_methylation_%s_split_bymethgroup.pdf", outputdir_meth_clustering, cur_sample, subfam), w = 6, h = 9)
        }
        p <- base::eval(base::parse(text = paste0("hms[['", names(hms), "']]", collapse = " + ")))
        mysaveandstore(sprintf("%s/%s_methylation_%s_split_bymethgroup.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 9)
        rm(p)

        # the same call made independently within each sample, so an element can change group
        # between samples
        sample_meth_groups <- assign_meth_group(region_meth_per_sample)
        group_counts_per_sample <- sample_meth_groups %>%
            count(sample, condition, meth_group) %>%
            complete(nesting(sample, condition), meth_group, fill = list(n = 0)) %>%
            group_by(sample) %>%
            mutate(frac = n / sum(n)) %>%
            ungroup() %>%
            mutate(sample = factor(sample, levels = sample_order))

        p <- group_counts_per_sample %>%
            ggplot(aes(x = sample, y = n, fill = meth_group)) +
            geom_col(color = "black", linewidth = 0.2) +
            facet_wrap(~condition, scales = "free_x", nrow = 1) +
            labs(x = "", y = "Elements", fill = "Region state") +
            mtclosed +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        mysaveandstore(sprintf("%s/methgroup_counts_per_sample_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 4)

        # samples differ in how many elements are detected at all, so the same composition as a
        # fraction is what actually compares across conditions
        p <- group_counts_per_sample %>%
            ggplot(aes(x = sample, y = frac, fill = meth_group)) +
            geom_col(color = "black", linewidth = 0.2) +
            facet_wrap(~condition, scales = "free_x", nrow = 1) +
            labs(x = "", y = "Fraction of elements", fill = "Region state") +
            mtclosed +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        mysaveandstore(sprintf("%s/methgroup_fraction_per_sample_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 10, h = 4)

        # per condition, the mean over its samples rather than the total, so conditions with more
        # samples are not simply taller. one point per sample behind each bar
        p <- group_counts_per_sample %>%
            group_by(condition, meth_group) %>%
            summarise(elements_per_sample = mean(n), .groups = "drop") %>%
            ggplot(aes(x = meth_group, y = elements_per_sample, fill = condition)) +
            geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.75) +
            geom_quasirandom(
                data = group_counts_per_sample, aes(x = meth_group, y = n, group = condition),
                dodge.width = 0.8, width = 0.12, size = 0.7, color = "black", inherit.aes = FALSE
            ) +
            scale_conditions +
            labs(x = "Region state", y = "Elements per sample") +
            mtclosed
        mysaveandstore(sprintf("%s/methgroup_elements_per_sample_by_condition_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 8, h = 4)

        # the same restricted to groups whose first region (1-329) is unmethylated
        u_counts <- group_counts_per_sample %>%
            filter(str_starts(as.character(meth_group), "U")) %>%
            mutate(meth_group = fct_drop(meth_group))
        p <- u_counts %>%
            group_by(condition, meth_group) %>%
            summarise(elements_per_sample = mean(n), .groups = "drop") %>%
            ggplot(aes(x = meth_group, y = elements_per_sample, fill = condition)) +
            geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.75) +
            geom_quasirandom(
                data = u_counts, aes(x = meth_group, y = n, group = condition),
                dodge.width = 0.8, width = 0.12, size = 0.7, color = "black", inherit.aes = FALSE
            ) +
            scale_conditions +
            labs(x = "Region state", y = "Elements per sample") +
            mtclosed
        mysaveandstore(sprintf("%s/methgroup_elements_per_sample_by_condition_r1unmeth_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 6, h = 4)

        # all four of those groups collapsed, so one count per sample of elements whose 1-329 is
        # unmethylated whatever the other two regions do
        u_collapsed <- group_counts_per_sample %>%
            filter(str_starts(as.character(meth_group), "U")) %>%
            group_by(sample, condition) %>%
            summarise(n = sum(n), .groups = "drop")
        # the count per sample against condition, adjusted the same way as the per element models
        u_lm_df <- u_collapsed %>%
            mutate(sample = as.character(sample)) %>%
            dplyr::select(-condition) %>%
            left_join(sample_table %>% mutate(sample = as.character(sample)), by = "sample")
        u_fit <- lm(reformulate(lm_right_hand_side, "n"), u_lm_df)
        u_stats <- emmeans(u_fit, ~condition) %>%
            contrast(method = "pairwise") %>%
            as.data.frame() %>%
            mutate(formula = paste("n ~", lm_right_hand_side), n_samples = nobs(u_fit), .before = 1)
        print(u_stats)

        p <- u_collapsed %>%
            group_by(condition) %>%
            summarise(elements_per_sample = mean(n), .groups = "drop") %>%
            ggplot(aes(x = condition, y = elements_per_sample, fill = condition)) +
            geom_col(color = "black", width = 0.75, show.legend = FALSE) +
            geom_quasirandom(
                data = u_collapsed, aes(x = condition, y = n),
                width = 0.12, size = 0.9, color = "black", inherit.aes = FALSE
            ) +
            scale_conditions +
            labs(x = "", y = "Elements per sample with 1-329 unmethylated") +
            mtclosed
        mysaveandstore(sprintf("%s/methgroup_r1unmeth_collapsed_by_condition_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 4, h = 4, sf = u_stats)

        # which samples the unmethylated calls are shared between. one set per sample, and only
        # the largest intersections, since 22 sets have far more combinations than are worth drawing
        # split rather than purrr, because library(mclust) above masks map with its own
        u_sets <- sample_meth_groups %>%
            filter(str_starts(as.character(meth_group), "U")) %>%
            distinct(sample, gene_id) %>%
            mutate(sample = as.character(sample)) %$%
            split(gene_id, sample)
        if (length(u_sets) > 1) {
            comb <- make_comb_mat(u_sets)
            comb <- comb[head(order(comb_size(comb), decreasing = TRUE), 30)]
            p <- UpSet(comb,
                set_order = intersect(as.character(sample_order), set_name(comb)),
                comb_order = order(comb_size(comb), decreasing = TRUE)
            )
            mysaveandstore(sprintf("%s/methgroup_r1unmeth_upset_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 12, h = 6)
        }

        # the same sharing as a spectrum: an element called unmethylated in one sample only sits at
        # the left, one called in every sample at the right
        u_recurrence <- sample_meth_groups %>%
            filter(str_starts(as.character(meth_group), "U")) %>%
            count(gene_id, name = "n_samples_u")
        p <- u_recurrence %>%
            ggplot(aes(x = n_samples_u)) +
            geom_histogram(binwidth = 1, color = "black", fill = "grey70") +
            labs(x = "Samples calling the element unmethylated over 1-329", y = "Elements") +
            mtclosed
        mysaveandstore(sprintf("%s/methgroup_r1unmeth_recurrence_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 5, h = 4)

        # and split by condition, as a fraction of that condition's samples, which separates
        # elements unmethylated everywhere from ones only unmethylated in one group
        p <- sample_meth_groups %>%
            filter(str_starts(as.character(meth_group), "U")) %>%
            count(condition, gene_id) %>%
            left_join(count(sample_table, condition, name = "n_condition_samples"), by = "condition") %>%
            mutate(frac_samples_u = n / n_condition_samples) %>%
            complete(gene_id, condition, fill = list(frac_samples_u = 0)) %>%
            ggplot(aes(x = frac_samples_u, fill = condition)) +
            geom_histogram(binwidth = 0.1, color = "black", linewidth = 0.2) +
            facet_wrap(~condition, nrow = 1) +
            scale_conditions +
            labs(x = "Fraction of the condition's samples calling it unmethylated", y = "Elements") +
            mtclosed
        mysaveandstore(sprintf("%s/methgroup_r1unmeth_recurrence_by_condition_%s.pdf", outputdir_meth_clustering, subfam), pl = p, w = 9, h = 3.5)

        library(viridis)


        pf <- merged %>%
            group_by(loc_superlowres_integrative_stranded, consensus_pos) %>%
            summarise(pctM = mean(pctM)) %>%
            mutate(pos_num = as.numeric(as.character(consensus_pos)))
        pf909 <- merged %>%
            group_by(loc_superlowres_integrative_stranded, consensus_pos) %>%
            summarise(pctM = mean(pctM)) %>%
            mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
            filter(pos_num < 1000)
        pabs <- ggplot(pf, aes(x = pos_num, y = pctM, color = loc_superlowres_integrative_stranded)) +
            geom_ribbon(
                data = pf %>% group_by(pos_num) %>%
                    summarise(ymin = min(pctM), ymax = max(pctM)),
                aes(x = pos_num, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                alpha = 0.3,
                fill = "grey70"
            ) +
            geom_line(
                aes(color = loc_superlowres_integrative_stranded),
                linewidth = 0.75
            ) +
            mtclosed
        mysaveandstore(sprintf("%s/signal_wpoint_byloc_%s1.pdf", outputdir_meth_clustering, subfam), pl = pabs, w = 10, h = 2)

        pabs909 <- ggplot(pf909, aes(x = pos_num, y = pctM, color = loc_superlowres_integrative_stranded)) +
            geom_ribbon(
                data = pf909 %>% group_by(pos_num) %>%
                    summarise(ymin = min(pctM), ymax = max(pctM)),
                aes(x = pos_num, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                alpha = 0.3,
                fill = "grey70"
            ) +
            geom_line(
                aes(color = loc_superlowres_integrative_stranded),
                linewidth = 0.75
            ) +
            mtclosed
        mss(sprintf("%s/signal_909wpoint_byloc_%s.pdf", outputdir_meth_clustering, subfam), pl = pabs909, fw = 1.5, fh = 0.5, plus_void = TRUE)


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


        # with genic loc clustering
        pf <- merged %>%
            group_by(condition, consensus_pos, loc_superlowres_integrative_stranded) %>%
            summarise(pctM = mean(pctM)) %>%
            mutate(pos_num = as.numeric(as.character(consensus_pos)))
        pf909 <- merged %>%
            group_by(condition, consensus_pos, loc_superlowres_integrative_stranded) %>%
            summarise(pctM = mean(pctM)) %>%
            mutate(pos_num = as.numeric(as.character(consensus_pos))) %>%
            filter(pos_num < 1000)
        pabs_clust <- ggplot(pf, aes(x = pos_num, y = pctM, color = condition)) +
            geom_ribbon(
                data = pf %>% group_by(pos_num, loc_superlowres_integrative_stranded) %>%
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
            facet_wrap(~loc_superlowres_integrative_stranded, nrow = 1) +
            scale_conditions +
            mtclosed
        mysaveandstore(sprintf("%s/signal_wpoint_%s_clustered.pdf", outputdir_meth_clustering, subfam), pl = pabs_clust, w = 12, h = 2)
        pabs909_clust <- ggplot(pf909, aes(x = pos_num, y = pctM, color = condition)) +
            geom_ribbon(
                data = pf909 %>% group_by(pos_num, loc_superlowres_integrative_stranded) %>%
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
            facet_wrap(~loc_superlowres_integrative_stranded, nrow = 1) +
            scale_conditions +
            mtclosed
        mysaveandstore(sprintf("%s/signal_909wpoint_%s_clustered.pdf", outputdir_meth_clustering, subfam), pl = pabs909_clust, w = 20, h = 2)

        pf_wide_clust <- merged %>%
            group_by(condition, consensus_pos, loc_superlowres_integrative_stranded) %>%
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
                facet_wrap(~loc_superlowres_integrative_stranded, nrow = 1) +
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
                facet_wrap(~loc_superlowres_integrative_stranded, nrow = 1) +
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
} # end analysis_grid loop








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
mysaveandstore(sprintf("%s/nonref_%s_methylation1_%s.pdf", outputdir_meth_clustering, "all", subfam), w = 7, h = 10)


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
                seqval <- dfs %>%
                    filter(consensus_pos == start_pos) %$% sequence_pos %>%
                    pluck(1)
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
            new_width <- mapping %>%
                filter(gene_id == mcols(x)$gene_id) %$% filter_pos %>%
                pluck(1)
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
