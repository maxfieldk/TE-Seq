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


rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)
rtedf$sample <- factor(rtedf$sample, levels = conf$samples)
rtedf$condition <- factor(rtedf$condition, levels = conf$levels)
flRTEpromoter <- read_delim(sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)
RMdf <- read_delim(sprintf("ldna/Rintermediates/%s/RMdf.tsv", params$mod_code), col_names = TRUE)

perelementdf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf_promoters.tsv", params$mod_code), col_names = TRUE)
perelementdf_promoters$sample <- factor(perelementdf_promoters$sample, levels = conf$samples)
perelementdf_promoters$condition <- factor(perelementdf_promoters$condition, levels = conf$levels)


dmrs_per_contrast <- list()
dmls_per_contrast <- list()
dmrsgr_per_contrast <- list()
dmlsgr_per_contrast <- list()
dmrsannot_per_contrast <- list()
dmrsgr_split_per_contrast <- list()
for (contrast in contrasts) {
    dmr_path <- sprintf("ldna/results/%s/tables/%s/dmrs.tsv", params$mod_code, contrast)
    dml_path <- sprintf("ldna/results/%s/tables/%s/dmls.tsv", params$mod_code, contrast)
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
###########################


#################
{
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

        for (contrast in contrasts) {
            cp <- parse_contrast(contrast)
            condition1 <- cp$condition1
            condition2 <- cp$condition2

            t05_col <- paste0("t05_", contrast)
            dm_intact_l1hs_elements <- flRTEpromoter %>%
                filter(rte_subfamily == "L1HS") %>%
                filter(intactness_req == "Intact") %>%
                filter(!!sym(t05_col) == "Hypo")
            dm_fl_l1hs_elements <- flRTEpromoter %>%
                filter(rte_subfamily == "L1HS") %>%
                filter(!!sym(t05_col) == "Hypo")
            # Recompute top_l1hs_movers for this contrast
            pfl1 <- perelementdf_promoters %>%
                filter(grepl("^L1", rte_subfamily))
            top_l1hs_movers_contrast <- pfl1 %>%
                filter(condition %in% c(condition1, condition2)) %>%
                group_by(gene_id, rte_subfamily, condition) %>%
                summarize(mean_meth = mean(mean_meth), .groups = "drop") %>%
                pivot_wider(names_from = condition, values_from = mean_meth) %>%
                mutate(dif = !!sym(condition1) - !!sym(condition2)) %>%
                mutate(abs_dif = abs(dif)) %>%
                arrange(-abs_dif) %>%
                group_by(rte_subfamily) %>%
                mutate(rank_change = row_number()) %>%
                ungroup() %>%
                filter(rte_subfamily == "L1HS") %$% gene_id %>%
                head(n = 10)
            topmovers_l1hs_elements <- flRTEpromoter %>%
                filter(gene_id %in% top_l1hs_movers_contrast)
            allfl_l1hs_elements <- flRTEpromoter %>%
                filter(rte_subfamily == "L1HS")
            element_sets_of_interst <- list("dm_fl_l1hs" = dm_fl_l1hs_elements, "dm_intact_l1hs" = dm_intact_l1hs_elements, "Top_Movers" = topmovers_l1hs_elements, "all_elements" = allfl_l1hs_elements)

            for (element_type in names(element_sets_of_interst)) {
                df <- element_sets_of_interst[[element_type]]
                dir.create(sprintf("ldna/Rintermediates/%s/%s/l1hs/", params$mod_code, contrast), recursive = TRUE)
                write_delim(df %>% dplyr::select(gene_id), sprintf("ldna/Rintermediates/%s/%s/l1hs/%s_gene_id.tsv", params$mod_code, contrast, element_type), col_names = FALSE)
                write_delim(df %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/%s/%s/l1hs/%s_promoters.bed", params$mod_code, contrast, element_type), col_names = FALSE, delim = "\t")
                write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ] %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/%s/%s/l1hs/%s_full_elements.bed", params$mod_code, contrast, element_type), col_names = FALSE, delim = "\t")
                write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ], sprintf("ldna/Rintermediates/%s/%s/l1hs/%s_full_elements.tsv", params$mod_code, contrast, element_type), col_names = TRUE, delim = "\t")

                for (element in df$gene_id) {
                    browser()
                    y_lim_lower <- 50
                    y_lim_upper <- 100
                    y_valmin <- y_lim_lower
                    y_valmax <- y_lim_lower + ((y_lim_upper - y_lim_lower) / 10)

                    if (rmannextended %>% filter(gene_id == element) %$% strand == "+") {
                        modifier <- rmannextended %>% filter(gene_id == element) %$% start
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

                        p1line <- pf %>%
                            filter(gene_id == element) %>%
                            mutate(start = start - offset) %>%
                            ggplot() +
                            geom_line(aes(x = start, y = rM, color = sample)) +
                            geom_line(data = . %>%
                                group_by(start, condition) %>%
                                summarise(rM = mean(rM)), aes(x = start, y = rM), color = "black") +
                            scale_samples_unique +
                            labs(y = "Methylation Rolling Mean") +
                            mtclosed +
                            theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
                        mysaveandstore(pl = p1line, sprintf("ldna/results/%s/plots/rte/%s/%s/%s_methylation_line.pdf", params$mod_code, contrast, element_type, element), 5, 5)
                    } else {
                        modifier <- rmannextended %>% filter(gene_id == element) %$% end
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

                        p1line <- pf %>%
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
                        mysaveandstore(pl = p1line, sprintf("ldna/results/%s/plots/rte/%s/%s/%s_methylation_line.pdf", params$mod_code, contrast, element_type, element), 5, 5)
                    }

                    p <- p1 + plot_layout(heights = c(1))

                    mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/%s/%s_methylation.pdf", params$mod_code, contrast, element_type, element), 5, 5)
                    mysaveandstore(pl = p1 + ggtitle(element) + mtclosed, sprintf("ldna/results/%s/plots/rte/%s/%s/%s_methylation_nc.pdf", params$mod_code, contrast, element_type, element), 5, 4)
                }

                for (element in df$gene_id) {
                    y_lim_lower <- 50
                    y_lim_upper <- 100
                    y_valmin <- y_lim_lower
                    y_valmax <- y_lim_lower + ((y_lim_upper - y_lim_lower) / 10)

                    if (rmannextended %>% filter(gene_id == element) %$% strand == "+") {
                        modifier <- rmannextended %>% filter(gene_id == element) %$% start
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
                        modifier <- rmannextended %>% filter(gene_id == element) %$% end
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

                    mysaveandstore(sprintf("ldna/results/%s/plots/rte/%s/%s/%s_methylation_conditionaveraged.pdf", params$mod_code, contrast, element_type, element), 5, 5)
                }
            }
        } # end contrast loop for element-level plots
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
}

{
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
        conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
        topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))

        for (contrast in contrasts) {
            t05_col <- paste0("t05_", contrast)

            l1hsintactdf <- l1hsintactmethdf %>%
                group_by(gene_id) %>%
                summarise(dm_direction = dplyr::first(!!sym(t05_col)), genic_loc = dplyr::first(genic_loc))
            dm_status <- l1hsintactdf %>%
                arrange(gene_id) %>%
                filter(gene_id %in% rownames(m)) %$% dm_direction
            is_sig <- !is.na(dm_status)
            pch <- rep("*", length(dm_status))
            pch[!is_sig] <- NA
            genic_locs <- l1hsintactdf %>%
                arrange(gene_id) %>%
                filter(gene_id %in% rownames(m)) %$% genic_loc
            row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), genic_loc = genic_locs, col = list(genic_loc = c("Genic" = "brown", "Intergenic" = "tan")))

            col_fun <- colorRamp2(c(50, 75, 100), c("red", "white", "blue"))
            heatmapL1UTR <- m %>%
                Heatmap(
                    name = "CpG Methylation",
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    column_names_rot = 45,
                    col = col_fun,
                    split = dm_status,
                    top_annotation = topAnn,
                    right_annotation = row_ha,
                    row_title = "Intact L1HS"
                )

            col_fun2 <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
            heatmapL1UTR2 <- m %>%
                Heatmap(
                    name = "CpG Methylation",
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    column_names_rot = 45,
                    split = dm_status,
                    col = col_fun2,
                    top_annotation = topAnn,
                    right_annotation = row_ha,
                    row_title = "Intact L1HS"
                )

            p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR, heatmap_legend_side = "right", annotation_legend_side = "right")))
            mysaveandstore(sprintf("ldna/results/%s/plots/%s/l1intactheatmap_5utr.pdf", params$mod_code, contrast), 7, 14)

            p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR2, heatmap_legend_side = "right", annotation_legend_side = "right")))
            mysaveandstore(sprintf("ldna/results/%s/plots/%s/l1intactheatmap_5utr_fullrange.pdf", params$mod_code, contrast), 7, 14)
        }
    } else {
        col_fun <- colorRamp2(c(50, 75, 100), c("red", "white", "blue"))
        heatmapL1UTR <- m %>%
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

        col_fun2 <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
        heatmapL1UTR2 <- m %>%
            Heatmap(
                name = "CpG Methylation",
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_names_rot = 45,
                col = col_fun2,
                row_title = "Intact L1HS"
            )

        p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR, heatmap_legend_side = "right", annotation_legend_side = "right")))
        mysaveandstore(sprintf("ldna/results/%s/plots/l1intactheatmap_5utr.pdf", params$mod_code), 7, 14)

        p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR2, heatmap_legend_side = "right", annotation_legend_side = "right")))
        mysaveandstore(sprintf("ldna/results/%s/plots/l1intactheatmap_5utr_fullrange.pdf", params$mod_code), 7, 14)
    }
}


file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "15", params$mod_code))



{
    # Per-element methylation: full-length L1HS, raw CpG methylation per sample
    # X-axis in element-relative coordinates (0 = 5' end, ~6000 = 3' end)
    # Each sample gets its own points + connecting line; no rolling mean

    fl_l1hs_df <- rtedf %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(rte_length_req == "FL") %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        mutate(
            elem_pos = if_else(rte_strand == "+", start - rte_start, rte_end - start)
        ) %>%
        filter(elem_pos >= 0, elem_pos <= 6500)

    fl_l1hs_elements <- fl_l1hs_df %>%
        pull(gene_id) %>%
        unique() %>%
        sort()

    outdir <- sprintf("ldna/results/%s/plots/rte/l1hs_fl_persite", params$mod_code)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    for (element in fl_l1hs_elements) {
        tryCatch(
            {
                pf <- fl_l1hs_df %>% filter(gene_id == element)
                if (nrow(pf) < 5) next

                # --- Heatmap: CpG index x-axis, magma palette scaled to data range ---
                cpg_index_map <- pf %>%
                    distinct(elem_pos) %>%
                    arrange(elem_pos) %>%
                    mutate(cpg_index = row_number())
                last_utr_cpg_idx <- cpg_index_map %>%
                    filter(elem_pos <= 909) %>%
                    pull(cpg_index) %>%
                    max(na.rm = TRUE)
                obs_min <- min(pf$pctM, na.rm = TRUE)
                obs_max <- max(pf$pctM, na.rm = TRUE)
                BIN_SIZE <- 100
                p_heatmap2 <- pf %>%
                    left_join(cpg_index_map, by = "elem_pos") %>%
                    mutate(
                        sample = factor(sample, levels = rev(
                            sample_table %>%
                                arrange(factor(condition, levels = conf$levels)) %>%
                                pull(sample_name)
                        )),
                        condition = factor(condition, levels = rev(conf$levels))
                    ) %>%
                    ggplot(aes(x = cpg_index, y = sample, fill = pctM)) +
                    geom_tile() +
                    geom_vline(xintercept = last_utr_cpg_idx, linetype = "dashed", color = "red", linewidth = 0.5) +
                    facet_grid(condition ~ ., scales = "free_y", space = "free_y") +
                    scale_fill_viridis_c(
                        option = "magma", direction = -1,
                        limits = c(obs_min, obs_max),
                        na.value = "grey85",
                        name = "Meth (%)"
                    ) +
                    scale_x_continuous(expand = expansion(mult = 0)) +
                    labs(x = "CpG index (5'→3')", y = NULL) +
                    mtclosed +
                    theme(
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.key.height = unit(0.5, "cm"),
                        legend.key.width = unit(0.3, "cm"),
                        strip.text.y = element_text(size = 11, angle = 0, hjust = 0),
                        panel.spacing.y = unit(0.5, "mm")
                    )
                # --- Spaghetti: per-sample lines colored by condition, bold condition mean on top ---
                p_spaghetti <- pf %>%
                    mutate(
                        bin = floor(elem_pos / BIN_SIZE) * BIN_SIZE + BIN_SIZE / 2,
                        condition = factor(condition, levels = conf$levels)
                    ) %>%
                    group_by(sample, condition, bin) %>%
                    summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
                    ggplot(aes(x = bin, y = mean_meth, color = condition, group = sample)) +
                    # Individual sample lines: thin, transparent
                    geom_line(alpha = 0.35, linewidth = 0.35) +
                    # Bold condition mean on top, grouped explicitly
                    stat_summary(
                        aes(group = condition),
                        fun = mean, geom = "line",
                        linewidth = 1.4, alpha = 0.95
                    ) +
                    geom_vline(xintercept = 909, linetype = "dashed", color = "grey40", linewidth = 0.4) +
                    annotate("text", x = 909, y = 103, label = "5'UTR", size = 2.5, color = "grey40", hjust = -0.1) +
                    scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100)) +
                    scale_x_continuous(
                        limits = c(0, 6000),
                        breaks = c(0, 909, 2000, 4000, 6000),
                        expand = expansion(mult = c(0.01, 0.02))
                    ) +
                    scale_conditions +
                    labs(
                        x = "Position in element (bp, 5'→3')",
                        y = "CpG Methylation (%)"
                    ) +
                    mtclosed

                p_combined3 <- p_heatmap2 / p_spaghetti +
                    plot_layout(heights = c(0.5, 1)) +
                    plot_annotation(
                        title = element,
                        theme = theme(plot.title = element_text(size = 11, face = "italic", hjust = 0.5))
                    )
                mysaveandstore(
                    pl = p_combined3,
                    sprintf("%s/%s_combined3.pdf", outdir, element),
                    7, 6
                )
            },
            error = function(e) {
                message(sprintf("Skipping %s: %s", element, conditionMessage(e)))
            }
        )
    }

    # --- All FL L1HS elements: spaghetti faceted by gene_id ---
    n_elements <- fl_l1hs_df %>%
        pull(gene_id) %>%
        n_distinct()
    ncols <- 10
    nrows <- ceiling(n_elements / ncols)

    p_all_spaghetti <- fl_l1hs_df %>%
        mutate(
            bin = floor(elem_pos / BIN_SIZE) * BIN_SIZE + BIN_SIZE / 2,
            condition = factor(condition, levels = conf$levels)
        ) %>%
        group_by(gene_id, sample, condition, bin) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
        ggplot(aes(x = bin, y = mean_meth, color = condition, group = sample)) +
        geom_line(alpha = 0.35, linewidth = 0.35) +
        stat_summary(
            aes(group = condition),
            fun = mean, geom = "line",
            linewidth = 1.1, alpha = 0.95
        ) +
        geom_vline(xintercept = 909, linetype = "dashed", color = "grey40", linewidth = 0.3) +
        facet_wrap(~gene_id, ncol = ncols) +
        scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
        scale_x_continuous(
            limits = c(0, 6000),
            breaks = c(0, 3000, 6000),
            labels = c("0", "3k", "6k")
        ) +
        scale_conditions +
        labs(
            x = "Position in element (bp, 5'→3')",
            y = "CpG Methylation (%)"
        ) +
        mtclosed +
        theme(
            strip.text = element_text(size = 5.5, face = "italic"),
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            panel.spacing = unit(1, "mm"),
            legend.position = "bottom"
        )

    png(
        sprintf("%s/all_elements_spaghetti.png", outdir),
        width = min(ncols * 2.5, 36), height = min(nrows * 1.8, 50),
        units = "in", res = 175
    )
    print(p_all_spaghetti)
    dev.off()

    mysaveandstore(
        pl = p_all_spaghetti,
        sprintf("%s/all_elements_spaghetti_highres.pdf", outdir),
        min(10), min(20),
        raster = TRUE
    )

    # --- Heatmap: mean of first 8 CpGs per element, rows = samples, columns = elements ---
    first8_avg <- fl_l1hs_df %>%
        group_by(gene_id, sample, condition) %>%
        arrange(elem_pos, .by_group = TRUE) %>%
        slice_head(n = 8) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop")

    # Order elements by mean methylation across all samples (low → high)
    element_order <- first8_avg %>%
        group_by(gene_id) %>%
        summarise(overall_mean = mean(mean_meth, na.rm = TRUE), .groups = "drop") %>%
        arrange(overall_mean) %>%
        pull(gene_id)

    sample_order <- sample_table %>%
        arrange(factor(condition, levels = conf$levels)) %>%
        pull(sample_name)

    p_first8_heatmap <- first8_avg %>%
        mutate(
            gene_id = factor(gene_id, levels = rev(element_order)),
            sample = factor(sample, levels = sample_order),
            condition = factor(condition, levels = conf$levels)
        ) %>%
        ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
        geom_tile() +
        facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
        scale_fill_viridis_c(
            option = "magma", direction = -1,
            limits = c(0, 100), na.value = "grey85",
            name = "Meth (%)"
        ) +
        labs(
            x = NULL, y = "Element (ordered by mean methylation)",
            title = "Mean methylation: first 8 CpGs per element"
        ) +
        mtclosed +
        theme(
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 5, face = "italic"),
            strip.text.x = element_text(size = 11),
            panel.spacing.x = unit(0.5, "mm"),
            plot.title = element_text(size = 9, hjust = 0.5)
        )

    n_el <- length(element_order)
    n_samp <- length(sample_order)
    mysaveandstore(
        pl = p_first8_heatmap,
        fn = sprintf("%s/all_elements_first8cpg_heatmap.pdf", outdir),
        w = max(4, n_samp * 0.3), h = max(6, n_el * 0.15)
    )

    # --- Two-column version for slides: split elements into two side-by-side panels ---
    half <- ceiling(n_el / 2)
    el_left <- element_order[(half + 1):n_el] # higher methylation
    el_right <- element_order[1:half] # lower methylation

    make_half_hm <- function(els) {
        first8_avg %>%
            filter(gene_id %in% els) %>%
            mutate(
                gene_id = factor(gene_id, levels = rev(els)),
                sample = factor(sample, levels = sample_order),
                condition = factor(condition, levels = conf$levels)
            ) %>%
            ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
            geom_tile() +
            facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
            scale_fill_viridis_c(
                option = "magma", direction = -1,
                limits = c(0, 100), na.value = "grey85",
                name = "Meth (%)"
            ) +
            labs(x = NULL, y = NULL) +
            mtclosed +
            theme(
                axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 5, face = "italic"),
                strip.text.x = element_text(size = 11),
                panel.spacing.x = unit(0.5, "mm")
            )
    }

    p_first8_heatmap_2col <- (make_half_hm(el_left) | make_half_hm(el_right)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = "Mean methylation: first 8 CpGs per element",
            theme = theme(plot.title = element_text(size = 9, hjust = 0.5))
        )

    mysaveandstore(
        pl = p_first8_heatmap_2col,
        fn = sprintf("%s/all_elements_first8cpg_heatmap_2col.pdf", outdir),
        w = max(8, n_samp * 0.6), h = max(4, (n_el / 2) * 0.15)
    )

    # --- Heatmap: mean of first 25 CpGs per element ---
    first25_avg <- fl_l1hs_df %>%
        group_by(gene_id, sample, condition) %>%
        arrange(elem_pos, .by_group = TRUE) %>%
        slice_head(n = 25) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop")

    element_order_25 <- first25_avg %>%
        group_by(gene_id) %>%
        summarise(overall_mean = mean(mean_meth, na.rm = TRUE), .groups = "drop") %>%
        arrange(overall_mean) %>%
        pull(gene_id)

    p_first25_heatmap <- first25_avg %>%
        mutate(
            gene_id = factor(gene_id, levels = rev(element_order_25)),
            sample = factor(sample, levels = sample_order),
            condition = factor(condition, levels = conf$levels)
        ) %>%
        ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
        geom_tile() +
        facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
        scale_fill_viridis_c(
            option = "magma", direction = -1,
            limits = c(0, 100), na.value = "grey85",
            name = "Meth (%)"
        ) +
        labs(
            x = NULL, y = "Element (ordered by mean methylation)",
            title = "Mean methylation: first 25 CpGs per element"
        ) +
        mtclosed +
        theme(
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 5, face = "italic"),
            strip.text.x = element_text(size = 11),
            panel.spacing.x = unit(0.5, "mm"),
            plot.title = element_text(size = 9, hjust = 0.5)
        )

    n_el_25 <- length(element_order_25)
    mysaveandstore(
        pl = p_first25_heatmap,
        fn = sprintf("%s/all_elements_first25cpg_heatmap.pdf", outdir),
        w = max(4, n_samp * 0.3), h = max(6, n_el_25 * 0.15)
    )

    # Two-column version: first 25 CpGs
    half_25 <- ceiling(n_el_25 / 2)
    el_left_25 <- element_order_25[(half_25 + 1):n_el_25]
    el_right_25 <- element_order_25[1:half_25]

    make_half_hm_25 <- function(els) {
        first25_avg %>%
            filter(gene_id %in% els) %>%
            mutate(
                gene_id = factor(gene_id, levels = rev(els)),
                sample = factor(sample, levels = sample_order),
                condition = factor(condition, levels = conf$levels)
            ) %>%
            ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
            geom_tile() +
            facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
            scale_fill_viridis_c(
                option = "magma", direction = -1,
                limits = c(0, 100), na.value = "grey85",
                name = "Meth (%)"
            ) +
            labs(x = NULL, y = NULL) +
            mtclosed +
            theme(
                axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 5, face = "italic"),
                strip.text.x = element_text(size = 11),
                panel.spacing.x = unit(0.5, "mm")
            )
    }

    p_first25_heatmap_2col <- (make_half_hm_25(el_left_25) | make_half_hm_25(el_right_25)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = "Mean methylation: first 25 CpGs per element",
            theme = theme(plot.title = element_text(size = 9, hjust = 0.5))
        )

    mysaveandstore(
        pl = p_first25_heatmap_2col,
        fn = sprintf("%s/all_elements_first25cpg_heatmap_2col.pdf", outdir),
        w = max(8, n_samp * 0.6), h = max(4, (n_el_25 / 2) * 0.15)
    )

    # --- Heatmap: first 8 CpGs, faceted by loc_superlowres_integrative_stranded_incl_pred ---
    first8_avg_loc <- fl_l1hs_df %>%
        group_by(gene_id, sample, condition, loc_superlowres_integrative_stranded_incl_pred) %>%
        arrange(elem_pos, .by_group = TRUE) %>%
        slice_head(n = 8) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop")

    p_first8_heatmap_loc <- first8_avg_loc %>%
        mutate(
            gene_id = factor(gene_id, levels = rev(element_order)),
            sample = factor(sample, levels = sample_order),
            condition = factor(condition, levels = conf$levels)
        ) %>%
        ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
        geom_tile() +
        facet_grid(loc_superlowres_integrative_stranded_incl_pred ~ condition, scales = "free", space = "free") +
        scale_fill_viridis_c(
            option = "magma", direction = -1,
            limits = c(0, 100), na.value = "grey85",
            name = "Meth (%)"
        ) +
        labs(
            x = NULL, y = "Element (ordered by mean methylation)",
            title = "Mean methylation: first 8 CpGs per element\n(faceted by genomic location)"
        ) +
        mtclosed +
        theme(
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 5, face = "italic"),
            strip.text = element_text(size = 9),
            panel.spacing = unit(0.5, "mm"),
            plot.title = element_text(size = 9, hjust = 0.5)
        )

    n_locs <- length(unique(first8_avg_loc$loc_superlowres_integrative_stranded_incl_pred))
    mysaveandstore(
        pl = p_first8_heatmap_loc,
        fn = sprintf("%s/all_elements_first8cpg_heatmap_byloc.pdf", outdir),
        w = max(5, n_samp * 0.35), h = max(8, n_el * 0.15 * n_locs / 2)
    )

    # Two-column version: first 8 CpGs, faceted by loc
    make_half_hm_loc <- function(els, avg_df) {
        avg_df %>%
            filter(gene_id %in% els) %>%
            mutate(
                gene_id = factor(gene_id, levels = rev(els)),
                sample = factor(sample, levels = sample_order),
                condition = factor(condition, levels = conf$levels)
            ) %>%
            ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
            geom_tile() +
            facet_grid(loc_superlowres_integrative_stranded_incl_pred ~ condition, scales = "free", space = "free") +
            scale_fill_viridis_c(
                option = "magma", direction = -1,
                limits = c(0, 100), na.value = "grey85",
                name = "Meth (%)"
            ) +
            labs(x = NULL, y = NULL) +
            mtclosed +
            theme(
                axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 5, face = "italic"),
                strip.text = element_text(size = 9),
                panel.spacing = unit(0.5, "mm")
            )
    }

    p_first8_heatmap_2col_loc <- (make_half_hm_loc(el_left, first8_avg_loc) | make_half_hm_loc(el_right, first8_avg_loc)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = "Mean methylation: first 8 CpGs per element (by genomic location)",
            theme = theme(plot.title = element_text(size = 9, hjust = 0.5))
        )

    mysaveandstore(
        pl = p_first8_heatmap_2col_loc,
        fn = sprintf("%s/all_elements_first8cpg_heatmap_2col_byloc.pdf", outdir),
        w = max(8, n_samp * 0.6), h = max(6, (n_el / 2) * 0.15 * n_locs / 2)
    )

    # --- Heatmap: first 25 CpGs, faceted by loc_superlowres_integrative_stranded_incl_pred ---
    first25_avg_loc <- fl_l1hs_df %>%
        group_by(gene_id, sample, condition, loc_superlowres_integrative_stranded_incl_pred) %>%
        arrange(elem_pos, .by_group = TRUE) %>%
        slice_head(n = 25) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), .groups = "drop")

    p_first25_heatmap_loc <- first25_avg_loc %>%
        mutate(
            gene_id = factor(gene_id, levels = rev(element_order_25)),
            sample = factor(sample, levels = sample_order),
            condition = factor(condition, levels = conf$levels)
        ) %>%
        ggplot(aes(x = sample, y = gene_id, fill = mean_meth)) +
        geom_tile() +
        facet_grid(loc_superlowres_integrative_stranded_incl_pred ~ condition, scales = "free", space = "free") +
        scale_fill_viridis_c(
            option = "magma", direction = -1,
            limits = c(0, 100), na.value = "grey85",
            name = "Meth (%)"
        ) +
        labs(
            x = NULL, y = "Element (ordered by mean methylation)",
            title = "Mean methylation: first 25 CpGs per element\n(faceted by genomic location)"
        ) +
        mtclosed +
        theme(
            axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 5, face = "italic"),
            strip.text = element_text(size = 9),
            panel.spacing = unit(0.5, "mm"),
            plot.title = element_text(size = 9, hjust = 0.5)
        )

    mysaveandstore(
        pl = p_first25_heatmap_loc,
        fn = sprintf("%s/all_elements_first25cpg_heatmap_byloc.pdf", outdir),
        w = max(5, n_samp * 0.35), h = max(8, n_el_25 * 0.15 * n_locs / 2)
    )

    # Two-column version: first 25 CpGs, faceted by loc
    p_first25_heatmap_2col_loc <- (make_half_hm_loc(el_left_25, first25_avg_loc) | make_half_hm_loc(el_right_25, first25_avg_loc)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = "Mean methylation: first 25 CpGs per element (by genomic location)",
            theme = theme(plot.title = element_text(size = 9, hjust = 0.5))
        )

    mysaveandstore(
        pl = p_first25_heatmap_2col_loc,
        fn = sprintf("%s/all_elements_first25cpg_heatmap_2col_byloc.pdf", outdir),
        w = max(8, n_samp * 0.6), h = max(6, (n_el_25 / 2) * 0.15 * n_locs / 2)
    )
}
