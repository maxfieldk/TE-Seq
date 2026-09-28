module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
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
# library(msigdbr)
library(Biostrings)
library(zoo)

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
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethylpaths = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_bedMethyl_unfiltered_hardtomapregions.bed", samples, samples),
            data = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss_unfiltered_hardtomapregions.tsv", sample_table$sample_name, sample_table$sample_name)
            # dmrs = "ldna/results/m/tables/dmrs.tsv",
            # dmls = "ldna/results/m/tables/dmls.tsv",
            # read_mods = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            # read_mods_cg = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("params", list(mod_code = "m"), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)


##########################
# PREP DATA FOR ANALYSIS
sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE), col_names = FALSE)
    df_m <- df %>% filter(X4 == params$mod_code)
    rm(df)
    gr <- GRanges(
        seqnames = df_m$X1,
        ranges = IRanges(start = df_m$X2, end = df_m$X2),
        cov = df_m$X10,
        pctM = as.double(df_m$X11)
    )
    gr$sample <- sample_name
    gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_grs[[sample_name]] <- gr
}

grs <- Reduce(c, sample_grs)
rm(sample_grs)
# filter out low coverage and ensure that all samples have the same cpgs NOTE THAT WE DON'T DO THIS NOW TO RETAIN RDNA LOCI
# grs <- grsunfiltered[grsunfiltered$cov > MINIMUMCOVERAGE]
grsdf <- tibble(as.data.frame(grs))
# grsdf %$% seqnames %>% unique()
# dir.create("ldna/Rintermediates", recursive = TRUE)
# grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
# seqnames <- grsdf$seqnames
# start <- grsdf$start
# end <- grsdf$end
# pos <- paste0(seqnames, "_", start, "_", end)
# grsdf$pos <- pos

# grsdfuntidy <- grsdf %>%
#     filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
#     pivot_wider(id_cols = c("pos", "seqnames"), names_from = "sample", values_from = "pctM", names_prefix = "pctM") %>%
#     drop_na()

# grsinboth <- grsdfuntidy %>% pull(pos)
# rm(grsdfuntidy)

# grsdffiltered <- grsdf %>%
#     filter(pos %in% grsinboth)

# grsdf <- grsdffiltered
# rm(grsdffiltered)

# cpg_islands <- rtracklayer::import(conf$cpg_islands)
# cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
# cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
# cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)


# grs <- GRanges(grsdf)
# grs_cpg_islands <- grs %>% subsetByOverlaps(cpg_islands)
# grs_cpg_islands$islandStatus <- "island"
# grs_cpgi_shelves <- grs %>% subsetByOverlaps(cpgi_shelves)
# grs_cpgi_shelves$islandStatus <- "shelf"
# grs_cpgi_shores <- grs %>% subsetByOverlaps(cpgi_shores)
# grs_cpgi_shores$islandStatus <- "shore"
# grs_cpg_opensea <- grs %>% subsetByOverlaps(cpgi_features, invert = TRUE)
# grs_cpg_opensea$islandStatus <- "opensea"

# grs <- c(grs_cpg_islands, grs_cpgi_shelves, grs_cpgi_shores, grs_cpg_opensea)
# grsdf <- tibble(as.data.frame(grs))

dir.create("ldna/Rintermediates", showWarnings = FALSE)
write_delim(grsdf, sprintf("ldna/Rintermediates/%s/grsdf_unfiltered_hardtomapregions.tsv", params$mod_code), col_names = TRUE)
# grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
# grs <- GRanges(grsdf)

rmann <- get_repeat_annotations(
    default_or_extended = "extended",
    keep_non_central = FALSE
)


merge_with_grs <- function(grs, region_frame) {
    mbo <- mergeByOverlaps(grs, region_frame)
    methdf <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    region_only_frame <- mbo$region_frame %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(region_seqnames = seqnames, region_start = start, region_end = end, region_strand = strand, region_width = width)
    regiondf <- bind_cols(methdf, region_only_frame)
    return(regiondf)
}


catlift <- rtracklayer::import("/users/mkelsey/data/Nanopore/alz/RTE/catLiftOffGenesV1.bed")

rdna_transcripts <- catlift[grepl("RNA45", mcols(catlift)$name)]
catlift[grepl("RNA45", mcols(catlift)$name)] %>%
    as.data.frame() %>%
    tibble() %$% name %>%
    table()

rdna_models <- rtracklayer::import("/users/mkelsey/data/Nanopore/alz/RTE/rdnaModel.bed")


plot_region <- function(region, meth_grs, highlights = NULL, group_var = "sample", smoothing = 0, type = "point", highlight_name_var = "name") {
    pf <- meth_grs %>%
        filter(seqnames == as.character(runValue(seqnames(region)))) %>%
        filter(start > start(region)) %>%
        filter(end < end(region)) %>%
        group_by(start, end, !!sym(group_var)) %>%
        summarise(pctM = mean(pctM))

    if (smoothing != 0) {
        pf <- pf %>%
            group_by(!!sym(group_var)) %>%
            mutate(plot_val = rollmean(pctM, smoothing, na.pad = TRUE, align = "center")) %>%
            filter(!is.na(plot_val)) %>%
            ungroup()
    } else {
        pf <- pf %>% mutate(plot_val = pctM)
    }


    if (type == "point") {
        p <- pf %>% ggplot() +
            geom_point(aes(x = start, y = plot_val, color = !!sym(group_var))) +
            scale_x_continuous(breaks = scales::breaks_pretty(3)) +
            ylim(c(0, 100)) +
            mtclosed
    } else if (type == "line") {
        p <- pf %>% ggplot() +
            geom_line(aes(x = start, y = plot_val, color = !!sym(group_var))) +
            scale_x_continuous(breaks = scales::breaks_pretty(3)) +
            ylim(c(0, 100)) +
            scale_samples +
            mtclosed
    }

    if (!is.null(highlights)) {
        pf_transcripts <- highlights %>%
            filter(start > start(region), end < end(region)) %>%
            mutate(middle = (start + end) / 2)
        p <- p +
            new_scale_colour() +
            new_scale_fill() +
            geom_rect(data = pf_transcripts, aes(xmin = start, ymin = 0, xmax = end, ymax = 100, fill = !!sym(highlight_name_var)), alpha = 0.2)
    }
    return(p)
}

# region <- GRanges(seqnames = "chr13", ranges = IRanges(start = 6051123, end = 6170503))
# meth_grs <- temp
# highlights <- temp_transcripts


# p <- plot_region(region, meth_grs, highlights, group_var = "sample", smoothing = 100, type = "line")
# mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/test1111.pdf", params$mod_code), 20, 4)



# models_grs <- merge_with_grs(grs, rdna_models)
# temp <- models_grs %>% filter(seqnames == "chr13")
# temp_transcripts <- rdna_transcripts %>%
#     as.data.frame() %>%
#     tibble() %>%
#     filter(seqnames == "chr13")
# library(zoo)
# pf <- temp %>%
#     filter(start > 6051123, end < 6070503) %>%
#     group_by(sample) %>%
#     mutate(rM = rollmean(pctM, 100, na.pad = TRUE, align = "center")) %>%
#     filter(!is.na(rM)) %>%
#     ungroup()
# pf_transcripts <- temp_transcripts %>% filter(start > 6051123, end < 6070503)

# pf <- temp %>%
#     filter(end < 7000000) %>%
#     group_by(seqnames, start, end, condition) %>%
#     summarise(pctM = mean(pctM)) %>%
#     group_by(condition) %>%
#     mutate(rM = rollmean(pctM, 100, na.pad = TRUE, align = "center")) %>%
#     filter(!is.na(rM)) %>%
#     ungroup()
# pf_transcripts <- temp_transcripts %>% filter(start > 5964549, end < 7000000)
# p <- pf %>% ggplot() +
#     geom_point(aes(x = start, y = rM, color = condition)) +
#     scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#     geom_rect(data = pf_transcripts, aes(xmin = start, ymin = 0, xmax = end, ymax = 100), alpha = 0.2) +
#     ylim(c(0, 100)) +
#     mtclosed +
#     scale_conditions
# mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/test_condition11.pdf", params$mod_code), 40, 4, raster = TRUE)

library(ggnewscale)
library(ggpubr)
library(ggrepel)
library(purrr)

##########################
# rDNA: methylation aggregated over each region definition
# 45S promoter (-186/+20 around TSS), the 45S transcript itself, the intergenic
# spacer (model region minus transcripts), and the full rDNA model region.
rdna_igs <- GenomicRanges::setdiff(rdna_models, rdna_transcripts, ignore.strand = TRUE)
mcols(rdna_igs)$name <- paste0("IGS_", as.character(seqnames(rdna_igs)))

rdna_regions <- list(
    "45S promoter" = promoters(rdna_transcripts, upstream = 186, downstream = 20),
    "45S transcript" = rdna_transcripts,
    "IGS" = rdna_igs,
    "rDNA model" = rdna_models
)


grs %>% subsetByOverlaps(promoters(rdna_transcripts, upstream = 1000, downstream = 1000))

rdna_grs <- imap(rdna_regions, function(region, region_type) {
    merge_with_grs(grs, region) %>%
        dplyr::select(seqnames, start, end, cov, pctM, sample, condition, name, region_start, region_end) %>%
        mutate(region_type = region_type)
}) %>%
    bind_rows() %>%
    mutate(region_type = factor(region_type, levels = names(rdna_regions)))

# per sample x region type
rdna_sample_meth <- rdna_grs %>%
    group_by(region_type, sample) %>%
    summarise(mean_meth = mean(pctM), n_cpgs = n(), .groups = "drop") %>%
    dplyr::rename(sample_name = sample) %>%
    left_join(sample_table)
dir.create(sprintf("ldna/results/%s/tables/hardtomapregions", params$mod_code), showWarnings = FALSE, recursive = TRUE)
write_delim(rdna_sample_meth, sprintf("ldna/results/%s/tables/hardtomapregions/rdna_region_meth_by_sample.tsv", params$mod_code), delim = "\t")

# per sample x region type x chromosome (each acrocentric arm separately)
rdna_sample_chr_meth <- rdna_grs %>%
    group_by(region_type, seqnames, sample, condition) %>%
    summarise(mean_meth = mean(pctM), n_cpgs = n(), .groups = "drop")
write_delim(rdna_sample_chr_meth, sprintf("ldna/results/%s/tables/hardtomapregions/rdna_region_meth_by_sample_chr.tsv", params$mod_code), delim = "\t")

# per individual region copy (each 45S transcript / promoter / IGS / model element)
rdna_copy_meth <- rdna_grs %>%
    group_by(region_type, seqnames, region_start, region_end, name, condition) %>%
    summarise(mean_meth = mean(pctM), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_meth)
write_delim(rdna_copy_meth, sprintf("ldna/results/%s/tables/hardtomapregions/rdna_region_meth_by_copy.tsv", params$mod_code), delim = "\t")

rdna_condition_stats <- function(df, group = "region_type") {
    df %>%
        group_by(!!sym(group)) %>%
        group_modify(~ summary(lm(mean_meth ~ condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z, .x)) %>% broom::tidy()) %>%
        ungroup()
}

library(ggpubr)
# one violin per region type
rdna_sample_meth %$% region_type
for (rt in names(rdna_regions)) {
    pf <- rdna_sample_meth %>% filter(region_type == rt)
    p <- pf %>%
        ggviolin(x = "condition", y = "mean_meth", fill = "condition", add = c("dotplot")) +
        scale_conditions +
        labs(y = "Mean methylation (%)", title = sprintf("45S rDNA: %s", rt)) +
        mtclosed
    stats <- summary(lm(mean_meth ~ condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z, pf)) %>% broom::tidy()
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_%s_mean_meth_violin.pdf", params$mod_code, gsub(" ", "_", rt)), 4, 4, sf = stats)
}

# all region types side by side
p <- rdna_sample_meth %>%
    ggplot(aes(x = region_type, y = mean_meth)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    scale_conditions +
    labs(x = "", y = "Mean methylation (%)", title = "rDNA methylation by region") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_regions_mean_meth_bar.pdf", params$mod_code), 5, 4, sf = rdna_condition_stats(rdna_sample_meth))

p <- rdna_sample_meth %>%
    ggviolin(x = "condition", y = "mean_meth", fill = "condition", add = c("dotplot"), facet.by = "region_type", nrow = 1) +
    scale_conditions +
    labs(y = "Mean methylation (%)") +
    mtclosed
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_regions_mean_meth_violin.pdf", params$mod_code), 10, 4, sf = rdna_condition_stats(rdna_sample_meth))

# broken out by chromosome
p <- rdna_sample_chr_meth %>%
    ggplot(aes(x = seqnames, y = mean_meth)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    facet_wrap(~region_type, nrow = 1) +
    scale_conditions +
    labs(x = "", y = "Mean methylation (%)") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_regions_mean_meth_bar_by_chr.pdf", params$mod_code), 12, 4)

# meta-profile across all 45S copies: CpG position relative to the TSS, with 1 kb of upstream IGS
rdna_45S_flank <- resize(rdna_transcripts, width = width(rdna_transcripts) + 1000, fix = "end") # all copies are + strand
rdna_profile <- merge_with_grs(grs, rdna_45S_flank) %>%
    mutate(rel_pos = start - (region_start + 1000)) %>%
    mutate(bin = floor(rel_pos / 100) * 100)
p <- rdna_profile %>%
    group_by(condition, bin) %>%
    summarise(pctM = mean(pctM), .groups = "drop") %>%
    ggplot(aes(x = bin, y = pctM, color = condition)) +
    geom_line() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_conditions +
    ylim(c(0, 100)) +
    labs(x = "Position relative to 45S TSS (bp)", y = "Methylation (%)", title = "45S rDNA meta-profile") +
    mtclosed
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_45S_metaprofile_condition.pdf", params$mod_code), 7, 4)
p <- rdna_profile %>%
    group_by(sample, bin) %>%
    summarise(pctM = mean(pctM), .groups = "drop") %>%
    ggplot(aes(x = bin, y = pctM, color = sample)) +
    geom_line() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_samples +
    ylim(c(0, 100)) +
    labs(x = "Position relative to 45S TSS (bp)", y = "Methylation (%)", title = "45S rDNA meta-profile") +
    mtclosed
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomapregions/rdna_45S_metaprofile_sample.pdf", params$mod_code), 7, 4)


##########################
# centromeres and centromeric satellites
cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
centromere <- cytobandsdf %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
cent <- reduce(centromere)

censat <- import.bed("resources/genomes/hs1/annotations/censat.bed")
mcols(censat)$types <- gsub("_.*", "", mcols(censat)$name)

grs_cent <- merge_with_grs(grs, cent)

# per-chromosome centromere plots with satellite annotation highlighted
for (chrofinterest in c("chr3", "chr4")) {
    roi <- cent[seqnames(cent) == chrofinterest]
    highlightstemp <- censat %>%
        subsetByOverlaps(roi) %>%
        as.data.frame()
    p <- plot_region(roi, grs_cent, highlights = highlightstemp, group_var = "sample", smoothing = 500, type = "line") + ggtitle(chrofinterest)
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight.pdf", params$mod_code, chrofinterest), 9, 4)
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight_tall.pdf", params$mod_code, chrofinterest), 9, 8)
}

p <- grs_cent %>%
    group_by(sample, seqnames) %>%
    mutate(plot_val = rollmean(pctM, 500, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(plot_val)) %>%
    ungroup() %>%
    ggplot() +
    geom_line(aes(x = start, y = plot_val, color = sample)) +
    facet_wrap(~seqnames, ncol = 4, scales = "free_x") +
    ylim(c(0, 100)) +
    scale_samples +
    mtclosed +
    theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
    )
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_together.pdf", params$mod_code), 12, 12)

grs_censat <- merge_with_grs(grs, censat)
censat_sample_meth <- grs_censat %>%
    group_by(sample, condition, types) %>%
    summarise(mean_meth = mean(pctM), .groups = "drop")
p <- censat_sample_meth %>%
    ggplot(aes(x = types, y = mean_meth)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    scale_conditions +
    labs(x = "", y = "Mean methylation (%)") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/censat_bar.pdf", params$mod_code), 8, 4, sf = rdna_condition_stats(censat_sample_meth, group = "types"))

mcols(censat)$types %>% unique()


##########################
# full length RTEs inside centromeres
cent_l1hs <- rmann %>%
    filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "HERVK-int")) %>%
    filter(rte_length_req == "FL") %>%
    GRanges() %>%
    subsetByOverlaps(cent)

for (gene_id in mcols(cent_l1hs)$gene_id) {
    element <- cent_l1hs[mcols(cent_l1hs)$gene_id == gene_id]
    p <- plot_region(element, grs_cent, group_var = "sample", smoothing = 5, type = "point") +
        scale_samples
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s.pdf", params$mod_code, gene_id), 5, 4)
    p <- plot_region(element, grs_cent, group_var = "condition", smoothing = 5, type = "point") +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_condition.pdf", params$mod_code, gene_id), 5, 4)
}


##########################
# satellites
sat <- rmann %>%
    filter(repeat_superfamily == "SAT") %>%
    GRanges()
satgrs <- merge_with_grs(grs, sat)

satdf <- satgrs %>%
    group_by(family, sample, condition) %>%
    summarise(mean_meth = mean(pctM), .groups = "drop") %>%
    mutate(family = fct_reorder(family, mean_meth, mean))

p <- satdf %>% ggviolin(x = "family", y = "mean_meth", fill = "condition", add = c("dotplot")) +
    scale_conditions +
    labs(x = "", y = "Mean methylation (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/satellites_violin.pdf", params$mod_code), 10, 4)

p <- satdf %>%
    ggplot(aes(x = family, y = mean_meth)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    scale_conditions +
    labs(x = "", y = "Mean methylation (%)") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/satellites_bar.pdf", params$mod_code), 8, 4, sf = rdna_condition_stats(satdf, group = "family"))

satgene <- satgrs %>%
    group_by(gene_id, family, length, condition) %>%
    summarise(mean_meth = mean(pctM), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_meth)
write_delim(satgene, sprintf("ldna/results/%s/tables/hardtomapregions/satellite_meth_by_element.tsv", params$mod_code), delim = "\t")

satgene




aaa <- grsdf %>%
    filter(seqnames == "chr22") %>%
    filter(start > 349817) %>%
    filter(end < 402540)
aaa %>% arrange(-start)
# interesting L1PA3 region.  chr22:336,762-385,857 #chr22:349,817-382,54
L1PA3_22p13_1_region <- GRanges(seqnames = "chr22", ranges = IRanges(start = 349817, end = 382540))
p <- plot_region(L1PA3_22p13_1_region, grsdf, group_var = "condition", smoothing = 5, type = "point") +
    scale_conditions
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/L1PA3_22p13_1_region_condition.pdf", params$mod_code), 5, 4)
highlightsdf <- rmann %>% filter(gene_id == "L1PA3_22p13_1")
p <- plot_region(L1PA3_22p13_1_region, grsdf, highlights = highlightsdf, highlight_name_var = "gene_id", group_var = "sample", smoothing = 5, type = "line") +
    scale_samples
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/L1PA3_22p13_1_region_sample11111.pdf", params$mod_code), 5, 4)
