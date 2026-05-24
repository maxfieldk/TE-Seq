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

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
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

grsunfiltered <- Reduce(c, sample_grs)
rm(sample_grs)
# filter out low coverage and ensure that all samples have the same cpgs
grs <- grsunfiltered[grsunfiltered$cov > 3]
grsdf <- tibble(as.data.frame(grs))
grsdf %$% seqnames %>% unique()
dir.create("ldna/Rintermediates", recursive = TRUE)
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
seqnames <- grsdf$seqnames
start <- grsdf$start
end <- grsdf$end
pos <- paste0(seqnames, "_", start, "_", end)
grsdf$pos <- pos

grsdfuntidy <- grsdf %>%
    filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
    pivot_wider(id_cols = c("pos", "seqnames"), names_from = "sample", values_from = "pctM", names_prefix = "pctM")

grsdfuntidy1 <- grsdfuntidy %>%
    mutate(sum_na_AD = case_when(is.na(pctMSEN1) + is.na(pctMSEN2) + is.na(pctMSEN3) + is.na(pctMSEN4) > 2 ~ "FAIL", TRUE ~ "PASS")) %>%
    mutate(sum_na_CTRL = case_when(is.na(pctMPRO1) + is.na(pctMPRO2) + is.na(pctMPRO3) + is.na(pctMPRO4) > 2 ~ "FAIL", TRUE ~ "PASS")) %>%
    mutate(keep = case_when(sum_na_AD == "PASS" & sum_na_CTRL == "PASS" ~ TRUE, TRUE ~ FALSE))

grsinboth <- grsdfuntidy1 %>%
    filter(keep == TRUE) %>%
    pull(pos)
# TODO pivot wider introduces sample names as variables, this needs to be addressed
grsdffiltered <- grsdf %>%
    filter(pos %in% grsinboth)


grsdf <- grsdffiltered
rm(grsdffiltered)
rm(grsdfuntidy)

cpg_islands <- rtracklayer::import(conf$cpg_islands)
cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)


grs <- GRanges(grsdf)
grs_cpg_islands <- grs %>% subsetByOverlaps(cpg_islands)
grs_cpg_islands$islandStatus <- "island"
grs_cpgi_shelves <- grs %>% subsetByOverlaps(cpgi_shelves)
grs_cpgi_shelves$islandStatus <- "shelf"
grs_cpgi_shores <- grs %>% subsetByOverlaps(cpgi_shores)
grs_cpgi_shores$islandStatus <- "shore"
grs_cpg_opensea <- grs %>% subsetByOverlaps(cpgi_features, invert = TRUE)
grs_cpg_opensea$islandStatus <- "opensea"

grs <- c(grs_cpg_islands, grs_cpgi_shelves, grs_cpgi_shores, grs_cpg_opensea)
grsdf <- tibble(as.data.frame(grs))

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


rdna_transcripts_promoters <- promoters(rdna_transcripts, upstream = 186, downstream = 20)
rdna_transcripts_promoters_grs <- merge_with_grs(grs, rdna_transcripts_promoters)

rdna_transcripts_promoters_grs %>%
    group_by(sample, name) %>%
    summarise(pctM = mean(pctM)) %>%
    ungroup() %>%
    pivot_wider(names_from = sample, values_from = pctM)

pf <- rdna_transcripts_promoters_grs %>%
    group_by(sample) %>%
    summarise(mean_meth = mean(pctM)) %>%
    ungroup() %>%
    dplyr::rename(sample_name = sample) %>%
    left_join(sample_table) %>%
    mutate(sample_name = fct_reorder(sample_name, mean_meth))

library(ggnewscale)
library(ggpubr)
p <- pf %>%
    ggviolin(x = "condition", y = "mean_meth", fill = "condition", add = c("dotplot")) +
    scale_conditions +
    ggtitle("45S rDNA Promoter Methylation") +
    mtclosed
stats <- summary(lm(mean_meth ~ condition, pf)) %>% broom::tidy()
mysaveandstore(fn = sprintf("ldna/results/%s/plots/hardtomapregions/rdna_mean_meth_violin.pdf", params$mod_code), 4, 4, sf = stats)

library(ggrepel)
# p <- pf %>%
#     mutate(sample_name = fct_reorder(paste0(sample_name, "_", age), age)) %>%
#     ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = sex)) +
#     geom_point(size = 3) +
#     scale_conditions +
#     geom_vline(xintercept = pf %>% filter(condition == "AD") %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
#     geom_vline(xintercept = pf %>% filter(condition == "CTRL") %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
#     geom_text_repel(aes(label = apoe)) +
#     # new_scale_fill() +
#     # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
#     scale_conditions +
#     mtopen
# stats <- summary(lm(mean_meth ~ condition + sex + age, pf)) %>% broom::tidy()
# mysaveandstore(fn = sprintf("ldna/results/%s/plots/rte/pca/mean_meth_point_withage_ordered.pdf", params$mod_code), 5, 4, sf = stats)



## rDNA
cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
cytobands <- cytobandsdf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
centromere <- cytobandsdf %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()



censat <- import.bed("resources/genomes/hs1/annotations/censat.bed")
censatdf <- censat %>%
    as.data.frame() %>%
    tibble() %>%
    mutate(types = gsub("_.*", "", name))

censat <- GRanges(censatdf)

censatdf %$% types %>% table()
censatdf %>% filter(types == "rDNA")

censatdf %>% filter(seqnames == "chr4")

grs_cent <- merge_with_grs(grs, reduce(centromere))
cent <- reduce(centromere)
grs_cent %>% GRanges()
grs_cent %>%
    GRanges() %>%
    subsetByOverlaps(cent)

# region <- GRanges(seqnames = "chr13", ranges = IRanges(start = 6051123, end = 6170503))

chrofinterest <- "chr3"
roi <- cent[seqnames(cent) == chrofinterest]
rmann %>%
    filter(rte_subfamily != "Other") %>%
    filter(rte_length_req == "FL") %>%
    GRanges() %>%
    subsetByOverlaps(roi) %>%
    as.data.frame() %>%
    tibble()
highlightstemp <- censat %>%
    subsetByOverlaps(roi) %>%
    as.data.frame()
p <- plot_region(roi, grs_cent, highlights = highlightstemp, group_var = "sample", smoothing = 500, type = "line") + ggtitle(chrofinterest)
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight.pdf", params$mod_code, chrofinterest), 9, 4)
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight_tall.pdf", params$mod_code, chrofinterest), 9, 8)

chrofinterest <- "chr4"
roi <- cent[seqnames(cent) == chrofinterest]
rmann %>%
    filter(rte_subfamily != "Other") %>%
    filter(rte_length_req == "FL") %>%
    GRanges() %>%
    subsetByOverlaps(roi) %>%
    as.data.frame() %>%
    tibble()
highlightstemp <- censat %>%
    subsetByOverlaps(roi) %>%
    as.data.frame()
p <- plot_region(roi, grs_cent, highlights = highlightstemp, group_var = "sample", smoothing = 500, type = "line") + ggtitle(chrofinterest)
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight11.pdf", params$mod_code, chrofinterest), 9, 4)
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_by_sample_whighlight_tall.pdf", params$mod_code, chrofinterest), 9, 8)

centplots <- list()
for (chr in cent %$% seqnames %>% unique()) {
    p <- plot_region(cent[seqnames(cent) == chr], grs_cent, group_var = "sample", smoothing = 500, type = "line") +
        scale_samples
    centplots[[chr]] <- p
}
# library(patchwork)
# patch <- wrap_plots(centplots, ncol = 4, guides = "collect")
# mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_together.pdf", params$mod_code, chr), pl= patch, 12, 12)

p <- grs_cent %>%
    group_by(!!sym("sample"), seqnames) %>%
    mutate(plot_val = rollmean(pctM, 500, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(plot_val)) %>%
    ungroup() %>%
    ggplot() +
    geom_line(aes(x = start, y = plot_val, color = !!sym("sample"))) +
    facet_wrap(~seqnames, ncol = 4, scales = "free_x") +
    ylim(c(0, 100)) +
    scale_samples +
    mtclosed +
    theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
    )

mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_together.pdf", params$mod_code, chr), 12, 12)


grs_censat <- merge_with_grs(grs, censat)



p <- grs_censat %>%
    group_by(sample, condition, types) %>%
    summarise(pctM = mean(pctM)) %>%
    # mutate(family = fct_reorder(family, pctM, mean)) %>%
    ggplot(aes(x = types, y = pctM)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    scale_conditions +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/censat_bar.pdf", params$mod_code), 8, 4)




# next repeats in the centromere
# r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
# r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
# rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
flRTEpromoter <- read_delim(sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)



cent_l1hs <- rmann %>%
    filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "HERVK-int")) %>%
    filter(rte_length_req == "FL") %>%
    GRanges() %>%
    subsetByOverlaps(cent)

# telomeredf <- read_delim(confALL$aref$ref_telomere, col_names = FALSE, delim = "\t")
# telomere <- telomeredf %>%
#     dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
#     GRanges()

# telo_l1hs <- rmann %>%
#     filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "HERVK-int")) %>%
#     filter(rte_length_req == "FL") %>%
#     GRanges() %>%
#     subsetByOverlaps(telomere)

for (gene_id in mcols(cent_l1hs)$gene_id) {
    p <- plot_region(cent_l1hs[mcols(cent_l1hs)$gene_id == gene_id], grs_cent, group_var = "sample", smoothing = 5, type = "point") +
        scale_samples
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s.pdf", params$mod_code, gene_id), 5, 4)
    p <- plot_region(cent_l1hs[mcols(cent_l1hs)$gene_id == gene_id], grs_cent, group_var = "condition", smoothing = 5, type = "point") +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/cent_%s_condition.pdf", params$mod_code, gene_id), 5, 4)
}






rmann %>%
    filter(repeat_superfamily == "SAT") %>%
    filter(length > 10000) %>%
    arrange(-length)

sat <- rmann %>%
    filter(repeat_superfamily == "SAT") %>%
    GRanges()
satgrs <- merge_with_grs(grs, sat)

satgrs %>% pw()

satdf <- satgrs %>%
    group_by(family, sample, condition) %>%
    summarise(pctM = mean(pctM))
satdfwgeneid <- satgrs %>%
    group_by(gene_id, family, sample, condition, length) %>%
    summarise(pctM = mean(pctM))

p <- satdf %>% ggviolin(x = "family", y = "pctM", fill = "condition", add = c("dotplot")) +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/satellites1.pdf", params$mod_code, gene_id), 10, 4)

p <- satdf %>% ggbarplot(x = "family", y = "pctM", fill = "condition", add = c("dotplot"), position = position_dodge()) +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/satellites_bar1.pdf", params$mod_code, gene_id), 8, 4)


p <- satdf %>%
    group_by(family) %>%
    mutate(family = fct_reorder(family, pctM, mean)) %>%
    ggplot(aes(x = family, y = pctM)) +
    stat_summary(aes(group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
    geom_point(aes(group = condition), position = position_dodge(width = 0.9)) +
    scale_conditions +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("ldna/results/%s/plots/hardtomap/satellites_bar1111.pdf", params$mod_code, gene_id), 8, 4)

satgene <- satdfwgeneid %>%
    group_by(condition, gene_id, length) %>%
    summarise(pctM = mean(pctM)) %>%
    pivot_wider(names_from = condition, values_from = pctM) %>%
    mutate(dif = SEN - PRO, pct_change = 100 * (dif / PRO))


satgene %>% arrange(dif)




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

library(ggnewscale)
