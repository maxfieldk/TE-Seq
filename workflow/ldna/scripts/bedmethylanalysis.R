source("workflow/scripts/defaults.R")
module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

library(rtracklayer)
library(Biostrings)
library(cowplot)
# library(zoo)
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

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

{
genome_lengths <- fasta.seqlengths(conf$reference)
chromosomesAll <- names(genome_lengths)
nonrefchromosomes <- grep("nonref", chromosomesAll, value = TRUE)
refchromosomes <- grep("^chr", chromosomesAll, value = TRUE)
autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE)
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

BMAtables <- list()

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", samples, samples),
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv",
            read_mods = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

dmlspath <- inputs$dmls
dmrspath <- inputs$dmrs

cdrpath <- conf$cdr
HORpath <- conf$HOR
cenSatpath <- conf$cenSat

ccrespath <- conf$ccres
refseqgenespath <- conf$refseq



ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

# RUN IF RESUMING
if (interactive()) {
    conditions <- conf$levels
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

    grsdf <- read_delim("Rintermediates/grsdf.tsv", col_names = TRUE)
    grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
    grs <- GRanges(grsdf)
    cpg_islands <- rtracklayer::import(conf$cpg_islands)
    cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
    cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
    cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)
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

    dmrs <- read_delim(dmrspath, delim = "\t", col_names = TRUE)
    dmls <- read_delim(dmlspath, delim = "\t", col_names = TRUE)
    dmls <- dmls %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
    dmrs <- dmrs %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
    dmrs <- dmrs %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
    dmrs$direction <- factor(dmrs$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))

    dmls <- dmls %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
    dmls$direction <- factor(dmls$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))


    dmrsgr <- GRanges(dmrs)
    dmlsgr <- GRanges(
        seqnames = dmls$chr,
        ranges = IRanges(start = dmls$pos, end = dmls$pos),
        mu_c2 = dmls$mu_c2,
        mu_c1 = dmls$mu_c1,
        diff_c2_minus_c1 = dmls$diff_c2_minus_c1,
        diff_c2_minus_c1.se = dmls$diff_c2_minus_c1.se,
        stat = dmls$stat,
        phi_c2 = dmls$phi_c2,
        phi_c1 = dmls$phi_c1,
        pval = dmls$pval,
        fdr = dmls$fdr,
        postprob.overThreshold = dmls$postprob.overThreshold,
        direction = dmls$direction
    )
}

##########################
if (!interactive()) {
    # PREP DATA FOR ANALYSIS
    conditions <- conf$levels
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

    sample_grs <- list()
    for (sample_name in samples[c(1,6)]) {
        df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethlpaths, value = TRUE), col_names = FALSE)
        df_m <- df %>% filter(X4 == "m")
        df_h <- df %>% filter(X4 == "h")
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
    # filter out low coverage and ensure that all samples have the same cpgs
    grs <- grs[grs$cov > MINIMUMCOVERAGE]
    grsdf <- tibble(as.data.frame(grs))
    grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
    seqnames <- grsdf$seqnames
    start <- grsdf$start
    end <- grsdf$end
    pos <- paste0(seqnames, "_", start, "_", end)
    grsdf$pos <- pos

    grsdfuntidy <- grsdf %>%
        filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
        pivot_wider(id_cols = c("pos", "seqnames"), names_from = "sample", values_from = "pctM", names_prefix = "pctM") %>%
        drop_na()

    grsinboth <- grsdfuntidy %>% pull(pos)
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

    dir.create("Rintermediates", showWarnings = FALSE)
    write_delim(grsdf, "Rintermediates/grsdf.tsv", col_names = TRUE)

    # SETTING UP SOME SUBSETS FOR EXPLORATION
    set.seed(75)
    grsdfs <- grsdf %>%
        group_by(sample, seqnames, islandStatus) %>%
        slice_sample(n = 1000)
    grss <- GRanges(grsdfs)
}

dmrs <- read_delim(dmrspath, delim = "\t", col_names = TRUE)
dmls <- read_delim(dmlspath, delim = "\t", col_names = TRUE)
dmls <- dmls %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
dmrs <- dmrs %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
dmrs <- dmrs %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
dmrs$direction <- factor(dmrs$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))

dmls <- dmls %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
dmls$direction <- factor(dmls$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))


dmrsgr <- GRanges(dmrs)
dmlsgr <- GRanges(
    seqnames = dmls$chr,
    ranges = IRanges(start = dmls$pos, end = dmls$pos),
    mu_c2 = dmls$mu_c2,
    mu_c1 = dmls$mu_c1,
    diff_c2_minus_c1 = dmls$diff_c2_minus_c1,
    diff_c2_minus_c1.se = dmls$diff_c2_minus_c1.se,
    stat = dmls$stat,
    phi_c2 = dmls$phi_c2,
    phi_c1 = dmls$phi_c1,
    pval = dmls$pval,
    fdr = dmls$fdr,
    postprob.overThreshold = dmls$postprob.overThreshold,
    direction = dmls$direction
)

write_delim(dmrs %>% dplyr::select(chr, start, end), "results/tables/dmrs.bed", delim = "\t", col_names = FALSE)
write_delim(dmrs %>% filter(direction == "SEN Hypo") %>% dplyr::select(chr, start, end), "results/tables/dmrs_hypo.bed", delim = "\t", col_names = FALSE)
write_delim(dmrs %>% filter(direction == "SEN Hyper") %>% dplyr::select(chr, start, end), "results/tables/dmrs_hyper.bed", delim = "\t", col_names = FALSE)

#########################################################

############
# GLOBAL
dir.create("results/plots/genomewide", showWarnings = FALSE)
a <- grsdf %>%
    group_by(sample) %>%
    summarize(mean = mean(pctM), sd = sd(pctM), n = n())
t.test(pctM ~ condition, data = grsdf, var.equal = TRUE)
t.test(pctM ~ condition, data = grsdf %>% filter(seqnames %in% chromosomesNoX), var.equal = TRUE)

p <- grsdf %>% ggplot() +
    geom_boxplot(aes(x = islandStatus, y = pctM, fill = condition)) +
    mtopen + scale_conditions
mysave(fn = "results/plots/genomewide/cpgislandstatusbox.png", w = 4, h = 4, res = 300, pl = p)
plots[["cpgislandstatusbox"]] <- p

p <- grsdf %>%
    group_by(islandStatus, condition) %>%
    summarize(pctM = mean(pctM)) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = pctM, fill = condition), position = "dodge", color = "black") +
    mtopen + scale_conditions +
    anchorbar
mysave(fn = "results/plots/genomewide/cpgislandstatusbar.png", w = 4, h = 4, res = 300, pl = p)
plots[["cpgislandstatusbar"]] <- p

##################################### DMR analysis
# how many dmrs
p <- dmrs %>%
    ggplot() +
    geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
    labs(x = "") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ggtitle("DMR Counts") +
    mtopen + scale_contrasts
mysave(fn = "results/plots/genomewide/dmr_number.png", 4, 4)
plots[["dmr_number"]] <- p

# what is their average length
p <- ggplot(data = dmrs) +
    geom_histogram(aes(length), fill = "#00A087FF") +
    ggtitle("DMR Lengths") +
    labs(x = "length (bp)") +
    xlim(0, 3000) +
    mtopen + scale_samples
mysave(fn = "results/plots/genomewide/dmr_length.png", w = 4, h = 4)
plots[["dmr_length"]] <- p

p <- ggplot(data = dmrs) +
    geom_histogram(aes(length, fill = direction), alpha = 0.7) +
    ggtitle("DMR Lengths") +
    labs(x = "length (bp)") +
    xlim(0, 3000) +
    mtopen + scale_samples
mysave(fn = "results/plots/genomewide/dmr_length_stratified.png", w = 4, h = 4)
plots[["dmr_length_stratified"]] <- p

# a positive difference means that sample 1 has more methylation than sample 2
meandiff <- mean(dmrs$diff_c2_minus_c1)

p <- ggplot() +
    geom_density(
        data = dmrs,
        aes(x = diff_c2_minus_c1), fill = mycolor
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = sprintf("DMR Methylation (%s-%s)", condition1, condition2), y = "Density") +
    ggtitle("DMR Methylation Density") +
    annotate("label", x = -Inf, y = Inf, label = "SEN Hyper", hjust = 0, vjust = 1) +
    annotate("label", x = Inf, y = Inf, label = "SEN Hypo", hjust = 1, vjust = 1) +
    mtopen + scale_samples
mysave(fn = "results/plots/genomewide/dmr_delta.png", w = 4, h = 4)
plots[["dmr_delta"]] <- p

## where are they?

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

p <- dmrsgrislandStatusdf %>%
    group_by(islandStatus, direction) %>%
    summarize(mean_diff = mean(diff_c2_minus_c1), n = n()) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = n, fill = direction), position = "dodge", color = "black") +
    mtopen + scale_contrasts +
    anchorbar +
    labs(x = "", y = "Count") +
    ggtitle("DMR")
mysave(fn = "results/plots/genomewide/dmr_count_islandstatus.png", 4, 4)
plots[["dmr_count_islandstatus"]] <- p


dmrlocdf <- dmrs %>%
    group_by(chr, direction) %>%
    summarize(n = n())
dmrlocdf$chr <- factor(dmrlocdf$chr, levels = chromosomes)
p <- ggplot(data = dmrlocdf) +
    geom_col(aes(y = chr, x = n, fill = direction), position = "dodge", color = "black") +
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
    labs(x = "count", y = "") +
    ggtitle("DMR Location") +
    mtopen + scale_contrasts
mysave(fn = "results/plots/genomewide/dmr_count.png", 5, 5)
plots[["dmr_count"]] <- p

p <- dmrs %>%
    ggplot(aes(x = meanMethy_c1, y = meanMethy_c2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(palette = "Spectral", direction = 1) +
    xlab(sprintf("CpG Methylation %s", condition1)) +
    ylab(sprintf("CpG Methylation %s", condition2)) +
    ggtitle("DMR Density") +
    mtopen
mysave(fn = "results/plots/genomewide/dmrdensity.png", 5, 5)
plots[["dmrdensity"]] <- p

p <- dmls %>%
    ggplot() +
    geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
    ggtitle("DML Counts") +
    labs(x = "", y = "Count") +
    anchorbar +
    mtopen + scale_contrasts
mysave(fn = "results/plots/genomewide/dml_count.png", 4, 4)
plots[["dml_count"]] <- p
# dmls
p <- ggplot() +
    geom_density(data = dmls, aes(x = diff_c2_minus_c1), fill = mycolor) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "DML Methylation (PRO-SEN)", y = "Density") +
    ggtitle("DML Methylation Density") +
    annotate("label", x = -Inf, y = Inf, label = paste0(condition2, " Hypo"), hjust = 0, vjust = 1) +
    annotate("label", x = Inf, y = Inf, label = paste0(condition2, " Hyper"), hjust = 1, vjust = 1) +
    mtopen + scale_contrasts

mysave(fn = "results/plots/genomewide/dml_delta.png", 4, 4)
plots[["dml_delta"]] <- p

dmllocdf <- dmls %>%
    group_by(chr, direction) %>%
    summarize(n = n())
dmllocdf$chr <- factor(dmllocdf$chr, levels = chromosomes)

p <- dmls %>%
    filter(chr %in% chromosomes) %>%
    ggplot(aes(x = mu_c1, y = mu_c2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_fill_distiller(palette = "Spectral", direction = 1) +
    xlab(sprintf("CpG Methylation %s", condition1)) +
    ylab(sprintf("CpG Methylation %s", condition2)) +
    ggtitle("DML Density") +
    mtopen
mysave(fn = "results/plots/genomewide/dmldensity.png", w = 5, h = 5)
plots[["dmldensity"]] <- p



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
    summarize(mean_diff = mean(diff_c2_minus_c1), n = n()) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = n, fill = direction), position = "dodge", color = "black") +
    mtopen + scale_contrasts +
    anchorbar +
    labs(x = "", y = "Count") +
    ggtitle("DML Island Status")
mysave(fn = "results/plots/genomewide/dml_count_islandstatus.png", 4, 4)
plots[["dml_count_islandstatus"]] <- p




dir.create("results/plots/figs")

# patch1 <- (pl_meth_by_chromosome)
# patch2 <- (pl_ndml + pl_dml_raster + pl_cdml)
# patch3 <- (pl_ndmr + pl_dmr_raster + pl_cdmr)
# patch <- patch1 + patch2 + patch3 + plot_layout(nrow = 3, ggene_ides = "collect") + plot_annotation(tag_levels = "A")
# png("results/plots/figs/global.png", height = 20, width = 20, res = 300, units = "in")
# print(patch & theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 14),
#     plot.title = element_text(size = 18),
#     plot.tag = element_text(size = 32)
# ))
# dev.off()

# patch1 <- (pl_cdml_island + pl_cdmr_island)
# patch2 <- (pl_ndml + pl_dml_raster + pl_cdml)
# patch3 <- (pl_ndmr + pl_dmr_raster + pl_cdmr)
# patch <- patch1 + patch2 + patch3 + plot_layout(nrow = 3, ggene_ides = "collect") + plot_annotation(tag_levels = "A")
# p <- patch & theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 14),
#     plot.title = element_text(size = 18),
#     plot.tag = element_text(size = 32)
# )
# mysave(fn = "results/plots/figs/global_2.png", w = 20, h = 20, res = 300, pl = p)


# rm(pl_meth_by_chromosome, patch, patch1, patch2, patch3, pl_ndml, pl_ddml, pl_cdml, pl_ndmr, pl_ddmr, pl_cdmr)
##########



### gimmemotifs
gimmehyper <- read_delim("results/gimme/hyper/gimme.roc.report.txt", delim = "\t")
gimmehypo <- read_delim("results/gimme/hypo/gimme.roc.report.txt", delim = "\t")

gimmehyper$Motif <- gsub("_HUMAN.H11MO.0.*", "", gimmehyper$Motif)
pl_gimme_hyper <- gimmehyper %>%
    arrange(`P-value`) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Motif, `log10 P-value`), y = `log10 P-value`)) +
    geom_col(color = "black", fill = mycolor) +
    scale_x_discrete(labels = scales::label_wrap(60)) +
    coord_flip() +
    ggtitle("Hyper Top 10 Motifs") +
    scale_color_continuous(trans = "reverse") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "") +
    theme(
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
png("results/plots/gimme/gimmebarplothyper.png", 3, 4, units = "in", res = 300)
print(pl_gimme_hyper)
dev.off()

gimmehypo %>%
    arrange(`P-value`) %>%
    head(10)
gimmehypo$Motif <- gsub("_HUMAN.H11MO.0.*", "", gimmehypo$Motif)
pl_gimme_hypo <- gimmehypo %>%
    arrange(`P-value`) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Motif, `log10 P-value`), y = `log10 P-value`)) +
    geom_col(color = "black", fill = mycolor) +
    scale_x_discrete(labels = scales::label_wrap(60)) +
    coord_flip() +
    scale_color_continuous(trans = "reverse") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "") +
    ggtitle("Hypo Top 10 Motifs") +
    theme(
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
png("results/plots/gimme/gimmebarplothypo.png", 3, 4, units = "in", res = 300)
print(pl_gimme_hypo)
dev.off()

p <- ggplot(gimmehypo, aes(x = `Recall at 10% FDR`)) +
    geom_histogram()
png("results/plots/gimme/histogramPval.png", 4, 4, units = "in", res = 300)
print(p)
dev.off()

####################
## RTEs
dir.create("results/plots/rte", showWarnings = FALSE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
r_repeatmasker_annotation %$% ltr_viral_status_req %>% unique()
RMdf <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
RM <- GRanges(RMdf)
### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        big_ontology_groups <- c(big_ontology_groups, RMdf %>%
            pull(!!sym(ontology)) %>%
            unique())
    }
    big_ontology_groups <- big_ontology_groups %>% unique()

    modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
    region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
    element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
}

# FULL ELEMENTS
# annotate whether full length elements promoters overlap DMRs
mbo <- mergeByOverlaps(RM, dmrsgr)
mergeddf <- tibble(as.data.frame(mbo))
mm <- mergeddf %>%
    group_by(gene_id) %>%
    mutate(max_val = max(dmrsgr.diff_c2_minus_c1), min_val = min(dmrsgr.diff_c2_minus_c1))
dups <- mm$gene_id %>% duplicated()
mm <- mm[!dups, ]
mm[(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "discordant"
mm[!(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "concordant"
mm <- mm %>%
    ungroup() %>%
    dplyr::rename(seqnames = RM.seqnames, start = RM.start, end = RM.end, strand = RM.strand)
# drop columns that start with RM.
mm <- mm %>% dplyr::select((!starts_with("RM.")))
merged <- GRanges(mm)

sboinvert <- subsetByOverlaps(RM, dmrsgr, invert = TRUE)
sboinvert$max_val <- NaN
sboinvert$min_val <- NaN
sboinvert$concordance <- NA_character_

RMfinal <- c(merged, sboinvert)
RMdffinal <- tibble(as.data.frame(RMfinal))
RMdfGOOD <- RMdffinal %>%
    mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))))
RMdmrs <- RMdfGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)


grouping_var <- "rte_length_req"
classes <- RMdmrs %>%
    pull(!!sym(grouping_var)) %>%
    unique() %>%
    na.omit()
classes <- classes[!str_detect(classes, "Other")]
rtegrl <- GRangesList()
for (rte in classes) {
    print(rte)
    mbo <- mergeByOverlaps(grs, GRanges(RMdmrs %>% filter(!!sym(grouping_var) == rte)))
    methgr <- mbo$grs
    methgr$type <- rte
    methgr$rte_family <- mbo$rte_family
    methgr$rte_subfamily_limited <- mbo$rte_subfamily_limited
    methgr$rte_subfamily <- mbo$rte_subfamily
    methgr$gene_id <- mbo$gene_id
    methgr$genic_loc <- mbo$genic_loc
    methgr$concordance <- mbo$concordance
    methgr$max_val <- mbo$max_val
    methgr$min_val <- mbo$min_val
    methgr$direction <- mbo$direction
    methgr$rte_length_req <- mbo$rte_length_req
    methgr$l1_intactness_req <- mbo$l1_intactness_req
    methgr$ltr_viral_status_req <- mbo$ltr_viral_status_req
    rtegrl[[rte]] <- methgr
}

rtegr <- unlist(rtegrl)
names(rtegr) <- NULL
rtedf <- tibble(as.data.frame(rtegr))

write_delim(rtedf, "Rintermediates/rtedf.tsv", col_names = TRUE)
# rtedf <- read_delim("Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf <- rtedf %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, rte_length_req, type, l1_intactness_req, ltr_viral_status_req) %>%
    summarize(mean_meth = mean(pctM))

perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))
# mask <- perelementdf$condition == "healthy"
# newC <- ifelse(mask, "ctrl", "alz")
# perelementdf$condition <- newC

write_delim(perelementdf, "Rintermediates/perelementdf.tsv", col_names = TRUE)
# perelementdf <- read_delim("Rintermediates/perelementdf.tsv", col_names = TRUE)


# PROMOTERS
# annotate whether full length elements promoters overlap DMRs
RMdf %$% ltr_viral_status_req %>% unique()
flelement <- RMdf %>% filter(str_detect(rte_length_req, ">"))
flSINE <- flelement %>% filter(rte_superfamily == "SINE")
flLINE <- flelement %>% filter(rte_superfamily == "LINE")
flFl_Provirus_5LTR <- flelement %>%
    filter(str_detect(gene_id, "LTR")) %>%
    filter(ltr_viral_status_req == "Fl_Provirus_5LTR")

flSINEgrs <- GRanges(flSINE)
flLINE5UTRgrs <- GRanges(flLINE) %>% resize(909)
flFl_Provirus_5LTRgrs <- GRanges(flFl_Provirus_5LTR)
flRTEpromotergrs <- c(c(flSINEgrs, flLINE5UTRgrs), flFl_Provirus_5LTRgrs)

mbo <- mergeByOverlaps(flRTEpromotergrs, dmrsgr)
mergeddf <- tibble(as.data.frame(mbo))
mm <- mergeddf %>%
    group_by(gene_id) %>%
    mutate(max_val = max(dmrsgr.diff_c2_minus_c1), min_val = min(dmrsgr.diff_c2_minus_c1))
dups <- mm$gene_id %>% duplicated()
mm <- mm[!dups, ]
mm[(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "discordant"
mm[!(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "concordant"
mm <- mm %>%
    ungroup() %>%
    dplyr::rename(seqnames = flRTEpromotergrs.seqnames, start = flRTEpromotergrs.start, end = flRTEpromotergrs.end, strand = flRTEpromotergrs.strand)
# drop columns that start with RM.
mm <- mm %>% dplyr::select((!starts_with("flRTEpromotergrs.")))
merged <- GRanges(mm)

sboinvert <- subsetByOverlaps(flRTEpromotergrs, dmrsgr, invert = TRUE)
sboinvert$max_val <- NaN
sboinvert$min_val <- NaN
sboinvert$concordance <- NA_character_

flRTEpromotergrsfinal <- c(merged, sboinvert)
flRTEpromoterdfinal <- tibble(as.data.frame(flRTEpromotergrsfinal))
flRTEpromoterdfGOOD <- flRTEpromoterdfinal %>%
    mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))))
flRTEpromoter <- flRTEpromoterdfGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
flRTEpromoter %$% rte_subfamily %>% unique()


pff <- flRTEpromoter %>%
    group_by(rte_length_req, genic_loc) %>%
    mutate(group_n = n()) %>%
    group_by(rte_length_req, direction, genic_loc) %>%
    summarise(n = n(), group_n = dplyr::first(group_n)) %>%
    mutate(frac_dm = n / group_n) %>%
    filter(!is.na(rte_length_req)) %>%
    filter(!is.na(direction)) %>%
    filter(direction != "discordant") %>%
    ungroup() %>%
    complete(rte_length_req, direction, genic_loc, fill = list(n = 0, group_n = 0, frac_dm = 0))
numdf <- flRTEpromoter %>%
    group_by(rte_length_req, genic_loc) %>%
    summarise(group_n_accurate = n())
p <- pff %>%
    left_join(numdf) %>%
    mutate(ann_axis = paste0(rte_length_req, "\n", "n=", group_n_accurate)) %>%
    ggplot() +
    geom_col(aes(x = ann_axis, y = frac_dm, fill = direction), position = "dodge", color = "black") +
    facet_wrap(~genic_loc, scales = "free_x", nrow = 2) +
    labs(x = "", y = "Fraction Differentially Methylated") +
    ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
    mtopen + scale_contrasts
mysave(sprintf("results/plots/rte/dmfl%s_promoter.png", "all"), 10, 7)
plots[[sprintf("dmfl%s_promoter", "all")]] <- p

classes <- flRTEpromoter$rte_length_req %>%
    unique() %>%
    na.omit()
classes <- classes[!str_detect(classes, "Other")]
classes <- classes[str_detect(classes, ">")]
rtegrl_promoters <- GRangesList()
for (rte in classes) {
    print(rte)
    mbo <- mergeByOverlaps(grs, GRanges(flRTEpromoter %>% filter(rte_length_req == rte)))
    methgr <- mbo$grs
    methgr$type <- rte
    methgr$rte_family <- mbo$rte_family
    methgr$rte_subfamily_limited <- mbo$rte_subfamily_limited
    methgr$rte_subfamily <- mbo$rte_subfamily
    methgr$gene_id <- mbo$gene_id
    methgr$genic_loc <- mbo$genic_loc
    methgr$concordance <- mbo$concordance
    methgr$max_val <- mbo$max_val
    methgr$min_val <- mbo$min_val
    methgr$direction <- mbo$direction
    methgr$rte_length_req <- mbo$rte_length_req
    methgr$l1_intactness_req <- mbo$l1_intactness_req
    methgr$ltr_viral_status_req <- mbo$ltr_viral_status_req
    rtegrl_promoters[[rte]] <- methgr
}

rtegr_promoters <- unlist(rtegrl_promoters)
names(rtegr_promoters) <- NULL
rtedf_promoters <- tibble(as.data.frame(rtegr_promoters))

write_delim(rtedf_promoters, "Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
# rtedf <- read_delim("Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, rte_length_req, type, l1_intactness_req, ltr_viral_status_req) %>%
    summarize(mean_meth = mean(pctM))

perelementdf_promoters <- perelementdf_promoters %>% filter(!is.na(rte_length_req))
# mask <- perelementdf$condition == "healthy"
# newC <- ifelse(mask, "ctrl", "alz")
# perelementdf$condition <- newC

write_delim(perelementdf_promoters, "Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
# perelementdf <- read_delim("Rintermediates/perelementdf.tsv", col_names = TRUE)

p <- perelementdf_promoters %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_length_req, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_length_req, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen + scale_conditions
mysave(fn = "results/plots/rte/repmasker_boxplot_promoters.png", 14, 6)
plots[["repmasker_boxplot_promoters"]] <- p

#################


l1hsintactmethgr <- rtedf %>%
    filter(str_detect(l1_intactness_req, "Intact")) %>%
    left_join(r_annotation_fragmentsjoined %>% dplyr::select(gene_id, start, end, strand) %>% dplyr::rename(element_strand = strand, element_start = start, element_end = end))
l1hsintactmethgr <- l1hsintactmethgr %>%
    mutate(rel_start = start - element_start) %>%
    mutate(rel_end = end - element_start)
write_delim(l1hsintactmethgr, "Rintermediates/l1hsintactdf.tsv", col_names = TRUE)

pf_pos <- l1hsintactmethgr %>%
    filter(element_strand == "+") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 5, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()
pf_neg <- l1hsintactmethgr %>%
    filter(element_strand == "-") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 5, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()

p <- pf_pos %>% ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtopen + scale_conditions
png("results/plots/rte/l1intact_Lines_pos_strand.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

p <- pf_pos %>%
    filter(rel_start < 910) %>%
    ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    xlim(c(1, 910)) +
    ylim(c(0, 100)) +
    mtopen + scale_conditions
png("results/plots/rte/l1intact_Lines_pos_strand_promoter.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

p <- pf_neg %>% ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtopen + scale_conditions
png("results/plots/rte/l1intact_Lines_neg_strand.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

p <- pf_neg %>%
    filter(rel_start < 910) %>%
    ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    xlim(c(1, 910)) +
    ylim(c(0, 100)) +
    mtopen + scale_conditions
png("results/plots/rte/l1intact_Lines_neg_strand_promoter.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

# chr4_83273990_83280141_+
# chr6_44705376_44711532_+
# for the fig
elements_of_interest <- c("chr6_44705376_44711532_+", "chr4_83273990_83280141_+")
### nice L1 annot
y_valmin <- 0
y_valmax <- 10
modifier <- 44705376
`5UTR` <- c(1, 909) + modifier
ORF1 <- c(911, 1927) + modifier
ORF2 <- c(1991, 5818) + modifier
EN <- c(1991, 2707) + modifier
RT <- c(3482, 4309) + modifier
`3UTR` <- c(5818, 6050) + modifier
l1.2.coord <- c(1, 6050) + modifier
data <- list(
    `5UTR` = `5UTR`,
    ORF1 = ORF1,
    ORF2 = ORF2,
    `3UTR` = `3UTR`
)
color_intervals <- tibble(
    start = sapply(data, "[", 1),
    end = sapply(data, "[", 2),
    Region = factor(names(data), levels = names(data))
)
###

pl_elementofinterest1 <- pf %>%
    filter(gene_id %in% elements_of_interest[1]) %>%
    ggplot() +
    geom_line(aes(x = start, y = rM, color = condition)) +
    ylim(c(0, 100)) +
    ggtitle(paste0("L1HS chr1_113508780_113514958_+")) +
    labs(y = "Intra-Element Methylation Rolling Mean") +
    geom_rect(aes(xmin = l1.2.coord[1], ymin = y_valmin, xmax = l1.2.coord[2], ymax = y_valmax), fill = "grey") +
    geom_rect(data = color_intervals, aes(xmin = start, xmax = end, ymin = y_valmin, ymax = y_valmax, fill = Region), alpha = 1) +
    geom_text(data = color_intervals, aes(x = (start + end) / 2, y = y_valmax, label = Region), vjust = -0.5) +
    mtopen + scale_conditions +
    scale_fill_manual(values = c(my_palette[5:4], paletteer::paletteer_d("dutchmasters::pearl_earring", 4)))
png("results/plots/rte/element_of_interest1.png", 5, 5, units = "in", res = 300)
print(pl_elementofinterest1)
dev.off()

y_valmin <- 0
y_valmax <- 10
modifier <- 83273990
`5UTR` <- c(1, 909) + modifier
ORF1 <- c(911, 1927) + modifier
ORF2 <- c(1991, 5818) + modifier
EN <- c(1991, 2707) + modifier
RT <- c(3482, 4309) + modifier
`3UTR` <- c(5818, 6050) + modifier
l1.2.coord2 <- c(1, 6050) + modifier
data <- list(
    `5UTR` = `5UTR`,
    ORF1 = ORF1,
    ORF2 = ORF2,
    `3UTR` = `3UTR`
)
color_intervals2 <- tibble(
    start = sapply(data, "[", 1),
    end = sapply(data, "[", 2),
    Region = factor(names(data), levels = names(data))
)
###
pl_elementofinterest2 <- pf %>%
    filter(gene_id %in% elements_of_interest[2]) %>%
    ggplot() +
    geom_line(aes(x = start, y = rM, color = condition)) +
    ylim(c(0, 100)) +
    ggtitle(paste0("L1HS chr4_83273990_83280141_+")) +
    labs(y = "Intra-Element Methylation Rolling Mean") +
    geom_rect(aes(xmin = l1.2.coord2[1], ymin = y_valmin, xmax = l1.2.coord2[2], ymax = y_valmax), fill = "grey") +
    geom_rect(data = color_intervals2, aes(xmin = start, xmax = end, ymin = y_valmin, ymax = y_valmax, fill = Region), alpha = 1) +
    geom_text(data = color_intervals2, aes(x = (start + end) / 2, y = y_valmax, label = Region), vjust = -0.5) +
    mtopen + scale_conditions +
    scale_fill_manual(values = c(my_palette[5:4], paletteer::paletteer_d("dutchmasters::pearl_earring", 4)))
png("results/plots/rte/element_of_interest2.png", 5, 5, units = "in", res = 300)
print(pl_elementofinterest2)
dev.off()


{ # now LTR5s
    LTR5Adf <- rtedf %>%
        filter(type == "LTR5A") %>%
        separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
        filter(element_stop - element_start > 600) %>%
        filter(cov > MINIMUMCOVERAGE)
    LTR5Bdf <- rtedf %>%
        filter(type == "LTR5B") %>%
        separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
        filter(element_stop - element_start > 600) %>%
        filter(cov > MINIMUMCOVERAGE)
    LTR5_Hsdf <- rtedf %>%
        filter(type == "LTR5_Hs") %>%
        separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
        filter(element_stop - element_start > 600) %>%
        filter(cov > MINIMUMCOVERAGE)
    LTR5df <- rtedf %>%
        filter(type == "LTR5") %>%
        separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
        filter(element_stop - element_start > 600) %>%
        filter(cov > MINIMUMCOVERAGE)

    pf <- LTR5_Hsdf %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        group_by(gene_id, condition) %>%
        mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
        filter(!is.na(rM)) %>%
        ungroup()
    p <- pf %>% ggplot() +
        geom_line(aes(x = start, y = rM, color = condition)) +
        scale_x_continuous(breaks = scales::breaks_pretty(3)) +
        facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
        ylim(c(0, 100)) +
        mtopen + scale_conditions
    png("results/plots/rte/LTR5_Hs_Lines.png", 12, 120, units = "in", res = 200)
    print(p)
    dev.off()

    pf <- LTR5Bdf %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        group_by(gene_id, condition) %>%
        mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
        filter(!is.na(rM)) %>%
        ungroup()
    p <- pf %>% ggplot() +
        geom_line(aes(x = start, y = rM, color = condition)) +
        scale_x_continuous(breaks = scales::breaks_pretty(3)) +
        facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
        ylim(c(60, 100)) +
        mtopen + scale_conditions
    png("results/plots/rte/LTR5B_Lines.png", 12, 120, units = "in", res = 200)
    print(p)
    dev.off()

    pf <- LTR5df %>%
        filter(cov > MINIMUMCOVERAGE) %>%
        group_by(gene_id, condition) %>%
        mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
        filter(!is.na(rM)) %>%
        ungroup()
    p <- pf %>% ggplot() +
        geom_line(aes(x = start, y = rM, color = condition)) +
        scale_x_continuous(breaks = scales::breaks_pretty(3)) +
        facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
        ylim(c(60, 100)) +
        mtopen + scale_conditions
    png("results/plots/rte/LTR5_Lines.png", 12, 120, units = "in", res = 200)
    print(p)
    dev.off()
}


# ####



l1hsintactmethdf <- l1hsintactmethgr %>%
    as.data.frame() %>%
    tibble()
for (gene_id in l1hsintactmethdf %$% gene_id %>% unique()) {
    pf <- l1hsintactmethdf %>%
        filter(gene_id == gene_id) %>%
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
        mtopen + scale_conditions
    dir.create("results/plots/rte/l1hsintact")
    png(paste0("results/plots/rte/l1hsintact/", gene_id, ".png"), 8, 3, units = "in", res = 300)
    print(p)
    dev.off()
}



# heatmap 5UTR
heatmapprep <- l1hsintactmethdf %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
        element_strand == "-" ~ (start > element_end - 909) & (start < element_end)
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
l1hsintactdf <- l1hsintactmethdf %>% group_by(gene_id) %>% summarise(concordance = dplyr::first(concordance), genic_loc = dplyr::first(genic_loc))
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
row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), genic_loc = genic_locs, col = list(genic_loc = c("Genic" = "brown", "NonGenic" = "tan")))
conditions <- c(sample_table %>% filter(condition == condition1) %$% condition, sample_table %>% filter(condition == condition2) %$% condition)
topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))
library(circlize)
col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
col_fun(seq(0, 50, by = 25))
heatmapL1UTR <- m %>%
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

png("results/plots/l1intactheatmap_5utr.png", 9, 14, units = "in", res = 300)
ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right")
dev.off()

plgrob <- grid.grabExpr(ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right"))
plots[["l1intactheatmap_5utr"]] <- plgrob

############## Read density
rm(grs)
mem_used()

dir.create("results/plots/reads")
library(Biostrings)

readslist <- list()
for (region in conf$rte_subfamily_read_level_analysis) {
    for (sample_name in samples) {
        df <- read_delim(
            grep(region,
                grep(sprintf("/%s/", sample_name),
                    inputs$read_mods,
                    value = TRUE
                ),
                value = TRUE
            )
        )
        df$region <- region
        df$sample <- sample_name
        df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grs <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0(conf$reference_annotation_dir, "/annotations/rte_beds/", region, ".bed"))
        mbo <- mergeByOverlaps(grs, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grs"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grs.", "", colnames(df1))
        readslist <- c(readslist, list(df1))
    }
}
reads <- Reduce(rbind, readslist)

readslistcg <- list()
for (region in conf$rte_subfamily_read_level_analysis) {
    for (sample_name in samples) {
        df <- read_delim(
            grep(region,
                grep(sprintf("/%s/", sample_name),
                    inputs$read_mods_cg,
                    value = TRUE
                ),
                value = TRUE
            )
        )
        df$region <- region
        df$sample <- sample_name
        df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grs <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0(conf$reference_annotation_dir, "/annotations/rte_beds/", region, ".bed"))
        mbo <- mergeByOverlaps(grs, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grs"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grs.", "", colnames(df1))
        readslistcg <- c(readslistcg, list(df1))
    }
}

readscg <- Reduce(rbind, readslistcg)
utr1 %>% filter(mod_code == "m") %>%
        mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0)) %$% mod_indicator %>% table()


read_analysis <- function(readsdf,
    region = "L1HS_l1_intactness_req_ALL",
    mod_code_var = "m",
    context = "CG") {

    readsdf1 <- readsdf %>% left_join(r_annotation_fragmentsjoined %>% dplyr::select(gene_id, start, end, strand) %>% dplyr::rename(element_strand = strand, element_start = start, element_end = end))
    utr1 <- readsdf1 %>%
        filter(mod_code == mod_code_var) %>%
        filter(case_when(
            element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
            element_strand == "-" ~ (start > element_end - 909) & (start < element_end)
        )) %>% dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))
    
    utr <- utr1 %>%
        group_by(gene_id, read_id, condition) %>%
        mutate(read_span = max(start) - min(start)) %>%
        mutate(num_cpgs_in_read = n()) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        relocate(gene_id) %>%
        ungroup()

    utr2 <- utr1 %>%
        group_by(gene_id, read_id, condition) %>%
        summarise(read_span = max(start) - min(start), num_cpgs_in_read = n(), fraction_meth = mean(mod_indicator)) 

    write_delim(utr, sprintf("Rintermediates/%s_%s_%s_reads.tsv", region, mod_code_var, context), delim = "\t")

    p <- utr %>% ggplot() +
        geom_density(aes(x = mod_qual, fill = condition), alpha = 0.3) +
        facet_wrap(vars(region))
    mysave(sprintf("results/plots/reads/modbase_score_dist_%s_%s_%s.png", region, mod_code_var, context), 12, 4, pl = p)


    p <- utr %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region))
    mysave(sprintf("results/plots/reads/read_span_distribution_%s_%s_%s.png", region, mod_code_var, context), 5, 5, pl = p)


    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(num_cpgs_in_read), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region))
    mysave(sprintf("results/plots/reads/read_num_cpg_distribution_%s_%s_%s.png", region, mod_code_var, context), 5, 5, pl = p)



    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(read_id, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_point(aes(x = read_span, y = fraction_meth, color = condition)) +
        facet_wrap(vars(region))
    mysave(sprintf("results/plots/reads/fraction_meth_distribution_%s_%s_%s.png", region, mod_code_var, context), 5, 5, pl = p)


    aa <- utr %>%
        filter(read_span > 600) %>%
        group_by(read_id, sample, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand))

    ab <- aa %>%
        mutate(unmeth = ifelse(fraction_meth > 0.5, 0, 1)) %>%
        group_by(sample, region, condition) %>%
        summarise(propUnmeth = mean(unmeth)) %>%
        group_by(condition, region) %>%
        mutate(meanProp = mean(propUnmeth))

    p <- ab %>%
        ggplot(aes(x = condition)) +
        stat_summary(aes(y = propUnmeth, fill = condition), color = "black", fun = "mean", geom = "bar") +
        geom_point(aes(y = propUnmeth)) +
        facet_wrap(vars(region)) +
        labs(x = "", y = "Proportion of Reads less than 50% methylated") +
        ggtitle("L1 5'UTR Methylation") +
        mtopen + scale_conditions +
        anchorbar
    mysave(sprintf("results/plots/reads/barplot_50pct_%s_%s_%s.png", region, mod_code_var, context), 4, 4, pl = p)
    plots[["barplot_50pct"]][[region]][[mod_code_var]][[context]] <- p

    # note that this isn't the right test for this

    # result <- ab %>%
    #     filter(region == "l1hsintact") %>%
    #     t.test(propUnmeth ~ condition, data = .)
    # result$p.value
}

read_analysis(reads, "L1HS_l1_intactness_req_ALL", "m", "NoContext")
read_analysis(reads, "L1HS_l1_intactness_req_ALL", "h", "NoContext")
read_analysis(reads, "L1HS_l1_intactness_req_ALL", "a", "NoContext")
read_analysis(readscg, "L1HS_l1_intactness_req_ALL", "m", "CpG")
read_analysis(readscg, "L1HS_l1_intactness_req_ALL", "h", "CpG")




utr %>%
    filter(read_span > 600) %>%
    group_by(read_id, sample, condition, region) %>%
    summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand))
p <- aa %>%
    ggplot() +
    geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.3) +
    labs(x = "Pct CpG Methylation per Read") +
    ggtitle("L1 5'UTR Methylation") +
    facet_wrap(vars(region)) +
    mtopen + scale_conditions
mysave("results/plots/reads/fraction_meth_density_distribution.png", 5, 7.5)
plots[["fraction_meth_density_distribution"]] <- p

p <- utr %>%
    filter(region == "L1HS_l1_intactness_req_ALL") %>%
    filter(read_span > 600) %>%
    group_by(read_id, condition, region) %>%
    summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
    ggplot() +
    geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.3) +
    labs(x = "Pct CpG Methylation per Read") +
    ggtitle("L1HS Intact 5'UTR Methylation") +
    mtopen + scale_conditions +
    anchorbar

mysave("results/plots/reads/fraction_meth_density_distribution_l1hsintact.png", 5, 7.5)
plots[["fraction_meth_density_distribution_l1hsintact"]] <- pl_l1hsintact_density


######### GENES


{
    directions <- c("Hypo", "Hyper", "Dif")
    mydir <- "results/plots/great"
    mydirtables <- "results/tables/great"
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

    et <- extendTSS(genes_gr, genome_lengths, gene_id_type = "SYMBOL")
}

## reactome
# {
#     dir.create(paste(mydir, "reactome", sep = "/"))
#     dir.create(paste(mydirtables, "reactome", sep = "/"))
#     gs <- as.list(reactomePATHID2EXTID)
#     tablesReactome <- list()
#     resultsReactome <- list()
#     for (direction in directions) {
#         if (direction == "Hypo") {
#             regions <- dmrsgr[grepl("Hypo", dmrsgr$direction)]
#         } else if (direction == "Hyper") {
#             regions <- dmrsgr[grepl("Hyper", dmrsgr$direction)]
#         } else {
#             regions <- dmrsgr
#         }

#         res <- great(regions, gene_sets = gs, extended_tss = et, background = chromosomes)
#         resultsReactome[[direction]] <- res

#         tb <- getEnrichmentTable(res)
#         anno.result <- mapIds(reactome.db,
#             keys = tb$id,
#             column = "PATHNAME", keytype = "PATHID", multiVals = "first"
#         )
#         tb <- cbind(tb, anno.result) %>%
#             dplyr::arrange(p_adjust) %>%
#             dplyr::relocate(anno.result)
#         tablesReactome[[direction]] <- tb
#         write_delim(tb, paste(mydirtables, "reactome", paste0(direction, "great_enrichment.tsv"), sep = "/"), delim = "\t")

#         png(paste(mydir, "reactome", paste0(direction, "volcano.png"), sep = "/"), height = 5, width = 5, res = 300, units = "in")
#         plotVolcano(res)
#         dev.off()

#         png(paste(mydir, "reactome", paste0(direction, "associations.png"), sep = "/"), height = 5, width = 10, res = 300, units = "in")
#         plotRegionGeneAssociations(res)
#         dev.off()
#     }

#     reactomePlots <- list()
#     for (direction in directions) {
#         p <- tablesReactome[[direction]] %>%
#             head(n = 20) %>%
#             mutate(anno.result = gsub("Homo sapiens: ", "", anno.result)) %>%
#             mutate(anno.result = fct_reorder(anno.result, fold_enrichment)) %>%
#             ggplot(aes(x = anno.result, y = fold_enrichment)) +
#             geom_segment(aes(x = anno.result, xend = anno.result, y = 0, yend = fold_enrichment, color = p_adjust)) +
#             geom_point(aes(color = p_adjust), size = 4, alpha = 0.6) +
#             scale_x_discrete(labels = scales::label_wrap(60)) +
#             coord_flip() +
#             scale_color_continuous(trans = "reverse") +
#             labs(x = "") +
#             theme(
#                 axis.text.y = element_text(color = "black"),
#                 panel.grid.major.y = element_blank(),
#                 panel.border = element_rect(color = "black", fill = NA, size = 1),
#                 axis.ticks.y = element_blank()
#             )
#         reactomePlots[[direction]] <- p
#         png(paste(mydir, "reactome", paste0(direction, "lollipop.png"), sep = "/"), height = 6, width = 6, res = 300, units = "in")
#         print(p)
#         dev.off()
#     }
# }
library(clusterProfiler)

{
    greatplots <- list()
    genecollections <- names(conf[["genesets_for_great"]])
    for (collection in genecollections) {
        tryCatch({
            dir.create(paste(mydir, collection, sep = "/"))
            dir.create(paste(mydirtables, collection, sep = "/"))
            genesets <- read.gmt(conf[["genesets_for_great"]][[collection]])
            tablesMsigdb <- list()
            results <- list()
            for (direction in directions) {
                if (direction == "Hypo") {
                    regions <- dmrsgr[grepl("Hypo", dmrsgr$direction)]
                } else if (direction == "Hyper") {
                    regions <- dmrsgr[grepl("Hyper", dmrsgr$direction)]
                } else {
                    regions <- dmrsgr
                }


                res <- great(regions, gene_sets = genesets, extended_tss = et, background = CHROMOSOMESINCLUDEDINANALYSIS_REF)
                results[[collection]][[direction]] <- res

                tb <- getEnrichmentTable(res)

                tb <- tb %>% dplyr::arrange(p_adjust)
                tablesMsigdb[[collection]][[direction]] <- tb
                write_delim(tb, paste(mydirtables, collection, paste0(direction, "great_enrichment.tsv"), sep = "/"))
                
                png(paste(mydir, collection, paste0(direction, "volcano.png"), sep = "/"), height = 5, width = 5, res = 300, units = "in")
                plotVolcano(res)
                dev.off()

                png(paste(mydir, collection, paste0(direction, "associations.png"), sep = "/"), height = 5, width = 10, res = 300, units = "in")
                plotRegionGeneAssociations(res)
                dev.off()

                p <- tablesMsigdb[[collection]][[direction]] %>%
                    head(n = 20) %>%
                    mutate(id = fct_reorder(id, fold_enrichment)) %>%
                    ggplot(aes(x = id, y = fold_enrichment)) +
                    geom_segment(aes(x = id, xend = id, y = 0, yend = fold_enrichment, color = p_adjust)) +
                    geom_point(aes(color = p_adjust), size = 4, alpha = 0.6) +
                    coord_flip() +
                    scale_color_continuous(trans = "reverse") +
                    labs(x = "") +
                    theme(
                        panel.grid.major.y = element_blank(),
                        panel.border = element_rect(color = "black", fill = NA, size = 1),
                        axis.ticks.y = element_blank()
                    )
                greatplots[[collection]][[direction]] <- p
                mysave(paste(mydir, collection, paste0(direction, "lollipop.png"), sep = "/"), 5, 10)
            }
        })
    }
    save(greatplots, file = "results/plots/objects/greatplots.rda")
}


{
    promotersplots <- list()
    promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
    write_delim(tibble(as.data.frame(promoters)) %>% mutate(score = 1000) %>% dplyr::select(seqnames, start, end, gene_id, score, strand), "Rintermediates/promoters.bed", col_names = FALSE, delim = "\t")

    hyporegions <- dmrsgr[grepl("Hypo", dmrsgr$direction)]
    hyperregions <- dmrsgr[grepl("Hyper", dmrsgr$direction)]

    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, dmrsgr, ignore.strand = TRUE))), "Rintermediates/promoters_dmregions.bed", col_names = FALSE, delim = "\t")
    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyporegions, ignore.strand = TRUE))), "Rintermediates/promoters_dmhyporegions.bed", col_names = FALSE, delim = "\t")
    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyperregions, ignore.strand = TRUE))), "Rintermediates/promoters_dmhyperregions.bed", col_names = FALSE, delim = "\t")

    mbo <- mergeByOverlaps(promoters, dmrsgr)


    mergeddf <- tibble(as.data.frame(mbo))
    # having to deal with elements matching multiple dmrs, and marking them as discordant in the event that the dmrs go in different directions.
    mm <- mergeddf %>%
        group_by(promoters.seqnames, promoters.start, promoters.end, promoters.width, promoters.strand, promoters.gene_id) %>%
        summarise(max_val = max(dmrsgr.diff_c2_minus_c1), min_val = min(dmrsgr.diff_c2_minus_c1))
    mm[(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "discordant"
    mm[!(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "concordant"

    mergeddf <- mm %>% mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "Discordant", ifelse(max_val > 0, "SEN Hypo", "SEN Hyper"))))

    merged <- GRanges(
        seqnames = mergeddf$promoters.seqnames,
        ranges = IRanges(start = mergeddf$promoters.start, end = mergeddf$promoters.end),
        strand = mergeddf$promoters.strand,
        gene_id = mergeddf$promoters.gene_id,
        max_val = mergeddf$max_val,
        min_val = mergeddf$min_val,
        concordance = mergeddf$concordance,
        direction = mergeddf$direction
    )

    genes_contrast_colors <- setNames(my_palette[c(1, 2, 9)], c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper"), "Discordant"))
    p <- mergeddf %>%
        group_by(direction) %>%
        summarise(n = n()) %>%
        ggplot() +
        geom_col(aes(x = direction, y = n, fill = direction), color = "black") +
        ggtitle("Promoter Methylation") +
        labs(x = "", y = "count") +
        theme(legend.position = "none") +
        mtopen + scale_contrasts +
        scale_fill_manual(values = genes_contrast_colors) +
        anchorbar

    mysave(pl = p, "results/plots/genes/genes_concordance.png", 3, 4)
    promotersplots[["genes_concordance"]] <- p

    pl_genes_density <- mergeddf %>%
        filter(concordance == "concordant") %>%
        ggplot() +
        geom_density(aes(x = max_val), fill = mycolor) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(x = "DMR Methylation (PRO-SEN)", y = "Density") +
        ggtitle("DM Promoters") +
        annotate("label", x = -Inf, y = Inf, label = "SEN Hyper", hjust = 0, vjust = 1) +
        annotate("label", x = Inf, y = Inf, label = "SEN Hypo", hjust = 1, vjust = 1) +
        theme(legend.position = "none") +
        mtopen + scale_conditions +
        anchorbar

    mysave(pl = p, "results/plots/genes/genes_density.png", 4, 4)
    promotersplots[["genes_density"]] <- p

    save(promotersplots, file = "results/plots/objects/promotersplots.rda")
    # highly_enriched_notch_sets <- tablesReactome[["Higher_in_Alz"]] %>% head(n = 15) %>% filter(grepl("NOTCH", anno.result) | grepl("LFNG", anno.result))
    # highly_enriched_notch_sets$id
    # gs[highly_enriched_notch_sets$id]

    # library(ComplexHeatmap)
    # m1 = make_comb_mat(gs[highly_enriched_notch_sets$id])
    # m1
    # p = UpSet(m1)
    # png("results/plots/genes/notch_upset_plot.png", 4, 4, units = "in", res = 300)
    # print(p)
    # dev.off()
}



## cCRES
{
    ccresplots <- list()
    # where are DMRs in cCREs?
    # cCREs
    ccresdf <- read_delim(ccrespath, col_names = FALSE)
    ccresgr <- GRanges(
        seqnames = ccresdf$X1,
        ranges = IRanges(start = ccresdf$X2, end = ccresdf$X3),
        type = ccresdf$X10,
        name = ccresdf$X4,
        closest_gene = ccresdf$X15
    )

    mbo <- mergeByOverlaps(ccresgr, dmrsgr)
    mbodf <- tibble(as.data.frame(mbo))

    p <- mbodf %>% ggplot() +
        geom_bar(aes(x = ccresgr.type, fill = direction), position = "dodge", color = "black") +
        coord_flip() +
        labs(x = "") +
        ggtitle("cCRE Methylation") +
        scale_y_continuous(expand = expansion(mult = c(0, .1))) +
        mtopen + scale_contrasts
    mysave(pl = p, "results/plots/ccres/dmrs_in_ccres.png", 5, 6)
    ccresplots[["dmrs_in_ccres"]] <- p

    totalccres <- ccresdf %>%
        group_by(X10) %>%
        summarize(n = n())
    totaldmccres <- mbodf %>%
        group_by(ccresgr.type, direction) %>%
        summarize(n = n())
    pctdmccres <- left_join(totaldmccres, totalccres, by = c("ccresgr.type" = "X10")) %>%
        mutate(pct = 100 * n.x / n.y) %>%
        mutate(myaxis = paste0(ccresgr.type, "\n", "n=", n.y)) %>%
        drop_na()

    p <- pctdmccres %>% ggplot() +
        geom_col(aes(x = myaxis, y = pct, fill = direction), position = "dodge", color = "black") +
        labs(x = "", y = "Pct Differentially Methylated") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        ggtitle("cCRE Methylation") +
        scale_y_continuous(expand = expansion(mult = c(0, .1))) +
        mtopen + scale_contrasts
    mysave(pl = p, "results/plots/ccres/dmrs_in_ccres_pct.png", 6, 6)
    ccresplots[["dmrs_in_ccres_pct"]] <- p

    save(ccresplots, file = "results/plots/objects/ccresplots.rda")
}







{ # GENES figures

    # reactomePlots <- list()
    # for (direction in directions) {
    #     p <- tablesReactome[[direction]] %>%
    #         head(n = 15) %>%
    #         mutate(anno.result = gsub("Homo sapiens: ", "", anno.result)) %>%
    #         mutate(anno.result = fct_reorder(anno.result, fold_enrichment)) %>%
    #         ggplot(aes(x = anno.result, y = fold_enrichment)) +
    #         geom_segment(aes(x = anno.result, xend = anno.result, y = 0, yend = fold_enrichment, color = p_adjust)) +
    #         geom_point(aes(color = p_adjust), size = 5) +
    #         scale_x_discrete(labels = scales::label_wrap(40)) +
    #         coord_flip() +
    #         scale_color_continuous(trans = "reverse") +
    #         labs(x = "") +
    #         theme(
    #             axis.text.y = element_text(color = "black"),
    #             panel.grid.major.y = element_blank(),
    #             panel.border = element_rect(color = "black", fill = NA, size = 1),
    #             axis.ticks.y = element_blank()
    #         )
    #     reactomePlots[[direction]] <- p
    # }

    # png("results/plots/genes/fig_genes.png", 10, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Dif_in_SEN"]] + ggtitle("Enriched Gene Sets")
    # patch <- patch + plot_annotation(tag_levels = "A")
    # print(patch)
    # dev.off()


    # png("results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched") | reactomePlots[["Higher_in_SEN"]] + ggtitle("Hyper Enrichment")
    # patch <- patch + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
    # print(patch)
    # dev.off()

    # png("results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched") | reactomePlots[["Higher_in_SEN"]] + ggtitle("SEN Hyper Enrichment")
    # patch <- patch + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
    # print(patch)
    # dev.off()



    # png("results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Lower_in_SEN"]] + ggtitle("SEN Hypo Enriched Sets") | reactomePlots[["Higher_in_SEN"]] + ggtitle("SEN Hyper Enriched Sets")
    # patch <- patch + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
    # print(patch)
    # dev.off()


#     layout <- ("
# AAACCDD
# AAACCDD
# BBBCCDD
# BBBCCDD
# GGGGGGG
# GGGGGGG
# GGGGGGG
# GGGGGGG
# ")
#     plist <- list(
#         pl_genes_bar, pl_genes_density, reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched"),
#         reactomePlots[["Higher_in_SEN"]] + ggtitle("Hyper Enriched"), pl_per_element_meth
#     )
#     patch <- wrap_plots(plist, design = layout) + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
#     png("results/plots/figs/fig_genes_both_directions_withRTE.png", 22, 20, units = "in", res = 300)
#     print(patch & theme(
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         plot.tag = element_text(size = 32)
#     ))
#     dev.off()

#     layout <- ("
# AAACCDDEEE
# AAACCDDEEE
# BBBCCDDFFF
# BBBCCDDFFF
# GGGGGGGGGG
# GGGGGGGGGG
# GGGGGGGGGG
# GGGGGGGGGG
# ")

#     plist <- list(
#         pl_genes_bar, pl_genes_density, reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched"),
#         reactomePlots[["Higher_in_SEN"]] + ggtitle("Hyper Enriched"), pl_gimme_hypo, pl_gimme_hyper, pl_per_element_meth
#     )
#     patch <- wrap_plots(plist, design = layout) + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
#     png("results/plots/figs/fig_genes_both_directions_withRTEandmotifs.png", 22, 20, units = "in", res = 300)
#     print(patch & theme(
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 18),
#         plot.tag = element_text(size = 32)
#     ))
#     dev.off()
# }


# # heatmap full element
# heatmapprep <- l1hsintactmethdf %>%
#     separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#     filter(case_when(
#         element_strand == "+" ~ (start > element_start) & (start < element_stop),
#         element_strand == "-" ~ (start > element_start) & (start < element_stop)
#     )) %>%
#     group_by(gene_id, sample) %>%
#     summarise(mean = mean(pctM)) %>%
#     pivot_wider(names_from = sample, values_from = mean)
# m <- as.matrix(heatmapprep %>% ungroup() %>% select(-gene_id))
# rownames(m) <- heatmapprep %$% gene_id
# m <- na.omit(m)
# pvals <- l1hsintactdf %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% direction
# length(pvals)
# is_sig <- !is.na(pvals)
# pch <- rep("*", length(pvals))
# pch[!is_sig] <- NA
# regions <- l1hsintactdf %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% region
# row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), region = regions, col = list(region = c("intronic" = "brown", "nonGenic" = "tan")))
# colsum <- colSums(m)
# conditions <- rep(c("control", "alz"), each = 3)
# topAnn <- ComplexHeatmap::HeatmapAnnotation(Sum = anno_barplot(colsum, axis = TRUE, axis_param = default_axis_param("column")), Condition = conditions, col = list(Condition = c("control" = "green", "alz" = "blue")))
# col_fun <- colorRamp2(c(80, 90, 100), c("red", "white", "blue"))
# col_fun(seq(80, 100, by = 10))
# heatmap <- m %>%
#     Heatmap(
#         name = "CpG Methylation",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = TRUE,
#         show_column_names = TRUE,
#         column_names_rot = 45,
#         col = col_fun,
#         top_annotation = topAnn,
#         right_annotation = row_ha,
#         row_title = "Intact L1HS"
#     )

# png("results/plots/l1intactheatmap_fullElement.png", 9, 14, units = "in", res = 300)
# draw(heatmap, heatmap_legend_side = "right")
# dev.off()



# # heatmap LTR5
# heatmapprep <- LTR5_Hsdf %>%
#     group_by(gene_id, sample) %>%
#     summarise(mean = mean(pctM)) %>%
#     pivot_wider(names_from = sample, values_from = mean)


# m <- as.matrix(heatmapprep %>% ungroup() %>% dplyr::select(-gene_id))
# rownames(m) <- heatmapprep %$% gene_id
# m <- na.omit(m)
# pvals <- LTR5_Hsdf %>%
#     group_by(gene_id, concordance) %>%
#     summarise(n = n()) %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% concordance
# length(pvals)
# is_sig <- !is.na(pvals)
# pch <- rep("*", length(pvals))
# pch[!is_sig] <- NA
# row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch))
# conditions <- rep(c("ctrl", "alz"), each = 3)
# topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = c("ctrl" = "#3C5488FF", "alz" = "#F39B7FFF")))
# col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
# col_fun(seq(0, 100, by = 25))
# heatmap <- m %>%
#     Heatmap(
#         name = "CpG Methylation",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = FALSE,
#         show_column_names = TRUE,
#         column_names_rot = 45,
#         col = col_fun,
#         top_annotation = topAnn,
#         right_annotation = row_ha,
#         row_title = "LTR5_Hs"
#     )

# png("results/plots/LTR5_Hs_heatmap.png", 9, 14, units = "in", res = 300)
# draw(heatmap, heatmap_legend_side = "right")
# dev.off()

# # heatmap LTR5A
# heatmapprep <- LTR5Adf %>%
#     group_by(gene_id, sample) %>%
#     summarise(mean = mean(pctM)) %>%
#     pivot_wider(names_from = sample, values_from = mean)


# m <- as.matrix(heatmapprep %>% ungroup() %>% dplyr::select(-gene_id))
# rownames(m) <- heatmapprep %$% gene_id
# m <- na.omit(m)
# pvals <- LTR5Adf %>%
#     group_by(gene_id, concordance) %>%
#     summarise(n = n()) %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% concordance
# length(pvals)
# is_sig <- !is.na(pvals)
# pch <- rep("*", length(pvals))
# pch[!is_sig] <- NA
# row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch))
# conditions <- rep(c("ctrl", "alz"), each = 3)
# topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = c("ctrl" = "#3C5488FF", "alz" = "#F39B7FFF")))
# col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
# col_fun(seq(0, 100, by = 25))
# heatmap <- m %>%
#     Heatmap(
#         name = "CpG Methylation",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = FALSE,
#         show_column_names = TRUE,
#         column_names_rot = 45,
#         col = col_fun,
#         top_annotation = topAnn,
#         right_annotation = row_ha,
#         row_title = "LTR5A"
#     )

# png("results/plots/LTR5A_heatmap.png", 9, 14, units = "in", res = 300)
# draw(heatmap, heatmap_legend_side = "right")
# dev.off()


# # heatmap LTR5B
# heatmapprep <- LTR5Bdf %>%
#     group_by(gene_id, sample) %>%
#     summarise(mean = mean(pctM)) %>%
#     pivot_wider(names_from = sample, values_from = mean)


# m <- as.matrix(heatmapprep %>% ungroup() %>% dplyr::select(-gene_id))
# rownames(m) <- heatmapprep %$% gene_id
# m <- na.omit(m)
# pvals <- LTR5Bdf %>%
#     group_by(gene_id, concordance) %>%
#     summarise(n = n()) %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% concordance
# length(pvals)
# is_sig <- !is.na(pvals)
# pch <- rep("*", length(pvals))
# pch[!is_sig] <- NA
# row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch))
# conditions <- rep(c("ctrl", "alz"), each = 3)
# topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = c("ctrl" = "#3C5488FF", "alz" = "#F39B7FFF")))
# col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
# col_fun(seq(0, 100, by = 25))
# heatmap <- m %>%
#     Heatmap(
#         name = "CpG Methylation",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = FALSE,
#         show_column_names = TRUE,
#         column_names_rot = 45,
#         col = col_fun,
#         top_annotation = topAnn,
#         right_annotation = row_ha,
#         row_title = "LTR5B"
#     )

# png("results/plots/LTR5B_heatmap.png", 9, 14, units = "in", res = 300)
# draw(heatmap, heatmap_legend_side = "right")
# dev.off()





# # heatmap LTR5
# heatmapprep <- LTR5df %>%
#     group_by(gene_id, sample) %>%
#     summarise(mean = mean(pctM)) %>%
#     pivot_wider(names_from = sample, values_from = mean)


# m <- as.matrix(heatmapprep %>% ungroup() %>% dplyr::select(-gene_id))
# rownames(m) <- heatmapprep %$% gene_id
# m <- na.omit(m)
# pvals <- LTR5df %>%
#     group_by(gene_id, concordance) %>%
#     summarise(n = n()) %>%
#     arrange(gene_id) %>%
#     filter(gene_id %in% rownames(m)) %$% concordance
# length(pvals)
# is_sig <- !is.na(pvals)
# pch <- rep("*", length(pvals))
# pch[!is_sig] <- NA
# row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch))
# conditions <- rep(c("ctrl", "alz"), each = 3)
# topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = c("ctrl" = "#3C5488FF", "alz" = "#F39B7FFF")))
# col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
# col_fun(seq(0, 100, by = 25))
# heatmap <- m %>%
#     Heatmap(
#         name = "CpG Methylation",
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         show_row_names = FALSE,
#         show_column_names = TRUE,
#         column_names_rot = 45,
#         col = col_fun,
#         top_annotation = topAnn,
#         right_annotation = row_ha,
#         row_title = "LTR5"
#     )

# png("results/plots/LTR5_heatmap.png", 9, 14, units = "in", res = 300)
# draw(heatmap, heatmap_legend_side = "right")
# dev.off()



# ##############################
# # cCREs
# dir.create("results/plots/ccres")
# ccrespath <- "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations3/cCREs/hs1-imr90-cCREsCuratedWithClosestGene.bed"
# refseqgenespath <- "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations2/ncbiRefSeqGenes.gtf"

# ccresdf <- read_delim(ccrespath, col_names = FALSE)
# ccresgr <- GRanges(
#     seqnames = ccresdf$X1,
#     ranges = IRanges(start = ccresdf$X2, end = ccresdf$X3),
#     type = ccresdf$X10,
#     name = ccresdf$X4,
#     closest_gene = ccresdf$X15
# )


# mbo <- mergeByOverlaps(grs, ccresgr)
# ccresmeth <- mbo$grs
# ccresmeth$type <- mbo$type
# ccresmethdf <- tibble(as.data.frame(ccresmeth))
# p <- grsdf %>% ggplot() +
#     geom_histogram(aes(x = pctM))
# png("results/plots/ccres/test.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()

# # group Av
# pal <- paletteer::paletteer_d("LaCroixColoR::Coconut", n = 5)
# palr <- rep(pal, times = c(1, 2, 2, 2, 2))
# p <- ccresmethdf %>%
#     filter(cov > 4) %>%
#     ggplot() +
#     geom_boxplot(aes(x = type, y = pctM, fill = condition), outlier.shape = NA) +
#     theme_cowplot() +
#     scale_fill_manual(values = palr) +
#     theme(aspect.ratio = 0.33) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     xlab("") +
#     ylab("CpG Fraction Methylated") +
#     ggtitle("", ) +
#     theme(plot.title = element_text(hjust = 0.5))
# png("results/plots/ccres/boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()





###################

# ## CENTROMERE
# # ANNOTATIONS
# cdr <- import.bed(cdrpath)

# HORdf <- read_delim(HORpath)
# HOR <- GRanges(
#     seqnames = HORdf$chr,
#     ranges = IRanges(start = HORdf$start, end = HORdf$end)
# )
# HOR$type <- HORdf[[1]]
# activeHOR <- HOR[HOR$type == "Active", ]





# cdrmeth <- subsetByOverlaps(grs, cdr)
# cdrmethdf <- tibble(as.data.frame(cdrmeth))
# cdrmethdf %>%
#     group_by(condition, seqnames) %>%
#     summarise(avMeth = mean(pctM), avCov = mean(cov))

# # group Av
# p <- cdrmethdf %>%
#     filter(cov > 4) %>%
#     ggplot() +
#     geom_boxplot(aes(x = seqnames, y = pctM, fill = condition), outlier.shape = NA) +
#     theme_cowplot() +
#     scale_fill_manual(values = c("#C4CED4FF", "#E13A3EFF")) +
#     theme(aspect.ratio = 0.33) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     xlab("") +
#     ylab("CpG Fraction Methylated") +
#     ggtitle("") +
#     theme(plot.title = element_text(hjust = 0.5))
# png("results/plots/centromere/boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()

# # group rolling av
# pf <- cdrmethdf %>%
#     filter(cov > 4) %>%
#     group_by(seqnames, condition) %>%
#     mutate(rM = rollmean(pctM, 50, na.pad = TRUE, align = "center")) %>%
#     filter(!is.na(rM)) %>%
#     ungroup()
# p <- pf %>% ggplot() +
#     geom_line(aes(x = start, y = rM, color = condition)) +
#     scale_color_manual(values = c("#484848", "#E13A3EFF")) +
#     scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#     facet_wrap(~seqnames, scales = "free_x")
# png("results/plots/centromere/cdrLines.png", 12, 12, units = "in", res = 300)
# print(p)
# dev.off()

# ###########################################

# cenRegion <- reduce(cenSat)
# lengths(cenRegion) %>% sum()

# cenRegionMeth <- subsetByOverlaps(grs, cenRegion)
# cenRegionMethdf <- tibble(as.data.frame(cenRegionMeth))

# # group Av
# p <- cenRegionMethdf %>%
#     filter(cov > 4) %>%
#     ggplot() +
#     geom_boxplot(aes(x = seqnames, y = pctM, fill = condition), outlier.shape = NA) +
#     theme_cowplot() +
#     scale_fill_manual(values = c("#C4CED4FF", "#E13A3EFF")) +
#     theme(aspect.ratio = 0.33) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     xlab("") +
#     ylab("CpG Fraction Methylated") +
#     ggtitle("", ) +
#     theme(plot.title = element_text(hjust = 0.5))
# png("results/plots/centromere/cenRegion_boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()

# # group rolling av
# pf <- cenRegionMethdf %>%
#     filter(cov > 4) %>%
#     group_by(seqnames, condition) %>%
#     mutate(rM = rollmean(pctM, 100, na.pad = TRUE, align = "center")) %>%
#     filter(!is.na(rM)) %>%
#     ungroup()
# p <- pf %>% ggplot() +
#     geom_line(aes(x = start, y = rM, color = condition), alpha = 0.5) +
#     scale_color_manual(values = c("#484848", "#E13A3EFF")) +
#     scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#     geom_rect(data = tibble(as.data.frame(activeHOR)), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "yellow", alpha = 0.3) +
#     facet_wrap(~seqnames, scales = "free_x")
# png("results/plots/centromere/cenRegion_Lines.png", 12, 12, units = "in", res = 300)
# p
# dev.off()


# ###########

# activeHOR %>% reduce()


# activeHORMeth <- subsetByOverlaps(grs, activeHOR)
# activeHORMethdf <- tibble(as.data.frame(activeHORMeth))

# # group Av
# p <- activeHORMethdf %>%
#     filter(cov > 4) %>%
#     ggplot() +
#     geom_boxplot(aes(x = seqnames, y = pctM, fill = condition), outlier.shape = NA) +
#     theme_cowplot() +
#     scale_fill_manual(values = c("#C4CED4FF", "#E13A3EFF")) +
#     theme(aspect.ratio = 0.33) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     xlab("") +
#     ylab("CpG Fraction Methylated") +
#     ggtitle("", ) +
#     theme(plot.title = element_text(hjust = 0.5))
# png("results/plots/centromere/activeHOR_boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()

# # group rolling av
# pf <- activeHORMethdf %>%
#     filter(cov > 4) %>%
#     group_by(seqnames, condition) %>%
#     mutate(rM = rollmean(pctM, 50, na.pad = TRUE, align = "center")) %>%
#     filter(!is.na(rM)) %>%
#     ungroup()
# p <- pf %>% ggplot() +
#     geom_line(aes(x = start, y = rM, color = condition), alpha = 0.5) +
#     scale_color_manual(values = c("#484848", "#E13A3EFF")) +
#     scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#     geom_rect(data = tibble(as.data.frame(cdr)), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "yellow", alpha = 0.3) +
#     facet_wrap(~seqnames, scales = "free_x")
# png("results/plots/centromere/activeHOR_Lines.png", 12, 12, units = "in", res = 300)
# print(p)
# dev.off()



# ####################




# cenSatMeth <- subsetByOverlaps(grs, cenSat)
# cenSatMethOL <- findOverlaps(grs, cenSat)
# typesOL <- cenSat$types[cenSatMethOL@to]
# cenSatMeth$types <- typesOL
# cenSatMethdf <- tibble(as.data.frame(cenSatMeth))
# cenSatMethdf %>%
#     group_by(types, condition) %>%
#     filter(cov > 4) %>%
#     summarise(avMeth = mean(pctM))


# sample_size <- cenSatMethdf %>%
#     group_by(types) %>%
#     summarize(num = n())
# # Plot

# p <- cenSatMethdf %>%
#     left_join(sample_size) %>%
#     mutate(myaxis = paste0(types, "\n", "n=", num)) %>%
#     filter(cov > 4) %>%
#     ggplot() +
#     geom_boxplot(aes(x = myaxis, y = pctM, fill = condition), outlier.shape = NA) +
#     theme_cowplot() +
#     scale_fill_manual(values = c("#C4CED4FF", "#E13A3EFF")) +
#     theme(aspect.ratio = 0.33) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     xlab("") +
#     ylab("CpG Fraction Methylated") +
#     ggtitle("", ) +
#     theme(plot.title = element_text(hjust = 0.5))
# png("results/plots/centromere/cenSatTypes_boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()
