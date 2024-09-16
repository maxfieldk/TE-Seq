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
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/outfiles/bedmethylanalysis.txt",
            promoters_bed = "ldna/Rintermediates/promoters.bed",
            dmrpromoterhyper_bed = "ldna/Rintermediates/promoters_dmhyperregions.bed",
            dmrpromoterhypo_bed = "ldna/Rintermediates/promoters_dmhyporegions.bed"
        ), env = globalenv())
    }
)

dmlspath <- inputs$dmls
dmrspath <- inputs$dmrs

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)

grsdf <- read_delim("ldna/Rintermediates/grsdf.tsv", col_names = TRUE)
grsdf %$% sample %>% unique()
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

RMdf <- read_delim("ldna/Rintermediates/RMdf.tsv", col_names = TRUE)
rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf <- read_delim("ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)
flRTEpromoter <- read_delim("ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
flRTEpromoter <- read_delim("ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
perelementdf_promoters <- read_delim("ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
reads <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
readscg <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
dmrs <- read_delim(inputs$dmrs, delim = "\t", col_names = TRUE)
dmls <- read_delim(inputs$dmls, delim = "\t", col_names = TRUE)

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

############
# GLOBAL
dir.create("ldna/results/plots/genomewide", showWarnings = FALSE)
a <- grsdf %>%
    group_by(sample) %>%
    summarize(mean = mean(pctM), sd = sd(pctM), n = n())
# t.test(pctM ~ condition, data = grsdf, var.equal = TRUE)
# t.test(pctM ~ condition, data = grsdf %>% filter(seqnames %in% chromosomesNoX), var.equal = TRUE)

p <- grsdfs %>% ggplot() +
    geom_boxplot(aes(x = islandStatus, y = pctM, fill = condition)) +
    mtopen +
    scale_conditions
mysaveandstore(fn = "ldna/results/plots/genomewide/cpgislandstatusbox.pdf", w = 6, h = 5, res = 300, pl = p)


p <- grsdfs %>%
    group_by(islandStatus, condition) %>%
    summarize(pctM = mean(pctM)) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = pctM, fill = condition), position = "dodge", color = "black") +
    mtopen +
    scale_conditions +
    anchorbar
mysaveandstore(fn = "ldna/results/plots/genomewide/cpgislandstatusbar_1000.pdf", w = 4, h = 4, res = 300, pl = p)

##################################### DMR analysis
# how many dmrs
p <- dmrs %>%
    ggplot() +
    geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
    labs(x = "") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ggtitle("DMR Counts") +
    mtopen +
    scale_methylation
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_number.pdf", 4, 4)

# what is their average length
p <- ggplot(data = dmrs) +
    geom_histogram(aes(length), fill = mycolor, color = "black") +
    ggtitle("DMR Lengths") +
    labs(x = "length (bp)") +
    xlim(0, 3000) +
    mtopen +
    anchorbar
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_length.pdf", w = 4, h = 4)

p <- ggplot(data = dmrs) +
    geom_histogram(aes(length, fill = direction), alpha = 0.7, color = "black") +
    ggtitle("DMR Lengths") +
    labs(x = "length (bp)") +
    xlim(0, 3000) +
    mtopen +
    scale_methylation +
    anchorbar
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_length_stratified.pdf", w = 4, h = 4)

# a positive difference means that sample 1 has more methylation than sample 2
meandiff <- mean(dmrs$diff_c2_minus_c1)

p <- ggplot() +
    geom_density(
        data = dmrs,
        aes(x = diff_c2_minus_c1), fill = mycolor
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = sprintf("DMR Methylation (%s-%s)", condition2, condition1), y = "Density") +
    ggtitle("DMR Methylation Density") +
    annotate("label", x = -Inf, y = Inf, label = "Hypo", hjust = 0, vjust = 1) +
    annotate("label", x = Inf, y = Inf, label = "Hyper", hjust = 1, vjust = 1) +
    mtopen
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_delta.pdf", w = 4, h = 4)

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
    mtopen +
    scale_methylation +
    anchorbar +
    labs(x = "", y = "Count") +
    ggtitle("DMR") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_count_islandstatus.pdf", 4, 4)


dmrlocdf <- dmrs %>%
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
mysaveandstore(fn = "ldna/results/plots/genomewide/dmr_count.pdf", 5, 5)

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
    mtclosed
mysaveandstore(fn = "ldna/results/plots/genomewide/dmrdensity.pdf", 5, 5)

p <- dmls %>%
    ggplot() +
    geom_bar(aes(x = direction, fill = direction), show.legend = FALSE, color = "black") +
    ggtitle("DML Counts") +
    labs(x = "", y = "Count") +
    anchorbar +
    mtclosed +
    scale_methylation
mysaveandstore(fn = "ldna/results/plots/genomewide/dml_count.pdf", 4, 4)
# dmls
p <- ggplot() +
    geom_density(data = dmls, aes(x = diff_c2_minus_c1), fill = mycolor) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = sprintf("DML Methylation (%s-%s)", condition2, condition1), y = "Density") +
    ggtitle("DML Methylation Density") +
    annotate("label", x = -Inf, y = Inf, label = "Hypo", hjust = 0, vjust = 1) +
    annotate("label", x = Inf, y = Inf, label = "Hyper", hjust = 1, vjust = 1) +
    mtopen +
    scale_contrasts +
    anchorbar

mysaveandstore(fn = "ldna/results/plots/genomewide/dml_delta.pdf", 4, 4)

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
    mtclosed
mysaveandstore(fn = "ldna/results/plots/genomewide/dmldensity.pdf", w = 5, h = 5)



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
    mtopen +
    scale_methylation +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "Count") +
    ggtitle("DML Island Status")
mysaveandstore(fn = "ldna/results/plots/genomewide/dml_count_islandstatus.pdf", 4, 4)




dir.create("ldna/results/plots/figs")
pff <- flRTEpromoter %>%
    group_by(rte_subfamily, genic_loc) %>%
    mutate(group_n = n()) %>%
    group_by(rte_subfamily, direction, genic_loc) %>%
    summarise(n = n(), group_n = dplyr::first(group_n)) %>%
    mutate(frac_dm = n / group_n) %>%
    filter(!is.na(rte_subfamily)) %>%
    filter(!is.na(direction)) %>%
    filter(direction != "discordant") %>%
    ungroup() %>%
    complete(rte_subfamily, direction, genic_loc, fill = list(n = 0, group_n = 0, frac_dm = 0))
numdf <- flRTEpromoter %>%
    group_by(rte_subfamily, genic_loc) %>%
    summarise(group_n_accurate = n())
p <- pff %>%
    left_join(numdf) %>%
    mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
    ggplot() +
    geom_col(aes(x = ann_axis, y = frac_dm, fill = direction), position = "dodge", color = "black") +
    facet_wrap(~genic_loc, scales = "free_x", nrow = 2) +
    labs(x = "", y = "Fraction Differentially Methylated") +
    ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
    mtclosed +
    scale_methylation
mysaveandstore(sprintf("ldna/results/plots/rte/dmfl%s_promoter_regionstrat.pdf", "all"), 12, 7)

pff <- flRTEpromoter %>%
    group_by(rte_subfamily) %>%
    mutate(group_n = n()) %>%
    group_by(rte_subfamily, direction) %>%
    summarise(n = n(), group_n = dplyr::first(group_n)) %>%
    mutate(frac_dm = n / group_n) %>%
    filter(!is.na(rte_subfamily)) %>%
    filter(!is.na(direction)) %>%
    filter(direction != "discordant") %>%
    ungroup() %>%
    complete(rte_subfamily, direction, fill = list(n = 0, group_n = 0, frac_dm = 0))
numdf <- flRTEpromoter %>%
    group_by(rte_subfamily) %>%
    summarise(group_n_accurate = n())
p <- pff %>%
    left_join(numdf) %>%
    mutate(ann_axis = paste0(rte_subfamily, "\n", "n=", group_n_accurate)) %>%
    ggplot() +
    geom_col(aes(x = ann_axis, y = frac_dm, fill = direction), position = "dodge", color = "black") +
    labs(x = "", y = "Fraction Differentially Methylated") +
    ggtitle(sprintf("Full Length %s Promoter Differential Methylation", "RTE")) +
    mtclosed +
    scale_methylation
mysaveandstore(sprintf("ldna/results/plots/rte/dmfl%s_promoter.pdf", "all"), 12, 4)

pf <- perelementdf_promoters
p <- pf %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_subfamily_n, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_subfamily_n, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
stats <- pf %>%
    group_by(sample, condition, rte_subfamily) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup() %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters.pdf", 14, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters.pdf", raster = TRUE, 14, 6)

p <- pf %>%
    filter(!grepl("HERVL", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_subfamily_n, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_subfamily_n, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
stats <- pf %>%
    group_by(sample, condition, rte_subfamily) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup() %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_1.pdf", 14, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_1.pdf", raster = TRUE, 14, 6)

p <- pf %>%
    filter(grepl("^L1", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot(aes(x = rte_subfamily_n, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_pwc(aes(group = condition)) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions
stats <- pf %>%
    group_by(sample, condition, rte_subfamily) %>%
    summarise(mean_meth = mean(mean_meth)) %>%
    ungroup() %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_L1s.pdf", 14, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_L1s.pdf", raster = TRUE, 12, 4)



p <- pf %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily) %>%
    mutate(n = n()) %>%
    mutate(rte_subfamily_n = paste0(rte_subfamily, "\nn=", n)) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = mean_meth, color = condition)) +
    geom_quasirandom(dodge.width = 0.75) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions

stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_L1s.pdf", 10, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_L1s.pdf", raster = TRUE, 12, 4)


p <- pf %>%
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
    ggtitle("RTE CpG Methylation") +
    mtopen +
    scale_conditions

stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_intactL1s.pdf", 10, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_L1s.pdf", raster = TRUE, 12, 4)

p <- pf %>%
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

stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", 10, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_braak_intactL1s.pdf", raster = TRUE, 12, 4)


p <- pf %>%
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
pf %>%
    filter(grepl("^L1HS", rte_subfamily)) %>%
    group_by(sample) %>%
    summarize(median = median(mean_meth)) %>%
    left_join(sample_table %>% dplyr::rename(sample = sample_name)) %>%
    ungroup() %>%
    group_by(condition) %>%
    summarize(mean_of_median = mean(median))
stats <- pf %>%
    compare_means(mean_meth ~ condition, data = ., group.by = "rte_subfamily", p.adjust.method = "fdr")
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_braak_L1s.pdf", 10, 6, sf = stats)
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters_by_sample_braak_L1s.pdf", raster = TRUE, 12, 4)



pfl1 <- pf %>% filter(grepl("^L1", rte_subfamily))
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
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_paired_promoters_L1s.pdf", 14, 6, raster = TRUE)

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
        mysaveandstore(fn = "ldna/results/plots/rte/l1hs_boxplot_promoters.pdf", 5, 4, sf = stats)
    },
    error = function(e) {
        mysaveandstore(fn = "ldna/results/plots/rte/l1hs_boxplot_promoters.pdf", 5, 4)
    }
)

#################

l1hsintactmethgr <- rtedf %>%
    filter(intactness_req == "Intact")
l1hsintactmethgr <- l1hsintactmethgr %>%
    mutate(rel_start = start - rte_start) %>%
    mutate(rel_end = end - rte_start)
write_delim(l1hsintactmethgr, "ldna/Rintermediates/l1hsintactdf.tsv", col_names = TRUE)

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

mysaveandstore("ldna/results/plots/rte/l1intact_Lines_pos_strand.pdf", 12, 30)


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
mysaveandstore("ldna/results/plots/rte/l1intact_Lines_pos_strand_promoter.pdf", 12, 30)

p <- pf_neg %>% ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtclosed +
    scale_conditions
mysaveandstore("ldna/results/plots/rte/l1intact_Lines_neg_strand.pdf", 12, 30)

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
mysaveandstore("ldna/results/plots/rte/l1intact_Lines_neg_strand_promoter.pdf", 12, 30)




# chr4_83273990_83280141_+
# chr6_44705376_44711532_+
# for the fig


element_anatomy <- read_delim("aref/default/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")

dm_intact_l1hs_elements <- flRTEpromoter %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(intactness_req == "Intact") %>%
    filter(direction == "Hypo")
dm_fl_l1hs_elements <- flRTEpromoter %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(direction == "Hypo")
topmovers_l1hs_elements <- flRTEpromoter %>%
    filter(gene_id %in% top_l1hs_movers_intact)
element_sets_of_interst <- list("dm_fl_l1hs" = dm_fl_l1hs_elements, "dm_intact_l1hs" = dm_intact_l1hs_elements, "Top_Movers" = topmovers_l1hs_elements)

library(ggnewscale)
library(patchwork)
l1hsflmethgr <- rtedf %>%
    filter(rte_subfamily == "L1HS")
l1hsflmethgr <- l1hsflmethgr %>%
    mutate(rel_start = start - rte_start) %>%
    mutate(rel_end = end - rte_start)
write_delim(l1hsflmethgr, "ldna/Rintermediates/l1hsfldf.tsv", col_names = TRUE)

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
    dir.create("ldna/Rintermediates/l1hs/", recursive = TRUE)
    write_delim(df %>% dplyr::select(gene_id), sprintf("ldna/Rintermediates/l1hs/%s_gene_id.tsv", element_type), col_names = FALSE)
    write_delim(df %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/l1hs/%s_promoters.bed", element_type), col_names = FALSE, delim = "\t")
    write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ] %>% dplyr::select(seqnames, start, end, strand, gene_id), sprintf("ldna/Rintermediates/l1hs/%s_full_elements.bed", element_type), col_names = FALSE, delim = "\t")
    write_delim(RMdf[match(df %$% gene_id, RMdf$gene_id), ], sprintf("ldna/Rintermediates/l1hs/%s_full_elements.tsv", element_type), col_names = TRUE, delim = "\t")

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
                geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
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

        mysaveandstore(sprintf("ldna/results/plots/rte/%s/%s_methylation.pdf", element_type, element), 5, 5)
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

        mysaveandstore(sprintf("ldna/results/plots/rte/%s/%s_methylation_conditionaveraged.pdf", element_type, element), 5, 5)
    }
}




# rtedf %$% ltr_viral_status %>% unique()


# LTR5Adf <- rtedf %>%
#     filter(ltr_viral_status == "5'LTR (FL Int)") %>%
#     separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#     filter(element_stop - element_start > 600) %>%
#     filter(cov > MINIMUMCOVERAGE)

# { # now LTR5s
#     LTR5Adf <- rtedf %>%
#         filter(type == "LTR5A") %>%
#         separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#         filter(element_stop - element_start > 600) %>%
#         filter(cov > MINIMUMCOVERAGE)
#     LTR5Bdf <- rtedf %>%
#         filter(type == "LTR5B") %>%
#         separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#         filter(element_stop - element_start > 600) %>%
#         filter(cov > MINIMUMCOVERAGE)
#     LTR5_Hsdf <- rtedf %>%
#         filter(type == "LTR5_Hs") %>%
#         separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#         filter(element_stop - element_start > 600) %>%
#         filter(cov > MINIMUMCOVERAGE)
#     LTR5df <- rtedf %>%
#         filter(type == "LTR5") %>%
#         separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#         filter(element_stop - element_start > 600) %>%
#         filter(cov > MINIMUMCOVERAGE)

#     pf <- LTR5_Hsdf %>%
#         filter(cov > MINIMUMCOVERAGE) %>%
#         group_by(gene_id, condition) %>%
#         mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
#         filter(!is.na(rM)) %>%
#         ungroup()
#     p <- pf %>% ggplot() +
#         geom_line(aes(x = start, y = rM, color = condition)) +
#         scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#         facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
#         ylim(c(0, 100)) +
#         mtopen + scale_conditions
#     png("ldna/results/plots/rte/LTR5_Hs_Lines.png", 12, 120, units = "in", res = 200)
#     print(p)
#     dev.off()

#     pf <- LTR5Bdf %>%
#         filter(cov > MINIMUMCOVERAGE) %>%
#         group_by(gene_id, condition) %>%
#         mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
#         filter(!is.na(rM)) %>%
#         ungroup()
#     p <- pf %>% ggplot() +
#         geom_line(aes(x = start, y = rM, color = condition)) +
#         scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#         facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
#         ylim(c(60, 100)) +
#         mtopen + scale_conditions
#     png("ldna/results/plots/rte/LTR5B_Lines.png", 12, 120, units = "in", res = 200)
#     print(p)
#     dev.off()

#     pf <- LTR5df %>%
#         filter(cov > MINIMUMCOVERAGE) %>%
#         group_by(gene_id, condition) %>%
#         mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
#         filter(!is.na(rM)) %>%
#         ungroup()
#     p <- pf %>% ggplot() +
#         geom_line(aes(x = start, y = rM, color = condition)) +
#         scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#         facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
#         ylim(c(60, 100)) +
#         mtopen + scale_conditions
#     png("ldna/results/plots/rte/LTR5_Lines.png", 12, 120, units = "in", res = 200)
#     print(p)
#     dev.off()
# }


# ####



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
            dir.create("ldna/results/plots/rte/l1hsintact")
            png(paste0("ldna/results/plots/rte/l1hsintact/", gene_id, ".png"), 8, 3, units = "in", res = 300)
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
library(circlize)
col_fun <- colorRamp2(c(50, 75, 100), c("red", "white", "blue"))
col_fun(seq(50, 100, by = 12.5))
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

p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR, heatmap_legend_side = "right", annotation_legend_side = "right")))
mysaveandstore("ldna/results/plots/l1intactheatmap_5utr.pdf", 7, 14)

col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
col_fun(seq(50, 100, by = 12.5))
heatmapL1UTR <- m %>%
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

p <- wrap_elements(grid.grabExpr(draw(heatmapL1UTR, heatmap_legend_side = "right", annotation_legend_side = "right")))
mysaveandstore("ldna/results/plots/l1intactheatmap_5utr_fullrange.pdf", 7, 14)
# plgrob <- grid.grabExpr(ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right"))
# plots[["l1intactheatmap_5utr"]] <- plgrob



read_analysis <- function(
    readsdf,
    region = "L1HS_intactness_req_ALL",
    mod_code_var = "m",
    context = "CG") {
    readsdf1 <- readsdf %>% left_join(r_annotation_fragmentsjoined %>% dplyr::select(gene_id, start, end, strand) %>% dplyr::rename(element_strand = strand, element_start = start, element_end = end))
    utr1 <- readsdf1 %>%
        filter(mod_code == mod_code_var) %>%
        filter(case_when(
            element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
            element_strand == "-" ~ (start > element_end - 909) & (start < element_end)
        )) %>%
        dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

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

    write_delim(utr, sprintf("ldna/Rintermediates/%s_%s_%s_reads.tsv", region, mod_code_var, context), delim = "\t")

    p <- utr %>% ggplot() +
        geom_density(aes(x = mod_qual, fill = condition), alpha = 0.3) +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/modbase_score_dist_%s_%s_%s.pdf", region, mod_code_var, context), 12, 4, pl = p)


    p <- utr %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/read_span_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)


    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(num_cpgs_in_read), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/read_num_cpg_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)



    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(read_id, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_point(aes(x = read_span, y = fraction_meth, color = condition)) +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)


    aa <- utr %>%
        filter(read_span > 600) %>%
        group_by(read_id, sample, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ungroup()

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
        geom_pwc(aes(x = condition, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "", y = "Pct Reads < 50% methylated") +
        ggtitle("5'UTR Methylation") +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/plots/reads/barplot_50pct_%s_%s_%s.pdf", region, mod_code_var, context), 4, 4, pl = p)




    p <- aa %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        labs(x = "Pct CpG Methylation per Read") +
        ggtitle("5'UTR Methylation") +
        facet_wrap(vars(region)) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_density_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 4, 4, pl = p)


    p <- utr %>%
        filter(region == "L1HS_intactness_req_ALL") %>%
        filter(read_span > 600) %>%
        group_by(read_id, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        labs(x = "Pct CpG Methylation per Read") +
        ggtitle("L1HS Intact 5'UTR Methylation") +
        mtopen +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_density_distribution_l1hsintact_%s_%s.pdf", mod_code_var, context), 4, 4, pl = p)
}

tryCatch(
    {
        read_analysis(readscg, "L1HS_intactness_req_ALL", "m", "CpG")
        read_analysis(reads, "L1HS_intactness_req_ALL", "m", "NoContext")
        # read_analysis(readscg, "L1HS_intactness_req_ALL", "h", "CpG")
        # read_analysis(reads, "L1HS_intactness_req_ALL", "h", "NoContext")
        # read_analysis(reads, "L1HS_intactness_req_ALL", "a", "NoContext")
    },
    error = function(e) {
        print(e)
    }
)



######### GENES


{
    directions <- c("Hypo", "Hyper", "Dif")
    mydir <- "ldna/results/plots/great"
    mydirtables <- "ldna/results/tables/great"
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

library(clusterProfiler)
library(scales)
library(msigdbr)

{
    gs <- msigdbr("human")
    tablesMsigdb <- list()
    results <- list()
    genecollections <- gs$gs_cat %>% unique()
    for (collection in genecollections) {
        tryCatch({
            dir.create(paste(mydir, collection, sep = "/"))
            dir.create(paste(mydirtables, collection, sep = "/"))
            genesets <- gs %>%
                filter(gs_cat == collection) %>%
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
            }
        })
    }
    save(file = "ldna/Rintermediates/tablesMsigdb.rds", tablesMsigdb)
    for (collection in genecollections) {
        tryCatch({
            dir.create(paste(mydir, collection, sep = "/"))
            dir.create(paste(mydirtables, collection, sep = "/"))
            for (direction in directions) {
                p <- tablesMsigdb[[collection]][[direction]] %>%
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
                mysaveandstore(paste(mydir, collection, paste0(direction, "lollipop.pdf"), sep = "/"), 7, 7)
                if (collection == "msigdbC5_GO") {
                    p <- tablesMsigdb[[collection]][[direction]] %>%
                        mutate(goterm = str_extract(id, "GO[a-zA-Z][a-zA-Z]")) %>%
                        mutate(id = gsub("GO[a-zA-Z][a-zA-Z]_", "", id)) %>%
                        mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                        mutate(id = fct_reorder(id, fold_enrichment)) %>%
                        group_by(goterm) %>%
                        head(n = 10) %>%
                        ggplot(aes(x = id, y = fold_enrichment)) +
                        geom_col(aes(fill = p_adjust), color = "black") +
                        facet_grid(goterm ~ ., scales = "free", space = "free") +
                        coord_flip() +
                        scale_color_continuous(trans = "reverse") +
                        labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                        mtclosed +
                        anchorbar
                    mysaveandstore(paste(mydir, collection, paste0(direction, "lollipop_faceted.pdf"), sep = "/"), 7, 7)
                }
            }
        })
    }
}




{
    promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
    write_delim(tibble(as.data.frame(promoters)) %>% mutate(score = 1000) %>% dplyr::select(seqnames, start, end, gene_id, score, strand), "ldna/Rintermediates/promoters.bed", col_names = FALSE, delim = "\t")

    hyporegions <- dmrsgr[grepl("Hypo", dmrsgr$direction)]
    hyperregions <- dmrsgr[grepl("Hyper", dmrsgr$direction)]

    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, dmrsgr, ignore.strand = TRUE))), outputs$promoters_bed, col_names = FALSE, delim = "\t")
    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyporegions, ignore.strand = TRUE))), outputs$dmrpromoterhypo_bed, col_names = FALSE, delim = "\t")
    write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyperregions, ignore.strand = TRUE))), outputs$dmrpromoterhyper_bed, col_names = FALSE, delim = "\t")

    mbo <- mergeByOverlaps(promoters, dmrsgr)


    mergeddf <- tibble(as.data.frame(mbo))
    # having to deal with elements matching multiple dmrs, and marking them as discordant in the event that the dmrs go in different directions.
    mm <- mergeddf %>%
        group_by(promoters.seqnames, promoters.start, promoters.end, promoters.width, promoters.strand, promoters.gene_id) %>%
        summarise(max_val = max(dmrsgr.diff_c2_minus_c1), min_val = min(dmrsgr.diff_c2_minus_c1))
    mm[(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "discordant"
    mm[!(mm$max_val > 0 & mm$min_val < 0), "concordance"] <- "concordant"

    mergeddf <- mm %>% mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "Discordant", ifelse(max_val > 0, "Hyper", "Hypo"))))

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

    p <- mergeddf %>%
        group_by(direction) %>%
        summarise(n = n()) %>%
        ggplot() +
        geom_col(aes(x = direction, y = n, fill = direction), color = "black") +
        ggtitle("Promoter Methylation") +
        labs(x = "", y = "count") +
        theme(legend.position = "none") +
        mtopen +
        scale_methylation +
        anchorbar

    mysaveandstore(pl = p, "ldna/results/plots/genes/genes_concordance.pdf", 4, 4)

    p <- mergeddf %>%
        filter(concordance == "concordant") %>%
        ggplot() +
        geom_density(aes(x = max_val), fill = mycolor) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(x = sprintf("DMR Methylation (%s-%s)", condition2, condition1), y = "Density") +
        ggtitle("DM Promoters") +
        annotate("label", x = -Inf, y = Inf, label = "Hyper", hjust = 0, vjust = 1) +
        annotate("label", x = Inf, y = Inf, label = "Hypo", hjust = 1, vjust = 1) +
        theme(legend.position = "none") +
        mtopen +
        anchorbar

    mysaveandstore(pl = p, "ldna/results/plots/genes/genes_density.pdf", 4, 4)

    # highly_enriched_notch_sets <- tablesReactome[["Higher_in_Alz"]] %>% head(n = 15) %>% filter(grepl("NOTCH", anno.result) | grepl("LFNG", anno.result))
    # highly_enriched_notch_sets$id
    # gs[highly_enriched_notch_sets$id]

    # library(ComplexHeatmap)
    # m1 = make_comb_mat(gs[highly_enriched_notch_sets$id])
    # m1
    # p = UpSet(m1)
    # png("ldna/results/plots/genes/notch_upset_plot.png", 4, 4, units = "in", res = 300)
    # print(p)
    # dev.off()
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
        closest_gene = ccresdf$X15
    )

    chromHMMgr <- import(conf$chromHMM)

    annotations_of_interest <- list(chromHMM = chromHMMgr, cCREs = ccresgr)
    # note I should be testing for enrichment of these sets too
    for (annotation_of_interest in names(annotations_of_interest)) {
        annot <- annotations_of_interest[[annotation_of_interest]]
        annotdf <- as.data.frame(annot) %>% tibble()
        mbo <- mergeByOverlaps(annot, dmrsgr)
        mbodf <- tibble(as.data.frame(mbo))

        p <- mbodf %>% ggplot() +
            geom_bar(aes(x = annot.name, fill = direction), position = "dodge", color = "black") +
            labs(x = "") +
            coord_flip() +
            ggtitle("cCRE Methylation") +
            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
            mtopen +
            scale_methylation
        mysaveandstore(pl = p, sprintf("ldna/results/plots/%s/dmrs_in_%s.pdf", annotation_of_interest, annotation_of_interest), 5, 6)

        total <- annotdf %>%
            group_by(name) %>%
            summarize(n = n())
        totaldm <- mbodf %>%
            group_by(annot.name, direction) %>%
            summarize(n = n())
        pctdm <- left_join(totaldm, total, by = c("annot.name" = "name")) %>%
            mutate(pct = 100 * n.x / n.y) %>%
            mutate(myaxis = paste0(annot.name, "\n", "n=", n.y)) %>%
            drop_na()

        p <- pctdm %>% ggplot() +
            geom_col(aes(x = myaxis, y = pct, fill = direction), position = "dodge", color = "black") +
            labs(x = "", y = "Pct Differentially Methylated") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
            coord_flip() +
            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
            mtopen +
            scale_methylation
        mysaveandstore(pl = p, sprintf("ldna/results/plots/%s/dmrs_in_%s_pct.pdf", annotation_of_interest, annotation_of_interest), 6, 6)
        p <- pctdm %>% ggplot() +
            geom_col(aes(x = annot.name, y = pct, fill = direction), position = "dodge", color = "black") +
            labs(x = "", y = "Pct Differentially Methylated") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            ggtitle(sprintf("%s Methylation", annotation_of_interest)) +
            coord_flip() +
            scale_y_continuous(expand = expansion(mult = c(0, .1))) +
            mtopen +
            scale_methylation
        mysaveandstore(pl = p, sprintf("ldna/results/plots/%s/dmrs_in_%s_pct_clean.pdf", annotation_of_interest, annotation_of_interest), 6, 6)
    }
}

library(regioneR)





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

    # png("ldna/results/plots/genes/fig_genes.png", 10, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Dif_in_SEN"]] + ggtitle("Enriched Gene Sets")
    # patch <- patch + plot_annotation(tag_levels = "A")
    # print(patch)
    # dev.off()


    # png("ldna/results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched") | reactomePlots[["Higher_in_SEN"]] + ggtitle("Hyper Enrichment")
    # patch <- patch + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
    # print(patch)
    # dev.off()

    # png("ldna/results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
    # patch <- (pl_genes_bar / pl_genes_density) | reactomePlots[["Lower_in_SEN"]] + ggtitle("Hypo Enriched") | reactomePlots[["Higher_in_SEN"]] + ggtitle("SEN Hyper Enrichment")
    # patch <- patch + plot_annotation(tag_levels = "A") + plot_layout(ggene_ides = "collect")
    # print(patch)
    # dev.off()



    # png("ldna/results/plots/genes/fig_genes_both_directions.png", 16, 7, units = "in", res = 300)
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
    #     png("ldna/results/plots/figs/fig_genes_both_directions_withRTE.png", 22, 20, units = "in", res = 300)
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
    #     png("ldna/results/plots/figs/fig_genes_both_directions_withRTEandmotifs.png", 22, 20, units = "in", res = 300)
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
    # row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), region = regions, col = list(region = c("intronic" = "brown", "Intergenic" = "tan")))
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

    # png("ldna/results/plots/l1intactheatmap_fullElement.png", 9, 14, units = "in", res = 300)
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

    # png("ldna/results/plots/LTR5_Hs_heatmap.png", 9, 14, units = "in", res = 300)
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

    # png("ldna/results/plots/LTR5A_heatmap.png", 9, 14, units = "in", res = 300)
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

    # png("ldna/results/plots/LTR5B_heatmap.png", 9, 14, units = "in", res = 300)
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

    # png("ldna/results/plots/LTR5_heatmap.png", 9, 14, units = "in", res = 300)
    # draw(heatmap, heatmap_legend_side = "right")
    # dev.off()



    # ##############################
    # # cCREs
    # dir.create("ldna/results/plots/ccres")
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
    # png("ldna/results/plots/ccres/test.png", 12, 6, units = "in", res = 300)
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
    # png("ldna/results/plots/ccres/boxplot.png", 12, 6, units = "in", res = 300)
    # print(p)
    # dev.off()
}


##################

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
mysaveandstore("ldna/results/plots/centromere/cenRegion_boxplot.pdf", 12, 5)

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
mysaveandstore("ldna/results/plots/centromere/cenRegion_Lines.pdf", 12, 12)
cenRegionMethdf %>% filter(seqnames == "chr22")
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
mysaveandstore("ldna/results/plots/centromere/censatTypes_boxplot.pdf", 12, 6)
