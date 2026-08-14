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


###########################

cpg_islands <- rtracklayer::import(conf$cpg_islands)
cpgi_shores <- rtracklayer::import(conf$cpgi_shores)
cpgi_shelves <- rtracklayer::import(conf$cpgi_shelves)
cpgi_features <- c(cpg_islands, cpgi_shelves, cpgi_shores)




refseq_gr <- import(conf$refseq_unaltered)
genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
mcols(genes_gr) %>% colnames()
mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)

grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
grsdf %$% sample %>% unique()
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
grs <- GRanges(grsdf)

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
possample <- sample(grsdf$pos, size = 1000, replace = FALSE)
possample <- grsdf %>%
    group_by(islandStatus) %>%
    slice_sample(n = 1000) %>%
    ungroup() %>%
    pull(pos)
grsdfs <- grsdf %>% filter(pos %in% possample)
grss <- GRanges(grsdfs)


grsdf %>%
    ungroup() %>%
    group_by(pos, condition) %>%
    summarise(variance = var(pctM)) %>%
    group_by(condition) %>%
    summarise(meanvar = mean(variance))



grsdfs %>%
    ungroup() %>%
    dplyr::select(pos, pctM, sample) %>%
    pivot_wider(id_cols = pos, names_from = sample, values_from = pctM) %>%
    na.omit()
############
# GLOBAL
dir.create(sprintf("ldna/results/%s/plots/genomewide", params$mod_code), showWarnings = FALSE)

grsdfsummary <- grsdf %>%
    mutate(methylated_sites = round(cov * pctM / 100)) %>%
    group_by(sample) %>%
    summarise(total_cov = sum(cov), total_meth = sum(methylated_sites))

pf <- grsdfsummary %>%
    mutate(mean_meth = total_meth / total_cov) %>%
    left_join(sample_table)
pf %>%
    group_by(condition) %>%
    summarise(mean_meth = mean(mean_meth))

p <- pf %>%
    ggplot(aes(y = sample, x = mean_meth, color = condition)) +
    geom_point(size = 3) +
    scale_conditions +
    geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
    # new_scale_fill() +
    # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
    scale_conditions +
    mtopen
if ("apoe" %in% names(pf)) {
    p <- p + geom_text_repel(aes(label = apoe))
}

library(broom)
if (enough_samples_per_condition_for_stats) {
    # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
    genomewide_model <- glmmTMB(
        cbind(total_meth, total_cov - total_meth) ~
            condition + sex + age_z + (1 | sample),
        data = grsdfsummary %>% left_join(sample_table %>%
            mutate(age_z = as.numeric(scale(age)))),
        family = binomial()
    )
    stats <- broom::tidy(genomewide_model)
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
} else {
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4)
}

tryCatch(
    {
        p <- pf %>%
            mutate(sample = fct_reorder(paste0(sample, "_", age), age)) %>%
            ggplot(aes(y = sample, x = mean_meth, color = condition, shape = !!sym(asc))) +
            geom_point(size = 3) +
            scale_conditions +
            geom_vline(xintercept = pf %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
            geom_vline(xintercept = pf %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
            geom_text_repel(aes(label = apoe)) +
            # new_scale_fill() +
            # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
            scale_conditions +
            mtopen
        library(broom)
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_point_withage_ordered_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
    },
    error = function(e) {

    }
)

# t.test(pctM ~ condition, data = grsdf, var.equal = TRUE)
# t.test(pctM ~ condition, data = grsdf %>% filter(seqnames %in% chromosomesNoX), var.equal = TRUE)

p <- grsdfs %>% ggplot() +
    geom_boxplot(aes(x = islandStatus, y = pctM, fill = condition), outliers = FALSE) +
    mtclosed +
    scale_conditions
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/cpgislandstatusbox.pdf", params$mod_code), w = 4, h = 4, res = 300, pl = p)


p <- grsdfs %>%
    group_by(islandStatus, condition) %>%
    summarize(pctM = mean(pctM)) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = pctM, fill = condition), position = "dodge", color = "black") +
    mtclosed +
    scale_conditions +
    anchorbar
mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/cpgislandstatusbar_1000.pdf", params$mod_code), w = 4, h = 4, res = 300, pl = p)





#########
for (sample in confALL$ldna$samples) {
    mean_cov <- read_delim(str_glue("ldna/intermediates/{sample}/coverage/analysis_default/{sample}.sorted.bg.mosdepth.mosdepth.summary.txt"))
    coverage <- rtracklayer::import(str_glue("ldna/intermediates/{sample}/coverage/analysis_default/{sample}.sorted.bg.mosdepth.regions.bed.gz"))


    mcols(coverage)$coverage <- as.numeric(mcols(coverage)$name)
    sample_fa <- Rsamtools::FaFile(sprintf("aref/extended/%s.fa", "A.REF"))

    # Define 100kb bins across the genome
    bins100k <- tileGenome(seqlengths(sample_fa),
        tilewidth = 100000,
        cut.last.tile.in.chrom = TRUE
    )

    hits <- findOverlaps(coverage, bins100k)

    # Aggregate coverage: average coverage of 1kb bins per 100kb bin
    mean_coverage <- tapply(mcols(coverage)$coverage[queryHits(hits)], subjectHits(hits), mean)

    # Store average coverage in the 100kb bins
    mcols(bins100k)$mean_coverage <- NA
    mcols(bins100k)$mean_coverage[as.integer(names(mean_coverage))] <- mean_coverage

    # Define 100kb bins across the genome
    bins10k <- tileGenome(seqlengths(sample_fa),
        tilewidth = 10000,
        cut.last.tile.in.chrom = TRUE
    )

    hits <- findOverlaps(coverage, bins10k)

    # Aggregate coverage: average coverage of 1kb bins per 100kb bin
    mean_coverage <- tapply(mcols(coverage)$coverage[queryHits(hits)], subjectHits(hits), mean)

    # Store average coverage in the 100kb bins
    mcols(bins10k)$mean_coverage <- NA
    mcols(bins10k)$mean_coverage[as.integer(names(mean_coverage))] <- mean_coverage


    segdupsdf <- read_delim("/users/mkelsey/data/Nanopore/alz/segdups.bed", col_names = TRUE, delim = "\t")
    segdups <- segdupsdf %>%
        dplyr::rename(seqnames = `#chrom`, start = chromStart, end = chromEnd) %>%
        GRanges()
    cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
    cytobands <- cytobandsdf %>%
        dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
        GRanges()
    centromere <- cytobandsdf %>%
        filter(X5 == "acen") %>%
        dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
        GRanges()
    telomeredf <- read_delim(confALL$aref$ref_telomere, col_names = FALSE, delim = "\t")
    telomere <- telomeredf %>%
        dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
        GRanges()
    mappability <- import("/users/mkelsey/data/Nanopore/alz/k50.Unique.Mappability.bb")



    bins100k_norep <- bins100k %>%
        subsetByOverlaps(centromere, invert = TRUE) %>%
        subsetByOverlaps(segdups, invert = TRUE) %>%
        subsetByOverlaps(telomere, invert = TRUE)
    df <- bins100k_norep %>%
        as.data.frame() %>%
        tibble() %>%
        mutate(refstatus = if_else(seqnames %in% nonrefchromosomes, "NonRef", "Ref")) %>%
        filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS_REF) %>%
        mutate(seqnames = factor(seqnames, levels = CHROMOSOMESINCLUDEDINANALYSIS_REF))
    df <- bins100k %>%
        as.data.frame() %>%
        tibble() %>%
        mutate(refstatus = if_else(seqnames %in% nonrefchromosomes, "NonRef", "Ref")) %>%
        filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS_REF) %>%
        mutate(seqnames = factor(seqnames, levels = CHROMOSOMESINCLUDEDINANALYSIS_REF))


    global_mean <- df %$% mean_coverage %>% mean()
    p <- df %>% ggplot() +
        geom_point(aes(x = start, y = mean_coverage, alpha = 0.2)) +
        facet_wrap(~seqnames, scales = "free") +
        geom_hline(yintercept = global_mean, color = "darkgreen", linewidth = 1.5) +
        ylim(c(0, 90)) +
        mtclosedgridh +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    mysaveandstore(fn = str_glue("RTE/ldna/results/m/plots/{sample}_coverage.pdf"), w = 12, h = 12, raster = TRUE)

    p <- df %>% ggplot() +
        geom_point(aes(x = start, y = mean_coverage, alpha = 0.2)) +
        facet_wrap(~seqnames, scales = "free") +
        ylim(c(0, 200)) +
        mtclosedgridh +
        theme(axis.text.x = element_blank())
    mysaveandstore(fn = str_glue("RTE/ldna/results/m/plots/{sample}_coverage_notext.pdf"), w = 12, h = 12, raster = TRUE)


    perchrom_mean <- df %>%
        group_by(seqnames) %>%
        summarise(mean_coverage = mean(mean_coverage))
    p <- df %>% ggplot() +
        geom_point(aes(x = start, y = mean_coverage, alpha = 0.2)) +
        facet_wrap(~seqnames, scales = "free") +
        geom_hline(data = perchrom_mean, aes(yintercept = mean_coverage), color = "darkgreen", linewidth = 1.5) +
        geom_hline(yintercept = global_mean, color = "blue", linewidth = 1.5) +
        ylim(c(0, 90)) +
        mtclosedgridh +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    mysaveandstore(fn = str_glue("RTE/ldna/results/m/plots/{sample}_coverage_individualmeans.pdf"), w = 12, h = 12, raster = TRUE)

    p <- df %>% ggplot() +
        geom_point(aes(x = start, y = mean_coverage, alpha = 0.2)) +
        facet_wrap(~seqnames, scales = "free") +
        geom_hline(data = perchrom_mean, aes(yintercept = mean_coverage), color = "darkgreen", linewidth = 1.5) +
        geom_hline(yintercept = global_mean, color = "blue", linewidth = 1.5) +
        ylim(c(0, 90)) +
        mtclosedgridh +
        theme(axis.text.x = element_blank())
    mysaveandstore(fn = str_glue("RTE/ldna/results/m/plots/{sample}_coverage_individualmeans_notext.pdf"), w = 12, h = 12, raster = TRUE)

    ####

    bins10k %>%
        as.data.frame() %$% mean_coverage %>%
        quantile(probs = seq(0, 1, 0.05), na.rm = TRUE)
    lowcovbins <- bins10k %>%
        as.data.frame() %>%
        tibble() %>%
        filter(mean_coverage < 5)
    lowcovbinsgrs <- GRanges(lowcovbins)

    lowcovbins %$% seqnames %>% table()

    lowcovbins %>% GRanges()

    lowcovgenes <- genes_gr %>%
        subsetByOverlaps(lowcovbinsgrs) %>%
        as.data.frame() %>%
        tibble()
}

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "1genomewide", params$mod_code))
