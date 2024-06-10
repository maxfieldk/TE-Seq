module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
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
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv",
            read_mods = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis


####
# RUN IF RESUMING
if (interactive()) {
    conditions <- conf$levels
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

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
}

##########################
if (!interactive()) {
    sample_grs <- list()
    for (sample_name in samples) {
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
    grsdf %$% seqnames %>% unique()
    dir.create("ldna/Rintermediates", showWarnings = FALSE)
    write_delim(grsdf %>% filter(grepl("*nonref*", seqnames)), "ldna/Rintermediates/grsdf_nonref.tsv", col_names = TRUE)
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

    dir.create("ldna/Rintermediates", showWarnings = FALSE)
    write_delim(grsdf, "ldna/Rintermediates/grsdf.tsv", col_names = TRUE)

    # SETTING UP SOME SUBSETS FOR EXPLORATION
    set.seed(75)
    grsdfs <- grsdf %>%
        group_by(sample, seqnames, islandStatus) %>%
        slice_sample(n = 1000)
    grss <- GRanges(grsdfs)
    write_delim(grsdfs, "ldna/Rintermediates/grsdfsmall.tsv", col_names = TRUE)
}



############
# GLOBAL
dir.create("ldna/results/plots/genomewide", showWarnings = FALSE)
a <- grsdf %>%
    group_by(sample) %>%
    summarize(mean = mean(pctM), sd = sd(pctM), n = n())

p <- grsdfs %>% ggplot() +
    geom_boxplot(aes(x = islandStatus, y = pctM, fill = condition)) +
    labs(x = "CpG") +
    mtopen + scale_conditions
mysaveandstore(fn = "ldna/results/plots/genomewide/cpgislandstatusbox.pdf", w = 4, h = 4, res = 300, pl = p)


p <- grsdfs %>%
    group_by(islandStatus, condition) %>%
    summarize(pctM = mean(pctM)) %>%
    ggplot() +
    geom_col(aes(x = islandStatus, y = pctM, fill = condition), position = "dodge", color = "black") +
    labs(x = "CpG") +
    mtopen + scale_conditions +
    anchorbar
mysaveandstore(fn = "ldna/results/plots/genomewide/cpgislandstatusbar_1000.pdf", w = 4, h = 4, res = 300, pl = p)


dir.create("ldna/results/plots/figs")

# patch1 <- (pl_meth_by_chromosome)
# patch2 <- (pl_ndml + pl_dml_raster + pl_cdml)
# patch3 <- (pl_ndmr + pl_dmr_raster + pl_cdmr)
# patch <- patch1 + patch2 + patch3 + plot_layout(nrow = 3, ggene_ides = "collect") + plot_annotation(tag_levels = "A")
# png("ldna/results/plots/figs/global.png", height = 20, width = 20, res = 300, units = "in")
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
# mysaveandstore(fn = "ldna/results/plots/figs/global_2.pdf", w = 20, h = 20, res = 300, pl = p)


# rm(pl_meth_by_chromosome, patch, patch1, patch2, patch3, pl_ndml, pl_ddml, pl_cdml, pl_ndmr, pl_ddmr, pl_cdmr)
##########



####################
## RTEs
dir.create("ldna/results/plots/rte", showWarnings = FALSE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
r_repeatmasker_annotation %$% ltr_viral_status_req %>% unique()
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
RM <- GRanges(rmann)
### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        big_ontology_groups <- c(big_ontology_groups, rmann %>%
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

rmann <- rmann %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)


grouping_var <- "rte_length_req"
classes <- rmann %>%
    pull(!!sym(grouping_var)) %>%
    unique() %>%
    na.omit()
classes <- classes[!str_detect(classes, "Other")]
rtegrl <- GRangesList()
for (rte in classes) {
    print(rte)
    mbo <- mergeByOverlaps(grs, GRanges(rmann %>% filter(!!sym(grouping_var) == rte)))
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

write_delim(rtedf, "ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
# rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf <- rtedf %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, rte_length_req, type, l1_intactness_req, ltr_viral_status_req) %>%
    summarize(mean_meth = mean(pctM))

perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))
# mask <- perelementdf$condition == "healthy"
# newC <- ifelse(mask, "ctrl", "alz")
# perelementdf$condition <- newC

write_delim(perelementdf, "ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)
# perelementdf <- read_delim("ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)


# PROMOTERS
# annotate whether full length elements promoters overlap DMRs
rmann %$% ltr_viral_status_req %>% unique()
flelement <- rmann %>% filter(str_detect(rte_length_req, ">"))
flSINE <- flelement %>% filter(rte_superfamily == "SINE")
flLINE <- flelement %>% filter(rte_superfamily == "LINE")
flFl_Provirus_5LTR <- flelement %>%
    filter(str_detect(gene_id, "LTR")) %>%
    filter(ltr_viral_status_req == "Fl_Provirus_5LTR")

flSINEgrs <- GRanges(flSINE)
flLINE5UTRgrs <- GRanges(flLINE) %>% resize(909)
flFl_Provirus_5LTRgrs <- GRanges(flFl_Provirus_5LTR)
flRTEpromotergrs <- c(c(flSINEgrs, flLINE5UTRgrs), flFl_Provirus_5LTRgrs)


flRTEpromoter <- flRTEpromotergrs %>% as.data.frame() %>% tibble() %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
flRTEpromoter %$% rte_subfamily %>% unique()

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

write_delim(rtedf_promoters, "ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
# rtedf_promoters <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, rte_length_req, type, l1_intactness_req, ltr_viral_status_req) %>%
    summarize(mean_meth = mean(pctM))

perelementdf_promoters <- perelementdf_promoters %>% filter(!is.na(rte_length_req))
# mask <- perelementdf$condition == "healthy"
# newC <- ifelse(mask, "ctrl", "alz")
# perelementdf$condition <- newC

write_delim(perelementdf_promoters, "ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
perelementdf <- read_delim("ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)
perelementdf_promoters <- read_delim("ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)

library(ggpubr)
perelementdf_promoters %$% type %>% table()
p <- perelementdf_promoters %>% filter(type == "L1HS >6kb") %>%
    gghistogram(x = "mean_meth",add = "mean", rug = TRUE)
mysaveandstore(fn = "ldna/results/plots/rte/l1hs_density_promoters.pdf", 4, 4)
perelementdf_promoters
lt25methl1 <- perelementdf_promoters %>% filter(type == "L1HS >6kb") %>% filter(mean_meth <= 25) %>% dplyr::select(gene_id) %>% mutate(MethStatus = "demethylated") %>% left_join(rmann) %>% mutate(position_string = paste0(seqnames, ":", start, "-", end))
mt75methl1 <- perelementdf_promoters %>% filter(type == "L1HS >6kb") %>% filter(mean_meth >= 75) %>% dplyr::select(gene_id) %>% mutate(MethStatus = "methylated") %>% left_join(rmann) %>% mutate(position_string = paste0(seqnames, ":", start, "-", end))
rbind(lt25methl1, mt75methl1) %>% dplyr::select(position_string, MethStatus) %>% dplyr::rename(loc = position_string, cluster = MethStatus) %>% write_delim("ldna/Rintermediates/l1_dif_meth_for_motif_analysis.tsv", col_names = TRUE, delim = "\t")
library(Rsamtools)
fa <- FaFile(conf$reference)
lt25methl1_ss <- getSeq(fa, lt25methl1 %>% GRanges())
names(lt25methl1_ss) <- lt25methl1$gene_id
mt75methl1_ss <- getSeq(fa, mt75methl1 %>% GRanges())
names(mt75methl1_ss) <- mt75methl1$gene_id
writeXStringSet(lt25methl1_ss, "ldna/Rintermediates/lt25methl1_ss.fa")
writeXStringSet(mt75methl1_ss, "ldna/Rintermediates/mt75methl1_ss.fa")


p <- perelementdf_promoters %>%
    ggplot() +
    geom_quasirandom(aes(x = rte_length_req, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = rte_length_req, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen + scale_conditions
mysaveandstore(fn = "ldna/results/plots/rte/repmasker_boxplot_promoters.pdf", 14, 6)

p <- perelementdf_promoters %>%
    filter(l1_intactness_req != "Other") %>%
    ggplot() +
    geom_quasirandom(aes(x = l1_intactness_req, y = mean_meth, color = condition), dodge.width = 0.75) +
    geom_boxplot(aes(x = l1_intactness_req, y = mean_meth, color = condition), alpha = 0.5, outlier.shape = NA) +
    xlab("") +
    ylab("Average CpG Methylation Per Element") +
    ggtitle("RTE CpG Methylation") +
    mtopen + scale_conditions
mysaveandstore(fn = "ldna/results/plots/rte/l1_intact_boxplot_promoters.pdf", 14, 6)

#################


l1hsintactmethgr <- rtedf %>%
    filter(str_detect(l1_intactness_req, "Intact")) %>%
    left_join(r_annotation_fragmentsjoined %>% dplyr::select(gene_id, start, end, strand) %>% dplyr::rename(element_strand = strand, element_start = start, element_end = end))
l1hsintactmethgr <- l1hsintactmethgr %>%
    mutate(rel_start = start - element_start) %>%
    mutate(rel_end = end - element_start)
write_delim(l1hsintactmethgr, "ldna/Rintermediates/l1hsintactdf.tsv", col_names = TRUE)

library(zoo)
pf_pos <- l1hsintactmethgr %>%
    filter(element_strand == "+") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()
pf_neg <- l1hsintactmethgr %>%
    filter(element_strand == "-") %>%
    as.data.frame() %>%
    tibble() %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, condition) %>%
    mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
    filter(!is.na(rM)) %>%
    ungroup()

p <- pf_pos %>% ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtclosed + scale_conditions
png("ldna/results/plots/rte/l1intact_Lines_pos_strand.png", 12, 30, units = "in", res = 300)
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
    mtclosed + scale_conditions
png("ldna/results/plots/rte/l1intact_Lines_pos_strand_promoter.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

p <- pf_neg %>% ggplot() +
    geom_line(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    ylim(c(0, 100)) +
    mtclosed + scale_conditions
png("ldna/results/plots/rte/l1intact_Lines_neg_strand.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

p <- pf_neg %>%
    filter(rel_start < 910) %>%
    ggplot() +
    geom_point(aes(x = rel_start, y = rM, color = condition)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
    xlim(c(1, 910)) +
    ylim(c(0, 100)) +
    mtclosed + scale_conditions
png("ldna/results/plots/rte/l1intact_Lines_neg_strand_promoter.png", 12, 30, units = "in", res = 300)
print(p)
dev.off()

# chr4_83273990_83280141_+
# chr6_44705376_44711532_+
# for the fig
element_anatomy <- read_delim("aref/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")

elements_of_interest <- c("L1HS_3q23_1", "L1HS_5p13.1_1")
### nice L1 annot
library(ggnewscale)
library(patchwork)
for (element in elements_of_interest) {
    if (rmann %>% filter(gene_id == element) %$% strand == "+") {
    modifier <- rmann %>% filter(gene_id == element) %$% start
    color_intervals <- element_anatomy %>% filter(!(feature %in% c("EN", "RT"))) %>% filter(gene_id == element) %>% mutate(across(where(is.numeric), ~ . + modifier))
    p1 <- pf_pos %>%
        filter(gene_id  == element) %>%
        ggplot() +
        geom_line(aes(x = start, y = rM, fill = sample, color = sample)) +
        labs(y = "Methylation Rolling Mean") +
        mtclosed + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

    } else {
    modifier <- rmann %>% filter(gene_id == element) %$% end
    color_intervals <- element_anatomy %>% filter(!(feature %in% c("EN", "RT"))) %>% filter(gene_id == element) %>% mutate(across(where(is.numeric), ~ modifier - .))
    p1 <- pf_neg %>%
        filter(gene_id  == element) %>%
        ggplot() +
        geom_line(aes(x = start, y = rM, fill = sample, color = sample)) +
        labs(y = "Methylation Rolling Mean") +
        mtclosed + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
    }
    p2 <- color_intervals %>% 
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
        geom_text(aes(x = (start + end) / 2, y = 1.5, label = feature)) +
        coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
        ggtitle(element) +
        scale_fill_paletteer_d("dutchmasters::milkmaid") + mtclosed + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(expand = c(0, 0.4)) + theme(legend.position = "none")
    
    p <- p2 / p1 + plot_layout(heights = c(0.2, 1))

    mysaveandstore(sprintf("ldna/results/plots/rte/%s_methylation_line.pdf", element), 5, 5)
}

for (element in elements_of_interest) {
    if (rmann %>% filter(gene_id == element) %$% strand == "+") {
    modifier <- rmann %>% filter(gene_id == element) %$% start
    color_intervals <- element_anatomy %>% filter(!(feature %in% c("EN", "RT"))) %>% filter(gene_id == element) %>% mutate(across(where(is.numeric), ~ . + modifier))
    p1 <- pf_pos %>%
        filter(gene_id  == element) %>%
        ggplot() +
        geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
        labs(y = "Methylation Rolling Mean") +
        mtclosed + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

    } else {
    modifier <- rmann %>% filter(gene_id == element) %$% end
    color_intervals <- element_anatomy %>% filter(!(feature %in% c("EN", "RT"))) %>% filter(gene_id == element) %>% mutate(across(where(is.numeric), ~ modifier - .))
    p1 <- pf_neg %>%
        filter(gene_id  == element) %>%
        ggplot() +
        geom_point(aes(x = start, y = rM, fill = sample, color = sample)) +
        labs(y = "Methylation Rolling Mean") +
        mtclosed + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

    }


    p2 <- color_intervals %>% 
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
        geom_text(aes(x = (start + end) / 2, y = 1.5, label = feature)) +
        coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
        ggtitle(element) +
        scale_fill_paletteer_d("dutchmasters::milkmaid") + mtclosed + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
        scale_y_continuous(expand = c(0, 0.4)) + theme(legend.position = "none")
    
    p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
    
    mysaveandstore(sprintf("ldna/results/plots/rte/%s_methylation.pdf", element), 5, 5)

}


rtedf %$% ltr_viral_status_req %>% unique()


#     LTR5Adf <- rtedf %>%
#         filter(ltr_viral_status_req == "Fl_Provirus_5LTR") %>%
#         separate(gene_id, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
#         filter(element_stop - element_start > 600) %>%
#         filter(cov > MINIMUMCOVERAGE)

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



# l1hsintactmethdf <- l1hsintactmethgr %>%
#     as.data.frame() %>%
#     tibble()
# for (gene_id in l1hsintactmethdf %$% gene_id %>% unique()) {
#     pf <- l1hsintactmethdf %>%
#         filter(gene_id == gene_id) %>%
#         filter(cov > MINIMUMCOVERAGE) %>%
#         group_by(sample) %>%
#         mutate(rM = rollmean(pctM, 15, na.pad = TRUE, align = "center")) %>%
#         filter(!is.na(rM)) %>%
#         ungroup()
#     p <- pf %>% ggplot() +
#         geom_line(aes(x = start, y = rM, color = condition)) +
#         scale_x_continuous(breaks = scales::breaks_pretty(3)) +
#         facet_wrap(~gene_id, ncol = 5, scales = "free_x") +
#         ylim(c(0, 100)) +
#         mtopen + scale_conditions
#     dir.create("ldna/results/plots/rte/l1hsintact")
#     png(paste0("ldna/results/plots/rte/l1hsintact/", gene_id, ".png"), 8, 3, units = "in", res = 300)
#     print(p)
#     dev.off()
# }



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
m <- m[!remove_rows, ] %>% as.matrix()

rownames(m)
l1hsintactdf <- l1hsintactmethdf %>% group_by(gene_id) %>% summarise(genic_loc = dplyr::first(genic_loc))

genic_locs <- l1hsintactdf %>%
    arrange(gene_id) %>%
    filter(gene_id %in% rownames(m)) %$% genic_loc
row_ha <- rowAnnotation(genic_loc = genic_locs, col = list(genic_loc = c("Genic" = "brown", "Intergenic" = "tan")))
conditions <- c(sample_table %>% filter(condition == condition1) %$% condition, sample_table %>% filter(condition == condition2) %$% condition)

conditions <- sample_table[match(colnames(m), sample_table$sample_name),]$condition
topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))
library(circlize)
col_fun <- colorRamp2(c(0, 50, 100), c("red", "white", "blue"))
col_fun(seq(0, 100, by = 12.5))
heatmapL1UTR <- m %>%
    Heatmap(
        name = "CpG Methylation",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_rot = 45,
        col = col_fun,
        top_annotation = topAnn,
        right_annotation = row_ha,
        row_title = "Intact L1HS"
    )

png("ldna/results/plots/l1intactheatmap_5utr.png", 5, 14, units = "in", res = 300)
ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right")
dev.off()

# plgrob <- grid.grabExpr(ComplexHeatmap::draw(heatmapL1UTR, heatmap_legend_side = "right"))
# plots[["l1intactheatmap_5utr"]] <- plgrob

############## Read density
rm(grs)
mem_used()

dir.create("ldna/results/plots/reads")
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
        eoi <- import(paste0(conf$reference_annotation_dir, "/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
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
        eoi <- import(paste0(conf$reference_annotation_dir, "/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
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
# utr1 %>% filter(mod_code == "m") %>%
#         mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0)) %$% mod_indicator %>% table()


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

    write_delim(utr, sprintf("ldna/Rintermediates/%s_%s_%s_reads.tsv", region, mod_code_var, context), delim = "\t")

    p <- utr %>% ggplot() +
        geom_density(aes(x = mod_qual, fill = condition), alpha = 0.3) +
        facet_wrap(vars(region)) + 
        mtclosed + scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/modbase_score_dist_%s_%s_%s.pdf", region, mod_code_var, context), 12, 4, pl = p)


    p <- utr %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region)) + 
        mtclosed + scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/read_span_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)


    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(gene_id, read_id, condition, region) %>%
        summarise(nc = max(num_cpgs_in_read), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = nc, fill = condition), alpha = 0.5) +
        facet_wrap(vars(region)) +
        mtclosed + scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/read_num_cpg_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)



    p <- utr %>%
        filter(read_span > 250) %>%
        group_by(read_id, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_point(aes(x = read_span, y = fraction_meth, color = condition)) +
        facet_wrap(vars(region)) +
        mtclosed + scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 5, 5, pl = p)


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
        labs(x = "", y = "Pct Reads < 50% methylated") +
        ggtitle("5'UTR Methylation") +
        mtclosedgridh + scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/plots/reads/barplot_50pct_%s_%s_%s.pdf", region, mod_code_var, context), 4, 4, pl = p)

    p <- aa %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        labs(x = "Pct CpG Methylation per Read") +
        ggtitle("5'UTR Methylation") +
        facet_wrap(vars(region)) +
        mtclosed + scale_conditions
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_density_distribution_%s_%s_%s.pdf", region, mod_code_var, context), 4, 4, pl = p)

    p <- utr %>%
        filter(region == "L1HS_l1_intactness_req_ALL") %>%
        filter(read_span > 600) %>%
        group_by(read_id, condition, region) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), read_span = max(read_span), strand = dplyr::first(element_strand)) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        labs(x = "Pct CpG Methylation per Read") +
        ggtitle("L1HS Intact 5'UTR Methylation") +
        mtopen + scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/plots/reads/fraction_meth_density_distribution_l1hsintact_%s_%s.pdf", mod_code_var, context), 4, 4, pl = p)


    # note that this isn't the right test for this

}

tryCatch({
    read_analysis(readscg, "L1HS_l1_intactness_req_ALL", "m", "CpG")
    read_analysis(reads, "L1HS_l1_intactness_req_ALL", "m", "NoContext")
    read_analysis(readscg, "L1HS_l1_intactness_req_ALL", "h", "CpG")
    read_analysis(reads, "L1HS_l1_intactness_req_ALL", "h", "NoContext")
    read_analysis(reads, "L1HS_l1_intactness_req_ALL", "a", "NoContext")
}, error = function(e) {
    print(e)
})




# ##################

# ## CENTROMERE
# # ANNOTATIONS
# cdr <- import.bed(conf$cdr)

# HORdf <- read_delim(conf$HOR)
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
# png("ldna/results/plots/centromere/boxplot.png", 12, 6, units = "in", res = 300)
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
# png("ldna/results/plots/centromere/cdrLines.png", 12, 12, units = "in", res = 300)
# print(p)
# dev.off()

# # ###########################################
# cenSat <- import.bed(conf$censat)
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
# png("ldna/results/plots/centromere/cenRegion_boxplot.png", 12, 6, units = "in", res = 300)
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
# png("ldna/results/plots/centromere/cenRegion_Lines.png", 12, 12, units = "in", res = 300)
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
# png("ldna/results/plots/centromere/activeHOR_boxplot.png", 12, 6, units = "in", res = 300)
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
# png("ldna/results/plots/centromere/activeHOR_Lines.png", 12, 12, units = "in", res = 300)
# print(p)
# dev.off()



# ####################



# grs <- GRanges(grsdf)
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
# png("ldna/results/plots/centromere/cenSatTypes_boxplot.png", 12, 6, units = "in", res = 300)
# print(p)
# dev.off()
