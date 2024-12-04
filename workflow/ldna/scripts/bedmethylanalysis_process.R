module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

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
            bedmethylpaths = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_bedMethyl.bed", samples, samples),
            data = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv",
            read_mods = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

dmlspath <- inputs$dmls
dmrspath <- inputs$dmrs

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis


##########################
# PREP DATA FOR ANALYSIS
sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE), col_names = FALSE)
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
dir.create("ldna/Rintermediates", recursive = TRUE)
write_delim(grsdf %>% filter(grepl("^NI", seqnames)), "ldna/Rintermediates/grsdf_nonref.tsv", col_names = TRUE)
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

if (conf$single_condition == "no") {
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



####################
## RTEs
dir.create("ldna/results/plots/rte", showWarnings = FALSE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
RM <- GRanges(rmann)
### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "_.*family")]
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
# annotate whether repeats overlap DMRs

if (conf$single_condition == "no") {
    conditions <- conf$levels
    condition1 <- conditions[1]
    condition2 <- conditions[2]
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
    rmannfinal <- tibble(as.data.frame(RMfinal))
    rmannGOOD <- rmannfinal %>%
        mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))))
    RMdf <- rmannGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
} else {
    RMdf <- rmann %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
}
write_delim(RMdf, "ldna/Rintermediates/RMdf.tsv", col_names = TRUE)
# RMdf <- read_delim("ldna/Rintermediates/RMdf.tsv", col_names = TRUE)


grouping_var <- "rte_subfamily"
rte_frame <- GRanges(RMdf %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))
mbo <- mergeByOverlaps(grs, rte_frame)
methdf <- mbo$grs %>%
    as.data.frame() %>%
    tibble()
rte_only_frame <- mbo$rte_frame %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
rtedf <- bind_cols(methdf, rte_only_frame)
write_delim(rtedf, "ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
# rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
if ("width" %in% colnames(RMdf)) {
    joinframe <<- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
} else {
    joinframe <<- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand)
}
perelementdf <- rtedf %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(joinframe)

perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))


write_delim(perelementdf, "ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)
# perelementdf <- read_delim("ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)


# RTE PROMOTERS
# annotate whether full length elements promoters overlap DMRs
flelement <- rmann %>% filter(rte_length_req == "FL")
rmann %$% rte_length_req %>% table()

flSINE <- flelement %>% filter(rte_superfamily == "SINE")
flLINE <- flelement %>% filter(rte_superfamily == "LINE")
rmann %$% ltr_viral_status %>% unique()
flFl_Provirus_5LTR <- flelement %>%
    filter(str_detect(gene_id, "LTR")) %>%
    filter(ltr_viral_status == "5'LTR (FL Int)")

flSINEgrs <- GRanges(flSINE)
flLINE5UTRgrs <- GRanges(flLINE) %>% resize(909)
flFl_Provirus_5LTRgrs <- GRanges(flFl_Provirus_5LTR)
flRTEpromotergrs <- c(c(flSINEgrs, flLINE5UTRgrs), flFl_Provirus_5LTRgrs)

if (conf$single_condition == "no") {
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
        mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, "Hyper", "Hypo"))))
} else {
    flRTEpromoterdfGOOD <- tibble(as.data.frame(flRTEpromotergrs))
}


flRTEpromoter <- flRTEpromoterdfGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
write_delim(flRTEpromoter, "ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
# flRTEpromoter <- read_delim("ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)

grouping_var <- "rte_subfamily"
rte_frame <- GRanges(flRTEpromoter %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))
mbo <- mergeByOverlaps(grs, rte_frame)
methdf <- mbo$grs %>%
    as.data.frame() %>%
    tibble()
rte_only_frame <- mbo$rte_frame %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
rtedf_promoters <- bind_cols(methdf, rte_only_frame)

write_delim(rtedf_promoters, "ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
# rtedf_promoters <- read_delim("ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(flRTEpromoter %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width))


write_delim(perelementdf_promoters, "ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
# perelementdf_promoters <- read_delim("ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)


# Genes
refseq_gr <- import(conf$refseq_unaltered)
genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
mcols(genes_gr) %>% colnames()
mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
mbo <- mergeByOverlaps(grs, promoters)
promoter_methdf <- mbo$grs %>%
    as.data.frame() %>%
    tibble()
promoter_methdf$gene_id <- mbo$promoters %>%
    as.data.frame() %>%
    tibble() %$% gene_id
write_delim(promoter_methdf, "ldna/Rintermediates/refseq_gene_promoter_methylation.tsv", col_names = TRUE)


############## Read density
## NOTE GOT UP TO HERE IN MAKING THIS WORK FOR SINGLE SAMPLE
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
        eoi <- import(paste0("aref/extended/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
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
write_delim(reads, "ldna/Rintermediates/reads_context_all.tsv", col_names = TRUE)
# reads <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)

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
        eoi <- import(paste0("aref/extended/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
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
write_delim(readscg, "ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
# readscg <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
