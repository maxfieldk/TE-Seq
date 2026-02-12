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
            bedmethylpaths = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_bedMethyl.bed", samples, samples),
            data = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/m/tables/dmrs.tsv",
            dmls = "ldna/results/m/tables/dmls.tsv",
            read_mods = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis),
            read_mods_cg_islands = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", "cpgI")
        ), env = globalenv())
        assign("params", list(mod_code = "m"), env = globalenv())
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
grs <- grsunfiltered[grsunfiltered$cov > MINIMUMCOVERAGE]
grsdf <- tibble(as.data.frame(grs))
grsdf %$% seqnames %>% unique()
dir.create("ldna/Rintermediates", recursive = TRUE)
write_delim(grsdf %>% filter(grepl("^NI_", seqnames)), sprintf("ldna/Rintermediates/%s/grsdf_nonref.tsv", params$mod_code), col_names = TRUE)
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
rm(grsdfuntidy)

# TODO pivot wider introduces sample names as variables, this needs to be addressed
grsdffiltered <- grsdf %>%
    filter(pos %in% grsinboth)
grsdf <- grsdffiltered
rm(grsdffiltered)

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
rm(grs)
grs <- c(grs_cpg_islands, grs_cpgi_shelves, grs_cpgi_shores, grs_cpg_opensea)
grsdf <- tibble(as.data.frame(grs))

dir.create("ldna/Rintermediates", showWarnings = FALSE)
write_delim(grsdf, sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
# grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
# grs <- GRanges(grsdf)

# SETTING UP SOME SUBSETS FOR EXPLORATION
set.seed(75)
grsdfs <- grsdf %>%
    group_by(sample, seqnames, islandStatus) %>%
    slice_sample(n = 1000)
grss <- GRanges(grsdfs)
write_delim(grsdfs, sprintf("ldna/Rintermediates/%s/grsdfsmall.tsv", params$mod_code), col_names = TRUE)
rm(grsdf)

if (conf$single_condition == "no") {
    dmrs <- read_delim(inputs$dmrs, delim = "\t", col_names = TRUE)
    dmls <- read_delim(inputs$dmls, delim = "\t", col_names = TRUE) %>% filter(fdrs <= 0.05)
    dmrsgr <- GRanges(dmrs)
    dmlsgr <- GRanges(
        seqnames = dmls$chr,
        ranges = IRanges(start = dmls$pos, end = dmls$pos),
        stat = dmls$stat,
        pval = dmls$pvals,
        fdr = dmls$fdrs,
        direction = dmls$direction
    )
}



####################
## RTEs
dir.create(sprintf("ldna/results/%s/plots/rte", params$mod_code), showWarnings = FALSE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
r_repeatmasker_annotation %>% filter(intactness_req == "Intact")
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
rm(r_annotation_fragmentsjoined)
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
    dmrtypes <- dmrs$dmr_type %>% unique()
    threshold_dfs <- list()
    for (dmrtype in dmrtypes) {
        dmrsgr_temp <- dmrs %>%
            filter(dmr_type == dmrtype) %>%
            GRanges()
        conditions <- conf$levels
        condition1 <- conditions[1]
        condition2 <- conditions[2]
        mbo <- mergeByOverlaps(RM, dmrsgr_temp)
        mergeddf <- tibble(as.data.frame(mbo))
        mm <- mergeddf %>%
            group_by(gene_id) %>%
            mutate(max_val = max(areaStat), min_val = min(areaStat))
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

        sboinvert <- subsetByOverlaps(RM, dmrsgr_temp, invert = TRUE)
        sboinvert$max_val <- NaN
        sboinvert$min_val <- NaN
        sboinvert$concordance <- NA_character_

        RMfinal <- c(merged, sboinvert)
        rmannfinal <- tibble(as.data.frame(RMfinal))
        rmannGOOD <- rmannfinal %>%
            mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))))
        threshold_df <- rmannGOOD %>%
            dplyr::select(gene_id, direction) %>%
            dplyr::rename(!!sym(dmrtype) := direction)
        threshold_dfs[[dmrtype]] <- threshold_df
    }
    RMdf <- full_join(rmann, purrr::reduce(threshold_dfs, full_join)) %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
} else {
    RMdf <- rmann %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
}
write_delim(RMdf, sprintf("ldna/Rintermediates/%s/RMdf.tsv", params$mod_code), col_names = TRUE)
# RMdf <- read_delim(sprintf("ldna/Rintermediates/%s/RMdf.tsv", params$mod_code), col_names = TRUE)


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

rte_frame <- GRanges(RMdf)
alltedf <- merge_with_grs(grs, rte_frame)
write_delim(alltedf, sprintf("ldna/Rintermediates/%s/alltedf.tsv", params$mod_code), col_names = TRUE)

rte_frame <- GRanges(RMdf %>% filter(rte_subfamily != "Other"))
rtedf <- merge_with_grs(grs, rte_frame)
write_delim(rtedf, sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

# rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)
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


write_delim(perelementdf, sprintf("ldna/Rintermediates/%s/perelementdf.tsv", params$mod_code), col_names = TRUE)
# perelementdf <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf.tsv", params$mod_code), col_names = TRUE)


# RTE PROMOTERS
# annotate whether full length elements promoters overlap DMRs
flelement <- rmann %>% filter(rte_length_req == "FL")
flSINE <- flelement %>% filter(rte_superfamily == "SINE")
flLINE <- flelement %>% filter(rte_superfamily == "LINE")
rmann %$% ltr_viral_status %>% unique()
flFl_Provirus_5LTR <- flelement %>%
    filter(str_detect(gene_id, "LTR")) %>%
    filter(ltr_viral_status == "5'LTR (FL Int)" | ltr_viral_status == "5'LTR (Trnc Int)")
flSINEgrs <- GRanges(flSINE)
flFl_Provirus_5LTRgrs <- GRanges(flFl_Provirus_5LTR)

# for L1HS I get several region of the 5UTR based off of consensus alignment
l1_aln_dir <- "ldna/results/m/plots/l1_alignment_meth"
fll1hs_consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", l1_aln_dir, "L1HS"))

filter_by_consensus_pos <- function(fl_grs, pos_mapping, pos_vec) {
    pos_start <- pos_vec[1]
    pos_end <- pos_vec[2]
    pos_genes <- pos_mapping %$% gene_id %>% unique()
    grs_genes <- mcols(fl_grs)$gene_id %>% unique()
    genes_to_map <- intersect(pos_genes, grs_genes)

    filter_pos_end_list <- list()
    i <- 0
    for (element in genes_to_map) {
        i <- i + 1
        print(i)
        print(element)
        dfs <- pos_mapping %>% filter(gene_id == element)
        seqval <- dfs %>%
            filter(consensus_pos == pos_end) %$% sequence_pos %>%
            pluck(1)
        if (!is.na(seqval)) {
            filter_pos <- seqval
        } else {
            start_pos <- pos_end - 1
            match <- FALSE
            while (match == FALSE) {
                seqval <- dfs %>% filter(consensus_pos == start_pos) %$% sequence_pos
                if (length(seqval) != 0) {
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
        filter_pos_end_list[[element]] <- filter_pos
    }
    print("end loop done")

    filter_pos_start_list <- list()
    i <- 0
    for (element in genes_to_map) {
        i <- i + 1
        print(i)
        print(element)
        dfs <- pos_mapping %>% filter(gene_id == element)
        if (pos_start == 0) {
            filter_pos <- 1
        } else {
            seqval <- dfs %>%
                filter(consensus_pos == pos_start) %$% sequence_pos %>%
                pluck(1)
            if (!is.na(seqval)) {
                filter_pos <- seqval
            } else {
                start_pos <- pos_start
                match <- FALSE
                while (match == FALSE) {
                    seqval <- dfs %>% filter(consensus_pos == start_pos) %$% sequence_pos
                    if (length(seqval) != 0) {
                        if (!is.na(seqval)) {
                            filter_pos <- seqval
                            match <- TRUE
                        } else {
                            start_pos <- start_pos + 1
                        }
                    } else {
                        filter_pos <- "NoRegionHomology"
                        match <- TRUE
                    }
                }
            }
        }
        filter_pos_start_list[[element]] <- filter_pos
    }
    print("start loop done")
    elements_keep <- names(filter_pos_start_list) %>% intersect(names(filter_pos_end_list))
    filter_pos_start <- filter_pos_start_list[elements_keep]
    filter_pos_end <- filter_pos_end_list[elements_keep]
    mapping <- tibble(gene_id = elements_keep, filter_pos_start = unlist(filter_pos_start), filter_pos_end = unlist(filter_pos_end))
    gene_ids_with_homology <- mapping %>%
        filter(filter_pos_start != "NoRegionHomology") %>%
        filter(filter_pos_end != "NoRegionHomology") %$% gene_id
    fl_grs_with_homology <- fl_grs[mcols(fl_grs)$gene_id %in% gene_ids_with_homology]
    grlist <- map(seq_along(fl_grs_with_homology), function(x) fl_grs_with_homology[x])
    grlistresized <- map(grlist, function(x) {
        if (mcols(x)$gene_id %in% mapping$gene_id) {
            tempdf <- mapping %>% filter(gene_id == mcols(x)$gene_id)
            final_width <- tempdf$filter_pos_end - tempdf$filter_pos_start + 1
            fix_end <- resize(x, width = (tempdf %$% filter_pos_end))
            fix_start <- resize(fix_end, width = final_width, fix = "end")
            return(fix_start)
        }
    })
    number_omitted <- length(grlist) - nrow(mapping)
    l1hs_resized <- purrr::reduce(grlistresized, c)
    return(l1hs_resized)
}

flL1HS5UTR <- filter_by_consensus_pos(fl_grs = flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), pos_mapping = fll1hs_consensus_index_long, pos_vec = c(0, 909))
flL1HS500 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 500))
flL1HS328 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 328))
flL1HSASP <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(400, 600))

flLINENOTL1HS <- flLINE %>%
    filter(rte_subfamily != "L1HS") %>%
    GRanges() %>%
    resize(909)


flRTEpromotergrs <- c(c(c(flSINEgrs, flLINENOTL1HS), flFl_Provirus_5LTRgrs), flL1HS5UTR)

mcols(flL1HS5UTR)$region <- "909"
mcols(flL1HS500)$region <- "500"
mcols(flL1HS328)$region <- "328"
mcols(flL1HSASP)$region <- "ASP"
l1hs_intra_utr_grs <- c(c(c(flL1HS328, flL1HS500), flL1HS5UTR), flL1HSASP)


if (conf$single_condition == "no") {
    dmrtypes <- dmrs$dmr_type %>% unique()
    threshold_dfs <- list()
    for (dmrtype in dmrtypes) {
        dmrsgr_temp <- dmrs %>%
            filter(dmr_type == dmrtype) %>%
            GRanges()
        mbo <- mergeByOverlaps(flRTEpromotergrs, dmrsgr_temp)
        mergeddf <- tibble(as.data.frame(mbo))
        mm <- mergeddf %>%
            group_by(gene_id) %>%
            mutate(max_val = max(areaStat), min_val = min(areaStat))
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

        sboinvert <- subsetByOverlaps(flRTEpromotergrs, dmrsgr_temp, invert = TRUE)
        sboinvert$max_val <- NaN
        sboinvert$min_val <- NaN
        sboinvert$concordance <- NA_character_

        flRTEpromotergrsfinal <- c(merged, sboinvert)
        flRTEpromoterdfinal <- tibble(as.data.frame(flRTEpromotergrsfinal))
        flRTEpromoterdfGOOD <- flRTEpromoterdfinal %>%
            mutate(direction = ifelse(is.na(concordance), NA_character_, ifelse(concordance == "discordant", "discordant", ifelse(max_val > 0, "Hyper", "Hypo"))))

        threshold_df <- flRTEpromoterdfGOOD %>%
            dplyr::select(gene_id, direction) %>%
            dplyr::rename(!!sym(dmrtype) := direction)
        threshold_dfs[[dmrtype]] <- threshold_df
    }
    flRTEpromoterdfGOOD <- full_join(tibble(as.data.frame(flRTEpromotergrs)), purrr::reduce(threshold_dfs, full_join)) %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
} else {
    flRTEpromoterdfGOOD <- tibble(as.data.frame(flRTEpromotergrs)) %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
}

flRTEpromoter <- flRTEpromoterdfGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
write_delim(flRTEpromoter, sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)
# flRTEpromoter <- read_delim(sprintf("ldna/Rintermediates/%s/flRTEpromoter.tsv", params$mod_code), col_names = TRUE)


grouping_var <- "rte_subfamily"
rte_frame <- GRanges(flRTEpromoter %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))

rtedf_promoters <- merge_with_grs(grs, rte_frame)
write_delim(rtedf_promoters, sprintf("ldna/Rintermediates/%s/rtedf_promoters.tsv", params$mod_code), col_names = TRUE)
# rtedf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf_promoters.tsv", params$mod_code), col_names = TRUE)

l1hs_intrautr <- merge_with_grs(grs, l1hs_intra_utr_grs)
write_delim(l1hs_intrautr, sprintf("ldna/Rintermediates/%s/l1hs_intrautr.tsv", params$mod_code), col_names = TRUE)

# alll1hsflids <- flRTEpromoter %>%
#     filter(rte_subfamily == "L1HS") %$% gene_id %>%
#     unique()
# coveredl1hsflids <- l1hs_intrautr %$% gene_id %>% unique()
# flL1HS_not_covered <- setdiff(alll1hsflids, coveredl1hsflids)


perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(flRTEpromoter %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand))
write_delim(perelementdf_promoters, sprintf("ldna/Rintermediates/%s/perelementdf_promoters.tsv", params$mod_code), col_names = TRUE)
rm(perelementdf_promoters)
# perelementdf_promoters <- read_delim(sprintf("ldna/Rintermediates/%s/perelementdf_promoters.tsv", params$mod_code), col_names = TRUE)

perl1hs_5utr_region <- l1hs_intrautr %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, region) %>%
    summarise(mean_meth = mean(pctM)) %>%
    ungroup()
write_delim(perl1hs_5utr_region, sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE)
# perl1hs_5utr_region <- read_delim(sprintf("ldna/Rintermediates/%s/perl1hs_5utr_region.tsv", params$mod_code), col_names = TRUE)



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
write_delim(promoter_methdf, sprintf("ldna/Rintermediates/%s/refseq_gene_promoter_methylation.tsv", params$mod_code), col_names = TRUE)


############## Read density
## NOTE GOT UP TO HERE IN MAKING THIS WORK FOR SINGLE SAMPLE
rm(grsdf)
rm(grs)
mem_used()

dir.create("ldna/results/plots/reads")
library(Biostrings)

# readslist <- list()
# for (region in conf$rte_subfamily_read_level_analysis) {
#     for (sample_name in samples) {
#         df <- read_delim(
#             grep(region,
#                 grep(sprintf("/%s/", sample_name),
#                     inputs$read_mods,
#                     value = TRUE
#                 ),
#                 value = TRUE
#             )
#         )
#         df$region <- region
#         df$sample <- sample_name
#         df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
#         grsx <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
#         eoi <- import(paste0("aref/extended/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
#         mbo <- mergeByOverlaps(grsx, eoi)
#         df1 <- as.data.frame(mbo) %>%
#             tibble() %>%
#             dplyr::select(starts_with("grsx"), name) %>%
#             dplyr::rename(gene_id = name)
#         colnames(df1) <- gsub("grsx.", "", colnames(df1))
#         readslist <- c(readslist, list(df1))
#     }
# }
# reads <- Reduce(rbind, readslist) %>% filter(mod_code == params$mod_code)
# write_delim(reads, sprintf("ldna/Rintermediates/%s/reads_context_all.tsv", params$mod_code), col_names = TRUE)
# reads <- read_delim(sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE)

readslistcg <- list()
for (region in conf$rte_subfamily_read_level_analysis) {
    for (sample_name in samples) {
        df <- read_delim(
            grep(
                "CpG",
                grep(region,
                    grep(sprintf("/%s/", sample_name),
                        inputs$read_mods_cg,
                        value = TRUE
                    ),
                    value = TRUE
                ),
                value = TRUE
            ),
        )
        df$region <- region
        df$sample <- sample_name
        df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grsx <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0("aref/extended/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
        strand(eoi) <- "*"
        mbo <- mergeByOverlaps(grsx, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grsx"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grsx.", "", colnames(df1))
        readslistcg <- c(readslistcg, list(df1))
    }
}

readscg <- Reduce(rbind, readslistcg)
write_delim(readscg, sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE)
# readscg <- read_delim(sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE)
