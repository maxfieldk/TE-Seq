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
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)

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
            read_mods = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis),
            read_mods_cg_islands = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", "cpgI")
        ), env = globalenv())
        assign("params", list(mod_code = "m"), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

contrasts <- conf$contrasts

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis


##########################
# PREP DATA FOR ANALYSIS
sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE), col_names = FALSE)
    df_m_sub <- df %>% filter(X4 == "m") %>% select(X1, X2, cov = X10, pctM = X11, Nmod_m = X12)
    df_h_sub <- df %>% filter(X4 == "h") %>% select(X1, X2, pctH = X11, Nmod_h = X12)
    rm(df)
    df_joined <- inner_join(df_m_sub, df_h_sub, by = c("X1", "X2"))
    rm(df_m_sub, df_h_sub)
    df_long <- bind_rows(
        transmute(df_joined, X1, X2, cov, mod = "m", pctM = as.double(pctM)),
        transmute(df_joined, X1, X2, cov, mod = "h", pctM = as.double(pctH)),
        transmute(df_joined, X1, X2, cov, mod = "u", pctM = pmax(0, (cov - Nmod_m - Nmod_h) / cov * 100))
    )
    rm(df_joined)
    gr <- GRanges(
        seqnames = df_long$X1,
        ranges = IRanges(start = df_long$X2, end = df_long$X2),
        cov = df_long$cov,
        mod = df_long$mod,
        pctM = df_long$pctM
    )
    gr$sample <- sample_name
    gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_grs[[sample_name]] <- gr
}

grsunfiltered <- Reduce(c, sample_grs)
rm(sample_grs)
# filter out low coverage and ensure that all samples have the same cpgs
grs <- grsunfiltered[grsunfiltered$cov > MINIMUMCOVERAGE]
rm(grsunfiltered)
grsdf <- tibble(as.data.frame(grs))
grsdf %$% seqnames %>% unique()

dir.create("ldna/Rintermediates", recursive = TRUE)
write_delim(grsdf %>% filter(grepl("^NI_", seqnames)), "ldna/Rintermediates/grsdf_nonref.tsv", col_names = TRUE)
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
seqnames <- grsdf$seqnames
start <- grsdf$start
end <- grsdf$end
pos <- paste0(seqnames, "_", start, "_", end)
grsdf$pos <- pos

grsdfuntidy <- grsdf %>%
    filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
    filter(mod == "m") %>%
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
rm(grs_cpg_islands, grs_cpgi_shelves, grs_cpgi_shores, grs_cpg_opensea, cpg_islands, cpgi_shores, cpgi_shelves, cpgi_features)
grsdf <- tibble(as.data.frame(grs))

write_delim(grsdf, "ldna/Rintermediates/grsdf.tsv", col_names = TRUE)
# grsdf <- read_delim("ldna/Rintermediates/grsdf.tsv", col_names = TRUE)
# grs <- GRanges(grsdf)

# SETTING UP SOME SUBSETS FOR EXPLORATION
set.seed(75)
grsdfs <- grsdf %>%
    group_by(sample, mod, seqnames, islandStatus) %>%
    slice_sample(n = 1000)
grss <- GRanges(grsdfs)
write_delim(grsdfs, "ldna/Rintermediates/grsdfsmall.tsv", col_names = TRUE)
rm(grsdf, grsdfs, grss)



####################
## RTEs
dir.create("ldna/results/plots/rte", showWarnings = FALSE, recursive = TRUE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
r_repeatmasker_annotation %>% filter(intactness_req == "Intact")
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
rm(r_annotation_fragmentsjoined)
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
rm(r_repeatmasker_annotation)

# FULL ELEMENTS
# annotate whether repeats overlap DMRs

RMdf <- rmann %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)

write_delim(RMdf, "ldna/Rintermediates/RMdf.tsv", col_names = TRUE)
# RMdf <- read_delim("ldna/Rintermediates/RMdf.tsv", col_names = TRUE)


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

rte_frame <- GRanges(RMdf %>% filter(rte_subfamily != "Other"))
rtedf <- merge_with_grs(grs, rte_frame)
rm(rte_frame)
write_delim(rtedf, "ldna/Rintermediates/rtedf.tsv", col_names = TRUE)

# rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
if ("width" %in% colnames(RMdf)) {
    joinframe <<- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
} else {
    joinframe <<- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand)
}
perelementdf <- rtedf %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, mod) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(joinframe)

perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))

write_delim(perelementdf, "ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)
rm(rtedf, perelementdf, RMdf, joinframe)
# perelementdf <- read_delim("ldna/Rintermediates/perelementdf.tsv", col_names = TRUE)


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
fll1hs_consensus_index_long %$% gene_id %>% unique()
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
                    seqval <- dfs %>%
                        filter(consensus_pos == start_pos) %$% sequence_pos %>%
                        pluck(1)
                    if (!is.null(seqval)) {
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
    if (length(elements_keep) == 0) {
        print("No overlapping gene_ids found between GRanges and consensus mapping")
        return(GRanges())
    }
    filter_pos_start <- filter_pos_start_list[elements_keep]
    filter_pos_end <- filter_pos_end_list[elements_keep]
    mapping <- tibble(gene_id = elements_keep, filter_pos_start = unlist(filter_pos_start), filter_pos_end = unlist(filter_pos_end)) %>%
        filter(filter_pos_start != "NoRegionHomology") %>%
        filter(filter_pos_end != "NoRegionHomology") %>%
        mutate(filter_pos_start = as.numeric(filter_pos_start), filter_pos_end = as.numeric(filter_pos_end))
    gene_ids_with_homology <- mapping %$% gene_id
    if (length(gene_ids_with_homology) == 0) {
        print("No gene_ids with homology found")
        return(GRanges())
    }
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
    grlistresized <- purrr::compact(grlistresized)
    number_omitted <- length(grlist) - length(grlistresized)
    if (length(grlistresized) == 0) {
        print("No resized GRanges produced")
        return(GRanges())
    }
    l1hs_resized <- purrr::reduce(grlistresized, c)
    return(l1hs_resized)
}
flL1HS5UTR <- filter_by_consensus_pos(fl_grs = flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), pos_mapping = fll1hs_consensus_index_long, pos_vec = c(0, 909))
flL1HS500 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 500))
flL1HS328 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 328))
flL1HS110 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 110))
flL1HSASP <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(400, 600))

flLINENOTL1HS <- flLINE %>%
    filter(rte_subfamily != "L1HS") %>%
    GRanges() %>%
    resize(909)


flRTEpromotergrs <- c(c(c(flSINEgrs, flLINENOTL1HS), flFl_Provirus_5LTRgrs), flL1HS5UTR)

mcols(flL1HS5UTR)$region <- "909"
mcols(flL1HS500)$region <- "500"
mcols(flL1HS328)$region <- "328"
mcols(flL1HS110)$region <- "110"
mcols(flL1HSASP)$region <- "ASP"
l1hs_intra_utr_grs <- c(c(c(c(flL1HS110, flL1HS328), flL1HS500), flL1HS5UTR), flL1HSASP)


flRTEpromoterdfGOOD <- tibble(as.data.frame(flRTEpromotergrs)) %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)


flRTEpromoter <- flRTEpromoterdfGOOD %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
write_delim(flRTEpromoter, "ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
# flRTEpromoter <- read_delim("ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)


grouping_var <- "rte_subfamily"
rte_frame <- GRanges(flRTEpromoter %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))

rtedf_promoters <- merge_with_grs(grs, rte_frame)
write_delim(rtedf_promoters, "ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
# rtedf_promoters <- read_delim("ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)

l1hs_intrautr <- merge_with_grs(grs, l1hs_intra_utr_grs)
write_delim(l1hs_intrautr, "ldna/Rintermediates/l1hs_intrautr.tsv", col_names = TRUE)

# # grsdf
# pf <- l1hs_intrautr %>%
#     filter(region == "500") %>%
#     group_by(sample, condition) %>%
#     summarise(mm = mean(pctM))
# pfmean <- pf %>% group_by(condition) %>% summarise(mm = mean(mm))
# p <- ggplot(pf) +
#     geom_point(aes(x = condition, y = mm), size = 1.5, shape = 16, alpha = 0.75) +
#     mtopen
# mysaveandstore(fn = "zz11111.pdf")
# alll1hsflids <- flRTEpromoter %>%
#     filter(rte_subfamily == "L1HS") %$% gene_id %>%
#     unique()
# coveredl1hsflids <- l1hs_intrautr %$% gene_id %>% unique()
# flL1HS_not_covered <- setdiff(alll1hsflids, coveredl1hsflids)


perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, mod) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(flRTEpromoter %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand))
write_delim(perelementdf_promoters, "ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
rm(perelementdf_promoters)
# perelementdf_promoters <- read_delim("ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)

perl1hs_5utr_region <- l1hs_intrautr %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition, mod, region) %>%
    summarise(mean_meth = mean(pctM)) %>%
    ungroup()
write_delim(perl1hs_5utr_region, "ldna/Rintermediates/perl1hs_5utr_region.tsv", col_names = TRUE)
# perl1hs_5utr_region <- read_delim("ldna/Rintermediates/perl1hs_5utr_region.tsv", col_names = TRUE)



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
        # Join m and h to one row per CpG per read, then compute u
        df_m <- df %>% filter(mod_code == "m") %>% rename(mod_qual_m = mod_qual) %>% select(-mod_code)
        df_h <- df %>% filter(mod_code == "h") %>% select(read_id, ref_position, mod_qual_h = mod_qual)
        rm(df)
        df_joined <- inner_join(df_m, df_h, by = c("read_id", "ref_position"))
        rm(df_m, df_h)

        df_joined$region <- region
        df_joined$sample <- sample_name
        df_joined$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grsx <- GRanges(df_joined %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0("aref/extended/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
        strand(eoi) <- "*"
        mbo <- mergeByOverlaps(grsx, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grsx"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grsx.", "", colnames(df1))

        # Pivot to long with m, h, u
        df1_long <- bind_rows(
            df1 %>% mutate(mod_code = "m", mod_qual = mod_qual_m) %>% select(-mod_qual_m, -mod_qual_h),
            df1 %>% mutate(mod_code = "h", mod_qual = mod_qual_h) %>% select(-mod_qual_m, -mod_qual_h),
            df1 %>% mutate(mod_code = "u", mod_qual = pmax(0, 1 - mod_qual_m - mod_qual_h)) %>% select(-mod_qual_m, -mod_qual_h)
        )
        readslistcg <- c(readslistcg, list(df1_long))
    }
}

readscg <- Reduce(rbind, readslistcg)
write_delim(readscg, "ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
# readscg <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)






###################### get bed files for L1 promoters

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
                    seqval <- dfs %>%
                        filter(consensus_pos == start_pos) %$% sequence_pos %>%
                        pluck(1)
                    if (!is.null(seqval)) {
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
    if (length(elements_keep) == 0) {
        print("No overlapping gene_ids found between GRanges and consensus mapping")
        return(GRanges())
    }
    filter_pos_start <- filter_pos_start_list[elements_keep]
    filter_pos_end <- filter_pos_end_list[elements_keep]
    mapping <- tibble(gene_id = elements_keep, filter_pos_start = unlist(filter_pos_start), filter_pos_end = unlist(filter_pos_end)) %>%
        filter(filter_pos_start != "NoRegionHomology") %>%
        filter(filter_pos_end != "NoRegionHomology") %>%
        mutate(filter_pos_start = as.numeric(filter_pos_start), filter_pos_end = as.numeric(filter_pos_end))
    gene_ids_with_homology <- mapping %$% gene_id
    if (length(gene_ids_with_homology) == 0) {
        print("No gene_ids with homology found")
        return(GRanges())
    }
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
    grlistresized <- purrr::compact(grlistresized)
    number_omitted <- length(grlist) - length(grlistresized)
    if (length(grlistresized) == 0) {
        print("No resized GRanges produced")
        return(GRanges())
    }
    l1hs_resized <- purrr::reduce(grlistresized, c)
    return(l1hs_resized)
}

flL1HS5UTRnonref <- filter_by_consensus_pos(fl_grs = flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(grepl("^NI", seqnames)) %>% GRanges(), pos_mapping = fll1hs_consensus_index_long, pos_vec = c(0, 909))

flL1HS5UTR <- filter_by_consensus_pos(fl_grs = flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), pos_mapping = fll1hs_consensus_index_long, pos_vec = c(0, 909))
flL1HS500 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 500))
flL1HS328 <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(0, 328))
flL1HSASP <- filter_by_consensus_pos(flLINE %>% filter(rte_subfamily == "L1HS") %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>% GRanges(), fll1hs_consensus_index_long, c(400, 600))

flL1HS5UTR_bed4 <- flL1HS5UTR %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::select(seqnames, start, end, gene_id)

flL1HS500_bed4 <- flL1HS500 %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::select(seqnames, start, end, gene_id)

flL1HS328_bed4 <- flL1HS328 %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::select(seqnames, start, end, gene_id)
write_delim(flL1HS5UTR_bed4, "ldna/Rintermediates/FL_L1HS_Ref_909.bed", delim = "\t", col_names = FALSE)
write_delim(flL1HS500_bed4, "ldna/Rintermediates/FL_L1HS_Ref_500.bed", delim = "\t", col_names = FALSE)
write_delim(flL1HS328_bed4, "ldna/Rintermediates/FL_L1HS_Ref_328.bed", delim = "\t", col_names = FALSE)
