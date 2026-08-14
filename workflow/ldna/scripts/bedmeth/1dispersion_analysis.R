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





##################

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


###########################



## READ DISPERSION ANALYSIS
{
    flyngl1 <- rmannextended %>%
        filter(rte_length_req == "FL") %>%
        filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2")

    flyngl1grs <- GRanges(flyngl1)
    flyngl1grspromoter <- promoters(flyngl1grs, upstream = 0, downstream = 909)


    l1cpggrs <- mergeByOverlaps(cpg_islands, flyngl1grspromoter)
    l1cpggrs1 <- l1cpggrs$cpg_islands
    mcols(l1cpggrs1)$gene_id <- mcols(l1cpggrs$flyngl1grspromoter)$gene_id

    promoterscpggrs <- mergeByOverlaps(cpg_islands, promoters)
    promoterscpggrs1 <- promoterscpggrs$cpg_islands
    mcols(promoterscpggrs1)$gene_id <- mcols(promoterscpggrs$promoters)$gene

    flanked_promoters <- flank(promoterscpggrs1, width = 10000, both = TRUE)

    merged <- mergeByOverlaps(flanked_promoters, l1cpggrs1)
    merged[!(mcols(merged$flanked_promoters)$name == mcols(merged$l1cpggrs1)$name), ]


    rmannextended %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        pw()


    reads_cgI_chr1 <- read_delim("ldna/intermediates/merged/merged_cpgI_reads1.tsv")
    reads_cgI_chr1_grs <- GRanges(reads_cgI_chr1 %>% dplyr::rename(seqnames = chrom, start = ref_position) %>% mutate(end = start))

    mbo <- mergeByOverlaps(reads_cgI_chr1_grs, promoters)
    mboreadsdf <- mbo$reads_cgI_chr1_grs %>%
        as.data.frame() %>%
        tibble()
    mbopromoters <- mbo$promoters %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width)

    readscg_genepromoters <- bind_cols(mboreadsdf, mbopromoters)
    rm(mbopromoters)
    rm(readsdf)
    rm(mbo)

    chr1rtes_grs <- rmannextended %>%
        filter(seqnames == "chr1") %>%
        filter(rte_subfamily != "Other") %>%
        GRanges()
    mbo <- mergeByOverlaps(reads_cgI_chr1_grs, chr1rtes_grs)
    mboreadsdf <- mbo$reads_cgI_chr1_grs %>%
        as.data.frame() %>%
        tibble()
    mbortes <- mbo$chr1rtes_grs %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width)

    readscg_rtes <- bind_cols(mboreadsdf, mbortes)
    rm(mboreadsdf)
    rm(mbortes)
    rm(mbo)


    by_cpg_genepromoters <- readscg_genepromoters %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    numCGneeded_genes <- by_cpg_genepromoters %>%
        group_by(gene_id, read_id) %>%
        summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        group_by(gene_id) %>%
        summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    by_cpg_genepromoters_ds <- by_cpg_genepromoters %>%
        filter(num_cpgs_in_read >= 25) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = 25, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()

    by_read_genepromoters_ds <- by_cpg_genepromoters_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id) %>%
        mutate(mean_fm = mean(fraction_meth))

    # by_read_genepromoters_ds %>%
    #     group_by(gene_id) %>%
    #     summarise(mean_fm = mean(fraction_meth)) %$% mean_fm %>%
    #     quantile(., probs = seq(0, 1, 0.05))
    # by_read_genepromoters_ds %$% gene_id %>% unique()


    # p <- by_read_genepromoters_ds %>%
    #     group_by(gene_id) %>%
    #     summarise(mean_fm = mean(fraction_meth)) %>%
    #     ggplot() +
    #     geom_histogram(aes(x = mean_fm))
    # mysaveandstore()

    # by_read_genepromoters_ds %>%
    #     filter(mean_fm > 0.7) %>%
    #     filter(mean_fm < 0.98) %$% gene_id %>%
    #     unique()


    dispersion_model_genes_high_meth <- glmmTMB(
        cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
        dispformula = ~condition,
        family = glmmTMB::betabinomial(),
        data = by_read_genepromoters_ds %>% filter(mean_fm > 0.7) %>% filter(mean_fm < 0.98) %>%
            mutate(meth = fraction_meth * 25, unmeth = 25 - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
    )

    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, "genes", "cpg25req")
    dir.create(outputdirtables, recursive = TRUE)
    summary(dispersion_model_genes_high_meth)
    res_text <- capture.output(dispersion_model_genes_high_meth %>% summary())
    writeLines(res_text, sprintf("%s/gene_high_meth_dispersion_model_summary.txt", outputdirtables))




    by_cpg_rtes <- readscg_rtes %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    numCGneeded_rtes <- by_cpg_rtes %>%
        group_by(gene_id, read_id) %>%
        summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        group_by(gene_id) %>%
        summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    by_cpg_rtes_ds <- by_cpg_rtes %>%
        filter(num_cpgs_in_read >= 25) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = 25, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()

    by_read_rtes_ds <- by_cpg_rtes_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id) %>%
        mutate(mean_fm = mean(fraction_meth)) %>%
        left_join(rmannextended %>% dplyr::select(gene_id, rte_subfamily, rte_length_req))




    rtemodels <- list()
    for (subfam in c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")) {
        dispersion_model_high_meth_fl <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = by_read_rtes_ds %>% filter(mean_fm > 0.7) %>% filter(mean_fm < 0.98) %>%
                filter(rte_subfamily == subfam) %>%
                filter(rte_length_req == "FL") %>%
                mutate(meth = fraction_meth * 25, unmeth = 25 - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
        )
        rtemodels[[subfam]] <- dispersion_model_high_meth_fl
    }

    res_text <- capture.output(map(dispersion_model_genes_high_meth, summary))
    writeLines(res_text, sprintf("%s/gene_high_meth_dispersion_model_summary.txt", outputdirtables))



    ################### reads to extract
    imprinted_genes_df <- read_csv("/users/mkelsey/data/Nanopore/alz/imprinted_genes.csv", comment = "#")
    imprinted_genes_df %$% Aliases %>%
        str_split(., ", ") %>%
        map(., ~ trimws(.x))

    temppromoteranalysis <- ip %>%
        group_by(gene_id) %>%
        summarise(mm = mean(mean_meth)) %>%
        filter(!grepl("-AS", gene_id))
    highmethgenes <- temppromoteranalysis %>% filter(mm > 75) %$% gene_id
    lowmethgenes <- temppromoteranalysis %>% filter(mm < 15) %$% gene_id
    imprinted_genes <- temppromoteranalysis %>% filter(gene_id %in% (imprinted_genes_df %$% gene_id)) %$% gene_id

    sex_chromosome_genes <- mcols(promoters[seqnames(promoters) == "chrX" | seqnames(promoters) == "chrY"])$gene_id

    # non interesting cpgIs
    # filter out any related to genes, enhancers, ccres, L1s
    gencode <- import("/users/mkelsey/data/Nanopore/alz/gencodeV35hs1.bed", format = "BED")

    boring_islands <- cpg_islands %>%
        subsetByOverlaps(gencode, invert = TRUE) %>%
        subsetByOverlaps(promoters, invert = TRUE) %>%
        subsetByOverlaps(refseq_gr[mcols(refseq_gr)$type == "gene"], invert = TRUE) %>%
        subsetByOverlaps(rmannextended %>% GRanges(), invert = TRUE) %>%
        subsetByOverlaps(ccresgr, invert = TRUE) %>%
        subsetByOverlaps(chromHMMgr[mcols(chromHMMgr)$name == "Quies"], invert = FALSE)

    set.seed(72)
    boring_islands_to_extract_reads_from <- boring_islands[width(boring_islands) < 1000 & !(seqnames(boring_islands) %in% c("chrX", "chrY"))] %>%
        sample(., size = 500)

    mcols(boring_islands_to_extract_reads_from)$score <- NULL
    mcols(boring_islands_to_extract_reads_from)$gene_id <- paste0(mcols(boring_islands_to_extract_reads_from)$name, "_", start(boring_islands_to_extract_reads_from))
    mcols(boring_islands_to_extract_reads_from)$name <- NULL
    boringgrs <- merge_with_grs(grs, boring_islands_to_extract_reads_from)
    boringgrsHMids <- boringgrs %>%
        as.data.frame() %>%
        tibble() %>%
        group_by(gene_id) %>%
        summarise(pctM = mean(pctM)) %>%
        filter(pctM > 85) %$% gene_id
    boringgrsHM %>% print(n = 50)
    as.data.frame(boring_islands_to_extract_reads_from[mcols(boring_islands_to_extract_reads_from)$gene_id %in% boringgrsHMids]) %>%
        tibble() %>%
        write_delim("ldna/Rintermediates/m/boringcpgi_highmeth.tsv", delim = "\t")

    set.seed(73)
    genes_to_extract_reads_from <- c(
        highmethgenes[!(highmethgenes %in% sex_chromosome_genes)] %>% sample(., size = 500),
        lowmethgenes[!(lowmethgenes %in% sex_chromosome_genes)] %>% sample(., size = 500),
        imprinted_genes[!(imprinted_genes %in% sex_chromosome_genes)]
    )

    genes_to_extract_reads_from_grs <- promoters[mcols(promoters)$gene_id %in% genes_to_extract_reads_from]

    set.seed(74)
    rtes_to_extract_reads_from_grs <- rmannextended %>%
        filter(!(seqnames %in% c("chrX", "chrY"))) %>%
        filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")) %>%
        filter(rte_length_req == "FL") %>%
        filter(element_start < 150) %>%
        GRanges() %>%
        promoters(., upstream = 0, downstream = 909) %>%
        as.data.frame() %>%
        tibble() %>%
        group_by(rte_subfamily) %>%
        slice_sample(n = 500, replace = FALSE) %>%
        GRanges()

    regions_to_extract_reads_from_with_mcols <- c(boring_islands_to_extract_reads_from, genes_to_extract_reads_from_grs, rtes_to_extract_reads_from_grs)
    strand(regions_to_extract_reads_from_with_mcols) <- "*"
    mcols(regions_to_extract_reads_from_with_mcols)$name <- mcols(regions_to_extract_reads_from_with_mcols)$gene_id
    rtracklayer::export.bed(regions_to_extract_reads_from_with_mcols, "ldna/Rintermediates/m/regions_to_extract_reads_from.bed")


    ###
    ad1ogreads <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/ldna/intermediates/AD1/methylation/analysis_default/AD1_readmods_CpG_L1HS_rte_length_req_ALL.tsv")
    ad1ogreadsnocontext <- read_delim("/users/mkelsey/data/Nanopore/alz/RTE/ldna/intermediates/AD1/methylation/analysis_default/AD1_readmods_NoContext_L1HS_rte_length_req_ALL.tsv")

    readscg %>% filter(read_id == "07d6cf27-a6cd-4378-bbbf-8daf468382c5")
    ad1ogreads %>% filter(read_id == "98ec6f08-1086-4c41-9a9b-65c07b2c637c")
    ad1ogreads %>% filter(read_id == "fcfc619a-3b42-4b3a-9b78-0798cd999f8b")
    ad1ogreads %>% filter(read_id == "c098ac16-64a6-4df3-a294-7fac505968d0")
    readscg %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")

    ad1ogreads %>% filter(read_id == "38c0e779-bcab-4df3-a7b6-f1a27c6bc335")

    ad1ogreadsnocontext %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")

    reads_rois %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")
    reads_rois <- read_delim("ldna/intermediates/merged/merged_rois_reads.tsv")
    reads_rois_grs <- GRanges(reads_rois %>% dplyr::rename(seqnames = chrom, start = ref_position) %>% mutate(end = start))


    mbo <- mergeByOverlaps(reads_rois_grs, regions_to_extract_reads_from_with_mcols)
    mboreadsdf <- mbo$reads_rois_grs %>%
        as.data.frame() %>%
        tibble()
    mborois <- mbo$regions_to_extract_reads_from_with_mcols %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(gene_seqnames = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width, gene_id)

    readscg_rois <- bind_cols(mboreadsdf, mborois)
    rm(mboreadsdf)
    rm(mborois)
    rm(mbo)


    by_cpg_readscg_rois <- readscg_rois %>%
        group_by(gene_id, read_id, sample_name) %>%
        mutate(num_cpgs_in_read = n()) %>%
        relocate(gene_id) %>%
        ungroup()


    # numCGneeded_rois <- by_cpg_readscg_rois %>%
    #     group_by(gene_id, read_id) %>%
    #     summarise(num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
    #     group_by(gene_id) %>%
    #     summarise(mncg = mean(num_cpgs_in_read), q10 = quantile(num_cpgs_in_read, probs = c(0.1)))

    # will downsample to 25
    downsample_to_n <- 20
    set.seed(10)
    by_cpg_rois_ds <- by_cpg_readscg_rois %>%
        filter(num_cpgs_in_read >= downsample_to_n) %>%
        group_by(sample_name, gene_id) %>%
        mutate(n_reads_pass_filter = n()) %>%
        filter(n_reads_pass_filter >= 10) %>%
        group_by(gene_id, read_id) %>%
        slice_sample(n = downsample_to_n, replace = FALSE) %>%
        mutate(fraction_meth = mean(mod_indicator)) %>%
        ungroup()
    by_cpg_rois_ds %>% write_csv("ldna/Rintermediates/m/by_cpg_rois_ds")

    by_read_rois_ds <- by_cpg_rois_ds %>%
        group_by(read_id, gene_id, sample_name) %>%
        summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read)) %>%
        ungroup() %>%
        group_by(gene_id, sample_name) %>%
        mutate(mean_fm = mean(fraction_meth)) # can use this to remove high meth genes that are actually low in one sample

    by_read_imprinted <- by_read_rois_ds %>%
        filter(gene_id %in% imprinted_genes) %>%
        mutate(roi = "imprinted_gene") %>%
        filter(mean_fm > .80)
    by_read_highmeth_gene <- by_read_rois_ds %>%
        filter(!(gene_id %in% imprinted_genes)) %>%
        filter(gene_id %in% highmethgenes) %>%
        mutate(roi = "high_meth_gene") %>%
        filter(mean_fm > .80)
    by_read_lowmeth_gene <- by_read_rois_ds %>%
        filter(!(gene_id %in% imprinted_genes)) %>%
        filter(gene_id %in% lowmethgenes) %>%
        mutate(roi = "low_meth_gene") %>%
        filter(mean_fm < .20)
    by_read_boringcpgi <- by_read_rois_ds %>%
        filter(grepl("CpG:", gene_id)) %>%
        mutate(roi = "isolated_cpgI") %>%
        filter(mean_fm > .80)
    by_read_yngl1s <- by_read_rois_ds %>%
        filter(gene_id %in% mcols(rtes_to_extract_reads_from_grs)$gene_id) %>%
        left_join(rmannextended) %>%
        mutate(roi = rte_subfamily) %>%
        filter(mean_fm > .80)

    dfs <- list(
        "imprinted" = by_read_imprinted,
        "highmethgene" = by_read_highmeth_gene,
        "lowmethgene" = by_read_lowmeth_gene,
        "boringcpgi" = by_read_boringcpgi
    )

    dfsl1s <- split(by_read_yngl1s %>% ungroup() %>% dplyr::select(-colnames(rmannextended)[!colnames(rmannextended) %in% c("gene_id", "rte_subfamily")]), by_read_yngl1s$rte_subfamily)

    dfsall <- c(dfs, dfsl1s)
    dfsallbound <- bind_rows(dfsall)


    cpg_fraction_reads_highly_demeth <- dfsallbound %>%
        mutate(fraction_meth_lt50 = case_when(
            fraction_meth < 0.50 ~ 1,
            TRUE ~ 0,
        )) %>%
        mutate(fraction_meth_lt25 = case_when(
            fraction_meth < 0.25 ~ 1,
            TRUE ~ 0,
        )) %>%
        mutate(fraction_meth_lt10 = case_when(
            fraction_meth < 0.1 ~ 1,
            TRUE ~ 0,
        )) %>%
        group_by(roi) %>%
        summarise(
            nreads = n(),
            mlt50 = mean(fraction_meth_lt50),
            mlt25 = mean(fraction_meth_lt25),
            mtl10 = mean(fraction_meth_lt10)
        )

    p <- dfsallbound %>% ggplot() +
        geom_density(aes(x = fraction_meth)) +
        facet_wrap(~roi) +
        mtclosed
    mysaveandstore(fn = "zzzte3.pdf", w = 8, h = 6)


    p <- dfsallbound %>% ggplot() +
        geom_histogram(
            aes(x = fraction_meth, y = after_stat(count / sum(count))),
            position = "identity"
        ) +
        facet_wrap(~roi) +
        mtclosed
    mysaveandstore(fn = "zzzte4.pdf", w = 8, h = 6)


    getdispmodel <- function(df) {
        dispersion_model_genes_high_meth <- glmmTMB(
            cbind(meth, unmeth) ~ 1 + (1 | sample_name) + (1 | gene_id) + (1 | sample_name:gene_id),
            dispformula = ~condition,
            family = glmmTMB::betabinomial(),
            data = df %>%
                mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>% left_join(sample_table %>% dplyr::select(sample_name, condition))
        )
    }

    models_dfs <- map(dfs, getdispmodel)
    models_l1s <- map(dfsl1s, getdispmodel)

    modelsall <- c(models_dfs, models_l1s)


    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, "allrois_types", sprintf("cpg%sreq", downsample_to_n))
    dir.create(outputdirtables, recursive = TRUE)
    res_text <- capture.output(map(models_l1s, summary))
    writeLines(res_text, sprintf("%s/model_summary.txt", outputdirtables))

    modelsall[[1]] %>% broom::tidy()
    fixef(modelsall[[1]])$disp
    modelsall[[1]]$fit$par["thetaf"] # Overdispersion theta (log scale)
    sigma(modelsall[[1]], type = "dispersion")
    VarCorr(modelsall[[1]])
    summary(modelsall[[1]])
    summary(modelsall[[1]]) %>% as.data.frame()
    summary(modelsall[[1]])$coefficients$disp

    disp_df <- map(modelsall, ~ summary(.x)$coefficients$disp %>%
        as.data.frame() %>%
        rownames_to_column("term")) %>%
        bind_rows(., .id = "roi")
    write_csv(disp_df, sprintf("%s/model_disp_coef.csv", outputdirtables))



    ### mixture model
    library(brms)

    # Define mixture of two binomial components
    mix <- mixture(binomial(), binomial())

    fit <- brm(
        bf(meth | trials(20) ~ 1), # Intercept-only for now
        family = mix,
        data = dfsall[["L1HS"]] %>%
            mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
            left_join(sample_table %>%
                dplyr::select(sample_name, condition)) %>%
            dplyr::select(gene_id, sample_name, meth),
        chains = 4,
        iter = 4000,
        #   control = list(adapt_delta = 0.95),
        cores = 4
    )

    library(flexmix)
    fitflex <- flexmix(
        cbind(meth, 20 - meth) ~ 1,
        data = dfsall[["L1HS"]] %>%
            mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
            left_join(sample_table %>%
                dplyr::select(sample_name, condition)) %>%
            dplyr::select(gene_id, sample_name, meth) %>%
            filter(gene_id == "L1HS_4q28.3_9"),
        k = 3, # number of components
        model = FLXglm(family = "binomial")
    )

    summary(fitflex)

    parameters(fitflex)

    posterior_probs <- posterior(fitflex)

    mixturemodelk2 <- function(df) {
        fit <- flexmix(
            cbind(meth, 20 - meth) ~ 1,
            data = df %>%
                ungroup() %>%
                mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
                left_join(sample_table %>%
                    dplyr::select(sample_name, condition)) %>%
                dplyr::select(gene_id, sample_name, meth),
            k = 2, # number of components
            model = FLXglm(family = "binomial")
        )
        tibble(proportions = prior(fit), mean_meth = unname(plogis(parameters(fit))))
    }


    mix_res <- bind_rows(map(dfsall, mixturemodelk2), .id = "type")

    dattemp <- dfsall[["imprinted"]] %>%
        mutate(meth = fraction_meth * downsample_to_n, unmeth = downsample_to_n - meth) %>%
        left_join(sample_table %>%
            dplyr::select(sample_name, condition)) %>%
        dplyr::select(gene_id, sample_name, meth) %>%
        filter(gene_id == "L1HS_4q28.3_9")


    mixturemodelk2(dfsall[["imprinted"]])
}






################ Overdispersion plots

library(ggplot2)
library(dplyr)

# Parameters

pf <- tibble(x = numeric(), prob = numeric(), condition = character(), ncpg = numeric())
for (numcpg in c(25, 34, 39)) {
    p <- 0.75
    phi_ctrl <- 13.61
    phi_ad <- 9.68

    # Beta-binomial parameters
    alpha_ctrl <- p * phi_ctrl
    beta_ctrl <- (1 - p) * phi_ctrl

    alpha_ad <- p * phi_ad
    beta_ad <- (1 - p) * phi_ad

    # Function for beta-binomial PMF
    dbetabinom <- function(x, n, alpha, beta) {
        choose(n, x) * beta(x + alpha, n - x + beta) / beta(alpha, beta)
    }
    x_vals <- 0:numcpg

    df <- bind_rows(
        data.frame(
            x = x_vals,
            prob = dbinom(x_vals, numcpg, p),
            condition = "Bin",
            ncpg = numcpg
        ),
        data.frame(
            x = x_vals,
            prob = dbetabinom(x_vals, numcpg, alpha_ctrl, beta_ctrl),
            condition = "BbC",
            ncpg = numcpg
        ),
        data.frame(
            x = x_vals,
            prob = dbetabinom(x_vals, numcpg, alpha_ad, beta_ad),
            condition = "BbAD",
            ncpg = numcpg
        )
    )
    pf <- pf %>% bind_rows(df)
}
pf <- pf %>% mutate(perc_meth = x / ncpg)
p <- ggplot(pf, aes(x = perc_meth, y = prob, color = condition)) +
    geom_line(size = 1.2) +
    geom_point() +
    labs(
        x = "Read Percent Methylated",
        y = "Probability",
        title = "Binomial vs. Beta-binomial (different φ)"
    ) +
    facet_wrap(~ncpg) +
    mtclosed +
    scale_palette
mysaveandstore(str_glue("overdispersion119.pdf"), w = 8.9, h = 4)

p <- ggplot(pf, aes(x = x, y = prob, color = condition)) +
    geom_line(size = 1.2) +
    geom_point() +
    labs(
        x = "Number of methylated CpGs per read",
        y = "Probability",
        title = "Binomial vs. Beta-binomial (different φ)"
    ) +
    facet_wrap(~ncpg, scales = "free_x") +
    mtclosed +
    scale_palette
mysaveandstore(str_glue("overdispersion1192.pdf"), w = 8.9, h = 4)

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "1dispersion_analysis", params$mod_code))
