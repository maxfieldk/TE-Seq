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

grsdf <- read_delim(sprintf("ldna/Rintermediates/%s/grsdf.tsv", params$mod_code), col_names = TRUE)
grsdf %$% sample %>% unique()
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
grs <- GRanges(grsdf)


###########################


refseq_gr <- import(conf$refseq_unaltered)
genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
mcols(genes_gr) %>% colnames()
mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
write_delim(tibble(as.data.frame(promoters)) %>% mutate(score = 1000) %>% dplyr::select(seqnames, start, end, gene_id, score, strand), sprintf("ldna/Rintermediates/%s/promoters.bed", params$mod_code), col_names = FALSE, delim = "\t")


dmrs_per_contrast <- list()
dmls_per_contrast <- list()
dmrsgr_per_contrast <- list()
dmlsgr_per_contrast <- list()
dmrsannot_per_contrast <- list()
dmrsgr_split_per_contrast <- list()
for (contrast in contrasts) {
    dmr_path <- sprintf("ldna/results/%s/tables/%s/dmrs.tsv", params$mod_code, contrast)
    dml_path <- sprintf("ldna/results/%s/tables/%s/dmls.tsv", params$mod_code, contrast)
    dmrs_per_contrast[[contrast]] <- read_delim(dmr_path, delim = "\t", col_names = TRUE) %>% filter(dmr_type %in% c("t01", "t05"))
    dmls_per_contrast[[contrast]] <- read_delim(dml_path, delim = "\t", col_names = TRUE) %>% filter(fdrs <= 0.2)

    dmrsgr_per_contrast[[contrast]] <- GRanges(dmrs_per_contrast[[contrast]])
    dmlsgr_per_contrast[[contrast]] <- GRanges(
        seqnames = dmls_per_contrast[[contrast]]$chr,
        ranges = IRanges(start = dmls_per_contrast[[contrast]]$pos, end = dmls_per_contrast[[contrast]]$pos),
        stat = dmls_per_contrast[[contrast]]$stat,
        pval = dmls_per_contrast[[contrast]]$pvals,
        fdr = dmls_per_contrast[[contrast]]$fdrs,
        direction = dmls_per_contrast[[contrast]]$direction
    )
    dmrsannot_per_contrast[[contrast]] <- dmrs_per_contrast[[contrast]] %>%
        mutate(direction_threshold = paste(direction, gsub("t", "", dmr_type), sep = "_")) %>%
        GRanges()
    dmrsgr_split_per_contrast[[contrast]] <- split(dmrsannot_per_contrast[[contrast]], dmrsannot_per_contrast[[contrast]]$direction_threshold)
}




######### GENES
{
    directions <- c("Hypo", "Hyper", "Dif")
    mydir <- sprintf("ldna/results/%s/plots/great", params$mod_code)
    mydirtables <- sprintf("ldna/results/%s/tables/great", params$mod_code)
    dir.create(mydir, recursive = TRUE, showWarnings = FALSE)
    dir.create(mydirtables, recursive = TRUE, showWarnings = FALSE)
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
    et <- extendTSS(genes_gr %>% resize(width = 1), genome_lengths, gene_id_type = "SYMBOL")
    et_noextension <- extendTSS(genes_gr %>% resize(width = 1), genome_lengths, gene_id_type = "SYMBOL", extension = 0)
}



{
    library(clusterProfiler)
    library(msigdbr)
    gs <- msigdbr("hs")
    # gs %>% filter(grepl("GOBP_HISTONE_H3_K27_TRI", gs_name))
    get_gs_enrichments <- function(gs, gs_ontology_level, outputdir, regionsgrs, etparam, et_mode_string, directions = c("Hypo", "Hyper", "Dif"), background = NULL) {
        outputdirplots <- file.path(outputdir, gs_ontology_level, et_mode_string)
        outputdirtables <- file.path(gsub("plots/great", "tables/great", outputdir), gs_ontology_level, et_mode_string)
        print(outputdirplots)
        dir.create(outputdirplots, recursive = TRUE)
        dir.create(outputdirtables, recursive = TRUE)
        tablesMsigdb <- list()
        genecollections <- gs %>%
            pluck(gs_ontology_level) %>%
            unique()
        for (collection in genecollections) {
            print(collection)
            tryCatch(
                {
                    dir.create(paste(outputdirplots, collection, sep = "/"), recursive = TRUE)
                    dir.create(paste(outputdirtables, collection, sep = "/"), recursive = TRUE)
                    genesets <- gs %>%
                        filter(!!sym(gs_ontology_level) == collection) %>%
                        dplyr::select(gs_name, gene_symbol) %>%
                        dplyr::rename(term = gs_name, gene = gene_symbol) %>%
                        as.data.frame()
                    for (direction in directions) {
                        if (direction == "Hypo") {
                            regionstemp <<- regionsgrs[grepl("Hypo", regionsgrs$direction)]
                        }
                        if (direction == "Hyper") {
                            regionstemp <<- regionsgrs[grepl("Hyper", regionsgrs$direction)]
                        }
                        if (direction == "Dif") {
                            regionstemp <<- regionsgrs
                        }

                        if (is.null(background)) {
                            res <- great(regionstemp, gene_sets = genesets, extended_tss = etparam, background = CHROMOSOMESINCLUDEDINANALYSIS_REF)
                        } else {
                            res <- great(regionstemp, gene_sets = genesets, extended_tss = etparam, background = background)
                        }
                        tb <- getEnrichmentTable(res)
                        if (nrow(tb) != 0) {
                            tb <- tb %>% dplyr::arrange(p_adjust)
                            tablesMsigdb[[collection]][[direction]] <- tb
                            write_delim(tb, paste(outputdirtables, collection, paste0(direction, "great_enrichment.tsv"), sep = "/"))

                            png(paste(outputdirplots, collection, paste0(direction, "volcano.png"), sep = "/"), height = 5, width = 5, res = 300, units = "in")
                            plotVolcano(res)
                            dev.off()

                            png(paste(outputdirplots, collection, paste0(direction, "associations.png"), sep = "/"), height = 5, width = 10, res = 300, units = "in")
                            plotRegionGeneAssociations(res)
                            dev.off()
                        }
                    }
                },
                error = function(e) {
                    print(e)
                }
            )
        }
        # save(file = sprintf("ldna/Rintermediates/%s/tablesMsigdb_%s_%s.rds", params$mod_code, gs_ontology_level, et_mode_string), tablesMsigdb)
        save(file = sprintf("%s/tablesMsigdb.rds", outputdirtables), tablesMsigdb)
        return(tablesMsigdb)
    }

    make_enrich_plots <- function(tempdirectory) {
        tsv_files <- list.files(
            path = tempdirectory,
            pattern = "\\.tsv$",
            recursive = TRUE,
            full.names = TRUE
        )
        for (tempfile in tsv_files) {
            tb <- read_delim(tempfile)
            outputdirplots <- gsub("tables", "plots", dirname(dirname(tempfile)))
            collection <- basename(dirname(tempfile))
            direction <- gsub("great_enrichment.tsv", "", basename(tempfile))

            tbnames <- tb %>%
                tibble() %>%
                mutate(id_nchar = nchar(id)) %>%
                mutate(id = case_when(
                    id_nchar < 40 ~ paste0(strrep("-", pmax(0, 40 - id_nchar)), id),
                    TRUE ~ id
                )) %>%
                mutate(mean_padj = (p_adjust + p_adjust_hyper) / 2)

            binom_df <- tbnames %>%
                dplyr::select(id, fold_enrichment, p_adjust, mean_padj) %>%
                mutate(type = "Binom")
            hyper_df <- tbnames %>%
                dplyr::select(id, fold_enrichment = fold_enrichment_hyper, p_adjust = p_adjust_hyper, mean_padj) %>%
                mutate(type = "Hyper")
            combined_df <- bind_rows(
                binom_df,
                hyper_df %>% mutate(fold_enrichment)
            )

            terms_to_plot <- tbnames %>%
                arrange(mean_padj) %>%
                head(n = 5) %$% id
            p <- combined_df %>%
                filter(id %in% terms_to_plot) %>%
                mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(id = fct_reorder(id, -mean_padj)) %>%
                ggplot(aes(x = id)) +
                geom_col(data = . %>% filter(type == "Binom"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position = position_nudge(x = 0.45 / 2), width = 0.45) +
                coord_flip() +
                scale_fill_distiller(
                    name = "BinomP",
                    palette = "Blues",
                    direction = -1,
                    limits = c(0, 1),
                    oob = scales::squish
                ) +
                new_scale_fill() +
                geom_col(data = . %>% filter(type == "Hyper"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position = position_nudge(x = -0.45 / 2), width = 0.45) +
                scale_fill_distiller(
                    name = "HyperP",
                    palette = "Greens",
                    direction = -1,
                    limits = c(0, 1),
                    oob = scales::squish
                ) +
                geom_vline(xintercept = 0, color = "black") +
                labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                mtclosedgridv +
                anchorbar +
                geom_hline(yintercept = 0, color = "black") +
                theme(axis.text.y = element_text(family = "mono"))
            mysaveandstore(pl = p, fn = paste(outputdirplots, collection, paste0(direction, "lollipop5.pdf"), sep = "/"), 6, 4)

            terms_to_plot <- tbnames %>%
                arrange(mean_padj) %>%
                head(n = 10) %$% id
            p <- combined_df %>%
                filter(id %in% terms_to_plot) %>%
                mutate(id = str_wrap(as.character(id) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(id = fct_reorder(id, -mean_padj)) %>%
                ggplot(aes(x = id)) +
                geom_col(data = . %>% filter(type == "Binom"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position = position_nudge(x = 0.45 / 2), width = 0.45) +
                coord_flip() +
                scale_fill_distiller(
                    name = "BinomP",
                    palette = "Blues",
                    direction = -1,
                    limits = c(0, 1),
                    oob = scales::squish
                ) +
                new_scale_fill() +
                geom_col(data = . %>% filter(type == "Hyper"), aes(y = fold_enrichment, fill = p_adjust, group = type), color = "black", position = position_nudge(x = -0.45 / 2), width = 0.45) +
                scale_fill_distiller(
                    name = "HyperP",
                    palette = "Greens",
                    direction = -1,
                    limits = c(0, 1),
                    oob = scales::squish
                ) +
                geom_vline(xintercept = 0, color = "black") +
                labs(x = "", title = sprintf("%s", collection), subtitle = sprintf("Direction: %s", ifelse(direction == "Dif", "Hypo|Hyper", direction))) +
                mtclosedgridv +
                anchorbar +
                geom_hline(yintercept = 0, color = "black") +
                theme(axis.text.y = element_text(family = "mono"))

            mysaveandstore(pl = p, fn = paste(outputdirplots, collection, paste0(direction, "lollipop.pdf"), sep = "/"), 6.5, 6)
        }
    }

    #####
    for (contrast in contrasts) {
        dmrs <- dmrs_per_contrast[[contrast]]
        dmls <- dmls_per_contrast[[contrast]]
        dmrsgr <- dmrsgr_per_contrast[[contrast]]
        dmlsgr <- dmlsgr_per_contrast[[contrast]]
        dmrsannot <- dmrsannot_per_contrast[[contrast]]
        dmrsgr_split <- dmrsgr_split_per_contrast[[contrast]]

        # # promoters but with background
        #     mydir <- sprintf("ldna/results/%s/plots/great_promoters/%s", params$mod_code, contrast)
        #     mydirtables <- sprintf("ldna/results/%s/tables/great_promoters/%s", params$mod_code, contrast)
        #     for (dmrtype in dmrs$dmr_type %>% unique()) {
        #         regions1 <- mergeByOverlaps(promoters, dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype])
        #         regions2 <- regions1$promoters
        #         mcols(regions2)$direction <- as.data.frame(regions1)$direction
        #         regions <- as.data.frame(regions2) %>%
        #             tibble() %>%
        #             distinct() %>%
        #             GRanges()


        #         tryCatch(
        #             {
        #                 get_gs_enrichments(
        #                     gs = gs,
        #                     gs_ontology_level = "gs_collection",
        #                     outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                     regionsgrs = regions,
        #                     etparam = et_noextension,
        #                     et_mode_string = "et_noextension",
        #                     directions = c("Hyper", "Hypo", "Dif"),
        #                     background = promoters
        #                 )
        #             },
        #             error = function(e) {}
        #         )
        #         tryCatch(
        #             {
        #                 get_gs_enrichments(
        #                     gs = gs,
        #                     gs_ontology_level = "gs_subcollection",
        #                     outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                     regionsgrs = regions,
        #                     etparam = et_noextension,
        #                     et_mode_string = "et_noextension",
        #                     directions = c("Hyper", "Hypo", "Dif"),
        #                     background = promoters
        #                 )
        #             },
        #             error = function(e) {}
        #         )
        #     }

        #     ##### cpg island background
        #     mydir <- sprintf("ldna/results/%s/plots/great_cpgislands/%s", params$mod_code, contrast)
        #     mydirtables <- sprintf("ldna/results/%s/tables/great_cpgislands/%s", params$mod_code, contrast)
        #     for (dmrtype in dmrs$dmr_type %>% unique()) {
        #         regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
        #         tryCatch(
        #             {
        #                 get_gs_enrichments(
        #                     gs = gs,
        #                     gs_ontology_level = "gs_collection",
        #                     outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                     regionsgrs = regions,
        #                     etparam = et,
        #                     et_mode_string = "et_withextension",
        #                     directions = c("Hyper", "Hypo", "Dif"),
        #                     background = cpg_islands
        #                 )
        #             },
        #             error = function(e) {}
        #         )
        #         tryCatch(
        #             {
        #                 get_gs_enrichments(
        #                     gs = gs,
        #                     gs_ontology_level = "gs_subcollection",
        #                     outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                     regionsgrs = regions,
        #                     etparam = et,
        #                     et_mode_string = "et_withextension",
        #                     directions = c("Hyper", "Hypo", "Dif"),
        #                     background = cpg_islands
        #                 )
        #             },
        #             error = function(e) {}
        #         )
        #     }

        ####

        ##### promoters and enhancers background
        chromHMMgr <- import(conf$chromHMM)
        chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]
        prom_no_mcols <- promoters
        mcols(prom_no_mcols) <- NULL
        enh_no_mcols <- chromHMM_enhancers_grs
        mcols(enh_no_mcols) <- NULL
        background <- c(enh_no_mcols, prom_no_mcols)

        mydir <- sprintf("ldna/results/%s/plots/great_prom_enh/%s", params$mod_code, contrast)
        mydirtables <- sprintf("ldna/results/%s/tables/great_prom_enh/%s", params$mod_code, contrast)

        for (dmrtype in dmrs$dmr_type %>% unique()) {
            regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_collection",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
            tryCatch(
                {
                    get_gs_enrichments(
                        gs = gs,
                        gs_ontology_level = "gs_subcollection",
                        outputdir = sprintf("%s/%s", mydir, dmrtype),
                        regionsgrs = regions,
                        etparam = et,
                        et_mode_string = "et_withextension",
                        directions = c("Hyper", "Hypo", "Dif"),
                        background = background
                    )
                },
                error = function(e) {}
            )
        }

        make_enrich_plots(sprintf("ldna/results/m/tables/great_prom_enh/%s", contrast))



        # ##### promoters and enhancers island background
        # chromHMMgr <- import(conf$chromHMM)
        # chromHMM_enhancers_grs <- chromHMMgr[grepl("Enh*", mcols(chromHMMgr)$name)]
        # prom_no_mcols <- promoters
        # mcols(prom_no_mcols) <- NULL
        # enh_no_mcols <- chromHMM_enhancers_grs
        # mcols(enh_no_mcols) <- NULL
        # prom_enh <- c(enh_no_mcols, prom_no_mcols)

        # background <- cpg_islands %>% subsetByOverlaps(prom_enh)
        # mydir <- sprintf("ldna/results/%s/plots/great_prom_enh_intersect_cpgI/%s", params$mod_code, contrast)
        # mydirtables <- sprintf("ldna/results/%s/tables/great_prom_enh_intersect_cpgI/%s", params$mod_code, contrast)
        # for (dmrtype in dmrs$dmr_type %>% unique()) {
        #     regions <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
        #     tryCatch(
        #         {
        #             get_gs_enrichments(
        #                 gs = gs,
        #                 gs_ontology_level = "gs_collection",
        #                 outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                 regionsgrs = regions,
        #                 etparam = et,
        #                 et_mode_string = "et_withextension",
        #                 directions = c("Hyper", "Hypo", "Dif"),
        #                 background = background
        #             )
        #         },
        #         error = function(e) {}
        #     )
        #     tryCatch(
        #         {
        #             get_gs_enrichments(
        #                 gs = gs,
        #                 gs_ontology_level = "gs_subcollection",
        #                 outputdir = sprintf("%s/%s", mydir, dmrtype),
        #                 regionsgrs = regions,
        #                 etparam = et,
        #                 et_mode_string = "et_withextension",
        #                 directions = c("Hyper", "Hypo", "Dif"),
        #                 background = background
        #             )
        #         },
        #         error = function(e) {}
        #     )
        # }
    }
}

####


{
    grsdf_island_status <- grsdf %>%
        group_by(sample, islandStatus) %>%
        mutate(methylated_sites = round(cov * pctM / 100)) %>%
        summarise(methylated_sites = sum(methylated_sites), cov = sum(cov)) %>%
        mutate(pctM = methylated_sites / cov)

    dat_tmp <- grsdf_island_status %>%
        ungroup() %>%
        dplyr::rename(sample_name = sample) %>%
        left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))
    genomewide_model_islandstatus <- glmmTMB(
        cbind(methylated_sites, cov - methylated_sites) ~
            condition * islandStatus + sex + age_z + (1 | sample_name),
        data = dat_tmp,
        family = binomial()
    )
    stats <- broom::tidy(genomewide_model_islandstatus)

    p <- grsdf_island_status %>%
        mutate(pctM = methylated_sites / cov) %>%
        left_join(sample_table) %>%
        ggplot() +
        geom_point(aes(x = islandStatus, y = pctM, color = condition), position = position_dodge(width = 0.75)) +
        mtopen +
        scale_conditions +
        anchorbar
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genomewide/mean_meth_islandstatus_genomewide.pdf", params$mod_code), 5, 4, sf = stats, pl = p)
}

{
    island_promoters <- merge_with_grs(grs[mcols(grs)$islandStatus == "island"], promoters)

    ip <- island_promoters %>%
        group_by(gene_id, sample, condition) %>%
        summarise(mean_meth = mean(pctM))


    ip_mean <- island_promoters %>%
        group_by(sample, condition) %>%
        summarise(mean_meth = mean(pctM)) %>%
        left_join(sample_table %>% mutate(sample = sample_name, age_z = as.numeric(scale(age))))


    #### multilevel STATS
    dat_tmp <- island_promoters %>%
        mutate(methylated_sites = round(cov * pctM / 100), total_sites = cov) %>%
        dplyr::rename(sample_name = sample) %>%
        group_by(gene_id, sample_name) %>%
        summarise(methylated_sites = sum(methylated_sites), total_sites = sum(total_sites)) %>%
        left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))


    p <- ip_mean %>%
        mutate(condition = factor(condition, levels = conf$levels)) %>%
        ggviolin(x = "condition", y = "mean_meth", fill = "condition", add = c("mean_se", "dotplot")) +
        scale_conditions +
        mtopen
    mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_violin.pdf", params$mod_code), 4.5, 3.75)

    p <- ip_mean %>%
        ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
        geom_point(size = 3) +
        scale_conditions +
        geom_vline(xintercept = ip_mean %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
        geom_vline(xintercept = ip_mean %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
        geom_text_repel(aes(label = apoe)) +

        # new_scale_fill() +
        # geom_tile(aes(x = -1, fill = sex), width = 2) + # Add metadata strip
        scale_conditions +
        mtopen

    library(broom)
    if (enough_samples_per_condition_for_stats) {
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()

        gene_model <- glmmTMB(
            cbind(methylated_sites, total_sites - methylated_sites) ~
                condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
            data = dat_tmp,
            family = binomial()
        )
        stats <- broom::tidy(gene_model)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide.pdf", params$mod_code), 5, 4, sf = stats)
    } else {
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide.pdf", params$mod_code), 5, 3.7)
    }

    island_promotersnoX <- island_promoters %>% filter(seqnames != "chrX")

    ipnoX <- island_promotersnoX %>%
        group_by(gene_id, sample, condition) %>%
        summarise(mean_meth = mean(pctM))


    ip_meannoX <- island_promotersnoX %>%
        group_by(sample, condition) %>%
        summarise(mean_meth = mean(pctM)) %>%
        left_join(sample_table %>% mutate(sample = sample_name, age_z = as.numeric(scale(age))))


    #### multilevel STATS THIS UP
    dat_tmp <- island_promotersnoX %>%
        mutate(methylated_sites = round(cov * pctM / 100), total_sites = cov) %>%
        dplyr::rename(sample_name = sample) %>%
        group_by(gene_id, sample_name) %>%
        summarise(methylated_sites = sum(methylated_sites), total_sites = sum(total_sites)) %>%
        left_join(sample_table %>% mutate(age_z = as.numeric(scale(age))))


    p <- ip_meannoX %>%
        ggplot(aes(y = sample_name, x = mean_meth, color = condition, shape = !!sym(asc))) +
        geom_point(size = 3) +
        scale_conditions +
        geom_vline(xintercept = ip_mean %>% filter(condition == condition2) %$% mean_meth %>% mean(), color = "blue", linetype = "dashed") +
        geom_vline(xintercept = ip_mean %>% filter(condition == condition1) %$% mean_meth %>% mean(), color = "grey", linetype = "dashed") +
        geom_text_repel(aes(label = apoe)) +
        scale_conditions +
        mtopen
    library(broom)
    if (enough_samples_per_condition_for_stats) {
        # stats <- summary(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), pf)) %>% tidy()

        gene_model <- glmmTMB(
            cbind(methylated_sites, total_sites - methylated_sites) ~
                condition + sex + age_z + (1 | sample_name) + (1 | gene_id),
            data = dat_tmp,
            family = binomial()
        )
        stats <- broom::tidy(gene_model)
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_noX.pdf", params$mod_code), 5, 4, sf = stats)
    } else {
        mysaveandstore(fn = sprintf("ldna/results/%s/plots/genes/mean_meth_point_genomewide_noX.pdf", params$mod_code), 5, 3.7)
    }

    #####

    results <- ip %>%
        left_join(sample_table) %>%
        group_by(gene_id) %>%
        summarise(
            model = list(lm(formula(sprintf("%s ~ %s", "mean_meth", lm_right_hand_side)), data = cur_data())),
            .groups = "drop"
        ) %>%
        mutate(
            tidied = map(model, broom::tidy)
        ) %>%
        unnest(tidied) %>%
        filter(term == paste0("condition", condition1)) %>%
        dplyr::select(gene_id, estimate, p.value) %>%
        mutate(
            padj = p.adjust(p.value, method = "fdr") # Adjust p-values for multiple testing
        )

    res <- results %>%
        mutate(signed_log10p = sign(estimate) * abs(log10(p.value))) %>%
        arrange(-signed_log10p)
    ordered_by_stat <- setNames(res[["signed_log10p"]], res$gene_id) %>% na.omit()
    # GSEA untargeted
    tryCatch(
        {
            gene_sets <- msigdbr(species = confALL$aref$species)
        },
        error = function(e) {
            gene_sets <<- msigdbr(species = "human")
        }
    )
    library(clusterProfiler)
    rm(gse_df)

    for (category in gene_sets %$% gs_collection %>% unique()) {
        cat(category, "\n")
        tryCatch({
            collection <- category
            msigdbr_df <- gene_sets %>% filter(gs_collection == category)
            msigdbr_t2g <- msigdbr_df %>%
                dplyr::distinct(gs_name, gene_symbol) %>%
                as.data.frame()
            gse <- GSEA(ordered_by_stat, TERM2GENE = msigdbr_t2g, maxGSSize = 100000, minGSSize = 1)
            df <- gse@result %>% tibble()
            df$collection <- collection
            df$contrast <- str_glue("condition_{condition2}_vs_{condition1}")
            # gse_results[[contrast]][[collection]] <- as.data.frame(df) %>% tibble()
            if (!exists("gse_df")) {
                gse_df <<- df
            } else {
                gse_df <<- rbind(gse_df, df)
            }
            genesettheme <- theme_gray() + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

            tryCatch(
                {
                    for (num in c(5, 10, 15)) {
                        dftemp <- arrange(df, -abs(NES)) %>%
                            group_by(sign(NES)) %>%
                            slice(1:num) %>%
                            mutate(log10padj = -log10(p.adjust))
                        dftemp <- dftemp %>% mutate(Description = str_wrap(as.character(Description) %>% gsub("_", " ", .), width = 40))
                        dftemp <- dftemp %>% mutate(`Gene Ratio` = 0.01 * as.numeric(gsub("%", "", gsub(",.*", "", gsub("tags=", "", leading_edge)))))
                        p <- ggplot(dftemp, aes(NES, fct_reorder(Description, NES), fill = log10padj)) +
                            geom_col(orientation = "y") +
                            scale_fill_gradient(high = "red", low = "white", limits = c(0, 5), oob = scales::squish) +
                            mtopen +
                            theme(axis.text.y = element_text(colour = "black")) +
                            ylab(NULL) +
                            labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                        mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)

                        p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust) * `sign(NES)`)) +
                            geom_col(orientation = "y") +
                            scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                            mtopen +
                            theme(axis.text.y = element_text(colour = "black")) +
                            ylab(NULL) +
                            guides(fill = guide_legend(title = "Signed \n-log10(p.adjust)")) +
                            labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                        mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                    }
                },
                error = function(e) {
                    print("")
                }
            )
        })
    }

    # ipcondition <- ip %>%
    #     group_by(gene_id, condition) %>%
    #     mutate(mean_meth = mean(mean_meth)) %>%
    #     group_by(gene_id, condition) %>%
    #     mutate(nrow = row_number()) %>%
    #     filter(nrow == 1) %>%
    #     dplyr::select(-nrow, -sample) %>%
    #     pivot_wider(id_cols = gene_id, names_from = condition, values_from = mean_meth) %>%
    #     mutate(dif = !!sym(condition2) - !!sym(condition1))


    # ipcondition %>%
    #     arrange(dif) %>%
    #     head(n = 30) %>%
    #     print(n = 30)
    # ipcondition %>%
    #     arrange(-dif) %>%
    #     head(n = 30) %>%
    #     print(n = 30)


    # ipwide <- ip %>%
    #     dplyr::select(-condition) %>%
    #     pivot_wider(names_from = sample, values_from = mean_meth)

    # ipwide %>% filter(gene_id == "LOC112268283")

    # library(dplyr)
    # library(broom)

    # results <- ip %>%
    #     group_by(gene_id) %>%
    #     summarise(
    #         model = list(lm(mean_meth ~ condition, data = cur_data())),
    #         .groups = "drop"
    #     ) %>%
    #     mutate(
    #         tidied = map(model, broom::tidy)
    #     ) %>%
    #     unnest(tidied) %>%
    #     filter(term == paste0("condition", condition1)) %>%
    #     dplyr::select(gene_id, estimate, p.value) %>%
    #     mutate(
    #         log2FC = log2(exp(estimate)), # Convert from natural log to log2 fold-change
    #         padj = p.adjust(p.value, method = "BH") # Adjust p-values
    #     )

    # results %>% filter(p.value < 0.05)

    # ipwide %>% filter(gene_id == "AGBL4-AS1")

    # # View results
    # head(results)
    ## END NEW CODE



    # threshold_dfs <- list()
    # for (dmrtype in dmrs$dmr_type %>% unique()) {
    #     dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
    #     hyporegions <- dmrsgr_temp[grepl("Hypo", dmrsgr_temp$direction)]
    #     hyperregions <- dmrsgr_temp[grepl("Hyper", dmrsgr_temp$direction)]

    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, dmrsgr_temp, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyporegions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hypo_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(promoters, hyperregions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/promoters_hyper_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")

    #     hypo <- mcols(subsetByOverlaps(promoters, hyporegions))$gene_id
    #     hyper <- mcols(subsetByOverlaps(promoters, hyperregions))$gene_id
    #     disc <- intersect(hypo, hyper)
    #     hypo_nd <- setdiff(hypo, disc)
    #     hyper_nd <- setdiff(hyper, disc)



    #     threshold_df <- bind_rows(
    #         tibble(gene_id = hypo_nd, !!sym(dmrtype) := "Hypo"),
    #         tibble(gene_id = hyper_nd, !!sym(dmrtype) := "Hyper"),
    #         tibble(gene_id = disc, !!sym(dmrtype) := "Discordant")
    #     )
    #     threshold_dfs[[dmrtype]] <- threshold_df
    # }
    # thresholddf <- purrr::reduce(threshold_dfs, full_join)
    # promoters_df <- full_join(tibble(as.data.frame(promoters)), thresholddf) %>%
    #     filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
    #     pivot_longer(cols = c("t05", "t01"), names_to = "dmr_type", values_to = "direction") %>%
    #     mutate(direction = factor(direction, levels = c("Discordant", "Hyper", "Hypo")))

    # total_possible <- nrow(mcols(promoters))

    # p <- promoters_df %>%
    #     group_by(dmr_type, direction) %>%
    #     summarise(n = n()) %>%
    #     ungroup() %>%
    #     tidyr::complete(dmr_type, direction, fill = list(n = 0)) %>%
    #     filter(!is.na(direction)) %>%
    #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
    #     ggplot() +
    #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
    #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
    #     ggtitle("Promoter Methylation") +
    #     labs(x = "Direction", y = str_glue("Count (out of {total_possible})")) +
    #     theme(legend.position = "none") +
    #     annotation_logticks(sides = "l") +
    #     scale_y_log10(
    #         breaks = scales::trans_breaks("log10", function(x) 10^x),
    #         labels = scales::trans_format("log10", math_format(10^.x))
    #     ) +
    #     mtopen +
    #     scale_methylation_thresholds
    # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/genes_concordance.pdf", params$mod_code), 5, 4)

    # threshold_dfs <- list()
    # for (dmrtype in dmrs$dmr_type %>% unique()) {
    #     dmrsgr_temp <- dmrsgr[mcols(dmrsgr)$dmr_type == dmrtype]
    #     hyporegions <- dmrsgr_temp[grepl("Hypo", dmrsgr_temp$direction)]
    #     hyperregions <- dmrsgr_temp[grepl("Hyper", dmrsgr_temp$direction)]
    #     chromHMM_enhancers_grs_with_id <- chromHMM_enhancers_grs
    #     mcols(chromHMM_enhancers_grs_with_id)$ID <- paste0("ID", 1:nrow(mcols(chromHMM_enhancers_grs)))
    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, dmrsgr_temp, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, hyporegions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_hypo_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")
    #     write_delim(tibble(as.data.frame(GenomicRanges::intersect(chromHMM_enhancers_grs_with_id, hyperregions, ignore.strand = TRUE))), sprintf("ldna/Rintermediates/%s/chromHMM_enhancers_hyper_%s.bed", params$mod_code, dmrtype), col_names = FALSE, delim = "\t")

    #     hypo <- mcols(subsetByOverlaps(chromHMM_enhancers_grs_with_id, hyporegions))$ID
    #     hyper <- mcols(subsetByOverlaps(chromHMM_enhancers_grs_with_id, hyperregions))$ID
    #     disc <- intersect(hypo, hyper)
    #     hypo_nd <- setdiff(hypo, disc)
    #     hyper_nd <- setdiff(hyper, disc)

    #     threshold_df <- bind_rows(
    #         tibble(ID = hypo_nd, !!sym(dmrtype) := "Hypo"),
    #         tibble(ID = hyper_nd, !!sym(dmrtype) := "Hyper"),
    #         tibble(ID = disc, !!sym(dmrtype) := "Discordant")
    #     )
    #     threshold_dfs[[dmrtype]] <- threshold_df
    # }
    # thresholddf <- purrr::reduce(threshold_dfs, full_join)
    # enh_df <- full_join(tibble(as.data.frame(chromHMM_enhancers_grs_with_id)), thresholddf) %>%
    #     filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS) %>%
    #     pivot_longer(cols = c("t05", "t01"), names_to = "dmr_type", values_to = "direction") %>%
    #     mutate(direction = factor(direction, levels = c("Discordant", "Hyper", "Hypo")))

    # total_possible <- nrow(mcols(chromHMM_enhancers_grs))

    # p <- enh_df %>%
    #     group_by(dmr_type, direction) %>%
    #     summarise(n = n()) %>%
    #     ungroup() %>%
    #     tidyr::complete(dmr_type, direction, fill = list(n = 0)) %>%
    #     filter(!is.na(direction)) %>%
    #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
    #     ggplot() +
    #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
    #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = direction, y = n, fill = direction_threshold), color = "black") +
    #     ggtitle("Enhancer Methylation") +
    #     labs(x = "Direction", y = str_glue("Count (out of {total_possible})")) +
    #     theme(legend.position = "none") +
    #     annotation_logticks(sides = "l") +
    #     scale_y_log10(
    #         breaks = scales::trans_breaks("log10", function(x) 10^x),
    #         labels = scales::trans_format("log10", math_format(10^.x))
    #     ) +
    #     mtopen +
    #     scale_methylation_thresholds
    # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/enhancer_concordance.pdf", params$mod_code), 5, 4)
    # #####
    # enh_df %>% bind_rows()
    # prom_enh_df <- promoters_df %>%
    #     mutate(biotype = "Promoter") %>%
    #     filter(!is.na(direction)) %>%
    #     dplyr::select(direction, biotype, dmr_type, ID) %>%
    #     bind_rows(enh_df %>% mutate(biotype = "Enhancer") %>% filter(!is.na(direction)) %>% dplyr::select(direction, dmr_type, biotype, ID))
    # p <- prom_enh_df %>%
    #     group_by(dmr_type, direction, biotype) %>%
    #     summarise(n = n()) %>%
    #     ungroup() %>%
    #     tidyr::complete(dmr_type, direction, biotype, fill = list(n = 0)) %>%
    #     filter(direction != "Discordant") %>%
    #     filter(!is.na(direction)) %>%
    #     mutate(direction_threshold = paste0(direction, "_", gsub("t", "", dmr_type))) %>%
    #     ggplot() +
    #     geom_col(data = . %>% filter(dmr_type == "t05"), aes(x = biotype, y = n, fill = direction_threshold), color = "black", position = "dodge2") +
    #     geom_col(data = . %>% filter(dmr_type == "t01"), aes(x = biotype, y = n, fill = direction_threshold), color = "black", position = "dodge2") +
    #     ggtitle("Enhancer Methylation") +
    #     labs(x = "Direction", y = str_glue("Count")) +
    #     theme(legend.position = "none") +
    #     annotation_logticks(sides = "l") +
    #     scale_y_log10(
    #         breaks = scales::trans_breaks("log10", function(x) 10^x),
    #         labels = scales::trans_format("log10", math_format(10^.x))
    #     ) +
    #     mtopen +
    #     scale_methylation_thresholds
    # mysaveandstore(pl = p, sprintf("ldna/results/%s/plots/genes/prom_enhancer_concordance.pdf", params$mod_code), 5, 4)

    # # dmrsgr_enh <- dmrsgr %>% subsetByOverlaps(chromHMM_enhancers_grs)
    # # dmrsgr_prom <- dmrsgr %>% subsetByOverlaps(promoters)
    # # length_overlap <- dmrsgr_enh %>% subsetByOverlaps(dmrsgr_prom) %>% mcols() %>% nrow()
    # # perc_prom_in_enh <- length_overlap/nrow(mcols(dmrsgr_prom))
    # # perc_enh_in_prom <- length_overlap/nrow(mcols(dmrsgr_enh))
    # enh_in_dmrs <- chromHMM_enhancers_grs %>% subsetByOverlaps(dmrsgr)
    # prom_in_dmrs <- promoters %>% subsetByOverlaps(dmrsgr)
    # enh_intersect <- enh_in_dmrs %>%
    #     subsetByOverlaps(prom_in_dmrs) %>%
    #     mcols() %>%
    #     nrow()
    # perc_prom_in_enh <- enh_intersect / nrow(mcols(dmrsgr_prom))
    # prom_intersect <- prom_in_dmrs %>%
    #     subsetByOverlaps(enh_in_dmrs) %>%
    #     mcols() %>%
    #     nrow()
    # perc_enh_in_prom <- prom_intersect / nrow(mcols(dmrsgr_enh))

    #####


    # for (dmrtype in dmrs$dmr_type %>% unique()) {
    #     library(clusterProfiler)
    #     hypo_genes <- promoters_df %>% filter(direction == "Hypo") %$% gene_id
    #     hyper_genes <- promoters_df %>% filter(direction == "Hyper") %$% gene_id
    #     background <- mcols(promoters)$gene_id

    #     gs <- msigdbr("human")
    #     tablesORA_subcollection <- list()
    #     results_ORA_hypo <- list()
    #     results_ORA_hyper <- list()
    #     genesubcollections <- gs$gs_subcollection %>% unique()
    #     for (collection in genesubcollections) {
    #         term2gene <- gs %>%
    #             filter(gs_subcollection == collection) %>%
    #             dplyr::rename(term = gs_name, gene = gene_symbol) %>%
    #             select(term, gene)
    #         res_hypo <- enricher(hypo_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

    #         results_ORA_hypo[[collection]] <- res_hypo %>% as.data.frame()

    #         res_hyper <- enricher(hyper_genes, universe = background, TERM2GENE = term2gene, pAdjustMethod = "fdr")

    #         results_ORA_hyper[[collection]] <- res_hyper %>% as.data.frame()
    #     }
}

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "16_great", params$mod_code))
