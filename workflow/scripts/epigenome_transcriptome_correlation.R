library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(purrr)
library(GenomicRanges)
library(AnnotationDbi)
library(zoo)
library(rtracklayer)
library(Biostrings)
library(GGally)


conf <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "outputdir" = "integrated/genomebrowserplots/dorado",
            "regions_of_interest" = "conf/integrated_regions_of_interest.bed",
            "r_annotation_fragmentsjoined" = conf[["lrna"]]$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf[["lrna"]]$r_repeatmasker_annotation,
            "txdb" = conf[["lrna"]]$txdb
        ), env = globalenv())
        assign("inputs", list(
            "srna_results" = "srna/results/agg/deseq/resultsdf.tsv",
            "lrna_results" = "lrna/results/agg/deseq/dorado/relaxed/resultsdf.tsv",
            "ldna_methylation" = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", conf[["ldna"]][["samples"]], conf[["ldna"]][["samples"]]),
            "rteprommeth" = "ldna/Rintermediates/perelementdf_promoters.tsv"
            "dmrs" = "ldna/results/tables/dmrs.CG_m.tsv",
            "dmls" = "ldna/results/tables/dmrs.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "plots" = "integrated/epigenome_transcriptome_correlation/objects/plots.rda"
        ), env = globalenv())
    }
)
## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
r_repeatmasker_annotation %$% ltr_viral_status_req %>% unique()
RMdf <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
RM <- GRanges(RMdf)



srna_samples <- conf[["srna"]]$samples
srna_conditions <- conf[["srna"]]$levels
srna_contrasts <- conf[["srna"]]$contrasts
srna_sample_table <- read_csv("conf/sample_table_srna.csv")
srna_df <- read_delim(inputs$srna_results, delim = "\t")
colnames(srna_df) <- paste0("srna_", colnames(srna_df)) %>% gsub("srna_gene_id", "gene_id", .)
srna_genes <- srna_df %>% filter(srna_gene_or_te == "gene")
srna_repeats <- srna_df %>% filter(srna_gene_or_te == "repeat")


lrna_samples <- conf[["lrna"]]$samples
lrna_conditions <- conf[["lrna"]]$levels
lrna_contrasts <- conf[["lrna"]]$contrasts
lrna_sample_table <- read_csv("conf/sample_table_lrna.csv")
lrna_df <- read_delim(inputs$lrna_results, delim = "\t")
colnames(lrna_df) <- paste0("lrna_", colnames(lrna_df)) %>%
    gsub("pro", "PRO", .) %>%
    gsub("sen", "SEN", .) %>%
    gsub("lrna_gene_id", "gene_id", .)
lrna_genes <- lrna_df %>% filter(lrna_gene_or_te == "gene")
lrna_genes$gene_id <- lrna_genes$gene_id %>% gsub("gene-", "", .)
lrna_repeats <- lrna_df %>% filter(lrna_gene_or_te == "repeat")

ldna_samples <- conf[["ldna"]]$samples
ldna_conditions <- conf[["ldna"]]$conditions
ldna_contrasts <- conf[["ldna"]]$contrasts
ldna_sample_table <- read_csv("conf/sample_table_ldna.csv")


join1 <- full_join(RMdf, lrna_repeats, by = "gene_id")
repeats <- full_join(join1, srna_repeats, by = "gene_id")
genes <- full_join(srna_genes, lrna_genes, by = "gene_id")

conditions_to_plot <- conf[["lrna"]]$levels
dnamethylationlist <- list()
for (condition in conditions_to_plot) {
    path <- grep(sprintf("%s", condition), inputs$ldna_methylation, value = TRUE)
    df <- read_delim(path, col_names = TRUE, delim = "\t") %>%
        mutate(start = pos) %>%
        mutate(end = pos + 1) %>%
        mutate(pctM = X / N) %>%
        dplyr::select(chr, start, end, pctM)
    dnamethylationlist[[condition]] <- df
}


methdif1 <- left_join(dnamethylationlist[[condition1]] %>% dplyr::rename(condition1 = pctM), dnamethylationlist[[condition2]] %>% dplyr::rename(condition2 = pctM))
methdif <- methdif1 %>%
    mutate(dif = condition2 - condition1)




conditions_to_plot <- conf[["lrna"]]$levels
condition1 <- conditions_to_plot[1]
condition2 <- conditions_to_plot[2]



shared_conditions <- intersect(intersect(srna_conditions, lrna_conditions), ldna_conditions)
shared_contrasts <- intersect(intersect(srna_contrasts, lrna_contrasts), ldna_contrasts)

srna_samples <- srna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
lrna_samples <- lrna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name
ldna_samples <- ldna_sample_table %>% filter(condition %in% shared_conditions) %$% sample_name

dir.create("results/plots/rte", showWarnings = FALSE)

plots <- list()



### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "_.*family")]
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

outputdir <- dirname(outputs$outfile)
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)



# txdb <- loadDb(params$txdb)
# columns(txdb)
# keys <- keys(txdb) %>% head()
# AnnotationDbi::select(txdb, keys = keys, columns = "TXNAME", keytype = "GENEID")
# library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
# assembly <- assembly(Genome = "custom", TxDb = txdb, OrgDb = "org.Hs.eg.db", BSgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0, gene.id.column = "GENEID", display.column = "GENEID")

#### srna to lrna correlation

for (gene_or_te in c("gene", "repeat")) {
    srna_s <- srna_sample_table %>%
        dplyr::filter(condition %in% shared_conditions) %>%
        pull(sample_name)
    lrna_s <- lrna_sample_table %>%
        dplyr::filter(condition %in% shared_conditions) %>%
        pull(sample_name)
    ldna_s <- ldna_sample_table %>%
        dplyr::filter(condition %in% shared_conditions) %>%
        pull(sample_name)

    if (gene_or_te == "gene") {
        tdf <- genes
    } else {
        tdf <- repeats
        matdf <- tdf %>% dplyr::select(gene_id, starts_with(paste0("lrna_", lrna_s)))
        mat <- matdf %>%
            dplyr::select(-gene_id) %>%
            as.matrix()
        rownames(mat) <- matdf$gene_id
        sum(rowSums(mat) > length(lrna_s) / length(shared_conditions), na.rm = TRUE)
        mask <- rowSums(mat) > length(lrna_s) / length(shared_conditions)
        # turn na into FALSE
        mask[is.na(mask)] <- FALSE
        matfiltered <- mat[mask, ]
        keep <- rownames(matfiltered)
        tdf <- tdf %>% filter(gene_id %in% keep)
    }


    # All genes
    pf <- tdf %>%
        dplyr::select(starts_with(paste0("srna_", srna_s)), starts_with(paste0("lrna_", lrna_s))) %>%
        mutate(across(where(is.numeric), ~ log2(.x + 1)))

    corr <- round(cor(pf, use = "pairwise.complete.obs"), 2)
    p <- ggcorrplot(corr,
        hc.order = TRUE, type = "lower",
        lab = TRUE
    )
    mysave(sprintf("%s/%s/cormatrix.png", outputdir, gene_or_te), 7, 7)
    plots[["cormatrix"]][[gene_or_te]][["all"]] <- p
    # DEGs
    for (contrast in shared_contrasts) {
        contrast_of_interest <- contrast
        contrast_level_2 <- contrast_of_interest %>%
            gsub("condition_", "", .) %>%
            gsub("_vs_.*", "", .)
        contrast_level_1 <- contrast_of_interest %>%
            gsub(".*_vs_", "", .)
        contrast_stat <- paste0("stat_", contrast_of_interest)
        contrast_padj <- paste0("padj_", contrast_of_interest)
        contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)

        pf <- tdf %>%
            filter((!!sym(paste0("srna_", contrast_padj)) < 0.05) | (!!sym(paste0("lrna_", contrast_padj)) < 0.05)) %>%
            dplyr::select(starts_with(paste0("srna_", srna_s)), starts_with(paste0("lrna_", lrna_s))) %>%
            dplyr::select(contains(contrast_level_1) | contains(contrast_level_2)) %>%
            mutate(across(where(is.numeric), ~ log2(.x + 1)))

        corr <- round(cor(pf, use = "pairwise.complete.obs"), 2)
        p <- ggcorrplot(corr,
            hc.order = FALSE, type = "lower",
            lab = TRUE
        )
        mysave(sprintf("%s/%s/cormatrix_degs_%s.png", outputdir, gene_or_te, contrast), 7, 7)
        plots[["cormatrix_degs"]][[gene_or_te]][[contrast]] <- p
    }
    # On a per condition basis
    for (cond in shared_conditions) {
        srna_s <- srna_sample_table %>%
            dplyr::filter(condition == cond) %>%
            pull(sample_name)
        lrna_s <- lrna_sample_table %>%
            dplyr::filter(condition == cond) %>%
            pull(sample_name)
        ldna_s <- ldna_sample_table %>%
            dplyr::filter(condition == cond) %>%
            pull(sample_name)

        pf <- tdf %>%
            dplyr::select(starts_with(paste0("srna_", srna_s)), starts_with(paste0("lrna_", lrna_s))) %>%
            mutate(across(where(is.numeric), ~ log2(.x + 1)))
        p <- ggpairs(pf)
        mysave(sprintf("%s/%s/correlation_%s.png", outputdir, gene_or_te, cond))
        plots[["correlation"]][[gene_or_te]][[cond]] <- p

        corr <- round(cor(pf, use = "pairwise.complete.obs"), 1)
        p <- ggcorrplot(corr)
        mysave(sprintf("%s/%s/cormatrix_%s.png", outputdir, gene_or_te, cond))
        plots[["cormatrix"]][[gene_or_te]][[cond]] <- p
    }
}





# dna methylation




# Initialize bedmethanalysis data
{
    confglobal <- conf
    conf <- conf$ldna
    sample_table <- ldna_sample_table

{
genome_lengths <- fasta.seqlengths(conf$reference)
chromosomesAll <- names(genome_lengths)
nonrefchromosomes <- grep("^NI", chromosomesAll, value = TRUE)
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

    conditions <- conf$conditions
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

    grsdf <- read_delim("ldna/Rintermediates/grsdf.tsv", col_names = TRUE)
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
    dmls <- dmls %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
    dmrs <- dmrs %>% filter(chr %in% CHROMOSOMESINCLUDEDINANALYSIS)
    dmrs <- dmrs %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
    dmrs$direction <- factor(dmrs$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))

    dmls <- dmls %>% mutate(direction = ifelse(diff_c2_minus_c1 > 0, paste0(condition2, " Hyper"), paste0(condition2, " Hypo")))
    dmls$direction <- factor(dmls$direction, levels = c(paste0(condition2, " Hypo"), paste0(condition2, " Hyper")))


    dmrsgr <- GRanges(dmrs)

    refseq_gr <- import(conf$refseq_unaltered)
    genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
    genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
    genes_gr <- genes_gr[mcols(genes_gr)[, "source"] %in% c("BestRefSeq", "Curated Genomic", "Gnomon"), ]
    mcols(genes_gr)$gene_id <- mcols(genes_gr)$Name
    mcols(genes_gr) %>% colnames()
    mcols(genes_gr) <- mcols(genes_gr)[, c("gene_id", "ID", "gene_biotype", "source")]
    promoters <- promoters(genes_gr, upstream = 5000, downstream = 1000)
    conf <- confglobal
}

# load gene promoter methylation data
mbo <- mergeByOverlaps(grs, promoters)
methgr <- mbo$grs
methgr$source <- mbo$source
methgr$gene_id <- mbo$Name
methgr$islandStatus <- mbo$islandStatus
promotergrs <- tibble(as.data.frame(methgr))
pergenedf <- promotergrs %>%
    group_by(gene_id, condition) %>%
    summarise(meanmeth = mean(pctM)) %>%
    pivot_wider(names_from = condition, values_from = meanmeth) %>%
    ungroup()
colnames(pergenedf) <- paste0("ldna_", colnames(pergenedf)) %>% gsub("ldna_gene_id", "gene_id", .)
pgdf <- left_join(pergenedf, genes, by = "gene_id")

# load rte promoter methylation data
rteprommeth <- read_delim(inputs$rteprommeth)
perrepeatdf <- rteprommeth %>%
    group_by(gene_id, condition) %>%
    summarise(meanmeth = mean(mean_meth)) %>%
    pivot_wider(names_from = condition, values_from = meanmeth) %>%
    ungroup()
colnames(perrepeatdf) <- paste0("ldna_", colnames(perrepeatdf)) %>% gsub("ldna_gene_id", "gene_id", .)
prdf <- left_join(perrepeatdf, repeats, by = "gene_id")

frames <- list("RefSeq Genes" = pgdf)
for (element_type in prdf %$% rte_subfamily_limited %>% unique()) {
    if (element_type != "Other") {
        frames[[element_type]] <- prdf %>% filter(rte_subfamily_limited == element_type)
    }
}

for (framename in names(frames)) {
    tryCatch(
        {
            frame <- frames[[framename]]

            pf <- frame %>%
                dplyr::select(starts_with(paste0("ldna_")), starts_with(paste0("srna_", srna_s))) %>%
                mutate(across(starts_with("srna"), ~ log2(.x + 1)))

            p <- ggpairs(pf) + ggtitle(framename)
            mysave(sprintf("%s/%s/ggpairs_%s.png", outputdir, "srna_ldna", framename), 10, 10)
            plots[["srna_ldna"]][["ggpairs"]][[framename]] <- p

            pf <- frame %>%
                dplyr::select(starts_with(paste0("ldna_")), starts_with(paste0("lrna_", lrna_s))) %>%
                mutate(across(starts_with("lrna"), ~ log2(.x + 1)))

            p <- ggpairs(pf) + ggtitle(framename)
            mysave(sprintf("%s/%s/ggpairs_%s.png", outputdir, "lrna_ldna", framename), 10, 10)
            plots[["lrna_ldna"]][["ggpairs"]][[framename]] <- p

            pf <- frame %>%
                mutate(ldna_cond2mcond1 = !!sym(paste0("ldna_", condition2)) - !!sym(paste0("ldna_", condition1))) %>%
                dplyr::select(ldna_cond2mcond1, paste0("srna_", contrast_log2FoldChange))
            p <- ggpairs(pf) + ggtitle(framename)
            mysave(sprintf("%s/%s/ggpairs_dif_l2fc_%s.png", outputdir, "srna_ldna", framename), 4, 4)
            plots[["srna_ldna"]][["ggpairs_dif_l2fc"]][[framename]] <- p

            pf <- frame %>%
                mutate(ldna_cond2mcond1 = !!sym(paste0("ldna_", condition2)) - !!sym(paste0("ldna_", condition1))) %>%
                dplyr::select(ldna_cond2mcond1, paste0("lrna_", contrast_log2FoldChange))
            p <- ggpairs(pf) + ggtitle(framename)
            mysave(sprintf("%s/%s/ggpairs_dif_l2fc_%s.png", outputdir, "lrna_ldna", framename), 4, 4)
            plots[["lrna_ldna"]][["ggpairs_dif_l2fc"]][[framename]] <- p
        },
        error = function(e) {
            print(e)
        }
    )
}

dir.create(dirname(outputs$plots), showWarnings = FALSE, recursive = TRUE)
save(plots, file = outputs$plots)