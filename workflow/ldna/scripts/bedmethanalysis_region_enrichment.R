module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(msigdbr)
library(regioneR)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            filtered_tldr = "aref/A.REF_tldr/A.REF.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/A.REF.fa",
            blast_njs = "aref/A.REF.njs",
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/A.REF_Analysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/A.REF_Analysis/transduction_df.csv"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "aref/A.REF"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "aref/A.REF"
        ), env = globalenv())
    }
)



dmrs <- read_delim(inputs$dmrs, delim = "\t", col_names = TRUE)

dmrsgr <- GRanges(dmrs)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)


library(BiocParallel)
register(MulticoreParam(12))

fa <- Rsamtools::FaFile("aref/A.REF.fa")
chr_lengths <- seqinfo(fa) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tibble() %>%
    dplyr::rename(seqnames = rowname, end = seqlengths) %>%
    mutate(start = 1) %>%
    dplyr::select(seqnames, start, end) %>%
    GRanges()


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
    width(chromHMMgr) %>% sum()
    annotations_of_interest <- list(chromHMM = chromHMMgr, cCREs = ccresgr)
    # note I should be testing for enrichment of these sets too
    permTestResults <- list()
    rm(resultsframe)
    rm(ptresults)
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
        permTestResults[[annotation_of_interest]] <- list()
        for (annotation_subtype in mcols(annot)$name %>% unique()) {
            soigr <- annot[mcols(annot)$name == annotation_subtype]
            for (direction in c("Hypo", "Hyper", "Both")) {
                if (direction == "Both") {
                    dmrs_to_test <- sample(dmrsgr, 5000)
                } else {
                    dmrs_to_test <- sample(dmrsgr[mcols(dmrsgr)$direction == direction], 5000)
                }
                pt <- permTest(
                    A = soigr, B = sample(dmrsgr, 5000), genome = chr_lengths, randomize.function = randomizeRegions,
                    evaluate.function = numOverlaps, ntimes = 30
                )
                row <- tibble(direction = direction, annotation_of_interest = annotation_of_interest, annotation_subtype = annotation_subtype, pvalue = pt$numOverlaps$pval, zscore = pt$numOverlaps$zscore)
                if (!exists("resultsframe")) {
                    resultsframe <- tibble(direction = character(), annotation_of_interest = character(), annotation_subtype = character(), pvalue = numeric(), zscore = numeric())
                }
                resultsframe <- resultsframe %>% add_row(row)
                if (!exists("ptresults")) {
                    ptresults <- list()
                }
                ptresults[[direction]][[annotation_subtype]] <- pt
            }
        }
    }
}
for (direction in c("Hypo", "Hyper", "Both")) {
    for (group in resultsframe %$% annotation_of_interest %>% unique()) {
        p <- resultsframe %>%
            filter(direction == {{ direction }}) %>%
            filter(annotation_of_interest == {{ group }}) %>%
            ggplot(aes(x = zscore, fill = pvalue, y = annotation_subtype)) +
            geom_col() +
            mtclosedgridv +
            labs(y = "", title = "DMR Enrichment")
        mysaveandstore(str_glue("ldna/results/plots/{group}/{group}_{direction}_zscore.pdf"))
    }
}
for (group in resultsframe %$% annotation_of_interest %>% unique()) {
    p <- resultsframe %>%
        filter(direction != "Both") %>%
        filter(annotation_of_interest == {{ group }}) %>%
        ggplot(aes(x = zscore, fill = direction, y = annotation_subtype)) +
        geom_col(position = "dodge") +
        mtclosedgridv +
        labs(y = "", title = "DMR Enrichment")
    mysaveandstore(str_glue("ldna/results/plots/{group}/{group}_bothdirections_zscore.pdf"))
}
