module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggVennDiagram")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)
library(patchwork)
library(zoo)
library("GenomicRanges")
library(rtracklayer)
library(genomation)


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
            "inputdir" = "srna/results/agg/deseq",
            "outputdir" = "srna/results/agg/repeatanalysis",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            txdbrefseq = "aref/default/A.REF_annotations/refseq.sqlite"
        ), env = globalenv())
        assign("inputs", list(
            bwF = sprintf("srna/outs/%s/star_output/%s.unique.F.bw", conf$samples, conf$samples),
            bwR = sprintf("srna/outs/%s/star_output/%s.unique.R.bw", conf$samples, conf$samples)
        ), env = globalenv())
        assign("outputs", list(
            outfile = "srna/results/agg/bigwig_plots/unique/bigwigplots.txt"
        ), env = globalenv())
    }
)

paths_bwF <- inputs$bwF
paths_bwR <- inputs$bwR


outputdir <- dirname(outputs$outfile)
contrasts <- conf$contrasts

r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
rmann <- r_annotation_fragmentsjoined %>%
    left_join(r_repeatmasker_annotation)

refseq <- import(conf$annotation_genes)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
temptxs <- coding_transcripts[!is.na(mcols(coding_transcripts)$tag)]
refseq_select <- temptxs[mcols(temptxs)$tag == "RefSeq Select"]
noncoding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NR", mcols(refseq)$transcript_id))]
transcripts <- c(coding_transcripts, noncoding_transcripts)
table(mcols(coding_transcripts)$tag)

tryCatch(
    {
        assign("norm_by_aligned_reads", read_delim("srna/qc/multiqc/multiqc_data/samtools-stats-dp.txt") %>%
            filter(Sample %in% sample_table$sample_name) %>%
            mutate(meanmandp = mean(`Mapped &amp; paired`)) %>%
            mutate(scale_factor = `Mapped &amp; paired` / meanmandp) %>%
            dplyr::select(Sample, scale_factor, `Mapped &amp; paired`, meanmandp) %>%
            dplyr::rename(sample_name = Sample), env = globalenv())
    },
    error = function(e) {
        assign("norm_by_aligned_reads", read_delim("srna/qc/multiqc/multiqc_data/multiqc_samtools_stats.txt") %>%
            filter(Sample %in% sample_table$sample_name) %>%
            mutate(meanmandp = mean(reads_mapped_and_paired)) %>%
            mutate(scale_factor = reads_mapped_and_paired / meanmandp) %>%
            dplyr::select(Sample, scale_factor, reads_mapped_and_paired, meanmandp) %>%
            dplyr::rename(sample_name = Sample), env = globalenv())
    }
)



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




grs_list <- list()
grs_total_score <- list()
for (sample in conf$samples) {
    bwF <- import(grep(sprintf("/%s/", sample), paths_bwF, value = TRUE))
    strand(bwF) <- "+"
    bwR <- import(grep(sprintf("/%s/", sample), inputs$bwR, value = TRUE))
    strand(bwR) <- "-"
    grstemp <- c(bwF, bwR)
    grstemp <- grstemp[!grepl("^NI_", seqnames(grstemp))]
    seqlevels(grstemp, pruning.mode = "coarse") <- seqlevelsInUse(grstemp)
    mcols(grstemp)$sample_name <- sample
    score <- grstemp$score %>% sum()
    grs_total_score[[sample]] <- score
    mcols(grstemp)$score <- mcols(grstemp)$score / score
    # mcols(grstemp)$score <- mcols(grstemp)$score / norm_by_aligned_reads$scale_factor[norm_by_aligned_reads$sample_name == sample]
    grs_list[[sample]] <- grstemp
}
grs <- Reduce(c, grs_list)


# total_score_df <- tibble(sample_name = names(grs_total_score), score = unlist(grs_total_score)) %>% mutate(scale_factor = score / 1000000)
# total_score_df$score / max(total_score_df$score)
# norm_by_aligned_reads$`Mapped &amp; paired` / max(norm_by_aligned_reads$`Mapped &amp; paired`)

txdb <- loadDb(params$txdbrefseq)
rmgrs <- rmann %>% GRanges()
introns <- intronsByTranscript(txdb, use.names = TRUE)
# introns <- introns[grepl("^NM", names(introns))]
introns <- introns[names(introns) %in% mcols(refseq_select)$transcript_id]
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
# exons <- exons[grepl("^NM", names(exons))]
exons <- exons[names(exons) %in% mcols(refseq_select)$transcript_id]
fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
# fiveUTRs <- fiveUTRs[grepl("^NM", names(fiveUTRs))]
fiveUTRs <- fiveUTRs[names(fiveUTRs) %in% mcols(refseq_select)$transcript_id]
threeUTRs <- threeUTRsByTranscript(txdb, use.names = TRUE)
# threeUTRs <- threeUTRs[grepl("^NM", names(threeUTRs))]
threeUTRs <- threeUTRs[names(threeUTRs) %in% mcols(refseq_select)$transcript_id]

getannotation <- function(to_be_annotated, regions_of_interest, property, name_in, name_out) {
    inregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = FALSE, ignore.strand = FALSE)
    mcols(inregions)[[property]] <- name_in
    outregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = TRUE, ignore.strand = FALSE)
    mcols(outregions)[[property]] <- name_out
    merged <- c(inregions, outregions)
    return(merged)
}

exonic_annot <- getannotation(grs, exons, "exonic", "Exonic", "NonExonic") %>%
    as.data.frame() %>%
    tibble()
intron_annot <- getannotation(grs, introns, "intronic", "Intronic", "NonIntronic") %>%
    as.data.frame() %>%
    tibble()
region_annot <- left_join(exonic_annot, intron_annot)
rm(exonic_annot)
rm(intron_annot)
rm_annot <- getannotation(grs, rmgrs, "repetitive", "Repetitive", "NonRepetitive") %>%
    as.data.frame() %>%
    tibble()
rm(grs)
region_annot1 <- region_annot %>% mutate(loc_integrative = case_when(
    exonic == "Exonic" ~ "Exonic",
    intronic == "Intronic" ~ "Intronic",
    TRUE ~ "Intergenic"
))
region_annot_rm <- left_join(region_annot, rm_annot)
region_annot_rm <- region_annot_rm %>% mutate(loc_integrative = case_when(
    exonic == "Exonic" ~ "Exonic",
    intronic == "Intronic" ~ "Intronic",
    repetitive == "Repetitive" ~ "Repetitive",
    TRUE ~ "Intergenic"
))

pf <- region_annot1 %>%
    group_by(sample_name, loc_integrative) %>%
    summarise(score_sum = sum(score)) %>%
    left_join(sample_table) %>%
    ungroup()

RINcol <- colnames(pf)[colnames(pf) == "RIN" | colnames(pf) == "batchCon_RIN"]
if (length(RINcol) > 0) {
    RINcoltmp <- RINcol[1]
    p <- pf %>% ggplot(aes(x = score_sum, y = !!sym(RINcoltmp))) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition)) +
        facet_wrap(~loc_integrative, scales = "free_x") +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_RIN_cor.pdf"), 8, 3.8)
}

barframe <- pf %>%
    group_by(condition, loc_integrative) %>%
    summarise(score_condition_mean = mean(score_sum))

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "loc_integrative", y = "score_sum", fill = "condition", scales = "free_y", add = c("mean_se"), position = position_dodge()) +
    stat_pwc(
        aes(group = condition),
        method = "t_test", label = "p.adj.format", ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_faceted.pdf"), 8, 3.8)

pf <- region_annot_rm %>%
    group_by(sample_name, loc_integrative) %>%
    summarise(score_sum = sum(score)) %>%
    left_join(sample_table) %>%
    ungroup()

RINcol <- colnames(pf)[colnames(pf) == "RIN" | colnames(pf) == "batchCon_RIN"]
if (length(RINcol) > 0) {
    RINcoltmp <- RINcol[1]
    p <- pf %>% ggplot(aes(x = score_sum, y = !!sym(RINcoltmp))) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition)) +
        facet_wrap(~loc_integrative, scales = "free_x") +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_RM_RIN_cor.pdf"), 7, 5)

    p <- pf %>%
        mutate(condition = factor(condition, levels = conf$levels)) %>%
        ggboxplot(x = "condition", y = "RIN", color = "condition") +
        geom_point() +
        geom_pwc(aes(group = condition),
            tip.length = 0,
            method = "t.test", label = "{p.adj.format}",
            p.adjust.method = "fdr", p.adjust.by = "panel",
            hide.ns = FALSE
        ) +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/RIN_by_condition.pdf"), 4, 4)
}

barframe <- pf %>%
    group_by(condition, loc_integrative) %>%
    summarise(score_condition_mean = mean(score_sum))

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "loc_integrative", y = "score_sum", fill = "condition", scales = "free_y", add = c("mean_se"), position = position_dodge()) +
    stat_pwc(
        aes(group = condition),
        method = "t_test", label = "p.adj.format", ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_rm.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_faceted_rm.pdf"), 8, 3.8)

region_annot_rm_rmpriority <- region_annot_rm %>% mutate(loc_integrative = case_when(
    exonic == "Exonic" ~ "Exonic",
    repetitive == "Repetitive" ~ "Repetitive",
    intronic == "Intronic" ~ "Intronic",
    TRUE ~ "Intergenic"
))

pf <- region_annot_rm_rmpriority %>%
    group_by(sample_name, loc_integrative) %>%
    summarise(score_sum = sum(score)) %>%
    left_join(sample_table) %>%
    ungroup()

RINcol <- colnames(pf)[colnames(pf) == "RIN" | colnames(pf) == "batchCon_RIN"]
if (length(RINcol) > 0) {
    RINcoltmp <- RINcol[1]
    p <- pf %>% ggplot(aes(x = score_sum, y = !!sym(RINcoltmp))) +
        stat_smooth(method = "lm", formula = y ~ x, color = "green", se = TRUE) +
        ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~~")),
            formula = y ~ x, parse = TRUE, label.y.npc = "top"
        ) +
        geom_point(aes(color = condition)) +
        facet_wrap(~loc_integrative) +
        scale_conditions +
        mtclosed
    mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_rmpriority_RIN_cor.pdf"), 6, 3.8)
}


barframe <- pf %>%
    group_by(condition, loc_integrative) %>%
    summarise(score_condition_mean = mean(score_sum))

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "loc_integrative", y = "score_sum", fill = "condition", scales = "free_y", add = c("mean_se"), position = position_dodge()) +
    stat_pwc(
        aes(group = condition),
        method = "t_test", label = "p.adj.format", ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_rmpriority.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "loc_integrative", y = "score_sum", fill = "sample_name", scales = "free_y", position = position_dodge()) +
    scale_palette + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")
mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_rmpriority_by_sample.pdf"), 6, 3.8)
p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = factor(sample, levels = conf$samples)) %>%
    ggbarplot(x = "sample_name", y = "score_sum", fill = "loc_integrative", scales = "free_y") +
    scale_palette + mtclosed + anchorbar + labs(x = "", y = "Read Fraction") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_rmpriority_by_sample_flip.pdf"), 1 + length(conf$samples) / 2.4, h = 5)


p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(str_glue("{outputdir}/genomic_context/gene_oriented_signal_faceted_rmpriority.pdf"), 8, 3.8)





groups_that_have_been_run <- c()
groups_not_to_run <- c()
for (ontology in c("rte_family", "rte_subfamily_limited")) {
    ontology_groups <- r_repeatmasker_annotation %>%
        pull(!!sym(ontology)) %>%
        unique()
    ontology_groups <- ontology_groups[ontology_groups != "Other"]
    for (group in ontology_groups) {
        if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
            print(group)
            groups_that_have_been_run <- c(groups_that_have_been_run, group)
            # maybe change to !!group
            groupframe <- rmann %>% dplyr::filter(!!sym(ontology) == group)
            eligible_modifiers <- c()
            for (modifier in modifiers) {
                values_present <- rmann %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
            }
            for (modifier in eligible_modifiers) {
                values_present <- rmann %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
                # requiring that elements be ful length / intact - else you get FL and truncated elements of varying lengths merged together - no good
                eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)])
                eligible_facet_modifiers <- c("ALL")
                eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE) %>% dplyr::distinct()
            }
            for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                filter_var <- eligible_modifier_combinations[i, ]$filter_var
                facet_var <- eligible_modifier_combinations[i, ]$facet_var
                print(filter_var)
                print(facet_var)
                tryCatch(
                    {
                        if (filter_var != "ALL") {
                            elements_of_interest <- rmann %>%
                                filter(!!sym(ontology) == group) %>%
                                filter(str_detect(!!sym(filter_var), "^Intact|FL$|^LTR"))
                        } else {
                            elements_of_interest <- rmann %>%
                                filter(!!sym(ontology) == group)
                        }

                        if (length(rownames(elements_of_interest)) > 5000) {
                            elements_of_interest <- elements_of_interest %>% sample_n(5000)
                        }
                        signal_group <- paste0(group, "_", filter_var)
                        windowsF <- elements_of_interest %>%
                            filter(strand == "+") %>%
                            GRanges()
                        windowsR <- elements_of_interest %>%
                            filter(strand == "-") %>%
                            GRanges()

                        {
                            nbin <- 100
                            smlF <- genomation::ScoreMatrixList(targets = paths_bwF, windows = windowsF, strand.aware = TRUE, bin.num = nbin)
                            smlR <- genomation::ScoreMatrixList(targets = paths_bwR, windows = windowsR, strand.aware = TRUE, bin.num = nbin)
                            mF <- genomation::plotMeta(smlF, plot = FALSE) %>% t()
                            mR <- genomation::plotMeta(smlR, plot = FALSE) %>% t()
                            m <- mF + mR
                            df <- m %>%
                                as.data.frame() %>%
                                tibble()
                            dfF <- mF %>%
                                as.data.frame() %>%
                                tibble()
                            dfR <- mR %>%
                                as.data.frame() %>%
                                tibble()
                            colnames(df) <- conf$samples
                            colnames(dfF) <- conf$samples
                            colnames(dfR) <- conf$samples

                            df$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
                            dfF$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
                            dfR$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
                            dfF$strand <- "+"
                            dfR$strand <- "-"
                            dfStranded <- rbind(dfF, dfR)
                            pf <- df %>%
                                pivot_longer(cols = c(-x), names_to = "sample_name", values_to = "value") %>%
                                left_join(sample_table) %>%
                                group_by(sample_name) %>%
                                mutate(smoothed_value = zoo::rollmean(value, 5, fill = NA, align = "left")) %>%
                                left_join(total_score_df) %>%
                                mutate(value = value / scale_factor) %>%
                                mutate(smoothed_value = smoothed_value / scale_factor)
                            pf1 <- pf %>%
                                ungroup() %>%
                                group_by(x, condition) %>%
                                mutate(condition_value = mean(value, na.rm = TRUE)) %>%
                                mutate(smoothed_condition_value = zoo::rollmean(condition_value, 5, fill = NA, align = "left")) %>%
                                ungroup()

                            pfStranded <- dfStranded %>%
                                pivot_longer(cols = c(-x, -strand), names_to = "sample_name", values_to = "value") %>%
                                left_join(sample_table) %>%
                                group_by(sample_name, strand) %>%
                                mutate(smoothed_value = zoo::rollmean(value, 5, fill = NA, align = "left")) %>%
                                left_join(total_score_df) %>%
                                mutate(value = value / scale_factor) %>%
                                mutate(smoothed_value = smoothed_value / scale_factor)
                            pfStranded1 <- pfStranded %>%
                                ungroup() %>%
                                group_by(x, strand, condition) %>%
                                mutate(condition_value = mean(value, na.rm = TRUE)) %>%
                                mutate(smoothed_condition_value = zoo::rollmean(condition_value, 5, fill = NA, align = "left")) %>%
                                ungroup()

                            element_anatomy <- read_delim("aref/default/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")
                            representative_element <- element_anatomy %>%
                                filter(gene_id == element_anatomy$gene_id[1]) %>%
                                filter(!(feature %in% c("EN", "RT")))
                            p2 <- representative_element %>%
                                ggplot() +
                                geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
                                geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
                                geom_text(aes(x = (start + end) / 2, y = 1.5, label = feature)) +
                                ggtitle(signal_group) +
                                scale_fill_paletteer_d("dutchmasters::milkmaid") +
                                mtclosed +
                                theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
                                scale_y_continuous(expand = c(0, 0.4)) +
                                theme(legend.position = "none")
                        }

                        {
                            p1 <- pf %>%
                                ggplot(aes(x = x, y = smoothed_value, color = sample_name)) +
                                geom_line() +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }

                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_by_sample.pdf"), 8, 6)

                            p1 <- pf %>%
                                group_by(x) %>%
                                summarise(value = mean(value, na.rm = TRUE)) %>%
                                ungroup() %>%
                                ggplot(aes(x = x, y = value)) +
                                geom_line() +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_all.pdf"))

                            p1 <- pf1 %>%
                                ggplot(aes(x = x, y = condition_value, color = condition)) +
                                geom_line() +
                                mtclosed +
                                scale_conditions +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_by_condition.pdf"), 6, 4)


                            p1 <- pfStranded %>%
                                ggplot(aes(x = x, y = smoothed_value, color = sample_name)) +
                                geom_line() +
                                facet_wrap(~strand, nrow = 2) +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_by_sample_stranded.pdf"), 8, 6)

                            p1 <- pfStranded %>%
                                group_by(x, strand) %>%
                                summarise(value = mean(value, na.rm = TRUE)) %>%
                                ungroup() %>%
                                ggplot(aes(x = x, y = value)) +
                                geom_line() +
                                facet_wrap(~strand, nrow = 2) +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_all_stranded.pdf"), 6, 6)

                            p1 <- pfStranded1 %>%
                                ggplot(aes(x = x, y = condition_value, color = condition)) +
                                geom_line() +
                                facet_wrap(~strand, nrow = 2) +
                                mtclosed +
                                scale_conditions +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if ((elements_of_interest %$% rte_family %>% unique() == "L1") & ((filter_var == "rte_length_req") | (filter_var == "intactness_req"))) {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(str_glue("{outputdir}/rte/{signal_group}_profile_by_condition_stranded.pdf"), 6, 6)
                        }
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        }
    }
}





#######
tryCatch(
    {
        draw_heatmaps <- function(elements_of_interest, element_set_string) {
            average_element_length <- elements_of_interest$length %>% median()
            # Define windows for sense transcription
            windows_F <- elements_of_interest %>%
                filter(strand == "+") %>%
                GRanges() %>%
                promoters(upstream = average_element_length, downstream = average_element_length * 2)

            windows_R <- elements_of_interest %>%
                filter(strand == "-") %>%
                GRanges() %>%
                promoters(upstream = average_element_length, downstream = average_element_length * 2)

            # Define number of bins
            nbin <- 36
            sml_senseF <- ScoreMatrixList(targets = paths_bwF, windows = windows_F, strand.aware = TRUE, bin.num = nbin)
            sml_antisenseF <- ScoreMatrixList(targets = paths_bwR, windows = windows_F, strand.aware = TRUE, bin.num = nbin)
            sml_senseR <- ScoreMatrixList(targets = paths_bwR, windows = windows_R, strand.aware = TRUE, bin.num = nbin)
            sml_antisenseR <- ScoreMatrixList(targets = paths_bwF, windows = windows_R, strand.aware = TRUE, bin.num = nbin)

            # no capping or scaling - just 99 percentile - often fails due to unresolved error with the winsorize parameter
            tryCatch(
                {
                    sense_matrices <- purrr::map2(
                        purrr::map(sml_senseF, ~ .x@.Data),
                        purrr::map(sml_senseR, ~ .x@.Data),
                        rbind
                    )
                    sm <- as(sense_matrices, "ScoreMatrixList")
                    sm@names <- conf$samples

                    antisense_matrices <- purrr::map2(
                        purrr::map(sml_antisenseF, ~ .x@.Data),
                        purrr::map(sml_antisenseR, ~ .x@.Data),
                        rbind
                    )
                    asm <- as(antisense_matrices, "ScoreMatrixList")
                    asm@names <- conf$samples

                    p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(sm, xcoords = c(-1, 2), winsorize = c(0, 99), common.scale = TRUE, order = TRUE, col = c("white", "blue"))))
                    mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_s.pdf"), w = 2 * length(conf$samples), h = 5)
                    p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(asm, xcoords = c(-1, 2), winsorize = c(0, 99), common.scale = TRUE, order = TRUE, col = c("white", "blue"))))
                    mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_a.pdf"), w = 2 * length(conf$samples), h = 5)
                },
                error = function(e) {

                }
            )

            # cap at 1
            cap_at_1 <- function(mat) {
                mat[mat > 1] <- 1
                return(mat)
            }
            # Apply scaling to ScoreMatrixList
            sml_senseF_capped_1 <- scaleScoreMatrixList(sml_senseF, scalefun = cap_at_1)
            sml_senseR_capped_1 <- scaleScoreMatrixList(sml_senseR, scalefun = cap_at_1)
            sml_antisenseF_capped_1 <- scaleScoreMatrixList(sml_antisenseF, scalefun = cap_at_1)
            sml_antisenseR_capped_1 <- scaleScoreMatrixList(sml_antisenseR, scalefun = cap_at_1)

            sense_matrices_capped_1 <- purrr::map2(
                purrr::map(sml_senseF_capped_1, ~ .x@.Data),
                purrr::map(sml_senseR_capped_1, ~ .x@.Data),
                rbind
            )
            sm_capped_1 <- as(sense_matrices_capped_1, "ScoreMatrixList")
            sm_capped_1@names <- conf$samples

            antisense_matrices_capped_1 <- purrr::map2(
                purrr::map(sml_antisenseF_capped_1, ~ .x@.Data),
                purrr::map(sml_antisenseR_capped_1, ~ .x@.Data),
                rbind
            )
            asm_capped_1 <- as(antisense_matrices_capped_1, "ScoreMatrixList")
            asm_capped_1@names <- conf$samples

            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(sm_capped_1, xcoords = c(-1, 2), grid = TRUE, order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_s_capped_1.pdf"), w = 2 * length(conf$samples), h = 5)
            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(asm_capped_1, xcoords = c(-1, 2), order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_a_capped_1.pdf"), w = 2 * length(conf$samples), h = 5)

            ## capped at 5
            cap_at_5 <- function(mat) {
                mat[mat > 5] <- 5
                return(mat)
            }
            # Apply scaling to ScoreMatrixList
            sml_senseF_capped_5 <- scaleScoreMatrixList(sml_senseF, scalefun = cap_at_5)
            sml_senseR_capped_5 <- scaleScoreMatrixList(sml_senseR, scalefun = cap_at_5)
            sml_antisenseF_capped_5 <- scaleScoreMatrixList(sml_antisenseF, scalefun = cap_at_5)
            sml_antisenseR_capped_5 <- scaleScoreMatrixList(sml_antisenseR, scalefun = cap_at_5)

            sense_matrices_capped_5 <- purrr::map2(
                purrr::map(sml_senseF_capped_5, ~ .x@.Data),
                purrr::map(sml_senseR_capped_5, ~ .x@.Data),
                rbind
            )
            sm_capped_5 <- as(sense_matrices_capped_5, "ScoreMatrixList")
            sm_capped_5@names <- conf$samples

            antisense_matrices_capped_5 <- purrr::map2(
                purrr::map(sml_antisenseF_capped_5, ~ .x@.Data),
                purrr::map(sml_antisenseR_capped_5, ~ .x@.Data),
                rbind
            )
            asm_capped_5 <- as(antisense_matrices_capped_5, "ScoreMatrixList")
            asm_capped_5@names <- conf$samples

            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(sm_capped_5, xcoords = c(-1, 2), grid = TRUE, order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_s_capped_5.pdf"), w = 2 * length(conf$samples), h = 5)
            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(asm_capped_5, xcoords = c(-1, 2), order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_a_capped_5.pdf"), w = 2 * length(conf$samples), h = 5)

            ## capped at 10
            cap_at_10 <- function(mat) {
                mat[mat > 10] <- 10
                return(mat)
            }
            # Apply scaling to ScoreMatrixList
            sml_senseF_capped_10 <- scaleScoreMatrixList(sml_senseF, scalefun = cap_at_10)
            sml_senseR_capped_10 <- scaleScoreMatrixList(sml_senseR, scalefun = cap_at_10)
            sml_antisenseF_capped_10 <- scaleScoreMatrixList(sml_antisenseF, scalefun = cap_at_10)
            sml_antisenseR_capped_10 <- scaleScoreMatrixList(sml_antisenseR, scalefun = cap_at_10)

            sense_matrices_capped_10 <- purrr::map2(
                purrr::map(sml_senseF_capped_10, ~ .x@.Data),
                purrr::map(sml_senseR_capped_10, ~ .x@.Data),
                rbind
            )
            sm_capped_10 <- as(sense_matrices_capped_10, "ScoreMatrixList")
            sm_capped_10@names <- conf$samples

            antisense_matrices_capped_10 <- purrr::map2(
                purrr::map(sml_antisenseF_capped_10, ~ .x@.Data),
                purrr::map(sml_antisenseR_capped_10, ~ .x@.Data),
                rbind
            )
            asm_capped_10 <- as(antisense_matrices_capped_10, "ScoreMatrixList")
            asm_capped_10@names <- conf$samples

            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(sm_capped_10, xcoords = c(-1, 2), grid = TRUE, order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_s_capped_10.pdf"), w = 2 * length(conf$samples), h = 5)
            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(asm_capped_10, xcoords = c(-1, 2), order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_a_capped_10.pdf"), w = 2 * length(conf$samples), h = 5)

            ## SCALED
            sml_senseF_scaled <- scaleScoreMatrixList(sml_senseF, row = TRUE)
            sml_senseR_scaled <- scaleScoreMatrixList(sml_senseR, row = TRUE)
            sml_antisenseF_scaled <- scaleScoreMatrixList(sml_antisenseF, row = TRUE)
            sml_antisenseR_scaled <- scaleScoreMatrixList(sml_antisenseR, row = TRUE)

            sense_matrices_scaled <- purrr::map2(
                purrr::map(sml_senseF_scaled, ~ .x@.Data),
                purrr::map(sml_senseR_scaled, ~ .x@.Data),
                rbind
            )
            sm_scaled <- as(sense_matrices_scaled, "ScoreMatrixList")
            sm_scaled@names <- conf$samples

            antisense_matrices_scaled <- purrr::map2(
                purrr::map(sml_antisenseF_scaled, ~ .x@.Data),
                purrr::map(sml_antisenseR_scaled, ~ .x@.Data),
                rbind
            )
            asm_scaled <- as(antisense_matrices_scaled, "ScoreMatrixList")
            asm_scaled@names <- conf$samples

            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(sm_scaled, xcoords = c(-1, 2), grid = TRUE, order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_s_scaled.pdf"), w = 2 * length(conf$samples), h = 5)
            p <- wrap_elements(grid.grabExpr(genomation::multiHeatMatrix(asm_scaled, xcoords = c(-1, 2), order = TRUE, col = c("white", "blue"))))
            mysaveandstore(pl = p, fn = str_glue("{outputdir}/{element_set_string}_a_scaled.pdf"), w = 2 * length(conf$samples), h = 5)
        }




        elements_of_interest <- rmann %>%
            filter(rte_family == "L1") %>%
            filter(req_integrative %in% c("Yng FL", "Yng Intact")) %>%
            filter(rte_length_req == "FL")
        included <- elements_of_interest %$% rte_subfamily %>%
            unique() %>%
            paste(collapse = "_")
        draw_heatmaps(elements_of_interest, str_glue("YngL1s_{included}"))



        for (subfam in rmann %>%
            filter(rte_subfamily != "Other") %$% rte_subfamily %>%
            unique()) {
            print("drawing heatmaps for")
            print(subfam)
            elements_of_interest <- rmann %>%
                filter(rte_subfamily == subfam) %>%
                filter(rte_length_req == "FL")
            draw_heatmaps(elements_of_interest, subfam)
        }
    },
    error = function(e) {
        print("heatmaps did not work")
    }
)



x <- tibble(OUT = "")
write_tsv(x, file = outputs$outfile)
