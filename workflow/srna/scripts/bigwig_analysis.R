module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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
            txdbrefseq = "aref/A.REF_annotations/refseq.sqlite"
        ), env = globalenv())
        assign("inputs", list(
            bwF = sprintf("srna/outs/%s/star_output/%s.F.bw", conf$samples, conf$samples),
            bwR = sprintf("srna/outs/%s/star_output/%s.R.bw", conf$samples, conf$samples)
        ), env = globalenv())
        assign("outputs", list(
            "environment" = "srna/results/agg/bigwig_plots/bigwigplots_environment.RData"
        ), env = globalenv())
    }
)

outputdir <- params$outputdir
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


paths_bwF <- inputs$bwF
paths_bwR <- inputs$bwR

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
                eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                eligible_facet_modifiers <- c("ALL")
                eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
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
                                filter(str_detect(!!sym(filter_var), "Intact|FL$|^LTR"))
                        } else {
                            elements_of_interest <- rmann %>%
                                filter(!!sym(ontology) == group)
                        }

                        if (length(rownames(elements_of_interest)) > 10000) {
                            elements_of_interest <- elements_of_interest %>% sample_n(10000)
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
                                left_join(norm_by_aligned_reads) %>%
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
                                left_join(norm_by_aligned_reads) %>%
                                mutate(value = value / scale_factor) %>%
                                mutate(smoothed_value = smoothed_value / scale_factor)
                            pfStranded1 <- pfStranded %>%
                                ungroup() %>%
                                group_by(x, strand, condition) %>%
                                mutate(condition_value = mean(value, na.rm = TRUE)) %>%
                                mutate(smoothed_condition_value = zoo::rollmean(condition_value, 5, fill = NA, align = "left")) %>%
                                ungroup()

                            element_anatomy <- read_delim("aref/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")
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
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }

                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_by_sample.pdf", signal_group), 8, 6)

                            p1 <- pf %>%
                                group_by(x) %>%
                                summarise(value = mean(value, na.rm = TRUE)) %>%
                                ungroup() %>%
                                ggplot(aes(x = x, y = value)) +
                                geom_line() +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_all.pdf", signal_group))

                            p1 <- pf1 %>%
                                ggplot(aes(x = x, y = condition_value, color = condition)) +
                                geom_line() +
                                mtclosed +
                                scale_conditions +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_by_condition.pdf", signal_group), 6, 4)


                            p1 <- pfStranded %>%
                                ggplot(aes(x = x, y = smoothed_value, color = sample_name)) +
                                geom_line() +
                                facet_wrap(~strand, nrow = 2) +
                                mtclosed +
                                scale_samples +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_by_sample_stranded.pdf", signal_group), 8, 6)

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
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_all_stranded.pdf", signal_group), 6, 6)

                            p1 <- pfStranded1 %>%
                                ggplot(aes(x = x, y = condition_value, color = condition)) +
                                geom_line() +
                                facet_wrap(~strand, nrow = 2) +
                                mtclosed +
                                scale_conditions +
                                labs(x = "Position (bp)", y = "Read Density", caption = "Multi")
                            if (elements_of_interest %$% rte_family %>% unique() == "L1") {
                                p2temp <- p2 + coord_cartesian(xlim = layer_scales(p1)$x$range$range)
                                p <- p2temp / p1 + plot_layout(heights = c(0.2, 1))
                            } else {
                                p <- p1
                            }
                            mysaveandstore(sprintf("srna/results/agg/bigwig_plots/rte/%s_profile_by_condition_stranded.pdf", signal_group), 6, 6)
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

grs_list <- list()
for (sample in sample_table$sample_name) {
    bwF <- import(grep(sprintf("/%s/", sample), paths_bwF, value = TRUE))
    strand(bwF) <- "+"
    bwR <- import(grep(sprintf("/%s/", sample), inputs$bwR, value = TRUE))
    strand(bwR) <- "-"
    grstemp <- c(bwF, bwR)
    grstemp <- grstemp[!grepl("^NI", seqnames(grstemp))]
    seqlevels(grstemp, pruning.mode = "coarse") <- seqlevelsInUse(grstemp)
    mcols(grstemp)$sample_name <- sample
    mcols(grstemp)$score <- mcols(grstemp)$score / norm_by_aligned_reads$scale_factor[norm_by_aligned_reads$sample_name == sample]
    grs_list[[sample]] <- grstemp
}
grs <- Reduce(c, grs_list)
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

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_faceted.pdf"), 8, 3.8)

pf <- region_annot_rm %>%
    group_by(sample_name, loc_integrative) %>%
    summarise(score_sum = sum(score)) %>%
    left_join(sample_table) %>%
    ungroup()

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

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_rm.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_faceted_rm.pdf"), 8, 3.8)

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

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_rmpriority.pdf"), 6, 3.8)

p <- pf %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    ggbarplot(x = "condition", y = "score_sum", facet.by = "loc_integrative", fill = "condition", scales = "free_y", add = c("mean_se", "dotplot")) +
    geom_pwc(
        method = "t_test", label = "p.adj.signif",
        ref.group = conf$levels[1],
        p.adjust.method = "fdr", hide.ns = TRUE
    ) + scale_conditions + mtclosed + anchorbar + labs(x = "", y = "Normalized Read Count")

mysaveandstore(sprintf("srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_faceted_rmpriority.pdf"), 8, 3.8)

if (conf$store_env_as_rds == "yes") {
    save.image(file = outputs$environment)
} else {
    x <- tibble(Env_file = "Opted not to store environment. If this is not desired, change 'store_plots_as_rds' to 'yes' in the relevant config file and rerun this rule.")
    write_tsv(x, file = outputs$environment)
}

# figures: modify plot compositions at will!
# load(outputs$environment)
tryCatch(
    {
        library(patchwork)
        p1 <- mysaveandstoreplots[["srna/results/agg/bigwig_plots/genomic_context/gene_oriented_signal_faceted.pdf"]]
        p2 <- mysaveandstoreplots[["srna/results/agg/bigwig_plots/rte/L1HS_rte_length_req_profile_by_sample_stranded.pdf"]]
        names(mysaveandstoreplots)
        ptch <- p1 + p2 + plot_layout(ncol = 2, guides = "collect")
        mysaveandstore(pl = ptch, fn = "srna/results/agg/bigwig_plots/composition.pdf", w = 15, h = 7)
    },
    error = function(e) {

    }
)
