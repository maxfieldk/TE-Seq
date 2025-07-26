if (interactive()) {
    module_name <<- "srna"
} else {
    module_name <<- snakemake@params$module_name
}
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)

set.seed(123)

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
# library("ggVennDiagram")
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
library(rtracklayer)
library(ComplexUpset)
library(patchwork)
library(scales)
library(ggbeeswarm)


# tryCatch(
#     {
#         params <- snakemake@params
#         inputs <- snakemake@input
#         outputs <- snakemake@output
#         print("sourced snake variables")
#     },
#     error = function(e) {
#         print("not sourced snake variables")
#         print(module_name)
#         if (module_name == "srna") {
#             assign("params", list(
#                 "inputdir" = "srna/results/agg/deseq",
#                 "outputdir" = "srna/results/agg/repeatanalysis",
#                 "counttype" = "telescope_multi",
#                 "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
#                 "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
#             ), env = globalenv())
#             assign("inputs", list(
#                 "resultsdf" = "srna/results/agg/deseq/resultsdf.tsv",
#                 "tetranscripts_resultsdf" = "srna/results/agg/deseq_tetranscripts/resultsdf.tsv"
#             ), env = globalenv())
#             assign("outputs", list(
#                 "environment" = "srna/results/agg/repeatanalysis/telescope_multi/repeatanalysisplots_environment.RData"
#             ), env = globalenv())
#         } else if (module_name == "lrna") {
#             assign("params", list(
#                 "inputdir" = "lrna/results/agg/deseq",
#                 "outputdir" = "lrna/results/agg/repeatanalysis/relaxed",
#                 "counttype" = "relaxed",
#                 "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
#                 "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
#             ), env = globalenv())
#             assign("inputs", list(
#                 "resultsdf" = "lrna/results/agg/deseq/resultsdf.tsv",
#                 "tetranscripts_resultsdf" = "lrna/results/agg/tetranscripts/resultsdf.tsv"
#             ), env = globalenv())
#             assign("outputs", list(
#                 "environment" = "lrna/results/agg/repeatanalysis/relaxed/repeatanalysisplots_environment.RData"
#             ), env = globalenv())
#         } else {
#             print("varible assignment failed")
#         }
#     }
# )



countspaths <- list.files(
    path = "lrna/intermediates",
    pattern = "\\.counts\\.txt$", # match files ending in .counts.txt
    recursive = TRUE,
    full.names = TRUE
)
relaxed_paths <- countspaths[grepl("relaxed/", countspaths)]
stringent_paths <- countspaths[grepl("stringent/", countspaths)]

names(relaxed_paths) <- basename(dirname(dirname(dirname(relaxed_paths))))
names(stringent_paths) <- basename(dirname(dirname(dirname(stringent_paths))))



relaxed_list <- list()
for (path in relaxed_paths) {
    sample_name <- basename(dirname(dirname(dirname(path))))
    df <- read_delim(path, delim = "\t")
    colnames(df) <- c("gene_id", sample_name)
    relaxed_list[[path]] <- df
}

relaxed <- purrr::reduce(relaxed_list, full_join, by = "gene_id")



rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE
)
resdf <- relaxed %>% full_join(rmann)

l1hsdf <- resdf %>% filter(rte_subfamily == "L1HS")

l1hsdf %>% arrange(-ORF1_1)
l1hsdf %>% arrange(-ORF2_1)

l1hsdf %>%
    arrange(-ORF2_1) %>%
    filter(rte_length_req == "FL") %>%
    filter(genic_loc == "Intergenic") %>%
    pw()

l1hsdf %>%
    group_by(rte_length_req, intactness_req) %>%
    summarise(morf1 = mean(ORF1_1), morf2 = mean(ORF2_1), mtot = mean(TOT_1))
l1hsdf %>% pw()

l1hsdf %>%
    filter(rte_length_req == "FL") %>%
    group_by(grepl("1017", orf_passes_mutation_threshold)) %>%
    summarise(morf1 = mean(ORF1_1), morf2 = mean(ORF2_1), mtot = mean(TOT_1), n = n())

l1hsdf %>%
    filter(rte_length_req == "FL") %>%
    group_by(grepl("3828", orf_passes_mutation_threshold)) %>%
    summarise(morf1 = mean(ORF1_1), morf2 = mean(ORF2_1), mtot = mean(TOT_1), n = n())

l1hsdf %>%
    filter(rte_length_req == "FL") %>%
    filter(!grepl("1017", orf_passes_mutation_threshold)) %>%
    group_by(grepl("1017", partial_orf_passes_mutation_threshold)) %>%
    summarise(morf1 = mean(ORF1_1), morf2 = mean(ORF2_1), mtot = mean(TOT_1), n = n())

l1hsdf %>%
    filter(rte_length_req == "FL") %>%
    filter(!grepl("3828", orf_passes_mutation_threshold)) %>%
    group_by(grepl("3828", partial_orf_passes_mutation_threshold)) %>%
    summarise(morf1 = mean(ORF1_1), morf2 = mean(ORF2_1), mtot = mean(TOT_1), n = n())

partial_orf_passes_mutation_threshold

l1hsdfintact <- l1hsdf %>% filter(intactness_req == "Intact")

library(magrittr)

totnormsum <- l1hsdf %$% TOT_1 %>%
    sum() %>%
    `/`(27963779)

l1hsdf %$% TOT_1 %>%
    sum() %>%
    `/`(27963779) %>%
    `/`(totnormsum)
l1hsdf %$% ORF1_1 %>%
    sum() %>%
    `/`(27632897) %>%
    `/`(totnormsum)
l1hsdf %$% ORF2_1 %>%
    sum() %>%
    `/`(6359779) %>%
    `/`(totnormsum)
l1hsdf %$% GFP_1 %>%
    sum() %>%
    `/`(3066184) %>%
    `/`(totnormsum)


totnormsum <- l1hsdfintact %$% TOT_1 %>%
    sum() %>%
    `/`(27963779)
l1hsdfintact %$% TOT_1 %>%
    sum() %>%
    `/`(27963779) %>%
    `/`(totnormsum)
l1hsdfintact %$% ORF1_1 %>%
    sum() %>%
    `/`(27632897) %>%
    `/`(totnormsum)
l1hsdfintact %$% ORF2_1 %>%
    sum() %>%
    `/`(6359779) %>%
    `/`(totnormsum)
l1hsdfintact %$% GFP_1 %>%
    sum() %>%
    `/`(3066184) %>%
    `/`(totnormsum)


short_read <- read_delim("srna/results/agg/deseq/resultsdf.tsv", delim = "\t") %>% left_join(rmann)

sl1hs <- short_read %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(counttype == "telescope_multi")



# rank cor
idforcor <- sl1hs %>% filter(refstatus == "Ref") %$% gene_id
ranks1 <- sl1hs %>%
    filter(gene_id %in% idforcor) %>%
    mutate(rank1 = rank(-N21_ORF1p_500_b1_l1)) %>%
    select(gene_id, rank1)

ranks2 <- l1hsdf %>%
    filter(gene_id %in% idforcor) %>%
    mutate(rank2 = rank(-ORF1_1)) %>%
    select(gene_id, rank2)
ranked <- inner_join(ranks1, ranks2, by = "gene_id")
cor.test(ranked$rank1, ranked$rank2, method = "spearman")


idforcor <- sl1hs %>% filter(refstatus == "Ref") %$% gene_id
ranks1 <- sl1hs %>%
    filter(gene_id %in% idforcor) %>%
    mutate(rank1 = rank(-N21_total_500_b1_l1)) %>%
    select(gene_id, rank1)

ranks2 <- l1hsdf %>%
    filter(gene_id %in% idforcor) %>%
    mutate(rank2 = rank(-TOT_1)) %>%
    select(gene_id, rank2)
ranked <- inner_join(ranks1, ranks2, by = "gene_id")
cor.test(ranked$rank1, ranked$rank2, method = "spearman")



refseq <- import(confALL$aref$ref_refseq_gtf)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
codinggene <- mcols(coding_transcripts) %>%
    as.data.frame() %$% gene_id %>%
    unique()

short_read %>%
    filter(gene_id %in% codinggene) %>%
    dplyr::select(gene_id, N21_total_500_b1_l1)
resdf %>%
    mutate(gene_id = gsub("gene-", "", gene_id)) %>%
    filter(gene_id %in% codinggene) %>%
    dplyr::select(gene_id, TOT_1)


t2 <- short_read %>%
    filter(gene_id %in% codinggene) %>%
    dplyr::select(gene_id, N21_total_500_b1_l1) %>%
    mutate(rank2 = rank(-N21_total_500_b1_l1)) %>%
    select(gene_id, rank2)
t1 <- resdf %>%
    mutate(gene_id = gsub("gene-", "", gene_id)) %>%
    filter(gene_id %in% codinggene) %>%
    mutate(rank1 = rank(-TOT_1)) %>%
    select(gene_id, rank1)
ranked <- inner_join(t1, t2, by = "gene_id")
cor.test(ranked$rank1, ranked$rank2, method = "spearman")
