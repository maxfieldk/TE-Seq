source("workflow/scripts/defaults.R")
module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

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
library(plotly)
library(DT)
library(ggExtra)
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
            "inputdir" = "srna/results/agg/deseq_telescope",
            "outputdir" = "srna/results/agg/repeatanalysis_telescope",
            "tecounttypes" = c("telescope_multi"),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
            bwF = sprintf("srna/outs/%s/star_output/%s.F.bw", conf$samples, conf$samples),
            bwR = sprintf("srna/outs/%s/star_output/%s.R.bw", conf$samples, conf$samples),
            txdbrefseq = "aref/A.REF_annotations/refseq.sqlite"

        ), env = globalenv())
        assign("outputs", list(
            "plots" = "srna/results/agg/bigwig_plots/plots.rds"
        ), env = globalenv())
    }
)

outputdir <- params$outputdir
contrasts <- conf$contrasts

r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
rmann <-r_annotation_fragmentsjoined %>%
    left_join(r_repeatmasker_annotation)

refseq <- import(conf$annotation_genes)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
temptxs <- coding_transcripts[!is.na(mcols(coding_transcripts)$tag)]
refseq_select <- temptxs[mcols(temptxs)$tag == "RefSeq Select"]
noncoding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NR", mcols(refseq)$transcript_id))]
transcripts <- c(coding_transcripts, noncoding_transcripts)
table(mcols(coding_transcripts)$tag)

### ONTOLOGY DEFINITION
# {
#     annot_colnames <- colnames(r_repeatmasker_annotation)
#     annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
#     ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
#     small_ontologies <- ontologies[grepl("subfamily", ontologies)]

#     big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
#     big_ontology_groups <- c()
#     for (ontology in big_ontologies) {
#         big_ontology_groups <- c(big_ontology_groups, resultsdf %>%
#             pull(!!sym(ontology)) %>%
#             unique())
#     }
#     big_ontology_groups <- big_ontology_groups %>% unique()

#     modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
#     region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
#     element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
# }


# windowsR <- coding_transcripts %>% as.data.frame() %>% tibble() %>% filter(gene_id %in% c("BRD1", "CPT1B", "ACE2", "STAU1")) %>% GRanges()
# windowsF <- windowsR
# strand(windowsF) <- "+"
for (group in c("L1HS")) {
windowsF <- rmann %>% filter(rte_subfamily == group) %>% filter(req_integrative != "Truncated") %>% filter(strand == "+") %>% GRanges()
windowsR <- rmann %>% filter(rte_subfamily == group) %>% filter(req_integrative != "Truncated") %>% filter(strand == "-") %>% GRanges()

{
smlF = genomation::ScoreMatrixList(targets = inputs$bwF, windows = windowsF, strand.aware = TRUE, bin.num = 100)
smlR = genomation::ScoreMatrixList(targets = inputs$bwR, windows = windowsR, strand.aware = TRUE, bin.num = 100)
mF = plotMeta(smlF, plot = FALSE) %>% t()
mR = plotMeta(smlR, plot = FALSE) %>% t()
m <- mF + mR
df = m %>% as.data.frame() %>% tibble()
dfF = mF %>% as.data.frame() %>% tibble()
dfR = mR %>% as.data.frame() %>% tibble()
colnames(df) <- conf$samples
colnames(dfF) <- conf$samples
colnames(dfR) <- conf$samples

df$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
dfF$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
dfR$x <- seq(0, width(windowsF) %>% mean(), length.out = nbin)
dfF$strand = "+"
dfR$strand = "-"
dfStranded = rbind(dfF, dfR)
pf <- df%>% pivot_longer(cols = c(-x), names_to = "sample_name", values_to = "value") %>% 
    left_join(sample_table) %>% group_by(sample_name) %>%
    mutate(smoothed_value = zoo::rollmean(value, 5, fill = NA, align = "left"))
pf1 <- pf %>% ungroup() %>% group_by(x, condition) %>% mutate(condition_value = mean(value, na.rm = TRUE)) %>% 
    mutate(smoothed_condition_value = zoo::rollmean(condition_value, 5, fill = NA, align = "left")) %>% ungroup()

pfStranded <- dfStranded %>% pivot_longer(cols = c(-x, -strand), names_to = "sample_name", values_to = "value") %>% 
    left_join(sample_table) %>% group_by(sample_name, strand) %>%
    mutate(smoothed_value = zoo::rollmean(value, 5, fill = NA, align = "left"))
pfStranded1 <- pfStranded %>% ungroup() %>% group_by(x, strand, condition) %>% mutate(condition_value = mean(value, na.rm = TRUE)) %>%
    mutate(smoothed_condition_value = zoo::rollmean(condition_value, 5, fill = NA, align = "left")) %>% ungroup()

element_anatomy <- read_delim("aref/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")
representative_element <- element_anatomy %>% filter(gene_id == element_anatomy$gene_id[1]) %>% filter(!(feature %in% c("EN", "RT")))
p2 <- representative_element %>% 
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
    geom_text(aes(x = (start + end) / 2, y = 1.5, label = feature)) +
    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
    ggtitle(group) +
    scale_fill_paletteer_d("dutchmasters::milkmaid") + mtclosed + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0, 0.4)) + theme(legend.position = "none")

p1 <- pf %>%
    ggplot(aes(x = x, y = smoothed_value, color = sample_name)) + geom_line() + mtclosed + scale_samples+ labs(x = "Position (bp)", y = "Read Density (cpm)", caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_by_sample.png", 8,4, group))

p1 <- pf %>% group_by(x) %>% summarise(value = mean(value, na.rm = TRUE)) %>% ungroup() %>%
    ggplot(aes(x = x, y = value)) + geom_line() + mtclosed + scale_samples + labs(x = "Position (bp)", y = "Read Density (cpm)", caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_all.png", group))

p1 <- pf1%>%
    ggplot(aes(x = x, y = condition_value, color = condition)) + geom_line() + mtclosed + scale_conditions + labs(x = "Position (bp)", y = "Read Density (cpm)", title = group, caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_by_condition.png", group))


p1 <- pfStranded %>%
    ggplot(aes(x = x, y = smoothed_value, color = sample_name)) + geom_line() + facet_wrap(~strand, nrow = 2) + mtclosed + scale_samples +
    labs(x = "Position (bp)", y = "Read Density (cpm)", caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_by_sample_stranded.png", 8,4, group))

p1 <- pfStranded %>% group_by(x, strand) %>% summarise(value = mean(value, na.rm = TRUE)) %>% ungroup() %>%
    ggplot(aes(x = x, y = value)) + geom_line() + facet_wrap(~strand, nrow = 2) + mtclosed + scale_samples + labs(x = "Position (bp)", y = "Read Density (cpm)", caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_all_stranded.png", group))

p1 <- pfStranded1 %>%
    ggplot(aes(x = x, y = condition_value, color = condition)) + geom_line() + facet_wrap(~strand, nrow = 2) + mtclosed + scale_conditions + labs(x = "Position (bp)", y = "Read Density (cpm)", title = group, caption = "Multi")
#p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
p <- p1 
mysaveandstore(sprintf("srna/results/agg/bigwig_plots/%s_profile_by_condition_stranded.png", group))
}
}

grs_list <- list()
for (sample in sample_table$sample_name) {
    bwF <- import(grep(sprintf("/%s/", sample), inputs$bwF, value = TRUE))
    strand(bwF) <- "+"
    bwR <- import(grep(sprintf("/%s/", sample), inputs$bwR, value = TRUE))
    strand(bwR) <- "-"
    grstemp <- c(bwF, bwR)
    mcols(grstemp)$sample_name <- sample
    grs_list[[sample]] <- grstemp
}
grs <- Reduce(c, grs_list)



txdb <- loadDb(inputs$txdbrefseq)
introns <- intronsByTranscript(txdb, use.names = TRUE)
introns <- introns[grepl("^NM", names(introns))]
introns <- introns[names(introns) %in% mcols(refseq_select)$transcript_id]
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
exons <- exons[grepl("^NM", names(exons))]
exons <- exons[names(exons) %in% mcols(refseq_select)$transcript_id]
fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
fiveUTRs <- fiveUTRs[grepl("^NM", names(fiveUTRs))]
fiveUTRs <- fiveUTRs[names(fiveUTRs) %in% mcols(refseq_select)$transcript_id]
threeUTRs <- threeUTRsByTranscript(txdb, use.names = TRUE)
threeUTRs <- threeUTRs[grepl("^NM", names(threeUTRs))]
threeUTRs <- threeUTRs[names(threeUTRs) %in% mcols(refseq_select)$transcript_id]

getannotation <- function(to_be_annotated, regions_of_interest, property, name_in, name_out) {
    inregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = FALSE, ignore.strand = FALSE)
    mcols(inregions)[[property]] <- name_in
    outregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = TRUE, ignore.strand = FALSE)
    mcols(outregions)[[property]] <- name_out
    merged <- c(inregions, outregions)
    return(merged)
}

exonic_annot <- getannotation(grs, exons, "exonic", "Exonic", "NonExonic")%>% as.data.frame() %>% tibble()
intron_annot <- getannotation(grs, introns, "intronic", "Intronic", "NonIntronic")%>% as.data.frame() %>% tibble()

region_annot <- left_join(exonic_annot, intron_annot)

region_annot1 <- region_annot %>% mutate(loc_integrative = case_when(
    exonic == "Exonic" ~ "Exonic",
    intronic == "Intronic" ~ "Intronic",
    TRUE ~ "NonGenic"
))

pf <- region_annot1 %>% group_by(sample_name, loc_integrative) %>% summarise(score_sum = sum(score)) %>% left_join(sample_table) %>% ungroup()

barframe <- pf %>% group_by(condition, loc_integrative) %>% summarise(score_condition_mean = mean(score_sum)) 

p<- pf %>% ggplot(aes(x = loc_integrative, y = score_sum, fill = sample_name)) + geom_col(position = position_dodge(preserve = "single")) + scale_samples_unique + mtopen + anchorbar
mysaveandstore("srna/results/agg/bigwig_plots/gene_oriented_signal.png", 12,5)

p<- pf %>% ggplot(aes(x = sample_name, y = score_sum, fill = sample_name)) + geom_col(position = position_dodge(preserve = "single")) +
    facet_wrap(~loc_integrative, scale = "free_y") + scale_samples_unique + mtclosed + anchorbar + theme(axis.text.x = element_blank())
mysaveandstore("srna/results/agg/bigwig_plots/gene_oriented_signal_faceted.png", 12,5)

save(mysaveandstoreplots, file = outputs$plots)
