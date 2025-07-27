module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
set.seed(123)

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(vcfR)
library(GenomicRanges)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            vcf = sprintf("aref/intermediates/%s/sniffles/sniffles.vcf", ifelse(conf$update_ref_with_tldr$per_sample == "yes", sample_table$sample_name, "A.REF")),
            r_annotation_fragmentsjoined = "aref/default/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv"
        ), env = globalenv())
        assign("outputs", list(
            absent_tes_to_mask = "aref/intermediates/A.REF/sniffles/tes_to_mask.bed"
        ), env = globalenv())
    }
)

rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)

SV <- read.vcfR(inputs$vcf)
tidyvcf_SV <- vcfR2tidy(SV)


dfSV <- tidyvcf_SV$fix %>%
    bind_cols(tidyvcf_SV$gt %>% dplyr::select(-ChromKey, -POS)) %>%
    filter(FILTER == "PASS") %>%
    filter(SVTYPE == "DEL") %>%
    mutate(seqnames = CHROM, start = POS, end = start - SVLEN) %>%
    dplyr::select(-CHROM, -POS) %>%
    dplyr::relocate(seqnames, start, end, SVLEN, SUPPORT, COVERAGE, VAF, gt_GT) %>%
    filter(SVLEN > -15000) %>%
    filter(gt_GT == "1/1") %>%
    filter(SUPPORT >= conf$rte_mask_supporting_read_count_threshold)

dfSVgrs <- GRanges(
    seqnames = dfSV$seqnames,
    ranges = IRanges(start = dfSV$start, end = dfSV$end),
    strand = "*", # or "+" / "-" if you have strand info
    mcols = dfSV[, !(names(dfSV) %in% c("seqnames", "start", "end"))]
)
dfSVgrs %>%
    width() %>%
    sum()

rmgrs <- rmfragments %>%
    filter(!grepl("__AS", gene_id)) %>%
    GRanges()

hits <- findOverlaps(rmgrs, dfSVgrs)

rm_hits <- rmgrs[queryHits(hits)]
sv_hits <- dfSVgrs[subjectHits(hits)]

overlaps <- pintersect(rm_hits, sv_hits)
percent_overlap <- width(overlaps) / width(rm_hits)
keep <- percent_overlap >= 0.05
fraction_deleted <- percent_overlap[percent_overlap >= 0.05]
tomask <- overlaps[keep]
mcols(tomask)$fraction_deleted <- fraction_deleted

maskbeddf <- as.data.frame(tomask) %>%
    tibble() %>%
    mutate(score = fraction_deleted, strand = ".") %>%
    dplyr::select(seqnames, start, end, gene_id, score, strand)

dir.create(dirname(outputs$absent_tes_to_mask), recursive = TRUE)
maskbeddf %>%
    dplyr::rename(name = gene_id) %>%
    GRanges() %>%
    rtracklayer::export(., con = outputs$absent_tes_to_mask, format = "bed")

##### NOW WRITE A DF ANNOTATING ZYGOSITY

dfSV <- tidyvcf_SV$fix %>%
    bind_cols(tidyvcf_SV$gt %>% dplyr::select(-ChromKey, -POS)) %>%
    filter(FILTER == "PASS") %>%
    filter(SVTYPE == "DEL") %>%
    mutate(seqnames = CHROM, start = POS, end = start - SVLEN) %>%
    dplyr::select(-CHROM, -POS) %>%
    dplyr::relocate(seqnames, start, end, SVLEN, SUPPORT, COVERAGE, VAF, gt_GT) %>%
    filter(SVLEN > -15000) %>%
    filter(
        case_when(
            gt_GT == "1/1" ~ SUPPORT >= conf$rte_mask_supporting_read_count_threshold,
            gt_GT == "0/1" ~ SUPPORT >= conf$rte_mask_supporting_read_count_threshold / 2
        )
    )

dfSVgrs <- GRanges(
    seqnames = dfSV$seqnames,
    ranges = IRanges(start = dfSV$start, end = dfSV$end),
    strand = "*", # or "+" / "-" if you have strand info
    mcols = dfSV[, !(names(dfSV) %in% c("seqnames", "start", "end"))]
)
dfSVgrs %>%
    width() %>%
    sum()

hits <- findOverlaps(rmgrs, dfSVgrs)

rm_hits <- rmgrs[queryHits(hits)]
sv_hits <- dfSVgrs[subjectHits(hits)]

overlaps <- pintersect(rm_hits, sv_hits)
mcols(overlaps)$gt <- mcols(sv_hits)$mcols.gt_GT
percent_overlap <- width(overlaps) / width(rm_hits)
keep <- percent_overlap >= 0.05
fraction_deleted <- percent_overlap[percent_overlap >= 0.05]
zygosite_annot <- overlaps[keep]
mcols(zygosite_annot)$fraction_deleted <- fraction_deleted

zygosite_annotdf <- zygosite_annot %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::select(seqnames, start, end, gene_id, fraction_deleted, gt)
write_csv(zygosite_annotdf, outputs$te_zygosity_annot)
