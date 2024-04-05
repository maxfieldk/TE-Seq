source("workflow/scripts/defaults.R")
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
library(rBLAST)

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")[["aref"]]
)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            filtered_tldr = "aref/tldr/tldr.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/annotations/repeatmasker_annotation.csv",
            ref = "aref/A.REF.fa",
            blast_njs = "aref/A.REF.njs"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/RefAnalysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/RefAnalysis/transduction_df.csv"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

df_filtered <- read_delim(inputs$filtered_tldr)

trsd <- df_filtered %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd %>% head(n = 4) %$% Transduction_3p


trsd_ss <- DNAStringSet(trsd$Transduction_3p)
names(trsd_ss) <- trsd$UUID


at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)

trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]
bl <- blast(db = "aref/ref")
bres <- tibble(predict(bl, trsd_ss_for_blast)) %>% left_join(rmfragments, by = c("QueryID" = "gene_id"))


bres %$% QueryID %>%
    unique() %>%
    length()

bres_hits <- bres %>%
    filter(!grepl("nonref", SubjectID)) %>%
    group_by(QueryID) %>%
    filter(Bits == max(Bits)) %>%
    filter(Perc.Ident == max(Perc.Ident)) %>%
    filter(Gap.Openings == min(Gap.Openings)) %>%
    filter(Alignment.Length == max(Alignment.Length))


l1grs <- GRanges(rmann %>% filter(grepl("L1HS|L1PA[2]", rte_subfamily)))

bres_hits_grs_prep <- bres_hits %>%
    dplyr::select(SubjectID, S.start, S.end) %>%
    group_by(QueryID) %>%
    mutate(seqnames = SubjectID) %>%
    mutate(start = min(S.start, S.end), end = max(S.start, S.end)) %>%
    dplyr::select(seqnames, start, end, QueryID)
bresgrs <- GRanges(bres_hits_grs_prep)
# extend by 500 bp on either side
bresgrs <- resize(bresgrs, width = width(bresgrs) + 2000, fix = "center")

subsetByOverlaps(bresgrs, l1grs)


save(mysaveandstoreplots, file = outputs$plots)
