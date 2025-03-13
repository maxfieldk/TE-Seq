module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
set.seed(123)

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(Biostrings)
library(Rsamtools)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            repmask_gff3 = "aref/default/A.REF_annotations/A.REF_repeatmasker.gff3",
            r_repeatmasker_annotation = "aref/default/A.REF_annotations/A.REF_repeatmasker_annotation.csv"
        ), env = globalenv())
        assign("outputs", list(
            repmask_gtf = "aref/default/A.REF_annotations/A.REF_repeatmasker_tetranscripts.gtf"
        ), env = globalenv())
    }
)


gff3df <- import(inputs$repmask_gff3) %>%
    as.data.frame() %>%
    tibble()

r_repeatmasker_annotation <- read_csv(inputs$r_repeatmasker_annotation) %>% dplyr::select(gene_id, rte_family, rte_superfamily)
gff3dfann <- gff3df %>%
    dplyr::select(-Target, -pctdiv, -pctdel) %>%
    left_join(r_repeatmasker_annotation)
tetranscriptsgr <- gff3dfann %>%
    mutate(
        split_ids = str_split(old_id, pattern = "/"),
        gene_id = map_chr(split_ids, ~ .x[length(.x) - 1]),
        family_id = map_chr(split_ids, ~ ifelse(length(.x) >= 2, .x[length(.x) - 2], .x[length(.x) - 1])),
        class_id = map_chr(split_ids, ~ ifelse(length(.x) >= 3, .x[length(.x) - 3],
            ifelse(length(.x) >= 2, .x[length(.x) - 2], .x[length(.x) - 1])
        ))
    ) %>%
    mutate(gene_name = paste0(map_chr(str_split(old_id, pattern = "/"), ~ .x[length(.x) - 1]), ":TE")) %>%
    dplyr::select(-old_id, -split_ids) %>%
    GRanges()

tetranscriptsgr <- gff3dfann %>%
    mutate(gene_id = map_chr(str_split(old_id, pattern = "/"), ~ .x[length(.x) - 1])) %>%
    mutate(family_id = rte_family) %>%
    mutate(class_id = rte_superfamily) %>%
    mutate(gene_name = paste0(gene_id, ":TE")) %>%
    dplyr::select(-rte_family, -rte_superfamily, -old_id) %>%
    filter(family_id != "Other") %>%
    GRanges()


rtracklayer::export(tetranscriptsgr, con = outputs$repmask_gtf, format = "GFF2")
