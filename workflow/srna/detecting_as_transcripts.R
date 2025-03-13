module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)


gtf <- import("srna/outs/agg/stringtie/merged.gtf")

r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- r_annotation_fragmentsjoined %>%
    left_join(r_repeatmasker_annotation)


l1hs <- rmann %>%
    filter(rte_subfamily %in% c("L1HS", "L1PA2", "L1PA3")) %>%
    filter(rte_length_req == "FL")
l1hs_asp_region <- GRanges(l1hs) %>% resize(width = 500)
l1hs_plus <- l1hs_asp_region[strand(l1hs_asp_region) == "+"]
l1hs_minus <- l1hs_asp_region[strand(l1hs_asp_region) == "-"]

strand(l1hs_plus) <- "-"
strand(l1hs_minus) <- "+"

aspgrs <- c(l1hs_plus, l1hs_minus)

gtf_start_500 <- resize(gtf, width = 500)

put_as_txid <- gtf_start_500 %>%
    subsetByOverlaps(aspgrs, ignore.strand = FALSE, minoverlap = 300) %>%
    as.data.frame() %>%
    tibble() %$% transcript_id
astxs <- gtf %>%
    as.data.frame() %>%
    tibble() %>%
    filter(transcript_id %in% put_as_txid)

astxs %>%
    filter(source == "StringTie") %>%
    filter(seqnames == "chr1")
