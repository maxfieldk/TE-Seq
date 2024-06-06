module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(ORFik)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            ref = "aref/A.REF.fa",
        txdbrefseq = "aref/A.REF_annotations/refseq.sqlite"
        ), env = globalenv())
        assign("outputs", list(
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            rmann_nonref = "aref/A.REF_annotations/A.REF_repeatmasker_rmann_nonref.csv",
            rmann = "aref/A.REF_annotations/A.REF_repeatmasker_rmann.csv"
        ), env = globalenv())
    }
)

rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)

rmfamilies <- rmfragments %>%
    dplyr::select(gene_id, family) %>%
    mutate(repeat_superfamily = case_when(
        str_detect(family, "^LINE") ~ "LINE",
        str_detect(family, "^SINE") ~ "SINE",
        str_detect(family, "^LTR") ~ "LTR",
        str_detect(family, "^Retroposon") ~ "Retroposon",
        str_detect(family, "^DNA") ~ "DNA",
        str_detect(family, "^Satellite") ~ "SAT",
    )) %>%
    mutate(repeat_superfamily = replace_na(repeat_superfamily, "Other")) %>%
    mutate(rte_superfamily = case_when(
        str_detect(family, "^LINE") ~ "LINE",
        str_detect(family, "^SINE") ~ "SINE",
        str_detect(family, "^LTR") ~ "LTR",
    )) %>%
    mutate(rte_superfamily = replace_na(rte_superfamily, "Other")) %>%
    mutate(rte_family = case_when(
        str_detect(family, "^LINE/L1") ~ "L1",
        str_detect(family, "^SINE/Alu") ~ "Alu",
        str_detect(family, "^LTR/ERV") ~ "ERV",
        str_detect(family, "^Retroposon/SVA") ~ "SVA",
    )) %>%
    mutate(rte_family = replace_na(rte_family, "Other")) %>%
    mutate(rte_subfamily = case_when(
        grepl("HERVK(.)*int$", family, perl = TRUE) ~ "HERVK_INT",
        grepl("LTR/ERVK/LTR(.)*$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5#$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5A.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5B.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5_Hs.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("HERVL(.)*int$", family, perl = TRUE) ~ "HERVL_INT",
        grepl("LTR/ERVL/LTR(.)*$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5#$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5A.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5B.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5_Hs.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("^SINE/Alu/AluY.*$", family, perl = TRUE) ~ "AluY",
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS",
        grepl("^LINE/L1/L1PA2$", family, perl = TRUE) ~ "L1PA2",
        grepl("^LINE/L1/L1PA3$", family, perl = TRUE) ~ "L1PA3",
        grepl("^LINE/L1/L1PA4$", family, perl = TRUE) ~ "L1PA4",
        grepl("^LINE/L1/L1PA5$", family, perl = TRUE) ~ "L1PA5",
        grepl("^LINE/L1/L1PA6$", family, perl = TRUE) ~ "L1PA6",
        grepl("^Retroposon/SVA/SVA_A$", family, perl = TRUE) ~ "SVA_A",
        grepl("^Retroposon/SVA/SVA_B$", family, perl = TRUE) ~ "SVA_B",
        grepl("^Retroposon/SVA/SVA_C$", family, perl = TRUE) ~ "SVA_C",
        grepl("^Retroposon/SVA/SVA_D$", family, perl = TRUE) ~ "SVA_D",
        grepl("^Retroposon/SVA/SVA_E$", family, perl = TRUE) ~ "SVA_E",
        grepl("^Retroposon/SVA/SVA_F$", family, perl = TRUE) ~ "SVA_F",
    )) %>%
    mutate(rte_subfamily = replace_na(rte_subfamily, "Other")) %>%
    mutate(rte_subfamily_limited = case_when(
        grepl("HERVK(.)*int$", family, perl = TRUE) ~ "HERVK_INT",
        grepl("LTR/ERVK/LTR(.)*$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5#$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5A.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5B.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5_Hs.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("^SINE/Alu/AluY.*$", family, perl = TRUE) ~ "AluY",
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS"
    )) %>%
    mutate(rte_subfamily_limited = replace_na(rte_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily_limited = case_when(
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS",
        grepl("^LINE/L1/L1PA2$", family, perl = TRUE) ~ "L1PA2",
        grepl("^LINE/L1/L1PA3$", family, perl = TRUE) ~ "L1PA3",
        grepl("^LINE/L1/L1PA4$", family, perl = TRUE) ~ "L1PA4",
        grepl("^LINE/L1/L1PA5$", family, perl = TRUE) ~ "L1PA5",
        grepl("^LINE/L1/L1PA6$", family, perl = TRUE) ~ "L1PA6"
    )) %>%
    mutate(l1_subfamily_limited = replace_na(l1_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily = ifelse(grepl("^LINE/L1", family, perl = TRUE), sapply(str_split(family, "/"), tail, n = 1), "Other")) %>%
    mutate(l1_subfamily = replace_na(l1_subfamily, "Other")) %>%
    mutate(herv_subfamily_limited = case_when(
        grepl("HERVK(.)*int$", family, perl = TRUE) ~ "HERVK_INT",
        grepl("LTR/ERVK/LTR(.)*$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5#$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5A.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5B.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("LTR/ERVK/LTR5_Hs.*ERVK$", family, perl = TRUE) ~ "HERVK_LTR",
        grepl("HERVL(.)*int$", family, perl = TRUE) ~ "HERVL_INT",
        grepl("LTR/ERVL/LTR(.)*$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5#$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5A.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5B.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
        grepl("LTR/ERVL/LTR5_Hs.*ERVL$", family, perl = TRUE) ~ "HERVL_LTR",
    )) %>%
    mutate(herv_subfamily_limited = replace_na(herv_subfamily_limited, "Other"))





rmlengthreq <- rmfragments %>%
    full_join(rmfamilies) %>%
    dplyr::select(gene_id, rte_subfamily, length) %>%
    mutate(rte_length = case_when(
        str_detect(rte_subfamily, "L1HS") & length > 5999 ~ "L1HS >6kb",
        str_detect(rte_subfamily, "L1HS") & length <= 5999 ~ "L1HS <6kb",
        str_detect(rte_subfamily, "L1PA2") & length > 5999 ~ "L1PA2 >6kb",
        str_detect(rte_subfamily, "L1PA2") & length <= 5999 ~ "L1PA2 <6kb",
        str_detect(rte_subfamily, "L1PA3") & length > 5999 ~ "L1PA3 >6kb",
        str_detect(rte_subfamily, "L1PA3") & length <= 5999 ~ "L1PA3 <6kb",
        str_detect(rte_subfamily, "L1PA4") & length > 5999 ~ "L1PA4 >6kb",
        str_detect(rte_subfamily, "L1PA4") & length <= 5999 ~ "L1PA4 <6kb",
        str_detect(rte_subfamily, "L1PA5") & length > 5999 ~ "L1PA5 >6kb",
        str_detect(rte_subfamily, "L1PA5") & length <= 5999 ~ "L1PA5 <6kb",
        str_detect(rte_subfamily, "L1PA6") & length > 5999 ~ "L1PA6 >6kb",
        str_detect(rte_subfamily, "L1PA6") & length <= 5999 ~ "L1PA6 <6kb",
        str_detect(rte_subfamily, "HERVL_LTR") & length > 850 ~ "HERVL_LTR >850bp",
        str_detect(rte_subfamily, "HERVL_LTR") & length <= 850 ~ "HERVL_LTR <850bp",
        str_detect(rte_subfamily, "HERVL_INT") & length > 6999 ~ "HERVL_INT >7kb",
        str_detect(rte_subfamily, "HERVL_INT") & length <= 6999 ~ "HERVL_INT <7kb",
        str_detect(rte_subfamily, "HERVK_LTR") & length > 850 ~ "HERVK_LTR >850bp",
        str_detect(rte_subfamily, "HERVK_LTR") & length <= 850 ~ "HERVK_LTR <850bp",
        str_detect(rte_subfamily, "HERVK_INT") & length > 6999 ~ "HERVK_INT >7kb",
        str_detect(rte_subfamily, "HERVK_INT") & length <= 6999 ~ "HERVK_INT <7kb",
        str_detect(rte_subfamily, "AluY") & length > 299 ~ "AluY >300bp",
        str_detect(rte_subfamily, "AluY") & length <= 299 ~ "AluY <300bp",
        str_detect(rte_subfamily, "SVA_A") & length > 1999 ~ "SVA_A >2kb",
        str_detect(rte_subfamily, "SVA_A") & length <= 1999 ~ "SVA_A <2kb",
        str_detect(rte_subfamily, "SVA_B") & length > 1999 ~ "SVA_B >2kb",
        str_detect(rte_subfamily, "SVA_B") & length <= 1999 ~ "SVA_B <2kb",
        str_detect(rte_subfamily, "SVA_C") & length > 1999 ~ "SVA_C >2kb",
        str_detect(rte_subfamily, "SVA_C") & length <= 1999 ~ "SVA_C <2kb",
        str_detect(rte_subfamily, "SVA_D") & length > 1999 ~ "SVA_D >2kb",
        str_detect(rte_subfamily, "SVA_D") & length <= 1999 ~ "SVA_D <2kb",
        str_detect(rte_subfamily, "SVA_E") & length > 1999 ~ "SVA_E >2kb",
        str_detect(rte_subfamily, "SVA_E") & length <= 1999 ~ "SVA_E <2kb",
        str_detect(rte_subfamily, "SVA_F") & length > 1999 ~ "SVA_F >2kb",
        str_detect(rte_subfamily, "SVA_F") & length <= 1999 ~ "SVA_F <2kb",
        str_detect(rte_subfamily, "Other") ~ "Other",
    )) %>%
    dplyr::select(gene_id, rte_length)

# Annotate Intactness
fa <- Rsamtools::FaFile(inputs$ref)

rmfragments %$% refstatus %>% unique()
ranges <- GRanges(rmfragments)
l1hsranges <- ranges[grepl("L1HS", ranges$gene_id)]
l1pa2ranges <- ranges[grepl("L1PA2", ranges$gene_id)]
l1hspa2ranges <- c(l1hsranges, l1pa2ranges)
gene_ids <- l1hspa2ranges$gene_id
l1hspa2ss <- getSeq(fa, l1hspa2ranges)
mcols(l1hspa2ss) <- mcols(l1hspa2ranges)
names(l1hspa2ss) <- mcols(l1hspa2ss)$gene_id
flss <- l1hspa2ss[width(l1hspa2ss) > 5999]

# orf analysis
pos <- ORFik::findORFs(flss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
names(pos) <- names(flss[as.integer(names(pos))])
# FILTERING ON THE BASIS OF ORF LENGTHS
pass <- c()
for (element in names(pos)) {
    if (3828 %in% width(pos[[element]]) & 1017 %in% width(pos[[element]])) {
        pass <- c(pass, element)
    }
}
l1hspass <- pass[pass %>% str_detect("L1HS")]
l1pa2pass <- pass[pass %>% str_detect("L1PA2")]

intactness <- rmfragments %>%
    dplyr::select(gene_id) %>%
    mutate(l1_intactness = case_when(
        gene_id %in% l1hspass ~ "L1HS ORFs Intact",
        gene_id %in% l1pa2pass ~ "L1PA2 ORFs Intact",
        TRUE ~ "Other"
    ))


req_annot <- full_join(intactness, rmlengthreq)
req_annot <- req_annot %>% mutate(req_integrative = case_when(
    l1_intactness == "L1HS ORFs Intact" ~ "Young ORFs Intact",
    l1_intactness == "L1PA2 ORFs Intact" ~ "Young ORFs Intact",
    str_detect(rte_length, ">") ~ "Young Full Length",
    str_detect(rte_length, "<") ~ "Young Truncated",  
    TRUE ~ "Old"
))

ltrs <- rmfamilies %>%
    filter(rte_subfamily == "HERVK_LTR") %>%
    left_join(rmfragments)
ltrsgrs <- GRanges(ltrs)
ltrsgrsextended <- resize(ltrsgrs, width = width(ltrsgrs) + (200 * 2), fix = "center")

ints <- rmfamilies %>%
    filter(rte_subfamily == "HERVK_INT") %>%
    left_join(rmfragments)
intsfl <- ints %>% filter(length > 6999)
intsgrs <- GRanges(ints)
intsflgrs <- GRanges(intsfl)
intsnotflgrs <- intsgrs %>% subsetByOverlaps(intsflgrs, invert = TRUE)

Solo_LTR <- ltrsgrsextended %>%
    subsetByOverlaps(intsgrs, invert = TRUE) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Fl_Provirus_5LTR <- ltrsgrs %>%
    flank(200) %>%
    subsetByOverlaps(intsflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Fl_Provirus_3LTR <- ltrsgrs %>%
    flank(200, start = FALSE) %>%
    subsetByOverlaps(intsflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Provirus_5LTR <- ltrsgrs %>%
    flank(200) %>%
    subsetByOverlaps(intsnotflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Provirus_3LTR <- ltrsgrs %>%
    flank(200, start = FALSE) %>%
    subsetByOverlaps(intsnotflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id

NOT_LTR_5Flanked_INT <- intsgrs %>%
    flank(200) %>%
    subsetByOverlaps(ltrsgrs, invert = TRUE) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
LTR_5Flanked_INT <- intsgrs %>%
    flank(200) %>%
    subsetByOverlaps(ltrsgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id

ltr_viral_status <- rmfragments %>%
    dplyr::select(gene_id) %>%
    mutate(ltr_viral_status = case_when(
        gene_id %in% Solo_LTR ~ "Solo_LTR",
        gene_id %in% Provirus_5LTR ~ "NOT_Fl_Provirus_5LTR",
        gene_id %in% Provirus_3LTR ~ "NOT_FL_Provirus_3LTR",
        gene_id %in% Fl_Provirus_5LTR ~ "Fl_Provirus_5LTR",
        gene_id %in% Fl_Provirus_3LTR ~ "Fl_Provirus_3LTR",
        gene_id %in% NOT_LTR_5Flanked_INT ~ "NOT_LTR_5Flanked_INT",
        gene_id %in% LTR_5Flanked_INT ~ "LTR_5Flanked_INT",
        TRUE ~ "Other"
    ))

# rmfamilies <- read_csv(outputs$r_annotation_families, col_names = TRUE)

# for annotation purposes, I will have to have the location of nonreference inserts be their insertion site
rmfragments_ref <- rmfragments %>% filter(!str_detect(seqnames, "nonref"))
rmfragments_nonref <- rmfragments %>% filter(str_detect(seqnames, "nonref"))
seqnamesNonref <- rmfragments_nonref %$% seqnames
insert_seqnames <- c()
insert_start <- c()
insert_end <- c()
for (seqnames in seqnamesNonref) {
    split <- str_split(seqnames, "_") %>% unlist()
    splitlen <- length(split)
    insert_seqnames <- c(insert_seqnames, split[splitlen - 2])
    insert_start <- c(insert_start, split[splitlen - 1])
    insert_end <- c(insert_end, split[splitlen])
}
rmfragments_nonref$insert_seqnames <- insert_seqnames
rmfragments_nonref$insert_start <- insert_start
rmfragments_nonref$insert_end <- insert_end

rmfragments_refgr <- GRanges(rmfragments_ref)
rmfragmentsgr_properinsertloc <- rmfragments_refgr
tryCatch(
    {
        rmfragments_nonrefgr <- GRanges(rmfragments_nonref %>% dplyr::select(-seqnames, -start, -end) %>% dplyr::relocate(gene_id, insert_seqnames, source, type, insert_start, insert_end, strand) %>% dplyr::rename(seqnames = insert_seqnames, start = insert_start, end = insert_end))
        rmfragmentsgr_properinsertloc <- c(rmfragments_refgr, rmfragments_nonrefgr)
    },
    error = function(e) {
    }
)

library(GenomicFeatures)
library(genomation)

# genes, centromere and telomere
refseq <- import(conf$ref_refseq_gtf)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
noncoding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NR", mcols(refseq)$transcript_id))]
transcripts <- c(coding_transcripts, noncoding_transcripts)

coding_transcript_upstream <- coding_transcripts %>% flank(10000)
noncoding_transcript_upstream <- noncoding_transcripts %>% flank(10000)
coding_transcript_downstream <- coding_transcripts %>% flank(10000,start=FALSE)
noncoding_transcript_downstream <- noncoding_transcripts %>% flank(10000,start=FALSE)

coding_transcript_adjacent <- c(coding_transcript_upstream, coding_transcript_downstream)
noncoding_transcript_adjacent <- c(noncoding_transcript_upstream, noncoding_transcript_downstream)
 

txdb <- loadDb(inputs$txdbrefseq)

introns <- intronsByTranscript(txdb, use.names = TRUE)
introns <- introns[grepl("^NM", names(introns))]
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
exons <- exons[grepl("^NM", names(exons))]
fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
fiveUTRs <- fiveUTRs[grepl("^NM", names(fiveUTRs))]
threeUTRs <- threeUTRsByTranscript(txdb, use.names = TRUE)
threeUTRs <- threeUTRs[grepl("^NM", names(threeUTRs))]


cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
cytobands <- cytobandsdf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
centromere <- cytobandsdf %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
telomeredf <- read_delim(conf$ref_telomere, col_names = FALSE, delim = "\t")
telomere <- telomeredf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()

getannotation <- function(to_be_annotated, regions_of_interest, property, name_in, name_out) {
    inregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = FALSE)
    inregions$prop <- name_in
    outregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = TRUE)
    outregions$prop <- name_out
    merged <- c(inregions, outregions)
    annot <- merged %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(gene_id, prop)
    annot[[property]] <- annot$prop
    annot <- annot %>% dplyr::select(-prop)
    return(annot)
}
genic_annot <- getannotation(rmfragmentsgr_properinsertloc, transcripts, "genic", "Genic", "Intergenic")
coding_tx_annot <- getannotation(rmfragmentsgr_properinsertloc, coding_transcripts, "coding_tx", "CodingTx", "NonCodingTx")
noncoding_tx_annot <- getannotation(rmfragmentsgr_properinsertloc, noncoding_transcripts, "noncoding_tx", "NoncodingTx", "NonNonCodingTx")
coding_tx_adjacent_annot <- getannotation(rmfragmentsgr_properinsertloc, coding_transcript_adjacent, "coding_tx_adjacent", "CodingTxAdjacent", "NonCodingTxAdjacent")
noncoding_tx_adjacent_annot <- getannotation(rmfragmentsgr_properinsertloc, noncoding_transcript_adjacent, "noncoding_tx_adjacent", "NoncodingTxAdjacent", "NonNonCodingTxAdjacent")
exonic_annot <- getannotation(rmfragmentsgr_properinsertloc, exons, "exonic", "Exonic", "NonExonic")
intron_annot <- getannotation(rmfragmentsgr_properinsertloc, introns, "intronic", "Intronic", "NonIntronic")
utr5_annot <- getannotation(rmfragmentsgr_properinsertloc, fiveUTRs, "utr5", "5UTR", "Non5UTR")
utr3_annot <- getannotation(rmfragmentsgr_properinsertloc, threeUTRs, "utr3", "3UTR", "Non3UTR")

cent_annot <- getannotation(rmfragmentsgr_properinsertloc, centromere, "centromeric", "Centromeric", "NonCentromeric")
telo_annot <- getannotation(rmfragmentsgr_properinsertloc, telomere, "telomeric", "Telomeric", "NonTelomeric")
genic_annot %$% genic %>% table()
cent_annot %$% centromeric %>% table()
telo_annot %$% telomeric %>% table()



region_annot <- full_join(genic_annot, coding_tx_annot) %>%
    full_join(noncoding_tx_annot) %>%
    full_join(coding_tx_adjacent_annot) %>%
    full_join(noncoding_tx_adjacent_annot) %>%
    full_join(cent_annot) %>%
    full_join(telo_annot) %>%
    full_join(exonic_annot) %>%
    full_join(intron_annot) %>%
    full_join(utr5_annot) %>%
    full_join(utr3_annot)


region_annot <- region_annot %>% mutate(loc_integrative = case_when(
    exonic == "Exonic" ~ "Exonic",
    utr5 == "5utr" ~ "5utr",
    utr3 == "3utr" ~ "3utr",
    intronic == "Intronic" ~ "Intronic",
    noncoding_tx == "NoncodingTx" ~ "NoncodingTx",
    coding_tx == "CodingTx" ~ "CodingTxOther",
    coding_tx_adjacent == "CodingTxAdjacent" ~ "CodingTxAdjacent",
    noncoding_tx_adjacent == "NoncodingTxAdjacent" ~ "NoncodingTxAdjacent",
    centromeric == "Centromeric" ~ "Centromeric",
    telomeric == "Telomeric" ~ "Telomeric",
    genic == "Intergenic" ~ "Intergenic",
    TRUE ~ "Other"
)) %>% mutate(loc_lowres_integrative = case_when(
    loc_integrative == "Centromeric" ~ "Intergenic",
    loc_integrative == "CodingTxAdjacent" ~ "Gene Adjacent",
    loc_integrative == "Intron" ~ "Genic",
    loc_integrative == "Exonic" ~ "Genic",
    loc_integrative == "Intronic" ~ "Genic",
    loc_integrative == "NoncodingTx" ~ "Genic",
    loc_integrative == "NoncodingTxAdjacent" ~ "Gene Adjacent",
    loc_integrative == "Intergenic" ~ "Intergenic",
    loc_integrative == "Telomeric" ~ "Intergenic",
    TRUE ~ "Other"))

# HERV ltr status
herv_ltr_map = tibble(
    name = c("HML2","HML2","HML2","HML2",
    "HML3",
    "HERVW",
    "HERVH","HERVH","HERVH","HERVH","HERVH",
    "HERV9","HERV9","HERV9","HERV9","HERV9","HERV9","HERV9",
    "HERV-FC1","HERV-FC1","HERV-FC1",
    "HERV-FC2"
    ),
    int = c("HERVK-int", "HERVK-int", "HERVK-int", "HERVK-int",
    "HERVK9-int",
    "HERV17-int",
    "HERVH-int","HERVH-int","HERVH-int","HERVH-int","HERVH-int",
    "HERV9-int","HERV9-int","HERV9-int","HERV9-int","HERV9-int","HERV9-int","HERV9-int",
    "HERV-Fc1-int","HERV-Fc1-int","HERV-Fc1-int",
    "HERV-Fc2-int"),
    ltr = c("LTR5", "LTR5A", "LTR5B", "LTR5_Hs",
    "MER9", 
    "LTR17",
    "LTR7", "LTR7A", "LTR7B", "LTR7C", "LTR7Y",
    "LTR12", "LTR12B", "LTR12C", "LTR12D", "LTR12E", "LTR12F", "LTR12_",
    "HERV-Fc1_LTR1", "HERV-Fc1_LTR2", "HERV-Fc1_LTR3",
    "HERV-Fc2_LTR")
)
group_frame <- tibble()
group_frame_any_ltr <- tibble()
for (provirus_type in herv_ltr_map %$% name %>% unique()) {
    int <- herv_ltr_map %>% filter(name == provirus_type) %$% int %>% unique()
    compatible_ltrs <- herv_ltr_map %>% filter(name == provirus_type) %$% ltr %>% unique()
    int_grs <- rmfragments %>% filter(grepl(sprintf("%s_.*", int), gene_id)) %>% GRanges()
    ltr_grs <- rmfragments %>% filter(grepl(paste0(compatible_ltrs,"_.*") %>% str_c(collapse = "|"), gene_id)) %>% GRanges()
    ltrs_all <- rmfragments %>% filter(grepl("LTR.*", gene_id)) %>% GRanges()
    for (int_id in mcols(int_grs)$gene_id) {
        print(int_id)
        group_id <- int_id %>% gsub(".*-int", provirus_type, .)
        int_gr <- int_grs[int_grs$gene_id == int_id]
        int_gr_extended <- resize(int_gr, width = width(int_gr) + (1000 * 2), fix = "center")
        ltr_ol <- subsetByOverlaps(ltr_grs, int_gr_extended, type = "any", minoverlap = 1)
        ltr_any_ol <- subsetByOverlaps(ltrs_all, int_gr_extended, type = "any", minoverlap = 1)
        if (length(ltr_ol) != 0) {
            gdf <- tibble(gene_id = c(int_id, mcols(ltr_ol)$gene_id))
            gdf$ltr_proviral_group_id <- group_id
            group_frame <<- bind_rows(group_frame, gdf)
        }
        if (length(ltr_any_ol) != 0) {
            gdf <- tibble(gene_id = c(int_id, mcols(ltr_ol)$gene_id))
            gdf$ltr_proviral_group_id <- group_id
            group_frame_any_ltr <<- bind_rows(group_frame_any_ltr, gdf)
        }
    }
}

ltr_proviral_groups <- group_frame
ltr_proviral_groups <- ltr_proviral_groups[!(ltr_proviral_groups %$% gene_id %>% duplicated()),] #for ltrs which overlap 2 ints, in order not to duplicate the feature I remove them from their second proviral group affiliation


annots <- rmfamilies %>%
    full_join(req_annot %>% rename_at(vars(-gene_id, -req_integrative), ~ paste0(., "_req"))) %>%
    full_join(ltr_viral_status %>% rename_at(vars(-gene_id), ~ paste0(., "_req"))) %>%
    full_join(ltr_proviral_groups) %>%
    full_join(region_annot %>% rename_at(vars(-gene_id, -loc_integrative), ~ paste0(., "_loc")))



region_annot %$% gene_id %>% duplicated() %>% sum()
ltr_proviral_groups %$% gene_id %>% duplicated() %>% sum()
ltr_viral_status %$% gene_id %>% duplicated() %>% sum()
req_annot %$% gene_id %>% duplicated() %>% sum()

rmfamilies %$% gene_id %>% duplicated() %>% sum()
annots %$% gene_id %>% duplicated() %>% sum()

rmann <- left_join(rmfragments, annots) 

write_csv(annots, outputs$r_repeatmasker_annotation)
write_csv(rmann, outputs$rmann)
write_csv(rmann %>% filter(refstatus == "NonRef"), outputs$rmann_nonref)
