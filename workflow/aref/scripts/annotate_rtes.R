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
    )) %>% mutate(repeat_superfamily = replace_na(repeat_superfamily, "Other")) %>%
    mutate(rte_superfamily = case_when(
        str_detect(family, "^LINE") ~ "LINE",
        str_detect(family, "^SINE") ~ "SINE",
        str_detect(family, "^LTR") ~ "LTR",
    )) %>% mutate(rte_superfamily = replace_na(rte_superfamily, "Other")) %>%
    mutate(rte_family = case_when(
        str_detect(family, "^LINE/L1") ~ "L1",
        str_detect(family, "^SINE/Alu") ~ "Alu",
        str_detect(family, "^LTR/ERV") ~ "ERV",
        str_detect(family, "^Retroposon/SVA") ~ "SVA",
    )) %>% mutate(rte_family = replace_na(rte_family, "Other"))

if (conf$species == "human") {
    rmfamilies <- rmfamilies %>%
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
} else if (conf$species == "mouse") {
    rmfamilies <- rmfamilies %>%
    mutate(rte_subfamily = case_when(
        rte_superfamily == "LINE" & grepl("L1MdTf", gene_id, perl = TRUE) ~ "L1MdTf",
        rte_superfamily == "LINE" & grepl("L1MdGf", gene_id, perl = TRUE) ~ "L1MdGf",
        rte_superfamily == "LINE" & grepl("L1MdA", gene_id, perl = TRUE) ~ "L1MdA",
        rte_superfamily == "LINE" & grepl("L1MdF", gene_id, perl = TRUE) ~ "L1MdF",
        rte_superfamily == "SINE" & grepl("B1", gene_id, perl = TRUE) ~ "B1",
        rte_superfamily == "SINE" & grepl("B2", gene_id, perl = TRUE) ~ "B2",
        rte_superfamily == "SINE" & grepl("ID", gene_id, perl = TRUE) ~ "ID",
        rte_superfamily == "SINE" & grepl("B4", gene_id, perl = TRUE) ~ "B4",
        rte_superfamily == "LTR" & grepl("MMTV-int", gene_id, perl = TRUE) ~ "MMTV-int",
        rte_superfamily == "LTR" & grepl("ETn.*int", gene_id, perl = TRUE) ~ "ETn-int",
        rte_superfamily == "LTR" & grepl("IAP.*int", gene_id, perl = TRUE) ~ "IAP-int",
        rte_superfamily == "LTR" & grepl("MMVL.*int", gene_id, perl = TRUE) ~ "MMVL-int",
        rte_superfamily == "LTR" & grepl("MMERGLN.*int", gene_id, perl = TRUE) ~ "MMERGLN-int",
        rte_superfamily == "LTR" & grepl("MuRRS.*int", gene_id, perl = TRUE) ~ "MuRRS-int",
        rte_superfamily == "LTR" & grepl("MuLV.*int", gene_id, perl = TRUE) ~ "MuLV-int",
        rte_superfamily == "LTR" & grepl("MERVL.*int", gene_id, perl = TRUE) ~ "MERVL-int",
    )) %>% mutate(rte_subfamily = replace_na(rte_subfamily, "Other")) %>%
    mutate(rte_subfamily_limited = case_when(
        rte_superfamily == "LINE" & grepl("L1MdTf_I", gene_id, perl = TRUE) ~ "L1MdTf_I",
        rte_superfamily == "LINE" & grepl("L1MdTf_II", gene_id, perl = TRUE) ~ "L1MdTf_II",
        rte_superfamily == "LINE" & grepl("L1MdTf_III", gene_id, perl = TRUE) ~ "L1MdTf_III",
        rte_superfamily == "LINE" & grepl("L1MdGf_I", gene_id, perl = TRUE) ~ "L1MdGf_I",
        rte_superfamily == "LINE" & grepl("L1MdGf_II", gene_id, perl = TRUE) ~ "L1MdGf_II",
        rte_superfamily == "LINE" & grepl("L1MdGf_III", gene_id, perl = TRUE) ~ "L1MdGf_III",
        rte_superfamily == "LINE" & grepl("L1MdA_I", gene_id, perl = TRUE) ~ "L1MdA_I",
        rte_superfamily == "LINE" & grepl("L1MdA_II", gene_id, perl = TRUE) ~ "L1MdA_II",
        rte_superfamily == "LINE" & grepl("L1MdA_III", gene_id, perl = TRUE) ~ "L1MdA_III",
        rte_superfamily == "SINE" & grepl("B1", gene_id, perl = TRUE) ~ "B1",
        rte_superfamily == "SINE" & grepl("B2", gene_id, perl = TRUE) ~ "B2",
        rte_superfamily == "SINE" & grepl("ID", gene_id, perl = TRUE) ~ "ID",
        rte_superfamily == "SINE" & grepl("B4", gene_id, perl = TRUE) ~ "B4",
        rte_superfamily == "LTR" & grepl("MMTV-int", gene_id, perl = TRUE) ~ "MMTV-int",
        rte_superfamily == "LTR" & grepl("ETn.*int", gene_id, perl = TRUE) ~ "ETn-int",
        rte_superfamily == "LTR" & grepl("IAP.*int", gene_id, perl = TRUE) ~ "IAP-int",
        rte_superfamily == "LTR" & grepl("MMVL.*int", gene_id, perl = TRUE) ~ "MMVL-int",
        rte_superfamily == "LTR" & grepl("MMERGLN.*int", gene_id, perl = TRUE) ~ "MMERGLN-int",
        rte_superfamily == "LTR" & grepl("MuRRS.*int", gene_id, perl = TRUE) ~ "MuRRS-int",
        rte_superfamily == "LTR" & grepl("MuLV.*int", gene_id, perl = TRUE) ~ "MuLV-int",
        rte_superfamily == "LTR" & grepl("MERVL.*int", gene_id, perl = TRUE) ~ "MERVL-int",
    )) %>% mutate(rte_subfamily_limited = replace_na(rte_subfamily_limited, "Other")) %>%
        mutate(l1_subfamily_limited = case_when(
        rte_superfamily == "LINE" & grepl("L1MdTf", gene_id, perl = TRUE) ~  gsub(".*/", "", family),
        rte_superfamily == "LINE" & grepl("L1MdGf", gene_id, perl = TRUE) ~ gsub(".*/", "", family),
        rte_superfamily == "LINE" & grepl("L1MdA", gene_id, perl = TRUE) ~ gsub(".*/", "", family),
        rte_superfamily == "LINE" & grepl("L1MdF", gene_id, perl = TRUE) ~ gsub(".*/", "", family),
    )) %>% mutate(l1_subfamily_limited = replace_na(l1_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily = ifelse(grepl("^LINE/L1", family, perl = TRUE), sapply(str_split(family, "/"), tail, n = 1), "Other")) %>%
    mutate(l1_subfamily = replace_na(l1_subfamily, "Other"))
} else {
    tempdf <- rmfamilies %>% left_join(rmfragments) %>% 
        group_by(rte_superfamily, family) %>% 
        summarise(pctdiv = mean(pctdiv)) %>% 
        filter(rte_superfamily != "Other") %>% 
        arrange(pctdiv)
    least_diverged_families <- tempdf %>% group_by(rte_superfamily) %>%
        slice_head(n = 7) %$% family
    least_diverged_families_stringent <- tempdf %>% group_by(rte_superfamily) %>%
        slice_head(n = 3) %$% family
    least_diverged_families_l1 <- tempdf %>% group_by(rte_superfamily) %>% filter(grepl("LINE", rte_superfamily)) %>%
        slice_head(n = 7) %$% family

    rmfamilies <- rmfamilies %>%
    mutate(rte_subfamily = ifelse(family %in% least_diverged_families, gsub(".*/", "", family), "Other")) %>% 
    mutate(rte_subfamily = replace_na(rte_subfamily, "Other")) %>%
    mutate(rte_subfamily_limited = ifelse(family %in% least_diverged_families_stringent, gsub(".*/", "", family), "Other")) %>% 
    mutate(rte_subfamily_limited = replace_na(rte_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily_limited = ifelse(family %in% least_diverged_families_l1, gsub(".*/", "", family), "Other")) %>% 
    mutate(l1_subfamily_limited = replace_na(rte_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily = ifelse(grepl("^LINE/L1", family, perl = TRUE), sapply(str_split(family, "/"), tail, n = 1), "Other")) %>%
    mutate(l1_subfamily = replace_na(l1_subfamily, "Other"))
}

# Annotate Intactness
fa <- Rsamtools::FaFile(inputs$ref)

rmfragments %$% refstatus %>% unique()
if (conf$species == "human") {
element_to_annotate <- c("L1HS","L1PA2")
} else if (conf$species == "mouse") {
element_to_annotate <- c("L1MdTf_I","L1MdTf_II","L1MdTf_III",
    "L1MdGf_I","L1MdGf_II",
    "L1MdA_I", "L1MdA_II", "L1MdA_III")
}
#trycatch needed if you are not using human/mouse
tryCatch({

element_info_list <- list()
for (element in element_to_annotate) {
    active_family_ranges <- GRanges(rmfragments %>% filter(grepl(paste0(element, "$"), family)))
    gene_ids <- active_family_ranges$gene_id
    active_family_ss <- getSeq(fa, active_family_ranges)
    mcols(active_family_ss) <- mcols(active_family_ranges)
    names(active_family_ss) <- mcols(active_family_ss)$gene_id
    flss <- active_family_ss[as.numeric(mcols(active_family_ss)$pctconsensuscovered) >= 95]

    # orf analysis
    z_score_cutoff <- 4
    pos <- ORFik::findORFs(flss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
    names(pos) <- names(flss[as.integer(names(pos))])

    orf_frame <- as.data.frame(pos) %>% tibble() %>% 
        mutate(bin_10 = cut(width, breaks = seq(min(width), max(width), by = 1)))

    binwidth_df <- orf_frame %>%
        group_by(width) %>%
        mutate(average_start = mean(start)) %>%
        summarise(n = n(), average_start = dplyr::first(average_start)) %>% 
        mutate(orf_length_zscore = (n - mean(n)) / sd(n))

    p <- ggplot(binwidth_df) + geom_point(aes(x = width, y = orf_length_zscore, color = orf_length_zscore > z_score_cutoff)) +
        geom_hline(yintercept = z_score_cutoff, color = "red") + mtopen +
        theme(legend.position = "none") + labs(x = "ORF Length", y = "Z-Score", title = element) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))
    mysave(sprintf("aref/A.REF_annotations/figures/%s_orf_lengths.pdf", element), 4, 4)
    modal_widths <- binwidth_df %>% filter(orf_length_zscore > z_score_cutoff) %>% arrange(-average_start) %$% width

    element_orf_intactness_df <- tibble(gene_id = names(flss))
    i = 1
    for (modal_width in modal_widths) {
        tdf <- orf_frame %>% filter(width == modal_width)
        tss <- flss[names(flss) %in% (tdf %$% group_name)]
        tdf <- tdf[match(names(tss), tdf %$% group_name),]
        names(tss) == tdf %$% group_name
        orf_ss <- subseq(tss, start = tdf$start, end = tdf$end)
        #need to remove all ambiguity prior to tranlsation
        consensus <- chartr("N", "A", replaceAmbiguities(DNAString(consensusString(orf_ss)), new="N"))
        consensus_ss <- DNAStringSet(consensus)
        names(consensus_ss) <- c("consensus")
        consensus_aa <- Biostrings::translate(consensus)
        consensus_aa_ss <- AAStringSet(consensus_aa)
        names(consensus_aa_ss) <- c("consensus")
        orf_consensus_path <- sprintf("aref/A.REF_annotations/figures/%s_orf_length_%s_consensus.fa", element, modal_width)
        orf_aa_consensus_path <- sprintf("aref/A.REF_annotations/figures/%s_orf_length_%s_aa_consensus.fa", element, modal_width)
        writeXStringSet(consensus_ss, file = orf_consensus_path)
        writeXStringSet(consensus_aa_ss, file = orf_aa_consensus_path)

        #get all orfs that are +/- 10% of the modal width
        tdf <- orf_frame %>% filter(width > modal_width*0.9 & width < modal_width*1.1)
        tss <- flss[names(flss) %in% (tdf %$% group_name)]
        tdf <- tdf[match(names(tss), tdf %$% group_name),]
        names(tss) == tdf %$% group_name
        orf_ss <- subseq(tss, start = tdf$start, end = tdf$end)


        library(seqinr)
        orf_fa_path <- sprintf("aref/A.REF_annotations/figures/%s_orf_length_around_%s_aa.fa", element, modal_width)
        writeXStringSet(c(Biostrings::translate(orf_ss), consensus_aa_ss), file = orf_fa_path)
        system(sprintf("echo $(pwd); mafft --auto %s > %s.aln.fa", orf_fa_path, orf_fa_path))
        aln <- read.alignment(sprintf("%s.aln.fa", orf_fa_path), format = "fasta")
        d <- dist.alignment(aln, "identity", gap = FALSE) %>% as.matrix()
        d_gapped <- dist.alignment(aln, "identity", gap = TRUE) %>% as.matrix()

        #pct identity to aa consensus
        dsqaured <- d["consensus",]**2
        d_gappedsquared <- d_gapped["consensus",]**2
        # tgz_file <- rBLAST::blast_db_get("pdbaa.tar.gz")
        # untar(tgz_file, exdir = "aref/blastdb/pdbaa")
        # library(rBLAST)
        # bl <- blast(db = "./aref/blastdb/pdbaa/pdbaa")
        # bres <- tibble(predict(bl, consensus_aa_ss, )) 

        tdf <- as.data.frame(d_gappedsquared)
        colnames(tdf) <- c(paste0("orf",i))
        tdf$gene_id <- rownames(tdf)
        tdf <- tibble(tdf)
        tdf <- full_join(element_orf_intactness_df, tdf) %>% filter(gene_id != "consensus")
        element_orf_intactness_df <- tdf
        i = i + 1
}
element_info_list[[element]] <- element_orf_intactness_df
}
#select all columns that start with orf and include their values in a tuple

rm(intactness_ann)
for (i in 1:length(element_info_list)) {
    tdf <- element_info_list[[i]] %>% mutate(across(starts_with("orf"), ~ replace_na(., 1))) %>% mutate(orf_distance_tuple = pmap(dplyr::select(., starts_with("orf")), c))  %>% dplyr::select(gene_id, orf_distance_tuple) 
    if (!exists("intactness_ann")) {
        intactness_ann <- tdf
    } else {
        intactness_ann <- full_join(intactness_ann, tdf)
    }
}

intactness_ann <- intactness_ann %>% mutate(intactness_req = ifelse(sapply(orf_distance_tuple, function(x) all(x < 0.05)), "Intact", "Not Intact"))%>%
  mutate(orf_passes_distance_threshold = map(orf_distance_tuple, ~ {
    threshold_indices <- which(.x < 0.05)
    if(length(threshold_indices) == 0) {
      "none"
    } else {
      paste0("orf", threshold_indices)
    }
  }))

}, error = function(e) {
    print("no elements annotated for intactness")
    intactness_ann <- rmfragments %>% dplyr::select(gene_id) %>% mutate(intactness_req = "Other")
})

length_ann <- rmfragments %>% dplyr::select(gene_id, pctconsensuscovered) %>% full_join(rmfamilies %>% dplyr::select(gene_id, rte_subfamily)) %>%
    mutate(rte_length_req = ifelse(pctconsensuscovered >= 95, "FL", "Trnc"))
divergence_ann <- rmfragments %>% dplyr::select(gene_id, family, pctdiv, pctconsensuscovered) %>% group_by(family) %>% mutate(family_av_pctdiv = mean(pctdiv, na.rm=TRUE), family_av_coverage = mean(pctconsensuscovered, na.rm=TRUE)) %>% 
    ungroup() %>% 
    mutate(divergence_age = ifelse(
        family_av_pctdiv < 15, "Yng","Old"))

req_annot <- left_join(length_ann, divergence_ann) %>% left_join(intactness_ann) %>% mutate(intactness_req = replace_na(intactness_req, "Other"))
req_annot <- req_annot %>% mutate(req_integrative = case_when(
    (intactness_req == "Intact") & (divergence_age == "Yng") ~ "Yng Intact",
    (intactness_req == "Intact") & (divergence_age == "Old") ~ "Old Intact",
    (pctconsensuscovered >= 95) & (divergence_age == "Yng") ~ "Yng FL",
    (pctconsensuscovered < 95) & (divergence_age == "Yng") ~ "Yng Trnc",
    (pctconsensuscovered >= 95) & (divergence_age == "Old") ~ "Old FL",
    (pctconsensuscovered < 95) & (divergence_age == "Old") ~ "Old Trnc",
    TRUE ~ "Unclassified"
))


#annotate LTR/Int relationship
capture_distance <- 500
ltrs <- rmfamilies %>%
    filter(rte_superfamily == "LTR") %>%
    filter(grepl("LTR", gene_id)) %>%
    filter(!grepl("-int", gene_id)) %>%
    left_join(rmfragments)
ltrsgrs <- GRanges(ltrs)
ltrsgrsextended <- resize(ltrsgrs, width = width(ltrsgrs) + (capture_distance * 2), fix = "center")

ints <- rmfamilies %>%
    filter(rte_superfamily == "LTR") %>%
    filter(grepl("-int$", family)) %>%
    left_join(rmfragments)
intsgrs <- GRanges(ints)
intsfl <- ints %>% filter(pctconsensuscovered >= 95)
intsflgrs <- GRanges(intsfl)
intsnotfl <- ints %>% filter(pctconsensuscovered < 95)
intsnotflgrs <- GRanges(intsnotfl)

Solo_LTR <- ltrsgrsextended %>%
    subsetByOverlaps(intsgrs, invert = TRUE) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Fl_Provirus_5LTR <- ltrsgrs %>%
    flank(capture_distance) %>%
    subsetByOverlaps(intsflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Fl_Provirus_3LTR <- ltrsgrs %>%
    flank(capture_distance, start = FALSE) %>%
    subsetByOverlaps(intsflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Provirus_5LTR <- ltrsgrs %>%
    flank(capture_distance) %>%
    subsetByOverlaps(intsnotflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
Provirus_3LTR <- ltrsgrs %>%
    flank(capture_distance, start = FALSE) %>%
    subsetByOverlaps(intsnotflgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id

NOT_LTR_5Flanked_INT <- intsgrs %>%
    flank(capture_distance) %>%
    subsetByOverlaps(ltrsgrs, invert = TRUE) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id
LTR_5Flanked_INT <- intsgrs %>%
    flank(capture_distance) %>%
    subsetByOverlaps(ltrsgrs) %>%
    mcols() %>%
    as.data.frame() %>%
    tibble() %$% gene_id

ltr_viral_status <- rmfragments %>%
    dplyr::select(gene_id) %>%
    mutate(ltr_viral_status = case_when(
        gene_id %in% Fl_Provirus_5LTR ~ "5'LTR (FL Int)",
        gene_id %in% Fl_Provirus_3LTR ~ "3'LTR (FL Int)",
        gene_id %in% Provirus_5LTR ~ "5'LTR (Trnc Int)",
        gene_id %in% Provirus_3LTR ~ "3'LTR (Trnc Int)",
        gene_id %in% Solo_LTR ~ "LTR (Solo)",
        gene_id %in% NOT_LTR_5Flanked_INT ~ "Int (No 5'LTR)",
        gene_id %in% LTR_5Flanked_INT ~ "Int (Has 5LTR)",
        TRUE ~ "Other"
    ))
ltr_viral_status$ltr_viral_status %>% table()

#ltr status of potentially active, >75% covered ltrs
group_frame <- tibble()
group_frame_any_ltr <- tibble()
element_types_to_scrutinize <- rmfamilies %>% filter(rte_subfamily != "Other") %>% filter(rte_superfamily == "LTR") %$% rte_subfamily %>% unique() %>% grep("int|INT", ., value = TRUE)

for (int in element_types_to_scrutinize) {
    element_family_df <- rmfamilies %>% filter(rte_subfamily == int) %$% family %>% unique() %>% str_split("/", simplify = TRUE) %>% unlist() %>% as.data.frame() %>% tibble()
    element_family_df <- element_family_df %>% tibble() %>% mutate(family_grep_string = paste0(V1, "/",V2, "/"))
    compatible_ltrs <- rmfamilies %>% filter(grepl(paste0(element_family_df %$% family_grep_string, collapse = "|"), family)) %>% filter(grepl("LTR", gene_id)) %$% family %>% unique() 
    int_grs <- rmfamilies %>% filter(rte_subfamily == int) %>% left_join(rmfragments) %>% filter(pctconsensuscovered >= 75) %>% GRanges()
    ltr_grs <- rmfragments %>% filter(family %in% compatible_ltrs) %>% GRanges()
    ltrs_all <- rmfragments %>% filter(grepl("LTR.*", gene_id)) %>% GRanges()
    for (int_id in mcols(int_grs)$gene_id) {
        print(int_id)
        provirus_type <- int_id %>% gsub("-int.*", "", .)
        group_id <- int_id %>% gsub(".*-int", provirus_type, .)
        int_gr <- int_grs[int_grs$gene_id == int_id]
        int_gr_extended <- resize(int_gr, width = width(int_gr) + (500 * 2), fix = "center")
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
    tryCatch({
        inregions$prop <- name_in
    }, error = function(e) {
        print("")
    })
    outregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = TRUE)
    tryCatch({
        outregions$prop <- name_out    
    }, error = function(e) {
        print("")
    })
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

dist_to_nearest_coding_tx <- distanceToNearest(rmfragmentsgr_properinsertloc, coding_transcripts, ignore.strand = TRUE) %>% mcols() %>% as.data.frame() %>% tibble() %$% distance
dist_to_nearest_noncoding_tx <- distanceToNearest(rmfragmentsgr_properinsertloc, noncoding_transcripts, ignore.strand = TRUE) %>% mcols() %>% as.data.frame() %>% tibble() %$% distance
dist_to_nearest_tx <- distanceToNearest(rmfragmentsgr_properinsertloc, transcripts, ignore.strand = TRUE) %>% mcols() %>% as.data.frame() %>% tibble() %$% distance
dist_to_nearest_txs_df <- tibble(gene_id = mcols(rmfragmentsgr_properinsertloc)$gene_id, dist_to_nearest_coding_tx = dist_to_nearest_coding_tx, dist_to_nearest_noncoding_tx = dist_to_nearest_noncoding_tx, dist_to_nearest_tx = dist_to_nearest_tx)



annots <- rmfamilies %>%
    full_join(req_annot) %>%
    full_join(ltr_viral_status) %>%
    full_join(ltr_proviral_groups) %>%
    full_join(region_annot %>% rename_at(vars(-gene_id, -loc_integrative), ~ paste0(., "_loc"))) %>%
    full_join(dist_to_nearest_txs_df)


region_annot %$% gene_id %>% duplicated() %>% sum()
ltr_proviral_groups %$% gene_id %>% duplicated() %>% sum()
req_annot %$% gene_id %>% duplicated() %>% sum()
rmfamilies %$% gene_id %>% duplicated() %>% sum()
annots %$% gene_id %>% duplicated() %>% sum()

rmann <- left_join(rmfragments, annots) 

write_csv(annots, outputs$r_repeatmasker_annotation)
write_csv(rmann, outputs$rmann)
write_csv(rmann %>% filter(refstatus == "NonRef"), outputs$rmann_nonref)

