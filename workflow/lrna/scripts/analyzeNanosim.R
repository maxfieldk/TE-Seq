source("workflow/scripts/defaults.R")
module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined
        ), env = globalenv())
        assign("inputs", list(
            bamgenome = "nanosim/nanosimAlignment_genome.primary.bam",
            bampgenome = "nanosim/nanosimAlignment_genome.perfect.primary.bam"
        ), env = globalenv())
        assign("outputs", list(plot = "nanosim/plots/mapping_accuracy_by_read_length.png"), env = globalenv())
    }
)

rm <- read_csv(params$r_annotation_fragmentsjoined)
rmgrs <- rm %>%
    filter(grepl("L1HS", gene_id)) %>%
    GRanges()


for (type in c("Perfect", "ONT_dRNA")) {
    if (type == "Perfect") {
        bamfile <- BamFile(inputs$bampgenome)
    } else {
        bamfile <- BamFile(inputs$bamgenome)
    }
    aln1 <- scanBam(bamfile)
    aln <- as.data.frame(aln1) %>% tibble()
    alndf <- aln %>%
        mutate(qsplitname = str_split(qname, "_", simplify = TRUE)[, 1]) %>%
        mutate(qsplitname = str_replace(qsplitname, "-", "_")) %>%
        mutate(qsplitname = str_replace(qsplitname, "-", ".")) %>%
        mutate(qsplitname = str_replace(qsplitname, "-", "_")) %>%
        mutate(qsplitname = ifelse(str_count(qsplitname, "_") + str_count(qsplitname, "\\.") > 2, qsplitname, str_replace(qsplitname, "\\.", "_"))) %>%
        mutate(start = pos, end = pos + qwidth) %>%
        mutate(seqnames = rname)

    p <- alndf %>% ggplot() +
        geom_histogram(aes(x = qwidth), binwidth = 100) +
        mtopen
    mysave(sprintf("nanosim/plots/read_length_distribution_%s.png", type), 4, 4)

    alndf1 <- alndf %>%
        select(-c(mrnm, mpos)) %>%
        mutate(readlengthgroup = ifelse(qwidth > 6000, "6kb", ifelse(qwidth > 5000, "5kb", ifelse(qwidth > 4000, "4kb", ifelse(qwidth > 3000, "3kb", ifelse(qwidth > 2000, "2kb", ifelse(qwidth > 1000, "1kb", ifelse(qwidth > 500, "0.5kb", "<0.5kb")))))))) %>%
        drop_na()

    grs <- GRanges(alndf1)
    mbo <- mergeByOverlaps(grs, rmgrs)
    pf <- mbo %>%
        as.data.frame() %>%
        tibble() %>%
        select(gene_id, qsplitname, readlengthgroup) %>%
        mutate(match = ifelse(gene_id == qsplitname, "Proper Alignment", "Misalignment")) %>%
        group_by(readlengthgroup, match) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        pivot_wider(names_from = match, values_from = n) %>%
        mutate(AlignmentAccuracy = `Proper Alignment` / (Misalignment + `Proper Alignment`))
    p <- pf %>% ggplot(aes(x = readlengthgroup, y = AlignmentAccuracy)) +
        geom_col() +
        mtopen +
        ggtitle(sprintf("Simulated %s L1HS Intact Reads", type)) +
        labs(y = "Mapping Accuracy", x = "Read Length") +
        anchorbar
    mysave(sprintf("nanosim/plots/mapping_accuracy_by_read_length_%s.png", type), 5, 5)

    p <- aln %>%
        select(qname, rname, qwidth, mapq) %>%
        drop_na() %>%
        ggplot() +
        geom_histogram(aes(x = mapq)) +
        mtopen +
        ggtitle(sprintf("Simulated %s L1HS Intact 6kb Reads", type)) +
        labs(y = "Count", x = "Mapq")
    mysave(sprintf("nanosim/plots/mapqhistogram_%s.png", type), 4, 4)
}



bamfileperfect <- BamFile(inputs$bampgenome)
bamfile <- BamFile(inputs$bamgenome)
aln1perfect <- scanBam(bamfileperfect)
aln1 <- scanBam(bamfile)
alnperfect <- as.data.frame(aln1perfect) %>% tibble()
aln <- as.data.frame(aln1) %>% tibble()
alnperfect$readtype <- "Perfect"
aln$readtype <- "ONT_dRNA"
aln2 <- bind_rows(aln, alnperfect)
alndf <- aln2 %>%
    mutate(qsplitname = str_split(qname, "_", simplify = TRUE)[, 1]) %>%
    mutate(qsplitname = str_replace(qsplitname, "-", "_")) %>%
    mutate(qsplitname = str_replace(qsplitname, "-", ".")) %>%
    mutate(qsplitname = str_replace(qsplitname, "-", "_")) %>%
    mutate(qsplitname = ifelse(str_count(qsplitname, "_") + str_count(qsplitname, "\\.") > 2, qsplitname, str_replace(qsplitname, "\\.", "_"))) %>%
    mutate(start = pos, end = pos + qwidth) %>%
    mutate(seqnames = rname)

alndf1 <- alndf %>%
    select(-c(mrnm, mpos)) %>%
    mutate(readlengthgroup = ifelse(qwidth > 6000, "6kb", ifelse(qwidth > 5000, "5kb", ifelse(qwidth > 4000, "4kb", ifelse(qwidth > 3000, "3kb", ifelse(qwidth > 2000, "2kb", ifelse(qwidth > 1000, "1kb", ifelse(qwidth > 500, "0.5kb", "<0.5kb")))))))) %>%
    drop_na()

grs <- GRanges(alndf1)
mbo <- mergeByOverlaps(grs, rmgrs)
pf <- mbo %>%
    as.data.frame() %>%
    tibble() %>%
    select(gene_id, qsplitname, readlengthgroup, readtype) %>%
    mutate(match = ifelse(gene_id == qsplitname, "Proper Alignment", "Misalignment")) %>%
    group_by(readtype, readlengthgroup, match) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = match, values_from = n) %>%
    mutate(AlignmentAccuracy = `Proper Alignment` / (Misalignment + `Proper Alignment`)) %>%
    dplyr::rename(`Error Profile` = readtype)
p <- pf %>% ggplot(aes(x = readlengthgroup, y = AlignmentAccuracy, fill = `Error Profile`)) +
    geom_col(pos = "dodge") +
    mtopen + scale_contrasts +
    ggtitle(sprintf("Simulated L1HS Intact Reads")) +
    labs(y = "Mapping Accuracy", x = "Read Length") +
    anchorbar
mysave(sprintf(outputs$plot), 6, 5)

pf <- mbo %>%
    as.data.frame() %>%
    tibble() %>%
    select(gene_id, qsplitname, readlengthgroup, readtype, mapq) %>%
    mutate(mapqGroup = ifelse(mapq > 0, "mapq > 0", "mapq = 0")) %>%
    mutate(match = ifelse(gene_id == qsplitname, "Proper Alignment", "Misalignment")) %>%
    group_by(readtype, readlengthgroup, match, mapqGroup) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = match, values_from = n) %>%
    replace_na(list(`Misalignment` = 0)) %>%
    mutate(AlignmentAccuracy = `Proper Alignment` / (Misalignment + `Proper Alignment`)) %>%
    dplyr::rename(`Error Profile` = readtype)


p <- pf %>% ggplot(aes(x = readlengthgroup, y = AlignmentAccuracy, fill = `Error Profile`)) +
    geom_col(pos = "dodge") +
    facet_wrap(~mapqGroup) +
    mtopen + scale_contrasts +
    ggtitle(sprintf("Simulated L1HS Intact Reads")) +
    labs(y = "Mapping Accuracy", x = "Read Length") +
    anchorbar
mysave(sprintf("nanosim/plots/mapping_accuracy_by_read_length_mapqfacet.png"), 8, 5)
