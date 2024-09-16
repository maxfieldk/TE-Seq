module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
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

sample <- "sen1"
rm <- read_csv(params$r_annotation_fragmentsjoined)
rmgrs <- rm %>%
    filter(grepl("L1HS|L1PA", gene_id)) %>%
    GRanges()

maskedreadsgr <- import(inputs$gtf)
maskedreadsdf <- as.data.frame(maskedreadsgr) %>% tibble()


maskedreadsdf %>%
    filter(str_detect(gene_id, "L1HS")) %>%
    filter(width > 1000) %>%
    select(-seqnames)

maskedreadsdffilterapplied <- maskedreadsdf %>%
    filter(str_detect(gene_id, "L1HS|L1PA")) %>%
    filter(width > 3000)

seqnamesforfilter <- maskedreadsdffilterapplied %>%
    select(seqnames)

write_delim(seqnamesforfilter, sprintf("intermediates/%s/repeatmasker/seqnamesforfilter.txt", sample), delim = "\n", col_names = FALSE)

system(sprintf("samtools view -F 256 -N intermediates/%s/repeatmasker/seqnamesforfilter.txt -o /users/mkelsey/data/Nanopore/dRNALF1/intermediates/%s/alignments/genome/%s.sorted.filtered_for_l1_rm.bam /users/mkelsey/data/Nanopore/dRNALF1/intermediates/%s/alignments/genome/%s.sorted.bam", sample, sample, sample, sample, sample))

bamFile <- BamFile("/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen1/alignments/genome/sen1.sorted.filtered_for_l1_rm.bam")
aln1 <- scanBam(bamFile)
aln <- as.data.frame(aln1) %>% tibble()
alndf <- aln %>%
    mutate(start = pos, end = pos + qwidth) %>%
    mutate(seqnames = rname)

p <- alndf %>% ggplot() +
    geom_histogram(aes(x = qwidth), binwidth = 100) +
    mtopen
mysave(sprintf("maskedread_length_distribution.png"), 4, 4)

alndf1 <- alndf %>%
    select(-c(mrnm, mpos)) %>%
    mutate(readlengthgroup = ifelse(qwidth > 6000, "6kb", ifelse(qwidth > 5000, "5kb", ifelse(qwidth > 4000, "4kb", ifelse(qwidth > 3000, "3kb", ifelse(qwidth > 2000, "2kb", ifelse(qwidth > 1000, "1kb", ifelse(qwidth > 500, "0.5kb", "<0.5kb")))))))) %>%
    drop_na()

alndf1 %>%
    filter(!is.na(gene_id)) %>%
    select(gene_id, qname, qwidth, readlengthgroup) %>%
    mutate(match = ifelse(gene_id == qname, "Proper Alignment", "Misalignment")) %>%
    write_delim("maskedread_alignment_summary.txt", delim = "\t")

maskedreadsforjoin <- maskedreadsdffilterapplied %>%
    dplyr::select(seqnames, width, pctdiv, gene_id) %>%
    dplyr::rename(qname = seqnames, masked_element_length = width)
maskedreadaln <- GRanges(alndf1 %>% left_join(maskedreadsforjoin))


mbo <- mergeByOverlaps(maskedreadaln, rmgrs)
mbodf <- mbo %>%
    as.data.frame() %>%
    tibble()
colnames(mbodf)
mbodf %>% select(qname, width, transcript_id)

mbodf %>%
    dplyr::select(masked_element_length, maskedreadaln.qwidth, maskedreadaln.gene_id, rmgrs.gene_id) %>%
    print(n = 100)

interestingpositions <- mbodf %>%
    select(maskedreadaln.seqnames, maskedreadaln.start, maskedreadaln.end, masked_element_length, maskedreadaln.qwidth, maskedreadaln.gene_id, rmgrs.gene_id)
write_delim(interestingpositions, sprintf("intermediates/%s/repeatmasker/masked_l1_coordinates.txt", sample), delim = "\t", col_names = TRUE)
