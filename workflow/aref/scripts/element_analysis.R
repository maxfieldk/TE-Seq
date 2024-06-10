module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

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
library(msa)
library(seqinr)
library(ape)
library(ggtree)
library(treeio)
library(ggmsa)
library(gggenes)
library(ggnewscale)
library(RColorBrewer)
library(ORFik)


tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/AD1_annotations/AD1_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/AD1_annotations/AD1_repeatmasker_annotation.csv",
            ref = "aref/AD1.fa"
        ), env = globalenv())
        assign("params", list(l13 = conf$l13fasta), env = globalenv())
        assign("outputs", list(
            outfile = "aref/AD1_Analysis/l1element_analysis.outfile",
            plots = "aref/AD1_Analysis/plots.Rda"
        ), env = globalenv())
    }
)

fa <- FaFile(inputs$ref)
rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies) %>%
    filter(refstatus != "NonCentral")
outputdir <- dirname(outputs$outfile)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

rmann %>% filter(l1_intactness_req == "L1HS ORFs Intact") %>% head() %>% write_delim(sprintf("%s/rmann_head.csv", outputdir), delim = "\t", col_names = TRUE)
rmann %$% l1_intactness_req %>% table()

options(scipen = 999)


rmanngr <- GRanges(rmann)
l1hsgr <- rmanngr[grepl("L1HS", rmanngr$gene_id)]
l1pa2gr <- rmanngr[grepl("L1PA2", rmanngr$gene_id)]
l1pa3gr <- rmanngr[grepl("L1PA3", rmanngr$gene_id)]

hs_pa2_pa3_gr <- c(c(l1hsgr, l1pa2gr), l1pa3gr)
hs_pa2_pa3 <- as.data.frame(hs_pa2_pa3_gr) %>% tibble()

# plot number of l1 sequences
pf <- rmann %>%
    filter(refstatus != "NonCentral") %>%
    filter(grepl("L1PA|L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily, refstatus) %>%
    summarise(n = n()) %>%
    ungroup()
p <- pf %>% ggplot() +
    geom_bar(aes(x = rte_subfamily, y = n), color = "black", stat = "identity") +
    labs(title = "Number of L1 Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
    facet_grid(~refstatus, scales = "free" , space = "free_x") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

mysaveandstore(sprintf("%s/l1familycount.pdf", outputdir), w = 8, h = 5)

p <- rmann %>%
    filter(grepl("L1PA|L1HS", rte_subfamily)) %>%
    filter(length > 5999) %>%
    group_by(rte_subfamily, refstatus, l1_intactness_req) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    ggplot() +
    geom_bar(aes(x = rte_subfamily, y = n, fill = l1_intactness_req), color = "black", stat = "identity") +
    labs(title = "Number of Full Length L1 Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
    facet_grid(~refstatus, scales = "free" , space = "free_x") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

mysaveandstore(sprintf("%s/l1FLfamilycount.pdf", outputdir), w = 8, h = 5)


gene_ids <- hs_pa2_pa3_gr$gene_id
hs_pa2_pa3_ss <- getSeq(fa, hs_pa2_pa3_gr)
mcols(hs_pa2_pa3_ss) <- mcols(hs_pa2_pa3_gr)
names(hs_pa2_pa3_ss) <- mcols(hs_pa2_pa3_gr)$gene_id
hs_pa2_pa3_fl_ss <- hs_pa2_pa3_ss[width(hs_pa2_pa3_ss) > 5999]

hs_pa2_pa3_intact_ss <- hs_pa2_pa3_fl_ss[grepl("Intact", mcols(hs_pa2_pa3_fl_ss)$l1_intactness_req)]

# last_n <- 20
# A_freq <- letterFrequency(subseq(flss, -last_n, -1), "A") / last_n
# png("l1hsAfreq.png", width = 6, height = 4, units = "in", res = 300)
# hist(A_freq)
# dev.off()

p <- rmann %>%
    filter(grepl("Intact", l1_intactness_req)) %>%
    ggplot() +
    geom_bar(aes(x = rte_subfamily, fill = refstatus), color = "black") +
    mtopen +
    anchorbar
mysaveandstore(sprintf("%s/l1_intact_refstatus.pdf", outputdir), 4, 4)

# ORF1 AND ORF2 SEQUENCE ANALYSES
# cannonical locations of RT and EN in ORF2 protein
RTcoords <- c(498, 773)
ENcoords <- c(1, 239)
RTcoordsnt <- 3 * RTcoords
ENcoordsnt <- 3 * ENcoords

pos <- ORFik::findORFs(hs_pa2_pa3_intact_ss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
names(pos) <- names(hs_pa2_pa3_intact_ss)
orf1 <- list()
orf2 <- list()
element_ann_df <- tibble()
for (element in names(hs_pa2_pa3_intact_ss)) {
    print(element)
    firstorf <- pos[[element]][1]
    secondorf <- pos[[element]][2]
    firstorfrow <- firstorf %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(start, end) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "ORF1")
    secondorfrow <- secondorf %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(start, end) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "ORF2")
    enrow <- tibble(start = secondorfrow$start, end = secondorfrow$start + ENcoordsnt[2]) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "EN")
    rtrow <- tibble(start = secondorfrow$start + RTcoordsnt[1], end = secondorfrow$start + RTcoordsnt[2]) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "RT")
    X5UTRrow <- tibble(start = 1, end = firstorfrow$start - 1) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "5UTR")
    X3UTRrow <- tibble(start = secondorfrow$end + 1, end = length(hs_pa2_pa3_intact_ss[[element]])) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "3UTR")

    element_ann_df <- rbind(element_ann_df, X5UTRrow, firstorfrow, enrow, rtrow, secondorfrow, X3UTRrow)
    orf1[[element]] <- hs_pa2_pa3_intact_ss[[element]][firstorf]
    orf2[[element]] <- hs_pa2_pa3_intact_ss[[element]][secondorf]
}

write_delim(element_ann_df %>% dplyr::select(gene_id, start, end, feature), sprintf("%s/intact_l1_anatomy_coordinates.tsv", outputdir), delim = "\t", col_names = TRUE)

orf1ss <- DNAStringSet(orf1)
orf2ss <- DNAStringSet(orf2)

orf1ps <- Biostrings::translate(orf1ss)
orf2ps <- Biostrings::translate(orf2ss)


ENps <- subseq(orf2ps, start = ENcoords[1], end = ENcoords[2])
RTps <- subseq(orf2ps, start = RTcoords[1], end = RTcoords[2])

writeXStringSet(orf1ps, sprintf("%s/orf1ps.fa", outputdir))
system(sprintf("echo $(pwd); mafft --auto %s/orf1ps.fa > %s/orf1ps.aln.fa", outputdir, outputdir))
aln <- read.alignment(sprintf("%s/orf1ps.aln.fa", outputdir), format = "fasta")
d <- dist.alignment(aln, "identity")
orf1Tree <- nj(d)
png(paste0(sprintf("%s/orf1Tree.png", outputdir)), width = 10, height = 30, units = "in", res = 300)
plot(orf1Tree, main = "Phylogenetic Tree of ORF1 Sequences")
dev.off()


######## TIME TO BLAST!
l1.3 <- readDNAStringSet(params$l13)
names(l1.3) <- "l1.3"
l1.3.orfs <- findORFs(l1.3, startCodon = "ATG", minimumLength = 300)
l1.3.orf1 <- DNAStringSet(l1.3[[1]][l1.3.orfs[[1]][1]])
names(l1.3.orf1) <- "orf1"
l1.3.orf2 <- DNAStringSet(l1.3[[1]][l1.3.orfs[[1]][2]])
names(l1.3.orf2) <- "orf2"


l1.3.seqs <- c(l1.3, l1.3.orf1, l1.3.orf2)

writeXStringSet(l1.3.seqs, sprintf("%s/l1.3.seqs.fa", outputdir))

system(sprintf("mkdir -p %s/blastdb; cd %s/blastdb; makeblastdb -in ../l1.3.seqs.fa -dbtype nucl -out l1.3.seqs", outputdir, outputdir))

## load a BLAST database (replace db with the location + name of the BLAST DB
## without the extension)

# will fail if you don't have any non-ref elements
library(ggbio)

bl <- blast(db = sprintf("%s/blastdb/l1.3.seqs", outputdir))
bres <- tibble(predict(bl, hs_pa2_pa3_ss)) %>% left_join(rmfragments, by = c("QueryID" = "gene_id"))

tryCatch({
nonref_aln_l13 <- bres %>%
    filter(SubjectID == "l1.3") %>%
    filter(refstatus == "NonRef") %>%
    rowwise() %>%
    mutate(minIget = min(S.start, S.end)) %>%
    mutate(maxIget = max(S.start, S.end))
# needed but cannot be loaded with orfik
{
    notsplitnonrefl1hs <- nonref_aln_l13 %>%
        filter(SubjectID == "l1.3") %>%
        filter(family == "LINE/L1/L1HS") %>%
        group_by(QueryID) %>%
        mutate(n = n()) %>%
        filter(n < 2)
    notsplitnonrefl1hsgrs <- GRanges(seqnames = "l1.3", ranges = IRanges(start = notsplitnonrefl1hs$minIget, end = notsplitnonrefl1hs$maxIget))
    p <- autoplot(notsplitnonrefl1hsgrs, aes(color = end)) + mtopen + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "L1HS Blast to L1.3")
    mysaveandstore(sprintf("%s/l13alnnotsplit.pdf", outputdir))
}


}, error = function(e) {
    print("fails if you don't have non-ref elements")
})

p <- bres %>%
    filter(SubjectID == "l1.3") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mtopen +
    anchorbar
mysaveandstore(sprintf("%s/l13identHist.pdf", outputdir))

p <- bres %>%
    filter(SubjectID == "orf1") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mtopen +
    anchorbar
mysaveandstore(sprintf("%s/orf1identHist.pdf", outputdir))

p <- bres %>%
    filter(SubjectID == "orf2") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mtopen +
    anchorbar
mysaveandstore(sprintf("%s/orf2identHist.pdf", outputdir))


# trying to make gene diagrams
# need to remove any (in this case it was just 1, full length elements which have multiple alignments to L1.3. I keep the longest one.
# brestidy <- bres %>%
#     group_by(QueryID, SubjectID) %>%
#     filter(Alignment.Length == max(Alignment.Length)) %>%
#     ungroup() %>%
#     pivot_wider(id_cols = QueryID, names_from = SubjectID, values_from = c(Q.start, Q.end))

# p <- ggplot(brestidy[1:5, ], aes(xmin = Q.start_l1.3, xmax = Q.end_l1.3, y = QueryID), fill = "white") +
#     geom_gene_arrow() +
#     geom_subgene_arrow(aes(xsubmin = Q.start_orf1, xsubmax = Q.end_orf1, fill = "orf1")) +
#     geom_subgene_arrow(aes(xsubmin = Q.start_orf2, xsubmax = Q.end_orf2, fill = "orf2")) +
#     facet_wrap(~QueryID, scales = "free" , space = "free_x", ncol = 1) +
#     scale_fill_brewer(palette = "Set3") +
#     theme_genes()
# mysaveandstore("geneDiagram", "blastRes")




##### ggtree
writeXStringSet(hs_pa2_pa3_intact_ss, sprintf("%s/hs_pa2_pa3_intact.fa", outputdir))
l1hs_intact <- hs_pa2_pa3_intact_ss[grepl("L1HS", names(hs_pa2_pa3_intact_ss))]
writeXStringSet(l1hs_intact, sprintf("%s/l1hs_intact.fa", outputdir))
l1pa2_intact <- hs_pa2_pa3_intact_ss[grepl("L1PA2", names(hs_pa2_pa3_intact_ss))]
writeXStringSet(l1pa2_intact, sprintf("%s/l1pa2_intact.fa", outputdir))
l1pa3_intact <- hs_pa2_pa3_intact_ss[grepl("L1PA3", names(hs_pa2_pa3_intact_ss))]
writeXStringSet(l1pa3_intact, sprintf("%s/l1pa3_intact.fa", outputdir))
system(sprintf("echo $(pwd); mafft --auto %s/hs_pa2_pa3_intact.fa > %s/hs_pa2_pa3_intact.aln.fa", outputdir, outputdir))
system(sprintf("echo $(pwd); mafft --auto %s/l1hs_intact.fa > %s/l1hs_intact.aln.fa", outputdir, outputdir))

l1hs_intact_aln <- read.alignment(sprintf("%s/l1hs_intact.aln.fa", outputdir), format = "fasta")
d <- dist.alignment(l1hs_intact_aln, "identity", gap = FALSE)
d_gapped <- dist.alignment(l1hs_intact_aln, "identity", gap = TRUE)

# this is the pct difference in identity between sequences
dsqaured <- 100 * d**2
d_gappedsquared <- 100 * d_gapped**2

tree_raw <- as.treedata(fastme.bal(dsqaured))
treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
# to df for further modifications
ttibble_raw <- as_tibble(tree_raw)
ttibblegapped_raw <- as_tibble(treegapped_raw)

ttibble_raw <- ttibble_raw %>%
    left_join(as.data.frame(mcols(hs_pa2_pa3_intact_ss)), by = c("label" = "gene_id")) %>%
    mutate(group = str_extract(label, "L1HS|L1PA2")) %>%
    mutate(group = factor(group, levels = c("L1HS", "L1PA2")))

ttibblegapped_raw <- ttibblegapped_raw %>%
    left_join(as.data.frame(mcols(hs_pa2_pa3_intact_ss)), by = c("label" = "gene_id")) %>%
    mutate(group = str_extract(label, "L1HS|L1PA2")) %>%
    mutate(group = factor(group, levels = c("L1HS", "L1PA2")))

l1hs_intact_tree <- as.treedata(ttibble_raw)
l1hs_intact_tree_gapped <- as.treedata(ttibblegapped_raw)

p <- ggtree(l1hs_intact_tree, layout = "fan") +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysaveandstore(sprintf("%s/l1hs_intact_tree_fan.pdf", outputdir), w = 7, h = 7)

p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysaveandstore(sprintf("%s/l1hs_intact_tree.pdf", outputdir), w = 7.5, h = 20)
# make the tree horizontal
p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label, angle = 90)) +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny") +
    coord_flip()
mysaveandstore(sprintf("%s/l1hs_intact_tree_horizontal.pdf", outputdir), w = 20, h = 8)


p <- ggtree(l1hs_intact_tree_gapped) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")

mysaveandstore(sprintf("%s/l1hs_intact_tree_gapped.pdf", outputdir), w = 7.5, h = 20)

# tree with msa
l1hs_intact_aln <- readDNAMultipleAlignment(sprintf("%s/l1hs_intact.aln.fa", outputdir), format = "fasta")
consensus <- consensusString(l1hs_intact_aln)
l1hs_intact_aln_ss <- as(l1hs_intact_aln, "DNAStringSet")
temp <- c(l1hs_intact_aln_ss %>% as.character(), consensus %>% as.character() %>% setNames("consensus"))
l1hs_intact_aln_c <- DNAMultipleAlignment(temp)
l1hs_intact_aln_ss_c <- as(l1hs_intact_aln_c, "DNAStringSet")
writeXStringSet(l1hs_intact_aln_ss_c, file = sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir))
p <- simplot(sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
    ggtitle("L1HS ORF1&2 Intact Consensus Accuracy") +
    mtopen
mysaveandstore(sprintf("%s/l1hs_intact_alnWconensus_similarity.pdf", outputdir), w = 10, h = 5)

p <- simplot(sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
    ggtitle("L1HS ORF1&2 Intact Consensus Accuracy") +
    geom_hline(yintercept = 0.95, col = "red", lty = 2) +
    mtopen
mysaveandstore(sprintf("%s/l1hs_intact_alnWconensus_similarity_ONTdRNAerror.pdf", outputdir), w = 10, h = 5)


data <- tidy_msa(l1hs_intact_aln_ss, 1, 100)

p <- ggtree(tr = l1hs_intact_tree_gapped) +
    geom_treescale(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    ggtitle("L1HS and L1PA2 ORF1&2 Intact Phylogeny") +
    geom_facet(geom = geom_msa, data = data, panel = "msa", font = NULL)
mysaveandstore(sprintf("%s/l1hs_intact_msa_first_100nt_treemsa.pdf", outputdir), w = 10, h = 20)

####
save(mysaveandstoreplots, file = outputs$plots)

x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)
