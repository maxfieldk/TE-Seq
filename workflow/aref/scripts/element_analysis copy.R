source("~/data/common/myDefaults.r")
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

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")[["aref"]]
)

tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/annotations/repeatmasker_annotation.csv",
            ref = "aref/ref.fa"
        ), env = globalenv())
        assign("params", list(l13 = conf$l13fasta), env = globalenv())
        assign("outputs", list(outfile = "aref/RefAnalysis/element_analysis.outfile"), env = globalenv())
    }
)

fa <- FaFile(inputs$ref)
rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)
outputdir <- dirname(outputs$outfile)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)




options(scipen = 999)


gtf <- GRanges(rmann)
l1hsgtf <- gtf[grepl("L1HS", gtf$gene_id)]
l1pa2gtf <- gtf[grepl("L1PA2", gtf$gene_id)]
l1pa3gtf <- gtf[grepl("L1PA3", gtf$gene_id)]

l1hspa23gtf <- c(c(l1hsgtf, l1pa2gtf), l1pa3gtf)
l1hspa23df <- as.data.frame(l1hspa23gtf) %>% tibble()

# plot number of l1 sequences
pf <- rmann %>%
    filter(grepl("L1PA|L1HS", rte_subfamily)) %>%
    group_by(rte_subfamily, refstatus) %>%
    summarise(n = n()) %>%
    ungroup()
p <- pf %>% ggplot() +
    geom_bar(aes(x = rte_subfamily, y = n), color = "black", stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of L1 Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
    facet_wrap(~refstatus, scales = "free") +
    mytheme +
    anchorbar
mysave(sprintf("%s/l1familycount.png", outputdir), w = 8, h = 5)

p <- rmann %>%
    filter(grepl("L1PA|L1HS", rte_subfamily)) %>%
    filter(length > 5999) %>%
    group_by(rte_subfamily, refstatus, l1_intactness_req) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    ggplot() +
    geom_bar(aes(x = rte_subfamily, y = n, fill = l1_intactness_req), color = "black", stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of Full Length L1 Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
    facet_wrap(~refstatus, scales = "free") +
    mytheme +
    anchorbar
mysave(sprintf("%s/l1FLfamilycount.png", outputdir), w = 8, h = 5)


gene_ids <- l1hspa23gtf$gene_id
l1hspa23ss <- getSeq(fa, l1hspa23gtf)
mcols(l1hspa23ss) <- mcols(l1hspa23gtf)
names(l1hspa23ss) <- mcols(l1hspa23ss)$gene_id
flss <- l1hspa23ss[width(l1hspa23ss) > 5999]


# last_n <- 20
# A_freq <- letterFrequency(subseq(flss, -last_n, -1), "A") / last_n
# png("l1hsAfreq.png", width = 6, height = 4, units = "in", res = 300)
# hist(A_freq)
# dev.off()


# orf analysis
pos <- ORFik::findORFs(flss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
str(pos)
names(pos) <- names(flss[as.integer(names(pos))])
tpm <- width(pos)


# toplot <- unlist(tpm)
# myvals <- seq(1010, 1027, by = 1)
# p <- hist(toplot[toplot %in% myvals], breaks = myvals, warn.unused = FALSE)
# mysave(sprintf("%s/orflength_distribution_orf1adjacent.png", outputdir))

# png("orflength_distribution_orf2adjacent.png", width = 6, height = 4, units = "in", res = 300)
# toplot <- unlist(tpm)
# myvals <- seq(3828 - 10, 3828 + 10, by = 1)
# hist(toplot[toplot %in% myvals], breaks = myvals, warn.unused = FALSE)
# dev.off()
# so we don't seem to be filtering out many L1s on the basis of their having one or so AAs inserted or deleted in ORF1 or ORF2


# FILTERING ON THE BASIS OF ORF LENGTHS
pass <- c()
for (element in names(pos)) {
    if (3828 %in% width(pos[[element]]) & 1017 %in% width(pos[[element]])) {
        pass <- c(pass, element)
    }
}
gtfpass <- l1hspa23gtf[l1hspa23gtf$gene_id %in% pass]
flsspass <- flss[pass]

pos[["L1HS_17p12_2"]]
pos[["L1PA3_22p13_1"]]
pos[["L1PA3_15p13_2"]]


gtfpassdf <- as.data.frame(gtfpass) %>% tibble()
gtfpassdf %$% refstatus %>% table()
gtfpassdf %>%
    filter(refstatus == "NonRef") %$% family %>%
    table()
p <- gtfpassdf %>% ggplot() +
    geom_bar(aes(x = family, fill = refstatus))
mysave(sprintf("%s/l1_intact_refstatus.png", outputdir))
###### WRITE PASS ELEMENTS IN VARIOUS FILE FORMATS
# gtf
export(gtfpass, sprintf("%s/l1hspa23_orfsintact.gtf", outputdir))
# dataframe
element_annotation <- tibble(gene_id = pass, intactness = "ORFs_Intact")

# location genic vs nongenic for all elements
# centromeric vs non centromeriuc for all elements

writeXStringSet(flsspass, sprintf("%s/l1hspa23intact.fa", outputdir))


# A_freq <- letterFrequency(subseq(flsspass, -last_n, -1), "A") / last_n
# png("l1hsAfreq_pass.png", width = 6, height = 4, units = "in", res = 300)
# hist(A_freq)
# dev.off()



# ORF1 AND ORF2 SEQUENCE ANALYSES
# cannonical locations of RT and EN in ORF2 protein
RTcoords <- c(498, 773)
ENcoords <- c(1, 239)
RTcoordsnt <- 3 * RTcoords
ENcoordsnt <- 3 * ENcoords


orf1 <- list()
orf2 <- list()
element_ann_df <- tibble()
for (element in pass) {
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
    X3UTRrow <- tibble(start = secondorfrow$end + 1, end = length(flsspass[[element]])) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "3UTR")

    element_ann_df <- rbind(element_ann_df, X5UTRrow, firstorfrow, enrow, rtrow, secondorfrow, X3UTRrow)
    orf1[[element]] <- flss[[element]][firstorf]
    orf2[[element]] <- flss[[element]][secondorf]
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
outputdir
system(sprintf("mkdir -p %s/blastdb; cd %s/blastdb; makeblastdb -in ../l1.3.seqs.fa -dbtype nucl -out l1.3.seqs", outputdir, outputdir))

## load a BLAST database (replace db with the location + name of the BLAST DB
## without the extension)

bl <- blast(db = sprintf("%s/blastdb/l1.3.seqs", outputdir))
short <- l1hspa23ss[width(l1hspa23ss) < 3000]
long <- l1hspa23ss[width(l1hspa23ss) > 5999]
bres <- tibble(predict(bl, l1hspa23ss)) %>% left_join(rmfragments, by = c("QueryID" = "gene_id"))

l1.3aln <- bres %>%
    filter(SubjectID == "l1.3") %>%
    filter(refstatus == "NonRef") %>%
    rowwise() %>%
    mutate(minIget = min(S.start, S.end)) %>%
    mutate(maxIget = max(S.start, S.end))

l13aln <- bres %>%
    filter(SubjectID == "l1.3") %>%
    rowwise() %>%
    mutate(minIget = min(S.start, S.end)) %>%
    mutate(maxIget = max(S.start, S.end))



grs <- GRanges(seqnames = "l1.3", ranges = IRanges(start = l1.3aln$minIget, end = l1.3aln$maxIget, name = l1.3aln$QueryID))

library(ggbio)
# needed but cannot be loaded with orfik
{
    notsplitnonrefl1hs <- bres %>%
        filter(SubjectID == "l1.3") %>%
        filter(refstatus == "NonRef") %>%
        filter(family == "LINE/L1/L1HS") %>%
        group_by(QueryID) %>%
        mutate(n = n()) %>%
        filter(n < 2)

    # non full length cov
    nflgrs <- grs[width(grs) < 6000]
    p <- autoplot(coverage(nflgrs)$l1.3) + mytheme + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "L1HS(>6000 bp) Blast to L1.3")
    mysave(sprintf("%s/l13covnotfulllength.png", outputdir))
    p <- autoplot(coverage(grs)$l1.3) + mytheme + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "NonRef L1HS Blast to L1.3")
    mysave(sprintf("%s/nonref_l13cov.png", outputdir))
    p <- autoplot(grs, aes(color = width)) + mytheme + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "L1HS Blast to L1.3")
    mysave(sprintf("%s/l13aln.png", outputdir))
    p <- autoplot(nflgrs, aes(color = width)) + mytheme + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "L1HS Blast to L1.3")
    mysave(sprintf("%s/l13alnnotfulllength.png", outputdir))
    # filtering out any split elements
    notsplit <- bres %>%
        filter(SubjectID == "l1.3") %>%
        group_by(QueryID) %>%
        summarise(n = n()) %>%
        filter(n < 2) %>%
        pull(QueryID)
    l1.3alnns <- bres %>%
        filter(SubjectID == "l1.3") %>%
        filter(QueryID %in% notsplit) %>%
        rowwise() %>%
        mutate(minIget = min(S.start, S.end)) %>%
        mutate(maxIget = max(S.start, S.end))
    nsgrs <- GRanges(seqnames = "l1.3", ranges = IRanges(start = l1.3alnns$minIget, end = l1.3alnns$maxIget))
    nsnflgrs <- nsgrs[width(nsgrs) < 6000]
    p <- autoplot(nsnflgrs, aes(color = end)) + mytheme + labs(x = "Position in L1.3 (bp)", y = "Coverage", title = "L1HS Blast to L1.3")
    mysave(sprintf("%s/l13alnnotsplitnotfulllength.png", outputdir))
}

bres <- tibble(predict(bl, flsspass))
p <- bres %>%
    filter(SubjectID == "l1.3") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mytheme
mysave(sprintf("%s/l13identHist.png", outputdir))

p <- bres %>%
    filter(SubjectID == "orf1") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mytheme
mysave(sprintf("%s/orf1identHist.png", outputdir))

p <- bres %>%
    filter(SubjectID == "orf2") %>%
    ggplot() +
    geom_histogram(aes(x = Perc.Ident)) +
    mytheme
mysave(sprintf("%s/orf2identHist.png", outputdir))


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
#     facet_wrap(~QueryID, scales = "free", ncol = 1) +
#     scale_fill_brewer(palette = "Set3") +
#     theme_genes()
# mysave("geneDiagram", "blastRes")




##### ggtree
writeXStringSet(flsspass, sprintf("%s/flsspass.fa", outputdir))
l1hs_intact <- flsspass[grepl("L1HS", names(flsspass))]
writeXStringSet(l1hs_intact, sprintf("%s/l1hs_intact.fa", outputdir))
l1pa2_intact <- flsspass[grepl("L1PA2", names(flsspass))]
writeXStringSet(l1pa2_intact, sprintf("%s/l1pa2_intact.fa", outputdir))
l1pa3_intact <- flsspass[grepl("L1PA3", names(flsspass))]
writeXStringSet(l1pa3_intact, sprintf("%s/l1pa3_intact.fa", outputdir))
system(sprintf("echo $(pwd); mafft --auto %s/flsspass.fa > %s/flss.aln.fa", outputdir, outputdir))
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
    left_join(gtfpassdf, by = c("label" = "gene_id")) %>%
    mutate(group = str_extract(label, "L1HS|L1PA2")) %>%
    mutate(group = factor(group, levels = c("L1HS", "L1PA2")))

ttibblegapped_raw <- ttibblegapped_raw %>%
    left_join(gtfpassdf, by = c("label" = "gene_id")) %>%
    mutate(group = str_extract(label, "L1HS|L1PA2")) %>%
    mutate(group = factor(group, levels = c("L1HS", "L1PA2")))

l1hs_intact_tree <- as.treedata(ttibble_raw)
l1hs_intact_tree_gapped <- as.treedata(ttibblegapped_raw)

p <- ggtree(l1hs_intact_tree, layout = "fan") +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysave(sprintf("%s/l1hs_intact_tree_fan.png", outputdir), w = 7, h = 7)

p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysave(sprintf("%s/l1hs_intact_tree.png", outputdir), w = 5, h = 20)

p <- ggtree(l1hs_intact_tree_gapped) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")

mysave(sprintf("%s/l1hs_intact_tree_gapped.png", outputdir), w = 5, h = 20)

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
    mytheme
mysave(sprintf("%s/l1hs_intact_alnWconensus_similarity.png", outputdir), w = 10, h = 5)

p <- simplot(sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
    ggtitle("L1HS ORF1&2 Intact Consensus Accuracy") +
    geom_hline(yintercept = 0.95, col = "red", lty = 2) +
    mytheme
mysave(sprintf("%s/l1hs_intact_alnWconensus_similarity_ONTdRNAerror.png", outputdir), w = 10, h = 5)


data <- tidy_msa(l1hs_intact_aln_ss, 1, 100)

p <- ggtree(tr = l1hs_intact_tree_gapped) +
    geom_treescale(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    ggtitle("L1HS and L1PA2 ORF1&2 Intact Phylogeny") +
    geom_facet(geom = geom_msa, data = data, panel = "msa", font = NULL)
mysave(sprintf("%s/flsspasstreemsa.png", outputdir), w = 10, h = 20)



### tree with genes didn't work.
# brestidy
# mapping <- aes(xmin = start, xmax = end, fill = gene, forward = direction)
# my_pal <- colorRampPalette(rev(brewer.pal(n = 10, name = "Set3")))


# datagenes <- brestidy %>%
#     mutate(molecule = QueryID, start = Q.start_l1.3, end = Q.end_l1.3) %>%
#     dplyr::select(molecule, start, end)
# p <- ggtree(tree1, inherit.aes = FALSE) +
#     geom_treescale() +
#     geom_tippoint(aes(color = group)) +
#     geom_tiplab(aes(label = ref_status)) +
#     ggtitle("L1HS and L1PA2 ORF1&2 Intact Phylogeny") +
#     geom_facet(geom = geom_msa, data = data, panel = "msa", font = NULL, color = "Chemistry_AA") +
#     new_scale_fill() +
#     scale_fill_manual(values = my_pal(10)) +
#     geom_facet(
#         geom = geom_motif,
#         mapping = aes(xmin = start, xmax = end),
#         data = datagenes,
#         panel = "Alignment",
#         on = "LINE/L1/L1HS/4323308",
#         arrowhead_height = unit(1, "mm"),
#         arrowhead_width = unit(1, "mm")
#     )

# mysave("flsspasstreemsagene", width = 10, height = 10)



# p <- ggtree(tree2) + theme_tree() +
#     geom_tiplab()
# leaf_order <- p$data %>%
#     filter(isTip) %>%
#     arrange(y) %>%
#     pull(label)

# genesp <- ggplot(brestidy %>% mutate(a = factor(QueryID, levels = leaf_order)), aes(xmin = Q.start_l1.3, xmax = Q.end_l1.3, y = a), fill = "white") +
#     geom_gene_arrow() +
#     geom_subgene_arrow(aes(xsubmin = Q.start_orf1, xsubmax = Q.end_orf1, fill = "orf1")) +
#     geom_subgene_arrow(aes(xsubmin = Q.start_orf2, xsubmax = Q.end_orf2, fill = "orf2")) +
#     facet_wrap(~QueryID, scales = "free", ncol = 1) +
#     scale_fill_brewer(palette = "Set3") +
#     theme_genes() +
#     theme(axis.line = element_blank(), axis.text = element_blank()) +
#     theme_void()





# # how many old intact do we pick up?
# l1hsintactpath <- "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations3/RTE/l1hsintact.bed"

# l1hsintactdf <- read_delim(l1hsintactpath, col_names = FALSE)

# l1hsintact <- GRanges(
#     seqnames = l1hsintactdf$X1,
#     ranges = IRanges(start = l1hsintactdf$X2, end = l1hsintactdf$X3),
#     strand = l1hsintactdf$X6,
#     name = l1hsintactdf$X4,
#     uid = paste(l1hsintactdf$X1, l1hsintactdf$X2, l1hsintactdf$X3, l1hsintactdf$X6, sep = "_"),
#     region = l1hsintactdf$X11
# )

# # see how many l1hsintact pass our filtering
# # first just how many are in gtf
# subsetByOverlaps(l1hsintact, l1hspa2gtf)
# # ok all of them
# subsetByOverlaps(l1hsintact, gtfpass, invert = TRUE)
# # huh like 30 didn't make it
# # let's investigate why
# dmi <- subsetByOverlaps(l1hsintact, gtfpass, invert = TRUE)

# dmi
# width(dmi)
# dmiss <- getSeq(fa, dmi)
# mcols(dmiss) <- mcols(dmi)
# names(dmiss) <- mcols(dmi)$uid
# # orf analysis
# pos <- findORFs(dmiss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
# names(pos) <- names(dmiss[as.integer(names(pos))])
# tpm <- width(pos)

# dmipass <- c()
# for (element in names(pos)) {
#     if (3828 %in% width(pos[[element]]) & 1017 %in% width(pos[[element]])) {
#         dmipass <- c(dmipass, element)
#     }
# }
# dmisspass <- dmiss[dmipass]
# # so most of the l1base l1hs intact elements have orfs which are disrupted. of the three which don't, two are less than 6000bp, not sure why the third didn't pass scrutiny.


####
x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)
