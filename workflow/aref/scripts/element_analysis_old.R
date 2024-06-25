

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


