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
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/A.REF.fa"
        ), env = globalenv())
        assign("params", list(l13 = conf$l13fasta), env = globalenv())
        assign("outputs", list(
            outfile = "aref/A.REF_Analysis/l1element_analysis.outfile",
            plots = "aref/A.REF_Analysis/l1element_analysis.rds"
        ), env = globalenv())
    }
)


library(Rsamtools)
fa <- Rsamtools::FaFile(inputs$ref)
rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies) %>%
    filter(refstatus != "NonCentral")
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

rmann %>%
    filter(str_detect(intactness_req, "Intact")) %>%
    head() %>%
    write_delim(sprintf("%s/rmann_head.csv", outputdir), delim = "\t", col_names = TRUE)
rmann %$% intactness_req %>% table()
rmann %>%
    filter(refstatus == "NonRef") %$% rte_subfamily %>%
    table()

options(scipen = 999)

genome_length <- seqlengths(fa) %>% sum()
repeat_lengths <- rmann %>%
    group_by(repeat_superfamily) %>%
    summarise(n = n(), basepairs = sum(length)) %>%
    ungroup()
rls <- repeat_lengths %$% basepairs %>% sum()
nonrepeat_genome_length <- genome_length - rls
repeat_lengths <- add_row(repeat_lengths, repeat_superfamily = "Non Repeat", n = 1, basepairs = nonrepeat_genome_length)
repeat_lengths <- mutate(repeat_lengths, pct = basepairs / genome_length * 100) %>% filter(pct > 1)

library(ggpubr)
p <- ggpie(repeat_lengths, "pct",
    label = "repeat_superfamily",
    lab.pos = "out", fill = "repeat_superfamily", lab.adjust = 1
) + scale_palette + theme(legend.position = "none")
mysaveandstore(sprintf("%s/repeat_pie.pdf", outputdir), w = 3, h = 3)


p <- rmann %>%
    filter(refstatus == "NonRef") %>%
    mutate(rte_subfamily = gsub("SVA.*", "SVA", rte_subfamily)) %>%
    group_by(rte_subfamily, refstatus, req_integrative) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    ggplot() +
    geom_bar(aes(x = rte_subfamily, y = n, fill = req_integrative), color = "black", stat = "identity") +
    labs(title = "NonRef Elements per Subfamily", x = "Family", y = "Number of NonRef Elements", fill = "req_integrative") +
    mtclosed +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/nonref_familycount.pdf", outputdir), w = 6, h = 5)

# Summary Statistics
families_of_interest <- rmann %>%
    filter(rte_subfamily != "Other") %$% rte_family %>%
    unique()
for (fam in families_of_interest) {
    p <- rmann %>%
        filter(rte_family == fam) %>%
        filter(rte_subfamily != "Other") %>%
        group_by(rte_subfamily, refstatus, req_integrative) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        ggplot() +
        geom_bar(aes(x = rte_subfamily, y = n, fill = req_integrative), color = "black", stat = "identity") +
        labs(title = "Full Length Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
        facet_grid(~refstatus, scales = "free", space = "free_x") +
        mtclosed +
        anchorbar +
        scale_palette +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(sprintf("%s/%s_familycount.pdf", outputdir, fam), w = 6, h = 5)

    p <- rmann %>%
        filter(rte_family == fam) %>%
        filter(rte_subfamily != "Other") %>%
        filter(rte_length_req == "FL") %>%
        group_by(rte_subfamily, refstatus, req_integrative) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        ggplot() +
        geom_bar(aes(x = rte_subfamily, y = n, fill = req_integrative), color = "black", stat = "identity") +
        labs(title = "Number of Full Length L1 Elements per Subfamily", x = "Family", y = "Number of Elements", fill = "Reference Status") +
        facet_grid(~refstatus, scales = "free", space = "free_x") +
        mtclosed +
        anchorbar +
        scale_palette +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(sprintf("%s/%s_FL_familycount.pdf", outputdir, fam), w = 6, h = 5)


    intact_fams <- rmann %>%
        filter(intactness_req != "Other") %$% rte_subfamily %>%
        unique()
    p <- rmann %>%
        filter(rte_family == fam) %>%
        filter(rte_subfamily != "Other") %>%
        filter(rte_subfamily %in% intact_fams) %>%
        filter(rte_length_req == "FL") %>%
        group_by(rte_subfamily, refstatus, req_integrative) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        ggplot() +
        geom_bar(aes(x = rte_subfamily, y = n, fill = req_integrative), color = "black", stat = "identity") +
        labs(title = "FL Elements per Subfamily", x = "Family", y = "Number of FL Elements", fill = "Reference Status") +
        facet_grid(~refstatus, scales = "free_x", space = "free_x") +
        mtclosed +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_palette
    mysaveandstore(sprintf("%s/%s_FL_active_families_familycount.pdf", outputdir, fam), w = 4, h = 5)
}

# human_per_ncl_per_year_mut_rate <- 0.33e-9

# 0.02/human_per_ncl_per_year_mut_rate
# l1hs <- rmann %>%
#     filter(grepl("L1HS", rte_subfamily))
# l1hs %>% dplyr::select(gene_id, pctdiv, pctconsensuscovered, refstatus, family_av_pctdiv)
# l1hs %>% filter(pctconsensuscovered > 90) %>% dplyr::select(gene_id, pctdiv, pctconsensuscovered, refstatus, family_av_pctdiv) %>% summarise(mean(pctdiv))


# sequence analyses
library(Rsamtools)
intact_fams <- rmann %>%
    filter(rte_family == "L1") %>%
    filter(intactness_req != "Other") %$% rte_subfamily %>%
    unique()
grs <- GRanges(rmann %>% filter(rte_subfamily %in% intact_fams))
grs_df <- as.data.frame(grs) %>% tibble()


gene_ids <- grs_df$gene_id
grs_ss <- Rsamtools::getSeq(fa, grs)
mcols(grs_ss) <- mcols(grs)
names(grs_ss) <- mcols(grs)$gene_id
grs_fl_ss <- grs_ss[mcols(grs_ss)$rte_length_req == "FL"]

grs_intact_ss <- grs_fl_ss[grepl("^Intact", mcols(grs_fl_ss)$intactness_req)]

# ORF1 AND ORF2 SEQUENCE ANALYSES
# cannonical locations of RT and EN in ORF2 protein

pos <- ORFik::findORFs(grs_intact_ss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
names(pos) <- names(grs_intact_ss)
orf1 <- list()
orf2 <- list()
element_ann_df <- tibble()
for (element in names(grs_intact_ss)) {
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
    X5UTRrow <- tibble(start = 1, end = firstorfrow$start - 1) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "5UTR")
    X3UTRrow <- tibble(start = secondorfrow$end + 1, end = length(grs_intact_ss[[element]])) %>%
        mutate(gene_id = element) %>%
        mutate(feature = "3UTR")

    element_ann_df <- rbind(element_ann_df, X5UTRrow, firstorfrow, secondorfrow, X3UTRrow)
    orf1[[element]] <- grs_intact_ss[[element]][firstorf]
    orf2[[element]] <- grs_intact_ss[[element]][secondorf]
}

write_delim(element_ann_df %>% dplyr::select(gene_id, start, end, feature), sprintf("%s/intact_l1_anatomy_coordinates.tsv", outputdir), delim = "\t", col_names = TRUE)

orf1ss <- DNAStringSet(orf1)
orf2ss <- DNAStringSet(orf2)

orf1ps <- Biostrings::translate(orf1ss)
orf2ps <- Biostrings::translate(orf2ss)

writeXStringSet(orf1ps, sprintf("%s/orf1ps.fa", outputdir))
system(sprintf("echo $(pwd); mafft --auto %s/orf1ps.fa > %s/orf1ps.aln.fa", outputdir, outputdir))
aln <- read.alignment(sprintf("%s/orf1ps.aln.fa", outputdir), format = "fasta")
d <- dist.alignment(aln, "identity")
orf1Tree <- nj(d)
png(paste0(sprintf("%s/orf1Tree.png", outputdir)), width = 10, height = 30, units = "in", res = 300)
plot(orf1Tree, main = "Phylogenetic Tree of ORF1 Sequences")
dev.off()


##### ggtree
# writeXStringSet(grs_intact_ss, sprintf("%s/l1_intact.fa", outputdir))
# system(sprintf("echo $(pwd); mafft --auto %s/l1_intact.fa > %s/l1_intact.aln.fa", outputdir, outputdir))

for (subfam in mcols(grs_intact_ss)$rte_subfamily %>% unique()) {
    subfam_ss <- grs_intact_ss[mcols(grs_intact_ss)$rte_subfamily == subfam]
    writeXStringSet(subfam_ss, sprintf("%s/%s_intact.fa", outputdir, subfam))
    system(sprintf("echo $(pwd); mafft --auto %s/%s_intact.fa > %s/%s_intact.aln.fa", outputdir, subfam, outputdir, subfam))

    subfam_intact_aln <- read.alignment(sprintf("%s/%s_intact.aln.fa", outputdir, subfam), format = "fasta")
    d <- dist.alignment(subfam_intact_aln, "identity", gap = FALSE)
    d_gapped <- dist.alignment(subfam_intact_aln, "identity", gap = TRUE)

    # this is the pct difference in identity between sequences
    dsqaured <- 100 * d**2
    d_gappedsquared <- 100 * d_gapped**2

    tree_raw <- as.treedata(fastme.bal(dsqaured))
    treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
    # to df for further modifications
    ttibble_raw <- as_tibble(tree_raw)
    ttibblegapped_raw <- as_tibble(treegapped_raw)

    ttibble_raw <- ttibble_raw %>%
        left_join(as.data.frame(mcols(grs_intact_ss)), by = c("label" = "gene_id")) %>%
        mutate(group = str_extract(label, "L1HS|L1PA2")) %>%
        mutate(group = factor(group, levels = c("L1HS", "L1PA2")))

    ttibblegapped_raw <- ttibblegapped_raw %>%
        left_join(as.data.frame(mcols(grs_intact_ss)), by = c("label" = "gene_id"))
    subfam_intact_tree <- as.treedata(ttibble_raw)
    subfam_intact_tree_gapped <- as.treedata(ttibblegapped_raw)

    p <- ggtree(subfam_intact_tree, layout = "fan") +
        geom_treescale() +
        geom_tippoint(aes(color = refstatus)) +
        ggtitle(sprintf("%s ORF1&2 Intact Phylogeny", subfam)) + scale_palette
    mysaveandstore(sprintf("%s/%s_intact_tree_fan.pdf", outputdir, subfam), w = 7, h = 7)

    p <- ggtree(subfam_intact_tree) +
        geom_treescale() +
        geom_tippoint(aes(color = refstatus)) +
        geom_tiplab(aes(label = label)) +
        theme_tree2() +
        xlab("% Mutation") +
        ggtitle(sprintf("%s ORF1&2 Intact Phylogeny", subfam)) + scale_palette
    mysaveandstore(sprintf("%s/%s_intact_tree.pdf", outputdir, subfam), w = 7.5, h = 20)
    # make the tree horizontal
    p <- ggtree(subfam_intact_tree) +
        geom_treescale() +
        geom_tippoint(aes(color = refstatus)) +
        geom_tiplab(aes(label = label, angle = 90)) +
        theme_tree2() +
        xlab("% Mutation") +
        ggtitle(sprintf("%s ORF1&2 Intact Phylogeny", subfam)) +
        coord_flip() + scale_palette
    mysaveandstore(sprintf("%s/%s_intact_tree_horizontal.pdf", outputdir, subfam), w = 20, h = 8)


    p <- ggtree(subfam_intact_tree_gapped) +
        geom_treescale() +
        geom_tippoint(aes(color = refstatus)) +
        geom_tiplab(aes(label = label)) +
        theme_tree2() +
        xlab("% Mutation") +
        ggtitle(sprintf("%s ORF1&2 Intact Phylogeny", subfam)) + scale_palette

    mysaveandstore(sprintf("%s/%s_intact_tree_gapped.pdf", outputdir, subfam), w = 7.5, h = 20)


    subfam_intact_aln <- readDNAMultipleAlignment(sprintf("%s/%s_intact.aln.fa", outputdir, subfam), format = "fasta")
    consensus <- consensusString(subfam_intact_aln)
    subfam_intact_aln_ss <- as(subfam_intact_aln, "DNAStringSet")
    temp <- c(subfam_intact_aln_ss %>% as.character(), consensus %>% as.character() %>% setNames("consensus"))
    subfam_intact_aln_c <- DNAMultipleAlignment(temp)
    subfam_intact_aln_ss_c <- as(subfam_intact_aln_c, "DNAStringSet")
    writeXStringSet(subfam_intact_aln_ss_c, file = sprintf("%s/%s_intact_alnWconensus.fa", outputdir, subfam))

    p <- simplot(sprintf("%s/%s_intact_alnWconensus.fa", outputdir, subfam), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
        ggtitle("%s Intact Consensus Accuracy") +
        mtopen
    mysaveandstore(sprintf("%s/%s_intact_alnWconensus_similarity.pdf", outputdir, subfam), w = 10, h = 5)

    p <- simplot(sprintf("%s/%s_intact_alnWconensus.fa", outputdir, subfam), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
        ggtitle(sprintf("%s Intact Consensus Accuracy", subfam)) +
        geom_hline(yintercept = 0.95, col = "red", lty = 2) +
        mtopen
    mysaveandstore(sprintf("%s/%s_intact_alnWconensus_similarity_ONTdRNAerror.pdf", outputdir, subfam), w = 10, h = 5)
}

####
# save(mysaveandstoreplots, file = outputs$plots)
if (conf$store_env_as_rds == "yes") {
    save.image(file = outputs$plots)
} else {
    x <- tibble(Env_file = "Opted not to store environment. If this is not desired, change 'store_plots_as_rds' to 'yes' in the relevant config file and rerun this rule.")
    write_tsv(x, file = outputs$plots)
}
