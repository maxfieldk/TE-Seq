module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

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
library(ggsankey)
library(circlize)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = sprintf("aref/%s_tldr/tldr.table.txt", sample = sample_table$sample_name),
            r_annotation_fragmentsjoined = sprintf("aref/%s_annotations/%s_repeatmasker.gtf.rformatted.fragmentsjoined.csv", sample_table$sample_name, sample_table$sample_name),
            r_repeatmasker_annotation = sprintf("aref/%s_annotations/%s_repeatmasker_annotation.csv", sample_table$sample_name, sample_table$sample_name),
            filtered_tldr = sprintf("aref/%s_tldr/%s.table.kept_in_updated_ref.txt", sample_table$sample_name, sample_table$sample_name),
            sample_refs = sprintf("aref/%s.fa", sample_table$sample_name),
            blast_njs = sprintf("aref/%s.njs", sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(
            plots = "ldna/results/insertion_analysis/insertion_analysis.rds",
            transduction_df = "ldna/results/insertion_analysis/transduction_df.csv",
            circlize_transduction = "ldna/results/insertion_analysis/circlize_transduction.png",
            circlize_phylogeny = "ldna/results/insertion_analysis/circlize_phylogeny.png"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

dfs_filtered <- list()
for (sample in sample_table$sample_name) {
    df <- read.table(grep(sprintf("%s_tldr", sample), inputs$filtered_tldr, value = TRUE), header = TRUE)
    df$sample_name <- sample
    df <- df %>% left_join(sample_table)
    dfs_filtered[[sample]] <- df
}

dff <- do.call(rbind, dfs_filtered) %>% tibble()

anns <- list()
for (sample in sample_table$sample_name) {
    df1 <- read_csv(grep(sprintf("%s_annotations", sample), inputs$r_annotation_fragmentsjoined, value = TRUE))
    df1$sample_name <- sample
    df2 <- read_csv(grep(sprintf("%s_annotations", sample), inputs$r_repeatmasker_annotation, value = TRUE))
    df <- df1 %>%
        left_join(sample_table) %>%
        left_join(df2)
    anns[[sample]] <- df
}
rmann <- do.call(rbind, anns) %>% tibble()
ann <- rmann %>% filter(refstatus == "NonRef")



p <- ann %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_family, scales = "free", space = "free_x", ncol = 3) +
    geom_bar(aes(fill = req_integrative)) +
    labs(x = "", y = "Count") +
    scale_palette +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/insertions_by_family_req.pdf", outputdir), 10, 4)

p <- ann %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_family, scales = "free", space = "free_x", ncol = 3) +
    geom_bar(aes(fill = loc_integrative)) +
    labs(x = "", y = "Count") +
    scale_palette +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/insertions_by_family_loc.pdf", outputdir), 10, 4)

p <- ann %>%
    filter(rte_family == "L1") %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_subfamily, scales = "free", space = "free_x", ncol = 3) +
    geom_bar(aes(fill = loc_integrative)) +
    labs(x = "", y = "Count") +
    scale_palette +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/l1_insertions_by_subfamily_loc.pdf", outputdir), 10, 4)

p <- ann %>%
    filter(rte_family == "L1") %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_subfamily, scales = "free", space = "free_x", ncol = 3) +
    geom_bar(aes(fill = req_integrative)) +
    labs(x = "", y = "Count") +
    scale_palette +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/l1_insertions_by_subfamily_req.pdf", outputdir), 10, 4)


dff %>% tibble()
p <- dff %>%
    filter(Family == "L1") %>%
    group_by(Chrom, Start, End, Subfamily) %>%
    summarize(n = n()) %>%
    arrange(n) %$% n %>%
    table() %>%
    as.data.frame() %>%
    mutate(Private = . == 1) %>%
    ggplot(aes(x = ., y = Freq)) +
    geom_col(aes(fill = Private)) +
    labs(x = "Insertion found in # samples", y = "Count") +
    ggtitle("L1 Insertions") +
    anchorbar +
    mtopengridh +
    scale_palette
mysaveandstore(sprintf("%s/l1_insertions_per_sample.pdf", outputdir), 4, 4)


insertions_filtered <- dff %>%
    filter(Family == "L1") %>%
    group_by(Chrom, Start, End, Subfamily) %>%
    arrange(desc(length(Transduction_3p))) %>%
    filter(row_number() == 1) %>%
    filter(Filter == "PASS")
# l1ta_sankey <- insertions_filtered %$% UUID
trsd_list <- list()
for (sample in sample_table$sample_name) {
    trsd <- insertions_filtered %>%
        filter(sample_name == sample) %>%
        mutate(TransductionLen = nchar(Transduction_3p)) %>%
        filter(TransductionLen > 20)
    trsd_sankey <- trsd %$% UUID
    trsd %>% head(n = 4) %$% Transduction_3p
    trsd_ss <- DNAStringSet(trsd$Transduction_3p)
    names(trsd_ss) <- trsd$UUID
    at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)
    trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]
    trsd_ss_for_blast_sankey <- names(trsd_ss_for_blast)

    bl <- blast(db = sub("\\.[^.]*$", "", grep(sprintf("%s.njs", sample), inputs$blast_njs, value = TRUE)))
    tryCatch(
        {
            bres <- tibble(predict(bl, trsd_ss_for_blast)) %>% left_join(insertions_filtered, by = c("QueryID" = "UUID"))
            bres %$% QueryID %>%
                unique() %>%
                length()

            bres_hits <- bres %>%
                filter(SubjectID != faName) %>%
                group_by(QueryID) %>%
                arrange(desc(Bits), desc(Perc.Ident), desc(Gap.Openings), desc(Alignment.Length)) %>%
                ungroup()
            # filter(Bits == max(Bits)) %>%
            # filter(Perc.Ident == max(Perc.Ident)) %>%
            # filter(Gap.Openings == min(Gap.Openings)) %>%
            # filter(Alignment.Length == max(Alignment.Length))


            l1grs <- GRanges(rmann %>% filter(sample_name == sample) %>% filter(grepl("L1HS|L1PA[2]", rte_subfamily)))

            bres_hits_grs_prep <- bres_hits %>%
                mutate(seqnames = SubjectID) %>%
                rowwise() %>%
                mutate(start = min(S.start, S.end), end = max(S.start, S.end)) %>%
                relocate(seqnames, start, end) %>%
                dplyr::select(seqnames, start, end, QueryID, Bits, Perc.Ident, Gap.Openings, Alignment.Length)
            bresgrs <- GRanges(bres_hits_grs_prep)
            # extend by 10000 bp on either side
            bresgrs <- resize(bresgrs, width = width(bresgrs) + 1000, fix = "center")

            trsd_adjacent_l1 <- subsetByOverlaps(bresgrs, l1grs)
            trsd_adjacent_l1 <- mergeByOverlaps(bresgrs, l1grs)

            trsd_adjacent_l1_df <- as.data.frame(trsd_adjacent_l1) %>%
                tibble() %>%
                dplyr::select(starts_with("bresgrs."), gene_id) %>%
                dplyr::rename(from_gene_id = gene_id)
            colnames(trsd_adjacent_l1_df) <- gsub("bresgrs.", "", colnames(trsd_adjacent_l1_df))
            trsd_adjacent_l1_df_final <- trsd_adjacent_l1_df %>%
                group_by(QueryID) %>%
                mutate(number_productive_hits = n()) %>%
                filter(Bits == max(Bits)) %>%
                ungroup()
            trsd_adjacent_l1_sankey <- trsd_adjacent_l1$QueryID
            trsd_list[[sample]] <- trsd_adjacent_l1_df_final
        },
        error = function(e) {
            print("Likely no blast hits")
        }
    )
}

TSDcalls <- do.call(rbind, trsd_list) %>% tibble()
from_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
to_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
for (index in 1:nrow(TSDcalls)) {
    from_gene_id <- TSDcalls[index, ]$from_gene_id
    if (grepl("nonref", TSDcalls[index, ]$seqnames)) {
        from_chr <- TSDcalls[index, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(3)
        from_start <- TSDcalls[index, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(4) %>%
            as.numeric()
        from_end <- TSDcalls[index, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(5) %>%
            as.numeric()
    } else {
        from_chr <- TSDcalls[index, ]$seqnames
        from_start <- TSDcalls[index, ]$start %>% as.numeric()
        from_end <- TSDcalls[index, ]$end %>% as.numeric()
    }

    from_df_row <- tibble(gene_id = from_gene_id, chr = from_chr, start = from_start, end = from_end)
    from_df <- bind_rows(from_df, from_df_row)


    to <- rmann %>%
        filter(seqnames == insertions_filtered %>%
            ungroup() %>%
            filter(UUID == TSDcalls[index, ]$QueryID) %$% faName) %>%
        filter(refstatus == "NonRef")
    to_gene_id <- to[1, ] %$% gene_id

    if (grepl("nonref", to[1, ]$seqnames)) {
        to_chr <- to[1, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(3)
        to_start <- to[1, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(4) %>%
            as.numeric()
        to_end <- to[1, ]$seqnames %>%
            str_split("_") %>%
            unlist() %>%
            pluck(5) %>%
            as.numeric()
    }

    to_df_row <- tibble(chr = to_chr, start = to_start, end = to_end, gene_id = to_gene_id, )
    to_df <- bind_rows(to_df, to_df_row)
}


from_elements <- from_df %>%
    pull(gene_id) %>%
    unique()
from_elements_colors <- tibble(gene_id = from_elements, color = adjust_transparency(as.character(paletteer_c("grDevices::Roma", direction = 1, n = length(from_elements))), alpha = 1))
from_df <- left_join(from_df, from_elements_colors, by = "gene_id")
cytobands <- read.cytoband(
    cytoband = ,
    species = NULL,
    sort.chr = TRUE
)

{
    png(outputs$circlize_transduction, width = 6, height = 6, units = "in", res = 300)

    circos.initializeWithIdeogram(
        cytoband = "/users/mkelsey/data/LF1/RTE/aref/annotations/cytobands.bed",
        sort.chr = TRUE,
        major.by = NULL,
        plotType = c("ideogram", "axis", "labels"),
        track.height = NULL,
        ideogram.height = convert_height(5, "mm")
    )

    circos.genomicLink(from_df, to_df, col = from_df$color, lwd = 2, direction = 1)
    title("L1HS Transduction Mapping")
    dev.off()
}



number_of_offspring <- from_df %>%
    group_by(gene_id) %>%
    summarize(n = n())
p <- number_of_offspring %>% ggplot(aes(x = fct_reorder(gene_id, n), y = n)) +
    geom_col() +
    mtopen +
    coord_flip() +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    labs(x = "Source Element", y = "Number of Offspring", caption = "Transduction-based Source Element Mapping") +
    anchorbar
mysaveandstore(sprintf("%s/offspring_per_source_element_transduction.pdf", outputdir), 4, 6)


# `NON_REF L1TA` <- ifelse(l1ta_sankey %in% l1ta_sankey, "NON-REF L1TA", "FAIL")
# `3' Tsd` <- ifelse(l1ta_sankey %in% trsd_sankey, "PASS", "FAIL")
# `AT content < 0.9` <- ifelse(l1ta_sankey %in% trsd_ss_for_blast_sankey, "PASS", "FAIL")
# `Source Element Identified` <- ifelse(l1ta_sankey %in% trsd_adjacent_l1_sankey, "PASS", "FAIL")

# sankey_df <- tibble(
#     `Non_reference_L1TA` = `NON_REF L1TA`,
#     `3_Prime_Transduction` = `3' Tsd`,
#     `AT_content` = `AT content < 0.9`,
#     `Adjacent_Source_Element` = `Source Element Identified`
# ) %>%
#     make_long(`Non_reference_L1TA`, `3_Prime_Transduction`, `AT_content`, `Adjacent_Source_Element`)

# {
# p <- ggplot(sankey_df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
# geom_sankey(flow.alpha = .6,
#             node.color = "gray30") +
# geom_sankey_label(size = 3, color = "white", fill = "gray40") +
# theme_sankey(base_size = 12) +
# labs(x = NULL, y = NULL) +
# theme(legend.position = "none",
#         plot.title = element_text(hjust = .5)) +
# ggtitle("Non-Ref L1TA Transduction Mapping")

# flow_labels <- sankey_df %>% group_by(x, node, next_x, next_node) %>% tally() %>% drop_na()
# # get corresponding positions of flows from the sankey plot
# flow_info <- layer_data(p) %>% select(xmax, flow_start_ymax, flow_start_ymin) %>% distinct() # get flow related key positions related from the plot
# flow_info <- flow_info[with(flow_info, order(xmax, flow_start_ymax)), ] # order the positions to match the order of flow_labels
# rownames(flow_info) <- NULL # drop the row indexes
# flow_info <- cbind(as.data.frame(flow_info), as.data.frame(flow_labels)) # bind the flow positions and the corresponding labels

# # add labels to the flows
# for (i in 1:nrow(flow_info)){
# p <- p + annotate("text", x = flow_info$xmax[i] + 0.1,
#                     y = (flow_info$flow_start_ymin[i] + flow_info$flow_start_ymax[i])/2,
#                     label = sprintf("%d", flow_info$n[i]), hjust = -1,
#                     size = 3
#                     )
# }

# mysaveandstore(sprintf("%s/transduction_sankey.pdf", outputdir), 7, 5)
# }

rmann %$% intactness_req %>% unique()
all_intact_l1hs <- rmann %>% filter(intactness_req == "Intact")
all_intact_l1hs <- all_intact_l1hs %>%
    group_by(seqnames, start, end) %>%
    filter(row_number() == 1) %>%
    ungroup()
all_intact_l1hs <- all_intact_l1hs %>% mutate(unique_gene_id = paste0(sample_name, ":", gene_id))
all_seqs <- list()
for (sample in sample_table$sample_name) {
    intact_l1hs_df <- all_intact_l1hs %>% filter(sample_name == sample)
    intact_l1hs_grs <- GRanges(intact_l1hs_df)
    # get biostrings from fasta
    fa <- Rsamtools::FaFile(grep(sprintf("%s.fa", sample), inputs$sample_refs, value = TRUE))

    # get the sequences
    seqs <- getSeq(fa, intact_l1hs_grs)
    names(seqs) <- intact_l1hs_df$unique_gene_id
    all_seqs[[sample]] <- seqs
}
all_seqs_ss <- Reduce(c, all_seqs)
writeXStringSet(all_seqs_ss, sprintf("%s/intact_l1hs.fa", outputdir))


system(sprintf("echo $(pwd); mafft --auto %s/intact_l1hs.fa > %s/intact_l1hs.aln.fa", outputdir, outputdir))

l1hs_intact_aln <- read.alignment(sprintf("%s/intact_l1hs.aln.fa", outputdir), format = "fasta")
d <- dist.alignment(l1hs_intact_aln, "identity", gap = FALSE)
d_gapped <- dist.alignment(l1hs_intact_aln, "identity", gap = TRUE)

# this is the pct difference in identity between sequences
dsqaured <- 100 * d**2
d_gappedsquared <- 100 * d_gapped**2

str(dsqaured)
tree_raw <- as.treedata(fastme.bal(dsqaured))
treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
# to df for further modifications
ttibble_raw <- as_tibble(tree_raw)
ttibblegapped_raw <- as_tibble(treegapped_raw)



non_ref_tree <- as.phylo(ttibble_raw %>% left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")))

non_ref_nodes <- ttibble_raw %>%
    left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")) %>%
    filter(refstatus == "NonRef") %>%
    pull(node)
ref_nodes <- ttibble_raw %>%
    left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")) %>%
    filter(refstatus == "Ref") %>%
    pull(node)
treedf <- ttibble_raw %>%
    left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")) %>%
    filter(!is.na(refstatus))
nearest_tip_per_tip <- castor::find_nearest_tips(as.phylo(ttibble_raw), target_tips = ref_nodes)$nearest_tip_per_tip
treedf$nearest_ref_tip <- nearest_tip_per_tip

treedf_nonref <- treedf %>%
    filter(refstatus == "NonRef") %>%
    select(node, nearest_ref_tip) %>%
    dplyr::rename(from = nearest_ref_tip, to = node)
# initialize new data frame
from_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
to_df <- tibble(chr = character(), start = integer(), end = integer(), gene_id = character())
for (index in 1:nrow(treedf_nonref)) {
    from <- treedf_nonref[index, ]$from
    to <- treedf_nonref[index, ]$to
    from_gene_id <- treedf %>%
        filter(node == from) %>%
        pull(label)
    to_gene_id <- treedf %>%
        filter(node == to) %>%
        pull(label)
    from_chr <- all_intact_l1hs %>%
        filter(unique_gene_id == from_gene_id) %>%
        pull(seqnames)
    to_chr <- all_intact_l1hs %>%
        filter(unique_gene_id == to_gene_id) %>%
        pull(seqnames) %>%
        str_extract("chr[0-9XY]+")
    from_start <- all_intact_l1hs %>%
        filter(unique_gene_id == from_gene_id) %>%
        pull(start)
    to_start <- all_intact_l1hs %>%
        filter(unique_gene_id == to_gene_id) %>%
        pull(start)
    from_end <- all_intact_l1hs %>%
        filter(unique_gene_id == from_gene_id) %>%
        pull(end)
    to_end <- all_intact_l1hs %>%
        filter(unique_gene_id == to_gene_id) %>%
        pull(end)

    from_df_row <- tibble(gene_id = from_gene_id, chr = from_chr, start = from_start, end = from_end)
    to_df_row <- tibble(gene_id = to_gene_id, chr = to_chr, start = to_start, end = to_end)
    from_df <- bind_rows(from_df, from_df_row)
    to_df <- bind_rows(to_df, to_df_row)
}

from_elements <- from_df %>%
    pull(gene_id) %>%
    unique()
from_elements_colors <- tibble(gene_id = from_elements, color = adjust_transparency(as.character(paletteer_c("grDevices::Roma", direction = 1, n = length(from_elements))), alpha = 1))
from_df <- left_join(from_df, from_elements_colors, by = "gene_id")
cytobands <- read.cytoband(
    cytoband = ,
    species = NULL,
    sort.chr = TRUE
)

{
    png(outputs$circlize_phylogeny, width = 6, height = 6, units = "in", res = 300)

    circos.initializeWithIdeogram(
        cytoband = "/users/mkelsey/data/LF1/RTE/aref/annotations/cytobands.bed",
        sort.chr = TRUE,
        major.by = NULL,
        plotType = c("ideogram", "axis", "labels"),
        track.height = NULL,
        ideogram.height = convert_height(5, "mm")
    )

    circos.genomicLink(from_df, to_df, col = from_df$color, lwd = 2, direction = 1)
    title("L1HS Source Element Mapping")
    dev.off()
}



number_of_offspring <- from_df %>%
    group_by(gene_id) %>%
    summarize(n = n())
p <- number_of_offspring %>% ggplot(aes(x = fct_reorder(gene_id, n), y = n)) +
    geom_col() +
    mtopen +
    coord_flip() +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    labs(x = "Source Element", y = "Number of Offspring", caption = "Phylogeny-based Source Element Mapping") +
    anchorbar
mysaveandstore(sprintf("%s/offspring_per_source_element.pdf", outputdir), 4, 6)

####### TREES
tree_raw <- as.treedata(fastme.bal(dsqaured))
treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
# to df for further modifications
ttibble_raw <- as_tibble(tree_raw)
ttibblegapped_raw <- as_tibble(treegapped_raw)

ttibble_raw <- ttibble_raw %>%
    left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")) %>%
    mutate(label = gsub(".*:", "", label))


ttibblegapped_raw <- ttibblegapped_raw %>%
    left_join(all_intact_l1hs %>% dplyr::select(unique_gene_id, refstatus), by = c("label" = "unique_gene_id")) %>%
    mutate(label = gsub(".*:", "", label))


l1hs_intact_tree <- as.treedata(ttibble_raw)
l1hs_intact_tree_gapped <- as.treedata(ttibblegapped_raw)

p <- ggtree(l1hs_intact_tree, layout = "fan") +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    scale_palette +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysaveandstore(sprintf("%s/l1hs_intact_tree_fan.pdf", outputdir), w = 7, h = 7)

p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    scale_palette +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysaveandstore(sprintf("%s/l1hs_intact_tree.pdf", outputdir), w = 5, h = 20)

p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    scale_palette +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")
mysaveandstore(sprintf("%s/l1hs_intact_tree_nolabels.pdf", outputdir), w = 5, h = 20)

# make the tree horizontal
p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label, angle = 90)) +
    scale_palette +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny") +
    coord_flip() + theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
mysaveandstore(sprintf("%s/l1hs_intact_tree_horizontal.pdf", outputdir), w = 20, h = 8)

p <- ggtree(l1hs_intact_tree) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    scale_palette +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny") +
    layout_dendrogram() + theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
mysaveandstore(sprintf("%s/l1hs_intact_tree_horizontal_nolabels.pdf", outputdir), w = 20, h = 8)

p <- ggtree(l1hs_intact_tree_gapped) +
    geom_treescale() +
    geom_tippoint(aes(color = refstatus)) +
    geom_tiplab(aes(label = label)) +
    scale_palette +
    theme_tree2() +
    xlab("% Mutation") +
    ggtitle("L1HS ORF1&2 Intact Phylogeny")

mysaveandstore(sprintf("%s/l1hs_intact_tree_gapped.pdf", outputdir), w = 5, h = 20)

## tree with msa
# l1hs_intact_aln <- readDNAMultipleAlignment(sprintf("%s/l1hs_intact.aln.fa", outputdir), format = "fasta")
# consensus <- consensusString(l1hs_intact_aln)
# l1hs_intact_aln_ss <- as(l1hs_intact_aln, "DNAStringSet")
# temp <- c(l1hs_intact_aln_ss %>% as.character(), consensus %>% as.character() %>% setNames("consensus"))
# l1hs_intact_aln_c <- DNAMultipleAlignment(temp)
# l1hs_intact_aln_ss_c <- as(l1hs_intact_aln_c, "DNAStringSet")
# writeXStringSet(l1hs_intact_aln_ss_c, file = sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir))
# p <- simplot(sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
#     ggtitle("L1HS ORF1&2 Intact Consensus Accuracy") +
#     mtopen
# mysaveandstore(sprintf("%s/l1hs_intact_alnWconensus_similarity.pdf", outputdir), w = 10, h = 5)

# p <- simplot(sprintf("%s/l1hs_intact_alnWconensus.fa", outputdir), "consensus", window = 1, group = TRUE, id = 1, sep = "_") +
#     ggtitle("L1HS ORF1&2 Intact Consensus Accuracy") +
#     geom_hline(yintercept = 0.95, col = "red", lty = 2) +
#     mtopen
# mysaveandstore(sprintf("%s/l1hs_intact_alnWconensus_similarity_ONTdRNAerror.pdf", outputdir), w = 10, h = 5)


# data <- tidy_msa(l1hs_intact_aln_ss, 1, 100)

# p <- ggtree(tr = l1hs_intact_tree_gapped) +
#     geom_treescale(aes(color = refstatus)) +
#     geom_tiplab(aes(label = label)) +
#     ggtitle("L1HS and L1PA2 ORF1&2 Intact Phylogeny") +
#     geom_facet(geom = geom_msa, data = data, panel = "msa", font = NULL)
# mysaveandstore(sprintf("%s/l1hs_intact_msa_first_100nt_treemsa.pdf", outputdir), w = 10, h = 20)

####

x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)
