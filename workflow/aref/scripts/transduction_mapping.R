source("workflow/scripts/defaults.R")
module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            filtered_tldr = "aref/A.REF_tldr/A.REF.table.kept_in_updated_ref.txt",
            r_annotation_fragmentsjoined = "aref/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/A.REF.fa",
            blast_njs = "aref/A.REF.njs"
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/A.REF_Analysis/tldr_plots/transduction_mapping.rds",
            transduction_df = "aref/A.REF_Analysis/transduction_df.csv"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "aref/A.REF"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
rmann <- left_join(rmfragments, rmfamilies)

df_filtered <- read_delim(inputs$filtered_tldr)

l1ta <- df_filtered %>%
    filter(Subfamily == "L1Ta")
l1ta_sankey <- l1ta %$% UUID

trsd <- l1ta %>%
    mutate(TransductionLen = nchar(Transduction_3p)) %>%
    filter(TransductionLen > 29)
trsd_sankey <- trsd %$% UUID
trsd %>% head(n = 4) %$% Transduction_3p


trsd_ss <- DNAStringSet(trsd$Transduction_3p)
names(trsd_ss) <- trsd$UUID


at_content <- letterFrequency(trsd_ss, letters = "AT", as.prob = TRUE)

trsd_ss_for_blast <- trsd_ss[which(at_content < 0.9)]
trsd_ss_for_blast_sankey <-names(trsd_ss_for_blast)

bl <- blast(db = sprintf("aref/%s", params$sample_or_ref))
bres <- tibble(predict(bl, trsd_ss_for_blast)) %>% left_join(rmfragments, by = c("QueryID" = "gene_id"))


bres %$% QueryID %>%
    unique() %>%
    length()

bres_hits <- bres %>%
    filter(!grepl("nonref", SubjectID)) %>%
    group_by(QueryID) %>%
    filter(Bits == max(Bits)) %>%
    filter(Perc.Ident == max(Perc.Ident)) %>%
    filter(Gap.Openings == min(Gap.Openings)) %>%
    filter(Alignment.Length == max(Alignment.Length))


l1grs <- GRanges(rmann %>% filter(grepl("L1HS|L1PA[2]", rte_subfamily)))

bres_hits_grs_prep <- bres_hits %>%
    dplyr::select(SubjectID, S.start, S.end) %>%
    group_by(QueryID) %>%
    mutate(seqnames = SubjectID) %>%
    mutate(start = min(S.start, S.end), end = max(S.start, S.end)) %>%
    dplyr::select(seqnames, start, end, QueryID)
bresgrs <- GRanges(bres_hits_grs_prep)
# extend by 500 bp on either side
bresgrs <- resize(bresgrs, width = width(bresgrs) + 10000, fix = "center")

trsd_adjacent_l1 <- subsetByOverlaps(bresgrs, l1grs)
trsd_adjacent_l1_sankey <- trsd_adjacent_l1$QueryID



`NON_REF L1TA` <- ifelse(l1ta_sankey %in% l1ta_sankey, "NON-REF L1TA", "FAIL")
`3' Tsd` <- ifelse(l1ta_sankey %in% trsd_sankey, "PASS", "FAIL")
`AT content < 0.9` <- ifelse(l1ta_sankey %in% trsd_ss_for_blast_sankey, "PASS", "FAIL")
`Source Element Identified` <- ifelse(l1ta_sankey %in% trsd_adjacent_l1_sankey, "PASS", "FAIL")

sankey_df <- tibble(
    `Non_reference_L1TA` = `NON_REF L1TA`,
    `3_Prime_Transduction` = `3' Tsd`,
    `AT_content` = `AT content < 0.9`,
    `Adjacent_Source_Element` = `Source Element Identified`
) %>%
    make_long(`Non_reference_L1TA`, `3_Prime_Transduction`, `AT_content`, `Adjacent_Source_Element`)

{
p <- ggplot(sankey_df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 12) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Non-Ref L1TA Transduction Mapping") 

flow_labels <- sankey_df %>% group_by(x, node, next_x, next_node) %>% tally() %>% drop_na()
# get corresponding positions of flows from the sankey plot
flow_info <- layer_data(p) %>% select(xmax, flow_start_ymax, flow_start_ymin) %>% distinct() # get flow related key positions related from the plot
flow_info <- flow_info[with(flow_info, order(xmax, flow_start_ymax)), ] # order the positions to match the order of flow_labels
rownames(flow_info) <- NULL # drop the row indexes
flow_info <- cbind(as.data.frame(flow_info), as.data.frame(flow_labels)) # bind the flow positions and the corresponding labels

# add labels to the flows
for (i in 1:nrow(flow_info)){
   p <- p + annotate("text", x = flow_info$xmax[i] + 0.1,
                       y = (flow_info$flow_start_ymin[i] + flow_info$flow_start_ymax[i])/2,
                       label = sprintf("%d", flow_info$n[i]), hjust = -1,
                       size = 3
                       )
 }

mysaveandstore(sprintf("%s/transduction_sankey.png", outputdir), 7, 5)
}



# l1hs_intact_aln <- read.alignment("aref/A.REF_Analysis/l1hs_intact.aln.fa", format = "fasta")
# d <- dist.alignment(l1hs_intact_aln, "identity", gap = FALSE)
# d_gapped <- dist.alignment(l1hs_intact_aln, "identity", gap = TRUE)

# # this is the pct difference in identity between sequences
# dsqaured <- 100 * d**2
# d_gappedsquared <- 100 * d_gapped**2

# str(dsqaured)
# tree_raw <- as.treedata(fastme.bal(dsqaured))
# treegapped_raw <- as.treedata(fastme.bal(d_gappedsquared))
# # to df for further modifications
# ttibble_raw <- as_tibble(tree_raw)
# ttibblegapped_raw <- as_tibble(treegapped_raw)


# rmannl1 <- rmann %>% filter(grepl("L1HS|L1PA[2]", rte_subfamily))

# non_ref_tree <- as.phylo(ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")))

# non_ref_nodes <- ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% filter(refstatus == "NonRef") %>% pull(node)
# ref_nodes <- ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% filter(refstatus == "Ref") %>% pull(node)
# ttibble_raw %>% left_join(rmannl1 %>% dplyr::select(gene_id, refstatus), by = c("label" = "gene_id")) %>% head(n = 200) %>% print(n = 200)

# castor::find_nearest_tips(as.phylo(ttibble_raw), target_tips = ref_nodes)

# tip_to_tip <- castor::find_nearest_tips(as.phylo(ttibble_raw), target_tips = ref_nodes)$nearest_tip_per_tip
# data.frame(from = seq(1:length(tip_to_tip)), to = tip_to_tip)

library(circlize)
circos.par("track.height" = 0.1)

circos.initialize(df$sectors, x = df$x)

df <- tibble(sectors = )




save(mysaveandstoreplots, file = outputs$plots)
