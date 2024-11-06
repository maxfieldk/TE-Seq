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
library(ORFik)


tryCatch(
    {
        inputs <- snakemake@input
        params <- snakemake@params
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "aref/extended/A.REF_annotations/A.REF_repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            r_repeatmasker_annotation = "aref/extended/A.REF_annotations/A.REF_repeatmasker_annotation.csv",
            ref = "aref/extended/A.REF.fa"
        ), env = globalenv())
        assign("params", list(l13 = conf$l13fasta), env = globalenv())
        assign("outputs", list(), env = globalenv())
    }
)


library(Rsamtools)
fa <- Rsamtools::FaFile(inputs$ref)
rmann <- read_csv(conf$rmann) %>%
    filter(refstatus != "NonCentral")

rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)

##############################
outputdir_already_computed <- "/users/mkelsey/data/LF1/RTE/ldna/results/plots/l1_alignment_meth"
outputdir <- "ldna/results/plots/l1_alignment_meth"

subfam <- "L1HS"
grs_fl <- rmann %>%
    filter(rte_length_req == "FL") %>%
    filter(rte_subfamily == subfam) %>%
    group_by(rte_subfamily) %>%
    mutate(length99 = quantile(length, 0.99)) %>%
    ungroup() %>%
    filter(length < length99 * 1.05) %>%
    GRanges()
grs_fl_ss <- getSeq(fa, grs_fl)
names(grs_fl_ss) <- mcols(grs_fl)$gene_id
fl_fa_path <- sprintf("%s/alignments/%s_fl.fa", outputdir, subfam)
dir.create(dirname(fl_fa_path), recursive = TRUE)
writeXStringSet(grs_fl_ss, fl_fa_path)
system(sprintf("echo $(pwd); mafft --auto %s/alignments/%s_fl.fa > %s/alignments/%s_fl.aln.fa", outputdir, subfam, outputdir, subfam))


aln <- readDNAMultipleAlignment(sprintf("%s/alignments/%s_fl.aln.fa", outputdir, subfam))
alnss <- readDNAStringSet(sprintf("%s/alignments/%s_fl.aln.fa", outputdir, subfam))
names(alnss) %>%
    duplicated() %>%
    sum()

consensusMat <- consensusMatrix(aln)
base_names <- rownames(consensusMat)
max_indices <- apply(consensusMat, 2, which.max)
max_indices %>% table()
consensus_vec <- base_names[max_indices]
consensus <- paste(consensus_vec, collapse = "")

consensus_ss_with_gaps <- as.character(consensus) %>% DNAStringSet()
names(consensus_ss_with_gaps) <- "consensus"
consensus_with_gaps_path <- sprintf("%s/alignments/%s_fl_consensus_with_gaps.fa", outputdir, subfam)
writeXStringSet(consensus_ss_with_gaps, consensus_with_gaps_path)
consensus_ss <- as.character(consensus) %>%
    gsub("-", "", .) %>%
    DNAStringSet()
names(consensus_ss) <- sprintf("%s_fl_consensus", subfam)
consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir, subfam)
writeXStringSet(consensus_ss, consensus_path)
consensus_ss <- readDNAStringSet(consensus_path)

getIndexMap <- function(gapped_seq_ss) {
    name <- names(gapped_seq_ss)
    gapped_seq_split <- as.character(gapped_seq_ss) %>% str_split_1(pattern = "")
    alignment_index_vec <- c()
    seq_index_vec <- c()
    seq_begun <- FALSE
    for (alignment_index in seq(1:length(gapped_seq_split))) {
        alignment_index_vec <- c(alignment_index_vec, alignment_index)
        if (seq_begun) {
            if (gapped_seq_split[[alignment_index]] != "-") {
                seq_index <- seq_index + 1
                seq_index_vec <- c(seq_index_vec, seq_index)
            } else {
                seq_index_vec <- c(seq_index_vec, NA)
            }
        } else if (gapped_seq_split[[alignment_index]] != "-") {
            seq_begun <- TRUE
            seq_index <- 1
            seq_index_vec <- c(seq_index_vec, seq_index)
        } else {
            seq_index_vec <- c(seq_index_vec, NA)
        }
    }
    tt <- tibble(alignment_pos = alignment_index_vec, gene_id = name, sequence_pos = seq_index_vec, seq = gapped_seq_split)
    return(tt)
}

for (seq_name in names(alnss)) {
    print(seq_name)
    ss <- alnss[seq_name]
    element_map <- getIndexMap(ss)
    if (!exists("index_df")) {
        index_df <<- element_map
    } else {
        index_df <<- bind_rows(index_df, element_map)
    }
}

alignment_index_long <<- bind_rows(index_df, getIndexMap(consensus_ss_with_gaps))

consensus_frame <- getIndexMap(consensus_ss_with_gaps) %>%
    dplyr::rename(consensus_pos = sequence_pos, consensus_seq = seq) %>%
    select(-gene_id)
consensus_index_long <- index_df %>%
    left_join(consensus_frame) %>%
    filter(!(is.na(sequence_pos) & is.na(consensus_pos)))

write_csv(alignment_index_long, sprintf("%s/%s_fl_mapping_to_alignment_table.csv", outputdir, subfam))
alignment_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_alignment_table.csv", outputdir, subfam))
write_csv(consensus_index_long, sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir, subfam))
consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir, subfam))


# Combine the individual data frames into one tidy data frame
cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()
consensus_ss[[1]][5762:5763]
cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

methdf <- rtedf %>% filter(rte_subfamily == subfam)
mdf <- methdf %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

senseelement <- mdf %>%
    filter(rte_strand == "+") %$% gene_id %>%
    pluck(1)

antisenseelement <- mdf %>%
    filter(rte_strand == "-") %$% gene_id %>%
    pluck(1)

cpgmapping_check <- cg_positions_df %>%
    filter(gene_id == senseelement) %$% sequence_pos %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == senseelement) %>%
    relocate(sequence_pos) %$% sequence_pos %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)
cpgmapping_check <- cg_positions_df %>%
    filter(gene_id == antisenseelement) %$% sequence_pos %>%
    unique() %>%
    sort()
methdf_check <- mdf %>%
    filter(gene_id == antisenseelement) %>%
    relocate(sequence_pos) %$% sequence_pos %>%
    unique() %>%
    sort()
print(methdf_check)
print(cpgmapping_check)

merged <- left_join(cg_positions_df, mdf, by = c("gene_id", "sequence_pos"))
cpg_order <- merged %$% consensus_pos %>%
    unique() %>%
    sort()
merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))
library(tidyHeatmap)

dat <- merged %>%
    filter(!is.na(pctM))

write_csv(dat, "meth_for_bayes.csv")


# methylation heatmaps
hms <- list()
for (sample in conf$samples) {
    p <- merged %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
    hms[[sample]] <- p
    dir.create(outputdir, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir, sample, subfam), w = 6, h = 6)
}



# motif analysis of consensus
motif_pwms_path <- "/oscar/home/mkelsey/anaconda/gimme/lib/python3.10/site-packages/data/motif_databases/HOCOMOCOv11_HUMAN.pfm"
conda_base_path <- system("conda info --base", intern = TRUE)
motif_search_output_path <- sprintf("%s/alignments/%s_fl_consensus.motifs.bed", outputdir, subfam)

system(
    sprintf(
        "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme scan %s -g %s -p %s -c 0.9 -n 10 -b > %s",
        conda_base_path, consensus_path, inputs$ref, motif_pwms_path, motif_search_output_path
    )
)
motifs_df <- as.data.frame(rtracklayer::import(motif_search_output_path)) %>% tibble()


motif_search_output_path <- sprintf("%s/alignments/%s_fl.motifs.bed", outputdir, subfam)
system(
    sprintf(
        "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme scan %s -g %s -p %s -c 0.9 -n 10 -b > %s",
        conda_base_path, fl_fa_path, inputs$ref, motif_pwms_path, motif_search_output_path
    )
)

motifs_all_df <- as.data.frame(rtracklayer::import(motif_search_output_path)) %>%
    tibble() %>%
    dplyr::rename(gene_id = seqnames)
motifs_all_df

# analysis of methylation by yy1 status
yy1_pos_consensus <- motifs_df %>% filter(grepl("YY1", name))
yy1_pos <- motifs_all_df %>% filter(grepl("YY1", name))

yy1sitesdf <- tibble(gene_id = yy1_pos$gene_id %>% unique(), yy1status = TRUE)
yy1sitesdf <- tibble(gene_id = merged$gene_id %>% unique()) %>%
    left_join(yy1sitesdf) %>%
    mutate(yy1status = replace_na(yy1status, FALSE))


merged1 <- merged %>% left_join(yy1sitesdf)
merged1 %$% yy1status %>% table()

hms <- list()
for (sample in conf$samples) {
    p <- merged1 %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        group_by(yy1status) %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
    hms[[sample]] <- p
    dir.create(outputdir, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s_yy1stratified.pdf", outputdir, sample, subfam), w = 6, h = 6)
}



merged1 %>%
    filter(!is.na(sample)) %>%
    mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
    mutate(UTR = ifelse((consensus_pos) < 909, TRUE, FALSE)) %>%
    mutate(UTRcpgisland = ifelse((consensus_pos) < 500, TRUE, FALSE)) %>%
    filter(UTRcpgisland == TRUE) %>%
    group_by(gene_id, sample, condition, yy1status) %>%
    summarise(pctM = mean(pctM)) %>%
    group_by(condition, yy1status) %>%
    summarise(pctM = mean(pctM))




# CLUSTER ANALYSIS
matdf <- merged1 %>%
    dplyr::select(consensus_pos, gene_id, pctM) %>%
    pivot_wider(names_from = consensus_pos, values_from = pctM)
rownames(matdf) <- matdf$gene_id
mat <- as.matrix(matdf %>% dplyr::select(-gene_id))
rownames(mat) <- matdf$gene_id

library(impute)
# Impute missing values in pctM (default is K-nearest neighbors)
imputed_data <- impute.knn(mat)$data
# Create a distance matrix from the imputed data
distance_matrix <- dist(imputed_data, method = "euclidean")
# Perform clustering
clusters <- hclust(distance_matrix, method = "ward.D2")
cluster_group <- cutree(clusters, k = 3)
clusterdf <- tibble(gene_id = names(cluster_group), cluster = cluster_group)
clusterdf %>% print(n = 1000)

merged1 %$% gene_id %>% unique()
merged1 %>%
    left_join(clusterdf) %>%
    filter(sample == !!sample) %>%
    group_by(gene_id) %>%
    mutate(cpgs_detected_per_element = n()) %>%
    ungroup() %$% gene_id %>%
    unique()
merged1 %>%
    filter(sample == !!sample) %$% gene_id %>%
    unique()

for (sample in conf$samples) {
    p <- merged1 %>%
        left_join(clusterdf) %>%
        filter(sample == !!sample) %>%
        group_by(gene_id) %>%
        mutate(cpgs_detected_per_element = n()) %>%
        ungroup() %>%
        filter(cpgs_detected_per_element > 50) %>%
        group_by(cluster) %>%
        heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
    hms[[sample]] <- p
    dir.create(outputdir, recursive = TRUE)
    mysaveandstore(sprintf("%s/%s_methylation_%s_cluster_stratified.pdf", outputdir, sample, subfam), w = 6, h = 6)
}

clustergrs <- clusterdf %>%
    left_join(rmann) %>%
    GRanges()
cluster_for_gimme_maelstrom %>%
    filter(cluster == "1") %>%
    print(n = 200)
cluster_for_gimme_maelstrom <- as.data.frame(promoters(clustergrs, upstream = 300, downstream = 909)) %>%
    tibble() %>%
    mutate(position_string = paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(position_string, cluster) %>%
    dplyr::rename(loc = position_string, cluster = cluster)

## promoters is giving negative ranges
comb_mat <- combn(cluster_for_gimme_maelstrom$cluster %>% unique(), 2)
for (i in 1:ncol(comb_mat)) {
    cluster1 <- comb_mat[1, i]
    cluster2 <- comb_mat[2, i]
    dir.create("ldna/results/cluster_analysis", recursive = TRUE)
    cluster_path <- sprintf("ldna/results/cluster_analysis/l1_clusters_for_motif_analysis_%s_%s.tsv", cluster1, cluster2)
    cluster_for_gimme_maelstrom %>%
        filter(cluster %in% c(cluster1, cluster2)) %>%
        arrange(cluster) %>%
        write_delim(cluster_path, col_names = TRUE, delim = "\t")
    outputdir_clusters <- "ldna/results/cluster_analysis"
    # motif analysis of consensus
    conda_base_path <- system("conda info --base", intern = TRUE)
    cluster_tf_enrichment_outputdir_path <- sprintf("%s/cluster_analysis_%s_%s", outputdir_clusters, cluster1, cluster2)
    dir.create(cluster_tf_enrichment_outputdir_path, recursive = TRUE)
    motif_pwms_path <- "/oscar/home/mkelsey/anaconda/gimme/lib/python3.10/site-packages/data/motif_databases/HOCOMOCOv11_HUMAN.pfm"

    system(
        sprintf(
            "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme maelstrom  --nthreads 16  %s %s %s -p %s",
            conda_base_path, cluster_path, inputs$ref, cluster_tf_enrichment_outputdir_path, motif_pwms_path
        )
    )
    results <- read_table()
}
