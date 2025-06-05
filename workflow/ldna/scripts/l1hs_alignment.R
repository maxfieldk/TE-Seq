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
        assign("params", list(
            l13 = conf$l13fasta,
            mod_code = "m"
        ), env = globalenv())
        assign("outputs", list(), env = globalenv())
    }
)


library(Rsamtools)
fa <- Rsamtools::FaFile(inputs$ref)


if (confALL$aref$update_ref_with_tldr$response == "yes") {
    if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
        rmannShared <- read_csv(conf$rmann)
        rmannSamples <- list()
        for (sample in sample_table$sample_name) {
            df <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample, sample))
            df$sample_name <- sample
            df <- df %>% left_join(sample_table)
            rmannSamples[[sample]] <- df
        }
        rmannnonref <- do.call(rbind, rmannSamples) %>% tibble()
        rmann <- bind_rows(rmannShared, rmannnonref) %>%
            filter(refstatus != "NonCentral")
    } else if (confALL$aref$update_ref_with_tldr$per_sample == "no") {
        rmann <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann.csv", "A.REF", "A.REF")) %>% filter(refstatus != "NonCentral")
    }
} else {
    rmann <- read_csv(sprintf("aref/extended/%s_annotations/%s_rmann.csv", "A.REF", "A.REF")) %>% filter(refstatus != "NonCentral")
}
##############################
outputdir_already_computed <- "/users/mkelsey/data/LF1/RTE/ldna/results/m/plots/l1_alignment_meth"
outputdir <- "ldna/results/m/plots/l1_alignment_meth"

subfam <- "L1HS"
grs_fl <- rmann %>%
    filter(rte_length_req == "FL") %>%
    filter(rte_subfamily == subfam) %>%
    group_by(rte_subfamily) %>%
    mutate(length99 = quantile(length, 0.99)) %>%
    ungroup() %>%
    filter(length < length99 * 1.05)
grs_fl_ref <- grs_fl %>%
    filter(refstatus == "Ref") %>%
    GRanges()

grs_fl_ss_ref <- getSeq(fa, grs_fl_ref)
names(grs_fl_ss_ref) <- mcols(grs_fl_ref)$gene_id
grs_fl_ss_nonref_list <- list()
if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
    for (sample in sample_table$sample_name) {
        sample_grs_fl <- grs_fl %>%
            filter(sample_name == sample) %>%
            GRanges()
        sample_fa <- Rsamtools::FaFile(sprintf("aref/extended/%s.fa", sample))
        grs_fl_ss_sample <- getSeq(sample_fa, sample_grs_fl)
        names(grs_fl_ss_sample) <- paste0(sample, "___", mcols(sample_grs_fl)$gene_id)
        grs_fl_ss_nonref_list[[sample]] <- grs_fl_ss_sample
    }
    grs_fl_ss_nonref <- purrr::reduce(grs_fl_ss_nonref_list, c)
} else {
    grs_fl_nonref <- grs_fl %>%
        filter(refstatus == "NonRef") %>%
        GRanges()
    grs_fl_ss_nonref <- getSeq(fa, grs_fl_nonref)
    names(grs_fl_ss_nonref) <- mcols(grs_fl_nonref)$gene_id
}
grs_fl_ss <- c(grs_fl_ss_ref, grs_fl_ss_nonref)

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
# alignment_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_alignment_table.csv", outputdir, subfam))
write_csv(consensus_index_long, sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir, subfam))
# consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir, subfam))

##########
consensus_index_long
consensus_index_long %>% filter(sequence_pos <10)
consensus_index_long %>% filter(gene_id == "L1HS_12q21.1_1") %>% pl()
consensus_index_long %>% filter(consensus_pos>10) %>% filter(sequence_pos <10)
consensus_index_long %>% filter(!is.na(consensus_pos)) %>% group_by(gene_id) %>% summarise(minsp = min(consensus_pos, na.rm = TRUE)) %$% minsp %>% table()
consensus_index_long %>% filter(sequence_pos == 1)
consensus_index_long %>% filter(gene_id == "L1HS_10p11.21_4")

consensus_index_long %>% filter(!is.na(sequence_pos)) %>% group_by(gene_id) %>% summarise(minsp = min(consensus_pos, na.rm = TRUE)) %$% minsp %>% table()
consensus_index_long %$% gene_id %>% unique()
consensus_index_long %>% filter(!is.na(sequence_pos)) %>% group_by(gene_id) %>% summarise(minsp = min(consensus_pos, na.rm = TRUE)) %$% minsp %>% table() %>% sum()

consensus_index_long %>% 
    filter(!is.na(sequence_pos)) %>% 
    mutate(m_ind = if_else(seq == consensus_seq, 0, 1)) %>% 
    filter(consensus_pos <21) %>% 
    filter(consensus_pos >= 10) %>% 
    group_by(gene_id) %>% 
    summarise(n = n(), mm = sum(m_ind)) %$% 
    mm %>% table()
perfect_yy1_elements <- consensus_index_long %>% 
    filter(!is.na(sequence_pos)) %>% 
    mutate(m_ind = if_else(seq == consensus_seq, 0, 1)) %>% 
    filter(consensus_pos <=21) %>% 
    filter(consensus_pos >= 11) %>%
    group_by(gene_id) %>% 
    summarise(n = n(), mm = sum(m_ind)) %>%
    filter(mm == 0) %>% 
    filter(n == 11) %$% gene_id
ok_yy1_elements <- consensus_index_long %>% 
    filter(!is.na(sequence_pos)) %>% 
    mutate(m_ind = if_else(seq == consensus_seq, 0, 1)) %>% 
    filter(consensus_pos <=21) %>% 
    filter(consensus_pos >= 11) %>%
    group_by(gene_id) %>% 
    summarise(n = n(), mm = sum(m_ind)) %>%
    filter(mm < 4) %>% 
    filter(n >9) %$% gene_id

p <- perelementdf_promoters %>% filter(rte_subfamily == "L1HS", rte_length_req == "FL") %>% 
mutate(yy1 = if_else(gene_id %in% ok_yy1_elements, "yes", "no")) %>%
ggplot() +
    geom_density(aes(x = mean_meth, group = yy1, color = yy1))
mysaveandstore("zztemp3.pdf")

##############
# METHYLATION CLUSTERING

# outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
# subfam <- "L1HS"
# consensus_index_long <- read_csv(sprintf("%s/%s_fl_mapping_to_consensus_table.csv", outputdir_meth_clustering, subfam))

# rtedf <- read_delim(sprintf("ldna/Rintermediates/%s/rtedf.tsv", params$mod_code), col_names = TRUE)

# consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
# consensus_ss <- readDNAStringSet(consensus_path)

# cg_indices <- consensus_ss %>%
#     vmatchPattern(pattern = "CG") %>%
#     start() %>%
#     unlist() %>%
#     as.numeric()
# consensus_ss[[1]][5762:5763]
# cg_positions_df <- consensus_index_long %>% filter(consensus_pos %in% cg_indices)

# methdf <- rtedf %>% filter(rte_subfamily == subfam)
# mdf <- methdf %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end) - 1))

# senseelement <- mdf %>%
#     filter(rte_strand == "+") %$% gene_id %>%
#     pluck(1)

# antisenseelement <- mdf %>%
#     filter(rte_strand == "-") %$% gene_id %>%
#     pluck(1)

# cpgmapping_check <- cg_positions_df %>%
#     filter(gene_id == senseelement) %$% sequence_pos %>%
#     unique() %>%
#     sort()
# methdf_check <- mdf %>%
#     filter(gene_id == senseelement) %>%
#     relocate(sequence_pos) %$% sequence_pos %>%
#     unique() %>%
#     sort()
# print(methdf_check)
# print(cpgmapping_check)
# cpgmapping_check <- cg_positions_df %>%
#     filter(gene_id == antisenseelement) %$% sequence_pos %>%
#     unique() %>%
#     sort()
# methdf_check <- mdf %>%
#     filter(gene_id == antisenseelement) %>%
#     relocate(sequence_pos) %$% sequence_pos %>%
#     unique() %>%
#     sort()
# print(methdf_check)
# print(cpgmapping_check)

# merged <- left_join(cg_positions_df, mdf, by = c("gene_id", "sequence_pos"))
# cpg_order <- merged %$% consensus_pos %>%
#     unique() %>%
#     sort()
# merged <- merged %>% mutate(consensus_pos = factor(consensus_pos, levels = cpg_order))
# library(tidyHeatmap)

# dat <- merged %>%
#     filter(!is.na(pctM))
# dat %$% gene_id %>%
#     unique() %>%
#     length()
# # write_csv(dat, "meth_for_bayes.csv")


# # methylation heatmaps
# hms <- list()
# for (sample in conf$samples) {
#     p <- merged %>%
#         filter(sample == !!sample) %>%
#         group_by(gene_id) %>%
#         mutate(cpgs_detected_per_element = n()) %>%
#         ungroup() %>%
#         filter(cpgs_detected_per_element > 50) %>%
#         heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE) %>%
#         as_ComplexHeatmap()
#     merged %$% consensus_pos %>% unique()
#     num_pos <- as.numeric(as.character(merged$consensus_pos)) %>% unique()
#     annotvec <- c(
#         rep("5UTR", length(num_pos[num_pos <= 909])),
#         rep("Body", length(num_pos[num_pos > 909]))
#     )
#     ha <- HeatmapAnnotation(anatomy = annotvec)
#     dir.create(outputdir_meth_clustering, recursive = TRUE)
#     mysaveandstore(sprintf("%s/%s_methylation_%s.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 6)
#     hms[[sample]] <- p
#     p <- p %v% ha
#     mysaveandstore(sprintf("%s/%s_methylation_%s_with_annot.pdf", outputdir_meth_clustering, sample, subfam), w = 6, h = 6)
# }
# # Generate the expression as a string and parse it
# p <- base::eval(base::parse(text = paste0("hms[['", conf$samples, "']]", collapse = " + ")))
# mysaveandstore(sprintf("%s/%s_methylation_%s1.pdf", outputdir_meth_clustering, "all", subfam), w = 36, h = 6)
# rm(p)
# p <- Reduce(`+`, hms) # Add all elements in the list

# # motif analysis of consensus
# motif_pwms_path <- "/oscar/home/mkelsey/anaconda/gimme/lib/python3.10/site-packages/data/motif_databases/HOCOMOCOv11_HUMAN.pfm"
# conda_base_path <- system("conda info --base", intern = TRUE)
# motif_search_output_path <- sprintf("%s/alignments/%s_fl_consensus.motifs.bed", outputdir_meth_clustering, subfam)

# system(
#     sprintf(
#         "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme scan %s -g %s -p %s -c 0.9 -n 10 -b > %s",
#         conda_base_path, consensus_path, inputs$ref, motif_pwms_path, motif_search_output_path
#     )
# )
# motifs_df <- as.data.frame(rtracklayer::import(motif_search_output_path)) %>% tibble()


# motif_search_output_path <- sprintf("%s/alignments/%s_fl.motifs.bed", outputdir_meth_clustering, subfam)
# system(
#     sprintf(
#         "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme scan %s -g %s -p %s -c 0.9 -n 10 -b > %s",
#         conda_base_path, fl_fa_path, inputs$ref, motif_pwms_path, motif_search_output_path
#     )
# )

# motifs_all_df <- as.data.frame(rtracklayer::import(motif_search_output_path)) %>%
#     tibble() %>%
#     dplyr::rename(gene_id = seqnames)
# motifs_all_df

# # analysis of methylation by yy1 status
# yy1_pos_consensus <- motifs_df %>% filter(grepl("YY1", name))
# yy1_pos <- motifs_all_df %>% filter(grepl("YY1", name))

# yy1sitesdf <- tibble(gene_id = yy1_pos$gene_id %>% unique(), yy1status = TRUE)
# yy1sitesdf <- tibble(gene_id = merged$gene_id %>% unique()) %>%
#     left_join(yy1sitesdf) %>%
#     mutate(yy1status = replace_na(yy1status, FALSE))


# yy1sitesdfbasedonconsensus <- consensus_index_long %>%
#     filter(consensus_pos < yy1_pos_consensus$end) %>%
#     filter(consensus_pos > yy1_pos_consensus$start) %>%
#     group_by(gene_id) %>%
#     mutate(seq_pos_void = ifelse(is.na(sequence_pos), 1, 0)) %>%
#     summarise(sumna = sum(seq_pos_void)) %>%
#     mutate(yy1status_by_consensus = ifelse(sumna > 2, FALSE, TRUE)) %>%
#     dplyr::select(-sumna)


# merged1 <- merged %>%
#     left_join(yy1sitesdf) %>%
#     left_join(yy1sitesdfbasedonconsensus)
# merged1 %$% yy1status %>% table()
# merged1 %$% yy1status_by_consensus %>% table()

# hms <- list()
# for (sample in conf$samples) {
#     p <- merged1 %>%
#         filter(sample == !!sample) %>%
#         group_by(gene_id) %>%
#         mutate(cpgs_detected_per_element = n()) %>%
#         ungroup() %>%
#         filter(cpgs_detected_per_element > 50) %>%
#         group_by(yy1status) %>%
#         heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
#     hms[[sample]] <- p
#     dir.create(outputdir, recursive = TRUE)
#     mysaveandstore(sprintf("%s/%s_methylation_%s_yy1stratified.pdf", outputdir, sample, subfam), w = 6, h = 6)

#     p <- merged1 %>%
#         filter(sample == !!sample) %>%
#         group_by(gene_id) %>%
#         mutate(cpgs_detected_per_element = n()) %>%
#         ungroup() %>%
#         filter(cpgs_detected_per_element > 50) %>%
#         group_by(yy1status_by_consensus) %>%
#         heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
#     hms[[sample]] <- p
#     dir.create(outputdir, recursive = TRUE)
#     mysaveandstore(sprintf("%s/%s_methylation_%s_yy1consensusstratified.pdf", outputdir, sample, subfam), w = 6, h = 6)
# }



# merged1 %>%
#     filter(!is.na(sample)) %>%
#     mutate(consensus_pos = as.integer(as.character(consensus_pos))) %>%
#     mutate(UTR = ifelse((consensus_pos) < 909, TRUE, FALSE)) %>%
#     mutate(UTRcpgisland = ifelse((consensus_pos) < 500, TRUE, FALSE)) %>%
#     filter(UTRcpgisland == TRUE) %>%
#     group_by(gene_id, sample, condition, yy1status) %>%
#     summarise(pctM = mean(pctM)) %>%
#     group_by(condition, yy1status) %>%
#     summarise(pctM = mean(pctM))




# # CLUSTER ANALYSIS
# matdf <- merged1 %>%
#     dplyr::select(consensus_pos, gene_id, pctM) %>%
#     pivot_wider(names_from = consensus_pos, values_from = pctM)
# rownames(matdf) <- matdf$gene_id
# mat <- as.matrix(matdf %>% dplyr::select(-gene_id))
# rownames(mat) <- matdf$gene_id

# library(impute)
# # Impute missing values in pctM (default is K-nearest neighbors)
# imputed_data <- impute.knn(mat)$data
# # Create a distance matrix from the imputed data
# distance_matrix <- dist(imputed_data, method = "euclidean")
# # Perform clustering
# clusters <- hclust(distance_matrix, method = "ward.D2")
# cluster_group <- cutree(clusters, k = 3)
# clusterdf <- tibble(gene_id = names(cluster_group), cluster = cluster_group)
# clusterdf %>% print(n = 1000)

# merged1 %$% gene_id %>% unique()
# merged1 %>%
#     left_join(clusterdf) %>%
#     filter(sample == !!sample) %>%
#     group_by(gene_id) %>%
#     mutate(cpgs_detected_per_element = n()) %>%
#     ungroup() %$% gene_id %>%
#     unique()
# merged1 %>%
#     filter(sample == !!sample) %$% gene_id %>%
#     unique()

# for (sample in conf$samples) {
#     p <- merged1 %>%
#         left_join(clusterdf) %>%
#         filter(sample == !!sample) %>%
#         group_by(gene_id) %>%
#         mutate(cpgs_detected_per_element = n()) %>%
#         ungroup() %>%
#         filter(cpgs_detected_per_element > 50) %>%
#         group_by(cluster) %>%
#         heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
#     hms[[sample]] <- p
#     dir.create(outputdir, recursive = TRUE)
#     mysaveandstore(sprintf("%s/%s_methylation_%s_cluster_stratified.pdf", outputdir, sample, subfam), w = 6, h = 6)
# }


# # L1HS_10p15.1_2 why all NA??
# # CLUSTER ANALYSIS TOSSING OUT ELEMENTS WITH MORE than 25% missing values.
# gene_ids_to_keep <- merged1 %>%
#     dplyr::select(consensus_pos, gene_id, pctM) %>%
#     group_by(gene_id) %>%
#     summarise(prop_cpg_no_data = sum(is.na(pctM) / n()))
# gene_ids_to_keep %$% prop_cpg_no_data %>% table()
# matdf <- merged1 %>%
#     left_join(gene_ids_to_keep) %>%
#     filter(number_cpg_no_data < 0.3) %>%
#     dplyr::select(consensus_pos, gene_id, pctM) %>%
#     pivot_wider(names_from = consensus_pos, values_from = pctM)
# rownames(matdf) <- matdf$gene_id
# mat <- as.matrix(matdf %>% dplyr::select(-gene_id))
# rownames(mat) <- matdf$gene_id

# library(impute)
# # Impute missing values in pctM (default is K-nearest neighbors)
# imputed_data <- impute.knn(mat)$data
# # Create a distance matrix from the imputed data
# distance_matrix <- dist(imputed_data, method = "euclidean")
# # Perform clustering
# clusters <- hclust(distance_matrix, method = "ward.D2")
# cluster_group <- cutree(clusters, k = 3)
# clusterdf <- tibble(gene_id = names(cluster_group), cluster = cluster_group)
# clusterdf %>% print(n = 1000)


# for (sample in conf$samples) {
#     p <- clusterdf %>%
#         left_join(merged1) %>%
#         filter(sample == !!sample) %>%
#         group_by(gene_id) %>%
#         mutate(cpgs_detected_per_element = n()) %>%
#         ungroup() %>%
#         filter(cpgs_detected_per_element > 50) %>%
#         group_by(cluster) %>%
#         heatmap(gene_id, consensus_pos, pctM, cluster_rows = TRUE, cluster_columns = FALSE)
#     hms[[sample]] <- p
#     dir.create(outputdir, recursive = TRUE)
#     mysaveandstore(sprintf("%s/%s_methylation_%s_cluster_stratified_only_elements_used_in_DE_TF_analysis.pdf", outputdir, sample, subfam), w = 6, h = 6)
# }


# clustergrs <- clusterdf %>%
#     left_join(rmann) %>%
#     GRanges()
# cluster_for_gimme_maelstrom <- as.data.frame(promoters(clustergrs, upstream = 300, downstream = 909)) %>%
#     tibble() %>%
#     mutate(position_string = paste0(seqnames, ":", start, "-", end)) %>%
#     dplyr::select(position_string, cluster) %>%
#     dplyr::rename(loc = position_string, cluster = cluster)

# cluster_for_homer <- promoters(clustergrs, upstream = 300, downstream = 909)
# ## promoters is giving negative ranges
# comb_mat <- combn(cluster_for_gimme_maelstrom$cluster %>% unique(), 2)
# for (i in 1:ncol(comb_mat)) {
#     cluster1 <- comb_mat[1, i]
#     cluster2 <- comb_mat[2, i]
#     dir.create(sprintf("ldna/results/%s/cluster_analysis", params$mod_code), recursive = TRUE)
#     cluster_path <- sprintf("ldna/results/%s/cluster_analysis/l1_clusters_for_motif_analysis_%s_%s.tsv", params$mod_code, cluster1, cluster2)
#     cluster_for_gimme_maelstrom %>%
#         filter(cluster %in% c(cluster1, cluster2)) %>%
#         arrange(cluster) %>%
#         write_delim(cluster_path, col_names = TRUE, delim = "\t")
#     outputdir_clusters <- sprintf("ldna/results/%s/cluster_analysis", params$mod_code)
#     # motif analysis of consensus
#     conda_base_path <- system("conda info --base", intern = TRUE)
#     cluster_tf_enrichment_outputdir_path_gimme <- sprintf("%s/cluster_analysis_%s_%s/gimme", outputdir_clusters, cluster1, cluster2)
#     cluster_tf_enrichment_outputdir_path_homer <- sprintf("%s/cluster_analysis_%s_%s/homer", outputdir_clusters, cluster1, cluster2)

#     dir.create(cluster_tf_enrichment_outputdir_path_gimme, recursive = TRUE)
#     dir.create(cluster_tf_enrichment_outputdir_path_homer, recursive = TRUE)

#     motif_pwms_path <- "/oscar/home/mkelsey/anaconda/gimme/lib/python3.10/site-packages/data/motif_databases/HOCOMOCOv11_HUMAN.pfm"

#     system(
#         sprintf(
#             "source %s/etc/profile.d/conda.sh && conda activate gimme && gimme maelstrom  --nthreads 16  %s %s %s -p %s",
#             conda_base_path, cluster_path, inputs$ref, cluster_tf_enrichment_outputdir_path_gimme, motif_pwms_path
#         )
#     )
#     # results <- read_table(sprintf("%s/"))
#     c1fapath <- sprintf("%s/cluster%s.fa", cluster_tf_enrichment_outputdir_path_homer, cluster1)
#     c2fapath <- sprintf("%s/cluster%s.fa", cluster_tf_enrichment_outputdir_path_homer, cluster2)
#     writeXStringSet(
#         getSeq(fa, cluster_for_homer[mcols(cluster_for_homer)$cluster == cluster1]),
#         c1fapath
#     )
#     writeXStringSet(
#         getSeq(fa, cluster_for_homer[mcols(cluster_for_homer)$cluster == cluster2]),
#         c2fapath
#     )


#     system(
#         sprintf(
#             "findMotifs.pl %s fasta %s -fasta %s > %s",
#             c1fapath, sprintf("%s/foreground%s_background%s", cluster_tf_enrichment_outputdir_path_homer, cluster1, cluster2), c2fapath, sprintf("%s/outputfile.txt", cluster_tf_enrichment_outputdir_path_homer)
#         )
#     )
#     system(
#         sprintf(
#             "findMotifs.pl %s fasta %s -fasta %s > %s",
#             c2fapath, sprintf("%s/foreground%s_background%s", cluster_tf_enrichment_outputdir_path_homer, cluster2, cluster1), c1fapath, sprintf("%s/outputfile.txt", cluster_tf_enrichment_outputdir_path_homer)
#         )
#     )
# }
