source("workflow/scripts/defaults.R")
module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(vcfR)
library(msigdbr)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            SV = sprintf("ldna/intermediates/%s/snpeff/snpeff_sv.pass.vcf", sample_table$sample_name),
            SML = sprintf("ldna/intermediates/%s/snpeff/snpeff_sml.pass.vcf", sample_table$sample_name)        
        ), env = globalenv())
        assign("outputs", list(
            plots = "ldna/results/plots/variation/vcf_plots.rds"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

cdl1 <- read_delim(conf$candidate_l1_regulators) %>%
    filter((`Minimum Effect 95% Credible Interval Estimate` < 0 & `Maximum Effect 95% Credible Interval Estimate` < 0) | (`Minimum Effect 95% Credible Interval Estimate` > 0 & `Maximum Effect 95% Credible Interval Estimate` > 0)) %>%
    dplyr::rename(gene_id = Symbol) %>% 
    mutate(type = ifelse(`Minimum Effect 95% Credible Interval Estimate` > 0 , "suppresor", "activator"))
putative_l1_regulators <- cdl1 %>% pull(gene_id)
putative_l1_activators <- cdl1 %>% filter(type == "activator") %>% pull(gene_id)
putative_l1_suppressors <- cdl1 %>% filter(type == "suppresor") %>% pull(gene_id)

omim <- read_delim(conf$omim, comment = "#")
colnames(omim)
omim_genes_string <- omim %>% dplyr::select( "Gene/Locus And Other Related Symbols" , Phenotypes) %>% filter(!str_detect(Phenotypes, "^Chromosome ")) %$%  `Gene/Locus And Other Related Symbols` 
omim_genes <- str_split(omim_genes_string, ", ") %>% unlist() %>% unique()

#data preprocessing
dflist <- list()
vcfdf <- list()
dfSVlist <- list()
for (sample in sample_table$sample_name) {

SV <- read.vcfR(inputs$SV)
tidyvcf_SV <- vcfR2tidy(SV)
dfSV <- tidyvcf_SV$fix %>% filter(FILTER == "PASS") %>% mutate(sample_name = sample)

vcf <- read.vcfR(inputs$snpeff)
tidyvcf <- vcfR2tidy(vcf)
df <- tidyvcf$fix%>% mutate(REF_LEN = nchar(REF), ALT_LEN = nchar(ALT)) %>%
    mutate(TYPE = ifelse(REF_LEN == 1 & ALT_LEN == 1, "SNP", ifelse(REF_LEN != ALT_LEN, "INDEL", "MIXED"))) %>% 
    mutate(sample_name = sample)
df <- df %>% cbind(tidyvcf$gt %>% dplyr::select(-c("ChromKey", "POS"))) %>% tibble()

dflist[[sample]] <- df
dfSVlist[[sample]] <- dfSV
vcfdf[[sample]] <- tidyvcf
}

df <- Reduce(bind_rows, dflist)
dfSV <- Reduce(bind_rows, dfSVlist)

snp <- df %>% filter(QUAL >= 20) %>% filter(TYPE == "SNP" | TYPE == "MIXED") %>% mutate(from_to = paste(REF, ALT, sep = "->")) %>% mutate(from_to_compressed = case_when(
    from_to == "A->G" ~ "T->C",
    from_to == "A->C" ~ "T->G",
    from_to == "A->T" ~ "T->A",
    from_to == "G->A" ~ "C->T",
    from_to == "G->C" ~ "C->G",
    from_to == "G->T" ~ "C->A",
    TRUE ~ "Other"
    ) %>% factor(levels = c("T->C", "T->A", "T->G", "C->T", "C->A","C->G", "Other"))
    ) %>% mutate(snp_type = case_when(
        from_to_compressed %in% c("T->C", "C->T") ~ "Transition",
        from_to_compressed %in% c("T->G", "T->A", "C->G", "C->A") ~ "Transversion",
        TRUE ~ "Other"
    ) %>% factor(levels = c("Transition", "Transversion", "Other"))) %>%
    mutate(ANN = str_split(ANN, "\\|")) %>% 
    mutate(region = sapply(ANN, function(x) x[2])) %>%
    mutate(severity = sapply(ANN, function(x) x[3])) %>% 
    mutate(gene_id = ifelse(region != "intergenic_region", sapply(ANN, function(x) x[4]), NA))

indel <- df %>% filter(QUAL >= 20) %>% filter(TYPE == "INDEL") %>%
    mutate(indel_type = ifelse(REF_LEN > ALT_LEN, "Deletion", "Insertion")) %>%
    mutate(ANN = str_split(ANN, "\\|")) %>% 
    mutate(region = sapply(ANN, function(x) x[2])) %>%
    mutate(severity = sapply(ANN, function(x) x[3])) %>% 
    mutate(gene_id = ifelse(region != "intergenic_region", sapply(ANN, function(x) x[4]), NA))

#SNVs
library(ggpubr)
p <- snp %>% group_by(from_to_compressed, snp_type, sample_name) %>% summarise(n = n()) %>% filter(snp_type != "Other") %>% 
    ggbarplot(x = "from_to_compressed", y = "n", fill = "snp_type", xlab = "Variant Type", ylab = "Count", facet = "sample_name") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/sm_variant_type.png", outputdir))


p <- snp %>% group_by(snp_type, sample_name) %>% summarise(n = n()) %>% filter(snp_type != "Other") %>% ungroup() %>% pivot_wider(names_from = snp_type, values_from = n) %>% group_by(sample_name) %>% summarize(tstv = Transition / Transversion) %>%
    ggbarplot(x = "sample_name", y = "tstv", xlab = "", ylab = "Ts/Tv", fill = "sample_name") +
    scale_samples +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/snp_tstv.png", outputdir))


p <- snp %>%  group_by(region, sample_name) %>% summarise(log2n = log2(n())) %>% 
    ggscatter(x = "region", y = "log2n", xlab = "", ylab = "Log2 Count", color = "sample_name") +
    scale_samples +
    mtopen + coord_flip()
mysaveandstore(sprintf("%s/snp_region.png", outputdir), 8, 5)

p <- snp %>%  group_by(severity, sample_name) %>% summarise(log2n = log2(n())) %>% 
    ggscatter(x = "severity", y = "log2n", xlab = "", ylab = "Log2 Count", color = "sample_name") +
    scale_samples +
    mtopen + coord_flip()
mysaveandstore(sprintf("%s/snp_severity.png", outputdir))

p <- snp %>%  
    gghistogram(x = "QUAL", xlab = "", ylab = "Count", facet = "sample_name", fill = "sample_name") +
    scale_samples +
    anchorbar +
    mtclosed
mysaveandstore(sprintf("%s/snp_quality.png", outputdir))


p <- indel %>% group_by(indel_type, sample_name) %>% summarise(n = n()) %>% filter(indel_type != "Other") %>% 
    ggbarplot(x = "indel_type", y = "n", fill = "indel_type", xlab = "Variant Type", ylab = "Count", facet = "sample_name") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/indel_variant_type.png", outputdir))

p <- indel %>%  group_by(region, sample_name) %>% summarise(log2n = log2(n())) %>% 
    ggscatter(x = "region", y = "log2n", xlab = "", ylab = "Log2 Count", color = "sample_name") +
    scale_samples +
    mtopen + coord_flip()
mysaveandstore(sprintf("%s/indel_region.png", outputdir), 12, 5)

p <- indel %>%  group_by(severity, sample_name) %>% summarise(log2n = log2(n())) %>% 
    ggscatter(x = "severity", y = "log2n", xlab = "", ylab = "Log2 Count", color = "sample_name") +
    scale_samples +
    mtopen + coord_flip()
mysaveandstore(sprintf("%s/indel_severity.png", outputdir))

p <- indel %>%  
    gghistogram(x = "QUAL", xlab = "", ylab = "Count", facet = "sample_name", fill = "sample_name") +
    scale_samples +
    anchorbar +
    mtclosed 
mysaveandstore(sprintf("%s/indel_quality.png", outputdir))




snp_lof_genes <- snp %>% filter(!is.na(LOF)) %$% LOF %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
snp_nmd_genes <- snp %>% filter(!is.na(NMD)) %$% NMD %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
snp_lof_genes[snp_lof_genes %in% putative_l1_regulators]

indel_lof_genes <- indel %>% filter(!is.na(LOF)) %$% LOF %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
indel_nmd_genes <- indel %>% filter(!is.na(NMD)) %$% NMD %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
indel_lof_genes[indel_lof_genes %in% putative_l1_regulators]

sv_lof_genes <- dfSV %>% filter(!is.na(LOF)) %$% LOF %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
sv_nmd_genes <- dfSV %>% filter(!is.na(NMD)) %$% NMD %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
sv_lof_genes[sv_lof_genes %in% putative_l1_regulators]

lof_genes <- c(snp_lof_genes, indel_lof_genes, sv_lof_genes) %>% unique()



library(RClusterProfiler)





gs <- msigdbr("human")

sets_of_interest <- c("KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY",
"REACTOME_INTERFERON_SIGNALING",
"REACTOME_VIRAL_INFECTION_PATHWAYS",
"REACTOME_VIRAL_MESSENGER_RNA_SYNTHESIS",
"KEGG_P53_SIGNALING_PATHWAY",
"REACTOME_CELLULAR_SENESCENCE",
"KEGG_PATHWAYS_IN_CANCER",
"REACTOME_INNATE_IMMUNE_SYSTEM",
"REACTOME_CYTOSOLIC_SENSORS_OF_PATHOGEN_A
SSOCIATED_DNA",
"WP_PATHWAYS_OF_NUCLEIC_ACID_METABOLISM_AND_INNATE_IMMUNE_SENSING",
"REACTOME_METABOLISM_OF_CARBOHYDRATES")


snp_lof_intersection <- list()
indel_lof_intersection <- list()
sv_lof_intersection <- list()
snp_lof_intersection_pasted <- list()
indel_lof_intersection_pasted <- list()
sv_lof_intersection_pasted <- list()
for (soi in sets_of_interest) {
    print(soi)
    genes <- gs %>% filter(gs_name == soi) %$% gene_symbol
    snp_lof_intersection[[soi]] <- snp_lof_genes[snp_lof_genes %in% genes] %>% unique()
    snp_lof_intersection_pasted[[soi]] <- paste(snp_lof_intersection[[soi]], collapse = ", ")
    indel_lof_intersection[[soi]] <- indel_lof_genes[indel_lof_genes %in% genes] %>% unique()
    indel_lof_intersection_pasted[[soi]] <- paste(indel_lof_intersection[[soi]], collapse = ", ")
    sv_lof_intersection[[soi]] <- sv_lof_genes[sv_lof_genes %in% genes] %>% unique()
    sv_lof_intersection_pasted[[soi]] <- paste(sv_lof_intersection[[soi]], collapse = ", ")
}

snp_lof_df <- snp_lof_intersection_pasted %>% as.data.frame() %>% tibble() %>% pivot_longer(everything(), names_to = "Gene Set", values_to = "Geneid") %>% mutate(mut_type = "SNP")
indel_lof_df <- indel_lof_intersection_pasted %>% as.data.frame() %>% tibble() %>% pivot_longer(everything(), names_to = "Gene Set", values_to = "Geneid") %>% mutate(mut_type = "INDEL")
sv_lof_df <- sv_lof_intersection_pasted %>% as.data.frame() %>% tibble() %>% pivot_longer(everything(), names_to = "Gene Set", values_to = "Geneid") %>% mutate(mut_type = "SV")
lof_df <- bind_rows(snp_lof_df, indel_lof_df, sv_lof_df) %>% pivot_wider(names_from = mut_type, values_from = Geneid) 
table_theme <- ttheme(
    base_style = "light",
  tbody.style = tbody_style(
   hjust = 0,
   x = 0.01
 ))
subtitle <- c("Curated genesets")
title <- c("LOF Genes")
p <- ggtexttable(lof_df, rows = NULL, theme = table_theme ) %>% 
    tab_add_title(text = subtitle, face = "plain", size = 10) %>%
    tab_add_title(text = title, face = "bold", padding = unit(0.1, "line"))
mysaveandstore(sprintf("%s/lof_genes.png", outputdir), 24, 4)


