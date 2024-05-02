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
library(rtracklayer)
library(Biostrings)
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
            vcf = sprintf("ldna/intermediates/%s/snpeff/snpeff.pass.vcf", sample_table$sample_name)
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


vcf <- read.vcfR(inputs$vcf)
Z <- vcfR2tidy(vcf)
df <- Z$fix%>% mutate(REF_LEN = nchar(REF), ALT_LEN = nchar(ALT)) %>% mutate(TYPE = ifelse(REF_LEN == 1 & ALT_LEN == 1, "SNP", ifelse(REF_LEN != ALT_LEN, "INDEL", "MIXED")))
df  %>% group_by(TYPE) %>% summarise(n = n())

snp <- df %>% filter(TYPE == "SNP" | TYPE == "MIXED")
indel <- df %>% filter(TYPE == "INDEL")



snp_lof_genes <- snp %>% filter(!is.na(LOF)) %$% LOF %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
snp_nmd_genes <- snp %>% filter(!is.na(NMD)) %$% NMD %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
snp_lof_genes[snp_lof_genes %in% putative_l1_regulators]

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
for (soi in sets_of_interest) {
    print(soi)
    genes <- gs %>% filter(gs_name == soi) %$% gene_symbol
    snp_lof_intersection[[soi]] <- snp_lof_genes[snp_lof_genes %in% genes]
    indel_lof_intersection[[soi]] <- indel_lof_genes[indel_lof_genes %in% genes]

}



indel_lof_genes <- indel %>% filter(!is.na(LOF)) %$% LOF %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
indel_nmd_genes <- indel %>% filter(!is.na(NMD)) %$% NMD %>% gsub("\\(", "", .) %>% gsub("\\)", "", .) %>% str_split("\\|") %>% sapply(head, 1)
indel_lof_genes[indel_lof_genes %in% putative_l1_regulators]
#look for L1 expression modulator dataset


gdf <- read_table("/users/mkelsey/data/n2102ep/RTE/ldna/intermediates/N2102EP1/snpeff/snpeffStats.genes.txt", comment = "# ") %>% dplyr::rename(gene_id = "GeneId")
highimpact <- gdf %>% filter(variants_impact_HIGH > 0)



highimpact %>% filter(gene_id %in% putative_l1_regulators) %>% pull(gene_id) %>% unique()
gdf %>% filter(gene_id %in% putative_l1_regulators) 

act_mut <- gdf %>% filter(gene_id %in% putative_l1_activators)
sup_mut <- gdf %>% filter(gene_id %in% putative_l1_suppressors) 

colnames(sup_mut)
c("variants_effect_start_lost", "variants_effect_stop_gained", "variants_effect_stop_lost")
sup_mut %>% filter() %>% pull(gene_id) %>% unique()
act_mut %>% filter(variants_effect_start_lost > 0 | variants_effect_stop_gained > 0 | variants_effect_stop_lost > 0) %>% pull(gene_id) %>% unique()

sup_mut %>% filter(variants_impact_HIGH > 0 ) %>% pull(gene_id) %>% unique()
act_mut %>% filter(variants_impact_HIGH > 0 ) %>% pull(gene_id) %>% unique()

