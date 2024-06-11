module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")

library(knitr)
library(rmarkdown)
library(circlize)
library(ComplexHeatmap)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggVennDiagram")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)
library(rtracklayer)
library(ComplexUpset)                            


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("params", list(
            "outputdir" = "lrna/results/agg/repeatanalysis/relaxed",
            "counttype" = "relaxed"
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "lrna/results/agg/deseq/resultsdf.tsv",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
                assign("outputs", list(
            "plots" = "lrna/results/agg/repeatanalysis/relaxed/repeatanalysis_plots.RData"
        ), env = globalenv())
    }
)

outputdir <- params$outputdir
contrasts <- conf$contrasts
sample_table <- read.csv(conf$sample_table)
counttype <- params$counttype

## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t") %>% filter(counttype == params$counttype)
r_annotation_fragmentsjoined <- read_csv(inputs$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(inputs$r_repeatmasker_annotation)
resultsdfwithgenes <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)
resultsdf <- resultsdfwithgenes %>% filter(gene_or_te != "gene")
refseq <- import(conf$annotation_genes)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
noncoding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NR", mcols(refseq)$transcript_id))]
transcripts <- c(coding_transcripts, noncoding_transcripts)


### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
    # ontologies <- ontologies[ontologies %in% conf$repeat_ontologies_to_scrutinize]
# ontologies <- ontologies[ontologies %in% conf$repeat_ontologies_to_scrutinize]
    ontologies <- c("rte_subfamily_limited")


    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        group <- resultsdf %>%
            pull(!!sym(ontology)) %>%
            unique()
        group <- group[!(group == "Other")]
        big_ontology_groups <- c(big_ontology_groups, group)
    }
    big_ontology_groups <- big_ontology_groups %>% unique()

    modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
    region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
    element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]

}

# subset <- resultsdf %>% filter(grepl("*tact*", intactness_req))
for (ontology in c("rte_subfamily_limited", "l1_subfamily_limited", "rte_family")) {
    print(ontology)

subset <- resultsdfwithgenes %>% filter(!!sym(ontology) != "Other") %>% filter(!is.na(!!sym(ontology)))
query <- subset %>% GRanges()
subject <- coding_transcripts
hits <- GenomicRanges::nearest(query, subject)

queryIDs <- mcols(query)$gene_id
subjectIDs <- mcols(subject)$gene_id[hits]
subset$nearest_gene <- subjectIDs[queryIDs]


te_gene_matrix_list <- list()
for (contrast in contrasts) {
    contrast_of_interest <- contrast
    contrast_level_2 <- contrast_of_interest %>%
        gsub("condition_", "", .) %>%
        gsub("_vs_.*", "", .)
    contrast_level_1 <- contrast_of_interest %>%
        gsub(".*_vs_", "", .)
    contrast_stat <- paste0("stat_", contrast_of_interest)
    contrast_padj <- paste0("padj_", contrast_of_interest)
    contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)
    contrast_samples <- sample_table %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)

resultsdfwithgenes[match(subjectIDs, resultsdfwithgenes$gene_id),] %>% pull(!!sym(contrast_log2FoldChange)) %>% length()
te_gene_matrix <- subset %>% dplyr::select(!!sym(ontology), !!sym(contrast_log2FoldChange),loc_integrative,req_integrative ) %>% dplyr::rename(TE = !!sym(contrast_log2FoldChange)) 
te_gene_matrix$GENE <- resultsdfwithgenes[match(subjectIDs, resultsdfwithgenes$gene_id),] %>% pull(!!sym(contrast_log2FoldChange))

te_gene_matrix <- te_gene_matrix %>% drop_na()
te_gene_matrix_list[[contrast]] <- te_gene_matrix
cor.test(te_gene_matrix$TE, te_gene_matrix$GENE, method = "spearman", )$estimate
cor.test(te_gene_matrix$TE, te_gene_matrix$GENE, method = "spearman", )$p.value
te_gene_matrix <- te_gene_matrix %>% drop_na()
cor_df <- te_gene_matrix %>% mutate(req_integrative = gsub(".*Intact.*", "Full Length", req_integrative)) %>% group_by(!!sym(ontology), req_integrative, loc_integrative) %>% mutate(groupN = n()) %>% filter(groupN > 4) %>% summarise(cor = cor.test(TE, GENE, method = "spearman")$estimate, pval = cor.test(TE, GENE, method = "spearman")$p.value)
cor_df %$% req_integrative %>% unique()
pf <- cor_df %>% ungroup() %>% filter(loc_integrative != "Centromeric") %>% complete(!!sym(ontology), req_integrative, loc_integrative)
p<-pf %>% mutate(genicfacet = ifelse(loc_integrative == "Intergenic", "", "Genic")) %>% ggplot() + geom_tile(aes(x = !!sym(ontology), y = loc_integrative, fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/%s/rte_genic_cor_%s.png", outputdir, counttype, contrast, ontology), 6, 4)

p <- pf %>% mutate(genicfacet = ifelse(loc_integrative == "Intergenic", "", "Genic")) %>% ggplot(aes(x = !!sym(ontology), y = loc_integrative)) + geom_tile(aes(fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", "")), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/%s/rte_genic_cor_pval_%s.png", outputdir, counttype, contrast, ontology), 6, 4)
}

te_gene_matrix_all <- Reduce(bind_rows, te_gene_matrix_list)
cor_df <- te_gene_matrix_all %>% mutate(req_integrative = gsub(".*Intact.*", "Full Length", req_integrative)) %>% group_by(!!sym(ontology), req_integrative, loc_integrative) %>% mutate(groupN = n()) %>% filter(groupN > 4) %>% summarise(cor = cor.test(TE, GENE, method = "spearman")$estimate, pval = cor.test(TE, GENE, method = "spearman")$p.value)
cor_df %$% req_integrative %>% unique()
pf <- cor_df %>% ungroup() %>% filter(loc_integrative != "Centromeric") %>% complete(!!sym(ontology), req_integrative, loc_integrative)

p<-pf %>% mutate(genicfacet = ifelse(loc_integrative == "Intergenic", "", "Genic")) %>% ggplot() + geom_tile(aes(x = !!sym(ontology), y = loc_integrative, fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/%s/rte_genic_cor_%s.png", outputdir, counttype, "pan_contrast", ontology), 8, 4)

p <- pf %>% mutate(genicfacet = ifelse(loc_integrative == "Intergenic", "", "Genic")) %>% ggplot(aes(x = !!sym(ontology), y = loc_integrative)) + geom_tile(aes(fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", "")), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/%s/rte_genic_cor_pval_%s.png", outputdir, counttype, "pan_contrast", ontology), 8, 4)

}


#### PLOTTING FUNCTIONS

pvp <- function(df, facet_var = "ALL", filter_var = "ALL", labels = "no", scale_log2 = "no") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), "FL$|Intact|towards"))
    }
if (scale_log2 == "no") {
    pf <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) 
} else {
    pf <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
        mutate(across(conf$levels, ~ log2(.x + 1)))
    }

    top_sig <- pf %>% filter(!!sym(contrast_padj) < 0.05) %>% arrange(padj) %>% head(6) %>% pull(gene_id) 
if (labels != "no") {
    p <- pf %>%
        {
            ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
                geom_point(aes(color = loc_integrative, shape = ifelse(is.na(padj), FALSE, ifelse(padj < 0.05, TRUE, FALSE)))) +
                scale_shape_manual(values = c(1, 19, 3)) +
                scale_size_manual(values = c(1,2,1)) +
labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2), 
                    caption = counttype_label,
                    color = "Genomic Context",
                    shape = "padj < 0.05") +
                mtclosed + scale_palette_alt +
                theme(aspect.ratio = 1) +
                coord_cartesian(clip = "off") +
                geom_text_repel(data = . %>% mutate(label = ifelse(gene_id %in% top_sig, gene_id, "")), aes(label = label), max.overlaps = Inf, show.legend = FALSE) +
                coord_fixed(
                    xlim = range(c(.[[contrast_level_1]], .[[contrast_level_2]])),
                    ylim = range(c(.[contrast_level_1], .[contrast_level_2]))
                ) + guides(fill=guide_legend(title="Genomic Context"))
         }
    } else {

    p <- pf %>%
        {
            ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
                geom_point(aes(color = loc_integrative, shape = ifelse(is.na(padj), FALSE, ifelse(padj < 0.05, TRUE, FALSE)))) +
                scale_shape_manual(values = c(1, 19, 3)) +
                scale_size_manual(values = c(1,2,1)) +
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2), 
                    caption = counttype_label,
                    color = "Genomic Context",
                    shape = "padj < 0.05") +
                mtclosed + scale_palette_alt +
theme(aspect.ratio = 1) +
                coord_cartesian(clip = "off") +
                coord_fixed(
                    xlim = range(c(.[[contrast_level_1]], .[[contrast_level_2]])),
                    ylim = range(c(.[contrast_level_1], .[contrast_level_2]))
                )
}

        }
    if (facet_var != "ALL") {
        p <- p + facet_wrap(facet_var)
    }
if (scale_log2 == "yes") {
        p <- p + labs(x = sprintf("log2(%s Norm Counts + 1)", contrast_level_1), y = sprintf("log2(%s Norm Counts + 1)", contrast_level_2))
    }
    return(p)
}
# df <- tidydf %>% filter(rte_subfamily == "L1HS")
# p <- pvp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "ALL", facet_var = "genic_loc", scale_log2 = "yes") + ggtitle("L1HS")
# mysave("temp1.png", 8, 8)

dep <- function(df, facet_var = "ALL", filter_var = "ALL") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), "FL$|Intact"))
    }
    if (facet_var == "ALL") {
        facet_var_1 <- NULL
    } else {
        facet_var_1 <- facet_var
    }
    vars_for_facet <- c("direction", "n", facet_var_1)
    plotframe <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
        mutate(direction = ifelse(padj < 0.05, ifelse(!!sym(contrast_level_2) > !!sym(contrast_level_1), "UP", "DOWN"), "NS")) %>%
        mutate(direction = ifelse(is.na(direction), "NS", direction)) %>%
        mutate(direction = direction %>% factor(levels = c("UP", "DOWN", "NS"))) %>%
        ungroup() %>%
        mutate(n = n()) %>%
        group_by(across(all_of(vars_for_facet))) %>%
        summarise(count = n()) %>%
        mutate(fraction = count / n) %>%
        filter(direction != "NS")
    if (facet_var == "ALL") {
        p <- plotframe %>% ggplot() +
            labs(x = "", y = "Number DE", caption = counttype_label) +
            geom_col(aes(x = direction, group = direction, y = count)) +
mtclosed +
            scale_palette +
            anchorbar +
                        guides(fill = "none")
    } else {
        p <- plotframe %>% ggplot() +
            labs(x = "", y = "Number DE", caption = counttype_label) +
            geom_col(aes(x = direction, group = direction, y = count)) +
            facet_wrap(facet_var) +
mtclosed +
            scale_palette +
            anchorbar +
                        guides(fill = "none")
    }
    return(p)
}

# need to wrap functions cALLs in try catch since filtering and checking for DE status can mean there are zero elements
# p <- dep(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "rte_length_req", facet_var = "genic_loc") + ggtitle("L1HS")
# mysave("temp1.png")
# p <- dep(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "rte_length_req") + ggtitle("L1HS")
# mysave("temp1.png")

stripp <- function(df, stats = "no", extraGGoptions = NULL, facet_var = "ALL", filter_var = "ALL") {
        if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), "FL$|Intact"))
    }
        if (facet_var != "ALL") {
        pf <- df %>%
            group_by(sample, !!sym(facet_var)) %>%
            summarise(sample_sum = sum(counts), condition = dplyr::first(condition))
            } else {
        pf <- df %>%
            group_by(sample) %>%
            summarise(sample_sum = sum(counts), condition = dplyr::first(condition))
        }
    p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", add = c("mean_se", "dotplot")) +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label) +
                theme(legend.position = "none") +
        mtclosedgridh +
        scale_conditions +
anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
    if (facet_var != "ALL") {
        p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = paste0(facet_var), add = c("mean_se", "dotplot")) +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label) +
        theme(legend.position = "none") +
        mtclosedgridh +
scale_conditions +
        anchorbar +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
    } else {
        p <- pf %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", add = c("mean_se", "dotplot")) +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label) +
        theme(legend.position = "none") +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
    }
    if (stats == "yes") {
        p <- p + geom_pwc(
    method = "t_test", label = "p.adj.signif",
    ref.group = conf$levels[1],
    p.adjust.method = "fdr", hide.ns = TRUE
  )
    }
    return(p)
}

# df <- tidydf %>% filter(rte_subfamily == "L1HS")
# p <- stripp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "ALL", facet_var = "genic_loc", stats = "yes") + ggtitle("L1HS")
# mysave("temp1.png")


myheatmap <- function(df, facet_var = "ALL", filter_var = "ALL", DEvar = "ALL", scaled = "notscaled", contrast_samples, condition_vec) {
    set_title <- group
    if (filter_var != "ALL") {
        if (str_detect(filter_var, "length_req")) {
            df <- df %>% filter(str_detect(!!sym(filter_var), "FL$"))
            set_title <- df %>%
                pull(!!sym(filter_var)) %>%
                unique()
        }
        if (str_detect(filter_var, "intact")) {
            df <- df %>% filter(str_detect(!!sym(filter_var), "Intact"))
        }
        set_title <- df %>%
            pull(!!sym(filter_var)) %>%
            unique()
    }
    show_row_names <- FALSE

    group_res <- df

    if (DEvar == "DE") {
        group_res <- df %>%
            filter((!!sym(contrast_padj)) < 0.05)
if (nrow(group_res) > 50) {
            show_row_names <- FALSE
        } else {
        show_row_names <- TRUE
}
        
    }
    m <- group_res %>%
        dplyr::select(contrast_samples) %>%
        as.matrix()
    rownames(m) <- group_res %$% gene_id

    if (scaled == "scaled") {
        m <- t(scale(t(m))) %>% na.omit()
        group_res <- group_res %>% filter(gene_id %in% rownames(m))
    }
group_res %$% gene_id %>% duplicated()
    split_annot <- group_res %>%
        filter(gene_id %in% rownames(m)) %>%
        mutate(split_annot = ifelse(!!sym(contrast_padj) < 0.05, ifelse(!!sym(contrast_log2FoldChange) > 0, "UP DE", "DOWN DE"), "NOT DE")) %>%
        mutate(split_annot = factor(split_annot, levels = c("UP DE", "NOT DE", "DOWN DE"))) %>%
        pull(split_annot)
    split_annot[is.na(split_annot)] <- "NOT DE"

    topAnn <- HeatmapAnnotation(Condition = condition_vec, col = list(Condition = unlist(condition_palette[condition_vec])))
    colors_for_de <- c("UP DE" = "red", "DOWN DE" = "blue", "NOT DE" = "lightgray")
    if (facet_var != "ALL") {
        facet_var_values <- group_res %>%
            pull(!!sym(facet_var)) %>%
            unique()
        facet_var_values_len <- length(facet_var_values)
        colors_for_facet <- c("#393939", "lightgray", "yellow", "green")
        row_ha <- rowAnnotation(Loc = group_res %>% pull(!!sym(facet_var)), DE = split_annot, col = list(DE = colors_for_de, Loc = setNames(colors_for_facet[1:facet_var_values_len], facet_var_values)))
        split_annot_df <- data.frame(de = split_annot, facet_var = group_res %>% pull(!!sym(facet_var)))
        hm <- m %>%
            Heatmap(
                name = "Normalized Counts",
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = show_row_names,
                show_parent_dend_line = FALSE,
                show_row_dend = FALSE,
                show_column_names = TRUE,
                column_names_rot = 90,
                row_title = NULL,
                cluster_row_slices = FALSE,
                row_split = split_annot_df,
                top_annotation = topAnn,
                right_annotation = row_ha,
                column_title = set_title
            )
    } else {
        row_ha <- rowAnnotation(DE = split_annot, col = list(DE = colors_for_de))
        hm <- m %>%
            Heatmap(
                name = "Normalized Counts",
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = show_row_names,
                show_column_names = TRUE,
                show_parent_dend_line = FALSE,
                show_row_dend = FALSE,
                column_names_rot = 90,
                row_title = NULL,
                row_title_rot = 0,
                cluster_row_slices = FALSE,
                row_split = split_annot,
                top_annotation = topAnn,
                right_annotation = row_ha,
                column_title = set_title
            )
    }
    p <- draw(hm)
    return(p)
}

# group <- "L1HS"
# counttype <- "telescope_multi"

# groupframe <- resultsdf %>%
#     filter(rte_subfamily == group) %>%
#     filter(counttype == counttype)
# p <- myheatmap(groupframe, facet_var = "genic_loc", filter_var = "rte_length_req", DEvar = "DE", scaled = "notscaled", contrast_samples = contrast_samples, condition_vec = condition_vec)
# mysave("temp1.png", 8, 8)


counttype_label <- ifelse(grepl("multi|Multi", counttype), "Multi", "Unique")
#### GETTING TIDY DATA
map <- setNames(sample_table$condition, sample_table$sample_name)
pvals <- colnames(resultsdf)[str_detect(colnames(resultsdf), "padj_condition")]
l2fc <- colnames(resultsdf)[str_detect(colnames(resultsdf), "log2FoldChange_condition")]
annotations <- c("length", colnames(r_repeatmasker_annotation))
strictly_annotations <- annotations[!(annotations %in% c("gene_id", "family"))]
colsToKeep <- c("gene_id", "family", pvals, l2fc, strictly_annotations)
tidydf <- resultsdf %>%
    filter(counttype == counttype) %>%
    select(all_of(colnames(resultsdf)[(colnames(resultsdf) %in% sample_table$sample_name) | (colnames(resultsdf) %in% colsToKeep)])) %>%
    pivot_longer(cols = -colsToKeep) %>%
    dplyr::rename(sample = name, counts = value) %>%
    mutate(condition = map_chr(sample, ~ map[[.]]))
tidydf$condition <- factor(tidydf$condition, levels = conf$levels)


rte_fams <- tidydf %$% rte_family %>% unique() %>% na.omit()
rte_fams <- rte_fams[rte_fams != "Other"]
for (rte_fam in rte_fams) {

df <- tidydf %>% filter(rte_family == rte_fam)

pf <- df %>%
        group_by(sample, req_integrative, genic_loc, condition) %>%
        summarise(sample_sum = sum(counts), condition = dplyr::first(condition)) %>% ungroup() %>% 
        mutate(req_integrative = ifelse(req_integrative == "Other", "Trucated Ancient", req_integrative))

p <- pf %>% arrange(req_integrative) %>%
    ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
    labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
    mtclosedgridh +
    scale_conditions +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar.png", outputdir, counttype, rte_fam), 15, 12)

p <- pf %>% arrange(req_integrative) %>%
    ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, ref.group = conf$levels[1])
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_stats.png", outputdir, counttype, rte_fam), 15, 12)

    p <- pf %>% arrange(req_integrative) %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = FALSE, ref.group = conf$levels[1])
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_stats_allsigannot.png", outputdir, counttype, rte_fam), 15, 12)

    p <- pf %>% arrange(req_integrative) %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc"), add = c("mean_se"), scales = "free_y") +
    labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
    mtclosedgridh +
    scale_palette +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, ref.group = conf$levels[1])
    mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_reqintegrative_stats.png", outputdir, counttype, rte_fam), 7, 7)

    p <- pf %>% arrange(req_integrative) %>% group_by(sample, condition, genic_loc) %>% summarise(sample_pantype_sum = sum(sample_sum)) %>% ungroup() %>%
        ggbarplot(x = "condition", y = "sample_pantype_sum", facet.by = c("genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = FALSE, ref.group = conf$levels[1]) +
        ggbarplot(data = pf %>% arrange(req_integrative), x = "condition", y = "sample_sum", fill = "req_integrative", facet.by = c("genic_loc"), scales = "free_y") +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
        mtclosedgridh +
        scale_palette +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_reqintegrative_stats_allsigannot.png", outputdir, counttype, rte_fam), 15, 12)


}


rte_subfams <- tidydf %$% rte_subfamily %>% unique() %>% na.omit()
rte_subfams <- rte_subfams[rte_subfams != "Other"]
for (rte_subfam in rte_subfams) {

    df <- tidydf %>% filter(rte_subfamily == rte_subfam)

    pf <- df %>%
            group_by(sample, req_integrative, genic_loc, condition) %>%
            summarise(sample_sum = sum(counts), condition = dplyr::first(condition)) %>% ungroup()  %>%
            mutate(req_integrative = ifelse(req_integrative == "Other", "Trucated Ancient", req_integrative))
    p <- pf %>% arrange(req_integrative) %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_subfam) +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar.png", outputdir, counttype, rte_subfam), 15, 12)

    p <- pf %>% arrange(req_integrative) %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_subfam) +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = TRUE, ref.group = conf$levels[1])
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_stats.png", outputdir, counttype, rte_subfam), 15, 12)

    p <- pf %>% arrange(req_integrative) %>%
        ggbarplot(x = "condition", y = "sample_sum", fill = "condition", facet.by = c("req_integrative", "genic_loc"), add = c("mean_se", "dotplot"), scales = "free_y") +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_subfam) +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_pwc(method = "t_test", label = "p.adj.format", p.adjust.method = "fdr", hide.ns = FALSE, ref.group = conf$levels[1])
        mysaveandstore(sprintf("%s/%s/pan_contrast/%s_bar_stats_allsigannot.png", outputdir, counttype, rte_subfam), 15, 12)

}


#### PLOTTING
groups_that_have_been_run <- c()
groups_not_to_run <- c()
for (ontology in ontologies) {
    ontology_groups <- r_repeatmasker_annotation %>%
        pull(!!sym(ontology)) %>%
        unique()
    ontology_groups <- ontology_groups[ontology_groups != "Other"]
    for (group in ontology_groups) {
        if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
            groups_that_have_been_run <- c(groups_that_have_been_run, group)
            #maybe change to !!group
            groupframe <- tidydf %>% dplyr::filter(!!sym(ontology) == group)
            eligible_modifiers <- c()
            for (modifier in modifiers) {
                values_present <- groupframe %>%
                pull(!!sym(modifier)) %>%
                unique()
            if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
                            }
                eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
                eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
                for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                filter_var <- eligible_modifier_combinations[i, ]$filter_var
                facet_var <- eligible_modifier_combinations[i, ]$facet_var
                plotting_functions <- c("stripp", "pvp", "dep")
                for (function_name in plotting_functions) {
                    tryCatch(
                        {
                            function_current <- get(function_name)
                            plot_width <- 5
                            plot_height <- 5
                            if (facet_var != "ALL") {
                                plot_width <- 8
                            }
                            if (filter_var != "ALL") {
                                plot_title <- groupframe %>%
                                    pull(!!sym(filter_var)) %>%
                                    unique() %>%
                                    grep("FL$|Intact", ., value = TRUE)
                            } else {
                                plot_title <- group
                            }
                            p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var) + ggtitle(plot_title)
                            mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                            if (function_name == "pvp") {
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_labels.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "no", scale_log2 = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_log2.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width +2, plot_height +2)
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "yes", scale_log2 = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_log2labels.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width + 2, plot_height + 2)
                            }
                            if (function_name == "stripp") {
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, stats = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_stats.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height + 2)
                            }
                        },
                        error = function(e) {
                            print(sprintf("Error with  %s %s %s %s %s %s", counttype, contrast, group, function_name, filter_var, facet_var))
                            print(e)
                            tryCatch(
                                {
                                    dev.off()
                                },
                                error = function(e) {
                                    print(e)
                                }
                            )
                        }
                    )
                }
            }
        }
    }
}


for (contrast in contrasts) {
    contrast_of_interest <- contrast
        contrast_level_2 <- contrast_of_interest %>%
        gsub("condition_", "", .) %>%
        gsub("_vs_.*", "", .)
    contrast_level_1 <- contrast_of_interest %>%
        gsub(".*_vs_", "", .)
    contrast_stat <- paste0("stat_", contrast_of_interest)
    contrast_padj <- paste0("padj_", contrast_of_interest)
    contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)
    contrast_samples <- sample_table %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
    groups_that_have_been_run <- c()
    groups_not_to_run <- c()
    for (ontology in ontologies) {
        ontology_groups <- r_repeatmasker_annotation %>%
            pull(!!sym(ontology)) %>%
            unique()
        ontology_groups <- ontology_groups[ontology_groups != "Other"]
        for (group in ontology_groups) {
            if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
                groups_that_have_been_run <- c(groups_that_have_been_run, group)
#maybe change to !!group
                groupframe <- tidydf %>% dplyr::filter(!!sym(ontology) == group)
                eligible_modifiers <- c()
                for (modifier in modifiers) {
                    values_present <- groupframe %>%
                        pull(!!sym(modifier)) %>%
                        unique()
                    if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                        eligible_modifiers <- c(eligible_modifiers, modifier)
                    }
                                    }
                    eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                    eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
                    eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
                
                for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                    filter_var <- eligible_modifier_combinations[i, ]$filter_var
                    facet_var <- eligible_modifier_combinations[i, ]$facet_var
                    plotting_functions <- c("stripp", "pvp", "dep")
                    for (function_name in plotting_functions) {
                        tryCatch(
                            {
                                function_current <- get(function_name)
                                plot_width <- 5
                                plot_height <- 5
                                if (facet_var != "ALL") {
                                    plot_width <- 8
                                }
                                if (filter_var != "ALL") {
                                    plot_title <- groupframe %>%
                                        pull(!!sym(filter_var)) %>%
                                        unique() %>%
                                        grep("FL$|Intact", ., value = TRUE)
                                } else {
                                    plot_title <- group
                                }
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var) + ggtitle(plot_title)
                                mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                if (function_name == "pvp") {
                                    p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "yes") + ggtitle(plot_title)
                                    mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_labels.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                }
                                if (function_name == "stripp") {
                                    p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, stats = "yes") + ggtitle(plot_title)
                                    mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_stats.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height + 2)
                                }
                            },
                            error = function(e) {
                                print(sprintf("Error with  %s %s %s %s %s %s", counttype, contrast, group, function_name, filter_var, facet_var))
                                print(e)
                                tryCatch(
                                    {
                                        dev.off()
                                    },
                                    error = function(e) {
                                        print(e)
                                    }
                                )
                            }
                        )
                    }
                }
            }
        }
    }
}




# Heatmaps


for (contrast in contrasts) {
    contrast_of_interest <- contrast
        contrast_level_2 <- contrast_of_interest %>%
        gsub("condition_", "", .) %>%
        gsub("_vs_.*", "", .)
    contrast_level_1 <- contrast_of_interest %>%
        gsub(".*_vs_", "", .)
    contrast_stat <- paste0("stat_", contrast_of_interest)
    contrast_padj <- paste0("padj_", contrast_of_interest)
    contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)
    contrast_samples <- sample_table %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
    groups_that_have_been_run <- c()
    groups_not_to_run <- c("AluY")
    for (ontology in ontologies) {
        print(ontology)
        ontology_groups <- r_repeatmasker_annotation %>%
            pull(!!sym(ontology)) %>%
            unique()
        ontology_groups <- ontology_groups[ontology_groups != "Other"]
        for (group in ontology_groups) {
            if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
                groups_that_have_been_run <- c(groups_that_have_been_run, group)
groupframe <- resultsdf %>%
                        filter(!!sym(ontology) == group)
                eligible_modifiers <- c()
                for (modifier in modifiers) {
                    values_present <- groupframe %>%
                        pull(!!sym(modifier)) %>%
                        unique()
                    if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                        eligible_modifiers <- c(eligible_modifiers, modifier)
                    }
                    eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                    eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
                    eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
                }
                # first plots without any modifiers
                                plotting_functions <- c("myheatmap")

                for (function_name in plotting_functions) {
                    for (DEvar in c("ALL", "DE")) {
                        for (scaled in c("notscaled", "scaled")) {
                            for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                                filter_var <- eligible_modifier_combinations[i, ]$filter_var
                                facet_var <- eligible_modifier_combinations[i, ]$facet_var

                                plot_width <- 7
                                plot_height <- 7
                                if (DEvar != "ALL") {
                                    plot_width <- 8
                                }
                                tryCatch(
                                    {
                                        function_current <- get(function_name)
                                        p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, DEvar = DEvar, scaled = scaled, contrast_samples = contrast_samples, condition_vec = condition_vec)
                                        mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_%s_%s.png", outputdir, counttype, contrast, function_name, group, filter_var, facet_var, DEvar, scaled), plot_width, plot_height)
                                                                            },
                                    error = function(e) {
                                        print(sprintf("Error with  %s %s %s %s %s %s %s %s %s", counttype, contrast, group, function_name, ontology, filter_var, facet_var, DEvar, scaled))
                                        print(e)
                                        tryCatch(
                                            {
                                                dev.off()
                                            },
                                            error = function(e) {
                                                print(e)
                                            }
                                        )
                                    }
                                )
                                                    }
                    }
                }
            }
        }
    }
}
}




### VENN DIAGRAMS

tryCatch(
    {
                results <- resultsdf %>% filter(counttype == counttype)
            ggvenn <- list()
            for (ontology in ontologies) {
                ontology_groups <- r_repeatmasker_annotation %>%
                    pull(!!sym(ontology)) %>%
                    unique() %>% unlist()
                ontology_groups <- ontology_groups[ontology_groups != "Other"]
                for (group in ontology_groups) {
                    df <- results %>%
                        dplyr::filter((!!sym(ontology)) == group) 

                        vdf <- df %>% mutate(across(starts_with("padj"), ~ ifelse(is.na(.), FALSE, ifelse(. < 0.05, TRUE, FALSE)))) %>%
                                    mutate(req_integrative = ifelse(req_integrative == "Other", "Trucated Ancient", req_integrative))
                        categories <- colnames(vdf) %>% grep("^padj", ., value = TRUE) %>% gsub("padj_condition_", "", .) %>% gsub("_vs.*", "", .)
                        colnames(vdf) <- colnames(vdf) %>% gsub("padj_condition_", "", .) %>% gsub("_vs.*", "", .)
                        up_set <- vdf %>% dplyr::filter(if_all(starts_with("log2F"), ~ .x > 0))
                        down_set <- vdf %>% dplyr::filter(if_all(starts_with("log2F"), ~ .x < 0))
                        tryCatch(
                            {
                            p1 <- up_set %>%
                                    upset(categories, name='condition', width_ratio=0.3,  min_degree=1,
                                    base_annotations=list(
                                    'Intersection size'= intersection_size(counts=FALSE,mapping=aes(fill=req_integrative)) + scale_palette + guides(fill = guide_legend(title="Element Features")))) +
                                    labs(title = sprintf("%s","UP"), caption = counttype) + 
                                    scale_y_continuous(breaks= pretty_breaks())

                            p2 <- down_set %>%
                                    upset(categories, name='condition', width_ratio=0.3,  min_degree=1,
                                    base_annotations=list(
                                    'Intersection size'= intersection_size(counts=FALSE,mapping=aes(fill=req_integrative)) + scale_palette + guides(fill = guide_legend(title="Element Features")))) +
                                    labs(title = sprintf("%s","DOWN"), caption = counttype) + 
                                    scale_y_continuous(breaks= pretty_breaks())
                            ncondtions <- length(categories)
                            mysave(pl = p1, sprintf("%s/%s/pan_contrast/venn/upset_%s_UP.png", outputdir, counttype, group), ncondtions*2, 6)
                            mysave(pl = p2, sprintf("%s/%s/pan_contrast/venn/upset_%s_DOWN.png", outputdir, counttype, group), ncondtions*2, 6)
                            p <- ((p1)/(p2)) + plot_annotation(title = sprintf("%s DE", group))
                            mysave(sprintf("%s/%s/pan_contrast/venn/upset_%s.png", outputdir, counttype, group), ncondtions*2, 10)
                            },
                            error = function(e) {
                                print(e)
                            }
)
                }
            }
        },
    error = function(e) {
        print(e)
        message("Venn diagrams failed - do you only have one contrast?")
        tryCatch(
            {
                dev.off()
            },
            error = function(e) {
                print(e)
            }
        )
    }
)

save(mysaveandstoreplots, file = outputs$plots)




    ######## Genes
    genes <- tidydf %>% filter(gene_or_te == "gene")

    p <- cpmdf %>%
        head(n = 20) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1), y = sen1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene ID", y = "Counts", title = "Top 10 expressed genes") +
        mtopen +
        anchorbar
    path <- paste0(plotdir, "topcounts.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()
    p <- cpmdf %>%
        head(n = 20) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1cpm), y = sen1cpm)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene ID", y = "Counts Per Million (CPM)", title = "Top 10 expressed genes") +
        mtopen +
        anchorbar

    path <- paste0(plotdir, "topcpm.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- genes %>%
        head(n = 20) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1), y = sen1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene ID", y = "Counts", title = "Top 10 expressed genes") +
        mtopen +
        anchorbar
    path <- paste0(plotdir, "topgenecounts.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()
    p <- genes %>%
        head(n = 20) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1cpm), y = sen1cpm)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene ID", y = "Counts Per Million (CPM)", title = "Top 10 expressed genes") +
        mtopen +
        anchorbar

    path <- paste0(plotdir, "topgenecpm.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- genes %>%
        filter(sen1 > 0) %>%
        ggplot() +
        geom_histogram(aes(x = sen1)) +
        xlim(0, 100) +
        labs(x = "Transcript Counts", y = "n", title = "Expressed Gene Count Distribution") +
        mtopen

    path <- paste0(plotdir, "expressedgenecountsdistribution.png")
    png(path, height = 5, width = 5, res = 300, units = "in")
    print(p)
    dev.off()

    p <- genes %>%
        filter(sen1 > 0) %>%
        ggplot() +
        geom_histogram(aes(x = Length)) +
        labs(x = "Read Length (bp)", y = "Count", title = "expressed gene length distribution") +
        xlim(0, 20000) +
        mtopen

    path <- paste0(plotdir, "expressedgenelengthdistribution.png")
    png(path, height = 5, width = 5, res = 300, units = "in")
    print(p)
    dev.off()

    p <- genes %>%
        mutate(detected = ifelse(sen1 > 0, "Expressed", "Not Expressed")) %>%
        group_by(detected) %>%
        summarise(n = n()) %>%
        ggplot() +
        geom_col(aes(x = detected, y = n)) +
        mtopen +
        anchorbar
    path <- paste0(plotdir, "genesidentified.png")
    png(path, height = 5, width = 5, res = 300, units = "in")
    print(p)
    dev.off()


    ################# Repeats
    dfrte <- cpmdf %>% filter(!grepl("gene-", gene_id))

    p <- dfrte %>%
        head(n = 20) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1), y = sen1)) +
        labs(x = "Gene ID", y = "Counts", title = "Top 20 expressed Repeats") +
        mtopen +
        coord_flip() +
        anchorbar
    path <- paste0(plotdir, "toprtecounts.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- dfrte %>%
        filter(grepl("L1HS|L1PA", gene_id)) %>%
        head(n = 20) %>%
        select(gene_id, sen1, Chr) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(gene_id, -sen1), y = sen1)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Gene ID", y = "Counts", title = "Top 20 expressed young L1s") +
        mtopen +
        coord_flip() +
        anchorbar

    path <- paste0(plotdir, "topyoungLINEcounts.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()
    ######################## TODO make plots for top expressed regions
    cpmdf %>%
        filter(grepl("LINE", gene_id)) %>%
        head(n = 20)

    p <- dfrte %>%
        filter(grepl("L1HS", family) | grepl("L1PA", family)) %>%
        group_by(family) %>%
        summarise(cpm = sum(sen1cpm)) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(family, cpm), y = cpm, fill = family)) +
        labs(x = "Family", y = "Counts", title = "") +
        mtopen +
        coord_flip() +
        anchorbar

    path <- paste0(plotdir, "linefamilies.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- dfrte %>%
        filter(grepl("AluY", family)) %>%
        group_by(family) %>%
        summarise(cpm = sum(sen1cpm)) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(family, cpm), y = cpm, fill = family)) +
        labs(x = "Family", y = "Counts", title = "") +
        mtopen +
        coord_flip() +
        anchorbar

    path <- paste0(plotdir, "aluyfamilies.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- dfrte %>%
        filter(grepl("HERVK", family)) %>%
        group_by(family) %>%
        summarise(cpm = sum(sen1cpm)) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(family, cpm), y = cpm, fill = family)) +
        labs(x = "Family", y = "Counts", title = "") +
        mtopen +
        coord_flip() +
        anchorbar

    path <- paste0(plotdir, "hervkfamilies.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()

    p <- dfrte %>%
        filter(grepl("LINE", family) | grepl("SINE", family) | grepl("ERV", family)) %>%
        group_by(rte_superfamily) %>%
        summarise(cpm = sum(sen1cpm)) %>%
        ggplot() +
        geom_col(aes(x = fct_reorder(rte_superfamily, cpm), y = cpm, fill = rte_superfamily)) +
        labs(x = "Family", y = "CPM", title = "") +
        mtopen +
        anchorbar

    path <- paste0(plotdir, "rtefamilies.png")
    png(path, height = 5, width = 10, res = 300, units = "in")
    print(p)
    dev.off()


    # highlight any fl L1HS for scrutiny

    dfrte %>%
        filter(grepl("L1HS", family)) %>%
        filter(sen1cpm > 0) %>%
        arrange(-sen1cpm) %>%
        head(n = 50) %>%
        write_csv(paste0(plotdir, "expressedL1HS.csv"))
    dfrte %>%
        filter(grepl("L1HS", family)) %>%
        filter(sen1cpm > 0) %>%
        arrange(-sen1cpm) %>%
        select(Chr, Start, End) %>%
        head(n = 50) %>%
        write_delim(paste0(plotdir, "expressedL1HS.bed"), col_names = FALSE, delim = "\t")


    dfrte %>%
        filter(grepl("L1HS", family)) %>%
        filter(sen1cpm > 0) %>%
        filter(length > 6000) %>%
        arrange(-sen1cpm) %>%
        head(n = 50) %>%
        write_csv(paste0(plotdir, "expressedL1HS6kb.csv"))
    dfrte %>%
        filter(grepl("L1HS", family)) %>%
        filter(sen1cpm > 0) %>%
        filter(length > 6000) %>%
        arrange(-sen1cpm) %>%
        select(Chr, Start, End) %>%
        head(n = 50) %>%
        write_delim(paste0(plotdir, "expressedL1HS6kb.bed"), col_names = FALSE, delim = "\t")

    # highlight any fl L1PA for scrutiny
    dfrte %>%
        filter(grepl("L1PA", family)) %>%
        filter(sen1cpm > 0) %>%
        arrange(-sen1cpm) %>%
        head(n = 50) %>%
        write_csv(paste0(plotdir, "expressedL1PA.csv"))
    dfrte %>%
        filter(grepl("L1PA", family)) %>%
        filter(sen1cpm > 0) %>%
        arrange(-sen1cpm) %>%
        select(Chr, Start, End) %>%
        head(n = 50) %>%
        write_delim(paste0(plotdir, "expressedL1PA.bed"), col_names = FALSE, delim = "\t")
}

################













################ Alignment analysis

# flag1 <- scanBamFlag(isDuplicate = FALSE, isNotPassingQualityControls = FALSE)
# param1 <- ScanBamParam(flag = flag1, what = "seq")
# aln1 <- readGAlignments("/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen1/alignments/sen1.sorted.bam", use.names = TRUE, param = param1)

aln <- readGAlignments("/users/mkelsey/data/Nanopore/dRNALF1/intermediates/sen1/alignments/genome/sen1.sorted.bam")

# how many reads mapped to RCS?

alndf <- as.data.frame(aln) %>% tibble()

p <- alndf %>%
    group_by(seqnames) %>%
    summarise(n = n(), meanWidth = mean(qwidth)) %>%
    ggplot() +
    geom_col(aes(x = seqnames, y = n)) +
    coord_flip() +
    mtopen +
    anchorbar
png("results/sen1/alignmentsbycontig.png", height = 5, width = 10, res = 300, units = "in")
print(p)
dev.off()

alndf %$% seqnames %>% unique()
library("GenomicFeatures")
alndff <- alndf %>% filter(seqnames != "chromosome:R64-1-1:VIII:450727:453240:1")
p <- alndff %>% ggplot() +
    geom_histogram(aes(x = qwidth)) +
    geom_vline(aes(xintercept = mean(qwidth)), color = "red") +
    labs(x = "Alignment Length (bp)", y = "Count", title = "Alignment Length Distribution") +
    xlim(0, 7500) +
    mtopen +
    anchorbar
png("results/sen1/alignmentsLengthHistogram.png", height = 5, width = 10, res = 300, units = "in")
print(p)
dev.off()
