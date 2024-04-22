source("workflow/scripts/defaults.R")
module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

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
library(plotly)
library(DT)
library(ggExtra)
library(rstatix)
library(purrr)
library(ggpubr)
library(GenomicRanges)
library(paletteer)

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
            "inputdir" = "srna/results/agg/deseq_telescope",
            "outputdir" = "srna/results/agg/repeatanalysis_telescope",
            "tecounttypes" = c("telescope_multi"),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "srna/results/agg/deseq_telescope/resultsdf.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "srna/results/agg/repeatanalysis_telescope/plots.outfile.txt"
        ), env = globalenv())
    }
)

outputdir <- params$outputdir
contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap
peptable <- read.csv(conf$peptable)
tecounttype <- "telescope_multi"


## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
resultsdf <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)

refseq <- import(conf$annotation_genes)
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]
noncoding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NR", mcols(refseq)$transcript_id))]
transcripts <- c(coding_transcripts, noncoding_transcripts)


resultsdf %$% rte_family %>% unique()

# subset <- resultsdf %>% filter(grepl("*tact*", l1_intactness_req))
for (ontology in c("rte_subfamily_limited", "l1_subfamily_limited", "rte_family")) {
    print(ontology)

subset <- resultsdf %>% filter(!!sym(ontology) != "Other") %>% filter(!is.na(!!sym(ontology)))
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
    contrast_samples <- peptable %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
subset %>% pull(!!sym(contrast_log2FoldChange)) %>% length()
resultsdf[match(subjectIDs, resultsdf$gene_id),] %>% pull(!!sym(contrast_log2FoldChange)) %>% length()
te_gene_matrix <- subset %>% dplyr::select(!!sym(ontology), !!sym(contrast_log2FoldChange),loc_integrative,req_integrative ) %>% dplyr::rename(TE = !!sym(contrast_log2FoldChange)) 
te_gene_matrix$GENE <- resultsdf[match(subjectIDs, resultsdf$gene_id),] %>% pull(!!sym(contrast_log2FoldChange))

te_gene_matrix <- te_gene_matrix %>% drop_na()
te_gene_matrix_list[[contrast]] <- te_gene_matrix
cor.test(te_gene_matrix$TE, te_gene_matrix$GENE, method = "spearman", )$estimate
cor.test(te_gene_matrix$TE, te_gene_matrix$GENE, method = "spearman", )$p.value
te_gene_matrix <- te_gene_matrix %>% drop_na()
cor_df <- te_gene_matrix %>% mutate(req_integrative = gsub(".*Intact.*", "Full Length", req_integrative)) %>% group_by(!!sym(ontology), req_integrative, loc_integrative) %>% mutate(groupN = n()) %>% filter(groupN > 4) %>% summarise(cor = cor.test(TE, GENE, method = "spearman")$estimate, pval = cor.test(TE, GENE, method = "spearman")$p.value)
cor_df %$% req_integrative %>% unique()
pf <- cor_df %>% ungroup() %>% filter(loc_integrative != "Centromeric") %>% complete(!!sym(ontology), req_integrative, loc_integrative)
p<-pf %>% mutate(genicfacet = ifelse(loc_integrative == "NonGenic", "", "Genic")) %>% ggplot() + geom_tile(aes(x = !!sym(ontology), y = loc_integrative, fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/rte_genic_cor_%s_%s.png", outputdir, tecounttype, ontology, contrast), 6, 4)

p <- pf %>% mutate(genicfacet = ifelse(loc_integrative == "NonGenic", "", "Genic")) %>% ggplot(aes(x = !!sym(ontology), y = loc_integrative)) + geom_tile(aes(fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", "")), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/rte_genic_cor_pval_%s_%s.png", outputdir, tecounttype, ontology, contrast), 6, 4)
}
te_gene_matrix_all <- Reduce(bind_rows, te_gene_matrix_list)
cor_df <- te_gene_matrix_all %>% mutate(req_integrative = gsub(".*Intact.*", "Full Length", req_integrative)) %>% group_by(!!sym(ontology), req_integrative, loc_integrative) %>% mutate(groupN = n()) %>% filter(groupN > 4) %>% summarise(cor = cor.test(TE, GENE, method = "spearman")$estimate, pval = cor.test(TE, GENE, method = "spearman")$p.value)
cor_df %$% req_integrative %>% unique()
pf <- cor_df %>% ungroup() %>% filter(loc_integrative != "Centromeric") %>% complete(!!sym(ontology), req_integrative, loc_integrative)

p<-pf %>% mutate(genicfacet = ifelse(loc_integrative == "NonGenic", "", "Genic")) %>% ggplot() + geom_tile(aes(x = !!sym(ontology), y = loc_integrative, fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/rte_genic_cor_%s_%s.png", outputdir, tecounttype, ontology, "all"), 6, 4)

p <- pf %>% mutate(genicfacet = ifelse(loc_integrative == "NonGenic", "", "Genic")) %>% ggplot(aes(x = !!sym(ontology), y = loc_integrative)) + geom_tile(aes(fill = cor)) +
    facet_grid(genicfacet ~req_integrative, space = "free", scales = "free") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", "")), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(-0.8,0,0.8), na.value = 'dark grey') +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "")
mysaveandstore(sprintf("%s/%s/rte_genic_cor_pval_%s_%s.png", outputdir, tecounttype, ontology, "all"), 6, 4)
}

### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        big_ontology_groups <- c(big_ontology_groups, resultsdf %>%
            pull(!!sym(ontology)) %>%
            unique())
    }
    big_ontology_groups <- big_ontology_groups %>% unique()

    modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
    region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
    element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
}



#### PLOTTING FUNCTIONS

pvp <- function(df, facet_var = "ALL", filter_var = "ALL", labels = "no") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), ">|Intact|towards"))
    }
    pf <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) 
    top_sig <- pf %>% filter(!!sym(contrast_padj) < 0.05) %>% arrange(padj) %>% head(6) %>% pull(gene_id) 
    if (labels != "no") {
    p <- pf %>%
        {
            ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
                geom_point(aes(color = loc_integrative, size = padj < 0.05, shape = padj < 0.05)) +
                scale_size_manual(values = c(1,2,1)) +
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2), 
                    caption = counttype_label,
                    color = "Genomic Context") +
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
                geom_point(aes(color = loc_integrative, size = padj < 0.05, shape = padj < 0.05)) +
                scale_size_manual(values = c(1,2,1)) +
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2), 
                    caption = counttype_label,
                    color = "Genomic Context") +
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
    return(p)
}
df <- tidydf %>% filter(rte_subfamily == "L1HS")
p <- pvp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "ALL", facet_var = "genic_loc") + ggtitle("L1HS")
mysave("temp1.png", 8, 8)

dep <- function(df, facet_var = "ALL", filter_var = "ALL") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), ">|Intact"))
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
    if (facet_var != "ALL") {
        stats <- "no"
    }
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), ">|Intact"))
    }
    
    stat.test <- df %>%
        group_by(sample) %>%
        summarise(sum(counts), condition = dplyr::first(condition)) %>%
        t_test(`sum(counts)` ~ condition) %>%
        add_significance() %>%
        add_xy_position(x = "condition")
    stat.test$p.format <- p_format(
        stat.test$p,
        accuracy = 0.001,
        leading.zero = TRUE
    )
    if (facet_var != "ALL") {
        pf <- df %>%
            group_by(sample, !!sym(facet_var)) %>%
            summarise(sample_sum = sum(counts), condition = dplyr::first(condition))
        summarydf <- pf %>%
            group_by(condition, !!sym(facet_var)) %>%
            summarise(
                n = n(),
                mean = mean(sample_sum),
                sd = sd(sample_sum)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1))
    } else {
        pf <- df %>%
            group_by(sample) %>%
            summarise(sample_sum = sum(counts), condition = dplyr::first(condition))
        summarydf <- pf %>%
            group_by(condition) %>%
            summarise(
                n = n(),
                mean = mean(sample_sum),
                sd = sd(sample_sum)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1))
    }

    p <- pf %>%
        ggplot() +
        geom_bar(data = summarydf, aes(x = condition, y = mean, fill = condition), stat = "identity") +
        geom_jitter(aes(x = condition, y = sample_sum), width = 0.2, alpha = 0.4) +
        geom_errorbar(data = summarydf, aes(x = condition, ymin = mean - se, ymax = mean + se), width = 0.2) +
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label) +
        extraGGoptions +
        theme(legend.position = "none") +
        mtclosedgridh +
        scale_conditions +
        anchorbar +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = "none")
    if (facet_var != "ALL") {
        p <- p + facet_wrap(facet_var, scales = "free")
    }
    if (stats == "yes") {
        p <- p +
            stat_pvalue_manual(stat.test, label = "p.format", bracket.nudge.y = 0) +
            scale_y_continuous(expand = expansion(mult = c(0.0, 0.1)))
    }
    return(p)
}

# df <- tidydf %>% filter(rte_subfamily == "L1HS")
# p <- stripp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "ALL", facet_var = "genic_loc") + ggtitle("L1HS")
# mysave("temp1.png")


myheatmap <- function(df, facet_var = "ALL", filter_var = "ALL", DEvar = "ALL", scaled = "notscaled", contrast_samples, condition_vec) {
    set_title <- group
    if (filter_var != "ALL") {
        if (str_detect(filter_var, "length_req")) {
            df <- df %>% filter(str_detect(!!sym(filter_var), ">"))
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
        show_row_names <- TRUE
    }
    m <- group_res %>%
        dplyr::select(contrast_samples) %>%
        as.matrix()
    rownames(m) <- group_res %$% gene_id

    if (scaled == "scaled") {
        m <- t(scale(t(m))) %>% na.omit()
        group_res <- group_res %>% filter(gene_id %in% rownames(m))
    }
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
                column_names_rot = 45,
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
# tecounttype <- "telescope_multi"

# groupframe <- resultsdf %>%
#     filter(rte_subfamily == group) %>%
#     filter(tecounttype == tecounttype)
# p <- myheatmap(groupframe, facet_var = "genic_loc", filter_var = "rte_length_req", DEvar = "DE", scaled = "notscaled", contrast_samples = contrast_samples, condition_vec = condition_vec)
# mysave("temp1.png", 8, 8)


for (tecounttype in params$tecounttypes) {
counttype_label <- ifelse(grepl("multi|Multi", tecounttype), "Multi", "Unique")
#### GETTING TIDY DATA
map <- setNames(peptable$condition, peptable$sample_name)
pvals <- colnames(resultsdf)[str_detect(colnames(resultsdf), "padj_condition")]
l2fc <- colnames(resultsdf)[str_detect(colnames(resultsdf), "log2FoldChange_condition")]
annotations <- c("length", colnames(r_repeatmasker_annotation))
strictly_annotations <- annotations[!(annotations %in% c("gene_id", "family"))]
colsToKeep <- c("gene_id", "family", pvals, l2fc, strictly_annotations)
tidydf <- resultsdf %>%
    filter(tecounttype == tecounttype) %>%
    select(all_of(colnames(resultsdf)[(colnames(resultsdf) %in% peptable$sample_name) | (colnames(resultsdf) %in% colsToKeep)])) %>%
    pivot_longer(cols = -colsToKeep) %>%
    rename(sample = name, counts = value) %>%
    mutate(condition = map_chr(sample, ~ map[[.]]))
tidydf$condition <- factor(tidydf$condition, levels = conf$levels)


rte_fams <- tidydf %$% rte_family %>% unique() %>% na.omit()
rte_fams <- rte_fams[rte_fams != "Other"]
for (rte_fam in rte_fams) {

df <- tidydf %>% filter(rte_family == rte_fam) %>% filter(tecounttype == "telescope_multi")

pf <- df %>%
        group_by(sample, req_integrative, genic_loc, condition) %>%
        summarise(sample_sum = sum(counts), condition = dplyr::first(condition)) %>% ungroup() %>% 
        mutate(req_integrative = ifelse(req_integrative == "Other", "Trucated Ancient", req_integrative))
# pf2 <- pf %>% group_by(sample, condition) %>% summarise(sample_sum = sum(sample_sum)) 

p <- pf %>% arrange(req_integrative) %>%
    ggplot() +
    geom_bar(aes(x = sample, y = sample_sum, fill = req_integrative), stat = "identity", position = position_stack(reverse = TRUE)) +
    labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_fam) +
    mtclosedgridh +
    scale_palette +
    anchorbar +
    facet_grid(col = vars(condition), row = vars(genic_loc), scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(sprintf("%s/%s/%s_bar.png", outputdir, tecounttype, rte_fam), 15, 8)
}

rte_subfams <- tidydf %$% rte_subfamily %>% unique() %>% na.omit()
rte_subfams <- rte_subfams[rte_subfams != "Other"]
for (rte_subfam in rte_subfams) {

df <- tidydf %>% filter(rte_subfamily == rte_subfam) %>% filter(tecounttype == "telescope_multi")

pf <- df %>%
        group_by(sample, req_integrative, genic_loc, condition) %>%
        summarise(sample_sum = sum(counts), condition = dplyr::first(condition)) %>% ungroup()  %>%
        mutate(req_integrative = ifelse(req_integrative == "Other", "Trucated Ancient", req_integrative))

# pf2 <- pf %>% group_by(sample, condition) %>% summarise(sample_sum = sum(sample_sum)) 

p <- pf %>% arrange(req_integrative) %>%
    ggplot() +
    geom_bar(aes(x = sample, y = sample_sum, fill = req_integrative), stat = "identity", position = position_stack(reverse = TRUE)) +
    labs(x = "", y = "Sum Normalized Counts", caption = counttype_label, title = rte_subfam) +
    mtclosedgridh +
    scale_palette +
    anchorbar +
    facet_grid(col = vars(condition), row = vars(genic_loc), scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    mysaveandstore(sprintf("%s/%s/%s_bar.png", outputdir, tecounttype, rte_subfam), 15, 8)
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
                values_present <- tidydf %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
            }
            for (modifier in eligible_modifiers) {
                values_present <- resultsdf %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
                eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
                eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
            }
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
                                    grep(">|Intact", ., value = TRUE)
                            } else {
                                plot_title <- group
                            }
                            p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var) + ggtitle(plot_title)
                            mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                            if (function_name == "pvp") {
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_labels.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                            }
                            if (function_name == "stripp") {
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, stats = "yes") + ggtitle(plot_title)
                                mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_stats.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height + 2)
                            }
                        },
                        error = function(e) {
                            print(sprintf("Error with  %s %s %s %s %s %s", tecounttype, contrast, group, function_name, filter_var, facet_var))
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
    contrast_samples <- peptable %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- peptable %>% filter(sample_name %in% contrast_samples) %$% condition
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
                    values_present <- tidydf %>%
                        filter(!!sym(ontology) == group) %>%
                        pull(!!sym(modifier)) %>%
                        unique()
                    if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                        eligible_modifiers <- c(eligible_modifiers, modifier)
                    }
                }
                for (modifier in eligible_modifiers) {
                    values_present <- resultsdf %>%
                        filter(!!sym(ontology) == group) %>%
                        pull(!!sym(modifier)) %>%
                        unique()
                    if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                        eligible_modifiers <- c(eligible_modifiers, modifier)
                    }
                    eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                    eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
                    eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
                }
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
                                        grep(">|Intact", ., value = TRUE)
                                } else {
                                    plot_title <- group
                                }
                                p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var) + ggtitle(plot_title)
                                mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                if (function_name == "pvp") {
                                    p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, labels = "yes") + ggtitle(plot_title)
                                    mysaveandstore(sprintf("% s/%s/%s/%s/%s_%s_%s_labels.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                }
                                if (function_name == "stripp") {
                                    p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, stats = "yes") + ggtitle(plot_title)
                                    mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_stats.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height + 2)
                                }
                            },
                            error = function(e) {
                                print(sprintf("Error with  %s %s %s %s %s %s", tecounttype, contrast, group, function_name, filter_var, facet_var))
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



# Heatmaps


for (tecounttype in params$tecounttypes) {
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
    contrast_samples <- peptable %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- peptable %>% filter(sample_name %in% contrast_samples) %$% condition
    groups_that_have_been_run <- c()
    groups_not_to_run <- c("AluY")
    for (ontology in small_ontologies) {
        print(ontology)
        ontology_groups <- r_repeatmasker_annotation %>%
            pull(!!sym(ontology)) %>%
            unique()
        ontology_groups <- ontology_groups[ontology_groups != "Other"]
        for (group in ontology_groups) {
            if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
                groups_that_have_been_run <- c(groups_that_have_been_run, group)
                eligible_modifiers <- c()
                for (modifier in modifiers) {
                    values_present <- resultsdf %>%
                        filter(!!sym(ontology) == group) %>%
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
                groupframe <- resultsdf %>% filter(!!sym(ontology) == group)
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
                                        mysaveandstore(sprintf("%s/%s/%s/%s/%s_%s_%s_%s_%s.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var, DEvar, scaled), plot_width, plot_height)
                                    },
                                    error = function(e) {
                                        print(sprintf("Error with  %s %s %s %s %s %s %s %s %s", tecounttype, contrast, group, function_name, ontology, filter_var, facet_var, DEvar, scaled))
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
}



### VENN DIAGRAMS

tryCatch(
    {
        for (direction in c("UP", "DOWN")) {
            for (tecounttype in params$tecounttypes) {
                results <- resultsdf %>% filter(tecounttype == tecounttype)
                ggvenn <- list()
                for (ontology in ontologies) {
                    for (group in r_repeatmasker_annotation %>%
                        dplyr::select(!!sym(ontology)) %>%
                        unique() %>%
                        unlist()) {
                        contrastL <- list()
                        for (contrast in contrasts) {
                            contrast_of_interest <- contrast
                            contrast_level_2 <- contrast_of_interest %>%
                                gsub("condition_", "", .) %>%
                                gsub("_vs_.*", "", .)
                            contrast_level_1 <- contrast_of_interest %>%
                                gsub(".*_vs_", "", .)
                            contrast_padj <- paste0("padj_", contrast_of_interest)
                            contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)

                            if (direction == "UP") {
                                sigsearched <- results %>%
                                    dplyr::filter((!!sym(ontology)) == group) %>%
                                    dplyr::filter((!!sym(contrast_padj)) < 0.05) %>%
                                    dplyr::filter((!!sym(contrast_log2FoldChange)) > 0)
                            } else {
                                sigsearched <- results %>%
                                    dplyr::filter((!!sym(ontology)) == group) %>%
                                    dplyr::filter((!!sym(contrast_padj)) < 0.05) %>%
                                    dplyr::filter((!!sym(contrast_log2FoldChange)) < 0)
                            }
                            if (!is.null(filter_var)) {
                                sigsearched <- sigsearched %>%
                                    dplyr::filter(str_detect(!!sym(filter_var), ">"))
                            }
                            sigsearched <- sigsearched %>%
                                dplyr::select(gene_id) %>%
                                na.omit() %$%
                                unlist(gene_id) %>%
                                list()
                            names(sigsearched) <- gsub("condition_", "", contrast)
                            contrastL <- c(contrastL, sigsearched)
                        }
                        p <- ggVennDiagram(contrastL, label_color = "black", set_color = "black") +
                            scale_fill_distiller(palette = "RdBu") +
                            scale_x_continuous(expand = expansion(mult = .2)) +
                            ggtitle(paste(group, direction, sep = " "))
                        if (!is.null(filter_var)) {
                            mysaveandstore(sprintf("%s/%s/%s/ggVenn_%s_%s_%s.png", outputdir, tecounttype, modifier, group, modifier, direction), 6, 6)
                        }
                        mysaveandstore(sprintf("%s/%s/ggVenn_%s_%s.png", outputdir, tecounttype, group, direction), 6, 6)
                    }
                }
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
x <- data.frame()
write.table(x, file = outputs$outfile, col.names = FALSE)
