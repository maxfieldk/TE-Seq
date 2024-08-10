if (interactive()) {
    module_name <<- "srna"
} else {
    module_name <<- snakemake@params$module_name
}
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)
set.seed(123)

library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(stringr)
library(cowplot)
library(pathview)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(ggbeeswarm)
library(paletteer)
library(forcats)
library(ggstance)
library(enrichplot)
library(circlize)
library(ComplexHeatmap)
library(msigdbr)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(dplyr)



# analysis parameters
{
    tryCatch(
        {
            params <- snakemake@params
            inputs <- snakemake@input
            outputs <- snakemake@output
        },
        error = function(e) {
            assign("params", list(
                "counttypes" = conf$counttypes,
                "sample_table" = conf$sample_table,
                "contrasts" = conf$contrasts,
                "genesets_for_heatmaps" = conf$genesets_for_heatmaps,
                "collections_for_gsea" = conf$collections_for_gsea,
                "inputdir" = "srna/results/agg/deseq2/featurecounts_genes",
                "outputdir" = "srna/results/agg/enrichment_analysis"
            ), env = globalenv())
            assign("inputs", list(
                resultsdf = paste0("srna/results/agg/deseq/resultsdf.tsv")
            ), env = globalenv())
            assign("outputs", list(
                plots = "srna/results/agg/enrichment_analysis/enrichment_analysis_plots.RData",
                results_table_targetted = "srna/results/agg/enrichment_analysis/results_table_targetted.tsv",
                results_table_unbiased = "srna/results/agg/enrichment_analysis/results_table_unbiased.tsv"
            ), env = globalenv())
        }
    )
}


# load results
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
resultsdf1 <- resultsdf1[resultsdf1$gene_id != "__no_feature", ]
res <- resultsdf1 %>% filter(counttype == counttype[1])
res <- res %>% filter(gene_or_te == "gene")

if (conf$gtf_id_mapping$response == "yes") {
    mapping <- read_delim(conf$gtf_id_mapping$path)
    res <- res %>%
        left_join(mapping) %>%
        mutate(gene_id = gene_name)
}

#### GETTING TIDY DATA
map <- sample_table %>%
    dplyr::select(sample_name, condition) %>%
    dplyr::rename(sample = sample_name)
pvals <- colnames(res)[str_detect(colnames(res), "padj_condition")]
l2fc <- colnames(res)[str_detect(colnames(res), "log2FoldChange_condition")]
colsToKeep <- c("gene_id", pvals, l2fc)

tidydf <- res %>%
    dplyr::select(all_of(colnames(res)[(colnames(res) %in% sample_table$sample_name) | (colnames(res) %in% colsToKeep)])) %>%
    pivot_longer(cols = -colsToKeep) %>%
    dplyr::rename(sample = name, counts = value) %>%
    left_join(map)
tidydf$condition <- factor(tidydf$condition, levels = conf$levels)

contrast_label_map <- tibble(contrast = params[["contrasts"]], label = gsub("_", " ", gsub("condition_", "", params[["contrasts"]])))

pvp <- function(df, facet_var = "ALL", filter_var = "ALL", labels = "no", scale_log2 = "no") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), "FL$|^Intact|towards"))
    }
    if (scale_log2 == "no") {
        pf <- df %>%
            group_by(across(all_of(c(colsToKeep, "condition")))) %>%
            summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj)), l2fc = dplyr::first(!!sym(contrast_l2fc))) %>%
            pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
            mutate(color_direction = ifelse(is.na(padj), "NS", ifelse(padj > 0.05, "NS", ifelse(l2fc > 0, "UP", "DOWN"))))
    } else {
        pf <- df %>%
            group_by(across(all_of(c(colsToKeep, "condition")))) %>%
            summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj)), l2fc = dplyr::first(!!sym(contrast_l2fc))) %>%
            pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
            mutate(across(conf$levels, ~ log2(.x + 1))) %>%
            mutate(color_direction = ifelse(is.na(padj), "NS", ifelse(padj > 0.05, "NS", ifelse(l2fc > 0, "UP", "DOWN"))))
    }

    top_sig <- pf %>%
        filter(!!sym(contrast_padj) < 0.05) %>%
        arrange(padj) %>%
        head(6) %>%
        pull(gene_id)
    if (labels != "no") {
        p <- pf %>%
            {
                ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
                    geom_point(aes(color = color_direction)) +
                    labs(
                        x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2),
                        title = set_title,
                        color = "DE",
                        shape = "padj < 0.05"
                    ) +
                    mtclosed +
                    scale_directions +
                    theme(aspect.ratio = 1) +
                    coord_cartesian(clip = "off") +
                    geom_text_repel(data = . %>% mutate(label = ifelse(gene_id %in% top_sig, gene_id, "")), aes(label = label), max.overlaps = Inf, show.legend = FALSE) +
                    coord_fixed(
                        xlim = range(c(.[[contrast_level_1]], .[[contrast_level_2]])),
                        ylim = range(c(.[contrast_level_1], .[contrast_level_2]))
                    ) +
                    guides(fill = guide_legend(title = "Genomic Context"))
            }
    } else {
        p <- pf %>%
            {
                ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
                    geom_point(aes(color = color_direction)) +
                    labs(
                        x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2),
                        title = set_title,
                        color = "DE",
                        shape = "padj < 0.05"
                    ) +
                    mtclosed +
                    scale_directions +
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

# GSEA targeted
genecollections <- names(params[["collections_for_gsea"]])
rm(gse_df)
for (contrast in params[["contrasts"]]) {
    # PREP RESULTS FOR GSEA
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
    res <- res %>% arrange(-!!sym(contrast_stat))

    ordered_by_stat <- setNames(res[[contrast_stat]], res$gene_id) %>% na.omit()
    resl2fc <- res %>% arrange(-!!sym(contrast_l2fc))
    ordered_by_l2fc <- setNames(resl2fc[[contrast_l2fc]], resl2fc$gene_id) %>% na.omit()
    for (collection in genecollections) {
        tryCatch(
            {
                genesets <- read.gmt(params[["collections_for_gsea"]][[collection]])
                gse <- GSEA(ordered_by_stat, TERM2GENE = genesets, maxGSSize = 100000, minGSSize = 1)
                df <- gse@result %>% tibble()
                df$collection <- collection
                df$contrast <- contrast
                # gse_results[[contrast]][[collection]] <- as.data.frame(df) %>% tibble()
                if (!exists("gse_df")) {
                    gse_df <<- df
                } else {
                    gse_df <<- rbind(gse_df, df)
                }
                genesettheme <- theme_gray() + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

                tryCatch(
                    {
                        for (num in c(5, 10, 15)) {
                            dftemp <- arrange(df, -abs(NES)) %>%
                                group_by(sign(NES)) %>%
                                slice(1:num) %>%
                                mutate(log10padj = -log10(p.adjust))
                            dftemp <- dftemp %>% mutate(Description = str_wrap(as.character(Description) %>% gsub("_", " ", .), width = 40))
                            dftemp <- dftemp %>% mutate(`Gene Ratio` = 0.01 * as.numeric(gsub("%", "", gsub(",.*", "", gsub("tags=", "", leading_edge)))))
                            p <- ggplot(dftemp, aes(NES, fct_reorder(Description, NES), fill = log10padj)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient(high = "red", low = "white", limits = c(0, 5), oob = scales::oob_squish) +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)

                            p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust) * `sign(NES)`)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                guides(fill = guide_legend(title = "Signed \n-log10(p.adjust)")) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )

                tryCatch(
                    {
                        for (genesetid in conf$genesets_for_gseaplot) {
                            p <- gseaplot2(gse, geneSetID = genesetid, pvalue_table = TRUE, subplots = 1:2, ES_geom = "line")
                        }
                        mysaveandstore(sprintf("%s/%s/gsea/%s_gsea.pdf", params[["outputdir"]], contrast, genesetid), w = 7, h = 3, res = 300)
                    },
                    error = function(e) {
                        print("")
                    }
                )
            },
            error = function(e) {
                print("probably no enrichments in this collection")
            }
        )
    }
}

tryCatch(
    {
        gres <- gse_df %>% tibble()
        gres %>% write_delim(outputs$results_table_targetted, delim = "\t")
        for (collec in genecollections) {
            grestemp <- gres %>%
                filter(collection == collec) %>%
                left_join(contrast_label_map) %>%
                filter(grepl(paste0(conf$levels[1], "$"), label))
            sigIDs <- grestemp %>%
                mutate(direction = ifelse(NES > 0, "UP", "DOWN")) %>%
                group_by(contrast, direction) %>%
                arrange(p.adjust) %>%
                slice_head(n = 5) %$% ID %>%
                unique()
            p <- grestemp %>%
                dplyr::filter(ID %in% sigIDs) %>%
                mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
                mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(label = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
                mutate(label = factor(label, levels = conf$levels)) %>%
                mutate(ID = fct_reorder(ID, NES)) %>%
                ggplot(aes(x = label, y = ID)) +
                geom_tile(aes(fill = NES), color = "black") +
                theme(legend.position = "none") +
                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                mtclosed +
                theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
                labs(x = "", y = "", title = collec) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                coord_equal()
            mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/targetted/gsea_top_%s_grid_1.pdf", collec), 6, 0.5 * length(sigIDs))
            mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/targetted/gsea_top_%s_grid_2.pdf", collec), 6, 0.75 * length(sigIDs))
            mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/targetted/gsea_top_%s_grid_3.pdf", collec), 6, 1 * length(sigIDs))
        }



        # plot core enrichments

        n_top_sets <- 5
        for (contrast in params[["contrasts"]]) {
            contrast_of_interest <- contrast
            contrast_level_2 <- contrast_of_interest %>%
                gsub("condition_", "", .) %>%
                gsub("_vs_.*", "", .)
            contrast_level_1 <- contrast_of_interest %>%
                gsub(".*_vs_", "", .)
            contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
            contrast_stat <- paste0("stat_", contrast)
            contrast_l2fc <- paste0("log2FoldChange_", contrast)
            contrast_padj <- paste0("padj_", contrast)
            contrast_samples <- sample_table %>%
                filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
                pull(sample_name)
            condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
            res <- res %>% arrange(-!!sym(contrast_stat))

            for (collec in genecollections) {
                gse <- gse_df %>% filter(collection == collec)
                df <- arrange(gse, -abs(NES)) %>%
                    group_by(sign(NES)) %>%
                    slice(1:n_top_sets)
                ids <- df$ID
                core_enrichments <- df$core_enrichment
                core_enrichments <- str_split(core_enrichments, "/")
                names(core_enrichments) <- ids
                for (geneset in c(names(core_enrichments), conf$genesets_for_gseaplot)) {
                    all_genes_in_set <- read.gmt(params[["collections_for_gsea"]][[collec]]) %>% filter(term == geneset) %$% gene
                    if (length(all_genes_in_set) == 0) {
                        next
                    }
                    params[["collections_for_gsea"]][[collection]]
                    genestoplot <- core_enrichments[[geneset]]
                    set_title <- geneset

                    contrastconditions <- gsub("condition_", "", contrast) %>%
                        str_split("_vs_") %>%
                        pluck(1)
                    sample_vec <- sample_table %>%
                        filter(condition %in% contrastconditions) %>%
                        pull(sample_name)
                    condition_vec <- sample_table %>%
                        filter(sample_name %in% sample_vec) %>%
                        pull(condition)

                    ngenes <- length(genestoplot)
                    ncol <- 4
                    barplotpf <- res %>%
                        filter(gene_id %in% genestoplot) %>%
                        pivot_longer(cols = conf$samples, names_to = "sample_name", values_to = "counts") %>%
                        left_join(sample_table) %>%
                        mutate(sample_name = factor(sample_name, levels = conf$samples)) %>%
                        mutate(condition = factor(condition, levels = conf$levels))
                    p <- barplotpf %>% ggbarplot(
                        x = "condition", y = "counts", fill = "condition", add = c("mean_se", "dotplot"),
                        facet.by = "gene_id", scales = "free", ncol = 4
                    ) +
                        labs(title = set_title, x = element_blank()) +
                        mtclosedgridh + scale_conditions + anchorbar +
                        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                    tryCatch(
                        {
                            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/barplot_%s.pdf", params[["outputdir"]], contrast, collec, set_title), 8, 2 * round(length(genestoplot) / ncol))

                            p <- pvp(tidydf %>% filter(gene_id %in% all_genes_in_set), labels = "no", scale_log2 = "yes")
                            mysaveandstore(sprintf("%s/%s/gsea/%s/all_genes/scatter_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = 5, h = 5, res = 300)

                            p <- pvp(tidydf %>% filter(gene_id %in% all_genes_in_set), labels = "yes", scale_log2 = "yes")
                            mysaveandstore(sprintf("%s/%s/gsea/%s/all_genes/scatter_labs_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = 5, h = 5, res = 300)

                            heatmapprep <- res %>% filter(gene_id %in% all_genes_in_set)
                            m <- as.matrix(heatmapprep %>% dplyr::select(sample_vec))
                            rownames(m) <- heatmapprep %$% gene_id
                            scaledm <- t(scale(t(m))) %>% na.omit()

                            pvals <- heatmapprep %>%
                                filter(gene_id %in% rownames(scaledm)) %>%
                                pull(!!sym(contrast_padj))
                            pch <- ifelse(pvals < 0.05, "*", "ns")
                            row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), col = list(pvalue = c("ns" = "grey", "*" = "red")))

                            pvalue_adj <- heatmapprep %>%
                                filter(gene_id %in% rownames(scaledm)) %>%
                                pull(!!sym(contrast_padj))
                            is_sig <- pvalue_adj < 0.05
                            pch <- rep("*", length(pvalue_adj))
                            pch[!is_sig] <- NA
                            pvalue_adj_col_fun <- colorRamp2(c(0, 2, 3), c("white", "blue", "red"))
                            ha <- rowAnnotation(pvalue_adj = anno_simple(-log10(pvalue_adj), col = pvalue_adj_col_fun, pch = pch))
                            conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
                            topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))

                            hm <- scaledm %>%
                                Heatmap(
                                    name = "Scaled Normalized Counts",
                                    cluster_rows = FALSE,
                                    cluster_columns = FALSE,
                                    show_row_names = TRUE,
                                    show_column_names = TRUE,
                                    column_names_rot = 90,
                                    split = pch,
                                    top_annotation = topAnn,
                                    right_annotation = ha,
                                    row_title = set_title,
                                    # border_gp = gpar(col = "black")
                                )
                            lgd_pvalue_adj <- Legend(
                                title = "p-value", col_fun = pvalue_adj_col_fun, at = c(0, 1, 2, 3),
                                labels = c("1", "0.1", "0.01", "0.001")
                            )
                            # and one for the significant p-values
                            lgd_sig <- Legend(pch = "*", type = "points", labels = "< 0.05")
                            # these two self-defined legends are added to the plot by `annotation_legend_list`
                            p <- wrap_elements(grid.grabExpr(draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")))
                            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/heatmap_%s.pdf", params[["outputdir"]], contrast, collection, set_title), w = min(0.75 * dim(scaledm)[2], 12), h = min(dim(scaledm)[1], 12), res = 300)

                            # heatmap all samples
                            m <- as.matrix(heatmapprep %>% dplyr::select(sample_table$sample_name))
                            rownames(m) <- heatmapprep %$% gene_id
                            scaledm <- t(scale(t(m))) %>% na.omit()

                            conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
                            topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))

                            hm <- scaledm %>%
                                Heatmap(
                                    name = "Scaled Normalized Counts",
                                    cluster_rows = TRUE,
                                    cluster_columns = FALSE,
                                    show_row_names = TRUE,
                                    show_column_names = TRUE,
                                    column_names_rot = 90,
                                    top_annotation = topAnn,
                                    row_title = set_title,
                                    border_gp = gpar(col = "black")
                                )

                            for (contrast_padj in grep("padj", colnames(heatmapprep), value = TRUE)) {
                                pvals <- heatmapprep %>%
                                    filter(gene_id %in% rownames(scaledm)) %>%
                                    pull(!!sym(contrast_padj))
                                pch <- ifelse(pvals < 0.05, "*", "ns")
                                row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), col = list(pvalue = c("ns" = "grey", "*" = "red")))

                                pvalue_adj <- heatmapprep %>%
                                    filter(gene_id %in% rownames(scaledm)) %>%
                                    pull(!!sym(contrast_padj))
                                is_sig <- pvalue_adj < 0.05
                                pch <- rep("*", length(pvalue_adj))
                                pch[!is_sig] <- NA
                                pvalue_adj_col_fun <- colorRamp2(c(0, 2, 3), c("tan", "tan", "tan"))
                                ha <- rowAnnotation(contrast_padj = anno_simple(-log10(pvalue_adj), col = pvalue_adj_col_fun, pch = pch))
                                hm <- hm + ha
                            }

                            lgd_pvalue_adj <- Legend(
                                title = "p-value", col_fun = pvalue_adj_col_fun, at = c(0, 1, 2, 3),
                                labels = c("1", "0.1", "0.01", "0.001")
                            )
                            # and one for the significant p-values
                            lgd_sig <- Legend(pch = "*", type = "points", labels = "< 0.05")
                            # these two self-defined legends are added to the plot by `annotation_legend_list`
                            p <- wrap_elements(grid.grabExpr(draw(hm, auto_adjust = FALSE, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "right", annotation_legend_side = "right")))
                            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/heatmap_%s.pdf", params[["outputdir"]], contrast, collection, set_title), w = min(0.75 * dim(scaledm)[2], 12), h = min(dim(scaledm)[1] / 7.5, 20), res = 300)
                        },
                        error = function(e) {
                            print(e)
                        }
                    )
                }
            }
        }
    },
    error = function(e) {
        print("probably no enrichments in custom sets")
    }
)



for (geneset in names(params[["genesets_for_heatmaps"]])) {
    # geneset <- geneset
    genestoplot <- read_csv(params[["genesets_for_heatmaps"]][[geneset]], col_names = FALSE) %>% pull(X1)
    set_title <- geneset


    # first for all conditions
    # make data long using samples
    ngenes <- length(genestoplot)
    ncol <- 4
    barplotpf <- res %>%
        filter(gene_id %in% genestoplot) %>%
        pivot_longer(cols = conf$samples, names_to = "sample_name", values_to = "counts") %>%
        left_join(sample_table) %>%
        mutate(sample_name = factor(sample_name, levels = conf$samples)) %>%
        mutate(condition = factor(condition, levels = conf$levels))



    p <- barplotpf %>% ggbarplot(
        x = "condition", y = "counts", fill = "condition", add = c("mean_se", "dotplot"),
        facet.by = "gene_id", scales = "free", ncol = 4
    ) +
        labs(title = set_title, x = element_blank()) +
        mtclosedgridh + scale_conditions + anchorbar +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    mysaveandstore(sprintf("%s/barplot_%s.pdf", params[["outputdir"]], set_title), 8, 2 * round(ngenes / ncol))

    heatmapprep <- res %>% filter(gene_id %in% genestoplot)
    m <- as.matrix(heatmapprep %>% dplyr::select(conf$samples))
    rownames(m) <- heatmapprep %$% gene_id
    scaledm <- t(scale(t(m))) %>% na.omit()

    conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
    topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))

    hm <- scaledm %>%
        Heatmap(
            name = "Normalized Count Z-score",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_names = TRUE,
            show_column_names = TRUE,
            column_names_rot = 90,
            top_annotation = topAnn,
            row_title = set_title,
            border_gp = gpar(col = "black")
        )
    p <- wrap_elements(grid.grabExpr(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")))
    hmscalew <- 0.7
    hmscaleh <- 0.7
    mysaveandstore(sprintf("%s/heatmap_%s.pdf", params[["outputdir"]], set_title), w = min(dim(scaledm)[2] * hmscalew, 12), h = min(dim(scaledm)[1] * hmscaleh, 12), res = 300)

    # now for each contrast
    for (contrast in params[["contrasts"]]) {
        contrast_of_interest <- contrast
        contrast_level_2 <- contrast_of_interest %>%
            gsub("condition_", "", .) %>%
            gsub("_vs_.*", "", .)
        contrast_level_1 <- contrast_of_interest %>%
            gsub(".*_vs_", "", .)
        contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
        contrast_stat <- paste0("stat_", contrast)
        contrast_l2fc <- paste0("log2FoldChange_", contrast)
        contrast_padj <- paste0("padj_", contrast)
        contrast_samples <- sample_table %>%
            filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
            pull(sample_name)
        condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
        res <- res %>% arrange(-!!sym(contrast_stat))

        contrastconditions <- gsub("condition_", "", contrast) %>%
            str_split("_vs_") %>%
            pluck(1)
        sample_vec <- sample_table %>%
            filter(condition %in% contrastconditions) %>%
            pull(sample_name)
        condition_vec <- sample_table %>%
            filter(sample_name %in% sample_vec) %>%
            pull(condition)

        heatmapprep <- res %>% filter(gene_id %in% genestoplot)
        m <- as.matrix(heatmapprep %>% dplyr::select(sample_vec))
        rownames(m) <- heatmapprep %$% gene_id
        scaledm <- t(scale(t(m))) %>% na.omit()

        pvals <- heatmapprep %>%
            filter(gene_id %in% rownames(scaledm)) %>%
            pull(!!sym(contrast_padj))
        pch <- ifelse(pvals < 0.05, "*", "ns")
        row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), col = list(pvalue = c("ns" = "grey", "*" = "red")))

        pvalue_adj <- heatmapprep %>%
            filter(gene_id %in% rownames(scaledm)) %>%
            pull(!!sym(contrast_padj))
        is_sig <- pvalue_adj < 0.05
        pch <- rep("*", length(pvalue_adj))
        pch[!is_sig] <- NA
        pvalue_adj_col_fun <- colorRamp2(c(0, 2, 3), c("white", "blue", "red"))
        ha <- rowAnnotation(pvalue_adj = anno_simple(-log10(pvalue_adj), col = pvalue_adj_col_fun, pch = pch))
        conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
        topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))
        hm <- scaledm %>%
            Heatmap(
                name = "Normalized Counts Z-score",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_names_rot = 90,
                split = pch,
                top_annotation = topAnn,
                right_annotation = ha,
                row_title = set_title,
                border_gp = gpar(col = "black")
            )
        lgd_pvalue_adj <- Legend(
            title = "p-value", col_fun = pvalue_adj_col_fun, at = c(0, 1, 2, 3),
            labels = c("1", "0.1", "0.01", "0.001")
        )
        # and one for the significant p-values
        lgd_sig <- Legend(pch = "*", type = "points", labels = "< 0.05")
        # these two self-defined legends are added to the plot by `annotation_legend_list`
        p <- wrap_elements(grid.grabExpr(draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")))
        mysaveandstore(sprintf("%s/%s/heatmap_%s.pdf", params[["outputdir"]], contrast, set_title), w = min(dim(scaledm)[2] * hmscalew, 12), h = min(dim(scaledm)[1] * hmscaleh, 12), res = 300)
    }
}


# GSEA untargeted
tryCatch(
    {
        gene_sets <- msigdbr(species = confALL$aref$species)
    },
    error = function(e) {
        gene_sets <<- msigdbr(species = "human")
    }
)
rm(gse_df)
for (contrast in params[["contrasts"]]) {
    # PREP RESULTS FOR GSEA
    contrast_of_interest <- contrast
    contrast_level_2 <- contrast_of_interest %>%
        gsub("condition_", "", .) %>%
        gsub("_vs_.*", "", .)
    contrast_level_1 <- contrast_of_interest %>%
        gsub(".*_vs_", "", .)
    contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
    contrast_padj <- paste0("padj_", contrast)
    contrast_samples <- sample_table %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
    res <- res %>% arrange(-!!sym(contrast_stat))

    ordered_by_stat <- setNames(res[[contrast_stat]], res$gene_id) %>% na.omit()
    resl2fc <- res %>% arrange(-!!sym(contrast_l2fc))
    ordered_by_l2fc <- setNames(resl2fc[[contrast_l2fc]], resl2fc$gene_id) %>% na.omit()
    for (category in gene_sets %$% gs_cat %>% unique()) {
        cat(category, "\n")
        tryCatch(
            {
                collection <- category
                msigdbr_df <- gene_sets %>% filter(gs_cat == category)
                msigdbr_t2g <- msigdbr_df %>%
                    dplyr::distinct(gs_name, gene_symbol) %>%
                    as.data.frame()
                gse <- GSEA(ordered_by_stat, TERM2GENE = msigdbr_t2g, maxGSSize = 100000, minGSSize = 1)
                df <- gse@result %>% tibble()
                df$collection <- collection
                df$contrast <- contrast
                # gse_results[[contrast]][[collection]] <- as.data.frame(df) %>% tibble()
                if (!exists("gse_df")) {
                    gse_df <<- df
                } else {
                    gse_df <<- rbind(gse_df, df)
                }
                genesettheme <- theme_gray() + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

                tryCatch(
                    {
                        for (num in c(5, 10, 15)) {
                            dftemp <- arrange(df, -abs(NES)) %>%
                                group_by(sign(NES)) %>%
                                slice(1:num) %>%
                                mutate(log10padj = -log10(p.adjust))
                            dftemp <- dftemp %>% mutate(Description = str_wrap(as.character(Description) %>% gsub("_", " ", .), width = 40))
                            dftemp <- dftemp %>% mutate(`Gene Ratio` = 0.01 * as.numeric(gsub("%", "", gsub(",.*", "", gsub("tags=", "", leading_edge)))))
                            p <- ggplot(dftemp, aes(NES, fct_reorder(Description, NES), fill = log10padj)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient(high = "red", low = "white", limits = c(0, 5), oob = scales::squish) +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)

                            p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust) * `sign(NES)`)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL) +
                                guides(fill = guide_legend(title = "Signed \n-log10(p.adjust)")) +
                                labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )

                tryCatch(
                    {
                        for (genesetid in conf$genesets_for_gseaplot) {
                            p <- gseaplot2(gse, geneSetID = genesetid, pvalue_table = TRUE, subplots = 1:2, ES_geom = "line")
                        }
                        mysaveandstore(sprintf("%s/%s/gsea/%s_gsea.pdf", params[["outputdir"]], contrast, genesetid), w = 8, h = 3, res = 300)
                    },
                    error = function(e) {
                        print("")
                    }
                )
            },
            error = function(e) {
                print("probably no enrichments in this collection")
            }
        )
    }
}

n_top_sets <- 5
for (contrast in params[["contrasts"]]) {
    contrast_of_interest <- contrast
    contrast_level_2 <- contrast_of_interest %>%
        gsub("condition_", "", .) %>%
        gsub("_vs_.*", "", .)
    contrast_level_1 <- contrast_of_interest %>%
        gsub(".*_vs_", "", .)
    contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
    contrast_padj <- paste0("padj_", contrast)
    contrast_samples <- sample_table %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- sample_table %>% filter(sample_name %in% contrast_samples) %$% condition
    contrast_significance <- paste0("Significance_", contrast)

    res <- res %>% arrange(-!!sym(contrast_stat))
    for (collec in gse_df %$% collection %>% unique()) {
        gse <- gse_df %>% filter(collection == collec)
        df <- arrange(gse, -abs(NES)) %>%
            group_by(sign(NES)) %>%
            slice(1:n_top_sets)
        ids <- df$ID
        for (geneset in c(ids, conf$genesets_for_gseaplot)) {
            tryCatch(
                {
                    genestoplot <- gene_sets %>%
                        filter(gs_name == geneset) %>%
                        pull(gene_symbol) %>%
                        unique()
                    set_title <- geneset

                    contrastconditions <- gsub("condition_", "", contrast) %>%
                        str_split("_vs_") %>%
                        pluck(1)
                    sample_vec <- sample_table %>%
                        filter(condition %in% contrastconditions) %>%
                        pull(sample_name)
                    condition_vec <- sample_table %>%
                        filter(sample_name %in% sample_vec) %>%
                        pull(condition)

                    p <- pvp(tidydf %>% filter(gene_id %in% genestoplot), labels = "no", scale_log2 = "yes")
                    mysaveandstore(sprintf("%s/%s/gsea/%s/all_genes/scatter_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = 5, h = 5, res = 300)

                    p <- pvp(tidydf %>% filter(gene_id %in% genestoplot), labels = "yes", scale_log2 = "yes")
                    mysaveandstore(sprintf("%s/%s/gsea/%s/all_genes/scatter_labs_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = 5, h = 5, res = 300)

                    heatmapprep <- res %>% filter(gene_id %in% genestoplot)
                    m <- as.matrix(heatmapprep %>% dplyr::select(sample_vec))
                    rownames(m) <- heatmapprep %$% gene_id
                    scaledm <- t(scale(t(m))) %>% na.omit()

                    pvals <- heatmapprep %>%
                        filter(gene_id %in% rownames(scaledm)) %>%
                        pull(!!sym(contrast_padj))
                    pch <- ifelse(pvals < 0.05, "*", "ns")
                    row_ha <- rowAnnotation(pvalue = anno_simple(pch, pch = pch), col = list(pvalue = c("ns" = "grey", "*" = "red")))

                    pvalue_adj <- heatmapprep %>%
                        filter(gene_id %in% rownames(scaledm)) %>%
                        pull(!!sym(contrast_padj))
                    is_sig <- pvalue_adj < 0.05
                    pch <- rep("*", length(pvalue_adj))
                    pch[!is_sig] <- NA
                    pvalue_adj_col_fun <- colorRamp2(c(0, 2, 3), c("white", "blue", "red"))
                    ha <- rowAnnotation(pvalue_adj = anno_simple(-log10(pvalue_adj), col = pvalue_adj_col_fun, pch = pch))
                    conditions <- sample_table[match(colnames(m), sample_table$sample_name), ]$condition
                    topAnn <- ComplexHeatmap::HeatmapAnnotation(Condition = conditions, col = list(Condition = condition_palette))

                    hm <- scaledm %>%
                        Heatmap(
                            name = "Scaled Normalized Counts",
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            column_names_rot = 90,
                            split = pch,
                            top_annotation = topAnn,
                            right_annotation = ha,
                            row_title = set_title,
                            # border_gp = gpar(col = "black")
                        )
                    lgd_pvalue_adj <- Legend(
                        title = "p-value", col_fun = pvalue_adj_col_fun, at = c(0, 1, 2, 3),
                        labels = c("1", "0.1", "0.01", "0.001")
                    )
                    # and one for the significant p-values
                    lgd_sig <- Legend(pch = "*", type = "points", labels = "< 0.05")
                    # these two self-defined legends are added to the plot by `annotation_legend_list`
                    p <- wrap_elements(grid.grabExpr(draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")))
                    mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/heatmap_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = min(dim(scaledm)[2], 12), h = min(dim(scaledm)[1], 12), res = 300)
                },
                error = function(e) {
                    print("")
                }
            )
        }
    }
}

gres <- gse_df %>% tibble()
gres %>% write_delim(outputs$results_table_unbiased, delim = "\t")
for (collec in gse_df %$% collection %>% unique()) {
    grestemp <- gres %>%
        filter(collection == collec) %>%
        filter(grepl(paste0(conf$levels[1], "$"), contrast)) %>%
        left_join(contrast_label_map)
    sigIDs <- grestemp %>%
        mutate(direction = ifelse(NES > 0, "UP", "DOWN")) %>%
        group_by(contrast, direction) %>%
        arrange(p.adjust) %>%
        slice_head(n = 5) %$% ID %>%
        unique()
    p <- grestemp %>%
        dplyr::filter(ID %in% sigIDs) %>%
        mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
        mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
        mutate(label = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
        mutate(label = factor(label, levels = conf$levels)) %>%
        mutate(ID = fct_reorder(ID, NES)) %>%
        ggplot(aes(x = label, y = ID)) +
        geom_tile(aes(fill = NES), color = "black") +
        theme(legend.position = "none") +
        scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
        mtclosed +
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
        labs(x = "", y = "", title = collec) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_equal()
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_1.pdf", collec), 6, 0.5 * length(sigIDs))
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_2.pdf", collec), 6, 0.75 * length(sigIDs))
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_3.pdf", collec), 6, 1 * length(sigIDs))

    grestemp <- gres %>%
        filter(collection == collec) %>%
        filter(grepl(paste0(conf$levels[1], "$"), contrast)) %>%
        left_join(contrast_label_map)
    sigIDs <- grestemp %>%
        mutate(direction = ifelse(NES > 0, "UP", "DOWN")) %>%
        group_by(contrast, direction) %>%
        arrange(p.adjust) %>%
        slice_head(n = 10) %$% ID %>%
        unique()
    p <- grestemp %>%
        dplyr::filter(ID %in% sigIDs) %>%
        mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
        mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
        mutate(label = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
        mutate(label = factor(label, levels = conf$levels)) %>%
        mutate(ID = fct_reorder(ID, NES)) %>%
        ggplot(aes(x = label, y = ID)) +
        geom_tile(aes(fill = NES), color = "black") +
        theme(legend.position = "none") +
        scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
        mtclosed +
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
        labs(x = "", y = "", title = collec) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_equal()
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_1_2.pdf", collec), 6, 0.5 * length(sigIDs))
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_2_2.pdf", collec), 6, 0.75 * length(sigIDs))
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/unbiased/gsea_top_%s_grid_3_2.pdf", collec), 6, 1 * length(sigIDs))
}

if (conf$store_env_as_rds == "yes") {
    save.image(file = outputs$environment)
} else {
    x <- tibble(Env_file = "Opted not to store environment. If this is not desired, change 'store_plots_as_rds' to 'yes' in the relevant config file and rerun this rule.")
    write_tsv(x, file = outputs$environment)
}

# figures: modify plot compositions at will!
# load(outputs$environment)
tryCatch(
    {
        library(patchwork)
        paste0("RTE/", names(mysaveandstoreplots))

        ptch <- wrap_plots(mysaveandstoreplots[names(mysaveandstoreplots) %>%
            grep("targetted", ., value = TRUE) %>%
            grep("gsea_top_.*grid_1", ., value = TRUE) %>%
            sort()]) +
            plot_layout(nrow = 1, guides = "collect")
        mysave(pl = ptch, fn = "srna/results/agg/enrichment_analysis/targetted/gsea_top_all.pdf", w = 10 * length(ptch), h = 10)


        ptch <- wrap_plots(mysaveandstoreplots[names(mysaveandstoreplots) %>%
            grep("unbiased", ., value = TRUE) %>%
            grep("gsea_top_.*grid_1.pdf", ., value = TRUE) %>%
            sort()])
        mysave(pl = ptch, fn = "srna/results/agg/enrichment_analysis/unbiased/gsea_top_all.pdf", w = 30, h = 30)

        # ptch <- wrap_plots(mysaveandstoreplots[
        #     c(
        #     names(mysaveandstoreplots) %>% grep("targetted", ., value = TRUE) %>% grep("gsea_top_.*grid", ., value = TRUE) %>% sort(),
        #     names(mysaveandstoreplots) %>% grep("targetted", ., value = TRUE) %>% grep("gsea_top_.*grid", ., value = TRUE) %>% sort()
        #     )
        # ])
        # mysave(pl = ptch, fn = "srna/results/agg/enrichment_analysis/unbiased/gsea_top_all.pdf", w = 30, h = 30)



        p1 <- mysaveandstoreplots[["srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_rte_subfamily_ALL.pdf"]]
        p2 <- mysaveandstoreplots[["srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_rte_subfamily_rte_length_req.pdf"]]



        names(mysaveandstoreplots)
        ptch <- p1 + p2 + plot_layout(ncol = 2, guides = "collect")
        mysaveandstore(pl = ptch, fn = "srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_combined.pdf", w = 6, h = 10)
    },
    error = function(e) {

    }
)

# for (name in names(mysaveandstoreplots)) {
#     tryCatch({
#         print(mysaveandstoreplots[[name]])
#         print(name)
#     },
#     error = function(e) {
#         print(name)
#         print("FAIL")
#     }
#     )}

# print(mysaveandstoreplots[[1]])

# mysaveandstoreplots[["srna/results/agg/enrichment_analysis/condition_ESEN_vs_PRO/gsea/msigdbH/core_enrichments/barplot_HALLMARK_INFLAMMATORY_RESPONSE.pdf"]]
# names(mysaveandstoreplots)[[1]]
# str(mysaveandstoreplots[[1]])
# print(mysaveandstoreplots[[1]])
# mysave(pl = mysaveandstoreplots[[1]])

# str(mysaveandstoreplots[[2]])
# print(mysaveandstoreplots[[2]])
# mysave(pl = mysaveandstoreplots[[2]])

# str(mysaveandstoreplots[[100]])
# print(mysaveandstoreplots[[100]])
# mysave(pl = mysaveandstoreplots[[100]])
