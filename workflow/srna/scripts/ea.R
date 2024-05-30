module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")

source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")

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
library(dplyr)
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
                "SenMayoHuman" = conf$SenMayoHuman,
                "genesets_for_heatmaps" = conf$genesets_for_heatmaps,
                "collections_for_gsea" = conf$collections_for_gsea,
                "inputdir" = "srna/results/agg/deseq2/featurecounts_genes",
                "outputdir" = "srna/results/agg/enrichment_analysis"
            ), env = globalenv())
            assign("inputs", list(
                resultsdf = paste0("srna/results/agg/deseq/resultsdf.tsv")
            ), env = globalenv())
            assign("outputs", list(outfile = "srna/results/agg/enrichment_analysis/outfile.txt"), env = globalenv())
        }
    )

}


# load results
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
resultsdf1 <- resultsdf1[resultsdf1$gene_id != "__no_feature", ]
res <- resultsdf1 %>% filter(tecounttype == tecounttype[1])
res <- res %>% filter(gene_or_te == "gene")
res %>%
    filter(str_detect(gene_id, "CDKN1A")) %>%
    select(conf$samples)

contrast_label_map <- tibble(contrast = params[["contrasts"]], label = gsub("_", " ", gsub("condition_", "", params[["contrasts"]])))

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
                                scale_fill_gradient(high = "red", low = "white", limits = c(0,5), oob = scales::oob_squish) +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL)+ labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)
                        
                            p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust)*`sign(NES)`)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL)+ guides(fill=guide_legend(title="Signed \n-log10(p.adjust)")) + labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                        
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )

                tryCatch(
                    {
                        for (genesetid in conf$genesets_for_gseaplot)
                        p <- gseaplot2(gse, geneSetID = genesetid, pvalue_table = TRUE, subplots = 1:2, ES_geom = "line")
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

gres <- gse_df %>% tibble()
for (collec in genecollections) {
    grestemp <- gres %>% filter(collection == collec) %>% left_join(contrast_label_map)
    sigIDs <- grestemp %>% group_by(contrast) %>% arrange(p.adjust) %>% slice_head(n = 10) %$% ID %>% unique()
    p <- grestemp %>% dplyr::filter(ID %in% sigIDs) %>% mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
        mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
        mutate(label = str_wrap(as.character(label) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
        mutate(label = factor(label, levels = conf$levels)) %>%
        ggplot(aes(x = label, y = ID)) + 
        geom_tile(aes(fill = NES), color = "black") + 
        theme(legend.position = "none") + 
        scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1) + 
        mtclosed + 
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
        labs(x = "", y = "", title = collec) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_equal()
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/gsea_top_%s_grid.pdf", collec), 7,12)
}
for (collec in genecollections) {
    grestemp <- gres %>% filter(collection == collec)
    sigIDs <- grestemp %>% group_by(contrast) %>% arrange(p.adjust) %>% slice_head(n = 10) %$% ID %>% unique()
    p <- grestemp %>% dplyr::filter(ID %in% sigIDs) %>% mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
        mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
        mutate(label = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
        mutate(label = factor(label, levels = conf$levels)) %>%
        ggplot(aes(x = label, y = ID)) + 
        geom_point(aes(color = NES, size = abs(NES))) + 
        theme(legend.position = "none") + 
        scale_color_paletteer_c("grDevices::RdYlBu", direction = -1) + 
        mtclosed + 
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
        labs(x = "", y = "", title = collec) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_equal()
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/gsea_top_%s_dot.pdf", collec), 7,12)
}

# plot core enrichments 
n_top_sets <- 5
for (contrast in params[["contrasts"]]) {
    contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
    contrast_padj <- paste0("padj_", contrast)
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
        for (geneset in names(core_enrichments)) {
            genestoplot <- core_enrichments[[geneset]]
            set_title <- geneset

            contrastconditions <- gsub("condition_", "", contrast) %>% str_split("_vs_") %>% pluck(1)
            sample_vec <- sample_table %>%
                filter(condition %in% contrastconditions) %>%
                pull(sample_name)
            condition_vec <- sample_table %>%
                filter(sample_name %in% sample_vec) %>%
                pull(condition)

            ngenes <- length(genestoplot)
            ncol <- 4
            barplotpf <- res %>% filter(gene_id %in% genestoplot) %>% pivot_longer(cols = conf$samples, names_to = "sample_name", values_to = "counts") %>% left_join(sample_table) %>%
                mutate(sample_name = factor(sample_name, levels = conf$samples)) %>%
                mutate(condition = factor(condition, levels = conf$levels))
            p <- barplotpf %>% ggbarplot(x = "condition", y = "counts", fill = "condition", add = c("mean_se", "dotplot"),
                facet.by = "gene_id", scales = "free", ncol = 4) + 
            labs(title = set_title, x = element_blank()) + 
            mtclosedgridh + scale_conditions + anchorbar +
            theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
            tryCatch(
                {
                    mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/barplot_%s.pdf", params[["outputdir"]], contrast, collec, set_title),8, 2*round(length(genestoplot)/ncol))
                },
                error = function(e) {
                    print(e)
                }
            )

            all_genes_in_set <- read.gmt(params[["collections_for_gsea"]][[collection]]) %>% filter(term == geneset) %$% gene
            pvppres <- res %>% filter(gene_id %in% all_genes_in_set) %>% mutate(logpadj = ifelse(-log10(!!sym(contrast_padj)) > 100, 100,-log10(!!sym(contrast_padj)) ))
            pvppres %$% logpadj
            pvprestopsig <- pvppres %>% arrange(-abs(!!sym(contrast_l2fc))) %>% slice_head(n = 15) %$% gene_id
            p <- pvppres  %>% ggscatter(x = contrast_l2fc, y = "logpadj", label = "gene_id", shape = 1, label.select = pvprestopsig, repel = TRUE) + mtclosed +
                labs(title = set_title, x = "log2fc", y = "-log10(padj)", subtitle = contrast_string) + scale_palette + scale_y_continuous(expand = expansion(mult = c(0, .075)))
            mysaveandstore(sprintf("%s/%s/gsea/%s/all_genes/scatter_%s.pdf", params[["outputdir"]], contrast, collec, set_title), w = 6, h = 6, res = 300)


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
            conditions <- sample_table[match(colnames(m), sample_table$sample_name),]$condition
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
            p <- draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/heatmap_%s.pdf", params[["outputdir"]], contrast, collection, set_title), w = min(dim(scaledm)[2], 12), h = min(dim(scaledm)[1], 12), res = 300)
        }
    }
}


for (geneset in names(params[["genesets_for_heatmaps"]])) {
    # geneset <- geneset
    genestoplot <- read_csv(params[["genesets_for_heatmaps"]][[geneset]], col_names = FALSE) %>% pull(X1)
    set_title <- geneset

    
    # first for all conditions
    #make data long using samples
    ngenes <- length(genestoplot)
    ncol <- 4
    barplotpf <- res %>% filter(gene_id %in% genestoplot) %>% pivot_longer(cols = conf$samples, names_to = "sample_name", values_to = "counts") %>% left_join(sample_table) %>%
        mutate(sample_name = factor(sample_name, levels = conf$samples)) %>%
        mutate(condition = factor(condition, levels = conf$levels))
    p <- barplotpf %>% ggplot() + geom_bar(aes(x = sample_name, y = counts, fill = condition), stat = "identity") + 
    labs(title = set_title, x = element_blank()) + 
    facet_wrap(~gene_id, ncol = 4, scales = "free") +
    mtclosedgridh + scale_conditions + anchorbar +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    mysaveandstore(sprintf("%s/barplot_%s.pdf", params[["outputdir"]], set_title), 8, 2*round(ngenes/ncol))

    heatmapprep <- res %>% filter(gene_id %in% genestoplot)
    m <- as.matrix(heatmapprep %>% dplyr::select(conf$samples))
    rownames(m) <- heatmapprep %$% gene_id
    scaledm <- t(scale(t(m))) %>% na.omit()

    conditions <- sample_table[match(colnames(m), sample_table$sample_name),]$condition
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
    p <- draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    hmscalew <- 0.7
    hmscaleh <- 0.7
    mysaveandstore(sprintf("%s/heatmap_%s.pdf", params[["outputdir"]], set_title), w = min(dim(scaledm)[2]*hmscalew, 12), h = min(dim(scaledm)[1]*hmscaleh, 12), res = 300)

    # now for each contrast
    for (contrast in params[["contrasts"]]) {
        contrast_stat <- paste0("stat_", contrast)
        contrast_l2fc <- paste0("log2FoldChange_", contrast)
        contrast_padj <- paste0("padj_", contrast)
        res <- res %>% arrange(-!!sym(contrast_stat))

        contrastconditions <- gsub("condition_", "", contrast) %>% str_split("_vs_") %>% pluck(1)
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
        conditions <- sample_table[match(colnames(m), sample_table$sample_name),]$condition
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
        p <- draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        mysaveandstore(sprintf("%s/%s/heatmap_%s.pdf", params[["outputdir"]], contrast, set_title), w = min(dim(scaledm)[2]*hmscalew, 12), h = min(dim(scaledm)[1]*hmscaleh, 12), res = 300)
    }
}


# GSEA untargeted 
gene_sets = msigdbr(species = confALL$aref$species)
rm(gse_df)
for (contrast in params[["contrasts"]]) {
    # PREP RESULTS FOR GSEA
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
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
                msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
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
                                slice(1:num)
                            dftemp <- dftemp %>% mutate(Description = str_wrap(as.character(Description) %>% gsub("_", " ", .), width = 40)) 
                            dftemp <- dftemp %>% mutate(`Gene Ratio` = 0.01 * as.numeric(gsub("%", "", gsub(",.*", "", gsub("tags=", "", leading_edge)))))
                            p <- ggplot(dftemp, aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient(high = "red", low = "white", limits = c(0,5), oob = scales::squish) +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL)+ labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/nes%s.pdf", params[["outputdir"]], contrast, collection, num), w = 8, h = min(num, 7), res = 300)
                        
                            p <- ggplot(dftemp, aes(`Gene Ratio`, fct_reorder(Description, `Gene Ratio`), fill = -log10(p.adjust)*`sign(NES)`)) +
                                geom_col(orientation = "y") +
                                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                                mtopen +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL)+ guides(fill=guide_legend(title="Signed \n-log10(p.adjust)")) + labs(caption = contrast_label_map %>% filter(contrast == !!contrast) %$% label)

                            mysaveandstore(sprintf("%s/%s/gsea/%s/dot%s.pdf", params[["outputdir"]], contrast, collection, num), w = 7.5, h = min(num, 10), res = 300)
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )

                tryCatch(
                    {
                        for (genesetid in conf$genesets_for_gseaplot)
                        p <- gseaplot2(gse, geneSetID = genesetid, pvalue_table = TRUE, subplots = 1:2, ES_geom = "line")
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
    contrast_string <- gsub("condition_", "", contrast) %>% str_replace_all("_vs_", " vs ")
    contrast_stat <- paste0("stat_", contrast)
    contrast_l2fc <- paste0("log2FoldChange_", contrast)
    contrast_padj <- paste0("padj_", contrast)
    contrast_significance <- paste0("Significance_", contrast)
    res <- res %>% arrange(-!!sym(contrast_stat))
    for (collec in gse_df %$% collection %>% unique()) {
        gse <- gse_df %>% filter(collection == collec)
        df <- arrange(gse, -abs(NES)) %>%
            group_by(sign(NES)) %>%
            slice(1:n_top_sets)
        ids <- df$ID
        for (geneset in ids) {
            genestoplot <- gene_sets %>% filter(gs_name == geneset) %>% pull(gene_symbol) %>% unique()
            set_title <- geneset

            contrastconditions <- gsub("condition_", "", contrast) %>% str_split("_vs_") %>% pluck(1)
            sample_vec <- sample_table %>%
                filter(condition %in% contrastconditions) %>%
                pull(sample_name)
            condition_vec <- sample_table %>%
                filter(sample_name %in% sample_vec) %>%
                pull(condition)

            pvppres <- res %>% filter(gene_id %in% genestoplot) %>% mutate(logpadj = ifelse(-log10(!!sym(contrast_padj)) > 100, 100,-log10(!!sym(contrast_padj)) ))
            pvppres %$% logpadj
            pvprestopsig <- pvppres %>% arrange(-abs(!!sym(contrast_l2fc))) %>% slice_head(n = 15) %$% gene_id
            p <- pvppres %>% ggscatter(x = contrast_l2fc, y = "logpadj", label = "gene_id", shape = 1, label.select = pvprestopsig, repel = TRUE) + mtclosed +
                labs(title = set_title, x = "log2fc", y = "-log10(padj)", subtitle = contrast_string) + scale_palette + scale_y_continuous(expand = expansion(mult = c(0, .075)))
            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/scatter_%s.pdf", params[["outputdir"]], contrast, collection, set_title), w = 6, h = 6, res = 300)


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
            conditions <- sample_table[match(colnames(m), sample_table$sample_name),]$condition
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
            p <- draw(hm, annotation_legend_list = list(lgd_pvalue_adj, lgd_sig), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
            mysaveandstore(sprintf("%s/%s/gsea/%s/core_enrichments/heatmap_%s.pdf", params[["outputdir"]], contrast, collection, set_title), w = min(dim(scaledm)[2], 12), h = min(dim(scaledm)[1], 12), res = 300)
        }
    }
}

gres <- gse_df %>% tibble()
for (collec in gse_df %$% collection %>% unique()) {
    grestemp <- gres %>% filter(collection == collec)
    sigIDs <- grestemp %>% group_by(contrast) %>% arrange(p.adjust) %>% slice_head(n = 10) %$% ID %>% unique()
    p <- grestemp %>% dplyr::filter(ID %in% sigIDs) %>% mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
        mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
        mutate(label = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
        mutate(label = factor(label, levels = conf$levels)) %>%
        ggplot(aes(x = label, y = ID)) + 
        geom_tile(aes(fill = NES), color = "black") + 
        theme(legend.position = "none") + 
        scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1) + 
        mtclosed + 
        theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
        labs(x = "", y = "", title = collec) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_equal()
    mysaveandstore(sprintf("srna/results/agg/enrichment_analysis/gsea_top_%s_grid.pdf", collec), 7,12)
}

save(mysaveandstoreplots, file = outputs$plots)