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
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(paletteer)
library(forcats)
library(ggstance)
library(enrichplot)
library(circlize)
library(ComplexHeatmap)
library(patchwork)



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
                "counttype" = "telescope_multi",
                "sample_table" = conf$sample_table,
                "contrasts" = conf$contrasts,
                "SenMayoHuman" = conf$SenMayoHuman,
                "genesets_for_heatmaps" = conf$genesets_for_heatmaps,
                "collections_for_gsea" = conf$collections_for_gsea,
                "inputdir" = "srna/results/agg/deseq",
                "outputdir" = "srna/results/agg/enrichment_analysis_repeats/telescope_multi",
                "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
                "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
            ), env = globalenv())
            assign("inputs", list(
                resultsdf = paste0("srna/results/agg/deseq/resultsdf.tsv")
            ), env = globalenv())
            assign("outputs", list(
                "environment" = "srna/results/agg/enrichment_analysis_repeats/telescope_multi/enrichment_analysis_repeats_environment.RData",
                "results_table" = "srna/results/agg/enrichment_analysis_repeats/telescope_multi/results_table.tsv"
            ), env = globalenv())
        }
    )
}

counttype <- params$counttype

## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
resultsdf1 <- resultsdf1[resultsdf1$gene_id != "__no_feature", ]
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
resultsdf <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)


# resultsdf %>% filter(grepl("L1HS", gene_id))


### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "_.*family")]
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





rm(gse_df)
for (contrast in params[["contrasts"]]) {
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

    res <- resultsdf %>%
        filter(counttype == !!counttype) %>%
        arrange(-!!sym(contrast_stat)) %>%
        filter(gene_or_te == "repeat")
    ordered_by_stat <- setNames(res %>% pull(!!sym(contrast_stat)), res$gene_id) %>% na.omit()

    for (ontology in small_ontologies) {
        print(ontology)
        ontology_groups <- r_repeatmasker_annotation %>%
            pull(!!sym(ontology)) %>%
            unique()
        ontology_groups <- ontology_groups[ontology_groups != "Other"]

        eligible_modifiers <- c()
        for (modifier in modifiers) {
            values_present <- resultsdf %>%
                filter(!!sym(ontology) != "Other") %>%
                pull(!!sym(modifier)) %>%
                unique()

            if (!("Other" %in% values_present)) {
                eligible_modifiers <- c(eligible_modifiers, modifier)
            }
        }
        eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
        eligible_facet_modifiers <- c(eligible_modifiers[grepl("genic_loc$", eligible_modifiers)], "ALL")
        eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)

        for (filter_var in eligible_filter_modifiers) {
            tryCatch(
                {
                    if (filter_var != "ALL") {
                        genesets <- resultsdf %>%
                            filter(!!sym(ontology) != "Other") %>%
                            filter(!!sym(filter_var) == "Intact" | !!sym(filter_var) == "FL") %>%
                            select(!!sym(ontology), gene_id)
                    } else {
                        genesets <- resultsdf %>%
                            select(!!sym(ontology), gene_id) %>%
                            filter(!!sym(ontology) != "Other")
                    }


                    # res %>% mutate(nrow = row_number()) %>% filter(grepl("L1HS", gene_id)) %>% filter(rte_length_req == "FL") %>% dplyr::select(gene_id, nrow) %>% print(n = 800
                    # )

                    gsestd <- GSEA(ordered_by_stat, TERM2GENE = genesets, maxGSSize = 10000, minGSSize = 1, pvalueCutoff = 0.05, seed = TRUE, scoreType = "std")
                    gsepos <- GSEA(ordered_by_stat, TERM2GENE = genesets, maxGSSize = 10000, minGSSize = 1, pvalueCutoff = 0.05, seed = TRUE, scoreType = "pos")
                    gseneg <- GSEA(ordered_by_stat, TERM2GENE = genesets, maxGSSize = 10000, minGSSize = 1, pvalueCutoff = 0.05, seed = TRUE, scoreType = "neg")
                    dfstd <- gsestd@result %>% tibble()
                    dfstd$test_type <- "std"
                    dfpos <- gsepos@result %>% tibble()
                    dfpos$test_type <- "pos"
                    dfneg <- gseneg@result %>% tibble()
                    dfneg$test_type <- "neg"
                    df <- bind_rows(bind_rows(dfstd, dfpos), dfneg)
                    df$collection <- ontology
                    df$contrast <- contrast
                    df$filter_var <- filter_var
                    if (!exists("gse_df")) {
                        gse_df <<- df
                    } else {
                        gse_df <<- rbind(gse_df, df)
                    }

                    # p <- dotplot(gse, showCategory = 20) + ggtitle(paste("GSEA", contrast, sep = " ")) + genesettheme + mtopen
                    # mysaveandstore(sprintf("%s/%s/%s/gsea/%s/%s/dotplot.pdf", params[["outputdir"]], counttype, contrast, ontology, filter_var), w = 3, h = 4, res = 300)

                    # p <- ridgeplot(gse, core_enrichment = FALSE) + ggtitle(paste("GSEA", contrast, sep = " ")) + xlab("Log2 FC") + xlim(c(-4, 4)) + genesettheme + mtopen
                    # mysaveandstore(sprintf("%s/%s/%s/gsea/%s/%s/ridgeplot.pdf", params[["outputdir"]], counttype, contrast, ontology, filter_var), w = 3, h = 4, res = 300)


                    # tryCatch(
                    #     {
                    #         for (num in c(5, 10, 15, 30)) {
                    #             df <- arrange(gse, -abs(NES)) %>%
                    #                 group_by(sign(NES)) %>%
                    #                 slice(1:num)
                    #             df <- df@result

                    #             p <- ggplot(df, aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
                    #                 geom_col(orientation = "y") +
                    #                 scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
                    #                 mtopen +
                    #                 theme(axis.text.y = element_text(colour = "black")) +
                    #                 ylab(NULL)

                    #             mysaveandstore(sprintf("%s/%s/%s/gsea/%s/%s/nes%s.pdf", params[["outputdir"]], counttype, contrast, ontology, filter_var, num), w = 3, h = min(num/2, 7), res = 300)
                    #         }
                    #     },
                    #     error = function(e) {
                    #         print("")
                    #     }
                    # )
                },
                error = function(e) {
                    print(e)
                }
            )
        }
    }
}

contrast_label_map <- tibble(contrast = params[["contrasts"]], label = gsub("condition_", "", params[["contrasts"]]))
gres <- gse_df %>% tibble()
for (test_type in c("std", "pos", "neg")) {
    for (ontology in small_ontologies) {
        for (filter_var in gres %$% filter_var %>% unique()) {
            grestemp <- gres %>%
                filter(test_type == !!test_type) %>%
                filter(collection == ontology) %>%
                filter(filter_var == !!filter_var) %>%
                filter(grepl(paste0(conf$levels[1], "$"), contrast)) %>%
                left_join(contrast_label_map)
            sigIDs <- grestemp %>%
                mutate(direction = ifelse(NES > 0, "UP", "DOWN")) %>%
                group_by(contrast, direction) %>%
                arrange(p.adjust) %>%
                slice_head(n = 5) %$% ID %>%
                unique()

            # human_subfam_ordering <- c(
            #     "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6",
            #     "AluY", "HERVK_LTR", "HERVK_INT", "HERVL_LTR", "HERVL_INT",
            #     "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"
            # )
            # other_ordering <- resultsdf %>%
            #     filter(rte_subfamily != "Other") %>%
            #     group_by(rte_family, rte_subfamily) %>%
            #     summarise(family_av_pctdiv = mean(family_av_pctdiv)) %>%
            #     mutate(rte_family = factor(rte_family, levels = c("L1", "Alu", "ERV", "SVA"))) %>%
            #     arrange(rte_family, family_av_pctdiv) %$% rte_subfamily
            # if (confALL$aref$species != "human") {
            #     subfam_ordering <<- other_ordering
            # } else {
            #     subfam_ordering <<- human_subfam_ordering
            # }
            p <- grestemp %>%
                dplyr::filter(ID %in% sigIDs) %>%
                mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
                mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(contrast = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
                mutate(contrast = factor(contrast, levels = conf$levels)) %>%
                ggplot(aes(x = contrast, y = ID)) +
                geom_tile(aes(fill = NES), color = "black") +
                theme(legend.position = "none") +
                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                mtclosed +
                theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
                labs(x = "", y = "", title = paste0(ontology, " ", filter_var, " ", subtitle = test_type)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                coord_equal()
            mysaveandstore(sprintf("%s/%s/gsea_top_rtes_%s_%s_TOP5EITHERDIRECTIONSHOWN.pdf", params[["outputdir"]], test_type, ontology, filter_var), w = 3 + .3 * (length(conf$levels) - 1), h = 3 + .3 * length(sigIDs))

            sigIDs <- grestemp %>%
                mutate(direction = ifelse(NES > 0, "UP", "DOWN")) %>%
                group_by(contrast, direction) %>%
                filter(p.adjust <= 0.05) %>%
                arrange(p.adjust) %$%
                ID %>%
                unique()
            p <- grestemp %>%
                dplyr::filter(ID %in% sigIDs) %>%
                mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
                mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(contrast = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
                mutate(contrast = factor(contrast, levels = conf$levels)) %>%
                ggplot(aes(x = contrast, y = ID)) +
                geom_tile(aes(fill = NES), color = "black") +
                theme(legend.position = "none") +
                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                mtclosed +
                theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
                labs(x = "", y = "", title = paste0(ontology, " ", filter_var, " ", subtitle = test_type)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                coord_equal()
            mysaveandstore(sprintf("%s/%s/gsea_top_rtes_%s_%s_ALLSIGSHOWN.pdf", params[["outputdir"]], test_type, ontology, filter_var), w = 3 + .3 * (length(conf$levels) - 1), h = 3 + .3 * length(sigIDs))

            all_levels <- resultsdf[, ontology] %>%
                unlist() %>%
                unique()
            all_levels <- all_levels[!is.na(all_levels)]
            all_levels <- all_levels[!(all_levels == "Other")]
            all_levels_df <- tibble(Description := all_levels)
            p <- grestemp %>%
                mutate(sig = ifelse(p.adjust < 0.05, "*", "")) %>%
                mutate(ID = str_wrap(as.character(ID) %>% gsub("_", " ", .), width = 40)) %>%
                mutate(contrast = str_wrap(as.character(contrast) %>% gsub("condition_", "", .) %>% gsub("_vs_.*", "", .), width = 40)) %>%
                mutate(contrast = factor(contrast, levels = conf$levels)) %>%
                ggplot(aes(x = contrast, y = ID)) +
                geom_tile(aes(fill = NES), color = "black") +
                theme(legend.position = "none") +
                scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
                mtclosed +
                theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", hjust = 1)) +
                labs(x = "", y = "", title = paste0(ontology, " ", filter_var, " ", subtitle = test_type)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                coord_equal()
            mysaveandstore(sprintf("%s/%s/gsea_top_rtes_%s_%s_ALLSHOWN.pdf", params[["outputdir"]], test_type, ontology, filter_var), w = 3 + .3 * (length(conf$levels) - 1), h = 3 + .3 * length(sigIDs))
        }
    }
}

gres %>% write_tsv(outputs$results_table)

x <- tibble(OUT = "")
write_tsv(x, file = outputs$environment)

# figures: modify plot compositions at will!
tryCatch(
    {
        library(patchwork)
        p1 <- mysaveandstoreplots[["srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_rte_subfamily_ALL.pdf"]]
        p2 <- mysaveandstoreplots[["srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_rte_subfamily_rte_length_req.pdf"]]
        names(mysaveandstoreplots)
        ptch <- p1 + p2 + plot_layout(ncol = 2, guides = "collect")
        mysaveandstore(pl = ptch, fn = "srna/results/agg/enrichment_analysis_repeats/telescope_multi/gsea_top_rtes_combined.pdf", w = 6, h = 10)
        print(mysaveandstoreplots[[1]])
    },
    error = function(e) {

    }
)
