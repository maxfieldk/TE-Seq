source("~/data/common/myDefaults.r")
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

# analysis parameters
{
    conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]


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
                "inputdir" = "results/agg/repeatanalysis_telescope",
                "outputdir" = "results/agg/enrichment_analysis_repeats",
                "tecounttypes" = c("telescope_multi"),
                "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
                "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
            ), env = globalenv())
            assign("inputs", list(
                "resultsdf" = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
            ), env = globalenv())
            assign("outputs", list(outfile = "results/agg/enrichment_analysis_repeats/outfile.txt"), env = globalenv())
        }
    )

    sample_table <- read_csv(params[["sample_table"]])
    peptable <- read.csv(conf$peptable)
}



## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
resultsdf1 <- resultsdf1[resultsdf1$gene_id != "__no_feature", ]
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
resultsdf <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)



gse_results <- list()
core_enrichments_for_plot <- list()

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



gse_results <- list()
core_enrichments_for_plot <- list()
EARTEplots <- list()

for (contrast in params[["contrasts"]]) {
    contrast_of_interest <- contrast
    contrast_level_1 <- contrast_of_interest %>%
        str_split("_") %>%
        unlist() %>%
        .[4]
    contrast_level_2 <- contrast_of_interest %>%
        str_split("_") %>%
        unlist() %>%
        .[2]
    contrast_stat <- paste0("stat_", contrast_of_interest)
    contrast_padj <- paste0("padj_", contrast_of_interest)
    contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)
    contrast_samples <- peptable %>%
        filter(condition %in% c(contrast_level_1, contrast_level_2)) %>%
        pull(sample_name)
    condition_vec <- peptable %>% filter(sample_name %in% contrast_samples) %$% condition

    for (tecounttype in resultsdf$tecounttype %>% unique()) {
        res <- resultsdf %>%
            filter(tecounttype == tecounttype) %>%
            arrange(-!!sym(contrast_stat)) %>%
            filter(gene_or_te == "repeat")
        ordered_by_stat <- setNames(res %>% pull(!!sym(contrast_stat)), res$gene_id) %>% na.omit()

        for (ontology in ontologies) {
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
            eligible_facet_modifiers <- c(eligible_modifiers[grepl("_loc$", eligible_modifiers)], "ALL")
            eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)

            for (filter_var in eligible_filter_modifiers) {
                tryCatch(
                    {
                        if (filter_var != "ALL") {
                            genesets <- resultsdf %>%
                                filter(!!sym(ontology) != "Other") %>%
                                filter(str_detect(!!sym(filter_var), ">|Intact")) %>%
                                select(!!sym(filter_var), gene_id)
                        } else {
                            genesets <- resultsdf %>%
                                select(!!sym(ontology), gene_id) %>%
                                filter(!!sym(ontology) != "Other")
                        }

                        gse <- GSEA(ordered_by_stat, TERM2GENE = genesets, maxGSSize = 10000, minGSSize = 1)
                        gse_results[[contrast]][[tecounttype]][[ontology]][[filter_var]] <- gse
                        genesettheme <- theme_gray() + theme(axis.text.y = element_text(colour = "black"))
                        p <- dotplot(gse, showCategory = 20) + ggtitle(paste("GSEA", contrast, sep = " ")) + genesettheme + mytheme
                        mysave(sprintf("%s/%s/%s/gsea/%s/%s/dotplot.png", params[["outputdir"]], tecounttype, contrast, ontology, filter_var), w = 4, h = 6, res = 300)
                        EARTEplots[[contrast]][[tecounttype]][[ontology]][[filter_var]][["dot"]] <- p

                        p <- ridgeplot(gse, core_enrichment = FALSE) + ggtitle(paste("GSEA", contrast, sep = " ")) + xlab("Log2 FC") + xlim(c(-4, 4)) + genesettheme + mytheme
                        mysave(sprintf("%s/%s/%s/gsea/%s/%s/ridgeplot.png", params[["outputdir"]], tecounttype, contrast, ontology, filter_var), w = 4, h = 6, res = 300)
                        EARTEplots[[contrast]][[tecounttype]][[ontology]][[filter_var]][["ridge"]] <- p


                        tryCatch(
                            {
                                for (num in c(5, 10, 15, 30)) {
                                    df <- arrange(gse, -abs(NES)) %>%
                                        group_by(sign(NES)) %>%
                                        slice(1:num)
                                    df <- df@result

                                    p <- ggplot(df, aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
                                        geom_col(orientation = "y") +
                                        scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
                                        mytheme +
                                        theme(axis.text.y = element_text(colour = "black")) +
                                        ylab(NULL)

                                    mysave(sprintf("%s/%s/%s/gsea/%s/%s/nes%s.png", params[["outputdir"]], tecounttype, contrast, ontology, filter_var, num), w = 4, h = min(num, 7), res = 300)
                                    EARTEplots[[contrast]][[tecounttype]][[ontology]][[filter_var]][["nes"]][[num]] <- p
                                }
                            },
                            error = function(e) {
                                print("")
                            }
                        )
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        }
    }
}

save(EARTEplots, file = sprintf("%s/EARTEplots.RData", params[["outputdir"]]))
save(gse_results, file = sprintf("%s/gse_results.RData", params[["outputdir"]]))

x <- data.frame()
write.table(x, file = outputs[["outfile"]], col.names = FALSE)
