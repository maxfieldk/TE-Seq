source("~/data/common/myDefaults.r")
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

conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]


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
            "inputdir" = "results/agg/deseq_telescope",
            "outputdir" = "results/agg/repeatanalysis_telescope",
            "tecounttypes" = c("telescope_multi"),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "results/agg/repeatanalysis_telescope/plots.outfile.txt"
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

#### PLOTTING FUNCTIONS

pvp <- function(df, facet_var = "ALL", filter_var = "ALL") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), ">|Intact|towards"))
    }
    p <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
        {
            ggplot(data = .) +
                geom_point(aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2), color = padj < 0.05)) +
                scale_color_manual(values = c("black", "red", "lightgray")) +
                geom_abline(intercept = 0, slope = 1) +
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2)) +
                mytheme +
                coord_fixed(
                    xlim = range(c(.[[contrast_level_1]], .[[contrast_level_2]])),
                    ylim = range(c(.[contrast_level_1], .[contrast_level_2]))
                ) +
                theme(aspect.ratio = 1)
        }
    if (facet_var != "ALL") {
        p <- p + facet_wrap(facet_var)
    }
    return(p)
}

# df <- tidydf %>% filter(rte_subfamily == "L1HS")
# p <- pvp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "l1_intactness_req", facet_var = "genic_loc") + ggtitle("L1HS")
# mysave("temp1.png")

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
            labs(x = "", y = "Number DE") +
            geom_col(aes(x = direction, fill = direction, y = count)) +
            anchorbar +
            mythemecontrast +
            guides(fill = "none")
    } else {
        p <- plotframe %>% ggplot() +
            labs(x = "", y = "Number DE") +
            geom_col(aes(x = direction, fill = direction, y = count)) +
            facet_wrap(facet_var) +
            anchorbar +
            mythemecontrast +
            guides(fill = "none")
    }
    return(p)
}

# need to wrap functions cALLs in try catch since filtering and checking for DE status can mean there are zero elements
# p <- dep(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "rte_length_req", facet_var = "genic_loc") + ggtitle("L1HS")
# mysave("temp1.png")
# p <- dep(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "rte_length_req") + ggtitle("L1HS")
# mysave("temp1.png")

stripp <- function(df, stats = "yes", extraGGoptions = NULL, facet_var = "ALL", filter_var = "ALL") {
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
        labs(x = "", y = "Sum Normalized Counts") +
        mythemecontrastrev +
        extraGGoptions +
        theme(legend.position = "none") +
        mytheme +
        anchorbar +
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

    topAnn <- HeatmapAnnotation(Condition = condition_vec, col = list(Condition = unlist(conf$condition_colors[condition_vec])))
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
                column_names_rot = 45,
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



#### PLOTTING
plots <- list()
for (contrast in contrasts) {
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
                groupframe <- tidydf %>% filter(!!sym(ontology) == group)
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
                    eligible_facet_modifiers <- c(eligible_modifiers[grepl("_loc$", eligible_modifiers)], "ALL")
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
                                plot_height <- 4
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
                                mysave(sprintf("%s/%s/%s/%s/%s_%s_%s.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                plots[[tecounttype]][[contrast]][[group]][[function_name]][[filter_var]][[facet_var]] <- p
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

# Heatmaps



for (contrast in contrasts) {
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
                    eligible_facet_modifiers <- c(eligible_modifiers[grepl("_loc$", eligible_modifiers)], "ALL")
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
                                        mysave(sprintf("%s/%s/%s/%s/%s_%s_%s_%s_%s.png", outputdir, tecounttype, contrast, function_name, group, filter_var, facet_var, DEvar, scaled), plot_width, plot_height)
                                        plots[[tecounttype]][[contrast]][[group]][[function_name]][[filter_var]][[facet_var]][[DEvar]][[scaled]] <- p
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




### VENN DIAGRAMS

tryCatch(
    {
        vennplots <- list()
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
                            contrast_level_1 <- contrast_of_interest %>%
                                str_split("_") %>%
                                unlist() %>%
                                .[4]
                            contrast_level_2 <- contrast_of_interest %>%
                                str_split("_") %>%
                                unlist() %>%
                                .[2]
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
                            mysave(sprintf("%s/%s/%s/ggVenn_%s_%s_%s.png", outputdir, tecounttype, modifier, group, modifier, direction), 6, 6)
                            vennplots[[tecounttype]][[group]][[modifier]][[direction]] <- p
                        }
                        mysave(sprintf("%s/%s/ggVenn_%s_%s.png", outputdir, tecounttype, group, direction), 6, 6)
                        vennplots[[tecounttype]][[group]][["unmodified"]][[direction]] <- p
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

save(plots, file = "results/agg/repeatanalysis_telescope/repeatanalysisplots_plots.RData")
tryCatch(
    {
        save(vennplots, file = "results/agg/repeatanalysis_telescope/repeatanalysisplots_vennplots.RData")
    },
    error = function(e) {
        print(e)
        message("Venn diagrams failed - do you only have one contrast?")
    }
)

########
# PROJECT SPECIFIC CODE

tryCatch(
    {
        old_t2t_refnonref_element_annotations <- read_delim("/users/mkelsey/data/ref/genomes/hs1/annotations6TldrDerived/refnonrefl1hspa2intact.bed")
        old_gene_names_gr <- old_t2t_refnonref_element_annotations %>%
            rename(seqnames = chr) %>%
            dplyr::select(gene_id, seqnames, start, end) %>%
            GRanges()

        t2t_updated <- resultsdf %>%
            filter(tecounttype == "telescope_multi") %>%
            filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2") %>%
            GRanges()


        audrey_mapping <- mergeByOverlaps(old_gene_names_gr, t2t_updated) %>%
            as.data.frame() %>%
            tibble()

        write_csv(audrey_mapping, "/users/mkelsey/data/Audrey/forAudrey/refnonrefl1hspa2intact_RNASEQ_results.csv")

        mapper <- audrey_mapping %>%
            dplyr::select(old_gene_names_gr.gene_id, t2t_updated.gene_id) %>%
            rename(old_gene_id = old_gene_names_gr.gene_id, gene_id = t2t_updated.gene_id)

        compartment_rs1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/240129_refnonref_collapse_pro_RS.txt", col_names = FALSE)
        compartment_rs <- compartment_rs1 %>%
            dplyr::select(X4, X22) %>%
            rename(old_gene_id = X4, compartment_rs = X22) %>%
            left_join(mapper) %>%
            select(-old_gene_id)

        compartment_ois1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/240129_refnonref_collapse_grow_RIS.txt", col_names = FALSE)
        compartment_ois <- compartment_ois1 %>%
            dplyr::select(X4, X22) %>%
            rename(old_gene_id = X4, compartment_ois = X22) %>%
            left_join(mapper) %>%
            select(-old_gene_id)

        switches <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/switches_w_dir.txt", delim = "\t", col_names = FALSE) %>%
            rename(subcompartment = X1, switch = X2)
        swtices_score <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/switches_score.tsv", delim = "\t", col_names = FALSE) %>%
            rename(subcompartment = X1, switch_rs_score = X2)
        subcompartment_rs1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/240129_refnonref_pro_sen_sub_collapsed.txt", col_names = FALSE)
        subcompartment_rs <- subcompartment_rs1 %>%
            dplyr::select(X4, X22) %>%
            rename(old_gene_id = X4, subcompartment = X22) %>%
            left_join(mapper) %>%
            select(-old_gene_id) %>%
            left_join(switches) %>%
            left_join(swtices_score) %>%
            rename(switch_rs = switch) %>%
            rename(subcompartment_rs = subcompartment) %>%
            replace_na(list(switch_rs = "no change", hotspots_subcompartment = "not_hotspot_sub"))


        subcompartment_ois1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/240129_refnonref_grow_RIS_collapsed_sub.txt", col_names = FALSE)
        subcompartment_ois <- subcompartment_ois1 %>%
            dplyr::select(X4, X22) %>%
            rename(old_gene_id = X4, subcompartment = X22) %>%
            left_join(mapper) %>%
            select(-old_gene_id) %>%
            left_join(switches) %>%
            rename(switch_ois = switch) %>%
            rename(subcompartment_ois = subcompartment) %>%
            replace_na(list(switch_ois = "no change", hotspots_subcompartment = "not_hotspot_sub"))


        hotspots_compartment1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/l1_hotspots_ab.txt", delim = "\t", col_names = FALSE)
        hotspots_compartment <- hotspots_compartment1 %>%
            dplyr::select(X1) %>%
            rename(old_gene_id = X1) %>%
            left_join(mapper) %>%
            select(-old_gene_id) %>%
            mutate(hotspots_compartment = "hotspot")

        hotspots_subcompartment1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/l1_hotspot_sub.txt", delim = "\t", col_names = FALSE)
        hotspots_subcompartment <- hotspots_subcompartment1 %>%
            dplyr::select(X1) %>%
            rename(old_gene_id = X1) %>%
            left_join(mapper) %>%
            select(-old_gene_id) %>%
            mutate(hotspots_subcompartment = "hotspot_sub")


        loop_S1 <- read_delim("/users/mkelsey/data/Audrey/fromAudrey/l1_S_dist.txt", delim = "\t", col_names = FALSE)
        loop_S <- loop_S1 %>%
            rename(old_gene_id = X4, distance_to_loop = X24) %>%
            mutate(InSenUniqueLoop = ifelse(distance_to_loop > 0, "Outside", "Inside")) %>%
            dplyr::select(old_gene_id, InSenUniqueLoop) %>%
            left_join(mapper) %>%
            select(-old_gene_id)


        # full_join(hotspots_subcompartment) %>%


        audrey_annotations <- full_join(compartment_rs, compartment_ois) %>%
            full_join(hotspots_compartment) %>%
            full_join(subcompartment_rs) %>%
            full_join(subcompartment_ois) %>%
            full_join(loop_S) %>%
            replace_na(list(hotspots_compartment = "not_hotspot", hotspots_subcompartment = "not_hotspot_sub")) %>%
            mutate(motion_comp_rs = ifelse(str_detect(compartment_rs, "towards"), compartment_rs, "no change")) %>%
            mutate(motion_comp_ois = ifelse(str_detect(compartment_ois, "towards"), compartment_ois, "no change"))
        audrey_modifiers <- c("compartment_rs", "motion_comp_rs", "hotspots_compartment", "switch_rs", "switch_rs_score", "InSenUniqueLoop")



        select_plots <- list()
        for (vst in c("NORM")) {
            for (contrast in contrasts) {
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
                groups_that_have_been_run <- c()
                groups_not_to_run <- c("AluY")
                for (group in c("L1HSvL1PA2")) {
                    if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
                        groups_that_have_been_run <- c(groups_that_have_been_run, group)
                        eligible_filter_modifiers <- "ALL"
                        eligible_facet_modifiers <- c(audrey_modifiers, "ALL")
                        eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
                        # first plots without any modifiers
                        if (vst == "VST") {
                            audrey_res <- left_join(audrey_annotations, vstresultsdf %>%
                                filter(tecounttype == "telescope_multi") %>%
                                filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2"))
                            groupframe <- audrey_res
                        } else {
                            audrey_res <- left_join(audrey_annotations, resultsdf %>%
                                filter(tecounttype == "telescope_multi") %>%
                                filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2"))
                            groupframe <- audrey_res
                        }

                        plotting_functions <- c("myheatmap")

                        for (function_name in plotting_functions) {
                            for (DEvar in c("ALL", "DE")) {
                                for (scaled in c("notscaled", "scaled")) {
                                    for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                                        filter_var <- "ALL"
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
                                                mysave(sprintf("%s/%s/%s/%s/%s/%s/%s_%s_%s_%s_%s.png", outputdir, "select_elements", vst, tecounttype, contrast, function_name, group, filter_var, facet_var, DEvar, scaled), plot_width, plot_height)
                                                select_plots[[vst]][[tecounttype]][[contrast]][[group]][[function_name]][[filter_var]][[facet_var]][[DEvar]][[scaled]] <- p
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
                        plotting_functions <- c("stripp", "pvp", "dep")

                        if (vst == "VST") {
                            audrey_tidyres <- left_join(audrey_annotations, vsttidydf %>%
                                filter(tecounttype == "telescope_multi") %>%
                                filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2"))
                            groupframe <- audrey_tidyres
                        } else {
                            audrey_tidyres <- left_join(audrey_annotations, tidydf %>%
                                filter(tecounttype == "telescope_multi") %>%
                                filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2"))
                            groupframe <- audrey_tidyres
                        }
                        for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                            colsToKeep <- c("gene_id", "family", pvals, l2fc, strictly_annotations, audrey_modifiers)
                            facet_var <- eligible_modifier_combinations[i, ]$facet_var
                            if (!!sym(facet_var) != "ALL") {
                                groupframetemp <- groupframe %>% filter(!!sym(facet_var) != ".")
                            } else {
                                groupframetemp <- groupframe
                            }
                            for (function_name in plotting_functions) {
                                tryCatch(
                                    {
                                        function_current <- get(function_name)
                                        plot_width <- 5
                                        plot_height <- 4
                                        if (facet_var != "ALL") {
                                            plot_width <- 8
                                        }
                                        if (filter_var != "ALL") {
                                            plot_title <- groupframetemp %>%
                                                pull(!!sym(filter_var)) %>%
                                                unique() %>%
                                                grep(">|Intact", ., value = TRUE)
                                        } else {
                                            plot_title <- group
                                        }
                                        p <- function_current(groupframetemp, filter_var = filter_var, facet_var = facet_var) + ggtitle(plot_title)

                                        mysave(sprintf("%s/%s/%s/%s/%s/%s/%s_%s_%s.png", outputdir, "select_elements", vst, tecounttype, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                        select_plots[[vst]][[tecounttype]][[contrast]][[group]][[function_name]][[filter_var]][[facet_var]] <- p
                                    },
                                    error = function(e) {
                                        print(sprintf("Error with  %s %s %s %s %s %s", tecounttype, contrast, group, function_name, filter_var, facet_var))
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
    },
    error = function(e) {
        print(e)
        message("Didn't run project specific code")
    }
)

tryCatch(
    {
        contrast_of_interest <- contrasts[1]
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

        audrey_tidyres <- left_join(audrey_annotations, tidydf %>%
            filter(tecounttype == "telescope_multi") %>%
            filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2"))
        groupframe <- audrey_tidyres

        groupframe %$% compartment_rs
        colnames(groupframe)
        pf <- groupframe %>%
            group_by(gene_id, condition) %>%
            summarise(counts = sum(counts), pval = dplyr::first(!!sym(contrast_padj)), switch_rs_score = dplyr::first(switch_rs_score), subcompartment_rs = dplyr::first(subcompartment_rs), subcomp = dplyr::first(switch_rs), loop = dplyr::first(InSenUniqueLoop)) %>%
            ungroup() %>%
            pivot_wider(names_from = condition, values_from = counts) %>%
            mutate(difSmP = SEN - PRO) %>%
            mutate(group = ifelse(subcomp == "no change", str_extract(subcompartment_rs, "A|B"), subcomp)) %>%
            filter(group != "NA") %>%
            mutate(group = factor(group, levels = c("B", "down", "up", "A"), ordered = FALSE)) %>%
            mutate(log2fc = log2((SEN + 1) / (PRO + 1)))



        summarydf <- pf %>%
            group_by(group) %>%
            summarise(
                n = n(),
                mean = mean(difSmP),
                sd = sd(difSmP),
                meanlog2fc = mean(log2fc),
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1)) %>%
            ungroup()


        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = group, y = difSmP, color = loop, shape = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = group, y = mean), shape = 19) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme_bw() +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("strip1.png", 6, 4)

        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = group, y = log2fc, color = loop, shape = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = group, y = meanlog2fc), shape = 19) +
            labs(x = "", y = "log2FC") +
            theme_bw() +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("strip1.2.png", 6, 4)

        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = group, y = difSmP, color = loop, shape = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = group, y = mean), shape = 19) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme_bw() +
            mythemecontrastrev +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # theme to remove grid lines

        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("strip1.1.png", 6, 4)

        model <- lm(difSmP ~ group, data = pf)
        summary(model)
        amodel <- aov(difSmP ~ group, data = pf)
        TukeyHSD(amodel)


        p <- pf %>%
            ggplot() +
            geom_violin(aes(x = group, y = difSmP), draw_quantiles = .50) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp11.png", 4, 4)
        p <- pf %>%
            ggplot() +
            geom_violin(aes(x = group, y = difSmP, fill = loop), draw_quantiles = .50) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp12.png", 6, 4)


        summarydf <- pf %>%
            group_by(group, loop) %>%
            summarise(
                n = n(),
                mean = mean(difSmP),
                sd = sd(difSmP)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1)) %>%
            ungroup()

        p <- pf %>%
            ggplot() +
            geom_violin(aes(x = group, y = difSmP, fill = loop)) +
            geom_point(data = summarydf, aes(x = group, y = mean, group = loop), shape = 19, position = position_dodge(width = 0.9), color = "black") +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp13.png", 6, 4)

        p <- pf %>%
            mutate(loopgroup = paste0(group, "_", loop)) %>%
            ggplot() +
            geom_jitter(aes(x = loopgroup, y = difSmP, color = loop), width = 0.2) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp14.png", 6, 4)

        pf <- groupframe %>%
            group_by(gene_id, condition) %>%
            summarise(counts = sum(counts), pval = dplyr::first(!!sym(contrast_padj)), switch_rs_score = dplyr::first(switch_rs_score), compartment = dplyr::first(compartment_rs), subcomp = dplyr::first(switch_rs), loop = dplyr::first(InSenUniqueLoop)) %>%
            ungroup() %>%
            pivot_wider(names_from = condition, values_from = counts) %>%
            mutate(difSmP = SEN - PRO)

        summarydf <- pf %>%
            group_by(subcomp) %>%
            summarise(
                n = n(),
                mean = mean(difSmP),
                sd = sd(difSmP)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1)) %>%
            ungroup()


        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = subcomp, y = difSmP, color = loop, shape = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = subcomp, y = mean), shape = 19) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp1.png", 6, 4)

        p <- pf %>%
            ggplot() +
            geom_violin(aes(x = subcomp, y = difSmP, fill = loop), draw_quantiles = c(0.5)) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp15.png", 6, 4)
        p <- pf %>%
            ggplot() +
            geom_violin(aes(x = subcomp, y = difSmP), draw_quantiles = c(0.5)) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp16.png", 6, 4)



        summarydf <- pf %>%
            group_by(switch_rs_score) %>%
            summarise(
                n = n(),
                mean = mean(difSmP),
                sd = sd(difSmP)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1)) %>%
            ungroup()

        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = switch_rs_score, y = difSmP, color = loop, shape = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = switch_rs_score, y = mean), shape = 19) +
            geom_smooth(aes(x = switch_rs_score, y = difSmP), method = "lm", se = FALSE) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp1.1.png", 6, 4)

        p <- pf %>%
            mutate(Significance = ifelse(is.na(pval), "ns", ifelse(pval < 0.05, "Padj < 0.05", "ns"))) %>%
            mutate(size_var = ifelse(is.na(pval), 0.2, -log10(pval))) %>%
            ggplot() +
            geom_jitter(aes(x = switch_rs_score, y = difSmP, color = loop, shape = Significance, size = size_var), width = 0.1) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            labs(x = "Subcompartment switch magnitude", y = "SEN - PRO Counts") +
            mythemecolor1 +
            scale_x_continuous(breaks = seq(-3, 7, by = 1)) # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp1.2.png", 4, 4)

        p <- pf %>%
            mutate(shape_var = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "Padj < 0.05", "NS"))) %>%
            mutate(size_var = ifelse(is.na(pval), 0.2, -log10(pval))) %>%
            ggplot() +
            geom_jitter(aes(x = switch_rs_score, y = difSmP, color = shape_var), width = 0.1) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            labs(x = "Subcompartment switch magnitude", y = "SEN - PRO Counts") +
            mythemecolor +
            scale_x_continuous(breaks = seq(-3, 7, by = 1)) # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp1.2.png", 3, 3)

        p <- summarydf %>%
            ggplot() +
            geom_bar(aes(x = subcomp, y = mean), stat = "identity") +
            geom_errorbar(aes(x = subcomp, ymin = mean - se, ymax = mean + se), width = 0.2) +
            labs(x = "", y = "SEN - PRO Counts") +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp2.png", 2, 4)


        pf <- groupframe %>%
            group_by(gene_id, condition) %>%
            summarise(counts = sum(counts), pval = dplyr::first(!!sym(contrast_padj)), subcomp = dplyr::first(switch_rs), loop = dplyr::first(InSenUniqueLoop)) %>%
            ungroup() %>%
            pivot_wider(names_from = condition, values_from = counts) %>%
            mutate(difSmP = SEN - PRO) %>%
            mutate(group = sprintf("%s | %s", subcomp, loop))


        summarydf <- pf %>%
            group_by(group) %>%
            summarise(
                n = n(),
                mean = mean(difSmP),
                sd = sd(difSmP)
            ) %>%
            mutate(se = sd / sqrt(n)) %>%
            mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1)) %>%
            ungroup() %>%
            mutate(group = paste0(group, " (n = ", n, ")"))


        p <- pf %>%
            ggplot() +
            geom_jitter(aes(x = group, y = difSmP, color = ifelse(is.na(pval), "NS", ifelse(pval < 0.05, "S", "NS"))), width = 0.2) +
            geom_point(data = summarydf, aes(x = group, y = mean), shape = 19) +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp3.png", 8, 4)

        p <- summarydf %>%
            ggplot() +
            geom_errorbar(aes(x = group, ymin = 0, ymax = mean + se), width = 0.2) +
            geom_bar(aes(x = group, y = mean), stat = "identity") +
            labs(x = "", y = "SEN - PRO Counts") +
            theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
            anchorbar +
            mythemecontrastrev
        # mythemecontrastrev +
        # mytheme +
        # guides(fill = "none")
        mysave("temp4.png", 3, 4)



        p <- plots$telescope_multi$condition_SEN_vs_PRO$L1HS$stripp$l1_intactness_req$ALL
        p <- p + theme_bw() + mytheme + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        mysave("temp21.png", 4, 4)

        resultsdf %>% filter(grepl("L1HS", gene_id)) %$% padj_condition_SEN_vs_PRO
    },
    error = function(e) {
        print(e)
        message("Didn't save select plots")
    }
)



x <- data.frame()
write.table(x, file = outputs$outfile, col.names = FALSE)
