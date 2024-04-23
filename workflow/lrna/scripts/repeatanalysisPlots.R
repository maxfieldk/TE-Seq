<<<<<<< HEAD
source("~/data/common/myDefaults.r")
module_name <- "srna"
source("workflow/srna/scripts/generate_colors_to_source.R")
conf <- configr::read.config(file = "conf/config.yaml")[["srna"]]

=======
source("workflow/scripts/defaults.R")
module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/generate_colors_to_source.R")

library(igvR)
>>>>>>> 1b1ec573103095d289cf8d739e3dc9d44a33f71d
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
<<<<<<< HEAD
library(paletteer)
=======
>>>>>>> 1b1ec573103095d289cf8d739e3dc9d44a33f71d

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
            "outputdir" = "results/agg/repeatanalysis/dorado/relaxed",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
        ), env = globalenv())
        assign("inputs", list(
            "resultsdf" = "results/agg/deseq/dorado/relaxed/resultsdf.tsv"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "results/agg/deseq/dorado/relaxed/plots.outfile.txt"
        ), env = globalenv())
    }
)
samples <- conf$samples
sample_table <- read_csv("conf/sample_table.csv")
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

outputdir <- params$outputdir
contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap



## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
resultsdf <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)

# resultsdf <- resultsdf %>% filter(gene_or_te == "repeat") %>% filter(!is.na(rte_family))

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

pvp <- function(df, facet_var = "ALL", filter_var = "ALL") {
    if (filter_var != "ALL") {
        df <- df %>% filter(str_detect(!!sym(filter_var), ">|Intact|towards"))
    }
    pf <- df %>%
        group_by(across(all_of(c(colsToKeep, "condition")))) %>%
        summarise(mean(counts), padj = dplyr::first(!!sym(contrast_padj))) %>%
        pivot_wider(names_from = condition, values_from = `mean(counts)`) 
    top_sig <- pf %>% filter(!!sym(contrast_padj) < 0.05) %>% arrange(padj) %>% head(6) %>% pull(gene_id) 
    p <- pf %>%
        {
            ggplot(data = ., mapping = aes(x = !!sym(contrast_level_1), y = !!sym(contrast_level_2))) +
                geom_point(aes(color = padj < 0.05)) +
                scale_color_manual(values = c("black", "red", "lightgray")) +
                geom_abline(intercept = 0, slope = 1) +
<<<<<<< HEAD
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2), caption = counttype_label) +
                mtclosedgrid +
                theme(aspect.ratio = 1) +
                coord_cartesian(clip = "off") +
                geom_label_repel(data = . %>% 
                        mutate(label = ifelse(gene_id %in% top_sig, gene_id, "")),
                    aes(label = label),
                    size = 2,
                    force_pull   = 0,
                    direction = "both",
                    box.padding = 0.25, max.overlaps = Inf, 
                    hjust = 1,
                    show.legend = FALSE) +
=======
                labs(x = sprintf("%s Norm Counts", contrast_level_1), y = sprintf("%s Norm Counts", contrast_level_2)) +
                mtopen +
>>>>>>> 1b1ec573103095d289cf8d739e3dc9d44a33f71d
                coord_fixed(
                    xlim = range(c(.[[contrast_level_1]], .[[contrast_level_2]])),
                    ylim = range(c(.[contrast_level_1], .[contrast_level_2]))
                )
        }
    if (facet_var != "ALL") {
        p <- p + facet_wrap(facet_var)
    }
    return(p)
}
# df <- tidydf %>% filter(rte_subfamily == "L1HS")
# p <- pvp(tidydf %>% filter(rte_subfamily == "L1HS"), filter_var = "ALL", facet_var = "genic_loc") + ggtitle("L1HS")
# mysave("temp1.png", 6, 6)

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
            mypalette +
            anchorbar +
            mtopen + scale_contrasts +
            guides(fill = "none")
    } else {
        p <- plotframe %>% ggplot() +
            labs(x = "", y = "Number DE", caption = counttype_label) +
            geom_col(aes(x = direction, group = direction, y = count)) +
            facet_wrap(facet_var) +
mtclosed +
            mypalette +
            anchorbar +
            mtopen + scale_contrasts +
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
<<<<<<< HEAD
        labs(x = "", y = "Sum Normalized Counts", caption = counttype_label) +
        extraGGoptions +
        theme(legend.position = "none") +
        mtclosedgridh +
        scale_conditions +
=======
        labs(x = "", y = "Sum Normalized Counts") +
        mtopen + scale_contrasts +
        extraGGoptions +
        theme(legend.position = "none") +
        mtopen +
>>>>>>> 1b1ec573103095d289cf8d739e3dc9d44a33f71d
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

#### TABLES
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
            groupframe <- tidydf %>% dplyr::filter(!!sym(ontology) == !!group)
            nontidygroupframe <- resultsdf %>% dplyr::filter(!!sym(ontology) == !!group)
            # Write table of
            expressed_elements <- groupframe %>%
                filter(counts > 0) %$% gene_id %>%
                unique()
            expressed_elements_df <- nontidygroupframe %>% filter(gene_id %in% expressed_elements)
            dir.create(sprintf("%s/tables/expressed_elements", outputdir), recursive = TRUE, showWarnings = FALSE)
            write_delim(expressed_elements_df, sprintf("%s/tables/expressed_elements/%s.tsv", outputdir, group), delim = "\t")
            write_delim(expressed_elements_df %>% dplyr::select(seqnames, start, end, gene_id, pctdiv, strand), sprintf("%s/tables/expressed_elements/%s.bed", outputdir, group), delim = "\t")

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
                differentially_expressed_elements_df <- expressed_elements_df %>% filter(!!sym(contrast_padj) <= 0.05)
                dir.create(sprintf("%s/tables/differentially_expressed_elements/%s", outputdir, contrast), recursive = TRUE, showWarnings = FALSE)
                write_delim(differentially_expressed_elements_df %>% filter(), sprintf("%s/tables/differentially_expressed_elements/%s/%s.tsv", outputdir, contrast, group), delim = "\t")
                write_delim(differentially_expressed_elements_df %>% dplyr::select(seqnames, start, end, gene_id, pctdiv, strand), sprintf("%s/tables/differentially_expressed_elements/%s/%s.bed", outputdir, contrast, group), delim = "\t", col_names = FALSE)
            }
        }
    }
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


#### PLOTTING
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
                                mysaveandstore(sprintf("%s/%s/%s/%s_%s_%s.png", outputdir, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                
                                if (function_name == "stripp") {
                                    p <- function_current(groupframe, filter_var = filter_var, facet_var = facet_var, stats = "no") + ggtitle(plot_title)
                                    mysaveandstore(sprintf("%s/%s/%s/%s_%s_%s_no_stats.png", outputdir, contrast, function_name, group, filter_var, facet_var), plot_width, plot_height)
                                }
                            },
                            error = function(e) {
                                print(sprintf("Error with  %s %s %s %s %s", contrast, group, function_name, filter_var, facet_var))
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
                                        mysaveandstore(sprintf("%s/%s/%s/%s_%s_%s_%s_%s.png", outputdir, contrast, function_name, group, filter_var, facet_var, DEvar, scaled), plot_width, plot_height)
                                                                            },
                                    error = function(e) {
                                        print(sprintf("Error with  %s %s %s %s %s %s %s %s", contrast, group, function_name, ontology, filter_var, facet_var, DEvar, scaled))
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
                            sigsearched <- resultsdf %>%
                                dplyr::filter((!!sym(ontology)) == group) %>%
                                dplyr::filter((!!sym(contrast_padj)) < 0.05) %>%
                                dplyr::filter((!!sym(contrast_log2FoldChange)) > 0)
                        } else {
                            sigsearched <- resultsdf %>%
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
                        mysaveandstore(sprintf("%s/%s/ggVenn_%s_%s_%s.png", outputdir, modifier, group, modifier, direction), 6, 6)
                        }
                        mysaveandstore(sprintf("%s/ggVenn_%s_%s.png", outputdir, group, direction), 6, 6)
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
