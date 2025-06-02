module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)

set.seed(123)
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
library("emmeans")
library("glmmTMB")


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = sprintf("%s/results/agg/deseq", conf$prefix),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())
        assign("outputs", list(
            stats = sprintf("srna/results/agg/repeatanalysis/telescope_multi/te_group_stats.csv")
        ), env = globalenv())
        assign("inputs", list(
            rte_counts = sprintf("%s/outs/%s/telescope/telescope-run_stats.tsv", conf$prefix, conf$samples),
            sizefactors = "srna/results/agg/deseq/telescope_multi/sizefactors.csv"
        ), env = globalenv())
    }
)



counttype <- params[["counttype"]]

rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE,
    append_NI_samplename_modifier = if (conf$per_sample_ref == "yes") TRUE else if (conf$per_sample_ref == "yes") TRUE else FALSE
)

if (counttype == "telescope_multi") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- gsub("run_stats", "TE_counts", grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE))
        tdf <- read_delim(path, comment = "#", col_names = TRUE) %>%
            dplyr::rename(gene_id = transcript)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

if (counttype == "telescope_unique") {
    bounddf <- tibble(gene_id = as.character())
    for (sample in conf$samples) {
        path <- grep(paste0("outs/", sample, "/telescope"), inputs$rte_counts, value = TRUE)
        tdf <- read_delim(path, comment = "#", col_names = FALSE) %>%
            dplyr::select(X1, X6) %>%
            dplyr::rename(gene_id = X1)
        bounddf <- full_join(bounddf, tdf, by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", conf$samples)
}

bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]
cts <- as.data.frame(bounddf1 %>% replace(is.na(.), 0)) %>%
    tibble() %>%
    mutate(across(-gene_id, ~ as.integer(round(.))))
rm(bounddf1)

tidydf <- cts %>%
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "counts") %>%
    {
        if (conf$per_sample_ref == "yes") {
            mutate(.,
                gene_id = case_when(
                    grepl("_NI_", gene_id) ~ paste0(sample_name, "__", gene_id),
                    TRUE ~ gene_id
                )
            ) %>%
                left_join(rmann %>% dplyr::rename(nonrefinsert_sample_name = sample_name))
        } else {
            left_join(., rmann)
        }
    } %>%
    filter(!grepl("__AS$", gene_id))

size_factors <- read_csv(inputs$sizefactors)
sample_table
#### New stats freq

strat_1_options <- c("none", "rte_length_req", "req_integrative")
strat_2_options <- c("none", "genic_loc", "loc_lowres_integrative_stranded")
options <- expand_grid(strat_1_options, strat_2_options)


if (is.null(conf$rte_subfamilies_for_aggregate_rte_stats) | conf$rte_subfamilies_for_aggregate_rte_stats == "") {
    subfams_for_stat <- rmann %$% rte_subfamily %>%
        unique() %>%
        grep("Other", ., invert = TRUE, value = TRUE)
} else {
    subfams_for_stat <- conf$rte_subfamilies_for_aggregate_rte_stats
}


# ensure batch variables used in linear model have more than one level!
batch_vars_to_use <- c()
if (any(grepl("^batch", colnames(sample_table)) | grepl("^covariate", colnames(sample_table)))) {
    for (value in colnames(sample_table)[grepl("^batch", colnames(sample_table)) | grepl("^covariate", colnames(sample_table))]) {
        number_unique_vals <- sample_table %>%
            pull(value) %>%
            unique() %>%
            length()
        if (number_unique_vals > 1) {
            batch_vars_to_use <- c(batch_vars_to_use, value)
        }
    }
}


results_list <- list()
for (i in 1:nrow(options)) {
    strat_1 <- options[[i, "strat_1_options"]]
    strat_2 <- options[[i, "strat_2_options"]]
    strat_vars <- c(strat_1, strat_2) %>% grep("none", ., invert = TRUE, value = TRUE)

    pf1 <- tidydf %>%
        filter(rte_subfamily %in% subfams_for_stat) %>%
        group_by(rte_subfamily, sample_name) %>%
        group_by(across(all_of(strat_vars)), .add = TRUE) %>%
        summarise(counts = sum(counts)) %>%
        ungroup() %>%
        left_join(sample_table) %>%
        left_join(size_factors)

    if (is.null(batch_vars_to_use)) {
        models_df <- pf1 %>%
            group_by(rte_subfamily) %>%
            group_by(across(all_of(strat_vars)), .add = TRUE) %>%
            nest() %>%
            mutate(
                model = map(data, ~ glmmTMB(counts ~ condition + offset(log(sizefactor)), family = nbinom2, data = .x))
            )
    } else {
        models_df <- pf1 %>%
            group_by(rte_subfamily) %>%
            group_by(across(all_of(strat_vars)), .add = TRUE) %>%
            nest() %>%
            mutate(
                model = map(data, ~ glmmTMB(formula(sprintf("counts ~ condition + %s + offset(log(sizefactor))", paste0(batch_vars_to_use, collapse = " + "))), family = nbinom2, data = .x))
            )
    }

    desired_contrasts <- conf$contrasts %>%
        gsub("^condition_", "", .) %>%
        strsplit("_vs_") %>%
        lapply(function(x) tibble(level1 = x[2], level2 = x[1])) %>%
        bind_rows() %>%
        mutate(contrast_string = paste(level1, "-", level2))

    results_df <- models_df %>%
        mutate(
            emmeans_obj = map(model, ~ emmeans(.x, ~condition)),
            contrast_df = map(emmeans_obj, ~ contrast(.x, method = "pairwise") %>% as.data.frame())
        ) %>%
        select(all_of(strat_vars), contrast_df) %>%
        unnest(contrast_df) %>%
        filter(contrast %in% desired_contrasts$contrast_string) %>%
        mutate(stratification_variables = paste(strat_vars, collapse = ";")) %>%
        ungroup() %>%
        mutate(p_adj = p.adjust(p.value, method = "fdr"))

    results_list[[paste(strat_vars, collapse = ";")]] <- results_df
}

resdf <- purrr::reduce(results_list, bind_rows) %>%
    mutate(p_adj_across_all_stratification_types = p.adjust(p.value, method = "fdr"))

dir.create(dirname(outputs$stats), recursive = TRUE)
write_csv(resdf, outputs$stats)

######
