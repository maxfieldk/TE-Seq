library(paletteer)
library(colorspace)

sample_table <- read_csv(conf[["sample_table"]])
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]

# these are optionally used depending on values set in config.yaml
custom_descriptive <- c("#4b9c7a", "#ca6628", "#6d68b5", "#d43f88", "#0068f9", "#ffe100", "#7B3A96FF", "#CB2027FF", "#808180FF", "#1B1919FF")
custom_condition <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF")

tryCatch(
    {
        condition_palette_length <<- length(as.character(paletteer_d(confALL$shared$default_palette_condition)))
    },
    error = function(e) {
        condition_palette_length <<- length(custom_condition)
    }
)
# build color palettes:
if (length(conf$levels) == 1) {
    manual_color_vec <- confALL$shared$condition_colors %>% unlist()
    if (all(conf$levels %in% names(manual_color_vec))) {
        condition_palette <- manual_color_vec
    } else {
        condition_palette <- setNames(c("grey"), conf$levels)
    }

    color_table <- sample_table %>%
        group_by(condition) %>%
        mutate(replicate = row_number()) %>%
        ungroup() %>%
        left_join(tibble(condition = names(condition_palette), color = condition_palette))
    num_replicates <- color_table %$% replicate %>%
        unique() %>%
        length()
    if (nrow(color_table) == 1) {
        color_table <- color_table %>%
            mutate(sample_unique_color = color)
    } else {
        color_table <- color_table %>%
            left_join(tibble(replicate = 1:num_replicates, shade_modifier = seq(-.4, 0.4, length.out = num_replicates))) %>%
            mutate(sample_unique_color = darken(color, shade_modifier))
    }

    sample_palette <- setNames(color_table$color, color_table$sample_name)
    sample_unique_palette <- setNames(color_table$sample_unique_color, color_table$sample_name)
    direction_palette <- setNames(c("red", "blue", "grey"), c("UP", "DOWN", "NS"))
    methylation_palette <- setNames(c("#ede6d1", "#2b2a24", "grey"), c("Hypo", "Hyper", "NS"))
    methylation_palette_thresholds <- setNames(c("#ede6d1", "#4a493f", "#e8d18c", "#2b2a24", "grey"), c("Hypo_05", "Hyper_05", "Hypo_01", "Hyper_01", "NS"))

    contrasts <- conf$contrasts
    contrast_base <- gsub(".*_vs_", "", contrasts)
    contrast_level <- gsub("condition_", "", gsub("_vs_.*", "", contrasts))
    contrasts_with_same_base <- contrast_base == conf$levels[1]
    contrasts_without_same_base <- !contrasts_with_same_base
    tojoin <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_with_same_base, ]
    toappend <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_without_same_base, ]

    tryCatch(
        {
            contrastframe <- color_table %>%
                left_join(tojoin) %>%
                filter(!is.na(contrast)) %>%
                group_by(contrast) %>%
                summarise(contrast_color = dplyr::first(color)) %>%
                dplyr::select(contrast, contrast_color)
            contrast_palette <- setNames(contrastframe$contrast_color, contrastframe$contrast)
            if (length(rownames(toappend)) > 0) {
                palette_index_start <- length(conf$levels)
                contrasts_to_append <- toappend$contrast
                tryCatch(
                    {
                        to_append_palette <<- setNames(as.character(paletteer_d(confALL$shared$default_palette_condition)[palette_index_start:palette_index_start + length(contrasts_to_append)]), contrasts_to_append)
                    },
                    error = function(e) {
                        to_append_palette <<- setNames(as.character(paletteer_c("viridis::cividis", direction = 1, n = length(contrasts_to_append))), contrasts_to_append)
                    }
                )
                contrast_palette <- c(contrast_palette, to_append_palette)
            }
            scale_contrasts <- list(scale_fill_manual(values = contrast_palette), scale_color_manual(values = contrast_palette))
        },
        error = function(e) {

        }
    )

    {
        tryCatch(
            {
                paletteer_d(confALL$shared$default_palette_descriptive)
                scale_palette <<- list(paletteer::scale_fill_paletteer_d(confALL$shared$default_palette_descriptive), paletteer::scale_color_paletteer_d(confALL$shared$default_palette_descriptive))
            },
            error = function(e) {
                scale_palette <<- list(scale_fill_manual(values = custom_descriptive), scale_color_manual(values = custom_descriptive))
            }
        )
        scale_palette_alt <- list(paletteer::scale_fill_paletteer_d("ggthemes::few_Dark"), paletteer::scale_color_paletteer_d("ggthemes::few_Dark"))

        scale_samples_unique <- list(scale_fill_manual(values = sample_unique_palette), scale_color_manual(values = sample_unique_palette))
        scale_samples <- list(scale_fill_manual(values = sample_palette), scale_color_manual(values = sample_palette))
        scale_conditions <- list(scale_fill_manual(values = condition_palette), scale_color_manual(values = condition_palette))
        scale_directions <- list(scale_fill_manual(values = direction_palette), scale_color_manual(values = direction_palette))
        scale_methylation <- list(scale_fill_manual(values = methylation_palette), scale_color_manual(values = methylation_palette))
        scale_methylation_thresholds <- list(scale_fill_manual(values = methylation_palette_thresholds), scale_color_manual(values = methylation_palette_thresholds))
    }
    mycolor <- "grey"
} else if (length(conf$levels) <= 1 + condition_palette_length) {
    manual_color_vec <- confALL$shared$condition_colors %>% unlist()
    if (all(conf$levels %in% names(manual_color_vec))) {
        condition_palette <- manual_color_vec
    } else {
        tryCatch(
            {
                condition_palette <<- setNames(c("grey", as.character(paletteer_d(confALL$shared$default_palette_condition, direction = 1, n = length(conf$levels) - 1))), conf$levels)
            },
            error = function(e) {
                condition_palette <<- setNames(c("grey", custom_condition[1:length(conf$levels) - 1]), conf$levels)
            }
        )
    }
    color_table <- sample_table %>%
        group_by(condition) %>%
        mutate(replicate = row_number()) %>%
        ungroup() %>%
        left_join(tibble(condition = names(condition_palette), color = condition_palette))
    num_replicates <- color_table %$% replicate %>%
        unique() %>%
        length()
    color_table <- color_table %>%
        left_join(tibble(replicate = 1:num_replicates, shade_modifier = seq(-.4, 0.4, length.out = num_replicates))) %>%
        mutate(sample_unique_color = darken(color, shade_modifier))

    sample_palette <- setNames(color_table$color, color_table$sample_name)
    sample_unique_palette <- setNames(color_table$sample_unique_color, color_table$sample_name)
    direction_palette <- setNames(c("red", "blue", "grey"), c("UP", "DOWN", "NS"))
    methylation_palette <- setNames(c("#ede6d1", "#2b2a24", "grey"), c("Hypo", "Hyper", "NS"))
    methylation_palette_thresholds <- setNames(c("#ede6d1", "#4a493f", "#e8d18c", "#2b2a24", "grey"), c("Hypo_05", "Hyper_05", "Hypo_01", "Hyper_01", "NS"))

    contrasts <- conf$contrasts
    contrast_base <- gsub(".*_vs_", "", contrasts)
    contrast_level <- gsub("condition_", "", gsub("_vs_.*", "", contrasts))
    contrasts_with_same_base <- contrast_base == conf$levels[1]
    contrasts_without_same_base <- !contrasts_with_same_base
    tojoin <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_with_same_base, ]
    toappend <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_without_same_base, ]

    tryCatch(
        {
            contrastframe <- color_table %>%
                left_join(tojoin) %>%
                filter(!is.na(contrast)) %>%
                group_by(contrast) %>%
                summarise(contrast_color = dplyr::first(color)) %>%
                dplyr::select(contrast, contrast_color)
            contrast_palette <- setNames(contrastframe$contrast_color, contrastframe$contrast)
            if (length(rownames(toappend)) > 0) {
                palette_index_start <- length(conf$levels)
                contrasts_to_append <- toappend$contrast
                tryCatch(
                    {
                        to_append_palette <<- setNames(as.character(paletteer_d(confALL$shared$default_palette_condition)[palette_index_start:palette_index_start + length(contrasts_to_append)]), contrasts_to_append)
                    },
                    error = function(e) {
                        to_append_palette <<- setNames(as.character(paletteer_c("viridis::cividis", direction = 1, n = length(contrasts_to_append))), contrasts_to_append)
                    }
                )
                contrast_palette <- c(contrast_palette, to_append_palette)
            }
            scale_contrasts <- list(scale_fill_manual(values = contrast_palette), scale_color_manual(values = contrast_palette))
        },
        error = function(e) {

        }
    )

    {
        tryCatch(
            {
                paletteer_d(confALL$shared$default_palette_descriptive)
                scale_palette <<- list(paletteer::scale_fill_paletteer_d(confALL$shared$default_palette_descriptive), paletteer::scale_color_paletteer_d(confALL$shared$default_palette_descriptive))
            },
            error = function(e) {
                scale_palette <<- list(scale_fill_manual(values = custom_descriptive), scale_color_manual(values = custom_descriptive))
            }
        )
        scale_palette_alt <- list(paletteer::scale_fill_paletteer_d("ggthemes::few_Dark"), paletteer::scale_color_paletteer_d("ggthemes::few_Dark"))

        scale_samples_unique <- list(scale_fill_manual(values = sample_unique_palette), scale_color_manual(values = sample_unique_palette))
        scale_samples <- list(scale_fill_manual(values = sample_palette), scale_color_manual(values = sample_palette))
        scale_conditions <- list(scale_fill_manual(values = condition_palette), scale_color_manual(values = condition_palette))
        scale_directions <- list(scale_fill_manual(values = direction_palette), scale_color_manual(values = direction_palette))
        scale_methylation <- list(scale_fill_manual(values = methylation_palette), scale_color_manual(values = methylation_palette))
        scale_methylation_thresholds <- list(scale_fill_manual(values = methylation_palette_thresholds), scale_color_manual(values = methylation_palette_thresholds))
    }
    mycolor <- "grey"
} else if (length(conf$levels) > 1 + condition_palette_length) {
    manual_color_vec <- confALL$shared$condition_colors %>% unlist()
    if (all(conf$levels %in% names(manual_color_vec))) {
        condition_palette <- manual_color_vec
    } else {
        tryCatch(
            {
                condition_palette <<- setNames(c("grey", custom_condition[1:length(conf$levels) - 1]), conf$levels)
            },
            error = function(e) {
                condition_palette <<- setNames(c("grey", as.character(paletteer_c("viridis::inferno", direction = 1, n = length(conf$levels) - 1))), conf$levels)
            }
        )
    }

    color_table <- sample_table %>%
        group_by(condition) %>%
        mutate(replicate = row_number()) %>%
        ungroup() %>%
        left_join(tibble(condition = names(condition_palette), color = condition_palette))
    num_replicates <- color_table %$% replicate %>%
        unique() %>%
        length()
    color_table <- color_table %>%
        left_join(tibble(replicate = 1:num_replicates, shade_modifier = seq(-.4, 0.4, length.out = num_replicates))) %>%
        mutate(sample_unique_color = darken(color, shade_modifier))

    sample_palette <- setNames(color_table$color, color_table$sample_name)
    sample_unique_palette <- setNames(color_table$sample_unique_color, color_table$sample_name)
    direction_palette <- setNames(c("red", "blue", "grey"), c("UP", "DOWN", "NS"))
    methylation_palette <- setNames(c("#ede6d1", "#2b2a24", "grey"), c("Hypo", "Hyper", "NS"))
    methylation_palette_thresholds <- setNames(c("#ede6d1", "#4a493f", "#e8d18c", "#2b2a24", "grey"), c("Hypo_05", "Hyper_05", "Hypo_01", "Hyper_01", "NS"))

    contrasts <- conf$contrasts
    contrast_base <- gsub(".*_vs_", "", contrasts)
    contrast_level <- gsub("condition_", "", gsub("_vs_.*", "", contrasts))
    contrasts_with_same_base <- contrast_base == conf$levels[1]
    contrasts_without_same_base <- !contrasts_with_same_base
    tojoin <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_with_same_base, ]
    toappend <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_without_same_base, ]

    tryCatch(
        {
            contrastframe <- color_table %>%
                left_join(tojoin) %>%
                filter(!is.na(contrast)) %>%
                group_by(contrast) %>%
                summarise(contrast_color = dplyr::first(color)) %>%
                dplyr::select(contrast, contrast_color)
            contrast_palette <- setNames(contrastframe$contrast_color, contrastframe$contrast)
            if (length(rownames(toappend)) > 0) {
                palette_index_start <- length(conf$levels)
                contrasts_to_append <- toappend$contrast
                to_append_palette <- setNames(as.character(paletteer_c("viridis::cividis", direction = 1, n = length(contrasts_to_append))), contrasts_to_append)
                contrast_palette <- c(contrast_palette, to_append_palette)
            }
            scale_contrasts <- list(scale_fill_manual(values = contrast_palette), scale_color_manual(values = contrast_palette))
        },
        error = function(e) {

        }
    )

    {
        tryCatch(
            {
                paletteer_d(confALL$shared$default_palette_descriptive)
                scale_palette <<- list(paletteer::scale_fill_paletteer_d(confALL$shared$default_palette_descriptive), paletteer::scale_color_paletteer_d(confALL$shared$default_palette_descriptive))
            },
            error = function(e) {
                scale_palette <<- list(scale_fill_manual(values = custom_descriptive), scale_color_manual(values = custom_descriptive))
            }
        )
        scale_palette_alt <- list(paletteer::scale_fill_paletteer_d("ggthemes::few_Dark"), paletteer::scale_color_paletteer_d("ggthemes::few_Dark"))
        scale_samples_unique <- list(scale_fill_manual(values = sample_unique_palette), scale_color_manual(values = sample_unique_palette))
        scale_samples <- list(scale_fill_manual(values = sample_palette), scale_color_manual(values = sample_palette))
        scale_conditions <- list(scale_fill_manual(values = condition_palette), scale_color_manual(values = condition_palette))
        scale_directions <- list(scale_fill_manual(values = direction_palette), scale_color_manual(values = direction_palette))
        scale_methylation <- list(scale_fill_manual(values = methylation_palette), scale_color_manual(values = methylation_palette))
        scale_methylation_thresholds <- list(scale_fill_manual(values = methylation_palette_thresholds), scale_color_manual(values = methylation_palette_thresholds))
    }
    mycolor <- "grey"
}

if (file.exists("conf/generate_colors_to_source_override.R")) {
    source("conf/generate_colors_to_source_override.R")
}
