library(paletteer)
library(colorspace)

{
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]

sample_table <- read_csv(conf[["sample_table"]])
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]
#build color palettes:
condition_palette <- setNames(c("darkgrey", as.character(paletteer_d(conf$default_palette, direction = 1, n = length(conf$levels)-1))), conf$levels)

color_table <- sample_table %>% group_by(condition) %>% mutate(replicate = row_number()) %>% ungroup() %>%
    left_join(tibble(condition = names(condition_palette), color = condition_palette))
num_replicates <- color_table %$% replicate %>% unique() %>% length()
color_table <- color_table %>% left_join(tibble(replicate = 1:num_replicates, shade_modifier = seq(-.4, 0.4, length.out=num_replicates))) %>%
    mutate(sample_unique_color = darken(color, shade_modifier))

sample_palette <- setNames(color_table$color, color_table$sample_name)
sample_unique_palette <- setNames(color_table$sample_unique_color, color_table$sample_unique_color)
direction_palette <- setNames(c("red", "blue", "grey"), c("UP", "DOWN", "NS"))
methylation_palette <- setNames(c("red", "blue", "grey"), c("Hypo", "Hyper", "NS"))

contrasts <- conf$contrasts
contrast_base <- gsub(".*_vs_", "", contrasts)
contrast_level <- gsub("condition_", "", gsub("_vs_.*", "", contrasts))
contrasts_with_same_base <- contrast_base == conf$levels[1]
contrasts_without_same_base <- !contrasts_with_same_base
tojoin <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_with_same_base,]
toappend <- tibble(contrast = contrasts, base = contrast_base, condition = contrast_level)[contrasts_without_same_base,]
contrastframe <- color_table %>% left_join(tojoin) %>% filter(!is.na(contrast)) %>% group_by(contrast) %>% summarise(contrast_color = dplyr::first(color)) %>% dplyr::select(contrast, contrast_color)
contrast_palette <- setNames(contrastframe$contrast_color, contrastframe$contrast)
if (length(rownames(toappend)) > 0) {
    palette_index_start <- length(conf$levels)
    contrasts_to_append <- toappend$contrast
    to_append_palette <- setNames(as.character(paletteer_d(conf$default_palette)[palette_index_start:palette_index_start+length(contrasts_to_append)]), contrasts_to_append)
    contrast_palette <- c(contrast_palette, to_append_palette)
}


scale_samples_unique <- list(scale_fill_manual(values = sample_unique_palette), scale_color_manual(values = sample_unique_palette))
scale_samples <- list(scale_fill_manual(values = sample_palette), scale_color_manual(values = sample_palette))
scale_conditions <- list(scale_fill_manual(values = condition_palette), scale_color_manual(values = condition_palette))
scale_directions <- list(scale_fill_manual(values = direction_palette), scale_color_manual(values = direction_palette))
scale_contrasts <- list(scale_fill_manual(values = contrast_palette), scale_color_manual(values = contrast_palette))
}


save(condition_palette, contrast_palette, sample_palette,sample_unique_palette, direction_palette,
    scale_samples_unique, scale_samples, scale_directions, scale_contrasts, scale_conditions,
    file = sprintf("conf/colors_%s.RData", module_name))
