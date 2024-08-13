module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(circlize)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(Biostrings)
library(ggpubr)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]


tryCatch(
    {
        params <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            r_annotation_fragmentsjoined = sprintf("aref/%s_annotations/%s_repeatmasker.gtf.rformatted.fragmentsjoined.csv", sample_table$sample_name, sample_table$sample_name),
            r_repeatmasker_annotation = sprintf("aref/%s_annotations/%s_repeatmasker_annotation.csv", sample_table$sample_name, sample_table$sample_name),
            filtered_tldr = sprintf("aref/%s_tldr/%s.table.kept_in_updated_ref.txt", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(
            plots = "ldna/results/insertion_analysis/insertion_analysis.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

dfs_filtered <- list()
for (sample in sample_table$sample_name) {
    df <- read.table(grep(sprintf("%s_tldr", sample), params$filtered_tldr, value = TRUE), header = TRUE)
    df$sample_name <- sample
    df <- df %>% left_join(sample_table)
    dfs_filtered[[sample]] <- df
}

dff <- do.call(rbind, dfs_filtered) %>% tibble()


for (element_type in dff %$% Subfamily %>% unique()) {
    dftemp <- dff %>% filter(Subfamily == element_type)
    p <- dftemp %>%
        gghistogram(x = "LengthIns") +
        mtopen +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/length_distribution.pdf", outputdir, element_type))
    dftemp %$% LengthIns %>% quantile()
    p <- dftemp %>%
        mutate(tsdlen = nchar(TSD)) %>%
        gghistogram(x = "tsdlen", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/tsd_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd5 = nchar(Transduction_5p)) %>%
        gghistogram(x = "trsd5", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/trsd5_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd3 = nchar(Transduction_3p)) %>%
        gghistogram(x = "trsd3", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/trsd3_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        arrange(-EndTE) %>%
        mutate(nrow = row_number()) %>%
        ggplot() +
        geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
        mtclosed +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/te_body_distribution.pdf", outputdir, element_type))
}




anns <- list()
for (sample in sample_table$sample_name) {
    df1 <- read_csv(grep(sprintf("%s_annotations", sample), params$r_annotation_fragmentsjoined, value = TRUE))
    df1$sample_name <- sample
    df2 <- read_csv(grep(sprintf("%s_annotations", sample), params$r_repeatmasker_annotation, value = TRUE))
    df <- df1 %>%
        filter(refstatus == "NonRef") %>%
        left_join(sample_table) %>%
        left_join(df2)
    anns[[sample]] <- df
}
ann <- do.call(rbind, anns) %>% tibble()



p <- ann %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_subfamily, scales = "free", space = "free_x", ncol = 2) +
    geom_bar(aes(fill = rte_length_req)) +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(w = 8, h = 6)

p <- ann %>%
    ggplot(aes(x = sample_name)) +
    facet_grid(~rte_subfamily, scales = "free", space = "free_x", ncol = 2) +
    geom_bar(aes(fill = loc_integrative_loc)) +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions") +
    anchorbar +
    mtclosedgridh +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(w = 8, h = 6)



p <- df %>%
    filter(Family != "NA") %>%
    ggplot(aes(x = Subfamily)) +
    facet_grid(~Family, scales = "free", space = "free_x") +
    geom_bar(color = "black") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Count") +
    paletteer::scale_fill_paletteer_d(conf$default_palette) +
    ggtitle("Non-reference RTE Insertions")
mysaveandstore(sprintf("%s/insertions_subfamily.pdf", outputdir), 5, 5)

save(mysaveandstoreplots, file = outputs$plots)
