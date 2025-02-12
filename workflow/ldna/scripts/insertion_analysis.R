module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

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
library(ggh4x)

samples <- conf$samples
source("conf/sample_table_source.R")


tryCatch(
    {
        params <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            rmann_nonref = sprintf("aref/extended/%s_annotations/%s_rmann_nonref.csv", sample_table$sample_name, sample_table$sample_name),
            filtered_tldr = sprintf("aref/extended/%s.table.kept_in_updated_ref.txt", sample_table$sample_name)
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
    df <- read.table(grep(sprintf("%s.table", sample), params$filtered_tldr, value = TRUE), header = TRUE)
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
    df1 <- read_csv(grep(sprintf("%s_annotations", sample), params$rmann, value = TRUE))
    df1$sample_name <- sample
    df <- df1 %>%
        left_join(sample_table)
    anns[[sample]] <- df
}
ann <- do.call(rbind, anns) %>% tibble()

nrdf <- ann %>%
    filter(refstatus == "NonRef") %>%
    dplyr::rename(UUID = nonref_UUID) %>%
    left_join(dff) %>%
    mutate(zygosity = factor(ifelse(fraction_reads_count >= 0.95, "homozygous", "heterozygous"), levels = c("homozygous", "heterozygous"))) %>%
    mutate(
        known_nonref = factor(
            case_when(
                confALL$aref$update_ref_with_tldr$known_nonref$response == "no" ~ "NotCalled",
                is.na(NonRef) ~ "novel",
                TRUE ~ "known"
            ),
            levels = c("novel", "known", "NotCalled")
        )
    )

# removing the rmann seqnames, start, stop, leaving only the tldr (ref genome coordinate based ones)
nrdfgrs <- GRanges(nrdf %>% dplyr::rename(rmann_seqnames = seqnames) %>% dplyr::select(-start, -end, -strand))
reduced_ranges <- reduce(nrdfgrs)
overlaps <- findOverlaps(nrdfgrs, reduced_ranges)
group_id <- paste0("GID", as.factor(subjectHits(overlaps)))
nrdfgrs$group_id <- group_id
nrdf_new <- tibble(as.data.frame(nrdfgrs)) %>%
    group_by(group_id) %>%
    mutate(num_samples_detected = n()) %>%
    mutate(shared = ifelse(num_samples_detected > 1, "Shared", "NotShared")) %>%
    relocate(num_samples_detected, shared) %>%
    ungroup()

p1 <- nrdf %>%
    group_by(rte_family, loc_integrative, sample_name) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "sample_name", y = "counts", fill = "loc_integrative") +
    facet_grid2(cols = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_subfamily_loc_integrative.pdf", outputdir), 12, 4, 9, 4, plus_void = TRUE)

p1 <- nrdf %>%
    group_by(rte_family, req_integrative, sample_name) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "sample_name", y = "counts", fill = "req_integrative") +
    facet_grid2(cols = vars(rte_family), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_subfamily_req_integrative.pdf", outputdir), 12, 4, 9, 4, , plus_void = TRUE)

p1 <- nrdf %>%
    filter(rte_family == "L1") %>%
    group_by(rte_subfamily, req_integrative, sample_name) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "sample_name", y = "counts", fill = "req_integrative") +
    facet_grid2(cols = vars(rte_subfamily), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_l1_req_integrative.pdf", outputdir), 8, 4, 6, 4, plus_void = TRUE)

p1 <- nrdf %>%
    filter(rte_family == "L1") %>%
    group_by(rte_subfamily, req_integrative, sample_name) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "sample_name", y = "counts", fill = "req_integrative") +
    facet_grid2(cols = vars(rte_subfamily), scales = "free", space = "free_x", independent = "y") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_l1_req_integrative1.pdf", outputdir), 11, 4, 6, 4, plus_void = TRUE)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, known_nonref) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "known_nonref") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_subfamily_nofacet_known.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, zygosity) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "zygosity") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_subfamily_nofacet_zygosity.pdf", outputdir), 6, 4)

p1 <- nrdf %>%
    group_by(rte_family, rte_subfamily, req_integrative) %>%
    summarise(count = n()) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "rte_subfamily", y = "counts", fill = "req_integrative") +
    ggtitle("Non-reference RTE Insertions") +
    labs(x = "Subfamily", y = "Counts") +
    scale_palette +
    mtclosed + anchorbar + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mss(pl = p1, sprintf("%s/insertions_subfamily_nofacet.pdf", outputdir), 6, 4, plus_void = TRUE)


novel_frac_df <- nrdf %>%
    group_by(rte_subfamily, known_nonref, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(known_nonref == "novel") %>%
    mutate(frac_novel = frac_total, frac_known = 1 - frac_total)

shared_frac_df <- nrdf_new %>%
    group_by(rte_subfamily, shared, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(shared == "Shared") %>%
    mutate(frac_shared = frac_total, frac_unique = 1 - frac_total)

homozyg_frac_df <- nrdf %>%
    group_by(rte_subfamily, zygosity, .drop = FALSE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(rte_subfamily) %>%
    mutate(family_n = sum(n)) %>%
    mutate(frac_total = n / family_n) %>%
    ungroup() %>%
    dplyr::filter(zygosity == "heterozygous") %>%
    mutate(frac_heterozygous = frac_total, frac_homozygous = 1 - frac_total)

hp <- homozyg_frac_df %>%
    pivot_longer(cols = c(frac_heterozygous, frac_homozygous)) %>%
    ggbarplot(x = "rte_subfamily", y = "value", fill = "name") +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("tan", "lightblue")) +
    labs(x = "Subfamily")
sp <- shared_frac_df %>%
    pivot_longer(cols = c(frac_shared, frac_unique)) %>%
    ggbarplot(x = "rte_subfamily", y = "value", fill = "name") +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("#018458", "#ffc354")) +
    labs(x = "Subfamily")

np <- novel_frac_df %>%
    pivot_longer(cols = c(frac_novel, frac_known)) %>%
    ggbarplot(x = "rte_subfamily", y = "value", fill = "name") +
    mtclosed + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("brown", "grey")) +
    labs(x = "Subfamily")
library(patchwork)

mss(pl = hp, sprintf("%s/insertions_subfamily_zygosity.pdf", outputdir), 6, 2, 3, 2.25, plus_void = TRUE)
mss(pl = sp, sprintf("%s/insertions_subfamily_sharedprop.pdf", outputdir), 6, 2, 3, 2.25, plus_void = TRUE)
mss(pl = np, sprintf("%s/insertions_subfamily_novelprop.pdf", outputdir), 6, 2, 3, 2.25, plus_void = TRUE)
ptch <- np / hp + plot_layout(heights = c(0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch_only_annot.pdf", outputdir), 6, 4, )
ptch <- sp / hp + plot_layout(heights = c(0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch_only_annot_shared.pdf", outputdir), 6, 4, )


ptch <- p1 / np / hp + plot_layout(heights = c(0.4, 0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch.pdf", outputdir), 6, 6)

ptch <- p1 / sp / hp + plot_layout(heights = c(0.4, 0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch_shared.pdf", outputdir), 6, 6)

ptch <- p1 / np / sp / hp + plot_layout(heights = c(0.4, 0.4, 0.4, 0.4), axis_titles = "collect")
mysaveandstore(pl = ptch, sprintf("%s/insertions_subfamily_patch_all.pdf", outputdir), 6, 8)


for (rte in homozyg_frac_df %$% rte_subfamily) {
    pf <- homozyg_frac_df %>%
        filter(rte_subfamily == rte) %>%
        pivot_longer(cols = c(frac_heterozygous, frac_homozygous), names_to = "label")

    p <- ggpie(pf, "value", "label", fill = "label") +
        scale_fill_manual(values = c("yellow", "blue"))
    mysaveandstore(sprintf("%s/%s_zygosity_pie.pdf", outputdir, rte), 3, 3)
}



p2 <- nrdf %>%
    group_by(loc_integrative) %>%
    summarise(count = n()) %>%
    mutate(loc_integrative = fct_reorder(loc_integrative, count)) %>%
    mutate(counts = count) %>%
    ggbarplot(x = "loc_integrative", y = "counts", fill = "loc_integrative") +
    ggtitle("Non-reference RTE Insertions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Genomic Context", y = "Count") +
    scale_palette +
    coord_flip() +
    mtclosed + anchorbar
mysaveandstore(pl = p2, sprintf("%s/insertions_genomic_context.pdf", outputdir), 6, 4)
library(patchwork)
p <- p1 + p2 + plot_layout(widths = c(1, 1), guides = "collect")
mysaveandstore(sprintf("%s/insertions_sub_fig.pdf", outputdir), 12, 5)


p <- nrdf %>%
    filter(rte_subfamily == "L1HS") %>%
    ggplot() +
    geom_histogram(aes(x = LengthIns)) +
    geom_vline(xintercept = 6000, color = "red", linetype = 2) +
    ggtitle("L1HS Insertion Lengths") +
    labs(x = "Length (bp)", y = "Count") +
    mtopen +
    scale_palette +
    anchorbar
mysaveandstore(sprintf("%s/l1hs_length_in_updated_ref.pdf", outputdir), 5, 5)


x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)








x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)
