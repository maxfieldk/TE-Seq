module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(zoo)
library(pryr)
library(circlize)
library(rGREAT)
library(reactome.db)
library(msigdb)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)
# library(ReMapEnrich)
library(msigdbr)
library(Biostrings)
library(ggpubr)



sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]
conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")
    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
    if ("chrY" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, grep("_chrX_|_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes)
        } else {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, grep("_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX)
        }
    } else if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrY, grep("_chrX_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrY)
    } else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
    }
}
#################### functions and themes

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", sample_table$sample_name, sample_table$sample_name),
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv"
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/outfiles/bedmethylanalysis.txt",
            promoters_bed = "ldna/Rintermediates/promoters.bed",
            dmrpromoterhyper_bed = "ldna/Rintermediates/promoters_dmhyperregions.bed",
            dmrpromoterhypo_bed = "ldna/Rintermediates/promoters_dmhyporegions.bed"
        ), env = globalenv())
    }
)

r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)

rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)
flRTEpromoter <- read_delim("ldna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
perelementdf_promoters <- read_delim("ldna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
reads <- read_delim("ldna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)





library(patchwork)

aa <- left_join(reads %>% head(), rmann, by = c("gene_id"))

aa %>% dplyr::select(start.x, end.x, start.y, end.y)

anatomy <- read_delim("aref/A.REF_Analysis/intact_l1_anatomy_coordinates.tsv")


l1hs <- rtedf %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(intactness_req == "Intact")
l1hs %$% gene_id %>% table()




l1hsP <- l1hs %>%
    filter(rte_strand == "+") %>%
    mutate(rel_cpg_pos = (start + element_start) - rte_start)
l1hsN <- l1hs %>%
    filter(rte_strand == "-") %>%
    mutate(rel_cpg_pos = (rte_end + element_start) - start)
l1hs <- bind_rows(l1hsP, l1hsN)


p <- l1hs %>%
    ggplot(aes(x = rel_cpg_pos, y = pctM)) +
    geom_point(alpha = 0.1)
mysaveandstore()

outputdir <- "ldna/results/plots/l1hs_meth_patterns"

color_intervals <- anatomy %>%
    filter(!(feature %in% c("EN", "RT"))) %>%
    group_by(feature) %>%
    summarise(start = mean(start), end = mean(end))

p1 <- l1hs %>%
    group_by(rel_cpg_pos) %>%
    summarise(mm = mean(pctM)) %>%
    ungroup() %>%
    ggplot(aes(x = rel_cpg_pos, y = mm)) +
    geom_line() +
    mtopen
p2 <- color_intervals %>%
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
    geom_text(aes(x = -200 + ((start + end) / 2), y = 1.5, label = feature)) +
    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
    scale_fill_paletteer_d("dutchmasters::milkmaid") +
    theme_map() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0, 0.4)) +
    theme(legend.position = "none")
p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
mysaveandstore(sprintf("%s/l1hs_all.pdf", outputdir), w = 12, h = 6)

p1 <- l1hs %>%
    group_by(rel_cpg_pos, condition) %>%
    summarise(mm = mean(pctM)) %>%
    ungroup() %>%
    ggplot(aes(x = rel_cpg_pos, y = mm)) +
    geom_line(aes(color = condition), alpha = 0.5) +
    scale_conditions +
    mtopen
p2 <- color_intervals %>%
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
    geom_text(aes(x = -200 + ((start + end) / 2), y = 1.5, label = feature)) +
    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
    scale_fill_paletteer_d("dutchmasters::milkmaid") +
    theme_map() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0, 0.4)) +
    theme(legend.position = "none")
p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
mysaveandstore(sprintf("%s/l1hs_condition.pdf", outputdir), w = 12, h = 6)

library(zoo)
pf <- l1hs %>%
    group_by(rel_cpg_pos) %>%
    summarise(mm = mean(pctM)) %>%
    ungroup()
pf$rmm <- rollmean(pf$mm, k = 10, align = "center", na.pad = TRUE)

p1 <- pf %>%
    ggplot(aes(x = rel_cpg_pos, y = rmm)) +
    geom_line() +
    mtclosedgrid
p2 <- color_intervals %>%
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.25, ymax = 0.75), fill = "darkgrey") +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature), alpha = 1) +
    geom_text(aes(x = -200 + ((start + end) / 2), y = 1.5, label = feature)) +
    coord_cartesian(xlim = layer_scales(p1)$x$range$range) +
    scale_fill_paletteer_d("dutchmasters::milkmaid") +
    theme_map() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0, 0.4)) +
    theme(legend.position = "none")
p <- p2 / p1 + plot_layout(heights = c(0.2, 1))
mysaveandstore(sprintf("%s/rm_l1hs_all.pdf", outputdir), w = 6, h = 4)


# l1hsN %>%
#     filter(gene_id == "L1HS_3p12.3_3") %>%
#     mutate(rel_cpg_pos = (rte_end + element_start) - start) %>% dplyr::select(sample,start, end, rte_start, rte_end, element_start, element_end, rel_cpg_pos ) %>%
#     arrange(sample, rel_cpg_pos) %>% print(n = 1000)

# l1hsN %$% gene_id %>% table()
