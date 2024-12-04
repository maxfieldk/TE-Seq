module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
# library(zoo)
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
# library(msigdbr)
library(Biostrings)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", samples, samples),
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/tables/dmrs.CG_m.tsv",
            dmls = "ldna/results/tables/dmls.CG_m.tsv",
            read_mods = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

### gimmemotifs
gimmehyper <- read_delim("ldna/results/gimme/hyper/gimme.roc.report.txt", delim = "\t")
gimmehypo <- read_delim("ldna/results/gimme/hypo/gimme.roc.report.txt", delim = "\t")

gimmehyper$Motif <- gsub("_HUMAN.H11MO.0.*", "", gimmehyper$Motif)
pl_gimme_hyper <- gimmehyper %>%
    arrange(`P-value`) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Motif, `log10 P-value`), y = `log10 P-value`)) +
    geom_col(color = "black", fill = mycolor) +
    scale_x_discrete(labels = scales::label_wrap(60)) +
    coord_flip() +
    ggtitle("Hyper Top 10 Motifs") +
    scale_color_continuous(trans = "reverse") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "") +
    theme(
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
png("ldna/results/plots/gimme/gimmebarplothyper.png", 3, 4, units = "in", res = 300)
print(pl_gimme_hyper)
dev.off()

gimmehypo %>%
    arrange(`P-value`) %>%
    head(10)
gimmehypo$Motif <- gsub("_HUMAN.H11MO.0.*", "", gimmehypo$Motif)
pl_gimme_hypo <- gimmehypo %>%
    arrange(`P-value`) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Motif, `log10 P-value`), y = `log10 P-value`)) +
    geom_col(color = "black", fill = mycolor) +
    scale_x_discrete(labels = scales::label_wrap(60)) +
    coord_flip() +
    scale_color_continuous(trans = "reverse") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "") +
    ggtitle("Hypo Top 10 Motifs") +
    theme(
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
png("ldna/results/plots/gimme/gimmebarplothypo.png", 3, 4, units = "in", res = 300)
print(pl_gimme_hypo)
dev.off()

p <- ggplot(gimmehypo, aes(x = `Recall at 10% FDR`)) +
    geom_histogram()
png("ldna/results/plots/gimme/histogramPval.png", 4, 4, units = "in", res = 300)
print(p)
dev.off()
