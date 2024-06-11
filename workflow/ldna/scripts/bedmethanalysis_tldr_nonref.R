module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

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
library(Biostrings)
# library(karyoploteR)
# library(glob)
####################
tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        params <- snakemake@params
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = sprintf("aref/%s_annotations/%s_repeatmasker.gtf.rformatted.fragmentsjoined.csv", sample_table$sample_name, sample_table$sample_name),
            r_repeatmasker_annotation = sprintf("aref/%s_annotations/%s_repeatmasker_annotation.csv",sample_table$sample_name, sample_table$sample_name),
            sample_refs = sprintf("aref/%s.fa", sample_table$sample_name),
            blast_njs = sprintf("aref/%s.njs", sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(
            plots = "ldna/results/plots/nonref_insert_meth_analysis.rds"
        ), env = globalenv())
    }
)

{
    chromosomesAll <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    chromosomes <- c(paste0("chr", 1:22), "chrX")
    chromosomesNoX <- c(paste0("chr", 1:22))
}

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

conditions <- unique(sample_table$condition)
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

anns <- list()
for (sample in sample_table$sample_name) {
    df1 <- read_csv(grep(sprintf("%s_annotations", sample), inputs$r_annotation_fragmentsjoined, value = TRUE))
    df1$sample_name <- sample
    df2 <- read_csv(grep(sprintf("%s_annotations", sample), inputs$r_repeatmasker_annotation, value = TRUE))
    df <- df1 %>% left_join(sample_table) %>% left_join(df2) 
    anns[[sample]] <- df
}
rmann <- bind_rows(anns)
ann <- rmann %>% filter(refstatus == "NonRef")
grsdf <- read_delim("ldna/Rintermediates/grsdf_nonref.tsv", col_names = TRUE)


########
l1s <- ann %>% filter(rte_subfamily == "L1HS") %>% filter(str_detect(intactness_req, "Intact"))



df <- grsdf %>% left_join(l1s, by = c("seqnames", "sample" = "sample_name", "condition"))
df <- df %>% dplyr::select(-element_start, -element_end) %>%
    dplyr::rename(start = start.x, stop = end.x, strand = strand.x, element_start = start.y, element_stop = end.y, element_strand = strand.y, uid = seqnames) %>%
    mutate(element_length = element_stop - element_start)


# reference insertions
rtedf_promoters <- read_delim("ldna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
rtedf <- read_delim("ldna/Rintermediates/rtedf.tsv", col_names = TRUE)

# l1hs_ref <- rtedf %>% filter(type == "L1HS")
# l1hs_ref <- l1hs_ref %>% separate(uid, c("element_seqnames", "element_start", "element_end", "element_strand"), sep = "_", remove = FALSE, convert = TRUE)
# l1hs_ref <- l1hs_ref %>% mutate(element_length = element_end - element_start)
# l1hs_fl_ref <- l1hs_ref %>% filter(element_length > 6000)
# l1hs_fl_ref %$% uid %>% unique()

l1hsintactmethpromotersdf <- rtedf_promoters %>% filter(rte_subfamily == "L1HS") %>% filter(str_detect(intactness_req, "Intact"))

nonref <- df %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
        element_strand == "-" ~ (start > element_stop - 909) & (start < element_stop)
    )) %>%
    group_by(gene_id, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(refstatus = "NonRef")
ref <- l1hsintactmethpromotersdf %>%
    group_by(gene_id, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(refstatus = "Ref")
p <- rbind(nonref, ref) %>% ggplot() +
    geom_beeswarm(aes(x = refstatus, y = mean, color = condition), alpha = 0.5, cex = 1) +
    stat_summary(aes(x = refstatus, y = mean), fun = mean, geom = "point", cex = 2) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
mysave("ldna/results/plots/tldr/beeswarm_5utr.png", 6, 4)



p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = refstatus, y = mean, fill = condition), outlier.alpha = 0) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
mysave("ldna/results/plots/tldr/box_hideoutliers_5utr.png", 6, 4)

p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = refstatus, y = mean, fill = condition)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
mysave("ldna/results/plots/tldr/box_5utr.png", 6, 4)


# Full element methylation
nonref <- df %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_stop),
        element_strand == "-" ~ (start > element_start) & (start < element_stop)
    )) %>%
    group_by(gene_id, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(refstatus = "NonRef")

l1hsintactmethdf <- rtedf %>% filter(rte_subfamily == "L1HS") %>% filter(str_detect(intactness_req, "Intact"))

ref <- l1hsintactmethdf %>%
    group_by(gene_id, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(refstatus = "Ref")
    
p <- rbind(nonref, ref) %>% ggplot() +
    geom_beeswarm(aes(x = refstatus, y = mean, color = condition), alpha = 0.5, cex = 1.3) +
    stat_summary(aes(x = refstatus, y = mean), fun = mean, geom = "point", cex = 2) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
mysaveandstore("ldna/results/plots/tldr/beeswarm.pdf", 6, 4)

p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = refstatus, y = mean, fill = condition), outlier.alpha = 0) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
mysaveandstore("ldna/results/plots/tldr/box_hideoutliers.pdf", 6, 4)

p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = refstatus, y = mean, fill = condition)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions

mysaveandstore("ldna/results/plots/tldr/box.pdf", 6, 4)



# sample_tables <- list()
# for (sample_name in sample_table$sample_name) {
#     tbl <- read_delim(paste0("tldr/", sample_name, ".filtered.table.txt"))
#     tbl$sample <- sample_name
#     tbl$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
#     sample_tables[[sample_name]] <- tbl
# }
# idf <- Reduce(rbind, sample_tables)


# my_palette[c(1, 2)]
# library(BSgenome.Hsapiens.UCSC.hs1)

# il1df <- idf %>%
#     filter(Subfamily == "L1Ta") %>%
#     mutate(length_te = EndTE - StartTE) %>%
#     filter(length_te > 6000)

# shared_nonref_counts <- tibble(counts = countOverlaps(il1grs, il1grs))
# p <- ggplot(shared_nonref_counts, aes(x = counts)) +
#     geom_bar(aes(fill = counts > 1)) +
#     theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
#     labs(x = "Insertion found in N samples", y = "Count", fill = "Shared") +
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
#     ggtitle("Non-reference fl-L1HS Insertion") +
#     mtopen + scale_conditions
# png("results/plots/tldr/shared_nonref_hist.png", width = 4, height = 4, units = "in", res = 300)
# print(p)
# dev.off()

# il1grs <- GRanges(
#     seqnames = il1df$Chrom,
#     ranges = IRanges(start = il1df$Start, end = il1df$End),
#     strand = il1df$Strand,
#     sample = il1df$sample,
#     condition = il1df$condition
# )

# il1grs_alz <- il1grs[il1grs$condition == "alz"]
# il1grs_ctrl <- il1grs[il1grs$condition == "ctrl"]

# png("results/plots/tldr/l1hsideogram.png", height = 12, width = 12.5, res = 300, units = "in")
# plot.new()
# kp <- plotKaryotype(genome = BSgenome.Hsapiens.UCSC.hs1, plot.type = 2, chromosomes = chromosomes, cex = 2)
# kpPlotRegions(kp, il1grs_alz, col = my_palette[1], data.panel = 1)
# kpPlotRegions(kp, il1grs_ctrl, col = my_palette[2], data.panel = 2)
# legend(x = "bottomright", fill = my_palette[c(1, 2)], legend = c("Alz", "Ctrl"))
# dev.off()
