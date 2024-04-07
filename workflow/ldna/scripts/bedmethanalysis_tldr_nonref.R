source("workflow/scripts/defaults.R")
module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
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
    library(karyoploteR)
    library(glob)
####################
{
    chromosomesAll <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    chromosomes <- c(paste0("chr", 1:22), "chrX")
    chromosomesNoX <- c(paste0("chr", 1:22))
}

sample_table <- read_csv("conf/private/sample_table.csv")
conditions <- unique(sample_table$condition)
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

sample_grs <- list()
sample_beds <- list()
for (sample_name in sample_table$sample_name) {
    dir <- paste0("tldr/", sample_name)
    files <- list.files(path = dir, pattern = "*.bedmethyl.bed", full.names = TRUE)
    df_list <- lapply(files, read_table, col_names = FALSE, cols(
        X1 = col_character(),
        X2 = col_double(),
        X3 = col_double(),
        X4 = col_character(),
        X5 = col_double(),
        X6 = col_character(),
        X7 = col_double(),
        X8 = col_double(),
        X9 = col_number(),
        X10 = col_double(),
        X11 = col_double(),
        X12 = col_double(),
        X13 = col_double(),
        X14 = col_double(),
        X15 = col_double(),
        X16 = col_double(),
        X17 = col_double(),
        X18 = col_double()
    ))
    df <- bind_rows(df_list)
    gr <- GRanges(
        seqnames = df$X1,
        ranges = IRanges(start = df$X2, end = df$X2),
        cov = df$X10,
        pctM = as.double(df$X11)
    )
    gr$sample <- sample_name
    gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_grs[[sample_name]] <- gr

    files <- list.files(path = dir, pattern = ".bed$", full.names = TRUE)
    files <- files[!grepl("bedmethyl", files)]
    bed_list <- lapply(files, read_delim, col_names = FALSE, col_types = cols(X1 = "c", X2 = "d", X3 = "d", X4 = "c", X5 = "c"))
    bed <- bind_rows(bed_list)
    bed$sample <- sample_name
    bed$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_beds[[sample_name]] <- bed
}

grs <- Reduce(c, sample_grs)
beds <- Reduce(rbind, sample_beds)
colnames(beds) <- c("seqnames", "start", "end", "name", "strand", "sample", "condition")
######
    # PREP DATA FOR ANALYSIS
    conditions <- conf$conditions
    condition1 <- conditions[1]
    condition2 <- conditions[2]
    condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
    condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

    sample_grs <- list()
    for (sample_name in samples) {
        df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethlpaths, value = TRUE), col_names = FALSE)
        df_m <- df %>% filter(X4 == "m")
        df_h <- df %>% filter(X4 == "h")
        rm(df)
        gr <- GRanges(
            seqnames = df_m$X1,
            ranges = IRanges(start = df_m$X2, end = df_m$X2),
            cov = df_m$X10,
            pctM = as.double(df_m$X11)
        )
        gr$sample <- sample_name
        gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
        sample_grs[[sample_name]] <- gr
    }

    grs <- Reduce(c, sample_grs)
    rm(sample_grs)
    # filter out low coverage and ensure that all samples have the same cpgs
    grs <- grs[grs$cov > MINIMUMCOVERAGE]
    grsdf <- tibble(as.data.frame(grs))
    grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
    seqnames <- grsdf$seqnames
    start <- grsdf$start
    end <- grsdf$end
    pos <- paste0(seqnames, "_", start, "_", end)
    grsdf$pos <- pos
########

l1s <- beds %>% filter(grepl("L1", name)) %$% seqnames
l1bed <- beds %>% filter(seqnames %in% l1s)
grsdf <- as.data.frame(grs) %>% tibble()
df <- grsdf %>% left_join(beds, by = c("seqnames", "sample", "condition"))
df <- df %>%
    rename(type = name, start = start.x, stop = end.x, strand = strand.x, element_start = start.y, element_stop = end.y, element_strand = strand.y, uid = seqnames) %>%
    mutate(element_length = element_stop - element_start)

l1hs_fl_non_ref <- df %>%
    filter(grepl("L1Ta", type)) %>%
    filter(element_length > 6000)
l1hs_fl_non_ref %$% uid %>% unique()
# reference insertions
rtedf <- read_delim("Rintermediates/rtedf.tsv", col_names = TRUE)
# l1hs_ref <- rtedf %>% filter(type == "L1HS")
# l1hs_ref <- l1hs_ref %>% separate(uid, c("element_seqnames", "element_start", "element_end", "element_strand"), sep = "_", remove = FALSE, convert = TRUE)
# l1hs_ref <- l1hs_ref %>% mutate(element_length = element_end - element_start)
# l1hs_fl_ref <- l1hs_ref %>% filter(element_length > 6000)
# l1hs_fl_ref %$% uid %>% unique()

l1hsintactmethdf <- rtedf %>% filter(glob_col_axis == "L1HS Intact\nn=89")
l1hsintactmethgr %$% uid %>% unique()

nonref <- l1hs_fl_non_ref %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
        element_strand == "-" ~ (start > element_stop - 909) & (start < element_stop)
    )) %>%
    group_by(uid, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(ref = "non-ref")
ref <- l1hsintactmethdf %>%
    separate(uid, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_start + 909),
        element_strand == "-" ~ (start > element_stop - 909) & (start < element_stop)
    )) %>%
    group_by(uid, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(ref = "ref")
p <- rbind(nonref, ref) %>% ggplot() +
    geom_beeswarm(aes(x = ref, y = mean, color = condition), alpha = 0.5, cex = 1.3) +
    stat_summary(aes(x = ref, y = mean), fun = mean, geom = "point", cex = 2) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
png("results/plots/tldr/beeswarm_5utr.png", width = 8, height = 4, units = "in", res = 300)
print(p)
dev.off()
p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = ref, y = mean, fill = condition), outlier.alpha = 0) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
png("results/plots/tldr/box_hideoutliers_5utr.png", width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()

p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = ref, y = mean, fill = condition)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
png("results/plots/tldr/box_5utr.png", width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()

nonref <- l1hs_fl_non_ref %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_stop),
        element_strand == "-" ~ (start > element_start) & (start < element_stop)
    )) %>%
    group_by(uid, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(ref = "non-ref")
ref <- l1hsintactmethdf %>%
    separate(uid, sep = "_", into = c("element_chr", "element_start", "element_stop", "element_strand"), convert = TRUE, remove = FALSE) %>%
    filter(case_when(
        element_strand == "+" ~ (start > element_start) & (start < element_stop),
        element_strand == "-" ~ (start > element_start) & (start < element_stop)
    )) %>%
    group_by(uid, sample, condition) %>%
    summarise(mean = mean(pctM)) %>%
    mutate(ref = "ref")
p <- rbind(nonref, ref) %>% ggplot() +
    geom_beeswarm(aes(x = ref, y = mean, color = condition), alpha = 0.5, cex = 1.3) +
    stat_summary(aes(x = ref, y = mean), fun = mean, geom = "point", cex = 2) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions

png("results/plots/tldr/beeswarm.png", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()
p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = ref, y = mean, fill = condition), outlier.alpha = 0) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions

png("results/plots/tldr/box_hideoutliers.png", width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()
p <- rbind(nonref, ref) %>% ggplot() +
    geom_boxplot(aes(x = ref, y = mean, fill = condition)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "", y = "Mean element methylation") +
    ggtitle("Reference vs Non-reference L1HS") +
    mtopen + scale_conditions
png("results/plots/tldr/box.png", width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()


sample_tables <- list()
for (sample_name in sample_table$sample_name) {
    tbl <- read_delim(paste0("tldr/", sample_name, ".filtered.table.txt"))
    tbl$sample <- sample_name
    tbl$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_tables[[sample_name]] <- tbl
}
idf <- Reduce(rbind, sample_tables)


my_palette[c(1, 2)]
library(BSgenome.Hsapiens.UCSC.hs1)

il1df <- idf %>%
    filter(Subfamily == "L1Ta") %>%
    mutate(length_te = EndTE - StartTE) %>%
    filter(length_te > 6000)

shared_nonref_counts <- tibble(counts = countOverlaps(il1grs, il1grs))
p <- ggplot(shared_nonref_counts, aes(x = counts)) +
    geom_bar(aes(fill = counts > 1)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = "Insertion found in N samples", y = "Count", fill = "Shared") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggtitle("Non-reference fl-L1HS Insertion") +
    mtopen + scale_conditions
png("results/plots/tldr/shared_nonref_hist.png", width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()

il1grs <- GRanges(
    seqnames = il1df$Chrom,
    ranges = IRanges(start = il1df$Start, end = il1df$End),
    strand = il1df$Strand,
    sample = il1df$sample,
    condition = il1df$condition
)

il1grs_alz <- il1grs[il1grs$condition == "alz"]
il1grs_ctrl <- il1grs[il1grs$condition == "ctrl"]

png("results/plots/tldr/l1hsideogram.png", height = 12, width = 12.5, res = 300, units = "in")
plot.new()
kp <- plotKaryotype(genome = BSgenome.Hsapiens.UCSC.hs1, plot.type = 2, chromosomes = chromosomes, cex = 2)
kpPlotRegions(kp, il1grs_alz, col = my_palette[1], data.panel = 1)
kpPlotRegions(kp, il1grs_ctrl, col = my_palette[2], data.panel = 2)
legend(x = "bottomright", fill = my_palette[c(1, 2)], legend = c("Alz", "Ctrl"))
dev.off()


# Create a basic circular karyotype with 4 sectors
png("results/plots/tldr/circlize.png", height = 12, width = 12.5, res = 300, units = "in")
plot.new()
circos.par("start.degree" = 90)
circos.par("gap.degree" = 4)
circos.par("cell.padding" = c(0.1, 0.1))
circos.par("track.height" = 0.4)
circos.initialize(factors = unique(seqnames(il1grs)), xlim = c(0, 100))
circos.highlight(il1grs, col = il1grs$condition, border = "black", track.height = 0.2, alpha = 0.7)
dev.off()

# Create a basic circular karyotype with 4 sectors
circos.par("start.degree" = 90)
circos.par("gap.degree" = 4)
circos.par("cell.padding" = c(0.1, 0.1))
circos.par("track.height" = 0.4)

# Define the sectors
sector_names <- unique(seqnames(il1grs))
sector_df <- data.frame(
    chromosome = sector_names,
    start = seq(1, length(sector_names), by = 1),
    end = seq(2, length(sector_names) + 1, by = 1)
)
rownames(sector_df) <- sector_names

# Initialize the circular karyotype
circos.initialize(sector_df)

# Plot the features on the karyotype
highlight(il1grs, col = il1grs$condition, border = "black", track.height = 0.2, alpha = 0.7)
