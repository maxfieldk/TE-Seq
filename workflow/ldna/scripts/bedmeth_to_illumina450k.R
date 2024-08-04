module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

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
library()
# library(ReMapEnrich)
# library(msigdbr)
library(Biostrings)
library(PAMES)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethylpaths = sprintf("ldna/intermediates/%s/methylation/hg38/%s_CG_bedMethyl.bed", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(csv = "ldna/intermediates/agg/cpg_formatted_methylation_filtered.csv"), env = globalenv())
    }
)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE), col_names = FALSE)
    df_m <- df %>% filter(X4 == "m")
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
grsdf %$% seqnames %>% unique()

cnv <- illumina450k_hg38 %>% tibble()

cnvgr <- makeGRangesFromDataFrame(cnv, start.field = "Genomic_Coordinate", end.field = "Genomic_Coordinate", keep.extra.columns = TRUE)
start(cnvgr) <- start(cnvgr) - 1
end(cnvgr) <- end(cnvgr) - 1

mbo <- mergeByOverlaps(grs, cnvgr)
mbodf <- as.data.frame(mbo) %>% tibble()
mbodf

bdf <- mbodf %>%
    dplyr::select(cnvgr.Probe, grs.pctM, grs.cov, grs.sample) %>%
    dplyr::rename(probe = cnvgr.Probe, pctM = grs.pctM, cov = grs.cov, sample = grs.sample)

wdf <- bdf %>%
    dplyr::select(-cov) %>%
    mutate(pctM = pctM / 100) %>%
    pivot_wider(names_from = sample, values_from = pctM)
path <- outputs$csv
dir.create(path, recursive = TRUE)
write_csv(wdf %>% filter(!grepl("^ch", probe)) %>% drop_na(), path)


# library(methylclock)
# library(impute)
# library(MethylationData)


# missCpGs <- checkClocks(wdf %>% filter(!grepl("^ch",probe)) %>% drop_na())

# res <- DNAmAge(wdf %>% filter(!grepl("^ch",probe)) %>% drop_na(), fastImp=TRUE)

# reswa <- res %>% left_join(sample_table %>% dplyr::rename(id = sample_name))

# cor(reswa$Horvath, reswa$age)
