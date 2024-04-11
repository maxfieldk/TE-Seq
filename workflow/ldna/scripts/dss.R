module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]

library(DSS)
library(BiocParallel)
library(readr)
require(bsseq)
library(readr)


tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("inputs", list(
            "data" = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            "r_annotation_fragmentsjoined" = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            "r_repeatmasker_annotation" = "annotations/repeatmasker_annotation.csv"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "annotations/rte_beds/outfile.txt"
        ), env = globalenv())
    }
)





sample_dfs <- list()
for (sample in sample_table$sample_name) {
    sample_dfs[[sample]] <- read.table(grep(sprintf("/%s/", sample), inputs$data, value = TRUE), header = TRUE)
}

BSobj <- makeBSseqData(sample_dfs, names(sample_dfs))
mParam <- MulticoreParam(workers = 12, progressbar = TRUE)

conditions <- conf$levels
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name
# need to adjust numbering given idiosyncracies of dmltest
dmlTest <- DMLtest(BSobj, group1 = condition2samples, group2 = condition1samples, smoothing = TRUE)

tryCatch(
    {
        head(dmlTest)
        dmls <- callDML(dmlTest)
        head(dmls)
        dmlcols <- colnames(dmls)
        dmlcols1 <- gsub("1", "_ctwo", dmlcols)
        dmlcols2 <- gsub("2", "_cone", dmlcols1)
        dmlcols3 <- gsub("two", "2", dmlcols2)
        dmlcols4 <- gsub("one", "1", dmlcols3)
        dmlcols5 <- gsub("diff", "diff_c2_minus_c1", dmlcols4)
        colnames(dmls) <- dmlcols5
        head(dmls)
    },
    error = function(e) {
        print("no DMLs")
        dmls <- data.frame()
    }
)

tryCatch(
    {
        dmrs <- callDMR(dmlTest)
        dmrs <- as.data.frame(dmrs)
        head(dmrs)
        dmrcols <- colnames(dmrs)
        dmrcols1 <- gsub("1", "_ctwo", dmrcols)
        dmrcols2 <- gsub("2", "_cone", dmrcols1)
        dmrcols3 <- gsub("two", "2", dmrcols2)
        dmrcols4 <- gsub("one", "1", dmrcols3)
        dmrcols5 <- gsub("diff.Methy", "diff_c2_minus_c1", dmrcols4)
        colnames(dmrs) <- dmrcols5
    },
    error = function(e) {
        print("no DMRs")
        dmrs <- data.frame()
    }
)
# save results
options(scipen = 500)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write_delim(dmls, outputs$dmls, delim = "\t", col_names = TRUE)
write_delim(dmrs, outputs$dmrs, delim = "\t", col_names = TRUE)
