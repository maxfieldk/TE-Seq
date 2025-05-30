module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    arrange(condition)

set.seed(123)
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "counttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = sprintf("%s/results/agg/deseq", conf$prefix),
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())
        assign("outputs", list(
            stats = sprintf("srna/results/agg/repeatanalysis/telescope_multi/te_group_stats.csv")
        ), env = globalenv())
        assign("inputs", list(
            rte_vstcounts = sprintf("srna/results/agg/deseq/telescope_multi/rtes/batchRemoved_no/vst_counts.csv"),
            sizefactors = "srna/results/agg/deseq/telescope_multi/sizefactors.csv"
        ), env = globalenv())
    }
)

vst_counts <- read_csv(inputs$rte_vstcounts)

counttype <- params[["counttype"]]
outputdir <- "srna/results/agg/rtepca"
rmann <- get_repeat_annotations(
    default_or_extended = "default",
    keep_non_central = FALSE,
    append_NI_samplename_modifier = FALSE
)

l1hsfl <- rmann %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(rte_length_req == "FL")
l1hs <- rmann %>%
    filter(rte_subfamily == "L1HS")
vsttemp <- vst_counts %>%
    filter(gene_id %in% l1hs$gene_id) %>%
    dplyr::select(-gene_id)

coldata <- sample_table %>% column_to_rownames("sample_name")

## PCA plots
pcaObj <- pca(vsttemp, metadata = coldata, removeVar = 0.1)

p <- screeplot(pcaObj, title = "") + mtopen + anchorbar
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "screeplot.pdf", sep = "/"), 4, 4)
p <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 4
) + mtopen
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "loadings.pdf", sep = "/"), 10, 7)

if (any(grepl("batchCat", colnames(coldata)))) {
    categorical_batch_vars <- grep("batchCat", colnames(coldata), value = TRUE)
    for (catbatchvar in categorical_batch_vars) {
        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
            colby = "condition", legendPosition = "right", shape = catbatchvar,
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen + scale_conditions
        mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), str_glue("pca_shapeby_{catbatchvar}.pdf"), sep = "/"), 5, 5)

        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
            colby = catbatchvar, legendPosition = "right",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen
        mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), str_glue("pca_batch_colorby_{catbatchvar}.pdf"), sep = "/"), 5, 5)

        p <- biplot(pcaObj,
            showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
            colby = catbatchvar, legendPosition = "right",
            labSize = 5, pointSize = 5, sizeLoadingsNames = 5
        ) + mtopen
        mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), str_glue("pca_batch_large_colorby_{catbatchvar}.pdf"), sep = "/"), 16, 16)

        p <- pairsplot(pcaObj, colby = catbatchvar, title = "Batch", legendPosition = "right")
        mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), str_glue("pca_pairs_batch_colorby_{catbatchvar}.pdf"), sep = "/"), 15, 15)
    }
} else {
    print("no batchCat")
}
batchestemp <- grep("batch", colnames(coldatatemp), value = TRUE)
tryCatch(
    {
        p <- eigencorplot(pcaObj, metavars = c(batchestemp, "condition"))
        mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "eigencor.pdf", sep = "/"), 8, 4)
    },
    error = function(e) {
        print("eigencorplot fail")
    }
)

p <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) + mtopen + scale_conditions
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca.pdf", sep = "/"), 5, 5)

p <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) + mtopen + scale_conditions
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_large.pdf", sep = "/"), 16, 16)


p <- pairsplot(pcaObj, colby = "condition", title = "Condition", legendPosition = "right")
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca_pairs_condition.pdf", sep = "/"), 15, 15)


p <- biplot(pcaObj,
    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) + mtopen + scale_conditions
mysaveandstore(paste(outputdir, counttype, subset, sprintf("batchRemoved_%s", batchnormed), "pca34.pdf", sep = "/"), 4, 4)
