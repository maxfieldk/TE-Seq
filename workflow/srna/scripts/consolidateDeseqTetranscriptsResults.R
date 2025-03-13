set.seed(123)
library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")
# conf <- configr::read.config(file = "conf/config.yaml")[["lrna"]]
module_name <- snakemake@params$module_name
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "inputdir" = "results/agg/deseq",
            "outputdir" = "results/agg/repeatanalysis"
        ), env = globalenv())
        assign("outputs", list(
            "resultsdf" = "results/agg/repeatanalysis/resultsdf.tsv"
        ), env = globalenv())
    }
)

multi <- read_delim(grep("multi", inputs$resultsdf, value = TRUE), delim = "\t") %>% mutate(tecounttype = "multi")
uniq <- read_delim(grep("uniq", inputs$resultsdf, value = TRUE), delim = "\t") %>% mutate(tecounttype = "unique")
resultsdf <- full_join(multi, uniq)
write_delim(resultsdf, outputs$resultsdf, delim = "\t")
