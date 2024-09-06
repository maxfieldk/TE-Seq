module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(magrittr)
library(forcats)

library(configr)
library(calcPCBrainAge)
set.seed(123)



tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            csv = "ldna/intermediates/agg/cpg_formatted_methylation.csv"
        ), env = globalenv())
        assign("outputs", list(
            csv = "ldna/results/tables/clocks/R_clock_results.csv"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$csv)
samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name
dft <- df %>%
    column_to_rownames("probe") %>%
    t() %>%
    as.data.frame()
df <- read_csv(inputs$csv)
df %>%
    as.data.frame() %>%
    library(calcPCBrainAge)
## basic example code
myPCBrainAges <- calcPCBrainAge(DNAm = dft) # This gives you a vector
datMeth <- dft


calcPCBrainAge(DNAm = dft, CpGImputation = imputeMissingBrainCpGs)
as.matrix(df) %>% head()
aa <- df %>%
    column_to_rownames("probe") %>%
    as.matrix()
t(aa) %>% str()
###############################
phenodf <- sample_table %>%
    dplyr::select(sample_name, age) %>%
    dplyr::rename(ID = sample_name, Age = age)
source("/users/mkelsey/data/Nanopore/alz/CorticalClock/PredCorticalAge/CorticalClock.r")
CorticalClock(betas = aa, pheno = phenodf, dir = "/users/mkelsey/data/Nanopore/alz/CorticalClock/PredCorticalAge/", IDcol = "ID", Agecol = "Age")
