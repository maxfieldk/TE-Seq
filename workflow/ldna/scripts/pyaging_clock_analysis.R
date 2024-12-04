module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
library(rtracklayer)
library(Biostrings)
library(cowplot)
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
set.seed(123)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            csv = "ldna/results/tables/clocks/pyaging_clock_results.csv"
        ), env = globalenv())
        assign("outputs", list(
            outfile = "ldna/results/plots/clocks/pyaging_clock_plots.outfile"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$outfile)
samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

df <- read_csv(inputs$csv) %>% dplyr::rename(sample_name = `...1`)
df <- sample_table %>% full_join(df)



df <- df %>%
    mutate(age_acc_altum = altumage - age) %>%
    mutate(age_acc_grim = pcgrimage - age) %>%
    mutate(age_acc_grim2 = grimage2 - age)

library(ggpubr)
p <- df %>% ggscatter(x = "altumage", y = "age", color = "condition") +
    mtopengrid +
    scale_conditions
mysaveandstore(sprintf("%s/altum_scatter.pdf", outputdir))

p <- df %>% ggscatter(x = "altumage", y = "age", color = "condition") +
    lims(x = c(0, 100), y = c(0, 100)) +
    geom_abline(slope = 1, intercept = 0) +
    mtopen +
    scale_conditions
mysaveandstore(sprintf("%s/altum_scatter_lims_set.pdf", outputdir), 4, 4)




p <- df %>% ggscatter(x = "pcgrimage", y = "age", color = "condition") +
    mtopengrid +
    scale_conditions
mysaveandstore(sprintf("%s/grim_scatter.pdf", outputdir))

p <- df %>% ggscatter(x = "grimage2", y = "age", color = "condition") +
    mtopengrid +
    scale_conditions
mysaveandstore(sprintf("%s/grim2_scatter.pdf", outputdir))


p <- df %>% ggscatter(x = "braak", y = "age_acc_altum") +
    mtopengrid
mysaveandstore(sprintf("%s/altum_resid.pdf", outputdir))

model <- lm(age_acc_altum ~ braak, df)
sf <- broom::tidy(summary(model))
p <- df %>% ggplot(aes(x = braak, y = age_acc_altum)) +
    geom_point() +
    geom_smooth(aes(x = unclass(braak), color = "1"),
        formula = y ~ x,
        method = lm, se = FALSE
    ) +
    geom_smooth(aes(x = unclass(braak), color = "2"),
        formula = y ~ poly(x, 2),
        method = lm, se = FALSE
    ) +
    geom_smooth(aes(x = unclass(braak), color = "3"),
        formula = y ~ poly(x, 3),
        method = lm, se = FALSE
    ) +
    scale_color_discrete("Trend", labels = c("linear", "quadratic", "cubic")) +
    mtclosed
mysaveandstore(sprintf("%s/altum_resid_fit.pdf", outputdir), sf = sf)


p <- df %>% ggscatter(x = "braak", y = "age_acc_grim") +
    mtopengrid
mysaveandstore(sprintf("%s/grim_resid.pdf", outputdir))

p <- df %>% ggscatter(x = "braak", y = "age_acc_grim2") +
    mtopengrid
mysaveandstore(sprintf("%s/grim2_resid.pdf", outputdir))
