module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)
library(magrittr)
library(ggpubr)

files <- list.files("srna/benchmarks", pattern = ".*.tsv",recursive = TRUE, include.dirs = TRUE, full.names = TRUE)

rm(bf)
for (file in files) {
    if (!exists("bf")) {
        bf <- readr::read_delim(file, delim = "\t")
        bf$rule <- dirname(file) %>% basename()
    } else {
        data <- readr::read_delim(file, delim = "\t")
        data$rule <- dirname(file) %>% basename()
        bf <- dplyr::bind_rows(bf, data)
    }
}

p1 <- ggplot(bf, aes(x = s, y = rule)) + 
    geom_point() + 
    labs(y = "", x = "time(s)") + mtclosedgridv
mysaveandstore(pl = p1, file.path("srna", "benchmarks", "time.pdf"), 6,bf %$% rule %>% unique %>% length/3)

p2 <- ggplot(bf, aes(x = cpu_time, y = rule)) + 
    geom_point() + 
    labs(y = "", x = "cpu_time(s)") + mtclosedgridv
mysaveandstore(pl = p2, file.path("srna", "benchmarks", "cpu_time.pdf"), 6,bf %$% rule %>% unique %>% length/3)

p3 <- ggplot(bf, aes(x = max_rss, y = rule)) + 
    geom_point() + 
    labs(y = "", x = "max_rss(Gb)") + mtclosedgridv
mysaveandstore(pl = p3, file.path("srna", "benchmarks", "max_rss.pdf"), 6,bf %$% rule %>% unique %>% length/3)


library(patchwork)
ptch <- (p1 | p2 | p3 ) + plot_layout(axes = "collect")
mysaveandstore(pl = ptch, file.path("srna", "benchmarks", "all.pdf"), 10,bf %$% rule %>% unique %>% length/3)
