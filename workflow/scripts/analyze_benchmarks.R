module_name <- "srna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)
library(magrittr)
library(ggpubr)
library(patchwork)

for (module in c("aref", "srna", "ldna", "lrna")) {
    tryCatch(
        {
            files <- list.files(sprintf("%s/benchmarks", module), pattern = ".*.tsv", recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
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

            default_time <- 70
            default_mem <- 12.5

            p1 <- ggplot(bf, aes(x = s / 60, y = rule)) +
                geom_point() +
                geom_vline(xintercept = default_time) +
                labs(y = "", x = "time(min)") +
                mtclosedgridv
            mysaveandstore(pl = p1, file.path(module, "benchmarks", "time.pdf"), 6, bf %$% rule %>% unique() %>% length() / 3)

            p2 <- ggplot(bf, aes(x = cpu_time / 60, y = rule)) +
                geom_point() +
                labs(y = "", x = "cpu_time(min)") +
                mtclosedgridv
            mysaveandstore(pl = p2, file.path(module, "benchmarks", "cpu_time.pdf"), 6, bf %$% rule %>% unique() %>% length() / 3)

            p3 <- ggplot(bf, aes(x = max_rss / 1000, y = rule)) +
                geom_point() +
                geom_vline(xintercept = default_mem) +
                labs(y = "", x = "max_rss(Gb)") +
                mtclosedgridv
            mysaveandstore(pl = p3, file.path(module, "benchmarks", "max_rss.pdf"), 6, bf %$% rule %>% unique() %>% length() / 3)


            ptch <- (p1 | p2 | p3) + plot_layout(axes = "collect")
            mysaveandstore(pl = ptch, file.path(module, "benchmarks", "all.pdf"), 10, bf %$% rule %>% unique() %>% length() / 3)
        },
        error = function(e) {
        }
    )
}
