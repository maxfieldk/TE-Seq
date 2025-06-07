library(readr)
library(magrittr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
# First argument
data_path <- args[1]
print(data_path)
tempreads <- read_delim(data_path)
tempreads %>%
    dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0)) %>%
    filter(mod_code == "m") %>%
    group_by(read_id) %>%
    mutate(n = n()) %>%
    dplyr::relocate(n) %>%
    dplyr::select(read_id, ref_position, chrom, mod_indicator) %>%
    write_delim(gsub(".temp.tsv", "", data_path), delim = "\t")
