library(readr)
library(dplyr)
sample_table <- read_csv(conf[["sample_table"]])
sample_table <- sample_table[match(conf$samples, sample_table$sample_name), ]
# set any factors here
# e.g. sample_table <- sample_table %>% mutate(braak = as.factor(braak, levels = c("O", "I", "II", "III", "IV", "V","V/VI", "VI", ordered = TRUE)))
