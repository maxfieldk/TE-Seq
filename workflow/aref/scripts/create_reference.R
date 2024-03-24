source("~/data/common/myDefaults.r")
library(configr)
library(Biostrings)

conf <- configr::read.config(file = "conf/config.yaml")[["aref"]]
tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = "aref/tldr/tldr.table.txt",
            reference = "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa"
        ), env = globalenv())
        assign("outputs", list(updated_reference = "updated_reference/hs1_lf1.fa"), env = globalenv())
    }
)

df <- read_delim(inputs$tldroutput) %>%
    mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS")

df <- df[!(df %$% faName %>% duplicated()), ]

# fasta
ss <- DNAStringSet(df %>% dplyr::arrange(faName) %$% Consensus)
names(ss) <- df %>% dplyr::arrange(faName) %$% faName
# writeXStringSet(ss, paste0(dirname(outputs$updated_reference), "/nonrefcontigs.fa"), append = FALSE, format = "fasta")

dir.create(dirname(outputs$updated_reference), recursive = TRUE, showWarnings = FALSE)
system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference))
writeXStringSet(ss, outputs$updated_reference, append = TRUE, format = "fasta")
system(sprintf("samtools faidx %s", outputs$updated_reference))
