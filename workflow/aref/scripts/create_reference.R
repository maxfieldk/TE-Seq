module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
library(configr)
library(Biostrings)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = "aref/A.REF_tldr/A.REF.table.txt",
            reference = "../genome_files/reference.ucsc.fa"
        ), env = globalenv())
        assign("outputs", list(
            updated_reference = "aref/A.REF-pre-ins-filtering.fa",
            non_ref_contigs = "aref/A.REF-pre-ins-filtering_nonrefcontigs.fa"
        ), env = globalenv())
    }
)

df <- read_delim(inputs$tldroutput) 
er <- df %$% EmptyReads %>% str_extract_all("\\|[0-9]+")

emptyreadsnum <- c()
for (i in 1:length(er)) {
    if (length(er[[i]]) == 1) {
        val <- er[[i]] %>% gsub("\\|", "", .) %>% as.numeric()
    } else {
        val <- er[[i]] %>% gsub("\\|", "", .) %>% as.numeric() %>% sum()
    }
    emptyreadsnum <- c(emptyreadsnum, val)
}

df$emptyreadsnum <- emptyreadsnum
df <- df %>% mutate(fraction_reads_count = UsedReads / (UsedReads + emptyreadsnum)) %>% 
    filter(fraction_reads_count > 0.20) %>%
    mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS")

df <- df[!(df %$% faName %>% duplicated()), ]

# fasta
ss <- DNAStringSet(df %>% dplyr::arrange(faName) %$% Consensus)
names(ss) <- df %>% dplyr::arrange(faName) %$% faName
writeXStringSet(ss, outputs$non_ref_contigs, append = FALSE, format = "fasta")

dir.create(dirname(outputs$updated_reference), recursive = TRUE, showWarnings = FALSE)
system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference))
writeXStringSet(ss, outputs$updated_reference, append = TRUE, format = "fasta")
system(sprintf("samtools faidx %s", outputs$updated_reference))
