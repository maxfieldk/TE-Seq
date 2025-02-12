module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
library(configr)
library(Biostrings)
set.seed(123)

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
            updated_reference_plusflank = "aref/extended/A.REF-pre-ins-filtering.fa",
            non_ref_contigs_plusflank = "aref/extended/A.REF-pre-ins-filtering_nonrefcontigs.fa"
        ), env = globalenv())
    }
)

tldr_output_dir <- gsub(".table.txt", "", inputs$tldroutput)
df <- read_delim(inputs$tldroutput)
er <- df %$% EmptyReads %>% str_extract_all("\\|[0-9]+")

emptyreadsnum <- c()
for (i in 1:length(er)) {
    if (length(er[[i]]) == 1) {
        val <- er[[i]] %>%
            gsub("\\|", "", .) %>%
            as.numeric()
    } else {
        val <- er[[i]] %>%
            gsub("\\|", "", .) %>%
            as.numeric() %>%
            sum()
    }
    emptyreadsnum <- c(emptyreadsnum, val)
}

df$emptyreadsnum <- emptyreadsnum %>% replace_na(0)
df <- df %>%
    mutate(fraction_reads_count = UsedReads / (UsedReads + emptyreadsnum)) %>%
    filter(fraction_reads_count > 0.25) %>%
    filter(SpanReads > 2) %>%
    filter(UsedReads > 5) %>%
    mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS") %>%
    filter(!is.na(TSD)) %>%
    filter(MedianMapQ == 60)

df <- df[!(df %$% faName %>% duplicated()), ]

# fasta
ss_plusflank <- DNAStringSet()
for (i in 1:nrow(df)) {
    UUID <- df[i, ]$UUID
    ss <- readDNAStringSet(sprintf("%s/%s.cons.ref.fa", tldr_output_dir, UUID))
    names(ss) <- df[i, ]$faName
    ss_plusflank <- c(ss_plusflank, ss)
}

dir.create(dirname(outputs$non_ref_contigs_plusflank), recursive = TRUE)
writeXStringSet(ss_plusflank, outputs$non_ref_contigs_plusflank, append = FALSE, format = "fasta")

system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference_plusflank))
writeXStringSet(ss_plusflank, outputs$updated_reference_plusflank, append = TRUE, format = "fasta")
system(sprintf("samtools faidx %s", outputs$updated_reference_plusflank))
