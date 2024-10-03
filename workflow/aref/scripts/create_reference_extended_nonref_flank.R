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
            updated_reference = "aref/default/A.REF-pre-ins-filtering.fa",
            non_ref_contigs = "aref/default/A.REF-pre-ins-filtering_nonrefcontigs.fa"
        ), env = globalenv())
    }
)

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
    filter(fraction_reads_count > 0.30) %>%
    filter(SpanReads > 3) %>%
    filter(UsedReads > 10) %>%
    mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS") %>%
    filter(!is.na(TSD)) %>%
    filter(MedianMapQ == 60)

df <- df[!(df %$% faName %>% duplicated()), ]

df <- df %>%
    mutate(Consensus_lower = str_extract_all(Consensus, pattern = "[:lower:]+")) %>%
    rowwise() %>%
    mutate(insLenMatch = which(as.vector((nchar(Consensus_lower))) == LengthIns)) %>%
    mutate(ins_consensus_noflank = Consensus_lower[insLenMatch]) %>%
    mutate(ins_consensus_30flank = str_sub(Consensus, str_locate(Consensus, ins_consensus_noflank)[1] - 30, str_locate(Consensus, ins_consensus_noflank)[2] + 30))


# fasta
ss <- DNAStringSet(df %>% dplyr::arrange(faName) %$% ins_consensus_30flank)
names(ss) <- df %>% dplyr::arrange(faName) %$% faName
writeXStringSet(ss, outputs$non_ref_contigs, append = FALSE, format = "fasta")

ss_plusflank <- DNAStringSet(df %>% dplyr::arrange(faName) %$% Consensus)
names(ss_plusflank) <- df %>% dplyr::arrange(faName) %$% faName
writeXStringSet(ss_plusflank, outputs$non_ref_contigs_plusflank, append = FALSE, format = "fasta")

dir.create(dirname(outputs$updated_reference), recursive = TRUE, showWarnings = FALSE)
system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference))
writeXStringSet(ss, outputs$updated_reference, append = TRUE, format = "fasta")
system(sprintf("samtools faidx %s", outputs$updated_reference))

system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference_plusflank))
writeXStringSet(ss_plusflank, outputs$updated_reference_plusflank, append = TRUE, format = "fasta")
system(sprintf("samtools faidx %s", outputs$updated_reference_plusflank))
