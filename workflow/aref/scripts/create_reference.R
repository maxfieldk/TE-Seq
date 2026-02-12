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
            tldroutput = "aref/AD5_tldr/AD5.table.txt",
            reference = "../genome_files/reference.ucsc.fa"
        ), env = globalenv())
        assign("params", list(
            sample_or_ref = "AD5"
        ), env = globalenv())
        assign("outputs", list(
            updated_reference = "aref/default/AD5-pre-ins-filtering.fa",
            non_ref_contigs = "aref/default/AD5-pre-ins-filtering_nonrefcontigs.fa"
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
    filter(fraction_reads_count > 0.25) %>%
    filter(SpanReads > 2) %>%
    filter(UsedReads > 5) %>%
    mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
    filter(!is.na(Family)) %>%
    filter(!is.na(StartTE)) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ == 60)

if (conf$update_ref_with_tldr$require_tsd == "yes") {
    df <- df %>% filter(!is.na(TSD))
}

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

dir.create(dirname(outputs$updated_reference), recursive = TRUE, showWarnings = FALSE)
system(sprintf("cp %s %s", inputs$reference, outputs$updated_reference))
writeXStringSet(ss, outputs$updated_reference, append = TRUE, format = "fasta")

conda_base_path <- system("conda info --base", intern = TRUE)
system(
    sprintf(
        "bash -c 'source %s/etc/profile.d/conda.sh && conda activate omics && samtools faidx %s'",
        conda_base_path, outputs$updated_reference
    )
)
