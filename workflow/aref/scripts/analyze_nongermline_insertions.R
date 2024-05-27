module_name <- "aref"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")

library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)
library(magrittr)
library(forcats)
library(jsonlite)
library(ggpmisc)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            tldroutput = sprintf("aref/%s_tldr/%s.table.txt", conf$samples, conf$samples),
            json = sprintf("aref/qc/%s/%spycoQC.json", conf$samples, conf$samples)
        ), env = globalenv())
        assign("outputs", list(
            plots = "aref/results/somatic_insertions/analyze_nongermline_insertions.rds"
        ), env = globalenv())
    }
)
outputdir <- dirname(outputs$plots)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)


# rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)
# rmfamilies <- read_csv(inputs$r_repeatmasker_annotation, col_names = TRUE)
# rmann <- left_join(rmfragments, rmfamilies)
dflist <- list()
rm(sample_sequencing_data)
for (sample in sample_table$sample_name) {
    if (conf$update_ref_with_tldr$per_sample == "yes") {
    df <- read.table(grep(sprintf("%s_tldr", sample), inputs$tldroutput, value = TRUE), header = TRUE) %>%
        mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>% tibble()
    } else {
    df <- read.table(inputs$tldroutput, header = TRUE) %>%
        mutate(faName = paste0("nonrefins_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>% tibble()
    }

    df$sample_name <- sample
    dflist[[sample]] <- df

    json <- fromJSON(grep(sprintf("/%s/", sample), inputs$json, value = TRUE))
    reads_number <- json[["All Reads"]][["basecall"]]$reads_number
    N50 <- json[["All Reads"]][["basecall"]]$N50
    bases_number <- json[["All Reads"]][["basecall"]]$bases_number
    row <- tibble(sample_name = sample, reads_number = reads_number, N50 = N50, bases_number = bases_number)
    if (!exists("sample_sequencing_data")) {
        sample_sequencing_data <- row
    } else {
        sample_sequencing_data <- rbind(sample_sequencing_data, row)
    }
}

dfall <- do.call(bind_rows, dflist) %>% tibble() %>% left_join(sample_sequencing_data) %>% left_join(sample_table)
sample_sequencing_data

dffilt <- dfall %>%  separate_wider_delim(EmptyReads,delim = "|", names = c("bamname", "emptyreadsnum"))  %>% mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>% 
    filter(fraction_reads_count < 0.2) %>% filter(UsedReads <5)


for (sample in unique(dffilt$sample_name)) {
    tempoutputdir <- sprintf("aref/%s_Analysis/tldr_plots/nongermline", sample)
    p <- dfall %>% ggplot(aes(x = UsedReads, fill = Filter == "PASS")) + geom_histogram() + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/usedreads_hist.pdf", tempoutputdir), 5, 4)

    p <- dfall %>% ggplot(aes(x = Family, fill = Filter == "PASS")) + geom_bar() + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/usedreads_bar.pdf", tempoutputdir), 5, 4)

    p <- dffilt %>% filter(UsedReads == 1) %>% ggplot(aes(x = Family, fill = Filter == "PASS")) + geom_bar() + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/single_read_bar.pdf", tempoutputdir), 5, 4)

    p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = Family, fill = Filter == "PASS")) + geom_bar() + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/single_read_fillpass_bar.pdf", tempoutputdir), 5, 4)

    dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% write_delim(sprintf("%s/single_read_pass.tsv", tempoutputdir), delim = "\t")

    p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = Family, fill = is.na(TSD))) + geom_bar() + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", tempoutputdir), 5, 4)

    p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = LengthIns)) + geom_histogram() + labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions") + facet_wrap(~Family) + mtclosed + anchorbar + scale_palette
    mysaveandstore(sprintf("%s/single_read_pass_insertion_length.pdf", tempoutputdir), 5, 3)

}

p <- dffilt %>% filter(UsedReads == 1) %>% ggplot(aes(x = sample_name, fill = Filter == "PASS")) + geom_bar() + facet_wrap(~Family) + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") + mtopen + anchorbar + scale_palette+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_bar.pdf", outputdir), 8, 5)

p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = sample_name, fill = Filter == "PASS")) + geom_bar()+ facet_wrap(~Family) + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtopen + anchorbar + scale_palette+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", outputdir), 8, 5)

dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% write_delim(sprintf("%s/single_read_pass.tsv", outputdir), delim = "\t")

p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = sample_name, fill = condition)) + geom_bar()+ facet_wrap(~Family) + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtopen + anchorbar + scale_conditions + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", outputdir), 8, 5)

p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = LengthIns)) + geom_histogram() + labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + facet_wrap(sample_name~Family) + mtclosed + anchorbar + scale_palette+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_insertion_length.pdf", outputdir), 10, 10)


tryCatch(
    {
        metadata_vars <- colnames(sample_table)[!(colnames(sample_table) %in% c("sample_name", "condition")) ]
        sequencing_metadata_vars <- c("N50", "reads_number", "bases_number")


        grouping_vars <- c("SingleReadSupport", "sample_name", "condition", sequencing_metadata_vars, metadata_vars)
        pf <- dfall %>% mutate(SingleReadSupport = ifelse(UsedReads == 1, "SingleReadSupport", "MultiReadSupport")) %>% group_by(across(grouping_vars)) %>% summarise(nins = n()) %>% ungroup()
        pfsrs <- pf %>% filter(SingleReadSupport == "SingleReadSupport")
        pfpass <- dfall %>% filter(Filter == "PASS") %>% mutate(SingleReadSupport = ifelse(UsedReads == 1, "SingleReadSupport", "MultiReadSupport")) %>% group_by(across(grouping_vars)) %>% summarise(nins = n()) %>% ungroup()
        pfpasssrs <- pfpass %>% filter(SingleReadSupport == "SingleReadSupport")

        p<- pf %>% ggplot(aes(x = N50, y = nins, color = sample_name)) + geom_point() +
            facet_wrap(~SingleReadSupport) + labs(x = "N50", y = "Number of Insertions") + mtclosed + scale_samples_unique
        mysaveandstore(sprintf("%s/n50_vs_insertions_facet.pdf", outputdir), 7, 4)

        p<- pf %>% ggplot(aes(x = reads_number, y = nins, color = sample_name)) + geom_point() +
            facet_wrap(~SingleReadSupport) + labs(x = "Total Reads Count", y = "Number of Insertions") + mtclosed + scale_samples_unique
        mysaveandstore(sprintf("%s/reads_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        p<- pf %>% ggplot(aes(x = bases_number, y = nins, color = sample_name)) + geom_point() +
            facet_wrap(~SingleReadSupport) + labs(x = "Total Bases Count", y = "Number of Insertions") + mtclosed + scale_samples_unique
        mysaveandstore(sprintf("%s/bases_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        predictors <- c(sequencing_metadata_vars, metadata_vars, "condition")
        model_formula <- paste("nins ~ ",paste(predictors, collapse="+"),sep = "")
        sink(sprintf("%s/lm_allins.txt", outputdir))
        model <- lm(as.formula(model_formula), data = pf)
        model %>% summary()
        sink()

        sink(sprintf("%s/lm_srsins.txt", outputdir))
        model <- lm(as.formula(model_formula), data = pfsrs)
        model %>% summary()
        sink()


        for (mvar in metadata_vars) {
            tryCatch({
            if (sample_table[[mvar]] %>% is.numeric()) {
                p <- pf %>% ggplot(aes(x = !!sym(mvar), y = nins, color = sample_name)) + geom_point() + 
                    labs(y = "Number of Insertions", title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") + mtopen +
                    scale_samples_unique + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            } else {
                p <- pf %>% ggplot(aes(x = sample_name, y = nins,fill = condition)) + facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") + geom_col() + 
                    labs(x = "", y = "Number of Insertions",, title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") + mtclosed + anchorbar +
                    scale_conditions + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
            mysaveandstore(sprintf("%s/ins_by_%s.pdf", outputdir, mvar), 5, 4)
            }, error = function(e) {
                print(e)
            })
        }
        for (mvar in metadata_vars) {
            tryCatch({
            if (sample_table[[mvar]] %>% is.numeric()) {
                p <- pfsrs %>% ggplot(aes(x = !!sym(mvar), y = nins, color = sample_name)) + geom_point() + 
                    labs(y = "Number of Insertions", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtopen +
                    scale_samples_unique + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            } else {
                p <- pfsrs %>% ggplot(aes(x = sample_name, y = nins,fill = condition)) + facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") + geom_col() + 
                    labs(x = "", y = "Number of Insertions",, title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtclosed + anchorbar +
                    scale_conditions + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
            mysaveandstore(sprintf("%s/ins_by_%s_singleread.pdf", outputdir, mvar), 5, 4)
            }, error = function(e) {
                print(e)
            })
        }
        for (mvar in metadata_vars) {
            tryCatch({
            if (sample_table[[mvar]] %>% is.numeric()) {
                p <- pfpass %>% ggplot(aes(x = !!sym(mvar), y = nins, color = sample_name)) + geom_point() + 
                    labs(y = "Number of Insertions", title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") + mtopen +
                    scale_samples_unique + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            } else {
                p <- pfpass %>% ggplot(aes(x = sample_name, y = nins,fill = condition)) + facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") + geom_col() + 
                    labs(x = "", y = "Number of Insertions",, title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") + mtclosed + anchorbar +
                    scale_conditions + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
            mysaveandstore(sprintf("%s/ins_by_%s_pass.pdf", outputdir, mvar), 5, 4)
            }, error = function(e) {
                print(e)
            })
        }
        for (mvar in metadata_vars) {
            tryCatch({
            if (sample_table[[mvar]] %>% is.numeric()) {
                p <- pfpasssrs %>% ggplot(aes(x = !!sym(mvar), y = nins, color = sample_name)) + geom_point() + 
                    labs(y = "Number of Insertions", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtopen +
                    scale_samples_unique + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            } else {
                p <- pfpasssrs %>% ggplot(aes(x = sample_name, y = nins,fill = condition)) + facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") + geom_col() + 
                    labs(x = "", y = "Number of Insertions",, title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtclosed + anchorbar +
                    scale_conditions + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
            mysaveandstore(sprintf("%s/ins_by_%s_singleread_pass.pdf", outputdir, mvar), 5, 4)
            }, error = function(e) {
                print(e)
            })
        }

    },
    error = function(e) {
        message("Likely single condition")
    }
)

##########

save(mysaveandstoreplots, file = outputs$plots)
