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
library(ggpubr)
library(ggh4x)

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




dflist <- list()
for (sample in sample_table$sample_name) {
    if (conf$update_ref_with_tldr$per_sample == "yes") {
        df <- read.table(grep(sprintf("%s_tldr", sample), inputs$tldroutput, value = TRUE), header = TRUE) %>%
            mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
            tibble()
    } else {
        df <- read.table("aref/A.REF_tldr/A.REF.table.txt", header = TRUE) %>%
            mutate(faName = paste0("NI_", Subfamily, "_", Chrom, "_", Start, "_", End)) %>%
            tibble()
    }

    df$sample_name <- sample
    df <- df %>%
        filter(grepl(sample, SampleReads)) %>%
        filter(!grepl(paste(setdiff(sample_table$sample_name, sample), collapse = "|"), SampleReads))
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

dfall <- do.call(bind_rows, dflist) %>%
    tibble() %>%
    left_join(sample_sequencing_data) %>%
    left_join(sample_table)
sample_sequencing_data

dffilt <- dfall %>%
    separate_wider_delim(EmptyReads, delim = "|", names = c("bamname", "emptyreadsnum")) %>%
    mutate(fraction_reads_count = UsedReads / (UsedReads + as.numeric(emptyreadsnum))) %>%
    filter(fraction_reads_count < 0.1) %>%
    filter(MedianMapQ >= 60) %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(!is.na(TSD))

for (element_type in dffilt %$% Subfamily %>% unique()) {
    dftemp <- dffilt %>% filter(Subfamily == element_type)
    p <- dftemp %>%
        gghistogram(x = "LengthIns") +
        mtopen +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/length_distribution.pdf", outputdir, element_type))
    dftemp %$% LengthIns %>% quantile()
    p <- dftemp %>%
        mutate(tsdlen = nchar(TSD)) %>%
        gghistogram(x = "tsdlen", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/tsd_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd5 = nchar(Transduction_5p)) %>%
        gghistogram(x = "trsd5", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/trsd5_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        mutate(trsd3 = nchar(Transduction_3p)) %>%
        gghistogram(x = "trsd3", fill = "blue") +
        labs(title = element_type)
    mtopen
    mysaveandstore(sprintf("%s/%s/trsd3_length_distribution.pdf", outputdir, element_type))

    p <- dftemp %>%
        arrange(-EndTE) %>%
        mutate(nrow = row_number()) %>%
        ggplot() +
        geom_segment(aes(x = StartTE, xend = EndTE, y = nrow, yend = nrow)) +
        mtclosed +
        labs(title = element_type)
    mysaveandstore(sprintf("%s/%s/te_body_distribution.pdf", outputdir, element_type))
}


for (sample in unique(dffilt$sample_name)) {
    tempoutputdir <- sprintf("aref/%s_Analysis/tldr_plots/nongermline", sample)
    dfallsample <- dfall %>% filter(sample_name == sample)
    dffiltsample <- dffilt %>% filter(sample_name == sample)

    p <- dfallsample %>% ggplot(aes(x = UsedReads, fill = Filter == "PASS")) +
        geom_histogram() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/usedreads_hist.pdf", tempoutputdir), 5, 4)

    p <- dfallsample %>% ggplot(aes(x = Family, fill = Filter == "PASS")) +
        geom_bar() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/usedreads_bar.pdf", tempoutputdir), 5, 4)

    p <- dffiltsample %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        ggplot(aes(x = Family, fill = Filter == "PASS")) +
        geom_bar() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/single_read_bar.pdf", tempoutputdir), 5, 4)

    p <- dffiltsample %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(MedianMapQ >= 60) %>%
        ggplot(aes(x = Family, fill = Filter == "PASS")) +
        geom_bar() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/single_read_fillpass_bar.pdf", tempoutputdir), 5, 4)

    dffiltsample %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(MedianMapQ >= 60) %>%
        write_delim(sprintf("%s/single_read_pass.tsv", tempoutputdir), delim = "\t")

    p <- dffiltsample %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(MedianMapQ >= 60) %>%
        ggplot(aes(x = Family, fill = is.na(TSD))) +
        geom_bar() +
        labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
        mtopen +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", tempoutputdir), 5, 4)

    p <- dffiltsample %>%
        filter(UsedReads == 1) %>%
        filter(SpanReads == 1) %>%
        filter(Filter == "PASS") %>%
        filter(MedianMapQ >= 60) %>%
        ggplot(aes(x = LengthIns)) +
        geom_histogram() +
        labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions") +
        facet_wrap(~Family) +
        mtclosed +
        anchorbar +
        scale_palette
    mysaveandstore(sprintf("%s/single_read_pass_insertion_length.pdf", tempoutputdir), 5, 3)
}

p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    ggplot(aes(x = sample_name, fill = Filter == "PASS")) +
    geom_bar() +
    facet_wrap(~Family) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions") +
    mtopen +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    mtclosed
mysaveandstore(sprintf("%s/single_read_bar.pdf", outputdir), 8, 5)

# p <- dffilt %>% filter(UsedReads == 1) %>% filter(Filter == "PASS") %>% ggplot(aes(x = sample_name, fill = Filter == "PASS")) + geom_bar()+ facet_wrap(~Family) + labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") + mtopen + anchorbar + scale_palette+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", outputdir), 8, 5)

dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    write_delim(sprintf("%s/single_read_pass.tsv", outputdir), delim = "\t")

p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    ggplot(aes(x = sample_name, fill = condition)) +
    geom_bar() +
    facet_wrap(~Family) +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_bar.pdf", outputdir), 8, 5)

p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    group_by(sample_name, Subfamily, condition) %>%
    summarise(nins = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Subfamily, values_from = nins) %>%
    ggplot(aes(x = L1HS, y = AluY)) +
    geom_point(aes(color = condition)) +
    scale_conditions +
    stat_cor() +
    stat_smooth(method = "lm", formula = y ~ x, geom = "smooth") +
    labs(x = "L1HS", y = "AluY", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_l1hs_aluy_corr_bar.pdf", outputdir), 4, 4)

dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    filter(Subfamily %in% c("AluY")) %>%
    filter(sample_name == "AD1") %>%
    filter(UUID == "a72cba1d-f8c4-4dd4-8054-ee7c97a6c704") %>%
    print(width = Inf)

p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    filter(Subfamily %in% c("L1HS", "AluY")) %>%
    ggplot(aes(x = TEMatch, fill = condition)) +
    geom_histogram() +
    facet_grid2(rows = vars(sample_name), cols = vars(Subfamily), scale = "free_x", axes = "all", remove_labels = "y") +
    labs(x = "TE Match", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_histogram.pdf", outputdir), 5, 20)

p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    filter(Subfamily %in% c("AluY")) %>%
    ggplot(aes(x = sample_name, fill = condition)) +
    geom_bar() +
    labs(x = "Supporting Read Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    mtclosed +
    anchorbar +
    scale_conditions +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_bar_aluy.pdf", outputdir), 8, 5)



p <- dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    ggplot(aes(x = LengthIns)) +
    geom_histogram() +
    labs(x = "Insertion Length", y = "Count", title = "RTE Somatic Insertions", subtitle = "Single Read Supported") +
    facet_grid(rows = vars(sample_name), cols = vars(Family), scales = "free") +
    mtclosed +
    anchorbar +
    scale_palette +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(sprintf("%s/single_read_pass_insertion_length.pdf", outputdir), 10, 10)

dffilt %>%
    filter(UsedReads == 1) %>%
    filter(SpanReads == 1) %>%
    filter(Filter == "PASS") %>%
    filter(MedianMapQ >= 60) %>%
    filter(LengthIns > 5000) %>%
    filter(Subfamily == "L1HS") %>%
    print(width = Inf)

tryCatch(
    {
        metadata_vars <- colnames(sample_table)[!(colnames(sample_table) %in% c("sample_name", "condition"))]
        sequencing_metadata_vars <- c("N50", "reads_number", "bases_number")


        grouping_vars <- c("sample_name", "condition", sequencing_metadata_vars, metadata_vars)
        pfpass <- dfall %>%
            filter(Filter == "PASS") %>%
            filter(UsedReads > 4) %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(SingleReadSupport = "MultiReadSupport")
        pfpasssrs <- dfall %>%
            filter(Filter == "PASS") %>%
            filter(UsedReads == 1) %>%
            filter(SpanReads == 1) %>%
            filter(MedianMapQ >= 60) %>%
            filter(fraction_reads_count < 0.1) %>%
            group_by(across(grouping_vars)) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(SingleReadSupport = "SingleReadSupport")
        pf <- bind_rows(pfpass, pfpasssrs)
        p <- pf %>% ggplot(aes(x = N50, y = nins, color = sample_name)) +
            geom_point(size = 2) +
            facet_wrap(~SingleReadSupport) +
            labs(x = "N50", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/n50_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- pf %>% ggplot(aes(x = reads_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~SingleReadSupport) +
            labs(x = "Total Reads Count", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/reads_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        p <- pf %>% ggplot(aes(x = bases_number, y = nins, color = sample_name)) +
            geom_point() +
            facet_wrap(~SingleReadSupport) +
            labs(x = "Total Bases Count", y = "Number of Insertions") +
            mtclosed +
            scale_samples_unique
        mysaveandstore(sprintf("%s/bases_number_vs_insertions_facet.pdf", outputdir), 7, 4)

        predictors <- c(sequencing_metadata_vars, metadata_vars, "condition")
        model_formula <- paste("nins ~ ", paste(predictors, collapse = "+"), sep = "")

        model <- lm(as.formula(model_formula), data = pf)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_allins2.txt", outputdir), delim = "\t")

        model <- lm(as.formula(model_formula), data = pfpasssrs)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_srsins.txt", outputdir), delim = "\t")

        predictors <- c(sequencing_metadata_vars, "condition")
        model_formula <- paste("nins ~ ", paste(predictors, collapse = "+"), sep = "")
        model <- lm(as.formula(model_formula), data = pfpasssrs)
        output <- coef(summary(model)) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            tibble()
        write_delim(output, file = sprintf("%s/lm_srsins_model_limited.txt", outputdir), delim = "\t")


        pfpass <- dfall %>%
            filter(Filter == "PASS") %>%
            filter(UsedReads > 4) %>%
            group_by(across(c(grouping_vars, Subfamily))) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(SingleReadSupport = "MultiReadSupport")
        pfpasssrs <- dfall %>%
            filter(Filter == "PASS") %>%
            filter(UsedReads == 1) %>%
            filter(SpanReads == 1) %>%
            filter(MedianMapQ >= 60) %>%
            filter(fraction_reads_count < 0.1) %>%
            group_by(across(c(grouping_vars, Subfamily))) %>%
            summarise(nins = n()) %>%
            ungroup() %>%
            mutate(SingleReadSupport = "SingleReadSupport")
        pf <- bind_rows(pfpass, pfpasssrs)
        for (insertion_type in pfpasssrs %$% Subfamily %>% unique()) {
            pf <- pfpasssrs %>% filter(Subfamily == insertion_type)
            for (mvar in metadata_vars) {
                tryCatch(
                    {
                        if (sample_table[[mvar]] %>% is.numeric()) {
                            p <- pf %>% ggscatter(x = mvar, y = "nins", color = "condition", size = 3) +
                                scale_conditions + stat_cor() +
                                stat_smooth(method = "lm", formula = y ~ x, geom = "smooth")
                        } else {
                            p <- pf %>% ggplot(aes(x = sample_name, y = nins, fill = condition)) +
                                facet_grid(cols = vars(!!sym(mvar)), scales = "free_x", space = "free") +
                                geom_col() +
                                labs(x = "", y = "Number of Insertions", , title = "RTE Somatic Insertions", subtitle = "Multi and Single Read Supported") +
                                mtclosed +
                                anchorbar +
                                scale_conditions +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        }
                        mysaveandstore(sprintf("%s/%s_by_%s_somatic_pass.pdf", outputdir, insertion_type, mvar), 5, 4)
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        }
    },
    error = function(e) {
        message("Likely single condition")
    }
)

##########

save(mysaveandstoreplots, file = outputs$plots)
