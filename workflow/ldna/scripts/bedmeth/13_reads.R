module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
source("conf/sample_table_source.R")
sample_table <- sample_table %>%
    mutate(sample_name = factor(sample_name, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels)) %>%
    mutate(sample = sample_name)
samples <- conf$samples
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

cell_types_cols <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")

cell_fractions <- read_delim("ldna/results/m/tables/scMD_cell_type_fractions_all.csv")
cell_fractions_scaled <- cell_fractions %>%
    mutate(Neuron = Inh + Exc) %>%
    mutate(across(all_of(c(cell_types_cols, "Neuron")), ~ as.numeric(scale(.)), .names = "{.col}_z"))

sample_table <- sample_table %>% left_join(cell_fractions_scaled)

set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(zoo)
library(pryr)
library(circlize)
library(rGREAT)
library(reactome.db)
library(msigdb)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)
# library(ReMapEnrich)
library(msigdbr)
library(Biostrings)
library(ggpubr)
library(PCAtools)
library(patchwork)
library(ggh4x)
library(ggnewscale)
library(betareg)
library(scales)
library(ggnewscale)
library(glmmTMB)
library(broom.mixed)

conditions <- conf$levels
contrasts <- conf$contrasts

# Helper to parse a contrast string into condition1 (reference) and condition2 (test)
parse_contrast <- function(contrast_string) {
    parts <- str_match(contrast_string, "^condition_(.+)_vs_(.+)$")
    list(condition2 = parts[1, 2], condition1 = parts[1, 3])
}

# For backwards compatibility, set condition1/condition2 from first contrast
first_contrast <- parse_contrast(contrasts[[1]])
condition1 <- first_contrast$condition1
condition2 <- first_contrast$condition2
condition1samples <- sample_table[sample_table$condition == condition1, ]$sample_name
condition2samples <- sample_table[sample_table$condition == condition2, ]$sample_name
enough_samples_per_condition_for_stats <- ifelse(length(condition1samples) < 3 | length(condition2samples) < 3, FALSE, TRUE)

adjustment_set <- c(conf$linear_model_adjustment_set)
lm_right_hand_side <- ifelse(is.null(adjustment_set), "condition", paste0(c("condition", adjustment_set), collapse = " + "))
adjustment_set_categorical <- if (is.null(adjustment_set)) {
    NULL
} else {
    sample_table %>%
        dplyr::select(all_of(adjustment_set)) %>%
        map(class) %>%
        unlist() %>%
        grep(pattern = "character|factor", value = TRUE) %>%
        names()
}

asc <- ifelse(is.null(adjustment_set_categorical), "", adjustment_set_categorical[1])


{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")
    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
    if ("chrY" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, grep("_chrX_|_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes)
        } else {
            CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, grep("_chrY_", nonrefchromosomes, invert = TRUE, value = TRUE))
            CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX)
        }
    } else if ("chrX" %in% conf$SEX_CHROMOSOMES_NOT_INCLUDED_IN_ANALYSIS) {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrY, grep("_chrX_", nonrefchromosomes, invert = TRUE, value = TRUE))
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrY)
    } else {
        CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
        CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
    }
}

#################### functions and themes
named_group_split <- function(.tbl, ...) {
    grouped <- group_by(.tbl, ...)
    names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

    grouped %>%
        group_split() %>%
        rlang::set_names(names)
}
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethlpaths = sprintf("ldna/intermediates/%s/methylation/%s_CG_bedMethyl.bed", sample_table$sample_name, sample_table$sample_name),
            data = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = "ldna/results/m/tables/dmrs.tsv",
            dmls = "ldna/results/m/tables/dmls.tsv"
        ), env = globalenv())
        assign("params", list(
            contrasts = conf$contrasts,
            mod_code = "m"
        ), env = globalenv())
        assign("outputs", list(
            promoters_bed = "ldna/Rintermediates/m/promoters_t05.bed",
            dmrpromoterhyper_bed = "ldna/Rintermediates/m/promoters_dmhyperregions_t05.bed",
            dmrpromoterhypo_bed = "ldna/Rintermediates/m/promoters_dmhyporegions_t05.bed"
        ), env = globalenv())
    }
)


merge_with_grs <- function(grs, rte_frame) {
    mbo <- mergeByOverlaps(grs, rte_frame)
    methdf <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    rte_only_frame <- mbo$rte_frame %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
    rtedf_promoters <- bind_cols(methdf, rte_only_frame)
    return(rtedf_promoters)
}

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

# rmannextended <- get_repeat_annotations(
#     default_or_extended = "default",
#     keep_non_central = FALSE
# )

# rmannextended %>% filter(rte_subfamily == "L1HS") %>% filter(refstatus == "Ref") %>% filter(intactness_req == "Intact")
rmannextended <- get_repeat_annotations(
    default_or_extended = "extended",
    keep_non_central = FALSE
)





##################


readscg <- read_delim(sprintf("ldna/Rintermediates/%s/reads_context_cpg.tsv", params$mod_code), col_names = TRUE) %>%
    mutate(sample = factor(sample, levels = sample_table$sample_name)) %>%
    mutate(condition = factor(condition, levels = conf$levels))



readscg_endfilt <- readscg %>%
    group_by(read_id) %>%
    mutate(read_n_total_mod = n()) %>%
    arrange(desc(forward_read_position), .by_group = TRUE) %>%
    group_split() %>%
    map_dfr(function(df) {
        is_mod <- df$mod_qual > 0.5
        mod_ix <- which(is_mod)

        if (length(mod_ix) < 2) {
            return(tibble())
        }

        between_unmod <- (mod_ix[2] - mod_ix[1]) - 1

        if (between_unmod < 4) {
            return(df[mod_ix[2]:nrow(df), ])
        } else if (length(mod_ix) >= 3) {
            return(df[mod_ix[3]:nrow(df), ])
        } else {
            return(tibble())
        }
    })

readscg_5endfilt <- readscg %>%
    group_by(read_id) %>%
    mutate(read_n_total_mod = n()) %>%
    arrange(forward_read_position, .by_group = TRUE) %>%
    group_split() %>%
    map_dfr(function(df) {
        is_mod <- df$mod_qual > 0.5
        mod_ix <- which(is_mod)

        if (length(mod_ix) < 2) {
            return(tibble())
        }

        between_unmod <- (mod_ix[2] - mod_ix[1]) - 1

        if (between_unmod < 4) {
            return(df[mod_ix[2]:nrow(df), ])
        } else if (length(mod_ix) >= 3) {
            return(df[mod_ix[3]:nrow(df), ])
        } else {
            return(tibble())
        }
    })


readscg_endfiltnew <- readscg %>%
    group_by(read_id) %>%
    mutate(read_n_total_mod = n()) %>%
    group_by(read_id, start) %>%
    mutate(mod_score = sum(mod_qual)) %>%
    group_by(read_id) %>%
    arrange(desc(forward_read_position), .by_group = TRUE) %>%
    group_split() %>%
    map_dfr(function(df) {
        mod_score <- df %$% mod_score

        is_mod <- mod_score > 0.5
        mod_ix <- which(is_mod)

        if (length(mod_ix) < 1) {
            return(tibble())
        } else {
            return(df[mod_ix[1]:nrow(df), ])
        }
    })

write_csv(readscg_endfiltnew, sprintf("ldna/Rintermediates/%s/reads_context_cpg_endfiltnew.tsv", params$mod_code))
# readscg_endfiltnew <- read_csv(sprintf("ldna/Rintermediates/%s/reads_context_cpg_endfiltnew.tsv", params$mod_code))



readscg_endfiltnew5 <- readscg %>%
    group_by(read_id) %>%
    mutate(read_n_total_mod = n()) %>%
    group_by(read_id, start) %>%
    mutate(mod_score = sum(mod_qual)) %>%
    group_by(read_id) %>%
    arrange(forward_read_position, .by_group = TRUE) %>%
    group_split() %>%
    map_dfr(function(df) {
        mod_score <- df %$% mod_score

        is_mod <- mod_score > 0.5
        mod_ix <- which(is_mod)

        if (length(mod_ix) < 1) {
            return(tibble())
        } else {
            return(df[mod_ix[1]:nrow(df), ])
        }
    })



{ # Get readIDs to filter out from bam for visualization purpose.
    # Will only remove reads that lose CpGs to filtering in the promoter

    flyngl1 <- rmannextended %>%
        filter(rte_length_req == "FL") %>%
        filter(rte_subfamily == "L1HS" | rte_subfamily == "L1PA2")

    flyngl1grs <- GRanges(flyngl1)
    flyngl1grspromoter <- promoters(flyngl1grs, upstream = 0, downstream = 909)

    unfilt_l1_reads <- readscg %>%
        filter(mod_code == params$mod_code) %>%
        GRanges() %>%
        subsetByOverlaps(flyngl1grspromoter, ignore.strand = TRUE) %>%
        as.data.frame() %>%
        tibble()
    filt_l1_reads <- readscg_endfilt %>%
        filter(mod_code == params$mod_code) %>%
        GRanges() %>%
        subsetByOverlaps(flyngl1grspromoter, ignore.strand = TRUE) %>%
        as.data.frame() %>%
        tibble()
    nrow(filt_l1_reads) / nrow(unfilt_l1_reads)
    l1_reads_ncpgs <- full_join(unfilt_l1_reads %>% group_by(read_id) %>% summarise(n_unfilt = n()), filt_l1_reads %>% group_by(read_id) %>% summarise(n_filt = n())) %>%
        mutate(dif = n_unfilt - n_filt)

    # readscg %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")
    # unfilt_l1_reads %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")
    # filt_l1_reads %>% filter(read_id == "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9")

    # readscg %>% filter(gene_id == "L1HS_4q28.3_9") %>% filter(sample == "AD1") %$% read_id %>% unique() %in% c(unfilt_l1_reads %$% read_id)

    # readscg %>% filter(gene_id == "L1HS_4q28.3_9") %>% filter(sample == "AD1") %$% read_id %>% unique() %in% c(filt_l1_reads %$% read_id)

    totally_lost_reads <- l1_reads_ncpgs %>% filter(is.na(dif))
    too_many_cpgs_potentially_compromised <- l1_reads_ncpgs %>%
        filter(!is.na(dif)) %>%
        filter(dif > 5)
    read_to_filter_out <- c(totally_lost_reads$read_id, too_many_cpgs_potentially_compromised$read_id)
    # "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9" %in% read_to_filter_out
    # "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9" %in% c(too_many_cpgs_potentially_compromised %$% read_id)
    # "d5021e72-c8d5-4e94-9bc2-28c05f4d35d9" %in% c(totally_lost_reads %$% read_id)

    write_lines(read_to_filter_out, "ldna/Rintermediates/yng_l1_promoter_potentially_compromised_reads.txt")
    l1_reads_ncpgs %>% pl()
}



########### READ ANALYSES
outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)
cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()


# //ANCHOR - READ ANALYSIS 2


# read_analysis2(readscg_endfiltnew %>% filter(mod_code == params$mod_code), cg_indices)
inputreaddf <- readscg_endfiltnew %>% filter(mod_code == params$mod_code)
mod_code_var <- params$mod_code
regions_of_interest <- list(c(0, 328), c(0, 500), c(0, 909))
required_fraction_of_total_cg <- 0.75
meth_bins <- c(0.5, 0.75)
context <- "CpG"
region <- "L1HS_FL"


read_analysis2 <- function(
    inputreaddf,
    cg_indices,
    mod_code_var = params$mod_code,
    regions_of_interest = list(c(0, 328), c(0, 500), c(0, 909), c(400, 600)),
    required_fraction_of_total_cg = 0.75,
    meth_bins = c(0.5, 0.75), # meth_bins = c(0.25, 0.5, 0.75) #need the 0,.5 bin for locus demeth analysis
    context = "CpG") {
    region <- "L1HS_FL"
    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new_withendfilter/%s_%s", mod_code_var, region, required_fraction_of_total_cg)
    dir.create(outputdirtables, recursive = TRUE)

    breaks <- c(0, meth_bins, 1)
    labels <- paste0("[", head(breaks, -1), ", ", tail(breaks, -1), ")")
    labels[length(labels)] <- sub("\\)$", "]", labels[length(labels)])


    readsdf1 <- inputreaddf %>%
        left_join(rmannextended %>%
            dplyr::select(gene_id, start, end, strand, rte_length_req, intactness_req) %>%
            dplyr::rename(element_strand = strand, element_start = start, element_end = end)) %>%
        filter(rte_length_req == "FL")


    by_cpg_l <- list()
    by_read_l <- list()
    by_sample_l <- list()
    by_gene_id_l <- list()


    for (region_of_interest in regions_of_interest) {
        roistart <- region_of_interest[1]
        roiend <- region_of_interest[2]
        roistring <- paste0(roistart, "to", roiend)

        numCGneeded <- ceiling(length(cg_indices[(cg_indices <= roiend) & (cg_indices >= roistart)]) * required_fraction_of_total_cg)

        utr1 <- readsdf1 %>%
            filter(mod_code == mod_code_var) %>%
            filter(case_when(
                element_strand == "+" ~ (start > element_start + roistart) & (start < element_start + roiend),
                element_strand == "-" ~ (start > element_end - roiend) & (start < element_end - roistart)
            )) %>%
            dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

        by_cpg_temp <- utr1 %>%
            group_by(gene_id, read_id, condition, sample) %>%
            mutate(num_cpgs_in_read = n()) %>%
            mutate(fraction_meth = mean(mod_indicator)) %>%
            relocate(gene_id) %>%
            ungroup() # %>%
        # filter(read_length - forward_read_position > 2000)

        by_read_temp <- by_cpg_temp %>%
            filter(num_cpgs_in_read >= numCGneeded) %>%
            group_by(read_id, gene_id, sample, condition, region) %>%
            summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read), strand = dplyr::first(element_strand), numCGneeded = dplyr::first(numCGneeded)) %>%
            ungroup()

        by_cpg_temp$subset <- as.character(roistring)
        by_read_temp$subset <- as.character(roistring)

        by_cpg_l[[as.character(roistring)]] <- by_cpg_temp
        by_read_l[[as.character(roistring)]] <- by_read_temp



        by_sample_temp <- by_read_temp %>%
            mutate(meth_bin = cut(fraction_meth, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)) %>%
            group_by(sample, region, condition, meth_bin) %>%
            summarise(n = n(), .groups = "drop_last") %>%
            mutate(prop_in_bin = n / sum(n)) %>%
            ungroup()

        by_sample_temp$subset <- as.character(roistring)
        by_sample_temp <- by_sample_temp %>% mutate(subset_bin = paste0(subset, "_", meth_bin))

        by_sample_l[[as.character(roistring)]] <- by_sample_temp


        by_gene_id_temp <- by_read_temp %>%
            mutate(meth_bin = cut(fraction_meth, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)) %>%
            group_by(sample, gene_id, region, meth_bin) %>%
            summarise(n = n(), .groups = "drop_last") %>%
            mutate(prop_in_bin = n / sum(n)) %>%
            ungroup() %>%
            complete(sample, gene_id, region, meth_bin, fill = list(n = 0, prop_in_bin = 0)) %>%
            group_by(sample, gene_id) %>%
            mutate(sample_pres = case_when(mean(n) > 0 ~ 1, TRUE ~ 0)) %>%
            group_by(gene_id) %>%
            mutate(group_size = sum(sample_pres) / (length(breaks) - 1)) %>%
            ungroup() %>%
            left_join(sample_table)
        by_gene_id_temp$subset <- as.character(roistring)
        by_gene_id_temp <- by_gene_id_temp %>% mutate(subset_bin = paste0(subset, "_", meth_bin))

        by_gene_id_l[[as.character(roistring)]] <- by_gene_id_temp
    }

    by_cpg <- purrr::reduce(by_cpg_l, bind_rows)
    by_read <- purrr::reduce(by_read_l, bind_rows)
    by_gene_id <- purrr::reduce(by_gene_id_l, bind_rows)
    by_sample <- purrr::reduce(by_sample_l, bind_rows)

    dir.create(dirname(sprintf("ldna/Rintermediates/%s/%s/cut_highly_demethylated_reads_by_cpg.csv", mod_code_var, region)), recursive = TRUE)
    by_cpg %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/cut_highly_demethylated_reads_by_cpg.csv", mod_code_var, region))
    by_read %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/cut_highly_demethylated_reads_by_read.csv", mod_code_var, region))
    by_gene_id %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/cut_highly_demethylated_reads_by_gene_id.csv", mod_code_var, region))
    by_sample %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/cut_highly_demethylated_reads_by_sample.csv", mod_code_var, region))

    # by_cpg <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", mod_code_var, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_read <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", mod_code_var, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_gene_id <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", mod_code_var, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_sample <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", mod_code_var, region)) %>%
    #         mutate(condition = factor(condition, levels = conf$levels)) %>%
    #         mutate(sample = factor(sample, levels = conf$samples))

    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        # stats <- by_sample %>%
        #     filter(subset != "400to600") %>%
        #     ungroup() %>%
        #     left_join(sample_table) %>%
        #     group_split(subset_bin) %>%
        #     set_names(unique(by_sample %>% filter(subset != "400to600") %$% subset_bin)) %>%
        #     map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "prop_in_bin", lm_right_hand_side)), .x)))) %>%
        #     imap_dfr(~ .x %>% mutate(subset = .y))
        stats_list <- list()
        i <- 1
        for (subset in unique(by_read$subset)) {
            for (bin in unique(by_sample$meth_bin)) {
                by_read_tmp <- by_read %>%
                    filter(subset == !!subset) %>%
                    mutate(meth_bin = cut(fraction_meth, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)) %>%
                    mutate(unmeth = as.integer(ifelse(meth_bin != bin, 0, 1))) %>%
                    dplyr::rename(sample_name = sample) %>%
                    group_by(sample_name, condition, gene_id) %>%
                    summarise(unmeth = sum(unmeth), total = n()) %>%
                    ungroup() %>%
                    left_join(sample_table) %>%
                    mutate(age_z = as.numeric(scale(age))) %>%
                    mutate(condition = factor(condition, levels = conf$levels))
                model_tmp <- glmmTMB(
                    cbind(unmeth, total - unmeth) ~
                        condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name) + (1 | gene_id),
                    data = by_read_tmp,
                    family = binomial()
                )
                library(broom.mixed)
                res_tmp <- broom::tidy(model_tmp) %>%
                    mutate(bin = !!bin) %>%
                    mutate(subset = !!subset)
                stats_list[[i]] <- res_tmp
                i <- i + 1
            }
        }
        stats <- purrr::reduce(stats_list, bind_rows)
        stats_padj <- stats %>%
            filter(grepl("condition", term)) %>%
            filter(subset != "400to600") %>%
            mutate(padj = p.adjust(p.value, method = "fdr"))
        statswpadj <- stats %>% left_join(stats_padj)
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = statswpadj)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }

    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        filter(meth_bin %in% labels[1:2]) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_first2bins.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 8, 4, pl = p)

    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        filter(meth_bin %in% labels[1:3]) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_first3bins.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    # p <- by_sample %>%
    #     filter(subset == "400to600") %>%
    #     mutate(meth_bin = as.character(meth_bin)) %>%
    #     left_join(sample_table) %>%
    #     ggplot(aes(x = Neuron, y = prop_in_bin, color = condition)) +
    #     geom_point(size = 2.5, alpha = 0.8) +
    #     # geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
    #     facet_wrap(vars(meth_bin), scales = "free_y") +
    #     labs(x = "Neuron Fraction (Exc + Inh)", y = "Proportion of Reads in Bin") +
    #     ggtitle("Read Methylation vs Neuron Fraction") +
    #     mtclosedgridh +
    #     scale_conditions
    # mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/scatter_neuron_vs_propinbin1.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 12, 8, pl = p)


    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        filter(meth_bin %in% labels[1:3]) %>%
        left_join(sample_table) %>%
        mutate(sample_label = sprintf("%s_%s_%.2f", sample_name, ancestry, Neuron)) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        ggrepel::geom_text_repel(aes(y = prop_in_bin, group = condition, label = sample_label), position = position_dodge(width = 0.9), size = 1.8, max.overlaps = 20, segment.size = 0.3) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation (labeled)")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_first3bins_labeled.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 12, 5, pl = p)


    p1 <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        filter(meth_bin %in% labels[1:2]) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    p2 <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        filter(meth_bin %in% labels[3:4]) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    p <- p1 / p2
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_first_highlowstack1.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 7, 6, pl = p)
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_first_highlowstack11.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 8, 7, pl = p)

    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        mutate(methbintype = ifelse(meth_bin %in% labels[1:2], "low", "high")) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_grid(cols = vars(methbintype, subset), scales = "free") +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_binsplit131.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    p <- by_sample %>%
        mutate(meth_bin = as.character(meth_bin)) %>%
        ggplot(aes(x = meth_bin)) +
        stat_summary(aes(y = prop_in_bin, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = prop_in_bin, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_bin, y = prop_in_bin, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation bin", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/barplot_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    p <- by_read %>%
        mutate(condition = factor(condition, levels = ifelse(is.na(condition2), c(condition1), c(condition2, condition1)))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/density_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 8, 4, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = ifelse(is.na(condition2), c(condition1), c(condition2, condition1)))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/density.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 8, 4, pl = p)


    # Get group-wise mean fraction_meth
    # beta_overlays <- by_read %>%
    # filter(subset != "400to600") %>%
    # group_by(condition, subset) %>%
    # summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop")

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(if (is.na(condition2)) condition1 else c(condition2, condition1)))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram_betaoverlay1.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(if (is.na(condition2)) condition1 else c(condition2, condition1)))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.1f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(if (is.na(condition2)) condition1 else c(condition2, condition1)))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram_smaller.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 7, 3, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = c(if (is.na(condition2)) condition1 else c(condition2, condition1)))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram_smaller_2.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 7, 4, pl = p)

    by_read %$% fraction_meth %>% quantile()

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))
    p <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram_betaoverlay_L1HS_4q28.3_9.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_7q11.22_1") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))
    p <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_7q11.22_1") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/readhistogram_betaoverlay_L1HS_7q11.22_1_.pdf", mod_code_var, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    by_gene_id %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = max(prop_in_bin)) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        filter(subset == "328") %$% condition %>%
        table()

    tryCatch(
        {
            p <- by_gene_id %>%
                mutate(meth_bin = as.character(meth_bin)) %>%
                filter(meth_bin == "0.5") %>%
                filter(subset == "0to909") %>%
                group_by(gene_id) %>%
                mutate(samples_detected_per_element = n()) %>%
                ungroup() %>%
                filter(samples_detected_per_element > 11) %>%
                distinct() %>%
                # group_by(direction) %>%
                tidyHeatmap::heatmap(gene_id, sample, prop_in_bin,
                    cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE,
                    show_row_dend = FALSE,
                    palette_value = circlize::colorRamp2(
                        c(0, 0.00001, seq(0.1, 1, length.out = 3)),
                        c("black", rev(RColorBrewer::brewer.pal(4, "Oranges")))
                    )
                )
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/heatmap_prop_in_bin_gene_id_consistent_across_samples.pdf", mod_code_var, region, required_fraction_of_total_cg), 5, 5, pl = p)

            #
            df_filtered <- by_gene_id %>%
                mutate(meth_bin = as.character(meth_bin)) %>%
                filter(meth_bin == "0.5", subset == "0to909") %>%
                group_by(gene_id) %>%
                mutate(samples_detected_per_element = n()) %>%
                ungroup() %>%
                filter(samples_detected_per_element > 11) %>%
                distinct()
            df_sorted <- df_filtered %>%
                group_by(sample) %>%
                arrange(prop_in_bin) %>%
                mutate(gene_rank = row_number()) %>%
                ungroup()
            p <- ggplot(df_sorted, aes(x = sample, y = gene_rank, fill = prop_in_bin)) +
                geom_tile() +
                scale_fill_gradientn(
                    colours = c("lightblue", rev(RColorBrewer::brewer.pal(4, "Oranges"))),
                    values = rescale(c(0, 0.00001, seq(0.1, 1, length.out = 3))),
                    name = "prop_in_bin"
                ) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_reverse(expand = c(0, 0)) + # No padding; bottom = high rank
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                labs(x = "sample", y = "genes (sorted by prop_in_bin per sample)") +
                mtclosed +
                theme(
                    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank()
                )
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/heatmap_prop_in_bin_noNA.pdf", mod_code_var, region, required_fraction_of_total_cg), 5, 5, pl = p)
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/heatmap_prop_in_bin_noNA.pdf", mod_code_var, region, required_fraction_of_total_cg), 5, 5, pl = p, raster = TRUE)
        },
        error = function(e) {

        }
    )

    #

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = max(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = max(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = max(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_xlim.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = max(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_xlim_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = mean(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_mean_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = mean(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_mean.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = mean(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_mean_xlim_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_bin, condition, subset) %>%
        summarise(max_frac = mean(prop_in_bin)) %>%
        left_join(rmannextended) %>%
        group_by(meth_bin, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_bin, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/roc_mean_xlim.pdf", mod_code_var, region, required_fraction_of_total_cg), 9, 4, pl = p)

    named_group_split <- function(.tbl, ...) {
        grouped <- group_by(.tbl, ...)
        names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

        grouped %>%
            group_split() %>%
            rlang::set_names(names)
    }

    if (conf$single_condition == "no") {
        groups_needed <- length(condition1samples) + 1
        groups_needed <- 10
        n_under_consideration <- by_gene_id %>%
            ungroup() %>%
            mutate(split_var = paste0(gene_id, "/", subset, "/", meth_bin)) %>%
            left_join(sample_table) %>%
            filter(group_size >= groups_needed) %>%
            group_by(subset_bin) %>%
            dplyr::select(subset_bin, gene_id) %>%
            distinct() %>%
            summarise(n_under_consideration = n())
        group_size_df <- by_gene_id %>% dplyr::select(gene_id, group_size, subset, meth_bin)

        # dat_tmp <- labels %>%
        #     map(~ by_read %>%
        #         mutate(meth_bin = cut(fraction_meth, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)) %>%
        #         mutate(unmeth = as.integer(ifelse(meth_bin == .x, 1, 0))) %>%
        #         filter(subset != "400to600") %>%
        #         group_by(sample, gene_id, subset) %>%
        #         mutate(total_read = n()) %>%
        #         ungroup() %>%
        #         group_by(sample, gene_id, subset, meth_bin) %>%
        #         summarise(unmeth = sum(unmeth), total_read = dplyr::first(total_read)) %>%
        #         ungroup()) %>%
        #     list_rbind() %>%
        #     dplyr::rename(sample_name = sample) %>%
        #     mutate(split_var = paste0(gene_id, "/", subset, "/", meth_bin)) %>%
        #     left_join(group_size_df %>% distinct()) %>%
        #     left_join(sample_table) %>%
        #     filter(group_size >= groups_needed) %>%
        #     mutate(age_z = as.numeric(scale(age)))


        dat_tmp <- by_read %>%
            mutate(unmeth = as.integer(ifelse(fraction_meth <= 0.5, 1, 0))) %>%
            filter(subset != "400to600") %>%
            group_by(sample, gene_id, subset) %>%
            mutate(total_read = n()) %>%
            ungroup() %>%
            group_by(sample, gene_id, subset) %>%
            summarise(unmeth = sum(unmeth), total_read = dplyr::first(total_read)) %>%
            ungroup() %>%
            dplyr::rename(sample_name = sample) %>%
            mutate(split_var = paste0(gene_id, "/", subset)) %>%
            left_join(group_size_df %>% distinct() %>% filter(meth_bin == "[0, 0.5)")) %>%
            left_join(sample_table) %>%
            filter(group_size >= groups_needed) %>%
            mutate(age_z = as.numeric(scale(age)))


        dat_tmp <- dat_tmp %>% filter(meth_bin == "[0, 0.5)")

        # Create an empty list to hold tidy model results
        stats_list <- list()
        # Split data by group
        dat_groups <- named_group_split(dat_tmp, split_var)
        library(broom.mixed)
        # Loop through each group
        for (i in seq_along(dat_groups)) {
            print(i)
            group_data <- dat_groups[[i]]
            group_name <- names(dat_groups)[i]

            print(group_name)

            # Check if unmeth values are constant (e.g. all 0s or all the same)
            if (length(unique(group_data$unmeth)) == 1 | sum(group_data$unmeth > 0) < 2) {
                print("no info")
                tidy_model <- tibble(
                    term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                    subset = group_name,
                    note = "No variation in unmeth"
                )
            } else {
                tryCatch({
                    model <- glmmTMB(
                        cbind(unmeth, total_read - unmeth) ~
                            condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name),
                        data = group_data,
                        family = binomial()
                    )
                    tidy_model <- broom::tidy(model) %>%
                        mutate(subset = group_name)
                }, error = function(e) {
                    print(paste("Model failed for", group_name, ":", e$message))
                    tidy_model <<- tibble(
                        term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                        subset = group_name, note = paste("Model error:", e$message)
                    )
                })
            }
            stats_list[[i]] <- tidy_model
        }

        # Combine the results into a single data frame
        stats <- purrr::reduce(stats_list, bind_rows) %>%
            tidyr::separate(subset, into = c("gene_id", "subset", "meth_bin"), sep = "/", convert = TRUE)


        write_csv(stats, sprintf("ldna/results/%s/tables/reads_new_withendfilter/%s_%s/by_gene_new.csv", mod_code_var, region, required_fraction_of_total_cg))

        gene_condition_stats <- stats %>%
            filter(term == paste0("condition", condition2)) %>%
            mutate(p.value = case_when(
                is.nan(p.value) ~ 1,
                TRUE ~ p.value
            )) %>%
            mutate(statistic = case_when(
                is.nan(statistic) ~ 0,
                TRUE ~ statistic
            ))

        tryCatch(
            {
                p <- stats %>%
                    mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                    filter(p.value <= 0.05) %>%
                    filter(grepl("condition", term)) %>%
                    mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                    count(meth_bin, subset, dir_stat) %>%
                    complete(meth_bin, subset, dir_stat, fill = list(n = 0)) %>%
                    ggplot(aes(x = as.character(meth_bin), y = n, fill = dir_stat)) +
                    geom_col(position = "dodge", color = "black") +
                    facet_wrap(vars(subset), nrow = 1) +
                    ggtitle(sprintf("Significant Loci")) +
                    labs(x = "Meth bin", y = sprintf("Number DM")) +
                    mtclosed +
                    anchorbar +
                    scale_methylation
                mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/num_de_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 5.5, 4, pl = p)
                tryCatch(
                    {
                        p <- stats %>%
                            mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                            filter(p.value <= 0.05) %>%
                            filter(grepl("condition", term)) %>%
                            mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                            count(meth_bin, subset, dir_stat) %>%
                            complete(meth_bin, subset, dir_stat, fill = list(n = 0)) %>%
                            ggplot(aes(x = as.character(meth_bin), y = n, fill = dir_stat)) +
                            geom_col(position = "dodge", color = "black") +
                            facet_wrap(vars(subset), nrow = 1) +
                            ggtitle(sprintf("Significant Loci")) +
                            labs(x = "Meth bin", y = sprintf("Number DM")) +
                            mtclosed +
                            anchorbar +
                            scale_methylation
                        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/num_de_wasp.pdf", mod_code_var, region, required_fraction_of_total_cg), 5.5, 4, pl = p)

                        p <- stats %>%
                            filter(subset != "400to600") %>%
                            mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                            filter(p.value <= 0.05) %>%
                            filter(grepl("condition", term)) %>%
                            mutate(meth_bin = factor(as.character(meth_bin), levels = c("0.1", "0.25", "0.5"))) %>%
                            mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                            count(meth_bin, subset, dir_stat) %>%
                            complete(meth_bin, subset, dir_stat, fill = list(n = 0)) %>%
                            ggplot(aes(x = as.character(meth_bin), y = n, fill = dir_stat)) +
                            geom_col(position = "dodge", color = "black") +
                            facet_wrap(vars(subset), nrow = 1) +
                            ggtitle(sprintf("Significant Loci")) +
                            labs(x = "Meth bin", y = sprintf("Number DM")) +
                            mtclosed +
                            anchorbar +
                            scale_methylation
                        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new_withendfilter/%s_%s/num_de.pdf", mod_code_var, region, required_fraction_of_total_cg), 5, 4, pl = p)
                    },
                    error = function(e) {
                        print("no sig de")
                    }
                )

                dispersion_models <- list()
                for (subsetofinterest in by_cpg$subset %>%
                    unique() %>%
                    grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
                    dispersion_model <- glmmTMB(
                        cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                        dispformula = ~condition,
                        family = glmmTMB::betabinomial(),
                        data = by_cpg %>% filter(subset == subsetofinterest) %>%
                            mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                    )
                    dispersion_models[[subsetofinterest]] <- dispersion_model
                }

                res_text <- capture.output(map(dispersion_models, summary))
                writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))
                # summary_disp <- summary(mod)$dispersion

                # # # Extract estimates and standard errors
                # log_phi_control <- 2.233794
                # se_log_phi_control <- 0.002305

                # log_phi_diff_AD <- -0.337098
                # se_log_phi_diff_AD <- 0.003029

                # # # Compute log(φ) for AD and its SE
                # log_phi_AD <- log_phi_control + log_phi_diff_AD
                # se_log_phi_AD <- sqrt(se_log_phi_control^2 + se_log_phi_diff_AD^2)

                # # # 95% CI on log-scale
                # z <- 1.96
                # ci_log_phi_control <- log_phi_control + c(-1, 1) * z * se_log_phi_control
                # ci_log_phi_AD <- log_phi_AD + c(-1, 1) * z * se_log_phi_AD

                # # # Exponentiate to get CIs for φ
                # phi_control <- exp(log_phi_control)
                # phi_AD <- exp(log_phi_AD)

                # ci_phi_control <- exp(ci_log_phi_control)
                # ci_phi_AD <- exp(ci_log_phi_AD)


                # dispersion_models_withreadid <- list()
                # for (subsetofinterest in by_cpg$subset %>%
                #     unique() %>%
                #     grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
                #     dispersion_model <- glmmTMB(
                #         cbind(meth, unmeth) ~ 1 + (1 | read_id),
                #         dispformula = ~condition,
                #         family = glmmTMB::betabinomial(),
                #         data = by_cpg %>% filter(subset == subsetofinterest) %>%
                #             mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                #     )
                #     dispersion_models[[subsetofinterest]] <- dispersion_model
                # }

                # res_text <- capture.output(map(dispersion_models, summary))
                # writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))


                # dat <- by_cpg %>%
                #     filter(subset == "0to909") %>%
                #     mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                # mod0 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                #     data = dat, family = binomial()
                # )

                # mod05 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = dat, family = glmmTMB::betabinomial()
                # )
                # mod05disp <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     dispformula = ~1,
                #     data = dat, family = glmmTMB::betabinomial()
                # )
                # datmorebinom <- dat %>%
                #     filter(fraction_meth < 0.95) %>%
                #     filter(fraction_meth > 0.35)
                # mod05morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = datmorebinom, family = glmmTMB::betabinomial()
                # )
                # mod05binom_morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = datmorebinom, family = binomial()
                # )


                # dispersion_model_morebinomdat <- glmmTMB(
                #     cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                #     dispformula = ~condition,
                #     family = glmmTMB::betabinomial(),
                #     data = datmorebinom
                # )

                # anova(mod0, mod05)
                # anova(mod05binom_morebinom, mod05morebinom)
            },
            error = function(e) {}
        )
    }
}


read_analysis2(readscg_endfiltnew %>% filter(mod_code == params$mod_code), cg_indices)
# read_analysis2(readscg_endfiltnew %>% filter(mod_code == params$mod_code), cg_indices)




# //ANCHOR - READ ANALYSIS 1
read_analysis1 <- function(
    df,
    cg_indices,
    mod_code_var = params$mod_code,
    regions_of_interest = list(c(0, 328), c(0, 500), c(0, 909), c(400, 600)),
    required_fraction_of_total_cg = 0.75,
    meth_thresholds = c(0.25, 0.5, 0.75),
    context = "CpG") {
    region <- "L1HS_FL"
    outputdirtables <- sprintf("ldna/results/%s/tables/reads_new/%s_%s", params$mod_code, region, required_fraction_of_total_cg)
    dir.create(outputdirtables, recursive = TRUE)

    readsdf1 <- df %>%
        filter(mod_code == mod_code_var) %>%
        left_join(rmannextended %>%
            dplyr::select(gene_id, start, end, strand, rte_length_req, intactness_req) %>%
            dplyr::rename(element_strand = strand, element_start = start, element_end = end)) %>%
        filter(rte_length_req == "FL")


    by_cpg_l <- list()
    by_read_l <- list()
    by_sample_l <- list()
    by_gene_id_l <- list()

    for (region_of_interest in regions_of_interest) {
        roistart <- region_of_interest[1]
        roiend <- region_of_interest[2]
        roistring <- paste0(roistart, "to", roiend)

        numCGneeded <- ceiling(length(cg_indices[(cg_indices <= roiend) & (cg_indices >= roistart)]) * required_fraction_of_total_cg)

        utr1 <- readsdf1 %>%
            filter(case_when(
                element_strand == "+" ~ (start > element_start + roistart) & (start < element_start + roiend),
                element_strand == "-" ~ (start > element_end - roiend) & (start < element_end - roistart)
            )) %>%
            dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0))

        by_cpg_temp <- utr1 %>%
            group_by(gene_id, read_id, condition, sample) %>%
            mutate(num_cpgs_in_read = n()) %>%
            mutate(fraction_meth = mean(mod_indicator)) %>%
            relocate(gene_id) %>%
            ungroup()

        by_read_temp <- by_cpg_temp %>%
            filter(num_cpgs_in_read >= numCGneeded) %>%
            group_by(read_id, gene_id, sample, condition, region) %>%
            summarise(fraction_meth = dplyr::first(fraction_meth), num_cpgs_in_read = dplyr::first(num_cpgs_in_read), strand = dplyr::first(element_strand), numCGneeded = dplyr::first(numCGneeded)) %>%
            ungroup()

        by_cpg_temp$subset <- as.character(roistring)
        by_read_temp$subset <- as.character(roistring)

        by_cpg_l[[as.character(roistring)]] <- by_cpg_temp
        by_read_l[[as.character(roistring)]] <- by_read_temp

        for (meth_threshold in meth_thresholds) {
            subset_threshold <- paste0(roistring, "_", meth_threshold)
            by_sample_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth > meth_threshold, 0, 1)) %>%
                group_by(sample, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                group_by(condition, region) %>%
                mutate(meanProp = mean(propUnmeth))
            by_sample_temp$meth_threshold <- meth_threshold
            by_sample_temp$subset <- as.character(roistring)
            by_sample_temp$subset_threshold <- subset_threshold

            by_gene_id_temp <- by_read_temp %>%
                mutate(unmeth = ifelse(fraction_meth >= meth_threshold, 0, 1)) %>%
                group_by(sample, gene_id, region, condition) %>%
                summarise(propUnmeth = mean(unmeth)) %>%
                ungroup() %>%
                group_by(gene_id) %>%
                mutate(group_size = n()) %>%
                ungroup()
            by_gene_id_temp$meth_threshold <- meth_threshold
            by_gene_id_temp$subset <- as.character(roistring)
            by_gene_id_temp$subset_threshold <- subset_threshold

            by_sample_l[[subset_threshold]] <- by_sample_temp
            by_gene_id_l[[subset_threshold]] <- by_gene_id_temp
        }
    }

    by_cpg <- purrr::reduce(by_cpg_l, bind_rows)
    by_read <- purrr::reduce(by_read_l, bind_rows)
    by_gene_id <- purrr::reduce(by_gene_id_l, bind_rows)
    by_sample <- purrr::reduce(by_sample_l, bind_rows)

    dir.create(dirname(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)), recursive = TRUE)
    by_cpg %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region))
    by_read %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region))
    by_gene_id %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region))
    by_sample %>% write_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region))

    # by_cpg <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_cpg.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_read <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_read.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_gene_id <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_gene_id.csv", params$mod_code, region)) %>%
    #     mutate(condition = factor(condition, levels = conf$levels)) %>%
    #     mutate(sample = factor(sample, levels = conf$samples))
    # by_sample <- read_csv(sprintf("ldna/Rintermediates/%s/%s/highly_demethylated_reads_by_sample.csv", params$mod_code, region)) %>%
    #         mutate(condition = factor(condition, levels = conf$levels)) %>%
    #         mutate(sample = factor(sample, levels = conf$samples))

    get_read_ecdf <- function(df, subset_val, breakpoints, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by specified variable
            summarise(
                percent_below = list(ecdf(fraction_meth)(breakpoints) * 100), # Compute ECDF
                .groups = "drop"
            ) %>%
            unnest_longer(percent_below) %>%
            mutate(
                threshold = rep(breakpoints, times = n() / length(breakpoints)), # Expand breakpoints
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), threshold, percent_below) # Keep relevant columns
    }


    get_read_quantiles <- function(df, subset_val, probs, group_var = NULL) {
        df %>%
            filter(subset == subset_val) %>%
            group_by(across(all_of(group_var))) %>% # Group by sample_name or another variable
            summarise(
                quantiles = list(quantile(fraction_meth, probs = probs)), # Store as list
                .groups = "drop"
            ) %>%
            unnest_longer(quantiles) %>%
            mutate(
                quantile = rep(probs, times = n() / length(probs)), # Expand probs for each group
                subset = subset_val
            ) %>%
            dplyr::select(subset, all_of(group_var), quantile, mean_meth = quantiles) # Keep relevant columns
    }
    subsets <- c("0to909", "0to500", "0to328")
    breakpoints <- seq(0, 1, 0.05)

    ecdf_reads <- map_dfr(subsets, ~ get_read_ecdf(by_read, .x, breakpoints, "sample"))
    quantile_reads <- map_dfr(subsets, ~ get_read_quantiles(by_read, .x, breakpoints, "sample"))


    ecdf_reads %>% write_mycsv(sprintf("%s/read_ecdf.csv", outputdirtables))
    quantile_reads %>% write_mycsv(sprintf("%s/read_quantiles.csv", outputdirtables))

    ecdf_reads_acrosssamplemean <- ecdf_reads %>%
        group_by(subset, threshold) %>%
        summarise(percent_below = mean(percent_below))
    quantile_reads_acrosssamplemean <- quantile_reads %>%
        group_by(subset, quantile) %>%
        summarise(mean_meth = mean(mean_meth))

    ecdf_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_ecdf_acrosssamplemean.csv", outputdirtables))
    quantile_reads_acrosssamplemean %>% write_mycsv(sprintf("%s/read_quantiles_acrosssamplemean.csv", outputdirtables))


    # Example ECDF plot function
    plot_ecdf <- function(ecdf_data) {
        ggplot(ecdf_data, aes(x = threshold, y = percent_below, color = as.factor(subset))) +
            geom_line() + # ECDF curve
            geom_point() + # Add points for clarity
            labs(
                title = "Empirical CDF of Fraction Methylation",
                x = "Fraction Methylation Threshold",
                y = "Percent Below Threshold",
                color = "Subset"
            ) +
            theme_minimal()
    }


    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        ggplot(aes(x = threshold, y = percent_below, color = sample)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_samples_unique +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)

    p <- ecdf_reads %>%
        filter(subset == "0to909") %>%
        left_join(sample_table) %>%
        group_by(condition, threshold) %>%
        summarise(percent_below = mean(percent_below)) %>%
        ggplot(aes(x = threshold, y = percent_below, color = condition)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Read Methylation Empirical CDF",
            x = "Methylation Threshold",
            y = "Percent Below Threshold",
            color = "Subset"
        ) +
        geom_vline(xintercept = c(0.1, 0.25, 0.5), linetype = "dashed") +
        coord_cartesian(ylim = c(0, 20)) +
        scale_conditions +
        mtclosed
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/ecdf_909_by_condition_zoom.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 5, 4, pl = p)


    p <- by_sample %>%
        filter(subset != "400to600") %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        # stats <- by_sample %>%
        #     filter(subset != "400to600") %>%
        #     ungroup() %>%
        #     left_join(sample_table) %>%
        #     group_split(subset_threshold) %>%
        #     set_names(unique(by_sample %>% filter(subset != "400to600") %$% subset_threshold)) %>%
        #     map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
        #     imap_dfr(~ .x %>% mutate(subset = .y))
        stats_list <- list()
        i <- 1
        for (subset in unique(by_cpg$subset)) {
            for (threshold in unique(by_sample$meth_threshold)) {
                by_cpg_tmp <- by_cpg %>%
                    filter(subset == !!subset) %>%
                    mutate(unmeth = ifelse(fraction_meth > threshold, 0, 1)) %>%
                    dplyr::rename(sample_name = sample) %>%
                    group_by(sample_name, condition, gene_id) %>%
                    summarise(unmeth = sum(unmeth), total = n()) %>%
                    ungroup() %>%
                    left_join(sample_table) %>%
                    mutate(age_z = as.numeric(scale(age)))
                model_tmp <- glmmTMB(
                    cbind(unmeth, total - unmeth) ~
                        condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name) + (1 | gene_id),
                    data = by_cpg_tmp,
                    family = binomial()
                )
                res_tmp <- broom::tidy(model_tmp) %>%
                    mutate(threshold = !!threshold) %>%
                    mutate(subset = !!subset)
                stats_list[[i]] <- res_tmp
                i <- i + 1
            }
        }
        stats <- purrr::reduce(stats_list, bind_rows)
        stats_padj <- stats %>%
            filter(term == paste0("condition", condition2)) %>%
            filter(subset != "400to600") %>%
            mutate(padj = p.adjust(p.value, method = "fdr"))
        statswpadj <- stats %>% left_join(stats_padj)
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = statswpadj)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }
    p <- by_sample %>%
        mutate(meth_threshold = as.character(meth_threshold)) %>%
        ggplot(aes(x = meth_threshold)) +
        stat_summary(aes(y = propUnmeth, group = condition, fill = condition), color = "black", fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
        geom_point(aes(y = propUnmeth, group = condition), position = position_dodge(width = 0.9)) +
        facet_wrap(vars(subset), nrow = 1) +
        geom_pwc(aes(x = meth_threshold, y = propUnmeth, group = condition), tip.length = 0, method = "wilcox_test") +
        labs(x = "Methylation Threshold", y = sprintf("Reads Fraction < # methylated")) +
        ggtitle(sprintf("Read Methylation")) +
        mtclosedgridh +
        scale_conditions +
        anchorbar
    if (enough_samples_per_condition_for_stats) {
        stats <- by_sample %>%
            ungroup() %>%
            left_join(sample_table) %>%
            group_split(subset_threshold) %>%
            set_names(unique(by_sample$subset_threshold)) %>%
            map(~ broom::tidy(summary(lm(formula(sprintf("%s ~ %s", "propUnmeth", lm_right_hand_side)), .x)))) %>%
            imap_dfr(~ .x %>% mutate(subset = .y))
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p, sf = stats)
    } else {
        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/barplot_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)
    }

    p <- by_read %>%
        mutate(condition = factor(condition, levels = ifelse(is.na(condition2), c(condition1), c(condition2, condition1)))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = ifelse(is.na(condition2), c(condition1), c(condition2, condition1)))) %>%
        ggplot() +
        geom_density(aes(x = fraction_meth, fill = condition), alpha = 0.7) +
        facet_wrap(vars(subset), nrow = 1) +
        ggtitle(sprintf("Read Density")) +
        labs(x = "", y = sprintf("Read Density")) +
        mtclosed +
        scale_conditions
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/density.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 8, 4, pl = p)


    # Get group-wise mean fraction_meth
    # beta_overlays <- by_read %>%
    # filter(subset != "400to600") %>%
    # group_by(condition, subset) %>%
    # summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop")

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))

    p <- by_read %>%
        filter(subset != "400to600") %>%
        mutate(condition = factor(condition, levels = ifelse(is.na(condition2), c(condition1), c(condition2, condition1)))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)

    beta_overlays <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        group_by(subset) %>%
        summarise(mean_meth = mean(fraction_meth, na.rm = TRUE), numCGneeded = dplyr::first(numCGneeded), .groups = "drop") %>%
        mutate(
            alpha = mean_meth * numCGneeded,
            beta = (1 - mean_meth) * numCGneeded
        ) %>%
        rowwise() %>%
        mutate(
            x = list(seq(0, 1, length.out = 200)),
            y = list(dbeta(x, alpha, beta))
        ) %>%
        unnest(c(x, y))
    p <- by_read %>%
        filter(subset != "400to600") %>%
        filter(gene_id == "L1HS_4q28.3_9") %>%
        mutate(condition = factor(condition, levels = c(condition2, condition1))) %>%
        ggplot() +
        geom_histogram(
            data = . %>% filter(condition == condition2),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_histogram(
            data = . %>% filter(condition == condition1),
            aes(x = fraction_meth, fill = condition, y = after_stat(count / sum(count))),
            alpha = 0.7,
            bins = 30,
            position = "identity"
        ) +
        geom_line(
            data = beta_overlays,
            aes(x = x, y = 7.5 * y / sum(y)), # normalize for proportion scale
            size = 0.8,
            color = "green",
            inherit.aes = FALSE
        ) +
        facet_wrap(vars(subset), nrow = 1) + # allows y-axis to vary if needed
        ggtitle("Read Methylation Histogram (Proportional)") +
        labs(x = "Read Methylation", y = "Proportion of Reads") +
        mtclosedgridh +
        scale_conditions +
        scale_y_continuous(expand = expansion(mult = c(0, .075)), labels = function(x) sprintf("%.3f", x))
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/readhistogram_betaoverlay_L1HS_4q28.3_9.pdf", params$mod_code, region, required_fraction_of_total_cg, context), 9, 4, pl = p)


    by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        filter(subset == "328") %$% condition %>%
        table()

    tryCatch(
        {
            p <- by_gene_id %>%
                mutate(meth_threshold = as.character(meth_threshold)) %>%
                filter(meth_threshold == "0.5") %>%
                filter(subset == "0to909") %>%
                group_by(gene_id) %>%
                mutate(samples_detected_per_element = n()) %>%
                ungroup() %>%
                filter(samples_detected_per_element > 11) %>%
                distinct() %>%
                # group_by(direction) %>%
                tidyHeatmap::heatmap(gene_id, sample, propUnmeth,
                    cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE,
                    show_row_dend = FALSE,
                    palette_value = circlize::colorRamp2(
                        c(0, 0.00001, seq(0.1, 1, length.out = 3)),
                        c("black", rev(RColorBrewer::brewer.pal(4, "Oranges")))
                    )
                )
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_gene_id_consistent_across_samples.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)

            #
            df_filtered <- by_gene_id %>%
                mutate(meth_threshold = as.character(meth_threshold)) %>%
                filter(meth_threshold == "0.5", subset == "0to909") %>%
                group_by(gene_id) %>%
                mutate(samples_detected_per_element = n()) %>%
                ungroup() %>%
                filter(samples_detected_per_element > 11) %>%
                distinct()
            df_sorted <- df_filtered %>%
                group_by(sample) %>%
                arrange(propUnmeth) %>%
                mutate(gene_rank = row_number()) %>%
                ungroup()
            p <- ggplot(df_sorted, aes(x = sample, y = gene_rank, fill = propUnmeth)) +
                geom_tile() +
                scale_fill_gradientn(
                    colours = c("lightblue", rev(RColorBrewer::brewer.pal(4, "Oranges"))),
                    values = rescale(c(0, 0.00001, seq(0.1, 1, length.out = 3))),
                    name = "propUnmeth"
                ) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_reverse(expand = c(0, 0)) + # No padding; bottom = high rank
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                labs(x = "sample", y = "genes (sorted by propUnmeth per sample)") +
                mtclosed +
                theme(
                    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank()
                )
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p)
            mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/heatmap_propunmeth_noNA.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 5, pl = p, raster = TRUE)
        },
        error = function(e) {

        }
    )

    #

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = max(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < %s methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    p <- by_gene_id %>%
        filter(subset != "400to600") %>%
        group_by(gene_id, meth_threshold, condition, subset) %>%
        summarise(max_frac = mean(propUnmeth)) %>%
        left_join(rmannextended) %>%
        group_by(meth_threshold, condition, subset) %>%
        arrange(max_frac) %>%
        mutate(ranked_row = row_number()) %>%
        mutate(condition = paste0(meth_threshold, "\n", condition)) %>%
        ggplot(aes(x = max_frac, y = ranked_row, color = condition)) +
        geom_point() +
        xlim(c(0, 0.25)) +
        facet_wrap(vars(subset), nrow = 1) +
        labs(title = sprintf("Read Methylation", roistring), y = "Unique Locus Rank", x = sprintf("Reads Fraction < # methylated", "#")) +
        mtclosedgridh +
        scale_color_brewer(palette = "Paired") +
        anchorbar
    mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/roc_mean_xlim.pdf", params$mod_code, region, required_fraction_of_total_cg), 9, 4, pl = p)

    named_group_split <- function(.tbl, ...) {
        grouped <- group_by(.tbl, ...)
        names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

        grouped %>%
            group_split() %>%
            rlang::set_names(names)
    }

    if (conf$single_condition == "no") {
        groups_needed <- length(condition1samples) + 1
        groups_needed <- 10
        n_under_consideration <- by_gene_id %>%
            ungroup() %>%
            mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
            left_join(sample_table) %>%
            filter(group_size >= groups_needed) %>%
            group_by(subset_threshold) %>%
            dplyr::select(subset_threshold, gene_id) %>%
            distinct() %>%
            summarise(n_under_consideration = n())
        group_size_df <- by_gene_id %>% dplyr::select(gene_id, group_size, subset, meth_threshold)

        dat_tmp <- meth_thresholds %>%
            map(~ by_read %>%
                mutate(unmeth = ifelse(fraction_meth > .x, 0, 1), meth_threshold = .x) %>%
                filter(subset != "400to600") %>%
                group_by(sample, gene_id, subset) %>%
                mutate(total_read = n()) %>%
                ungroup() %>%
                group_by(sample, gene_id, subset, meth_threshold) %>%
                summarise(unmeth = sum(unmeth), total_read = dplyr::first(total_read)) %>%
                ungroup()) %>%
            list_rbind() %>%
            dplyr::rename(sample_name = sample) %>%
            mutate(split_var = paste0(gene_id, "/", subset, "/", meth_threshold)) %>%
            left_join(group_size_df %>% distinct()) %>%
            left_join(sample_table) %>%
            filter(group_size >= groups_needed) %>%
            mutate(age_z = as.numeric(scale(age)))




        # Create an empty list to hold tidy model results
        stats_list <- list()

        # Split data by group
        dat_groups <- named_group_split(dat_tmp, split_var)
        library(broom.mixed)
        # Loop through each group
        for (i in seq_along(dat_groups)) {
            print(i)
            group_data <- dat_groups[[i]]
            group_name <- names(dat_groups)[i]

            print(group_name)

            # Check if unmeth values are constant (e.g. all 0s or all the same)
            if (length(unique(group_data$unmeth)) == 1 | sum(group_data$unmeth > 0) < 2) {
                print("no info")
                tidy_model <- tibble(
                    term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                    subset = group_name,
                    note = "No variation in unmeth"
                )
            } else {
                tryCatch({
                    model <- glmmTMB(
                        cbind(unmeth, total_read - unmeth) ~
                            condition + sex + ancestry + Oligo_z + Astro_z + Micro_z + Inh_z + Exc_z + OPC_z + (1 | sample_name),
                        data = group_data,
                        family = binomial()
                    )
                    tidy_model <- broom::tidy(model) %>%
                        mutate(subset = group_name)
                }, error = function(e) {
                    print(paste("Model failed for", group_name, ":", e$message))
                    tidy_model <<- tibble(
                        term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                        subset = group_name, note = paste("Model error:", e$message)
                    )
                })
            }
            stats_list[[i]] <- tidy_model
        }

        # Combine the results into a single data frame
        stats <- purrr::reduce(stats_list, bind_rows) %>%
            tidyr::separate(subset, into = c("gene_id", "subset", "meth_threshold"), sep = "/", convert = TRUE)


        write_csv(stats, sprintf("ldna/results/%s/tables/reads_new/%s_%s/by_gene_new.csv", params$mod_code, region, required_fraction_of_total_cg))

        gene_condition_stats <- stats %>%
            filter(term == paste0("condition", condition2)) %>%
            mutate(p.value = case_when(
                is.nan(p.value) ~ 1,
                TRUE ~ p.value
            )) %>%
            mutate(statistic = case_when(
                is.nan(statistic) ~ 0,
                TRUE ~ statistic
            ))

        tryCatch(
            {
                p <- stats %>%
                    mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                    filter(p.value <= 0.05) %>%
                    filter(grepl("condition", term)) %>%
                    mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                    count(meth_threshold, subset, dir_stat) %>%
                    complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                    ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                    geom_col(position = "dodge", color = "black") +
                    facet_wrap(vars(subset), nrow = 1) +
                    ggtitle(sprintf("Significant Loci")) +
                    labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                    mtclosed +
                    anchorbar +
                    scale_methylation
                mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 5.5, 4, pl = p)
                tryCatch(
                    {
                        p <- stats %>%
                            mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                            filter(p.value <= 0.05) %>%
                            filter(grepl("condition", term)) %>%
                            mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                            count(meth_threshold, subset, dir_stat) %>%
                            complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                            ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                            geom_col(position = "dodge", color = "black") +
                            facet_wrap(vars(subset), nrow = 1) +
                            ggtitle(sprintf("Significant Loci")) +
                            labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                            mtclosed +
                            anchorbar +
                            scale_methylation
                        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de_wasp.pdf", params$mod_code, region, required_fraction_of_total_cg), 5.5, 4, pl = p)

                        p <- stats %>%
                            filter(subset != "400to600") %>%
                            mutate(p.value = ifelse(is.nan(p.value), 1, p.value)) %>%
                            filter(p.value <= 0.05) %>%
                            filter(grepl("condition", term)) %>%
                            mutate(meth_threshold = factor(as.character(meth_threshold), levels = c("0.1", "0.25", "0.5"))) %>%
                            mutate(dir_stat = factor(ifelse(statistic > 0, "Hypo", "Hyper"), levels = c("Hyper", "Hypo"))) %>%
                            count(meth_threshold, subset, dir_stat) %>%
                            complete(meth_threshold, subset, dir_stat, fill = list(n = 0)) %>%
                            ggplot(aes(x = as.character(meth_threshold), y = n, fill = dir_stat)) +
                            geom_col(position = "dodge", color = "black") +
                            facet_wrap(vars(subset), nrow = 1) +
                            ggtitle(sprintf("Significant Loci")) +
                            labs(x = "Meth Threshold", y = sprintf("Number DM")) +
                            mtclosed +
                            anchorbar +
                            scale_methylation
                        mysaveandstore(sprintf("ldna/results/%s/plots/reads_new/%s_%s/num_de.pdf", params$mod_code, region, required_fraction_of_total_cg), 5, 4, pl = p)
                    },
                    error = function(e) {
                        print("no sig de")
                    }
                )

                dispersion_models <- list()
                for (subsetofinterest in by_cpg$subset %>%
                    unique() %>%
                    grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
                    dispersion_model <- glmmTMB(
                        cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                        dispformula = ~condition,
                        family = glmmTMB::betabinomial(),
                        data = by_cpg %>% filter(subset == subsetofinterest) %>%
                            mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                    )
                    dispersion_models[[subsetofinterest]] <- dispersion_model
                }

                res_text <- capture.output(map(dispersion_models, summary))
                writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))
                summary_disp <- summary(mod)$dispersion


                # could it be because read indels mess it up?
                by_cpg_withgene_id <- by_cpg %>%
                    filter(element_strand == "+", strand == "+") %>%
                    mutate(sequence_pos = (start - element_start) + 2) %>%
                    left_join(cg_positions_df)
                by_cpg_withgene_id %>%
                    filter(!is.na(consensus_pos)) %$% element_strand %>%
                    table()

                by_cpg_withgene_id <- by_cpg %>%
                    filter(element_strand == "+", strand == "-") %>%
                    mutate(sequence_pos = (start - element_start) + 1) %>%
                    left_join(cg_positions_df)
                by_cpg_withgene_id %>%
                    filter(!is.na(consensus_pos)) %$% element_strand %>%
                    table()

                by_cpg_withgene_id <- by_cpg %>%
                    filter(element_strand == "-", strand == "+") %>%
                    mutate(sequence_pos = (element_end - start) - 1) %>%
                    left_join(cg_positions_df)
                by_cpg_withgene_id %>%
                    filter(!is.na(consensus_pos)) %$% element_strand %>%
                    table()

                by_cpg_withgene_id <- by_cpg %>%
                    filter(element_strand == "-", strand == "-") %>%
                    mutate(sequence_pos = (element_end - start)) %>%
                    left_join(cg_positions_df)
                by_cpg_withgene_id %>%
                    filter(!is.na(consensus_pos)) %$% element_strand %>%
                    table()


                by_cpg_withgene_id <- by_cpg %>%
                    mutate(sequence_pos = case_when(
                        element_strand == "+" & strand == "+" ~ (start - element_start) + 2,
                        element_strand == "+" & strand == "-" ~ (start - element_start) + 1,
                        element_strand == "-" & strand == "+" ~ (element_end - start) - 1,
                        element_strand == "-" & strand == "-" ~ (element_end - start)
                    )) %>%
                    left_join(cg_positions_df)
                by_cpg_withgene_id <- by_cpg_withgene_id %>% filter(!is.na(consensus_pos))

                dispersion_models_withcpgid <- list()
                for (subsetofinterest in by_cpg$subset %>%
                    unique() %>%
                    grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
                    dispersion_model <- glmmTMB(
                        cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | consensus_pos) + (1 | sample:gene_id),
                        dispformula = ~condition,
                        family = glmmTMB::betabinomial(),
                        data = by_cpg_withgene_id %>% filter(subset == subsetofinterest) %>%
                            mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                    )
                    dispersion_models_withcpgid[[subsetofinterest]] <- dispersion_model
                }

                res_text_withcpgid <- capture.output(map(dispersion_models_withcpgid, summary))
                writeLines(res_text_withcpgid, sprintf("%s/dispersion_model_withcpgid_summary.txt", outputdirtables))
                summary_disp_withcpgid <- summary(mod)$dispersion
                summary(dispersion_model)
                # # Extract estimates and standard errors
                # log_phi_control <- 2.610989
                # se_log_phi_control <- 0.001710

                # log_phi_diff_AD <- -0.340575
                # se_log_phi_diff_AD <- 0.002263

                # # Compute log(φ) for AD and its SE
                # log_phi_AD <- log_phi_control + log_phi_diff_AD
                # se_log_phi_AD <- sqrt(se_log_phi_control^2 + se_log_phi_diff_AD^2)

                # # 95% CI on log-scale
                # z <- 1.96
                # ci_log_phi_control <- log_phi_control + c(-1, 1) * z * se_log_phi_control
                # ci_log_phi_AD <- log_phi_AD + c(-1, 1) * z * se_log_phi_AD

                # # Exponentiate to get CIs for φ
                # phi_control <- exp(log_phi_control)
                # phi_AD <- exp(log_phi_AD)

                # ci_phi_control <- exp(ci_log_phi_control)
                # ci_phi_AD <- exp(ci_log_phi_AD)


                # dispersion_models_withreadid <- list()
                # for (subsetofinterest in by_cpg$subset %>%
                #     unique() %>%
                #     grep(pattern = "400to600", ., invert = TRUE, value = TRUE)) {
                #     dispersion_model <- glmmTMB(
                #         cbind(meth, unmeth) ~ 1 + (1 | read_id),
                #         dispformula = ~condition,
                #         family = glmmTMB::betabinomial(),
                #         data = by_cpg %>% filter(subset == subsetofinterest) %>%
                #             mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                #     )
                #     dispersion_models[[subsetofinterest]] <- dispersion_model
                # }

                # res_text <- capture.output(map(dispersion_models, summary))
                # writeLines(res_text, sprintf("%s/dispersion_model_summary.txt", outputdirtables))


                # dat <- by_cpg %>%
                #     filter(subset == "0to909") %>%
                #     mutate(meth = fraction_meth * num_cpgs_in_read, unmeth = num_cpgs_in_read - meth)
                # mod0 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                #     data = dat, family = binomial()
                # )

                # mod05 <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = dat, family = glmmTMB::betabinomial()
                # )
                # mod05disp <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     dispformula = ~1,
                #     data = dat, family = glmmTMB::betabinomial()
                # )
                # datmorebinom <- dat %>%
                #     filter(fraction_meth < 0.95) %>%
                #     filter(fraction_meth > 0.35)
                # mod05morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = datmorebinom, family = glmmTMB::betabinomial()
                # )
                # mod05binom_morebinom <- glmmTMB(cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id),
                #     data = datmorebinom, family = binomial()
                # )


                # dispersion_model_morebinomdat <- glmmTMB(
                #     cbind(meth, unmeth) ~ 1 + (1 | sample) + (1 | gene_id) + (1 | sample:gene_id),
                #     dispformula = ~condition,
                #     family = glmmTMB::betabinomial(),
                #     data = datmorebinom
                # )

                # anova(mod0, mod05)
                # anova(mod05binom_morebinom, mod05morebinom)
            },
            error = function(e) {}
        )
    }
}

read_analysis1(readscg_endfiltnew, cg_indices)

file.create(sprintf("ldna/outfiles/bedmeth_%s_%s.done", "13_reads", params$mod_code))
