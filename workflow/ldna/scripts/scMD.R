module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
# source("workflow/scripts/defaults.R")
# source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(magrittr)
library(GenomicRanges)
library(configr)
library(Biostrings)
library(readr)
library(scMD)
library(dplyr)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

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
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            bedmethylpaths = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_bedMethyl.bed", samples, samples),
            data = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name),
            dmrs = sprintf("ldna/results/m/tables/%s/dmrs.tsv", conf$contrasts),
            dmls = sprintf("ldna/results/m/tables/%s/dmls.tsv", conf$contrasts),
            read_mods = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis),
            read_mods_cg = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", conf$rte_subfamily_read_level_analysis),
            read_mods_cg_islands = sprintf("ldna/intermediates/%s/methylation/analysis_default/%s_readmods_%s_%s.tsv", samples, samples, "CpG", "cpgI")
        ), env = globalenv())
        assign("params", list(mod_code = "m"), env = globalenv())
        assign("outputs", list(outfile = "ldna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

contrasts <- conf$contrasts

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis


##########################
# PREP DATA FOR ANALYSIS
sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE), col_names = FALSE)
    # Sum 5mC ("m") and 5hmC ("h") to match bisulfite references which cannot distinguish them
    df_mh <- df %>%
        filter(X4 %in% c("m", "h")) %>%
        group_by(X1, X2) %>%
        summarise(cov = dplyr::first(X10), pctM = sum(as.double(X11)), .groups = "drop")
    rm(df)
    gr <- GRanges(
        seqnames = df_mh$X1,
        ranges = IRanges(start = df_mh$X2 + 1, end = df_mh$X2 + 1), # bedMethyl is 0-based; GRanges is 1-based
        cov = df_mh$cov,
        pctM = df_mh$pctM
    )
    gr$sample <- sample_name
    gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_grs[[sample_name]] <- gr
}


#

##########################
# scMD cell type deconvolution
# Liftover T2T → hg19 (scMD references are hg19-based)
chain_t2t_hg38 <- import.chain("/users/mkelsey/data/ref/genomes/hs1/liftOver/hs1ToHg38.over.chain")
chain_hg38_hg19 <- import.chain("/users/mkelsey/data/ref/genomes/hg38/hg38ToHg19.over.chain")

bulk_list <- list()
for (sample_name in samples) {
    gr <- sample_grs[[sample_name]]
    # Filter for minimum coverage and autosomes only
    gr <- gr[gr$cov >= MINIMUMCOVERAGE & seqnames(gr) %in% autosomes]

    # Two-step liftover: T2T → hg38 → hg19
    n_t2t <- length(gr)
    gr_hg38 <- unlist(liftOver(gr, chain_t2t_hg38))
    n_hg38 <- length(gr_hg38)
    gr_hg19 <- unlist(liftOver(gr_hg38, chain_hg38_hg19))
    n_hg19 <- length(gr_hg19)
    cat(sprintf(
        "%s: T2T=%d → hg38=%d (lost %.1f%%) → hg19=%d (lost %.1f%%) | total lost: %.1f%%\n",
        sample_name, n_t2t, n_hg38, 100 * (1 - n_hg38 / n_t2t),
        n_hg19, 100 * (1 - n_hg19 / n_hg38),
        100 * (1 - n_hg19 / n_t2t)
    ))

    # Create TargetID with chr prefix (matches Tian reference)
    target_ids <- paste0(seqnames(gr_hg19), ":", start(gr_hg19))
    betas <- gr_hg19$pctM / 100
    bulk_list[[sample_name]] <- setNames(betas, target_ids)
}
# rm(sample_grs)

# Keep sites present in all samples
common_sites <- Reduce(intersect, lapply(bulk_list, names))
print(sprintf("Number of common CpG sites across all samples: %d", length(common_sites)))
bulk_mat <- do.call(cbind, lapply(bulk_list, function(x) x[common_sites]))
colnames(bulk_mat) <- names(bulk_list)
rownames(bulk_mat) <- common_sites
rm(bulk_list)

# scMD references use two TargetID formats:
#   Lee: "5:34189635" (no chr prefix)
#   Tian: "chr5:34189635" (with chr prefix)
# Provide both so intersect() matches against each reference
nochr_sites <- gsub("^chr", "", common_sites)
bulk_mat_both <- rbind(bulk_mat, bulk_mat)
rownames(bulk_mat_both) <- c(common_sites, nochr_sites)
rm(bulk_mat)

# If only 1 sample, add a dummy column to prevent R from dropping matrix dimensions
# inside scMD (single-column matrix gets coerced to vector on subset)
added_dummy <- FALSE
if (ncol(bulk_mat_both) == 1) {
    bulk_mat_both <- cbind(bulk_mat_both, dummy_ = bulk_mat_both[, 1])
    added_dummy <- TRUE
}

# Diagnostics: check overlap with scMD references and marker coverage
{
    e <- new.env()
    lazyLoad("/users/mkelsey/anaconda/scMD/lib/R/library/scMD/data/Rdata", envir = e)

    bulk_sites_chr <- common_sites
    bulk_sites_nochr <- gsub("^chr", "", common_sites)

    # Check overlap with each reference
    lee_sites <- rownames(e$Lee_sig_all_WGBS)
    tian_sites <- rownames(e$Tian_sig_all_WGBS)

    lee_overlap <- intersect(bulk_sites_nochr, lee_sites)
    tian_overlap <- intersect(bulk_sites_chr, tian_sites)

    cat(sprintf(
        "Lee reference: %d total signature sites, %d overlap with bulk (%.1f%%)\n",
        length(lee_sites), length(lee_overlap), 100 * length(lee_overlap) / length(lee_sites)
    ))
    cat(sprintf(
        "Tian reference: %d total signature sites, %d overlap with bulk (%.1f%%)\n",
        length(tian_sites), length(tian_overlap), 100 * length(tian_overlap) / length(tian_sites)
    ))

    # Check which marker CpGs are used after get_sig selects top 100 per cell type
    for (ref_name in c("Lee", "Tian")) {
        if (ref_name == "Lee") {
            sig <- e$Lee_sig_all_WGBS
            dm <- e$Lee_DF_WGBS
            overlap <- lee_overlap
        } else {
            sig <- e$Tian_sig_all_WGBS
            dm <- e$Tian_DF_WGBS
            overlap <- tian_overlap
        }
        # Reproduce get_sig logic: top 100 markers per cell type
        ct_ind <- colnames(sig)
        dm_filtered <- dm %>% filter(TargetID %in% rownames(sig))
        markers <- lapply(ct_ind, function(ct) {
            dm_filtered %>%
                arrange(.data[[ct]]) %>%
                dplyr::slice(1:100) %>%
                pull(TargetID)
        })
        names(markers) <- ct_ind
        all_markers <- unique(unlist(markers))
        markers_in_bulk <- intersect(all_markers, if (ref_name == "Lee") bulk_sites_nochr else bulk_sites_chr)

        cat(sprintf(
            "\n%s: %d total markers selected, %d present in bulk (%.1f%%)\n",
            ref_name, length(all_markers), length(markers_in_bulk),
            100 * length(markers_in_bulk) / length(all_markers)
        ))

        # Per cell type marker coverage
        for (ct in ct_ind) {
            ct_markers <- markers[[ct]]
            ct_in_bulk <- sum(ct_markers %in% (if (ref_name == "Lee") bulk_sites_nochr else bulk_sites_chr))
            cat(sprintf("  %s: %d/100 markers covered\n", ct, ct_in_bulk))
        }

        # Check mean coverage of marker sites vs non-marker sites in first sample
        first_sample <- samples[1]
        if (ref_name == "Lee") {
            marker_betas <- bulk_mat_both[bulk_sites_nochr %in% markers_in_bulk, first_sample]
        } else {
            marker_betas <- bulk_mat_both[bulk_sites_chr %in% markers_in_bulk, first_sample]
        }
        cat(sprintf(
            "  Marker beta stats in %s: mean=%.3f, median=%.3f, NAs=%d\n",
            first_sample, mean(marker_betas, na.rm = TRUE),
            median(marker_betas, na.rm = TRUE), sum(is.na(marker_betas))
        ))
    }

    # Test ±1 offset: bedMethyl is 0-based; scMD references may be 1-based
    # Use numeric arithmetic on pre-split parts to avoid slow regex on 30M sites
    bulk_chr_parts <- do.call(rbind, strsplit(bulk_sites_chr, ":", fixed = TRUE))
    bulk_pos <- as.numeric(bulk_chr_parts[, 2])
    bulk_prefix_chr <- paste0(bulk_chr_parts[, 1], ":")
    bulk_prefix_nochr <- gsub("^chr", "", bulk_prefix_chr)
    bulk_sites_chr_p1 <- paste0(bulk_prefix_chr, bulk_pos + 1)
    bulk_sites_chr_m1 <- paste0(bulk_prefix_chr, bulk_pos - 1)
    bulk_sites_nochr_p1 <- paste0(bulk_prefix_nochr, bulk_pos + 1)
    bulk_sites_nochr_m1 <- paste0(bulk_prefix_nochr, bulk_pos - 1)

    cat("\n=== OFFSET TEST (all autosomes) ===\n")
    cat(sprintf("Lee  exact:  %d / %d\n", length(intersect(bulk_sites_nochr, lee_sites)), length(lee_sites)))
    cat(sprintf("Lee  +1:     %d / %d\n", length(intersect(bulk_sites_nochr_p1, lee_sites)), length(lee_sites)))
    cat(sprintf("Lee  -1:     %d / %d\n", length(intersect(bulk_sites_nochr_m1, lee_sites)), length(lee_sites)))
    cat(sprintf("Tian exact:  %d / %d\n", length(intersect(bulk_sites_chr, tian_sites)), length(tian_sites)))
    cat(sprintf("Tian +1:     %d / %d\n", length(intersect(bulk_sites_chr_p1, tian_sites)), length(tian_sites)))
    cat(sprintf("Tian -1:     %d / %d\n", length(intersect(bulk_sites_chr_m1, tian_sites)), length(tian_sites)))

    # Check whether references are hg38 or hg19:
    # Sample 200 Tian positions on chr1 and check if they sit on CG dinucleotides in hg38
    hg38_fa <- "/users/mkelsey/data/ref/genomes/hg38/hg38.p13.sorted.fa"
    if (file.exists(hg38_fa)) {
        library(Rsamtools)
        hg38 <- FaFile(hg38_fa)
        test_sites <- head(tian_sites[grep("^chr1:", tian_sites)], 200)
        test_pos <- as.numeric(sub(".*:", "", test_sites))
        test_gr <- GRanges("chr1", IRanges(start = test_pos, width = 2))
        test_seq <- as.character(getSeq(hg38, test_gr))
        n_cg <- sum(test_seq == "CG")
        cat(sprintf("\n=== GENOME BUILD TEST ===\n"))
        cat(sprintf(
            "Tian chr1: %d/200 sites are CG dinucleotides in hg38 → %s\n",
            n_cg, ifelse(n_cg > 150, "LIKELY hg38", "NOT hg38")
        ))
        cat("First 10 sequences at Tian sites in hg38:\n")
        print(head(test_seq, 10))
    } else {
        cat("\nhg38 fasta not found at", hg38_fa, "- skipping genome build test\n")
    }

    rm(e)
}

# Run scMD with pre-filtered references so marker selection only considers sites in bulk
# This ensures all 100 markers per cell type are present in the data
dmet_list <- c("CIBERSORT", "EPIC", "FARDEEP", "DCQ", "ICeDT", "NNLS", "RPC")
data("Lee_7ct_WGBS")
data("Tian_7ct_WGBS")

bulk_sites <- rownames(bulk_mat_both)

# Filter DM_df to only include sites present in bulk
Lee_DF_filtered <- Lee_DF_WGBS %>% filter(TargetID %in% bulk_sites)
Tian_DF_filtered <- Tian_DF_WGBS %>% filter(TargetID %in% bulk_sites)
cat(sprintf(
    "Lee DM_df: %d → %d sites after filtering to bulk\n",
    nrow(Lee_DF_WGBS), nrow(Lee_DF_filtered)
))
cat(sprintf(
    "Tian DM_df: %d → %d sites after filtering to bulk\n",
    nrow(Tian_DF_WGBS), nrow(Tian_DF_filtered)
))

# Patch sc_MD_deconv to handle methods that fail internally
# (scMD bug: sapply over res_Ens crashes if a method returns no p_hat)
sc_MD_deconv_safe <- function(bulk, sc_mtx, DM_df, dmet_list, nmrk = 100, ncluster = 5) {
    overlap_features <- intersect(rownames(sc_mtx), rownames(bulk))
    sc_mtx <- sc_mtx[overlap_features, ]
    bulk <- bulk[overlap_features, ]
    celltype_ind <- colnames(sc_mtx)
    sig <- scMD:::get_sig(nmrk = nmrk, sc_mtx, DM_df, ct_ind = celltype_ind)
    bulk <- bulk[intersect(rownames(sig), rownames(bulk)), ]
    phat_all <- list()
    if ("NNLS" %in% dmet_list) {
        suppressMessages(frac_est <- MIND::est_frac(sig, bulk))
        phat_all <- append(phat_all, list(NNLS = frac_est))
    }
    if ("RPC" %in% dmet_list) {
        suppressMessages(frac_rpc <- EpiDISH::epidish(bulk, sig, method = "RPC")$estF)
        phat_all <- append(phat_all, list(RPC = frac_rpc))
    }
    dmet_list_ens <- setdiff(dmet_list, c("Houseman", "NNLS", "RPC"))
    ref_list_beta <- list(beta = list(
        ref_matrix = sig,
        meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
        data_name = "beta"
    ))
    suppressMessages(res_Ens <- EnsDeconv::gen_all_res_list(
        count_bulk = bulk,
        ref_list = ref_list_beta,
        params = EnsDeconv::get_params(
            data_type = "singlecell-rna", data_name = "beta",
            n_markers = 50, Marker.Method = "none",
            TNormalization = c("none"), CNormalization = c("none"),
            Scale = c("linear", "log"), dmeths = dmet_list_ens
        ),
        enableFileSaving = FALSE, parallel_comp = TRUE, ncore = ncluster
    ))
    ind <- sapply(res_Ens, function(x) tryCatch(length(x[["a"]][["p_hat"]][[1]]), error = function(e) 0))
    res_Ens <- res_Ens[which(ind == 1)]
    phat_list1 <- unlist(lapply(res_Ens, EnsDeconv:::getphat), recursive = FALSE)
    ref_list_Mval <- list(Mval = list(
        ref_matrix = scMD:::BetaToMvalue(sig),
        meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
        data_name = "Mval"
    ))
    suppressMessages(res_Ens_mval <- EnsDeconv::gen_all_res_list(
        count_bulk = scMD:::BetaToMvalue(bulk),
        ref_list = ref_list_Mval,
        params = EnsDeconv::get_params(
            data_type = "singlecell-rna", data_name = "Mval",
            n_markers = 50, Marker.Method = "none",
            TNormalization = c("none"), CNormalization = c("none"),
            Scale = c("linear"), dmeths = dmet_list_ens
        ),
        enableFileSaving = FALSE, parallel_comp = TRUE, ncore = ncluster
    ))
    ind <- sapply(res_Ens_mval, function(x) tryCatch(length(x[["a"]][["p_hat"]][[1]]), error = function(e) 0))
    res_Ens_mval <- res_Ens_mval[which(ind == 1)]
    phat_list2 <- unlist(lapply(res_Ens_mval, EnsDeconv:::getphat), recursive = FALSE)
    phat_all <- append(phat_all, append(phat_list1, phat_list2))
    return(phat_all)
}

Lee_res <- sc_MD_deconv_safe(
    bulk = bulk_mat_both,
    sc_mtx = Lee_sig_all_WGBS,
    DM_df = Lee_DF_filtered, dmet_list = dmet_list, nmrk = 100
)
Tian_res <- sc_MD_deconv_safe(bulk_mat_both,
    sc_mtx = Tian_sig_all_WGBS,
    DM_df = Tian_DF_filtered, dmet_list = dmet_list, nmrk = 100
)

# Ensemble across methods, excluding outlier methods
phat_all <- c(Lee_res, Tian_res)
saveRDS(phat_all, "ldna/results/m/tables/phat_all.rds")
exclude_methods <- c("EPIC_Mval_linear", "DCQ_beta_log")
keep_idx <- !grepl(paste(exclude_methods, collapse = "|"), names(phat_all))

cat(sprintf(
    "Excluding %d/%d method variants: %s\n",
    sum(!keep_idx), length(keep_idx), paste(exclude_methods, collapse = ", ")
))

lee_keep <- !grepl(paste(exclude_methods, collapse = "|"), names(Lee_res))
tian_keep <- !grepl(paste(exclude_methods, collapse = "|"), names(Tian_res))

make_fractions <- function(phat_list) {
    recomp <- scMD:::CTS_EnsDeconv_wrapper(phat_all = phat_list)
    ctf <- as.data.frame(recomp$ensemble_p)
    if (added_dummy) {
        ctf <- ctf[!grepl("^dummy_", rownames(ctf)), , drop = FALSE]
    }
    ctf$sample_name <- rownames(ctf)
    ctf %>% left_join(sample_table %>% select(sample_name, condition), by = "sample_name")
}

cell_type_fractions <- make_fractions(phat_all[keep_idx])
cell_type_fractions_lee <- make_fractions(Lee_res[lee_keep])
cell_type_fractions_tian <- make_fractions(Tian_res[tian_keep])

dir.create("ldna/results/m/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(cell_type_fractions, "ldna/results/m/tables/scMD_cell_type_fractions.csv")
write_csv(cell_type_fractions_lee, "ldna/results/m/tables/scMD_cell_type_fractions_Lee.csv")
write_csv(cell_type_fractions_tian, "ldna/results/m/tables/scMD_cell_type_fractions_Tian.csv")
# cell_type_fractions <- read_csv("ldna/results/m/tables/scMD_cell_type_fractions.csv")
# cell_type_fractions_lee <- read_csv("ldna/results/m/tables/scMD_cell_type_fractions_Lee.csv")
# cell_type_fractions_tian <- read_csv("ldna/results/m/tables/scMD_cell_type_fractions_Tian.csv")

print("scMD cell type deconvolution complete")
print(cell_type_fractions)


# Plotting helper
cell_types <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")
dir.create("ldna/results/m/plots/scMD", recursive = TRUE, showWarnings = FALSE)

plot_scmd <- function(ctf, label) {
    ctl <- ctf %>%
        pivot_longer(cols = all_of(cell_types), names_to = "cell_type", values_to = "fraction")

    p <- ggplot(ctl, aes(x = cell_type, y = fraction, fill = condition)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1.5, alpha = 0.7) +
        labs(x = "", y = "Estimated Fraction", title = sprintf("scMD Cell Type Deconvolution (%s)", label)) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        )
    ggsave(sprintf("ldna/results/m/plots/scMD/cell_type_fractions_%s.pdf", label), p, width = 6, height = 5)

    # Statistical tests: cell type fractions by condition
    cat(sprintf("\n=== Cell type fraction stats (%s) ===\n", label))
    ctf_test <- ctf %>% mutate(Neuron = Inh + Exc)
    ctf_test$condition <- relevel(factor(ctf_test$condition), ref = "YNG")
    test_types <- c(cell_types, "Neuron")
    stats_results <- list()
    for (ct in test_types) {
        fit <- aov(as.formula(paste(ct, "~ condition")), data = ctf_test)
        aov_p <- summary(fit)[[1]][["Pr(>F)"]][1]
        # Pairwise vs YNG (Dunnett-like using linear model contrasts)
        lm_fit <- lm(as.formula(paste(ct, "~ condition")), data = ctf_test)
        coefs <- summary(lm_fit)$coefficients
        pairwise_rows <- grep("^condition", rownames(coefs))
        for (r in pairwise_rows) {
            contrast_name <- gsub("condition", "", rownames(coefs)[r])
            stats_results[[length(stats_results) + 1]] <- data.frame(
                cell_type = ct,
                contrast = sprintf("%s vs YNG", contrast_name),
                estimate = coefs[r, "Estimate"],
                p_value = coefs[r, "Pr(>|t|)"]
            )
        }
        stats_results[[length(stats_results) + 1]] <- data.frame(
            cell_type = ct, contrast = "ANOVA (overall)", estimate = NA, p_value = aov_p
        )
    }
    stats_df <- bind_rows(stats_results)
    cat(sprintf("\n"))
    print(stats_df %>% arrange(p_value))
    write_csv(stats_df, sprintf("ldna/results/m/tables/scMD_celltype_stats_%s.csv", label))

    p <- ggplot(ctl, aes(x = sample_name, y = fraction, fill = cell_type)) +
        geom_col() +
        labs(x = "", y = "Estimated Fraction", title = sprintf("scMD Cell Type Composition (%s)", label)) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) +
        facet_wrap(~condition, scales = "free_x")
    ggsave(sprintf("ldna/results/m/plots/scMD/cell_type_composition_stacked_%s.pdf", label), p, width = 10, height = 5)

    p <- ctl %>%
        filter(cell_type %in% c("Exc", "Inh", "Oligo")) %>%
        pivot_wider(names_from = cell_type, values_from = fraction) %>%
        mutate(Neuron = Exc + Inh) %>%
        ggplot(aes(x = Oligo, y = Neuron, color = condition)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "grey50") +
        labs(
            x = "Oligodendrocyte Fraction", y = "Neuron Fraction (Exc + Inh)",
            title = sprintf("Neuron vs Oligodendrocyte (%s)", label)
        ) +
        theme_minimal() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
    ggsave(sprintf("ldna/results/m/plots/scMD/neuron_vs_oligo_%s.pdf", label), p, width = 6, height = 5)
}

plot_scmd(cell_type_fractions, "ensemble")
plot_scmd(cell_type_fractions_lee, "Lee")
plot_scmd(cell_type_fractions_tian, "Tian")

# Neuron fraction by condition: plot and test
plot_neu_by_condition <- function(ctf, label) {
    ctf_neu <- ctf %>% mutate(Neuron = Inh + Exc)
    conditions <- unique(ctf_neu$condition)

    p <- ggplot(ctf_neu, aes(x = condition, y = Neuron, fill = condition)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
        labs(
            x = "", y = "Neuron Fraction (Exc + Inh)",
            title = sprintf("Neuron Fraction by Condition (%s)", label)
        ) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            legend.position = "none"
        )
    ggsave(sprintf("ldna/results/m/plots/scMD/neuron_by_condition_%s.pdf", label), p, width = 5, height = 5)

    # Statistical tests
    cat(sprintf("\n=== Neuron fraction by condition (%s) ===\n", label))
    cat("Group means:\n")
    print(ctf_neu %>% group_by(condition) %>% summarise(mean_Neu = mean(Neuron), sd_Neu = sd(Neuron), n = n()))
    fit <- aov(Neuron ~ condition, data = ctf_neu)
    cat("ANOVA:\n")
    print(summary(fit))
    if (length(conditions) > 2) {
        cat("Tukey HSD:\n")
        print(TukeyHSD(fit))
    }
}

plot_neu_by_condition(cell_type_fractions, "ensemble")
plot_neu_by_condition(cell_type_fractions_lee, "Lee")
plot_neu_by_condition(cell_type_fractions_tian, "Tian")

# Per-method stacked plots split by reference
phat <- phat_all
n_methods <- length(phat)
n_per_ref <- n_methods / 2
method_dfs <- list()
for (i in seq_along(phat)) {
    ref <- if (i <= n_per_ref) "Lee" else "Tian"
    mat <- phat[[i]]
    method_name <- names(phat)[i]
    # Extract short method name
    short_name <- gsub(".*TRUE_([A-Za-z]+)_.*", "\\1", method_name)
    if (short_name == method_name) short_name <- method_name # fallback for NNLS
    # Extract scale (beta vs Mval)
    scale <- ifelse(grepl("Mval", method_name), "Mval", "beta")
    if (short_name == "NNLS") scale <- "beta"
    # Extract linear vs log
    lin <- ifelse(grepl("_log_", method_name), "log", "linear")
    if (short_name == "NNLS") lin <- ""
    label <- paste0(short_name, "_", scale, ifelse(lin != "", paste0("_", lin), ""))

    df <- as.data.frame(mat)
    df$sample_name <- rownames(mat)
    df$method <- label
    df$reference <- ref
    method_dfs[[i]] <- df
}
all_methods_df <- bind_rows(method_dfs) %>%
    left_join(sample_table %>% select(sample_name, condition), by = "sample_name")

cell_types <- c("Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc")
all_methods_long <- all_methods_df %>%
    filter(!(method %in% c("DCQ_beta_log", "EPIC_Mval_linear"))) %>%
    pivot_longer(cols = any_of(cell_types), names_to = "cell_type", values_to = "fraction")

for (ref in c("Lee", "Tian")) {
    p <- all_methods_long %>%
        filter(reference == ref) %>%
        ggplot(aes(x = sample_name, y = fraction, fill = cell_type)) +
        geom_col() +
        facet_wrap(~method, ncol = 2) +
        labs(x = "", y = "Estimated Fraction", title = sprintf("scMD Per-Method Estimates (%s reference)", ref)) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)
        )
    ggsave(sprintf("ldna/results/m/plots/scMD/cell_type_by_method_%s2.pdf", ref), p, width = 14, height = 16)
}
