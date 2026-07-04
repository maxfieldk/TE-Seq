outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
subfam <- "L1HS"
consensus_path <- sprintf("%s/alignments/%s_fl_consensus.fa", outputdir_meth_clustering, subfam)
consensus_ss <- readDNAStringSet(consensus_path)
cg_indices <- consensus_ss %>%
    vmatchPattern(pattern = "CG") %>%
    start() %>%
    unlist() %>%
    as.numeric()

mod_code_var <- "m"
regions_of_interest <- list(c(0, 145), c(0, 328), c(0, 500), c(0, 909))
required_fraction_of_total_cg <- 0.75
meth_bins <- c(0.5, 0.7, 0.9)
context <- "CpG"
region <- "L1HS_FL"
outputdirtables <- sprintf("ldna/results/%s/tables/reads_custom/%s_%s", mod_code_var, region, required_fraction_of_total_cg)
dir.create(outputdirtables, recursive = TRUE)

breaks <- c(0, meth_bins, 1)
labels <- paste0("[", head(breaks, -1), ", ", tail(breaks, -1), ")")
labels[length(labels)] <- sub("\\)$", "]", labels[length(labels)])


readsdf1 <- readscg %>%
    left_join(rmann %>%
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


by_read %>%
    filter(subset == "0to145") %$% num_cpgs_in_read %>%
    table()


mdf <- by_read %>%
    # filter(num_cpgs_in_read == 8) %>%
    filter(subset == "0to145") %>%
    group_by(gene_id, condition) %>%
    summarise(fm = mean(fraction_meth, na.rm = TRUE)) %>%
    pivot_wider(names_from = condition, values_from = fm) %>%
    ungroup() %>%
    mutate(dif = SEN - PRO) %>%
    filter(is.finite(dif)) %>%
    mutate(dif_quantile = ecdf(dif)(dif))
mdf %>%
    arrange(dif) %>%
    print(n = 20)



mdf <- by_read %>%
    # filter(num_cpgs_in_read == 8) %>%
    filter(subset == "0to145") %>%
    mutate(lt15 = ifelse(fraction_meth < 0.15, 1, 0)) %>%
    group_by(gene_id, condition) %>%
    summarise(lt15 = mean(lt15)) %>%
    pivot_wider(names_from = condition, values_from = lt15) %>%
    mutate(dif = SEN - PRO)
mdf %>%
    left_join(rmann) %>%
    filter(intactness_req == "Intact") %>%
    arrange(-dif)
mdf %>%
    filter(PRO == 0) %>%
    arrange(-dif)

mdf %>% filter(gene_id == "L1HS_7p14.3_1")
