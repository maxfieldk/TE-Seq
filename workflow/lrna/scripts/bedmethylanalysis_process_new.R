module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
set.seed(123)

library(rtracklayer)
library(Biostrings)
library(cowplot)
library(pryr)
library(magrittr)
library(forcats)
library(ComplexHeatmap)
library(GenomicRanges)
library(configr)
library(ggbeeswarm)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

conf_lr <- confALL[["lrna"]]
conf_sr <- confALL[["srna"]]
sample_table_lrna <- read_csv("conf/sample_table_lrna.csv")
sample_table_srna <- read_csv("conf/sample_table_srna.csv")

{
    genome_lengths <- fasta.seqlengths(conf$reference)
    genome_fa <- readDNAStringSet(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("^NI_", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE) %>% str_sort(numeric = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE) %>% str_sort(numeric = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")
    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
    CHROMOSOMESINCLUDEDINANALYSIS <- c(autosomes, chrX, chrY, nonrefchromosomes)
    CHROMOSOMESINCLUDEDINANALYSIS_REF <- c(autosomes, chrX, chrY)
}

#################### Modification code mapping
MOD_CODE_MAP <- tibble::tibble(
    mod_code = c("a", "m", "17596", "17802", "19227", "19228", "69426"),
    mod_name = c("m6A", "m5C", "Am", "pseU", "Um", "Cm", "inosine"),
    mod_long = c(
        "N6-methyladenosine",
        "5-methylcytosine",
        "2'-O-methyladenosine",
        "pseudouridine",
        "2'-O-methyluridine",
        "2'-O-methylcytosine",
        "inosine"
    )
)

MOD_CANONICAL_BASE <- c(
    m6A = "A", m5C = "C", Am = "A", pseU = "T", Um = "T", Cm = "C", inosine = "A"
)

#################### functions and themes
tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        mod_name_val <- "m6A"
        regions_val <- unique(c(conf$rte_subfamily_read_level_analysis, conf$rte_subfamily_extra_modification_analysis))
        assign("params", list(
            mod_name = mod_name_val,
            regions = regions_val
        ), env = globalenv())
        assign("inputs", list(
            bedmethylpaths = sprintf("lrna/intermediates/%s/methylation/bymod/%s_%s_bedMethyl.bed", samples, samples, mod_name_val),
            read_mods = sprintf("lrna/intermediates/%s/methylation/bymod/%s_readmods_NoContext_%s_%s.tsv", samples, samples, conf$rte_subfamily_read_level_analysis, mod_name_val),
            region_bedmethylpaths = as.vector(outer(
                sprintf("lrna/intermediates/%s/methylation/bymod/%s", samples, samples),
                sprintf("_%s_%s_bedMethyl.bed", regions_val, mod_name_val),
                paste0
            ))
        ), env = globalenv())
        assign("outputs", list(grsdf = sprintf("lrna/Rintermediates/%s/grsdf.tsv", mod_name_val)), env = globalenv())
    }
)

current_mod_name <- params$mod_name
cat(sprintf("\n=== Processing modification: %s ===\n", current_mod_name))

ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

# //ANCHOR - DATA FIRST PROCESS

if (!file.exists(sprintf("lrna/Rintermediates/%s/grsdf.tsv", current_mod_name))) {
    ##########################
    # PREP DATA FOR ANALYSIS
    # Input files are already split by modification type

    dir.create(sprintf("lrna/Rintermediates/%s", current_mod_name), recursive = TRUE, showWarnings = FALSE)
    dir.create(sprintf("lrna/results/%s/plots/rte", current_mod_name), recursive = TRUE, showWarnings = FALSE)

    sample_grs <- list()
    for (sample_name in samples) {
        cat(sprintf("Loading bedMethyl for %s (%s)...\n", sample_name, current_mod_name))
        filepath <- grep(sprintf("/%s/", sample_name), inputs$bedmethylpaths, value = TRUE)
        # Skip if file is empty (mod not present for this sample)
        if (file.size(filepath) == 0) next
        df <- read_table(filepath, col_names = FALSE)
        if (nrow(df) == 0) next
        gr <- GRanges(
            seqnames = df$X1,
            ranges = IRanges(start = df$X2, end = df$X2),
            cov = df$X10,
            pctM = as.double(df$X11)
        )
        gr$sample <- sample_name
        gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
        sample_grs[[sample_name]] <- gr
        rm(df, gr)
        gc()
    }

    grs <- Reduce(c, sample_grs)
    rm(sample_grs)

    grsdf <- tibble(as.data.frame(grs))
    write_delim(grsdf %>% filter(grepl("^NI_", seqnames)), sprintf("lrna/Rintermediates/%s/grsdf_nonref.tsv", current_mod_name), col_names = TRUE)
    grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
    pos <- paste0(grsdf$seqnames, "_", grsdf$start, "_", grsdf$end)
    grsdf$pos <- pos

    write_delim(grsdf, sprintf("lrna/Rintermediates/%s/grsdf.tsv", current_mod_name), col_names = TRUE)



    ####################
    ## RTEs

    ## Load Data and add annotations
    r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
    r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
    rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
    rm(r_annotation_fragmentsjoined)
    RM <- GRanges(rmann)

    ### ONTOLOGY DEFINITION
    {
        annot_colnames <- colnames(r_repeatmasker_annotation)
        annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
        ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "_.*family")]
        small_ontologies <- ontologies[grepl("subfamily", ontologies)]

        big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
        big_ontology_groups <- c()
        for (ontology in big_ontologies) {
            big_ontology_groups <- c(big_ontology_groups, rmann %>%
                pull(!!sym(ontology)) %>%
                unique())
        }
        big_ontology_groups <- big_ontology_groups %>% unique()

        modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
        region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
        element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
    }

    RMdf <- rmann %>% filter(seqnames %in% CHROMOSOMESINCLUDEDINANALYSIS)
    write_delim(RMdf, sprintf("lrna/Rintermediates/%s/RMdf.tsv", current_mod_name), col_names = TRUE)

    merge_with_grs <- function(grs, rte_frame) {
        mbo <- mergeByOverlaps(grs, rte_frame)
        methdf <- mbo$grs %>%
            as.data.frame() %>%
            tibble()
        rte_only_frame <- mbo$rte_frame %>%
            as.data.frame() %>%
            tibble() %>%
            dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
        rtedf_merged <- bind_cols(methdf, rte_only_frame)
        return(rtedf_merged)
    }

    if ("width" %in% colnames(RMdf)) {
        joinframe <- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
    } else {
        joinframe <- RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand)
    }

    # All RTEs
    rte_frame <- GRanges(RMdf)
    alltedf <- merge_with_grs(grs, rte_frame)
    write_delim(alltedf, sprintf("lrna/Rintermediates/%s/alltedf.tsv", current_mod_name), col_names = TRUE)

    # RTEs excluding "Other"
    rte_frame <- GRanges(RMdf %>% filter(rte_subfamily != "Other"))
    grs <- GRanges(grsdf)
    rtedf <- merge_with_grs(grs, rte_frame)
    write_delim(rtedf, sprintf("lrna/Rintermediates/%s/rtedf.tsv", current_mod_name), col_names = TRUE)

    perelementdf <- rtedf %>%
        group_by(gene_id, sample, condition) %>%
        summarize(mean_meth = mean(pctM)) %>%
        left_join(joinframe)
    perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))
    write_delim(perelementdf, sprintf("lrna/Rintermediates/%s/perelementdf.tsv", current_mod_name), col_names = TRUE)

    # FL RTE body modification analysis
    rte_frame_fl <- GRanges(RMdf %>% filter(rte_subfamily != "Other") %>% filter(rte_length_req == "FL"))
    rtedf_fl <- merge_with_grs(grs, rte_frame_fl)
    write_delim(rtedf_fl, sprintf("lrna/Rintermediates/%s/rtedf_fl.tsv", current_mod_name), col_names = TRUE)

    perelementdf_fl <- rtedf_fl %>%
        group_by(gene_id, sample, condition) %>%
        summarize(mean_meth = mean(pctM)) %>%
        left_join(joinframe)
    perelementdf_fl <- perelementdf_fl %>% filter(!is.na(rte_length_req))
    write_delim(perelementdf_fl, sprintf("lrna/Rintermediates/%s/perelementdf_fl.tsv", current_mod_name), col_names = TRUE)

    # Gene body modification analysis
    refseq_gr <- import(conf$annotation_genes)
    genes_gr <- refseq_gr[mcols(refseq_gr)[, "type"] == "gene", ]
    genes_gr <- genes_gr[seqnames(genes_gr) %in% CHROMOSOMESINCLUDEDINANALYSIS, ]
    rm(refseq_gr)

    mbo <- mergeByOverlaps(grs, genes_gr)
    gene_methdf <- mbo$grs %>%
        as.data.frame() %>%
        tibble()
    gene_methdf$gene_id <- mbo$genes_gr %>%
        as.data.frame() %>%
        tibble() %$% gene_id
    write_delim(gene_methdf, sprintf("lrna/Rintermediates/%s/gene_body_modification.tsv", current_mod_name), col_names = TRUE)

    # rm(grs, alltedf, rtedf, perelementdf, rtedf_fl, perelementdf_fl, gene_methdf)


    ############## Read-level modification analysis
    dir.create(sprintf("lrna/results/%s/plots/reads", current_mod_name), recursive = TRUE, showWarnings = FALSE)

    readslist <- list()
    for (region in conf$rte_subfamily_read_level_analysis) {
        for (sample_name in samples) {
            filepath <- grep(region,
                grep(sprintf("/%s/", sample_name),
                    inputs$read_mods,
                    value = TRUE
                ),
                value = TRUE
            )
            if (length(filepath) == 0) next
            if (file.size(filepath) == 0) next
            df <- read_delim(filepath)
            if (nrow(df) == 0) next
            df$region <- region
            df$sample <- sample_name
            df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
            grsx <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
            eoi <- import(paste0("aref/default/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
            strand(eoi) <- "*"
            mbo <- mergeByOverlaps(grsx, eoi)
            df1 <- as.data.frame(mbo) %>%
                tibble() %>%
                dplyr::select(starts_with("grsx"), name) %>%
                dplyr::rename(gene_id = name)
            colnames(df1) <- gsub("grsx.", "", colnames(df1))
            readslist <- c(readslist, list(df1))
        }
    }
    if (length(readslist) > 0) {
        reads <- Reduce(rbind, readslist)
        write_delim(reads, sprintf("lrna/Rintermediates/%s/reads_context_all.tsv", current_mod_name), col_names = TRUE)
    }
} else {
    grsdf <- read_delim(sprintf("lrna/Rintermediates/%s/grsdf.tsv", current_mod_name), col_names = TRUE)
    r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
    rmann <- left_join(read_csv(conf$r_annotation_fragmentsjoined), r_repeatmasker_annotation)
    # alltedf <- read_delim(sprintf("lrna/Rintermediates/%s/alltedf.tsv", current_mod_name), col_names = TRUE)
    rtedf <- read_delim(sprintf("lrna/Rintermediates/%s/rtedf.tsv", current_mod_name), col_names = TRUE)
    # perelementdf <- read_delim(sprintf("lrna/Rintermediates/%s/perelementdf.tsv", current_mod_name), col_names = TRUE)
    rtedf_fl <- read_delim(sprintf("lrna/Rintermediates/%s/rtedf_fl.tsv", current_mod_name), col_names = TRUE)
    perelementdf_fl <- read_delim(sprintf("lrna/Rintermediates/%s/perelementdf_fl.tsv", current_mod_name), col_names = TRUE)
    # gene_methdf <- read_delim(sprintf("lrna/Rintermediates/%s/gene_body_modification.tsv", current_mod_name), col_names = TRUE)
    reads_path <- sprintf("lrna/Rintermediates/%s/reads_context_all.tsv", current_mod_name)
    if (file.exists(reads_path)) {
        reads <- read_delim(reads_path, col_names = TRUE)
    }
}

############################

# //ANCHOR - PAN-RTE SUBFAMILY SUMMARY
dir.create(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily"), recursive = TRUE, showWarnings = FALSE)

# Per-element mean modification across all sites
per_element <- rtedf_fl %>%
    filter(cov >= 5) %>%
    group_by(gene_id, sample, condition, rte_superfamily, rte_family, rte_subfamily) %>%
    summarise(
        mean_pctM = mean(pctM, na.rm = TRUE),
        n_sites = n(),
        n_high = sum(pctM >= 50),
        pct_high = 100 * n_high / n_sites,
        .groups = "drop"
    )

# Per-subfamily summary
per_subfamily <- per_element %>%
    group_by(rte_superfamily, rte_family, rte_subfamily, condition) %>%
    summarise(
        n_elements = n_distinct(gene_id),
        mean_mod = mean(mean_pctM, na.rm = TRUE),
        median_mod = median(mean_pctM, na.rm = TRUE),
        mean_n_sites = mean(n_sites),
        mean_pct_high = mean(pct_high, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(n_elements >= 5)

# 1. Barplot: mean modification by subfamily, colored by superfamily
p <- per_subfamily %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, mean_mod)) %>%
    ggplot(aes(x = rte_subfamily, y = mean_mod, fill = rte_superfamily)) +
    geom_col(position = "dodge") +
    facet_wrap(~condition) +
    coord_flip() +
    labs(x = NULL, y = str_glue("Mean % {current_mod_name}"), fill = "Superfamily") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/mean_mod_by_subfamily.pdf"),
    w = 10, h = max(4, nrow(per_subfamily %>% distinct(rte_subfamily)) * 0.3)
)

# 2. Violin/boxplot: per-element modification distribution by superfamily
p <- per_element %>%
    ggplot(aes(x = rte_superfamily, y = mean_pctM, fill = rte_superfamily)) +
    geom_violin(alpha = 0.3, scale = "width") +
    geom_boxplot(width = 0.2, outlier.size = 0.3) +
    facet_wrap(~condition) +
    labs(x = NULL, y = str_glue("Mean % {current_mod_name} per element")) +
    mtopen +
    theme(legend.position = "none")
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/mod_distribution_by_superfamily.pdf"), w = 8, h = 5)

# 3. Violin/boxplot: per-element modification by top subfamilies (top 15 by element count)
top_subfams <- per_element %>%
    count(rte_subfamily) %>%
    slice_max(n, n = 15) %>%
    pull(rte_subfamily)

p <- per_element %>%
    filter(rte_subfamily %in% top_subfams) %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, mean_pctM, .fun = median)) %>%
    ggplot(aes(x = rte_subfamily, y = mean_pctM, fill = rte_superfamily)) +
    geom_violin(alpha = 0.3, scale = "width") +
    geom_boxplot(width = 0.2, outlier.size = 0.3) +
    facet_wrap(~condition) +
    coord_flip() +
    labs(x = NULL, y = str_glue("Mean % {current_mod_name} per element"), fill = "Superfamily") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/mod_distribution_top_subfamilies.pdf"), w = 10, h = 8)

# 4. Scatter: mean sites per element vs mean modification (bubble = n_elements)
p <- per_subfamily %>%
    ggplot(aes(x = mean_n_sites, y = mean_mod, color = rte_superfamily, size = n_elements)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(aes(label = ifelse(n_elements >= 20, rte_subfamily, NA)), size = 2.5, max.overlaps = 15, na.rm = TRUE) +
    facet_wrap(~condition) +
    scale_size_continuous(range = c(1, 8)) +
    labs(
        x = "Mean # modified sites per element", y = str_glue("Mean % {current_mod_name}"),
        color = "Superfamily", size = "# elements"
    ) +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/sites_vs_mod_bubble.pdf"), w = 10, h = 6)

# 5. Heatmap: % of highly modified sites (>=50%) per subfamily × condition
p <- per_subfamily %>%
    filter(n_elements >= 10) %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, mean_pct_high)) %>%
    ggplot(aes(x = condition, y = rte_subfamily, fill = mean_pct_high)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.0f%%", mean_pct_high)), size = 2.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 25, name = "% sites\n>=50% mod") +
    facet_wrap(~rte_superfamily, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    mtopen +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/pct_high_sites_heatmap.pdf"),
    w = 8, h = max(4, nrow(per_subfamily %>% filter(n_elements >= 10) %>% distinct(rte_subfamily)) * 0.25)
)

# 6. Element count and coverage summary table
p <- per_subfamily %>%
    filter(n_elements >= 10) %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, n_elements)) %>%
    ggplot(aes(x = rte_subfamily, y = n_elements, fill = rte_superfamily)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "# FL elements with modification data", fill = "Superfamily") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/element_counts.pdf"),
    w = 8, h = max(4, nrow(per_subfamily %>% filter(n_elements >= 10) %>% distinct(rte_subfamily)) * 0.25)
)

# 7. Condition comparison: paired dot plot per subfamily
p <- per_subfamily %>%
    filter(n_elements >= 10) %>%
    dplyr::select(rte_subfamily, rte_superfamily, condition, mean_mod) %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, mean_mod)) %>%
    ggplot(aes(x = mean_mod, y = rte_subfamily, color = condition)) +
    geom_line(aes(group = rte_subfamily), color = "grey70") +
    geom_point(size = 2.5) +
    facet_wrap(~rte_superfamily, scales = "free_y") +
    labs(x = str_glue("Mean % {current_mod_name}"), y = NULL, color = "Condition") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/condition_comparison_dotplot.pdf"),
    w = 10, h = max(4, nrow(per_subfamily %>% filter(n_elements >= 10) %>% distinct(rte_subfamily)) * 0.25)
)

# 8. Robustly modified sites: positions >50% modified in >=50% of elements
# Robustly modified sites: use element-relative positions (binned to 1bp resolution)
# Compute relative position within each element, then bin across elements
site_by_element <- rtedf_fl %>%
    filter(cov >= 5) %>%
    mutate(rel_pos = ifelse(rte_strand == "+",
        start - rte_start,
        rte_end - start
    )) %>%
    group_by(rte_superfamily, rte_subfamily, condition, gene_id, rel_pos) %>%
    summarise(site_pctM = mean(pctM), .groups = "drop")

# For each relative position within a subfamily, count elements and % highly modified
site_robust <- site_by_element %>%
    mutate(is_high = site_pctM >= 50) %>%
    group_by(rte_superfamily, rte_subfamily, condition, rel_pos) %>%
    summarise(n_elements = n(), n_high = sum(is_high), pct_high = 100 * n_high / n_elements, .groups = "drop") %>%
    filter(n_elements >= 5, pct_high >= 50)

# Count robust sites per subfamily
robust_counts <- site_robust %>%
    group_by(rte_superfamily, rte_subfamily, condition) %>%
    summarise(n_robust_sites = n(), .groups = "drop")

# Get mean element length per subfamily
mean_lengths <- rmann %>%
    filter(rte_length_req == "FL") %>%
    group_by(rte_subfamily) %>%
    summarise(mean_length = mean(as.numeric(end - start), na.rm = TRUE), .groups = "drop")

robust_counts <- robust_counts %>%
    left_join(mean_lengths, by = "rte_subfamily") %>%
    mutate(robust_per_kb = 1000 * n_robust_sites / mean_length) %>%
    filter(!is.na(mean_length))

# Raw count
p <- robust_counts %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, n_robust_sites)) %>%
    ggplot(aes(x = rte_subfamily, y = n_robust_sites, fill = rte_superfamily)) +
    geom_col(position = "dodge") +
    facet_wrap(~condition) +
    coord_flip() +
    labs(x = NULL, y = str_glue("# sites >50% {current_mod_name} in >=50% of elements"), fill = "Superfamily") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/robust_sites_raw_count.pdf"),
    w = 10, h = max(4, n_distinct(robust_counts$rte_subfamily) * 0.3)
)

# Normalized to mean element length
p <- robust_counts %>%
    mutate(rte_subfamily = fct_reorder(rte_subfamily, robust_per_kb)) %>%
    ggplot(aes(x = rte_subfamily, y = robust_per_kb, fill = rte_superfamily)) +
    geom_col(position = "dodge") +
    facet_wrap(~condition) +
    coord_flip() +
    labs(x = NULL, y = str_glue("Robust {current_mod_name} sites per kb"), fill = "Superfamily") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/pan_subfamily/robust_sites_per_kb.pdf"),
    w = 10, h = max(4, n_distinct(robust_counts$rte_subfamily) * 0.3)
)

############################

# //ANCHOR - DATA HARMONIZATION

# yl1 <- rtedf_fl %>%
#     filter(rte_family == "L1") %>%
#     filter(rte_subfamily != "Other")
yl1_alllength <- rtedf %>%
    filter(rte_family == "L1") %>%
    filter(rte_subfamily != "Other")

outputdir_meth_clustering <- "ldna/results/m/plots/l1_alignment_meth"
mapping_files <- list.files(outputdir_meth_clustering, pattern = "_fl_mapping_to_consensus_table\\.csv$", full.names = TRUE)
# consensus_index_long <- map_dfr(mapping_files, read_csv)
all_combined <- read_csv(sprintf("ldna/results/m/plots/l1_alignment/all_subfamilies_mapping_to_consensus_harmonized.csv"))
consensus_index_long <- all_combined

l1hs_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1HS_consensus.fa")
l1pa2_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1PA2_consensus.fa")
l1pa3_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1PA3_consensus.fa")
l1pa4_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1PA4_consensus.fa")
l1pa5_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1PA5_consensus.fa")
l1pa6_cs <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/L1PA6_consensus.fa")
aligned <- readDNAStringSet("ldna/results/m/plots/l1_alignment/consensus/all_subfamilies_aligned.fa")

# Get all alignment columns (harmonized positions) where any subfamily has an A
modbase_indices <- which(Reduce("|", lapply(as.character(aligned), function(s) strsplit(s, "")[[1]] == unname(MOD_CANONICAL_BASE[current_mod_name]))))

# Add harmonized_pos to modbase_positions_df by detecting subfamily from gene_id
# Then filter to harmonized positions where any subfamily has an A
modbase_positions_df <- consensus_index_long %>%
    mutate(subfamily_consensus = str_extract(gene_id, "^[^_]+")) %>%
    filter(harmonized_pos %in% modbase_indices) %>%
    distinct(gene_id, sequence_pos, .keep_all = TRUE) %>%
    left_join(rmann %>% mutate(element_start_stranded = case_when(strand == "+" ~ start, strand == "-" ~ end)) %>%
        dplyr::select(gene_id, element_start_stranded, strand)) %>%
    mutate(genome_pos = case_when(
        strand == "+" ~ element_start_stranded + sequence_pos - 2,
        strand == "-" ~ element_start_stranded - sequence_pos
    )) %>%
    dplyr::select(-element_start_stranded, -strand)

# yl1ann <- yl1 %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end)))
# yl1annm <- left_join(modbase_positions_df, yl1ann, by = c("gene_id", "sequence_pos")) %>% filter(!is.na(seqnames))
# yl1annm %>% write_delim("lrna/Rintermediates/yl1annm.tsv")

yl1ann_alllength <- yl1_alllength %>% mutate(sequence_pos = ifelse(rte_strand == "+", (start - rte_start) + 2, (rte_end - end)))
yl1annm_alllength <- left_join(modbase_positions_df, yl1ann_alllength, by = c("gene_id", "sequence_pos")) %>% filter(!is.na(seqnames))
yl1annm_alllength %>% write_delim("lrna/Rintermediates/yl1annm_alllength.tsv")
yl1ann <- yl1annm_alllength %>% filter(rte_length_req == "FL")

tas <- unique(consensus_index_long %$% gene_id)
tas %>% grep("L1HS", .)

# //ANCHOR - POSITION PLOTS

aaa <- yl1annm %>%
    group_by(harmonized_pos, rte_subfamily) %>%
    summarise(sc = sum(cov), n = n(), mm = mean(pctM), .groups = "drop")

# Principled filtering: require positions to be in top quartile of observation count
# and have enough total coverage to be reliable
n_samples <- length(samples)
min_n <- quantile(aaa$n, 0.25) # at least 25th percentile of observations
min_cov <- quantile(aaa$sc, 0.25) # at least 25th percentile of total coverage
min_pctM <- 10 # minimum mean modification % to be considered "modified"

cat(sprintf("Position filters: min_n=%d (Q25), min_cov=%d (Q25), min_pctM=%d\n", min_n, min_cov, min_pctM))
cat(sprintf("Total positions: %d\n", length(unique(aaa$harmonized_pos))))

poskeep <- aaa %>%
    filter(n >= min_n, sc >= min_cov) %$% harmonized_pos %>%
    unique()

posall <- aaa %$% harmonized_pos %>% unique()

posmeth <- aaa %>%
    filter(harmonized_pos %in% poskeep) %>%
    filter(mm >= 10) %$% harmonized_pos %>%
    unique()

poskeep_l1hs <- aaa %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(n >= min_n, sc >= min_cov) %$% harmonized_pos %>%
    unique()
posmeth_l1hs <- aaa %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(harmonized_pos %in% poskeep_l1hs) %>%
    filter(mm >= 10) %$% harmonized_pos %>%
    unique()

cat(sprintf(
    "poskeep: %d, posmeth: %d, poskeep_l1hs: %d, posmeth_l1hs: %d\n",
    length(poskeep), length(posmeth), length(poskeep_l1hs), length(posmeth_l1hs)
))
p <- aaa %>%
    ggplot(aes(x = harmonized_pos, y = sc)) +
    geom_point() +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    mtclosed
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_cov1.pdf"), w = 12, h = 8)

p <- aaa %>%
    ggplot(aes(x = harmonized_pos, y = n)) +
    geom_point() +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    mtclosed
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_n1.pdf"), w = 12, h = 8)

p <- aaa %>%
    filter(harmonized_pos %in% poskeep) %>%
    ggplot(aes(x = harmonized_pos, y = sc)) +
    geom_point() +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    mtclosed
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posfilt_cov1.pdf"), w = 12, h = 8)

p <- aaa %>%
    filter(harmonized_pos %in% poskeep) %>%
    ggplot(aes(x = harmonized_pos, y = n)) +
    geom_point() +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    mtclosed
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posfilt_n1.pdf"), w = 12, h = 8)

library(tidyHeatmap)
subfam_order <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6")
# Heatmap of mean modification by subfamily x position (poskeep)
p <- aaa %>%
    filter(harmonized_pos %in% poskeep) %>%
    filter(n >= 30, sc >= 50) %>%
    mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order)) %>%
    group_by(rte_subfamily) %>%
    heatmap(rte_subfamily, harmonized_pos, mm,
        .scale = "none",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "Harmonized position",
        palette_value = circlize::colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
    )
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posfilt_mm_heatmap12.pdf"), w = 14, h = 5)

# Same heatmap filtered to posmeth
p <- aaa %>%
    filter(harmonized_pos %in% posmeth) %>%
    filter(n >= 50, sc >= 100) %>%
    mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order)) %>%
    group_by(rte_subfamily) %>%
    heatmap(rte_subfamily, harmonized_pos, mm,
        .scale = "none",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "Harmonized position (modified)"
    )
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posmeth_mm_heatmap1111.pdf"), w = 14, h = 5)

# Per-element heatmap: N_ORF1 condition, elements as rows grouped by subfamily
element_hm_data <- yl1annm %>%
    filter(condition == "N_ORF1") %>%
    filter(harmonized_pos %in% posmeth) %>%
    filter(cov >= 5) %>%
    group_by(gene_id, harmonized_pos, rte_subfamily) %>%
    summarise(mean_pctM = mean(pctM), .groups = "drop") %>%
    # Keep elements with data at enough positions
    group_by(gene_id) %>%
    filter(n() >= 0.5 * length(posmeth)) %>%
    ungroup()

# Remove positions where >70% of elements have no data
n_elements <- n_distinct(element_hm_data$gene_id)
pos_with_data <- element_hm_data %>%
    group_by(harmonized_pos) %>%
    summarise(n_with_data = n_distinct(gene_id), .groups = "drop") %>%
    filter(n_with_data >= 0.3 * n_elements)
element_hm_data <- element_hm_data %>%
    filter(harmonized_pos %in% pos_with_data$harmonized_pos) %>%
    mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order))

p <- element_hm_data %>%
    group_by(rte_subfamily) %>%
    heatmap(gene_id, harmonized_pos, mean_pctM,
        .scale = "none",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        column_title = "Harmonized position (modified)"
    )
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/element_posmeth_heatmap_orf11.pdf"), w = 14, h = 10)
####





subfam_count <- aaa %>%
    filter(harmonized_pos %in% posmeth, mm >= 25) %>%
    group_by(harmonized_pos) %>%
    summarise(n_subfams_above25 = n_distinct(rte_subfamily), .groups = "drop")

p <- aaa %>%
    filter(harmonized_pos %in% posmeth) %>%
    left_join(subfam_count, by = "harmonized_pos") %>%
    mutate(n_subfams_above25 = replace_na(n_subfams_above25, 0)) %>%
    ggplot(aes(x = harmonized_pos, y = mm, color = factor(n_subfams_above25))) +
    geom_point() +
    geom_hline(yintercept = 25) +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    labs(color = "# subfamilies\n>25% modified") +
    mtclosedgrid
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posmeth_mm1111.pdf"), w = 12, h = 8)

p <- aaa %>%
    filter(harmonized_pos %in% posmeth) %>%
    left_join(subfam_count, by = "harmonized_pos") %>%
    mutate(n_subfams_above25 = replace_na(n_subfams_above25, 0)) %>%
    ggplot(aes(x = harmonized_pos, y = mm, color = factor(n_subfams_above25), size = log2(sc))) +
    geom_point() +
    geom_hline(yintercept = 25) +
    facet_wrap(~rte_subfamily, scales = "free_y") +
    labs(color = "# subfamilies\n>25% modified") +
    mtclosedgrid
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/position_posmeth_mm11111.pdf"), w = 12, h = 8)


# //ANCHOR - HEATMAPs and PCAs
filters <- list("meth" = posmeth_l1hs, "all" = posall)
for (pos_filter_name in names(filters)) {
    pos_filter <- filters[[pos_filter_name]]
    tryCatch(
        {
            heatmap_l1hs <- yl1annm %>%
                filter(rte_subfamily == "L1HS") %>%
                filter(harmonized_pos %in% pos_filter) %>%
                group_by(gene_id, harmonized_pos, condition) %>%
                summarise(mean_pctM = mean(pctM, na.rm = TRUE), cov = sum(cov), .groups = "drop") %>%
                filter(cov >= 10) %>%
                mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))))

            # Filter to elements with >= 50% position coverage
            n_pos_hm <- n_distinct(heatmap_l1hs$harmonized_pos)
            good_elements <- heatmap_l1hs %>%
                group_by(gene_id, condition) %>%
                summarise(n_pos = n(), .groups = "drop") %>%
                filter(n_pos >= n_pos_hm * 0.5) %>%
                pull(gene_id) %>%
                unique()

            # Average across samples (collapse condition), pivot manually to ensure numeric matrix
            heatmap_mat <- heatmap_l1hs %>%
                filter(gene_id %in% good_elements) %>%
                pivot_wider(id_cols = gene_id, names_from = harmonized_pos, values_from = mean_pctM, values_fn = mean) %>%
                column_to_rownames("gene_id") %>%
                as.matrix()
            # Order columns numerically
            heatmap_mat <- heatmap_mat[, order(as.numeric(colnames(heatmap_mat)))]
            # Remove rows that are all NA
            heatmap_mat <- heatmap_mat[rowSums(!is.na(heatmap_mat)) > 0, ]
            # Remove columns where fewer than 20% of elements have data
            col_coverage <- colMeans(!is.na(heatmap_mat))
            heatmap_mat <- heatmap_mat[, col_coverage >= 0.2]

            p <- Heatmap(heatmap_mat,
                name = "% modified",
                col = circlize::colorRamp2(c(0, 50, 100), c("blue", "white", "red")),
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                na_col = "grey90",
                use_raster = TRUE,
                column_title = str_glue("L1HS {current_mod_name} modification at top positions")
            )
            mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/element_heatmap11111.pdf"), w = 10, h = 8)

            # PCA: elements as points (rows of heatmap_mat, positions as features)
            pca_mat_elements <- heatmap_mat
            pca_mat_elements[is.na(pca_mat_elements)] <- 0
            pca_mat_elements <- pca_mat_elements[, apply(pca_mat_elements, 2, var) > 0, drop = FALSE]
            pca_elements <- prcomp(pca_mat_elements, center = TRUE, scale. = FALSE)
            var_el <- summary(pca_elements)$importance[2, 1:2] * 100
            pca_el_df <- as.data.frame(pca_elements$x[, 1:2]) %>%
                rownames_to_column("gene_id") %>%
                left_join(rmann %>% dplyr::select(gene_id, rte_subfamily) %>% distinct(), by = "gene_id")

            p <- pca_el_df %>%
                left_join(rmann) %>%
                ggplot(aes(x = PC1, y = PC2, color = intactness_req)) +
                geom_point(alpha = 0.5, size = 1) +
                labs(
                    x = sprintf("PC1 (%.1f%%)", var_el[1]),
                    y = sprintf("PC2 (%.1f%%)", var_el[2]),
                    title = "PCA of L1HS elements by modification profile"
                ) +
                mtopen
            mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/pca_elements11.pdf"), w = 7, h = 5)

            p <- pca_el_df %>%
                left_join(rmann) %>%
                ggplot(aes(x = PC1, y = PC2, color = loc_superlowres_integrative_stranded)) +
                geom_point(alpha = 0.5, size = 1) +
                labs(
                    x = sprintf("PC1 (%.1f%%)", var_el[1]),
                    y = sprintf("PC2 (%.1f%%)", var_el[2]),
                    title = "PCA of L1HS elements by modification profile"
                ) +
                mtopen
            mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/pca_elements_loc.pdf"), w = 7, h = 5)

            # PCA: positions as points (transpose, elements as features)
            pca_mat_pos <- t(heatmap_mat)
            pca_mat_pos[is.na(pca_mat_pos)] <- 0
            pca_mat_pos <- pca_mat_pos[, apply(pca_mat_pos, 2, var) > 0, drop = FALSE]
            pca_pos <- prcomp(pca_mat_pos, center = TRUE, scale. = FALSE)
            var_pos <- summary(pca_pos)$importance[2, 1:2] * 100
            pca_pos_df <- as.data.frame(pca_pos$x[, 1:2]) %>%
                rownames_to_column("harmonized_pos") %>%
                mutate(harmonized_pos_num = as.numeric(harmonized_pos))
            library(ggrepel)
            p <- pca_pos_df %>%
                ggplot(aes(x = PC1, y = PC2, label = harmonized_pos)) +
                geom_point(alpha = 0.7, size = 2) +
                geom_text_repel(size = 2.5, max.overlaps = 15) +
                labs(
                    x = sprintf("PC1 (%.1f%%)", var_pos[1]),
                    y = sprintf("PC2 (%.1f%%)", var_pos[2]),
                    title = "PCA of positions by element modification profile"
                ) +
                mtopen
            mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/pca_positions.pdf"), w = 7, h = 5)

            # Element PCA loadings (features = positions)
            el_loadings <- as.data.frame(pca_elements$rotation[, 1:2]) %>%
                rownames_to_column("position") %>%
                mutate(position = factor(position, levels = position[order(as.numeric(position))]))

            p <- el_loadings %>%
                pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
                ggplot(aes(x = reorder(position, as.numeric(as.character(position))), y = loading)) +
                geom_col() +
                facet_wrap(~PC, ncol = 1, scales = "free_y") +
                labs(x = "Harmonized position", y = "Loading", title = "Element PCA: position contributions") +
                mtopen +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
            mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/pca_elements_loadings1.pdf"), w = 10, h = 6)

            # Position PCA loadings (features = elements, show top 20 by abs loading)
            pos_loadings <- as.data.frame(pca_pos$rotation[, 1:2]) %>%
                rownames_to_column("gene_id")

            for (pc in c("PC1", "PC2")) {
                top_genes <- pos_loadings %>%
                    arrange(desc(abs(!!sym(pc)))) %>%
                    head(20)
                p <- top_genes %>%
                    ggplot(aes(x = reorder(gene_id, !!sym(pc)), y = !!sym(pc))) +
                    geom_col() +
                    coord_flip() +
                    labs(x = NULL, y = "Loading", title = sprintf("Position PCA: top 20 element contributions to %s", pc)) +
                    mtopen
                mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/{pos_filter_name}/pca_positions_loadings_{pc}.pdf"), w = 8, h = 6)
            }
        },
        error = function(e) {
            print(str_glue("issue with heatmap or PCA"))
        }
    )
}

# //ANCHOR - DRACH EVO
if (current_mod_name == "m6A") {
    library(ggrepel)
    l1hs_cs
    aligned


    # detect DRACH motifs in consensus
    # DRACH = [D=A/G/T][R=A/G][A][C][H=A/C/T], m6A target is the A at position 3
    drach_matches <- vmatchPattern("DRACH", l1hs_cs, fixed = FALSE)
    # Get the A position (3rd base in each DRACH match)
    drach_a_positions <- start(drach_matches[[1]]) + 2
    cat(sprintf("Found %d DRACH motifs in L1HS consensus\n", length(drach_a_positions)))

    yl1annm %>%
        filter(condition == "N_ORF1") %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(consensus_pos == 2036) %>%
        pw()

    library(ggrepel)

    p <- yl1annm %>%
        filter(condition == "N_ORF1") %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(!is.na(seqnames)) %>%
        group_by(harmonized_pos, condition, rte_subfamily) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            label = ifelse(mm > 50, as.character(harmonized_pos), NA_character_),
            is_drach = harmonized_pos %in% drach_a_positions
        ) %>%
        ggplot(aes(x = harmonized_pos, y = mm)) +
        geom_point(data = . %>% filter(is_drach), color = "black", size = 2.5) +
        geom_point(aes(color = condition)) +
        facet_wrap(~rte_subfamily) +
        geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_along_element1111111.pdf"), w = 12, h = 8)


    p <- yl1annm %>%
        filter(condition == "N_ORF1") %>%
        filter(rte_subfamily == "L1HS") %>%
        filter(!is.na(seqnames)) %>%
        group_by(harmonized_pos, condition, rte_subfamily) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            label = ifelse(mm > 50, as.character(harmonized_pos), NA_character_),
            is_drach = harmonized_pos %in% drach_a_positions
        ) %>%
        ggplot(aes(x = harmonized_pos, y = mm)) +
        geom_point(data = . %>% filter(is_drach), color = "black", size = 2.5) +
        geom_point(aes(color = condition)) +
        facet_wrap(~rte_subfamily) +
        geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/mod_along_element1111111.pdf"), w = 12, h = 8)

    p <- yl1annm %>%
        filter(!is.na(seqnames)) %>%
        mutate(motif = ifelse(consensus_pos %in% drach_a_positions, "DRACH", "non-DRACH")) %>%
        group_by(motif, condition) %>%
        summarise(
            mean_meth = mean(pctM, na.rm = TRUE),
            se_meth = sd(pctM, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        ) %>%
        ggplot(aes(x = motif, y = mean_meth, fill = motif)) +
        geom_col() +
        geom_errorbar(aes(ymin = mean_meth - se_meth, ymax = mean_meth + se_meth), width = 0.2) +
        facet_wrap(~condition) +
        labs(x = NULL, y = "Mean % modified") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_vs_nondrach_barplot.pdf"), w = 6, h = 4)


    p <- yl1annm %>%
        filter(!is.na(seqnames)) %>%
        filter(consensus_pos %in% drach_a_positions) %>%
        group_by(consensus_pos, condition) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
        mutate(label = ifelse(mm > 50, as.character(consensus_pos), NA_character_)) %>%
        ggplot(aes(x = consensus_pos, y = mm)) +
        geom_point(aes(color = condition)) +
        facet_wrap(~condition) +
        geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_along_element_drach_only1111.pdf"), w = 8, h = 4)


    top_pos <- yl1annm %>%
        filter(!is.na(seqnames)) %>%
        filter(consensus_pos %in% drach_a_positions) %>%
        group_by(consensus_pos, condition) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
        filter(mm > 30) %$% consensus_pos %>%
        unique()


    # Build mapping from L1HS ungapped position to alignment column
    l1hs_aligned <- as.character(aligned[["L1HS"]])
    l1hs_chars <- strsplit(l1hs_aligned, "")[[1]]
    ungapped_pos <- 0
    l1hs_pos_to_aln_col <- integer(0)
    for (col in seq_along(l1hs_chars)) {
        if (l1hs_chars[col] != "-") {
            ungapped_pos <- ungapped_pos + 1
            l1hs_pos_to_aln_col[ungapped_pos] <- col
        }
    }

    # For each top_pos, extract flanking window from all subfamilies
    flank <- 2
    subfams <- names(aligned)
    window_data <- list()

    for (pos in top_pos) {
        # +1 to center on the A (pos points to base before the A)
        aln_col <- l1hs_pos_to_aln_col[pos] + 1
        col_start <- max(1, aln_col - flank)
        col_end <- min(width(aligned)[1], aln_col + flank)

        for (sf in subfams) {
            aln_seq <- as.character(subseq(aligned[[sf]], start = col_start, end = col_end))
            bases <- strsplit(aln_seq, "")[[1]]
            # Shift so the target position is at 0 instead of 1
            relative_pos <- seq(col_start - aln_col, col_end - aln_col)
            window_data <- c(window_data, list(tibble(
                l1hs_pos = pos,
                subfamily = sf,
                rel_pos = relative_pos,
                base = bases
            )))
        }
    }
    window_df <- bind_rows(window_data)

    # Check whether each subfamily's 5-base window fits the DRACH motif
    # DRACH: D=[A/G/T], R=[A/G], A, C, H=[A/C/T]
    drach_check <- window_df %>%
        arrange(l1hs_pos, subfamily, rel_pos) %>%
        group_by(l1hs_pos, subfamily) %>%
        summarise(motif = paste0(base, collapse = ""), .groups = "drop") %>%
        mutate(is_drach = grepl("^[AGT][AG]AC[ACT]$", motif))

    # Plot: tile heatmap colored by nucleotide, faceted by position
    base_colors <- c("A" = "#4DAF4A", "T" = "#E41A1C", "C" = "#377EB8", "G" = "#FF7F00", "-" = "grey90")

    plot_df <- window_df %>%
        mutate(
            subfamily = factor(subfamily, levels = rev(subfams)),
            l1hs_pos = factor(l1hs_pos, levels = sort(unique(l1hs_pos)))
        ) %>%
        left_join(drach_check %>% mutate(
            subfamily = factor(subfamily, levels = rev(subfams)),
            l1hs_pos = factor(l1hs_pos, levels = sort(unique(l1hs_pos)))
        ), by = c("l1hs_pos", "subfamily"))

    # Star marker for non-DRACH rows
    star_df <- plot_df %>%
        filter(!is_drach) %>%
        group_by(l1hs_pos, subfamily) %>%
        dplyr::slice(1)

    p <- plot_df %>%
        ggplot(aes(x = rel_pos, y = subfamily)) +
        geom_tile(aes(fill = base), color = "white", linewidth = 0.3) +
        geom_text(aes(label = base), size = 3) +
        geom_text(
            data = star_df, aes(x = flank + 0.8, y = subfamily, label = "*"),
            size = 9, color = "black", inherit.aes = FALSE
        ) +
        scale_fill_manual(values = base_colors) +
        facet_wrap(~l1hs_pos, ncol = 4, labeller = labeller(l1hs_pos = function(x) paste0("L1HS pos ", x))) +
        scale_x_continuous(breaks = seq(-flank, flank)) +
        coord_cartesian(xlim = c(-flank - 0.5, flank + 1.2), clip = "off") +
        labs(x = "Position relative to site", y = NULL) +
        mtopen +
        theme(strip.text = element_text(size = 8))
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_site_subfamily_alignment212.pdf"),
        w = 4 * min(length(top_pos), 4) + 1, h = 2 + length(subfams) * 0.4 * ceiling(length(top_pos) / 4)
    )



    # Detect DRACH motifs in all subfamily consensus sequences, normalize to element length
    all_cs_list <- list(L1HS = l1hs_cs, L1PA2 = l1pa2_cs, L1PA3 = l1pa3_cs, L1PA4 = l1pa4_cs, L1PA5 = l1pa5_cs)

    drach_all <- list()
    for (sf_name in names(all_cs_list)) {
        cs <- all_cs_list[[sf_name]]
        dm <- vmatchPattern("DRACH", cs, fixed = FALSE)
        a_pos <- start(dm[[1]]) + 2
        cs_len <- width(cs)[1]
        drach_all <- c(drach_all, list(tibble(
            subfamily = sf_name,
            position = a_pos,
            norm_position = a_pos / cs_len,
            consensus_length = cs_len
        )))
        cat(sprintf("Found %d DRACH motifs in %s consensus (length %d)\n", length(a_pos), sf_name, cs_len))
    }
    drach_all_df <- bind_rows(drach_all) %>%
        mutate(subfamily = factor(subfamily, levels = names(all_cs_list)))

    # Compare DRACH density along normalized element length
    p <- drach_all_df %>%
        ggplot(aes(x = norm_position, fill = subfamily)) +
        geom_histogram(bins = 50, alpha = 0.7) +
        facet_wrap(~subfamily, ncol = 1) +
        labs(x = "Normalized position along element", y = "DRACH motif count") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_density_by_subfamily.pdf"), w = 8, h = 8)

    # Overlay density comparison
    p <- drach_all_df %>%
        ggplot(aes(x = norm_position, color = subfamily)) +
        geom_density(linewidth = 0.8) +
        labs(x = "Normalized position along element", y = "DRACH motif density") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_density_overlay.pdf"), w = 8, h = 4)

    # Total DRACH count per subfamily
    p <- drach_all_df %>%
        group_by(subfamily) %>%
        summarise(n_drach = n(), consensus_length = first(consensus_length), .groups = "drop") %>%
        mutate(drach_per_kb = n_drach / (consensus_length / 1000)) %>%
        ggplot(aes(x = subfamily, y = drach_per_kb, fill = subfamily)) +
        geom_col() +
        labs(x = NULL, y = "DRACH motifs per kb") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_count_per_kb_by_subfamily.pdf"), w = 6, h = 4)

    # //ANCHOR - DRACH ALL SUBFAMILIES (harmonized positions)
    # Detect DRACH motifs in ALL subfamily consensus sequences and map to harmonized positions
    all_cs_list_full <- list(L1HS = l1hs_cs, L1PA2 = l1pa2_cs, L1PA3 = l1pa3_cs, L1PA4 = l1pa4_cs, L1PA5 = l1pa5_cs, L1PA6 = l1pa6_cs)

    drach_harmonized <- list()
    for (sf_name in names(all_cs_list_full)) {
        cs <- all_cs_list_full[[sf_name]]
        dm <- vmatchPattern("DRACH", cs, fixed = FALSE)
        a_pos <- start(dm[[1]]) + 2 # A is 3rd base in DRACH
        # Map consensus_pos to harmonized_pos
        sf_mapping <- pos_to_aln_df %>% filter(subfamily_consensus == sf_name)
        drach_harmonized <- c(drach_harmonized, list(
            tibble(subfamily = sf_name, consensus_pos = a_pos) %>%
                left_join(sf_mapping, by = c("subfamily" = "subfamily_consensus", "consensus_pos")) %>%
                filter(!is.na(harmonized_pos))
        ))
        cat(sprintf("%s: %d DRACH A-positions mapped to harmonized coords\n", sf_name, length(a_pos)))
    }
    drach_harmonized_df <- bind_rows(drach_harmonized) %>%
        mutate(subfamily = factor(subfamily, levels = names(all_cs_list_full)))

    # All unique harmonized positions where ANY subfamily has a DRACH
    drach_harmonized_positions <- unique(drach_harmonized_df$harmonized_pos)
    cat(sprintf("Total unique DRACH harmonized positions (any subfamily): %d\n", length(drach_harmonized_positions)))

    # How many subfamilies share each DRACH position
    drach_sharing <- drach_harmonized_df %>%
        group_by(harmonized_pos) %>%
        summarise(n_subfams = n_distinct(subfamily), subfams = paste(sort(unique(subfamily)), collapse = ","), .groups = "drop")

    p <- drach_sharing %>%
        ggplot(aes(x = factor(n_subfams))) +
        geom_bar() +
        labs(x = "# subfamilies with DRACH at position", y = "# harmonized positions") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_sharing_barplot.pdf"), w = 6, h = 4)

    # Modification at DRACH positions across all subfamilies
    drach_meth_summary <- yl1annm %>%
        filter(!is.na(seqnames)) %>%
        filter(harmonized_pos %in% drach_harmonized_positions) %>%
        group_by(harmonized_pos, rte_subfamily, condition) %>%
        summarise(mm = mean(pctM, na.rm = TRUE), n = n(), .groups = "drop")

    # DRACH vs non-DRACH barplot using harmonized positions (all subfamilies)
    p <- yl1annm %>%
        filter(!is.na(seqnames)) %>%
        mutate(motif = ifelse(harmonized_pos %in% drach_harmonized_positions, "DRACH", "non-DRACH")) %>%
        group_by(motif, condition, rte_subfamily) %>%
        summarise(mean_meth = mean(pctM, na.rm = TRUE), se_meth = sd(pctM, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
        mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order)) %>%
        ggplot(aes(x = motif, y = mean_meth, fill = motif)) +
        geom_col(position = "dodge") +
        geom_errorbar(aes(ymin = mean_meth - se_meth, ymax = mean_meth + se_meth), width = 0.2) +
        facet_grid(rte_subfamily ~ condition) +
        labs(x = NULL, y = "Mean % modified") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_vs_nondrach_all_subfams.pdf"), w = 8, h = 10)

    # Modification along element at DRACH positions, all subfamilies
    p <- drach_meth_summary %>%
        left_join(drach_sharing, by = "harmonized_pos") %>%
        mutate(
            rte_subfamily = factor(rte_subfamily, levels = subfam_order),
            label = ifelse(mm > 30, as.character(harmonized_pos), NA_character_)
        ) %>%
        ggplot(aes(x = harmonized_pos, y = mm, color = condition)) +
        geom_point(alpha = 0.7) +
        geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, na.rm = TRUE) +
        facet_wrap(~rte_subfamily, ncol = 1) +
        labs(x = "Harmonized position", y = "Mean % modified (DRACH sites only)") +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_mod_along_element_all_subfams.pdf"), w = 12, h = 12)

    # Heatmap: mean modification at DRACH harmonized positions, subfamilies as rows
    drach_hm_data <- drach_meth_summary %>%
        filter(n >= 100) %>%
        group_by(harmonized_pos, rte_subfamily) %>%
        summarise(mm = mean(mm), .groups = "drop") %>%
        mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order))

    p <- drach_hm_data %>%
        group_by(rte_subfamily) %>%
        heatmap(rte_subfamily, harmonized_pos, mm,
            .scale = "none",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_title = "DRACH harmonized positions"
        )
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_heatmap_all_subfams111.pdf"), w = 14, h = 5)

    # Per-condition heatmaps
    for (cond in unique(drach_meth_summary$condition)) {
        drach_hm_cond <- drach_meth_summary %>%
            filter(condition == cond, n >= 100) %>%
            mutate(rte_subfamily = factor(rte_subfamily, levels = subfam_order))

        p <- drach_hm_cond %>%
            group_by(rte_subfamily) %>%
            heatmap(rte_subfamily, harmonized_pos, mm,
                .scale = "none",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                column_title = str_glue("DRACH positions — {cond}"),
            )
        mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_heatmap_{cond}.pdf"), w = 14, h = 5)
    }

    # Alignment window visualization for top DRACH positions in at least 1 subfamily
    top_drach_harmonized <- drach_meth_summary %>%
        filter(n > 100) %>%
        group_by(harmonized_pos, rte_subfamily) %>%
        summarise(mm = mean(mm, na.rm = TRUE), .groups = "drop") %>%
        filter(mm > 30) %$% harmonized_pos %>%
        unique()

    flank <- 2
    subfams <- names(aligned)
    window_data_all <- list()

    for (hpos in top_drach_harmonized) {
        col_start <- max(1, hpos - flank)
        col_end <- min(width(aligned)[1], hpos + flank)

        for (sf in subfams) {
            aln_seq <- as.character(subseq(aligned[[sf]], start = col_start, end = col_end))
            bases <- strsplit(aln_seq, "")[[1]]
            relative_pos <- seq(col_start - hpos, col_end - hpos)
            window_data_all <- c(window_data_all, list(tibble(
                harmonized_pos = hpos,
                subfamily = sf,
                rel_pos = relative_pos,
                base = bases
            )))
        }
    }
    window_df_all <- bind_rows(window_data_all)

    drach_check_all <- window_df_all %>%
        arrange(harmonized_pos, subfamily, rel_pos) %>%
        group_by(harmonized_pos, subfamily) %>%
        summarise(motif = paste0(base, collapse = ""), .groups = "drop") %>%
        mutate(is_drach = grepl("^[AGT][AG]AC[ACT]$", motif))

    base_colors <- c("A" = "#4DAF4A", "T" = "#E41A1C", "C" = "#377EB8", "G" = "#FF7F00", "-" = "grey90")

    plot_df_all <- window_df_all %>%
        mutate(
            subfamily = factor(subfamily, levels = rev(subfams)),
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))
        ) %>%
        left_join(drach_check_all %>% mutate(
            subfamily = factor(subfamily, levels = rev(subfams)),
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos)))
        ), by = c("harmonized_pos", "subfamily"))

    star_df_all <- plot_df_all %>%
        filter(!is_drach) %>%
        group_by(harmonized_pos, subfamily) %>%
        dplyr::slice(1)

    # Mean modification per subfamily at each position (averaged across conditions)
    meth_annot <- drach_meth_summary %>%
        group_by(harmonized_pos, rte_subfamily) %>%
        summarise(mm = mean(mm, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            subfamily = factor(rte_subfamily, levels = rev(subfams)),
            harmonized_pos = factor(harmonized_pos, levels = sort(unique(as.numeric(as.character(plot_df_all$harmonized_pos))))),
            mm_label = sprintf("%.0f%%", mm)
        ) %>%
        filter(!is.na(harmonized_pos))

    p <- plot_df_all %>%
        ggplot(aes(x = rel_pos, y = subfamily)) +
        geom_tile(aes(fill = base), color = "white", linewidth = 0.3) +
        geom_text(aes(label = base), size = 3) +
        geom_text(
            data = star_df_all, aes(x = flank + 0.8, y = subfamily, label = "*"),
            size = 9, color = "black", inherit.aes = FALSE
        ) +
        geom_text(
            data = meth_annot, aes(x = flank + 1.8, y = subfamily, label = mm_label, color = mm),
            size = 5, fontface = "bold", inherit.aes = FALSE
        ) +
        scale_color_gradient2(low = "blue", mid = "grey40", high = "red", midpoint = 25, guide = "none") +
        scale_fill_manual(values = base_colors) +
        facet_wrap(~harmonized_pos, ncol = 4, labeller = labeller(harmonized_pos = function(x) paste0("pos ", x))) +
        scale_x_continuous(breaks = seq(-flank, flank)) +
        coord_cartesian(xlim = c(-flank - 0.5, flank + 2.5), clip = "off") +
        labs(x = "Position relative to site", y = NULL) +
        mtopen +
        theme(strip.text = element_text(size = 8))
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_alignment_all_subfams1.pdf"),
        w = 4 * min(length(top_drach_harmonized), 4) + 2,
        h = 2 + length(subfams) * 0.4 * ceiling(length(top_drach_harmonized) / 4)
    )

    # Simplified: only positions where DRACH status CHANGES across subfamilies
    drach_variable <- drach_check_all %>%
        group_by(harmonized_pos) %>%
        summarise(n_drach = sum(is_drach), n_nondrach = sum(!is_drach), .groups = "drop") %>%
        filter(n_drach > 0, n_nondrach > 0) %$% harmonized_pos
    cat(sprintf("Positions with variable DRACH status across subfamilies: %d\n", length(drach_variable)))

    plot_df_var <- plot_df_all %>% filter(as.numeric(as.character(harmonized_pos)) %in% drach_variable)
    star_df_var <- star_df_all %>% filter(as.numeric(as.character(harmonized_pos)) %in% drach_variable)
    meth_annot_var <- meth_annot %>% filter(as.numeric(as.character(harmonized_pos)) %in% drach_variable)

    p <- plot_df_var %>%
        ggplot(aes(x = rel_pos, y = subfamily)) +
        geom_tile(aes(fill = base), color = "white", linewidth = 0.3) +
        geom_text(aes(label = base), size = 3) +
        geom_text(
            data = star_df_var, aes(x = flank + 0.8, y = subfamily, label = "*"),
            size = 9, color = "black", inherit.aes = FALSE
        ) +
        geom_text(
            data = meth_annot_var, aes(x = flank + 1.8, y = subfamily, label = mm_label, color = mm),
            size = 5, fontface = "bold", inherit.aes = FALSE
        ) +
        scale_color_gradient2(low = "blue", mid = "grey40", high = "red", midpoint = 25, guide = "none") +
        scale_fill_manual(values = base_colors) +
        facet_wrap(~harmonized_pos, ncol = 4, labeller = labeller(harmonized_pos = function(x) paste0("pos ", x))) +
        scale_x_continuous(breaks = seq(-flank, flank)) +
        coord_cartesian(xlim = c(-flank - 0.5, flank + 2.5), clip = "off") +
        labs(x = "Position relative to site", y = NULL, title = "Positions with variable DRACH status") +
        mtopen +
        theme(strip.text = element_text(size = 8))
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1s/drach_alignment_variable_only1.pdf"),
        w = 4 * min(length(drach_variable), 4) + 2,
        h = 2 + length(subfams) * 0.4 * ceiling(length(drach_variable) / 4)
    )

    # //ANCHOR - DRACH POLYMORPHISM WITHIN L1HS ELEMENTS
    # Scan each FL L1HS element's actual genomic sequence for DRACH motifs,
    # then map back to harmonized positions to find polymorphic sites

    # Get FL L1HS element info
    l1hs_fl_ids <- rmann %>%
        filter(rte_subfamily == "L1HS", rte_length_req == "FL") %>%
        pull(gene_id)
    l1hs_elements <- rmann %>%
        filter(gene_id %in% l1hs_fl_ids) %>%
        dplyr::select(gene_id, seqnames, start, end, strand) %>%
        distinct()
    cat(sprintf("FL L1HS elements: %d\n", nrow(l1hs_elements)))

    # Extract full element sequences in one batch
    el_gr <- GRanges(
        seqnames = l1hs_elements$seqnames,
        ranges = IRanges(start = l1hs_elements$start, end = l1hs_elements$end)
    )
    el_seqs <- getSeq(Rsamtools::FaFile(conf$reference), el_gr)
    # Reverse complement minus-strand elements so sequence is in consensus orientation
    minus_idx <- as.character(l1hs_elements$strand) == "-"
    el_seqs[minus_idx] <- reverseComplement(el_seqs[minus_idx])
    names(el_seqs) <- l1hs_elements$gene_id

    # Scan each element for DRACH motifs
    l1hs_mapping <- pos_to_aln_df %>% filter(subfamily_consensus == "L1HS")

    drach_per_element <- list()
    for (i in seq_along(el_seqs)) {
        gid <- names(el_seqs)[i]
        seq_i <- el_seqs[[i]]
        dm <- matchPattern("DRACH", seq_i, fixed = FALSE)
        if (length(dm) == 0) next
        a_positions <- start(dm) + 2L # A is 3rd base in DRACH
        motifs <- as.character(dm)

        drach_per_element <- c(drach_per_element, list(tibble(
            gene_id = gid,
            element_pos = a_positions,
            motif = motifs,
            is_drach = TRUE
        )))
    }
    drach_element_df <- bind_rows(drach_per_element)
    cat(sprintf(
        "Total DRACH hits across %d elements: %d\n",
        n_distinct(drach_element_df$gene_id), nrow(drach_element_df)
    ))

    # Map element positions to consensus positions using consensus_index_long
    # consensus_index_long has: gene_id, sequence_pos, consensus_pos
    drach_element_mapped <- drach_element_df %>%
        left_join(consensus_index_long %>% dplyr::select(gene_id, sequence_pos, consensus_pos),
            by = c("gene_id", "element_pos" = "sequence_pos")
        ) %>%
        filter(!is.na(consensus_pos)) %>%
        left_join(l1hs_mapping %>% dplyr::select(consensus_pos, harmonized_pos),
            by = "consensus_pos"
        ) %>%
        filter(!is.na(harmonized_pos))

    # For each element × harmonized_pos, mark whether it has DRACH
    # First, get ALL harmonized positions that are A in any subfamily (these are the positions we track)
    all_element_positions <- modbase_positions_df %>%
        filter(gene_id %in% l1hs_fl_ids) %>%
        dplyr::select(gene_id, harmonized_pos) %>%
        distinct()

    # Join: for each element × position, TRUE if DRACH found there, FALSE otherwise
    drach_poly_df <- all_element_positions %>%
        left_join(
            drach_element_mapped %>%
                dplyr::select(gene_id, harmonized_pos, motif) %>%
                mutate(is_drach = TRUE) %>%
                distinct(gene_id, harmonized_pos, .keep_all = TRUE),
            by = c("gene_id", "harmonized_pos")
        ) %>%
        mutate(
            is_drach = replace_na(is_drach, FALSE),
            motif = replace_na(motif, "")
        )

    # Identify polymorphic DRACH sites: positions where >=20% of elements differ from majority
    drach_poly_summary <- drach_poly_df %>%
        group_by(harmonized_pos) %>%
        summarise(
            n_elements = n(),
            n_drach = sum(is_drach),
            n_nondrach = sum(!is_drach),
            pct_drach = 100 * n_drach / n_elements,
            .groups = "drop"
        ) %>%
        mutate(is_polymorphic = pmin(pct_drach, 100 - pct_drach) >= 10)

    polymorphic_positions <- drach_poly_summary %>%
        filter(is_polymorphic) %>%
        pull(harmonized_pos)
    cat(sprintf(
        "Polymorphic DRACH sites (>=20%% minor allele): %d / %d\n",
        length(polymorphic_positions), nrow(drach_poly_summary)
    ))

    # Summary: how many polymorphic sites
    p <- drach_poly_summary %>%
        ggplot(aes(x = pct_drach)) +
        geom_histogram(bins = 20, fill = "steelblue", color = "white") +
        geom_vline(xintercept = c(20, 80), linetype = "dashed", color = "red") +
        labs(
            x = "% elements with DRACH motif", y = "# positions",
            title = sprintf("L1HS DRACH polymorphism (%d positions)", nrow(drach_poly_summary))
        ) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_polymorphism_histogram.pdf"), w = 6, h = 4)

    # Where are polymorphic sites along the element
    p <- drach_poly_summary %>%
        mutate(status = case_when(
            is_polymorphic ~ "Polymorphic",
            pct_drach >= 80 ~ "Conserved DRACH",
            TRUE ~ "Conserved non-DRACH"
        )) %>%
        ggplot(aes(x = harmonized_pos, y = pct_drach, color = status)) +
        geom_point(size = 2) +
        geom_hline(yintercept = c(20, 80), linetype = "dashed", alpha = 0.5) +
        labs(x = "Harmonized position", y = "% elements with DRACH", color = NULL) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_polymorphism_along_element.pdf"), w = 12, h = 5)

    # For polymorphic positions: compare methylation between DRACH and non-DRACH elements
    if (length(polymorphic_positions) > 0) {
        # Join DRACH status with methylation data
        poly_meth <- yl1annm %>%
            filter(rte_subfamily == "L1HS") %>%
            filter(harmonized_pos %in% polymorphic_positions) %>%
            left_join(drach_poly_df %>% dplyr::select(gene_id, harmonized_pos, is_drach, motif),
                by = c("gene_id", "harmonized_pos")
            )

        # Boxplot: methylation by DRACH status at each polymorphic position
        p <- poly_meth %>%
            filter(!is.na(is_drach)) %>%
            mutate(drach_status = ifelse(is_drach, "DRACH", "non-DRACH")) %>%
            group_by(gene_id, harmonized_pos, condition, drach_status) %>%
            summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop") %>%
            ggplot(aes(x = drach_status, y = mean_pctM, fill = drach_status)) +
            geom_violin(alpha = 0.3, scale = "width") +
            geom_boxplot(width = 0.2, outlier.shape = NA) +
            facet_wrap(~harmonized_pos, scales = "free_y") +
            labs(x = NULL, y = "Mean % modified", title = "Methylation at polymorphic DRACH sites") +
            mtopen +
            theme(legend.position = "none")
        mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_poly_meth_boxplot.pdf"),
            w = 3 * min(length(polymorphic_positions), 4), h = 3 * ceiling(length(polymorphic_positions) / 4)
        )

        # Alignment window + methylation for polymorphic positions
        # Show the actual 5-base motif distribution per position
        poly_motif_summary <- drach_poly_df %>%
            filter(harmonized_pos %in% polymorphic_positions) %>%
            group_by(harmonized_pos, motif, is_drach) %>%
            summarise(n = n(), .groups = "drop") %>%
            group_by(harmonized_pos) %>%
            mutate(pct = 100 * n / sum(n)) %>%
            arrange(harmonized_pos, -n)

        # Get mean methylation per motif variant
        poly_motif_meth <- drach_poly_df %>%
            filter(harmonized_pos %in% polymorphic_positions) %>%
            left_join(
                yl1annm %>%
                    filter(rte_subfamily == "L1HS") %>%
                    group_by(gene_id, harmonized_pos) %>%
                    summarise(mean_pctM = mean(pctM, na.rm = TRUE), .groups = "drop"),
                by = c("gene_id", "harmonized_pos")
            ) %>%
            group_by(harmonized_pos, motif, is_drach) %>%
            summarise(mm = mean(mean_pctM, na.rm = TRUE), n = n(), .groups = "drop") %>%
            filter(!is.na(mm)) %>%
            group_by(harmonized_pos) %>%
            mutate(pct = 100 * n / sum(n)) %>%
            arrange(harmonized_pos, -n)

        # Tile plot: show top motif variants per position with their frequency and methylation
        p <- poly_motif_meth %>%
            group_by(harmonized_pos) %>%
            slice_head(n = 5) %>%
            ungroup() %>%
            mutate(
                drach_label = ifelse(is_drach, "DRACH", "*"),
                motif_label = sprintf("%s (%d%%, %.0f%%mod)", motif, round(pct), mm),
                harmonized_pos = factor(harmonized_pos)
            ) %>%
            ggplot(aes(x = 1, y = reorder(motif_label, mm), fill = mm)) +
            geom_tile(color = "white") +
            geom_text(aes(label = drach_label), size = 3) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 25, name = "% mod") +
            facet_wrap(~harmonized_pos, scales = "free_y", ncol = 3) +
            labs(x = NULL, y = NULL, title = "Motif variants at polymorphic DRACH sites") +
            mtopen +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/drach_poly_motif_variants.pdf"),
            w = 4 * min(length(polymorphic_positions), 3) + 1, h = 3 * ceiling(length(polymorphic_positions) / 3) + 1
        )
    }
}





# //ANCHOR - Element Expression load

# Summarize mean modification per element across positions
heatmap_l1hs <- yl1annm %>%
    filter(rte_subfamily == "L1HS") %>%
    filter(harmonized_pos %in% posall) %>%
    group_by(gene_id, harmonized_pos, condition) %>%
    summarise(mean_pctM = mean(pctM, na.rm = TRUE), cov = sum(cov), .groups = "drop") %>%
    filter(cov >= 10) %>%
    mutate(harmonized_pos = factor(harmonized_pos, levels = sort(unique(harmonized_pos))))

element_meth <- heatmap_l1hs %>%
    group_by(gene_id, condition) %>%
    summarise(mean_mod = mean(mean_pctM, na.rm = TRUE), n_pos = n(), sc = sum(cov), .groups = "drop") %>%
    filter(n_pos >= length(posall) / 1.5, sc >= 3 * length(posall))

elements_to_analyze <- unique(element_meth %$% gene_id)

contrasts_lr <- c("condition_N_ORF1_vs_N_TOT")

deseq_lr <- tibble(gene_id = elements_to_analyze)
for (contrast in contrasts_lr) {
    filepath <- str_glue("lrna/results/agg/deseq/relaxed/{contrast}/results_rtes.csv")
    df <- read_csv(filepath, show_col_types = FALSE)
    colnames(df)[1] <- "gene_id"
    df <- df %>%
        filter(gene_id %in% elements_to_analyze) %>%
        dplyr::select(gene_id, baseMean, log2FoldChange, lfcSE, padj) %>%
        dplyr::rename(
            !!str_glue("baseMean_{contrast}") := baseMean,
            !!str_glue("l2fc_{contrast}") := log2FoldChange,
            !!str_glue("lfcSE_{contrast}") := lfcSE,
            !!str_glue("padj_{contrast}") := padj
        )
    deseq_lr <- deseq_lr %>% left_join(df, by = "gene_id")
}
# Fill NAs: L2FC -> 0, padj -> 1
for (contrast in contrasts_lr) {
    deseq_lr <- deseq_lr %>%
        mutate(!!str_glue("l2fc_{contrast}") := replace_na(!!sym(str_glue("l2fc_{contrast}")), 0)) %>%
        mutate(!!str_glue("padj_{contrast}") := replace_na(!!sym(str_glue("padj_{contrast}")), 1))
}

normed_lr_raw <- read_csv("lrna/results/agg/deseq/relaxed/counttablesizenormed.csv", show_col_types = FALSE)
colnames(normed_lr_raw)[1] <- "gene_id"
normed_lr <- normed_lr_raw %>% filter(gene_id %in% elements_to_analyze)

# Pivot to long, join with sample metadata, compute per-condition means
normed_lr_long <- normed_lr %>%
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "normed_count") %>%
    left_join(sample_table_lrna %>% dplyr::select(sample_name, condition, bcondition, libtype), by = "sample_name")

normed_lr_mean <- normed_lr_long %>%
    group_by(gene_id, condition) %>%
    summarise(mean_normed = mean(normed_count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_normed, names_prefix = "mean_lr_")


enrichment_lr <- deseq_lr %>%
    left_join(normed_lr_mean, by = "gene_id") %>%
    mutate(across(starts_with("mean_lr_"), ~ replace_na(., 0)))



# Join with enrichment data
element_combined <- element_meth %>%
    left_join(enrichment_lr, by = "gene_id") %>%
    filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT))

# 1. Scatter: mean modification vs ORF1 enrichment (log2FC)
p <- element_combined %>%
    ggplot(aes(x = mean_mod, y = l2fc_condition_N_ORF1_vs_N_TOT)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~condition) +
    labs(x = "Mean %  modification", y = "ORF1 enrichment (log2FC)") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_vs_orf1_enrichment_scatter5.pdf"), w = 8, h = 4)

# 2. Scatter: mean modification vs expression (mean_lr_N_TOT and mean_lr_N_ORF1)
p <- element_combined %>%
    ggplot(aes(x = mean_mod, y = log2(mean_lr_N_TOT + 1))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~condition) +
    labs(x = "Mean %  modification", y = "log2(mean_lr_N_TOT + 1)") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_vs_expression_tot_scatter.pdf"), w = 8, h = 4)

p <- element_combined %>%
    ggplot(aes(x = mean_mod, y = log2(mean_lr_N_ORF1 + 1))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~condition) +
    labs(x = "Mean %  modification", y = "log2(mean_lr_N_ORF1 + 1)") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_vs_expression_orf1_scatter.pdf"), w = 8, h = 4)

p <- element_combined %>%
    ggplot(aes(x = mean_mod, y = mean_lr_N_ORF1)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~condition) +
    labs(x = "Mean %  modification", y = "mean_lr_N_ORF1)") +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_vs_expression_linear_orf1_scatter.pdf"), w = 8, h = 4)

# 3. Per-position correlation with enrichment
pos_enrichment_cor <- heatmap_l1hs %>%
    left_join(enrichment_lr, by = "gene_id") %>%
    filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT)) %>%
    group_by(harmonized_pos, condition) %>%
    summarise(
        cor_enrichment = cor(mean_pctM, l2fc_condition_N_ORF1_vs_N_TOT, use = "pairwise.complete.obs"),
        cor_expression = cor(mean_pctM, log2(mean_lr_N_TOT + 1), use = "pairwise.complete.obs"),
        n = n(),
        .groups = "drop"
    ) %>%
    filter(n >= 10)

p <- pos_enrichment_cor %>%
    pivot_longer(c(cor_enrichment, cor_expression), names_to = "metric", values_to = "correlation") %>%
    mutate(
        metric = recode(metric, cor_enrichment = "ORF1 enrichment", cor_expression = "Expression"),
        harmonized_pos = as.numeric(as.character(harmonized_pos)),
        label = ifelse(abs(correlation) > 0.2, as.character(harmonized_pos), NA_character_)
    ) %>%
    ggplot(aes(x = harmonized_pos, y = correlation, color = metric)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    facet_wrap(~condition) +
    labs(x = "Harmonized position", y = "Pearson correlation", color = NULL) +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_position_cor_with_enrichment141.pdf"), w = 10, h = 5)

p <- pos_enrichment_cor %>%
    pivot_longer(c(cor_enrichment, cor_expression), names_to = "metric", values_to = "correlation") %>%
    mutate(
        metric = recode(metric, cor_enrichment = "ORF1 enrichment", cor_expression = "Expression"),
        harmonized_pos = as.numeric(as.character(harmonized_pos)),
        label = ifelse(abs(correlation) > 0.2, as.character(harmonized_pos), NA_character_)
    ) %>%
    ggplot(aes(x = harmonized_pos, y = correlation, color = metric)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    facet_grid(rows = vars(condition), cols = vars(metric)) +
    labs(x = "Harmonized position", y = "Pearson correlation", color = NULL) +
    mtopen
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_position_cor_with_enrichmentfacet_14.pdf"), w = 10, h = 5)


# 4. Binned heatmap: split elements into modification quartiles, show enrichment
p <- element_combined %>%
    mutate(mod_bin = cut(mean_mod, breaks = quantile(mean_mod, probs = seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE)) %>%
    filter(!is.na(mod_bin)) %>%
    ggplot(aes(x = mod_bin, y = l2fc_condition_N_ORF1_vs_N_TOT)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~condition) +
    labs(x = "Modification quartile", y = "ORF1 enrichment (log2FC)") +
    mtopen +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/mod_quartile_vs_enrichment_boxplot.pdf"), w = 8, h = 5)

dir.create(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_pos"), recursive = TRUE, showWarnings = FALSE)


positions_of_interest <- pos_enrichment_cor %>%
    pivot_longer(c(cor_enrichment, cor_expression), names_to = "metric", values_to = "correlation") %>%
    mutate(
        metric = recode(metric, cor_enrichment = "ORF1 enrichment", cor_expression = "Expression"),
        harmonized_pos = as.numeric(as.character(harmonized_pos)),
        label = ifelse(abs(correlation) > 0.2, as.character(harmonized_pos), NA_character_)
    ) %>%
    group_by(metric) %>%
    arrange(abs(correlation)) %>%
    slice_tail(n = 10) %$% harmonized_pos %>%
    unique()
pos_enrichment_cor %>% filter(harmonized_pos %in% positions_of_interest)
# positions_of_interest <- c(5664, 602)

pos_df <- pos_df %>% left_join(rmann)
for (pos_i in positions_of_interest) {
    pos_df <- heatmap_l1hs %>%
        filter(harmonized_pos == pos_i) %>%
        left_join(enrichment_lr, by = "gene_id") %>%
        filter(!is.na(l2fc_condition_N_ORF1_vs_N_TOT))
    if (nrow(pos_df) < 10) next

    p <- pos_df %>%
        ggplot(aes(x = mean_pctM, y = l2fc_condition_N_ORF1_vs_N_TOT)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE) +
        facet_wrap(~condition) +
        labs(x = sprintf("%%  at position %s", pos_i), y = "ORF1 enrichment (log2FC)", title = sprintf("Position %s", pos_i)) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_pos/pos{pos_i}_vs_orf1_enrichment.pdf"), w = 8, h = 4)

    p <- pos_df %>%
        ggplot(aes(x = mean_pctM, y = log2(mean_lr_N_TOT + 1))) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE) +
        facet_wrap(~condition) +
        labs(x = sprintf("%%  at position %s", pos_i), y = "log2(mean_lr_N_TOT + 1)", title = sprintf("Position %s", pos_i)) +
        mtopen
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_pos/pos{pos_i}_vs_expression.pdf"), w = 8, h = 4)


    p <- pos_df %>%
        ggplot(aes(y = mean_pctM, x = req_integrative)) +
        geom_violin(aes(fill = req_integrative), alpha = 0.3, scale = "width") +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        geom_jitter(alpha = 0.3, size = 0.8, width = 0.1) +
        facet_wrap(~condition) +
        labs(y = sprintf("%%  at position %s", pos_i), x = "Intactness", title = sprintf("Position %s", pos_i)) +
        mtopen +
        theme(legend.position = "none")
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_pos/pos{pos_i}_vs_intactness.pdf"), w = 8, h = 4)

    p <- pos_df %>%
        filter(!is.na(orf1_intactness)) %>%
        ggplot(aes(y = mean_pctM, x = orf1_intactness)) +
        geom_violin(aes(fill = orf1_intactness), alpha = 0.3, scale = "width") +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        geom_jitter(alpha = 0.3, size = 0.8, width = 0.1) +
        facet_wrap(~condition) +
        labs(y = sprintf("%%  at position %s", pos_i), x = "ORF1 intactness", title = sprintf("Position %s", pos_i)) +
        mtopen +
        theme(legend.position = "none")
    mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_pos/pos{pos_i}_vs_orf1_intactness11.pdf"), w = 8, h = 4)
}


#### Overall read modification plots


flreads1 <- reads %>%
    group_by(read_id) %>%
    mutate(read_element_span = abs(min(start) - max(end)))
flreads <- flreads1 %>% mutate(read_fl = read_element_span >= 5000)
flreads %$% rte_subfamily %>% table()


readsdf1 <- flreads %>%
    left_join(rmann %>%
        dplyr::select(gene_id, start, end, strand, rte_length_req, intactness_req, rte_subfamily) %>%
        dplyr::rename(element_strand = strand, element_start = start, element_end = end), by = c("gene_id")) %>%
    filter(strand == element_strand) %>%
    dplyr::mutate(mod_indicator = ifelse(mod_qual > 0.5, 1, 0)) %>%
    mutate(genome_pos = start) %>%
    left_join(modbase_positions_df %>% dplyr::select(gene_id, genome_pos, harmonized_pos))


by_read <- readsdf1 %>%
    group_by(read_id, gene_id, condition, sample) %>%
    summarise(mm = mean(mod_indicator))


by_pos <- readsdf1 %>%
    group_by(harmonized_pos, gene_id, condition, sample) %>%
    summarise(mm = mean(mod_indicator))

p <- by_read %>% ggplot(aes(x = mm)) +
    geom_histogram() +
    mtclosed
mysaveandstore(str_glue("lrna/results/plots/modification_analysis/{current_mod_name}/rtes/l1hs/per_read/density.pdf"), w = 4, h = 4)




by_read %>% filter(mm >= 0.05)
