module_name <- "lrna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")

library(rtracklayer)
library(Biostrings)
library(cowplot)
# library(zoo)
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
# library(msigdbr)
library(Biostrings)

samples <- conf$samples
sample_table <- read_csv(sprintf("conf/sample_table_%s.csv", conf$prefix))
sample_table <- sample_table[match(samples, sample_table$sample_name), ]

{
    genome_lengths <- fasta.seqlengths(conf$reference)
    chromosomesAll <- names(genome_lengths)
    nonrefchromosomes <- grep("nonref", chromosomesAll, value = TRUE)
    refchromosomes <- grep("^chr", chromosomesAll, value = TRUE)
    autosomes <- grep("^chr[1-9]", refchromosomes, value = TRUE)
    chrX <- c("chrX")
    chrY <- c("chrY")

    MINIMUMCOVERAGE <- conf$MINIMUM_COVERAGE_FOR_METHYLATION_ANALYSIS
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
            bedmethlpaths = sprintf("lrna/intermediates/%s/methylation/%s_m6a_bedMethyl.bed", samples, samples),
            dmrs = "lrna/results/tables/dmrs.CG_m.tsv",
            dmls = "lrna/results/tables/dmls.CG_m.tsv",
            read_mods = sprintf("lrna/intermediates/%s/methylation/%s_readmods_%s_%s.tsv", samples, samples, "NoContext", conf$rte_subfamily_read_level_analysis)
        ), env = globalenv())
        assign("outputs", list(outfile = "lrna/outfiles/bedmethylanalysis.txt"), env = globalenv())
    }
)

dir.create("lrna/Rintermediates", showWarnings = FALSE)


ref_annotation_dir <- conf$reference_annotation_dir
rte_subfamily_read_level_analysis <- conf$rte_subfamily_read_level_analysis

##########################
# PREP DATA FOR ANALYSIS
conditions <- conf$levels
condition1 <- conditions[1]
condition2 <- conditions[2]
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name

sample_grs <- list()
for (sample_name in samples) {
    df <- read_table(grep(sprintf("/%s/", sample_name), inputs$bedmethlpaths, value = TRUE), col_names = FALSE)
    df_m <- df %>% filter(X4 == "a")
    rm(df)
    gr <- GRanges(
        seqnames = df_m$X1,
        ranges = IRanges(start = df_m$X2, end = df_m$X2),
        cov = df_m$X10,
        pctM = as.double(df_m$X11)
    )
    gr$sample <- sample_name
    gr$condition <- sample_table[sample_table$sample_name == sample_name, ]$condition
    sample_grs[[sample_name]] <- gr
}

grs <- Reduce(c, sample_grs)
rm(sample_grs)
# filter out low coverage and ensure that all samples have the same cpgs
grs <- grs[grs$cov > MINIMUMCOVERAGE]
grsdf <- tibble(as.data.frame(grs))
grsdf %$% seqnames %>% unique()
write_delim(grsdf %>% filter(grepl("*nonref*", seqnames)), "lrna/Rintermediates/grsdf_nonref.tsv", col_names = TRUE)
grsdf$seqnames <- factor(grsdf$seqnames, levels = chromosomesAll)
seqnames <- grsdf$seqnames
start <- grsdf$start
end <- grsdf$end
pos <- paste0(seqnames, "_", start, "_", end)
grsdf$pos <- pos

grsdfuntidy <- grsdf %>%
    pivot_wider(id_cols = c("pos", "seqnames"), names_from = "sample", values_from = "pctM", names_prefix = "pctM") %>%
    drop_na()

grsinboth <- grsdfuntidy %>% pull(pos)
# TODO pivot wider introduces sample names as variables, this needs to be addressed
grsdffiltered <- grsdf %>%
    filter(pos %in% grsinboth)
grsdf <- grsdffiltered
rm(grsdffiltered)
rm(grsdfuntidy)

write_delim(grsdf, "lrna/Rintermediates/grsdf.tsv", col_names = TRUE)


# SETTING UP SOME SUBSETS FOR EXPLORATION
set.seed(75)
grsdfs <- grsdf %>%
    group_by(sample, seqnames) %>%
    slice_sample(n = 1000)
grss <- GRanges(grsdfs)
write_delim(grsdfs, "lrna/Rintermediates/grsdfsmall.tsv", col_names = TRUE)


####################
## RTEs
dir.create("lrna/results/plots/rte", showWarnings = FALSE)

## Load Data and add annotations
r_annotation_fragmentsjoined <- read_csv(conf$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(conf$r_repeatmasker_annotation)
r_repeatmasker_annotation %$% ltr_viral_status %>% unique()
rmann <- left_join(r_annotation_fragmentsjoined, r_repeatmasker_annotation)
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

# FULL ELEMENTS


rmannfinal <- tibble(as.data.frame(RM))
RMdf <- rmannfinal
write_delim(RMdf, "lrna/Rintermediates/RMdf.tsv", col_names = TRUE)
# RMdf <- read_delim("lrna/Rintermediates/RMdf.tsv", col_names = TRUE)

grouping_var <- "rte_subfamily"
rte_frame <- GRanges(RMdf %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))
mbo <- mergeByOverlaps(grs, rte_frame)
methdf <- mbo$grs %>%
    as.data.frame() %>%
    tibble()
rte_only_frame <- mbo$rte_frame %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
rtedf <- bind_cols(methdf, rte_only_frame)
write_delim(rtedf, "lrna/Rintermediates/rtedf.tsv", col_names = TRUE)
# rtedf <- read_delim("lrna/Rintermediates/rtedf.tsv", col_names = TRUE)
perelementdf <- rtedf %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(RMdf %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width))

perelementdf <- perelementdf %>% filter(!is.na(rte_length_req))


write_delim(perelementdf, "lrna/Rintermediates/perelementdf.tsv", col_names = TRUE)
# perelementdf <- read_delim("lrna/Rintermediates/perelementdf.tsv", col_names = TRUE)


# RTE PROMOTERS
# annotate whether full length elements promoters overlap DMRs
flelement <- rmann %>% filter(rte_length_req == "FL")
rmann %$% rte_length_req %>% table()

flSINE <- flelement %>% filter(rte_superfamily == "SINE")
flLINE <- flelement %>% filter(rte_superfamily == "LINE")
rmann %$% ltr_viral_status %>% unique()
flFl_Provirus_5LTR <- flelement %>%
    filter(str_detect(gene_id, "LTR")) %>%
    filter(ltr_viral_status == "5'LTR (FL Int)")

flSINEgrs <- GRanges(flSINE)
flLINE5UTRgrs <- GRanges(flLINE) %>% resize(909)
flFl_Provirus_5LTRgrs <- GRanges(flFl_Provirus_5LTR)
flRTEpromotergrs <- c(c(flSINEgrs, flLINE5UTRgrs), flFl_Provirus_5LTRgrs)

flRTEpromotergrsfinal <- flRTEpromotergrs
flRTEpromoterdfinal <- tibble(as.data.frame(flRTEpromotergrsfinal))
flRTEpromoterdfGOOD <- flRTEpromoterdfinal
flRTEpromoter <- flRTEpromoterdfGOOD
write_delim(flRTEpromoter, "lrna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)
# flRTEpromoter <- read_delim("lrna/Rintermediates/flRTEpromoter.tsv", col_names = TRUE)

grouping_var <- "rte_subfamily"
rte_frame <- GRanges(flRTEpromoter %>% filter(!!sym(grouping_var) != "Other") %>% filter(rte_length_req == "FL"))
mbo <- mergeByOverlaps(grs, rte_frame)
methdf <- mbo$grs %>%
    as.data.frame() %>%
    tibble()
rte_only_frame <- mbo$rte_frame %>%
    as.data.frame() %>%
    tibble() %>%
    dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width)
rtedf_promoters <- bind_cols(methdf, rte_only_frame)

write_delim(rtedf_promoters, "lrna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
# rtedf_promoters <- read_delim("lrna/Rintermediates/rtedf_promoters.tsv", col_names = TRUE)
perelementdf_promoters <- rtedf_promoters %>%
    filter(cov > MINIMUMCOVERAGE) %>%
    group_by(gene_id, sample, condition) %>%
    summarize(mean_meth = mean(pctM)) %>%
    left_join(flRTEpromoter %>% dplyr::rename(rte_seqnames = seqnames, rte_start = start, rte_end = end, rte_strand = strand, rte_width = width))

write_delim(perelementdf_promoters, "lrna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)
# perelementdf_promoters <- read_delim("lrna/Rintermediates/perelementdf_promoters.tsv", col_names = TRUE)


############## Read density
rm(grs)
mem_used()

dir.create("lrna/results/plots/reads")
library(Biostrings)

readslist <- list()
for (region in conf$rte_subfamily_read_level_analysis) {
    for (sample_name in samples) {
        df <- read_delim(
            grep(region,
                grep(sprintf("/%s/", sample_name),
                    inputs$read_mods,
                    value = TRUE
                ),
                value = TRUE
            )
        )
        df$region <- region
        df$sample <- sample_name
        df$condition <- sample_table[sample_table$sample_name == sample_name, "condition"]
        grs <- GRanges(df %>% dplyr::rename(seqnames = chrom, start = ref_position, strand = ref_strand) %>% mutate(end = start))
        eoi <- import(paste0("aref/A.REF_annotations/A.REF_rte_beds/", region, ".bed"))
        mbo <- mergeByOverlaps(grs, eoi)
        df1 <- as.data.frame(mbo) %>%
            tibble() %>%
            dplyr::select(starts_with("grs"), name) %>%
            dplyr::rename(gene_id = name)
        colnames(df1) <- gsub("grs.", "", colnames(df1))
        readslist <- c(readslist, list(df1))
    }
}
reads <- Reduce(rbind, readslist)
write_delim(reads, "lrna/Rintermediates/reads_context_all.tsv", col_names = TRUE)
# reads <- read_delim("lrna/Rintermediates/reads_context_cpg.tsv", col_names = TRUE)
