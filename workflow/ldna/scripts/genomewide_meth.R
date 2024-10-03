module_name <- "ldna"
conf <- configr::read.config(file = "conf/config.yaml")[[module_name]]
confALL <- configr::read.config(file = "conf/config.yaml")
source("workflow/scripts/defaults.R")
source("workflow/scripts/generate_colors_to_source.R")
### BSSEQ
library(readr)
library(bsseq)
library(BiocParallel)
set.seed(123)

tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            dss = sprintf("ldna/intermediates/%s/methylation/%s_CG_m_dss.tsv", sample_table$sample_name, sample_table$sample_name)
        ), env = globalenv())
        assign("outputs", list(plots = "ldna/results/plots/genomewide/genomewide_meth_plots.rds"), env = globalenv())
    }
)
# couldn't load dss, here is a function from that package
makeBSseqData <- function(dat, sampleNames) {
    n0 <- length(dat)

    if (missing(sampleNames)) {
        sampleNames <- paste("sample", 1:n0, sep = "")
    }

    alldat <- dat[[1]]
    if (any(alldat[, "N"] < alldat[, "X"], na.rm = TRUE)) {
        stop("Some methylation counts are greater than coverage.\n")
    }

    ix.X <- which(colnames(alldat) == "X")
    ix.N <- which(colnames(alldat) == "N")
    colnames(alldat)[ix.X] <- "X1"
    colnames(alldat)[ix.N] <- "N1"

    if (n0 > 1) { ## multiple replicates, merge data
        for (i in 2:n0) {
            thisdat <- dat[[i]]
            if (any(thisdat[, "N"] < thisdat[, "X"], na.rm = TRUE)) {
                stop("Some methylation counts are greater than coverage.\n")
            }

            ix.X <- which(colnames(thisdat) == "X")
            ix.N <- which(colnames(thisdat) == "N")
            colnames(thisdat)[c(ix.X, ix.N)] <- paste(c("X", "N"), i, sep = "")
            alldat <- merge(alldat, thisdat, all = TRUE)
        }
    }

    ## make BSseq object
    ix.X <- grep("X", colnames(alldat))
    ix.N <- grep("N", colnames(alldat))
    alldat[is.na(alldat)] <- 0
    M <- as.matrix(alldat[, ix.X, drop = FALSE])
    Cov <- as.matrix(alldat[, ix.N, drop = FALSE])
    colnames(M) <- colnames(Cov) <- sampleNames

    ## order CG sites according to positions
    idx <- split(1:length(alldat$chr), alldat$chr)
    M.ordered <- M
    Cov.ordered <- Cov
    pos.ordered <- alldat$pos

    for (i in seq(along = idx)) {
        thisidx <- idx[[i]]
        thispos <- alldat$pos[thisidx]
        dd <- diff(thispos)
        if (min(dd) < 0) { # not ordered
            warning(paste0("CG positions in chromosome ", names(idx)[i], " is not ordered. Reorder CG sites.\n"))
            iii <- order(thispos)
            M.ordered[thisidx, ] <- M[thisidx, ][iii, ]
            Cov.ordered[thisidx, ] <- Cov[thisidx, ][iii, ]
            pos.ordered[thisidx] <- alldat$pos[thisidx][iii]
        }
    }

    result <- BSseq(chr = alldat$chr, pos = pos.ordered, M = M.ordered, Cov = Cov.ordered)

    ##    result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M, Cov=Cov)

    result
}


sample_dfs <- list()
for (sample in sample_table$sample_name) {
    sample_dfs[[sample]] <- read_delim(grep(sprintf("/%s/", sample), inputs$dss, value = TRUE), col_names = TRUE)
}

BSobj <- makeBSseqData(sample_dfs, names(sample_dfs))
mParam <- MulticoreParam(workers = 12, progressbar = TRUE)

conditions <- conf$levels
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name
condition1 <- conditions[1]
condition2 <- conditions[2]
# condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name[1]
# condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name[1]
# need to adjust numbering given idiosyncracies of dmltest
# BSobj <- realize(BSobj, "HDF5Array")
# BSobj <- BSmooth(BSobj, BPPARAM = mParam)
# BSobj

chrSizesdf <- read_delim("aref/A.REF.fa.fai", delim = "\t", col_names = FALSE) %>%
    dplyr::select(X1, X2) %>%
    dplyr::rename(chr = X1, seqlengths = X2)
# centromere
chrSizesCSdf <- chrSizesdf %>% mutate(cumsum = cumsum(seqlengths))
cumsum_zerostart <- c(0, chrSizesCSdf$cumsum[1:length(chrSizesCSdf$cumsum) - 1])
chrSizesCSdf$cumsum_zerostart <- cumsum_zerostart

chrSizes <- setNames(chrSizesdf$seqlengths, chrSizesdf$chr)
bins <- tileGenome(chrSizes, tilewidth = 50000, cut.last.tile.in.chrom = T)
globalmeth <- getMeth(BSobj, bins, type = "raw", what = "perRegion")
head(globalmeth)
mcols(bins) <- globalmeth
globalmeth <- as.data.frame(bins) %>%
    tibble() %>%
    mutate(csum = cumsum(as.numeric(width)))

cytobands <- read.delim(conf$ref_cytobands, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "gieStain")) %>%
    tibble() %>%
    mutate(chr = factor(chr, levels = levels(globalmeth$seqnames)))

centromeres <- cytobands %>%
    filter(gieStain == "acen") %>%
    group_by(chr) %>%
    summarize(start = min(start), end = max(end)) %>%
    as.data.frame() %>%
    tibble() %>%
    left_join(chrSizesCSdf) %>%
    mutate(start_cumsum = cumsum_zerostart + start) %>%
    mutate(end_cumsum = cumsum_zerostart + end)

axis_set <- globalmeth %>%
    group_by(seqnames) %>%
    summarize(center = mean(csum))

for (sample in sample_table$sample_name) {
    p <- globalmeth %>%
        ggplot() +
        geom_point(aes(x = csum, y = !!sym(sample), color = factor(seqnames)), alpha = 0.2) +
        geom_segment(data = centromeres, aes(x = start_cumsum, xend = end_cumsum, y = 1.05, yend = 1.05), color = "red", size = 3) +
        scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$seqnames)))) +
        scale_x_continuous(label = axis_set$seqnames, breaks = axis_set$center) +
        labs(x = "", y = "Pct CpG Methylation", title = "Global Methylation") +
        mtopengridh +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

    mysaveandstore(pl = p, fn = sprintf("ldna/results/plots/genomewide/globalmeth_%s.pdf", sample), w = 14, h = 4, res = 300)
}


tryCatch(
    {
        globalmeth1 <- globalmeth %>%
            pivot_longer(cols = sample_table$sample_name, names_to = "sample_name", values_to = "methylation") %>%
            left_join(sample_table) %>%
            group_by(condition, csum, seqnames) %>%
            summarize(methylation = mean(methylation, na.rm = TRUE)) %>%
            pivot_wider(names_from = condition, values_from = methylation) %>%
            mutate(diff = !!sym(condition2) - !!sym(condition1)) %>%
            mutate(direction = factor(ifelse(diff > 0, "Hyper", "Hypo"), levels = c("Hyper", "Hypo")))

        p <- globalmeth1 %>%
            ggplot() +
            geom_point(aes(x = csum, y = diff, color = factor(seqnames), alpha = 0.2)) +
            geom_segment(data = centromeres, aes(x = start_cumsum, xend = end_cumsum, y = 1.05, yend = 1.05), color = "red", size = 3) +
            scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$seqnames)))) +
            scale_x_continuous(label = axis_set$seqnames, breaks = axis_set$center) +
            labs(x = "", y = sprintf("%s - %s CpG Methylation", condition2, condition1), title = "Global Methylation") +
            mtopengridh +
            theme(legend.position = "none") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        mysaveandstore(pl = p, fn = sprintf("ldna/results/plots/genomewide/globalmeth_diff.pdf", sample), w = 14, h = 4, res = 300)

        p <- globalmeth %>% ggplot(aes(x = direction)) +
            geom_bar(aes(fill = direction), position = "dodge") +
            labs(x = "", y = "Number of Regions", title = "Global Methylation") +
            mtopen +
            scale_methylation +
            theme(legend.position = "none")
        mysaveandstore(pl = p, fn = sprintf("ldna/results/plots/genomewide/globalmeth_diff_bar.pdf", sample), w = 4, h = 4, res = 300)
    },
    error = function(e) {
        message("No difference plot generated")
    }
)

x <- tibble(OUT = "")
write_tsv(x, file = outputs$plots)
