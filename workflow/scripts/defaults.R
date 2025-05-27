# PACKAGES
library("stringr")
library("dplyr")
library("ggplot2")
library("readr")
library("magrittr")
library("forcats")
library("purrr")
library("tibble")
library("tidyr")
library("cowplot")
library("paletteer")
# library("GenomicRanges")

# PATHS

# FUNCTIONS

get_repeat_annotations <- function(
    default_or_extended = "default",
    keep_non_central = TRUE,
    append_NI_samplename_modifier = FALSE) {
    rmannShared <- read_csv(confALL$aref$rmann_shared)
    if (confALL$aref$update_ref_with_tldr$response == "yes") {
        if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
            rmannSamples <- list()
            for (sample in confALL$aref$samples) {
                df <- read_csv(sprintf("aref/%s/%s_annotations/%s_rmann_nonref.csv", default_or_extended, sample, sample))
                df$sample_name <- sample
                if (append_NI_samplename_modifier) {
                    df <- df %>%
                        left_join(sample_table) %>%
                        mutate(gene_id = paste0(sample_name, "__", gene_id))
                }
                rmannSamples[[sample]] <- df
            }
            rmannnonref <- do.call(rbind, rmannSamples) %>% tibble()
            rmann <- bind_rows(rmannShared, rmannnonref)
            if (!keep_non_central) {
                rmann <<- rmann %>% filter(refstatus != "NonCentral")
            }
        } else if (confALL$aref$update_ref_with_tldr$per_sample == "no") {
            rmann <- rmannShared
            if (!keep_non_central) {
                rmann <<- rmann %>% filter(refstatus != "NonCentral")
            }
        }
    } else {
        rmann <- rmannShared
    }

    rmann <- rmann %>%
        mutate(req_integrative = factor(req_integrative, levels = c("Old Trnc", "Old FL", "Yng Trnc", "Yng FL", "Yng Intact"))) %>%
        mutate(ltr_viral_status = factor(ltr_viral_status, levels = c("Int (Has 5LTR)", "Int (No 5'LTR)", "5'LTR (FL Int)", "3'LTR (FL Int)", "5'LTR (Trnc Int)", "3'LTR (Trnc Int)", "LTR (Solo)", "Other")))

    return(rmann)
}

rm2granges <- function(regex, rm, filter = c("length > 0", "length > 0")) {
    tmp <- rm %>% filter(str_detect(gene_id, regex))
    tmp <- tmp %>% filter(!!!rlang::parse_exprs(filter))
    gr <- GRanges(seqnames = tmp$chr, ranges = IRanges(start = tmp$start, end = tmp$end), strand = tmp$strand, gene_id = tmp$gene_id)
    mcols(gr) <- tmp %>% select(-chr, -start, -end, -strand, -gene_id)
    return(gr)
}

pw <- function(to_be_printed) {
    print(to_be_printed, width = Inf)
}
pl <- function(to_be_printed) {
    print(to_be_printed, n = 200)
}
pwl <- function(to_be_printed) {
    print(to_be_printed, n = 200, width = Inf)
}
mysave <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, pl = p, store = store_var, raster = FALSE) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    if (raster == TRUE) {
        tryCatch(
            {
                png(gsub(".pdf", ".png", fn), width = w, height = h, units = "in", res = res)
                print(pl)
                dev.off()
                print(paste(getwd(), gsub(".pdf", ".png", fn), sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )
    } else {
        tryCatch(
            {
                cairo_pdf(fn, width = w, height = h, family = "Helvetica")
                print(pl)
                dev.off()
                print(paste(getwd(), fn, sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )
    }
}
store_var <- "no"
mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, pl = p, store = store_var, raster = FALSE, sf = NULL, sfm = NULL) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    if (raster == TRUE) {
        tryCatch(
            {
                png(gsub(".pdf", ".png", fn), width = w, height = h, units = "in", res = res)
                print(pl)
                dev.off()
                print(paste(getwd(), gsub(".pdf", ".png", fn), sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )
    } else {
        tryCatch(
            {
                cairo_pdf(fn, width = w, height = h, family = "Helvetica")
                print(pl)
                dev.off()
                print(paste(getwd(), fn, sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )

        if (!exists("mysaveandstoreplots")) {
            mysaveandstoreplots <<- list()
        }
        if (store == "yes") {
            mysaveandstoreplots[[fn]] <<- pl
            print("plot_stored!")
        }
    }
    if (!is.null(sf)) {
        write_delim(sf, gsub(".pdf", "_stats.tsv", fn), delim = "\t", col_names = TRUE)
        print(paste(getwd(), gsub(".pdf", "_stats.tsv", fn), sep = "/"))
    }
    if (!is.null(sfm)) {
        write_delim(sfm, gsub(".pdf", "_modelstats.tsv", fn), delim = "\t", col_names = TRUE)
        print(paste(getwd(), gsub(".pdf", "_modelstats.tsv", fn), sep = "/"))
    }
}


mss <- function(fn = "ztmp.pdf", w = 5, h = 5, wv = 4, hv = 4, res = 600, pl = p, store = store_var, raster = FALSE, sf = NULL, plus_void = FALSE) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    if (raster == TRUE) {
        tryCatch(
            {
                png(gsub(".pdf", ".png", fn), width = w, height = h, units = "in", res = res)
                print(pl)
                dev.off()
                print(paste(getwd(), gsub(".pdf", ".png", fn), sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )
        if (plus_void == TRUE) {
            tryCatch(
                {
                    pl_no_title_or_legend <- pl + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = "none")
                    fn_no_title_or_legend <- gsub(".pdf", "_void.pdf", fn)
                    png(gsub(".pdf", ".png", fn_no_title_or_legend), width = wv, height = hv, units = "in", res = res)
                    print(pl_no_title_or_legend)
                    dev.off()
                    print(paste(getwd(), gsub(".pdf", ".png", fn_no_title_or_legend), sep = "/"))
                },
                error = function(e) {
                    print("plot not saved")
                    print(e)
                    tryCatch(
                        {
                            dev.off()
                        },
                        error = function(e) {
                            print(e)
                        }
                    )
                }
            )
        }
    } else {
        tryCatch(
            {
                cairo_pdf(fn, width = w, height = h, family = "Helvetica")
                print(pl)
                dev.off()
                print(paste(getwd(), fn, sep = "/"))
            },
            error = function(e) {
                print("plot not saved")
                print(e)
                tryCatch(
                    {
                        dev.off()
                    },
                    error = function(e) {
                        print(e)
                    }
                )
            }
        )
        if (plus_void == TRUE) {
            tryCatch(
                {
                    pl_no_title_or_legend <- pl + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = "none")
                    fn_no_title_or_legend <- gsub(".pdf", "_void.pdf", fn)
                    cairo_pdf(fn_no_title_or_legend, width = wv, height = hv, family = "Helvetica")
                    print(pl_no_title_or_legend)
                    dev.off()
                    print(paste(getwd(), fn_no_title_or_legend, sep = "/"))
                },
                error = function(e) {
                    print("plot not saved")
                    print(e)
                    tryCatch(
                        {
                            dev.off()
                        },
                        error = function(e) {
                            print(e)
                        }
                    )
                }
            )
        }

        if (!exists("mysaveandstoreplots")) {
            mysaveandstoreplots <<- list()
        }
        if (store == "yes") {
            mysaveandstoreplots[[fn]] <<- pl
            print("plot_stored!")
        }
    }
    if (!is.null(sf)) {
        write_delim(sf, gsub(".pdf", "_stats.tsv", fn), delim = "\t", col_names = TRUE)
        print(paste(getwd(), gsub(".pdf", "_stats.tsv", fn), sep = "/"))
    }
}


write_mycsv <- function(table, fn) {
    dir.create(dirname(fn), recursive = TRUE)
    write_csv(table, fn)
}

# VARIABLES
chromosomesAll <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
chromosomes <- c(paste0("chr", 1:22), "chrX")
chromosomesNoX <- c(paste0("chr", 1:22))

# GGPLOT
# mypalette <- list(
#     scale_colour_paletteer_d("ggsci::default_aaas", direction = 1), scale_fill_paletteer_d("ggsci::default_aaas", direction = 1)
#     )


mtopen <- theme_cowplot(font_family = "helvetica")
mtopengrid <- theme_cowplot(font_family = "helvetica") + background_grid(minor = "none")
mtopengridh <- theme_cowplot(font_family = "helvetica") + background_grid(major = "y", minor = "none")
mtopengridv <- theme_cowplot(font_family = "helvetica") + background_grid(major = "x", minor = "none")

mtclosed <- theme_cowplot(font_family = "helvetica") +
    panel_border(color = "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm"))
mtclosedgrid <- theme_cowplot(font_family = "helvetica") +
    panel_border(color = "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm")) + background_grid(minor = "none")
mtclosedgridh <- theme_cowplot(font_family = "helvetica") +
    panel_border(color = "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm")) + background_grid(major = "y", minor = "none")
mtclosedgridv <- theme_cowplot(font_family = "helvetica") +
    panel_border(color = "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm")) + background_grid(major = "x", minor = "none")



anchorbar <- list(
    scale_y_continuous(expand = expansion(mult = c(0, .075)))
)
anchorh <- list(
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)),
    scale_x_continuous(expand = expansion(mult = c(0, .075)))
)
anchorv <- list(
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)),
    scale_y_continuous(expand = expansion(mult = c(0, .075)))
)
anchorb <- list(
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)),
    scale_y_continuous(expand = expansion(mult = c(0, .075))),
    scale_x_continuous(expand = expansion(mult = c(0, .075)))
)


# mytheme <- list(
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# my_palette <- paletteer::paletteer_d("ggsci::nrc_npg")

# mythemecolor <- list(
#     scale_fill_manual(values = my_palette),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemecolor1 <- list(
#     scale_color_manual(values = my_palette),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemefill <- list(
#     scale_fill_manual(values = my_palette),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )

# bluelight <- "#8399c9"
# blue <- "#3C5488FF"
# bluedark <- "#273759"
# orangelight <- "#F39B7FFF"
# orange <- "#e94716"
# orangedark <- "#75240b"
# mycolor <- "#00A087FF"
# mythemesamples <- list(
#     scale_fill_manual(values = rep(my_palette[4:5], each = 3)),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemesamplesdif <- list(
#     scale_fill_manual(values = c(bluelight, blue, bluedark, orangelight, orange, orangedark)),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemeconditions <- list(
#     scale_color_manual(values = my_palette[5:4]),
#     scale_fill_manual(values = my_palette[5:4]),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemecontrast <- list(
#     scale_fill_manual(values = my_palette[c(1, 2, 9)]),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
# mythemecontrastrev <- list(
#     scale_fill_manual(values = my_palette[c(2, 1, 9)]),
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
# )
