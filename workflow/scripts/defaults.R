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
    keep_non_central = TRUE) {
    if (default_or_extended == "default") {
        rmannShared <- read_csv("aref/default/A.REF_annotations/A.REF_rmann.csv")
    } else if (default_or_extended == "extended") {
        rmannShared <- read_csv("aref/extended/A.REF_annotations/A.REF_rmann.csv")
    } else {
        print("FAIL LOAD RMANN")
    }
    if (confALL$aref$update_ref_with_tldr$response == "yes") {
        if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
            rmannSamples <- list()
            for (sample in confALL$aref$samples) {
                df <- read_csv(sprintf("aref/%s/%s_annotations/%s_rmann_nonref.csv", default_or_extended, sample, sample))
                df$nonref_insert_sample_name <- sample
                rmannSamples[[sample]] <- df
            }
            rmannnonref <- bind_rows(rmannSamples) %>% tibble()
            rmann <- bind_rows(rmannShared, rmannnonref)
            if (!keep_non_central) {
                rmann <- rmann %>% filter(refstatus != "NonCentral")
            }
        } else if (confALL$aref$update_ref_with_tldr$per_sample == "no") {
            rmann <- rmannShared
            if (!keep_non_central) {
                rmann <- rmann %>% filter(refstatus != "NonCentral")
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

inspect_plot_structure <- function(plot_obj) {
    p_build <- ggplot_build(plot_obj)

    # Number of facets
    nfacets <- nrow(p_build$layout$layout)
    ncol_facets <- length(unique(p_build$layout$layout$COL))
    nrow_facets <- length(unique(p_build$layout$layout$ROW))

    # Number of bars per facet
    plot_data <- p_build$data[[1]]
    bars_per_facet <- length(unique(plot_data[["group"]]))

    # Extract titles and subtitles
    plot_labels <- plot_obj$labels
    title_present_t_or_f <- !is.null(plot_labels$title) && plot_labels$title != ""
    subtitle_present_t_or_f <- !is.null(plot_labels$subtitle) && plot_labels$subtitle != ""

    # legend titles
    legend_present_t_or_f <- !is.null(p_build$plot$labels$color) || !is.null(p_build$plot$labels$fill) || !is.null(p_build$plot$labels$shape)
    if (legend_present_t_or_f) {
        vars_to_check <- c(
            p_build$plot$labels$fill,
            p_build$plot$labels$color,
            p_build$plot$labels$shape
        ) %>% as.character()

        vals <- list()
        for (var in vars_to_check) {
            vals[[var]] <- plot_obj$data %>%
                pull(var) %>%
                unique() %>%
                as.character()
        }
        maxnchar_vals <- vals %>%
            map(~ max(nchar(.x))) %>%
            unlist()
        ## assuming ggplot puts the legends side by side
        max_legend_text_len <- sum(maxnchar_vals)
    }

    # Check if facet titles are present
    facet_title_present_t_or_f <- !is.null(plot_obj$facet$params$facets) ||
        !is.null(plot_obj$facet$params$rows) ||
        !is.null(plot_obj$facet$params$cols)

    # Extract axis text (x-axis and y-axis) and compute max length
    x_axis_labels <- p_build$layout$panel_params[[1]]$x$get_labels()
    y_axis_labels <- p_build$layout$panel_params[[1]]$y$get_labels()

    max_xaxis_text_len <- if (!is.null(x_axis_labels)) max(nchar(as.character(x_axis_labels)), na.rm = TRUE) else 0
    max_yaxis_text_len <- if (!is.null(y_axis_labels)) max(nchar(as.character(y_axis_labels)), na.rm = TRUE) else 0

    # detemine whether x axis is rotated
    theme_info <- ggplot2:::plot_theme(plot_obj)
    x_text_element <- theme_info$axis.text.x
    x_angle <- x_text_element$angle

    # Check for axis titles
    x_axis_title_present_t_or_f <- !is.null(plot_labels$x) && plot_labels$x != ""
    y_axis_title_present_t_or_f <- !is.null(plot_labels$y) && plot_labels$y != ""

    # Return extracted information as a list
    list(
        nfacets = nfacets,
        ncol_facets = ncol_facets,
        nrow_facets = nrow_facets,
        bars_per_facet = bars_per_facet,
        title_present = title_present_t_or_f,
        subtitle_present = subtitle_present_t_or_f,
        facet_title_present = facet_title_present_t_or_f,
        max_xaxis_text_len = max_xaxis_text_len,
        max_yaxis_text_len = max_yaxis_text_len,
        x_axis_title_present = x_axis_title_present_t_or_f,
        y_axis_title_present = y_axis_title_present_t_or_f,
        x_angle = x_angle,
        max_legend_text_len = max_legend_text_len,
        legend_present_t_or_f = legend_present_t_or_f
    )
}

get_print_dims <- function(plot_inspection_res, room_for_stats = FALSE) {
    nfacets <- plot_inspection_res$nfacets
    ncol_facets <- plot_inspection_res$ncol_facets
    nrow_facets <- plot_inspection_res$nrow_facets
    bars_per_facet <- plot_inspection_res$bars_per_facet
    title_present <- plot_inspection_res$title_present
    subtitle_present <- plot_inspection_res$subtitle_present
    facet_title_present <- plot_inspection_res$facet_title_present
    max_xaxis_text_len <- plot_inspection_res$max_xaxis_text_len
    max_yaxis_text_len <- plot_inspection_res$max_yaxis_text_len
    x_angle <- plot_inspection_res$x_angle
    x_axis_title_present <- plot_inspection_res$x_axis_title_present
    y_axis_title_present <- plot_inspection_res$y_axis_title_present
    legend_present_t_or_f <- plot_inspection_res$legend_present_t_or_f
    max_legend_text_len <- plot_inspection_res$max_legend_text_len

    width_intercept <- 1 + (y_axis_title_present * 2 / 8) + (max_yaxis_text_len / 10) + (!is.null(nrow_facets)) * 2 / 8
    width_legend <- legend_present_t_or_f * ((max_legend_text_len / 10))
    width_plot <- ncol_facets * bars_per_facet * 2 / 8
    total_width <- width_intercept + width_legend + width_plot

    height_intercept <- 1 + (x_axis_title_present * 2 / 8) + (title_present * 2 / 8) + (subtitle_present * 2 / 8) + ifelse(!(is.null(x_angle)), sin((x_angle) * pi / 180) * max_xaxis_text_len / 10, 2 / 8) + (!is.null(ncol_facets)) * 2 / 8
    height_plot <- ifelse(room_for_stats, nrow_facets * 2.25, nrow_facets * 1.5)
    total_height <- height_intercept + height_plot

    list(
        w = total_width,
        h = total_height
    )
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
mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, pl = p, store = store_var, raster = FALSE, sf = NULL, sfm = NULL, auto_width = FALSE, auto_height = FALSE, auto_dims_stats = FALSE) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)
    if (auto_width || auto_height) {
        auto_dims <- tryCatch(
            {
                get_print_dims(plot_inspection_res = inspect_plot_structure(plot_obj = pl), room_for_stats = auto_dims_stats)
            },
            error = function(e) {
                print("autodims error")
                NULL # Return NULL if an error occurs
            }
        )

        if (!is.null(auto_dims)) {
            if (auto_width) {
                w <- auto_dims$w
            }
            if (auto_height) {
                h <- auto_dims$h
            }
        }
    }
    # print("getprintdims")
    # print(sprintf("%s", w))
    # print(sprintf("%s", h))
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


make_void_plot <- function(pl, target_inches = 3, axis_text_size = 12) {
    # Calculate number of mono characters needed to fill target_inches
    # Mono char width at a given pt size: ~0.6 * pt_size in points, convert to inches (/72)
    char_width_in <- 0.6 * axis_text_size / 72
    pad_width <- ceiling(target_inches / char_width_in)

    pad_fn <- function(x) {
        formatted <- as.character(x)
        stringr::str_pad(formatted, width = pad_width, side = "right", pad = "$")
    }

    # Strip all non-data elements, set mono font, rotate x labels 90 degrees
    pl <- pl + theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14, family = "mono", margin = margin(0, 0, 1, 0, "mm")),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = "mono", size = axis_text_size),
        axis.text.y = element_text(family = "mono", size = axis_text_size)
    )

    # Modify existing scale label functions in-place to preserve breaks/limits/transforms
    x_modified <- FALSE
    y_modified <- FALSE
    for (sc in pl$scales$scales) {
        if ("x" %in% sc$aesthetics) {
            sc$labels <- pad_fn
            x_modified <- TRUE
        }
        if ("y" %in% sc$aesthetics) {
            sc$labels <- pad_fn
            y_modified <- TRUE
        }
    }

    # If no explicit scale exists, detect type and add one with padding
    if (!x_modified || !y_modified) {
        built <- ggplot_build(pl)
        if (!x_modified) {
            if (inherits(built$layout$panel_scales_x[[1]], "ScaleDiscretePosition")) {
                pl <- pl + scale_x_discrete(labels = pad_fn)
            } else {
                pl <- pl + scale_x_continuous(labels = pad_fn)
            }
        }
        if (!y_modified) {
            if (inherits(built$layout$panel_scales_y[[1]], "ScaleDiscretePosition")) {
                pl <- pl + scale_y_discrete(labels = pad_fn)
            } else {
                pl <- pl + scale_y_continuous(labels = pad_fn)
            }
        }
    }

    pl
}

mss <- function(fn = "ztmp.pdf", w = 5, h = 5, fw = NULL, fh = NULL, res = 600, pl = p, store = store_var, raster = FALSE, sf = NULL, plus_void = FALSE, void_target_inches = 4, axis_text_size = 12) {
    # fw/fh = per-panel data area (inches). mss detects facet layout
    # from the ggplot object and scales total width/height accordingly.
    # Void plots account for $-padded axis labels per panel (free scales)
    # or once (shared scales).

    # Detect facet layout if fw or fh provided and pl is a ggplot
    ncol_facet <- 1L
    nrow_facet <- 1L
    free_x <- FALSE
    free_y <- FALSE
    if ((!is.null(fw) || !is.null(fh)) && inherits(pl, "gg")) {
        built <- ggplot_build(pl)
        layout <- built$layout$layout
        ncol_facet <- max(layout$COL, 1L)
        nrow_facet <- max(layout$ROW, 1L)
        free_x <- isTRUE(pl$facet$params$free$x)
        free_y <- isTRUE(pl$facet$params$free$y)
    }

    if (!is.null(fw)) {
        panel_w <- fw * 1.5
        w <- ncol_facet * panel_w + 3
        if (free_y) {
            # Each column gets its own y-axis with $ padding
            wv <- ncol_facet * (panel_w + void_target_inches)
        } else {
            # Single y-axis padding + all panel columns
            wv <- void_target_inches + ncol_facet * panel_w
        }
    } else {
        wv <- w
    }
    if (!is.null(fh)) {
        panel_h <- fh * 1.5
        h <- nrow_facet * panel_h + 3
        if (free_x) {
            # Each row gets its own x-axis with $ padding
            hv <- nrow_facet * (panel_h + void_target_inches)
        } else {
            # Single x-axis padding + all panel rows
            hv <- void_target_inches + nrow_facet * panel_h
        }
    } else {
        hv <- h
    }
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    # Enforce consistent axis text size on both axes
    if (!is.null(axis_text_size) && inherits(pl, "gg")) {
        pl <- pl + theme(
            axis.text.x = element_text(size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size)
        )
    }

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
                    pl_void <- make_void_plot(pl, target_inches = void_target_inches, axis_text_size = axis_text_size)
                    fn_void <- gsub(".pdf", "_void.pdf", fn)
                    png(gsub(".pdf", ".png", fn_void), width = wv, height = hv, units = "in", res = res)
                    print(pl_void)
                    dev.off()
                    print(paste(getwd(), gsub(".pdf", ".png", fn_void), sep = "/"))
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
                    pl_void <- make_void_plot(pl, target_inches = void_target_inches, axis_text_size = axis_text_size)
                    fn_void <- gsub(".pdf", "_void.pdf", fn)
                    cairo_pdf(fn_void, width = wv, height = hv, family = "Helvetica")
                    print(pl_void)
                    dev.off()
                    print(paste(getwd(), fn_void, sep = "/"))
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


standardize_axis_labels_num <- function(digits_before = 3, digits_after = 3) {
    function(x) {
        # Construct a format string like "%06.3f"
        # Width = digits_before + digits_after + 1 (for decimal point)
        total_width <- digits_before + digits_after + 1
        fmt <- paste0("%0", total_width, ".", digits_after, "f")
        sprintf(fmt, x)
    }
}
standardize_axis_labels_chr <- function(width = NULL, side = c("left", "right"), pad = " ") {
    side <- match.arg(side)

    function(x) {
        x <- as.character(x)
        # If width not specified, use the longest string
        if (is.null(width)) width <- max(nchar(x), na.rm = TRUE)

        if (side == "left") {
            stringr::str_pad(x, width = width, side = "left", pad = pad)
        } else {
            stringr::str_pad(x, width = width, side = "right", pad = pad)
        }
    }
}

snN <- function(f, s = 0) {
    list(
        scale_y_continuous(labels = standardize_axis_labels_num(digits_before = f, digits_after = s)),
        theme(text = element_text(family = "mono"))
    )
}
snNab <- function(f, s = 0) {
    list(
        scale_y_continuous(labels = standardize_axis_labels_num(digits_before = f, digits_after = s), expand = expansion(mult = c(0, .075))),
        theme(text = element_text(family = "mono"))
    )
}

scN <- function(f) {
    list(
        scale_y_continuous(labels = standardize_axis_labels_chr(width = f, side = "left", pad = " ")),
        theme(text = element_text(family = "mono"))
    )
}

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
