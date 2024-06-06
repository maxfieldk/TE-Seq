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
rm2granges <- function(regex, rm, filter = c("length > 0", "length > 0")) {
    tmp <- rm %>% filter(str_detect(gene_id, regex))
    tmp <- tmp %>% filter(!!!rlang::parse_exprs(filter))
    gr <- GRanges(seqnames = tmp$chr, ranges = IRanges(start = tmp$start, end = tmp$end), strand = tmp$strand, gene_id = tmp$gene_id)
    mcols(gr) <- tmp %>% select(-chr, -start, -end, -strand, -gene_id)
    return(gr)
}

mysave <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, pl = p, store = store_var, raster = FALSE) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    if (raster == TRUE) {
            tryCatch({
            png(gsub(".pdf", ".png", fn), width = w, height = h, units = "in", res = res)
            print(pl)
            dev.off()
            print(paste(getwd(),gsub(".pdf", ".png", fn), sep = "/"))
            }
            , error = function(e) {
                print("plot not saved")
                print(e)
            }
            )
    } else {
    tryCatch(
        {
            cairo_pdf(fn, width = w, height = h, family = "Helvetica")
            print(pl)
            dev.off()
            print(paste(getwd(),fn, sep = "/"))
        },
        error = function(e) {
            print("plot not saved")
            print(e)
        }
    )
}
}

mysaveandstore <- function(fn = "ztmp.pdf", w = 5, h = 5, res = 600, pl = p, store = store_var, raster = FALSE) {
    dn <- dirname(fn)
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)

    if (raster == TRUE) {
            tryCatch({
            png(gsub(".pdf", ".png", fn), width = w, height = h, units = "in", res = res)
            print(pl)
            dev.off()
            print(paste(getwd(),gsub(".pdf", ".png", fn), sep = "/"))
            }
            , error = function(e) {
                print("plot not saved")
                print(e)
            }
            )
    } else {
    tryCatch(
        {
            cairo_pdf(fn, width = w, height = h, family = "Helvetica")
            print(pl)
            dev.off()
            print(paste(getwd(),fn, sep = "/"))
        },
        error = function(e) {
            print("plot not saved")
            print(e)
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
    panel_border(color= "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm"))
mtclosedgrid <- theme_cowplot(font_family = "helvetica") +
    panel_border(color= "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm")) + background_grid(minor = "none")
mtclosedgridh <- theme_cowplot(font_family = "helvetica") +
    panel_border(color= "black") +
    theme(axis.line = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(4, "mm")) + background_grid(major = "y", minor = "none")
mtclosedgridv <- theme_cowplot(font_family = "helvetica") +
    panel_border(color= "black") +
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

