#
# _common.R — shared settings for plot.R (style after Claus O. Wilke, Fundamentals of Data Visualization)
#   v3 (2026-08-26) — when editing, bump the version by 1 and update the date.
#     v3: register FONT in the PostScript/PDF font DB — removes metric-lookup warnings under Rscript.
#     v2: save_fig() background changed from `background =` to `bg =` (transparent-PNG bug, see ★ below).
#
#   Usage: place this file next to plot.R and, at the top of plot.R,
#             source("_common.R")
#
#   Provides
#     - path constants : ROOT / DATA(= data/) / FIGDIR(= figs/) / INPUT(= input/)
#     - palette        : OKABE (color-blind safe)
#     - theme          : th_default() / th_grid_h() / th_grid_hv() / th_top() / th_bottom()
#     - saving         : save_fig()  → figs/<name>.pdf (vector, main) + figs/<name>.png (preview)
#
#   Required packages: ggplot2, cowplot, ragg  (optional: readr, dplyr, latex2exp)
#

suppressMessages({
  library(ggplot2)
  library(cowplot)
})

# ── paths (resolved relative to the script location, so it works from any working directory) ──
.args <- commandArgs(FALSE)
.fa   <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT   <- if (length(.fa)) dirname(normalizePath(.fa)) else normalizePath(getwd())
DATA   <- file.path(ROOT, "data")     # csv files written by run.jl
INPUT  <- file.path(ROOT, "input")    # input data
FIGDIR <- file.path(ROOT, "figs")     # figure output
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# helper that reads data/<stem>.csv
read_data <- function(stem, ...) {
  p <- if (grepl("\\.csv$", stem)) file.path(DATA, stem) else file.path(DATA, paste0(stem, ".csv"))
  if (!file.exists(p)) stop("File not found: ", p, "  — check that run.jl has been run first", call. = FALSE)
  readr::read_csv(p, show_col_types = FALSE, ...)
}

# ── colour (Okabe-Ito) ───────────────────────────────────────────────────────
# Where colour carries no meaning, use black / grey instead.
OKABE <- c(orange     = "#E69F00", skyblue = "#56B4E9", green      = "#009E73",
           yellow     = "#F0E442", blue    = "#0072B2", vermillion = "#D55E00",
           purple     = "#CC79A7", black   = "#000000")
GREY_REF    <- "grey50"   # reference line (e.g. R = 1)
GREY_RIBBON <- "grey85"   # confidence-interval ribbon
GREY_BAR    <- "grey35"   # neutral bar / needle

okabe <- function(...) unname(OKABE[c(...)])
scale_colour_okabe <- function(...) scale_colour_manual(values = unname(OKABE), ...)
scale_fill_okabe   <- function(...) scale_fill_manual(values   = unname(OKABE), ...)

# ── font ─────────────────────────────────────────────────────────────────────
# macOS default. If a journal requires a specific font, change it only here.
FONT <- "Arial"
BASE_SIZE <- 11

# cairo_pdf · cairo_ps round the glyph advance of small text, so letters such as t·e in "commute" touch or drift apart
# (disappears when the same figure is drawn at 2× size). Vector output (pdf · eps) therefore draws text as FreeType
# outlines via showtext — size unchanged, spacing exact. png (ragg) is already correct and does not use it.
.ST_OK <- requireNamespace("showtext", quietly = TRUE)
if (.ST_OK) {
  .AF <- "/System/Library/Fonts/Supplemental/"
  sysfonts::font_add(FONT, regular = paste0(.AF, "Arial.ttf"), bold = paste0(.AF, "Arial Bold.ttf"),
                     italic = paste0(.AF, "Arial Italic.ttf"), bolditalic = paste0(.AF, "Arial Bold Italic.ttf"),
                     symbol = paste0(.AF, "STIXGeneralItalic.otf"))   # Greek letters in plotmath (italic, as in the text equations)
}
with_outlined_text <- function(expr) {
  if (.ST_OK) { showtext::showtext_auto(TRUE); on.exit(showtext::showtext_auto(FALSE)) }
  force(expr)
}

# Under Rscript the default device is classic postscript-based, so FONT metrics are looked up in its font DB.
# Arial is not registered there, so at every step where ggplot measures grob sizes
# ("cowplot::get_legend()", "ggplotGrob()", plot_grid alignment, etc.)
#   font family 'Arial' not found in PostScript font database
# is emitted hundreds of times. **Actual output (cairo_pdf · ragg) is unaffected** — it is purely a metric-lookup
# warning. Registering Helvetica (practically identical metrics) under the FONT name silences it.
local({
  reg <- function(f, db) {
    if (FONT %in% names(db())) return(invisible(NULL))
    a <- list(db()$Helvetica); names(a) <- FONT
    try(suppressWarnings(do.call(f, a)), silent = TRUE)
  }
  reg(grDevices::postscriptFonts, grDevices::postscriptFonts)
  reg(grDevices::pdfFonts,        grDevices::pdfFonts)
})

# ── theme ────────────────────────────────────────────────────────────────────
# Wilke: keep grid lines minimal. Time series / continuous y values get horizontal grid only (th_grid_h).
.th_common <- function(size = BASE_SIZE) {
  theme(axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.margin  = margin(6, 10, 6, 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
}
th_grid_h  <- function(size = BASE_SIZE)
  theme_minimal_hgrid(font_size = size, font_family = FONT) +
  theme(axis.line.x = element_line(linewidth = 0.4)) + .th_common(size)
th_grid_hv <- function(size = BASE_SIZE)
  theme_minimal_grid(font_size = size, font_family = FONT) + .th_common(size)
th_default <- function(size = BASE_SIZE)
  theme_half_open(font_size = size, font_family = FONT) + .th_common(size)

# For vertically stacked multi-panel figures: top panels drop x-axis labels, only the bottom panel keeps them.
th_top    <- function(size = BASE_SIZE) th_grid_h(size) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(4, 10, 2, 6))
th_bottom <- function(size = BASE_SIZE) th_grid_h(size) +
  theme(plot.margin = margin(2, 10, 6, 6))

# Common math labels (plotmath). latex2exp::TeX("$R_t$") also works if needed.
lab_Rt <- expression(italic(R)[t])

# ── saving ───────────────────────────────────────────────────────────────────
# Default: PDF (vector, main) + PNG (preview).
#   PDF  cairo_pdf — LaTeX math · subscripts · italics stay vector and remain sharp when zoomed.
#   PNG  ragg 400 dpi, white background — needed because markdown cannot render pdf inline.
#   TIFF for journal submission. ggsave produces an alpha channel (RGBA), which violates the spec, so
#        the device is opened directly and flattened onto a white background.
#   EPS  cairo_ps, fonts subset-embedded. Only when a journal requires EPS.
#
#   ★ The background must be given as `bg =` (`background =` does not work).  [v2, 2026-08-26]
#     ggsave() computes its own `bg` from the theme's `plot.background` fill and passes it to the device;
#     the cowplot themes used here (theme_minimal_hgrid · theme_half_open) have
#     `plot.background = element_blank()`, so that value becomes NA (transparent) and **overrides** the device's
#     `background =` argument. The result is an RGBA(0,0,0,0) transparent PNG, and **on a dark-mode viewer the
#     background shows black and black text · axes vanish**.
#     `bg = "white"` replaces the value ggsave computes, so it reliably takes effect.
#     Check: PNG header colortype 2 (RGB, no alpha) is fine; 6 (RGBA) means transparency crept in.
#
#   width/height in inches. single column ≈ 3.5, double column ≈ 7.0–7.5.
save_fig <- function(plot, name, width = 5.8, height = 3.5,
                     pdf = TRUE, png = TRUE, tiff = FALSE, eps = FALSE,
                     dpi_png = 400, dpi_tiff = 300, bg = "white") {
  stem <- if (grepl("/", name)) name else file.path(FIGDIR, name)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  made <- character(0)

  if (pdf) {
    with_outlined_text(ggsave(paste0(stem, ".pdf"), plot, device = cairo_pdf,
                              width = width, height = height, family = FONT, bg = bg))
    made <- c(made, "pdf")
  }
  if (png) {
    ggsave(paste0(stem, ".png"), plot, device = ragg::agg_png,
           width = width, height = height, dpi = dpi_png, bg = bg)
    made <- c(made, "png")
  }
  if (tiff) {
    # the device is opened directly, bypassing ggsave's bg computation → pass background =.
    ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height,
                   units = "in", res = dpi_tiff, compression = "lzw",
                   background = bg)
    print(plot); grDevices::dev.off()
    made <- c(made, "tiff")
  }
  if (eps) {
    with_outlined_text(ggsave(paste0(stem, ".eps"), plot, device = cairo_ps,
                              width = width, height = height, family = FONT,
                              fallback_resolution = 1200, bg = bg))
    made <- c(made, "eps")
  }
  cat("Saved: ", basename(stem), " [", paste(made, collapse = ", "), "]  ",
      width, "x", height, " in\n", sep = "")
  invisible(stem)
}
