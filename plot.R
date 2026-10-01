#!/usr/bin/env Rscript
#
# plot.R — 2018-01 … 2023-12 weekday · weekend · holiday population movement
#
#   Part 1  2022-09 — calibration (prior vs posterior) · 3 components · OD heatmaps · daytime presence-ratio maps · validation
#   Part 2  2018–2023 — monthly time series of the 3 weekday components (all ages · not age-standardised)
#
#   Input  : data/*.csv · data/od/od_radiation_202209_*.csv.gz  (produced by run.jl; this script does no computation)
#   Output : figs/<name>.pdf (main) + figs/<name>.png (preview)
#   Run    : Rscript plot.R
#
#   Figures use **English labels for the paper** and carry no in-figure captions — captions are in the paper.
#   **All result figures use the main seed (`radiation`).** The mobile-phone OD seed (`mobile`) appears only in the
#   validation figures (seed_check_*).

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(viridisLite)
})

.args <- commandArgs(FALSE); .fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
source(file.path(if (length(.fa)) dirname(normalizePath(.fa)) else getwd(), "_common.R"))

CANON    <- "canon"                                      # must match CANON in run.jl
SEED     <- "radiation"                                  # main seed
DTLEV    <- c("weekday", "weekend", "holiday")
AGE_LAB  <- c("00"="0-9","10"="10-14","15"="15-19","20"="20-24","25"="25-29","30"="30-34",
              "35"="35-39","40"="40-44","45"="45-49","50"="50-54","55"="55-59","60"="60-64",
              "65"="65-69","70"="70-74","75"="75+")
AGE3LEV  <- c("0-19", "20-64", "65+")
COMP     <- c("commute", "noncommute", "nonmove")
COMP_LAB <- c(commute = "commute", noncommute = "non-commute", nonmove = "non-move")
COMP_COL <- c(commute = okabe("blue"), noncommute = okabe("orange"), nonmove = "grey60")
# Component colours are **this single set in every figure** (decomp_by_daytype · decomp_by_age · cnc_inner_outer).
# inner/outer are separated by tone of the same colour — outer = full colour, inner = light tone mixed 55% with white.
tint <- function(col, f) { m <- grDevices::col2rgb(col) / 255; grDevices::rgb(t(m * (1 - f) + f)) }
SEED_LAB <- c(radiation = "Radiation seed (main)", mobile = "Mobile-phone OD seed (validation)")
PTOT     <- 51466658                                     # registered population 2022-09 (run.jl §4)

# facet strip — no grey box. Next to a heatmap panel a grey fill becomes a mid-brightness band that is
# "neither data nor margin" and draws the eye.
strip_plain <- function(size = 9)
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = size))
fdt <- function(x) factor(x, levels = DTLEV)
YLAB_PCT <- "Share of registered population (%)"
# Paper figures: drawn at manuscript width (6.25 in) and included in LaTeX **without scaling**,
# so text prints at the same size in every figure — axis titles · panel titles SD_FS pt, tick labels SD_FS × 12/14 pt (cowplot).
SD_W  <- 6.25
SD_FS <- 9


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ Part 1 — 2022-09                                                            ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# ══ Figure 6 — prior vs posterior (φ, ψ, χ) ══════════════════════════════════
pri  <- read_data("prior_marginals") %>% mutate(kind = "prior")
pos  <- read_data("posterior_marginals") %>% filter(variant == CANON) %>%
  select(param, value, weight) %>% mutate(kind = "posterior")
# panel title = parameter symbol (plotmath, label_parsed). Meaning is given in the figure legend.
PLAB <- c(phi = "italic(phi)", psi = "italic(psi)", chi = "italic(chi)")   # phi = curly φ (same glyph as \varphi in the text)
d1 <- bind_rows(pri, pos) %>%
  mutate(param = factor(PLAB[param], levels = PLAB),
         kind  = factor(kind, levels = c("prior", "posterior")))

# For psi·chi the posterior **equals** the prior (flat likelihood). A thin line over a thick one smears a
# grey edge that reads as a non-existent shift, so the posterior is drawn as a filled area with the prior as a **dotted** line
# on top — where they coincide the grey dots trace the blue area's outline, and only on the phi axis do they separate.
# Layer order matters — draw the posterior first (thick), then the prior dotted line **on top**.
save_fig(ggplot(d1, aes(value, weight)) +
           geom_area(data = filter(d1, kind == "posterior"),
                     fill = okabe("blue"), alpha = 0.18) +
           geom_line(data = filter(d1, kind == "posterior"),
                     aes(colour = kind, linetype = kind), linewidth = 1.2) +
           geom_line(data = filter(d1, kind == "prior"),
                     aes(colour = kind, linetype = kind), linewidth = 0.65) +
           facet_wrap(~ param, scales = "free", ncol = 3, labeller = label_parsed) +
           scale_colour_manual(NULL, values = c(prior = "grey25", posterior = okabe("blue"))) +
           scale_linetype_manual(NULL, values = c(prior = "21", posterior = "solid")) +
           scale_x_continuous(NULL) +
           scale_y_continuous("Probability per grid point", expand = expansion(mult = c(0, 0.05))) +
           th_grid_h(SD_FS) + strip_plain(SD_FS + 1) +
           theme(legend.position = "top"),
         "prior_posterior", width = SD_W, height = 2.7, eps = TRUE)


# ══ Figure 2 — 3 components for the three day types (95% CI) ═════════════════
ps  <- read_data("posterior_summary") %>% filter(variant == CANON)
fin <- read_data("decomp_share_final") %>% filter(seed == SEED) %>%
  group_by(daytype) %>%
  summarise(commute = 100 * sum(commute) / sum(registered),
            noncommute = 100 * sum(noncommute) / sum(registered),
            nonmove = 100 * sum(nonmove) / sum(registered), .groups = "drop") %>%
  pivot_longer(-daytype, names_to = "component", values_to = "pct")
# Whiskers are drawn **only where theta actually enters**. On weekdays all three components are independent of theta
# (lambda = 1 . omega = 1 . phi not applied), and weekend commute is fixed by construction via sum_i Chat = P lambda kappa,
# so its CI width is 0.001%p — such a hairline whisker would be misread as "measured very precisely",
# so it is not drawn at all. The remaining 5 whiskers are exactly where theta moves.
# Rows are hidden with alpha = 0 **rather than removed**. position_dodge splits slots by the number of rows actually at
# that x, so removing a row would shift the remaining whiskers onto the wrong bars.
ci <- ps %>% filter(grepl("^(weekend|holiday)_(commute|noncommute|nonmove)$", quantity)) %>%
  separate(quantity, c("daytype", "component"), sep = "_") %>%
  mutate(across(c(q025, median, q975), ~ 100 * .x / PTOT),
         show = ifelse(daytype == "weekend" & component == "commute", 0, 1))

fx <- function(d) d %>%
  mutate(daytype = fdt(daytype),
         component = factor(component, levels = COMP, labels = unname(COMP_LAB)))

save_fig(ggplot(fx(fin), aes(daytype, pct, fill = component)) +
           geom_col(width = 0.6, position = position_dodge(0.66)) +
           geom_errorbar(data = fx(ci),
                         aes(y = median, ymin = q025, ymax = q975, alpha = show),
                         width = 0.16, linewidth = 0.4, colour = "grey20",
                         position = position_dodge(0.66)) +
           scale_alpha_identity() +
           scale_fill_manual(NULL, values = unname(COMP_COL)) +
           scale_x_discrete(NULL) +
           scale_y_continuous(YLAB_PCT, expand = expansion(mult = c(0, 0.05))) +
           th_grid_h(SD_FS) +
           theme(legend.position = "top"),
         "decomp_by_daytype", width = SD_W, height = 3.3, eps = TRUE)


# ══ Figure 3 — 15 age groups × 3 components × 3 day types ═════════════════════
d3 <- read_data("decomp_share_long", col_types = cols(age_group = col_character())) %>%
  filter(seed == SEED) %>%
  mutate(age = factor(AGE_LAB[age_group], levels = unname(AGE_LAB)),
         component = factor(component, levels = COMP),
         daytype = fdt(daytype))

save_fig(ggplot(d3, aes(age, pct, fill = component)) +
           geom_col(width = 0.8) +
           facet_wrap(~ daytype, ncol = 1) +
           scale_fill_manual(values = COMP_COL, labels = COMP_LAB) +
           scale_x_discrete("Age group") +
           scale_y_continuous(YLAB_PCT,
                              expand = expansion(mult = c(0, 0.02))) +
           th_grid_h(SD_FS) + strip_plain(SD_FS) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top"),
         "decomp_by_age", width = SD_W, height = 6.0, eps = TRUE)


# ══ C · NC split into inner / outer (supplementary) ══════════════════════════
cnc <- read_data("cnc_share_daytype", col_types = cols(age_group = col_character())) %>%
  filter(seed == SEED) %>%
  group_by(daytype) %>%
  summarise(across(c(registered, commute_inner, commute_outer,
                     noncommute_inner, noncommute_outer, nonmove), sum), .groups = "drop") %>%
  pivot_longer(-c(daytype, registered)) %>%
  mutate(pct = 100 * value / registered,
         comp = factor(ifelse(grepl("^commute", name), "commute",
                       ifelse(grepl("^noncommute", name), "non-commute", "non-move")),
                       levels = unname(COMP_LAB)),
         daytype = fdt(daytype))

CNC_KEY <- c(commute_outer = "commute — to another district (outer)",
             commute_inner = "commute — within own district (inner)",
             noncommute_outer = "non-commute — to another district (outer)",
             noncommute_inner = "non-commute — within own district (inner)",
             nonmove = "non-move")
CNC_COL <- c(commute_outer = COMP_COL[["commute"]], commute_inner = tint(COMP_COL[["commute"]], .55),
             noncommute_outer = COMP_COL[["noncommute"]], noncommute_inner = tint(COMP_COL[["noncommute"]], .55),
             nonmove = COMP_COL[["nonmove"]])
cnc <- cnc %>% mutate(key = factor(name, levels = names(CNC_KEY)))

save_fig(ggplot(cnc, aes(comp, pct, fill = key)) +
           geom_col(width = 0.65) + facet_wrap(~ daytype, ncol = 3) +
           scale_fill_manual(NULL, values = CNC_COL, labels = CNC_KEY) +
           scale_x_discrete(NULL) +
           scale_y_continuous(YLAB_PCT, expand = expansion(mult = c(0, 0.05))) +
           guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
           th_grid_h() + strip_plain() +
           theme(legend.position = "top", legend.text = element_text(size = 8),
                 axis.text.x = element_text(angle = 30, hjust = 1)),
         "cnc_inner_outer", width = 7.0, height = 4.4)


# ═══════════════════════════════════════════════════════════════════════════
# OD heatmaps — shared colour conventions
# ═══════════════════════════════════════════════════════════════════════════
# Values are persons spanning many orders of magnitude, so fill is log10. No grey (= missing) is used — exact 0,
# floating-point residuals and values below 0.0005 persons are all "effectively 0 persons", so they get the colour-bar floor colour.
ADM1_EN <- c("11"="Seoul","26"="Busan","27"="Daegu","28"="Incheon","29"="Gwangju","30"="Daejeon",
             "31"="Ulsan","36"="Sejong","41"="Gyeonggi","42"="Gangwon","43"="Chungbuk",
             "44"="Chungnam","45"="Jeonbuk","46"="Jeonnam","47"="Gyeongbuk","48"="Gyeongnam",
             "50"="Jeju")
ticks <- read_data("adm1_ticks", col_types = cols(adm1_code = col_character())) %>%
  mutate(name = ADM1_EN[adm1_code], span = i_max - i_min + 1) %>%
  arrange(i_min) %>% mutate(alt = rep(c("a", "b"), length.out = n()))
tick_lab   <- filter(ticks, span >= 14 | adm1_code == "28")   # 11 labels (large panels)
tick_small <- filter(ticks, span >= 20)                       # 5 labels (15-panel figures)

# adm1 ruler band — boundaries are shown as a band outside the panel rather than lines over the data.
BAND_COL <- c(a = "grey85", b = "grey45")
adm1_band <- function(pad) unlist(lapply(names(BAND_COL), function(g) {
  r <- filter(ticks, alt == g)
  list(annotate("rect", xmin = r$i_min - .5, xmax = r$i_max + .5,
                ymin = 0.5 - pad, ymax = -0.5, fill = BAND_COL[[g]], colour = NA),
       annotate("rect", ymin = r$i_min - .5, ymax = r$i_max + .5,
                xmin = 0.5 - pad, xmax = -0.5, fill = BAND_COL[[g]], colour = NA))
}), recursive = FALSE)

NSGG      <- 250
# heatmap palette — viridisLite rocket reversed so that 0 is paper-coloured (perceptually uniform).
HM_COLS   <- rev(rocket(256))
FLOOR_COL <- HM_COLS[1]
log_labels <- function(brk) parse(text = c(paste0("'\u2264 '*10^", brk[1]),   # ≤ as a plain character (avoids the symbol font)
                                           paste0("10^", brk[-1])))
fill_log <- function(lim)
  scale_fill_gradientn("Persons", colours = HM_COLS, limits = lim, oob = squish,
                       na.value = FLOOR_COL, breaks = seq(lim[1], floor(lim[2])),
                       labels = log_labels(seq(lim[1], floor(lim[2]))))
hm_axes <- function(pad, tk, sz)
  list(adm1_band(pad),
       scale_x_continuous("Origin  j", breaks = tk$i_mid, labels = tk$name,
                          expand = c(0, 0), limits = c(0.5 - pad, NSGG + 0.5)),
       scale_y_continuous("Destination  i", breaks = tk$i_mid, labels = tk$name,
                          expand = c(0, 0), limits = c(0.5 - pad, NSGG + 0.5)),
       coord_fixed(), th_default(8),
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = sz),
             axis.text.y = element_text(size = sz),
             panel.background = element_rect(fill = FLOOR_COL, colour = NA),
             panel.spacing = unit(2, "pt"),
             legend.position = "right", legend.key.height = unit(1.1, "cm")),
       strip_plain(9))
hm_common <- function(d, lim = c(0, 5.5), pad = 6, tk = tick_lab, sz = 6)
  ggplot(d, aes(jx, iy, fill = log10(pmax(value, 10^lim[1])))) +
    geom_raster() + fill_log(lim) + hm_axes(pad, tk, sz)

ccol <- cols(from_cd = col_character(), to_cd = col_character())
# raster index — ascending adm2 region code (same order as `S` in run.jl). Seoul is at the lower left.
lev  <- sort(unique(read_data("map_daytime_ratio",
                              col_types = cols(sgg_cd = col_character(),
                                               age = col_character()))$sgg_cd))


# ══ age groups 3 × 3 day types (one figure per kind; supplementary) ═══════════
d5 <- read_data("heatmap_agegroup3_daytype", col_types = ccol) %>%
  mutate(jx = match(from_cd, lev), iy = match(to_cd, lev),
         age3 = factor(age3, levels = AGE3LEV), daytype = fdt(daytype))
for (k in COMP[1:2]) {
  save_fig(hm_common(filter(d5, kind == k)) + facet_grid(daytype ~ age3),
           paste0("heatmap_agegroup3_daytype_", k), width = 7.5, height = 8.0)
}
# commute · non-commute in one figure (Figure 4). 3 day types × (2 kinds × 3 age groups),
# fitted to page width (7.2 in) — colour scale identical to the two figures above.
d5b <- filter(d5, kind %in% COMP[1:2]) %>%
  mutate(kind = factor(COMP_LAB[kind], levels = unname(COMP_LAB[1:2])))
# Repeating "non-commute" over six panels makes titles collide at manuscript width. So panels carry only the age group,
# and a group title (Commute · Non-commute) is placed once over each set of three. The two sets are slightly spaced apart.
d5c <- d5b %>% mutate(col = factor(paste(kind, age3), levels = as.vector(outer(unname(COMP_LAB[1:2]), AGE3LEV, paste)[, ]) |>
                                     matrix(nrow = 2) |> t() |> as.vector()))
col_lab <- setNames(rep(AGE3LEV, 2), levels(d5c$col))
p4 <- hm_common(d5c, sz = 5.5, tk = filter(tick_lab, adm1_code != "43")) +   # Chungbuk label dropped: it overlaps Chungnam
  facet_grid(daytype ~ col, labeller = labeller(col = col_lab)) +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.9)) +
  theme(legend.position = "bottom", legend.key.height = unit(0.25, "cm"), legend.key.width = unit(1.4, "cm"),
        panel.spacing.x = unit(c(2, 2, 9, 2, 2), "pt"),
        strip.text = element_text(size = SD_FS, face = "bold"), axis.title = element_text(size = SD_FS),
        legend.title = element_text(size = SD_FS), legend.text = element_text(size = SD_FS * 12 / 14))
g4 <- ggplotGrob(p4)
pc <- g4$layout[grepl("^panel-1-", g4$layout$name), ]; pc <- pc[order(pc$l), ]
g4 <- gtable::gtable_add_rows(g4, unit(1.5, "lines"), pos = 0)
for (k in 1:2) {
  cols <- pc$l[(3 * k - 2):(3 * k)]
  g4 <- gtable::gtable_add_grob(g4, grid::textGrob(c("Commute", "Non-commute")[k],
                                gp = grid::gpar(fontsize = SD_FS + 1, fontface = "bold", fontfamily = FONT)),
                                t = 1, l = min(cols), r = max(cols), clip = "off", name = paste0("group-title-", k))
}
save_fig(g4, "heatmap_agegroup3_daytype_both", width = SD_W, height = 4.4, eps = TRUE)
rm(d5, d5b, d5c, p4, g4); invisible(gc())


# ══ all 15 age groups (one figure per day type × kind; supplementary) ═════════
# Input is `data/od/od_radiation_202209_<daytype>.csv.gz` — it already holds the main-seed Ĉ · N̂C by age group in this
# form, so only the raster index is attached, without aggregation. **All six figures share the same colour scale.**
LOG_LIM15 <- c(0, 5)
hm15 <- function(d)
  ggplot(d, aes(jx, iy, fill = log10(pmax(value, 10^LOG_LIM15[1])))) +
    geom_raster() + fill_log(LOG_LIM15) + facet_wrap(~ age, ncol = 5) +
    hm_axes(10, tick_small, 5) + theme(strip.text = element_text(face = "bold", size = 8))

for (dt in DTLEV) {
  od <- read_csv(file.path(DATA, "od", paste0("od_", SEED, "_202209_", dt, ".csv.gz")), show_col_types = FALSE,
                  col_types = cols(kind = col_character(), from_cd = col_character(),
                                   to_cd = col_character(), age_group = col_character(),
                                   value = col_double())) %>%
    mutate(jx = match(from_cd, lev), iy = match(to_cd, lev),
           age = factor(AGE_LAB[age_group], levels = unname(AGE_LAB)))
  for (k in COMP[1:2]) {
    save_fig(hm15(filter(od, kind == k)),
             paste0("heatmap_agegroup15_", k, "_", dt), width = 7.5, height = 5.4)
  }
  rm(od); invisible(gc())
}


# ═══════════════════════════════════════════════════════════════════════════
# Daytime population / registered population maps (IPF destination marginal)
# ═══════════════════════════════════════════════════════════════════════════
# The IPF **origin** marginal satisfies Σ_i M̂^g_ij = P^g_j by definition, so its ratio is always 1. The meaningful one is
# the **destination (daytime) marginal** Σ_j M̂^g_ij; divided by the registered population it gives the day type's
# 15 h presence ratio. **It does not depend on seed or θ** — it is the very target that IPF fits.
CAP_ADM1 <- c("11", "28", "41")          # capital region = Seoul · Incheon · Gyeonggi
mpoly <- read_data("map_sgg_polygons",
                   col_types = cols(sgg_cd = col_character(), gid = col_character(),
                                    adm1 = col_character(), cap_outline = col_logical()))
mext  <- read_data("map_extent")
dayr  <- read_data("map_daytime_ratio",
                   col_types = cols(sgg_cd = col_character(), level = col_character(),
                                    age = col_character())) %>% mutate(daytype = fdt(daytype))

# Colour — diverging palette. Stops are dense around 1 so near-1 values still split blue/red. Values outside [0.70, 1.30] are squished.
RATIO_LIM   <- c(0.70, 1.30)
RATIO_COLS  <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#FDBBA1", "#E8613C", "#B2182B")
RATIO_STOPS <- c(0.70, 0.88, 0.9625, 1.00, 1.03, 1.096, 1.30)
ratio_scale <- scale_fill_gradientn(
  "Daytime ÷ registered", colours = RATIO_COLS, values = rescale(RATIO_STOPS, from = RATIO_LIM),
  limits = RATIO_LIM, oob = squish,
  breaks = c(0.70, 0.85, 1.00, 1.15, 1.30),
  labels = c("≤0.70", "0.85", "1.00", "1.15", "≥1.30"),
  guide = guide_colourbar(barwidth = 14, barheight = 0.5, title.position = "left", title.vjust = 1))

# Capital-region outline on the national panel — obtained by **drawing order** instead of computing a union.
cap_ring <- function(d) list(
  geom_polygon(data = filter(d, cap_outline), fill = NA, colour = "black", linewidth = 0.8),
  geom_polygon(data = filter(d, cap_outline), colour = NA))

# scope = national | capital · lvl = all | age3 | age15 · fct = facet
ratio_panel <- function(scope, lvl, fct, ttl, lw = 0.05, strip = 9, dts = DTLEV) {
  ex  <- filter(mext, scope == !!scope)
  pol <- if (scope == "capital") filter(mpoly, adm1 %in% CAP_ADM1) else mpoly
  sub <- filter(dayr, level == lvl, daytype %in% dts)
  lab <- if (lvl == "age15") unname(AGE_LAB) else sort(unique(sub$age))
  d <- pol %>%
    left_join(select(sub, sgg_cd, age, daytype, ratio), by = "sgg_cd",
              relationship = "many-to-many") %>%
    mutate(age = factor(if (lvl == "age15") AGE_LAB[age] else age, levels = lab))
  ggplot(d, aes(x, y, group = gid, fill = ratio)) +
    geom_polygon(colour = "white", linewidth = lw) +
    (if (scope == "national") cap_ring(d) else NULL) +
    ratio_scale + fct +
    coord_equal(xlim = c(ex$xmin, ex$xmax), ylim = c(ex$ymin, ex$ymax), expand = FALSE) +
    labs(title = ttl) + theme_void(base_size = 8) + strip_plain(strip) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 9),
          panel.spacing = unit(2, "pt"),
          plot.background = element_rect(fill = "white", colour = NA))
}
ratio_legend <- get_plot_component(
  ratio_panel("capital", "all", facet_grid(. ~ daytype), "") + theme(legend.position = "bottom"),
  "guide-box-bottom", return_all = TRUE)


# ── 3 age groups × 3 day types (one figure per scope; supplementary) ─────────
save_fig(plot_grid(ratio_panel("national", "age3", facet_grid(daytype ~ age),
                               "National — capital region outlined", lw = 0.04),
                   ratio_legend, ncol = 1, rel_heights = c(1, 0.05)),
         "map_daytime_ratio_agegroup3_national", width = 7.0, height = 8.6)
save_fig(plot_grid(ratio_panel("capital", "age3", facet_grid(daytype ~ age),
                               "Capital region — Seoul · Incheon · Gyeonggi", lw = 0.04),
                   ratio_legend, ncol = 1, rel_heights = c(1, 0.05)),
         "map_daytime_ratio_agegroup3_capital", width = 7.0, height = 7.0)

# ── all 15 age groups (one figure per day type, national; supplementary) ─────
for (dt in DTLEV) {
  save_fig(plot_grid(ratio_panel("national", "age15", facet_wrap(~ age, nrow = 3),
                                 paste0("National — ", dt), lw = 0.03, strip = 8, dts = dt),
                     ratio_legend, ncol = 1, rel_heights = c(1, 0.05)),
           paste0("map_daytime_ratio_agegroup15_", dt), width = 7.5, height = 7.8)
}


# ═══════════════════════════════════════════════════════════════════════════
# Seed comparison (validation) — radiation (main) vs mobile-phone OD (mobile)
# ═══════════════════════════════════════════════════════════════════════════
# ══ Figure 7 — district scatter ══════════════════════════════════════════════
sc <- read_data("sgg_scatter", col_types = cols(sgg_cd = col_character())) %>%
  filter(scope == "all") %>%
  pivot_longer(ends_with("_pct"), names_to = "comp", values_to = "pct") %>%
  mutate(comp = factor(sub("_pct$", "", comp), levels = COMP, labels = unname(COMP_LAB))) %>%
  select(daytype, seed, sgg_cd, comp, pct) %>%
  pivot_wider(names_from = seed, values_from = pct) %>%
  mutate(daytype = fdt(daytype))

# Both axes use **the same scale** — every panel 0–80% (holiday non-commute max 77.6%), coord_fixed. Points on the line mean agreement.
save_fig(ggplot(sc, aes(radiation, mobile)) +
           geom_abline(slope = 1, intercept = 0, colour = GREY_REF, linewidth = 0.4) +
           geom_point(alpha = 0.35, size = 0.7, colour = okabe("blue")) +
           facet_grid(daytype ~ comp) +
           scale_x_continuous("Radiation seed — share by district (%)",
                              limits = c(0, 80), breaks = seq(0, 80, 20)) +
           scale_y_continuous("Mobile-phone OD seed — share by district (%)",
                              limits = c(0, 80), breaks = seq(0, 80, 20)) +
           coord_fixed() + th_grid_hv(SD_FS) + strip_plain(SD_FS),
         "seed_check_scatter", width = SD_W, height = 6.4, eps = TRUE)

# ══ national 3 components: two seeds side by side (supplementary) ════════════
sd19 <- read_data("seed_comparison") %>%
  filter(level == "national") %>%
  mutate(across(c(radiation, mobile, diff_pp), as.numeric)) %>%
  mutate(comp = factor(sub("_pct$", "", metric), levels = COMP, labels = unname(COMP_LAB)),
         daytype = fdt(daytype)) %>%
  pivot_longer(c(radiation, mobile), names_to = "seed", values_to = "pct") %>%
  mutate(seed = factor(SEED_LAB[seed], levels = unname(SEED_LAB)))

save_fig(ggplot(sd19, aes(comp, pct, fill = seed)) +
           geom_col(width = 0.66, position = position_dodge(0.72)) +
           facet_wrap(~ daytype, ncol = 3) +
           scale_fill_manual(values = okabe("blue", "vermillion")) +
           scale_x_discrete(NULL) +
           scale_y_continuous(YLAB_PCT, expand = expansion(mult = c(0, 0.05))) +
           th_grid_h() + strip_plain() +
           theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1)),
         "seed_check_national", width = 7.5, height = 4.0)


# ═══════════════════════════════════════════════════════════════════════════
# Held-out comparison with tourism visitor data (Korea Tourism Data Lab visitor data)
# ═══════════════════════════════════════════════════════════════════════════
# ⚠️ **This is held-out.** The tourism data enter none of the likelihood, seed or marginals.
# rows = visited place (destination) · columns = residence (origin); diagonal set to 0, then row-normalised — normalisation
# removes the destination marginal that IPF fixes, leaving only the origin composition we actually model.

# ══ Figure 8 — adm1 17×17 row-normalised heatmap (tourism vs main NC) ═══════
vh <- read_data("visitor_heatmap_adm1",
                col_types = cols(visit_code = col_character(), res_code = col_character())) %>%
  mutate(daytype = fdt(daytype),
         source  = factor(source, levels = c("tourism", "model_NC"),
                          labels = c("Tourism visits (observed)", "Model non-commute")))
A1LAB <- ADM1_EN
A1ORD <- names(A1LAB)
vh <- vh %>% mutate(vx = match(res_code, A1ORD), vy = match(visit_code, A1ORD))

# The 17 × 17 cells are drawn with geom_tile (vector rectangles) — geom_raster becomes a small image that viewers smooth when enlarged, so it looks blurry.
save_fig(ggplot(vh, aes(vx, vy, fill = ifelse(prob > 0, prob, NA))) +
           geom_tile() + facet_grid(daytype ~ source) +
           scale_fill_gradientn("Probability", colours = HM_COLS, na.value = "grey55",
                                limits = c(0, max(vh$prob)), trans = "sqrt",
                                breaks = c(0.01, 0.05, 0.15, 0.3, 0.5, 0.8)) +
           scale_x_continuous("Origin (province of residence)", breaks = seq_along(A1ORD), labels = A1LAB,
                              expand = c(0, 0)) +
           scale_y_continuous("Destination (visited province)", breaks = seq_along(A1ORD), labels = A1LAB,
                              expand = c(0, 0)) +
           coord_fixed() +
           th_default(SD_FS) + strip_plain(SD_FS) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
                 axis.text.y = element_text(size = 6.5),
                 panel.spacing = unit(4, "pt"),
                 legend.position = "right", legend.key.height = unit(1.1, "cm"),
                 legend.title = element_text(size = SD_FS),
                 legend.text = element_text(size = SD_FS * 12 / 14)),
         "visitor_adm1_heatmap", width = SD_W, height = 8.0, eps = TRUE)


# ══ JSD by province (17 adm1): main NC vs tourism data (supplementary) ═══════
ja <- read_data("visitor_jsd_rows_adm1",
                col_types = cols(adm1_code = col_character())) %>%
  mutate(daytype = fdt(daytype),
         nm = unname(ADM1_EN[adm1_code]))
# Mass (visit share) goes into the axis label instead of a separate panel — here mass is auxiliary information
# telling how much weight to give "this province's JSD".
lab1 <- ja %>% filter(daytype == "weekday") %>%
  transmute(nm, lab = sprintf("%s  (%.1f%%)", nm, 100 * visit_share))
ord  <- ja %>% group_by(nm) %>% summarise(m = mean(jsd), .groups = "drop") %>%
  arrange(m) %>% left_join(lab1, by = "nm")
ja   <- ja %>% left_join(lab1, by = "nm") %>% mutate(lab = factor(lab, levels = ord$lab))
med  <- ja %>% group_by(daytype) %>% summarise(m = median(jsd), .groups = "drop")

save_fig(ggplot(ja, aes(jsd, lab, colour = daytype)) +
           geom_vline(data = med, aes(xintercept = m, colour = daytype),
                      linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
           geom_point(size = 2.4, alpha = 0.9) +
           scale_colour_manual(NULL, values = okabe("blue", "orange", "vermillion")) +
           scale_x_continuous("Row-wise Jensen–Shannon divergence (base 2)",
                              limits = c(0, max(ja$jsd) * 1.04),
                              expand = expansion(mult = c(0, 0.02))) +
           scale_y_discrete(NULL) +
           th_grid_hv() +
           theme(legend.position = "top", legend.key.size = unit(10, "pt")),
         "visitor_jsd_adm1", width = 6.0, height = 4.6)

cat("\n2022-09 figures done\n")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ Part 2 — 2018-01 … 2023-12                                                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# ══ Figure 5 — weekday all-age 3 components (not age-standardised) · 1 × 3 ═══
# Grey band = first wave (2020-02-23) – end of social distancing (2022-04-18). Solid lines = edges of the grey interval, dotted = in between.
# Purple = tightened, green = relaxed (two Okabe-Ito colours not overlapping the day-type colours blue · orange · vermillion).
# y axis — the three components sit at 60 · 25 · 14%, so a common axis flattens the changes, while free axes per panel make
# change sizes incomparable. Hence **the range width is fixed at 2.5%p** (only the centre differs by series) with ticks every 0.5%p.
# Monthly means are drawn as **steps covering the whole month** — a line joining mid-month points blurs which months precede an intervention date.
ev <- read_data("covid_events") %>% mutate(date = as.Date(date))
EV_COL <- c(tighten = okabe("purple"), relax = okabe("green"))
EV_LAB <- c(tighten = "Social distancing tightened", relax = "Social distancing relaxed")
Y_SPAN <- 2.5
mstart <- function(ym) as.Date(paste0(substr(ym, 1, 4), "-", substr(ym, 5, 6), "-01"))
wc <- read_data("weekday_components_crude", col_types = cols(ym = col_character())) %>%
  mutate(m0 = mstart(ym), comp = factor(COMP_LAB[component], levels = unname(COMP_LAB)))
# Row extending the last step month (2023-12) to 2024-01-01
wc_end <- wc %>% group_by(comp) %>% filter(m0 == max(m0)) %>% ungroup() %>% mutate(m0 = as.Date("2024-01-01"))
wlim <- wc %>% group_by(comp) %>%
  summarise(mid = (min(pct) + max(pct)) / 2, .groups = "drop") %>%
  mutate(lo = round((mid - Y_SPAN / 2) / 0.25) * 0.25, hi = lo + Y_SPAN) %>%
  select(comp, lo, hi) %>% pivot_longer(c(lo, hi), values_to = "pct") %>%
  mutate(m0 = as.Date("2018-01-01"))
ev_layers <- function() list(
  annotate("rect", xmin = min(ev$date), xmax = max(ev$date), ymin = -Inf, ymax = Inf, fill = "grey90"),
  geom_vline(data = ev, aes(xintercept = date, colour = kind,
                            linetype = ifelse(date %in% range(ev$date), "edge", "inner")),
             linewidth = 0.5, key_glyph = "path"),
  scale_colour_manual(NULL, values = EV_COL, labels = EV_LAB, breaks = c("tighten", "relax")),
  scale_linetype_manual(values = c(edge = "solid", inner = "dotted"), guide = "none"))

save_fig(ggplot(wc, aes(m0, pct)) + ev_layers() + geom_blank(data = wlim) +
           geom_step(data = bind_rows(wc, wc_end), direction = "hv", colour = "black", linewidth = 0.5) +
           facet_wrap(~ comp, nrow = 1, scales = "free_y") +
           scale_x_date(NULL, date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
           scale_y_continuous(YLAB_PCT, expand = c(0, 0),
                              breaks = function(l) seq(ceiling(l[1] * 2) / 2, floor(l[2] * 2) / 2, by = 0.5),
                              labels = label_number(accuracy = 0.1)) +
           th_grid_h(SD_FS) + strip_plain(SD_FS) +
           theme(legend.position = "top", legend.key.width = unit(18, "pt"),
                 legend.text = element_text(margin = margin(l = 4, r = 30))),
         "weekday_components_crude", width = SD_W, height = 2.6, eps = TRUE)

cat("\nDone\n")
