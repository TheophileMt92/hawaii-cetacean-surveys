#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 06 — Figures for the humpback density surface model
#
# Palette continues 02: dark ocean ground, amber for whales, the blue ramp
# reused for uncertainty. Everything is drawn in the projected coordinates of
# 05 (Lambert azimuthal equal-area, km) with coord_equal, so geom_raster stays
# on a regular grid and no layer needs reprojecting at draw time.
#
# Produces, in outputs/figures/:
#   fig_effort.png  fig_density.png  fig_cv.png  fig_depth.png  fig_panel.png
#
# Needs: ggplot2, patchwork, rnaturalearth
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(here); library(sf); library(ggplot2); library(patchwork)

P  <- readRDS(here("data", "derived", "dsm_prediction.rds"))
DF <- readRDS(here("data", "derived", "df_humpback.rds"))

FIG <- here("outputs", "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ── palette ──────────────────────────────────────────────────────────────────
BG    <- "#06101A"   # ground
LAND  <- "#1B2733"
INK   <- "#E9F1F5"   # titles
LABEL <- "#D2DEE6"
SOFT  <- "#B6C6D2"
RULE  <- "#9FD9E8"   # study-area boundary
PT    <- "#FFC24D"   # whales, as in 02
TRACK <- "#3E5A6B"

# Sequential, single hue, lightness monotonic (OKLab L 0.30 → 0.91).
DENS <- c("#111C26", "#3A2A12", "#5C3F14", "#84591A",
          "#AD7620", "#D49529", "#FFC24D", "#FFDE96")
# 02's bathymetry ramp, reused here for uncertainty.
CVR  <- c("#0E2739", "#163E55", "#1F5875", "#2E7796",
          "#4E9BB5", "#7CBFD4", "#B9E2EE")

# ── geometry ─────────────────────────────────────────────────────────────────
sf_to_df <- function(g) {
  cc <- as.data.frame(st_coordinates(g))
  idc <- setdiff(names(cc), c("X", "Y"))
  data.frame(x = cc$X, y = cc$Y,
             grp = do.call(paste, c(cc[idc], sep = "_")))
}

poly_df <- sf_to_df(P$poly)

# Coastline. s2 is switched off for the crop: the bbox is a planar rectangle and
# s2 treats it as a spherical polygon, which fails on some sf/s2 combinations.
old_s2 <- sf::sf_use_s2(FALSE)
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid() |>
  st_crop(st_bbox(c(xmin = -163, xmax = -152, ymin = 16, ymax = 25), crs = 4326)) |>
  st_transform(P$prj)
sf::sf_use_s2(old_s2)
land_df <- sf_to_df(land)

prj_xy <- function(lon, lat) {
  st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326) |>
    st_transform(P$prj) |> st_coordinates()
}

# Tracklines from the segment endpoints, and detections at their own positions.
seg <- DF$segments
e1 <- prj_xy(seg$lon1, seg$lat1); e2 <- prj_xy(seg$lon2, seg$lat2)
trk <- data.frame(x = e1[, 1], y = e1[, 2], xend = e2[, 1], yend = e2[, 2])

if (is.null(DF$obs$lon) || is.null(DF$obs$lat)) {
  stop("df_humpback.rds carries no detection positions. Re-run ",
       "R/04_detection_function.R, whose obs block now keeps lon/lat.")
}

keep   <- DF$obs$distance <= DF$truncation
det_xy <- prj_xy(DF$obs$lon[keep], DF$obs$lat[keep])
det    <- DF$obs[keep, ] |> mutate(x = det_xy[, 1], y = det_xy[, 2])

# ── per-cell uncertainty ─────────────────────────────────────────────────────
# On a log link the standard error of the linear predictor is, to first order,
# the coefficient of variation of the prediction. The detection-function and
# g(0) components are spatially constant, so the pattern is the GAM's; adding
# them in quadrature shows the floor they impose everywhere.
lp <- predict(P$model, P$grid, type = "link", se.fit = TRUE)
grid <- P$grid |>
  mutate(cv = sqrt(as.numeric(lp$se.fit)^2 +
                     P$cv_parts[["detection"]]^2 + P$cv_parts[["g0"]]^2))

cv_floor <- sqrt(P$cv_parts[["detection"]]^2 + P$cv_parts[["g0"]]^2)
cat("per-cell CV range:", round(range(grid$cv), 3),
    "| floor from detection + g(0):", round(cv_floor, 3), "\n")

# ── theme ────────────────────────────────────────────────────────────────────
th <- theme_void(base_family = "Helvetica") +
  theme(
    plot.background   = element_rect(fill = BG, colour = NA),
    panel.background  = element_rect(fill = BG, colour = NA),
    legend.background = element_rect(fill = BG, colour = NA),
    plot.title    = element_text(colour = INK,  size = 13, face = "bold",
                                 hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(colour = SOFT, size = 9.2, hjust = 0,
                                 margin = margin(b = 8), lineheight = 1.25),
    plot.caption  = element_text(colour = SOFT, size = 7.6, hjust = 0,
                                 margin = margin(t = 8), lineheight = 1.3),
    legend.title  = element_text(colour = LABEL, size = 8.4),
    legend.text   = element_text(colour = SOFT,  size = 7.8),
    legend.position = "bottom",
    legend.key.height = unit(6, "pt"),
    legend.key.width  = unit(52, "pt"),
    legend.margin = margin(t = 2),
    plot.margin   = margin(12, 14, 12, 14))

# One extent for every map, so effort, density and uncertainty are directly
# comparable. Without this each panel takes its limits from its own layers and
# the three polygons come out at three different scales.
bb  <- st_bbox(P$poly)
PAD <- 25                                     # km
XL  <- c(bb[["xmin"]] - PAD, bb[["xmax"]] + PAD)
YL  <- c(bb[["ymin"]] - PAD, bb[["ymax"]] + PAD)
ASP <- diff(YL) / diff(XL)                    # panel aspect every plot adopts

frame <- list(
  geom_polygon(data = land_df, aes(x, y, group = grp), fill = LAND, colour = NA),
  geom_path(data = poly_df, aes(x, y, group = grp),
            colour = RULE, linewidth = 0.35, alpha = 0.75),
  coord_equal(xlim = XL, ylim = YL, expand = FALSE), th)

guide_bar <- guide_colourbar(title.position = "top", ticks.colour = NA,
                             frame.colour = NA)

# ── A. effort and detections ─────────────────────────────────────────────────
p_eff <- ggplot() +
  geom_segment(data = trk, aes(x, y, xend = xend, yend = yend),
               colour = TRACK, linewidth = 0.28) +
  frame +
  geom_point(data = det, aes(x, y, size = size),
             colour = PT, alpha = 0.85, shape = 16) +
  scale_size_area(max_size = 4.6, name = "group size", breaks = c(1, 2, 4, 6)) +
  labs(title = "Survey effort and humpback detections",
       subtitle = sprintf("%s systematic segments, %s km on effort, %d detections",
                          nrow(seg), format(round(sum(seg$dist)), big.mark = ","),
                          nrow(det))) +
  guides(size = guide_legend(title.position = "top",
                             override.aes = list(colour = PT)))

# ── B. predicted density ─────────────────────────────────────────────────────
# Square-root colour scaling: density is strongly right-skewed and a linear
# ramp would render all but the peak as background.
pk <- grid |> slice_max(dens, n = 1)

p_dens <- ggplot() +
  geom_raster(data = grid, aes(x, y, fill = dens)) +
  frame +
  scale_fill_gradientn(colours = DENS, trans = "sqrt",
                       name = expression("animals km"^-2*" (sqrt scale)"),
                       labels = scales::label_number(accuracy = 0.01)) +
  annotate("point", x = pk$x, y = pk$y, shape = 21, size = 7,
           colour = INK, fill = NA, stroke = 0.35, alpha = 0.8) +
  annotate("segment", x = pk$x + 38, xend = pk$x + 105, y = pk$y + 38,
           yend = pk$y + 95, colour = INK, linewidth = 0.3, alpha = 0.8) +
  annotate("text", x = pk$x + 112, y = pk$y + 108, hjust = 0, vjust = 0,
           label = sprintf("peak %.2f animals km\u207b\u00b2", max(grid$dens)),
           colour = INK, size = 2.9, family = "Helvetica") +
  labs(title = "Predicted humpback whale density",
       subtitle = sprintf("%s animals (95%% CI %s–%s) over %s km²",
                          format(round(P$N), big.mark = ","),
                          format(round(P$ci[1]), big.mark = ","),
                          format(round(P$ci[2]), big.mark = ","),
                          format(round(nrow(grid) * P$cell^2), big.mark = ","))) +
  guides(fill = guide_bar)

# ── C. uncertainty ───────────────────────────────────────────────────────────
p_cv <- ggplot() +
  geom_raster(data = grid, aes(x, y, fill = cv)) +
  frame +
  scale_fill_gradientn(colours = CVR, name = "coefficient of variation") +
  labs(title = "Where the estimate is least certain",
       subtitle = sprintf("never below %.2f anywhere: g(0) alone contributes %.2f",
                          min(grid$cv), P$cv_parts[["g0"]])) +
  guides(fill = guide_bar)

# ── D. depth response ────────────────────────────────────────────────────────
d_rng <- P$depth_range
pd <- data.frame(depth = seq(d_rng[1], d_rng[2], length.out = 250),
                 x = median(P$segments$x),
                 y = median(P$segments$y), off.set = 1)
tm <- predict(P$model, pd, type = "terms", se.fit = TRUE)
k  <- grep("depth", colnames(tm$fit), fixed = TRUE)[1]
pd <- pd |> mutate(f = tm$fit[, k], se = tm$se.fit[, k])

p_depth <- ggplot(pd, aes(depth, f)) +
  geom_ribbon(aes(ymin = f - 2 * se, ymax = f + 2 * se),
              fill = PT, alpha = 0.16) +
  geom_line(colour = PT, linewidth = 0.8) +
  geom_rug(data = P$segments, aes(x = depth), inherit.aes = FALSE,
           sides = "b", colour = SOFT, alpha = 0.35, length = unit(3, "pt")) +
  scale_x_continuous(labels = function(v) format(-v, big.mark = ",")) +
  labs(title = "Depth response",
       subtitle = "partial effect on log density, ±2 SE; rug shows segment depths",
       x = "Depth (m)", y = "s(depth)") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    plot.background  = element_rect(fill = BG, colour = NA),
    panel.background = element_rect(fill = BG, colour = NA),
    panel.grid.major = element_line(colour = "#16232E", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text  = element_text(colour = SOFT, size = 7.8),
    axis.title = element_text(colour = LABEL, size = 8.4),
    plot.title = element_text(colour = INK, size = 13, face = "bold",
                              margin = margin(b = 3)),
    plot.subtitle = element_text(colour = SOFT, size = 9.2,
                                 margin = margin(b = 8)),
    aspect.ratio = ASP,
    plot.margin = margin(12, 14, 12, 14))

# ── write ────────────────────────────────────────────────────────────────────
# Height is set from the panel aspect plus room for titles and legend. Too
# short and coord_equal letterboxes the map, so each panel ends up sized by
# whatever vertical space its own legend happens to leave — which is what made
# the four panels disagree.
PW <- 9                                        # canvas width, inches
sv <- function(p, file, w = PW, h = (w - 1.6) * ASP + 2.1) {
  ggsave(file.path(FIG, file), p, width = w, height = h, dpi = 200, bg = BG)
}
sv(p_eff,   "fig_effort.png")
sv(p_dens,  "fig_density.png")
sv(p_cv,    "fig_cv.png")
sv(p_depth, "fig_depth.png")

CAP <- paste(
  "Density surface model of humpback whale (Megaptera novaeangliae),",
  "WHICEAS winter 2020. Tweedie GAM on 10 km segments, half-normal detection",
  "function truncated at 5.5 km, corrected for g(0) = 0.68 (Barlow 2015).",
  "Data: NOAA PIFSC, processed with LTabundR.")

panel <- (p_eff | p_dens) / (p_cv | p_depth) +
  plot_annotation(
    caption = CAP,
    theme = theme(plot.background = element_rect(fill = BG, colour = NA),
                  plot.caption = element_text(colour = SOFT, size = 8,
                                              hjust = 0, family = "Helvetica",
                                              margin = margin(t = 6, b = 4)),
                  plot.margin = margin(6, 8, 6, 8)))

PANEL_W <- 16
ggsave(file.path(FIG, "fig_panel.png"), panel,
       width = PANEL_W,
       height = 2 * ((PANEL_W / 2 - 1.6) * ASP + 2.1) + 0.5,
       dpi = 170, bg = BG)

cat("written to", FIG, "\n")
