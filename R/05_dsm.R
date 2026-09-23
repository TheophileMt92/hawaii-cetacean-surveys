#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 05 — Density surface model, humpback whale, WHICEAS winter 2020
#
# Fits a count-response GAM over 10 km segments with a Tweedie error structure,
# predicts onto a 5 km grid covering the WHICEAS study area, and corrects for
# incomplete trackline detection with g(0) = 0.68 (CV 0.36), the Beaufort-
# weighted value Bradford, Yano & Oleson (2022) take from Barlow (2015).
#
# Their design-based estimate for this survey — 2,975 (CV 0.40, 95% CI
# 1,407-6,291) over 403,822 km² — is an independent check, not a target.
#
# Produces:
#   data/derived/dsm_prediction.rds   grid, predictions, fitted model
#   outputs/dsm_models.csv            GAM comparison
#   outputs/abundance.csv             abundance with the uncertainty budget
#   outputs/figures/dsm_check.png     GAM diagnostics
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(here); library(sf); library(terra)
library(dsm); library(mgcv)

G0    <- 0.68   # trackline detection probability, Barlow (2015)
G0_CV <- 0.36
CELL  <- 5      # km, prediction grid resolution

# Equal-area projection centred on the study area. +ellps rather than +datum
# avoids a proj.db datum lookup, which some Anaconda installs shadow.
PRJ <- "+proj=laea +lat_0=20.5 +lon_0=-157.5 +ellps=WGS84 +units=km +no_defs"
stopifnot(!is.na(st_crs(PRJ)$proj4string))

fit <- readRDS(here("data", "derived", "df_humpback.rds"))
ddf <- fit$model
seg <- fit$segments
obs <- fit$obs |> filter(distance <= fit$truncation)

# ── study area ───────────────────────────────────────────────────────────────
cruz10 <- readRDS(here("data", "derived", "cruz_10km.rds"))
whi <- cruz10$settings$strata[["WHICEAS"]]

# Close the ring if the stratum table does not already repeat its first vertex;
# st_polygon() requires closure and produces invalid geometry without it.
ring <- as.matrix(whi[, c("Lon", "Lat")])
if (!isTRUE(all.equal(ring[1, ], ring[nrow(ring), ], check.attributes = FALSE))) {
  ring <- rbind(ring, ring[1, ])
}
stopifnot(nrow(ring) >= 4)

poly <- st_sfc(st_polygon(list(ring)), crs = 4326) |> st_transform(PRJ)
stopifnot(st_is_valid(poly))

area_km2 <- as.numeric(st_area(poly))
cat("WHICEAS polygon area:", round(area_km2),
    "km² | Bradford et al.: 403,822 km²\n")

# ── bathymetry ───────────────────────────────────────────────────────────────
# Cached to disk: marmap re-downloads otherwise, and the file is gitignored.
BATH <- here("data", "raw", "bathy_whiceas.rds")
if (!file.exists(BATH)) {
  bath <- marmap::getNOAA.bathy(lon1 = -162.5, lon2 = -152.5,
                                lat1 = 16.5,   lat2 = 24.5,
                                resolution = 2, keep = FALSE)
  saveRDS(bath, BATH)
} else bath <- readRDS(BATH)

bxyz <- marmap::as.xyz(bath); names(bxyz) <- c("x", "y", "z")
bath_r <- rast(bxyz, type = "xyz")            # lon/lat, no CRS set deliberately

depth_at <- function(lon, lat) {
  as.numeric(terra::extract(bath_r, cbind(lon, lat))[, 1])
}

# ── segment covariates ───────────────────────────────────────────────────────
seg_xy <- st_as_sf(seg, coords = c("mlon", "mlat"), crs = 4326, remove = FALSE) |>
  st_transform(PRJ) |> st_coordinates()

seg_dsm <- seg |>
  transmute(Sample.Label = seg_id,
            Effort = dist,
            x = seg_xy[, 1], y = seg_xy[, 2],
            lon = mlon, lat = mlat,
            depth = depth_at(mlon, mlat))

cat("segments with missing depth:", sum(is.na(seg_dsm$depth)),
    "| positive depth:", sum(seg_dsm$depth >= 0, na.rm = TRUE), "\n")
stopifnot(sum(is.na(seg_dsm$depth)) == 0)

# A handful of segments may clip an island shelf; depth is used as a smooth
# covariate so a positive value is meaningless rather than merely extreme.
seg_dsm <- seg_dsm |> filter(depth < 0)
obs     <- obs |> filter(Sample.Label %in% seg_dsm$Sample.Label)

cat("modelling on", nrow(seg_dsm), "segments,", nrow(obs), "detections,",
    sum(obs$size), "individuals\n")

# ── prediction grid ──────────────────────────────────────────────────────────
grid_pts <- st_make_grid(poly, cellsize = CELL, what = "centers") |>
  st_as_sf() |> st_filter(poly)

ll <- st_transform(grid_pts, 4326) |> st_coordinates()
xy <- st_coordinates(grid_pts)

grid <- tibble(x = xy[, 1], y = xy[, 2],
               lon = ll[, 1], lat = ll[, 2],
               depth = depth_at(ll[, 1], ll[, 2])) |>
  filter(!is.na(depth), depth < 0) |>
  mutate(off.set = CELL^2)

cat("grid cells:", nrow(grid),
    "| water area", round(nrow(grid) * CELL^2), "km²\n")

# Depth range must cover the grid, or the smooth is extrapolating.
cat("segment depth range:", round(range(seg_dsm$depth)), "\n")
cat("grid depth range:   ", round(range(grid$depth)), "\n")

# ── candidate models ─────────────────────────────────────────────────────────
# 72 detections over ~450 segments is sparse, so basis dimensions are kept low
# deliberately: k is an upper bound the REML penalty shrinks from, but with data
# this thin a generous k invites the smooth to chase individual encounters.
mk <- function(f) dsm(f, ddf.obj = ddf, segment.data = seg_dsm,
                      observation.data = obs, family = tw(), method = "REML")

mods <- list(
  depth        = mk(count ~ s(depth, k = 5)),
  space        = mk(count ~ s(x, y, k = 12)),
  `depth+space`= mk(count ~ s(depth, k = 5) + s(x, y, k = 12))
)

gcmp <- tibble(model = names(mods),
               aic = vapply(mods, AIC, numeric(1)),
               dev_expl = vapply(mods, function(m)
                 round(100 * summary(m)$dev.expl, 1), numeric(1)),
               edf = vapply(mods, function(m) round(sum(m$edf), 2), numeric(1))) |>
  mutate(delta_aic = round(aic - min(aic), 2), aic = round(aic, 2)) |>
  arrange(aic)

print(as.data.frame(gcmp))
readr::write_csv(gcmp, here("outputs", "dsm_models.csv"))

mod <- mods[[gcmp$model[1]]]
cat("selected:", gcmp$model[1], "\n")
print(summary(mod))

png(here("outputs", "figures", "dsm_check.png"), 1100, 900, res = 120)
par(mfrow = c(2, 2)); gam.check(mod)
dev.off()

# ── prediction ───────────────────────────────────────────────────────────────
# The shallowest segment sits in ~99 m of water, but the grid reaches the
# shoreline, and s(depth) is log-linear: extended unchecked it keeps raising
# density all the way inshore, placing ~14% of the animals in 0.9% of the area
# that was never surveyed. Predictions are made with depth clamped to the
# surveyed range; the unclamped surface is retained as a sensitivity.
d_rng <- range(seg_dsm$depth)
grid$depth_obs <- grid$depth
grid$depth     <- pmin(pmax(grid$depth, d_rng[1]), d_rng[2])
grid$clamped   <- grid$depth != grid$depth_obs

cat("surveyed depth range:", round(d_rng), "m | cells clamped:",
    sum(grid$clamped), sprintf("(%.1f%% of area)\n", 100 * mean(grid$clamped)))

grid$abund <- as.numeric(predict(mod, grid, grid$off.set))
grid$dens  <- grid$abund / grid$off.set          # individuals per km²

g_un <- grid; g_un$depth <- g_un$depth_obs
N_unclamped <- sum(as.numeric(predict(mod, g_un, g_un$off.set)))

N_uncorrected <- sum(grid$abund)
cat("uncorrected N:", round(N_uncorrected), "clamped |",
    round(N_unclamped), "unclamped\n")

# GAM uncertainty. dsm_varprop propagates uncertainty in the detection
# function's covariate parameters; with a covariate-free half-normal there is
# nothing for it to propagate, so the GAM and detection components are combined
# by the delta method instead.
vg <- dsm_var_gam(mod, grid, off.set = grid$off.set)
sv <- summary(vg)
cv_gam <- sv$cv

sd <- summary(ddf)$ds
cv_df <- sd$average.p.se / sd$average.p          # 0.109 from 04's summary

cv_total <- sqrt(cv_gam^2 + cv_df^2 + G0_CV^2)

N <- N_uncorrected / G0
se <- N * cv_total
C  <- exp(1.96 * sqrt(log(1 + cv_total^2)))      # lognormal CI

res <- tibble(
  quantity = c("N (uncorrected)", "N (g0-corrected)", "density /km²",
               "CV gam", "CV detection", "CV g(0)", "CV total",
               "lower 95%", "upper 95%",
               "N unclamped (sensitivity)"),
  value = c(round(N_uncorrected), round(N), signif(N / area_km2, 3),
            round(cv_gam, 3), round(cv_df, 3), G0_CV, round(cv_total, 3),
            round(N / C), round(N * C),
            round(N_unclamped / G0)))

print(as.data.frame(res))
readr::write_csv(res, here("outputs", "abundance.csv"))

cat("\nBradford et al. (2022): 2,975 (CV 0.40, 95% CI 1,407-6,291)\n")

saveRDS(list(model = mod, models = mods, comparison = gcmp, grid = grid,
             poly = poly, prj = PRJ, cell = CELL, area_km2 = area_km2,
             N = N, cv = cv_total, ci = c(N / C, N * C), g0 = G0,
             cv_parts = c(gam = cv_gam, detection = cv_df, g0 = G0_CV),
             N_unclamped = N_unclamped / G0, depth_range = d_rng,
             ddf = ddf, segments = seg_dsm, obs = obs),
        here("data", "derived", "dsm_prediction.rds"))
