#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 04 — Detection function for humpback whale, WHICEAS winter 2020
#
# Follows Bradford, Yano & Oleson (2022), NOAA-TM-NMFS-PIFSC-135, which is
# PIFSC's own analysis of this survey: 5.5 km truncation, Beaufort as the
# candidate covariate, no adjustment terms. Their mean ESW was 3.72 km, so
# that figure doubles as a check on this fit.
#
# Produces:
#   data/derived/df_humpback.rds     fitted detection function + the data it saw
#   outputs/detection_models.csv     model comparison table
#   outputs/figures/df_diagnostic.png
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(here); library(Distance)

dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

SPP   <- "076"    # Megaptera novaeangliae, asserted in 03
TRUNC <- 5.5      # km, Bradford et al. (2022)

cruz10 <- readRDS(here("data", "derived", "cruz_10km.rds"))

# ── effort ───────────────────────────────────────────────────────────────────
# Systematic on-effort segments inside the WHICEAS study area, winter 2020.
# `use` is LTabundR's own flag for effort it considers analysable.
seg0 <- cruz10$cohorts$all$segments
stopifnot(is.logical(seg0$OnEffort), is.logical(seg0$use))

seg <- seg0 |>
  filter(year == 2020, EffType == "S", OnEffort, use, stratum == "WHICEAS")

cat("segments:", nrow(seg),
    "| effort", round(sum(seg$dist)), "km",
    "| Beaufort", round(mean(seg$avgBft, na.rm = TRUE), 2),
    "mean,", sum(is.na(seg$avgBft)), "missing\n")
stopifnot(nrow(seg) > 300, sum(is.na(seg$avgBft)) == 0)

# ── detections ───────────────────────────────────────────────────────────────
# Keyed to the segments above by seg_id rather than re-filtered on stratum, so
# the detections cannot disagree with the effort they will be offset against.
sit <- cruz10$cohorts$all$sightings |>
  filter(species == SPP, included, seg_id %in% seg$seg_id)

cat("humpback detections on those segments:", nrow(sit), "\n")
stopifnot(nrow(sit) > 30)

# Distance and Beaufort must both be present; ds() would drop rows silently.
miss <- sit |> summarise(no_dist = sum(is.na(PerpDistKm)),
                         no_bft  = sum(is.na(Bft)),
                         no_size = sum(is.na(best)))
print(miss)

sit <- sit |> filter(!is.na(PerpDistKm), !is.na(Bft))

beyond <- sum(sit$PerpDistKm > TRUNC)
cat("beyond", TRUNC, "km truncation:", beyond,
    sprintf("(%.0f%%)\n", 100 * beyond / nrow(sit)))

df_dat <- sit |>
  transmute(object = row_number(),
            distance = PerpDistKm,
            Bft,
            size = best)

# Segment linkage for the DSM, kept out of the data ds() sees: a partial set of
# Region.Label / Area / Sample.Label / Effort columns changes how ds() behaves.
obs <- sit |>
  transmute(object = row_number(),
            Sample.Label = seg_id,
            size = best,
            distance = PerpDistKm,
            lon = Lon, lat = Lat)
stopifnot(identical(obs$object, df_dat$object))

cat("fitting on", sum(df_dat$distance <= TRUNC), "detections",
    "| mean group size", round(mean(df_dat$size[df_dat$distance <= TRUNC],
                                    na.rm = TRUE), 2), "\n")

# ── candidate models ─────────────────────────────────────────────────────────
# Half-normal and hazard-rate, each with and without Beaufort. Adjustment terms
# are switched off throughout: they are not permitted alongside covariates, and
# allowing them for the null models only would make the AIC comparison unfair.
fit <- function(key, formula) {
  ds(df_dat, truncation = TRUNC, key = key, formula = formula,
     adjustment = NULL, quiet = TRUE)
}

m <- list(
  `hn ~1`   = fit("hn", ~1),
  `hn ~Bft` = fit("hn", ~Bft),
  `hr ~1`   = fit("hr", ~1),
  `hr ~Bft` = fit("hr", ~Bft)
)

# Mean effective strip width, averaged over the covariate values actually seen.
esw <- vapply(m, function(x) mean(predict(x, esw = TRUE)$fitted), numeric(1))

cmp <- tibble(model = names(m),
              aic   = vapply(m, function(x) AIC(x)$AIC, numeric(1)),
              esw   = round(esw, 3),
              p     = round(esw / TRUNC, 3),
              cvm_p = vapply(m, function(x) gof_ds(x, plot = FALSE)$dsgof$CvM$p,
                             numeric(1))) |>
  mutate(delta_aic = round(aic - min(aic), 2), aic = round(aic, 2)) |>
  arrange(aic)

print(as.data.frame(cmp))
readr::write_csv(cmp, here("outputs", "detection_models.csv"))

cat("\nBradford et al. (2022) mean ESW: 3.72 km\n")

# ── keep the best by AIC ─────────────────────────────────────────────────────
best_name <- cmp$model[1]
best <- m[[best_name]]
cat("selected:", best_name, "\n")
print(summary(best))

png(here("outputs", "figures", "df_diagnostic.png"), 1200, 500, res = 120)
par(mfrow = c(1, 2))
plot(best, main = best_name, xlab = "Perpendicular distance (km)")
plot(df_dat$Bft[df_dat$distance <= TRUNC],
     df_dat$distance[df_dat$distance <= TRUNC],
     xlab = "Beaufort", ylab = "Distance (km)", pch = 16, col = "#00000060")
dev.off()

saveRDS(list(model = best, all_models = m, comparison = cmp,
             data = df_dat, obs = obs, truncation = TRUNC, segments = seg),
        here("data", "derived", "df_humpback.rds"))
