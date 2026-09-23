#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 03 — Reprocess PIFSC WinCruz DAS data into short segments for density
#      surface modelling
#
# The public WinCruz sighting archive used in 01/02 carries no effort and no
# perpendicular distances, so it cannot support a DSM. The LTabundR vignette
# distributes the underlying DAS files and PIFSC's own processing settings,
# which can.
#
# Produces:
#   data/derived/cruz_10km.rds       processed surveys, ~10 km segments
#   outputs/species_by_year.csv      detections by year and species
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(here); library(LTabundR)

dir.create(here("data", "derived"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs"),         recursive = TRUE, showWarnings = FALSE)

VIGNETTE <- here("data", "raw", "LTabundR-vignette")

# ── source data ──────────────────────────────────────────────────────────────
# DAS files and PIFSC's processing settings live in the vignette repository.
if (!dir.exists(VIGNETTE)) {
  system2("git", c("clone",
                   "https://github.com/PIFSC-Protected-Species-Division/LTabundR-vignette.git",
                   shQuote(VIGNETTE)))
}

# cruz_1720: PIFSC's processed object for HICEAS 2017 + WHICEAS 2020.
# We use it for two things: its settings, and as the reference to check
# our reprocessing against.
nm <- load(file.path(VIGNETTE, "whiceas_cruz_1720.RData"))
stopifnot("cruz_1720" %in% nm)

seg_ref <- cruz_1720$cohorts$all$segments
sit_ref <- cruz_1720$cohorts$all$sightings

# ── what surveys are in here ─────────────────────────────────────────────────
# 2017 months 7-11 is HICEAS; 2020 months 1-3 is WHICEAS, the winter survey.
# Note "WHICEAS" in the `stratum` column is a study-area name, not a survey.
print(table(seg_ref$year, seg_ref$month))

# ── why the supplied segments will not do ────────────────────────────────────
# LTabundR segments effort into ~150 km blocks because, for design-based
# abundance, segments are the bootstrap resampling unit. In a DSM segments are
# the rows of the GAM, and a spatial smooth learns nothing from 106 rows
# averaging 100 km. The other prepared objects in the vignette
# (whiceas_cruz.RData, whiceas_cruz_segment.RData) are segmented the same way,
# so there is no shorter-segmented version to borrow.
cat("supplied segments:", nrow(seg_ref),
    "| median length", round(median(seg_ref$dist), 1), "km",
    "| total effort", round(sum(seg_ref$dist)), "km\n")

# ── which species can carry a detection function ─────────────────────────────
# ~60-80 detections is the conventional minimum. Only humpback whales in the
# 2020 winter survey come close; nothing in 2017 is within a factor of seven.
species_by_year <- sit_ref |>
  filter(OnEffort, EffType == "S", included) |>
  count(year, species, name = "detections") |>
  left_join(cruz_1720$settings$survey$species_codes |>
              select(species = code, common),
            by = "species") |>
  arrange(desc(detections))

print(head(as_tibble(species_by_year), 15))
readr::write_csv(species_by_year, here("outputs", "species_by_year.csv"))

# ── reprocess with short segments ────────────────────────────────────────────
# PIFSC's own settings, with one parameter changed, so the difference from the
# published object is exactly one line and anyone can see what it is.
OUT <- here("data", "derived", "cruz_10km.rds")

if (!file.exists(OUT)) {
  set <- cruz_1720$settings
  set$survey$segment_target_km <- 10        # was 150

  das <- file.path(VIGNETTE, "data", "surveys",
                   "CenPac1986-2020_Final_alb_edited.das")
  stopifnot(file.exists(das))

  # The edited file is the one cruz_1720 was built from (see `file_das` in its
  # sightings table). Covers 1986-2020, so this takes a few minutes.
  cruz10 <- process_surveys(das, settings = set)
  saveRDS(cruz10, OUT)
} else {
  cruz10 <- readRDS(OUT)
}

# ── check the reprocessing ───────────────────────────────────────────────────
s10 <- cruz10$cohorts$all$segments |> filter(year == 2020)

cat("2020 segments:", nrow(s10),
    "| median length", round(median(s10$dist), 2), "km\n")
stopifnot(nrow(s10) > 400, median(s10$dist) > 8, median(s10$dist) < 12)

# Humpback is species code 076 in LTabundR's own table; the whole analysis
# hangs on this, so it is asserted rather than assumed.
hump <- cruz10$settings$survey$species_codes |> filter(code == "076")
print(hump)
stopifnot(grepl("Megaptera", hump$scientific_name))
