# Cetaceans in the Hawaiian Islands EEZ

Two analyses of cetacean survey data from the NOAA Pacific Islands Fisheries
Science Center (PIFSC), written in R and wrapped in Quarto reports.

**[Read the reports →](https://theophilemt92.github.io/hawaii-cetacean-seasonality/)**

![](outputs/figures/fig_panel.png)

---

## 1. Seasonal cetacean sighting records

An animated map of where cetaceans were recorded across the Hawaiian Islands
EEZ, month by month, over bathymetry.

370 sightings of 24 identified species, from nine survey cruises between 2009
and 2017, inside the 2,474,715 km² EEZ. The report is as much about what the
data cannot say as what it can: the public archive records sightings without
recording effort, so it can describe where animals were seen but not how many
there are, and apparent seasonal pattern is partly the pattern of when ships
were at sea.

Source: [NOAA InPort item 18141](https://www.fisheries.noaa.gov/inport/item/18141).

## 2. A density surface model for humpback whales

What the sighting archive cannot support, line-transect data can. This report
reprocesses PIFSC's own WinCruz DAS files into 10 km segments, fits a detection
function and a spatial model, and estimates abundance for the winter 2020
survey (WHICEAS).

448 systematic on-effort segments, 4,440 km of trackline, 72 humpback whale
detections within a 5.5 km truncation. A half-normal detection function and a
Tweedie GAM over depth and position, corrected for incomplete trackline
detection with g(0) = 0.68.

**Estimate: 2,790 animals (95% CI 1,193–6,529)**, against the published
design-based figure of 2,975 (CI 1,407–6,291) for the same survey — 6% apart,
from a different modelling framework.

Source: DAS files and processing settings from the
[LTabundR vignette repository](https://github.com/PIFSC-Protected-Species-Division/LTabundR-vignette).

---

## Repository

```
R/01_prepare_data.R        download and filter the public sighting archive
R/02_seasonal_map.R        twelve-frame animated map, MP4
R/03_process_das.R         DAS → 10 km segments with LTabundR
R/04_detection_function.R  detection function, model selection
R/05_dsm.R                 density surface model, prediction, variance
R/06_figures.R             DSM figures

index.qmd                  landing page
cetacean-seasonal-map.qmd  report 1
humpback-dsm.qmd           report 2
outputs/                   model comparison tables, abundance, figures
```

Scripts run in order within each analysis; `01–02` and `03–06` are independent
of each other. `03` clones the vignette repository and reprocesses the full
1986–2020 DAS file, which takes a few minutes. Everything downstream is
seconds.

`data/` is not tracked — both analyses download or clone what they need.

### Packages

`LTabundR`, `Distance`, `dsm`, `mgcv`, `sf`, `terra`, `marmap`, `ggplot2`,
`patchwork`, `rnaturalearth`, `rphylopic`, `av`.

## Data and attribution

Survey data collected by NOAA Pacific Islands Fisheries Science Center,
Protected Species Division, under the Hawaiian Islands Cetacean and Ecosystem
Assessment Survey programme, and distributed publicly. Processing for the
second analysis uses
[LTabundR](https://github.com/PIFSC-Protected-Species-Division/LTabundR), the
Division's own R package.

Reference estimate: Bradford, A.L., Yano, K.M. & Oleson, E.M. (2022).
*Abundance estimates of cetaceans from a 2020 survey of the Hawaiian Islands
EEZ*. NOAA Technical Memorandum NMFS-PIFSC-135. Trackline detection
probability from Barlow, J. (2015), *Marine Mammal Science* 31: 923–943.

Bathymetry from ETOPO via `marmap`.

## Author

Théophile L. Mouton — quantitative marine ecologist.
