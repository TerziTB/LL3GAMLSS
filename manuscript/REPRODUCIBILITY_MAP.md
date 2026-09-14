# Manuscript reproducibility map

Run commands from the repository root. `manuscript/build_manuscript_figures.R`
recreates Figs. 3–5. Fig. 6 is recreated by
`analysis/plot_era5_all_scale_maps.R` from the archived compact ERA5-Land
selection tables. Figs. 1 and 2 are author-prepared assets and are intentionally
not distributed in this interim release. The Word manuscript is maintained
manually; this submission-ready repository does not include a manuscript
generator.

| Item | Archived input | Generator |
|---|---|---|
| Table 1 | Candidate definitions in the manuscript | Maintained manually in the Word manuscript |
| Table 2 | `analysis/manuscript_inputs/station_metadata.csv` | Maintained manually in the Word manuscript |
| Table 3 | Stationary-null, scale-trend, joint-trend and interval-coverage summaries in `validation/results/` | Maintained manually in the Word manuscript |
| Table 4 | `validation/results/seasonal_selection/seasonal_selection_summary.csv` | Maintained manually in the Word manuscript |
| Table 5 | `analysis/manuscript_inputs/station_summary.csv` | Maintained manually in the Word manuscript |
| Table 6 | `analysis/manuscript_inputs/era_area.csv` and `era_basin.csv` | Maintained manually in the Word manuscript |
| Fig. 1 | Author-prepared asset to be supplied | Not included in the interim release |
| Fig. 2 | Author-prepared asset to be supplied | Not included in the interim release |
| Fig. 3 | Stress, stationary-null and interval-coverage summaries in `validation/results/` | `manuscript/build_manuscript_figures.R` |
| Fig. 4 | `analysis/manuscript_inputs/station_selection.csv` | `manuscript/build_manuscript_figures.R` |
| Fig. 5 | `analysis/manuscript_inputs/figure5_stationary_nonstationary.csv` | `manuscript/build_manuscript_figures.R` |
| Fig. 6 | Compact final ERA5-Land cell selections and basin geometry-derived plotting inputs | `analysis/plot_era5_all_scale_maps.R` |

`analysis/build_manuscript_inputs.R` rebuilds the compact manuscript inputs.
It deliberately excludes station 17934 after the quality-control decision and
does not copy precipitation, temperature, PET, or climatic-water-balance values
from the access-controlled station workbook.
