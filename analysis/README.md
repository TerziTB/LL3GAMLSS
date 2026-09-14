# Analysis workflows

The repository retains scripts and compact derived outputs for the Seyhan
station application and the independent ERA5-Land gridded application. It does
not contain restricted station observations, raw ERA5-Land downloads, fitted
object caches, or temporary processing files.

## Restricted station application

`seyhan_real_data_analysis.R` expects the provider-controlled station data in
`artifact_work/inspection/workbook_values.json`. That path is ignored by Git.
The public repository contains only compact derived summaries needed to audit
model selection and manuscript results.

Main scripts:

1. `seyhan_real_data_analysis.R`
2. `summarize_seyhan_results.R`
3. `empirical_dependence_bootstrap.R`
4. `build_manuscript_inputs.R`

The missing-input guard intentionally stops the station workflow when the
restricted input is unavailable.

## ERA5-Land application

The main input is the Copernicus ERA5-Land monthly-averaged reanalysis for
January 1950 through December 2025 over `[39.3, 34.3, 36.6, 37.0]` in
`[north, west, south, east]` order. Required NetCDF variables are `t2m` and
`tp`; optional process-consistency variables are `sd`, `swvl1`, `swvl2`, and
`swvl3`. Raw NetCDF files are not distributed.

Supply the main NetCDF as the first argument and, when available, the
2022–2024 hourly-monthly precipitation correction as the second argument:

```sh
Rscript analysis/prepare_era5_basin_data.R path/to/main.nc path/to/precipitation_correction.nc
Rscript analysis/run_era5_basin_ll3.R --full --workers=4 --refresh-cache
Rscript analysis/validate_era5_nspei.R analysis/era5_results/full_195001_202512
Rscript analysis/plot_era5_all_scale_maps.R analysis/era5_results/full_195001_202512
```

The first script creates `analysis/era5_results/era5_basin_monthly.rds`, which
is intentionally ignored. The full fit creates task caches and fitted-object
RDS files that are also ignored. Compact CSV summaries, figures, settings,
package versions, and session information are retained.

ERA5-Land soil moisture is used only for an internal process-consistency check.
It is not a predictor in the LL3 models and does not affect PET, water balance,
model selection, or nSPEI values.

## HydroBASINS boundary and redistribution

The analysis uses the HydroBASINS v1c level-8 polygon with
`MAIN_BAS = 2080001450`. The converted one-feature GeoPackage is not included
in this public repository. Obtain the European HydroBASINS level-8 data from
the [official HydroBASINS page](https://www.hydrosheds.org/products/hydrobasins),
select that feature, and save it as
`analysis/era5_input/seyhan_hydrobasins_lev08.gpkg`.

HydroBASINS is available for scientific, educational, and commercial use but
is governed by the HydroSHEDS v1 license, including attribution and
redistribution conditions. Excluding the converted subset is a conservative
repository decision, not a legal determination. Before redistributing it,
confirm that the intended distribution method complies with the current
[HydroSHEDS terms](https://www.hydrosheds.org/products/hydrobasins).

Derived maps and summaries use HydroBASINS v1 data, © WWF (2006–2022), under
the HydroSHEDS v1 license. WWF has not evaluated the derived outputs and gives
no warranty regarding their suitability. Cite:

Lehner B, Grill G (2013) Global river hydrography and network routing: baseline
data and new approaches to study the world's large river systems.
*Hydrological Processes* 27:2171–2186. https://doi.org/10.1002/hyp.9740

## Archived public outputs

- `era5_results/final_full_195001_202512/`: basin and cell model selections,
  calibrated nSPEI summaries, soil-moisture consistency results, audit files,
  figures, versions, and session information.
- `manuscript_inputs/`: compact table and figure inputs.
- `seyhan_results/`: station model-selection summaries and dependence/bootstrap
  diagnostics without the underlying meteorological observations.
