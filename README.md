# LL3GAMLSS

`LL3GAMLSS` is an R package implementing a support-safe three-parameter
shifted log-logistic distribution for classic `gamlss`. Its validated
nonstationary formulation estimates a constant threshold and permits
covariates in the scale and/or shape parameters.


## Repository guide

- [Installation](#installation)
- [Contribution and validated scope](#contribution-and-validated-scope)
- [Validation design, commands, and archived outputs](validation/README.md)
- [Station and ERA5-Land analyses](analysis/README.md)
- [Manuscript reproducibility map](manuscript/REPRODUCIBILITY_MAP.md)
- [Machine-readable validation summary](VALIDATION_SUMMARY.json)
- [Citation](#citation)
- [License](#license)

## Contribution and validated scope

This project does **not** claim to introduce the first nonstationary SPEI.
Nonstationary log-logistic SPEI formulations, including GAMLSS models with
covariate-dependent location parameters, have already been published.

The intended contribution is narrower and reproducible:

- an installable classic-`gamlss` family for the shifted three-parameter
  log-logistic distribution, which is not supplied by `gamlss.dist`;
- a constant-threshold formulation with nonstationary scale and shape;
- explicit support-boundary safeguards and an inference-readiness gate;
- validation against `FAdist`, numerical derivatives, and a separate direct
  optimizer;
- an operational calendar-month, multi-timescale LL3-SPEI workflow; and
- iid parametric and moving-block residual bootstrap procedures.

“Three-parameter” means that `mu`, `sigma`, and `nu` are estimated. In the
validated nonstationary scope, `mu` is constant in time. Covariate-dependent
`mu` remains experimental because a changing threshold changes conditional
support and requires separate theoretical and simulation validation.

## Installation

```r
install.packages(c("gamlss", "gamlss.dist", "FAdist", "numDeriv"))
install.packages(".", repos = NULL, type = "source")
library(LL3GAMLSS)
```

For development checks, install `testthat` as well.

The package is also installable from a source archive attached to a future
GitHub Release. Source archives are intentionally ignored in normal Git
history.

## Minimal reproducible example

```r
library(LL3GAMLSS)

set.seed(2026)
n <- 66L * 12L
climate <- data.frame(
  date = seq(as.Date("1960-01-01"), by = "month", length.out = n),
  time_scaled = as.numeric(scale(seq_len(n)))
)
climate$water_balance <- rLL3(
  n,
  mu = -50,
  sigma = exp(log(70) + 0.20 * climate$time_scaled),
  nu = 2.5
)

fit <- fit_LL3_spei(
  climate,
  response = "water_balance",
  date = "date",
  scale = 3,
  sigma.formula = ~ time_scaled,
  nu.formula = ~ 1
)

head(fit$results[c("date", "accumulated_balance", "LL3_SPEI")])
```

## Distribution functions

```r
p <- c(0.01, 0.5, 0.99)
x <- qLL3(p, mu = -50, sigma = 70, nu = 1.8)
pLL3(x, mu = -50, sigma = 70, nu = 1.8)

LL3_to_FAdist(mu = -50, sigma = 70, nu = 1.8)
# shape = 1 / nu, scale = log(sigma), thres = mu
```

## Fitting the validated candidate set

```r
models <- fit_LL3_candidate_models(
  data = dat,
  response = "balance",
  covariate = "time_scaled"
)

LL3_model_table(
  models[c("stationary", "sigma", "nu", "joint")],
  n = nrow(dat)
)
```

The four models are stationary; covariate-dependent scale; covariate-dependent
shape; and joint scale-shape nonstationarity, always with a constant threshold.

Minimum AIC across this four-model set has a substantial family-wise
false-selection rate in the supplied stationary-null simulations, especially
for short records. Pre-specify scientifically justified predictors, report BIC
and out-of-sample calibration, and do not describe minimum-AIC selection as a
hypothesis test.

## Operational LL3-SPEI

Canonical monthly SPEI accumulates climatic water balance at a selected scale,
fits a separate distribution for each calendar month, and transforms the
conditional CDF to a standard normal variate. `fit_LL3_spei()` implements that
workflow for LL3 models:

```r
spei12 <- fit_LL3_spei(
  data = climate,
  response = "water_balance",
  date = "date",
  scale = 12,
  sigma.formula = ~ time_scaled,
  nu.formula = ~ 1,
  reference_start = as.Date("1961-01-01"),
  reference_end = as.Date("2020-12-01")
)

head(spei12$results)
plot(spei12$results$date, spei12$index, type = "l")
```

Dates must be unique and increasing. The current operational prediction path
supports parametric formula terms. The underlying family can fit GAMLSS
smoothers, but smooth-term out-of-sample prediction needs a separately
validated prediction method before it is used in this workflow.

The original SPEI uses the three-parameter log-logistic distribution by
default, but other distributions are admissible and implemented by the
official `SPEI` package. An application paper should compare LL3 with credible
alternatives rather than treating LL3 as a mathematical requirement.

## Boundary policy

The threshold link is

\[
\eta_\mu=\log(lower-\mu), \qquad \mu=lower-\exp(\eta_\mu),
\]

where `lower` lies just below the minimum fitting response. Only contact within
floating-point tolerance is projected. Material violations raise an error.

```r
diagnostic <- LL3_boundary_diagnostic(fit, y)
LL3_assert_inference_ready(fit, y)
```

Boundary contact is not a routine successful solution. It indicates
non-regular threshold estimation; ordinary Hessian, Wald, AIC, and BIC
interpretations may not apply. The operational SPEI and uncertainty helpers
therefore reject boundary-contact fits by default. Report a pre-specified
threshold sensitivity analysis whenever fitted thresholds approach the
boundary.

## Independent likelihood comparison

```r
comparison <- compare_LL3_gamlss_direct(
  fit,
  data = dat,
  mu.formula = balance ~ 1,
  sigma.formula = ~ time_scaled,
  nu.formula = ~ time_scaled,
  n_starts = 20
)
```

Agreement validates optimization. It is not itself an uncertainty method.

## Uncertainty and dependence

The OPG curvature is a stable working curvature for GAMLSS optimization, not
the reported covariance matrix. For iid conditional observations:

```r
boot_iid <- LL3_parametric_bootstrap(
  fit,
  data = dat,
  response = "balance",
  mu.formula = balance ~ 1,
  sigma.formula = ~ time_scaled,
  nu.formula = ~ time_scaled,
  B = 999,
  seed = 2026
)
```

If normalized residuals retain temporal dependence:

```r
boot_block <- LL3_moving_block_bootstrap(
  fit,
  data = dat,
  response = "balance",
  mu.formula = balance ~ 1,
  sigma.formula = ~ time_scaled,
  nu.formula = ~ time_scaled,
  B = 999,
  block_length = 12,
  seed = 2026
)
```

The moving-block procedure resamples circular blocks of normalized residuals
and maps them through the fitted conditional LL3 quantile function.

## Validation and package checks

The complete validation inventory, execution order, archived outputs, and
interpretation limits are documented in [validation/README.md](validation/README.md).

The research-validation scripts below are provided in the full repository,
not in the installable package tarball. From a full repository checkout, run:

```r
source(file.path("validation", "validate_probability_functions.R"))
source(file.path("validation", "validate_derivatives.R"))
source(file.path("validation", "validate_stationary_MLE.R"))
```

Package tests and checks can be run from a full repository checkout or an
unpacked source-package directory:

```sh
Rscript -e "testthat::test_local('.')"
R CMD build .
R CMD check --no-manual LL3GAMLSS_0.1.1.tar.gz
```

In the full repository, the scripts in `validation/` cover stationary-null
selection, scale and joint
nonstationarity, interval coverage, optimizer discrepancies, lower-bound
sensitivity, short-record/dependence/contamination stress scenarios, and the
full seasonal selection workflow used by the empirical analysis. Run
`simulate_seasonal_model_selection.R` to evaluate all twelve calendar-month
fits at 1-, 3-, 6-, and 12-month scales under iid and AR(1) scenarios.
The default 30 replicates per design cell are a compact workflow validation,
not a high-precision power study. Longer-scale sums are not exactly LL3, so
the archived summary reports both operational selection completion and the
stricter rate at which all four candidates remain inference-ready.

`validate_spei_maxlik_equivalence.R` also exercises the maximum-likelihood
generalized-logistic engine supplied with SPEI 1.8.1. The script records an
installed-version limitation: the public `spei(..., fit = "max-lik")` path
returns missing coefficients because its starting-value switch omits that
method. The benchmark therefore calls SPEI's own `parglo.maxlik` engine with
the package's unbiased-PWM start. It flags fitted generalized-logistic shapes
with non-negative kappa because those have an upper finite endpoint and are
not parameter-equivalent to the lower-threshold LL3 family.

## Reproducibility

Seeds, requested and completed replicate counts, Monte Carlo standard errors,
fit failures, boundary contacts, settings, runtimes, and session information
are written beside each validation result. In the full repository (not the
installable package tarball), the `analysis/manuscript_inputs/` directory
contains compact derived inputs needed to rebuild manuscript tables and figures
without redistributing the restricted station observations.

The empirical and gridded workflows, required external inputs, and output
provenance are documented in [analysis/README.md](analysis/README.md). The
[manuscript reproducibility map](manuscript/REPRODUCIBILITY_MAP.md) connects
each retained table and figure to its archived inputs and generator.

## Citation

Use `citation("LL3GAMLSS")` after installation. Citation metadata are also
provided in [CITATION.cff](CITATION.cff) and [inst/CITATION](inst/CITATION).
Repository and archive identifiers should be added only after they exist.

## License

LL3GAMLSS is distributed under the GNU General Public License version 3. See
[LICENSE](LICENSE). Third-party datasets remain subject to their own terms;
see [analysis/README.md](analysis/README.md) before redistributing analysis
inputs.

## Key references

- Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
  for location, scale and shape. *Applied Statistics*, 54, 507–554.
- Vicente-Serrano, S. M., Beguería, S. and López-Moreno, J. I. (2010). A
  multiscalar drought index sensitive to global warming: the SPEI. *Journal of
  Climate*, 23, 1696–1718.
- Bazrafshan, J., Cheraghalizadeh, M. and Shahgholian, K. (2022). Development
  of a non-stationary SPEI for drought monitoring in a changing climate.
  *Water Resources Management*, 36, 3523–3543.
- Masanta, S. K. and Srinivas, V. V. (2022). Proposal and evaluation of
  nonstationary versions of SPEI and SDDI based on climate covariates.
  *Journal of Hydrology*, 610, 127808.
