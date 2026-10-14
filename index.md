## neutr R package

### Parameter inference and fast generation for three classes of neutral species assemblage models.

### Jerome Chave

<jerome.chave@cnrs.fr>

Documentation website: <https://jeromechave.github.io/neutr/>

## History

- **First release**: May 2024
- **Second release**: October 2026
  - Fixes:
    1.  Package without `usethis` dependency
    2.  Fixed incorrect stick-breaking equation
    3.  Fixed random number generator bug for large values
    4.  Fixed comparative maximum likelihood value
    5.  Theta no longer restricted to be ≥ 1
    6.  Now Ewen’s and Pitman likelihoods are comparable
    7.  Package no longer calls `nloptr` and `pracma`
    8.  Cleaner management of boundary values
    9.  Now includes Etienne Sampling Formula (from untb package, with
        improvements)

## Description

The neutr R package contains a number of functions that perform the
following tasks. \* Estimation of the model parameter for four neutral
models: the Ewens model, the Etienne model (a local community with
dispersal limitation), the multideme model and the Pitman models. \*
Parameter inference is based on the maximization of the likelihood
functions. \* The generation of typical abundance distributions of the
Ewens and Pitman models given model parameter(s).

The above tasks can be performed with very large sample sizes, on the
order of up to 10^12 individuals (up to about 10^5 individuals for the
Etienne model).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("jeromechave/neutr")
```

The package needs only `dqrng` (fast random numbers); everything else is
base R.

## Quick start

``` r

library(neutr)

# Generate a community: Ewens model, theta = 11.3, 234,373 individuals
set.seed(1)
com <- generate_hoppe_urn(theta = 11.3, J = 234373)
length(com$abundance)      # number of species in this draw (expected: com$k)

# Estimate the parameters
ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
ew <- optim_ewens(ab)                        # Ewens theta
pm <- optim_pitman(ab, c(10, 0.1))           # Pitman theta and sigma

# Compare the nested models (Ewens is Pitman with sigma = 0): the log-likelihoods are comparable
2 * (pm$logl - ew$logl)                      # likelihood-ratio statistic

# Etienne model: a local community with dispersal limitation (Ewens is Etienne with m = 1)
et <- optim_etienne(c(8, 5, 3, 2, 1, 1))     # theta = 7.05, m = 0.226 (Etienne 2005)

# Multideme model: one row per deme, one column per species
optim_multideme(rbind(c(44, 37, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1),
                      c(240, 20, 48, 2, 21, 1, 3, 2, 5, 2, 0, 1, 1)))
```

See
[`vignette("Using-the-neutr-package")`](https://jeromechave.github.io/neutr/articles/Using-the-neutr-package.md)
for a longer tour.

## Functions

| Task | Function |
|----|----|
| ML estimate of Ewens’ $`\theta`$ | [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md) |
| ML estimate of Pitman’s $`(\theta, \sigma)`$ | [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md) |
| ML estimate of Etienne’s $`(\theta, m)`$ (local community with dispersal limitation) | [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md) |
| Log-likelihood of Etienne’s sampling formula | [`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md) |
| Etienne’s coefficients $`K(D,A)`$, on the log scale | [`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md) |
| ML estimate of the local immigration rates (K demes) | [`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md) |
| Exact Hoppe urn, individual by individual | [`generate_hoppe_urn0()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn0.md) |
| Ewens community via the GEM representation (any sample size) | [`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md) |
| Pitman community via the GEM representation (any sample size) | [`generate_pitman_urn()`](https://jeromechave.github.io/neutr/reference/generate_pitman_urn.md) |
| Expected number of species | [`kest()`](https://jeromechave.github.io/neutr/reference/kest.md) |
| Expected number of species above an abundance threshold | [`kest_gt()`](https://jeromechave.github.io/neutr/reference/kest_gt.md) (integral), [`kest_gt0()`](https://jeromechave.github.io/neutr/reference/kest_gt0.md) (simulation) |
| Drop zeros and `NA`s | [`remove_zero()`](https://jeromechave.github.io/neutr/reference/remove_zero.md) |

## Good to know

- [`generate_hoppe_urn0()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn0.md)
  uses the `dqrng` random stream, which
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not control:
  pass its `seed` argument for reproducibility. The other generators
  follow [`set.seed()`](https://rdrr.io/r/base/Random.html).
- The `logl` returned by
  [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md),
  [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)
  and
  [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
  is on the same scale. The Ewens model is nested in both the Pitman
  model (`sigma = 0`) and the Etienne model (`m = 1`), so the models can
  be compared by likelihood ratio or AIC.
- [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
  and
  [`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)
  return standard errors, and likelihood-ratio confidence intervals with
  `ci = TRUE`;
  [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
  also lists all the local maxima it finds. Together with the vignette
  section on TeTame data files, this covers the functionality of the
  TeTame software.
- For repeated evaluations of
  [`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md)
  on the same data, compute
  [`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md)
  once and pass it as `log_kda`.
- Estimates that fall on the edge of their search interval trigger a
  warning: the likelihood then has no interior maximum for those data
  (e.g. all species are singletons).
- [`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)
  expects **rows = demes, columns = species**.

See [NEWS.md](https://jeromechave.github.io/neutr/NEWS.md) for the list
of changes.
