## neutr R package
### Parameter inference and fast generation for three classes of neutral species assemblage models.
### Jerome Chave
jerome.chave@cnrs.fr

## History

- **First release**: May 2024
- **Second release**: October 2026
  - Fixes:
    1. Package without `usethis` dependency
    2. Fixed incorrect stick-breaking equation
    3. Fixed random number generator bug for large values
    4. Fixed comparative maximum likelihood value
    5. Theta no longer restricted to be ≥ 1
    6. Now Ewen’s and Pitman likelihoods are comparable
    7. Package no longer calls `nloptr` and `pracma`
    8. Cleaner management of boundary values


## Description

The neutr R package contains a number of functions that perform the following tasks.
* Estimation of the model parameter for three neutral models: the Ewens model, the multideme model and the Pitman models. 
* Parameter inference is based on the maximization of the likelihood functions. 
* The generation of typical abundance distributions of the Ewens and Pitman models given model parameter(s). 

The above tasks can be performed with very large sample sizes, on the order of up to 10^12 individuals. 

## Installation

```r
# install.packages("remotes")
remotes::install_github("jeromechave/neutr")
```

The package needs only `dqrng` (fast random numbers); everything else is base R.

## Quick start

```r
library(neutr)

# Generate a community: Ewens model, theta = 11.3, 234,373 individuals
set.seed(1)
com <- generate.hoppe.urn(theta = 11.3, J = 234373)
length(com$abundance)      # number of species in this draw (expected: com$k)

# Estimate the parameters
ab <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
ew <- optim.ewens(ab)                        # Ewens theta
pm <- optim.pitman(ab, c(10, 0.1))           # Pitman theta and sigma

# Compare the nested models (Ewens is Pitman with sigma = 0): the log-likelihoods are comparable
2 * (pm$logl - ew$logl)                      # likelihood-ratio statistic

# Multideme model: one row per deme, one column per species
optim.multideme(rbind(c(44, 37, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1),
                      c(240, 20, 48, 2, 21, 1, 3, 2, 5, 2, 0, 1, 1)))
```

See `vignette("Using-the-neutr-package")` for a longer tour.

## Functions

| Task | Function |
|---|---|
| ML estimate of Ewens' $\theta$ | `optim.ewens()` |
| ML estimate of Pitman's $(\theta, \sigma)$ | `optim.pitman()` |
| ML estimate of the local immigration rates (K demes) | `optim.multideme()` |
| Exact Hoppe urn, individual by individual | `generate.hoppe.urn0()` |
| Ewens community via the GEM representation (any sample size) | `generate.hoppe.urn()` |
| Pitman community via the GEM representation (any sample size) | `generate.pitman.urn()` |
| Expected number of species | `kest()` |
| Expected number of species above an abundance threshold | `kest.gt()` (integral), `kest.gt0()` (simulation) |
| Drop zeros and `NA`s | `zeroremove()` |

## Good to know

* `generate.hoppe.urn0()` uses the `dqrng` random stream, which `set.seed()` does not
  control: pass its `seed` argument for reproducibility. The other generators follow `set.seed()`.
* The `logl` returned by `optim.ewens()` and `optim.pitman()` is on the same scale, so the two
  (nested) models can be compared by likelihood ratio or AIC.
* Estimates that fall on the edge of their search interval trigger a warning: the likelihood then has
  no interior maximum for those data (e.g. all species are singletons).
* `optim.multideme()` expects **rows = demes, columns = species**.

See [NEWS.md](NEWS.md) for the list of changes.
