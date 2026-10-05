# neutr 0.1.1

## Estimation

* `optim.ewens()` and `optim.pitman()` return the log-likelihood (`logl`) on the same scale
  (the log-probability of the partition), so the Ewens model, which is the Pitman model with
  `sigma = 0`, can be compared directly with the Pitman model by likelihood ratio or AIC.
* `optim.ewens()` solves the score equation `k = theta (digamma(theta + J) - digamma(theta))`
  directly, to full numerical precision.
* `optim.pitman()` maximises the likelihood with `stats::optim(method = "L-BFGS-B")` and the exact
  gradient, from two starting points. New arguments `lower` and `upper` set the search interval,
  and the result has a `converged` element.
* `optim.multideme()` scans each deme likelihood on a log grid before refining the maximum.
  New `logl` element in the result, new `verbose` and `I_max` arguments; a data frame of counts is
  accepted. A deme with a single individual has `NA` immigration rate (with a warning).
* The estimators warn when an estimate lies on the edge of its search interval, i.e. when the
  likelihood has no interior maximum for the data.

## Generation

* `generate.hoppe.urn0()` is vectorised and much faster. New `seed` argument (the `dqrng` generator
  is not controlled by `set.seed()`).
* `generate.hoppe.urn()` and `generate.pitman.urn()` are exact for any sample size up to 2^53,
  and the total number of individuals is always exactly `J`.
* `kest.gt0()` gains an `nrep` argument (default 100).

## General

* Arguments are validated, with error messages that name the calling function.
* New vignette, expanded README, and a `testthat` test suite.
* `neutr` now imports only `dqrng` and `stats`.
