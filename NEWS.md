# neutr 0.2.0

## Etienne sampling formula

* New `optim.etienne()`: maximum likelihood estimation of the fundamental biodiversity number
  `theta` and the immigration rate `m` of Etienne's (2005) sampling formula, for a local
  community with dispersal limitation. The likelihood surface can be multimodal: it is profiled
  over both parameters before being refined with `stats::optim(method = "L-BFGS-B")` and the
  exact gradient. Returns `theta`, `m`, `I`, `logl` and `converged`, and warns when an estimate
  lies on the edge of its search interval (for `m`, the upper bound means no evidence of
  dispersal limitation).
* `optim.etienne()` returns standard errors (`se_theta`, `se_m`, `se_I`, from the inverse
  Hessian), all the distinct local maxima found (`local_maxima`, as TeTame reports a second
  maximum), and with `ci = TRUE` likelihood-ratio confidence intervals (`ci_theta`, `ci_m`,
  `ci_I`).
* New `logl.etienne()`: log-likelihood of Etienne's sampling formula, vectorised over `theta` and
  `m`. It is on the same scale as the `logl` of `optim.ewens()` and `optim.pitman()`; with
  `m = 1` it equals the Ewens log-likelihood, so the Ewens model is nested in the Etienne model.
  `full = TRUE` gives the probability of the species abundance distribution itself.
* `optim.multideme()` returns standard errors of `I` and `m` (`se_I`, `se_m`, from the exact
  curvature of the likelihood) and, with `ci = TRUE`, likelihood-ratio confidence intervals
  (`ci_I`, `ci_m`). These replace the `Std_I` and `Std_m` columns of TeTame, whose curvature
  lacks a factor `0.01 I`.
* With these additions, `neutr` covers the functionality of the TeTame 2.1 software (Chave and
  Jabot): see the vignette for reading TeTame data files and plotting the likelihood surface.
* New `logkda()`: the coefficients `K(D,A)` of the formula, computed entirely on the log scale.
  This code comes from package `untb` (R. K. S. Hankin), where it replaces earlier
  implementations that required PARI/GP or the packages `Brobdingnag`, `partitions` and
  `polynom`: `neutr` gains no new dependency. Results agree with PARI/GP to about 1e-15.

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
