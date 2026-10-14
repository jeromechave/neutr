# Maximum likelihood estimation of \\\theta\\ and \\m\\ based on Etienne sampling formula

Etienne sampling formula is an exact solution of Hubbell's (2001) model
where a single deme is connected to a single regional pool with
dispersal limitation. In Hubbell's (2001) model, a dead individual is
immediately replaced by offspring of other local individuals or by
immigrants from the regional pool with probability \\m\\. The diversity
of the regional pool is described by \\\theta\\. Etienne sampling
formula infers \\m\\ and \\\theta\\ based on one sample of the deme. If
there is no dispersal limitation, the sample can be considered as
representative of the regional pool, corresponding to the limiting case
\\m=1\\, where the Etienne sampling formula becomes equivalent to the
Ewens sampling formula.

`optim_etienne()` computes the maximum likelihood estimates of the
fundamental biodiversity number \\\theta\\ and the immigration rate
\\m\\ of Etienne's sampling formula (Etienne 2005) for a species
abundance distribution \\(n_1, ..., n_k)\\ sampled from a single local
community that receives immigrants from a neutral metacommunity. The
input should be a vector and initial values of \\(\theta, m)\\ can be
supplied with the argument `init_vals = c(theta_init, m_init)` with the
condition `theta_init > 0` and `0 < m_init < 1`.

The coefficients \\K(D,A)\\ are computed once with
[`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md).
The log-likelihood is then maximised with the box-constrained
quasi-Newton method L-BFGS-B
([`stats::optim()`](https://rdrr.io/r/stats/optim.html)) using the exact
gradient, on the scale \\(\log\theta, \log I)\\ with
\\I=\frac{m}{1-m}(J-1)\\. The optimisation is The likelihood surface can
be multimodal, with narrow peaks next to ridges that run to the edge of
the search domain, so it is first profiled in both directions: on a grid
of 61 values of \\\log I\\ the likelihood is maximised over \\\theta\\,
and on a grid of 61 values of \\\log\theta\\ it is maximised over \\I\\.
The optimisation is then run from every local maximum of the two
profiles, from `init_vals`, and from the Ewens estimate with \\m\\ close
to 1, and the best result is kept.

## Usage

``` r
optim_etienne(
  input_abundances,
  init_vals = c(10, 0.5),
  lower = c(1e-06, 1e-08),
  upper = NULL,
  ci = FALSE,
  level = 0.95
)
```

## Arguments

- input_abundances:

  a vector of integers (species abundances); species with zero abundance
  are ignored

- init_vals:

  initial values of (theta, m), a vector of two real numbers

- lower, upper:

  bounds of the search interval for `c(theta, m)`

- ci:

  if `TRUE`, also return likelihood-ratio confidence intervals (default
  `FALSE`)

- level:

  confidence level of the intervals (default 0.95)

## Value

A list with

- theta:

  fundamental biodiversity number \\\theta\\, a real value

- m:

  immigration rate, a real value

- I:

  rescaled immigration rate, \\I=\frac{m}{1-m}(J-1)\\, a real value

- logl:

  maximal log-likelihood, a real value; it is comparable with the `logl`
  value of
  [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)
  and
  [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)

- converged:

  `TRUE` if the optimiser reported convergence, or stopped on a
  line-search failure at a point where the (projected) gradient vanishes

- se_theta, se_m, se_I:

  standard errors of `theta`, `m` and `I`, real values

- local_maxima:

  a data frame of the distinct local maxima, best first, with columns
  `theta`, `m`, `I`, `logl` and `boundary`

- ci_theta, ci_m, ci_I:

  only if `ci = TRUE`: confidence intervals, vectors `c(lower, upper)`

## Search interval and boundary warnings

By default \\\theta\\ is searched in `[1e-6, 100 J]` and \\m\\ in
`[1e-8, 1 - 1e-8]`, where `J` is the number of individuals. A warning is
issued if an estimate lies on the edge of its interval; use `lower` and
`upper` to change it. An estimate of \\m\\ at its upper bound means that
the data show no evidence of dispersal limitation: the Ewens model
([`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md))
is then adequate.

## Uncertainty

The standard errors `se_theta`, `se_m` and `se_I` are obtained from the
inverse of the Hessian of the log-likelihood (on the scale
\\(\log\theta, \log I)\\, then the delta method). They account for the
correlation between the two parameters, unlike the `Std_` columns of the
TeTame software, and are `NA` when an estimate lies on a bound of the
search interval. As the likelihood is often strongly asymmetric, the
likelihood-ratio confidence intervals returned with `ci = TRUE` are
preferable: for each parameter, the range over which the profile
log-likelihood (maximised over the other parameter) stays within
`qchisq(level, 1) / 2` of its maximum. A limit equal to a bound of the
search interval means that the interval extends at least to that bound.
When the likelihood runs along a ridge towards \\m = 1\\, \\\theta\\ and
\\m\\ are weakly identified, the chi-squared approximation is
conservative, and the intervals are wider than their nominal level (in
200 simulations with \\\theta = 20\\, \\m = 0.1\\, \\J = 1000\\, 95%
intervals covered the true values in 99.5% of the replicates). When the
likelihood runs along a ridge towards \\m = 1\\, the Ewens sampling
formula should be used.

## Local maxima

`local_maxima` lists the distinct local maxima of the likelihood found
by the optimisation (TeTame reports a second maximum in the same way).
The first row is the global maximum. Points on a bound of the search
interval are flagged (`boundary`); they are kept only when they are
local maxima of the profile likelihood there. Two flagged points with
the same likelihood typically mark the two ends of a flat ridge, along
which the data cannot discriminate between \\\theta\\ and \\m\\.

## Log-likelihood

`logl` is on the same scale as the `logl` of
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)
and
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)
(see
[`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md)).
The Ewens model is the Etienne model with \\m = 1\\, so the two models
are nested and can be compared by the likelihood-ratio statistic
\\2(\mathrm{logl}\_{Etienne} - \mathrm{logl}\_{Ewens})\\ or by AIC.

## Sample size

The computing time of \\K(D,A)\\ grows roughly as \\J^2\\: about 2 s for
\\10^4\\ individuals and 30 s for \\4 \times 10^4\\ on a typical
computer. Unlike the other estimators of the package, the Etienne model
is therefore practical up to about \\10^5\\ individuals. Please exercise
caution for large datasets. One practical recommendation for large
datasets is to consider them as a collection of demes and to use the
multi-deme model instead of the Hubbell model.

## References

Etienne, R. S. (2005) A new sampling formula for neutral biodiversity.
*Ecology Letters* 8: 253-260.

## See also

[`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md),
[`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md),
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md),
[`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)

## Examples

``` r
## Example of Etienne (2005), supplementary material: theta = 7.05, m = 0.226
zoo <- c(8, 5, 3, 2, 1, 1)
et <- optim_etienne(zoo)
et
#> $theta
#> [1] 7.047969
#> 
#> $m
#> [1] 0.2263591
#> 
#> $I
#> [1] 5.559197
#> 
#> $logl
#> [1] -28.56215
#> 
#> $converged
#> [1] TRUE
#> 
#> $se_theta
#> [1] 15.10683
#> 
#> $se_m
#> [1] 0.3349456
#> 
#> $se_I
#> [1] 10.63283
#> 
#> $local_maxima
#>      theta         m        I      logl boundary
#> 1 7.047969 0.2263591 5.559197 -28.56215    FALSE
#> 

## Is there evidence of dispersal limitation?
## Compare with the nested Ewens model
ew <- optim_ewens(zoo)
2 * (et$logl - ew$logl) # likelihood-ratio statistic
#> [1] 0.3577858

## Confidence intervals; a sample with two local maxima (TeTame test data)
optim_etienne(zoo, ci = TRUE)[c("ci_theta", "ci_m")]
#> $ci_theta
#>       lower       upper 
#>    0.811495 2000.000000 
#> 
#> $ci_m
#>      lower      upper 
#> 0.04102535 0.99999999 
#> 
optim_etienne(c(15, 81, 80, 2, 1, 2))$local_maxima
#>       theta           m          I      logl boundary
#> 1 21.747365 0.006867291    1.24466 -195.2637    FALSE
#> 2  1.056303 0.980993441 9290.41505 -195.2851    FALSE
```
