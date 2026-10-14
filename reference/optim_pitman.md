# Maximum likelihood estimation of \\\theta\\ and \\\sigma\\ based on Pitman sampling formula

`optim_pitman()` computes the maximum likelihood estimate of Pitman
\\(\theta,\sigma)\\ parameters for a given species abundance
distribution \\(n_1, ..., n_k)\\ where \\n_i\\ is the number of
organisms of species i in the sample and k is the total number of
species (for all \\i, n_i \> 0\\). The input should be a vector and
initial values of \\(\theta,\sigma)\\ can be supplied with the argument
`init_vals = c(theta_init, sigma_init)` with the condition
`theta_init > 0` and `0 < sigma_init < 1`.

The log-likelihood function is maximised with the box-constrained
quasi-Newton method L-BFGS-B
([`stats::optim()`](https://rdrr.io/r/stats/optim.html)) using the exact
gradient, on the scale \\(\log\theta,\sigma)\\. The optimisation is run
from `init_vals` and from a second starting point derived from the Ewens
estimate, and the better result is kept.

## Usage

``` r
optim_pitman(
  input_abundances,
  init_vals = c(10, 0.1),
  lower = c(1e-06, 1e-08),
  upper = NULL
)
```

## Arguments

- input_abundances:

  a vector of integers (species abundances); species with zero abundance
  are ignored

- init_vals:

  initial values of (theta, sigma), a vector of two real numbers

- lower, upper:

  bounds of the search interval for `c(theta, sigma)`

## Value

A list with

- theta:

  parameter \\\theta\\, a real value

- sigma:

  parameter \\\sigma\\, a real value

- logl:

  maximal log-likelihood (log-probability of the partition, see
  [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)),
  a real value; it is comparable with the `logl` value of
  [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)

- converged:

  `TRUE` if the optimiser reported convergence

## Search interval and boundary warnings

By default \\\theta\\ is searched in `[1e-6, J]` and \\\sigma\\ in
`[1e-8, 1 - 1e-8]`, where `J` is the number of individuals. A warning is
issued if an estimate lies on the edge of its interval (the likelihood
then has no interior maximum, or the interval is too narrow); use
`lower` and `upper` to change it.

## See also

[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md),
[`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)

## Examples

``` r
optim_pitman(c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1), c(10.0, 0.1))
#> $theta
#> [1] 0.7110667
#> 
#> $sigma
#> [1] 0.2781073
#> 
#> $logl
#> [1] -432.4741
#> 
#> $converged
#> [1] TRUE
#> 
```
