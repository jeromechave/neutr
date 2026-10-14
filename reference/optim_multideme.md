# Maximum likelihood estimation of local immigration rates \\(m_1,...,m_K)\\ based on the K-deme sampling formula

`optim_multideme()` computes the maximum likelihood estimate for the
multi-deme parameters for a given species abundance matrix \\(n\_{1j},
..., n\_{kj})\\ for deme j (\\j\in\\1,...,K\\\\) where \\n\_{ij}\\ is
the number of organisms of species i in deme j and k is the total number
of species (for all \\i, n_i \geq 0\\). The formulation in this
description is that of species abundance distributions but is valid for
any partition in a subdivided setting (e.g., word frequencies in
multiple books).

**The matrix must have one row per deme and one column per species.**
The regional species abundance distribution is obtained by summing the
demes (rows).

## Usage

``` r
optim_multideme(
  input_abundance_matrix,
  verbose = FALSE,
  I_max = NULL,
  ci = FALSE,
  level = 0.95
)
```

## Arguments

- input_abundance_matrix:

  a matrix \\n\_{ij}\\ of abundances (one row per deme (local site) j,
  one column per species i). A data frame of counts is accepted.

- verbose:

  if `TRUE`, report progress over the demes with
  [`message()`](https://rdrr.io/r/base/message.html) (default `FALSE`)

- I_max:

  upper bound of the search interval for `I`: `NULL` (default, the deme
  size `J`), or a positive number, or one number per deme

- ci:

  if `TRUE`, also return likelihood-ratio confidence intervals (default
  `FALSE`)

- level:

  confidence level of the intervals (default 0.95)

## Value

A list with

- I:

  rescaled immigration rate, \\I=\frac{m}{1-m}(J-1)\\, a vector of real
  numbers

- m:

  immigration rate, a vector of real numbers

- J:

  local community size, a vector of integers

- k:

  number of classes/species, a vector of integers

- logl:

  maximal log-likelihood of each deme (up to a constant that does not
  depend on `I`), a vector of real numbers

- se_I, se_m:

  standard errors of `I` and `m`, vectors of real numbers

- ci_I, ci_m:

  only if `ci = TRUE`: confidence intervals of `I` and `m`, matrices
  with one row per deme and columns `lower` and `upper`

## Boundary and degenerate cases

- A deme with a single individual carries no information on the
  immigration rate: `I` and `m` are `NA` and a warning is issued.

- A deme containing a single species has its likelihood maximal at
  `I = 0` (`m = 0`); this value is returned as is.

- The search interval for `I` is `[1e-6, I_max]` with `I_max` equal to
  the deme size `J` by default. If an estimate lies on the edge of this
  interval (for instance a deme that is essentially a random sample of
  the regional pool, for which the likelihood keeps increasing with
  `I`), a warning is issued; increase `I_max` to explore larger values.

## Uncertainty

The standard errors `se_I` and `se_m` are obtained from the exact second
derivative of the log-likelihood in `I` (and the delta method for `m`);
they are `NA` when the estimate lies on a bound of the search interval
or for a deme with a single species. With `ci = TRUE`, likelihood-ratio
confidence intervals are also returned: the range of `I` over which the
log-likelihood stays within `qchisq(level, 1) / 2` of its maximum. A
limit equal to a bound of the search interval means that the interval
extends at least to that bound. These replace the `Std_I` and `Std_m`
columns of the TeTame software.

## See also

[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md),
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md),
[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)

## Examples

``` r
input_abundances1 <- c(44, 37, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
input_abundances2 <- c(240, 20, 48, 2, 21, 1, 3, 2, 5, 2, 0, 1, 1)
input_abundance_matrix <- rbind(input_abundances1, input_abundances2)
optim_multideme(input_abundance_matrix)
#> $I
#> [1]  69.44552 197.91949
#> 
#> $m
#> [1] 0.3363867 0.3645467
#> 
#> $J
#> [1] 138 346
#> 
#> $k
#> [1] 13 12
#> 
#> $logl
#> [1] -255.5496 -392.8408
#> 
#> $se_I
#> [1]  34.80389 128.98987
#> 
#> $se_m
#> [1] 0.1118761 0.1509746
#> 
```
