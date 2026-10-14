# Log-likelihood of Etienne's sampling formula

`logl_etienne()` computes the log-likelihood of Etienne's sampling
formula for a species abundance distribution \\(n_1, ..., n_k)\\, given
the fundamental biodiversity number \\\theta\\ of the metacommunity and
the immigration rate \\m\\ into the local community. With the rescaled
immigration rate \\I=\frac{m}{1-m}(J-1)\\ and \\J=\sum n_i\\, the
log-likelihood is \$\$k\log\theta + \log\Gamma(\theta) +
\sum_i\log\Gamma(n_i) + \log\Gamma(I) - \log\Gamma(I+J) +
\log\sum\_{A=k}^{J} K(D,A)\frac{I^A}{\Gamma(\theta+A)}\$\$ with
\\K(D,A)\\ given by
[`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md).
The function is vectorised over `theta` and `m`.

## Usage

``` r
logl_etienne(input_abundances, theta, m, log_kda = NULL, full = FALSE)
```

## Arguments

- input_abundances:

  a vector of integers (species abundances); species with zero abundance
  are possible but are ignored in the inference

- theta:

  fundamental biodiversity number \\\theta\\, positive real values

- m:

  immigration rate, real values in \\(0, 1\]\\

- log_kda:

  optional, the output of `logkda(input_abundances)`. It depends on the
  data only and is the expensive part of the calculation: supply it when
  the likelihood is evaluated many times on the same data

- full:

  if `TRUE`, return the log-probability of the species abundance
  distribution rather than of the partition (see Details); default
  `FALSE`

## Value

the log-likelihood, a vector of real values (one per value of `theta`
and `m`, recycled to a common length)

## Log-likelihood

By default `logl_etienne()` returns the log-probability of the
partition, on the same scale as the `logl` of
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)
and
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md).
For \\m = 1\\ (no dispersal limitation) Etienne's formula reduces to
Ewens sampling formula and the two log-likelihoods are identical, so the
Ewens model is nested in the Etienne model and the two can be compared
by likelihood ratio or AIC.

With `full = TRUE` the log-probability of the species abundance
distribution itself (Etienne 2005, equation 6) is returned instead: it
adds the term \\\log J! - \sum_i \log n_i! - \sum_j \log \Phi_j!\\,
where \\\Phi_j\\ is the number of species with \\j\\ individuals, which
depends on the data only.

## Acknowledgement

This function is adapted from the `etienne()` function of package `untb`
by Robin K. S. Hankin (Hankin 2007). The present version computes the
coefficients \\K(D,A)\\ on the log scale in pure R (no PARI/GP),
evaluates the likelihood without loss of precision when \\m\\ is close
to 1, and puts it on the same scale as the Ewens likelihood (i.e.,
`logl_etienne` is comparable with `logl.ewens`).

## References

Etienne, R. S. (2005) A new sampling formula for neutral biodiversity.
*Ecology Letters* 8: 253-260.

Hankin, R. K. S. (2007) Introducing untb, an R package for simulating
ecological drift under the unified neutral theory of biodiversity.
*Journal of Statistical Software* 22(12).

## See also

[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md),
[`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md),
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)

## Examples

``` r
zoo <- c(8, 5, 3, 2, 1, 1)
logl_etienne(zoo, theta = 7, m = 0.2)
#> [1] -28.58198

## profile over m, computing K(D,A) once
lk <- logkda(zoo)
logl_etienne(zoo, theta = 7, m = c(0.01, 0.1, 0.2, 0.5, 1), log_kda = lk)
#> [1] -37.47753 -29.36160 -28.58198 -29.24036 -30.61063

## m = 1 is Ewens sampling formula
all.equal(
  logl_etienne(zoo, theta = 3, m = 1),
  length(zoo) * log(3) + lgamma(3) - lgamma(3 + sum(zoo)) + sum(lgamma(zoo))
)
#> [1] TRUE
```
