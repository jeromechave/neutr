# Maximum likelihood estimation of \\\theta\\ based on Ewens sampling formula

`optim_ewens()` computes the maximum likelihood estimate of Ewens
\\\theta\\ parameter for a given species abundance distribution \\(n_1,
..., n_k)\\ where \\n_i\\ is the number of organisms of species i in the
sample and k is the total number of species in the sample (for all \\i,
n_i \> 0\\). The input should be a vector; species with zero abundance
are ignored.

The maximum likelihood estimate is the root of \\k =
\theta\[\psi(\theta+J)-\psi(\theta)\]\\ (the expected number of species
equals the observed one), where \\\psi\\ is the digamma function and
\\J=\sum n_i\\. It is found by a one-dimensional root search to full
numerical precision, so no general-purpose optimiser is involved.

## Usage

``` r
optim_ewens(input_abundances)
```

## Arguments

- input_abundances:

  a vector of integers (species abundances)

## Value

A list with

- theta:

  the \\\theta\\ parameter of Ewens sampling formula, a real value

- logl:

  the maximal value of the log-likelihood (see Details), a real value;
  it is comparable with the `logl` value of
  [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)

## Boundary cases

The likelihood function has no interior maximum when all individuals
belong to a single species (`k = 1`, the estimate is `theta = 0`) or
when every species has a single individual (`k = J`, the estimate is
`theta = Inf`). A warning is issued and the boundary value is returned.

## Log-likelihood

`logl` is the log-probability of the partition (the exchangeable
partition probability function of Ewens sampling formula). Terms that
depend on the data only, and not on the model, are omitted. It is on the
same scale as the `logl` returned by
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md).
Note that the Ewens model is equal to the Pitman model when \\\sigma =
0\\, so the two models are nested and can be compared directly, for
instance with the likelihood-ratio statistic
\\2(\mathrm{logl}\_{Pitman} - \mathrm{logl}\_{Ewens})\\ or with AIC.

## See also

[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md),
[`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)

## Examples

``` r
input_abundances <- c(234, 87, 34, 5, 4, 3, 3, 2, 2, 1, 1, 1, 1)
optim_ewens(input_abundances)
#> $theta
#> [1] 2.476993
#> 
#> $logl
#> [1] -433.6187
#> 
```
