# Expected number of classes with more than N individuals given parameters \\\theta\\ and sample size J using the species frequency spectrum

The expected number of species \\\>N\\ is computed using the
approximation based on the species frequency spectrum \\\psi(x)=\theta
x^{-1}(1-x)^{\theta-1}\\ such that the number of species \\\>N\\ is the
integral from N/J to 1 of \\(1-(1-x)^J)\psi(x)\\. This is an
approximation of the exact formula, but it turns out to be accurate.

## Usage

``` r
kest_gt(theta, J, N)
```

## Arguments

- theta:

  parameter \\\theta\\, a positive real value

- J:

  number of individuals, a whole number

- N:

  minimal abundance, with `0 <= N < J`

## Value

An object of class `"integrate"` as returned by
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html): the
expected number of classes \\\>N\\ is its `$value` component, and
`$abs.error` is the estimated absolute error of the numerical
integration.

## Examples

``` r
kest_gt(11.3, 234373, 50)$value
#> [1] 62.12428
```
