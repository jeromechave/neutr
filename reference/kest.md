# Expected number of classes given parameter \\\theta\\ and sample size J for the Hoppe urn model

The expected number of species is computed using an exact expression
based on the digamma function,
\\\theta\[\psi(\theta+J)-\psi(\theta)\]\\. The function is vectorised
over `theta` and `J`.

## Usage

``` r
kest(theta, J)
```

## Arguments

- theta:

  parameter \\\theta\\, a positive real value

- J:

  number of individuals

## Value

the expected number of classes (species), a real number

## Examples

``` r
kest(11.3, 234373)
#> [1] 112.8284
```
