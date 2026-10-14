# Expected number of classes with more than N individuals given parameters \\\theta\\ and sample size J using the GEM representation

The expected number of species with more than `N` individuals is
estimated by averaging `nrep` independent draws of the neutral model
with the Griffiths-Engen-McCloskey representation
([`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)).

## Usage

``` r
kest_gt0(theta, J, N, nrep = 100)
```

## Arguments

- theta:

  parameter \\\theta\\, a positive real value

- J:

  number of individuals, a whole number

- N:

  minimal abundance (a species is counted if its abundance is `> N`)

- nrep:

  number of simulated communities, at least 2 (default 100)

## Value

A list with

- k:

  the expected number of classes with more than `N` individuals, a real
  number

- sigmak:

  the standard deviation of that number across the `nrep` draws, a real
  number (a measure of the spread of single communities, not the
  standard error of `k`)

## See also

[`kest_gt()`](https://jeromechave.github.io/neutr/reference/kest_gt.md)
for a deterministic approximation

## Examples

``` r
set.seed(1)
kest_gt0(11.3, 234373, 50)
#> $k
#> [1] 63.22
#> 
#> $sigmak
#> [1] 6.178457
#> 
```
