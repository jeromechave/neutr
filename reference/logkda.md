# Etienne's coefficients \\K(D,A)\\ of the sampling formula with dispersal limitation

`logkda()` computes the logarithm of the coefficients \\K(D,A)\\ of
Etienne's sampling formula (Etienne 2005, equation A11), for a species
abundance distribution \\D = (n_1, ..., n_k)\\ and all \\A\\ from k to
\\J=\sum n_i\\. These coefficients depend on the data only, not on the
parameters, and are the expensive part of
[`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md)
and
[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md).

\\K(D,A)\\ is the coefficient of \\x^A\\ in the polynomial
\\\prod\_{i=1}^k T\_{n_i}(x)\\, where \\T_n(x)=\sum_a T_n\[a\] x^a\\
with \\T_1\[1\]=1\\ and \\T_n\[a\] = T\_{n-1}\[a\] +
T\_{n-1}\[a-1\](a-1)/(n-1)\\, i.e. \\T_n\[a\] =
\|s(n,a)\|(a-1)!/(n-1)!\\ with \\s\\ the Stirling numbers of the first
kind. The recursion and the polynomial products are evaluated entirely
on the log scale, so that the calculation cannot overflow whatever the
sample size.

## Usage

``` r
logkda(input_abundances)
```

## Arguments

- input_abundances:

  a vector of integers (species abundances); species with zero abundance
  are ignored

## Value

\\\log K(D,A)\\ for \\A = k, k+1, ..., J\\, a vector of `J - k + 1` real
values (the first and the last are always 0)

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
[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)

## Examples

``` r
logkda(c(8, 5, 3, 2, 1, 1))
#>  [1] 0.000000 1.970769 3.322156 4.324023 5.073586 5.613912 5.963897 6.128735
#>  [9] 6.103853 5.876264 5.423493 4.709155 3.672072 2.197225 0.000000

## a single species of abundance n gives the Stirling numbers |s(n, a)|
exp(logkda(5)) * factorial(4) / factorial(0:4)
#> [1] 24 50 35 10  1
```
