# Generation of a neutral partition given parameter \\\theta\\ using Hoppe's urn

`generate_hoppe_urn0()` creates a rank-abundance distribution based on
Hoppe's urn scheme for a given parameter \\\theta\\ and sample size `J`.
The process starts at zero species, corresponding to the single black
ball with weight \\\theta\\. If the black ball is picked, a new class is
created, else the picked colored ball is duplicated (its abundance
increases by one unit).

This is an *exact* simulation of the urn, but it does not loop over the
individuals: the genealogy of the `J` individuals is resolved in a
handful of vectorised passes, and the random numbers come from the fast
`dqrng` generator. Memory use is proportional to `J`, so it is meant for
`J` up to a few hundred million; for larger samples use
[`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md).

## Usage

``` r
generate_hoppe_urn0(theta, J, seed = NULL)
```

## Arguments

- theta:

  parameter \\\theta\\, a positive real value

- J:

  number of individuals, a whole number

- seed:

  optional integer; if given, the `dqrng` generator is seeded with it
  (this changes the state of the `dqrng` generator for the whole
  session)

## Value

`generate_hoppe_urn0()` returns the species abundances, an integer
vector sorted in decreasing order and summing to `J`.

## Reproducibility

`generate_hoppe_urn0()` draws its random numbers from `dqrng`, which has
its own random number stream:
[`set.seed()`](https://rdrr.io/r/base/Random.html) has **no effect** on
it. Use the `seed` argument instead.
[`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)
and
[`generate_pitman_urn()`](https://jeromechave.github.io/neutr/reference/generate_pitman_urn.md)
use R's standard generator and follow
[`set.seed()`](https://rdrr.io/r/base/Random.html).

## See also

[`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)
for a faster approach on very large samples

## Examples

``` r
generate_hoppe_urn0(11.3, 234373, seed = 1)
#>   [1] 60996 24243 21303 13817 13297  9833  7847  7082  7008  6435  5523  5501
#>  [13]  5265  3799  3324  3117  3002  2915  2772  2584  2480  2240  1865  1816
#>  [25]  1747  1708  1383  1277  1032   918   845   793   744   668   562   447
#>  [37]   417   259   242   240   231   228   212   206   186   173   152   144
#>  [49]   130   117   109   105   102   102    92    79    70    63    61    57
#>  [61]    53    42    36    29    27    26    21    18    16    11    11    10
#>  [73]    10     9     8     8     6     6     6     5     5     4     3     3
#>  [85]     3     3     2     2     2     2     2     2     2     1     1     1
#>  [97]     1     1     1     1     1     1     1     1     1     1     1     1
```
