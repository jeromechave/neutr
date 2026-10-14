# Generation of a Pitman partition given parameters \\\theta\\ and \\\sigma\\ from the Griffiths-Engen-McCloskey representation

The Pitman urn scheme is as follows: starting at zero species,
corresponding to the single black ball with weight \\\theta\\, a species
i, of k species total, is selected with probability
\\(n_i-\sigma)/(n+\theta)\\, where \\\sigma\\ is strictly between 0.0
and 1.0 and n is the current number of balls; the black ball is picked
with probability \\(k\sigma+\theta)/(n+\theta)\\.
`generate_pitman_urn()` creates a rank-abundance distribution for given
\\\theta\\, \\\sigma\\ and `J`. Beta-distributed random variables are
drawn from \\Beta(1-\sigma,\theta+i\sigma)\\ (for the i-th stick) and
the partition is the multinomial draw of `J` individuals according to
the resulting Griffiths-Engen-McCloskey (GEM) weights.

The draw is exact: sticks are generated in blocks until every individual
has been assigned, and the few individuals left in the far tail are
drawn with the exact sequential Pitman urn. The total number of
individuals is always exactly `J`.

## Usage

``` r
generate_pitman_urn(theta, sigma, J)
```

## Arguments

- theta:

  parameter \\\theta\\, a real value \\\ge 0\\

- sigma:

  parameter \\\sigma\\, a real value strictly between 0 and 1 (use
  [`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)
  for \\\sigma = 0\\)

- J:

  number of individuals, a whole number (up to 2^53)

## Value

A list with

- abundance:

  the species abundances, sorted in decreasing order and summing to `J`

- k:

  the *expected* number of species for these parameters; the number of
  species in this particular draw is `length(abundance)`

## Examples

``` r
set.seed(1)
generate_pitman_urn(11.3, 0.1, 234373)$k
#> [1] 193.5583
```
