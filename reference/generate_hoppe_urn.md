# Generation of a neutral partition given parameter \\\theta\\ using the Griffiths-Engen-McCloskey representation

`generate_hoppe_urn()` creates a rank-abundance distribution for a given
\\\theta\\ and `J`. Beta-distributed random variables are drawn from
\\Beta(1,\theta)\\ to build the stick-breaking weights of the
Griffiths-Engen-McCloskey (GEM) distribution, and the partition is a
multinomial draw of `J` individuals according to these weights. Its cost
does not depend on `J` (up to the multinomial draw), which makes it
suitable for very large samples (10^12 individuals and beyond).

The draw is exact: sticks are generated in blocks until every individual
has been assigned. The total number of individuals is always exactly
`J`.

## Usage

``` r
generate_hoppe_urn(theta, J)
```

## Arguments

- theta:

  parameter \\\theta\\, a positive real value

- J:

  number of individuals, a whole number (up to 2^53)

## Value

A list with

- abundance:

  the species abundances, sorted in decreasing order and summing to `J`

- k:

  the *expected* number of species for these parameters (see
  [`kest()`](https://jeromechave.github.io/neutr/reference/kest.md));
  the number of species in this particular draw is `length(abundance)`

## Examples

``` r
set.seed(1)
out <- generate_hoppe_urn(11.3, 234373)
length(out$abundance) # close to out$k
#> [1] 123
```
