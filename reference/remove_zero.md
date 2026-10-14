# General use function

`remove_zero()` is used to remove zeros and NAs in a vector of reals or
of integers. This is a simple utility function.

## Usage

``` r
remove_zero(x)
```

## Arguments

- x:

  vector of real

## Value

the vector `x` without its zero and missing (`NA`/`NaN`) elements

## Examples

``` r
remove_zero(c(2.1, 0.0, 3.0, NA))
#> [1] 2.1 3.0
```
