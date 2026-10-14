# Deprecated functions in neutr

These functions were renamed in neutr 0.2.1 to follow snake_case naming.
The old names still work but are deprecated and will be removed in a
future release. Please use the new names instead:

|  |  |
|----|----|
| Deprecated | Replacement |
| `generate.hoppe.urn()` | [`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md) |
| `generate.hoppe.urn0()` | [`generate_hoppe_urn0()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn0.md) |
| `generate.pitman.urn()` | [`generate_pitman_urn()`](https://jeromechave.github.io/neutr/reference/generate_pitman_urn.md) |
| `kest.gt()` | [`kest_gt()`](https://jeromechave.github.io/neutr/reference/kest_gt.md) |
| `kest.gt0()` | [`kest_gt0()`](https://jeromechave.github.io/neutr/reference/kest_gt0.md) |
| `logl.etienne()` | [`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md) |
| `optim.etienne()` | [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md) |
| `optim.ewens()` | [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md) |
| `optim.multideme()` | [`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md) |
| `optim.pitman()` | [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md) |
| `zeroremove()` | [`remove_zero()`](https://jeromechave.github.io/neutr/reference/remove_zero.md) |

The deprecated functions take the same arguments and return the same
values as their replacements. They signal a soft deprecation warning
with
[`lifecycle::deprecate_soft()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html):
it is shown when a deprecated function is called from the global
environment or from tests, but not when it is called from another
package.

## Usage

``` r
generate.hoppe.urn(theta, J)

generate.hoppe.urn0(theta, J, seed = NULL)

generate.pitman.urn(theta, sigma, J)

kest.gt(theta, J, N)

kest.gt0(theta, J, N, nrep = 100)

logl.etienne(input_abundances, theta, m, log_kda = NULL, full = FALSE)

optim.etienne(
  input_abundances,
  init_vals = c(10, 0.5),
  lower = c(1e-06, 1e-08),
  upper = NULL,
  ci = FALSE,
  level = 0.95
)

optim.ewens(input_abundances)

optim.multideme(
  input_abundance_matrix,
  verbose = FALSE,
  I_max = NULL,
  ci = FALSE,
  level = 0.95
)

optim.pitman(
  input_abundances,
  init_vals = c(10, 0.1),
  lower = c(1e-06, 1e-08),
  upper = NULL
)

zeroremove(x)
```

## Arguments

- theta, sigma, J, seed, N, nrep, input_abundances,
  input_abundance_matrix:

  See the documentation of the replacement function.

- m, log_kda, full, init_vals, lower, upper, ci, level, verbose, I_max,
  x:

  See the documentation of the replacement function.

## Value

See the documentation of the replacement function.
