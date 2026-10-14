# Package index

## Parameter estimation

Maximum likelihood estimation of the parameters of the four neutral
models. The `logl` values of
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md),
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)
and
[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
are on the same scale, so nested models can be compared directly.

- [`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)
  : Maximum likelihood estimation of \\\theta\\ based on Ewens sampling
  formula
- [`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)
  : Maximum likelihood estimation of \\\theta\\ and \\\sigma\\ based on
  Pitman sampling formula
- [`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
  : Maximum likelihood estimation of \\\theta\\ and \\m\\ based on
  Etienne sampling formula
- [`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)
  : Maximum likelihood estimation of local immigration rates
  \\(m_1,...,m_K)\\ based on the K-deme sampling formula

## Etienne sampling formula

Likelihood of Etienne’s sampling formula and its coefficients.

- [`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md)
  : Log-likelihood of Etienne's sampling formula
- [`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md)
  : Etienne's coefficients \\K(D,A)\\ of the sampling formula with
  dispersal limitation

## Generation of species abundance distributions

Exact simulation of the Ewens and Pitman models.

- [`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)
  : Generation of a neutral partition given parameter \\\theta\\ using
  the Griffiths-Engen-McCloskey representation
- [`generate_hoppe_urn0()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn0.md)
  : Generation of a neutral partition given parameter \\\theta\\ using
  Hoppe's urn
- [`generate_pitman_urn()`](https://jeromechave.github.io/neutr/reference/generate_pitman_urn.md)
  : Generation of a Pitman partition given parameters \\\theta\\ and
  \\\sigma\\ from the Griffiths-Engen-McCloskey representation

## Number of species

Expected number of species, overall or above an abundance threshold.

- [`kest()`](https://jeromechave.github.io/neutr/reference/kest.md) :
  Expected number of classes given parameter \\\theta\\ and sample size
  J for the Hoppe urn model
- [`kest_gt()`](https://jeromechave.github.io/neutr/reference/kest_gt.md)
  : Expected number of classes with more than N individuals given
  parameters \\\theta\\ and sample size J using the species frequency
  spectrum
- [`kest_gt0()`](https://jeromechave.github.io/neutr/reference/kest_gt0.md)
  : Expected number of classes with more than N individuals given
  parameters \\\theta\\ and sample size J using the GEM representation

## Utilities

- [`remove_zero()`](https://jeromechave.github.io/neutr/reference/remove_zero.md)
  : General use function

## Deprecated functions

Dotted names replaced by snake_case names in neutr 0.2.1.

- [`generate.hoppe.urn()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`generate.hoppe.urn0()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`generate.pitman.urn()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`kest.gt()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`kest.gt0()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`logl.etienne()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`optim.etienne()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`optim.ewens()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`optim.multideme()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`optim.pitman()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  [`zeroremove()`](https://jeromechave.github.io/neutr/reference/neutr-deprecated.md)
  : Deprecated functions in neutr
