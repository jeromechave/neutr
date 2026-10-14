# neutr: Inference and Fast Generation of Neutral Species Assemblages

Parameter inference and fast generation for neutral species assemblage
models: the Ewens model, the Etienne model (a local community with
dispersal limitation), the multideme model and the Pitman model.

## Estimation

Parameters are estimated by maximising the likelihood:
[`optim_ewens()`](https://jeromechave.github.io/neutr/reference/optim_ewens.md)
for the Ewens model,
[`optim_etienne()`](https://jeromechave.github.io/neutr/reference/optim_etienne.md)
for the Etienne model,
[`optim_pitman()`](https://jeromechave.github.io/neutr/reference/optim_pitman.md)
for the Pitman model and
[`optim_multideme()`](https://jeromechave.github.io/neutr/reference/optim_multideme.md)
for the multideme model. The Ewens model is nested in both the Etienne
model (\\m = 1\\) and the Pitman model (\\\sigma = 0\\), and the three
functions return their log-likelihood on the same scale.

## Etienne sampling formula

[`logl_etienne()`](https://jeromechave.github.io/neutr/reference/logl_etienne.md)
evaluates the likelihood of Etienne's sampling formula and
[`logkda()`](https://jeromechave.github.io/neutr/reference/logkda.md)
its coefficients \\K(D,A)\\, computed on the log scale so that large
samples cannot overflow.

## Generation

[`generate_hoppe_urn0()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn0.md)
simulates Hoppe's urn exactly, individual by individual;
[`generate_hoppe_urn()`](https://jeromechave.github.io/neutr/reference/generate_hoppe_urn.md)
and
[`generate_pitman_urn()`](https://jeromechave.github.io/neutr/reference/generate_pitman_urn.md)
use the Griffiths-Engen-McCloskey representation, whose cost does not
depend on the sample size, and scale to 10^12 individuals and beyond.

## Expected richness

[`kest()`](https://jeromechave.github.io/neutr/reference/kest.md),
[`kest_gt()`](https://jeromechave.github.io/neutr/reference/kest_gt.md)
and
[`kest_gt0()`](https://jeromechave.github.io/neutr/reference/kest_gt0.md)
give the expected number of species (possibly above an abundance
threshold).

## See also

Useful links:

- <https://github.com/jeromechave/neutr>

- Report bugs at <https://github.com/jeromechave/neutr/issues>

## Author

Jerome Chave
