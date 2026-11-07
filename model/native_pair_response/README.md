# Spatial pair correlations

Two pairs of native repair records measure spatial correlations on one
fixed-occupation source. The exact response distinguishes exclusion from
an independent pair that retains collisions as failures. Complete final
readback also measures smooth quadratic occupation fields; under the
existing density-limit hypotheses these approach quadratic Gaussian
fluctuations. The [proof](PROOF.md) states the common state, clock and
measurement contracts. These source observables require a physical
spacetime and particle identification before they become laboratory claims.

Reproduce the exact checks with:

```sh
python -m model.native_pair_response.packet --check
python -m pytest -q -W error model/native_pair_response
```
