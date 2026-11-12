# Localized charged readouts and their spin phase

A prepared charged energy line of the interacting eight-cell source stays
localized near its fixed compensating reference for a certified finite
duration. The bound includes screened configurations and every preparation
failure. Native spin rotation gives a minus sign after one full turn and
the identity after two, with a quantum reference making the relative phase
observable. One electron calibration converts the source transition energy
to GeV; identifying an autonomous physical electron requires a separate
particle and reference comparison.

The [proof](PROOF.md) derives the localization and persistence bounds from
the actual charged-source projector estimate. The [receipt](results.json)
contains both native lepton channels, a million-source-time-unit horizon,
complete record probabilities and the finite rotation controls.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_particle_readouts --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_particle_readouts --check
python -m pytest -q model/openwave_particle_readouts
```
