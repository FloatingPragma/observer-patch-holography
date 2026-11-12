# Prepared charged readouts on one interacting source

A prepared charged energy line of the interacting eight-cell source stays
localized near its fixed compensating reference for a certified finite
duration. The bound includes screened configurations and every preparation
failure. Native spin rotation gives a minus sign after one full turn and
the identity after two, with a quantum reference making the relative phase
observable. One electron calibration converts the source transition energy
to GeV; identifying an autonomous physical electron requires a separate
particle and reference comparison.

The [localization proof](PROOF.md) derives the localization and persistence bounds from
the actual charged-source projector estimate. The [receipt](results.json)
contains both native lepton channels, a million-source-time-unit horizon,
complete record probabilities and the finite rotation controls.

The [motion theorem](MOTION.md) connects an actual full-H energy selection,
an explicit position preparation and subsequent charged arrivals under the
same clock. Both channels have a strictly positive arrival probability at
source time three that exceeds its upper bound at time six. This is a
finite, off-resonant flux excursion with retained position-measurement
disturbance, rather than a dispersion relation for a free electron.

The [static-response theorem](MAGNETIC.md) proves that this source absorbs
arbitrary central hypercharge link backgrounds by unitary conjugation.
Its spectrum and complete paired preparation/readout laws are unchanged.
The current operator is nonzero, and a fixed-state control responds to a
background quench. The result identifies why the separate free-envelope
magnetic response cannot be transferred to this interacting source.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_particle_readouts --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_particle_readouts --check
python -m pytest -q model/openwave_particle_readouts
```
