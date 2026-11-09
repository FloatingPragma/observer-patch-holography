# Nonzero charged spectral weight on the native208 source

The lowest paired charge thresholds of the unchanged finite source are
visible to its complete CAR detector. Exact spectral-window counts,
a rational ground trial and a complete rational pointer apparatus certify
an original weight close to0.498558 for each charge sign, before any
success normalization. The source excitation energy is approximately
1.107713183 in its declared clock units.

The trial's trace-distance error is below1.065e-7. All40 detector records
are evaluated as exact rational probabilities, including each first-failure
label and the charge/parity-compensating references. A separate source-clock
Fourier plan bounds a narrower energy readout while retaining its finite
resolution, clock error and full outcome space.

These are accessible finite source excitations. The calculation establishes
no relativistic mass shell, physical spin assignment, GeV mass or observed
particle identification. The [proof](PROOF.md) states the precise source,
state, reference and apparatus attachments. Inherited inputs are pinned
in `source_pins.json`; the exact outward certificate has no floating
eigensolver dependency.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -W error -m pytest -q model/particle_charged_residue
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_charged_residue --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_charged_residue --check
```
