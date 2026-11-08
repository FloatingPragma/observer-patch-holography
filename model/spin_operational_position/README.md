# Measuring a finite spin coordinate through source controls

The selected eight-site spin source has a complete coordinate measurement
compiled into the existing six-level control service. Its eight coordinate
outcomes and dark complement retain the corresponding state changes and
every source failure. The 49,228-pulse decoded core is executed with an
untouched reference; the full 49,324-pulse encoded detector has an analytic
all-input error bound. These finite labels are reflection coordinates.
Continuum detector convergence and particle identification require separate
proofs.

[PROOF.md](PROOF.md) also gives a constructive source implementation of every
computable finite quantum instrument under the same declared coherent
memory, pointer and refinable control interface. The complex rank-one
compiler uses exact rational angle algorithms. It does not infer these
services from bare occupation repair or claim efficient large-register
execution.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.spin_operational_position --check model/spin_operational_position/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --strict-markers --strict-config model/spin_operational_position
```
