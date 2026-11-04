# Searching every one-particle mode

The unchanged periodic Weyl/contact source has a complete energy detector on
all 1,024 one-particle states of its eight-site-per-axis cube. The search
finds sixteen zero modes and six pairs of signed energy shells, with exact
multiplicities. Native spin rotations complete the source's Fourier readout,
including every pointer and abort outcome. A stationary one-particle state
cannot occupy only six sites; the infinite lattice has no square-summable
stationary one-particle state. These results constrain this sector's particle
search. They identify neither observed species nor localized many-body defects.

The [proof](PROOF.md) gives the spectral instrument, its native error bound
and the localization theorems. The [receipt](results.json) binds the source
and retains both exact calculations and independent complete-column checks.
The source Hamiltonian, its energy origin, contact law and mode count are
unchanged. The observer reads their specification to construct the detector.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.one_particle_search --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/one_particle_search
```
