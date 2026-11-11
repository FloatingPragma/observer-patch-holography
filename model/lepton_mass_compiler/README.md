# From the lepton readout to a free particle mass

The complete lepton analyzer can drive a specified fermion Hamiltonian whose
three free rest energies equal the analyzer's squared amplitudes. This is an
explicit operator construction with a controlled limit, including every failed
reset and its quantum archive. The same kinetic law gives free Dirac propagation
with all eight lattice valleys retained. One empirical electron anchor converts
the three coefficients to 0.00051099895069, 0.1056594144959 and 1.7769849708134 GeV.

The [proof](PROOF.md) and [receipt](results.json) state the complete source,
reference, pulse clock and precision cost. Identification with observed
interacting leptons requires physical clock and field matching; the programmed
phase, neutral scalar background and ancillary preparation are declared inputs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.lepton_mass_compiler --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/lepton_mass_compiler
```
