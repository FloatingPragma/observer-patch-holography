# Quantum precession with complete field outputs

The cap-two quantum source supports finite-time Larmor dynamics with the
same rest energy, current derivative, field and clock as its static magnetic
coefficient. Exact angular-sector Ritz certificates enclose finite-field
splittings; a uniform propagator bound retains arbitrary reference systems.
Complete records include joint photon-number/angular outputs and every
success/failure state coefficient.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_quantum_larmor --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_quantum_larmor --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_quantum_larmor
```

`PROOF.md` states the analytic finite-time theorem and interval branch
certificate. Physical electron identification and magnetic-field calibration
are not inferred from these conditional confined-source records.
