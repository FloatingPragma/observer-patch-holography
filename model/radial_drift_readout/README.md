# Observing a defect while its Hamiltonian acts

The same interacting radial field can be read repeatedly during a positive
Hamiltonian-time interval. Only the origin's incident hopping and contact
enter the error bound, so the bound is independent of box volume. The
[proof](PROOF.md) retains the complete field, reference and native archive.
For a supplied timed local implementation lasting at most one millionth of
a source time unit, the one-preparation ground-versus-charged classifier
has error below one in 2,500 without compensating the Hamiltonian.

The [receipt](results.json) contains exact local norms, the finite pulse-rate
schedule and independent full-Fock controls. The timing and native-port
embedding are specified apparatus requirements; the integrated source word
alone does not establish physical bandwidth. All failures are scored and
retained, and the measured field need not remain in its initial energy sector.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.radial_drift_readout --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m pytest -q model/radial_drift_readout
```
