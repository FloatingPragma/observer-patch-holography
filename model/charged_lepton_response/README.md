# Charged fields and candidate particle readouts

Native event weights supply a complete coherent family generator with Koide
balance. Its full fermion realization and free relativistic limit retain
the same coefficients. Choosing one phase from the muon/electron ratio gives
a tau/electron ratio compatible with the stated experimental uncertainty.
This is a calibrated action test. The interacting chiral theory and its
measured-mass correspondence require their own physical estimates.

The separate C-duty memory construction realizes a supplied pixel-dependent
mass law. Its four coordinate branches and the candidate map preserve the
source inputs and distinguish action modes from identified particles.

See the [results](RESULTS.md), [proof and dependencies](PROOF.md), and
[candidate map](CANDIDATE_MAP.md).

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.charged_lepton_response --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.charged_lepton_response --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/charged_lepton_response
```
