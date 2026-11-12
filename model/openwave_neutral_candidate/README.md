# A neutral gapped sector and its boundary dependence

The native Higgs-lepton composite is neutral under the existing gauge
currents and changes conserved lepton number by one. On the complete
finite source, its charge sector has a positive energy threshold above
the actual interacting ground and cannot decay to that ground under the
unchanged Hamiltonian.

The original antiperiodic source has no zero-momentum odd state. An
explicit periodic-boundary variant changes twelve seam signs and admits
a neutral zero-momentum threshold with the same certified bounds. The
[proof](PROOF.md) gives the boundary transfer and a conditional common-scale
GeV interval. No physical dark-matter species, neutrino identity or
relativistic mass shell is assigned.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_neutral_candidate --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_neutral_candidate --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/openwave_neutral_candidate
```
