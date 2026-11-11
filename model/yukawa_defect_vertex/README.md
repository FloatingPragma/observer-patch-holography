# A Higgs vertex measured on the finite source

The existing six-mode lepton source has a complete 320-state quantum
Hamiltonian, including both charged Higgs components and their antiparticles.
Its primitive lepton interaction is one. This package derives the complete
fermion and scalar two-point records and reconstructs that interaction from
an anomalous three-point response. A separate charged-scalar pulse produces
a nonzero pair-pointer response while retaining its quantum charge reference.

The measured scalar spectral moment is `14/5`, whereas the fermion spectral
sum is the identity. Normalizing those complete responses gives the finite
contact coefficient `sqrt(14/5)`. Its normalization comes from the source
records. The primitive unit interaction and the finite source clock remain
declared construction data. This is a one-family onsite result, with a
separate [proof](PROOF.md) of what the coefficient measures; it supplies no
on-shell particle masses or three-family hierarchy.

Every source energy, charged-probe branch and scalar no-operation outcome is
retained. The original ground-preparation probability is `1/160`; all other
preparation outcomes remain in the complete records.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/yukawa_defect_vertex
OPENBLAS_NUM_THREADS=1 python -m model.yukawa_defect_vertex --write
OPENBLAS_NUM_THREADS=1 python -m model.yukawa_defect_vertex --check
```
