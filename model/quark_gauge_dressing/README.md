# Electric dressing of native quark records

An exact particle-hole symmetry of the full compact source exchanges its
up and down quantum interactions. A declared electric-field extension
breaks the equality of their charged energy records: exact Gauss constraints
on the committed cone produce a positive charge-squared field energy.
Every spectral line keeps its original probability, including zero weights
and the failed probe outcome. One common electric coefficient is used.

The [proof](PROOF.md) distinguishes this geometry-dependent electrostatic
energy from an intrinsic rest mass. The [receipt](results.json) retains the
complete graph, symmetry and charged records.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/quark_gauge_dressing
OPENBLAS_NUM_THREADS=1 python -m model.quark_gauge_dressing --write
OPENBLAS_NUM_THREADS=1 python -m model.quark_gauge_dressing --check
```
