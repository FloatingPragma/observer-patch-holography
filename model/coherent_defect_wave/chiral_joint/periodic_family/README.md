# Periodic three-axis quantum source

The full native Weyl/gauge/Higgs source has three commuting translations and
one exact finite Hamiltonian on a periodic cube. Every Fock occupation,
quantum link and scalar memory is retained. A physical eight-cell control
executes all three translations and a complete 64-outcome joint phase readout.
Its interacting ground, relativistic limit and particle masses are separate
requirements.

```sh
python -m model.coherent_defect_wave.chiral_joint.periodic_family --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/periodic_family
```

See [the proof](PROOF.md), [receipt](results.json) and
[exact source API](model.py). `Configuration` fixes the side, spacing and
uniform source coefficients. Default spacing pi gives the precise frozen
two-cell coefficients `+i/4,-i/4`. `kinetic_terms` exposes every actual
oriented path and its coefficient; `column`, `gauss_column` and
`translation_column` act on complete keys `(mask,links,scalars)`.

The default eight-cell logical dimension is `2^240 * 60^24 * 5^8`. Sparse
execution on a specified state is not a restriction of that carrier. The
source uses raw long-range quantum path products, and imports no free-field
or open-two-cell ground. Its schema is `oph.chiral_periodic_quantum_family.v1`.
