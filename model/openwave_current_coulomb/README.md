# Native current and a finite Coulomb energy choice

The actual PW119 Gauss flux determines a minimum-norm longitudinal
field through the finite Green operator. This is an exact readout of
the same charge, scalar, reference and link records. Each of the six
side-two unit-charge arrivals gives the quadratic flux value 7/48,
with its original positive probability and every failure retained.

The original unit-link energy gives a linear separation cost in its
complete single-unscreened-charge sector. Adding a positive Coulomb
term cannot cancel that cost. A separate executable Hamiltonian
replaces the link penalty by the quadratic Green energy. Its physical
static trial states have path-independent Coulomb-shaped energy
differences, with explicit unpenalized flux loops and no transferred
ground-state or particle-spectrum claims.

The [proof](PROOF.md) states both constructions and the exact native
current attachment. The [receipt](results.json) retains the complete
finite source data, independent field profile and free coupling.
No physical fine-structure constant is derived.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_current_coulomb --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_current_coulomb --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_current_coulomb
```
