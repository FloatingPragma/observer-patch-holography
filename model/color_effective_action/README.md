# Color response after eliminating Higgs states

The same finite gauge–Higgs Hamiltonian determines an exact six-state
effective action after all 432 occupied-scalar states are integrated out.
Exact intervals expose additional flux transitions and scale-dependent
electric and magnetic coefficients. The local raising links also fail a
direct Wilson-unitarity test. The [proof](PROOF.md) and source-bound
receipt retain the full action; identifying its controlled
continuum limit with a physical strong coupling requires an additional
source-to-field construction.

The separately implemented U(1) rotor completion has a
[cutoff-uniform resolvent and action limit](REGULATOR_PROOF.md) on its
fixed spatial graph. Its coordinate dynamics and source identity are
kept distinct from the finite raising-link Higgs model.
The [nonabelian coordinate theorem](NONABELIAN_REGULATOR.md) supplies
explicit representation-growth conditions and a conditional regulator
bound, together with a finite experiment distinguishing that completion
from the original source.

From the repository root:

```sh
python -m model.color_effective_action --check model/color_effective_action/results.json
python -m pytest -q model/color_effective_action
```
