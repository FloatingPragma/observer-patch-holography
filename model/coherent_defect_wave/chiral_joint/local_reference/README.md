# The ground of the quantum cell

The full native quantum cell has an interacting, gauge-invariant ground
doublet with spin one-half. Exact collective-algebra bounds cover every
local physical charge sector, and a complete 104-state reduction gives
the ground energy and positive gap, approximately 0.0895553 in the declared
Hamiltonian units. The rotation-invariant reference is the normalized
ground projector. Its mean native fermion occupation is fifteen.

The [proof](PROOF.md) and [receipt](results.json) apply to one cell with
zero Gauss charges. A stationary many-cell reference requires control of
the specified intercell interactions and their charged local sectors.

Run from the repository root with the pinned `model/requirements.txt`:

```sh
python -m model.coherent_defect_wave.chiral_joint.local_reference --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/local_reference
```

`--write` regenerates the receipt. Both paths use exact rational algebra
and root isolation; neither uses a numerical eigenbasis proposal.
