# Boundary records and regional information

The original selected Wilson state has strictly nonzero coherence between
regional boundary-charge sectors. Recording those sectors leaves the
regional marginal entropy unchanged while disturbing the global source.
The physical Schur decomposition gives the exact difference between
ambient and gauge-invariant regional-algebra mutual information: twice
the boundary representation term plus the lost center coherence.

The actual source action, polynomial modular slope and finite support
floor give positive bounds on all three terms. They deteriorate with
capacity and supply no geometric area law or Newton constant.
[Proof](PROOF.md) distinguishes these entropy domains and retains the
complete source, reference and record outputs.

```sh
python -m model.particle_primary_boundary_information --check
python -m pytest model/particle_primary_boundary_information/test_model.py -q -W error
```
