# Moving excitations between interacting blocks

A lattice of the complete eight-cell interacting sources has a unique
ground, a volume-independent gap, and an isolated moving excitation band.
The coupling uses the actual native spin current dressed by its isolated
spectral projectors. The three current polarizations keep their internal
cube momenta while the excitation propagates between blocks.

The [proof](PROOF.md) derives exact finite momentum and detector laws,
controlled preparation and boundary errors, and a Schrödinger spacing
limit. Spectral dressing is a declared architecture; the result supplies
no relativistic continuum or equivalence to the uniform original source.

From the repository root, with `model/requirements.txt` installed:

```sh
python -m model.interacting_block_propagation --write
python -m model.interacting_block_propagation --check
python -m pytest -q model/interacting_block_propagation
```

[Source ancestry](PROVENANCE.md) identifies the immutable parents.
The machine receipt contains exact rational bounds and complete finite
position distributions. No measured masses enter the construction.
