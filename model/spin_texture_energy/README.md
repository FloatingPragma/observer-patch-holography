# Energy of topological spin textures

The ordered quantum ground of the derived cubic spin system supplies an
exact energy formula for arbitrary local spin rotations. It gives a
positive lower bound on the direction-gradient energy and a controlled
static limit for smooth textures. Compact Hopf textures in this family
have energy proportional to their size; an isolated radial hedgehog has
infinite energy. These statements concern rotations of the actual quantum
ground, with no product-state approximation.

[The derivation](PROOF.md) and [certificate](results.json) distinguish this
family from energy minimization over all states. A finite-energy particle
carried by a topological defect requires a selected state, stabilizing
dynamics and a spectral identification beyond this static result.

Run `python -m model.spin_texture_energy --check` and
`python -m pytest -q model/spin_texture_energy` from the repository root.
