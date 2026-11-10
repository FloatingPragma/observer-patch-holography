# Preparing an interacting spin reference

A declared finite algebraic controller prepares the actual local quantum
ground, makes the eight-cell cube trial state, and compares the energies
of two copies under the full interacting Hamiltonian. Exact bounds give
high ground fidelity and retained failure probabilities. The same
collective spin readout has controlled probabilities after this preparation.
The schedule has an explicit finite length; an efficient native gate
implementation is outside its scope.

See the [proof](PROOF.md) and [receipt](results.json). Replay with
`python -m model.coherent_defect_wave.chiral_joint.interacting_ground_preparation --check`
and `pytest -q model/coherent_defect_wave/chiral_joint/interacting_ground_preparation`.
