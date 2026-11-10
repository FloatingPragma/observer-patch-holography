# Collective spin records

A spin intervention on the eight-cell interacting ground has a strictly
positive first excitation line and a strictly positive weight outside that
line. Both use the actual source Hamiltonian at hopping coefficient
`2^-34`; the three source momenta are the three X points. Exact algebra and
full complementary-space bounds control the probabilities. The preparation
filter's failure probability `1-1/D` is retained. The result concerns a
finite collective response, with no relativistic mass identification.

See the [proof](PROOF.md) and [receipt](results.json). Reproduce with
`python -m model.coherent_defect_wave.chiral_joint.collective_spin_response --check`
and `pytest -q model/coherent_defect_wave/chiral_joint/collective_spin_response`.
