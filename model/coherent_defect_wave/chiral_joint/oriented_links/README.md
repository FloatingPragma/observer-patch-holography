# Quantum links and cubic rotations

A 119-state gauge link permits spatial orientation reversal while preserving
the actual native matter currents. Its transport matrices have norm one, and
the complete three-dimensional periodic quantum Hamiltonian has all 24 proper
cubic rotations without conjugating matter charges. The executable retains
every fermion occupation, scalar state and link state. Its link architecture
is declared explicitly; a relativistic limit and measured-particle identification
require a physical reference and separate estimates.

```sh
python -m model.coherent_defect_wave.chiral_joint.oriented_links --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/oriented_links
```

See the [derivation](PROOF.md), [exact operator API](model.py) and
[algebraic receipt](results.json).
