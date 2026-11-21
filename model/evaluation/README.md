# Evaluation of candidate laws

This package scores every candidate law against the [requirements](../REQUIREMENTS.md): a law runs on the fixed
spherical echosahedral federation and fulfils the premises of the finite theorems together, with every
calibration declared and no outcome hardcoded.

```sh
python -m model evaluate              # print the evaluation of every candidate
python -m model evaluate --law NAME   # one candidate
python -m model evaluate --write      # write each EVALUATION.json, SCORECARD.md and the table of CANDIDATE.md
python -m model evaluate --check      # fail when a written file differs from a fresh evaluation
python -m pytest -q model/evaluation
```

## What is checked

| Part | How |
| --- | --- |
| Native cell | the candidate's labelling of one cell is isomorphic to the fixed echosahedron |
| Locality | every nonzero coupling of the candidate's generator joins incident registers of one cell |
| Rotations | the generator commutes with the sixty rotations of the cell, mapped into the candidate's labelling |
| Rank-three response | the slowest nonconstant level of the candidate's response operator has multiplicity three |
| Federations | the generator builds, Hermitian, on the federations of 20, 80 and 320 cells |
| Parameters | every parameter declares its provenance; a calibrated one names its observation; none declares a particle mass, mass ratio, mixing angle, number of dimensions, metric, lattice or gauge group |
| Premises | every row of [premise_register.json](../../tracking/premise_register.json), with a role from [roles.json](roles.json) and the candidate's status and evidence; undeclared rows are open, a failed row names its counterexample, and the evidence files must exist |
| Joint requirements | one law for all premises, the Einstein stress carried by the propagating excitations, gauge currents acting on the same excitations, one patch measuring both, and one dictionary for rest energy and inertia |

A premise of type `structural_rule` is supplied by the law. Selection rules, representation choices and numerical
inputs may be supplied by the designer as implementable calibration, except where [roles.json](roles.json) assigns
them to the law because they would encode an observed outcome. Empirical imports and external mathematics are
imports.

The carrier checks run on each candidate's own code, so a change to a law changes its evaluation. The premise
statuses are declarations backed by evidence files; a declaration becomes a check when the candidate's evidence
includes the test or certificate that proves it.
