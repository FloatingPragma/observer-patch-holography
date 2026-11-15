# Local-clock reduction and scrutiny of record gluing

[DERIVATION.md](DERIVATION.md) gives the universal analytic arguments and their
physical premises; [CONTRACT.md](CONTRACT.md) fixes objective, deliverables and
exit. This is a self-contained reduction of the proposed RG law stated in
`../source_selection_model/RECORD_GLUING.md`.

The new results are:

* Tetrahedral carrier symmetry makes a positive quadratic local clock round.
  Time-faithful refinement then excludes all faster coarse dependencies;
  directional coverage and finite composition attain every strict deadline.
  Exact scaling and continuous rotation symmetry are unnecessary.
* NOT/Toffoli with clean cofinal classical workspace implement every
  finite classical map. An explicit matrix-algebra witness makes classical
  copying compatible with twelve central ports and usable quantum seams.
  Actual gate-transition constraints also give a finite classical entropy
  minimum equal to the executed input ensemble, including defective outputs
  when the program is defective; feasibility is independent of target success.
* Lossless independent records impose a transcript-capacity bound. Literal
  dense M1 observation has at least order q^5 read incidences per fixed
  horizon and order q^2 fresh input bandwidth per interior receiver.

LC1--LC3 are proposed physical premises. The native all-level source
completion, response completeness and simulator binding are not certified.
Finite algebra compatibility is not full A1--A3 membership. The entropy
Hessian motivates a candidate quadratic form but supplies no time calibration.

Run from the repository root:

```sh
python code/rg_principle/receipt.py verify
python -m pytest -q code/rg_principle
cd Lean
lake build Geometry.RecordGluingPrincipleAxiomAudit
```

`receipt.py build` regenerates the compact receipt. It is a summary of full
recomputation, not an archive of all execution traces. The verifier checks
source bytes, replays every input of all 263 circuits through a separate
interpreter, checks the complete named catalog against independent target
functions, retained input, output truth table, clean work and inverse
execution, and runs 7,880 localized intermediate interventions. It executes
442,496 gates including inverse and intervention replays and 40 paired
intermediate-writer lowering probes. Gates and auxiliary
registers of the reference circuits are counted exactly. These are software
reference counts, not measurements of a native source compiler.

Exact symbolic controls check the carrier rotations, invariant quadratic
forms, a nonround invariant norm, nonlinear latency refinement, entropy
susceptibility and the classical-copy Kraus map on all 1,728 algebra basis
units. Integer controls exercise the dense read bounds. Analytic proofs
establish the all-size statements. Twelve Lean theorems cover quadratic
isotropy, path/limit bounds, finite strict deadlines, information capacity,
reversible evaluation and the interior-cube inequalities. The transitive
axiom audit permits only `propext`, `Classical.choice` and `Quot.sound`, and
tests rejection of `sorryAx` and compiler-trust proof shortcuts.

Adversarial tests include malformed types, incomplete truth tables, dropped
or retargeted gates, dirty work, cached-writer substitution, additive tuple
collisions, stale source pins, omitted evidence and forged clock/resource
claims. Coordinated substitutions of both a circuit and its truth table,
cached-writer substitutions during receipt recomputation, and corrupted
Kraus operators also fail. The CLI rejects junk under Python optimization as well. The suite
does not download a simulator, require cloud storage, or infer a theorem from
a continuum fit.
