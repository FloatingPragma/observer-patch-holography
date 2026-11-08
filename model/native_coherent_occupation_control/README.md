# Coherent occupation controls

The twelve-port source supplies exact rotations that preserve occupation,
local phases and an entangling conditional phase. Combining two rotations
with one native boundary repair creates a heralded Bell resource with
success above 0.4999, retaining every failure and record. The finite evaluator
checks the addressed source handoff without projecting or normalizing its
state. These explicitly declared coherent services preserve cell charge;
preparing a growing field reference requires a further coherent coupling.

The [source-control proof](DERIVED_CONTROL.md),
[phase construction](DIAGONAL_CONTROL.md) and
[complete receipt](results.json) give the definitions and bounds.
[PROOF.md](PROOF.md) establishes the exact native parity operation, bounded
ideal retries and cell-charge obstruction. The
[literal-control counterexample](LITERAL_CONTROL_AUDIT.md) preserves the
fixed-number handoff limitation of the 22-pulse beam.

Run from the repository root:

```sh
python -m model.native_coherent_occupation_control.packet --check
python -m pytest -q model/native_coherent_occupation_control
```

The source and test hashes are checked before and after each owner replay.
Receipts retain exact analytic bounds; floating successful-array values are
checked diagnostics rather than cross-platform equality commitments.
