# Build receipt

Verified on 2026-07-11 (arm64 macOS):

```text
Lean 4.29.1, commit f72c35b3f637c8c6571d353742168ab66cc22c00
Lake 5.0.0-src+f72c35b
Mathlib input revision v4.29.1
Mathlib commit 5e932f97dd25535344f80f9dd8da3aab83df0fe6
```

Both integration modes succeeded:

```text
cd <repository>/Lean
lake build ObservableNormalForms
Build completed successfully (8259 jobs).

cd <artifact-root>
lake build
Build completed successfully (8259 jobs).
```

The full parent project also succeeded:

```text
cd <repository>/Lean
lake build
Build completed successfully (8265 jobs).
```

That full build also replays the repository's concrete target. The former
`Primitives.lean` admissions (`localRepair`, `Repair`, and
`repair_respects_gauge`) are discharged. The `ObservableNormalForms`
sources remain independently admission-free, and their axiom audit reports no
`sorryAx`.

`ObservableNormalForms/AxiomAudit.lean` was included in both builds. Its
theorem-level reports contain only the standard axioms `propext`,
`Classical.choice`, and `Quot.sound` where required; several finite and exact
theorems report no axioms. No theorem reports `sorryAx`.

The source admission audit

```sh
rg -n '^\s*(sorry|admit)\b|:=\s*(sorry|admit)\b' \
  ObservableNormalForms --glob '*.lean'
```

returns no matches.
