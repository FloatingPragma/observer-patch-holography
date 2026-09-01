# Build receipt

Verified on 2026-08-14 (x86_64 Windows) from OPH base
`6c2a96c24a82f3ab73fd16f5dc75abad3505b745`:

```text
Lean 4.29.1, commit f72c35b3f637c8c6571d353742168ab66cc22c00
Lake 5.0.0-src+f72c35b
Mathlib input revision v4.29.1
Mathlib commit 5e932f97dd25535344f80f9dd8da3aab83df0fe6
```

The standalone artifact build succeeded:

```text
cd <short-path-exact-copy>/Lean/ObservableNormalForms
lake build
Build completed successfully (8264 jobs).
```

The full parent project also succeeded:

```text
cd <short-path-exact-copy>/Lean
lake build
Build completed successfully (8539 jobs).
```

Windows could not materialize all generated Lean filenames below the
repository worktree path because that absolute path exceeds the platform path
limit. The two
builds therefore used short-path copies. Before recording this receipt, SHA-256
comparison found 0 differences across all 26 standalone artifact files and 0
differences across all 292 parent Lean inputs (including toolchain and manifest
files, excluding generated `.lake` trees).

`ObservableNormalForms/AxiomAudit.lean` was also executed directly after the
builds and exited successfully. Its theorem-level reports contain only the
standard axioms `propext`, `Classical.choice`, and `Quot.sound` where required;
several finite and exact theorems report no axioms. No theorem reports
`sorryAx`.

The source admission audit

```sh
rg -n '^\s*(sorry|admit)\b|:=\s*(sorry|admit)\b' \
  ObservableNormalForms --glob '*.lean'
```

returns no matches.

## Refresh receipt: 2026-09-01

The archive boundary was revalidated while developing from frozen stacked OPH
base `bbe16c0b2697e6a4beb64454832af81f1a89822e`. The library root already
imported `ObservableNormalForms/SchedulerClassObstructions.lean`; the manifest
and README now list that tracked standalone module explicitly. Parent-workspace
OPH/Tower theorem rows were removed from this standalone proof index and remain
in the parent `Lean/docs/PROOF_INDEX.md`.

With Lean 4.29.1 and the pinned Mathlib revision, the standalone command

```text
lake build
Build completed successfully (8265 jobs).
```

exited 0. The fresh parent-workspace build also exited 0 after 8,668 jobs. The
standalone axiom audit reported only `propext`, `Classical.choice`, and
`Quot.sound` where required and no `sorryAx`. After the boundary repair, the
documented hash-generation recipe was rerun and `sha256sum -c HASHES.sha256`
verified every listed file.
