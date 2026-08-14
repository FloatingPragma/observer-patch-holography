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
