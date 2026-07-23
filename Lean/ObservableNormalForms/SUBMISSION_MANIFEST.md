# Submission artifact manifest

The directory `ObservableNormalForms/` is the archive root. Include exactly
the following tracked files:

```text
.gitignore
BUILD_RECEIPT.md
HASHES.sha256
ObservableNormalForms.lean
ObservableNormalForms/AxiomAudit.lean
ObservableNormalForms/ConditionalResampling.lean
ObservableNormalForms/Exact.lean
ObservableNormalForms/Examples/Rule90.lean
ObservableNormalForms/Functional.lean
ObservableNormalForms/ObserverConfluence.lean
ObservableNormalForms/Refinement.lean
ObservableNormalForms/Repair.lean
ObservableNormalForms/Stability.lean
ObservableNormalForms/Stochastic.lean
PROOF_INDEX.md
README.md
SUBMISSION_MANIFEST.md
lake-manifest.json
lakefile.lean
lean-toolchain
```

Do not include `.lake/`, editor metadata, generated object files, or build
logs.

## Reproduction recipe

```sh
lake exe cache get
lake build 2>&1 | tee lake-build.log
rg -n '^\s*(sorry|admit)\b|:=\s*(sorry|admit)\b' \
  ObservableNormalForms --glob '*.lean'
rg -n 'sorryAx' lake-build.log
```

Expected result: `lake build` succeeds; the source audit finds no admissions;
the theorem-level axiom output in `lake-build.log` contains no `sorryAx`.

Before creating a public archive, regenerate the file hashes with:

```sh
find . -type f \
  ! -path './.lake/*' \
  ! -name 'lake-build.log' \
  ! -name 'HASHES.sha256' \
  -print0 | sort -z | xargs -0 sha256sum > HASHES.sha256
```

Add `HASHES.sha256` to the archive after verifying it with
`sha256sum -c HASHES.sha256`.
