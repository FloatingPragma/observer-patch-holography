# Lean formalization workspace

This is the umbrella Lean 4 / Mathlib project for the repository. The paper
stack cites a sorry-free 111-theorem subset covering finite observer consensus,
public records, normal forms, and related coupling algebra. Lean checks those
formal statements exactly. The continuum gravity branch and physical
screen-to-matter carriers remain paper-level conditional constructions.

## Layout

```text
Lean/
├── ObserverPatchHolography/
│   ├── Source/                         concrete carrier and dynamics
│   ├── Proofs/ObservableNormalForms/   standalone neutral paper artifact
│   ├── README.md
│   └── PROOF_INDEX.md
├── Main.lean
├── lakefile.lean
├── lake-manifest.json
└── lean-toolchain
```

`ObserverPatchHolography/Source/ObserverPatchHolography.lean` is the public
umbrella module.  It retains Jonathan Hill's `OPH` development, re-exports the
separate `ObservableNormalForms` namespace, and imports a small bridge showing
how the generic boundary-identification theorem specializes to the concrete
local-repair interface.

The neutral submission project remains a single canonical source tree.  To
prepare its archive, zip the contents of
`ObserverPatchHolography/Proofs/ObservableNormalForms/`; the outer repository
path is not part of the archive.

## Build

```sh
cd Lean
lake exe cache get
lake build
```

The proof receipt is the library build above.  The tiny console entry point is
optional and requires a native executable build (`lake build oph:exe`).

The neutral submission artifact also builds independently:

```sh
cd Lean/ObserverPatchHolography/Proofs/ObservableNormalForms
lake build
```

See `ObserverPatchHolography/README.md` and
`ObserverPatchHolography/PROOF_INDEX.md` for the concrete proof status, and the
nested proof package's own `README.md` and `PROOF_INDEX.md` for its manuscript
coverage.
