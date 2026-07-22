# Lean formalization workspace

This is the umbrella Lean 4 / Mathlib project for the repository. It contains
a sorry-free proof subset covering finite observer consensus, public records,
normal forms, coupling algebra, and the exact algebraic/compositional kernel of
the corrected Einstein branch. Lean checks those formal statements exactly.
Continuum geometry, asymptotic tails, physical identification, and existence
of an Einstein-admissible source tower remain explicit premises rather than
proved facts.

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
