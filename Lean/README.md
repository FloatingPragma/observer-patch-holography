# Lean formalization workspace

This is the umbrella Lean 4 / Mathlib project for the repository. It contains
a sorry-free proof subset covering finite observer consensus, public records,
normal forms, coupling algebra, the screen/trichotomy arithmetic, and the
exact algebraic/compositional kernel of the corrected Einstein branch. Lean
checks those formal statements exactly. Continuum geometry, asymptotic tails,
physical identification, and existence of an Einstein-admissible source tower
remain explicit premises rather than proved facts.

## Layout

One Lake workspace, four libraries, one directory each:

```text
Lean/
├── ObserverPatchHolography.lean        umbrella module of the main library
├── ObserverPatchHolography/            main OPH library: carrier, repair,
│   ├── Bridges/                        consensus, coupling algebra, collar
│   └── EinsteinBranch/                 chain, Einstein-branch kernel
├── EventAlgebra.lean
├── EventAlgebra/                       neutral finite event algebras
│                                       (journal artifact, Mathlib-only)
├── Screen/                             OPHScreen library: icosahedral screen
│                                       arithmetic, A5 corpus, trichotomy
├── ObservableNormalForms/              standalone neutral submission package
│                                       (own lakefile; also built here)
├── docs/                               proof indices and application notes
├── Main.lean
├── lakefile.lean
├── lake-manifest.json
└── lean-toolchain
```

`ObserverPatchHolography.lean` is the public umbrella module.  It retains
Jonathan Hill's `OPH` development, re-exports the separate
`ObservableNormalForms` namespace, and imports a small bridge showing how the
generic boundary-identification theorem specializes to the concrete
local-repair interface.

The neutral submission project remains a single canonical source tree.  To
prepare its archive, zip the contents of `ObservableNormalForms/`; the outer
repository path is not part of the archive.

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
cd Lean/ObservableNormalForms
lake exe cache get
lake build
```

Its `.lake/packages` is a symlink into the umbrella project's package
directory, so the shared Mathlib checkout is reused when building locally;
continuous integration recreates it as a real directory.

## Documentation

- `docs/PROOF_INDEX.md` — proof-to-paper mapping and formalisation status
- `docs/LIBRARY_GUIDE.md` — scope and module guide for the main library
- `docs/EINSTEIN_BRANCH_INDEX.md` — Einstein-branch statement audit
- `docs/BRIDGE_BOUNDARY_INDEX.md` — cross-paper boundary map
- `docs/BOUNDARY_FIBER_APPLICATION.md` — #304 application note
- `ObservableNormalForms/README.md` and its `PROOF_INDEX.md` — manuscript
  coverage of the neutral submission package
