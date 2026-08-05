# Lean formalization workspace

This is the umbrella Lean 4 / Mathlib project for the repository. It contains
a sorry-free proof subset covering finite observer consensus, public records,
normal forms, coupling algebra, the screen/trichotomy arithmetic, and the
exact algebraic/compositional kernel of the typed Einstein branch. It also
checks the finite de Sitter capacity-transfer identities, the eigenvalue signs
of the declared analytic Hessian action, and the pure-de-Sitter shock
normalization. The source contains no admitted proofs. Twenty-three finite
proofs use `native_decide`, whose generated native-code evaluation axioms are
tracked by `tools/check_lean_native_decide_inventory.py`; the other audited
proofs state their axiom dependencies in their module reports. Continuum
geometry, asymptotic tails, physical identification, and existence of an
Einstein-admissible source tower are explicit premises rather than proved
facts.

## Layout

One Lake workspace, seven Lean libraries across their source directories:

```text
Lean/
├── ObserverPatchHolography.lean        umbrella module of the main library
├── ObserverPatchHolography/            main OPH library: carrier, repair,
│   ├── Bridges/                        consensus, coupling algebra, collar
│   └── EinsteinBranch/                 chain, Einstein-branch kernel
├── EventAlgebra.lean
├── EventAlgebra/                       neutral finite event algebras
│                                       (journal artifact, Mathlib-only)
├── Thermodynamics/                     finite repair thermodynamics,
│                                       stationary/nonreversible H-theorem
│                                       interfaces, reversible transport, graph diffusion,
│                                       and Einstein first-law premise links
├── Screen/                             OPHScreen library: icosahedral screen
│                                       arithmetic, A5 corpus, trichotomy
├── Dynamics.lean
├── Dynamics/                           OPHConstruction reusable dynamics
│                                       interfaces
├── InformationProjection.lean
├── InformationProjection/              conditional finite-history projection
├── Time/                               time/order type ledger and explicit
│                                       realization-map boundary
├── Tower/                              timeless consensus-tower interface and
│                                       constant finite-fiber adaptor
├── Locality.lean
├── ObserverPatchHolography/Locality/    fixed-word locality helpers
├── Variational.lean
├── Variational/                         scalar variational helpers and bridge
│                                       obstruction
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

- `docs/PROOF_INDEX.md`: proof-to-paper mapping and formalisation status
- `docs/LIBRARY_GUIDE.md`: scope and module guide for the main library
- `docs/EINSTEIN_BRANCH_INDEX.md`: Einstein-branch statement audit
- `docs/B4_LOCALITY_BOUNDARY.md`: fixed-word locality and physical boundary
- `docs/B5_WARD_BRIDGE.md`: finite continuity and Ward-premise boundary
- `docs/B7_HISTORY_BRIDGE.md`: conditional history helpers and interface no-go
- `docs/B8_TRANSPORT_KERNEL.md`: finite Green--Kubo and graph-transport boundary
- `docs/A1_TIME_ORDER_LEDGER.md`: distinct time/order types, explicit bridge API,
  and affine clock-gauge boundary
- `docs/A3_CONSENSUS_TOWER.md`: directed finite tower interface, constant
  projective-partition adaptor, and E1/E2 wiring boundary
- `docs/B1_PUBLIC_RECORD_ALGEBRA.md`: exact active-label record algebra,
  sharp no-cloning theorem, and mixed-state adapter boundary
- `docs/B2_PUBLICIZATION_DYNAMICS.md`: normalized Kraus data, solvable
  publicization semigroup, literal bounded-operator exponential, fixed algebra,
  and physical-channel boundary
- `docs/B3_PUBLIC_PRIVATE_DYNAMICS.md`: stochastic public maps, the continuous
  permutation-flow obstruction, full-private-block innerness, and the open
  public-automorphism classification, central-block, and converse-generator
  boundaries
- `docs/BRIDGE_BOUNDARY_INDEX.md`: cross-paper boundary map
- `docs/BOUNDARY_FIBER_APPLICATION.md`: #304 application note
- `ObservableNormalForms/README.md` and its `PROOF_INDEX.md`: manuscript
  coverage of the neutral submission package
