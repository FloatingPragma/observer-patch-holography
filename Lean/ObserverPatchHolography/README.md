# Observer-Patch Holography: Lean 4 Formalisation

Lean 4 / Mathlib formalisation effort for *Paradise as Fixed-Point Consensus*
(B. Müller, 2026; source in `paper/paradise_as_fixed_point_consensus.tex`),
with **Proposition 4.2** as the primary target.

## Scope

This project is an active Lean 4 / Mathlib formalisation and proof-audit
surface for the OPH consensus layer. Contents:

- Lake project with pinned `leanprover/lean4:v4.29.1` and Mathlib `v4.29.1`.
- An abstract-rewriting skeleton (Newman's lemma, normal-form uniqueness,
  descent termination, fixed-point zero-potential corollary) in
  `Source/ObserverPatchHolography/AbstractRewriting.lean`.
- A concrete OPH carrier layer in `Source/ObserverPatchHolography/Primitives.lean`
  for records, patch/interface data, observable overlap maps, the mismatch
  potential `Φ`, gauge equivalence, edge-consistency, and non-vacuity witnesses.
- Machine-checked proof-audit material for the consensus/reconstruction layer,
  including boundary-fiber observer uniqueness, commutation-based confluence,
  concrete countermodels separating confluence from observer-facing uniqueness,
  and axiom audits for the discharged reconstruction statements.
- Declared `sorry`-bearing signatures for the paper-incomplete
  asynchronous/transactional repair machinery: `localRepair`, `Repair`, and
  `repair_respects_gauge`.
- A sorry-free **#304 boundary-fiber carrier witness** in
  `Source/ObserverPatchHolography/Rule90.lean` (PR #385): the linear Rule 90 carrier
  discharges the `Hfib` binder of `boundary_fiber_observer_unique` on a proper
  information-set boundary, with a bad-boundary counterexample, a non-trivial
  gauge, and a local-repair no-go (`H1`–`H3` route only). A carrier-level
  witness; it does **not** advance the Prop 4.2 target. See `PROOF_INDEX.md`.
- A sorry-free **Part-A coupling-algebra layer** (13 lemmas, standard axioms
  only): `Source/ObserverPatchHolography/BridgeEquivalence.lean` (bridge
  count/tick equivalence, 5 lemmas),
  `Source/ObserverPatchHolography/CapacityFixedPoint.lean` (capacity
  fixed-point uniqueness schema, 4 lemmas), and
  `Source/ObserverPatchHolography/SeedPi.lean` (CAP-P seed statement,
  4 lemmas). These formalise the ALGEBRAIC layer of the coupling theorem and
  carry no physical-derivation content; the physical identities I1/I2 are
  outside the formalised set. Numeric interval enclosures stay in the Python
  certificates (`code/capacity_readback/`, `code/P_derivation/`); no
  floating-point numerics enter Lean. See `PROOF_INDEX.md`.
- A standalone, application-neutral proof package in
  `Proofs/ObservableNormalForms/`.  Its generic endpoint theorem is connected
  to the concrete local-repair interface by
  `Source/ObserverPatchHolography/Bridges/ObservableNormalForms.lean`.  This
  bridge characterizes the #304 application premise for arbitrary carriers;
  the premise itself is discharged on the declared domain by the
  `BoundaryFiber` module below.
- A sorry-free **#304 application theorem**
  (`Source/ObserverPatchHolography/BoundaryFiber.lean`): the verified
  rooted-tree packet-net domain (*Reality* Definition `def:tree-packet-domain`)
  as an `OPHCarrier` family (`TreePacketNet`), the declared physical
  boundary/sector map `B_OPH` as the concrete root-packet readback
  (`TreePacketNet.BOPH`), and the gate theorem
  `BOPH_injective_modulo_gauge`: consistent records with equal `B_OPH`
  readback are `gaugeEquiv`, for every net in the class. The generic
  hypotheses (`ObservationPreserving`, `CompleteFor`) are discharged for the
  domain's own tree repair, giving the endpoint payoff
  `BOPH_observerEndpointUnique` with no confluence input. Witnesses: the
  paper's exported four-vertex `ℤ₃`/`ℤ₂` instance with a proper gauge
  quotient, and a deficient-readback failure. Outside the declared class,
  boundary identifiability stays a named per-net premise (countermodels:
  `demoCarrier_Hfib_fails`, `rule90_Hfib_bad_fails`). See
  `BOUNDARY_FIBER_APPLICATION.md` for the hypothesis-by-hypothesis bridge.
- A sorry-free **finite event algebras** library (`Source/EventAlgebra/`,
  lake target `EventAlgebra`, 64 audited declarations, standard axioms
  only): events as Hermitian idempotents, states as positive trace-one
  matrices, Born weights (reality, nonnegativity, normalisation,
  additivity, complement bound, monotonicity), Lüders conditioning
  (state preservation, repeatability, idempotence, compatibility,
  classical restriction, fixed-point characterisation), the conditional
  expectation onto a commutative center (projector laws, state
  preservation, trace selfadjointness, Pythagoras, contractivity,
  uniqueness, compatibility with conditioning on central events), the
  expectation functional, and the Tsirelson bound `‖S‖ ≤ 2√2` in
  abstract unital C*-rings with a matrix instantiation. **This bundle is
  OPH-vocabulary-free by design** (namespace `EventAlgebra`, Mathlib-only
  imports, no repository vocabulary); it is the journal-neutral surface
  for submission; every lemma is tagged **algebra-only** or **consumes a
  tracial state** in its doc comment. Inventory in `PROOF_INDEX.md`
  ("Finite event algebras"); Mathlib friction log in
  `Source/EventAlgebra/MATHLIB_NOTES.md`. Not a Prop 4.2 / Def 4.1 item.

What is **not** yet present and is the target:

> **Proposition 4.2 (Fixed-point reading of reality).** Define the public
> world as the quotient-normal form `World = NF(x) / ∼_gauge`, where `NF(x)`
> is the terminal state reached by accepted repair and `∼_gauge` identifies
> hidden local presentations with the same declared observable overlap data
> (Definition 4.1). Then `World ∈ Fix(Repair)` and `Repair(World) = World`.
> When OPH confluence and completeness conditions hold, this terminal public
> state is independent of update schedule on the physical quotient.

A theorem-grade Lean statement matching Prop 4.2 requires:

- `World`, `Records`, `Repair`, `Patch`, `Obs`, and the local accepted
  repair-step relation as concrete Lean structures matching paper
  Definition 4.1 and the OPH preliminaries. Note: several names are TeX
  macros in *Paradise* (lines 28–31) whose structural content lives in the
  companion paper *Reality as a Consensus Protocol*; i.e. the target is
  **paper-incomplete** as well as Lean-incomplete.
- The mismatch potential `Φ : Records → NNReal` with the paper's concrete
  formula `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))`
  (*Paradise* line 300).
- The gauge equivalence `∼_gauge` on declared observable overlap data,
  proved an equivalence relation AND a `Repair`-congruence (required by
  Prop 4.2 sentence 2's "on the physical quotient" clause).
- `NF` defined as the terminal state of accepted repair built from local
  recovery moves (line 297), with those local repairs composed under the
  asynchronous update schedules used by the consensus companion.
- OPH-specific Lyapunov descent/termination, `Confluence`, and
  `Completeness` obligations as Lean definitions (Prop 4.2 hypothesis,
  line 326, with details supplied by the consensus companion).
- Schedule independence on the physical quotient, transferring the
  abstract-rewriting confluence result to the structured OPH setting.

The full quotient-normal-form theorem is not yet fully formalised as a single
`World` statement. The abstract-rewriting module is the generic skeleton, while
`Primitives` discharges concrete carrier and reconstruction subclaims and keeps
the remaining asynchronous repair obligations explicit so they cannot be
silently elided. See `PROOF_INDEX.md` for the proof-to-paper map and completion
tracker.

## Building

    cd Lean
    lake exe cache get        # fetch pre-built Mathlib oleans
    lake build                # build all three proof libraries
                              # (ObservableNormalForms, ObserverPatchHolography,
                              #  EventAlgebra)

The `Main` console entry point is optional and not part of the proof receipt;
build it separately with `lake build oph:exe` if needed.

Lean CI is manual. Use GitHub Actions, `Lean CI`, `Run workflow` when the
Lean formalisation changes need a hosted check. The workflow allows exactly
the 3 intentional `sorry` warnings in `Primitives.lean` and fails if new proof
debt appears elsewhere or the count changes.

## Provenance

- PR #299 (closed 2026-05-18 unmerged) shipped the abstract-rewriting
  skeleton as a claimed Proposition 4.2 formalisation. Audit verdict: the
  proofs are sorry-free but generic; they do not reach OPH-specific
  structure.
- This scaffold ports those proofs into a properly-built Lake project,
  applies accurate labels, and lays out the gap to be closed.
- Jonathan Hill contributed the substantive Lean formalisation and proof-audit
  work that closed concrete carrier/reconstruction subclaims, added non-vacuity
  witnesses, separated confluence from observer-facing uniqueness, and exposed
  the remaining asynchronous repair obligations explicitly.
- Coordination: "OPH LEAN Proofs" working group (Bernhard Mueller, Ben
  Cassie, Dula, Jonathan Hill). Cross-audit between auditors is required before
  PRs are merged.

## License And Patent Policy

This formalisation surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE) and [OPH Open Use And Anti-Patent Covenant](../../PATENTS.md).
