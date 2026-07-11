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
- A standalone, application-neutral proof package in
  `Proofs/ObservableNormalForms/`.  Its generic endpoint theorem is connected
  to the concrete local-repair interface by
  `Source/ObserverPatchHolography/Bridges/ObservableNormalForms.lean`.  This
  bridge characterizes the remaining #304 premise; it does not prove the
  declared physical boundary map injective.

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
    lake build                # build the library + Main executable
    ./.lake/build/bin/oph     # optional: run the entry point

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
