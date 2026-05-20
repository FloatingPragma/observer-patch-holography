# Observer-Patch Holography — Lean 4 Formalisation

Lean 4 / Mathlib formalisation effort for *Paradise as Fixed-Point Consensus*
(B. Müller, 2026; source in `paper/paradise_as_fixed_point_consensus.tex`),
with **Proposition 4.2** as the primary target.

## Scope and current state

This project is **early-stage scaffolding**. What is currently present:

- Lake project with pinned `leanprover/lean4:v4.29.1` and Mathlib `v4.29.1`.
- An abstract-rewriting skeleton (Newman's lemma, normal-form uniqueness,
  descent termination, fixed-point zero-potential corollary) in
  `ObserverPatchHolography/AbstractRewriting.lean`.
- Sorry-bearing primitive signatures in
  `ObserverPatchHolography/Primitives.lean` (Records, Patch, Obs, Repair,
  local accepted repair steps, Φ, gauge equivalence, OPH-Confluence, and
  repair-completeness).

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
  companion paper *Reality as a Consensus Protocol* — i.e. the target is
  **paper-incomplete**, not just Lean-incomplete.
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
  `Completeness` obligations as Lean definitions (Prop 4.2 hypothesis —
  line 326, with details supplied by the consensus companion).
- Schedule independence on the physical quotient, transferring the
  abstract-rewriting confluence result to the structured OPH setting.

None of this is currently formalised. The abstract-rewriting module is the
preliminary skeleton, not a statement about `World`. The sorry-bearing
`Primitives` module declares the obligations explicitly so they cannot be
silently elided. See `PROOF_INDEX.md` for the proof-to-paper map and
completion tracker.

## Building

    cd LEAN
    lake exe cache get        # fetch pre-built Mathlib oleans
    lake build                # build the library + Main executable
    ./.lake/build/bin/oph     # optional: run the entry point

CI runs on every push and pull request — see
`.github/workflows/lean-ci.yml`. The workflow currently allows exactly the
10 intentional `sorry` warnings in `Primitives.lean` and fails if new proof
debt appears elsewhere or the count changes.

## Provenance

- PR #299 (closed 2026-05-18 unmerged) shipped the abstract-rewriting
  skeleton as a claimed Proposition 4.2 formalisation. Audit verdict: the
  proofs are sorry-free but generic; they do not reach OPH-specific
  structure.
- This scaffold ports those proofs into a properly-built Lake project,
  relabels them honestly, and lays out the gap to be closed.
- Coordination: "OPH LEAN Proofs" working group (Bernhard Mueller, Ben
  Cassie, Dula). Cross-audit between auditors is required before PRs are
  merged.
