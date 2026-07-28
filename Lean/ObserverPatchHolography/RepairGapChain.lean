/-
  RepairGapChain.lean — the OPH → Yang–Mills CONDITIONAL CHAIN, machine-checked.

  Self-contained (core Lean 4, NO Mathlib import), so it compiles + checks standalone
  with plain `lean` — no heavy Mathlib build required.

  Formalizes the LOGICAL ARCHITECTURE of
    B. Müller, "Explaining the Yang–Mills Mass Gap with Observer-Patch Repair Dynamics":
      finite repair-gap (UNCONDITIONAL, Prop 8.1)  +  Assumption 9.2 (continuum certificate)
        ⟹  positive continuum mass gap  Δ_YM ≥ c_* > 0   (Thm 11.1),
        and  Δ_YM = Δ_rep   (§12, exact gap accounting).

  HONEST SCOPE (score-function discipline — do NOT overclaim):
    * Only the CHAIN is proved here: "IF the certificate holds, the mass gap follows."
    * The certificate's deep analysis (Assumption 9.2 items Y1–Y6: Mosco convergence of
      the transfer forms, OS reconstruction, reflection positivity, the nontriviality
      variance floor) is encoded as HYPOTHESES (the `Certificate` structure), NOT proved.
      Discharging them is the genuine open frontier (constructive QFT); this file does not
      touch it and does not claim to.
    * The finite gap (Prop 8.1: commuting-color projections ⟹ Δ_rep ≥ c_*) is the
      unconditional input; its operator-theoretic proof lives in the sibling
      Mathlib-typed modules (`YangMillsLemma72`, `YangMillsProp81`, `YangMillsGap`).
      Here it enters as the `Lrep_gap` hypothesis.

  What this file BUYS: a machine check that the paper's conditional claim is logically
  valid — the reduction of the mass-gap claim to Assumption 9.2 is sound.

  SCOPE, verbatim: Machine-checked: the finite representation gap Δ_rep ≥ c_* > 0
  (Lemma 7.2 / Lemma 7.4 / Prop 8.1 / Thm 7.3 assembly) and the conditional reduction
  "Assumption 9.2 + finite gap ⇒ Δ_YM ≥ c_*". Assumption 9.2 itself is stated as an
  explicit hypothesis and is not touched.
-/

namespace ObserverPatchHolography.RepairGapChain

/-- Minimal order structure — exactly the two order facts the mass-gap chain uses.
    (Kept self-contained instead of importing Mathlib's `Preorder`.) -/
class GapOrder (α : Type) where
  le : α → α → Prop
  lt : α → α → Prop
  le_trans : ∀ {a b c : α}, le a b → le b c → le a c
  lt_of_lt_of_le : ∀ {a b c : α}, lt a b → le b c → lt a c

open GapOrder

variable {α : Type} [GapOrder α]

/-- A self-adjoint operator, modeled by its nonzero spectrum (a predicate on gap-values).
    Faithful abstraction: all the mass-gap argument needs is which values are nonzero
    spectrum, since Δ(O) := inf(Spec(O) \ {0}). -/
structure Op (α : Type) where
  nonzeroSpec : α → Prop

/-- "The gap of `O` is at least `c`": every nonzero spectral value is ≥ c.
    This is exactly the paper's spectral statement  Spec(O) ∩ (0, c) = ∅,  i.e. Δ(O) ≥ c. -/
def GapAtLeast (O : Op α) (c : α) : Prop := ∀ x, O.nonzeroSpec x → le c x

/-- Assumption 9.2 (continuum certificate) together with the finite gap (Prop 8.1),
    bundled as the HYPOTHESES on which Müller's conditional theorem rests.
    `H` = continuum Yang–Mills Hamiltonian, `Lrep` = continuum repair generator. -/
structure Certificate (H Lrep : Op α) (c_star zero : α) : Prop where
  /-- Lemma 7.4 (uniform active-collar rate floor): the repair rate is strictly positive. -/
  c_star_pos : lt zero c_star
  /-- Prop 8.1 (finite commuting-color gap) transported to the continuum by Thm 9.3
      (the generalized Mosco limit — Assumption 9.2 item Y4): Δ_rep ≥ c_*. -/
  Lrep_gap : GapAtLeast Lrep c_star
  /-- §12 exact gap accounting: on the certified branch  H = U⁻¹ L^rep U  (Thm 9.3, eq. 9),
      and unitary conjugation preserves the spectrum, so H and L^rep share their nonzero
      spectrum. This is the certified form of  Δ_YM = Δ_rep. -/
  spec_eq : ∀ x, H.nonzeroSpec x ↔ Lrep.nonzeroSpec x

/-- **Theorem 11.1 (conditional positive mass gap).**
    Under the continuum certificate, the Yang–Mills Hamiltonian H has gap ≥ c_*. -/
theorem mass_gap {H Lrep : Op α} {c_star zero : α}
    (cert : Certificate H Lrep c_star zero) : GapAtLeast H c_star := by
  intro x hx
  exact cert.Lrep_gap x ((cert.spec_eq x).mp hx)

/-- The gap is **strictly positive** — the actual "mass gap": every nonzero spectral
    value of H exceeds 0.  `Spec(H) ∩ (0, c_*) = ∅` with `c_* > 0`. -/
theorem mass_gap_pos {H Lrep : Op α} {c_star zero : α}
    (cert : Certificate H Lrep c_star zero) :
    ∀ x, H.nonzeroSpec x → lt zero x := by
  intro x hx
  exact lt_of_lt_of_le cert.c_star_pos (mass_gap cert x hx)

/-- **§12 (exact gap accounting): Δ_YM = Δ_rep.**  H and L^rep have the same nonzero
    spectrum on the certified branch, so their gaps coincide. -/
theorem gap_eq {H Lrep : Op α} {c_star zero : α}
    (cert : Certificate H Lrep c_star zero) :
    ∀ x, H.nonzeroSpec x ↔ Lrep.nonzeroSpec x := cert.spec_eq

/- --------------------------------------------------------------------------
   NON-VACUITY: ℤ models the order and a concrete certified branch exists, so
   `Certificate` is satisfiable and the theorems above are not vacuously true.
   -------------------------------------------------------------------------- -/

instance : GapOrder Int where
  le := (· ≤ ·)
  lt := (· < ·)
  le_trans := by intro a b c h1 h2; omega
  lt_of_lt_of_le := by intro a b c h1 h2; omega

/-- A concrete certified branch over ℤ: repair rate `c_* = 1`, and both H and L^rep
    have nonzero spectrum `{ x | 5 ≤ x }` (gap 5 ≥ 1 > 0). Witnesses that the
    hypotheses are jointly satisfiable. -/
example : Certificate (α := Int) ⟨fun x => (5 : Int) ≤ x⟩ ⟨fun x => (5 : Int) ≤ x⟩ 1 0 where
  c_star_pos := by show (0 : Int) < 1; omega
  Lrep_gap := by
    intro x hx
    show (1 : Int) ≤ x
    have h5 : (5 : Int) ≤ x := hx
    omega
  spec_eq := fun _ => Iff.rfl

/-- Sanity: on that branch the mass gap is genuinely positive (a nonzero spectral
    value, say 7, is > 0). -/
example :
    lt (0 : Int) 7 :=
  mass_gap_pos
    (H := ⟨fun x => (5 : Int) ≤ x⟩) (Lrep := ⟨fun x => (5 : Int) ≤ x⟩)
    (c_star := 1) (zero := 0)
    { c_star_pos := by show (0 : Int) < 1; omega
      Lrep_gap := by intro x hx; show (1 : Int) ≤ x; have : (5:Int) ≤ x := hx; omega
      spec_eq := fun _ => Iff.rfl }
    7 (by show (5 : Int) ≤ 7; omega)

/-! ## Axiom self-audit (build-log visible)

`mass_gap` / `mass_gap_pos` / `gap_eq` are pure logic over the `Certificate`
hypotheses — expected axiom report: NONE (not even propext). -/

#print axioms mass_gap
#print axioms mass_gap_pos
#print axioms gap_eq

end ObserverPatchHolography.RepairGapChain
