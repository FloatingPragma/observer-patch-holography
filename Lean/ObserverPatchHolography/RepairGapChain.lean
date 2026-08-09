/-
  RepairGapChain.lean — the OPH → Yang–Mills CONDITIONAL CHAIN, machine-checked.

  Self-contained (core Lean 4, NO Mathlib import), so it compiles + checks standalone
  with plain `lean` — no heavy Mathlib build required.

  Formalizes only a LOGICAL IMPLICATION SCHEMA associated with
    B. Müller, "Explaining the Yang–Mills Mass Gap with Observer-Patch Repair Dynamics":
      assumed finite repair-gap  +  assumed continuum spectral identification
        ⟹  positive continuum mass gap  Δ_YM ≥ c_* > 0   (Thm 11.1),
        and  Δ_YM = Δ_rep   (§12, exact gap accounting).

  HONEST SCOPE (score-function discipline — do NOT overclaim):
    * Only the CHAIN is proved here: "IF the certificate holds, the mass gap follows."
    * The certificate's deep analysis (Assumption 9.2 items Y1–Y6: Mosco convergence of
      the transfer forms, OS reconstruction, reflection positivity, the nontriviality
      variance floor) is encoded as HYPOTHESES (the `Certificate` structure), NOT proved.
      Discharging them is the genuine open frontier (constructive QFT); this file does not
      touch it and does not claim to.
    * The finite gap is a HYPOTHESIS here (`Lrep_gap`), not an unconditional
      result of this file.  In the current paper it is intended to be supplied
      by the noncommuting Dobrushin/approximate-tensorisation theorem.  The
      sibling Mathlib modules (`YangMillsLemma72`, `YangMillsProp81`,
      `YangMillsGap`) instead prove a legacy/special pairwise-commuting branch.
    * `spec_eq` is likewise a HYPOTHESIS.  This file proves neither continuum
      convergence nor compatibility between any finite-gap theorem and that
      spectral identification.
    * Consequently the legacy commuting theorem cannot be inserted into the
      current continuum receipt automatically.  A separate source proof of
      commutation and certificate compatibility would be required, and would
      define a different conditional branch from the current Dobrushin one.

  What this file buys: a machine check of the elementary implication from the
  explicitly bundled hypotheses.  It does not machine-check that the current
  paper establishes, or can consistently compose, those hypotheses.

  SCOPE: machine-checked only: "assumed positive finite gap + assumed equality
  of nonzero spectra ⇒ positive continuum gap".  Neither the current
  noncommuting finite-stage theorem nor the continuum receipt is proved here.
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

/-- A continuum spectral certificate together with an assumed finite gap,
    bundled as the HYPOTHESES on which the conditional theorem rests.
    `H` = continuum Yang–Mills Hamiltonian, `Lrep` = continuum repair generator. -/
structure Certificate (H Lrep : Op α) (c_star zero : α) : Prop where
  /-- Lemma 7.4 (uniform active-collar rate floor): the repair rate is strictly positive. -/
  c_star_pos : lt zero c_star
  /-- Assumed repair-generator gap `Δ_rep ≥ c_*`.  On the current paper's
      branch this is intended to come from the noncommuting Dobrushin argument;
      this field does not prove that theorem or transport any legacy
      commuting-projection result to the continuum. -/
  Lrep_gap : GapAtLeast Lrep c_star
  /-- Assumed exact spectral accounting: `H` and `Lrep` share their nonzero
      spectrum.  The structure records this compatibility; it does not derive
      it from a finite-stage theorem or a continuum limit. -/
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
