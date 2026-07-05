import Mathlib

/-!
# Abstract Rewriting: Newman's Lemma and Fixed-Point Corollary

The abstract rewriting skeleton used by *Paradise as Fixed-Point Consensus*:

* Newman's lemma (terminating + locally confluent → confluent)
* Uniqueness of normal forms (schedule independence on a generic ARS)
* A deterministic-operator specialisation with descent potential
* A corollary that any fixed point of a repair operator has zero
  mismatch potential (motivational only — see doc-comment on
  `fixedPt_zero_potential`)

**Scope warning.** These are generic abstract-rewriting results stated over
an opaque type `X` and an opaque relation `r : X → X → Prop`. They are
**preliminary**: they do not commit to OPH-specific structure (`World`,
gauge quotient, observer patches, declared observable overlap data, local
recovery moves). The theorem-grade formalisation of Proposition 4.2 lives
downstream and is tracked in `PROOF_INDEX.md`.

Source: PR #299 (Aristotle output, verified sorry-free) — relabelled per the
2026-05-18 audit verdict that this is a skeleton, not a Prop 4.2 statement.
-/

open Function Relation

namespace OPH.AbstractRewriting

variable {X : Type*}

/-! ## Abstract Rewriting System -/

/-- A relation `r` is terminating if there are no infinite forward chains
    `x₀ → x₁ → x₂ → ⋯`. Formally, `WellFounded (flip r)` where `r x y`
    means "`x` rewrites to `y`". -/
def Terminating (r : X → X → Prop) : Prop := WellFounded (flip r)

/-- A relation is locally confluent if any two one-step rewrites from the
    same source can be joined by multi-step rewrites. -/
def LocallyConfluent (r : X → X → Prop) : Prop :=
  ∀ x y z, r x y → r x z → ∃ w, ReflTransGen r y w ∧ ReflTransGen r z w

/-- A relation is confluent (Church–Rosser) if any two multi-step rewrites
    from the same source can be joined. -/
def Confluent (r : X → X → Prop) : Prop :=
  ∀ x y z, ReflTransGen r x y → ReflTransGen r x z →
    ∃ w, ReflTransGen r y w ∧ ReflTransGen r z w

/-- A normal form is an element with no successors under `r`. -/
def IsNormalForm (r : X → X → Prop) (x : X) : Prop := ∀ y, ¬r x y

/-- `x` reduces to normal form `y`: `x →* y` and `y` is a normal form. -/
def ReducesToNF (r : X → X → Prop) (x y : X) : Prop :=
  ReflTransGen r x y ∧ IsNormalForm r y

/-! ## Newman's Lemma -/

/-- **Newman's Lemma.** A terminating, locally confluent relation is
    confluent. -/
theorem newman_lemma (r : X → X → Prop)
    (hterm : Terminating r) (hlc : LocallyConfluent r) :
    Confluent r := by
      have h_wf_ind : ∀ x, ∀ y z, ReflTransGen r x y → ReflTransGen r x z → ∃ w, ReflTransGen r y w ∧ ReflTransGen r z w := by
        intro x;
        induction' x using hterm.induction with x ih;
        intro y z hy hz
        by_cases hxy : y = x;
        · aesop;
        · by_cases hxz : z = x;
          · exact ⟨ y, by aesop ⟩;
          · obtain ⟨y', hy'⟩ : ∃ y', r x y' ∧ ReflTransGen r y' y := by
              have := hy.cases_head;
              tauto
            obtain ⟨z', hz'⟩ : ∃ z', r x z' ∧ ReflTransGen r z' z := by
              have := hz.cases_head; aesop;
            obtain ⟨ w, hw₁, hw₂ ⟩ := hlc x y' z' hy'.1 hz'.1;
            obtain ⟨ u, hu₁, hu₂ ⟩ := ih y' hy'.1 y w hy'.2 hw₁;
            obtain ⟨ v, hv₁, hv₂ ⟩ := ih z' hz'.1 z u hz'.2 ( hw₂.trans hu₂ );
            exact ⟨ v, hu₁.trans hv₂, hv₁ ⟩;
      exact h_wf_ind

/-- Confluence implies uniqueness of normal forms. -/
theorem unique_normal_form (r : X → X → Prop) (hconf : Confluent r)
    {x y z : X} (hy : ReducesToNF r x y) (hz : ReducesToNF r x z) : y = z := by
      obtain ⟨ w, hw₁, hw₂ ⟩ := hconf x y z hy.1 hz.1;
      have hy_normal : ∀ w, ReflTransGen r y w → y = w := by
        intro w hw; induction hw <;> simp_all +decide [ ReducesToNF ] ;
        exact False.elim ( hy.2 _ ‹_› );
      have hz_normal : ∀ w, ReflTransGen r z w → z = w := by
        intro w hw; induction hw <;> simp_all +decide [ ReducesToNF ] ;
        exact False.elim ( hz.2 _ ‹_› );
      rw [ hy_normal w hw₁, hz_normal w hw₂ ]

/-- Schedule independence (generic ARS): terminating + locally confluent
    implies normal forms are unique. -/
theorem newman_unique_nf (r : X → X → Prop)
    (hterm : Terminating r) (hlc : LocallyConfluent r)
    {x y z : X} (hy : ReducesToNF r x y) (hz : ReducesToNF r x z) : y = z :=
  unique_normal_form r (newman_lemma r hterm hlc) hy hz

/-! ## Deterministic Repair Operator -/

section DeterministicRepair

variable (T : X → X)

/-- The one-step rewriting relation induced by `T`: `x` steps to `T x`
    precisely when `x` is not a fixed point. -/
def stepRel (T : X → X) (x y : X) : Prop := y = T x ∧ T x ≠ x

/-- A normal form of `stepRel T` is a fixed point of `T`. -/
theorem nf_is_fixedPt (w : X) (hnf : IsNormalForm (stepRel T) w) :
    IsFixedPt T w := by
      exact Classical.not_not.1 fun h => hnf _ ⟨ rfl, h ⟩

/-- A descent potential `Φ : X → ℕ` that strictly decreases on
    non-fixed-points makes the induced rewriting terminating. -/
theorem descent_terminating
    (Φ : X → ℕ)
    (hdesc : ∀ x, T x ≠ x → Φ (T x) < Φ x) :
    Terminating (stepRel T) := by
      refine' ⟨ fun x => _ ⟩;
      have h_wf : WellFounded (fun x y : X => stepRel T y x) := by
        refine' ⟨ fun x => _ ⟩;
        induction' n : Φ x using Nat.strong_induction_on with n ih generalizing x;
        refine' Acc.intro x fun y hy => _;
        cases hy ; aesop;
      exact h_wf.apply x

/-- A deterministic operator is trivially locally confluent. -/
theorem deterministic_locally_confluent :
    LocallyConfluent (stepRel T) := by
      intro x y z hxy hxz;
      cases hxy ; cases hxz ; aesop

/-- Combined: for a terminating deterministic repair operator the normal
    form is a fixed point of `T` and is unique (schedule-independent). -/
theorem deterministic_full
    (Φ : X → ℕ)
    (hdesc : ∀ x, T x ≠ x → Φ (T x) < Φ x)
    {x y z : X}
    (hy : ReducesToNF (stepRel T) x y)
    (hz : ReducesToNF (stepRel T) x z) :
    y = z ∧ IsFixedPt T y :=
  ⟨newman_unique_nf _ (descent_terminating T Φ hdesc) (deterministic_locally_confluent T) hy hz,
   nf_is_fixedPt T y hy.2⟩

end DeterministicRepair

/-! ## Corollary: fixed points have zero mismatch -/

/-- If positive mismatch potential always triggers repair, any fixed point
    has `Φ = 0`. **Motivational only, no formal paper anchor.** The
    *Paradise* passage at lines 330–334 is prose (not a labelled
    proposition/corollary), and the paper makes no formal `Φ(W) = 0`
    claim. This theorem is abstract-rewriting decoration; it does not
    translate any statement of the paper. Kept in the preliminary skeleton
    on that understanding. -/
theorem fixedPt_zero_potential
    (T : X → X) (Φ : X → NNReal)
    (hrepair : ∀ x, Φ x > 0 → T x ≠ x)
    (w : X) (hfp : IsFixedPt T w) :
    Φ w = 0 := by
      exact le_antisymm ( le_of_not_gt fun h => hrepair w h hfp ) ( zero_le _ )

end OPH.AbstractRewriting
