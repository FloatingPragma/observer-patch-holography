import Mathlib

/-!
# Coupling-Algebra Part A: Capacity Fixed-Point Schema

The abstract uniqueness schema behind the capacity readback certificates:

* `averagingMap_unique_fixedPt` — the averaging map
  `C(x) = (1 − λ)·x + λ·K` with `0 < λ ≤ 1` and constant target `K` has
  exactly one fixed point, namely `K`. This pins the averaging-map argument
  used by the capacity certificates.
* `banach_unique_fixedPt` — a `ContractingWith`-map on a nonempty complete
  metric space has exactly one fixed point. Cited from Mathlib
  (`ContractingWith.fixedPoint` and its uniqueness lemma), not reproved.
* `capacity_Icc_unique_fixedPt` — the Banach schema applied to the capacity
  map's interval domain: a self-map of `Icc a b` that is Lipschitz on
  `Icc a b` with constant `K < 1` has exactly one fixed point in `Icc a b`.

**Scope warning.** This module formalises the ALGEBRAIC layer only: the
uniqueness schema is stated abstractly, with no floating-point numerics in
Lean. The numeric interval enclosures (which specific `a`, `b`, `K` the
capacity map satisfies) stay in the Python certificates
(`code/capacity_readback/`), and the physical identities I1/I2 are outside
the formalised set.
-/

namespace OPH.CapacityFixedPoint

open Function

open scoped NNReal

/-- The averaging map `C(x) = (1 − λ)·x + λ·K` toward the constant
    target `K`. -/
def averagingMap (lam K x : ℝ) : ℝ := (1 - lam) * x + lam * K

/-- For `λ ≠ 0`, `x` is a fixed point of the averaging map iff `x = K`. -/
theorem averagingMap_isFixedPt_iff {lam K : ℝ} (hlam : lam ≠ 0) (x : ℝ) :
    IsFixedPt (averagingMap lam K) x ↔ x = K := by
  simp only [IsFixedPt, averagingMap]
  constructor
  · intro h
    have h0 : lam * (K - x) = 0 := by linear_combination h
    rcases mul_eq_zero.mp h0 with h' | h'
    · exact absurd h' hlam
    · linarith
  · rintro rfl
    ring

/-- **Averaging-map uniqueness schema.** For `0 < λ ≤ 1` the averaging map
    `C(x) = (1 − λ)·x + λ·K` has exactly one fixed point, namely `K`.
    (Only `λ ≠ 0` is used; the upper bound is carried because the
    certificates instantiate the schema with `λ ≤ 1`.) -/
theorem averagingMap_unique_fixedPt {lam K : ℝ} (h0 : 0 < lam) (_h1 : lam ≤ 1) :
    ∃! x : ℝ, IsFixedPt (averagingMap lam K) x :=
  ⟨K, (averagingMap_isFixedPt_iff h0.ne' K).mpr rfl,
    fun y hy => (averagingMap_isFixedPt_iff h0.ne' y).mp hy⟩

/-- **Banach uniqueness schema (cited from Mathlib).** A map that is
    `ContractingWith K` (Lipschitz with constant `K < 1`) on a nonempty
    complete metric space has exactly one fixed point. Wraps
    `ContractingWith.fixedPoint` / `ContractingWith.fixedPoint_unique`. -/
theorem banach_unique_fixedPt {α : Type*} [MetricSpace α]
    [Nonempty α] [CompleteSpace α] {K : ℝ≥0} {f : α → α}
    (hf : ContractingWith K f) :
    ∃! x : α, IsFixedPt f x :=
  ⟨ContractingWith.fixedPoint f hf, hf.fixedPoint_isFixedPt,
    fun _ hy => hf.fixedPoint_unique hy⟩

/-- **Capacity interval schema.** A self-map of the interval `Icc a b`
    (`a ≤ b`) that is Lipschitz on `Icc a b` with constant `K < 1` has
    exactly one fixed point in `Icc a b`. This is the form in which the
    Banach schema applies to the capacity map's certified interval domain;
    the numeric `a`, `b`, `K` witnessing the hypotheses live in the Python
    certificates, not in Lean. -/
theorem capacity_Icc_unique_fixedPt {a b : ℝ} (hab : a ≤ b)
    {K : ℝ≥0} (hK : K < 1) {f : ℝ → ℝ}
    (hmaps : Set.MapsTo f (Set.Icc a b) (Set.Icc a b))
    (hlip : LipschitzOnWith K f (Set.Icc a b)) :
    ∃! x : ℝ, x ∈ Set.Icc a b ∧ f x = x := by
  haveI : Nonempty (Set.Icc a b) := (Set.nonempty_Icc.mpr hab).to_subtype
  haveI : CompleteSpace (Set.Icc a b) := isClosed_Icc.isComplete.completeSpace_coe
  have hg : ContractingWith K (hmaps.restrict f _ _) :=
    ⟨hK, hlip.mapsToRestrict hmaps⟩
  obtain ⟨x, hx, huniq⟩ := banach_unique_fixedPt hg
  refine ⟨(x : ℝ), ⟨x.2, ?_⟩, ?_⟩
  · have hval := congrArg Subtype.val hx
    simpa using hval
  · rintro y ⟨hy, hfy⟩
    have hyfix : IsFixedPt (hmaps.restrict f _ _) ⟨y, hy⟩ := by
      apply Subtype.ext
      simpa using hfy
    exact congrArg Subtype.val (huniq ⟨y, hy⟩ hyfix)

-- Axiom audit: these must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms averagingMap_isFixedPt_iff
#print axioms averagingMap_unique_fixedPt
#print axioms banach_unique_fixedPt
#print axioms capacity_Icc_unique_fixedPt

end OPH.CapacityFixedPoint
