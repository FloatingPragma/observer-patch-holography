import Mathlib

/-!
# Typed capacity-closure principle and reserve-branch boundary

This module separates three formal steps that are easy to conflate in prose.

1. An equivalence between the simulator-side and simulated-side state types,
   together with a commuting readback square, makes their readings equal.
2. A source return law plus one such closure identification makes the
   identified value a fixed point.
3. Uniqueness requires a separate hypothesis about the source return map.

The last section records the algebraic distinction between the finite
one-step presence factor `1 - P/24` and the Poisson/projective-limit factor
`exp (-P/24)`.  For nonzero `P` they are unequal.  A translation-invariant
normalized weight on the six `ZMod 6` classes is uniform.  Neither theorem
selects a physical reserve reading, identifies a global capacity, or supplies
a horizon attachment.
-/

namespace OPH.CapacityClosurePrinciple

open Function

/-! ## Typed commuting readback -/

/-- Two state descriptions, an explicit equivalence between their types, and
    readout maps into one shared quantity type.  The `commutes` field is the
    substantive same-quantity identification. -/
structure EquivalentReadback (Outer Inner Quantity : Type*) where
  identify : Outer ≃ Inner
  outerRead : Outer → Quantity
  innerRead : Inner → Quantity
  commutes : ∀ outer, innerRead (identify outer) = outerRead outer

/-- Once the typed readback square commutes, the two readings are equal. -/
theorem readings_eq_of_commuting_square
    {Outer Inner Quantity : Type*}
    (readback : EquivalentReadback Outer Inner Quantity)
    (outer : Outer) :
    readback.innerRead (readback.identify outer) =
      readback.outerRead outer :=
  readback.commutes outer

/-! ## Source return and fixed points -/

/-- A declared source return map.  Constructing this map from the source
    architecture is separate from the closure principle. -/
structure SourceReturnLaw (Quantity : Type*) where
  returnMap : Quantity → Quantity

/-- One same-quantity closure identification for a declared source return
    law.  The equality is the commuting readback obligation. -/
structure ClosureIdentification {Quantity : Type*}
    (law : SourceReturnLaw Quantity) where
  value : Quantity
  commutes : law.returnMap value = value

/-- The closure equality makes the identified value a fixed point of the
    declared source return map. -/
theorem source_return_fixed
    {Quantity : Type*}
    {law : SourceReturnLaw Quantity}
    (closure : ClosureIdentification law) :
    IsFixedPt law.returnMap closure.value :=
  closure.commutes

/-- A uniqueness hypothesis selects the closure value among the fixed points.
    Uniqueness is an input here, rather than a consequence of self-reference
    or of the commuting readback equality. -/
theorem closure_value_eq_of_unique
    {Quantity : Type*}
    {law : SourceReturnLaw Quantity}
    (closure : ClosureIdentification law)
    (unique : ∀ x, IsFixedPt law.returnMap x → x = closure.value)
    {candidate : Quantity}
    (fixed : IsFixedPt law.returnMap candidate) :
    candidate = closure.value :=
  unique candidate fixed

/-! ## Finite-presence and Poisson branches -/

/-- The finite one-step presence factor attached to reserve coordinate `P/24`. -/
noncomputable def finitePresenceFactor (P : ℝ) : ℝ := 1 - P / 24

/-- The Poisson or projective-limit factor attached to reserve coordinate
    `P/24`. -/
noncomputable def poissonFactor (P : ℝ) : ℝ := Real.exp (-P / 24)

/-- For nonzero `P`, the finite one-step presence factor is strictly below
    the Poisson/projective-limit factor.  They therefore cannot be exchanged
    without an additional semantic premise. -/
theorem finitePresenceFactor_lt_poissonFactor
    {P : ℝ} (hP : P ≠ 0) :
    finitePresenceFactor P < poissonFactor P := by
  have hne : -P / 24 ≠ 0 := div_ne_zero (neg_ne_zero.mpr hP) (by norm_num)
  change 1 - P / 24 < Real.exp (-P / 24)
  calc
    1 - P / 24 = -P / 24 + 1 := by ring
    _ < Real.exp (-P / 24) := Real.add_one_lt_exp hne

/-- Multiplication by a positive conditional baseline preserves the strict
    distinction between the two capacity rows. -/
theorem finitePresenceCapacity_lt_poissonCapacity
    {N0 P : ℝ} (hN0 : 0 < N0) (hP : P ≠ 0) :
    N0 * finitePresenceFactor P < N0 * poissonFactor P :=
  mul_lt_mul_of_pos_left (finitePresenceFactor_lt_poissonFactor hP) hN0

/-! ## Six-class uniformity -/

/-- Translation invariance on all six residue classes makes the weight
    constant.  No objective, grammar-completeness, or physical attachment
    premise is needed once this invariance is supplied. -/
theorem z6_translation_invariant_constant
    (weight : ZMod 6 → ℝ)
    (invariant : ∀ shift index, weight (index + shift) = weight index) :
    ∀ i j, weight i = weight j := by
  intro i j
  have h := invariant (j - i) i
  simpa [sub_eq_add_neg, add_assoc, add_comm, add_left_comm] using h.symm

/-- A normalized translation-invariant weight on the six `ZMod 6` classes is
    exactly uniform. -/
theorem z6_translation_invariant_uniform
    (weight : ZMod 6 → ℝ)
    (invariant : ∀ shift index, weight (index + shift) = weight index)
    (normalized : ∑ index, weight index = 1) :
    ∀ index, weight index = 1 / 6 := by
  intro index
  have constant := z6_translation_invariant_constant weight invariant
  have hsum : ∑ j : ZMod 6, weight j = 6 * weight index := by
    calc
      ∑ j : ZMod 6, weight j = ∑ _j : ZMod 6, weight index := by
        apply Finset.sum_congr rfl
        intro j _hj
        exact constant j index
      _ = 6 * weight index := by norm_num
  rw [normalized] at hsum
  linarith

-- Axiom audit: these must report only the standard Mathlib classical axioms.
#print axioms readings_eq_of_commuting_square
#print axioms source_return_fixed
#print axioms closure_value_eq_of_unique
#print axioms finitePresenceFactor_lt_poissonFactor
#print axioms finitePresenceCapacity_lt_poissonCapacity
#print axioms z6_translation_invariant_constant
#print axioms z6_translation_invariant_uniform

end OPH.CapacityClosurePrinciple
