import Mathlib

open scoped BigOperators

namespace OPH.UnitSplit12

/-! This file proves only the arithmetic consequence of a twelve-slot
positive-integer split with total weight twelve.  It does not derive the
number of slots, an Euler/cost functional, a source selector, an `A5` action,
or a physical current algebra.

Removability controls (issue #568): each premise of the unit-splitting
theorem is load-bearing.  Dropping integrality (rational weights), dropping
strict positivity, or changing the twelve-slot domain each admits an
explicit non-unit witness, so the theorem fails without them. -/

/-- Twelve positive natural-number weights whose sum is twelve are all one. -/
theorem unit_split_of_positive_sum
    (q : Fin 12 → ℕ)
    (hpos : ∀ i, 1 ≤ q i)
    (hsum : ∑ i, q i = 12) :
    ∀ i, q i = 1 := by
  intro i
  have hrest :
      ∑ j ∈ (Finset.univ.erase i), 1 ≤
        ∑ j ∈ (Finset.univ.erase i), q j := by
    exact Finset.sum_le_sum fun j _ ↦ hpos j
  have hrest_count : ∑ _j ∈ (Finset.univ.erase i), 1 = 11 := by
    simp
  have hdecomp :
      (∑ j ∈ (Finset.univ.erase i), q j) + q i = 12 := by
    rw [Finset.sum_erase_add Finset.univ q (Finset.mem_univ i)]
    exact hsum
  rw [hrest_count] at hrest
  have hqi := hpos i
  omega

/-- Removing integrality prevents the theorem: positive rational weights
summing to twelve need not be units. -/
theorem no_unit_split_without_integrality :
    ∃ q : Fin 12 → ℚ, (∀ i, 0 < q i) ∧ (∑ i, q i = 12) ∧ ¬(∀ i, q i = 1) := by
  refine ⟨fun i => if i = 0 then 1/2 else if i = 1 then 3/2 else 1, ?_, ?_, ?_⟩
  · intro i
    dsimp only
    split_ifs <;> norm_num
  · simp [Fin.sum_univ_succ]
    norm_num
  · intro h
    have h0 := h 0
    norm_num at h0

/-- Removing strict positivity prevents the theorem: nonnegative integer
weights summing to twelve need not be units. -/
theorem no_unit_split_without_positivity :
    ∃ q : Fin 12 → ℕ, (∑ i, q i = 12) ∧ ¬(∀ i, q i = 1) := by
  refine ⟨fun i => if i = 0 then 0 else if i = 1 then 2 else 1, ?_, ?_⟩
  · decide
  · intro h
    have h0 := h 0
    simp at h0

/-- Changing the twelve-slot domain prevents the theorem: eleven positive
integer weights summing to twelve need not be units. -/
theorem no_unit_split_without_twelve_slots :
    ∃ q : Fin 11 → ℕ, (∀ i, 1 ≤ q i) ∧ (∑ i, q i = 12) ∧ ¬(∀ i, q i = 1) := by
  refine ⟨fun i => if i = 0 then 2 else 1, ?_, ?_, ?_⟩
  · intro i
    dsimp only
    split_ifs <;> omega
  · decide
  · intro h
    have h0 := h 0
    simp at h0

#print axioms unit_split_of_positive_sum
#print axioms no_unit_split_without_integrality
#print axioms no_unit_split_without_positivity
#print axioms no_unit_split_without_twelve_slots

end OPH.UnitSplit12
