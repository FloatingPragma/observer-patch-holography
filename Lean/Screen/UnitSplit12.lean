import Mathlib

open scoped BigOperators

namespace OPH.UnitSplit12

/-! This file proves only the arithmetic consequence of a twelve-slot
positive-integer split with total weight twelve.  It does not derive the
number of slots, an Euler/cost functional, a source selector, an `A5` action,
or a physical current algebra. -/

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

#print axioms unit_split_of_positive_sum

end OPH.UnitSplit12
