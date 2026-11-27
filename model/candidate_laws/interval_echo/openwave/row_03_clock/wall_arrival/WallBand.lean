import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace OPH.WallArrival

/-- A signed pair of constituent energies is within two of a multiple
of 32 once each signed energy is within one of its signed centre. -/
theorem pair_radius (a b s t : ℝ) (m : ℤ)
    (ha : |a-s| ≤ 1) (hb : |b-t| ≤ 1) (hc : s+t = 32*(m:ℝ)) :
    |(a+b)-32*(m:ℝ)| ≤ 2 := by
  rw [← hc, show (a+b)-(s+t) = (a-s)+(b-t) by ring]
  exact (abs_add_le _ _).trans (by linarith)

/-- Adding r such pairs gives distance at most 2r from the same lattice. -/
theorem paired_radius (xs : List ℝ)
    (h : ∀ x ∈ xs, ∃ m : ℤ, |x-32*(m:ℝ)| ≤ 2) :
    ∃ k : ℤ, |xs.sum-32*(k:ℝ)| ≤ 2*(xs.length:ℝ) := by
  induction xs with
  | nil => exact ⟨0, by simp⟩
  | cons a rest ih =>
    obtain ⟨m, hm⟩ := h a (by simp)
    obtain ⟨k, hk⟩ := ih (by
      intro x hx
      exact h x (List.mem_cons_of_mem _ hx))
    refine ⟨m+k, ?_⟩
    simp only [List.sum_cons, List.length_cons, Nat.cast_add_one, Int.cast_add]
    rw [show a + rest.sum - 32*((m:ℝ)+(k:ℝ)) =
      (a-32*(m:ℝ))+(rest.sum-32*(k:ℝ)) by ring]
    exact (abs_add_le _ _).trans (by linarith)

/-- The radius-fourteen neighbourhood of 32Z misses (14,18). -/
theorem silent_interval (w : ℝ) (k : ℤ) (h : |w-32*(k:ℝ)| ≤ 14)
    (lo : 14 < w) (hi : w < 18) : False := by
  obtain ⟨hl,hu⟩ := abs_le.mp h
  by_cases hk : k ≤ 0
  · have hkr : (k:ℝ) ≤ 0 := by exact_mod_cast hk
    linarith
  · have hk1 : (1:ℤ) ≤ k := by omega
    have hkr : (1:ℝ) ≤ k := by exact_mod_cast hk1
    linarith

/-- Arbitrary fixed effects on the exact one-constituent subspace. -/
theorem one_constituent_width (a b : ℝ) (ha : 15 ≤ a ∧ a ≤ 17)
    (hb : 15 ≤ b ∧ b ≤ 17) : |a-b| ≤ 2 := by
  obtain ⟨al,au⟩ := ha
  obtain ⟨bl,bu⟩ := hb
  exact abs_le.mpr ⟨by linarith, by linarith⟩

#print axioms paired_radius
#print axioms silent_interval

end OPH.WallArrival
