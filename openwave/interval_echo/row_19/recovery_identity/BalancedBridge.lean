import Mathlib.Tactic

namespace BalancedBridge

/-- For three physical endpoint differences the cofactors balance the two
reflection-symmetric bridge outcomes. This identity is exact at population
probabilities even when the experiment only bounds those probabilities. -/
theorem cross_balances (s a : Fin 3 → ℝ) :
    let c : Fin 3 → ℝ := ![s 1 * a 2 - s 2 * a 1, s 2 * a 0 - s 0 * a 2,
                            s 0 * a 1 - s 1 * a 0]
    (∑ i, c i * s i) = 0 ∧ (∑ i, c i * a i) = 0 := by
  dsimp
  constructor <;> simp [Fin.sum_univ_three] <;> ring

/-- The signs of the cofactors exchange the two endpoints of each physical
preparation. The absolute coefficients are ordinary mixture weights. -/
theorem endpoint_mix {n : ℕ} (c plus minus : Fin n → ℝ)
    (hn : 0 < ∑ i, |c i|) :
    let z := ∑ i, |c i|
    (∀ i, 0 ≤ |c i|/z) ∧ (∑ i, |c i|/z) = 1 ∧
    (∑ i, |c i|/z * (if 0 ≤ c i then plus i else minus i)) -
    (∑ i, |c i|/z * (if 0 ≤ c i then minus i else plus i)) =
    (∑ i, c i * (plus i - minus i))/z := by
  dsimp
  refine ⟨fun i => div_nonneg (abs_nonneg _) hn.le, ?_, ?_⟩
  · rw [← Finset.sum_div, div_self hn.ne']
  · rw [← Finset.sum_sub_distrib, Finset.sum_div]
    apply Finset.sum_congr rfl
    intro i _
    by_cases h : 0 ≤ c i
    · simp [h, abs_of_nonneg h]
      ring
    · simp [h, abs_of_neg (lt_of_not_ge h)]
      ring

theorem no_range_function (u v : ℝ × ℝ) (x y : ℝ)
    (same : u = v) (different : x ≠ y) :
    ¬ ∃ f : (ℝ × ℝ) → ℝ, f u = x ∧ f v = y := by
  rintro ⟨f, hx, hy⟩
  exact different (hx.symm.trans ((congrArg f same).trans hy))

end BalancedBridge
