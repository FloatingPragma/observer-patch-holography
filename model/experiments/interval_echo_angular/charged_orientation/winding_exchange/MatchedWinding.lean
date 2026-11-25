import Mathlib.Tactic

namespace OPH.MatchedWinding

/-- Two flow records on an oriented path with the same charge divergence
differ by a constant circulation. A cycle restricts this identity to one lap. -/
theorem equal_divergence_constant (a b : ℕ → ℤ)
    (h : ∀ j, a (j+1) - a j = b (j+1) - b j) (k : ℕ) :
    a k - b k = a 0 - b 0 := by
  induction k with
  | zero => rfl
  | succ j ih =>
    change a (j+1) - b (j+1) = a 0 - b 0
    linarith [h j]

/-- A diagonal head rephasing that squares to the identity transfers every
power of a finite generator, including all terms of the free propagator. -/
theorem head_gauge_powers {R : Type*} [Monoid R] (D H : R)
    (hd : D * D = 1) (n : ℕ) : (D * H * D)^n = D * H^n * D := by
  induction n with
  | zero => simp [hd]
  | succ k ih =>
    rw [pow_succ, ih, pow_succ]
    calc
      D * H^k * D * (D * H * D) = D * H^k * (D * D) * H * D := by
        simp only [mul_assoc]
      _ = D * (H^k * H) * D := by simp [hd, mul_assoc]

end OPH.MatchedWinding
