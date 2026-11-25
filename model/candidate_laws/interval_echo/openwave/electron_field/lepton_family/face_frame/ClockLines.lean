import Mathlib

namespace ObserverRecords.FamilyClock

/-- The evolution-versus-idle amplitude for three stationary components. -/
def sample (a b c x y z : ℂ) (n : ℕ) : ℂ :=
  a * x ^ n + b * y ^ n + c * z ^ n

/-- Determinant of the three-by-three Hankel matrix of five consecutive readings. -/
def witness (s : ℕ → ℂ) : ℂ :=
  s 0 * (s 2 * s 4 - s 3 ^ 2) -
  s 1 * (s 1 * s 4 - s 2 * s 3) +
  s 2 * (s 1 * s 3 - s 2 ^ 2)

theorem factorizes (a b c x y z : ℂ) :
    witness (sample a b c x y z) = a * b * c * ((x-y)*(x-z)*(y-z)) ^ 2 := by
  simp only [witness, sample]
  ring

theorem three_visible_lines (a b c x y z : ℂ)
    (ha : a ≠ 0) (hb : b ≠ 0) (hc : c ≠ 0)
    (hxy : x ≠ y) (hxz : x ≠ z) (hyz : y ≠ z) :
    witness (sample a b c x y z) ≠ 0 := by
  rw [factorizes]
  exact mul_ne_zero (mul_ne_zero (mul_ne_zero ha hb) hc)
    (pow_ne_zero _ (mul_ne_zero (mul_ne_zero (sub_ne_zero.mpr hxy)
      (sub_ne_zero.mpr hxz)) (sub_ne_zero.mpr hyz)))

/-- Every two-tone clock has zero witness at the very same five delays. -/
theorem two_lines (a b x y : ℂ) : witness (fun n => a*x^n + b*y^n) = 0 := by
  unfold witness
  ring

end ObserverRecords.FamilyClock
