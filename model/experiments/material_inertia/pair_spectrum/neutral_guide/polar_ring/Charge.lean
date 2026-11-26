import Mathlib.Data.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

namespace IntervalEcho.PolarRing

/-- Each of the real and imaginary parts of a compressed charge effect is a
scalar matrix if it commutes with every local phase and orientation exchange.
A finite set names the dimers in the upper orientation. -/
theorem scalar_effect {ι : Type*} [DecidableEq ι]
    (K : Finset ι → Finset ι → ℝ)
    (phase : ∀ (x y : Finset ι) (i : ι),
      (if i ∈ x then -(K x y) else K x y) =
      (if i ∈ y then -(K x y) else K x y))
    (toggle : ∀ (i : ι) (x : Finset ι), i ∉ x →
      K (insert i x) (insert i x) = K x x)
    (x y : Finset ι) : K x y = if x = y then K ∅ ∅ else 0 := by
  classical
  have diagonal (s : Finset ι) : K s s = K ∅ ∅ := by
    induction s using Finset.induction_on with
    | empty => rfl
    | @insert i t hi ih => exact (toggle i t hi).trans ih
  by_cases same : x = y
  · subst y
    simp [diagonal]
  · rw [if_neg same]
    by_contra nonzero
    apply same
    ext i
    have h := phase x y i
    by_cases hx : i ∈ x <;> by_cases hy : i ∈ y
    · simp [hx, hy]
    · simp [hx, hy] at h
      exact False.elim (nonzero (by linarith))
    · simp [hx, hy] at h
      exact False.elim (nonzero (by linarith))
    · simp [hx, hy]

end IntervalEcho.PolarRing
