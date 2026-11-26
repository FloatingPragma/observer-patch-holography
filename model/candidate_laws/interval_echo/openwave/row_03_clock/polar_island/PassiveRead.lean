import Mathlib.Analysis.Complex.Exponential
import Mathlib.Tactic.Linarith

namespace IntervalEcho.PolarIslandClock

/-- On a band carrying one copy of each translation character, a translated
average of any terminal effect has only diagonal matrix elements. -/
theorem read_diagonal {ι : Type*} (phase : ι → ℂ) (hphase : Function.Injective phase)
    (read : ι → ι → ℂ) (commute : ∀ a b, phase a * read a b = phase b * read a b)
    {a b : ι} (hne : a ≠ b) : read a b = 0 := by
  have h : (phase a - phase b) * read a b = 0 := by
    rw [sub_mul, commute, sub_self]
  exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr (fun he => hne (hphase he)))

/-- Every translated average has constant probability for any incoming
density matrix and arbitrary free energies on this multiplicity-one band. -/
theorem read_stationary {ι : Type*} [Fintype ι] [DecidableEq ι]
    (phase : ι → ℂ) (hphase : Function.Injective phase)
    (read : ι → ι → ℂ) (commute : ∀ a b, phase a * read a b = phase b * read a b)
    (rho : ι → ι → ℂ) (energy : ι → ℝ) (t : ℝ) :
    (∑ a, ∑ b, rho b a * read a b *
      Complex.exp (Complex.I * ((energy a - energy b : ℝ) : ℂ) * (t : ℂ))) =
      ∑ a, rho a a * read a a := by
  apply Finset.sum_congr rfl
  intro a _
  rw [Finset.sum_eq_single a]
  · simp
  · intro b _ hb
    simp [read_diagonal phase hphase read commute (Ne.symm hb)]
  · simp

/-- With the min-max interval for the N-dimer band and a sharp all-lower
reference of zero mean energy, all band beats lie below the rest gap. -/
theorem line_below_rest {e₁ e₂ er g n κ : ℝ}
    (lo₁ : 5 * κ - n / 2 ≤ e₁) (hi₁ : e₁ ≤ 5 * κ + n / 2)
    (lo₂ : 5 * κ - n / 2 ≤ e₂) (hi₂ : e₂ ≤ 5 * κ + n / 2)
    (lor : 5 * κ - n / 2 ≤ er) (hg : g ≤ 0) (separate : 3 * n < 10 * κ) :
    |e₁ - e₂| < er - g := by
  apply abs_lt.mpr
  constructor <;> linarith

end IntervalEcho.PolarIslandClock
