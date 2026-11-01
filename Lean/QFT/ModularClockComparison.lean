import Mathlib.Data.Complex.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.Group
import Mathlib.Tactic.Ring

/-!
# Relative modular clock comparisons and finite resonance

The group identities apply to the global and embedded regional unitary flows
of one faithful finite field reference. They require no regional covariance.
The coefficient theorem is the finite spectral intertwiner test: a linear
map cannot transmit a target clock frequency absent from the source.

The actual dimer spectral attachment and rational exponential/logarithm
enclosures are checked separately in `model/modular_clock_comparison`.
They are not encoded as a Lean value here. No physical time, reference
selection, Einstein premise or continuum theorem is asserted.
-/

namespace OPH.QFT.ModularClockComparison

section Comparison

variable {G : Type*} [Group G]

def relative (target origin : G) : G := target * origin⁻¹

theorem relative_self (a : G) : relative a a = 1 := by
  simp [relative]

theorem relative_transitive (a b c : G) :
    relative a b * relative b c = relative a c := by
  simp [relative, mul_assoc]

theorem relative_corrects (a b : G) : relative a b * b = a := by
  simp [relative, mul_assoc]

/-- Comparison is a time cocycle for the origin action. It need not be a
one-parameter group on its own. -/
theorem relative_time_cocycle (target origin : ℝ → G)
    (ht : ∀ s t, target (s + t) = target s * target t)
    (ho : ∀ s t, origin (s + t) = origin s * origin t)
    (s t : ℝ) :
    relative (target (s+t)) (origin (s+t)) =
      relative (target s) (origin s) *
        (origin s * relative (target t) (origin t) * (origin s)⁻¹) := by
  rw [ht, ho]
  simp only [relative, mul_inv_rev]
  group

end Comparison

section Resonance

variable {I J : Type*}

/-- Matrix coefficients of a linear map intertwine diagonal generators.
For modular generators these are the differences of logarithmic weights. -/
def FrequencyIntertwines (source : I → ℂ) (target : J → ℂ)
    (coefficients : J → I → ℂ) : Prop :=
  ∀ j i, target j * coefficients j i = coefficients j i * source i

theorem nonresonant_coefficient_zero
    (source : I → ℂ) (target : J → ℂ) (coefficients : J → I → ℂ)
    (h : FrequencyIntertwines source target coefficients) (j : J) (i : I)
    (hne : target j ≠ source i) : coefficients j i = 0 := by
  have hz : (target j - source i) * coefficients j i = 0 := by
    calc
      (target j - source i) * coefficients j i =
          target j * coefficients j i - coefficients j i * source i := by ring
      _ = 0 := sub_eq_zero.mpr (h j i)
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr hne)

/-- If all nonzero target frequencies miss the source spectrum, every
intertwiner has stationary range. This is a statement about every linear
map; CP/unitality assumptions cannot restore the missing frequencies. -/
theorem nonresonant_output_stationary [Fintype I]
    (source : I → ℂ) (target : J → ℂ) (coefficients : J → I → ℂ)
    (h : FrequencyIntertwines source target coefficients)
    (hn : ∀ j, target j ≠ 0 → ∀ i, target j ≠ source i)
    (input : I → ℂ) (j : J) :
    target j * (∑ i, coefficients j i * input i) = 0 := by
  by_cases hj : target j = 0
  · simp [hj]
  · have hc : ∀ i, coefficients j i = 0 :=
      fun i => nonresonant_coefficient_zero source target coefficients h j i (hn j hj i)
    simp [hc]

end Resonance

#print axioms relative_transitive
#print axioms relative_time_cocycle
#print axioms nonresonant_coefficient_zero
#print axioms nonresonant_output_stationary

end OPH.QFT.ModularClockComparison
