import GibbsEntropyLimit

set_option autoImplicit false

/-!
# Gap-free finite Gibbs energy control

The entropy identity and log-cardinality bound refer to the same finite
Gibbs family. A single attained energy level lower-bounds its partition
function. Taking this level to be the ground energy gives an explicit
mean-excess bound without a nonzero spectral-gap premise.

For a finite Hermitian source Hamiltonian, its complete spectral list
instantiates the finite carrier. A physical interpretation of its energy,
and approximation or refinement of that state, require their own bounds.
-/

namespace OPH.Thermodynamics

open scoped BigOperators

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω] [Nonempty Ω]

theorem gibbsMean_sub_attained_le_log_card_div
    (E : Ω → ℝ) (β E₀ : ℝ) (hβ : 0 < β)
    (hattained : ∃ x, E x = E₀) :
    (∑ x, gibbs E β x * E x) - E₀ ≤ Real.log (Fintype.card Ω) / β := by
  obtain ⟨x₀, hx₀⟩ := hattained
  have hz : Real.exp (-β * E₀) ≤ partitionZ E β := by
    rw [← hx₀]
    exact Finset.single_le_sum
      (fun x _ => (Real.exp_pos (-β * E x)).le) (Finset.mem_univ x₀)
  have hl := Real.log_le_log (Real.exp_pos (-β * E₀)) hz
  rw [Real.log_exp] at hl
  have hs := gibbsEntropy_le_log_card E β
  rw [gibbsEntropy_normalization] at hs
  apply (le_div_iff₀ hβ).mpr
  nlinarith

/-- The energy budget uses normalized original weights, without a ground gap. -/
theorem gibbs_excess_sum_le_log_card_div
    (E : Ω → ℝ) (β E₀ : ℝ) (hβ : 0 < β)
    (hattained : ∃ x, E x = E₀) :
    (∑ x, gibbs E β x * (E x - E₀)) ≤ Real.log (Fintype.card Ω) / β := by
  have hsum : ∑ x, gibbs E β x = 1 := by
    simp only [gibbs, ← Finset.sum_div]
    exact div_self (partitionZ_pos E β).ne'
  simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul, hsum, one_mul]
  exact gibbsMean_sub_attained_le_log_card_div E β E₀ hβ hattained

end OPH.Thermodynamics
