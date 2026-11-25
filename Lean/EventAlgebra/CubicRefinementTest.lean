import EventAlgebra.OperationalAdditivityBoundary

set_option autoImplicit false

namespace EventAlgebra

/-!
# Resolving a quantitative ambiguity in three public record menus

A three-outcome instrument assigns probabilities a and b to two distinct
public records representing the same half effect. Their union has probability
a + b by ordinary counting. Separate readers of the half effect and the union
may be implemented differently. The committed transverse cubic has an exact
minimum disagreement across these three comparisons. No Born rule or matrix
state occurs in the lower bound.
-/

/-- For any calibrated binary readings of a half effect and its double, three
cross-context errors bounded by ε imply a discrepancy at most 3ε. -/
theorem record_refinement_discrepancy
    (a b half whole ε : ℝ)
    (ha : |a - half| ≤ ε) (hb : |b - half| ≤ ε)
    (hu : |a + b - whole| ≤ ε) : |2 * half - whole| ≤ 3 * ε := by
  rcases abs_le.mp ha with ⟨ha₁, ha₂⟩
  rcases abs_le.mp hb with ⟨hb₁, hb₂⟩
  rcases abs_le.mp hu with ⟨hu₁, hu₂⟩
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- The eight-context fixture admits the committed cubic valuation, but a
common refinement must disagree with one separate binary reading by at least
one part in 256. -/
theorem cubic_refinement_minimum
    (a b ε : ℝ)
    (ha : |a - producedCubicValuation halfWitness| ≤ ε)
    (hb : |b - producedCubicValuation halfWitness| ≤ ε)
    (hu : |a + b - producedCubicValuation witnessSum| ≤ ε) : 1/256 ≤ ε := by
  have h := record_refinement_discrepancy a b _ _ ε ha hb hu
  rw [producedCubicValuation_halfWitness, producedCubicValuation_witnessSum] at h
  norm_num at h
  linarith

/-- A normalized ternary law attains the lower bound simultaneously on all
three comparisons. -/
theorem cubic_refinement_minimum_attained :
    ∃ a b c : ℝ, 0 ≤ a ∧ 0 ≤ b ∧ 0 ≤ c ∧ a + b + c = 1 ∧
      |a - producedCubicValuation halfWitness| = 1/256 ∧
      |b - producedCubicValuation halfWitness| = 1/256 ∧
      |a + b - producedCubicValuation witnessSum| = 1/256 := by
  refine ⟨141/512, 141/512, 115/256, ?_⟩
  norm_num [producedCubicValuation_halfWitness, producedCubicValuation_witnessSum]

#print axioms record_refinement_discrepancy
#print axioms cubic_refinement_minimum
#print axioms cubic_refinement_minimum_attained

end EventAlgebra
