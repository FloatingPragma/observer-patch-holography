import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Data.Complex.Basic

open scoped BigOperators

namespace OPH.CompactEnergySum

/-- Once the material gateway terms annihilate the product, each surviving
local matrix acts on its own coordinate. Summing the local eigenvalue
equations yields the exact many-body eigenvalue, even at nodal coordinates. -/
theorem product_spectrum {I X : Type*} [Fintype I] [Fintype X] [DecidableEq I]
    (h : I → X → X → ℂ) (v : I → X → ℂ) (e : I → ℂ)
    (localEigen : ∀ i b, ∑ a, h i b a * v i a = e i * v i b) (x : I → X) :
    (∑ i, (∑ a, h i (x i) a * v i a) *
      ∏ j ∈ Finset.univ.erase i, v j (x j)) =
    (∑ i, e i) * ∏ j, v j (x j) := by
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro i hi
  rw [localEigen, mul_assoc, Finset.mul_prod_erase Finset.univ (fun j => v j (x j)) hi]

#print axioms product_spectrum
end OPH.CompactEnergySum
