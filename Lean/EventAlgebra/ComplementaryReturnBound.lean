import EventAlgebra.Basic

/-! # A lower bound from two compatible return tests

For two commuting return events, the weight of their intersection is
at least the sum of the two individual return weights minus one.
For the two complementary qubit bases the events act on a reference
and the returned system. Their intersection is the ideal Bell ray.
The identification of those projectors and the statistical application
are given in the executable pair-return experiment.
-/

namespace EventAlgebra

open Matrix

variable {n : ℕ}

/-- A normalized finite observer table gives a lower bound on the joint
return weight from the two marginal weights. Leakage is allowed. -/
theorem joint_return_weight_lower {J Z X : Matrix (Fin n) (Fin n) ℂ}
    (hJ : IsState J) (hZ : IsEvent Z) (hX : IsEvent X)
    (hZX : Z * X = X * Z) :
    (bornWeight J Z).re + (bornWeight J X).re - 1 ≤
      (bornWeight J (Z * X)).re := by
  have hc : (1 - Z) * (1 - X) = (1 - X) * (1 - Z) := by
    simp only [sub_mul, mul_sub, one_mul, mul_one, hZX]
    abel
  have h := bornWeight_re_nonneg hJ.1 (hZ.compl.mul_of_commute hX.compl hc)
  simp only [sub_mul, mul_sub, one_mul, mul_one, bornWeight_sub,
    bornWeight_one hJ, Complex.sub_re, Complex.one_re] at h
  linarith

#print axioms joint_return_weight_lower

end EventAlgebra
