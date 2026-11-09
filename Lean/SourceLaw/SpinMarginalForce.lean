import Mathlib.Data.Real.Basic
import Mathlib.Tactic

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A finite two-site bound on the mean-spin force replacement

For a normalized two-spin marginal with diagonal `(a,b,c,a)` and real
middle off-diagonal entry `z`, positivity bounds the squared ordered moment
by `1/4-a+a^2-z^2`. An antiferromagnetic bond-energy gain `δ` gives
`a+z ≤ -δ`. These explicit scalar premises imply a positive difference
between the exact twist coefficient `1/4-a-z/2` and the squared moment.
The source marginal, energy estimate, joint limit and force attachment
are separate executable/prose premises; this module proves the finite
inequality and its rational specialization.
-/

namespace OPH.SourceLaw.SpinMarginalForce

/-- The exact twist coefficient cannot equal the squared mean spin under these premises. -/
theorem excess_lower {a z m δ : ℝ}
    (ha : 0 ≤ a) (hδ : 0 ≤ δ) (hz : a + z ≤ -δ)
    (hm : m ^ 2 ≤ 1 / 4 - a + a ^ 2 - z ^ 2) :
    δ ^ 2 + δ / 2 ≤ (1 / 4 - a - z / 2) - m ^ 2 := by
  have hfirst : z + a + δ ≤ 0 := by linarith
  have hsecond : z - a - δ ≤ 0 := by linarith
  have hprod := mul_nonneg_of_nonpos_of_nonpos hfirst hsecond
  have had := mul_nonneg ha hδ
  nlinarith

/-- The zero-field source bond gain `1/120` forces the stated rational margin. -/
theorem source_margin {a z m : ℝ}
    (ha : 0 ≤ a) (hz : a + z ≤ -(1 / 120 : ℝ))
    (hm : m ^ 2 ≤ 1 / 4 - a + a ^ 2 - z ^ 2) :
    (61 / 14400 : ℝ) ≤ (1 / 4 - a - z / 2) - m ^ 2 := by
  have h := excess_lower ha (by norm_num : (0 : ℝ) ≤ 1 / 120) hz hm
  norm_num at h ⊢
  exact h

#print axioms excess_lower
#print axioms source_margin

end OPH.SourceLaw.SpinMarginalForce
