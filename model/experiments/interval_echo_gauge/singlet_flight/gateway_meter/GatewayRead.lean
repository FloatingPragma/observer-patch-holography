import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith

namespace IntervalEcho.GatewayMeter

/-- A coefficient of (sqrt(5) B + C)^2 equals the corresponding
coefficient of 4 sqrt(5) (sqrt(5) B + C), using the integer incidence data. -/
theorem projector_product (b c bb bc cb cc s : ℝ) (root : s^2 = 5)
    (hbb : bb = 2*b) (hbc : bc = 2*c) (hcb : cb = 2*c) (hcc : cc = 10*b) :
    s^2*bb + s*(bc+cb) + cc = 4*s*(s*b+c) := by
  rw [hbb,hbc,hcb,hcc]
  linear_combination (-2*b) * root

/-- The root-gated dark effect is bounded by the root weight even for
coherent superpositions entangled with an arbitrary finite reference. -/
theorem loop_dark_le_root {ι : Type*} (s : Finset ι)
    (u v x y : ι → ℝ) :
    (∑ j ∈ s, ((u j - v j)^2 + (x j - y j)^2)/2) ≤
    ∑ j ∈ s, (u j^2 + v j^2 + x j^2 + y j^2) := by
  apply Finset.sum_le_sum
  intro j _
  nlinarith [sq_nonneg (u j + v j), sq_nonneg (x j + y j)]

end IntervalEcho.GatewayMeter

