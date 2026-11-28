import ObserverPatchHolography.EinsteinBranch.RecordMass

/-! Finite observer comparisons for two serial contacts and a recorded walk.
The route, clock, and probe coefficients are supplied by the native experiment. -/
namespace OPH.EinsteinBranch.ReceiverGravity

/-- The extra arrival probability at two serial contacts with a common fair tag. -/
theorem serial_excess (x y : ℝ) :
    (x^2 + y^2)/2 - ((x+y)/2)^2 = ((x-y)/2)^2 := by ring

/-- First-hit decomposition and parity on a geometrically returning visit.
Here h is the hitting probability, d the expected visits from the receiver,
g the expected visits from the source, and u the odd-visit probability. -/
theorem odd_visits (g d h u r : ℝ)
    (hvisits : d*(1-r) = 1) (hgreen : g = h*d)
    (hodd : (1+r)*u = h) : (2*d-1)*u = g := by
  have hret : d*r = d-1 := by linarith
  rw [hgreen, ← hodd]
  linear_combination -u*hret

/-- For a common return factor the measured odd-visit field has the same
exterior screened balance as the calibrated Green function. -/
theorem exterior_balance (c α x y z w : ℝ) (hc : c ≠ 0)
    (balance : (1-3*α/4)*(c*x) = α*((c*y)+(c*z)+(c*w))/12) :
    (1-3*α/4)*x = α*(y+z+w)/12 := by
  apply mul_left_cancel₀ hc
  linear_combination balance

/-- The release readings on two material loads cannot be obtained from any
positive mixture of the two ordinary deposits at the enrolled ports. -/
theorem ordinary_cone (x y : ℝ) (hx : 0 ≤ x) :
    (65/4)*x + (65/2)*y ≤ 2*((17/2)*x+(65/4)*y) := by linarith

theorem inward_and_slow (d : ℝ) (hd : 0 < d) :
    0 < (17/4)*d ∧ 0 < (65/4)*d ∧ -24*d < 0 := by
  constructor <;> [linarith ; constructor <;> linarith]

theorem no_common_dose (d x y : ℝ) (hd : 0 < d) (hx : 0 ≤ x)
    (h₁ : (17/4)*d = (17/2)*x+(65/4)*y)
    (h₂ : (65/4)*d = (65/4)*x+(65/2)*y) : False := by
  linarith [ordinary_cone x y]

/-- At an open, triangle-free crossing the signed release curvature is half
of the drop in the measured onsite energy. Changing only passive intracell
marks changes each energy by the same coefficient per mark of a given load. -/
theorem passive_response (q near far : ℝ) :
    ((1/2+8*q^2)*near-(1/2+8*q^2)*far)/2 = (1/4+4*q^2)*(near-far) := by ring

/-- Summing the gateway coherences of a single extended launch preserves
its response coefficient when written contacts pair gateways to private ports.
The real edge weights can have either sign. -/
theorem shell_response {n : ℕ} (q : ℝ) (w d : Fin n → ℝ) :
    (∑ i, w i * ((1/4+4*q^2)*d i)) = (1/4+4*q^2)*(∑ i, w i*d i) := by
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- Mixing in ordinary matter leaves the same residual once its separately
measured contribution is subtracted. A pure passive correction may only
enter the ordinary response cone without subtraction above this dose. -/
theorem combined_cone_requires (d x y : ℝ)
    (h : (65/4)*d+(65/4)*x+(65/2)*y ≤ 2*((17/4)*d+(17/2)*x+(65/4)*y)) :
    31*d ≤ 3*x := by linarith

/-- A positive area response together with a positive bulk record-entropy
response precludes stationarity at this fixed volume for any positive gain. -/
theorem stationarity_obstruction (area entropy gain : ℝ)
    (ha : 0 < area) (hs : 0 < entropy) (hg : 0 ≤ gain) :
    0 < gain*area+entropy := by positivity

end OPH.EinsteinBranch.ReceiverGravity
