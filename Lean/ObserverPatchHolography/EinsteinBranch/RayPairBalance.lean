import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-! Null contractions of a symmetric stress-curvature residual on a time-space plane.
The time-space entry carries the oriented flux. Each read uses the same frame. -/
namespace OPH.EinsteinBranch.RayPairBalance

def contraction (d j p t x : ℝ) : ℝ := d*t^2 + 2*j*t*x + p*x^2

/-- A rest read and two future unit timelike reads give the outgoing null charge.
Reversing the spatial axis gives the incoming charge. -/
theorem from_three_observers (d j p : ℝ) :
    contraction d j p 1 1 =
      -(9/16)*contraction d j p 1 0 +
      (81/160)*contraction d j p (5/3) (4/3) +
      (9/160)*contraction d j p (5/3) (-4/3) := by
  unfold contraction
  ring

theorem null_error_bound (d j p ε : ℝ)
    (h₀ : |contraction d j p 1 0| ≤ ε)
    (h₁ : |contraction d j p (5/3) (4/3)| ≤ ε)
    (h₂ : |contraction d j p (5/3) (-4/3)| ≤ ε) :
    |contraction d j p 1 1| ≤ (9/8)*ε := by
  rw [from_three_observers]
  obtain ⟨a,b⟩ := abs_le.mp h₀
  obtain ⟨c,d'⟩ := abs_le.mp h₁
  obtain ⟨e,f⟩ := abs_le.mp h₂
  apply abs_le.mpr
  constructor <;> linarith

/-- Two states with equal geometric flux and opposite matter flux. Even
stress terms of the two states can differ arbitrarily. -/
theorem opposite_flux_lower_bound (evenA evenB geom current ε : ℝ)
    (ap : |evenA + 2*(geom-current)| ≤ ε) (am : |evenA - 2*(geom-current)| ≤ ε)
    (bp : |evenB + 2*(geom+current)| ≤ ε) (bm : |evenB - 2*(geom+current)| ≤ ε) :
    2*|current| ≤ ε := by
  obtain ⟨a,b⟩ := abs_le.mp ap
  obtain ⟨c,d⟩ := abs_le.mp am
  obtain ⟨e,f⟩ := abs_le.mp bp
  obtain ⟨g,h⟩ := abs_le.mp bm
  rcases le_total 0 current with pos | neg
  · rw [abs_of_nonneg pos]
    linarith
  · rw [abs_of_nonpos neg]
    linarith

/-- With equal number currents a common shift of the energy zero cancels
between the two sources. The opposite ray readings remove the even term
separately on each source. -/
theorem matched_number_current_bound (evenA evenB gA gB jA jB nA nB k c ε : ℝ)
    (matched : nA = nB)
    (ap : |evenA + 2 * (gA - k * (jA + c * nA))| ≤ ε)
    (am : |evenA - 2 * (gA - k * (jA + c * nA))| ≤ ε)
    (bp : |evenB + 2 * (gB - k * (jB + c * nB))| ≤ ε)
    (bm : |evenB - 2 * (gB - k * (jB + c * nB))| ≤ ε) :
    |k * (jA - jB)| ≤ ε + |gA - gB| := by
  subst nB
  obtain ⟨apLo, apHi⟩ := abs_le.mp ap
  obtain ⟨amLo, amHi⟩ := abs_le.mp am
  obtain ⟨bpLo, bpHi⟩ := abs_le.mp bp
  obtain ⟨bmLo, bmHi⟩ := abs_le.mp bm
  obtain ⟨gLo, gHi⟩ := abs_le.mp (le_refl |gA - gB|)
  apply abs_le.mpr
  constructor
  · linarith only [gLo, apHi, amLo, bpLo, bmHi]
  · linarith only [gHi, apLo, amHi, bpHi, bmLo]

/-- For a diagonal homogeneous comparison the all-sky null residual is a
convex combination of its three axial values. The axes attain its extrema. -/
theorem diagonal_all_sky (a b c x y z ε : ℝ)
    (unit : x^2+y^2+z^2=1) (ha : |a| ≤ ε) (hb : |b| ≤ ε) (hc : |c| ≤ ε) :
    |a*x^2+b*y^2+c*z^2| ≤ ε := by
  obtain ⟨aLo,aHi⟩ := abs_le.mp ha
  obtain ⟨bLo,bHi⟩ := abs_le.mp hb
  obtain ⟨cLo,cHi⟩ := abs_le.mp hc
  have s := congrArg (fun t : ℝ => ε*t) unit
  dsimp at s
  apply abs_le.mpr
  constructor
  · nlinarith [mul_le_mul_of_nonneg_right aLo (sq_nonneg x),
      mul_le_mul_of_nonneg_right bLo (sq_nonneg y), mul_le_mul_of_nonneg_right cLo (sq_nonneg z)]
  · nlinarith [mul_le_mul_of_nonneg_right aHi (sq_nonneg x),
      mul_le_mul_of_nonneg_right bHi (sq_nonneg y), mul_le_mul_of_nonneg_right cHi (sq_nonneg z)]

/-- Rate derivative obtained from the canonical momentum equations of the
finite-mode semiclassical box. C is the sum of the three pair products of
logarithmic expansion rates, P the sum of directional pressures. -/
noncomputable def canonicalSlope (θ C Λ κ h p P : ℝ) : ℝ :=
  (C+Λ)/2 + κ*(p-P/2) - θ*h

/-- Hamilton's momentum update and the derivative of its volume prefactor
give the declared rate slope. The three forces use the same area-free energy
and positive volume calibration as the source Hamiltonian. -/
theorem canonical_from_momenta (θ C Λ κ h p q r : ℝ) :
    (-C-Λ+κ*p) - ((-C-Λ+κ*p)+(-C-Λ+κ*q)+(-C-Λ+κ*r))/2 - θ*h =
      canonicalSlope θ C Λ κ h p (p+q+r) := by
  unfold canonicalSlope
  ring

def timeCurvature (x y z : ℝ) : ℝ := x*y+x*z+y*z

/-- The measured spatial curvature on the first axis follows from Hamilton's
equations alone. Cyclic permutation gives the other two components. -/
theorem canonical_spatial (x y z p q r Λ κ : ℝ) :
    -(canonicalSlope (x+y+z) (timeCurvature x y z) Λ κ y q (p+q+r) + y^2)
    -(canonicalSlope (x+y+z) (timeCurvature x y z) Λ κ z r (p+q+r) + z^2)
    -y*z = κ*p-Λ := by
  unfold canonicalSlope timeCurvature
  ring

/-- Derivative of the Hamiltonian constraint defect on the canonical flow,
using each mode's directional stress continuity equation. -/
theorem canonical_constraint_rate (x y z p q r ρ Λ κ : ℝ) :
    (y+z)*canonicalSlope (x+y+z) (timeCurvature x y z) Λ κ x p (p+q+r)
    + (x+z)*canonicalSlope (x+y+z) (timeCurvature x y z) Λ κ y q (p+q+r)
    + (x+y)*canonicalSlope (x+y+z) (timeCurvature x y z) Λ κ z r (p+q+r)
    + κ*(x*(ρ+p)+y*(ρ+q)+z*(ρ+r)) =
    -(x+y+z)*(timeCurvature x y z-Λ-κ*ρ) := by
  unfold canonicalSlope timeCurvature
  ring

/-- On each regular proper-time interval the product V D is constant.
Zero defect on the initial slice gives zero defect at the final slice. -/
theorem propagate_constraint (vol defect expansion : ℝ → ℝ) (a b : ℝ) (hab : a < b)
    (continuous : ContinuousOn (fun t => vol t * defect t) (Set.Icc a b))
    (rateV : ∀ t ∈ Set.Ioo a b, HasDerivAt vol (expansion t * vol t) t)
    (rateD : ∀ t ∈ Set.Ioo a b, HasDerivAt defect (-expansion t * defect t) t)
    (initial : defect a = 0) (positive : 0 < vol b) : defect b = 0 := by
  have rate : ∀ t ∈ Set.Ioo a b, HasDerivAt (fun t => vol t * defect t) 0 t := by
    intro t ht
    apply ((rateV t ht).mul (rateD t ht)).congr_deriv
    ring
  obtain ⟨_, _, hc⟩ := exists_hasDerivAt_eq_slope _ (fun _ => (0:ℝ)) hab continuous rate
  rw [initial, mul_zero, sub_zero] at hc
  exact (mul_eq_zero.mp ((div_eq_zero_iff.mp hc.symm).resolve_right
    (sub_ne_zero.mpr hab.ne'))).resolve_left positive.ne'

/-- With diagonal gradient and stress and initial constraint propagated, every
normalized null vector reads exactly the same coupled energy and pressure. -/
theorem canonical_null_balance (C ρ px py pz Λ κ x y z : ℝ)
    (constraint : C = Λ+κ*ρ) (unit : x^2+y^2+z^2=1) :
    C + (κ*px-Λ)*x^2 + (κ*py-Λ)*y^2 + (κ*pz-Λ)*z^2 =
      κ*(ρ+px*x^2+py*y^2+pz*z^2) := by
  rw [constraint]
  linear_combination -Λ*unit

/-- In an orthonormal frame with signature (-,+,+,+), the time component
of the Ricci tensor is half the sum of the four diagonal Einstein components.
Substitute the four components of E = G + Λ g − κ T before taking a local
Newtonian limit. The vacuum term has the opposite sign to ordinary density. -/
theorem time_trace_reversal (g₀ g₁ g₂ g₃ ρ p₁ p₂ p₃ Λ κ : ℝ) :
    (g₀+g₁+g₂+g₃)/2 = κ*(ρ+p₁+p₂+p₃)/2 - Λ +
      ((g₀-Λ-κ*ρ)+(g₁+Λ-κ*p₁)+(g₂+Λ-κ*p₂)+(g₃+Λ-κ*p₃))/2 := by
  ring

/-- The trace-reversed residual has sharp componentwise upper bound 2ε. -/
theorem time_trace_residual_bound (a b c d ε : ℝ) (ha : |a| ≤ ε) (hb : |b| ≤ ε)
    (hc : |c| ≤ ε) (hd : |d| ≤ ε) : |(a+b+c+d)/2| ≤ 2*ε := by
  obtain ⟨aLo,aHi⟩ := abs_le.mp ha
  obtain ⟨bLo,bHi⟩ := abs_le.mp hb
  obtain ⟨cLo,cHi⟩ := abs_le.mp hc
  obtain ⟨dLo,dHi⟩ := abs_le.mp hd
  apply abs_le.mpr
  constructor <;> linarith

end OPH.EinsteinBranch.RayPairBalance
