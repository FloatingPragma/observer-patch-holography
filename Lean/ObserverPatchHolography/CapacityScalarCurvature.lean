import Mathlib.Tactic

/-! Conditional canonical-scalar curvature bound for the declared capacity history.
The source-to-scalar physical map is not assumed proved. The polynomial identity
is the exact finite consequence of the conservation-based reconstruction. -/
namespace OPH.CapacityScalarCurvature

noncomputable def massRatio (x r : ℝ) : ℝ :=
    -3*r^2*x^4/16-r^2*x^3/8+5*r^2*x^2/16+r^2*x/4-r*x^5
    -9*r*x^4/8+r*x^3-r*x^2/8-r*x/4+r-5*x^6/4-2*x^5
    +21*x^4/16+9*x^3/8-31*x^2/16+5

theorem lower (x r : ℝ) (hx : 0 ≤ x) (hx1 : x ≤ 1)
    (hr : 0 ≤ r) (hr1 : r ≤ 1) : 2 ≤ massRatio x r := by
  have hxc : 0 ≤ 1-x := sub_nonneg.mpr hx1
  have hrc : 0 ≤ 1-r := sub_nonneg.mpr hr1
  have hid : massRatio x r - 2 =
    1 * x^0 * (1-x)^6 * ((3 / 1 : ℝ)*(1-r)^2 + 2*(7 / 2 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) +
    6 * x^1 * (1-x)^5 * ((3 / 1 : ℝ)*(1-r)^2 + 2*(167 / 48 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) +
    15 * x^2 * (1-x)^4 * ((689 / 240 : ℝ)*(1-r)^2 + 2*(133 / 40 : ℝ)*r*(1-r) + (233 / 60 : ℝ)*r^2) +
    20 * x^3 * (1-x)^3 * ((427 / 160 : ℝ)*(1-r)^2 + 2*(499 / 160 : ℝ)*r*(1-r) + (15 / 4 : ℝ)*r^2) +
    15 * x^4 * (1-x)^2 * ((203 / 80 : ℝ)*(1-r)^2 + 2*(359 / 120 : ℝ)*r*(1-r) + (37 / 10 : ℝ)*r^2) +
    6 * x^5 * (1-x)^1 * ((19 / 8 : ℝ)*(1-r)^2 + 2*(65 / 24 : ℝ)*r*(1-r) + (10 / 3 : ℝ)*r^2) +
    1 * x^6 * (1-x)^0 * ((1 / 4 : ℝ)*(1-r)^2 + 2*(0 / 1 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) := by
    unfold massRatio
    ring
  have hnonneg : 0 ≤
    1 * x^0 * (1-x)^6 * ((3 / 1 : ℝ)*(1-r)^2 + 2*(7 / 2 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) +
    6 * x^1 * (1-x)^5 * ((3 / 1 : ℝ)*(1-r)^2 + 2*(167 / 48 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) +
    15 * x^2 * (1-x)^4 * ((689 / 240 : ℝ)*(1-r)^2 + 2*(133 / 40 : ℝ)*r*(1-r) + (233 / 60 : ℝ)*r^2) +
    20 * x^3 * (1-x)^3 * ((427 / 160 : ℝ)*(1-r)^2 + 2*(499 / 160 : ℝ)*r*(1-r) + (15 / 4 : ℝ)*r^2) +
    15 * x^4 * (1-x)^2 * ((203 / 80 : ℝ)*(1-r)^2 + 2*(359 / 120 : ℝ)*r*(1-r) + (37 / 10 : ℝ)*r^2) +
    6 * x^5 * (1-x)^1 * ((19 / 8 : ℝ)*(1-r)^2 + 2*(65 / 24 : ℝ)*r*(1-r) + (10 / 3 : ℝ)*r^2) +
    1 * x^6 * (1-x)^0 * ((1 / 4 : ℝ)*(1-r)^2 + 2*(0 / 1 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) := by positivity
  linarith

theorem upper (x r : ℝ) (hx : 0 ≤ x) (hx1 : x ≤ 1)
    (hr : 0 ≤ r) (hr1 : r ≤ 1) : massRatio x r ≤ 6 := by
  have hxc : 0 ≤ 1-x := sub_nonneg.mpr hx1
  have hrc : 0 ≤ 1-r := sub_nonneg.mpr hr1
  have hid : 6 - massRatio x r =
    1 * x^0 * (1-x)^6 * ((1 / 1 : ℝ)*(1-r)^2 + 2*(1 / 2 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) +
    6 * x^1 * (1-x)^5 * ((1 / 1 : ℝ)*(1-r)^2 + 2*(25 / 48 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) +
    15 * x^2 * (1-x)^4 * ((271 / 240 : ℝ)*(1-r)^2 + 2*(27 / 40 : ℝ)*r*(1-r) + (7 / 60 : ℝ)*r^2) +
    20 * x^3 * (1-x)^3 * ((213 / 160 : ℝ)*(1-r)^2 + 2*(141 / 160 : ℝ)*r*(1-r) + (1 / 4 : ℝ)*r^2) +
    15 * x^4 * (1-x)^2 * ((117 / 80 : ℝ)*(1-r)^2 + 2*(121 / 120 : ℝ)*r*(1-r) + (3 / 10 : ℝ)*r^2) +
    6 * x^5 * (1-x)^1 * ((13 / 8 : ℝ)*(1-r)^2 + 2*(31 / 24 : ℝ)*r*(1-r) + (2 / 3 : ℝ)*r^2) +
    1 * x^6 * (1-x)^0 * ((15 / 4 : ℝ)*(1-r)^2 + 2*(4 / 1 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) := by
    unfold massRatio
    ring
  have hnonneg : 0 ≤
    1 * x^0 * (1-x)^6 * ((1 / 1 : ℝ)*(1-r)^2 + 2*(1 / 2 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) +
    6 * x^1 * (1-x)^5 * ((1 / 1 : ℝ)*(1-r)^2 + 2*(25 / 48 : ℝ)*r*(1-r) + (0 / 1 : ℝ)*r^2) +
    15 * x^2 * (1-x)^4 * ((271 / 240 : ℝ)*(1-r)^2 + 2*(27 / 40 : ℝ)*r*(1-r) + (7 / 60 : ℝ)*r^2) +
    20 * x^3 * (1-x)^3 * ((213 / 160 : ℝ)*(1-r)^2 + 2*(141 / 160 : ℝ)*r*(1-r) + (1 / 4 : ℝ)*r^2) +
    15 * x^4 * (1-x)^2 * ((117 / 80 : ℝ)*(1-r)^2 + 2*(121 / 120 : ℝ)*r*(1-r) + (3 / 10 : ℝ)*r^2) +
    6 * x^5 * (1-x)^1 * ((13 / 8 : ℝ)*(1-r)^2 + 2*(31 / 24 : ℝ)*r*(1-r) + (2 / 3 : ℝ)*r^2) +
    1 * x^6 * (1-x)^0 * ((15 / 4 : ℝ)*(1-r)^2 + 2*(4 / 1 : ℝ)*r*(1-r) + (4 / 1 : ℝ)*r^2) := by positivity
  linarith

#print axioms lower
#print axioms upper
end OPH.CapacityScalarCurvature

