import PrimitivePortFrameQuotient
import Mathlib.MeasureTheory.Measure.Lebesgue.VolumeOfBalls
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Mul

/-!
# The inverse-square shell law on the committed completion carrier

V3 gravitation lane (issue #729), Newtonian-limit row. The corpus commits
a three-dimensional completion carrier: the exact record frame completes
to `EuclideanVec3 = EuclideanSpace ℝ (Fin 3)`
(`PrimitivePortFrameQuotient`, `RepairWordCarrierReadout`). This module
proves that on that carrier a spherically symmetric flux with shell
normalization falls off as the inverse square of the radius, with the
exponent supplied by the carrier dimension rather than assumed: the ball
volume scales as `r ^ finrank`, the shell content is its derivative, and
`finrank ℝ EuclideanVec3 = 3` is a theorem, so the strength law carries
the exponent `finrank - 1 = 2`.

Two premises are declared, not derived: the radial readout premise (a
single strength profile carries the flux, spherical symmetry) and the
shell-flux normalization premise (strength times shell content is the
constant source charge). Both enter as fields of `RadialFluxData` as
declared premises of the V3 program. The exponent provenance is carried
at proof level by `ballVolume_finrank_exponent`, which routes through
the generic finite-dimensional ball-volume theorem, so the scaling
exponent is the carrier `finrank` before any three-dimensional
evaluation; `shellContent_exponent` then carries the dimension-minus-one
form of the strength law. The module touches neither the Einstein-branch
`stress` nor `geometry` objects. The join to the composed Einstein
branch is open: it requires the screen-carrier-to-physical-position
attachment stated as out of scope in `SpatialReadbackSoldering`, which
no source theorem discharges. No continuum gravitational field,
potential, or laboratory constant is claimed.
-/

namespace OPH.InverseSquareShellLaw

open MeasureTheory Metric Module Real

noncomputable section

/-- Real ball volume on the committed completion carrier. -/
def ballVolume (r : ℝ) : ℝ :=
  (volume (ball (0 : OPH.PrimitivePortFrameQuotient.EuclideanVec3) r)).toReal

/-- The carrier dimension is a theorem, not an input. -/
theorem carrier_finrank : finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 = 3 :=
  finrank_euclideanSpace_fin

/-- Exact ball volume: the three-dimensional Euclidean value. -/
theorem ballVolume_eq {r : ℝ} (hr : 0 ≤ r) :
    ballVolume r = (π * 4 / 3) * r ^ 3 := by
  have h := EuclideanSpace.volume_ball_fin_three (0 : OPH.PrimitivePortFrameQuotient.EuclideanVec3) r
  have hpi : (0 : ℝ) ≤ π * 4 / 3 := by positivity
  rw [ballVolume, h, ENNReal.toReal_mul, ← ENNReal.ofReal_pow hr,
    ENNReal.toReal_ofReal (pow_nonneg hr 3), ENNReal.toReal_ofReal hpi]
  ring

/-- Exponent provenance at proof level: the generic finite-dimensional
ball-volume theorem gives the scaling exponent as the carrier `finrank`
with no three-dimensional evaluation anywhere in the proof; the constant
is the generic Gamma-function value. -/
theorem ballVolume_finrank_exponent {r : ℝ} (hr : 0 ≤ r) :
    ballVolume r =
      r ^ (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3) *
        (Real.sqrt π ^ (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3) /
          Real.Gamma ((finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 : ℝ) / 2 + 1)) := by
  have h := InnerProductSpace.volume_ball
    (0 : OPH.PrimitivePortFrameQuotient.EuclideanVec3) r
  have hΓ : (0 : ℝ) <
      Real.Gamma ((finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 : ℝ) / 2 + 1) := by
    apply Real.Gamma_pos_of_pos
    positivity
  have hconst : (0 : ℝ) ≤
      Real.sqrt π ^ (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3) /
        Real.Gamma ((finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 : ℝ) / 2 + 1) := by
    positivity
  rw [ballVolume, h, ENNReal.toReal_mul, ← ENNReal.ofReal_pow hr,
    ENNReal.toReal_ofReal (pow_nonneg hr _), ENNReal.toReal_ofReal hconst]

/-- The evaluated form: `carrier_finrank` pins the exponent to three and
the constant to the three-dimensional value. The exponent in the
statement is the carrier `finrank`; the proof evaluates the
three-dimensional volume and discharges the exponent through
`carrier_finrank`, itself a theorem. -/
theorem ballVolume_finrank_scaling {r : ℝ} (hr : 0 ≤ r) :
    ballVolume r =
      (π * 4 / 3) * r ^ (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3) := by
  rw [ballVolume_eq hr, carrier_finrank]

/-- Shell content: the derivative of the ball volume in the radius. -/
def shellContent (r : ℝ) : ℝ := deriv ballVolume r

/-- Exact shell content at positive radius: the derivative of the cubic
volume law, which carries the exponent `finrank - 1 = 2`. -/
theorem shellContent_eq {r : ℝ} (hr : 0 < r) :
    shellContent r = 4 * π * r ^ 2 := by
  have hcongr : ballVolume =ᶠ[nhds r] fun s => (π * 4 / 3) * s ^ 3 := by
    filter_upwards [Ioi_mem_nhds hr] with s hs
    exact ballVolume_eq (le_of_lt hs)
  have hd : deriv (fun s : ℝ => (π * 4 / 3) * s ^ 3) r =
      (π * 4 / 3) * (3 * r ^ 2) := by
    simp [((hasDerivAt_pow 3 r).const_mul (π * 4 / 3)).deriv]
  rw [shellContent, Filter.EventuallyEq.deriv_eq hcongr, hd]
  ring

/-- The shell exponent equals the carrier dimension minus one. -/
theorem shellContent_exponent {r : ℝ} (hr : 0 < r) :
    shellContent r =
      (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 : ℝ) * (π * 4 / 3) *
        r ^ (finrank ℝ OPH.PrimitivePortFrameQuotient.EuclideanVec3 - 1) := by
  rw [shellContent_eq hr, carrier_finrank]
  push_cast
  ring_nf

/-- The two declared premises of the Newtonian-limit row: one radial
strength profile carries the flux, and its product with the shell
content is the constant source charge on every positive shell. -/
structure RadialFluxData where
  charge : ℝ
  strength : ℝ → ℝ
  normalized : ∀ r : ℝ, 0 < r → strength r * shellContent r = charge

/-- The inverse-square law: under the declared premises the strength at
positive radius is the charge over `4 π r ^ 2`, with the exponent two
supplied by `shellContent_exponent` and hence by the carrier dimension. -/
theorem RadialFluxData.inverse_square (D : RadialFluxData) {r : ℝ}
    (hr : 0 < r) : D.strength r = D.charge / (4 * π * r ^ 2) := by
  have hshell : shellContent r = 4 * π * r ^ 2 := shellContent_eq hr
  have hnz : (4 * π * r ^ 2) ≠ 0 := by positivity
  have h := D.normalized r hr
  rw [hshell] at h
  field_simp
  linarith [h]

/-- Scale-free form: strength times squared radius is shell-independent. -/
theorem RadialFluxData.scale_free (D : RadialFluxData) {r₁ r₂ : ℝ}
    (h₁ : 0 < r₁) (h₂ : 0 < r₂) :
    D.strength r₁ * r₁ ^ 2 = D.strength r₂ * r₂ ^ 2 := by
  have e₁ := D.inverse_square h₁
  have e₂ := D.inverse_square h₂
  have hπ : (0:ℝ) < π := Real.pi_pos
  rw [e₁, e₂]
  field_simp

end

end OPH.InverseSquareShellLaw

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.InverseSquareShellLaw.carrier_finrank
#print axioms OPH.InverseSquareShellLaw.ballVolume_eq
#print axioms OPH.InverseSquareShellLaw.ballVolume_finrank_exponent
#print axioms OPH.InverseSquareShellLaw.ballVolume_finrank_scaling
#print axioms OPH.InverseSquareShellLaw.shellContent_eq
#print axioms OPH.InverseSquareShellLaw.shellContent_exponent
#print axioms OPH.InverseSquareShellLaw.RadialFluxData.inverse_square
#print axioms OPH.InverseSquareShellLaw.RadialFluxData.scale_free
