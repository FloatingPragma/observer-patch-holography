import PrimitivePortMetricAttachment
import SeamCurrentDirichletGenerator

namespace OPH.SeamCurrentPhysicalMetricAttachment

open OPH.SeamCurrentHomogeneousAction (DirectedSeam)
open OPH.SeamCurrentDirichletGenerator

/-!
# Same-metric physical attachment for the two propagation actions

The source packet fixes a normalized port frame and the exact norm of every
unit-current seam event in that frame.  The port-dual packet separately fixes
normalized weight `1/12`.  Neither normalized statement identifies the UV
cut-area density `aCell` with a port-dual spherical sector.  In particular,
`aCell` is not a preferred spatial ruler merely because the two normalized
weights agree.

This file records the conditional algebra available if one joint physical
attachment identifies all three objects: the normalized port frame, the
primitive vertex action, and the seam-current action.  The same attachment
must also identify one port-dual sector with the physical UV cut-area
element.  Under that typed same-object premise, the vertex and seam scales
are no longer independent: both coefficients follow exactly from `aCell` and
the response-Gram seam norm.

The structure below is an explicit physical-identification record, not an
A1--A3 theorem or the output of a source theorem.  No current normalized
incidence, state-weight, or Gram theorem constructs an inhabitant.  It
therefore supplies a precise conditional route for issue #664 without
promoting a physical scale from the current source packet.
-/

/-- An explicit joint attachment of one physical similarity metric to the
normalized support frame and both propagation actions.

`scaleSq` is the squared physical length represented by one unit of the
normalized response metric.  `vertexActionSq` and `seamActionSq` are squared
lengths for the exact actions used by the vertex-12 and seam-edge-30 branches.
The port-dual area equation is an explicit physical attachment: it is not
inferred from equal normalized weights. -/
structure PortDualSameMetricAttachment
    (aCell P ell scaleSq vertexActionSq : ℝ)
    (seamActionSq : DirectedSeam → ℝ) : Prop where
  scaleSqPositive : 0 < scaleSq
  cellAreaPositive : 0 < aCell
  ellPositive : 0 < ell
  /-- One physical UV cut-area element is identified with one of the twelve
  equal port-dual sectors of the same round metric shell. -/
  portDualPhysicalArea : 12 * aCell = 4 * Real.pi * scaleSq
  /-- Definition of the dimensionless pixel coordinate on this attachment. -/
  cellScale : aCell = P * ell ^ 2
  /-- The vertex-12 action uses a unit port generator in this same metric. -/
  vertexUsesSameMetric : vertexActionSq = scaleSq
  /-- The edge-30 action uses the exact unit-current seam displacement in
  this same metric. -/
  seamUsesSameMetric : ∀ e,
    seamActionSq e = scaleSq *
      OPH.PrimitivePortTranslationBridge.dot
        (completionSeamStep e) (completionSeamStep e)

/-- The attached port-dual cell fixes the one common metric conversion.
This is the conditional `kappa_vertex = 3/pi` relation, expressed without
identifying an area with a length. -/
theorem common_scale_sq_eq_three_over_pi_mul_cell
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq) :
    scaleSq = (3 / Real.pi) * aCell := by
  have hpi : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  field_simp [hpi]
  nlinarith [h.portDualPhysicalArea]

/-- The primitive vertex action has the port-sector coefficient `3/pi` once
the joint physical attachment is supplied. -/
theorem vertex_action_sq_eq_three_over_pi_mul_cell
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq) :
    vertexActionSq = (3 / Real.pi) * aCell := by
  rw [h.vertexUsesSameMetric,
    common_scale_sq_eq_three_over_pi_mul_cell h]

/-- Pixel-coordinate form of the primitive vertex scale on the same
attachment. -/
theorem vertex_action_sq_eq_pixel_scale
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq) :
    vertexActionSq = (3 / Real.pi) * P * ell ^ 2 := by
  rw [vertex_action_sq_eq_three_over_pi_mul_cell h, h.cellScale]
  ring

/-- The exact response-Gram theorem fixes the seam action relative to the
same physical metric unit. -/
theorem seam_action_sq_eq_internal_factor_mul_scale
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (e : DirectedSeam) :
    seamActionSq e = (2 - 2 / Real.sqrt 5) * scaleSq := by
  rw [h.seamUsesSameMetric e,
    completionSeamStep_norm_sq_eq_two_sub_two_div_sqrt_five]
  ring

/-- On the conditional shared metric, the seam scale is fixed by the physical
UV cut area with no second metric coefficient. -/
theorem seam_action_sq_eq_cell
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (e : DirectedSeam) :
    seamActionSq e =
      ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) * aCell := by
  rw [seam_action_sq_eq_internal_factor_mul_scale h e,
    common_scale_sq_eq_three_over_pi_mul_cell h]
  ring

/-- Equivalent coefficient used by the simulator receipt:
`kappa_edge = 6 (1 - 1/sqrt(5)) / pi`. -/
theorem seam_metric_coefficient_equivalent :
    (3 / Real.pi) * (2 - 2 / Real.sqrt 5) =
      6 * (1 - 1 / Real.sqrt 5) / Real.pi := by
  have hpi : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  have hsqrt : Real.sqrt 5 ≠ 0 := by positivity
  field_simp [hpi, hsqrt]
  ring

/-- The exact ratio between the two attached branch scales is the internal
response-Gram seam norm.  A scale derived for either action transfers to the
other only through the same-object premise represented by `h`. -/
theorem seam_over_vertex_squared_ratio
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (e : DirectedSeam) :
    seamActionSq e =
      (2 - 2 / Real.sqrt 5) * vertexActionSq := by
  rw [h.vertexUsesSameMetric]
  exact seam_action_sq_eq_internal_factor_mul_scale h e

/-- The same-object attachment converts the pixel definition into the exact
dimensionless seam-to-reference-area relation. -/
theorem seam_action_sq_eq_pixel_scale
    {aCell P ell scaleSq vertexActionSq : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (e : DirectedSeam) :
    seamActionSq e =
      ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) * P * ell ^ 2 := by
  rw [seam_action_sq_eq_cell h e, h.cellScale]
  ring

/-- The exact seam metric coefficient is strictly positive. -/
theorem seam_metric_coefficient_pos :
    0 < (3 / Real.pi) * (2 - 2 / Real.sqrt 5) := by
  have hsqrtPos : 0 < Real.sqrt 5 := by positivity
  have hsqrtSq : Real.sqrt 5 ^ 2 = 5 := by norm_num
  have hsqrtGtOne : 1 < Real.sqrt 5 := by nlinarith
  have hdiv : 2 / Real.sqrt 5 < 2 := by
    rw [div_lt_iff₀ hsqrtPos]
    nlinarith
  exact mul_pos (div_pos (by norm_num) Real.pi_pos) (sub_pos.mpr hdiv)

/-- A nonnegative outer detuning coordinate gives the algebraic bound
`P >= goldenRatio` under its defining parametrization.  This statement does
not identify `alpha` with the inner coupling `1 / A_T(P)`, impose the
self-referential closure equation, or complete the source-side endpoint map. -/
theorem nonnegative_outer_detuning_ge_goldenRatio
    {P alpha : ℝ} (halpha : 0 ≤ alpha)
    (hdetuning : P = Real.goldenRatio + alpha * Real.sqrt Real.pi) :
    Real.goldenRatio ≤ P := by
  rw [hdetuning]
  exact le_add_of_nonneg_right
    (mul_nonneg halpha (Real.sqrt_nonneg Real.pi))

/-- Combining the nonnegative outer detuning parametrization with the joint
metric attachment produces a strictly positive algebraic lower bound in
`ell` units.  This is a physical same-action bound only after the attachment
premise is discharged, `alpha` is identified with the inner coupling through
the full closure, and `ell` is independently attached to a physical scale. -/
theorem seam_action_sq_goldenRatio_lower_bound
    {aCell P ell scaleSq vertexActionSq alpha : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (halpha : 0 ≤ alpha)
    (hdetuning : P = Real.goldenRatio + alpha * Real.sqrt Real.pi)
    (e : DirectedSeam) :
    0 < ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) *
        Real.goldenRatio * ell ^ 2 ∧
      ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) *
          Real.goldenRatio * ell ^ 2 ≤ seamActionSq e := by
  let kappa : ℝ := (3 / Real.pi) * (2 - 2 / Real.sqrt 5)
  have hkappa : 0 < kappa := seam_metric_coefficient_pos
  have hP : Real.goldenRatio ≤ P :=
    nonnegative_outer_detuning_ge_goldenRatio halpha hdetuning
  have hsqPos : 0 < ell ^ 2 := sq_pos_of_pos h.ellPositive
  have hscaled :
      Real.goldenRatio * ell ^ 2 ≤ P * ell ^ 2 :=
    mul_le_mul_of_nonneg_right hP (le_of_lt hsqPos)
  constructor
  · exact mul_pos (mul_pos hkappa Real.goldenRatio_pos) hsqPos
  · rw [seam_action_sq_eq_pixel_scale h e]
    simpa only [kappa, mul_assoc] using
      mul_le_mul_of_nonneg_left hscaled (le_of_lt hkappa)

/-- The same nonnegative detuning parametrization gives the corresponding
conditional lower bound for the primitive vertex action. -/
theorem vertex_action_sq_goldenRatio_lower_bound
    {aCell P ell scaleSq vertexActionSq alpha : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (halpha : 0 ≤ alpha)
    (hdetuning : P = Real.goldenRatio + alpha * Real.sqrt Real.pi) :
    0 < (3 / Real.pi) * Real.goldenRatio * ell ^ 2 ∧
      (3 / Real.pi) * Real.goldenRatio * ell ^ 2 ≤ vertexActionSq := by
  have hkappa : 0 < (3 / Real.pi : ℝ) := div_pos (by norm_num) Real.pi_pos
  have hP : Real.goldenRatio ≤ P :=
    nonnegative_outer_detuning_ge_goldenRatio halpha hdetuning
  have hsqPos : 0 < ell ^ 2 := sq_pos_of_pos h.ellPositive
  have hscaled :
      Real.goldenRatio * ell ^ 2 ≤ P * ell ^ 2 :=
    mul_le_mul_of_nonneg_right hP (le_of_lt hsqPos)
  constructor
  · exact mul_pos (mul_pos hkappa Real.goldenRatio_pos) hsqPos
  · rw [vertex_action_sq_eq_pixel_scale h]
    simpa only [mul_assoc] using
      mul_le_mul_of_nonneg_left hscaled (le_of_lt hkappa)

/-- A source-certified lower interval for `P`, together with a lower interval
for the independent physical scale `ell`, yields a positive same-action
squared-length bound.  The current corpus has not constructed the physical
attachment `h` or a source-native SI lower interval for `ell`. -/
theorem seam_action_sq_positive_lower_bound
    {aCell P ell scaleSq vertexActionSq ellMin PMin : ℝ}
    {seamActionSq : DirectedSeam → ℝ}
    (h : PortDualSameMetricAttachment
      aCell P ell scaleSq vertexActionSq seamActionSq)
    (hPMin : 0 < PMin) (hP : PMin ≤ P)
    (hEllMin : 0 < ellMin) (hEll : ellMin ≤ ell)
    (e : DirectedSeam) :
    0 < ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) *
        PMin * ellMin ^ 2 ∧
      ((3 / Real.pi) * (2 - 2 / Real.sqrt 5)) *
          PMin * ellMin ^ 2 ≤ seamActionSq e := by
  let kappa : ℝ := (3 / Real.pi) * (2 - 2 / Real.sqrt 5)
  have hkappa : 0 < kappa := seam_metric_coefficient_pos
  have hPPos : 0 < P := lt_of_lt_of_le hPMin hP
  have hsq : ellMin ^ 2 ≤ ell ^ 2 := by nlinarith
  have hprod : PMin * ellMin ^ 2 ≤ P * ell ^ 2 := by
    exact mul_le_mul hP hsq (sq_nonneg ellMin) (le_of_lt hPPos)
  constructor
  · positivity
  · rw [seam_action_sq_eq_pixel_scale h e]
    simpa only [kappa, mul_assoc] using
      mul_le_mul_of_nonneg_left hprod (le_of_lt hkappa)

#print axioms common_scale_sq_eq_three_over_pi_mul_cell
#print axioms vertex_action_sq_eq_three_over_pi_mul_cell
#print axioms vertex_action_sq_eq_pixel_scale
#print axioms seam_action_sq_eq_internal_factor_mul_scale
#print axioms seam_action_sq_eq_cell
#print axioms seam_metric_coefficient_equivalent
#print axioms seam_over_vertex_squared_ratio
#print axioms seam_action_sq_eq_pixel_scale
#print axioms seam_metric_coefficient_pos
#print axioms nonnegative_outer_detuning_ge_goldenRatio
#print axioms seam_action_sq_goldenRatio_lower_bound
#print axioms vertex_action_sq_goldenRatio_lower_bound
#print axioms seam_action_sq_positive_lower_bound

end OPH.SeamCurrentPhysicalMetricAttachment
