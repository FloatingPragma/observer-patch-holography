import Mathlib
import PrimitivePortScaleBoundary

namespace OPH.PrimitivePortMetricAttachment

/-!
# Primitive-port metric attachment boundary

The oriented icosahedral packet has three distinguished transitive incidence
orbits: twelve ports, twenty faces, and thirty edges.  If a physical metric
packet identifies one of those orbit elements with one equal-area cell on a
sphere whose areal radius is the primitive translation hop `a`, then the orbit
cardinality fixes an exact coefficient relating `a^2` to the cell area.

This file proves that conditional arithmetic and the resulting ambiguity.  It
does not identify a screen cell with a port, face, or edge sector; prove that
the translation hop is the areal radius; construct a quotient-visible metric
readout; or select any of the three attachments.  Equal trace or equal
probability on the port projections is not, by itself, an area-measure
identification.
-/

/-- Scalar shadow of the missing metric attachment.  `orbitCard` equal-area
cells cover a round shell of areal radius `a`.  The equation is a premise,
not a consequence of the finite incidence packet. -/
structure EqualOrbitShellAttachment
    (orbitCard : ℕ) (a cellArea : ℝ) : Prop where
  hopPositive : 0 < a
  cellAreaPositive : 0 < cellArea
  shellAreaPartition :
    (orbitCard : ℝ) * cellArea = 4 * Real.pi * a ^ 2

/-- An equal-area orbit attachment fixes the squared hop in terms of the cell
area and the orbit cardinality. -/
theorem hop_sq_of_equal_orbit_shell
    (orbitCard : ℕ) (a cellArea : ℝ)
    (h : EqualOrbitShellAttachment orbitCard a cellArea) :
    a ^ 2 = ((orbitCard : ℝ) / (4 * Real.pi)) * cellArea := by
  have hpi : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  field_simp [hpi]
  nlinarith [h.shellAreaPartition]

/-- If the cell area is `P * ell^2`, the orbit attachment gives the
corresponding conditional metric coefficient. -/
theorem metric_relation_of_equal_orbit_shell
    (orbitCard : ℕ) (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment orbitCard a cellArea)
    (hcell : cellArea = P * ell ^ 2) :
    a ^ 2 = ((orbitCard : ℝ) / (4 * Real.pi)) * P * ell ^ 2 := by
  rw [hop_sq_of_equal_orbit_shell orbitCard a cellArea h, hcell]
  ring

/-- Port-sector attachment: twelve equal cells give `kappa = 3/pi`. -/
theorem port_orbit_metric_relation
    (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment 12 a cellArea)
    (hcell : cellArea = P * ell ^ 2) :
    a ^ 2 = (3 / Real.pi) * P * ell ^ 2 := by
  rw [metric_relation_of_equal_orbit_shell 12 a cellArea P ell h hcell]
  ring

/-- On the port-sector attachment, substituting the selected metric relation
into the frozen coefficient packet gives exact `P`- and `ell`-dependent
coefficients.  The result remains conditional on the metric attachment and on
the primitive-port branch. -/
theorem port_orbit_frozen_coefficient_consequences
    (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment 12 a cellArea)
    (hcell : cellArea = P * ell ^ 2) :
    OPH.PrimitivePortScaleBoundary.scaleC4 a =
        -(3 / (20 * Real.pi)) * P * ell ^ 2 ∧
      OPH.PrimitivePortScaleBoundary.scaleB0 a =
        (3 / (280 * Real.pi ^ 2)) * P ^ 2 * ell ^ 4 ∧
      OPH.PrimitivePortScaleBoundary.scaleB6 a =
        (2 / (875 * Real.pi ^ 2)) * P ^ 2 * ell ^ 4 := by
  have hpi : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  have hscale := port_orbit_metric_relation a cellArea P ell h hcell
  constructor
  · unfold OPH.PrimitivePortScaleBoundary.scaleC4
    rw [hscale]
    field_simp [hpi]
  constructor
  · unfold OPH.PrimitivePortScaleBoundary.scaleB0
    rw [show a ^ 4 = (a ^ 2) ^ 2 by ring, hscale]
    field_simp [hpi]
    ring
  · unfold OPH.PrimitivePortScaleBoundary.scaleB6
    rw [show a ^ 4 = (a ^ 2) ^ 2 by ring, hscale]
    field_simp [hpi]
    ring

/-- Face-sector attachment: twenty equal cells give `kappa = 5/pi`. -/
theorem face_orbit_metric_relation
    (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment 20 a cellArea)
    (hcell : cellArea = P * ell ^ 2) :
    a ^ 2 = (5 / Real.pi) * P * ell ^ 2 := by
  rw [metric_relation_of_equal_orbit_shell 20 a cellArea P ell h hcell]
  ring

/-- Edge-sector attachment: thirty equal cells give `kappa = 15/(2*pi)`. -/
theorem edge_orbit_metric_relation
    (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment 30 a cellArea)
    (hcell : cellArea = P * ell ^ 2) :
    a ^ 2 = (15 / (2 * Real.pi)) * P * ell ^ 2 := by
  rw [metric_relation_of_equal_orbit_shell 30 a cellArea P ell h hcell]
  ring

/-- The three incidence-orbit metric coefficients are pairwise distinct.
Symmetry and the orbit census therefore do not choose one without a typed
cell-to-orbit attachment. -/
theorem incidence_orbit_metric_coefficients_pairwise_distinct :
    (3 / Real.pi : ℝ) ≠ 5 / Real.pi ∧
      (3 / Real.pi : ℝ) ≠ 15 / (2 * Real.pi) ∧
      (5 / Real.pi : ℝ) ≠ 15 / (2 * Real.pi) := by
  have hpi : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  constructor
  · intro h
    field_simp [hpi] at h
    norm_num at h
  constructor
  · intro h
    field_simp [hpi] at h
    norm_num at h
  · intro h
    field_simp [hpi] at h
    norm_num at h

/-- On the port-sector attachment, the dimensionless hop ratio is fixed by
the positive square root.  This remains conditional on the attachment. -/
theorem port_orbit_hop_ratio
    (a cellArea P ell : ℝ)
    (h : EqualOrbitShellAttachment 12 a cellArea)
    (hcell : cellArea = P * ell ^ 2)
    (hell : 0 < ell) (hP : 0 ≤ P) :
    a / ell = Real.sqrt ((3 / Real.pi) * P) := by
  apply OPH.PrimitivePortScaleBoundary.metric_ratio_eq_sqrt
      a (3 / Real.pi) P ell (le_of_lt h.hopPositive) hell
  · positivity
  · exact port_orbit_metric_relation a cellArea P ell h hcell

/-- A certified interval for `P` propagates to an exact squared-hop interval
on the port-sector attachment.  Without the attachment this is not a physical
bound on the primitive translation scale. -/
theorem port_orbit_squared_hop_ratio_interval
    (a cellArea P ell Plo Phi : ℝ)
    (h : EqualOrbitShellAttachment 12 a cellArea)
    (hcell : cellArea = P * ell ^ 2)
    (hell : ell ≠ 0) (hlo : Plo ≤ P) (hhi : P ≤ Phi) :
    (3 / Real.pi) * Plo ≤ (a / ell) ^ 2 ∧
      (a / ell) ^ 2 ≤ (3 / Real.pi) * Phi := by
  have hscale := port_orbit_metric_relation a cellArea P ell h hcell
  have hratio :=
    OPH.PrimitivePortScaleBoundary.metric_ratio_of_scale_relation
      a (3 / Real.pi) P ell hell hscale
  have hkappa : 0 ≤ (3 / Real.pi : ℝ) := by positivity
  rw [hratio]
  exact ⟨mul_le_mul_of_nonneg_left hlo hkappa,
    mul_le_mul_of_nonneg_left hhi hkappa⟩

#print axioms hop_sq_of_equal_orbit_shell
#print axioms metric_relation_of_equal_orbit_shell
#print axioms port_orbit_metric_relation
#print axioms port_orbit_frozen_coefficient_consequences
#print axioms face_orbit_metric_relation
#print axioms edge_orbit_metric_relation
#print axioms incidence_orbit_metric_coefficients_pairwise_distinct
#print axioms port_orbit_hop_ratio
#print axioms port_orbit_squared_hop_ratio_interval

end OPH.PrimitivePortMetricAttachment
