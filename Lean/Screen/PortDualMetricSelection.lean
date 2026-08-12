import PrimitivePortDualMeasure
import OrientedFaceInvariantMetric

open scoped BigOperators

/-!
# The pinned uniform port-dual measure point selects the compact family G

Proof-kernel companion to
`code/b14_jacobi/port_dual_metric_selector.certificate.json`.  Two
committed layers meet here.  The dual-measure module proves that equal
normalized weights on the twenty reconstructed faces and the declared
barycentric one-third share give port-dual weight exactly `1/12` at
every port.  The invariant-metric module proves that every
sector-balanced positive carrier metric selects `G` as the unique
nearest classified compact family.  The uniform weight induces, under
the declared diagonal measure-to-metric rule, the scalar carrier inner
product with every sector scale equal to `1/12`, which is
sector-balanced, so the selection theorem applies at the measure point.
This module states that instantiation, together with the exact scaled
distance values at the measure point.  The identification of the
uniform-measure metric with the balanced cone coordinates
`(1/12, 1/12, 1/12, 1/12)` is certified by the companion Python
producer and verifier, and enters here only through the substitution of
the common scale into all three metric slots.

**Boundary.**  The equal-face-weight rule and the barycentric share are
the declared premises of the dual-measure module, and the diagonal
measure-to-metric rule and the nearest-point repair rule are declared
discriminator choices; none is derived from the source.  The compact
locus is conditional on the three named textbook lemmas of the pinned
classification certificate.  The theorems upgrade the committed
metric-conditional selection to a measure-conditional selection whose
measure is pinned by a committed theorem.  No physical current,
coupling, laboratory gauge field, or source bracket selection is
claimed.
-/

namespace OPH.PortDualMetricSelection

open OPH.PrimitivePortDualMeasure
open OPH.OrientedFaceInvariantMetric
open OPH.OrientedFaceBracketSelector

/-- The common port-dual scale as a real number. -/
noncomputable def measureScale : ℝ := 1 / 12

lemma measureScale_pos : (0 : ℝ) < measureScale := by
  norm_num [measureScale]

/-- At the sector-balanced measure point every classified family other
than `G` is strictly farther from the pinned face bracket. -/
theorem measure_point_unique_nearest_G
    (family : CompactFamily) (h : family ≠ .G) :
    dG2 measureScale measureScale measureScale <
      squaredDistanceAt measureScale measureScale measureScale family :=
  balanced_unique_nearest_G measureScale_pos measureScale_pos family h

/-- The bridge from the committed dual-measure theorem: under its
premises the port-dual weight at any port is the measure scale, and the
selection theorem holds at that weight. -/
theorem port_dual_weight_selects_G
    (vertexShare : ℚ) (faceWeight : FaceTriple → ℚ)
    (hbarycentric : vertexShare = 1 / 3)
    (hconstant : ∀ f ∈ baseFaces, ∀ g ∈ baseFaces,
      faceWeight f = faceWeight g)
    (hnormalized : ∑ f ∈ baseFaces, faceWeight f = 1)
    (p : Fin 12) (family : CompactFamily) (h : family ≠ .G) :
    dG2 ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ)
        ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ)
        ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ) <
      squaredDistanceAt
        ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ)
        ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ)
        ((portDualWeight vertexShare faceWeight p : ℚ) : ℝ) family := by
  have hw := equal_face_weights_imply_port_dual_one_twelfth
    vertexShare faceWeight hbarycentric hconstant hnormalized p
  rw [hw]
  have hcast : (((1 : ℚ) / 12 : ℚ) : ℝ) = measureScale := by
    norm_num [measureScale]
  rw [hcast]
  exact measure_point_unique_nearest_G family h

/-- Exact distance receipt: the squared `G` distance at the measure
point is twelve times the pinned reference distance. -/
theorem measure_point_dG2 :
    dG2 measureScale measureScale measureScale = 12 * distanceG := by
  have h : dG2 measureScale measureScale measureScale =
      12 * dG2 1 1 1 := by
    simp only [dG2, measureScale]
    ring
  rw [h, reference_G]

/-- Exact distance receipt for the mirror family `F`. -/
theorem measure_point_dF2 :
    dF2 measureScale measureScale measureScale = 12 * distanceF := by
  have h : dF2 measureScale measureScale measureScale =
      12 * dF2 1 1 1 := by
    simp only [dF2, measureScale]
    ring
  rw [h, reference_F]

/-- Exact distance receipt for the excluded plane `P`. -/
theorem measure_point_dP2 :
    dP2 measureScale measureScale measureScale = 12 * distanceP := by
  have h : dP2 measureScale measureScale measureScale =
      12 * dP2 1 1 1 := by
    simp only [dP2, measureScale]
    ring
  rw [h, reference_P]

#print axioms measure_point_unique_nearest_G
#print axioms port_dual_weight_selects_G
#print axioms measure_point_dG2
#print axioms measure_point_dF2
#print axioms measure_point_dP2

end OPH.PortDualMetricSelection
