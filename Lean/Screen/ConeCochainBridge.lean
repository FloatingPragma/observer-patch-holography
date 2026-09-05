import LocalFaceMaxwellAction

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.ConeCochainBridge

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.SeamCurrentCarrierQuotient

noncomputable section

/-!
# A gauge-covariant cochain extension to the cone of the committed surface

The cone has one apex, twelve boundary vertices, thirty boundary edges,
twelve radial edges, twenty boundary triangles, thirty radial triangles,
and twenty tetrahedra. Radial edges are oriented apex-to-vertex; the radial
triangle above an oriented edge `(l,r)` is `(apex,l,r)`. These are explicit
finite cochain objects, not a geometric embedding or a physical volume.

The committed graph Green operator constructs the radial edge potential
`g(A)=G boundary(A)`. Scalar data extend with their mean at the apex, so
the extension commutes with the gradient. Curvature has boundary value `CA`
and radial value `A-dg(A)`, the committed cycle projector. A second extension
constructed with the committed dual-tree section makes the curl square
commute as well. Closed extension of a boundary two-form exists exactly
when its total flux vanishes; no arbitrary-flux extension is claimed.

No metric, Hodge operator, action equality, dynamics, refinement, source
selection, physical time, or continuum approximation is supplied here.
-/

abbrev ConeZero := ℝ × (Fin 12 → ℝ)
abbrev ConeOne := (Fin 12 → ℝ) × (Fin 30 → ℝ)
abbrev ConeTwo := (Fin 30 → ℝ) × (Fin 20 → ℝ)

/-- The first component is the apex value; the second is the boundary trace. -/
def coneGradient (u : ConeZero) : ConeOne :=
  (fun p ↦ u.2 p - u.1, realCoboundary u.2)

/-- The radial triangle value is `A(l,r)-g(r)+g(l)`. -/
def coneCurl (a : ConeOne) : ConeTwo :=
  (a.2 - realCoboundary a.1, faceCurvature a.2)

/-- Oriented tetrahedral boundary: outer face minus the three radial faces. -/
def coneDivergence (b : ConeTwo) : Fin 20 → ℝ :=
  b.2 - faceCurvature b.1

/-- Incidence commutation in degrees zero and one on the actual carrier. -/
theorem coneCurl_gradient (u : ConeZero) :
    coneCurl (coneGradient u) = 0 := by
  apply Prod.ext
  · funext e
    simp [coneCurl, coneGradient, realCoboundary_apply]
  · exact faceCurvature_coboundary u.2

/-- Incidence commutation in degrees one and two, at every cone tetrahedron. -/
theorem coneDivergence_curl (a : ConeOne) :
    coneDivergence (coneCurl a) = 0 := by
  simp [coneDivergence, coneCurl, map_sub, faceCurvature_coboundary]

/-- The scalar mean used at the apex, with no selected boundary vertex. -/
def boundaryMean (phi : Fin 12 → ℝ) : ℝ := (∑ p, phi p) / 12

/-- The radial potential is computed by the committed exact Green matrix. -/
def radialPotential : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 12 → ℝ) :=
  greenMatrixR.mulVecLin ∘ₗ realBoundary

theorem radialPotential_apply (A : Fin 30 → ℝ) :
    radialPotential A = greenMatrixR.mulVec (realBoundary A) := rfl

/-- The Green identity fixes the precise constant-gauge correction. -/
theorem radialPotential_gradient (phi : Fin 12 → ℝ) :
    radialPotential (realCoboundary phi) = fun p ↦ phi p - boundaryMean phi := by
  rw [radialPotential_apply, ← realLaplacian_apply, realLaplacian_eq_mulVec,
    Matrix.mulVec_mulVec, green_mul_laplacian_real, Matrix.sub_mulVec,
    Matrix.one_mulVec, Matrix.smul_mulVec, allOnesR_mulVec]
  funext p
  simp [boundaryMean]
  ring

/-- The radial gauge potential has zero mean. -/
theorem radialPotential_sum_zero (A : Fin 30 → ℝ) :
    (∑ p, radialPotential A p) = 0 := by
  simp only [radialPotential_apply, Matrix.mulVec, dotProduct]
  rw [Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro q _
  rw [← Finset.sum_mul]
  have hcol : (∑ p, greenMatrixR p q) = 0 := by
    have hs : ∀ p, greenMatrixR p q = greenMatrixR q p := by
      intro p
      exact congrFun (congrFun greenMatrixR_symm q) p
    simp_rw [hs]
    exact greenMatrixR_row_sum q
  rw [hcol, zero_mul]

def extendZero (phi : Fin 12 → ℝ) : ConeZero := (boundaryMean phi, phi)

def extendOne (A : Fin 30 → ℝ) : ConeOne := (radialPotential A, A)

/-- The extension preserves the boundary scalar and edge cochains literally. -/
theorem extension_trace (phi : Fin 12 → ℝ) (A : Fin 30 → ℝ) :
    (extendZero phi).2 = phi ∧ (extendOne A).2 = A := ⟨rfl, rfl⟩

/-- The scalar-to-edge square commutes, including the apex gauge value. -/
theorem extendOne_gradient (phi : Fin 12 → ℝ) :
    extendOne (realCoboundary phi) = coneGradient (extendZero phi) := by
  apply Prod.ext
  · exact radialPotential_gradient phi
  · rfl

/-- Boundary gauge changes lift to ordinary cone gauge changes. -/
theorem extendOne_gauge_covariant (A : Fin 30 → ℝ) (phi : Fin 12 → ℝ) :
    extendOne (A + realCoboundary phi) =
      extendOne A + coneGradient (extendZero phi) := by
  rw [← extendOne_gradient]
  exact Prod.ext (map_add radialPotential A (realCoboundary phi)) rfl

/-- The radial curvature is exactly the committed cycle component. -/
theorem coneCurl_extendOne (A : Fin 30 → ℝ) :
    coneCurl (extendOne A) = (fieldProjector A, faceCurvature A) := by
  apply Prod.ext
  · exact (fieldProjector_apply A).symm
  · rfl

theorem coneCurl_gauge_invariant (A : Fin 30 → ℝ) (phi : Fin 12 → ℝ) :
    coneCurl (extendOne (A + realCoboundary phi)) =
      coneCurl (extendOne A) := by
  rw [coneCurl_extendOne, coneCurl_extendOne]
  apply Prod.ext
  · change fieldProjector (A + realCoboundary phi) = fieldProjector A
    rw [map_add, fieldProjector_coboundary, add_zero]
  · change faceCurvature (A + realCoboundary phi) = faceCurvature A
    rw [map_add, faceCurvature_coboundary, add_zero]

/-- Outer curvature is preserved exactly, with no reindexing of the faces. -/
theorem coneCurl_boundary_trace (A : Fin 30 → ℝ) :
    (coneCurl (extendOne A)).2 = faceCurvature A := rfl

/-- Every extended potential has exactly zero tetrahedral divergence. -/
theorem extended_curvature_closed (A : Fin 30 → ℝ) :
    coneDivergence (coneCurl (extendOne A)) = 0 :=
  coneDivergence_curl (extendOne A)

/-- Removing a gradient does not alter face curvature. -/
theorem faceCurvature_projector (A : Fin 30 → ℝ) :
    faceCurvature (fieldProjector A) = faceCurvature A := by
  rw [fieldProjector_apply, map_sub, faceCurvature_coboundary, sub_zero]

/-- Equal curvatures have equal radial cycle readouts. -/
theorem projector_eq_of_curvature_eq (A B : Fin 30 → ℝ)
    (h : faceCurvature A = faceCurvature B) :
    fieldProjector A = fieldProjector B := by
  apply (fieldProjector_eq_iff A B).mpr
  rw [← ker_faceCurvature_eq_gradient, LinearMap.mem_ker, map_sub, h, sub_self]

/-- The radial component of a boundary two-form is obtained from the exact
dual-tree inverse, followed by the graph cycle projector. Closedness is
proved only for zero-total boundary forms. -/
def extendTwo (B : Fin 20 → ℝ) : ConeTwo :=
  (fieldProjector (faceBoundarySection B), B)

theorem extendTwo_trace (B : Fin 20 → ℝ) : (extendTwo B).2 = B := rfl

/-- The edge-to-face square commutes on every boundary potential. -/
theorem extendTwo_curvature (A : Fin 30 → ℝ) :
    extendTwo (faceCurvature A) = coneCurl (extendOne A) := by
  rw [coneCurl_extendOne]
  apply Prod.ext
  · exact projector_eq_of_curvature_eq _ A
      (faceBoundarySection_is_section _ (faceCurvature_total_zero A))
  · rfl

/-- A zero-flux boundary two-form extends to a closed cone two-form. -/
theorem extendTwo_closed (B : Fin 20 → ℝ) (hB : faceTotal B = 0) :
    coneDivergence (extendTwo B) = 0 := by
  unfold coneDivergence extendTwo
  rw [faceCurvature_projector, faceBoundarySection_is_section B hB, sub_self]

/-- On zero-total boundary forms the extension is additive. No linearity
claim on an unrestricted, nonzero-flux extension is needed. -/
theorem extendTwo_add (B F : Fin 20 → ℝ)
    (hB : faceTotal B = 0) (hF : faceTotal F = 0) :
    extendTwo (B + F) = extendTwo B + extendTwo F := by
  apply Prod.ext
  · change fieldProjector (faceBoundarySection (B + F)) =
      fieldProjector (faceBoundarySection B) + fieldProjector (faceBoundarySection F)
    rw [← map_add]
    apply projector_eq_of_curvature_eq
    rw [faceBoundarySection_is_section _ (by rw [map_add, hB, hF, add_zero]),
      map_add, faceBoundarySection_is_section B hB, faceBoundarySection_is_section F hF]
  · rfl

/-- The same zero-flux extension respects scalar multiplication. -/
theorem extendTwo_smul (c : ℝ) (B : Fin 20 → ℝ) (hB : faceTotal B = 0) :
    extendTwo (c • B) = c • extendTwo B := by
  apply Prod.ext
  · change fieldProjector (faceBoundarySection (c • B)) =
      c • fieldProjector (faceBoundarySection B)
    rw [← map_smul]
    apply projector_eq_of_curvature_eq
    rw [faceBoundarySection_is_section _ (by rw [map_smul, hB, smul_zero]),
      map_smul, faceBoundarySection_is_section B hB]
  · rfl

/-- Every co-closed radial solution with the right curl agrees with this
extension. This allows comparison to other exact inverses, including a
Moore--Penrose construction, on the zero-flux domain. -/
theorem extendTwo_radial_unique (B : Fin 20 → ℝ) (hB : faceTotal B = 0)
    (r : Fin 30 → ℝ) (hr : realBoundary r = 0) (hc : faceCurvature r = B) :
    (extendTwo B).1 = r := by
  change fieldProjector (faceBoundarySection B) = r
  rw [← fieldProjector_fixes r hr]
  exact projector_eq_of_curvature_eq _ _
    ((faceBoundarySection_is_section B hB).trans hc.symm)

/-- The exact obstruction and constructive converse: closed cone extension
exists precisely for zero total boundary flux. This quantifies over all
radial two-cochains, not only the constructed extension. -/
theorem closed_extension_iff_zero_flux (B : Fin 20 → ℝ) :
    (∃ b : ConeTwo, b.2 = B ∧ coneDivergence b = 0) ↔ faceTotal B = 0 := by
  constructor
  · rintro ⟨b, hb, hclosed⟩
    have h : b.2 = faceCurvature b.1 := sub_eq_zero.mp hclosed
    rw [← hb, h]
    exact faceCurvature_total_zero b.1
  · intro hB
    exact ⟨extendTwo B, rfl, extendTwo_closed B hB⟩

/-- A nonzero constant flux cannot extend to a closed cone two-form. -/
theorem constant_unit_flux_has_no_closed_extension :
    ¬ ∃ b : ConeTwo, b.2 = (fun _ ↦ 1) ∧ coneDivergence b = 0 := by
  rw [closed_extension_iff_zero_flux]
  norm_num [faceTotal]

/-- Fixing the apex to zero breaks covariance even for a constant boundary
gauge. The apex mean in `extendZero` is therefore load-bearing. -/
theorem zero_apex_constant_gauge_fails :
    extendOne (realCoboundary (fun _ ↦ 1)) ≠
      coneGradient ((0 : ℝ), fun _ ↦ 1) := by
  intro h
  have hd : realCoboundary (fun _ ↦ (1 : ℝ)) = 0 :=
    (realCoboundary_eq_zero_iff _).mpr ⟨1, rfl⟩
  rw [hd] at h
  have hh := congrFun (congrArg Prod.fst h) (0 : Fin 12)
  simp [extendOne, coneGradient] at hh

end

#print axioms radialPotential_gradient
#print axioms radialPotential_sum_zero
#print axioms extendOne_gradient
#print axioms extendOne_gauge_covariant
#print axioms coneCurl_extendOne
#print axioms coneCurl_gauge_invariant
#print axioms extended_curvature_closed
#print axioms extendTwo_curvature
#print axioms extendTwo_closed
#print axioms extendTwo_add
#print axioms extendTwo_smul
#print axioms extendTwo_radial_unique
#print axioms closed_extension_iff_zero_flux
#print axioms constant_unit_flux_has_no_closed_extension
#print axioms zero_apex_constant_gauge_fails

end OPH.ConeCochainBridge
