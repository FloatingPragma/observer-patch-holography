import PositionSpaceMaxwellAction
import ObserverPatchHolography.CoreAxioms

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.LocalFaceMaxwellAction

open OPH.CoreAxioms
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.SeamCurrentCarrierQuotient

/-!
# A genuinely local face-curvature action on the committed carrier

This module uses the twenty oriented triangular faces already committed by
`CoreAxioms.orientedFaces`.  A real potential lives on the thirty committed
seams.  Its curvature on a face is the signed sum around that face, and the
kinetic action is the sum of the twenty squared face curvatures.  Thus the
construction is local by definition; the exact Hessian has five nonzero
entries per seam row (150 of 900 entries), unlike the dense Green/Hodge
projector in `PositionSpaceMaxwellAction`.

The exact dual-tree section proves that the face-incidence map has rank 19.
Its kernel is exactly the image of the port coboundary.  Consequently the
local sourced action is gauge invariant exactly for conserved seam sources,
is solvable exactly for those sources, and its minimizers are unique modulo
port gauge.

The scalar-potential Coulomb action is kept as a separately typed sector.  A
neutral port load has stationary potential `G rho` and field `d(G rho)`.
`StaticSources` then forms an explicit direct-product composition of a port
charge and a conserved seam current.  It supplies no map identifying those
two source types and no unproved continuity equation.

Everything here is exact finite mathematics on the committed tables.  No
physical position, time, current, detector, continuum, or photon attachment
is claimed.
-/

noncomputable section

/-! ## Exact oriented-face incidence -/

/-- Indexed form of the twenty committed oriented faces. -/
def faceVertices : Fin 20 → Fin 12 × Fin 12 × Fin 12 :=
  ![(0, 1, 2), (0, 3, 1), (0, 2, 4), (0, 6, 3), (0, 4, 6),
    (1, 5, 2), (1, 3, 7), (1, 7, 5), (2, 8, 4), (2, 5, 8),
    (3, 6, 9), (3, 9, 7), (4, 10, 6), (4, 8, 10), (5, 7, 11),
    (5, 11, 8), (6, 10, 9), (7, 9, 11), (8, 11, 10), (9, 10, 11)]

/-- The indexed table is exactly, and not merely isomorphic to, the A1
oriented-face packet. -/
theorem faceVertices_matches_committed :
    List.ofFn faceVertices = orientedFaces := by decide

/-- The three directed boundary edges of indexed face `f`. -/
def indexedFaceEdges (f : Fin 20) : List (Fin 12 × Fin 12) :=
  let t := faceVertices f
  [(t.1, t.2.1), (t.2.1, t.2.2), (t.2.2, t.1)]

/-- Signed face--seam incidence.  `+1` means that the committed smaller-to-
larger seam orientation agrees with the face orientation, and `-1` means it
is opposite. -/
def faceIncidenceZ (f : Fin 20) (e : Fin 30) : ℤ :=
  if (seamLeft e, seamRight e) ∈ indexedFaceEdges f then 1
  else if (seamRight e, seamLeft e) ∈ indexedFaceEdges f then -1
  else 0

/-- Every face contains exactly three committed seams. -/
theorem faceIncidence_three_per_face :
    ∀ f : Fin 20,
      (Finset.univ.filter (fun e : Fin 30 ↦ faceIncidenceZ f e ≠ 0)).card = 3 := by
  decide

/-- Every seam bounds exactly two faces, with opposite induced signs. -/
theorem faceIncidence_two_opposite_per_seam :
    ∀ e : Fin 30,
      (Finset.univ.filter (fun f : Fin 20 ↦ faceIncidenceZ f e ≠ 0)).card = 2 ∧
        (∑ f : Fin 20, faceIncidenceZ f e) = 0 := by
  decide

/-- Exact incidence support: twenty rows with three nonzeros each. -/
theorem faceIncidence_total_support :
    (∑ f : Fin 20,
      (Finset.univ.filter (fun e : Fin 30 ↦ faceIncidenceZ f e ≠ 0)).card) =
        60 := by
  simp [faceIncidence_three_per_face]

/-- Boundary-of-boundary is zero for the exact committed tables. -/
theorem face_port_incidence_product_zero :
    ∀ (f : Fin 20) (p : Fin 12),
      (∑ e : Fin 30, faceIncidenceZ f e * incidenceZ e p) = 0 := by
  set_option maxRecDepth 16384 in
    decide

/-- Orientation-erasure mutation control.  Replacing signed incidence by
its absolute value already makes boundary-of-boundary nonzero on face zero
at port zero. -/
theorem unsigned_face_incidence_breaks_boundary :
    (∑ e : Fin 30, |faceIncidenceZ 0 e| * incidenceZ e 0) = -2 := by
  decide

/-- Integer Hessian of the local face energy. -/
def localKineticZ (e d : Fin 30) : ℤ :=
  ∑ f : Fin 20, faceIncidenceZ f e * faceIncidenceZ f d

/-- Each seam couples only to itself and the four seams sharing one of its
two faces. -/
theorem localKineticZ_five_per_row :
    ∀ e : Fin 30,
      (Finset.univ.filter (fun d : Fin 30 ↦ localKineticZ e d ≠ 0)).card = 5 := by
  set_option maxRecDepth 16384 in
  set_option maxHeartbeats 2000000 in
    decide

/-- Exact sparsity receipt: 150 nonzero coefficients out of 900. -/
theorem localKineticZ_total_support :
    (∑ e : Fin 30,
      (Finset.univ.filter (fun d : Fin 30 ↦ localKineticZ e d ≠ 0)).card) =
        150 := by
  simp [localKineticZ_five_per_row]

/-- The local Hessian is symmetric. -/
theorem localKineticZ_symm (e d : Fin 30) :
    localKineticZ e d = localKineticZ d e := by
  unfold localKineticZ
  exact Finset.sum_congr rfl fun f _ ↦ mul_comm _ _

/-- Every seam lies on exactly two faces, so the local kinetic diagonal is
exactly two. -/
theorem localKineticZ_diagonal_two :
    ∀ e : Fin 30, localKineticZ e e = 2 := by
  decide

/-- Real scalar extension of face incidence. -/
def faceIncidenceR (f : Fin 20) (e : Fin 30) : ℝ :=
  (faceIncidenceZ f e : ℝ)

def faceIncidenceMatrixR : Matrix (Fin 20) (Fin 30) ℝ :=
  Matrix.of faceIncidenceR

/-- Local curvature: the oriented sum of the three seam potentials on each
committed face. -/
def faceCurvature : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 20 → ℝ) :=
  faceIncidenceMatrixR.mulVecLin

theorem faceCurvature_apply (A : Fin 30 → ℝ) (f : Fin 20) :
    faceCurvature A f = ∑ e : Fin 30, faceIncidenceR f e * A e := rfl

/-- Face pairing and energy. -/
def faceInner (F H : Fin 20 → ℝ) : ℝ := ∑ f : Fin 20, F f * H f
def faceEnergy (F : Fin 20 → ℝ) : ℝ := ∑ f : Fin 20, F f ^ 2

theorem faceInner_comm (F H : Fin 20 → ℝ) : faceInner F H = faceInner H F :=
  Finset.sum_congr rfl fun _ _ ↦ mul_comm _ _

theorem faceInner_self_eq_energy (F : Fin 20 → ℝ) :
    faceInner F F = faceEnergy F :=
  Finset.sum_congr rfl fun f _ ↦ (pow_two (F f)).symm

theorem faceEnergy_nonneg (F : Fin 20 → ℝ) : 0 ≤ faceEnergy F :=
  Finset.sum_nonneg fun f _ ↦ sq_nonneg (F f)

theorem faceEnergy_eq_zero_iff (F : Fin 20 → ℝ) :
    faceEnergy F = 0 ↔ F = 0 := by
  unfold faceEnergy
  rw [Finset.sum_eq_zero_iff_of_nonneg fun f _ ↦ sq_nonneg (F f)]
  constructor
  · intro h
    funext f
    exact pow_eq_zero_iff two_ne_zero |>.mp (h f (Finset.mem_univ f))
  · rintro rfl f _
    simp

theorem face_port_incidence_product_zero_real (f : Fin 20) (p : Fin 12) :
    (∑ e : Fin 30, faceIncidenceR f e * incidenceR e p) = 0 := by
  rw [show (0 : ℝ) = ((0 : ℤ) : ℝ) by norm_num,
    ← face_port_incidence_product_zero f p, Int.cast_sum]
  exact Finset.sum_congr rfl fun e _ ↦ by
    rw [Int.cast_mul, incidenceR_eq_cast]
    unfold faceIncidenceR
    rfl

/-- Gauge invariance of local curvature: the curl of a gradient vanishes. -/
theorem faceCurvature_coboundary (χ : Fin 12 → ℝ) :
    faceCurvature (realCoboundary χ) = 0 := by
  funext f
  rw [faceCurvature_apply]
  simp only [realCoboundary_apply, ← incidenceR_sum_apply]
  calc
    (∑ e : Fin 30,
        faceIncidenceR f e * ∑ p : Fin 12, incidenceR e p * χ p)
        = ∑ e : Fin 30, ∑ p : Fin 12,
            faceIncidenceR f e * (incidenceR e p * χ p) := by
              refine Finset.sum_congr rfl fun e _ ↦ ?_
              rw [Finset.mul_sum]
    _ = ∑ p : Fin 12, ∑ e : Fin 30,
          faceIncidenceR f e * (incidenceR e p * χ p) := Finset.sum_comm
    _
        = ∑ p : Fin 12,
            (∑ e : Fin 30, faceIncidenceR f e * incidenceR e p) * χ p := by
              refine Finset.sum_congr rfl fun p _ ↦ ?_
              rw [Finset.sum_mul]
              exact Finset.sum_congr rfl fun e _ ↦ by ring
    _ = 0 := by simp [face_port_incidence_product_zero_real]

/-! ## Exact rank and `ker curl = im gradient` -/

/-- Sum of all oriented face loads. -/
def faceTotal : (Fin 20 → ℝ) →ₗ[ℝ] ℝ where
  toFun ρ := ∑ f : Fin 20, ρ f
  map_add' ρ σ := Finset.sum_add_distrib
  map_smul' a ρ := by simp [Finset.mul_sum]

/-- Curvatures have zero total because each seam occurs on two oppositely
oriented faces. -/
theorem faceCurvature_total_zero (A : Fin 30 → ℝ) :
    faceTotal (faceCurvature A) = 0 := by
  change (∑ f : Fin 20, faceCurvature A f) = 0
  simp_rw [faceCurvature_apply]
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun e _ ↦ ?_
  rw [← Finset.sum_mul]
  have h := (faceIncidence_two_opposite_per_seam e).2
  rw [show (∑ f : Fin 20, faceIncidenceR f e) = 0 by
    rw [show (0 : ℝ) = ((0 : ℤ) : ℝ) by norm_num, ← h, Int.cast_sum]
    rfl]
  simp

set_option maxHeartbeats 2000000 in
/-- Explicit dual-tree section.  Only the nineteen displayed tree seams are
used; the formula was derived from the committed face--seam incidence and is
checked below, not assumed. -/
def faceBoundarySection (ρ : Fin 20 → ℝ) : Fin 30 → ℝ := ![
    -(ρ 1 + ρ 3 + ρ 6 + ρ 10 + ρ 11 + ρ 16 + ρ 17 + ρ 19),
    ρ 2 + ρ 4 + ρ 8 + ρ 12 + ρ 13 + ρ 18,
    -(ρ 3 + ρ 10 + ρ 16 + ρ 19), ρ 4 + ρ 12, 0,
    -(ρ 5 + ρ 7 + ρ 9 + ρ 14 + ρ 15),
    ρ 6 + ρ 11 + ρ 17, -(ρ 7 + ρ 14), 0,
    -(ρ 8 + ρ 13 + ρ 18), ρ 9 + ρ 15, 0,
    ρ 10 + ρ 16 + ρ 19, -(ρ 11 + ρ 17), 0, -ρ 12,
    ρ 13 + ρ 18, 0, ρ 14, -ρ 15, 0, -(ρ 16 + ρ 19), 0,
    ρ 17, 0, -ρ 18, 0, ρ 19, 0, 0]

set_option maxHeartbeats 2000000 in
/-- The dual-tree section reconstructs every zero-total face load. -/
theorem faceBoundarySection_is_section (ρ : Fin 20 → ℝ)
    (hρ : faceTotal ρ = 0) :
    faceCurvature (faceBoundarySection ρ) = ρ := by
  funext f
  fin_cases f <;>
    simp [faceTotal, faceCurvature_apply, faceIncidenceR, faceIncidenceZ,
      indexedFaceEdges, faceVertices, faceBoundarySection,
      seamLeft, seamRight, Fin.sum_univ_succ] at hρ ⊢ <;>
    linarith

/-- The range of local curvature is exactly the zero-total face subspace. -/
theorem range_faceCurvature :
    LinearMap.range faceCurvature = LinearMap.ker faceTotal := by
  ext ρ
  constructor
  · rintro ⟨A, rfl⟩
    exact LinearMap.mem_ker.mpr (faceCurvature_total_zero A)
  · intro hρ
    rw [LinearMap.mem_ker] at hρ
    exact ⟨faceBoundarySection ρ, faceBoundarySection_is_section ρ hρ⟩

theorem faceTotal_surjective : Function.Surjective faceTotal := by
  intro x
  let ρ : Fin 20 → ℝ := fun f ↦ if f = 0 then x else 0
  refine ⟨ρ, ?_⟩
  simp [faceTotal, ρ]

/-- The local curvature has exact rank nineteen. -/
theorem faceCurvature_range_finrank :
    Module.finrank ℝ (LinearMap.range faceCurvature) = 19 := by
  rw [range_faceCurvature]
  have h := LinearMap.finrank_range_add_finrank_ker faceTotal
  have hrange : LinearMap.range faceTotal = ⊤ :=
    LinearMap.range_eq_top.mpr faceTotal_surjective
  rw [hrange] at h
  simp at h ⊢
  omega

/-- The curvature kernel has exact rank eleven. -/
theorem faceCurvature_ker_finrank :
    Module.finrank ℝ (LinearMap.ker faceCurvature) = 11 := by
  have h := LinearMap.finrank_range_add_finrank_ker faceCurvature
  rw [faceCurvature_range_finrank] at h
  simp at h ⊢
  omega

/-- **Exact discrete Poincare lemma on the committed sphere:** a seam
potential has zero face curvature exactly when it is a port gradient. -/
theorem ker_faceCurvature_eq_gradient :
    LinearMap.ker faceCurvature = LinearMap.range realCoboundary := by
  apply (Submodule.eq_of_le_of_finrank_eq ?_ ?_).symm
  · intro A hA
    rcases hA with ⟨χ, rfl⟩
    exact LinearMap.mem_ker.mpr (faceCurvature_coboundary χ)
  · rw [faceCurvature_ker_finrank, gauge_finrank]

/-! ## Local Maxwell operator and action -/

/-- Transpose face incidence. -/
def faceCodifferential : (Fin 20 → ℝ) →ₗ[ℝ] (Fin 30 → ℝ) :=
  faceIncidenceMatrixRᵀ.mulVecLin

theorem faceCodifferential_apply (F : Fin 20 → ℝ) (e : Fin 30) :
    faceCodifferential F e = ∑ f : Fin 20, faceIncidenceR f e * F f := rfl

theorem faceCurvature_codifferential_adjoint (A : Fin 30 → ℝ)
    (F : Fin 20 → ℝ) :
    faceInner (faceCurvature A) F = realSeamInner A (faceCodifferential F) := by
  unfold faceInner realSeamInner
  simp only [faceCurvature_apply, faceCodifferential_apply, Finset.sum_mul,
    Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun e _ ↦
    Finset.sum_congr rfl fun f _ ↦ by ring

/-- Sparse local Maxwell operator `C^T C`. -/
def localMaxwellOperator : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 30 → ℝ) :=
  faceCodifferential ∘ₗ faceCurvature

theorem localMaxwellOperator_apply (A : Fin 30 → ℝ) (e : Fin 30) :
    localMaxwellOperator A e =
      ∑ d : Fin 30, (localKineticZ e d : ℝ) * A d := by
  rw [show localMaxwellOperator A e = faceCodifferential (faceCurvature A) e by rfl,
    faceCodifferential_apply]
  simp only [faceCurvature_apply]
  calc
    (∑ f : Fin 20,
      faceIncidenceR f e * ∑ d : Fin 30, faceIncidenceR f d * A d)
        = ∑ f : Fin 20, ∑ d : Fin 30,
            faceIncidenceR f e * (faceIncidenceR f d * A d) := by
              refine Finset.sum_congr rfl fun f _ ↦ ?_
              rw [Finset.mul_sum]
    _ = ∑ d : Fin 30, ∑ f : Fin 20,
          faceIncidenceR f e * (faceIncidenceR f d * A d) := Finset.sum_comm
    _
        = ∑ d : Fin 30,
            (∑ f : Fin 20, faceIncidenceR f e * faceIncidenceR f d) * A d := by
              refine Finset.sum_congr rfl fun d _ ↦ ?_
              rw [Finset.sum_mul]
              exact Finset.sum_congr rfl fun f _ ↦ by ring
    _ = ∑ d : Fin 30, (localKineticZ e d : ℝ) * A d := by
      refine Finset.sum_congr rfl fun d _ ↦ ?_
      congr 1
      rw [localKineticZ, Int.cast_sum]
      exact Finset.sum_congr rfl fun f _ ↦ by
        rw [Int.cast_mul]
        rfl

theorem realSeamInner_localMaxwellOperator (A B : Fin 30 → ℝ) :
    realSeamInner A (localMaxwellOperator B) =
      faceInner (faceCurvature A) (faceCurvature B) := by
  rw [faceCurvature_codifferential_adjoint]
  rfl

theorem ker_localMaxwellOperator :
    LinearMap.ker localMaxwellOperator = LinearMap.ker faceCurvature := by
  ext A
  rw [LinearMap.mem_ker, LinearMap.mem_ker]
  constructor
  · intro hK
    have hzero : faceEnergy (faceCurvature A) = 0 := by
      rw [← faceInner_self_eq_energy, ← realSeamInner_localMaxwellOperator, hK]
      simp [realSeamInner]
    exact (faceEnergy_eq_zero_iff _).mp hzero
  · intro hC
    change faceCodifferential (faceCurvature A) = 0
    rw [hC, map_zero]

theorem localMaxwellOperator_range_le_conserved :
    LinearMap.range localMaxwellOperator ≤ LinearMap.ker realBoundary := by
  rintro J ⟨A, rfl⟩
  rw [LinearMap.mem_ker]
  have hpair : realPortInner (realBoundary (localMaxwellOperator A))
      (realBoundary (localMaxwellOperator A)) = 0 := by
    rw [← realCoboundary_boundary_adjoint]
    rw [realSeamInner_localMaxwellOperator,
      faceCurvature_coboundary]
    simp [faceInner]
  exact (realPortInner_self_eq_zero_iff _).mp hpair

theorem localMaxwellOperator_range_finrank :
    Module.finrank ℝ (LinearMap.range localMaxwellOperator) = 19 := by
  have h := LinearMap.finrank_range_add_finrank_ker localMaxwellOperator
  rw [ker_localMaxwellOperator, faceCurvature_ker_finrank] at h
  simp at h ⊢
  omega

/-- The local Maxwell operator reaches exactly the conserved seam-current
subspace. -/
theorem range_localMaxwellOperator :
    LinearMap.range localMaxwellOperator = LinearMap.ker realBoundary := by
  apply Submodule.eq_of_le_of_finrank_eq localMaxwellOperator_range_le_conserved
  rw [localMaxwellOperator_range_finrank, cycleSpace_finrank]

/-- The genuinely local sourced face-curvature action. -/
def localSourcedAction (J A : Fin 30 → ℝ) : ℝ :=
  (1 / 2) * faceEnergy (faceCurvature A) - realSeamInner J A

/-- Gauge invariance holds exactly when the separately typed seam source is
conserved. -/
theorem localSourcedAction_gauge_invariant_iff (J : Fin 30 → ℝ) :
    (∀ (A : Fin 30 → ℝ) (χ : Fin 12 → ℝ),
      localSourcedAction J (A + realCoboundary χ) = localSourcedAction J A) ↔
      realBoundary J = 0 := by
  constructor
  · intro h
    have hp : ∀ χ : Fin 12 → ℝ,
        realSeamInner J (realCoboundary χ) = 0 := by
      intro χ
      have h0 := h 0 χ
      unfold localSourcedAction at h0
      rw [zero_add, faceCurvature_coboundary, map_zero] at h0
      simp [faceEnergy, realSeamInner] at h0
      simpa [realSeamInner] using h0
    have hself := hp (realBoundary J)
    rw [realSeamInner_comm, realCoboundary_boundary_adjoint] at hself
    exact (realPortInner_self_eq_zero_iff _).mp hself
  · intro hJ A χ
    unfold localSourcedAction
    have hC : faceCurvature (A + realCoboundary χ) = faceCurvature A := by
      rw [map_add, faceCurvature_coboundary, add_zero]
    rw [hC, realSeamInner_add_right]
    have hp : realSeamInner J (realCoboundary χ) = 0 := by
      rw [realSeamInner_comm, realCoboundary_boundary_adjoint, hJ]
      simp [realPortInner]
    rw [hp, add_zero]

theorem localSourcedAction_expansion (J A h : Fin 30 → ℝ) :
    localSourcedAction J (A + h) = localSourcedAction J A +
      realSeamInner (localMaxwellOperator A - J) h +
      (1 / 2) * faceEnergy (faceCurvature h) := by
  unfold localSourcedAction
  rw [map_add]
  have hsquare : faceEnergy (faceCurvature A + faceCurvature h) =
      faceEnergy (faceCurvature A) +
        2 * faceInner (faceCurvature A) (faceCurvature h) +
        faceEnergy (faceCurvature h) := by
    unfold faceEnergy faceInner
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun f _ ↦ by
      simp only [Pi.add_apply]
      ring
  have hcross : faceInner (faceCurvature A) (faceCurvature h) =
      realSeamInner (localMaxwellOperator A) h := by
    rw [faceInner_comm, ← realSeamInner_localMaxwellOperator]
    exact realSeamInner_comm _ _
  have hsource := realSeamInner_add_right J A h
  have hsub := realSeamInner_sub_left (localMaxwellOperator A) J h
  rw [hsquare, hcross, hsource, hsub]
  ring

theorem localStationary_iff (J A : Fin 30 → ℝ) :
    (∀ h : Fin 30 → ℝ,
      localSourcedAction J (A + h) = localSourcedAction J A +
        (1 / 2) * faceEnergy (faceCurvature h)) ↔
      localMaxwellOperator A = J := by
  constructor
  · intro hs
    have hlin : ∀ h : Fin 30 → ℝ,
        realSeamInner (localMaxwellOperator A - J) h = 0 := by
      intro h
      have h1 := localSourcedAction_expansion J A h
      have h2 := hs h
      linarith
    have hself := hlin (localMaxwellOperator A - J)
    rw [realSeamInner_self_eq_energy] at hself
    exact sub_eq_zero.mp ((realSeamEnergy_eq_zero_iff _).mp hself)
  · intro hA h
    rw [localSourcedAction_expansion, hA, sub_self]
    simp [realSeamInner]

theorem localStationary_solvable_iff (J : Fin 30 → ℝ) :
    (∃ A : Fin 30 → ℝ, localMaxwellOperator A = J) ↔ realBoundary J = 0 := by
  rw [← LinearMap.mem_range, range_localMaxwellOperator,
    LinearMap.mem_ker]

theorem localSolutions_differ_by_gauge (J A B : Fin 30 → ℝ)
    (hA : localMaxwellOperator A = J) :
    localMaxwellOperator B = J ↔
      B - A ∈ LinearMap.range realCoboundary := by
  rw [← hA]
  constructor
  · intro hB
    have hK : localMaxwellOperator (B - A) = 0 := by rw [map_sub, hB, sub_self]
    rw [← ker_faceCurvature_eq_gradient, ← ker_localMaxwellOperator]
    exact LinearMap.mem_ker.mpr hK
  · intro hBA
    have hK : localMaxwellOperator (B - A) = 0 := by
      have hCmem : B - A ∈ LinearMap.ker faceCurvature := by
        rw [ker_faceCurvature_eq_gradient]
        exact hBA
      have hKmem : B - A ∈ LinearMap.ker localMaxwellOperator := by
        rw [ker_localMaxwellOperator]
        exact hCmem
      exact LinearMap.mem_ker.mp hKmem
    have : B = A + (B - A) := by abel
    rw [this, map_add, hK, add_zero]

theorem localSourced_global_minimum (J A : Fin 30 → ℝ)
    (hA : localMaxwellOperator A = J) (B : Fin 30 → ℝ) :
    localSourcedAction J A ≤ localSourcedAction J B ∧
      (localSourcedAction J B = localSourcedAction J A ↔
        B - A ∈ LinearMap.range realCoboundary) := by
  have hab : A + (B - A) = B := by abel
  have hgap : localSourcedAction J B = localSourcedAction J A +
      (1 / 2) * faceEnergy (faceCurvature (B - A)) := by
    have h := localSourcedAction_expansion J A (B - A)
    rw [hab, hA, sub_self] at h
    simpa [realSeamInner] using h
  constructor
  · rw [hgap]
    exact le_add_of_nonneg_right
      (mul_nonneg (by norm_num) (faceEnergy_nonneg _))
  · constructor
    · intro heq
      have hz : faceEnergy (faceCurvature (B - A)) = 0 := by
        rw [hgap] at heq
        linarith
      have hC := (faceEnergy_eq_zero_iff _).mp hz
      rw [← ker_faceCurvature_eq_gradient]
      exact LinearMap.mem_ker.mpr hC
    · intro hg
      rw [← ker_faceCurvature_eq_gradient] at hg
      rw [hgap, LinearMap.mem_ker.mp hg]
      simp [faceEnergy]

/-! ## The separately typed scalar Coulomb action -/

/-- Scalar-potential action for a port charge. -/
def scalarCoulombAction (ρ φ : Fin 12 → ℝ) : ℝ :=
  (1 / 2) * realSeamEnergy (realCoboundary φ) - realPortInner ρ φ

theorem scalarCoulombAction_expansion (ρ φ h : Fin 12 → ℝ) :
    scalarCoulombAction ρ (φ + h) = scalarCoulombAction ρ φ +
      realPortInner (realLaplacian φ - ρ) h +
      (1 / 2) * realSeamEnergy (realCoboundary h) := by
  unfold scalarCoulombAction
  rw [map_add]
  have hsquare : realSeamEnergy (realCoboundary φ + realCoboundary h) =
      realSeamEnergy (realCoboundary φ) +
        2 * realSeamInner (realCoboundary φ) (realCoboundary h) +
        realSeamEnergy (realCoboundary h) := by
    unfold realSeamEnergy realSeamInner
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun e _ ↦ by
      simp only [Pi.add_apply]
      ring
  have hcross : realSeamInner (realCoboundary φ) (realCoboundary h) =
      realPortInner (realLaplacian φ) h := by
    rw [realSeamInner_comm, realCoboundary_boundary_adjoint,
      ← realLaplacian_apply, realPortInner_comm]
  have hsource : realPortInner ρ (φ + h) =
      realPortInner ρ φ + realPortInner ρ h := by
    unfold realPortInner
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun p _ ↦ by simp [mul_add]
  have hsub : realPortInner (realLaplacian φ - ρ) h =
      realPortInner (realLaplacian φ) h - realPortInner ρ h := by
    unfold realPortInner
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun p _ ↦ by simp [sub_mul]
  rw [hsquare, hcross, hsource, hsub]
  ring

theorem scalarStationary_iff (ρ φ : Fin 12 → ℝ) :
    (∀ h : Fin 12 → ℝ,
      scalarCoulombAction ρ (φ + h) = scalarCoulombAction ρ φ +
        (1 / 2) * realSeamEnergy (realCoboundary h)) ↔
      realLaplacian φ = ρ := by
  constructor
  · intro hs
    have hlin : ∀ h : Fin 12 → ℝ,
        realPortInner (realLaplacian φ - ρ) h = 0 := by
      intro h
      have h1 := scalarCoulombAction_expansion ρ φ h
      have h2 := hs h
      linarith
    have hself := hlin (realLaplacian φ - ρ)
    exact sub_eq_zero.mp ((realPortInner_self_eq_zero_iff _).mp hself)
  · intro hφ h
    rw [scalarCoulombAction_expansion, hφ, sub_self]
    simp [realPortInner]

/-- Constant-shift invariance of the scalar action is equivalent to charge
neutrality. -/
theorem scalarCoulombAction_constant_gauge_iff (ρ : Fin 12 → ℝ) :
    (∀ (φ : Fin 12 → ℝ) (c : ℝ),
      scalarCoulombAction ρ (φ + fun _ ↦ c) = scalarCoulombAction ρ φ) ↔
      (∑ p : Fin 12, ρ p) = 0 := by
  constructor
  · intro h
    have h0 := h 0 1
    unfold scalarCoulombAction at h0
    have hc : realCoboundary (fun _ : Fin 12 ↦ (1 : ℝ)) = 0 :=
      (realCoboundary_eq_zero_iff _).mpr ⟨1, rfl⟩
    rw [zero_add, hc, map_zero] at h0
    simp [realSeamEnergy, realPortInner] at h0
    exact h0
  · intro hρ φ c
    unfold scalarCoulombAction
    have hc : realCoboundary (fun _ : Fin 12 ↦ c) = 0 :=
      (realCoboundary_eq_zero_iff _).mpr ⟨c, rfl⟩
    rw [map_add, hc, add_zero]
    unfold realPortInner
    simp only [Pi.add_apply, mul_add, Finset.sum_add_distrib]
    have hconst : (∑ p : Fin 12, ρ p * c) = 0 := by
      calc
        (∑ p : Fin 12, ρ p * c) = (∑ p : Fin 12, ρ p) * c :=
          (Finset.sum_mul Finset.univ (fun p ↦ ρ p) c).symm
        _ = 0 := by rw [hρ, zero_mul]
    rw [hconst, add_zero]

/-- Canonical Green potential; unlike `realCoulombField`, this is a port
scalar rather than a seam field. -/
def greenPotential (ρ : Fin 12 → ℝ) : Fin 12 → ℝ := greenMatrixR.mulVec ρ

theorem greenPotential_stationary (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) :
    realLaplacian (greenPotential ρ) = ρ := by
  change realBoundary (realCoulombField ρ) = ρ
  exact realCoulombField_gauss ρ hρ

theorem greenPotential_field (ρ : Fin 12 → ℝ) :
    realCoboundary (greenPotential ρ) = realCoulombField ρ := rfl

theorem scalarCoulomb_global_minimum (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) (ψ : Fin 12 → ℝ) :
    scalarCoulombAction ρ (greenPotential ρ) ≤ scalarCoulombAction ρ ψ ∧
      (scalarCoulombAction ρ ψ = scalarCoulombAction ρ (greenPotential ρ) ↔
        ∃ c : ℝ, ψ - greenPotential ρ = fun _ ↦ c) := by
  have hsol := greenPotential_stationary ρ hρ
  have hadd : greenPotential ρ + (ψ - greenPotential ρ) = ψ := by abel
  have hgap : scalarCoulombAction ρ ψ = scalarCoulombAction ρ (greenPotential ρ) +
      (1 / 2) * realSeamEnergy
        (realCoboundary (ψ - greenPotential ρ)) := by
    have h := scalarCoulombAction_expansion ρ (greenPotential ρ)
      (ψ - greenPotential ρ)
    rw [hadd, hsol, sub_self] at h
    simpa [realPortInner] using h
  constructor
  · rw [hgap]
    exact le_add_of_nonneg_right
      (mul_nonneg (by norm_num) (realSeamEnergy_nonneg _))
  · constructor
    · intro heq
      have hz : realSeamEnergy (realCoboundary (ψ - greenPotential ρ)) = 0 := by
        rw [hgap] at heq
        linarith
      exact (realCoboundary_eq_zero_iff _).mp
        ((realSeamEnergy_eq_zero_iff _).mp hz)
    · rintro ⟨c, hc⟩
      rw [hgap, hc]
      have hz : realCoboundary (fun _ : Fin 12 ↦ c) = 0 :=
        (realCoboundary_eq_zero_iff _).mpr ⟨c, rfl⟩
      rw [hz]
      simp [realSeamEnergy]

/-! ## Honest direct-product static composition -/

/-- The two source types are deliberately distinct fields.  This structure
asserts neutrality and conservation separately and contains no conversion
map between a port charge and a seam current. -/
structure StaticSources where
  charge : Fin 12 → ℝ
  seamCurrent : Fin 30 → ℝ
  charge_neutral : (∑ p : Fin 12, charge p) = 0
  current_conserved : realBoundary seamCurrent = 0

structure StaticPotentials where
  scalar : Fin 12 → ℝ
  vector : Fin 30 → ℝ

/-- Direct sum of the scalar Coulomb action and local face-Maxwell action. -/
def staticAction (S : StaticSources) (X : StaticPotentials) : ℝ :=
  scalarCoulombAction S.charge X.scalar +
    localSourcedAction S.seamCurrent X.vector

/-- Independent scalar-constant and vector-gradient gauge transformation. -/
def staticGaugeTransform (X : StaticPotentials) (c : ℝ)
    (χ : Fin 12 → ℝ) : StaticPotentials where
  scalar := X.scalar + fun _ ↦ c
  vector := X.vector + realCoboundary χ

/-- Gauge invariance of the typed direct product uses charge neutrality and
seam-current conservation separately. -/
theorem staticAction_gauge_invariant (S : StaticSources) (X : StaticPotentials)
    (c : ℝ) (χ : Fin 12 → ℝ) :
    staticAction S (staticGaugeTransform X c χ) = staticAction S X := by
  unfold staticAction staticGaugeTransform
  rw [(scalarCoulombAction_constant_gauge_iff S.charge).2 S.charge_neutral]
  rw [(localSourcedAction_gauge_invariant_iff S.seamCurrent).2
    S.current_conserved]

/-- A canonical scalar representative and an arbitrary chosen local vector
solution.  The latter exists exactly because the seam source is conserved. -/
def canonicalStaticPotentials (S : StaticSources) : StaticPotentials where
  scalar := greenPotential S.charge
  vector := Classical.choose ((localStationary_solvable_iff S.seamCurrent).2
    S.current_conserved)

theorem canonicalStaticPotentials_vector_stationary (S : StaticSources) :
    localMaxwellOperator (canonicalStaticPotentials S).vector = S.seamCurrent :=
  Classical.choose_spec ((localStationary_solvable_iff S.seamCurrent).2
    S.current_conserved)

/-- The direct-product action has a global minimizer, and equality is
exactly independent scalar-constant gauge plus vector-gradient gauge. -/
theorem staticAction_global_minimum (S : StaticSources) (Y : StaticPotentials) :
    staticAction S (canonicalStaticPotentials S) ≤ staticAction S Y ∧
      (staticAction S Y = staticAction S (canonicalStaticPotentials S) ↔
        (∃ c : ℝ,
          Y.scalar - (canonicalStaticPotentials S).scalar = fun _ ↦ c) ∧
        Y.vector - (canonicalStaticPotentials S).vector ∈
          LinearMap.range realCoboundary) := by
  have hs := scalarCoulomb_global_minimum S.charge S.charge_neutral Y.scalar
  have hv := localSourced_global_minimum S.seamCurrent
    (canonicalStaticPotentials S).vector
    (canonicalStaticPotentials_vector_stationary S) Y.vector
  constructor
  · unfold staticAction
    have hs' : scalarCoulombAction S.charge
        (canonicalStaticPotentials S).scalar ≤
        scalarCoulombAction S.charge Y.scalar := by
      simpa [canonicalStaticPotentials] using hs.1
    linarith [hs', hv.1]
  · constructor
    · intro heq
      have hse : scalarCoulombAction S.charge Y.scalar =
          scalarCoulombAction S.charge (canonicalStaticPotentials S).scalar := by
        unfold staticAction at heq
        have hsn : scalarCoulombAction S.charge
            (canonicalStaticPotentials S).scalar ≤
            scalarCoulombAction S.charge Y.scalar := by
          simpa [canonicalStaticPotentials] using hs.1
        have hvn := hv.1
        linarith
      have hve : localSourcedAction S.seamCurrent Y.vector =
          localSourcedAction S.seamCurrent
            (canonicalStaticPotentials S).vector := by
        unfold staticAction at heq
        have hsn : scalarCoulombAction S.charge
            (canonicalStaticPotentials S).scalar ≤
            scalarCoulombAction S.charge Y.scalar := by
          simpa [canonicalStaticPotentials] using hs.1
        have hvn := hv.1
        linarith
      have hscalar : ∃ c : ℝ,
          Y.scalar - (canonicalStaticPotentials S).scalar = fun _ ↦ c := by
        apply hs.2.mp
        simpa [canonicalStaticPotentials] using hse
      exact ⟨hscalar, hv.2.mp hve⟩
    · rintro ⟨hscalar, hvector⟩
      have hse := hs.2.mpr hscalar
      have hve := hv.2.mpr hvector
      unfold staticAction
      simpa [canonicalStaticPotentials] using congrArg₂ (· + ·) hse hve

end

end OPH.LocalFaceMaxwellAction

/- Axiom audit: exact finite tables and ordinary finite-dimensional real
linear algebra only. -/
#print axioms OPH.LocalFaceMaxwellAction.faceVertices_matches_committed
#print axioms OPH.LocalFaceMaxwellAction.face_port_incidence_product_zero
#print axioms OPH.LocalFaceMaxwellAction.unsigned_face_incidence_breaks_boundary
#print axioms OPH.LocalFaceMaxwellAction.localKineticZ_five_per_row
#print axioms OPH.LocalFaceMaxwellAction.localKineticZ_total_support
#print axioms OPH.LocalFaceMaxwellAction.localKineticZ_diagonal_two
#print axioms OPH.LocalFaceMaxwellAction.faceBoundarySection_is_section
#print axioms OPH.LocalFaceMaxwellAction.ker_faceCurvature_eq_gradient
#print axioms OPH.LocalFaceMaxwellAction.range_localMaxwellOperator
#print axioms OPH.LocalFaceMaxwellAction.localSourcedAction_gauge_invariant_iff
#print axioms OPH.LocalFaceMaxwellAction.localStationary_solvable_iff
#print axioms OPH.LocalFaceMaxwellAction.localSourced_global_minimum
#print axioms OPH.LocalFaceMaxwellAction.greenPotential_stationary
#print axioms OPH.LocalFaceMaxwellAction.scalarCoulombAction_constant_gauge_iff
#print axioms OPH.LocalFaceMaxwellAction.scalarCoulomb_global_minimum
#print axioms OPH.LocalFaceMaxwellAction.staticAction_gauge_invariant
#print axioms OPH.LocalFaceMaxwellAction.staticAction_global_minimum
