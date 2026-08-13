import Mathlib.LinearAlgebra.CrossProduct
import SeamCurrentFreePhotonLift

open scoped Matrix

namespace OPH.ModalMaxwellFactorizationBoundary

open OPH.PrimitivePortTranslationBridge
open OPH.SeamCurrentFreePhotonLift

/-!
# Exact modal cross-product factorization and its Maxwell boundary

The committed light-signal surface supplies a nonnegative scalar FZ-12
symbol on a rank-two transverse fibre and a conditional second-order
oscillator.  This file proves that, at every nonzero chart momentum, that
scalar transverse action has an explicit real cross-product square root.
After complexification, the standard Fourier curl multiplier `i k ×`
and opposite electric/magnetic signs give a first-order Maxwell-shaped system
whose second iterate is exactly the existing FZ-12 wave operator and whose
curl contribution has zero modal divergence.

The result is deliberately modal.  The multiplier
`photonModeFrequency a k / sqrt (dot k k)` depends on momentum, so no local
position-space curl operator is claimed.  The final negative control replaces
the required opposite electric/magnetic signs by equal signs and proves that
this mutation has the wrong second-order sign.

No theorem identifies either amplitude with a physical electric or magnetic
field, supplies the reality pairing between opposite momenta, connects this
complex modal divergence to the rational seam-incidence Gauss receipt,
constructs a U(1) potential, gauge quotient, or action, supplies a conserved
physical current, assembles the modes into a field, proves locality or
continuum control, or establishes Lorentz covariance or laboratory readout.
-/

/-- The real three-coordinate modal carrier already used by the free-vector
lift. -/
abbrev ModalVec3 := FreePhotonVec3

/-- A pair of modal amplitudes.  The names electric/magnetic are intentionally
not built into the type: physical identification remains open. -/
abbrev PairedModalState := ModalVec3 × ModalVec3

/-- The exact scalar needed to turn the coordinate cross product into a
square root of the FZ-12 transverse symbol at nonzero momentum. -/
noncomputable def modalCurlScale (a : ℝ) (k : ModalVec3) : ℝ :=
  photonModeFrequency a k / Real.sqrt (dot k k)

/-- Momentum-space curl candidate.  At `k = 0` Lean's totalized division
makes this zero; the square-root theorem is stated only for `k ≠ 0`. -/
noncomputable def modalCurl (a : ℝ) (k v : ModalVec3) : ModalVec3 :=
  modalCurlScale a k • crossProduct k v

/-- The normalized curl scale has exactly the coefficient needed to square
to the committed frequency. -/
theorem modalCurlScale_sq_mul_dot_self {a : ℝ} {k : ModalVec3}
    (hk : k ≠ 0) :
    modalCurlScale a k ^ 2 * dot k k = photonModeFrequency a k ^ 2 := by
  have hkk : 0 < dot k k := dot_self_pos hk
  have hsqrt : 0 < Real.sqrt (dot k k) := Real.sqrt_pos.2 hkk
  have hsqrt_sq : Real.sqrt (dot k k) ^ 2 = dot k k :=
    Real.sq_sqrt hkk.le
  unfold modalCurlScale
  rw [div_pow]
  field_simp [ne_of_gt hsqrt]
  nlinarith

/-- Modal `div curl = 0`: the cross-product square root always lands in the
transverse plane. -/
theorem dot_modalCurl_zero (a : ℝ) (k v : ModalVec3) :
    dot k (modalCurl a k v) = 0 := by
  unfold modalCurl
  rw [dot_smul_right]
  have hcross : dot k (crossProduct k v) = 0 := by
    simpa [dot, dotProduct] using dot_self_cross k v
  rw [hcross, mul_zero]

/-- The modal curl candidate maps every vector into the committed transverse
fibre. -/
theorem modalCurl_mem_transverse (a : ℝ) (k v : ModalVec3) :
    modalCurl a k v ∈ TransversePolarization k := by
  change dot k (modalCurl a k v) = 0
  exact dot_modalCurl_zero a k v

/-- The exact vector triple-product identity in the `dot` spelling used by
the OPH modal modules. -/
theorem cross_cross_eq_dot_smul_sub_dot_smul
    (k v : ModalVec3) :
    crossProduct k (crossProduct k v) =
      dot k v • k - dot k k • v := by
  simpa [dot, dotProduct] using cross_cross_eq_smul_sub_smul' k k v

/-- **Exact curl square root.**  On the transverse fibre at nonzero momentum,
applying the modal curl twice is minus the committed FZ-12 spatial action. -/
theorem modalCurl_sq_on_transverse {a : ℝ} {k v : ModalVec3}
    (hk : k ≠ 0) (hv : v ∈ TransversePolarization k) :
    modalCurl a k (modalCurl a k v) = -photonSpatialAction a k v := by
  have hvdot : dot k v = 0 := hv
  have hscale := modalCurlScale_sq_mul_dot_self (a := a) (k := k) hk
  rw [photonSpatialAction_eq_frequency_sq]
  unfold modalCurl
  simp only [map_smul, smul_smul]
  rw [cross_cross_eq_dot_smul_sub_dot_smul, hvdot]
  ext d
  simp only [zero_smul, zero_sub, Pi.smul_apply, Pi.neg_apply,
    smul_eq_mul]
  rw [← hscale]
  ring

/-- The modal curl is linear in its amplitude argument. -/
theorem modalCurl_smul (a : ℝ) (k : ModalVec3) (c : ℝ) (v : ModalVec3) :
    modalCurl a k (c • v) = c • modalCurl a k v := by
  unfold modalCurl
  simp only [map_smul, smul_smul]
  rw [mul_comm]

/-- Cross-product skewness in the OPH `dot` spelling. -/
theorem dot_modalCurl_skew (a : ℝ) (k u v : ModalVec3) :
    dot u (modalCurl a k v) = -dot (modalCurl a k u) v := by
  have htriple : dot u (crossProduct k v) =
      -dot (crossProduct k u) v := by
    change u ⬝ᵥ crossProduct k v =
      -(crossProduct k u ⬝ᵥ v)
    calc
      u ⬝ᵥ crossProduct k v = k ⬝ᵥ crossProduct v u :=
        triple_product_permutation u k v
      _ = k ⬝ᵥ -(crossProduct u v) := by
        rw [cross_anticomm]
      _ = -(k ⬝ᵥ crossProduct u v) := dotProduct_neg _ _
      _ = -(v ⬝ᵥ crossProduct k u) := by
        rw [triple_product_permutation v k u]
      _ = -(crossProduct k u ⬝ᵥ v) := by
        rw [dotProduct_comm]
  unfold modalCurl
  rw [dot_smul_right, htriple]
  rw [dot_comm (modalCurlScale a k • crossProduct k u) v,
    dot_smul_right, dot_comm v (crossProduct k u)]
  ring

/-! ## The exact complex Fourier-Maxwell-shaped pairing -/

/-- Complex modal amplitudes.  A physical real field would additionally
require a proved conjugate relation between the `k` and `-k` fibres. -/
abbrev ComplexModalVec3 := Fin 3 → ℂ

/-- A pair of complex modal amplitudes.  The type does not identify them with
physical electric or magnetic fields. -/
abbrev ComplexPairedModalState := ComplexModalVec3 × ComplexModalVec3

/-- Scalar extension of one real OPH modal vector. -/
noncomputable def complexifyModalVector (v : ModalVec3) : ComplexModalVec3 :=
  fun d ↦ (v d : ℂ)

/-- The real chart momentum embedded in the complex modal carrier. -/
noncomputable def complexMomentum (k : ModalVec3) : ComplexModalVec3 :=
  complexifyModalVector k

/-- Complex bilinear modal divergence against the real chart momentum. -/
noncomputable def complexMomentumDot
    (k : ModalVec3) (v : ComplexModalVec3) : ℂ :=
  complexMomentum k ⬝ᵥ v

/-- The complex transverse fibre at chart momentum `k`. -/
def ComplexModeIsTransverse (k : ModalVec3) (v : ComplexModalVec3) : Prop :=
  complexMomentumDot k v = 0

/-- The real squared momentum agrees with the bilinear complex dot product
after scalar extension. -/
theorem complexMomentum_dot_self (k : ModalVec3) :
    complexMomentum k ⬝ᵥ complexMomentum k = (dot k k : ℂ) := by
  simp only [complexMomentum, complexifyModalVector, dot, dotProduct]
  push_cast
  rfl

/-- Fourier curl multiplier `i (omega/|k|) k ×`.  The momentum-dependent
normalization makes this a modal/pseudodifferential construction, not a
position-space locality theorem. -/
noncomputable def fourierCurl
    (a : ℝ) (k : ModalVec3) (v : ComplexModalVec3) : ComplexModalVec3 :=
  (Complex.I * (modalCurlScale a k : ℂ)) •
    crossProduct (complexMomentum k) v

/-- Fourier `div curl = 0`. -/
theorem complexMomentumDot_fourierCurl_zero
    (a : ℝ) (k : ModalVec3) (v : ComplexModalVec3) :
    complexMomentumDot k (fourierCurl a k v) = 0 := by
  unfold complexMomentumDot fourierCurl
  rw [dotProduct_smul, dot_self_cross]
  simp

/-- Fourier curl is complex-linear in its amplitude argument. -/
theorem fourierCurl_smul
    (a : ℝ) (k : ModalVec3) (c : ℂ) (v : ComplexModalVec3) :
    fourierCurl a k (c • v) = c • fourierCurl a k v := by
  unfold fourierCurl
  simp only [map_smul, smul_smul]
  rw [mul_comm]

/-- Fourier curl commutes with amplitude negation. -/
theorem fourierCurl_neg
    (a : ℝ) (k : ModalVec3) (v : ComplexModalVec3) :
    fourierCurl a k (-v) = -fourierCurl a k v := by
  simpa using fourierCurl_smul a k (-1 : ℂ) v

/-- Complexification of the committed FZ-12 scalar spatial action. -/
noncomputable def complexPhotonSpatialAction
    (a : ℝ) (k : ModalVec3) (v : ComplexModalVec3) : ComplexModalVec3 :=
  ((photonModeFrequency a k ^ 2 : ℝ) : ℂ) • v

/-- The complex spatial action is exactly the scalar extension of the
committed real FZ-12 action. -/
theorem complexPhotonSpatialAction_complexifies
    (a : ℝ) (k v : ModalVec3) :
    complexPhotonSpatialAction a k (complexifyModalVector v) =
      complexifyModalVector (photonSpatialAction a k v) := by
  rw [photonSpatialAction_eq_frequency_sq]
  ext d
  simp [complexPhotonSpatialAction, complexifyModalVector]

/-- **Exact Fourier-curl square.**  On the complex transverse fibre, the
standard `i k ×` multiplier squares to the positive complexified FZ-12
spatial action. -/
theorem fourierCurl_sq_on_transverse
    {a : ℝ} {k : ModalVec3} {v : ComplexModalVec3}
    (hk : k ≠ 0) (hv : ComplexModeIsTransverse k v) :
    fourierCurl a k (fourierCurl a k v) =
      complexPhotonSpatialAction a k v := by
  have hscale := modalCurlScale_sq_mul_dot_self (a := a) (k := k) hk
  have htriple :
      crossProduct (complexMomentum k)
          (crossProduct (complexMomentum k) v) =
        complexMomentumDot k v • complexMomentum k -
          (dot k k : ℂ) • v := by
    rw [cross_cross_eq_smul_sub_smul']
    rw [complexMomentum_dot_self]
    rfl
  unfold fourierCurl complexPhotonSpatialAction
  simp only [map_smul, smul_smul]
  rw [htriple, hv]
  ext d
  simp only [zero_smul, zero_sub, Pi.smul_apply, Pi.neg_apply,
    smul_eq_mul]
  have hscaleComplex := congrArg (fun x : ℝ ↦ (x : ℂ)) hscale
  push_cast at hscaleComplex ⊢
  rw [← hscaleComplex]
  calc
    Complex.I * (modalCurlScale a k : ℂ) *
          (Complex.I * (modalCurlScale a k : ℂ)) *
          -((dot k k : ℂ) * v d) =
        -(Complex.I * Complex.I) *
          ((modalCurlScale a k : ℂ) ^ 2 * (dot k k : ℂ)) * v d := by
      ring
    _ = (modalCurlScale a k : ℂ) ^ 2 * (dot k k : ℂ) * v d := by
      rw [Complex.I_mul_I]
      ring

/-- The standard opposite-sign Fourier Maxwell-shaped generator:
`E' = curl B`, `B' = -curl E`.  The prime is still only the supplied
auxiliary modal parameter. -/
noncomputable def maxwellShapedModalGenerator
    (a : ℝ) (k : ModalVec3)
    (state : ComplexPairedModalState) : ComplexPairedModalState :=
  (fourierCurl a k state.2, -fourierCurl a k state.1)

/-- **Exact modal Maxwell-shaped factorization.**  The opposite-sign
first-order pairing squares to the existing second-order FZ-12 wave equation
on both complex transverse amplitudes. -/
theorem maxwellShapedModalGenerator_sq_wave
    {a : ℝ} {k : ModalVec3} (hk : k ≠ 0)
    (state : ComplexPairedModalState)
    (hstate : ComplexModeIsTransverse k state.1 ∧
      ComplexModeIsTransverse k state.2) :
    let second := maxwellShapedModalGenerator a k
      (maxwellShapedModalGenerator a k state)
    second.1 + complexPhotonSpatialAction a k state.1 = 0 ∧
      second.2 + complexPhotonSpatialAction a k state.2 = 0 := by
  dsimp [maxwellShapedModalGenerator]
  rw [fourierCurl_neg]
  rw [fourierCurl_sq_on_transverse hk hstate.1,
    fourierCurl_sq_on_transverse hk hstate.2]
  constructor <;> simp

/-- The Maxwell-shaped generator preserves both complex transverse
constraints algebraically. -/
theorem maxwellShapedModalGenerator_transverse
    (a : ℝ) (k : ModalVec3) (state : ComplexPairedModalState) :
    ComplexModeIsTransverse k
        (maxwellShapedModalGenerator a k state).1 ∧
      ComplexModeIsTransverse k
        (maxwellShapedModalGenerator a k state).2 := by
  constructor
  · exact complexMomentumDot_fourierCurl_zero a k state.2
  · unfold maxwellShapedModalGenerator ComplexModeIsTransverse
    change complexMomentumDot k (-fourierCurl a k state.1) = 0
    unfold complexMomentumDot
    rw [dotProduct_neg]
    change -(complexMomentumDot k (fourierCurl a k state.1)) = 0
    rw [complexMomentumDot_fourierCurl_zero, neg_zero]

/-! ## Same-sign adversarial control -/

/-- Mutation control: use the same sign in both off-diagonal curl blocks. -/
noncomputable def sameSignCurlMutation
    (a : ℝ) (k : ModalVec3)
    (state : ComplexPairedModalState) : ComplexPairedModalState :=
  (fourierCurl a k state.2, fourierCurl a k state.1)

/-- The same-sign mutation squares to the **positive** spatial action, the
wrong sign for the committed oscillator equation. -/
theorem sameSignCurlMutation_sq_positive
    {a : ℝ} {k : ModalVec3} (hk : k ≠ 0)
    (state : ComplexPairedModalState)
    (hstate : ComplexModeIsTransverse k state.1 ∧
      ComplexModeIsTransverse k state.2) :
    sameSignCurlMutation a k (sameSignCurlMutation a k state) =
      (complexPhotonSpatialAction a k state.1,
        complexPhotonSpatialAction a k state.2) := by
  apply Prod.ext
  · exact fourierCurl_sq_on_transverse hk hstate.1
  · exact fourierCurl_sq_on_transverse hk hstate.2

/-- **Adversarial rejection.**  On any transverse state whose first spatial
action is nonzero, the same-sign mutation cannot obey the committed
second-order wave equation. -/
theorem sameSignCurlMutation_fails_wave
    {a : ℝ} {k : ModalVec3} (hk : k ≠ 0)
    (state : ComplexPairedModalState)
    (hstate : ComplexModeIsTransverse k state.1 ∧
      ComplexModeIsTransverse k state.2)
    (haction : complexPhotonSpatialAction a k state.1 ≠ 0) :
    let second := sameSignCurlMutation a k
      (sameSignCurlMutation a k state)
    second.1 + complexPhotonSpatialAction a k state.1 ≠ 0 := by
  dsimp
  rw [show (sameSignCurlMutation a k
      (sameSignCurlMutation a k state)).1 =
      complexPhotonSpatialAction a k state.1 from
    congrArg Prod.fst (sameSignCurlMutation_sq_positive hk state hstate)]
  intro hsum
  apply haction
  funext d
  have hd := congrFun hsum d
  simp only [Pi.add_apply, Pi.zero_apply] at hd
  change complexPhotonSpatialAction a k state.1 d = 0
  exact add_self_eq_zero.mp hd

/-! ## Axiom audit -/

#print axioms modalCurlScale_sq_mul_dot_self
#print axioms dot_modalCurl_zero
#print axioms modalCurl_sq_on_transverse
#print axioms dot_modalCurl_skew
#print axioms complexMomentumDot_fourierCurl_zero
#print axioms complexPhotonSpatialAction_complexifies
#print axioms fourierCurl_sq_on_transverse
#print axioms maxwellShapedModalGenerator_sq_wave
#print axioms maxwellShapedModalGenerator_transverse
#print axioms sameSignCurlMutation_sq_positive
#print axioms sameSignCurlMutation_fails_wave

end OPH.ModalMaxwellFactorizationBoundary
