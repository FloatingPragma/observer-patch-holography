import Mathlib
import A5PortAction
import PortGramRepairCovariance

open scoped BigOperators Matrix

namespace OPH.PortGramA5Isometry

open OPH.PortGramRepairCovariance

/-!
# Proper port rotations act isometrically on the repair-selected carrier

The registered proper icosahedral action consists of sixty explicit
permutations of the twelve port labels.  This file puts those rows into exact
permutation matrices and checks that every row

* is a genuine, non-wrapping permutation of the twelve ports;
* commutes with the declared repair mean and its selected low-band projector;
* preserves the normalized port Gram matrix; and
* restricts to an isometric linear action on the intrinsic rank-three carrier
  `range(P3)`.

The exact finite matrix statements are closed by native evaluation.  The
range restriction and isometry statements are theorem proofs over the real
counting space.  The name `A5` follows the registered port-action module; the
formal object used here is the explicit order-sixty proper-rotation list.
Identifying it with abstract `A5`, selecting this local action as a physical
position action, extending it to the metric completion, and gluing actions
across observer overlaps or refinements are outside this file.
-/

/-- An index into the sixty registered proper-rotation rows. -/
abbrev ProperRotation := Fin 60

/-- The registered permutation row indexed by a proper rotation. -/
def rotationRow (r : ProperRotation) : List Nat :=
  OPH.A5PortAction.perms.get
    ⟨r.val, by rw [OPH.A5PortAction.perms_length]; exact r.isLt⟩

/-- The port map read from a registered row.  Reduction modulo twelve makes
the definition total; `portMap_no_wrap` proves that no registered entry is
actually reduced. -/
def portMap (r : ProperRotation) (i : Fin 12) : Fin 12 :=
  ⟨OPH.A5PortAction.app (rotationRow r) i.val % 12,
    Nat.mod_lt _ (by norm_num)⟩

set_option maxHeartbeats 4000000 in
/-- Every registered row stays in `0,...,11`; the total definition above does
not hide a wrapped label. -/
theorem portMap_no_wrap :
    ∀ r : ProperRotation, ∀ i : Fin 12,
      (portMap r i).val = OPH.A5PortAction.app (rotationRow r) i.val := by
  native_decide

set_option maxHeartbeats 4000000 in
/-- Every registered port map is bijective. -/
theorem portMap_bijective :
    ∀ r : ProperRotation, Function.Bijective (portMap r) := by
  native_decide

/-- Exact permutation matrix of a registered proper rotation.  Column `j`
is sent to the unit vector at `portMap r j`. -/
def rotationQ (r : ProperRotation) : PortMatrix := fun i j ↦
  if i = portMap r j then 1 else 0

set_option maxHeartbeats 8000000 in
/-- Every registered matrix is orthogonal over the exact coefficient field. -/
theorem rotationQ_orthogonal :
    ∀ r : ProperRotation,
      (rotationQ r).transpose * rotationQ r = 1 ∧
        rotationQ r * (rotationQ r).transpose = 1 := by
  native_decide

set_option maxHeartbeats 8000000 in
/-- Every registered proper rotation commutes with the selected low-band
projector. -/
theorem rotationQ_commutes_pLow :
    ∀ r : ProperRotation,
      rotationQ r * pLow = pLow * rotationQ r := by
  native_decide

set_option maxHeartbeats 8000000 in
/-- The complete declared repair mean is equivariant under every registered
proper port rotation. -/
theorem rotationQ_commutes_repairMean :
    ∀ r : ProperRotation,
      rotationQ r * repairMean = repairMean * rotationQ r := by
  native_decide

set_option maxHeartbeats 8000000 in
/-- Every registered proper rotation commutes with the exact normalized port
Gram matrix. -/
theorem rotationQ_commutes_portGram :
    ∀ r : ProperRotation,
      rotationQ r * portGram = portGram * rotationQ r := by
  native_decide

set_option maxHeartbeats 8000000 in
/-- Conjugation by every registered proper rotation preserves the exact
normalized port Gram matrix. -/
theorem rotationQ_preserves_portGram :
    ∀ r : ProperRotation,
      rotationQ r * portGram * (rotationQ r).transpose = portGram := by
  native_decide

set_option maxHeartbeats 12000000 in
/-- The sixty proper rotations remain pairwise distinct after restriction to
the selected rank-three band.  Thus the finite carrier representation loses
no element of the registered rotation group. -/
theorem selected_band_action_faithful :
    ∀ r s : ProperRotation,
      rotationQ r * pLow = rotationQ s * pLow → r = s := by
  native_decide

/-! ## Real action and restriction to the intrinsic carrier -/

/-- Real permutation matrix obtained from the exact registered matrix. -/
noncomputable def rotationR (r : ProperRotation) :
    Matrix (Fin 12) (Fin 12) ℝ :=
  realMatrix (rotationQ r)

/-- Linear action of a registered proper rotation on port counting space. -/
noncomputable def rotationLinear (r : ProperRotation) :
    PortVector →ₗ[ℝ] PortVector :=
  Matrix.toLin' (rotationR r)

/-- The matrix action sends each unit port record to the unit record at the
registered permuted label. -/
theorem rotationLinear_portBasis (r : ProperRotation) (p : Fin 12) :
    rotationLinear r (portBasis p) = portBasis (portMap r p) := by
  funext i
  by_cases h : i = portMap r p <;>
    simp [rotationLinear, rotationR, realMatrix, rotationQ, portBasis,
      Matrix.toLin'_apply, Matrix.mulVec, dotProduct, evalReal, h,
      QuadraticAlgebra.re_one, QuadraticAlgebra.im_one]

/-- Real orthogonality inherited from the exact permutation matrix. -/
theorem rotationR_transpose_mul (r : ProperRotation) :
    (rotationR r).transpose * rotationR r = 1 := by
  have h := congrArg realMatrix (rotationQ_orthogonal r).1
  simpa [rotationR, realMatrix] using h

/-- Real commutation with the repair-selected projector. -/
theorem rotationR_commutes_pLow (r : ProperRotation) :
    rotationR r * pLowR = pLowR * rotationR r := by
  have h := congrArg realMatrix (rotationQ_commutes_pLow r)
  simpa [rotationR, realMatrix] using h

/-- The proper-rotation action commutes with projection onto the intrinsic
carrier. -/
theorem pLowLinear_rotation_commute (r : ProperRotation) (x : PortVector) :
    pLowLinear (rotationLinear r x) =
      rotationLinear r (pLowLinear x) := by
  change pLowR *ᵥ (rotationR r *ᵥ x) =
    rotationR r *ᵥ (pLowR *ᵥ x)
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
    (rotationR_commutes_pLow r).symm]

/-- The induced action permutes the intrinsic port generators according to
the same registered port row. -/
theorem rotationLinear_intrinsicPortVector
    (r : ProperRotation) (p : Fin 12) :
    rotationLinear r (intrinsicPortVector p) =
      intrinsicPortVector (portMap r p) := by
  unfold intrinsicPortVector
  calc
    rotationLinear r (pLowLinear (2 • portBasis p)) =
        pLowLinear (rotationLinear r (2 • portBasis p)) :=
      (pLowLinear_rotation_commute r (2 • portBasis p)).symm
    _ = pLowLinear (2 • portBasis (portMap r p)) := by
      apply congrArg pLowLinear
      rw [map_nsmul, rotationLinear_portBasis]

/-- Every registered proper rotation preserves the counting-space scalar
product. -/
theorem rotationLinear_isometry (r : ProperRotation) (x y : PortVector) :
    portDot (rotationLinear r x) (rotationLinear r y) = portDot x y := by
  change (rotationR r *ᵥ x) ⬝ᵥ (rotationR r *ᵥ y) = x ⬝ᵥ y
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    Matrix.mulVec_mulVec, rotationR_transpose_mul, Matrix.one_mulVec]

/-- Restriction of a registered proper rotation to `range(P3)`. -/
noncomputable def carrierRotation (r : ProperRotation) :
    IntrinsicCarrier →ₗ[ℝ] IntrinsicCarrier where
  toFun x :=
    ⟨rotationLinear r x, by
      rcases x.property with ⟨y, hy⟩
      refine ⟨rotationLinear r y, ?_⟩
      rw [pLowLinear_rotation_commute]
      exact congrArg (rotationLinear r) hy⟩
  map_add' x y := by
    apply Subtype.ext
    exact (rotationLinear r).map_add x y
  map_smul' a x := by
    apply Subtype.ext
    exact (rotationLinear r).map_smul a x

/-- The restricted action carries each intrinsic generator to the generator
at the registered permuted port. -/
theorem carrierRotation_intrinsicPortGenerator
    (r : ProperRotation) (p : Fin 12) :
    carrierRotation r (intrinsicPortGenerator p) =
      intrinsicPortGenerator (portMap r p) := by
  apply Subtype.ext
  exact rotationLinear_intrinsicPortVector r p

/-- The restricted action preserves the intrinsic carrier Gram form exactly. -/
theorem carrierRotation_isometry
    (r : ProperRotation) (x y : IntrinsicCarrier) :
    portDot (carrierRotation r x) (carrierRotation r y) = portDot x y := by
  exact rotationLinear_isometry r x y

/-- Each restricted carrier rotation is injective. -/
theorem carrierRotation_injective (r : ProperRotation) :
    Function.Injective (carrierRotation r) := by
  intro x y hxy
  apply Subtype.ext
  have hambient : rotationLinear r x = rotationLinear r y :=
    congrArg Subtype.val hxy
  change rotationR r *ᵥ (x : PortVector) =
    rotationR r *ᵥ (y : PortVector) at hambient
  have hback := congrArg (fun z ↦ (rotationR r).transpose *ᵥ z) hambient
  simpa only [Matrix.mulVec_mulVec, rotationR_transpose_mul,
    Matrix.one_mulVec] using hback

/-- In particular, every restricted rotation preserves squared distance in
the intrinsic carrier. -/
theorem carrierRotation_dist_sq
    (r : ProperRotation) (x y : IntrinsicCarrier) :
    portDot (carrierRotation r x - carrierRotation r y)
        (carrierRotation r x - carrierRotation r y) =
      portDot (x - y) (x - y) := by
  simpa only [map_sub, Submodule.coe_sub] using
    carrierRotation_isometry r (x - y) (x - y)

#print axioms portMap_no_wrap
#print axioms portMap_bijective
#print axioms rotationQ_orthogonal
#print axioms rotationQ_commutes_pLow
#print axioms rotationQ_commutes_repairMean
#print axioms rotationQ_commutes_portGram
#print axioms rotationQ_preserves_portGram
#print axioms selected_band_action_faithful
#print axioms rotationLinear_portBasis
#print axioms pLowLinear_rotation_commute
#print axioms rotationLinear_intrinsicPortVector
#print axioms rotationLinear_isometry
#print axioms carrierRotation_intrinsicPortGenerator
#print axioms carrierRotation_isometry
#print axioms carrierRotation_injective
#print axioms carrierRotation_dist_sq

end OPH.PortGramA5Isometry
