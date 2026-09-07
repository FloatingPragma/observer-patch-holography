import Mathlib
import A5SixAxes
import A5PortSixAxesBridge

/-!
# The canonical `SL(2, F5)` cover and the committed port action

This module constructs the projective action from the generic Mathlib action,
descends it through the center quotient, identifies its image pointwise with
the sixty permutations committed in `A5SixAxes.L60`, and packages the sixty
committed twelve-port rotations as a typed subgroup isomorphic to that image.

The result stops at the typed finite-group interfaces. In particular, it does
not identify the abstract group with `A5`, identify `SL(2,5)` with the binary
icosahedral group, invoke McKay, transport the Golden sectors as typed
representations, derive or select `φ`, state a mass law, or identify physical
rotations.
-/

open scoped LinearAlgebra.Projectivization

namespace OPH.PSL2F5SixAxesBridge

local instance : Fact (Nat.Prime 5) := ⟨by decide⟩

abbrev F5 := ZMod 5
abbrev SL2F5 := Matrix.SpecialLinearGroup (Fin 2) F5
abbrev PSL2F5 := Matrix.ProjectiveSpecialLinearGroup (Fin 2) F5
abbrev P1F5 := ℙ F5 (Fin 2 → F5)

/-- A computational equality instance for the matrix-subtype presentation of
`SL(2, F5)`.  Mathlib deliberately does not export this instance through the
abbreviation. -/
local instance : DecidableEq SL2F5 :=
  inferInstanceAs (DecidableEq
    {A : Matrix (Fin 2) (Fin 2) F5 // A.det = 1})

/-- Decidable center membership, specialized through Mathlib's scalar-matrix
characterization. -/
local instance : DecidablePred (· ∈ Subgroup.center SL2F5) := fun A =>
  decidable_of_iff' (∃ r : F5,
    r ^ Fintype.card (Fin 2) = 1 ∧ Matrix.scalar (Fin 2) r = A)
    Matrix.SpecialLinearGroup.mem_center_iff

/-! ## The canonical quotient -/

/-- The defining quotient from `SL(2, F5)` by its center. -/
def slToPsl : SL2F5 →* PSL2F5 :=
  QuotientGroup.mk' (Subgroup.center SL2F5)

theorem slToPsl_surjective : Function.Surjective slToPsl :=
  QuotientGroup.mk'_surjective _

theorem slToPsl_ker_center : slToPsl.ker = Subgroup.center SL2F5 :=
  QuotientGroup.ker_mk' _

/-! ## Coordinate-preserving projective-line equivalences -/

/-- A linear equivalence induces an equivalence of projectivizations. -/
def projectivizationEquiv {V W : Type*} [AddCommGroup V] [Module F5 V]
    [AddCommGroup W] [Module F5 W] (e : V ≃ₗ[F5] W) :
    ℙ F5 V ≃ ℙ F5 W where
  toFun := Projectivization.map e.toLinearMap e.injective
  invFun := Projectivization.map e.symm.toLinearMap e.symm.injective
  left_inv p := by
    induction p using Projectivization.ind with
    | h v hv =>
        simp only [Projectivization.map_mk]
        congr 1
        exact e.symm_apply_apply v
  right_inv p := by
    induction p using Projectivization.ind with
    | h v hv =>
        simp only [Projectivization.map_mk]
        congr 1
        exact e.apply_symm_apply v

/-- The explicit coordinate identification from the vector presentation to
the usual affine-line-plus-infinity presentation. -/
def p1EquivOnePoint : P1F5 ≃ OnePoint F5 :=
  (projectivizationEquiv (LinearEquiv.finTwoArrow F5 F5)).trans
    (OnePoint.equivProjectivization F5).symm

/-- Put the five affine `ZMod 5` coordinates first and infinity at index `5`. -/
def onePointEquivSix : OnePoint F5 ≃ Fin 6 :=
  (Equiv.optionCongr (ZMod.finEquiv 5).symm.toEquiv).trans
    finSuccEquivLast.symm

/-- The load-bearing relabelling `[z:1] ↔ z`, `[1:0] ↔ 5`. -/
def p1EquivSix : P1F5 ≃ Fin 6 :=
  p1EquivOnePoint.trans onePointEquivSix

@[simp]
theorem onePointEquivSix_affine (z : F5) :
    onePointEquivSix (z : OnePoint F5) =
      Fin.castSucc ((ZMod.finEquiv 5).symm z) := by
  change finSuccEquivLast.symm (some ((ZMod.finEquiv 5).symm z)) = _
  simp

@[simp]
theorem onePointEquivSix_infinity :
    onePointEquivSix (OnePoint.infty : OnePoint F5) = Fin.last 5 := by
  change finSuccEquivLast.symm none = _
  simp

@[simp]
theorem p1EquivOnePoint_affine (z : F5) :
    p1EquivOnePoint (Projectivization.mk F5 ![z, 1] (by simp)) = z := by
  change (OnePoint.equivProjectivization F5).symm
      (Projectivization.map
        (LinearEquiv.finTwoArrow F5 F5).toLinearMap
        (LinearEquiv.finTwoArrow F5 F5).injective
        (Projectivization.mk F5 ![z, 1] (by simp))) = z
  rw [Projectivization.map_mk]
  simp

@[simp]
theorem p1EquivOnePoint_infinity :
    p1EquivOnePoint (Projectivization.mk F5 ![1, 0] (by simp)) =
      (OnePoint.infty : OnePoint F5) := by
  change (OnePoint.equivProjectivization F5).symm
      (Projectivization.map
        (LinearEquiv.finTwoArrow F5 F5).toLinearMap
        (LinearEquiv.finTwoArrow F5 F5).injective
        (Projectivization.mk F5 ![1, 0] (by simp))) = OnePoint.infty
  rw [Projectivization.map_mk]
  simp

@[simp]
theorem p1EquivSix_affine (z : F5) :
    p1EquivSix (Projectivization.mk F5 ![z, 1] (by simp)) =
      Fin.castSucc ((ZMod.finEquiv 5).symm z) := by
  simp [p1EquivSix]

@[simp]
theorem p1EquivSix_infinity :
    p1EquivSix (Projectivization.mk F5 ![1, 0] (by simp)) = Fin.last 5 := by
  simp [p1EquivSix]

/-! ## The locally constructed projective action -/

local instance : MulAction SL2F5 (OnePoint F5) :=
  MulAction.compHom (OnePoint F5) Matrix.SpecialLinearGroup.toGL

/-- The `SL(2, F5)` action on the affine line plus infinity, obtained by
restricting Mathlib's generic `GL(2)` projective action. -/
def slOnePointAction : SL2F5 →* Equiv.Perm (OnePoint F5) :=
  MulAction.toPermHom SL2F5 (OnePoint F5)

/-- The same locally constructed action on `P¹(F5)`. -/
def slProjectiveAction : SL2F5 →* Equiv.Perm P1F5 :=
  p1EquivOnePoint.symm.permCongrHom.toMonoidHom.comp slOnePointAction

/-- The concrete six-coordinate form of the projective action. -/
def slToSix : SL2F5 →* Equiv.Perm (Fin 6) :=
  p1EquivSix.permCongrHom.toMonoidHom.comp slProjectiveAction

@[simp]
theorem slToSix_apply (A : SL2F5) (x : P1F5) :
    slToSix A (p1EquivSix x) = p1EquivSix (slProjectiveAction A x) := by
  simp [slToSix]

/-! ## Kernel, center, and standard generators -/

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
/-- The concrete action has exactly the scalar center as kernel.  This is a
closed exhaustive check over the 120 determinant-one matrices, not an order
comparison. -/
theorem slToSix_ker : slToSix.ker = Subgroup.center SL2F5 := by
  ext A
  exact (by decide : ∀ A : SL2F5,
    A ∈ slToSix.ker ↔ A ∈ Subgroup.center SL2F5) A

/-- Faithful coordinate transport does not alter the kernel. -/
theorem slProjectiveAction_ker :
    slProjectiveAction.ker = Subgroup.center SL2F5 := by
  ext A
  have hconj : A ∈ slProjectiveAction.ker ↔ A ∈ slToSix.ker := by
    change slProjectiveAction A = 1 ↔
      p1EquivSix.permCongrHom (slProjectiveAction A) = 1
    constructor
    · intro hA
      rw [hA, map_one]
    · intro hA
      exact p1EquivSix.permCongrHom.injective
        (hA.trans (map_one p1EquivSix.permCongrHom).symm)
  rw [hconj, slToSix_ker]

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
/-- The center consists exactly of the two scalar matrices `+I` and `-I`. -/
theorem center_mem_iff_plus_minus_one (A : SL2F5) :
    A ∈ Subgroup.center SL2F5 ↔ A = 1 ∨ A = -1 := by
  exact (by decide : ∀ A : SL2F5,
    A ∈ Subgroup.center SL2F5 ↔ A = 1 ∨ A = -1) A

theorem center_eq_plus_minus_one :
    (Subgroup.center SL2F5 : Set SL2F5) = ({1, -1} : Set SL2F5) := by
  ext A
  simpa [Set.mem_insert_iff, Set.mem_singleton_iff] using
    center_mem_iff_plus_minus_one A

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
theorem center_card_two : Fintype.card (Subgroup.center SL2F5) = 2 := by
  decide

/-- The determinant-one lift of `z ↦ z + 1`. -/
def tLift : SL2F5 :=
  ⟨!![1, 1; 0, 1], by decide⟩

/-- The determinant-one lift of `z ↦ -1/z`. -/
def sLift : SL2F5 :=
  ⟨!![0, -1; 1, 0], by decide⟩

theorem tLift_action : slToSix tLift = OPH.A5SixAxes.t := by
  decide

theorem sLift_action : slToSix sLift = OPH.A5SixAxes.s := by
  decide

/-! ## Descent to `PSL(2, F5)` -/

/-- The projective action descended through the defining center quotient. -/
def pslProjectiveAction : PSL2F5 →* Equiv.Perm P1F5 :=
  QuotientGroup.lift (Subgroup.center SL2F5) slProjectiveAction
    slProjectiveAction_ker.symm.le

/-- The descended action in the committed six-coordinate convention. -/
def pslToSix : PSL2F5 →* Equiv.Perm (Fin 6) :=
  QuotientGroup.lift (Subgroup.center SL2F5) slToSix slToSix_ker.symm.le

@[simp]
theorem pslProjectiveAction_slToPsl (A : SL2F5) :
    pslProjectiveAction (slToPsl A) = slProjectiveAction A :=
  rfl

@[simp]
theorem pslToSix_slToPsl (A : SL2F5) :
    pslToSix (slToPsl A) = slToSix A :=
  rfl

theorem pslProjectiveAction_injective :
    Function.Injective pslProjectiveAction := by
  rw [← MonoidHom.ker_eq_bot_iff]
  rw [pslProjectiveAction, QuotientGroup.ker_lift,
    slProjectiveAction_ker, QuotientGroup.map_mk'_self]

theorem pslToSix_injective : Function.Injective pslToSix := by
  rw [← MonoidHom.ker_eq_bot_iff]
  rw [pslToSix, QuotientGroup.ker_lift, slToSix_ker,
    QuotientGroup.map_mk'_self]

theorem psl_actions_intertwine (q : PSL2F5) (x : P1F5) :
    pslToSix q (p1EquivSix x) =
      p1EquivSix (pslProjectiveAction q x) := by
  obtain ⟨A, rfl⟩ := slToPsl_surjective q
  exact slToSix_apply A x

theorem tClass_action : pslToSix (slToPsl tLift) = OPH.A5SixAxes.t := by
  simpa using tLift_action

theorem sClass_action : pslToSix (slToPsl sLift) = OPH.A5SixAxes.s := by
  simpa using sLift_action

/-! ## The committed subgroup and exact image -/

/-- The subgroup whose carrier is exactly the committed list `L60`. -/
def SixAxisGroup : Subgroup (Equiv.Perm (Fin 6)) where
  carrier := {g | g ∈ OPH.A5SixAxes.L60}
  one_mem' := OPH.A5SixAxes.one_mem
  mul_mem' hg hh := OPH.A5SixAxes.mul_closed _ hg _ hh
  inv_mem' hg := OPH.A5SixAxes.inv_closed _ hg

@[simp]
theorem mem_sixAxisGroup_iff (g : Equiv.Perm (Fin 6)) :
    g ∈ SixAxisGroup ↔ g ∈ OPH.A5SixAxes.L60 :=
  Iff.rfl

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
/-- Every determinant-one matrix acts by one of the committed sixty rows. -/
theorem slToSix_mem_L60 :
    ∀ A : SL2F5, slToSix A ∈ OPH.A5SixAxes.L60 := by
  decide

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
/-- Every committed row has an explicit determinant-one preimage.  The
existential witnesses are found by a closed exhaustive check, so this is a
concrete image theorem rather than a comparison of group orders. -/
theorem every_L60_row_has_sl_preimage :
    ∀ i : Fin 60, ∃ A : SL2F5, slToSix A = OPH.A5SixAxes.el i := by
  decide

/-- The image of the abstract projective group is exactly the committed
sixty-row subgroup. -/
theorem pslToSix_range : pslToSix.range = SixAxisGroup := by
  ext g
  constructor
  · rintro ⟨q, rfl⟩
    obtain ⟨A, rfl⟩ := slToPsl_surjective q
    exact slToSix_mem_L60 A
  · intro hg
    obtain ⟨i, rfl⟩ := OPH.A5SixAxes.mem_iff_el hg
    obtain ⟨A, hA⟩ := every_L60_row_has_sl_preimage i
    exact ⟨slToPsl A, hA⟩

/-- The descended action with codomain restricted to the exact committed
subgroup. -/
def pslToSixAxisGroup : PSL2F5 →* SixAxisGroup :=
  pslToSix.codRestrict SixAxisGroup fun q => by
    rw [← pslToSix_range]
    exact ⟨q, rfl⟩

theorem pslToSixAxisGroup_injective :
    Function.Injective pslToSixAxisGroup := by
  intro q r h
  apply pslToSix_injective
  exact congrArg Subtype.val h

theorem pslToSixAxisGroup_surjective :
    Function.Surjective pslToSixAxisGroup := by
  intro g
  have hg : (g : Equiv.Perm (Fin 6)) ∈ pslToSix.range := by
    rw [pslToSix_range]
    exact g.property
  obtain ⟨q, hq⟩ := hg
  exact ⟨q, Subtype.ext hq⟩

/-- The abstract Mathlib quotient is concretely isomorphic to the subgroup
of the sixty committed six-axis permutations. -/
noncomputable def psl_equiv_six_axis_group : PSL2F5 ≃* SixAxisGroup :=
  MulEquiv.ofBijective pslToSixAxisGroup
    ⟨pslToSixAxisGroup_injective, pslToSixAxisGroup_surjective⟩

#print axioms slToPsl_ker_center
#print axioms slProjectiveAction_ker
#print axioms pslProjectiveAction_injective
#print axioms center_eq_plus_minus_one
#print axioms pslToSix_range
#print axioms psl_equiv_six_axis_group

end OPH.PSL2F5SixAxesBridge

namespace OPH.A5PortGroupBridge

open OPH.A5PortSixAxesBridge
open OPH.PSL2F5SixAxesBridge

/-!
# The committed twelve-port rotations as a typed finite group

The sixty committed twelve-port rotations are packaged as a subgroup of
`Equiv.Perm (Fin 12)` and identified with the committed six-axis subgroup.

The load-bearing check is multiplication on the actual `Fin 12`
permutations. The port-row multiplication and inverse indices are transported
through the certified `rowEquiv` into the exact multiplication and inverse
tables of `A5SixAxes`; the resulting identities are then checked pointwise on
the committed port permutations. `PortGroup` is defined as the range of those
sixty rows, so it contains exactly the committed rotations.

Composing this group isomorphism with `psl_equiv_six_axis_group` gives

    PSL2F5 ≃* PortGroup.

BOUNDARY. This section does not prove `PSL(2,5) ≅ A5`, identify `SL(2,5)` with
the binary icosahedral group, construct an `SU(2)` lift, invoke McKay, transport
the golden sectors as typed representations, derive or select `φ`, state a
mass law, or identify physical rotations. -/

/-! ## 1. Transport the six-axis row law back to the committed port rows -/

/-- Multiplication index for the port rows, transported through the certified
row equivalence to the exact six-axis multiplication table. -/
def portMulIndex (i j : Fin 60) : Fin 60 :=
  rowEquiv.symm (OPH.A5SixAxes.mulT (rowEquiv i) (rowEquiv j))

/-- Inverse index for a port row, transported through the certified row
equivalence to the exact six-axis inverse table. -/
def portInvIndex (i : Fin 60) : Fin 60 :=
  rowEquiv.symm (OPH.A5SixAxes.invT (rowEquiv i))

@[simp]
theorem rowEquiv_portMulIndex (i j : Fin 60) :
    rowEquiv (portMulIndex i j) =
      OPH.A5SixAxes.mulT (rowEquiv i) (rowEquiv j) := by
  simp [portMulIndex]

@[simp]
theorem rowEquiv_portInvIndex (i : Fin 60) :
    rowEquiv (portInvIndex i) = OPH.A5SixAxes.invT (rowEquiv i) := by
  simp [portInvIndex]

@[simp]
theorem rowEquiv_zero : rowEquiv (0 : Fin 60) = 0 := by decide

/-- Row zero is the identity port permutation. -/
theorem portEl_zero : portEl 0 = 1 := by
  apply Equiv.ext
  intro k
  revert k
  decide

/- The quadratic check is split into the same five-row bands used by
`A5SixAxes.rowF_mul`. This keeps kernel reduction on raw action values rather
than asking `decide` to compare `Equiv.Perm` proof structures. -/

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band0 :
    ∀ (i j : Fin 60) (k : Fin 12), i.val < 5 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band1 :
    ∀ (i j : Fin 60) (k : Fin 12), 5 ≤ i.val → i.val < 10 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band2 :
    ∀ (i j : Fin 60) (k : Fin 12), 10 ≤ i.val → i.val < 15 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band3 :
    ∀ (i j : Fin 60) (k : Fin 12), 15 ≤ i.val → i.val < 20 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band4 :
    ∀ (i j : Fin 60) (k : Fin 12), 20 ≤ i.val → i.val < 25 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band5 :
    ∀ (i j : Fin 60) (k : Fin 12), 25 ≤ i.val → i.val < 30 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band6 :
    ∀ (i j : Fin 60) (k : Fin 12), 30 ≤ i.val → i.val < 35 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band7 :
    ∀ (i j : Fin 60) (k : Fin 12), 35 ≤ i.val → i.val < 40 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band8 :
    ∀ (i j : Fin 60) (k : Fin 12), 40 ≤ i.val → i.val < 45 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band9 :
    ∀ (i j : Fin 60) (k : Fin 12), 45 ≤ i.val → i.val < 50 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band10 :
    ∀ (i j : Fin 60) (k : Fin 12), 50 ≤ i.val → i.val < 55 →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem portEl_mul_band11 :
    ∀ (i j : Fin 60) (k : Fin 12), 55 ≤ i.val →
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by decide

/-- Pointwise multiplication law for all committed port rows. -/
theorem portEl_mul_apply :
    ∀ (i j : Fin 60) (k : Fin 12),
      (portEl i * portEl j) k = portEl (portMulIndex i j) k := by
  intro i j k
  rcases Nat.lt_or_ge i.val 5 with h | h0
  · exact portEl_mul_band0 i j k h
  rcases Nat.lt_or_ge i.val 10 with h | h1
  · exact portEl_mul_band1 i j k h0 h
  rcases Nat.lt_or_ge i.val 15 with h | h2
  · exact portEl_mul_band2 i j k h1 h
  rcases Nat.lt_or_ge i.val 20 with h | h3
  · exact portEl_mul_band3 i j k h2 h
  rcases Nat.lt_or_ge i.val 25 with h | h4
  · exact portEl_mul_band4 i j k h3 h
  rcases Nat.lt_or_ge i.val 30 with h | h5
  · exact portEl_mul_band5 i j k h4 h
  rcases Nat.lt_or_ge i.val 35 with h | h6
  · exact portEl_mul_band6 i j k h5 h
  rcases Nat.lt_or_ge i.val 40 with h | h7
  · exact portEl_mul_band7 i j k h6 h
  rcases Nat.lt_or_ge i.val 45 with h | h8
  · exact portEl_mul_band8 i j k h7 h
  rcases Nat.lt_or_ge i.val 50 with h | h9
  · exact portEl_mul_band9 i j k h8 h
  rcases Nat.lt_or_ge i.val 55 with h | h10
  · exact portEl_mul_band10 i j k h9 h
  exact portEl_mul_band11 i j k h10

/-- Multiplication of the actual twelve-port permutations is exactly the
transported six-axis multiplication table. -/
theorem portEl_mul (i j : Fin 60) :
    portEl i * portEl j = portEl (portMulIndex i j) := by
  apply Equiv.ext
  exact portEl_mul_apply i j

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Pointwise inverse law on the committed port rows. -/
theorem portEl_inv_apply :
    ∀ (i : Fin 60) (k : Fin 12),
      (portEl i)⁻¹ k = portEl (portInvIndex i) k := by decide

/-- Inversion of the actual twelve-port permutations is exactly the
transported six-axis inverse table. -/
theorem portEl_inv (i : Fin 60) :
    (portEl i)⁻¹ = portEl (portInvIndex i) := by
  apply Equiv.ext
  exact portEl_inv_apply i

/-! ## 2. The exact subgroup carried by the committed port rows -/

/-- The subgroup consisting exactly of the sixty committed twelve-port
rotations. Range membership makes the row ancestry explicit. -/
def PortGroup : Subgroup (Equiv.Perm (Fin 12)) where
  carrier := Set.range portEl
  one_mem' := ⟨0, portEl_zero⟩
  mul_mem' := by
    intro g h hg hh
    obtain ⟨i, rfl⟩ := hg
    obtain ⟨j, rfl⟩ := hh
    exact ⟨portMulIndex i j, (portEl_mul i j).symm⟩
  inv_mem' := by
    intro g hg
    obtain ⟨i, rfl⟩ := hg
    exact ⟨portInvIndex i, (portEl_inv i).symm⟩

/-- Distinct row indices give distinct committed port permutations. The proof
reuses the certified faithfulness of the antipodal quotient. -/
theorem portEl_injective : Function.Injective portEl := by
  intro i j hij
  apply quotient_action_faithful
  funext x
  simp only [quotientAxis]
  rw [hij]

/-- Choose the unique committed row underlying a `PortGroup` element. -/
noncomputable def portIndex (g : PortGroup) : Fin 60 :=
  Classical.choose (show ∃ i : Fin 60, portEl i = g.1 from g.property)

/-- The chosen row is the underlying port permutation. -/
theorem portIndex_spec (g : PortGroup) : portEl (portIndex g) = g.1 :=
  Classical.choose_spec (show ∃ i : Fin 60, portEl i = g.1 from g.property)

/-- The range witness is unique, so choosing the index of a committed row
returns that row's index. -/
theorem portIndex_portEl (i : Fin 60) :
    portIndex (⟨portEl i, ⟨i, rfl⟩⟩ : PortGroup) = i := by
  apply portEl_injective
  exact portIndex_spec _

/-- The identity element has row index zero. -/
theorem portIndex_one : portIndex (1 : PortGroup) = 0 := by
  apply portEl_injective
  calc
    portEl (portIndex (1 : PortGroup)) = (1 : PortGroup).1 := portIndex_spec _
    _ = 1 := rfl
    _ = portEl 0 := portEl_zero.symm

/-- The unique row index of a product is the transported multiplication
index. -/
theorem portIndex_mul (g h : PortGroup) :
    portIndex (g * h) = portMulIndex (portIndex g) (portIndex h) := by
  apply portEl_injective
  calc
    portEl (portIndex (g * h)) = (g * h).1 := portIndex_spec _
    _ = g.1 * h.1 := rfl
    _ = portEl (portIndex g) * portEl (portIndex h) := by
      rw [portIndex_spec g, portIndex_spec h]
    _ = portEl (portMulIndex (portIndex g) (portIndex h)) :=
      portEl_mul _ _

/-! ## 3. The typed homomorphism to the exact six-axis subgroup -/

/-- Row zero is the identity in the explicit six-axis list. -/
theorem sixEl_zero : OPH.A5SixAxes.el 0 = 1 := by
  apply Equiv.ext
  intro x
  rw [OPH.A5SixAxes.el_apply, Equiv.Perm.one_apply]
  revert x
  decide

/-- Multiplication of explicit six-axis rows, exposed from the certified raw
row table. -/
theorem sixEl_mul (i j : Fin 60) :
    OPH.A5SixAxes.el i * OPH.A5SixAxes.el j =
      OPH.A5SixAxes.el (OPH.A5SixAxes.mulT i j) := by
  apply Equiv.ext
  intro x
  rw [Equiv.Perm.mul_apply, OPH.A5SixAxes.el_apply,
    OPH.A5SixAxes.el_apply, OPH.A5SixAxes.el_apply]
  exact OPH.A5SixAxes.rowF_mul i j x

/-- Send a committed port rotation to the corresponding committed six-axis
rotation. -/
noncomputable def portToSix : PortGroup →*
    OPH.PSL2F5SixAxesBridge.SixAxisGroup where
  toFun g :=
    ⟨OPH.A5SixAxes.el (rowEquiv (portIndex g)), OPH.A5SixAxes.el_mem _⟩
  map_one' := by
    apply Subtype.ext
    change OPH.A5SixAxes.el (rowEquiv (portIndex (1 : PortGroup))) = 1
    rw [portIndex_one, rowEquiv_zero, sixEl_zero]
  map_mul' g h := by
    apply Subtype.ext
    change OPH.A5SixAxes.el (rowEquiv (portIndex (g * h))) =
      OPH.A5SixAxes.el (rowEquiv (portIndex g)) *
        OPH.A5SixAxes.el (rowEquiv (portIndex h))
    rw [portIndex_mul, rowEquiv_portMulIndex]
    exact (sixEl_mul _ _).symm

/-- The port-to-six-axis homomorphism is injective. -/
theorem portToSix_injective : Function.Injective portToSix := by
  intro g h hgh
  have hval :
      OPH.A5SixAxes.el (rowEquiv (portIndex g)) =
        OPH.A5SixAxes.el (rowEquiv (portIndex h)) :=
    congrArg Subtype.val hgh
  have hrow :
      OPH.A5SixAxes.rowF (rowEquiv (portIndex g)) =
        OPH.A5SixAxes.rowF (rowEquiv (portIndex h)) := by
    funext x
    simpa [OPH.A5SixAxes.el_apply] using
      congrArg (fun e : Equiv.Perm (Fin 6) => e x) hval
  have hidx : rowEquiv (portIndex g) = rowEquiv (portIndex h) :=
    six_axis_rows_injective hrow
  have hpidx : portIndex g = portIndex h := rowEquiv.injective hidx
  apply Subtype.ext
  calc
    g.1 = portEl (portIndex g) := (portIndex_spec g).symm
    _ = portEl (portIndex h) := congrArg portEl hpidx
    _ = h.1 := portIndex_spec h

/-- Every committed six-axis group element comes from one committed port
rotation. -/
theorem portToSix_surjective : Function.Surjective portToSix := by
  intro g
  obtain ⟨j, hj⟩ := OPH.A5SixAxes.mem_iff_el g.property
  let i : Fin 60 := rowEquiv.symm j
  let p : PortGroup := ⟨portEl i, ⟨i, rfl⟩⟩
  refine ⟨p, ?_⟩
  apply Subtype.ext
  change OPH.A5SixAxes.el (rowEquiv (portIndex p)) = g.1
  have hpi : portIndex p = i := by
    dsimp [p]
    exact portIndex_portEl i
  rw [hpi]
  change OPH.A5SixAxes.el (rowEquiv (rowEquiv.symm j)) = g.1
  rw [rowEquiv.apply_symm_apply]
  exact hj.symm

/-- The actual twelve-port rotation subgroup is isomorphic to the exact
committed six-axis subgroup. -/
noncomputable def portGroupEquivSixAxisGroup :
    PortGroup ≃* OPH.PSL2F5SixAxesBridge.SixAxisGroup :=
  MulEquiv.ofBijective portToSix
    ⟨portToSix_injective, portToSix_surjective⟩

/-- Mathlib's abstract projective group is isomorphic to the subgroup of the
sixty committed twelve-port rotations. -/
noncomputable def pslEquivPortGroup :
    OPH.PSL2F5SixAxesBridge.PSL2F5 ≃* PortGroup :=
  OPH.PSL2F5SixAxesBridge.psl_equiv_six_axis_group.trans
    portGroupEquivSixAxisGroup.symm

end OPH.A5PortGroupBridge

/- Axiom audit: no `sorry`, `admit`, new axioms, or `native_decide`. -/

#print axioms OPH.A5PortGroupBridge.portEl_zero
#print axioms OPH.A5PortGroupBridge.portEl_mul
#print axioms OPH.A5PortGroupBridge.portEl_inv
#print axioms OPH.A5PortGroupBridge.portEl_injective
#print axioms OPH.A5PortGroupBridge.portToSix_injective
#print axioms OPH.A5PortGroupBridge.portToSix_surjective
#print axioms OPH.A5PortGroupBridge.portGroupEquivSixAxisGroup
#print axioms OPH.A5PortGroupBridge.pslEquivPortGroup
