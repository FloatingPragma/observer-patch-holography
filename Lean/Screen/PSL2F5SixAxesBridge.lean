import Mathlib
import A5SixAxes

/-!
# The canonical `SL(2, F5)` cover of the committed six-axis action

This module constructs the projective action from the generic Mathlib action,
descends it through the center quotient, and identifies its image pointwise
with the sixty permutations committed in `A5SixAxes.L60`.

The result stops at the typed six-axis subgroup.  In particular, it does not
identify the abstract group with `A5`, append the pointwise port-row bridge as
a group homomorphism, or transport the Golden-sector anti-representations.
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
