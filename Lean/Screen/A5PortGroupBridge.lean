import A5PortSixAxesBridge
import PSL2F5SixAxesBridge

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

BOUNDARY. This file does not prove `PSL(2,5) ≅ A5`, identify `SL(2,5)` with
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
