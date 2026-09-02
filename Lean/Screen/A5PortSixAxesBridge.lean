import A5PortModule
import A5SixAxes

namespace OPH.A5PortSixAxesBridge

/-! # The antipodal port quotient and the explicit six-axis action

The sixty permutations in `A5PortAction` commute with the antipode on the
twelve ports.  This file quotients those ports into the six antipodal pairs
`{0,11}, {1,10}, ..., {5,6}` and checks the resulting finite action against
the repository's explicit six-axis `PSL(2,F5)` model in `A5SixAxes`.

The comparison is direct: one fixed relabeling of the six axes and one
explicit bijection of the sixty row indices make all `60 * 6` action values
equal.  The quotient action is also faithful, so no nontrivial listed port
rotation is hidden by passage to antipodal pairs.

This proves only an equivalence of two committed finite actions.  It does not
identify either action with Mathlib's abstract `PSL(2, ZMod 5)`, prove
`PSL(2,5) ≅ A5`, construct `2I ≅ SL(2,5)`, formalize McKay, derive `φ` from
`E8`, select `27^φ`, prove a charged-lepton mass law, identify physical
rotations, or introduce a physical observable. -/

/-- The antipodal-pair label of a port.  The literal table makes the six
committed pairs explicit. -/
def axisOfPort : Fin 12 → Fin 6 :=
  ![0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0]

/-- Antipode on the twelve port labels. -/
def antipode (k : Fin 12) : Fin 12 :=
  ⟨11 - k.val, by omega⟩

/-- Row `i` of the genuine port permutations supplied by `A5PortModule`.
That module proves row-for-row equality with `A5PortAction.perms`. -/
def portEl (i : Fin 60) : Equiv.Perm (Fin 12) :=
  OPH.A5PortModule.P60.get (Fin.cast (by decide) i)

/-- Every listed port rotation induces a well-defined map on antipodal
pairs: either representative has the same image-axis label. -/
theorem quotient_respects_antipode : ∀ (i : Fin 60) (k : Fin 12),
    axisOfPort (portEl i (antipode k)) = axisOfPort (portEl i k) := by
  decide

/-- The six-axis quotient action, using ports `0,...,5` as representatives. -/
def quotientAxis (i : Fin 60) (x : Fin 6) : Fin 6 :=
  axisOfPort (portEl i ⟨x.val, by omega⟩)

/-- The fixed relabeling `[0,1,2,4,5,3]` from antipodal-pair labels to the
coordinates used by `A5SixAxes`. -/
def axisRelabel : Equiv.Perm (Fin 6) :=
  OPH.A5SixAxes.perm ![0, 1, 2, 4, 5, 3] ![0, 1, 2, 5, 3, 4]

theorem axisRelabel_left_inv (x : Fin 6) :
    axisRelabel.symm (axisRelabel x) = x := by
  exact axisRelabel.symm_apply_apply x

theorem axisRelabel_right_inv (x : Fin 6) :
    axisRelabel (axisRelabel.symm x) = x := by
  exact axisRelabel.apply_symm_apply x

/-- Conjugate the antipodal quotient by the fixed six-axis relabeling. -/
def bridgedAxis (i : Fin 60) (x : Fin 6) : Fin 6 :=
  axisRelabel (quotientAxis i (axisRelabel.symm x))

set_option maxRecDepth 16384 in
/-- Explicit row matching from port-action rows to `A5SixAxes.rowF` rows.
The inverse table is supplied as part of the equivalence, and both inverse
laws are kernel checked. -/
def rowEquiv : Fin 60 ≃ Fin 60 :=
  ⟨![0, 3, 6, 9, 4, 11, 12, 17, 14, 18,
      20, 23, 28, 25, 27, 41, 42, 47, 48, 45,
      50, 55, 56, 58, 52, 33, 34, 38, 37, 30,
      31, 36, 39, 35, 32, 53, 59, 57, 54, 51,
      44, 49, 46, 43, 40, 26, 24, 29, 22, 21,
      19, 15, 16, 13, 10, 5, 8, 7, 2, 1],
    ![0, 59, 58, 1, 4, 55, 2, 57, 56, 3,
      54, 5, 6, 53, 8, 51, 52, 7, 9, 50,
      10, 49, 48, 11, 46, 13, 45, 14, 12, 47,
      29, 30, 34, 25, 26, 33, 31, 28, 27, 32,
      44, 15, 16, 43, 40, 19, 42, 17, 18, 41,
      20, 39, 24, 35, 38, 21, 22, 37, 23, 36],
    by decide,
    by decide⟩

/-- Headline bridge: after the explicit axis relabeling and row bijection,
the complete antipodal quotient action is pointwise the repository's
explicit six-axis `PSL(2,F5)` action. -/
theorem bridged_axis_eq_six_axis : ∀ (i : Fin 60) (x : Fin 6),
    bridgedAxis i x = OPH.A5SixAxes.rowF (rowEquiv i) x := by
  decide

/-- Every row of the repository's six-axis model is realized by a row of the
antipodal port quotient. -/
theorem every_six_axis_row_realized (j : Fin 60) :
    ∃ i : Fin 60, ∀ x : Fin 6,
      bridgedAxis i x = OPH.A5SixAxes.rowF j x := by
  refine ⟨rowEquiv.symm j, fun x => ?_⟩
  rw [bridged_axis_eq_six_axis, rowEquiv.apply_symm_apply]

set_option maxRecDepth 16384 in
/-- The sixty raw value rows in `A5SixAxes` are pairwise distinct. -/
theorem six_axis_rows_injective : Function.Injective OPH.A5SixAxes.rowF := by
  decide

/-- Distinct listed port rotations remain distinct on the six antipodal
axes; the quotient representation has no hidden kernel. -/
theorem quotient_action_faithful : Function.Injective quotientAxis := by
  intro i j hij
  apply rowEquiv.injective
  apply six_axis_rows_injective
  funext x
  rw [← bridged_axis_eq_six_axis, ← bridged_axis_eq_six_axis]
  simp only [bridgedAxis]
  rw [hij]

end OPH.A5PortSixAxesBridge

/- Axiom audit: standard axioms only; no `sorry`, `admit`, or new axioms. -/

#print axioms OPH.A5PortSixAxesBridge.quotient_respects_antipode
#print axioms OPH.A5PortSixAxesBridge.bridged_axis_eq_six_axis
#print axioms OPH.A5PortSixAxesBridge.every_six_axis_row_realized
#print axioms OPH.A5PortSixAxesBridge.quotient_action_faithful
