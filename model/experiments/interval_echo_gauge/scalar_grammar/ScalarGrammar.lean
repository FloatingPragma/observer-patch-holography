import ElectroweakBreakingComposition

/-! The three-address material control spans the exact four directions of the
finite scalar theorem. The third address is an external phase reference.
The source certificate checks the physical contact and degree generators. -/

namespace OPH.NativeScalarGrammar

open OPH.ElectroweakBreakingComposition

/-- A block is specified by its upper/lower diagonal and real/imaginary parts
of its upper off-diagonal. Each coefficient is independently measurable
against the parked reference; the imaginary off-diagonal is minus y. -/
noncomputable def coefficients (g g' : ℝ) (a : EWParam) : ℝ × ℝ × ℝ × ℝ :=
  ((g*a.a3 + g'*a.b)/2, (g'*a.b-g*a.a3)/2, g*a.a1/2, g*a.a2/2)

noncomputable def blockAction (d : ℝ × ℝ × ℝ × ℝ) (h : Doublet) : Doublet :=
  fun j => (![(((d.1 : ℂ) * h 0 +
                  ((d.2.2.1 : ℂ) - Complex.I * (d.2.2.2 : ℂ)) * h 1)),
                ((((d.2.2.1 : ℂ) + Complex.I * (d.2.2.2 : ℂ)) * h 0 +
                  (d.2.1 : ℂ) * h 1))] : Doublet) j * Complex.I

theorem registered_action (g g' : ℝ) (a : EWParam) (h : Doublet) :
    blockAction (coefficients g g' a) h = infAction g g' a h := by
  funext j
  fin_cases j <;> simp [blockAction, coefficients, infAction] <;> ring_nf <;>
    simp [Complex.I_sq] <;> ring

/-- The full real parameter space is reproduced, including the central phase. -/
theorem recover_parameters (g g' : ℝ) (hg : 0 < g) (hg' : 0 < g')
    (u v x y : ℝ) :
    coefficients g g' ⟨2*x/g, 2*y/g, (u-v)/g, (u+v)/g'⟩ = (u,v,x,y) := by
  simp only [coefficients, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> field_simp <;> ring

/-- The finite grammar includes one colorless doublet, positive couplings
and exactly the imported infinitesimal action with its normalized factors. -/
theorem finite_packet (D : EWBreakingPremiseData) :
    (OPH.ExteriorSelection.colorDim 1 = 1 ∧
     OPH.ExteriorSelection.weakDim 1 = 2 ∧ scalarChargeSix = 3) ∧
    (0 < D.gW ∧ 0 < D.gY) ∧
    ∀ (A : EWParam) (h : Doublet),
      blockAction (coefficients D.gW D.gY A) h = infAction D.gW D.gY A h :=
  ⟨⟨scalar_row_receipt.1, scalar_row_receipt.2.1, scalar_row_receipt.2.2.2.2.1⟩,
    ⟨D.gW_pos,D.gY_pos⟩, registered_action D.gW D.gY⟩

end OPH.NativeScalarGrammar
