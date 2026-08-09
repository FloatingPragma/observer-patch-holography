import QFT.SourceContextWeb
import EventAlgebra.Basic

/-!
# Tomography no-go for the current real S3 source context web

Every projector in `SourceContextWeb` is a real symmetric matrix.  The two
pure Pauli-Y states below are distinct, but their opposite imaginary
coherences cancel against every such effect.  Both states therefore assign
the exact weight `1/2` to both outcomes of every diagonal or S3-conjugated
binary context in the current web.

Boundary.  This theorem concerns informational completeness of the declared
real two-dimensional context web.  It does not deny that the web contains
noncommuting real projectors, and it does not apply after a genuinely complex
source context (for example a Pauli-Y projector) has been earned.  No missing
context statistics are imputed to the preregistered run.
-/

namespace OPH.QFT

open Matrix

noncomputable section

/-! ## Two exact states with opposite imaginary coherence -/

/-- Entrywise complexification of a real two-by-two matrix. -/
def complexifyRealMatrix (M : Matrix (Fin 2) (Fin 2) ℝ) :
    Matrix (Fin 2) (Fin 2) ℂ := fun i j => (M i j : ℂ)

/-- Pure `+Y` state. -/
def rhoYPlus : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(1 / 2 : ℂ), -(Complex.I / 2); Complex.I / 2, 1 / 2]

/-- Pure `-Y` state. -/
def rhoYMinus : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(1 / 2 : ℂ), Complex.I / 2; -(Complex.I / 2), 1 / 2]

theorem rhoYPlus_isEvent : EventAlgebra.IsEvent rhoYPlus := by
  constructor
  · show rhoYPlusᴴ = rhoYPlus
    ext i j
    fin_cases i <;> fin_cases j <;> simp [rhoYPlus] <;> ring
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [rhoYPlus, Matrix.mul_apply, Fin.sum_univ_two] <;>
      ring_nf <;> simp [Complex.I_sq] <;> norm_num

theorem rhoYMinus_isEvent : EventAlgebra.IsEvent rhoYMinus := by
  constructor
  · show rhoYMinusᴴ = rhoYMinus
    ext i j
    fin_cases i <;> fin_cases j <;> simp [rhoYMinus] <;> ring
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [rhoYMinus, Matrix.mul_apply, Fin.sum_univ_two] <;>
      ring_nf <;> simp [Complex.I_sq] <;> norm_num

theorem rhoYPlus_isState : EventAlgebra.IsState rhoYPlus := by
  constructor
  · exact rhoYPlus_isEvent.posSemidef
  · norm_num [rhoYPlus, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

theorem rhoYMinus_isState : EventAlgebra.IsState rhoYMinus := by
  constructor
  · exact rhoYMinus_isEvent.posSemidef
  · norm_num [rhoYMinus, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

/-- The states differ exactly in their imaginary off-diagonal entry. -/
theorem rhoYPlus_ne_rhoYMinus : rhoYPlus ≠ rhoYMinus := by
  intro h
  have h01 := congrFun (congrFun h (0 : Fin 2)) (1 : Fin 2)
  have him := congrArg Complex.im h01
  norm_num [rhoYPlus, rhoYMinus] at him

/-! ## Real symmetric effects erase the sign of the Y coherence -/

/-- The `+Y` state assigns weight one half to every real symmetric
trace-one two-by-two matrix. -/
theorem rhoYPlus_bornWeight_real_symmetric_trace_one
    (M : Matrix (Fin 2) (Fin 2) ℝ)
    (hsym : M 0 1 = M 1 0) (htrace : M 0 0 + M 1 1 = 1) :
    EventAlgebra.bornWeight rhoYPlus (complexifyRealMatrix M) = 1 / 2 := by
  simp [EventAlgebra.bornWeight, rhoYPlus, complexifyRealMatrix,
    Matrix.trace, Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]
  rw [hsym]
  have htraceC : (M 0 0 : ℂ) + (M 1 1 : ℂ) = 1 := by
    exact_mod_cast htrace
  linear_combination (1 / 2 : ℂ) * htraceC

/-- The same exact half-weight formula for the `-Y` state. -/
theorem rhoYMinus_bornWeight_real_symmetric_trace_one
    (M : Matrix (Fin 2) (Fin 2) ℝ)
    (hsym : M 0 1 = M 1 0) (htrace : M 0 0 + M 1 1 = 1) :
    EventAlgebra.bornWeight rhoYMinus (complexifyRealMatrix M) = 1 / 2 := by
  simp [EventAlgebra.bornWeight, rhoYMinus, complexifyRealMatrix,
    Matrix.trace, Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]
  rw [hsym]
  have htraceC : (M 0 0 : ℂ) + (M 1 1 : ℂ) = 1 := by
    exact_mod_cast htrace
  linear_combination (1 / 2 : ℂ) * htraceC

theorem conjProjector_zero_eq_record : conjProjector 0 = recordProjector := by
  simp [conjProjector, gaugeIrrep]

/-- Every record-side projector in the S3 web is real symmetric. -/
theorem conjProjector_offdiag_symmetric (g : Fin 6) :
    conjProjector g 0 1 = conjProjector g 1 0 := by
  have hRecord : recordProjectorᵀ = recordProjector := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [recordProjector]
  have hSymm : (conjProjector g)ᵀ = conjProjector g := by
    simp [conjProjector, Matrix.transpose_mul, hRecord, Matrix.mul_assoc]
  have h := congrFun (congrFun hSymm (0 : Fin 2)) (1 : Fin 2)
  simpa [Matrix.transpose_apply] using h.symm

/-- Every record-side projector in the S3 web has trace one. -/
theorem conjProjector_trace_one_real (g : Fin 6) :
    conjProjector g 0 0 + conjProjector g 1 1 = 1 := by
  have hTrace : Matrix.trace (conjProjector g) = 1 := by
    rw [conjProjector, Matrix.trace_mul_cycle,
      (gaugeIrrep_unitary g).2, one_mul]
    norm_num [recordProjector, Matrix.trace, Matrix.diag, Fin.sum_univ_two]
  simpa [Matrix.trace, Matrix.diag, Fin.sum_univ_two] using hTrace

/-- Both states assign exact weight one half to every realized conjugated
record-side effect. -/
theorem rhoYPlus_conjProjector_weight (g : Fin 6) :
    EventAlgebra.bornWeight rhoYPlus
      (complexifyRealMatrix (conjProjector g)) = 1 / 2 :=
  rhoYPlus_bornWeight_real_symmetric_trace_one _
    (conjProjector_offdiag_symmetric g) (conjProjector_trace_one_real g)

theorem rhoYMinus_conjProjector_weight (g : Fin 6) :
    EventAlgebra.bornWeight rhoYMinus
      (complexifyRealMatrix (conjProjector g)) = 1 / 2 :=
  rhoYMinus_bornWeight_real_symmetric_trace_one _
    (conjProjector_offdiag_symmetric g) (conjProjector_trace_one_real g)

/-- Exact equality of record-side weights on the entire S3 web. -/
theorem rhoY_sourceWeb_record_weights_equal (g : Fin 6) :
    EventAlgebra.bornWeight rhoYPlus
      (complexifyRealMatrix (conjProjector g)) =
    EventAlgebra.bornWeight rhoYMinus
      (complexifyRealMatrix (conjProjector g)) := by
  rw [rhoYPlus_conjProjector_weight, rhoYMinus_conjProjector_weight]

/-! ## Both outcomes and the tomography obstruction -/

theorem complexifyRealMatrix_one :
    complexifyRealMatrix (1 : Matrix (Fin 2) (Fin 2) ℝ) = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [complexifyRealMatrix]

theorem complexifyRealMatrix_sub (A B : Matrix (Fin 2) (Fin 2) ℝ) :
    complexifyRealMatrix (A - B) =
      complexifyRealMatrix A - complexifyRealMatrix B := by
  ext i j
  simp [complexifyRealMatrix]

theorem conjCompanionProjector_eq_complement (g : Fin 6) :
    conjCompanionProjector g = 1 - conjProjector g := by
  have h := conjContext_complete g
  apply (eq_sub_iff_add_eq).2
  simpa [add_comm] using h

theorem complexify_conjCompanionProjector (g : Fin 6) :
    complexifyRealMatrix (conjCompanionProjector g) =
      1 - complexifyRealMatrix (conjProjector g) := by
  rw [conjCompanionProjector_eq_complement,
    complexifyRealMatrix_sub, complexifyRealMatrix_one]

/-- Both states also assign exact weight one half to every companion-side
effect of the web. -/
theorem rhoYPlus_conjCompanion_weight (g : Fin 6) :
    EventAlgebra.bornWeight rhoYPlus
      (complexifyRealMatrix (conjCompanionProjector g)) = 1 / 2 := by
  rw [complexify_conjCompanionProjector,
    EventAlgebra.bornWeight_sub, EventAlgebra.bornWeight_one rhoYPlus_isState,
    rhoYPlus_conjProjector_weight]
  norm_num

theorem rhoYMinus_conjCompanion_weight (g : Fin 6) :
    EventAlgebra.bornWeight rhoYMinus
      (complexifyRealMatrix (conjCompanionProjector g)) = 1 / 2 := by
  rw [complexify_conjCompanionProjector,
    EventAlgebra.bornWeight_sub, EventAlgebra.bornWeight_one rhoYMinus_isState,
    rhoYMinus_conjProjector_weight]
  norm_num

theorem rhoY_sourceWeb_companion_weights_equal (g : Fin 6) :
    EventAlgebra.bornWeight rhoYPlus
      (complexifyRealMatrix (conjCompanionProjector g)) =
    EventAlgebra.bornWeight rhoYMinus
      (complexifyRealMatrix (conjCompanionProjector g)) := by
  rw [rhoYPlus_conjCompanion_weight, rhoYMinus_conjCompanion_weight]

/-- Exact half-weights on the separately declared diagonal context. -/
theorem rhoYPlus_diagonal_weights :
    EventAlgebra.bornWeight rhoYPlus
        (complexifyRealMatrix recordProjector) = 1 / 2 ∧
      EventAlgebra.bornWeight rhoYPlus
        (complexifyRealMatrix companionProjector) = 1 / 2 := by
  constructor
  · apply rhoYPlus_bornWeight_real_symmetric_trace_one
    · rfl
    · norm_num [recordProjector]
  · apply rhoYPlus_bornWeight_real_symmetric_trace_one
    · rfl
    · norm_num [companionProjector]

theorem rhoYMinus_diagonal_weights :
    EventAlgebra.bornWeight rhoYMinus
        (complexifyRealMatrix recordProjector) = 1 / 2 ∧
      EventAlgebra.bornWeight rhoYMinus
        (complexifyRealMatrix companionProjector) = 1 / 2 := by
  constructor
  · apply rhoYMinus_bornWeight_real_symmetric_trace_one
    · rfl
    · norm_num [recordProjector]
  · apply rhoYMinus_bornWeight_real_symmetric_trace_one
    · rfl
    · norm_num [companionProjector]

theorem rhoY_sourceWeb_diagonal_weights_equal :
    EventAlgebra.bornWeight rhoYPlus
        (complexifyRealMatrix recordProjector) =
      EventAlgebra.bornWeight rhoYMinus
        (complexifyRealMatrix recordProjector) ∧
    EventAlgebra.bornWeight rhoYPlus
        (complexifyRealMatrix companionProjector) =
      EventAlgebra.bornWeight rhoYMinus
        (complexifyRealMatrix companionProjector) := by
  constructor
  · rw [rhoYPlus_diagonal_weights.1, rhoYMinus_diagonal_weights.1]
  · rw [rhoYPlus_diagonal_weights.2, rhoYMinus_diagonal_weights.2]

/-- The current real S3 context web is not tomographically complete: two
certified, distinct states have identical exact weights on every outcome of
every declared binary context. -/
theorem current_source_context_web_not_tomographically_complete :
    ∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      EventAlgebra.IsState rho ∧ EventAlgebra.IsState sigma ∧ rho ≠ sigma ∧
      (EventAlgebra.bornWeight rho
            (complexifyRealMatrix recordProjector) =
          EventAlgebra.bornWeight sigma
            (complexifyRealMatrix recordProjector) ∧
        EventAlgebra.bornWeight rho
            (complexifyRealMatrix companionProjector) =
          EventAlgebra.bornWeight sigma
            (complexifyRealMatrix companionProjector)) ∧
      (∀ g : Fin 6,
        EventAlgebra.bornWeight rho
            (complexifyRealMatrix (conjProjector g)) =
          EventAlgebra.bornWeight sigma
            (complexifyRealMatrix (conjProjector g)) ∧
        EventAlgebra.bornWeight rho
            (complexifyRealMatrix (conjCompanionProjector g)) =
          EventAlgebra.bornWeight sigma
            (complexifyRealMatrix (conjCompanionProjector g))) := by
  exact ⟨rhoYPlus, rhoYMinus, rhoYPlus_isState, rhoYMinus_isState,
    rhoYPlus_ne_rhoYMinus, rhoY_sourceWeb_diagonal_weights_equal,
    fun g => ⟨rhoY_sourceWeb_record_weights_equal g,
      rhoY_sourceWeb_companion_weights_equal g⟩⟩

/-! ## Exact missing-direction control -/

/-- A genuinely complex `+Y` projector distinguishes the two states
perfectly, showing exactly which operator direction is absent from the real
web. -/
theorem pauliY_context_distinguishes_states :
    EventAlgebra.bornWeight rhoYPlus rhoYPlus = 1 ∧
      EventAlgebra.bornWeight rhoYMinus rhoYPlus = 0 := by
  constructor
  · rw [EventAlgebra.bornWeight, rhoYPlus_isEvent.2, rhoYPlus_isState.2]
  · norm_num [EventAlgebra.bornWeight, rhoYPlus, rhoYMinus, Matrix.trace,
      Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]
    ring_nf
    simp [Complex.I_sq]

#print axioms rhoYPlus_isState
#print axioms rhoYMinus_isState
#print axioms rhoYPlus_ne_rhoYMinus
#print axioms rhoY_sourceWeb_record_weights_equal
#print axioms rhoY_sourceWeb_companion_weights_equal
#print axioms rhoY_sourceWeb_diagonal_weights_equal
#print axioms current_source_context_web_not_tomographically_complete
#print axioms pauliY_context_distinguishes_states

end


end OPH.QFT
