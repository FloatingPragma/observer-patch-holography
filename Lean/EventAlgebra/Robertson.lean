import EventAlgebra.StateExpectation

/-!
# Robertson inequality from a supplied state pairing

This module isolates the Cauchy--Schwarz step in the finite-dimensional
Robertson uncertainty relation.  A consumer supplies the two centred
state-action vectors in a complex inner-product space.  Their squared norms
are the supplied variances and twice the imaginary part of their pairing is
the supplied commutator readout.  Nothing here selects a state, observable,
Hamiltonian, sector, or physical instrument.

The abstract Cauchy--Schwarz lemma is followed by a concrete finite-state
adapter.  For a supplied density matrix and two supplied Hermitian
observables, the positive-semidefinite matrix seminorm provides the pairing;
the adapter centres both observables and proves the state-level Robertson
bound.  Nothing in the adapter selects those inputs.
-/

namespace EventAlgebra.Robertson

variable {V : Type*} [SeminormedAddCommGroup V] [InnerProductSpace ℂ V]

/-- The variance supplied by a centred state-action vector. -/
noncomputable def suppliedVariance (x : V) : ℝ :=
  (inner ℂ x x).re

/-- The real Robertson commutator readout, `2 Im ⟨x,y⟩`.  In a concrete
state representation the consumer must separately prove that this equals
`-i` times the expectation of the observable commutator. -/
noncomputable def suppliedCommutatorReadout (x y : V) : ℝ :=
  2 * (inner ℂ x y).im

/-- Cauchy--Schwarz gives the Robertson inequality for any supplied pair of
centred state-action vectors. -/
theorem supplied_pairing_robertson (x y : V) :
    suppliedCommutatorReadout x y ^ 2 / 4 ≤
      suppliedVariance x * suppliedVariance y := by
  have hIm :
      (inner ℂ x y).im * (inner ℂ x y).im ≤
        (inner ℂ x x).re * (inner ℂ y y).re := by
    calc
      (inner ℂ x y).im * (inner ℂ x y).im ≤
          Complex.normSq (inner ℂ x y) :=
        Complex.im_sq_le_normSq (inner ℂ x y)
      _ = ‖inner ℂ x y‖ * ‖inner ℂ y x‖ := by
        rw [Complex.normSq_eq_norm_sq, pow_two, norm_inner_symm]
      _ ≤ (inner ℂ x x).re * (inner ℂ y y).re :=
        inner_mul_inner_self_le x y
  unfold suppliedCommutatorReadout suppliedVariance
  nlinarith

/-- A zero first variance forces the supplied commutator readout to vanish. -/
theorem commutatorReadout_eq_zero_of_variance_left_eq_zero
    {x y : V} (h : suppliedVariance x = 0) :
    suppliedCommutatorReadout x y = 0 := by
  have hRobertson := supplied_pairing_robertson x y
  have hsq : 0 ≤ suppliedCommutatorReadout x y ^ 2 := sq_nonneg _
  rw [h, zero_mul] at hRobertson
  nlinarith

/-- A zero second variance forces the supplied commutator readout to vanish. -/
theorem commutatorReadout_eq_zero_of_variance_right_eq_zero
    {x y : V} (h : suppliedVariance y = 0) :
    suppliedCommutatorReadout x y = 0 := by
  have hRobertson := supplied_pairing_robertson x y
  have hsq : 0 ≤ suppliedCommutatorReadout x y ^ 2 := sq_nonneg _
  rw [h, mul_zero] at hRobertson
  nlinarith

/-! ## Exact skew-pairing control -/

/-- The complex-line control has unit norms and nonzero skew pairing.  It is
an algebraic saturation control, not a control involving observables. -/
theorem complex_line_control_values :
    suppliedVariance (1 : ℂ) = 1 ∧
      suppliedVariance Complex.I = 1 ∧
      suppliedCommutatorReadout (1 : ℂ) Complex.I = 2 := by
  norm_num [suppliedVariance, suppliedCommutatorReadout, inner]

/-- The complex-line control saturates the Robertson inequality. -/
theorem complex_line_control_saturates :
    suppliedCommutatorReadout (1 : ℂ) Complex.I ^ 2 / 4 =
      suppliedVariance (1 : ℂ) * suppliedVariance Complex.I := by
  norm_num [suppliedVariance, suppliedCommutatorReadout, inner]

/-! ## Adapter for a supplied finite matrix state -/

open Matrix
open scoped ComplexOrder

variable {n : ℕ}

/-- Centre a supplied observable by the real part of its expectation in a
supplied state.  The theorem `expectation_eq_ofReal_re_of_isHermitian` below
identifies that real part with the full expectation for Hermitian inputs. -/
noncomputable def centeredObservable
    (ρ A : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  A - ((EventAlgebra.expectation ρ A).re : ℂ) • 1

/-- The variance of a supplied observable, defined by the squared centred
observable in the supplied state. -/
noncomputable def stateVariance
    (ρ A : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (EventAlgebra.expectation ρ
    (centeredObservable ρ A * centeredObservable ρ A)).re

/-- Twice the imaginary part of the centred ordered-product expectation.
For Hermitian inputs, `neg_I_mul_commutator_expectation_eq_readout` below
identifies it with `-i` times the ordinary commutator expectation. -/
noncomputable def stateCommutatorReadout
    (ρ A B : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  2 * (EventAlgebra.expectation ρ
    (centeredObservable ρ A * centeredObservable ρ B)).im

/-- Centring preserves Hermiticity. -/
theorem centeredObservable_isHermitian
    {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hA : A.IsHermitian) : (centeredObservable ρ A).IsHermitian := by
  unfold centeredObservable
  apply hA.sub
  exact Matrix.isHermitian_one.smul (by simp [isSelfAdjoint_iff])

/-- A Hermitian state and Hermitian observable have a real expectation. -/
theorem expectation_eq_ofReal_re_of_isHermitian
    {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian) :
    (((EventAlgebra.expectation ρ A).re : ℝ) : ℂ) =
      EventAlgebra.expectation ρ A := by
  exact EventAlgebra.bornWeight_eq_re hρ.1.1 hA

/-- The supplied-state expectation of the centred observable is zero. -/
theorem expectation_centeredObservable_eq_zero
    {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian) :
    EventAlgebra.expectation ρ (centeredObservable ρ A) = 0 := by
  change EventAlgebra.stateExpectationLinearMap ρ
    (A - ((EventAlgebra.expectation ρ A).re : ℂ) • 1) = 0
  rw [map_sub, map_smul]
  change EventAlgebra.expectation ρ A -
    ((EventAlgebra.expectation ρ A).re : ℂ) •
      EventAlgebra.expectation ρ 1 = 0
  rw [EventAlgebra.expectation_one hρ]
  simp only [smul_eq_mul, mul_one]
  exact sub_eq_zero.mpr
    (expectation_eq_ofReal_re_of_isHermitian hρ hA).symm

/-- Centring does not change the observable commutator. -/
theorem centeredObservable_commutator
    (ρ A B : Matrix (Fin n) (Fin n) ℂ) :
    centeredObservable ρ A * centeredObservable ρ B -
        centeredObservable ρ B * centeredObservable ρ A =
      A * B - B * A := by
  simp only [centeredObservable, sub_mul, mul_sub, smul_mul_assoc,
    mul_smul_comm, one_mul, mul_one]
  module

/-- The finite-state Robertson theorem.  The density matrix, Hermitian
observables, and their physical interpretation are supplied inputs. -/
theorem finite_state_pairing_robertson
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) :
    stateCommutatorReadout ρ A B ^ 2 / 4 ≤
      stateVariance ρ A * stateVariance ρ B := by
  letI : SeminormedAddCommGroup (Matrix (Fin n) (Fin n) ℂ) :=
    Matrix.toMatrixSeminormedAddCommGroup ρ hρ.1
  letI : InnerProductSpace ℂ (Matrix (Fin n) (Fin n) ℂ) :=
    Matrix.toMatrixInnerProductSpace ρ hρ.1
  let A₀ := centeredObservable ρ A
  let B₀ := centeredObservable ρ B
  have hA₀ : A₀.IsHermitian := centeredObservable_isHermitian hA
  have hB₀ : B₀.IsHermitian := centeredObservable_isHermitian hB
  have hRobertson := supplied_pairing_robertson A₀ B₀
  have hAA : inner ℂ A₀ A₀ = EventAlgebra.expectation ρ (A₀ * A₀) := by
    change (A₀ * ρ * A₀ᴴ).trace = (ρ * (A₀ * A₀)).trace
    rw [hA₀.eq]
    calc
      (A₀ * ρ * A₀).trace = (A₀ * (A₀ * ρ)).trace :=
        trace_mul_comm (A₀ * ρ) A₀
      _ = ((A₀ * A₀) * ρ).trace := by rw [mul_assoc]
      _ = (ρ * (A₀ * A₀)).trace := trace_mul_comm _ _
  have hBB : inner ℂ B₀ B₀ = EventAlgebra.expectation ρ (B₀ * B₀) := by
    change (B₀ * ρ * B₀ᴴ).trace = (ρ * (B₀ * B₀)).trace
    rw [hB₀.eq]
    calc
      (B₀ * ρ * B₀).trace = (B₀ * (B₀ * ρ)).trace :=
        trace_mul_comm (B₀ * ρ) B₀
      _ = ((B₀ * B₀) * ρ).trace := by rw [mul_assoc]
      _ = (ρ * (B₀ * B₀)).trace := trace_mul_comm _ _
  have hAB : inner ℂ A₀ B₀ = EventAlgebra.expectation ρ (A₀ * B₀) := by
    change (B₀ * ρ * A₀ᴴ).trace = (ρ * (A₀ * B₀)).trace
    rw [hA₀.eq]
    calc
      (B₀ * ρ * A₀).trace = (A₀ * (B₀ * ρ)).trace :=
        trace_mul_comm (B₀ * ρ) A₀
      _ = ((A₀ * B₀) * ρ).trace := by rw [mul_assoc]
      _ = (ρ * (A₀ * B₀)).trace := trace_mul_comm _ _
  unfold suppliedCommutatorReadout suppliedVariance at hRobertson
  rw [hAB, hAA, hBB] at hRobertson
  simpa [stateCommutatorReadout, stateVariance, A₀, B₀] using hRobertson

/-- A state variance is nonnegative for every supplied Hermitian
observable. -/
theorem stateVariance_nonnegative
    {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian) :
    0 ≤ stateVariance ρ A := by
  let A₀ := centeredObservable ρ A
  have hA₀ : A₀.IsHermitian := centeredObservable_isHermitian hA
  have hPSD : (A₀ * A₀).PosSemidef := by
    simpa [hA₀.eq] using Matrix.posSemidef_conjTranspose_mul_self A₀
  exact (Complex.le_def.mp (EventAlgebra.expectation_nonneg hρ.1 hPSD)).1

/-- For Hermitian supplied observables, the real readout used by the bound
is the real part of `-i` times the expectation of their centred commutator. -/
theorem commutatorReadout_eq_neg_I_mul_expectation
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) :
    (((-Complex.I) * EventAlgebra.expectation ρ
      (centeredObservable ρ A * centeredObservable ρ B -
        centeredObservable ρ B * centeredObservable ρ A)).re) =
      stateCommutatorReadout ρ A B := by
  let A₀ := centeredObservable ρ A
  let B₀ := centeredObservable ρ B
  have hA₀ : A₀.IsHermitian := centeredObservable_isHermitian hA
  have hB₀ : B₀.IsHermitian := centeredObservable_isHermitian hB
  have hConj :
      EventAlgebra.expectation ρ (B₀ * A₀) =
        star (EventAlgebra.expectation ρ (A₀ * B₀)) := by
    change (ρ * (B₀ * A₀)).trace = star ((ρ * (A₀ * B₀)).trace)
    rw [← trace_conjTranspose]
    simp only [conjTranspose_mul, hρ.1.1.eq, hA₀.eq, hB₀.eq,
      mul_assoc]
    exact (trace_mul_cycle' B₀ A₀ ρ).symm
  change (((-Complex.I) * EventAlgebra.expectation ρ
    (A₀ * B₀ - B₀ * A₀)).re) =
    2 * (EventAlgebra.expectation ρ (A₀ * B₀)).im
  change (((-Complex.I) * EventAlgebra.bornWeight ρ
    (A₀ * B₀ - B₀ * A₀)).re) =
    2 * (EventAlgebra.bornWeight ρ (A₀ * B₀)).im
  have hConj' : EventAlgebra.bornWeight ρ (B₀ * A₀) =
      star (EventAlgebra.bornWeight ρ (A₀ * B₀)) := hConj
  rw [EventAlgebra.bornWeight_sub, hConj']
  simp
  ring

/-- The same real-part bridge written with the ordinary, uncentred observable
commutator. -/
theorem commutatorReadout_eq_neg_I_mul_uncentered_expectation
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) :
    (((-Complex.I) * EventAlgebra.expectation ρ (A * B - B * A)).re) =
      stateCommutatorReadout ρ A B := by
  rw [← centeredObservable_commutator ρ A B]
  exact commutatorReadout_eq_neg_I_mul_expectation hρ hA hB

/-- For Hermitian supplied observables, `-i` times the expectation of the
ordinary commutator is real and equals the real Robertson readout as a
complex number. -/
theorem neg_I_mul_commutator_expectation_eq_readout
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) :
    (-Complex.I) * EventAlgebra.expectation ρ (A * B - B * A) =
      (stateCommutatorReadout ρ A B : ℂ) := by
  let A₀ := centeredObservable ρ A
  let B₀ := centeredObservable ρ B
  have hA₀ : A₀.IsHermitian := centeredObservable_isHermitian hA
  have hB₀ : B₀.IsHermitian := centeredObservable_isHermitian hB
  have hConj :
      EventAlgebra.expectation ρ (B₀ * A₀) =
        star (EventAlgebra.expectation ρ (A₀ * B₀)) := by
    change (ρ * (B₀ * A₀)).trace = star ((ρ * (A₀ * B₀)).trace)
    rw [← trace_conjTranspose]
    simp only [conjTranspose_mul, hρ.1.1.eq, hA₀.eq, hB₀.eq,
      mul_assoc]
    exact (trace_mul_cycle' B₀ A₀ ρ).symm
  rw [← centeredObservable_commutator ρ A B]
  change (-Complex.I) * EventAlgebra.expectation ρ
      (A₀ * B₀ - B₀ * A₀) =
    ((2 * (EventAlgebra.expectation ρ (A₀ * B₀)).im : ℝ) : ℂ)
  change (-Complex.I) * EventAlgebra.bornWeight ρ
      (A₀ * B₀ - B₀ * A₀) =
    ((2 * (EventAlgebra.bornWeight ρ (A₀ * B₀)).im : ℝ) : ℂ)
  have hConj' : EventAlgebra.bornWeight ρ (B₀ * A₀) =
      star (EventAlgebra.bornWeight ρ (A₀ * B₀)) := hConj
  rw [EventAlgebra.bornWeight_sub, hConj']
  apply Complex.ext <;> simp [two_mul]

/-- Conventional finite-state Robertson inequality written directly with
the ordinary observable commutator. -/
theorem finite_state_robertson_commutator
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : EventAlgebra.IsState ρ) (hA : A.IsHermitian)
    (hB : B.IsHermitian) :
    Complex.normSq (EventAlgebra.expectation ρ (A * B - B * A)) / 4 ≤
      stateVariance ρ A * stateVariance ρ B := by
  have hbridge := neg_I_mul_commutator_expectation_eq_readout hρ hA hB
  have hnorm := congrArg Complex.normSq hbridge
  have hnegI : Complex.normSq (-Complex.I) = 1 := by
    norm_num [Complex.normSq]
  have hreal :
      Complex.normSq ((stateCommutatorReadout ρ A B : ℝ) : ℂ) =
        stateCommutatorReadout ρ A B ^ 2 := by
    simp [Complex.normSq, pow_two]
  rw [Complex.normSq.map_mul, hnegI, one_mul, hreal] at hnorm
  rw [hnorm]
  exact finite_state_pairing_robertson hρ hA hB

/-! ## Two-level matrix controls -/

/-- The pure `+z` two-level state used by the exact controls. -/
def pauliState : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1, 0; 0, 0]

/-- Pauli `X`. -/
def pauliX : Matrix (Fin 2) (Fin 2) ℂ :=
  !![0, 1; 1, 0]

/-- Pauli `Y`. -/
def pauliY : Matrix (Fin 2) (Fin 2) ℂ :=
  !![0, -Complex.I; Complex.I, 0]

/-- Pauli `Z`. -/
def pauliZ : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1, 0; 0, -1]

theorem pauliState_isEvent : EventAlgebra.IsEvent pauliState := by
  constructor
  · show pauliStateᴴ = pauliState
    ext i j
    fin_cases i <;> fin_cases j <;> simp [pauliState]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [pauliState, Matrix.mul_apply, Fin.sum_univ_two]

theorem pauliState_isState : EventAlgebra.IsState pauliState := by
  constructor
  · exact pauliState_isEvent.posSemidef
  · norm_num [pauliState, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

theorem pauliX_isHermitian : pauliX.IsHermitian := by
  show pauliXᴴ = pauliX
  ext i j
  fin_cases i <;> fin_cases j <;> simp [pauliX]

theorem pauliY_isHermitian : pauliY.IsHermitian := by
  show pauliYᴴ = pauliY
  ext i j
  fin_cases i <;> fin_cases j <;> simp [pauliY]

theorem pauliZ_isHermitian : pauliZ.IsHermitian := by
  show pauliZᴴ = pauliZ
  ext i j
  fin_cases i <;> fin_cases j <;> simp [pauliZ]

/-- The `X,Y` control observables do not commute. -/
theorem pauliX_pauliY_ne_pauliY_pauliX : pauliX * pauliY ≠ pauliY * pauliX := by
  intro h
  have h00 := congrFun (congrFun h 0) 0
  norm_num [pauliX, pauliY, Matrix.mul_apply, Fin.sum_univ_two] at h00
  have him := congrArg Complex.im h00
  norm_num at him

/-- The zero-variance `Z,X` control observables also do not commute. -/
theorem pauliZ_pauliX_ne_pauliX_pauliZ : pauliZ * pauliX ≠ pauliX * pauliZ := by
  intro h
  have h01 := congrFun (congrFun h 0) 1
  norm_num [pauliZ, pauliX, Matrix.mul_apply, Fin.sum_univ_two] at h01

/-- The Pauli `X,Y` control has unit variances and a nonzero commutator
readout, saturating the finite-state Robertson inequality. -/
theorem pauli_xy_noncommuting_control :
    stateVariance pauliState pauliX = 1 ∧
      stateVariance pauliState pauliY = 1 ∧
      stateCommutatorReadout pauliState pauliX pauliY = 2 := by
  norm_num [stateVariance, stateCommutatorReadout, centeredObservable,
    EventAlgebra.expectation, EventAlgebra.stateExpectationLinearMap,
    EventAlgebra.bornWeight, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    pauliState, pauliX, pauliY, Fin.sum_univ_two]

/-- In the same state, `Z` has zero variance and its commutator readout with
`X` vanishes. -/
theorem pauli_z_zero_variance_control :
    stateVariance pauliState pauliZ = 0 ∧
      stateCommutatorReadout pauliState pauliZ pauliX = 0 := by
  norm_num [stateVariance, stateCommutatorReadout, centeredObservable,
    EventAlgebra.expectation, EventAlgebra.stateExpectationLinearMap,
    EventAlgebra.bornWeight, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    pauliState, pauliX, pauliZ, Fin.sum_univ_two]

#print axioms EventAlgebra.Robertson.supplied_pairing_robertson
#print axioms EventAlgebra.Robertson.commutatorReadout_eq_zero_of_variance_left_eq_zero
#print axioms EventAlgebra.Robertson.commutatorReadout_eq_zero_of_variance_right_eq_zero
#print axioms EventAlgebra.Robertson.complex_line_control_values
#print axioms EventAlgebra.Robertson.complex_line_control_saturates
#print axioms EventAlgebra.Robertson.centeredObservable_isHermitian
#print axioms EventAlgebra.Robertson.expectation_eq_ofReal_re_of_isHermitian
#print axioms EventAlgebra.Robertson.expectation_centeredObservable_eq_zero
#print axioms EventAlgebra.Robertson.centeredObservable_commutator
#print axioms EventAlgebra.Robertson.finite_state_pairing_robertson
#print axioms EventAlgebra.Robertson.stateVariance_nonnegative
#print axioms EventAlgebra.Robertson.commutatorReadout_eq_neg_I_mul_expectation
#print axioms EventAlgebra.Robertson.commutatorReadout_eq_neg_I_mul_uncentered_expectation
#print axioms EventAlgebra.Robertson.neg_I_mul_commutator_expectation_eq_readout
#print axioms EventAlgebra.Robertson.finite_state_robertson_commutator
#print axioms EventAlgebra.Robertson.pauliState_isState
#print axioms EventAlgebra.Robertson.pauliX_pauliY_ne_pauliY_pauliX
#print axioms EventAlgebra.Robertson.pauliZ_pauliX_ne_pauliX_pauliZ
#print axioms EventAlgebra.Robertson.pauli_xy_noncommuting_control
#print axioms EventAlgebra.Robertson.pauli_z_zero_variance_control

end EventAlgebra.Robertson
