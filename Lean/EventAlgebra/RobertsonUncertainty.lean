import EventAlgebra.StateExpectation

/-!
# Robertson uncertainty for the finite state pairing

For a state `ρ` (positive-semidefinite, trace one) on the finite event
algebra `Matrix (Fin n) (Fin n) ℂ`, the trace pairing

  `⟪X, Y⟫_ρ = Tr(ρ Xᴴ Y)`

is a positive-semidefinite sesquilinear form: exactly a
`PreInnerProductSpace.Core` over `ℂ`, not an inner product, because a
mixed or rank-deficient state annihilates part of the algebra.
Cauchy–Schwarz for this pre-inner product yields the finite **Robertson
uncertainty relation**

  `‖Tr(ρ [A, B])‖² ≤ 4 · Var_ρ(A) · Var_ρ(B)`

for Hermitian observables `A`, `B`, where `Var_ρ(A)` is the state
pairing of the centred observable `A - Tr(ρ A) • 1` with itself.

Controls:

* **zero-variance control** — `Var_ρ(A) = 0` forces `Tr(ρ [A, B]) = 0`
  for every Hermitian `B`: a dispersion-free observable has no
  uncertainty obstruction against any partner;
* **noncommuting witness** — the qubit ground state with `σ_x`, `σ_y`
  attains the bound exactly: `‖Tr(ρ [σ_x, σ_y])‖² = 4` and
  `Var(σ_x) · Var(σ_y) = 1`, so both sides equal `4` and neither
  vanishes; the inequality is sharp and nonvacuous;
* **noncommutativity control** — `σ_x σ_y ≠ σ_y σ_x`, so the witness
  commutator is not trivially zero.

## Tagging convention

As in `EventAlgebra.Basic`. Every statement here is trace-dependent:
the pairing *is* the trace pairing of `EventAlgebra.StateExpectation`.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder
open ComplexConjugate

variable {n : ℕ}

/-! ## The state pairing -/

/-- The **state pairing** of two observables under `ρ`: the trace pairing
`Tr(ρ Xᴴ Y)`, conjugate-linear in the first argument and linear in the
second. For a state `ρ` this is a positive-semidefinite sesquilinear
form — the finite GNS pairing of `ρ`. -/
noncomputable def stateInner (ρ X Y : Matrix (Fin n) (Fin n) ℂ) : ℂ :=
  expectation ρ (Xᴴ * Y)

/-- **Trace-dependent.** The state pairing evaluated as a trace. -/
theorem stateInner_eq_trace (ρ X Y : Matrix (Fin n) (Fin n) ℂ) :
    stateInner ρ X Y = (ρ * (Xᴴ * Y)).trace :=
  rfl

/-- **Trace-dependent.** Additivity in the first argument. -/
theorem stateInner_add_left (ρ X Y Z : Matrix (Fin n) (Fin n) ℂ) :
    stateInner ρ (X + Y) Z = stateInner ρ X Z + stateInner ρ Y Z := by
  simp [stateInner_eq_trace, conjTranspose_add, add_mul, mul_add, trace_add]

/-- **Trace-dependent.** Conjugate homogeneity in the first argument. -/
theorem stateInner_smul_left (ρ X Y : Matrix (Fin n) (Fin n) ℂ) (r : ℂ) :
    stateInner ρ (r • X) Y = conj r * stateInner ρ X Y := by
  rw [starRingEnd_apply, stateInner_eq_trace, stateInner_eq_trace,
    conjTranspose_smul, smul_mul_assoc, mul_smul_comm, trace_smul,
    smul_eq_mul]

/-- **Trace-dependent.** Hermitian symmetry of the state pairing: for a
Hermitian `ρ`, swapping the arguments conjugates the value. -/
theorem stateInner_conj_symm {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.IsHermitian) (X Y : Matrix (Fin n) (Fin n) ℂ) :
    conj (stateInner ρ Y X) = stateInner ρ X Y := by
  rw [starRingEnd_apply, stateInner_eq_trace, stateInner_eq_trace,
    ← trace_conjTranspose, conjTranspose_mul, conjTranspose_mul,
    conjTranspose_conjTranspose, hρ.eq, trace_mul_comm]

/-- **Trace-dependent.** Positivity of the diagonal of the state pairing
in the partial order of `ℂ`: `Tr(ρ Xᴴ X) ≥ 0` for positive-semidefinite
`ρ`. This is the positivity of the finite GNS pairing. -/
theorem stateInner_self_nonneg {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (X : Matrix (Fin n) (Fin n) ℂ) :
    0 ≤ stateInner ρ X X :=
  expectation_nonneg hρ (posSemidef_conjTranspose_mul_self X)

/-- **Trace-dependent.** The state pairing of a positive-semidefinite `ρ`
bundled as a `PreInnerProductSpace.Core` over `ℂ`. This is exactly the
right structure: the pairing is positive semidefinite but not definite
(a rank-deficient state annihilates the complement of its support), so
no `InnerProductSpace.Core` exists in general. -/
@[implicit_reducible]
noncomputable def stateInnerCore {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) :
    PreInnerProductSpace.Core ℂ (Matrix (Fin n) (Fin n) ℂ) where
  inner := stateInner ρ
  conj_inner_symm X Y := stateInner_conj_symm hρ.1 X Y
  re_inner_nonneg X := (Complex.nonneg_iff.mp (stateInner_self_nonneg hρ X)).1
  add_left X Y Z := stateInner_add_left ρ X Y Z
  smul_left X Y r := stateInner_smul_left ρ X Y r

/-- **Trace-dependent.** Cauchy–Schwarz for the state pairing: the
squared norm of the off-diagonal pairing is bounded by the product of the
real diagonal pairings. Inherited from the semi-definite Cauchy–Schwarz
inequality of `PreInnerProductSpace.Core`. -/
theorem stateInner_cauchySchwarz {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (X Y : Matrix (Fin n) (Fin n) ℂ) :
    ‖stateInner ρ X Y‖ ^ 2 ≤ (stateInner ρ X X).re * (stateInner ρ Y Y).re := by
  have h := InnerProductSpace.Core.inner_mul_inner_self_le
    (𝕜 := ℂ) (c := stateInnerCore hρ) X Y
  have hsymm : ‖stateInner ρ Y X‖ = ‖stateInner ρ X Y‖ := by
    rw [← stateInner_conj_symm hρ.1 X Y, RCLike.norm_conj]
  have h' : ‖stateInner ρ X Y‖ * ‖stateInner ρ Y X‖ ≤
      (stateInner ρ X X).re * (stateInner ρ Y Y).re := h
  rw [hsymm] at h'
  calc ‖stateInner ρ X Y‖ ^ 2
      = ‖stateInner ρ X Y‖ * ‖stateInner ρ X Y‖ := by rw [pow_two]
    _ ≤ (stateInner ρ X X).re * (stateInner ρ Y Y).re := h'

/-! ## Centred observables and variance -/

/-- The **centred observable**: `A` minus its expectation times the
identity. -/
noncomputable def deviation (ρ A : Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ :=
  A - expectation ρ A • 1

/-- The **variance** of an observable under `ρ`: the real part of the
state pairing of the centred observable with itself. For a Hermitian
observable and a state this is the usual `Tr(ρ (A - ⟨A⟩)²)`. -/
noncomputable def variance (ρ A : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (stateInner ρ (deviation ρ A) (deviation ρ A)).re

/-- **Trace-dependent.** Variance is nonnegative. -/
theorem variance_nonneg {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (A : Matrix (Fin n) (Fin n) ℂ) :
    0 ≤ variance ρ A :=
  (Complex.nonneg_iff.mp (stateInner_self_nonneg hρ (deviation ρ A))).1

/-- **Trace-dependent.** The expectation of a Hermitian observable under
a Hermitian `ρ` is self-conjugate (real). -/
theorem expectation_conj {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.IsHermitian) (hA : A.IsHermitian) :
    conj (expectation ρ A) = expectation ρ A := by
  rw [starRingEnd_apply]
  show star ((ρ * A).trace) = (ρ * A).trace
  rw [← trace_conjTranspose, conjTranspose_mul, hρ.eq, hA.eq, trace_mul_comm]

/-- **Trace-dependent.** The centred observable of a Hermitian observable
under a Hermitian `ρ` is Hermitian. -/
theorem deviation_isHermitian {ρ A : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.IsHermitian) (hA : A.IsHermitian) :
    (deviation ρ A).IsHermitian := by
  show (A - expectation ρ A • 1)ᴴ = A - expectation ρ A • 1
  rw [conjTranspose_sub, conjTranspose_smul, conjTranspose_one, hA.eq,
    ← starRingEnd_apply, expectation_conj hρ hA]

/-- **Algebra-only.** Centring drops out of the commutator: the
commutator of the centred observables is the commutator of the
observables. -/
theorem deviation_commutator (ρ A B : Matrix (Fin n) (Fin n) ℂ) :
    deviation ρ A * deviation ρ B - deviation ρ B * deviation ρ A =
      A * B - B * A := by
  simp only [deviation, sub_mul, mul_sub, smul_mul_assoc, mul_smul_comm,
    one_mul, mul_one]
  module

/-! ## The Robertson uncertainty relation -/

/-- **Trace-dependent.** The **finite Robertson uncertainty relation**:
for a state `ρ` and Hermitian observables `A`, `B`,

  `‖Tr(ρ (AB - BA))‖² ≤ 4 · Var_ρ(A) · Var_ρ(B)`.

The proof is the Cauchy–Schwarz inequality for the state pairing applied
to the centred observables, plus the observation that the commutator
expectation is twice the imaginary part of their pairing. -/
theorem robertson_uncertainty {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    ‖expectation ρ (A * B - B * A)‖ ^ 2 ≤
      4 * (variance ρ A * variance ρ B) := by
  obtain ⟨hpsd, -⟩ := hρ
  have hΔA : (deviation ρ A).IsHermitian := deviation_isHermitian hpsd.1 hA
  have hΔB : (deviation ρ B).IsHermitian := deviation_isHermitian hpsd.1 hB
  -- The commutator expectation is `c - conj c` for the pairing `c` of the
  -- centred observables.
  have h1 : stateInner ρ (deviation ρ A) (deviation ρ B) =
      expectation ρ (deviation ρ A * deviation ρ B) := by
    show expectation ρ ((deviation ρ A)ᴴ * deviation ρ B) = _
    rw [hΔA.eq]
  have h3 : stateInner ρ (deviation ρ B) (deviation ρ A) =
      expectation ρ (deviation ρ B * deviation ρ A) := by
    show expectation ρ ((deviation ρ B)ᴴ * deviation ρ A) = _
    rw [hΔB.eq]
  have hkey : expectation ρ (A * B - B * A) =
      stateInner ρ (deviation ρ A) (deviation ρ B) -
        conj (stateInner ρ (deviation ρ A) (deviation ρ B)) := by
    rw [stateInner_conj_symm hpsd.1, h1, h3, ← deviation_commutator ρ A B]
    exact bornWeight_sub ρ _ _
  -- `‖c - conj c‖² = 4 (im c)²`.
  have him : ∀ z : ℂ, ‖z - conj z‖ ^ 2 = 4 * z.im ^ 2 := by
    intro z
    rw [Complex.sub_conj]
    rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real]
    rw [Real.norm_eq_abs, sq_abs]
    ring
  have hcs := stateInner_cauchySchwarz hpsd (deviation ρ A) (deviation ρ B)
  have him_le : (stateInner ρ (deviation ρ A) (deviation ρ B)).im ^ 2 ≤
      ‖stateInner ρ (deviation ρ A) (deviation ρ B)‖ ^ 2 := by
    have h := Complex.abs_im_le_norm
      (stateInner ρ (deviation ρ A) (deviation ρ B))
    have h0 := norm_nonneg (stateInner ρ (deviation ρ A) (deviation ρ B))
    nlinarith [abs_nonneg (stateInner ρ (deviation ρ A) (deviation ρ B)).im,
      sq_abs (stateInner ρ (deviation ρ A) (deviation ρ B)).im]
  have hvar : (stateInner ρ (deviation ρ A) (deviation ρ A)).re *
      (stateInner ρ (deviation ρ B) (deviation ρ B)).re =
      variance ρ A * variance ρ B := rfl
  calc ‖expectation ρ (A * B - B * A)‖ ^ 2
      = 4 * (stateInner ρ (deviation ρ A) (deviation ρ B)).im ^ 2 := by
        rw [hkey, him]
    _ ≤ 4 * ‖stateInner ρ (deviation ρ A) (deviation ρ B)‖ ^ 2 := by
        linarith
    _ ≤ 4 * (variance ρ A * variance ρ B) := by
        rw [← hvar]; linarith

/-- **Trace-dependent.** Robertson uncertainty in the conventional
normalization `(1/4)‖⟨[A,B]⟩‖² ≤ Var(A)·Var(B)`. -/
theorem robertson_uncertainty_div {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    ‖expectation ρ (A * B - B * A)‖ ^ 2 / 4 ≤
      variance ρ A * variance ρ B := by
  have h := robertson_uncertainty hρ hA hB
  linarith

/-- **Trace-dependent.** Zero-variance control: a dispersion-free
observable has vanishing commutator expectation against every partner.
In particular the Robertson bound is never violated trivially — a zero
variance kills the left-hand side too. -/
theorem commutator_expectation_eq_zero_of_variance_eq_zero
    {ρ A B : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hA : A.IsHermitian) (hB : B.IsHermitian)
    (h0 : variance ρ A = 0) :
    expectation ρ (A * B - B * A) = 0 := by
  have hr := robertson_uncertainty hρ hA hB
  rw [h0, zero_mul, mul_zero] at hr
  have hnorm : ‖expectation ρ (A * B - B * A)‖ = 0 := by
    nlinarith [norm_nonneg (expectation ρ (A * B - B * A))]
  exact norm_eq_zero.mp hnorm

/-! ## The qubit witness: the bound is sharp and nonvacuous

The ground state `|0⟩⟨0|` of a qubit with the Pauli observables `σ_x`,
`σ_y` attains the Robertson bound with both sides equal to `4`: the
commutator expectation is `2i`, of squared norm `4`, and both variances
are `1`. -/

/-- The qubit ground-state density matrix `|0⟩⟨0|`. -/
def qubitGround : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, 0]

/-- The Pauli `σ_x` observable. -/
def pauliX : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; 1, 0]

/-- The Pauli `σ_y` observable. -/
def pauliY : Matrix (Fin 2) (Fin 2) ℂ := !![0, -Complex.I; Complex.I, 0]

/-- The qubit ground state is a state. -/
theorem qubitGround_isState : IsState qubitGround := by
  constructor
  · have hd : qubitGround = diagonal ![1, 0] := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [qubitGround, diagonal]
    rw [hd]
    refine posSemidef_diagonal_iff.mpr ?_
    intro i
    fin_cases i
    · exact zero_le_one
    · exact le_refl 0
  · simp [qubitGround, trace_fin_two]

/-- `σ_x` is Hermitian. -/
theorem pauliX_isHermitian : pauliX.IsHermitian := by
  show pauliXᴴ = pauliX
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pauliX, conjTranspose_apply]

/-- `σ_y` is Hermitian. -/
theorem pauliY_isHermitian : pauliY.IsHermitian := by
  show pauliYᴴ = pauliY
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pauliY, conjTranspose_apply]

/-- `σ_x² = 1`. -/
theorem pauliX_mul_self : pauliX * pauliX = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pauliX, Matrix.mul_apply, Fin.sum_univ_two]

/-- `σ_y² = 1`. -/
theorem pauliY_mul_self : pauliY * pauliY = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pauliY, Matrix.mul_apply, Fin.sum_univ_two, Complex.I_mul_I]

/-- **Noncommutativity control.** `σ_x σ_y ≠ σ_y σ_x`: the witness
commutator is not trivially zero. -/
theorem pauliX_pauliY_not_commute : pauliX * pauliY ≠ pauliY * pauliX := by
  intro h
  have h00 := congrFun (congrFun h 0) 0
  simp [pauliX, pauliY, Matrix.mul_apply, Fin.sum_univ_two,
    Complex.ext_iff] at h00
  norm_num at h00

/-- **Trace-dependent.** `⟨σ_x⟩ = 0` in the ground state. -/
theorem expectation_qubitGround_pauliX :
    expectation qubitGround pauliX = 0 := by
  show (qubitGround * pauliX).trace = 0
  rw [trace_fin_two]
  simp [qubitGround, pauliX, Matrix.mul_apply, Fin.sum_univ_two]

/-- **Trace-dependent.** `⟨σ_y⟩ = 0` in the ground state. -/
theorem expectation_qubitGround_pauliY :
    expectation qubitGround pauliY = 0 := by
  show (qubitGround * pauliY).trace = 0
  rw [trace_fin_two]
  simp [qubitGround, pauliY, Matrix.mul_apply, Fin.sum_univ_two]

/-- **Trace-dependent.** `Var(σ_x) = 1` in the ground state. -/
theorem variance_qubitGround_pauliX : variance qubitGround pauliX = 1 := by
  have hd : deviation qubitGround pauliX = pauliX := by
    rw [deviation, expectation_qubitGround_pauliX, zero_smul, sub_zero]
  rw [variance, hd]
  have h : stateInner qubitGround pauliX pauliX = 1 := by
    show expectation qubitGround (pauliXᴴ * pauliX) = 1
    rw [pauliX_isHermitian.eq, pauliX_mul_self]
    exact expectation_one qubitGround_isState
  rw [h]
  simp

/-- **Trace-dependent.** `Var(σ_y) = 1` in the ground state. -/
theorem variance_qubitGround_pauliY : variance qubitGround pauliY = 1 := by
  have hd : deviation qubitGround pauliY = pauliY := by
    rw [deviation, expectation_qubitGround_pauliY, zero_smul, sub_zero]
  rw [variance, hd]
  have h : stateInner qubitGround pauliY pauliY = 1 := by
    show expectation qubitGround (pauliYᴴ * pauliY) = 1
    rw [pauliY_isHermitian.eq, pauliY_mul_self]
    exact expectation_one qubitGround_isState
  rw [h]
  simp

/-- **Trace-dependent.** The ground-state commutator expectation of
`σ_x`, `σ_y` is `2i`. -/
theorem commutator_expectation_qubit :
    expectation qubitGround (pauliX * pauliY - pauliY * pauliX) =
      2 * Complex.I := by
  show (qubitGround * (pauliX * pauliY - pauliY * pauliX)).trace =
    2 * Complex.I
  rw [trace_fin_two]
  simp [qubitGround, pauliX, pauliY, Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- **Trace-dependent.** Nonvacuity of the witness: the commutator
expectation is nonzero, so the Robertson bound says something. -/
theorem robertson_qubit_nonvacuous :
    expectation qubitGround (pauliX * pauliY - pauliY * pauliX) ≠ 0 := by
  rw [commutator_expectation_qubit]
  simp [Complex.I_ne_zero]

/-- **Trace-dependent.** The qubit witness attains the Robertson bound
exactly: both sides equal `4`. The inequality is sharp and neither side
collapses to zero. -/
theorem robertson_qubit_witness :
    ‖expectation qubitGround (pauliX * pauliY - pauliY * pauliX)‖ ^ 2 = 4 ∧
      variance qubitGround pauliX * variance qubitGround pauliY = 1 := by
  constructor
  · rw [commutator_expectation_qubit]
    rw [norm_mul, Complex.norm_I, mul_one]
    norm_num
  · rw [variance_qubitGround_pauliX, variance_qubitGround_pauliY]
    norm_num

-- Axiom audit: all statements must remain within the standard Mathlib basis.
#print axioms stateInner_cauchySchwarz
#print axioms robertson_uncertainty
#print axioms robertson_uncertainty_div
#print axioms commutator_expectation_eq_zero_of_variance_eq_zero
#print axioms deviation_commutator
#print axioms qubitGround_isState
#print axioms pauliX_pauliY_not_commute
#print axioms variance_qubitGround_pauliX
#print axioms variance_qubitGround_pauliY
#print axioms commutator_expectation_qubit
#print axioms robertson_qubit_nonvacuous
#print axioms robertson_qubit_witness

end EventAlgebra
