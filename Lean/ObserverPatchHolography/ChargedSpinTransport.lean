import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Analysis.Complex.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Record comparison for charged spin transport

A half beat translates the two admitted energy components, including their
correlations with an idle reference. The observed claim compares the entire
receiver alphabet, which may contain rejected preparations and measurement
refusals. The Bell test uses raw probabilities for a common admitted event:
the two observers compare three signed Pauli readings with the admission
frequency. Neither implication identifies the underlying Hamiltonian, the
physical charges, or the apparatus implementing these readings.
-/

open scoped BigOperators ComplexConjugate

namespace OPH.ChargedSpinTransport

section Flight

variable {V : Type*} [AddCommGroup V] [Module ℂ V]

/-- Opposite translation characters on the two accepted components turn the
free half beat into the same translation, up to one common phase. The
ambient module can include an arbitrary retained reference. -/
theorem halfBeat_translation (flight shift : V →ₗ[ℂ] V) (phase : ℂ) (even odd : V)
    (flightEven : flight even = phase • even) (flightOdd : flight odd = -(phase • odd))
    (shiftEven : shift even = even) (shiftOdd : shift odd = -odd) :
    flight (even + odd) = phase • shift (even + odd) := by
  simp [map_add, flightEven, flightOdd, shiftEven, shiftOdd, smul_add]

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℂ W]

/-- The complete raw receiver law is transported under the exact half beat.
The outcome set can include every terminal refusal and may retain a reference.
The detector is the same linear amplitude map in the two comparisons. -/
theorem halfBeat_all_records {Record : Type*} (receiver : Record → V →ₗ[ℂ] W)
    (flight shift : V →ₗ[ℂ] V) (phase : ℂ) (even odd : V)
    (unitPhase : ‖phase‖ = 1)
    (flightEven : flight even = phase • even) (flightOdd : flight odd = -(phase • odd))
    (shiftEven : shift even = even) (shiftOdd : shift odd = -odd) (r : Record) :
    ‖receiver r (flight (even + odd))‖ ^ 2 = ‖receiver r (shift (even + odd))‖ ^ 2 := by
  rw [halfBeat_translation flight shift phase even odd flightEven flightOdd shiftEven shiftOdd,
    map_smul, norm_smul, unitPhase, one_mul]

/-- Summing the raw law over an admitted alphabet preserves the denominator
as well as every outcome. This identity precedes any conditioning. -/
theorem halfBeat_admission_weight {Record : Type*} (admitted : Finset Record)
    (receiver : Record → V →ₗ[ℂ] W) (flight shift : V →ₗ[ℂ] V)
    (phase : ℂ) (even odd : V) (unitPhase : ‖phase‖ = 1)
    (flightEven : flight even = phase • even) (flightOdd : flight odd = -(phase • odd))
    (shiftEven : shift even = even) (shiftOdd : shift odd = -odd) :
    ∑ r ∈ admitted, ‖receiver r (flight (even + odd))‖ ^ 2 =
      ∑ r ∈ admitted, ‖receiver r (shift (even + odd))‖ ^ 2 := by
  apply Finset.sum_congr rfl
  intro r _
  exact halfBeat_all_records receiver flight shift phase even odd unitPhase
    flightEven flightOdd shiftEven shiftOdd r

/-- A norm allowance on an admitted receiver amplitude propagates to the
record probability before dividing by the acceptance probability. -/
theorem record_error_of_amplitude {K : Type*} [NormedAddCommGroup K] (a b : K) {ε : ℝ}
    (aNorm : ‖a‖ ≤ 1) (bNorm : ‖b‖ ≤ 1) (error : ‖a - b‖ ≤ ε) :
    |‖a‖ ^ 2 - ‖b‖ ^ 2| ≤ 2 * ε := by
  have hdiff : |‖a‖ - ‖b‖| ≤ ε := by
    apply abs_le.mpr
    constructor
    · have h := norm_sub_norm_le b a
      rw [norm_sub_rev] at h
      linarith
    · exact (norm_sub_norm_le a b).trans error
  have hsum : 0 ≤ ‖a‖ + ‖b‖ := add_nonneg (norm_nonneg _) (norm_nonneg _)
  rw [show ‖a‖ ^ 2 - ‖b‖ ^ 2 = (‖a‖ - ‖b‖) * (‖a‖ + ‖b‖) by ring,
    abs_mul, abs_of_nonneg hsum, mul_comm 2]
  exact mul_le_mul hdiff (show ‖a‖ + ‖b‖ ≤ 2 by linarith) hsum
    ((abs_nonneg _).trans hdiff)

/-- Phase or approximation errors on the two components add in amplitude.
Contractive record maps propagate them to one uniform raw probability bound. -/
theorem halfBeat_record_error {A : Type*} [NormedAddCommGroup A] [NormedSpace ℂ A]
    (flight shift : A →ₗ[ℂ] A) (receiver : A →ₗ[ℂ] W) (phase : ℂ) (even odd : A)
    {εEven εOdd : ℝ} (unitPhase : ‖phase‖ = 1)
    (evenError : ‖flight even - phase • shift even‖ ≤ εEven)
    (oddError : ‖flight odd - phase • shift odd‖ ≤ εOdd)
    (contractive : ∀ v, ‖receiver v‖ ≤ ‖v‖)
    (actualNorm : ‖receiver (flight (even + odd))‖ ≤ 1)
    (idealNorm : ‖receiver (shift (even + odd))‖ ≤ 1) :
    |‖receiver (flight (even + odd))‖ ^ 2 - ‖receiver (shift (even + odd))‖ ^ 2| ≤
      2 * (εEven + εOdd) := by
  have rotatedNorm :
      ‖receiver (phase • shift (even + odd))‖ = ‖receiver (shift (even + odd))‖ := by
    rw [map_smul, norm_smul, unitPhase, one_mul]
  rw [← rotatedNorm]
  apply record_error_of_amplitude _ _ actualNorm (rotatedNorm.trans_le idealNorm)
  rw [← map_sub]
  refine (contractive _).trans ?_
  have hsum : flight (even + odd) - phase • shift (even + odd) =
      (flight even - phase • shift even) + (flight odd - phase • shift odd) := by
    simp only [map_add, smul_add]
    abel
  rw [hsum]
  exact (norm_add_le _ _).trans (add_le_add evenError oddError)

end Flight

section Bell

/-- Total raw probability for four joint spin outcomes of a pure component. -/
def spinAdmission (a b c d : ℂ) : ℝ :=
  Complex.normSq a + Complex.normSq b + Complex.normSq c + Complex.normSq d

/-- The unnormalized reference/receiver Bell probability for Φ⁺. -/
noncomputable def bellNumerator (a d : ℂ) : ℝ := Complex.normSq (a + d) / 2

/-- The three signed Pauli-pair numerators keep the original-start denominator. -/
def xxNumerator (a b c d : ℂ) : ℝ := 2 * ((conj a * d + conj b * c).re)
def yyNumerator (a b c d : ℂ) : ℝ := 2 * ((-(conj a * d) + conj b * c).re)
def zzNumerator (a b c d : ℂ) : ℝ :=
  Complex.normSq a - Complex.normSq b - Complex.normSq c + Complex.normSq d

/-- A common admitted spin event converts three signed Pauli records to the
Bell score, without an isotropy assumption on the received spin. -/
theorem pauli_bell_identity (a b c d : ℂ) :
    4 * bellNumerator a d =
      spinAdmission a b c d + xxNumerator a b c d - yyNumerator a b c d +
        zzNumerator a b c d := by
  simp only [bellNumerator, spinAdmission, xxNumerator, yyNumerator, zzNumerator,
    Complex.normSq_apply, Complex.add_re, Complex.add_im, Complex.mul_re,
    Complex.conj_re, Complex.conj_im, Complex.neg_re]
  ring

/-- One unnormalized product component cannot put more than half of its
admitted probability on the reference/receiver Bell report. -/
theorem product_bell_bound (a b c d : ℂ) :
    2 * bellNumerator (a * c) (b * d) ≤
      spinAdmission (a * c) (a * d) (b * c) (b * d) := by
  have h :
      2 * bellNumerator (a * c) (b * d) + Complex.normSq (a * conj d - b * conj c) =
        spinAdmission (a * c) (a * d) (b * c) (b * d) := by
    simp only [bellNumerator, spinAdmission, Complex.normSq_apply, Complex.add_re,
      Complex.add_im, Complex.mul_re, Complex.mul_im, Complex.conj_re, Complex.conj_im,
      Complex.sub_re, Complex.sub_im]
    ring
  linarith [Complex.normSq_nonneg (a * conj d - b * conj c)]

/-- The same upper bound holds for every finite separable mixture, with
arbitrary nonnegative weights and no normalization of the admitted event. -/
theorem separable_bell_bound {I : Type*} (terms : Finset I) (weight : I → ℝ)
    (nonnegative : ∀ i ∈ terms, 0 ≤ weight i) (a b c d : I → ℂ) :
    2 * (∑ i ∈ terms, weight i * bellNumerator (a i * c i) (b i * d i)) ≤
      ∑ i ∈ terms, weight i *
        spinAdmission (a i * c i) (a i * d i) (b i * c i) (b i * d i) := by
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro i hi
  rw [mul_left_comm]
  exact mul_le_mul_of_nonneg_left (product_bell_bound (a i) (b i) (c i) (d i))
    (nonnegative i hi)

/-- A finite separable description bounds the signed Pauli sum by the
joint admission frequency, in raw original-start units. -/
theorem separable_pauli_bound {I : Type*} (terms : Finset I) (weight : I → ℝ)
    (nonnegative : ∀ i ∈ terms, 0 ≤ weight i) (a b c d : I → ℂ) :
    ∑ i ∈ terms, weight i *
      (xxNumerator (a i * c i) (a i * d i) (b i * c i) (b i * d i) -
        yyNumerator (a i * c i) (a i * d i) (b i * c i) (b i * d i) +
        zzNumerator (a i * c i) (a i * d i) (b i * c i) (b i * d i)) ≤
    ∑ i ∈ terms, weight i *
      spinAdmission (a i * c i) (a i * d i) (b i * c i) (b i * d i) := by
  apply Finset.sum_le_sum
  intro i hi
  apply mul_le_mul_of_nonneg_left _ (nonnegative i hi)
  have product := product_bell_bound (a i) (b i) (c i) (d i)
  have read := pauli_bell_identity (a i * c i) (a i * d i) (b i * c i) (b i * d i)
  linarith

/-- Two observers can reject a finite separable description using raw
count intervals. The three Pauli reports must use the same admission event.
All errors may include finite clock, count and calibration uncertainty. -/
theorem pauli_count_interval_excludes_separable
    {p x y z bell pUpper xLower yUpper zLower : ℝ}
    (admission : p ≤ pUpper) (readX : xLower ≤ x) (readY : y ≤ yUpper)
    (readZ : zLower ≤ z) (identity : 4 * bell = p + x - y + z)
    (margin : pUpper < xLower - yUpper + zLower) :
    p < 2 * bell := by
  linarith

end Bell

end OPH.ChargedSpinTransport
