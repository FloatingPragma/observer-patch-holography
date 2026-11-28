import Mathlib.Data.Real.Sqrt
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
The four configurations are ordered by the three material hops from the
first to the second hub word. The face flip connects the two hub words.
After subtracting an energy, the hub diagonal is d and the rim diagonals are
e and f. All observations below are on this isolated, zero-flux component.
-/

namespace OPH.MaterialLoop

def Mode (t g d e f : ℝ) (u : Fin 4 → ℝ) : Prop :=
  d * u 0 - t * u 1 - g * u 3 = 0 ∧
  -t * u 0 + e * u 1 - t * u 2 = 0 ∧
  -t * u 1 + f * u 2 - t * u 3 = 0 ∧
  -g * u 0 - t * u 2 + d * u 3 = 0

private theorem modeFour (t g d e f a b c z : ℝ) :
    Mode t g d e f ![a,b,c,z] ↔
      d * a - t * b - g * z = 0 ∧ -t * a + e * b - t * c = 0 ∧
      -t * b + f * c - t * z = 0 ∧ -g * a - t * c + d * z = 0 := by
  rfl

/-- The two hub records determine the two rim amplitudes. -/
theorem rimFromHub {t g d e f : ℝ} {u : Fin 4 → ℝ}
    (hu : Mode t g d e f u) :
    t * u 1 = d * u 0 - g * u 3 ∧ t * u 2 = d * u 3 - g * u 0 := by
  exact ⟨by linarith [hu.1], by linarith [hu.2.2.2]⟩

theorem hubInjective {t g d e f : ℝ} {u v : Fin 4 → ℝ}
    (ht : t ≠ 0) (hu : Mode t g d e f u) (hv : Mode t g d e f v)
    (h0 : u 0 = v 0) (h3 : u 3 = v 3) : u = v := by
  have limitsU := rimFromHub hu
  have limitsV := rimFromHub hv
  have h1 : u 1 = v 1 := by
    apply (mul_left_cancel₀ ht)
    rw [limitsU.1, limitsV.1, h0, h3]
  have h2 : u 2 = v 2 := by
    apply (mul_left_cancel₀ ht)
    rw [limitsU.2, limitsV.2, h0, h3]
  funext i
  fin_cases i <;> assumption

theorem reducedEquations {t g d e f : ℝ} {u : Fin 4 → ℝ}
    (hu : Mode t g d e f u) :
    (e * d - t ^ 2 + t * g) * u 0 + (-e * g - t * d) * u 3 = 0 ∧
    (-f * g - t * d) * u 0 + (f * d - t ^ 2 + t * g) * u 3 = 0 := by
  rcases hu with ⟨a, b, c, d'⟩
  constructor
  · linear_combination t * b + e * a - t * d'
  · linear_combination t * c + f * d' - t * a

private theorem coefficientsZero {a b x y p q : ℝ}
    (h : x * q - y * p ≠ 0) (h1 : a * x + b * y = 0)
    (h2 : a * p + b * q = 0) : a = 0 ∧ b = 0 := by
  have hA : (x * q - y * p) * a = 0 := by
    linear_combination q * h1 - y * h2
  have hB : (x * q - y * p) * b = 0 := by
    linear_combination x * h2 - p * h1
  exact ⟨(mul_eq_zero.mp hA).resolve_left h, (mul_eq_zero.mp hB).resolve_left h⟩

/-- Two independent stationary modes require the sharp rate threshold.
    Hub injection makes the two-by-two hub determinant the independence test. -/
theorem stationaryThreshold {t g d e f : ℝ} {u v : Fin 4 → ℝ}
    (ht : 0 < t) (hg : 0 ≤ g)
    (hu : Mode t g d e f u) (hv : Mode t g d e f v)
    (h : u 0 * v 3 - u 3 * v 0 ≠ 0) :
    t ≤ g ∧ e = f ∧ d ^ 2 = g * (g - t) := by
  have equationsU := reducedEquations hu
  have equationsV := reducedEquations hv
  have ⟨hd, he⟩ := coefficientsZero h equationsU.1 equationsV.1
  have ⟨hf, _⟩ := coefficientsZero h equationsU.2 equationsV.2
  have prod : t * (d ^ 2 - g * (g - t)) = 0 := by
    linear_combination -g * hd - d * he
  have square := (mul_eq_zero.mp prod).resolve_left ht.ne'
  by_cases hzero : g = 0
  · subst g
    have dzero : d = 0 := by
      have prod' : t * d = 0 := by nlinarith [he]
      exact (mul_eq_zero.mp prod').resolve_left ht.ne'
    subst d
    nlinarith [hd, mul_pos ht ht]
  · have gpos : 0 < g := lt_of_le_of_ne hg (Ne.symm hzero)
    have hrate : t ≤ g := by
      by_contra! contrary
      have neg := mul_neg_of_pos_of_neg gpos (sub_neg.mpr contrary)
      nlinarith [sq_nonneg d, neg]
    have levels : e = f := by
      have prod' : (e - f) * g = 0 := by linear_combination hf - he
      exact sub_eq_zero.mp ((mul_eq_zero.mp prod').resolve_right hzero)
    exact ⟨hrate, levels, sub_eq_zero.mp square⟩

/-- Every rate at or above threshold admits a scalar choice with two modes. -/
theorem thresholdAttained {t g : ℝ} (ht : 0 < t) (hg : t ≤ g) :
    ∃ d e : ℝ, ∃ u v : Fin 4 → ℝ,
      Mode t g d e e u ∧ Mode t g d e e v ∧
      u 0 * v 3 - u 3 * v 0 = 1 := by
  have gpos := lt_of_lt_of_le ht hg
  let d := Real.sqrt (g * (g - t))
  have dsq : d ^ 2 = g * (g - t) :=
    Real.sq_sqrt (mul_nonneg gpos.le (sub_nonneg.mpr hg))
  let e := -t * d / g
  have he : e * g + t * d = 0 := by
    dsimp [e]
    field_simp
    ring
  have hd : e * d - t ^ 2 + t * g = 0 := by
    apply (mul_left_cancel₀ gpos.ne')
    linear_combination d * he - t * dsq
  refine ⟨d, e, ![1, d / t, -g / t, 0], ![0, -g / t, d / t, 1], ?_, ?_, ?_⟩
  · rw [modeFour]
    refine ⟨?_, ?_, ?_, ?_⟩
    · field_simp
      ring
    · field_simp
      linear_combination hd
    · field_simp
      linear_combination -he
    · field_simp
      ring
  · rw [modeFour]
    refine ⟨?_, ?_, ?_, ?_⟩
    · field_simp
      ring
    · field_simp
      linear_combination -he
    · field_simp
      linear_combination hd
    · field_simp
      ring
  · change (1 : ℝ) * 1 - 0 * 0 = 1
    norm_num

/-- The symmetric physical rim at 1045 has potentials (0,8,8).
    At t=4 its even reflection coordinates have generator (4,-4;-4,12). -/
def evenPolynomial (E : ℝ) : ℝ := (4 - E) * (12 - E) - 16

/-- The odd reflection coordinates have generator (4+2g,-4;-4,20). -/
def oddPolynomial (g E : ℝ) : ℝ := (4 + 2 * g - E) * (20 - E) - 16

theorem evenReflectionMode (g E x y : ℝ) :
    Mode 4 g (4 + g - E) (16 - E) (16 - E) ![x,y,y,x] ↔
      (4-E)*x - 4*y = 0 ∧ -4*x + (12-E)*y = 0 := by
  rw [modeFour]
  constructor <;> intro h
  · exact ⟨by linear_combination h.1, by linear_combination h.2.1⟩
  · exact ⟨by linear_combination h.1, by linear_combination h.2,
      by linear_combination h.2, by linear_combination h.1⟩

theorem oddReflectionMode (g E x y : ℝ) :
    Mode 4 g (4 + g - E) (16 - E) (16 - E) ![x,y,-y,-x] ↔
      (4+2*g-E)*x - 4*y = 0 ∧ -4*x + (20-E)*y = 0 := by
  rw [modeFour]
  constructor <;> intro h
  · exact ⟨by linear_combination h.1, by linear_combination h.2.1⟩
  · exact ⟨by linear_combination h.1, by linear_combination h.2,
      by linear_combination -h.2, by linear_combination -h.1⟩

def shiftedPhysical (g E : ℝ) : Matrix (Fin 4) (Fin 4) ℝ := fun i j =>
  match i.val, j.val with
  | 0, 0 | 3, 3 => 4 + g - E
  | 1, 1 | 2, 2 => 16 - E
  | 0, 1 | 1, 0 | 1, 2 | 2, 1 | 2, 3 | 3, 2 => -4
  | 0, 3 | 3, 0 => -g
  | _, _ => 0

/-- Even and odd parities factor the determinant of this same four-state matrix. -/
theorem physicalCharacteristic (g E : ℝ) :
    (shiftedPhysical g E).det = evenPolynomial E * oddPolynomial g E := by
  have r0 : Fin.succ (0 : Fin 3) = 1 := by decide
  have r1 : Fin.succ (1 : Fin 3) = 2 := by decide
  have r2 : Fin.succ (2 : Fin 3) = 3 := by decide
  have c10 : Fin.succAbove (1 : Fin 4) (0 : Fin 3) = 0 := by decide
  have c11 : Fin.succAbove (1 : Fin 4) (1 : Fin 3) = 2 := by decide
  have c12 : Fin.succAbove (1 : Fin 4) (2 : Fin 3) = 3 := by decide
  have c30 : Fin.succAbove (3 : Fin 4) (0 : Fin 3) = 0 := by decide
  have c31 : Fin.succAbove (3 : Fin 4) (1 : Fin 3) = 1 := by decide
  have c32 : Fin.succAbove (3 : Fin 4) (2 : Fin 3) = 2 := by decide
  rw [Matrix.det_succ_row_zero, Fin.sum_univ_four]
  simp only [Matrix.det_fin_three, Matrix.submatrix_apply, Fin.zero_succAbove,
    r0, r1, r2, c10, c11, c12, c30, c31, c32]
  norm_num [shiftedPhysical, evenPolynomial, oddPolynomial]
  ring

theorem nativeCharacteristic (E : ℝ) :
    (shiftedPhysical (1/2) E).det =
      (E^2 - 16*E + 32) * (E-4) * (E-21) := by
  rw [physicalCharacteristic]
  unfold evenPolynomial oddPolynomial
  ring

theorem tunedCrossing :
    evenPolynomial (8 + 4 * Real.sqrt 2) = 0 ∧
    oddPolynomial ((20 + 16 * Real.sqrt 2) / 7) (8 + 4 * Real.sqrt 2) = 0 := by
  have hs : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  constructor
  · unfold evenPolynomial
    nlinarith
  · unfold oddPolynomial
    nlinarith

/-- The physical crossing is the only reflection crossing at a nonnegative rate. -/
theorem crossingUnique {g E : ℝ} (hg : 0 ≤ g)
    (hEven : evenPolynomial E = 0) (hOdd : oddPolynomial g E = 0) :
    E = 8 + 4 * Real.sqrt 2 ∧ g = (20 + 16 * Real.sqrt 2) / 7 := by
  have hs := Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
  have s0 := Real.sqrt_nonneg (2 : ℝ)
  have s1 : 1 < Real.sqrt 2 := by nlinarith
  have s3 : Real.sqrt 2 < 3 := by nlinarith
  have factor : (E - (8 + 4 * Real.sqrt 2)) * (E - (8 - 4 * Real.sqrt 2)) = 0 := by
    unfold evenPolynomial at hEven
    nlinarith
  rcases mul_eq_zero.mp factor with high | low
  · have hE := sub_eq_zero.mp high
    refine ⟨hE, ?_⟩
    have hgt : 0 < 20 - E := by linarith
    have htuned := tunedCrossing.2
    rw [← hE] at htuned
    unfold oddPolynomial at hOdd htuned
    have final : (20 - E) * (g - (20 + 16 * Real.sqrt 2) / 7) = 0 := by
      linear_combination (hOdd - htuned) / 2
    exact sub_eq_zero.mp ((mul_eq_zero.mp final).resolve_left hgt.ne')
  · have hE := sub_eq_zero.mp low
    have a : 0 < 4 - E := by linarith
    have b : 0 < 20 - E := by linarith
    have bb := mul_nonneg hg b.le
    unfold evenPolynomial at hEven
    unfold oddPolynomial at hOdd
    nlinarith only [hEven, hOdd, a, bb]

#print axioms stationaryThreshold
#print axioms thresholdAttained
#print axioms nativeCharacteristic
#print axioms tunedCrossing
#print axioms crossingUnique

end OPH.MaterialLoop
