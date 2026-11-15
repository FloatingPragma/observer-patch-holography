import SourceLaw.TransferStructure

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Reflection positivity of the stationary chain, one mode at a time

The stationary chain of a mode has covariance `C(k) = c0 r^|k|`: `chain_covariance_row` and
`contact_identity_mode` show that this Toeplitz kernel inverts, row by row, the interior precision
of the register mode (slab block `2/Δτ + Δτλ`, time-neighbour block `−1/Δτ`); the finite open
register differs at its two end slabs. For the chain, `Δτ > 0` and `λ > 0`. Reflecting time
across a bond pairs a time `t ≥ 0` on one side with
`−t − 1` on the other, so the bond-reflection Gram form of coefficients `v` at times `t` is
`Σ_{a,b} v_a v_b c0 r^{t_a + t_b + 1}`; site reflection gives exponents `t_a + t_b`.

* `bondGram_eq_sq`, `siteGram_eq_sq`: both forms are `c0 r^ε (Σ_a v_a r^{t_a})²`.
* `bondGram_nonneg_iff`: bond reflection positivity of a mode holds exactly when `r ≥ 0`.
* `chain_bond_rp`: the stationary chain whose interior precision is that of the register mode
  satisfies it (the finite open register itself is not claimed reflection positive).
* `bondGram_neg_witness`: a negative transfer eigenvalue violates it.

At the level of Gaussian measures, reflection positivity of the measure is equivalent to that of
its covariance (standard; recorded in the companion note).
-/

namespace OPH.SourceLaw

open scoped BigOperators

/-- Bond-reflection Gram form of a mode with covariance `c0 r^|k|`. -/
def bondGram {n : ℕ} (c0 r : ℝ) (t : Fin n → ℕ) (v : Fin n → ℝ) : ℝ :=
  ∑ a, ∑ b, v a * v b * (c0 * r ^ (t a + t b + 1))

/-- Site-reflection Gram form of a mode with covariance `c0 r^|k|`. -/
def siteGram {n : ℕ} (c0 r : ℝ) (t : Fin n → ℕ) (v : Fin n → ℝ) : ℝ :=
  ∑ a, ∑ b, v a * v b * (c0 * r ^ (t a + t b))

theorem bondGram_eq_sq {n : ℕ} (c0 r : ℝ) (t : Fin n → ℕ) (v : Fin n → ℝ) :
    bondGram c0 r t v = c0 * r * (∑ a, v a * r ^ t a) ^ 2 := by
  rw [bondGram, sq, Finset.sum_mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [pow_succ, pow_add]
  ring

theorem siteGram_eq_sq {n : ℕ} (c0 r : ℝ) (t : Fin n → ℕ) (v : Fin n → ℝ) :
    siteGram c0 r t v = c0 * (∑ a, v a * r ^ t a) ^ 2 := by
  rw [siteGram, sq, Finset.sum_mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [pow_add]
  ring

theorem siteGram_nonneg {n : ℕ} (c0 r : ℝ) (t : Fin n → ℕ) (v : Fin n → ℝ) (hc0 : 0 ≤ c0) :
    0 ≤ siteGram c0 r t v := by
  rw [siteGram_eq_sq]
  positivity

/-- Bond reflection positivity of a mode holds exactly when its transfer eigenvalue is
nonnegative. -/
theorem bondGram_nonneg_iff (c0 r : ℝ) (hc0 : 0 < c0) :
    (∀ (n : ℕ) (t : Fin n → ℕ) (v : Fin n → ℝ), 0 ≤ bondGram c0 r t v) ↔ 0 ≤ r := by
  constructor
  · intro h
    have h1 := h 1 (fun _ => 0) (fun _ => 1)
    rw [bondGram_eq_sq] at h1
    have e : (∑ a : Fin 1, (fun _ => (1 : ℝ)) a * r ^ (fun _ => 0) a) = 1 := by simp
    rw [e, one_pow, mul_one] at h1
    exact (mul_nonneg_iff_of_pos_left hc0).1 h1
  · intro hr n t v
    rw [bondGram_eq_sq]
    positivity

/-- A negative transfer eigenvalue violates bond reflection positivity. -/
theorem bondGram_neg_witness (c0 r : ℝ) (hc0 : 0 < c0) (hr : r < 0) :
    bondGram c0 r (n := 1) (fun _ => 0) (fun _ => 1) < 0 := by
  rw [bondGram_eq_sq]
  have e : (∑ a : Fin 1, (fun _ => (1 : ℝ)) a * r ^ (fun _ => 0) a) = 1 := by simp
  rw [e, one_pow, mul_one]
  exact mul_neg_of_pos_of_neg hc0 hr

section Chain

variable {dτ lam : ℝ}

/-- For `Δτ > 0` and `λ > 0`, the stationary chain whose interior precision is that of the
register mode satisfies bond reflection positivity. -/
theorem chain_bond_rp (hd : 0 < dτ) (hl : 0 < lam) {n : ℕ} (t : Fin n → ℕ)
    (v : Fin n → ℝ) : 0 ≤ bondGram (equalTimeVar dτ lam) (transferEig dτ lam) t v :=
  ((bondGram_nonneg_iff _ _ (equalTimeVar_pos hd hl)).2 (transferEig_pos hd hl).le) n t v

/-- The Toeplitz kernel `c0 r^|k|` inverts the interior precision of the register mode away from
the diagonal: for `m ≥ 1`, `(2/Δτ + Δτλ) c0 r^m − (1/Δτ)(c0 r^(m−1) + c0 r^(m+1)) = 0`. Together with
`contact_identity_mode` (the diagonal row) this identifies `c0 r^|k|` as the covariance of the
stationary bi-infinite chain whose interior precision is that of the register mode (the finite
open register differs from it at its two end slabs). -/
theorem chain_covariance_row (hd : 0 < dτ) (hl : 0 < lam) (m : ℕ) :
    (2 / dτ + dτ * lam) * (equalTimeVar dτ lam * transferEig dτ lam ^ (m + 1))
      - 1 / dτ * (equalTimeVar dτ lam * transferEig dτ lam ^ m
        + equalTimeVar dτ lam * transferEig dτ lam ^ (m + 2)) = 0 := by
  have hr0 : transferEig dτ lam ≠ 0 := (transferEig_pos hd hl).ne'
  have hd0 : dτ ≠ 0 := hd.ne'
  have hsum := transferEig_add_inv hd hl
  have key : (2 + dτ ^ 2 * lam) * transferEig dτ lam = 1 + transferEig dτ lam ^ 2 := by
    rw [← hsum]
    field_simp
    try ring
  have e : 2 / dτ + dτ * lam = (2 + dτ ^ 2 * lam) / dτ := by
    field_simp
    try ring
  rw [e]
  linear_combination (equalTimeVar dτ lam * transferEig dτ lam ^ m / dτ) * key

end Chain

#print axioms bondGram_eq_sq
#print axioms siteGram_nonneg
#print axioms bondGram_nonneg_iff
#print axioms bondGram_neg_witness
#print axioms chain_bond_rp
#print axioms chain_covariance_row

end OPH.SourceLaw
