import Mathlib

/-!
# Einstein-branch tensor algebra

This module formalizes the finite linear-algebra and compositional part of the
corrected Einstein branch.  The continuum inputs are not asserted here:
timelike/null coverage, the Ward identity, and the contracted Bianchi identity
occur as theorem hypotheses.  In addition to the generic Minkowski lemmas, the
module gives an explicit nine-direction tomography frame in four dimensions
and a checked decoder for its trace-free coordinates.

There are no admissions, new axioms, or native-code proof shortcuts.
-/

namespace OPH.EinsteinBranch

noncomputable section

variable {n : ℕ}

/-- Vectors in Minkowski `ℝ^(1+n)`. -/
abbrev V (n : ℕ) := Fin (n + 1) → ℝ

/-- Bilinear forms, represented by their matrix entries. -/
abbrev Mat (n : ℕ) := Fin (n + 1) → Fin (n + 1) → ℝ

/-- Matrix evaluation as a bilinear form. -/
def bilinOf (B : Mat n) (x y : V n) : ℝ := ∑ i, ∑ j, B i j * x i * y j

/-- The associated quadratic form. -/
def quadOf (B : Mat n) (u : V n) : ℝ := bilinOf B u u

/-- Minkowski metric with signature `(-,+,…,+)`. -/
def eta (n : ℕ) : Mat n :=
  fun i j ↦ if i = j then (if i = 0 then -1 else 1) else 0

/-- A standard basis vector. -/
def e (i : Fin (n + 1)) : V n := Pi.single i 1

theorem bilinOf_add_left (B : Mat n) (x x' y : V n) :
    bilinOf B (x + x') y = bilinOf B x y + bilinOf B x' y := by
  unfold bilinOf
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [Pi.add_apply]
  ring

theorem bilinOf_add_right (B : Mat n) (x y y' : V n) :
    bilinOf B x (y + y') = bilinOf B x y + bilinOf B x y' := by
  unfold bilinOf
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [Pi.add_apply]
  ring

theorem bilinOf_smul_left (B : Mat n) (c : ℝ) (x y : V n) :
    bilinOf B (c • x) y = c * bilinOf B x y := by
  unfold bilinOf
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

theorem bilinOf_smul_right (B : Mat n) (c : ℝ) (x y : V n) :
    bilinOf B x (c • y) = c * bilinOf B x y := by
  unfold bilinOf
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

theorem bilinOf_single_single (B : Mat n) (p q : Fin (n + 1)) :
    bilinOf B (e p) (e q) = B p q := by
  unfold bilinOf e
  rw [Finset.sum_eq_single p]
  · rw [Finset.sum_eq_single q]
    · simp
    · intro l _ hl
      simp [Pi.single_eq_of_ne hl]
    · intro h
      exact absurd (Finset.mem_univ q) h
  · intro k _ hk
    apply Finset.sum_eq_zero
    intro l _
    simp [Pi.single_eq_of_ne hk]
  · intro h
    exact absurd (Finset.mem_univ p) h

theorem quadOf_sub_smul (A C : Mat n) (c : ℝ) (u : V n) :
    quadOf (fun i j ↦ A i j - c * C i j) u = quadOf A u - c * quadOf C u := by
  unfold quadOf bilinOf
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  ring

theorem quadOf_smul (B : Mat n) (c : ℝ) (x : V n) :
    quadOf B (c • x) = c ^ 2 * quadOf B x := by
  unfold quadOf
  rw [bilinOf_smul_left, bilinOf_smul_right]
  ring

theorem quadOf_double (B : Mat n) (hsymm : ∀ i j, B i j = B j i)
    (a b : ℝ) (i : Fin (n + 1)) :
    quadOf B (a • e 0 + b • e i) =
      a ^ 2 * B 0 0 + 2 * a * b * B 0 i + b ^ 2 * B i i := by
  unfold quadOf
  simp only [bilinOf_add_left, bilinOf_add_right, bilinOf_smul_left,
    bilinOf_smul_right, bilinOf_single_single]
  rw [hsymm i 0]
  ring

theorem quadOf_triple (B : Mat n) (hsymm : ∀ i j, B i j = B j i)
    (a b c : ℝ) (i j : Fin (n + 1)) :
    quadOf B (a • e 0 + b • e i + c • e j) =
      a ^ 2 * B 0 0 + b ^ 2 * B i i + c ^ 2 * B j j
        + 2 * a * b * B 0 i + 2 * a * c * B 0 j + 2 * b * c * B i j := by
  unfold quadOf
  simp only [bilinOf_add_left, bilinOf_add_right, bilinOf_smul_left,
    bilinOf_smul_right, bilinOf_single_single]
  rw [hsymm i 0, hsymm j 0, hsymm j i]
  ring

theorem eta_symm : ∀ i j : Fin (n + 1), eta n i j = eta n j i := by
  intro i j
  unfold eta
  by_cases h : i = j
  · subst h
    rfl
  · rw [if_neg h, if_neg (Ne.symm h)]

@[simp] theorem eta_zero_zero : eta n 0 0 = -1 := by simp [eta]

theorem eta_diag_spatial {i : Fin (n + 1)} (hi : i ≠ 0) : eta n i i = 1 := by
  simp [eta, hi]

theorem eta_off_diag {i j : Fin (n + 1)} (hij : i ≠ j) : eta n i j = 0 := by
  simp [eta, hij]

theorem e_zero_at_zero : (e (0 : Fin (n + 1))) 0 = 1 := Pi.single_eq_same 0 1

theorem double_at_zero (a b : ℝ) {i : Fin (n + 1)} (hi : i ≠ 0) :
    (a • e 0 + b • e i : V n) 0 = a := by
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, e]
  rw [Pi.single_eq_same, Pi.single_eq_of_ne (Ne.symm hi)]
  ring

theorem triple_at_zero (a b c : ℝ) {i j : Fin (n + 1)} (hi : i ≠ 0) (hj : j ≠ 0) :
    (a • e 0 + b • e i + c • e j : V n) 0 = a := by
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, e]
  rw [Pi.single_eq_same, Pi.single_eq_of_ne (Ne.symm hi), Pi.single_eq_of_ne (Ne.symm hj)]
  ring

/-! ## Timelike polarization and null ambiguity -/

/-- A symmetric form which vanishes on every future unit timelike vector is zero. -/
theorem unit_timelike_determines (B : Mat n) (hsymm : ∀ i j, B i j = B j i)
    (h : ∀ u : V n, quadOf (eta n) u = -1 → 0 < u 0 → quadOf B u = 0) :
    ∀ i j, B i j = 0 := by
  have cone : ∀ w : V n, quadOf (eta n) w < 0 → 0 < w 0 → quadOf B w = 0 := by
    intro w hw hw0
    set c : ℝ := Real.sqrt (-(quadOf (eta n) w)) with hc
    have hcpos : 0 < c := Real.sqrt_pos.mpr (by linarith)
    have hc2 : c ^ 2 = -(quadOf (eta n) w) := Real.sq_sqrt (by linarith)
    have hqu : quadOf (eta n) (c⁻¹ • w) = -1 := by
      rw [quadOf_smul]
      have hdiv : (c⁻¹) ^ 2 * quadOf (eta n) w = quadOf (eta n) w / c ^ 2 := by
        field_simp
      rw [hdiv, hc2, div_neg, div_self (ne_of_lt hw)]
    have hu0 : 0 < (c⁻¹ • w) 0 := by
      simp only [Pi.smul_apply, smul_eq_mul]
      positivity
    have hzero := h (c⁻¹ • w) hqu hu0
    have hw_eq : w = c • (c⁻¹ • w) := by
      rw [smul_smul, mul_inv_cancel₀ hcpos.ne', one_smul]
    calc
      quadOf B w = quadOf B (c • (c⁻¹ • w)) := by rw [← hw_eq]
      _ = c ^ 2 * quadOf B (c⁻¹ • w) := quadOf_smul B c _
      _ = 0 := by rw [hzero, mul_zero]
  have h00 : B 0 0 = 0 := by
    have he0 : quadOf B (e 0) = B 0 0 := bilinOf_single_single B 0 0
    have heta : quadOf (eta n) (e (0 : Fin (n + 1))) = -1 := by
      have hsingle := bilinOf_single_single (eta n) (0 : Fin (n + 1)) 0
      unfold quadOf
      rw [hsingle, eta_zero_zero]
    rw [← he0]
    exact h (e 0) heta (by rw [e_zero_at_zero]; norm_num)
  have hspace : ∀ i : Fin (n + 1), i ≠ 0 → B 0 i = 0 ∧ B i i = 0 := by
    intro i hi
    have hval : ∀ s : ℝ, s ^ 2 < 1 →
        B 0 0 + 2 * s * B 0 i + s ^ 2 * B i i = 0 := by
      intro s hs
      have heta : quadOf (eta n) ((1 : ℝ) • e 0 + s • e i) = -1 + s ^ 2 := by
        rw [quadOf_double (eta n) eta_symm]
        rw [eta_zero_zero, eta_diag_spatial hi, eta_off_diag (Ne.symm hi)]
        ring
      have hcone := cone ((1 : ℝ) • e 0 + s • e i)
        (by rw [heta]; linarith) (by rw [double_at_zero 1 s hi]; norm_num)
      rw [quadOf_double B hsymm] at hcone
      linarith
    have h12 := hval (1 / 2) (by norm_num)
    have h13 := hval (1 / 3) (by norm_num)
    constructor <;> linarith
  have hmixed : ∀ i j : Fin (n + 1), i ≠ 0 → j ≠ 0 → i ≠ j → B i j = 0 := by
    intro i j hi hj hij
    have heta : quadOf (eta n)
        ((1 : ℝ) • e 0 + (1 / 2 : ℝ) • e i + (1 / 2 : ℝ) • e j) = -(1 / 2) := by
      rw [quadOf_triple (eta n) eta_symm]
      rw [eta_zero_zero, eta_diag_spatial hi, eta_diag_spatial hj,
        eta_off_diag (Ne.symm hi), eta_off_diag (Ne.symm hj), eta_off_diag hij]
      ring
    have hcone := cone
      ((1 : ℝ) • e 0 + (1 / 2 : ℝ) • e i + (1 / 2 : ℝ) • e j)
      (by rw [heta]; norm_num)
      (by rw [triple_at_zero 1 (1 / 2) (1 / 2) hi hj]; norm_num)
    rw [quadOf_triple B hsymm] at hcone
    have hii := (hspace i hi).2
    have hjj := (hspace j hj).2
    have h0i := (hspace i hi).1
    have h0j := (hspace j hj).1
    linarith
  intro i j
  by_cases hi : i = 0
  · by_cases hj : j = 0
    · rw [hi, hj]
      exact h00
    · rw [hi]
      exact (hspace j hj).1
  · by_cases hj : j = 0
    · rw [hj, hsymm i 0]
      exact (hspace i hi).1
    · by_cases hij : i = j
      · rw [hij]
        exact (hspace j hj).2
      · exact hmixed i j hi hj hij

/-- A symmetric form vanishing on the null cone is a multiple of `eta`. -/
theorem null_cone_determines (B : Mat n) (hsymm : ∀ i j, B i j = B j i)
    (h : ∀ k : V n, quadOf (eta n) k = 0 → quadOf B k = 0) :
    ∀ i j, B i j = -(B 0 0) * eta n i j := by
  have hpm : ∀ i : Fin (n + 1), i ≠ 0 → B 0 i = 0 ∧ B i i = -(B 0 0) := by
    intro i hi
    have hval : ∀ s : ℝ, s ^ 2 = 1 →
        B 0 0 + 2 * s * B 0 i + s ^ 2 * B i i = 0 := by
      intro s hs
      have heta : quadOf (eta n) ((1 : ℝ) • e 0 + s • e i) = 0 := by
        rw [quadOf_double (eta n) eta_symm]
        rw [eta_zero_zero, eta_diag_spatial hi, eta_off_diag (Ne.symm hi)]
        linarith
      have hcone := h ((1 : ℝ) • e 0 + s • e i) heta
      rw [quadOf_double B hsymm] at hcone
      linarith
    have hp := hval 1 (by norm_num)
    have hm := hval (-1) (by norm_num)
    constructor <;> linarith
  have hmixed : ∀ i j : Fin (n + 1), i ≠ 0 → j ≠ 0 → i ≠ j → B i j = 0 := by
    intro i j hi hj hij
    set a : ℝ := Real.sqrt 2 with ha
    have ha2 : a ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    have heta : quadOf (eta n) (a • e 0 + (1 : ℝ) • e i + (1 : ℝ) • e j) = 0 := by
      rw [quadOf_triple (eta n) eta_symm]
      rw [eta_zero_zero, eta_diag_spatial hi, eta_diag_spatial hj,
        eta_off_diag (Ne.symm hi), eta_off_diag (Ne.symm hj), eta_off_diag hij]
      linarith
    have hcone := h (a • e 0 + (1 : ℝ) • e i + (1 : ℝ) • e j) heta
    rw [quadOf_triple B hsymm] at hcone
    have h0i := (hpm i hi).1
    have h0j := (hpm j hj).1
    have hii := (hpm i hi).2
    have hjj := (hpm j hj).2
    rw [h0i, h0j, hii, hjj, ha2] at hcone
    linarith
  intro i j
  by_cases hi : i = 0
  · subst hi
    by_cases hj : j = 0
    · subst hj
      rw [eta_zero_zero]
      ring
    · rw [eta_off_diag (Ne.symm hj), (hpm j hj).1]
      ring
  · by_cases hj : j = 0
    · subst hj
      rw [eta_off_diag hi, hsymm i 0, (hpm i hi).1]
      ring
    · by_cases hij : i = j
      · subst hij
        rw [eta_diag_spatial hi, (hpm i hi).2]
        ring
      · rw [eta_off_diag hij, hmixed i j hi hj hij]
        ring

/-- Null matching leaves exactly a metric term. -/
theorem jacobson_step (F T : Mat n)
    (hF : ∀ i j, F i j = F j i) (hT : ∀ i j, T i j = T j i) (κ : ℝ)
    (h : ∀ k : V n, quadOf (eta n) k = 0 → quadOf F k = κ * quadOf T k) :
    ∃ lam : ℝ, ∀ i j, F i j = κ * T i j + lam * eta n i j := by
  let B : Mat n := fun i j ↦ F i j - κ * T i j
  have hBsymm : ∀ i j, B i j = B j i := by
    intro i j
    dsimp [B]
    rw [hF i j, hT i j]
  have hBnull : ∀ k : V n, quadOf (eta n) k = 0 → quadOf B k = 0 := by
    intro k hk
    dsimp [B]
    rw [quadOf_sub_smul, h k hk, sub_self]
  have hmetric := null_cone_determines B hBsymm hBnull
  refine ⟨-(B 0 0), fun i j ↦ ?_⟩
  have hij := hmetric i j
  dsimp [B] at hij
  linarith

/-- Timelike scalar matching upgrades to entrywise tensor matching. -/
theorem tensor_upgrade (Gm gm T : Mat n) (Λ κ : ℝ)
    (hsymm : ∀ i j, Gm i j + Λ * gm i j - κ * T i j =
      Gm j i + Λ * gm j i - κ * T j i)
    (hall : ∀ u : V n, quadOf (eta n) u = -1 → 0 < u 0 →
      quadOf (fun i j ↦ Gm i j + Λ * gm i j) u = κ * quadOf T u) :
    ∀ i j, Gm i j + Λ * gm i j = κ * T i j := by
  let Y : Mat n := fun i j ↦ Gm i j + Λ * gm i j - κ * T i j
  have hYsymm : ∀ i j, Y i j = Y j i := by simpa [Y] using hsymm
  have hYvan : ∀ u : V n, quadOf (eta n) u = -1 → 0 < u 0 → quadOf Y u = 0 := by
    intro u h1 h2
    dsimp [Y]
    rw [quadOf_sub_smul, hall u h1 h2, sub_self]
  have hzero := unit_timelike_determines Y hYsymm hYvan
  intro i j
  have hij := hzero i j
  dsimp [Y] at hij
  linarith

/-! ## Explicit nine-direction null tomography in four dimensions -/

/-- Coordinates
`(T₀₀,T₀₁,T₀₂,T₀₃,T₁₁,T₁₂,T₁₃,T₂₂,T₂₃)`
for an `eta`-trace-free symmetric form; `T₃₃=T₀₀-T₁₁-T₂₂`. -/
abbrev TraceFreeCoords := Fin 9 → ℝ

/-- The trace-free symmetric matrix represented by nine coordinates. -/
def traceFreeMat (x : TraceFreeCoords) : Mat 3 :=
  ![![x 0, x 1, x 2, x 3],
    ![x 1, x 4, x 5, x 6],
    ![x 2, x 5, x 7, x 8],
    ![x 3, x 6, x 8, x 0 - x 4 - x 7]]

theorem traceFreeMat_symm (x : TraceFreeCoords) :
    ∀ i j, traceFreeMat x i j = traceFreeMat x j i := by
  intro i j
  fin_cases i <;> fin_cases j <;> rfl

/-- The normalized diagonal spatial coordinate `1/√3`, written as `√3/3`. -/
def invSqrtThree : ℝ := Real.sqrt 3 / 3

theorem invSqrtThree_sq : invSqrtThree ^ 2 = (1 / 3 : ℝ) := by
  rw [invSqrtThree, div_pow, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 3)]
  norm_num

/-- The paper's explicit frame: six axial directions and three body diagonals.
Every vector has the form `k=(1,n)` with `n∈S²`. -/
def tomographyDirections : Fin 9 → V 3 :=
  ![![1, 1, 0, 0],
    ![1, -1, 0, 0],
    ![1, 0, 1, 0],
    ![1, 0, -1, 0],
    ![1, 0, 0, 1],
    ![1, 0, 0, -1],
    ![1, invSqrtThree, invSqrtThree, invSqrtThree],
    ![1, invSqrtThree, invSqrtThree, -invSqrtThree],
    ![1, invSqrtThree, -invSqrtThree, invSqrtThree]]

/-- Every direction in the explicit tomography frame is Minkowski-null. -/
theorem tomographyDirections_null (r : Fin 9) :
    quadOf (eta 3) (tomographyDirections r) = 0 := by
  have ha := invSqrtThree_sq
  fin_cases r <;>
    simp [quadOf, bilinOf, eta, tomographyDirections, Fin.sum_univ_succ]
  all_goals ring_nf
  all_goals rw [ha]
  all_goals norm_num

/-- Evaluation of a trace-free form on the nine null directions, expanded as
the paper's nine-row design map. -/
def tomographyDesign (x : TraceFreeCoords) : Fin 9 → ℝ :=
  ![x 0 + 2 * x 1 + x 4,
    x 0 - 2 * x 1 + x 4,
    x 0 + 2 * x 2 + x 7,
    x 0 - 2 * x 2 + x 7,
    2 * x 0 + 2 * x 3 - x 4 - x 7,
    2 * x 0 - 2 * x 3 - x 4 - x 7,
    (4 / 3 : ℝ) * x 0 + 2 * invSqrtThree * (x 1 + x 2 + x 3) +
      (2 / 3 : ℝ) * (x 5 + x 6 + x 8),
    (4 / 3 : ℝ) * x 0 + 2 * invSqrtThree * (x 1 + x 2 - x 3) +
      (2 / 3 : ℝ) * (x 5 - x 6 - x 8),
    (4 / 3 : ℝ) * x 0 + 2 * invSqrtThree * (x 1 - x 2 + x 3) +
      (2 / 3 : ℝ) * (-x 5 + x 6 - x 8)]

/-- The actual null-null charge map. -/
def tomographyCharge (x : TraceFreeCoords) (r : Fin 9) : ℝ :=
  quadOf (traceFreeMat x) (tomographyDirections r)

theorem tomographyCharge_eq_design (x : TraceFreeCoords) :
    tomographyCharge x = tomographyDesign x := by
  have ha := invSqrtThree_sq
  funext r
  fin_cases r <;>
    simp [tomographyCharge, tomographyDesign, tomographyDirections, traceFreeMat,
      quadOf, bilinOf, Fin.sum_univ_succ]
  all_goals ring_nf
  all_goals rw [ha]
  all_goals ring

/-- An explicit inverse of the nine-row tomography design. -/
def tomographyDecoder (q : Fin 9 → ℝ) : TraceFreeCoords :=
  let s := Real.sqrt 3
  ![(q 0 + q 1 + q 2 + q 3 + q 4 + q 5) / 8,
    (q 0 - q 1) / 4,
    (q 2 - q 3) / 4,
    (q 4 - q 5) / 4,
    (3 * q 0 + 3 * q 1 - q 2 - q 3 - q 4 - q 5) / 8,
    (-(s + 1) * q 0 + (s - 1) * q 1 - (s + 1) * q 2 + (s - 1) * q 3 -
      q 4 - q 5 + 3 * q 6 + 3 * q 7) / 4,
    (-(s + 1) * q 0 + (s - 1) * q 1 - q 2 - q 3 - (s + 1) * q 4 +
      (s - 1) * q 5 + 3 * q 6 + 3 * q 8) / 4,
    (-q 0 - q 1 + 3 * q 2 + 3 * q 3 - q 4 - q 5) / 8,
    ((1 + s) * q 0 + (1 - s) * q 1 + q 2 + q 3 + q 4 + q 5 -
      3 * q 7 - 3 * q 8) / 4]

/-- The decoder is a left inverse of the explicit design. -/
theorem tomographyDecoder_design (x : TraceFreeCoords) :
    tomographyDecoder (tomographyDesign x) = x := by
  funext i
  have hs : (Real.sqrt 3) ^ 2 = (3 : ℝ) := Real.sq_sqrt (by norm_num)
  fin_cases i <;>
    simp [tomographyDecoder, tomographyDesign, invSqrtThree] <;>
    nlinarith

/-- The decoder reconstructs the trace-free tensor from its nine actual
null-null charges. -/
theorem tomographyDecoder_charge (x : TraceFreeCoords) :
    tomographyDecoder (tomographyCharge x) = x := by
  rw [tomographyCharge_eq_design]
  exact tomographyDecoder_design x

@[simp] theorem tomographyDecoder_zero :
    tomographyDecoder (0 : Fin 9 → ℝ) = 0 := by
  funext i
  fin_cases i <;> simp [tomographyDecoder]

/-- Nine explicit directions determine all nine trace-free components. -/
theorem tomographyCharge_injective : Function.Injective tomographyCharge := by
  intro x y h
  rw [← tomographyDecoder_charge x, ← tomographyDecoder_charge y, h]

/-- The coefficient of the metric part of a symmetric four-tensor. -/
def metricCoefficient (B : Mat 3) : ℝ :=
  (-B 0 0 + B 1 1 + B 2 2 + B 3 3) / 4

/-- Trace-free coordinates obtained after subtracting the metric component. -/
def traceFreeCoordsOf (B : Mat 3) : TraceFreeCoords :=
  let lam := metricCoefficient B
  ![B 0 0 + lam, B 0 1, B 0 2, B 0 3,
    B 1 1 - lam, B 1 2, B 1 3, B 2 2 - lam, B 2 3]

theorem traceFree_decomposition (B : Mat 3) (hsymm : ∀ i j, B i j = B j i) :
    ∀ i j, B i j = traceFreeMat (traceFreeCoordsOf B) i j +
      metricCoefficient B * eta 3 i j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [traceFreeMat, traceFreeCoordsOf, metricCoefficient, eta]
  all_goals first | ring | rw [hsymm]

/-- Vanishing on the nine displayed null directions already forces a
symmetric form to be a metric multiple. -/
theorem nine_null_directions_determine_mod_metric (B : Mat 3)
    (hsymm : ∀ i j, B i j = B j i)
    (hzero : ∀ r, quadOf B (tomographyDirections r) = 0) :
    ∀ i j, B i j = metricCoefficient B * eta 3 i j := by
  let x := traceFreeCoordsOf B
  have hmat : traceFreeMat x = fun i j ↦ B i j - metricCoefficient B * eta 3 i j := by
    funext i j
    rw [traceFree_decomposition B hsymm i j]
    ring
  have hcharge : tomographyCharge x = 0 := by
    funext r
    rw [tomographyCharge, hmat, quadOf_sub_smul, hzero r,
      tomographyDirections_null r]
    simp
  have hx : x = 0 := by
    rw [← tomographyDecoder_charge x, hcharge, tomographyDecoder_zero]
  have hx' : traceFreeCoordsOf B = 0 := by simpa [x] using hx
  intro i j
  rw [traceFree_decomposition B hsymm i j, hx']
  fin_cases i <;> fin_cases j <;> norm_num [traceFreeMat]

/-- Equality of the nine charges determines two symmetric tensors up to one
multiple of the Minkowski metric. -/
theorem nine_charge_metric_ambiguity (F T : Mat 3)
    (hF : ∀ i j, F i j = F j i) (hT : ∀ i j, T i j = T j i)
    (hcharge : ∀ r,
      quadOf F (tomographyDirections r) = quadOf T (tomographyDirections r)) :
    ∃ lam : ℝ, ∀ i j, F i j = T i j + lam * eta 3 i j := by
  let B : Mat 3 := fun i j ↦ F i j - 1 * T i j
  have hBsymm : ∀ i j, B i j = B j i := by
    intro i j
    dsimp [B]
    rw [hF i j, hT i j]
  have hBzero : ∀ r, quadOf B (tomographyDirections r) = 0 := by
    intro r
    dsimp [B]
    rw [quadOf_sub_smul, hcharge r]
    ring
  refine ⟨metricCoefficient B, fun i j ↦ ?_⟩
  have hij := nine_null_directions_determine_mod_metric B hBsymm hBzero i j
  dsimp [B] at hij
  linarith

/-! ## Ward/Bianchi composition and componentwise constancy -/

/-- Forward coordinate difference of a scalar field. -/
def dpar {P : Type*} (step : Fin (n + 1) → P → P) (f : P → ℝ)
    (i : Fin (n + 1)) (p : P) : ℝ :=
  f (step i p) - f p

/-- Discrete divergence, contracting forward differences in the first index. -/
def ddiv {P : Type*} (step : Fin (n + 1) → P → P) (M : P → Mat n)
    (p : P) (j : Fin (n + 1)) : ℝ :=
  ∑ i, (M (step i p) i j - M p i j)

/-- Reachability by forward coordinate steps. -/
def Reachable {P : Type*} (step : Fin (n + 1) → P → P) : P → P → Prop :=
  Relation.ReflTransGen (fun p q ↦ ∃ i, q = step i p)

/-- Connectivity generated by forward or backward coordinate steps. -/
def SymmReachable {P : Type*} (step : Fin (n + 1) → P → P) : P → P → Prop :=
  Relation.ReflTransGen
    (fun p q ↦ (∃ i, q = step i p) ∨ (∃ i, p = step i q))

/-- The divergence of `lam • eta` is the eta-contracted discrete gradient. -/
theorem ddiv_lam_eta {P : Type*} (step : Fin (n + 1) → P → P) (lam : P → ℝ)
    (p : P) (j : Fin (n + 1)) :
    ddiv step (fun q ↦ fun i j' ↦ lam q * eta n i j') p j =
      ∑ i, dpar step lam i p * eta n i j := by
  unfold ddiv dpar
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- Nondegeneracy of the diagonal Minkowski metric, in the exact form used by
the divergence argument. -/
theorem row_eta_cancel (r : Fin (n + 1) → ℝ)
    (h : ∀ j, ∑ i, r i * eta n i j = 0) : ∀ i, r i = 0 := by
  intro j
  have hsum := h j
  have hsingle : ∑ i, r i * eta n i j = r j * eta n j j := by
    apply Finset.sum_eq_single j
    · intro b _ hb
      rw [eta_off_diag hb, mul_zero]
    · intro habs
      exact absurd (Finset.mem_univ j) habs
  rw [hsingle] at hsum
  by_cases hj : j = 0
  · subst hj
    rw [eta_zero_zero] at hsum
    linarith
  · rw [eta_diag_spatial hj] at hsum
    linarith

/-- If the geometric and stress sides are respectively Bianchi- and
Ward-conserved, the pointwise metric coefficient is invariant under every
coordinate step.  The two conservation statements are explicit hypotheses. -/
theorem step_invariant_of_divergence_free {P : Type*}
    (step : Fin (n + 1) → P → P) (F T : P → Mat n) (lam : P → ℝ) (κ : ℝ)
    (hpoint : ∀ p, ∀ i j, F p i j = κ * T p i j + lam p * eta n i j)
    (hdivF : ∀ p j, ddiv step F p j = 0)
    (hdivT : ∀ p j, ddiv step T p j = 0) :
    ∀ (i : Fin (n + 1)) (p : P), lam (step i p) = lam p := by
  intro i p
  have htrace : ∀ j, ∑ i', dpar step lam i' p * eta n i' j = 0 := by
    intro j
    have hlin : ∑ i', dpar step lam i' p * eta n i' j =
        ddiv step F p j - κ * ddiv step T p j := by
      rw [← ddiv_lam_eta step lam p j]
      unfold ddiv
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro i' _
      rw [hpoint (step i' p) i' j, hpoint p i' j]
      ring
    rw [hlin, hdivF, hdivT, mul_zero, sub_zero]
  have hz := row_eta_cancel (fun i' ↦ dpar step lam i' p) htrace i
  unfold dpar at hz
  linarith

/-- A step-invariant scalar is constant along symmetric reachability. -/
theorem lam_eq_of_symmReachable {P : Type*}
    (step : Fin (n + 1) → P → P) (lam : P → ℝ)
    (hstep : ∀ (i : Fin (n + 1)) (p : P), lam (step i p) = lam p)
    {p q : P} (h : SymmReachable step p q) : lam q = lam p := by
  induction h with
  | refl => rfl
  | @tail b c _ hbc ih =>
      rcases hbc with ⟨i, rfl⟩ | ⟨i, hb⟩
      · rw [hstep, ih]
      · have h1 : lam b = lam c := by rw [hb, hstep]
        rw [← h1]
        exact ih

/-- Ward plus Bianchi conservation makes the metric term one constant on a
connected discrete component. -/
theorem lambda_constant_symm {P : Type*}
    (step : Fin (n + 1) → P → P) (p₀ : P)
    (hconn : ∀ q : P, SymmReachable step p₀ q)
    (F T : P → Mat n) (lam : P → ℝ) (κ : ℝ)
    (hpoint : ∀ p, ∀ i j, F p i j = κ * T p i j + lam p * eta n i j)
    (hdivF : ∀ p j, ddiv step F p j = 0)
    (hdivT : ∀ p j, ddiv step T p j = 0) :
    ∃ Λ : ℝ, ∀ p, lam p = Λ := by
  refine ⟨lam p₀, fun p ↦ ?_⟩
  exact lam_eq_of_symmReachable step lam
    (step_invariant_of_divergence_free step F T lam κ hpoint hdivF hdivT)
    (hconn p)

/-- Pointwise null matching, symmetry, explicit Ward/Bianchi premises, and
component connectivity compose to an Einstein-form equation with one constant
metric coefficient. -/
theorem einstein_equation_with_constant_symm {P : Type*}
    (step : Fin (n + 1) → P → P) (p₀ : P)
    (hconn : ∀ q : P, SymmReachable step p₀ q)
    (F T : P → Mat n) (κ : ℝ)
    (hFsymm : ∀ p, ∀ i j, F p i j = F p j i)
    (hTsymm : ∀ p, ∀ i j, T p i j = T p j i)
    (hnull : ∀ p, ∀ k : V n, quadOf (eta n) k = 0 →
      quadOf (F p) k = κ * quadOf (T p) k)
    (hdivF : ∀ p j, ddiv step F p j = 0)
    (hdivT : ∀ p j, ddiv step T p j = 0) :
    ∃ Λ : ℝ, ∀ p, ∀ i j, F p i j = κ * T p i j + Λ * eta n i j := by
  choose lam hlam using fun p ↦
    jacobson_step (F p) (T p) (hFsymm p) (hTsymm p) κ (hnull p)
  obtain ⟨Λ, hΛ⟩ := lambda_constant_symm step p₀ hconn F T lam κ
    (fun p i j ↦ hlam p i j) hdivF hdivT
  exact ⟨Λ, fun p i j ↦ by rw [hlam p i j, hΛ p]⟩

/-! ## Per-theorem axiom audit -/

#print axioms bilinOf_add_left
#print axioms bilinOf_add_right
#print axioms bilinOf_smul_left
#print axioms bilinOf_smul_right
#print axioms bilinOf_single_single
#print axioms quadOf_sub_smul
#print axioms quadOf_smul
#print axioms quadOf_double
#print axioms quadOf_triple
#print axioms eta_symm
#print axioms eta_zero_zero
#print axioms eta_diag_spatial
#print axioms eta_off_diag
#print axioms e_zero_at_zero
#print axioms double_at_zero
#print axioms triple_at_zero
#print axioms unit_timelike_determines
#print axioms null_cone_determines
#print axioms jacobson_step
#print axioms tensor_upgrade
#print axioms traceFreeMat_symm
#print axioms invSqrtThree_sq
#print axioms tomographyDirections_null
#print axioms tomographyCharge_eq_design
#print axioms tomographyDecoder_design
#print axioms tomographyDecoder_charge
#print axioms tomographyDecoder_zero
#print axioms tomographyCharge_injective
#print axioms traceFree_decomposition
#print axioms nine_null_directions_determine_mod_metric
#print axioms nine_charge_metric_ambiguity
#print axioms ddiv_lam_eta
#print axioms row_eta_cancel
#print axioms step_invariant_of_divergence_free
#print axioms lam_eq_of_symmReachable
#print axioms lambda_constant_symm
#print axioms einstein_equation_with_constant_symm

end

end OPH.EinsteinBranch
