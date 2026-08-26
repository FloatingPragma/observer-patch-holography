import CarrierEvolutionFlow

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.FieldSectorEnergyInnerProduct

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierModeOscillators
open OPH.CarrierEvolutionFlow

/-!
# The energy inner product of the field sector inside the strict Courant window

STATUS.  Exact finite real and complex linear algebra on `Fin 2` and on
`ι → Fin 2 → ℝ`, elementary trigonometry, and the committed mode-form
`modeForm`, step `stepMatrix`, conjugator, flow `rotFlow`/`modeFlow`/
`assembledFlow` and energy identity `fieldEnergyScaled_scalarHistory` of
`CarrierEvolutionFlow`.  The step `h` is a declared Courant number and the
real parameter `t` is a declared flow parameter; physical time enters only
through the declared source-clock row, which this module does not touch.
No Born rule, no measurement, no field-to-observation attachment, no
register row discharged.

WHAT IS PROVED.
1. Positive definiteness is exactly the strict window.  For `h ≠ 0` the
   per-mode energy form `modeForm h lam (a, b) = lam a² + h lam a b + b²`
   is positive definite on `ℝ²` if and only if `0 < h² lam < 4`
   (`modeForm_posDef_iff`).  The failures are explicit: for `lam ≤ 0` the
   state `(1, 0)` has form `lam ≤ 0` (`modeForm_nonpos_of_lam_nonpos`); on
   the boundary `h² lam = 4` the nonzero state `(-h/2, 1)` has form `0`
   and the form is the rank-one square `lam (a + h b / 2)²`
   (`modeForm_boundary_zero`, `modeForm_boundary_rank_one`); outside,
   `h² lam > 4`, the same state has strictly negative form
   (`modeForm_exterior_neg`).  This is the same window `h² Λ < 4` as the
   committed `stability_certificate`.
2. The energy inner product.  `energyInner h lam x y` is the polarization
   of `modeForm` (`energyInner_polarization`, `energyInner_self`); it is
   symmetric and bilinear, positive definite inside the window, and it is
   packaged as a Mathlib `InnerProductSpace.Core ℝ (Fin 2 → ℝ)`
   (`energyInnerCore`).  The step is an isometry of it for every `h, lam`
   (`stepMatrix_energy_isometry`, a polynomial identity), and the per-mode
   flow is a one-parameter group of isometries (`modeFlow_energy_isometry`,
   `rotFlow_energy_isometry`, `modeFlow_isometry_group`).
3. Uniqueness of the invariant form.  For `h ≠ 0` every symmetric bilinear
   form on `ℝ²` invariant under `stepMatrix h lam` is the real multiple
   `B (0,1) (0,1) · energyInner h lam` (`invariant_symmetric_form_eq_smul`),
   and conversely every real multiple is symmetric and invariant
   (`invariant_symmetric_forms_iff`).  Inside the window a positive
   definite invariant form is a positive multiple
   (`invariant_inner_product_unique_up_to_positive_scale`).  The multiple
   is fixed by the committed staggered form: the symmetric bilinear forms
   whose diagonal is the committed staggered energy of the mode history
   are exactly `(‖v‖² / 2) · energyInner` (`staggered_energy_fixes_scale`,
   `staggered_energy_inner_product_iff`).  Polarization alone recovers a
   symmetric form from its diagonal, so the forward direction of that
   equivalence uses symmetry and the diagonal only and invariance is a
   consequence; the role of invariance is that the polarized energy is the
   only invariant candidate up to scale.  Uniqueness needs
   only `h ≠ 0` (the step is then not scalar); `sin θ ≠ 0` is not needed.
4. Complex structure and unitarity.  `complexStructure h lam` is the
   explicit matrix `(1 / sin θ) !![h² lam / 2, h; -h lam, -h² lam / 2]`,
   equal to the conjugated quarter turn `Q · rot (π/2) · Q⁻¹`
   (`complexStructure_eq_conj`), with `J² = -1` (`complexStructure_sq`),
   `energyInner`-orthogonal (`complexStructure_energy_isometry`) and
   `energyInner`-skew (`complexStructure_skew`).  The flow is
   `rotFlow h lam t = cos (θ t / h) · 1 + sin (θ t / h) · J`
   (`rotFlow_eq_cos_add_sin`).  The coordinate
   `modeCoordinate h lam x = (Q⁻¹ x)₀ + (Q⁻¹ x)₁ i ∈ ℂ` is an `ℝ`-linear
   bijection `ℝ² → ℂ` (`modeCoordinate_bijective`) under which `J` is
   multiplication by `i` (`modeCoordinate_complexStructure`) and the flow
   is multiplication by the unit phase `exp (i θ t / h)`
   (`modeCoordinate_rotFlow`).  The sesquilinear form
   `modeHermitian x y = energyInner x y - i energyInner x (J y)` equals
   `sin² θ · conj (z x) · z y` (`modeHermitian_eq_coordinate`), is
   conjugate-symmetric, `i`-linear in the second slot, has diagonal
   `modeForm` and is preserved by the flow (`modeHermitian_rotFlow`).  The
   orbit of every coordinate is differentiable with derivative
   `i (θ / h) · z(t)` (`modeCoordinate_flow_hasDerivAt`): the generator is
   the real number `θ / h` on this one-dimensional complex line, with
   sign convention recorded below.  Stone's theorem is not invoked; the
   generator is exhibited (`mode_unitary_group_explicit_generator`).
5. Assembled sector.  On a finite family of pairwise orthogonal eigenvectors
   `assembledInner h lam v x y = ∑ i, (‖v i‖² / 2) energyInner (x i) (y i)`
   polarizes the committed assembled energy (`assembledInner_self`), is
   preserved by `assembledFlow` for admissible data (`assembledInner_flow`),
   and is positive definite exactly when every mode is inside the window and
   every `v i` is nonzero: a Mathlib `InnerProductSpace.Core` is built
   (`assembledInnerCore`), the componentwise coordinates evolve by the phases
   `exp (i θ_i t / h)` (`assembledCoordinate_flow`) with derivative
   `i (θ_i / h)` (`assembledCoordinate_hasDerivAt`), and the coefficient
   state is recovered from the seam potential and electric field
   (`coefficients_of_fields`).  The gradient sector is the exact radical:
   for `lam i = 0` the amplitude direction is a nonzero null vector
   (`gradient_direction_null`), the diagonal of `assembledInner` vanishes
   precisely on the states supported on gradient amplitudes
   (`assembledInner_radical`), and the form on a gradient mode is not
   definite (`gradient_sector_not_definite`).  The Hilbert reading lives
   on the stated finite orthogonal family of curl modes, not on the whole
   curl sector (`orthogonal_family_hilbert_reading`); any symmetric
   bilinear form on the assembled state space whose diagonal is the
   committed assembled energy is `assembledInner`
   (`assembled_form_of_diagonal`); the four committed modes instantiate
   the reading at every declared `h ≠ 0` with `h² (3 + √5) < 4`
   (`curl_family_hilbert_instance`), in particular at `h = 1/2`
   (`curl_family_hilbert_instance_half`).
6. Degenerate eigenspaces.  On any two-mode block the step-invariant
   symmetric forms are at least two-dimensional (independent weights on
   the two summands), so invariance never fixes relative scales between
   modes.  What a degenerate block adds is a mode-mixing invariant form:
   the cross form `crossForm h lam x y = energyInner (x 0) (y 1) +
   energyInner (x 1) (y 0)` is symmetric and invariant under the block
   step exactly when the two eigenvalues coincide (`crossForm_invariant_iff`,
   declared `h ≠ 0`, `lam ≠ 0`), and neither it nor the plain direct sum is
   a multiple of the other (`degenerate_block_invariant_forms_not_unique`).
   Invariance therefore does not fix the orthogonal mode decomposition
   inside a degenerate eigenspace; the committed assembled energy does, by
   its diagonal (`assembled_form_of_diagonal`).

SIGN AND FRAME CONVENTIONS.  `rot α = !![cos α, -sin α; sin α, cos α]`,
`conjugator h lam = !![h, 0; -h² lam / 2, -sin θ]`, `θ = modeAngle h lam =
arccos (1 - h² lam / 2) ∈ (0, π)` as in `CarrierEvolutionFlow`.  With `J =
Q rot (π/2) Q⁻¹` the flow is multiplication by `exp (+ i θ t / h)`; the
conjugate structure `-J` is equally `energyInner`-orthogonal and turns the
phase into `exp (- i θ t / h)`, which for `h > 0` is the Schrödinger sign
`exp (- i H t)` with `H = θ / h > 0`.  Every statement below holds for
every declared `h ≠ 0`; for `h < 0` the number `θ / h` is negative and the
two orientations exchange roles.  The orientation of `J` is a convention;
the generator magnitude `|θ / h|` and the modulus-one phase are not.

BOUNDARY.  `h` is a declared Courant number; `t` a declared flow
parameter; the identification with physical time is the source clock and
duration row (not touched); the flow is on the coefficient state space of
the stated family and on the field span through `coefficients_of_fields`.
Not claimed: any Born rule, measurement, observable, or attachment of the
complex line to an observer; any statement about the whole curl sector
beyond the stated family (the selection of nineteen orthogonal seam
eigenvectors remains open as in `CarrierEvolutionFlow`); any physical
unit for `θ / h`.  Rows cited: OL-C2 (partial, quantum dynamics surface:
the field-sector flow is now a one-parameter unitary group with explicit
generator on a Hilbert space whose inner product is derived from the
committed energy, not declared) and PR-15 (no discharge).  Issues #730,
#733, #736.

FALSIFIER.  A nonzero state with nonpositive `modeForm` inside the window,
a step that fails to preserve `energyInner`, a symmetric invariant form
not proportional to `energyInner` for some `h ≠ 0`, a `J` with `J² ≠ -1`,
or a coordinate orbit whose derivative differs from `i (θ / h) z` would
make the module wrong.

Axiom audit.  The `#print axioms` lines at the end show at most `propext`,
`Classical.choice`, `Quot.sound`.
-/

noncomputable section

/-! ## 1. Positive definiteness of the mode form is exactly the strict window -/

theorem fin2_ne_zero_iff (x : Fin 2 → ℝ) : x ≠ 0 ↔ x 0 ≠ 0 ∨ x 1 ≠ 0 := by
  constructor
  · intro hx
    by_contra hc
    push Not at hc
    apply hx
    funext i
    fin_cases i
    · exact hc.1
    · exact hc.2
  · rintro h rfl
    simp at h

theorem modeForm_completed_square (h lam : ℝ) (x : Fin 2 → ℝ) :
    modeForm h lam x = lam * (x 0 + h * x 1 / 2) ^ 2 + (1 - h ^ 2 * lam / 4) * x 1 ^ 2 := by
  unfold modeForm
  ring

theorem lam_pos_of_window (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) : 0 < lam := by
  by_contra hc
  push Not at hc
  have := mul_nonpos_of_nonneg_of_nonpos (sq_nonneg h) hc
  linarith

theorem modeForm_nonneg_of_window (h lam : ℝ) (h4 : h ^ 2 * lam ≤ 4) (hl : 0 ≤ lam) (x : Fin 2 → ℝ) : 0 ≤ modeForm h lam x := by
  rw [modeForm_completed_square]
  have hc : 0 ≤ 1 - h ^ 2 * lam / 4 := by linarith
  positivity

/-- Inside the strict window the mode form is positive on nonzero states. -/
theorem modeForm_pos (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4)
    (x : Fin 2 → ℝ) (hx : x ≠ 0) : 0 < modeForm h lam x := by
  rw [modeForm_completed_square]
  have hl := lam_pos_of_window h lam h0
  have hc : 0 < 1 - h ^ 2 * lam / 4 := by linarith
  have hA : 0 ≤ lam * (x 0 + h * x 1 / 2) ^ 2 := mul_nonneg hl.le (sq_nonneg _)
  have hB : 0 ≤ (1 - h ^ 2 * lam / 4) * x 1 ^ 2 := mul_nonneg hc.le (sq_nonneg _)
  by_cases h2 : x 1 = 0
  · have h1 : x 0 ≠ 0 := by
      rcases (fin2_ne_zero_iff x).mp hx with h1 | h1
      · exact h1
      · exact absurd h2 h1
    have hsq : 0 < (x 0 + h * x 1 / 2) ^ 2 := by
      rw [h2, mul_zero, zero_div, add_zero]
      positivity
    have := mul_pos hl hsq
    linarith
  · have hsq : 0 < x 1 ^ 2 := by positivity
    have := mul_pos hc hsq
    linarith

/-- For `lam ≤ 0` the amplitude state `(1, 0)` has nonpositive form. -/
theorem modeForm_nonpos_of_lam_nonpos (h lam : ℝ) (hl : lam ≤ 0) :
    modeForm h lam ![1, 0] ≤ 0 := by
  unfold modeForm
  simp
  exact hl

/-- On the boundary `h² lam = 4` the nonzero state `(-h/2, 1)` has form zero. -/
theorem modeForm_boundary_zero (h lam : ℝ) (hb : h ^ 2 * lam = 4) :
    modeForm h lam ![-(h / 2), 1] = 0 := by
  rw [modeForm_completed_square]
  have : (1 : ℝ) - h ^ 2 * lam / 4 = 0 := by rw [hb]; norm_num
  rw [this]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

/-- On the boundary the form is the rank-one square `lam (a + h b / 2)²`. -/
theorem modeForm_boundary_rank_one (h lam : ℝ) (hb : h ^ 2 * lam = 4) (x : Fin 2 → ℝ) :
    modeForm h lam x = lam * (x 0 + h * x 1 / 2) ^ 2 := by
  rw [modeForm_completed_square]
  have : (1 : ℝ) - h ^ 2 * lam / 4 = 0 := by rw [hb]; norm_num
  rw [this, zero_mul, add_zero]

/-- Outside the window, `h² lam > 4`, the state `(-h/2, 1)` has negative form. -/
theorem modeForm_exterior_neg (h lam : ℝ) (hb : 4 < h ^ 2 * lam) :
    modeForm h lam ![-(h / 2), 1] < 0 := by
  rw [modeForm_completed_square]
  simp
  linarith

theorem boundary_state_ne_zero (h : ℝ) : (![-(h / 2), 1] : Fin 2 → ℝ) ≠ 0 := by
  intro hc
  have := congrFun hc 1
  simp at this

theorem amplitude_state_ne_zero : (![1, 0] : Fin 2 → ℝ) ≠ 0 := by
  intro hc
  have := congrFun hc 0
  simp at this

/-- **(1) Positive definiteness is exactly the strict window.**  For a
declared `h ≠ 0`, the mode form is positive on every nonzero state if and
only if `0 < h² lam < 4`. -/
theorem modeForm_posDef_iff (h lam : ℝ) (hh : h ≠ 0) :
    (∀ x : Fin 2 → ℝ, x ≠ 0 → 0 < modeForm h lam x) ↔ (0 < h ^ 2 * lam ∧ h ^ 2 * lam < 4) := by
  constructor
  · intro hpos
    have hsq : 0 < h ^ 2 := by positivity
    constructor
    · have h1 := hpos _ amplitude_state_ne_zero
      have hl : 0 < lam := by
        unfold modeForm at h1
        simpa using h1
      positivity
    · have h2 := hpos _ (boundary_state_ne_zero h)
      rw [modeForm_completed_square] at h2
      simp at h2
      linarith
  · rintro ⟨h0, h4⟩ x hx
    exact modeForm_pos h lam h0 h4 x hx

/-! ## 2. The energy inner product by polarization -/

/-- The energy inner product `lam a a' + (h lam / 2)(a b' + b a') + b b'`. -/
def energyInner (h lam : ℝ) (x y : Fin 2 → ℝ) : ℝ :=
  lam * x 0 * y 0 + h * lam / 2 * (x 0 * y 1 + x 1 * y 0) + x 1 * y 1

/-- Gram matrix of the energy inner product. -/
def energyGram (h lam : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![lam, h * lam / 2; h * lam / 2, 1]

theorem energyInner_eq_gram (h lam : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam x y = x ⬝ᵥ (energyGram h lam).mulVec y := by
  unfold energyInner energyGram
  rw [mulVec_fin_two]
  simp [dotProduct, Fin.sum_univ_two]
  ring

theorem energyInner_self (h lam : ℝ) (x : Fin 2 → ℝ) :
    energyInner h lam x x = modeForm h lam x := by
  unfold energyInner modeForm
  ring

/-- **(2) Polarization.**  The inner product is the polarization of the mode form. -/
theorem energyInner_polarization (h lam : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam x y = (modeForm h lam (x + y) - modeForm h lam x - modeForm h lam y) / 2 := by
  unfold energyInner modeForm
  simp only [Pi.add_apply]
  ring

theorem energyInner_comm (h lam : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam x y = energyInner h lam y x := by
  unfold energyInner
  ring

theorem energyInner_add_left (h lam : ℝ) (x y z : Fin 2 → ℝ) :
    energyInner h lam (x + y) z = energyInner h lam x z + energyInner h lam y z := by
  unfold energyInner
  simp only [Pi.add_apply]
  ring

theorem energyInner_smul_left (h lam : ℝ) (c : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam (c • x) y = c * energyInner h lam x y := by
  unfold energyInner
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

theorem energyInner_add_right (h lam : ℝ) (x y z : Fin 2 → ℝ) :
    energyInner h lam x (y + z) = energyInner h lam x y + energyInner h lam x z := by
  rw [energyInner_comm, energyInner_add_left, energyInner_comm h lam y, energyInner_comm h lam z]

theorem energyInner_smul_right (h lam : ℝ) (c : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam x (c • y) = c * energyInner h lam x y := by
  rw [energyInner_comm, energyInner_smul_left, energyInner_comm]

theorem energyInner_pos (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4)
    (x : Fin 2 → ℝ) (hx : x ≠ 0) : 0 < energyInner h lam x x := by
  rw [energyInner_self]
  exact modeForm_pos h lam h0 h4 x hx

theorem energyInner_self_eq_zero_iff (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4)
    (x : Fin 2 → ℝ) : energyInner h lam x x = 0 ↔ x = 0 := by
  constructor
  · intro hz
    by_contra hx
    exact (energyInner_pos h lam h0 h4 x hx).ne' hz
  · rintro rfl
    unfold energyInner
    simp

/-- **(2) The energy form is an inner product** (Mathlib `InnerProductSpace.Core`)
on the mode state space inside the strict window. -/
@[reducible] def energyInnerCore (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4) :
    InnerProductSpace.Core ℝ (Fin 2 → ℝ) where
  inner := energyInner h lam
  conj_inner_symm x y := by
    simp only [RCLike.conj_to_real]
    exact energyInner_comm h lam y x
  re_inner_nonneg x := by
    simp only [RCLike.re_to_real]
    rw [energyInner_self]
    exact modeForm_nonneg_of_window h lam h4.le (lam_pos_of_window h lam h0).le x
  add_left := energyInner_add_left h lam
  smul_left x y r := by
    simp only [RCLike.conj_to_real]
    exact energyInner_smul_left h lam r x y
  definite x hx := (energyInner_self_eq_zero_iff h lam h0 h4 x).mp hx

/-- **(2) The step is an isometry of the energy inner product**, for every
`h, lam` (a polynomial identity; no window is needed). -/
theorem stepMatrix_energy_isometry (h lam : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam ((stepMatrix h lam).mulVec x) ((stepMatrix h lam).mulVec y) =
      energyInner h lam x y := by
  unfold energyInner
  rw [stepMatrix_mulVec, stepMatrix_mulVec]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

/-- The step preserves the mode form for every `h, lam`. -/
theorem modeForm_stepMatrix (h lam : ℝ) (x : Fin 2 → ℝ) :
    modeForm h lam ((stepMatrix h lam).mulVec x) = modeForm h lam x := by
  rw [← energyInner_self, stepMatrix_energy_isometry, energyInner_self]

theorem rotFlow_energy_isometry (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x y : Fin 2 → ℝ) :
    energyInner h lam ((rotFlow h lam t).mulVec x) ((rotFlow h lam t).mulVec y) =
      energyInner h lam x y := by
  rw [energyInner_polarization, energyInner_polarization, ← Matrix.mulVec_add,
    modeForm_rotFlow h lam hh h0 h4, modeForm_rotFlow h lam hh h0 h4,
    modeForm_rotFlow h lam hh h0 h4]

/-- **(2) The per-mode flow is a one-parameter group of isometries.** -/
theorem modeFlow_energy_isometry (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) (t : ℝ)
    (x y : Fin 2 → ℝ) :
    energyInner h lam ((modeFlow h lam t).mulVec x) ((modeFlow h lam t).mulVec y) =
      energyInner h lam x y := by
  rw [energyInner_polarization, energyInner_polarization, ← Matrix.mulVec_add,
    modeForm_modeFlow h lam hh hadm, modeForm_modeFlow h lam hh hadm,
    modeForm_modeFlow h lam hh hadm]

theorem modeFlow_isometry_group (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) :
    modeFlow h lam 0 = 1 ∧
    (∀ s t, modeFlow h lam s * modeFlow h lam t = modeFlow h lam (s + t)) ∧
    modeFlow h lam h = stepMatrix h lam ∧
    Continuous (modeFlow h lam) ∧
    (∀ t x y, energyInner h lam ((modeFlow h lam t).mulVec x) ((modeFlow h lam t).mulVec y) =
      energyInner h lam x y) :=
  ⟨modeFlow_zero h lam hh hadm, modeFlow_add h lam hh hadm, modeFlow_step h lam hh hadm,
    modeFlow_continuous h lam, modeFlow_energy_isometry h lam hh hadm⟩

/-! ## 3. Uniqueness of the invariant symmetric form -/

theorem fin2_decomp (x : Fin 2 → ℝ) : x = x 0 • ![1, 0] + x 1 • ![0, 1] := by
  ext i
  fin_cases i <;> simp

/-- A bilinear form on `ℝ²` is determined by its four basis values. -/
theorem bilin_expand (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ)) (x y : Fin 2 → ℝ) :
    B x y = x 0 * y 0 * B ![1, 0] ![1, 0] + x 0 * y 1 * B ![1, 0] ![0, 1] +
      x 1 * y 0 * B ![0, 1] ![1, 0] + x 1 * y 1 * B ![0, 1] ![0, 1] := by
  conv_lhs => rw [fin2_decomp x, fin2_decomp y]
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
  ring

theorem stepMatrix_mulVec_e0 (h lam : ℝ) :
    (stepMatrix h lam).mulVec ![1, 0] = ![1, -(h * lam)] := by
  rw [stepMatrix_mulVec]
  ext i; fin_cases i <;> simp

theorem stepMatrix_mulVec_e1 (h lam : ℝ) :
    (stepMatrix h lam).mulVec ![0, 1] = ![h, 1 - h ^ 2 * lam] := by
  rw [stepMatrix_mulVec]
  ext i; fin_cases i
  · simp
  · simp; ring

/-- **(3) Uniqueness of the invariant form.**  For a declared `h ≠ 0` every
symmetric bilinear form invariant under the step is the multiple
`B (0,1) (0,1) · energyInner h lam`.  No window is needed: the step is not
a scalar matrix as soon as `h ≠ 0`, so its commutant is spanned by `1` and
itself, and symmetry fixes the combination. -/
theorem invariant_symmetric_form_eq_smul (h lam : ℝ) (hh : h ≠ 0)
    (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ)) (hsymm : ∀ x y, B x y = B y x)
    (hinv : ∀ x y, B ((stepMatrix h lam).mulVec x) ((stepMatrix h lam).mulVec y) = B x y)
    (x y : Fin 2 → ℝ) : B x y = B ![0, 1] ![0, 1] * energyInner h lam x y := by
  have hq' : B ![0, 1] ![1, 0] = B ![1, 0] ![0, 1] := hsymm _ _
  have e00 := hinv ![1, 0] ![1, 0]
  have e01 := hinv ![1, 0] ![0, 1]
  have e11 := hinv ![0, 1] ![0, 1]
  rw [stepMatrix_mulVec_e0, bilin_expand] at e00
  rw [stepMatrix_mulVec_e0, stepMatrix_mulVec_e1, bilin_expand] at e01
  rw [stepMatrix_mulVec_e1, bilin_expand] at e11
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one] at e00 e01 e11
  rw [hq'] at e00 e01 e11
  rw [bilin_expand B x y, hq']
  set p := B ![1, 0] ![1, 0] with hp
  set q := B ![1, 0] ![0, 1] with hq
  set r := B ![0, 1] ![0, 1] with hr
  have hA : h * lam * (h * lam * r - 2 * q) = 0 := by linear_combination e00
  have hB : h * (p - 2 * h * lam * q - lam * r + h ^ 2 * lam ^ 2 * r) = 0 := by
    linear_combination e01
  have hB' : p - 2 * h * lam * q - lam * r + h ^ 2 * lam ^ 2 * r = 0 :=
    (mul_eq_zero.mp hB).resolve_left hh
  by_cases hl : lam = 0
  · subst hl
    have hp0 : p = 0 := by linear_combination hB'
    have hq0 : q = 0 := by
      have : h * (h * p + 2 * q) = 0 := by linear_combination e11
      have h2 := (mul_eq_zero.mp this).resolve_left hh
      linear_combination h2 / 2 - h / 2 * hp0
    unfold energyInner
    rw [hp0, hq0]
    ring
  · have hl' : h * lam ≠ 0 := mul_ne_zero hh hl
    have hq1 : q = h * lam * r / 2 := by
      have := (mul_eq_zero.mp hA).resolve_left hl'
      linear_combination -this / 2
    have hp1 : p = lam * r := by
      rw [hq1] at hB'
      linear_combination hB'
    unfold energyInner
    rw [hp1, hq1]
    ring

/-- **(3) The invariant symmetric forms are exactly the real multiples of
`energyInner`** (both directions), for a declared `h ≠ 0`. -/
theorem invariant_symmetric_forms_iff (h lam : ℝ) (hh : h ≠ 0)
    (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ)) :
    ((∀ x y, B x y = B y x) ∧
      (∀ x y, B ((stepMatrix h lam).mulVec x) ((stepMatrix h lam).mulVec y) = B x y)) ↔
    ∃ c : ℝ, ∀ x y, B x y = c * energyInner h lam x y := by
  constructor
  · rintro ⟨hsymm, hinv⟩
    exact ⟨B ![0, 1] ![0, 1], invariant_symmetric_form_eq_smul h lam hh B hsymm hinv⟩
  · rintro ⟨c, hc⟩
    refine ⟨fun x y ↦ ?_, fun x y ↦ ?_⟩
    · rw [hc, hc, energyInner_comm]
    · rw [hc, hc, stepMatrix_energy_isometry]

/-- **(3) Inside the window a positive definite invariant form is a positive
multiple of `energyInner`.** -/
theorem invariant_inner_product_unique_up_to_positive_scale (h lam : ℝ) (hh : h ≠ 0)
    (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ)) (hsymm : ∀ x y, B x y = B y x)
    (hinv : ∀ x y, B ((stepMatrix h lam).mulVec x) ((stepMatrix h lam).mulVec y) = B x y)
    (hpos : ∀ x, x ≠ 0 → 0 < B x x) :
    ∃ c : ℝ, 0 < c ∧ ∀ x y, B x y = c * energyInner h lam x y := by
  refine ⟨B ![0, 1] ![0, 1], ?_, invariant_symmetric_form_eq_smul h lam hh B hsymm hinv⟩
  apply hpos
  intro hc
  have := congrFun hc 1
  simp at this

/-- Every state is the initial state of a scalar profile. -/
theorem state_realized (h : ℝ) (hh : h ≠ 0) (x : Fin 2 → ℝ) :
    ∃ a : ℕ → ℝ, ![a 0, velocity h a 0] = x := by
  refine ⟨fun n ↦ x 0 + (n : ℝ) * h * x 1, ?_⟩
  ext i
  fin_cases i
  · simp
  · simp only [velocity, Fin.mk_one, Fin.isValue, Matrix.cons_val_one]
    push_cast
    field_simp
    ring

/-- **(3) The staggered energy fixes the scale.**  A symmetric bilinear form
whose diagonal on every initial state is the committed staggered form of
the corresponding mode history at step `0` is `(‖v‖² / 2) · energyInner`.
The hypothesis ranges over every scalar profile `a`, solution of the
evolution or not: the step-`0` value of the committed form depends only on
`(a 0, a 1)`, i.e. on the state `(a 0, velocity h a 0)`.  Polarization
alone gives this; invariance is not used here. -/
theorem staggered_energy_fixes_scale (h lam : ℝ) (hh : h ≠ 0) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ))
    (hsymm : ∀ x y, B x y = B y x)
    (hdiag : ∀ a : ℕ → ℝ, B ![a 0, velocity h a 0] ![a 0, velocity h a 0] =
      fieldEnergyScaled h (scalarHistory a v) (fun _ ↦ 0) 0)
    (x y : Fin 2 → ℝ) : B x y = (realSeamEnergy v / 2) * energyInner h lam x y := by
  have hquad : ∀ z, B z z = (realSeamEnergy v / 2) * modeForm h lam z := by
    intro z
    obtain ⟨a, ha⟩ := state_realized h hh z
    rw [← ha, hdiag a, fieldEnergyScaled_scalarHistory h hh a v lam hv 0]
  have hpol : B x y = (B (x + y) (x + y) - B x x - B y y) / 2 := by
    simp only [map_add, LinearMap.add_apply]
    rw [hsymm y x]
    ring
  rw [hpol, hquad, hquad, hquad, energyInner_polarization]
  ring

/-- **(3) The inner product with the committed staggered energy as diagonal.**
Symmetric, step-invariant, with the committed staggered energy as diagonal,
if and only if equal to `(‖v‖² / 2) · energyInner h lam`.  The forward
direction uses symmetry and the diagonal only (polarization); step
invariance is then a consequence, not an input. -/
theorem staggered_energy_inner_product_iff (h lam : ℝ) (hh : h ≠ 0) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) (B : LinearMap.BilinForm ℝ (Fin 2 → ℝ)) :
    ((∀ x y, B x y = B y x) ∧
      (∀ x y, B ((stepMatrix h lam).mulVec x) ((stepMatrix h lam).mulVec y) = B x y) ∧
      (∀ a : ℕ → ℝ, B ![a 0, velocity h a 0] ![a 0, velocity h a 0] =
        fieldEnergyScaled h (scalarHistory a v) (fun _ ↦ 0) 0)) ↔
    ∀ x y, B x y = (realSeamEnergy v / 2) * energyInner h lam x y := by
  constructor
  · rintro ⟨hsymm, _, hdiag⟩
    exact staggered_energy_fixes_scale h lam hh v hv B hsymm hdiag
  · intro hB
    refine ⟨fun x y ↦ ?_, fun x y ↦ ?_, fun a ↦ ?_⟩
    · rw [hB, hB, energyInner_comm]
    · rw [hB, hB, stepMatrix_energy_isometry]
    · rw [hB, energyInner_self, fieldEnergyScaled_scalarHistory h hh a v lam hv 0]

/-! ## 4. Complex structure, complex coordinate, unitarity, explicit generator -/

/-- The complex structure `J = (1 / sin θ) !![h² lam / 2, h; -h lam, -h² lam / 2]`. -/
def complexStructure (h lam : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  (Real.sin (modeAngle h lam))⁻¹ • !![h ^ 2 * lam / 2, h; -(h * lam), -(h ^ 2 * lam / 2)]

/-- The integer-coefficient core `!![h² lam / 2, h; -h lam, -h² lam / 2]` squares to
`-sin² θ` times the identity. -/
theorem complexStructure_core_sq (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4) :
    (!![h ^ 2 * lam / 2, h; -(h * lam), -(h ^ 2 * lam / 2)] : Matrix (Fin 2) (Fin 2) ℝ) *
      !![h ^ 2 * lam / 2, h; -(h * lam), -(h ^ 2 * lam / 2)] =
      (-(Real.sin (modeAngle h lam) ^ 2)) • (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
  have hsq := sin_sq_modeAngle h lam h0 h4
  rw [Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j
  · simp
    linear_combination (1 : ℝ) * hsq
  · simp
    ring
  · simp
    ring
  · simp
    linear_combination (1 : ℝ) * hsq

theorem rot_pi_div_two : rot (Real.pi / 2) = !![0, -1; 1, 0] := by
  unfold rot
  rw [Real.cos_pi_div_two, Real.sin_pi_div_two]

/-- **(4)** `J² = -1` inside the window. -/
theorem complexStructure_sq (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4) :
    complexStructure h lam * complexStructure h lam = -1 := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  unfold complexStructure
  rw [Matrix.smul_mul, Matrix.mul_smul, complexStructure_core_sq h lam h0 h4, smul_smul,
    smul_smul]
  rw [show (Real.sin (modeAngle h lam))⁻¹ * (Real.sin (modeAngle h lam))⁻¹ *
    -(Real.sin (modeAngle h lam) ^ 2) = -1 by field_simp]
  simp

/-- **(4)** `J` is the conjugated quarter turn `Q · rot (π/2) · Q⁻¹`. -/
theorem complexStructure_eq_conj (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) :
    complexStructure h lam =
      conjugator h lam * rot (Real.pi / 2) * conjugatorInv h lam := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  have hsq := sin_sq_modeAngle h lam h0 h4
  unfold complexStructure conjugator conjugatorInv
  rw [rot_pi_div_two, Matrix.mul_fin_two, Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j
  · simp
    field_simp
  · simp
    field_simp
  · simp
    field_simp
    linear_combination (4 : ℝ) * hsq
  · simp
    field_simp

theorem rot_eq_cos_add_sin (α : ℝ) :
    rot α = Real.cos α • (1 : Matrix (Fin 2) (Fin 2) ℝ) + Real.sin α • rot (Real.pi / 2) := by
  rw [rot_pi_div_two]
  unfold rot
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- **(4) The flow is `cos (θ t / h) · 1 + sin (θ t / h) · J`.** -/
theorem rotFlow_eq_cos_add_sin (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) :
    rotFlow h lam t = Real.cos (modeAngle h lam * t / h) • (1 : Matrix (Fin 2) (Fin 2) ℝ) +
      Real.sin (modeAngle h lam * t / h) • complexStructure h lam := by
  rw [complexStructure_eq_conj h lam hh h0 h4]
  unfold rotFlow
  rw [rot_eq_cos_add_sin, Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one, conjugator_mul_inv h lam hh h0 h4]

/-- `Q⁻¹ J = rot (π/2) Q⁻¹` on vectors. -/
theorem conjugatorInv_mulVec_complexStructure (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    (conjugatorInv h lam).mulVec ((complexStructure h lam).mulVec x) =
      (rot (Real.pi / 2)).mulVec ((conjugatorInv h lam).mulVec x) := by
  rw [complexStructure_eq_conj h lam hh h0 h4, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
    ← Matrix.mul_assoc, ← Matrix.mul_assoc, conjugatorInv_mul h lam hh h0 h4, Matrix.one_mul]

/-- `Q⁻¹ · rotFlow t = rot (θ t / h) · Q⁻¹` on vectors. -/
theorem conjugatorInv_mulVec_rotFlow (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x : Fin 2 → ℝ) :
    (conjugatorInv h lam).mulVec ((rotFlow h lam t).mulVec x) =
      (rot (modeAngle h lam * t / h)).mulVec ((conjugatorInv h lam).mulVec x) := by
  unfold rotFlow
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec, ← Matrix.mul_assoc, ← Matrix.mul_assoc,
    conjugatorInv_mul h lam hh h0 h4, Matrix.one_mul]

/-- The energy inner product is `sin² θ` times the Euclidean pairing of the
pulled-back coordinates. -/
theorem energyInner_eq_pulled (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    energyInner h lam x y = Real.sin (modeAngle h lam) ^ 2 *
      ((conjugatorInv h lam).mulVec x 0 * (conjugatorInv h lam).mulVec y 0 +
        (conjugatorInv h lam).mulVec x 1 * (conjugatorInv h lam).mulVec y 1) := by
  rw [energyInner_polarization, ← sin_sq_mul_pulledForm h lam hh h0 h4,
    ← sin_sq_mul_pulledForm h lam hh h0 h4, ← sin_sq_mul_pulledForm h lam hh h0 h4]
  unfold pulledForm sq2
  rw [Matrix.mulVec_add]
  simp only [Pi.add_apply]
  ring

/-- **(4) `J` is `energyInner`-orthogonal.** -/
theorem complexStructure_energy_isometry (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    energyInner h lam ((complexStructure h lam).mulVec x) ((complexStructure h lam).mulVec y) =
      energyInner h lam x y := by
  rw [energyInner_eq_pulled h lam hh h0 h4, energyInner_eq_pulled h lam hh h0 h4,
    conjugatorInv_mulVec_complexStructure h lam hh h0 h4,
    conjugatorInv_mulVec_complexStructure h lam hh h0 h4, rot_pi_div_two, mulVec_fin_two,
    mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

/-- **(4) `J` is `energyInner`-skew:** `⟨J x, y⟩ = -⟨x, J y⟩`. -/
theorem complexStructure_skew (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    energyInner h lam ((complexStructure h lam).mulVec x) y =
      -energyInner h lam x ((complexStructure h lam).mulVec y) := by
  rw [energyInner_eq_pulled h lam hh h0 h4, energyInner_eq_pulled h lam hh h0 h4,
    conjugatorInv_mulVec_complexStructure h lam hh h0 h4,
    conjugatorInv_mulVec_complexStructure h lam hh h0 h4, rot_pi_div_two, mulVec_fin_two,
    mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

theorem energyInner_complexStructure_self (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    energyInner h lam x ((complexStructure h lam).mulVec x) = 0 := by
  have h1 := complexStructure_skew h lam hh h0 h4 x x
  rw [energyInner_comm] at h1
  linarith

/-- The complex coordinate `z = (Q⁻¹ x)₀ + (Q⁻¹ x)₁ i`. -/
def modeCoordinate (h lam : ℝ) (x : Fin 2 → ℝ) : ℂ :=
  ((conjugatorInv h lam).mulVec x 0 : ℂ) + ((conjugatorInv h lam).mulVec x 1 : ℂ) * Complex.I

theorem modeCoordinate_add (h lam : ℝ) (x y : Fin 2 → ℝ) :
    modeCoordinate h lam (x + y) = modeCoordinate h lam x + modeCoordinate h lam y := by
  unfold modeCoordinate
  rw [Matrix.mulVec_add]
  simp only [Pi.add_apply]
  push_cast
  ring

theorem modeCoordinate_smul (h lam : ℝ) (c : ℝ) (x : Fin 2 → ℝ) :
    modeCoordinate h lam (c • x) = (c : ℂ) * modeCoordinate h lam x := by
  unfold modeCoordinate
  rw [Matrix.mulVec_smul]
  simp only [Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

/-- **(4) `J` is multiplication by `i` in the complex coordinate.** -/
theorem modeCoordinate_complexStructure (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    modeCoordinate h lam ((complexStructure h lam).mulVec x) =
      Complex.I * modeCoordinate h lam x := by
  unfold modeCoordinate
  rw [conjugatorInv_mulVec_complexStructure h lam hh h0 h4, rot_pi_div_two, mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  apply Complex.ext <;> simp

/-- **(4) The flow is multiplication by the unit phase `exp (i θ t / h)`.** -/
theorem modeCoordinate_rotFlow (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x : Fin 2 → ℝ) :
    modeCoordinate h lam ((rotFlow h lam t).mulVec x) =
      Complex.exp ((modeAngle h lam * t / h : ℝ) * Complex.I) * modeCoordinate h lam x := by
  unfold modeCoordinate
  rw [conjugatorInv_mulVec_rotFlow h lam hh h0 h4, Complex.exp_ofReal_mul_I]
  unfold rot
  rw [mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  apply Complex.ext <;> simp <;> ring

theorem modeCoordinate_eq_zero_iff (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) : modeCoordinate h lam x = 0 ↔ x = 0 := by
  constructor
  · intro hz
    have hre := congrArg Complex.re hz
    have him := congrArg Complex.im hz
    unfold modeCoordinate at hre him
    simp at hre him
    have hu : (conjugatorInv h lam).mulVec x = 0 := by
      ext i; fin_cases i
      · exact hre
      · exact him
    have : x = (conjugator h lam * conjugatorInv h lam).mulVec x := by
      rw [conjugator_mul_inv h lam hh h0 h4, Matrix.one_mulVec]
    rw [this, ← Matrix.mulVec_mulVec, hu, Matrix.mulVec_zero]
  · rintro rfl
    unfold modeCoordinate
    simp

theorem modeCoordinate_surjective (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (z : ℂ) : ∃ x, modeCoordinate h lam x = z := by
  refine ⟨(conjugator h lam).mulVec ![z.re, z.im], ?_⟩
  unfold modeCoordinate
  rw [Matrix.mulVec_mulVec, conjugatorInv_mul h lam hh h0 h4, Matrix.one_mulVec]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  exact Complex.re_add_im z

/-- **(4) The complex coordinate is an `ℝ`-linear bijection `ℝ² → ℂ`.** -/
theorem modeCoordinate_bijective (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : Function.Bijective (modeCoordinate h lam) := by
  constructor
  · intro x y hxy
    have : modeCoordinate h lam (x - y) = 0 := by
      rw [sub_eq_add_neg, modeCoordinate_add, ← neg_one_smul ℝ y, modeCoordinate_smul, hxy]
      push_cast
      ring
    exact sub_eq_zero.mp ((modeCoordinate_eq_zero_iff h lam hh h0 h4 _).mp this)
  · exact modeCoordinate_surjective h lam hh h0 h4

/-- The sesquilinear form `⟨x, y⟩ = energyInner x y - i · energyInner x (J y)`. -/
def modeHermitian (h lam : ℝ) (x y : Fin 2 → ℝ) : ℂ :=
  (energyInner h lam x y : ℂ) -
    (energyInner h lam x ((complexStructure h lam).mulVec y) : ℂ) * Complex.I

/-- **(4) In the coordinate the form is `sin² θ · conj z · w`.** -/
theorem modeHermitian_eq_coordinate (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    modeHermitian h lam x y = (Real.sin (modeAngle h lam) : ℂ) ^ 2 *
      (starRingEnd ℂ) (modeCoordinate h lam x) * modeCoordinate h lam y := by
  unfold modeHermitian modeCoordinate
  rw [energyInner_eq_pulled h lam hh h0 h4, energyInner_eq_pulled h lam hh h0 h4,
    conjugatorInv_mulVec_complexStructure h lam hh h0 h4, rot_pi_div_two, mulVec_fin_two,
    map_add, map_mul, Complex.conj_ofReal, Complex.conj_ofReal, Complex.conj_I]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  apply Complex.ext <;> simp <;> ring

theorem modeHermitian_self (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    modeHermitian h lam x x = (modeForm h lam x : ℂ) := by
  unfold modeHermitian
  rw [energyInner_complexStructure_self h lam hh h0 h4, energyInner_self]
  push_cast
  ring

theorem modeHermitian_conj_symm (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    (starRingEnd ℂ) (modeHermitian h lam y x) = modeHermitian h lam x y := by
  rw [modeHermitian_eq_coordinate h lam hh h0 h4, modeHermitian_eq_coordinate h lam hh h0 h4,
    map_mul, map_mul, Complex.conj_conj, map_pow, Complex.conj_ofReal]
  ring

theorem modeHermitian_complexStructure_right (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    modeHermitian h lam x ((complexStructure h lam).mulVec y) =
      Complex.I * modeHermitian h lam x y := by
  rw [modeHermitian_eq_coordinate h lam hh h0 h4, modeHermitian_eq_coordinate h lam hh h0 h4,
    modeCoordinate_complexStructure h lam hh h0 h4]
  ring

theorem modeHermitian_complexStructure_left (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x y : Fin 2 → ℝ) :
    modeHermitian h lam ((complexStructure h lam).mulVec x) y =
      -Complex.I * modeHermitian h lam x y := by
  rw [modeHermitian_eq_coordinate h lam hh h0 h4, modeHermitian_eq_coordinate h lam hh h0 h4,
    modeCoordinate_complexStructure h lam hh h0 h4, map_mul, Complex.conj_I]
  ring

theorem modeHermitian_add_right (h lam : ℝ) (x y z : Fin 2 → ℝ) :
    modeHermitian h lam x (y + z) = modeHermitian h lam x y + modeHermitian h lam x z := by
  unfold modeHermitian
  rw [Matrix.mulVec_add, energyInner_add_right, energyInner_add_right]
  push_cast
  ring

theorem modeHermitian_smul_right (h lam : ℝ) (c : ℝ) (x y : Fin 2 → ℝ) :
    modeHermitian h lam x (c • y) = (c : ℂ) * modeHermitian h lam x y := by
  unfold modeHermitian
  rw [Matrix.mulVec_smul, energyInner_smul_right, energyInner_smul_right]
  push_cast
  ring

theorem conj_phase_mul_phase (α : ℝ) :
    (starRingEnd ℂ) (Complex.exp ((α : ℂ) * Complex.I)) * Complex.exp ((α : ℂ) * Complex.I) = 1 := by
  rw [mul_comm, Complex.mul_conj, Complex.normSq_eq_norm_sq, Complex.norm_exp_ofReal_mul_I]
  norm_num

/-- **(4) Unitarity.**  The flow preserves the sesquilinear form. -/
theorem modeHermitian_rotFlow (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x y : Fin 2 → ℝ) :
    modeHermitian h lam ((rotFlow h lam t).mulVec x) ((rotFlow h lam t).mulVec y) =
      modeHermitian h lam x y := by
  rw [modeHermitian_eq_coordinate h lam hh h0 h4, modeHermitian_eq_coordinate h lam hh h0 h4,
    modeCoordinate_rotFlow h lam hh h0 h4, modeCoordinate_rotFlow h lam hh h0 h4, map_mul]
  have hu := conj_phase_mul_phase (modeAngle h lam * t / h)
  linear_combination ((Real.sin (modeAngle h lam) : ℂ) ^ 2 *
    (starRingEnd ℂ) (modeCoordinate h lam x) * modeCoordinate h lam y) * hu

/-- The phase has modulus one. -/
theorem phase_norm_one (h lam t : ℝ) :
    ‖Complex.exp ((modeAngle h lam * t / h : ℝ) * Complex.I)‖ = 1 :=
  Complex.norm_exp_ofReal_mul_I _

/-- **(4) Explicit generator.**  The orbit of every coordinate is differentiable
in the flow parameter with derivative `i (θ / h) · z(t)`. -/
theorem modeCoordinate_flow_hasDerivAt (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) (t : ℝ) :
    HasDerivAt (fun s : ℝ ↦ modeCoordinate h lam ((rotFlow h lam s).mulVec x))
      (Complex.I * (modeAngle h lam / h : ℝ) * modeCoordinate h lam ((rotFlow h lam t).mulVec x))
      t := by
  have hfun : (fun s : ℝ ↦ modeCoordinate h lam ((rotFlow h lam s).mulVec x)) =
      fun s : ℝ ↦ Complex.exp ((modeAngle h lam * s / h : ℝ) * Complex.I) *
        modeCoordinate h lam x :=
    funext fun s ↦ modeCoordinate_rotFlow h lam hh h0 h4 s x
  rw [hfun, modeCoordinate_rotFlow h lam hh h0 h4]
  have h1 : HasDerivAt (fun s : ℝ ↦ modeAngle h lam * s / h) (modeAngle h lam / h) t := by
    have := ((hasDerivAt_id t).const_mul (modeAngle h lam)).div_const h
    simpa using this
  have h2 : HasDerivAt (fun s : ℝ ↦ ((modeAngle h lam * s / h : ℝ) : ℂ) * Complex.I)
      (((modeAngle h lam / h : ℝ) : ℂ) * Complex.I) t :=
    h1.ofReal_comp.mul_const Complex.I
  have h3 := (h2.cexp).mul_const (modeCoordinate h lam x)
  convert h3 using 1
  ring

/-- **(4) One-dimensional complex Hilbert line with unitary flow and explicit
generator.**  Inside the window: the coordinate is an `ℝ`-linear bijection to
`ℂ` under which `J` is `i`; the sesquilinear form is `sin² θ · conj z · w`,
conjugate-symmetric, with real positive diagonal `modeForm`; the flow is the
unit phase `exp (i θ t / h)`, preserves the form, and every orbit has
derivative `i (θ / h) z`.  Stone's theorem is not invoked; the generator
`θ / h` is exhibited.  Sign: with `J = Q rot (π/2) Q⁻¹` the phase is
`exp (+ i θ t / h)`; the conjugate structure `-J` gives `exp (- i θ t / h)`;
`θ / h` is positive for `h > 0` and negative for `h < 0`. -/
theorem mode_unitary_group_explicit_generator (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) :
    Function.Bijective (modeCoordinate h lam) ∧
    (∀ x y, modeCoordinate h lam (x + y) = modeCoordinate h lam x + modeCoordinate h lam y) ∧
    (∀ (c : ℝ) x, modeCoordinate h lam (c • x) = (c : ℂ) * modeCoordinate h lam x) ∧
    (∀ x, modeCoordinate h lam ((complexStructure h lam).mulVec x) =
      Complex.I * modeCoordinate h lam x) ∧
    complexStructure h lam * complexStructure h lam = -1 ∧
    (∀ x y, modeHermitian h lam x y = (Real.sin (modeAngle h lam) : ℂ) ^ 2 *
      (starRingEnd ℂ) (modeCoordinate h lam x) * modeCoordinate h lam y) ∧
    (∀ x y, (starRingEnd ℂ) (modeHermitian h lam y x) = modeHermitian h lam x y) ∧
    (∀ x, modeHermitian h lam x x = (modeForm h lam x : ℂ)) ∧
    (∀ x, x ≠ 0 → 0 < modeForm h lam x) ∧
    (∀ t x, modeCoordinate h lam ((rotFlow h lam t).mulVec x) =
      Complex.exp ((modeAngle h lam * t / h : ℝ) * Complex.I) * modeCoordinate h lam x) ∧
    (∀ t, ‖Complex.exp ((modeAngle h lam * t / h : ℝ) * Complex.I)‖ = 1) ∧
    (∀ t x y, modeHermitian h lam ((rotFlow h lam t).mulVec x) ((rotFlow h lam t).mulVec y) =
      modeHermitian h lam x y) ∧
    (∀ x t, HasDerivAt (fun s : ℝ ↦ modeCoordinate h lam ((rotFlow h lam s).mulVec x))
      (Complex.I * (modeAngle h lam / h : ℝ) *
        modeCoordinate h lam ((rotFlow h lam t).mulVec x)) t) :=
  ⟨modeCoordinate_bijective h lam hh h0 h4, modeCoordinate_add h lam, modeCoordinate_smul h lam,
    modeCoordinate_complexStructure h lam hh h0 h4, complexStructure_sq h lam h0 h4,
    modeHermitian_eq_coordinate h lam hh h0 h4, modeHermitian_conj_symm h lam hh h0 h4,
    modeHermitian_self h lam hh h0 h4, modeForm_pos h lam h0 h4,
    modeCoordinate_rotFlow h lam hh h0 h4, phase_norm_one h lam,
    modeHermitian_rotFlow h lam hh h0 h4, modeCoordinate_flow_hasDerivAt h lam hh h0 h4⟩


/-! ## 5. The assembled sector: inner product on the curl sector, radical on the gradient sector -/

theorem modeForm_zero_lam (h : ℝ) (x : Fin 2 → ℝ) : modeForm h 0 x = x 1 ^ 2 := by
  unfold modeForm
  ring

/-- On a gradient mode the form vanishes on the whole amplitude direction. -/
theorem gradient_form_vanishes (h a : ℝ) : modeForm h 0 ![a, 0] = 0 := by
  rw [modeForm_zero_lam]
  simp

/-- **(5) The gradient sector is not definite.** -/
theorem gradient_sector_not_definite (h : ℝ) :
    ¬ ∀ x : Fin 2 → ℝ, x ≠ 0 → 0 < modeForm h 0 x := by
  intro hc
  have := hc _ amplitude_state_ne_zero
  rw [gradient_form_vanishes] at this
  exact lt_irrefl 0 this

theorem modeForm_nonneg_of_admissible (h lam : ℝ) (hadm : Admissible h lam) (x : Fin 2 → ℝ) :
    0 ≤ modeForm h lam x := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · rw [modeForm_zero_lam]
    positivity
  · exact modeForm_nonneg_of_window h lam h4.le (lam_pos_of_window h lam h0).le x

theorem realSeamEnergy_pos_of_ne_zero (v : Fin 30 → ℝ) (hv : v ≠ 0) : 0 < realSeamEnergy v :=
  lt_of_le_of_ne (realSeamEnergy_nonneg v) fun hz ↦ hv ((realSeamEnergy_eq_zero_iff v).mp hz.symm)

section Assembled

variable {ι : Type} [Fintype ι]

/-- The assembled inner product `∑ i, (‖v i‖² / 2) energyInner h (lam i) (x i) (y i)`. -/
def assembledInner (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x y : ι → Fin 2 → ℝ) : ℝ :=
  ∑ i, (realSeamEnergy (v i) / 2) * energyInner h (lam i) (x i) (y i)

/-- **(5) The assembled inner product polarizes the committed assembled energy.** -/
theorem assembledInner_self (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x : ι → Fin 2 → ℝ) :
    assembledInner h lam v x x = assembledEnergy h lam v x := by
  unfold assembledInner assembledEnergy
  simp only [energyInner_self]

theorem assembledInner_comm (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x y : ι → Fin 2 → ℝ) :
    assembledInner h lam v x y = assembledInner h lam v y x := by
  unfold assembledInner
  exact Finset.sum_congr rfl fun i _ ↦ by rw [energyInner_comm]

theorem assembledInner_add_left (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (x y z : ι → Fin 2 → ℝ) :
    assembledInner h lam v (x + y) z = assembledInner h lam v x z + assembledInner h lam v y z := by
  unfold assembledInner
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ ↦ by
    simp only [Pi.add_apply]
    rw [energyInner_add_left]
    ring

theorem assembledInner_smul_left (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (c : ℝ)
    (x y : ι → Fin 2 → ℝ) :
    assembledInner h lam v (c • x) y = c * assembledInner h lam v x y := by
  unfold assembledInner
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ ↦ by
    simp only [Pi.smul_apply]
    rw [energyInner_smul_left]
    ring

theorem assembledInner_add_right (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (x y z : ι → Fin 2 → ℝ) :
    assembledInner h lam v x (y + z) = assembledInner h lam v x y + assembledInner h lam v x z := by
  rw [assembledInner_comm, assembledInner_add_left, assembledInner_comm h lam v y,
    assembledInner_comm h lam v z]

/-- The assembled inner product is the polarization of its diagonal. -/
theorem assembledInner_polarization (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (x y : ι → Fin 2 → ℝ) :
    assembledInner h lam v x y = (assembledInner h lam v (x + y) (x + y) -
      assembledInner h lam v x x - assembledInner h lam v y y) / 2 := by
  rw [assembledInner_add_left, assembledInner_add_right, assembledInner_add_right,
    assembledInner_comm h lam v y x]
  ring

/-- **(5) The committed assembled energy determines the assembled inner product.**
Every symmetric bilinear form on the assembled state space whose diagonal is
the committed assembled energy is `assembledInner` (polarization; no
invariance hypothesis, no window). -/
theorem assembled_form_of_diagonal (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (B : LinearMap.BilinForm ℝ (ι → Fin 2 → ℝ)) (hsymm : ∀ x y, B x y = B y x)
    (hdiag : ∀ x, B x x = assembledEnergy h lam v x) (x y : ι → Fin 2 → ℝ) :
    B x y = assembledInner h lam v x y := by
  have hpol : B x y = (B (x + y) (x + y) - B x x - B y y) / 2 := by
    simp only [map_add, LinearMap.add_apply]
    rw [hsymm y x]
    ring
  rw [hpol, hdiag, hdiag, hdiag, ← assembledInner_self h lam v (x + y),
    ← assembledInner_self h lam v x, ← assembledInner_self h lam v y,
    ← assembledInner_polarization]

/-- **(5) The assembled flow preserves the assembled inner product.** -/
theorem assembledInner_flow (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (t : ℝ) (x y : ι → Fin 2 → ℝ) :
    assembledInner h lam v (assembledFlow h lam t x) (assembledFlow h lam t y) =
      assembledInner h lam v x y := by
  unfold assembledInner assembledFlow
  exact Finset.sum_congr rfl fun i _ ↦ by rw [modeFlow_energy_isometry h (lam i) hh (hadm i)]

theorem assembledInner_nonneg (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) : 0 ≤ assembledInner h lam v x x := by
  unfold assembledInner
  refine Finset.sum_nonneg fun i _ ↦ ?_
  rw [energyInner_self]
  exact mul_nonneg (by have := realSeamEnergy_nonneg (v i); linarith)
    (modeForm_nonneg_of_admissible h (lam i) (hadm i) (x i))

/-- **(5) Positive definiteness on the curl sector.**  With every mode inside the
window and every `v i` nonzero the assembled form is positive on nonzero states. -/
theorem assembledInner_pos (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (hv : ∀ i, v i ≠ 0)
    (x : ι → Fin 2 → ℝ) (hx : x ≠ 0) : 0 < assembledInner h lam v x x := by
  obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
  unfold assembledInner
  refine Finset.sum_pos' (fun j _ ↦ ?_) ⟨i, Finset.mem_univ i, ?_⟩
  · rw [energyInner_self]
    exact mul_nonneg (by have := realSeamEnergy_nonneg (v j); linarith)
      (modeForm_nonneg_of_window h (lam j) (hwin j).2.le (lam_pos_of_window h (lam j) (hwin j).1).le
        (x j))
  · rw [energyInner_self]
    exact mul_pos (by have := realSeamEnergy_pos_of_ne_zero (v i) (hv i); linarith)
      (modeForm_pos h (lam i) (hwin i).1 (hwin i).2 (x i) hi)

/-- **(5) The assembled inner product is an inner product** (Mathlib
`InnerProductSpace.Core`) on the coefficient state space of a curl-sector family. -/
@[reducible] def assembledInnerCore (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (hv : ∀ i, v i ≠ 0) :
    InnerProductSpace.Core ℝ (ι → Fin 2 → ℝ) where
  inner := assembledInner h lam v
  conj_inner_symm x y := by
    simp only [RCLike.conj_to_real]
    exact assembledInner_comm h lam v y x
  re_inner_nonneg x := by
    simp only [RCLike.re_to_real]
    exact assembledInner_nonneg h lam v (fun i ↦ Or.inr (hwin i)) x
  add_left := assembledInner_add_left h lam v
  smul_left x y r := by
    simp only [RCLike.conj_to_real]
    exact assembledInner_smul_left h lam v r x y
  definite x hx := by
    by_contra hne
    exact (assembledInner_pos h lam v hwin hv x hne).ne' hx

/-- **(5) The gradient amplitude direction is a nonzero null vector.** -/
theorem gradient_direction_null [DecidableEq ι] (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (i : ι) (hi : lam i = 0) :
    (Pi.single i ![1, 0] : ι → Fin 2 → ℝ) ≠ 0 ∧
      assembledInner h lam v (Pi.single i ![1, 0]) (Pi.single i ![1, 0]) = 0 := by
  constructor
  · intro hc
    have := congrFun hc i
    rw [Pi.single_eq_same] at this
    exact amplitude_state_ne_zero this
  · unfold assembledInner
    rw [Finset.sum_eq_single i]
    · rw [Pi.single_eq_same, energyInner_self, hi, gradient_form_vanishes, mul_zero]
    · intro j _ hji
      rw [Pi.single_eq_of_ne hji, energyInner_self]
      unfold modeForm
      simp
    · intro hi'
      exact absurd (Finset.mem_univ i) hi'

/-- **(5) The radical of the assembled form is exactly the gradient amplitudes.**
For admissible data with nonzero `v i`, the diagonal vanishes if and only if
every gradient mode has zero velocity and every curl mode has zero state. -/
theorem assembledInner_radical (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (hadm : ∀ i, Admissible h (lam i)) (hv : ∀ i, v i ≠ 0) (x : ι → Fin 2 → ℝ) :
    assembledInner h lam v x x = 0 ↔
      ∀ i, (lam i = 0 → x i 1 = 0) ∧ (lam i ≠ 0 → x i = 0) := by
  unfold assembledInner
  have hterm : ∀ i, 0 ≤ (realSeamEnergy (v i) / 2) * energyInner h (lam i) (x i) (x i) := by
    intro i
    rw [energyInner_self]
    exact mul_nonneg (by have := realSeamEnergy_nonneg (v i); linarith)
      (modeForm_nonneg_of_admissible h (lam i) (hadm i) (x i))
  rw [Finset.sum_eq_zero_iff_of_nonneg fun i _ ↦ hterm i]
  have hvpos : ∀ i, realSeamEnergy (v i) / 2 ≠ 0 := fun i ↦ by
    have := realSeamEnergy_pos_of_ne_zero (v i) (hv i)
    positivity
  constructor
  · intro hz i
    have hi := hz i (Finset.mem_univ i)
    rw [energyInner_self] at hi
    have hm : modeForm h (lam i) (x i) = 0 := (mul_eq_zero.mp hi).resolve_left (hvpos i)
    rcases hadm i with hl | ⟨h0, h4⟩
    · refine ⟨fun _ ↦ ?_, fun hne ↦ absurd hl hne⟩
      rw [hl, modeForm_zero_lam] at hm
      exact pow_eq_zero_iff two_ne_zero |>.mp hm
    · have hne : lam i ≠ 0 := by
        rintro hl
        rw [hl] at h0
        simp at h0
      refine ⟨fun hl ↦ absurd hl hne, fun _ ↦ ?_⟩
      by_contra hx
      exact (modeForm_pos h (lam i) h0 h4 (x i) hx).ne' hm
  · intro hz i _
    rw [energyInner_self]
    rcases hadm i with hl | ⟨h0, h4⟩
    · rw [hl, modeForm_zero_lam, (hz i).1 hl]
      ring
    · have hne : lam i ≠ 0 := by
        rintro hl
        rw [hl] at h0
        simp at h0
      rw [(hz i).2 hne]
      unfold modeForm
      simp

/-- Componentwise complex coordinates. -/
def assembledCoordinate (h : ℝ) (lam : ι → ℝ) (x : ι → Fin 2 → ℝ) (i : ι) : ℂ :=
  modeCoordinate h (lam i) (x i)

omit [Fintype ι] in
/-- **(5) The assembled flow is the diagonal phase `exp (i θ_i t / h)`.** -/
theorem assembledCoordinate_flow (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (t : ℝ) (x : ι → Fin 2 → ℝ) (i : ι) :
    assembledCoordinate h lam (assembledFlow h lam t x) i =
      Complex.exp ((modeAngle h (lam i) * t / h : ℝ) * Complex.I) *
        assembledCoordinate h lam x i := by
  unfold assembledCoordinate assembledFlow
  rw [modeFlow_of_pos h (lam i) t (hwin i).1, modeCoordinate_rotFlow h (lam i) hh (hwin i).1 (hwin i).2]

omit [Fintype ι] in
/-- **(5) Explicit diagonal generator `θ_i / h`.** -/
theorem assembledCoordinate_hasDerivAt (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (x : ι → Fin 2 → ℝ) (i : ι) (t : ℝ) :
    HasDerivAt (fun s : ℝ ↦ assembledCoordinate h lam (assembledFlow h lam s x) i)
      (Complex.I * (modeAngle h (lam i) / h : ℝ) *
        assembledCoordinate h lam (assembledFlow h lam t x) i) t := by
  have hfun : (fun s : ℝ ↦ assembledCoordinate h lam (assembledFlow h lam s x) i) =
      fun s : ℝ ↦ modeCoordinate h (lam i) ((rotFlow h (lam i) s).mulVec (x i)) := by
    funext s
    unfold assembledCoordinate assembledFlow
    rw [modeFlow_of_pos h (lam i) s (hwin i).1]
  rw [hfun]
  unfold assembledCoordinate assembledFlow
  rw [modeFlow_of_pos h (lam i) t (hwin i).1]
  exact modeCoordinate_flow_hasDerivAt h (lam i) hh (hwin i).1 (hwin i).2 (x i) t

/-- The assembled sesquilinear form. -/
def assembledHermitian (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x y : ι → Fin 2 → ℝ) : ℂ :=
  ∑ i, ((realSeamEnergy (v i) / 2 : ℝ) : ℂ) * modeHermitian h (lam i) (x i) (y i)

theorem assembledHermitian_eq_coordinate (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (x y : ι → Fin 2 → ℝ) :
    assembledHermitian h lam v x y = ∑ i,
      ((realSeamEnergy (v i) / 2 : ℝ) : ℂ) * (Real.sin (modeAngle h (lam i)) : ℂ) ^ 2 *
        (starRingEnd ℂ) (assembledCoordinate h lam x i) * assembledCoordinate h lam y i := by
  unfold assembledHermitian assembledCoordinate
  exact Finset.sum_congr rfl fun i _ ↦ by
    rw [modeHermitian_eq_coordinate h (lam i) hh (hwin i).1 (hwin i).2]
    ring

theorem assembledHermitian_self (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (x : ι → Fin 2 → ℝ) :
    assembledHermitian h lam v x x = (assembledEnergy h lam v x : ℂ) := by
  unfold assembledHermitian assembledEnergy
  push_cast
  exact Finset.sum_congr rfl fun i _ ↦ by
    rw [modeHermitian_self h (lam i) hh (hwin i).1 (hwin i).2]

/-- **(5) Unitarity of the assembled flow.** -/
theorem assembledHermitian_flow (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4) (t : ℝ) (x y : ι → Fin 2 → ℝ) :
    assembledHermitian h lam v (assembledFlow h lam t x) (assembledFlow h lam t y) =
      assembledHermitian h lam v x y := by
  unfold assembledHermitian assembledFlow
  exact Finset.sum_congr rfl fun i _ ↦ by
    rw [modeFlow_of_pos h (lam i) t (hwin i).1,
      modeHermitian_rotFlow h (lam i) hh (hwin i).1 (hwin i).2]

/-- Pairing of a combination of pairwise orthogonal seam vectors with one member. -/
theorem seamInner_sum_single (v : ι → Fin 30 → ℝ)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) (α : ι → ℝ) (j : ι) :
    realSeamInner (∑ i, α i • v i) (v j) = α j * realSeamEnergy (v j) := by
  rw [realSeamInner_sum_left, Finset.sum_eq_single j]
  · rw [seamInner_smul_left, realSeamInner_self_eq_energy]
  · intro i _ hij
    rw [seamInner_smul_left, horth i j hij, mul_zero]
  · intro hj
    exact absurd (Finset.mem_univ j) hj

/-- **(5) The coefficient state is recovered from the seam potential and the
seam electric field** when the family is nonzero and pairwise orthogonal, so
the assembled inner product is faithful on the field span. -/
theorem coefficients_of_fields (v : ι → Fin 30 → ℝ) (hv : ∀ i, v i ≠ 0)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) (x y : ι → Fin 2 → ℝ)
    (hpot : potentialOf v x = potentialOf v y) (hel : electricOf v x = electricOf v y) : x = y := by
  funext j k
  have hE := (realSeamEnergy_pos_of_ne_zero (v j) (hv j)).ne'
  fin_cases k
  · have h1 := congrArg (fun w ↦ realSeamInner w (v j)) hpot
    simp only [potentialOf] at h1
    rw [seamInner_sum_single v horth, seamInner_sum_single v horth] at h1
    exact mul_right_cancel₀ hE h1
  · have h1 := congrArg (fun w ↦ realSeamInner w (v j)) hel
    simp only [electricOf] at h1
    rw [realSeamInner_comm, seamInner_neg_right, realSeamInner_comm (-_), seamInner_neg_right,
      realSeamInner_comm, seamInner_sum_single v horth, realSeamInner_comm,
      seamInner_sum_single v horth, neg_inj] at h1
    exact mul_right_cancel₀ hE h1

/-- **(5) The Hilbert reading on a finite orthogonal family of curl modes.**
For a finite family of pairwise orthogonal nonzero eigenvectors with every
eigenvalue inside the strict window (hence a family of curl modes; the
statement is about the stated family, not about the whole curl sector): the assembled form is a positive definite symmetric bilinear
form preserved by the assembled flow; its diagonal is the committed staggered
form of the generated history at every step; the componentwise complex
coordinates evolve by the diagonal phases `exp (i θ_i t / h)` with derivative
`i (θ_i / h)`; the sesquilinear form is preserved; and the coefficient state
is recovered from the seam potential and electric field.  The gradient
sector is excluded by `gradient_direction_null` and `assembledInner_radical`. -/
theorem orthogonal_family_hilbert_reading (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hwin : ∀ i, 0 < h ^ 2 * lam i ∧ h ^ 2 * lam i < 4)
    (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i) (hnz : ∀ i, v i ≠ 0)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) :
    (∀ x, x ≠ 0 → 0 < assembledInner h lam v x x) ∧
    (∀ x y, assembledInner h lam v x y = assembledInner h lam v y x) ∧
    (∀ x y z, assembledInner h lam v (x + y) z =
      assembledInner h lam v x z + assembledInner h lam v y z) ∧
    (∀ (c : ℝ) x y, assembledInner h lam v (c • x) y = c * assembledInner h lam v x y) ∧
    (∀ t x y, assembledInner h lam v (assembledFlow h lam t x) (assembledFlow h lam t y) =
      assembledInner h lam v x y) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h lam v x) (fun _ ↦ 0) n =
      assembledInner h lam v x x) ∧
    (∀ t x i, assembledCoordinate h lam (assembledFlow h lam t x) i =
      Complex.exp ((modeAngle h (lam i) * t / h : ℝ) * Complex.I) *
        assembledCoordinate h lam x i) ∧
    (∀ x i t, HasDerivAt (fun s : ℝ ↦ assembledCoordinate h lam (assembledFlow h lam s x) i)
      (Complex.I * (modeAngle h (lam i) / h : ℝ) *
        assembledCoordinate h lam (assembledFlow h lam t x) i) t) ∧
    (∀ t x y, assembledHermitian h lam v (assembledFlow h lam t x) (assembledFlow h lam t y) =
      assembledHermitian h lam v x y) ∧
    (∀ x, assembledHermitian h lam v x x = (assembledInner h lam v x x : ℂ)) ∧
    (∀ x y, potentialOf v x = potentialOf v y → electricOf v x = electricOf v y → x = y) := by
  have hadm : ∀ i, Admissible h (lam i) := fun i ↦ Or.inr (hwin i)
  refine ⟨assembledInner_pos h lam v hwin hnz, assembledInner_comm h lam v,
    assembledInner_add_left h lam v, assembledInner_smul_left h lam v,
    assembledInner_flow h lam v hh hadm, fun x n ↦ ?_, assembledCoordinate_flow h lam hh hwin,
    assembledCoordinate_hasDerivAt h lam hh hwin, assembledHermitian_flow h lam v hh hwin,
    fun x ↦ ?_, coefficients_of_fields v hnz horth⟩
  · rw [assembledInner_self]
    exact assembledHistory_energy h lam v hh hadm hv horth x n
  · rw [assembledHermitian_self h lam v hh hwin, assembledInner_self]

end Assembled

/-! ### The four committed curl modes -/

/-- The four committed curl eigenvectors `twoMode, threeMode, fiveMode, goldenMode`. -/
def curlFamily : Fin 4 → Fin 30 → ℝ := fun i ↦ fiveFamily 0 (Fin.castSucc i)

/-- Their eigenvalues `2, 3, 5, 3 + √5`. -/
def curlLam : Fin 4 → ℝ := fun i ↦ fiveLam (Fin.castSucc i)

theorem curlFamily_eigen (i : Fin 4) :
    localMaxwellOperator (curlFamily i) = curlLam i • curlFamily i :=
  fiveFamily_eigen 0 (Fin.castSucc i)

theorem curlFamily_orth (i j : Fin 4) (hij : i ≠ j) :
    realSeamInner (curlFamily i) (curlFamily j) = 0 :=
  fiveFamily_orth 0 _ _ (fun hc ↦ hij (Fin.castSucc_injective 4 hc))

theorem curlFamily_ne_zero (i : Fin 4) : curlFamily i ≠ 0 := by
  fin_cases i
  · exact twoMode_ne_zero
  · exact threeMode_ne_zero
  · exact fiveMode_ne_zero
  · exact goldenMode_ne_zero

theorem curlLam_window (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (i : Fin 4) :
    0 < h ^ 2 * curlLam i ∧ h ^ 2 * curlLam i < 4 := by
  have hp : 0 < h ^ 2 := by positivity
  have h2 := two_lt_sqrt5
  fin_cases i
  · show 0 < h ^ 2 * 2 ∧ h ^ 2 * 2 < 4
    constructor <;> nlinarith
  · show 0 < h ^ 2 * 3 ∧ h ^ 2 * 3 < 4
    constructor <;> nlinarith
  · show 0 < h ^ 2 * 5 ∧ h ^ 2 * 5 < 4
    constructor <;> nlinarith
  · show 0 < h ^ 2 * (3 + Real.sqrt 5) ∧ h ^ 2 * (3 + Real.sqrt 5) < 4
    constructor <;> nlinarith

/-- **(5) Instance.**  The four committed curl modes at every declared `h` with
`h² (3 + √5) < 4` inhabit the curl-sector Hilbert reading. -/
theorem curl_family_hilbert_instance (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    (∀ x, x ≠ 0 → 0 < assembledInner h curlLam curlFamily x x) ∧
    (∀ t x y, assembledInner h curlLam curlFamily (assembledFlow h curlLam t x)
      (assembledFlow h curlLam t y) = assembledInner h curlLam curlFamily x y) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h curlLam curlFamily x) (fun _ ↦ 0) n =
      assembledInner h curlLam curlFamily x x) ∧
    (∀ t x i, assembledCoordinate h curlLam (assembledFlow h curlLam t x) i =
      Complex.exp ((modeAngle h (curlLam i) * t / h : ℝ) * Complex.I) *
        assembledCoordinate h curlLam x i) ∧
    (∀ x i t, HasDerivAt (fun s : ℝ ↦ assembledCoordinate h curlLam (assembledFlow h curlLam s x) i)
      (Complex.I * (modeAngle h (curlLam i) / h : ℝ) *
        assembledCoordinate h curlLam (assembledFlow h curlLam t x) i) t) := by
  obtain ⟨h1, _, _, _, h5, h6, h7, h8, _, _, _⟩ :=
    orthogonal_family_hilbert_reading h curlLam curlFamily hh (curlLam_window h hh h4) curlFamily_eigen
      curlFamily_ne_zero curlFamily_orth
  exact ⟨h1, h5, h6, h7, h8⟩

theorem sqrt5_lt_three : Real.sqrt 5 < 3 := by
  rw [Real.sqrt_lt' (by norm_num)]
  norm_num

/-- The declared value `h = 1/2` lies inside the window of all four committed modes. -/
theorem half_step_window : ((1 : ℝ) / 2) ^ 2 * (3 + Real.sqrt 5) < 4 := by
  have := sqrt5_lt_three
  nlinarith

/-- **(5) Concrete inhabitant.**  At the declared `h = 1/2` the assembled form
on the four committed curl modes is positive definite and flow-invariant. -/
theorem curl_family_hilbert_instance_half :
    (∀ x, x ≠ 0 → 0 < assembledInner (1 / 2) curlLam curlFamily x x) ∧
    (∀ t x y, assembledInner (1 / 2) curlLam curlFamily (assembledFlow (1 / 2) curlLam t x)
      (assembledFlow (1 / 2) curlLam t y) = assembledInner (1 / 2) curlLam curlFamily x y) :=
  let hi := curl_family_hilbert_instance (1 / 2) (by norm_num) half_step_window
  ⟨hi.1, hi.2.1⟩

/-! ## 6. Degenerate eigenspaces: invariance does not select the form -/

/-- The block step on a two-mode block with eigenvalues `lam, mu`. -/
def blockStep (h lam mu : ℝ) (x : Fin 2 → Fin 2 → ℝ) : Fin 2 → Fin 2 → ℝ :=
  ![(stepMatrix h lam).mulVec (x 0), (stepMatrix h mu).mulVec (x 1)]

theorem blockStep_zero (h lam mu : ℝ) (x : Fin 2 → Fin 2 → ℝ) :
    blockStep h lam mu x 0 = (stepMatrix h lam).mulVec (x 0) := rfl

theorem blockStep_one (h lam mu : ℝ) (x : Fin 2 → Fin 2 → ℝ) :
    blockStep h lam mu x 1 = (stepMatrix h mu).mulVec (x 1) := rfl

/-- The block step is the `h`-value of the assembled flow on the block. -/
theorem blockStep_eq_assembledFlow (h lam mu : ℝ) (hh : h ≠ 0) (hl : Admissible h lam)
    (hm : Admissible h mu) (x : Fin 2 → Fin 2 → ℝ) :
    blockStep h lam mu x = assembledFlow h ![lam, mu] h x := by
  funext i
  fin_cases i
  · show (stepMatrix h lam).mulVec (x 0) = (modeFlow h (![lam, mu] 0) h).mulVec (x 0)
    simp only [Fin.isValue, Matrix.cons_val_zero]
    rw [modeFlow_step h lam hh hl]
  · show (stepMatrix h mu).mulVec (x 1) = (modeFlow h mu h).mulVec (x 1)
    rw [modeFlow_step h mu hh hm]

/-- The plain direct sum `energyInner_lam (x 0) (y 0) + energyInner_mu (x 1) (y 1)`. -/
def directSumForm (h lam mu : ℝ) (x y : Fin 2 → Fin 2 → ℝ) : ℝ :=
  energyInner h lam (x 0) (y 0) + energyInner h mu (x 1) (y 1)

/-- The cross form `energyInner_lam (x 0) (y 1) + energyInner_lam (x 1) (y 0)`. -/
def crossForm (h lam : ℝ) (x y : Fin 2 → Fin 2 → ℝ) : ℝ :=
  energyInner h lam (x 0) (y 1) + energyInner h lam (x 1) (y 0)

theorem directSumForm_comm (h lam mu : ℝ) (x y : Fin 2 → Fin 2 → ℝ) :
    directSumForm h lam mu x y = directSumForm h lam mu y x := by
  unfold directSumForm
  rw [energyInner_comm, energyInner_comm h mu]

theorem crossForm_comm (h lam : ℝ) (x y : Fin 2 → Fin 2 → ℝ) :
    crossForm h lam x y = crossForm h lam y x := by
  unfold crossForm
  rw [energyInner_comm, energyInner_comm h lam (x 1)]
  ring

theorem directSumForm_blockStep (h lam mu : ℝ) (x y : Fin 2 → Fin 2 → ℝ) :
    directSumForm h lam mu (blockStep h lam mu x) (blockStep h lam mu y) =
      directSumForm h lam mu x y := by
  unfold directSumForm
  rw [blockStep_zero, blockStep_zero, blockStep_one, blockStep_one, stepMatrix_energy_isometry,
    stepMatrix_energy_isometry]

/-- **(6) On a degenerate block the cross form is invariant.** -/
theorem crossForm_blockStep_equal (h lam : ℝ) (x y : Fin 2 → Fin 2 → ℝ) :
    crossForm h lam (blockStep h lam lam x) (blockStep h lam lam y) = crossForm h lam x y := by
  unfold crossForm
  rw [blockStep_zero, blockStep_zero, blockStep_one, blockStep_one, stepMatrix_energy_isometry,
    stepMatrix_energy_isometry]

theorem energyInner_e0_e1 (h lam : ℝ) : energyInner h lam ![1, 0] ![0, 1] = h * lam / 2 := by
  unfold energyInner
  simp

theorem energyInner_e1_e1 (h lam : ℝ) : energyInner h lam ![0, 1] ![0, 1] = 1 := by
  unfold energyInner
  simp

/-- **(6) The cross form is invariant under the block step if and only if the
two eigenvalues coincide** (declared `h ≠ 0`, curl mode `lam ≠ 0`). -/
theorem crossForm_invariant_iff (h lam mu : ℝ) (hh : h ≠ 0) (hl : lam ≠ 0) :
    (∀ x y, crossForm h lam (blockStep h lam mu x) (blockStep h lam mu y) = crossForm h lam x y) ↔
      mu = lam := by
  constructor
  · intro hinv
    have h1 := hinv ![![1, 0], 0] ![0, ![0, 1]]
    unfold crossForm at h1
    rw [blockStep_zero, blockStep_zero, blockStep_one, blockStep_one] at h1
    simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.mulVec_zero] at h1
    rw [stepMatrix_mulVec_e0, stepMatrix_mulVec_e1, energyInner_e0_e1] at h1
    unfold energyInner at h1
    simp at h1
    have h2 : h ^ 3 * lam * (mu - lam) = 0 := by linear_combination 2 * h1
    have h3 : h ^ 3 * lam ≠ 0 := mul_ne_zero (pow_ne_zero 3 hh) hl
    exact sub_eq_zero.mp ((mul_eq_zero.mp h2).resolve_left h3)
  · rintro rfl
    exact crossForm_blockStep_equal h mu

/-- **(6) Inside a degenerate eigenspace the invariant symmetric forms are not
unique up to scale, and the extra freedom mixes the modes.**  Both the direct
sum and the cross form are symmetric and invariant under the block step, and
neither is a real multiple of the other.  Non-uniqueness up to scale alone is
not special to the degenerate case (independent weights on two summands are
invariant for any eigenvalues); the point of the block with equal eigenvalues
is that a mode-mixing form is invariant (`crossForm_invariant_iff`), so
invariance does not fix the orthogonal mode decomposition.  A symmetric form
whose diagonal is the committed assembled energy is `assembledInner`
(`assembled_form_of_diagonal`); the diagonal, not invariance, fixes it. -/
theorem degenerate_block_invariant_forms_not_unique (h lam : ℝ) :
    (∀ x y, directSumForm h lam lam x y = directSumForm h lam lam y x) ∧
    (∀ x y, directSumForm h lam lam (blockStep h lam lam x) (blockStep h lam lam y) =
      directSumForm h lam lam x y) ∧
    (∀ x y, crossForm h lam x y = crossForm h lam y x) ∧
    (∀ x y, crossForm h lam (blockStep h lam lam x) (blockStep h lam lam y) =
      crossForm h lam x y) ∧
    (¬ ∃ c : ℝ, ∀ x y, crossForm h lam x y = c * directSumForm h lam lam x y) ∧
    (¬ ∃ c : ℝ, ∀ x y, directSumForm h lam lam x y = c * crossForm h lam x y) := by
  refine ⟨directSumForm_comm h lam lam, directSumForm_blockStep h lam lam, crossForm_comm h lam,
    crossForm_blockStep_equal h lam, ?_, ?_⟩
  · rintro ⟨c, hc⟩
    have h1 := hc ![![0, 1], 0] ![![0, 1], 0]
    have h2 := hc ![![0, 1], ![0, 1]] ![![0, 1], ![0, 1]]
    unfold crossForm directSumForm at h1 h2
    simp [energyInner] at h1 h2
    have hc0 : c = 0 := by linarith
    rw [hc0] at h2
    norm_num at h2
  · rintro ⟨c, hc⟩
    have h1 := hc ![![0, 1], 0] ![![0, 1], 0]
    unfold crossForm directSumForm at h1
    simp [energyInner] at h1


end

end OPH.FieldSectorEnergyInnerProduct

#print axioms OPH.FieldSectorEnergyInnerProduct.modeForm_posDef_iff
#print axioms OPH.FieldSectorEnergyInnerProduct.modeForm_boundary_zero
#print axioms OPH.FieldSectorEnergyInnerProduct.modeForm_boundary_rank_one
#print axioms OPH.FieldSectorEnergyInnerProduct.modeForm_exterior_neg
#print axioms OPH.FieldSectorEnergyInnerProduct.modeForm_nonpos_of_lam_nonpos
#print axioms OPH.FieldSectorEnergyInnerProduct.energyInner_polarization
#print axioms OPH.FieldSectorEnergyInnerProduct.energyInnerCore
#print axioms OPH.FieldSectorEnergyInnerProduct.stepMatrix_energy_isometry
#print axioms OPH.FieldSectorEnergyInnerProduct.modeFlow_energy_isometry
#print axioms OPH.FieldSectorEnergyInnerProduct.modeFlow_isometry_group
#print axioms OPH.FieldSectorEnergyInnerProduct.invariant_symmetric_form_eq_smul
#print axioms OPH.FieldSectorEnergyInnerProduct.invariant_symmetric_forms_iff
#print axioms OPH.FieldSectorEnergyInnerProduct.invariant_inner_product_unique_up_to_positive_scale
#print axioms OPH.FieldSectorEnergyInnerProduct.staggered_energy_fixes_scale
#print axioms OPH.FieldSectorEnergyInnerProduct.staggered_energy_inner_product_iff
#print axioms OPH.FieldSectorEnergyInnerProduct.complexStructure_sq
#print axioms OPH.FieldSectorEnergyInnerProduct.complexStructure_eq_conj
#print axioms OPH.FieldSectorEnergyInnerProduct.rotFlow_eq_cos_add_sin
#print axioms OPH.FieldSectorEnergyInnerProduct.complexStructure_energy_isometry
#print axioms OPH.FieldSectorEnergyInnerProduct.complexStructure_skew
#print axioms OPH.FieldSectorEnergyInnerProduct.modeCoordinate_bijective
#print axioms OPH.FieldSectorEnergyInnerProduct.modeCoordinate_complexStructure
#print axioms OPH.FieldSectorEnergyInnerProduct.modeCoordinate_rotFlow
#print axioms OPH.FieldSectorEnergyInnerProduct.modeHermitian_eq_coordinate
#print axioms OPH.FieldSectorEnergyInnerProduct.modeHermitian_rotFlow
#print axioms OPH.FieldSectorEnergyInnerProduct.modeCoordinate_flow_hasDerivAt
#print axioms OPH.FieldSectorEnergyInnerProduct.mode_unitary_group_explicit_generator
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledInner_self
#print axioms OPH.FieldSectorEnergyInnerProduct.assembled_form_of_diagonal
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledInner_flow
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledInner_pos
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledInnerCore
#print axioms OPH.FieldSectorEnergyInnerProduct.gradient_direction_null
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledInner_radical
#print axioms OPH.FieldSectorEnergyInnerProduct.gradient_sector_not_definite
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledCoordinate_flow
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledCoordinate_hasDerivAt
#print axioms OPH.FieldSectorEnergyInnerProduct.assembledHermitian_flow
#print axioms OPH.FieldSectorEnergyInnerProduct.coefficients_of_fields
#print axioms OPH.FieldSectorEnergyInnerProduct.orthogonal_family_hilbert_reading
#print axioms OPH.FieldSectorEnergyInnerProduct.curl_family_hilbert_instance
#print axioms OPH.FieldSectorEnergyInnerProduct.curl_family_hilbert_instance_half
#print axioms OPH.FieldSectorEnergyInnerProduct.blockStep_eq_assembledFlow
#print axioms OPH.FieldSectorEnergyInnerProduct.crossForm_invariant_iff
#print axioms OPH.FieldSectorEnergyInnerProduct.degenerate_block_invariant_forms_not_unique
