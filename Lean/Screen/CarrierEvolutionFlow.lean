import CarrierModeOscillators

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.CarrierEvolutionFlow

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierModeOscillators

/-!
# A continuous one-parameter flow whose step is the scaled Maxwell evolution

STATUS.  Exact finite real linear algebra on `Fin 2`, `Fin 30`, `Fin 20`,
elementary trigonometry, and the committed integer eigenvectors of
`ScaledMaxwellStability` and `CarrierModeOscillators`.  Every history below
is a solution of the committed zero-current scaled Ampere evolution
`AmpereEvolutionScaled h A 0 0` in the temporal gauge `φ = 0`.  The step
`h` is a declared Courant number and the real parameter `t` is a declared
flow parameter.  No register row is discharged.

WHAT IS PROVED.
1. Per-mode state map.  For a history `A n = a n • v` on an eigenvector
   `CᵀC v = lam • v`, `v ≠ 0`, with velocity-like variable
   `b n = (a (n+1) - a n) / h`, the zero-current temporal-gauge evolution
   is equivalent to the state map `(a, b) ↦ (a + h b, b - h lam (a + h b))`
   at every step (`ampere_iff_stepMatrix`, through the scalar recurrence
   `ampere_iff_recurrence`); the map is the matrix
   `stepMatrix h lam = !![1, h; -h lam, 1 - h² lam]` with determinant one
   (`stepMatrix_det`) and trace `2 - h² lam` (`stepMatrix_trace`).  For
   `0 < h² lam < 4` the explicit matrix
   `conjugator h lam = !![h, 0; -h² lam / 2, -sin θ]`, `θ = modeAngle h lam`,
   with explicit inverse `conjugatorInv`, conjugates the rotation
   `rot θ = !![cos θ, -sin θ; sin θ, cos θ]` to the step
   (`stepMatrix_eq_conj`, `conjugatorInv_eq_inv`).
2. Continuous flow per mode.  `rotFlow h lam t = Q · rot (θ t / h) · Q⁻¹`
   obeys `rotFlow_zero`, the group law `rotFlow_add`, `rotFlow_step`
   (`rotFlow h lam h = stepMatrix h lam`), `rotFlow_nat_step`
   (`rotFlow h lam (n h) = stepMatrix h lam ^ n`), and is continuous in `t`
   (`rotFlow_continuous`).  The pulled-back Euclidean form
   `pulledForm h lam x = ‖Q⁻¹ x‖²` is preserved (`pulledForm_rotFlow`) and
   equals `modeForm h lam x / sin² θ` with
   `modeForm h lam (a, b) = lam a² + h lam a b + b²`
   (`pulledForm_eq_modeForm`), while the committed staggered form of the
   mode history is `(‖v‖² / 2) modeForm` (`fieldEnergyScaled_scalarHistory`);
   the flow conserves the committed energy (`modeForm_rotFlow`).
3. Gradient sector.  On the kernel of `CᵀC` the step is the shear
   `!![1, h; 0, 1]`, the `h`-step of the shear flow `!![1, t; 0, 1]`, a
   one-parameter group preserving `modeForm h 0 = b²`
   (`shearFlow_add`, `shearFlow_step`, `modeForm_shearFlow`).
4. Assembled flow.  `modeFlow h lam t` is the shear flow for `lam = 0` and
   the rotation conjugate otherwise; for a finite family of pairwise
   orthogonal eigenvectors `v i` with `lam i = 0 ∨ 0 < h² lam i < 4`, the
   componentwise coefficient flow `assembledFlow` on states `ι → Fin 2 → ℝ` is a
   one-parameter group, linear and continuous in `t`, whose `n h` values
   generate the history `A n = ∑ i, (Φ (n h) x) i 0 • v i` solving the
   committed evolution with electric field `-∑ i, (Φ (n h) x) i 1 • v i`
   and with committed staggered form `∑ i, (‖v i‖² / 2) modeForm (x i)` at
   every step, the same value the flow preserves at every real `t`
   (`assembled_flow`).  Eigenvectors of `CᵀC` with distinct eigenvalues are
   orthogonal (`eigen_orthogonal`), so the four committed modes
   `twoMode, threeMode, fiveMode, goldenMode` together with any gradient
   seam field `d χ` instantiate the family (`carrier_flow_five_modes`).
5. Stone reading at scope (`carrier_flow_stone_reading`): on the stated
   coefficient state space there is a one-parameter group of linear maps, continuous in
   `t`, preserving the committed energy, whose `h`-step is the committed
   evolution.  Not claimed: unitarity on a Hilbert space (no inner product
   on the state space is declared here beyond the committed seam pairing);
   the corpus Stone representation (`Dynamics/StoneConverse.lean`,
   `EventAlgebra/SchroedingerFrameFlow.lean`,
   `EventAlgebra/QuantumAdequacySurface.lean`, `unique_continuous_flow`,
   the OL-C2 row surface) starts from a supplied continuous
   star-automorphism flow on the private algebra, and this module supplies
   a flow on the field sector's stated subspace, with the identification
   of `t` with physical time left to the source clock and duration row.
   The theorem does not assume the listed vectors are nonzero, so the
   coefficient-to-field map need not be injective; a faithful flow on the
   actual field span requires nonzero linear independence or quotienting the
   redundant coefficient directions.
   Extension to the whole curl sector: the six face projectors
   `projTwoR, projThreeR, projFiveR` (`ScaledMaxwellStability`) and
   `goldenPlusR, goldenMinusR` (`GoldenSectorCharacters`, with
   `normalR_plus`, `normalR_minus`: face eigenvalues `3 + √5`, `3 - √5`,
   integer tables over `ℤ[φ]`) resolve the face space, and
   `codifferential_eigen` transports their columns to seam eigenvectors of
   `CᵀC`; what is missing is the selection, from those columns through
   `faceCodifferential`, of nineteen pairwise orthogonal seam vectors with
   their `CᵀC` eigen receipts.  The assembled theorem is stated for any
   such family, so the extension is that selection.

PRIOR WORK.  `ScaledMaxwellStability` (scaled packet, energy conservation,
Courant window, `fiveMode`, `goldenMode`, face projectors);
`CarrierModeOscillators` (`modeAngle`, `cosHistory`, `sinHistory`,
`scalarHistory_ampere`, `cos_recurrence`, `twoMode`, `threeMode`,
`codifferential_eigen`); `CarrierModeEquivariance` (`eigen_pull`, the
icosahedral transport of eigenvectors); `GoldenSectorCharacters`
(`goldenPlusR`, `goldenMinusR`, `normalR_plus`, `normalR_minus`: the
golden face sector split into the `3 + √5` and `3 - √5` eigenspaces with
integer tables, idempotent, orthogonal, trace `3` each);
`Dynamics/StoneConverse.lean` (`stonePropagator`,
`finiteStoneConverse_of_continuous`);
`EventAlgebra/SchroedingerFrameFlow.lean` (`FlowedPublicFrame`,
`schroedinger_frame_duality`); `EventAlgebra/QuantumAdequacySurface.lean`
(`unique_continuous_flow`, the OL-C2 row surface: a supplied continuous
star-automorphism group on a matrix block has a self-adjoint generator).
`LocalFaceMaxwellAction` (`ker_faceCurvature_eq_gradient`,
`ker_localMaxwellOperator`: the kernel of `CᵀC` is exactly the gradient
space, harmonic part zero).  The cosine and sine histories of
`CarrierModeOscillators` are the `n h` values of `rotFlow` on the initial
states `(1, (cos θ - 1) / h)` and `(0, sin θ / h)`; the flow here
interpolates them to real `t`.

ROWS TOUCHED (none discharged).  Quantum dynamics row: the flow is
constructed on the stated coefficient state space; its physical time identification is
declared, not proved.  Light-signal row: no propagation speed is attached.
Source clock and duration row: `h` and `t` are declared parameters with no
unit.  Physical spacetime attachment row: the carrier is the committed
combinatorial complex.  Coupled-action row: the kinetic term
`(h/2) ‖E‖²` is declared in `ScaledMaxwellStability`.  Laboratory clock
and energy calibration import: no tick is used.  Gravitation-route energy
identification: the staggered form is the committed field energy, with no
mass attached.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): the kinetic term
is a declared enrichment, and the flow inherits it.

CONVENTIONS.  Forward differences: `E n = -(h⁻¹ • (A (n+1) - A n))` in the
temporal gauge, `B n = C (A n)`.  The state of a mode is `![a, b]` with
`a` the amplitude and `b = (a (n+1) - a n) / h`, so `E n = -(b n) • v`.
The rotation `rot α = !![cos α, -sin α; sin α, cos α]`; the conjugator
sends the second basis vector to `-sin θ` times the second axis, which fixes
the orientation of the rotation.  Angles are in radians per step.

FALSIFIER.  A step matrix with determinant other than one, a conjugation
`Q · rot θ · Q⁻¹` differing from the step in some entry, a flow value at
`t = h` differing from the step, or a mode history whose committed
staggered form is not `(‖v‖² / 2) (lam a² + h lam a b + b²)` would make the
module wrong.

Axiom audit.  The `#print axioms` lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## 1. The per-mode step map, derived from the committed evolution -/

/-- Velocity-like variable of a scalar profile: `b n = (a (n+1) - a n) / h`. -/
def velocity (h : ℝ) (a : ℕ → ℝ) (n : ℕ) : ℝ := (a (n + 1) - a n) / h

/-- The per-mode step matrix `!![1, h; -h lam, 1 - h² lam]`. -/
def stepMatrix (h lam : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![1, h; -(h * lam), 1 - h ^ 2 * lam]

theorem mulVec_fin_two (p q r s : ℝ) (x : Fin 2 → ℝ) :
    (!![p, q; r, s]).mulVec x = ![p * x 0 + q * x 1, r * x 0 + s * x 1] := by
  ext i
  fin_cases i
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]

theorem stepMatrix_mulVec (h lam : ℝ) (x : Fin 2 → ℝ) :
    (stepMatrix h lam).mulVec x =
      ![x 0 + h * x 1, x 1 - h * lam * (x 0 + h * x 1)] := by
  unfold stepMatrix
  rw [mulVec_fin_two]
  ext i
  fin_cases i
  · simp
  · simp
    ring

/-- The step is symplectic: determinant one. -/
theorem stepMatrix_det (h lam : ℝ) : (stepMatrix h lam).det = 1 := by
  unfold stepMatrix
  rw [Matrix.det_fin_two_of]
  ring

/-- Trace `2 - h² lam`. -/
theorem stepMatrix_trace (h lam : ℝ) : (stepMatrix h lam).trace = 2 - h ^ 2 * lam := by
  unfold stepMatrix
  rw [Matrix.trace_fin_two_of]
  ring

theorem exists_ne_zero_of_ne_zero (v : Fin 30 → ℝ) (hv : v ≠ 0) : ∃ e, v e ≠ 0 := by
  by_contra hc
  push Not at hc
  exact hv (funext hc)

/-- **Derivation of the scalar recurrence.**  For an eigenvector `v ≠ 0` and
`h ≠ 0`, the history `a n • v` solves the zero-current temporal-gauge
evolution if and only if the amplitude obeys
`a (n+2) - 2 a (n+1) + a n + h² lam a (n+1) = 0` at every step. -/
theorem ampere_iff_recurrence (h : ℝ) (hh : h ≠ 0) (a : ℕ → ℝ) (v : Fin 30 → ℝ)
    (hv0 : v ≠ 0) (lam : ℝ) (hv : localMaxwellOperator v = lam • v) :
    AmpereEvolutionScaled h (scalarHistory a v) (fun _ ↦ 0) (fun _ ↦ 0) ↔
      ∀ n, a (n + 2) - 2 * a (n + 1) + a n + h ^ 2 * lam * a (n + 1) = 0 := by
  constructor
  · intro hA n
    obtain ⟨e, he⟩ := exists_ne_zero_of_ne_zero v hv0
    have hn := hA n
    have hlm : faceCodifferential (magneticField (scalarHistory a v) (n + 1)) =
        a (n + 1) • (lam • v) := by
      show localMaxwellOperator (a (n + 1) • v) = _
      rw [map_smul, hv]
    rw [hlm, scalarHistory_electricField, scalarHistory_electricField] at hn
    have hne := congrFun hn e
    simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul, sub_zero] at hne
    rw [show n + 1 + 1 = n + 2 from rfl] at hne
    have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
    have key : (a (n + 2) - 2 * a (n + 1) + a n + h ^ 2 * lam * a (n + 1)) * v e = 0 := by
      have h1 : h⁻¹ * ((a (n + 2) - 2 * a (n + 1) + a n + h ^ 2 * lam * a (n + 1)) * v e)
          = 0 := by
        linear_combination (-1 : ℝ) * hne + (h * lam * a (n + 1) * v e) * hinv
      rcases mul_eq_zero.mp h1 with h2 | h2
      · exact absurd h2 (inv_ne_zero hh)
      · exact h2
    rcases mul_eq_zero.mp key with h2 | h2
    · exact h2
    · exact absurd h2 he
  · intro hc
    exact scalarHistory_ampere h hh a v lam hv hc

/-- The scalar recurrence is the step map on the state `![a n, b n]`. -/
theorem recurrence_iff_stepMatrix (h : ℝ) (hh : h ≠ 0) (a : ℕ → ℝ) (lam : ℝ) :
    (∀ n, a (n + 2) - 2 * a (n + 1) + a n + h ^ 2 * lam * a (n + 1) = 0) ↔
      ∀ n, ![a (n + 1), velocity h a (n + 1)] =
        (stepMatrix h lam).mulVec ![a n, velocity h a n] := by
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  constructor
  · intro hc n
    rw [stepMatrix_mulVec]
    have hr := hc n
    ext i; fin_cases i
    · simp only [velocity, Fin.zero_eta, Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
      field_simp
      ring
    · simp only [velocity, Fin.mk_one, Fin.isValue, Matrix.cons_val_one,
        Matrix.cons_val_zero, div_eq_mul_inv]
      rw [show n + 1 + 1 = n + 2 from rfl]
      linear_combination h⁻¹ * hr - (h * lam * a n) * hinv
  · intro hs n
    have h1 := congrFun (hs n) 1
    rw [stepMatrix_mulVec] at h1
    simp only [velocity, Fin.isValue, Matrix.cons_val_one,
      Matrix.cons_val_zero, div_eq_mul_inv] at h1
    rw [show n + 1 + 1 = n + 2 from rfl] at h1
    linear_combination h * h1 -
      (a (n + 2) - 2 * a (n + 1) + a n + h ^ 2 * lam * (a (n + 1) - a n)) * hinv

/-- **(1) Per-mode state map.**  On an eigenvector `v ≠ 0` the committed
zero-current temporal-gauge evolution of `a n • v` is the state map
`stepMatrix h lam` on `![a n, (a (n+1) - a n) / h]`. -/
theorem ampere_iff_stepMatrix (h : ℝ) (hh : h ≠ 0) (a : ℕ → ℝ) (v : Fin 30 → ℝ)
    (hv0 : v ≠ 0) (lam : ℝ) (hv : localMaxwellOperator v = lam • v) :
    AmpereEvolutionScaled h (scalarHistory a v) (fun _ ↦ 0) (fun _ ↦ 0) ↔
      ∀ n, ![a (n + 1), velocity h a (n + 1)] =
        (stepMatrix h lam).mulVec ![a n, velocity h a n] :=
  (ampere_iff_recurrence h hh a v hv0 lam hv).trans (recurrence_iff_stepMatrix h hh a lam)

/-- The electric field of a scalar history is minus the velocity times `v`. -/
theorem electricField_scalarHistory_velocity (h : ℝ) (a : ℕ → ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    electricFieldScaled h (scalarHistory a v) (fun _ ↦ 0) n = -(velocity h a n) • v := by
  rw [scalarHistory_electricField]
  unfold velocity
  congr 1
  ring


/-! ## 2. Rotation, conjugator, and the continuous per-mode flow -/

/-- Rotation matrix `!![cos α, -sin α; sin α, cos α]`. -/
def rot (α : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![Real.cos α, -Real.sin α; Real.sin α, Real.cos α]

theorem rot_zero : rot 0 = 1 := by
  unfold rot
  ext i j; fin_cases i <;> fin_cases j <;> simp

theorem rot_mul (α β : ℝ) : rot α * rot β = rot (α + β) := by
  unfold rot
  rw [Matrix.mul_fin_two, Real.cos_add, Real.sin_add]
  ext i j; fin_cases i <;> fin_cases j <;> simp <;> ring

theorem rot_continuous : Continuous rot := by
  unfold rot
  refine continuous_matrix fun i j => ?_
  fin_cases i <;> fin_cases j <;> simp <;> fun_prop

/-- Euclidean square on `Fin 2`. -/
def sq2 (y : Fin 2 → ℝ) : ℝ := y 0 ^ 2 + y 1 ^ 2

theorem sq2_rot (α : ℝ) (y : Fin 2 → ℝ) : sq2 ((rot α).mulVec y) = sq2 y := by
  unfold rot sq2
  rw [mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  have := Real.sin_sq_add_cos_sq α
  linear_combination (y 0 ^ 2 + y 1 ^ 2) * this

theorem sin_modeAngle_pos (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4) :
    0 < Real.sin (modeAngle h lam) := by
  apply Real.sin_pos_of_pos_of_lt_pi
  · unfold modeAngle
    rw [Real.arccos_pos]
    linarith
  · exact modeAngle_lt_pi h lam h4

theorem sin_sq_modeAngle (h lam : ℝ) (h0 : 0 < h ^ 2 * lam) (h4 : h ^ 2 * lam < 4) :
    Real.sin (modeAngle h lam) ^ 2 = 1 - (1 - h ^ 2 * lam / 2) ^ 2 := by
  rw [Real.sin_sq, cos_modeAngle h lam h0.le h4.le]

/-- The conjugator `!![h, 0; -h² lam / 2, -sin θ]`, `θ = modeAngle h lam`. -/
def conjugator (h lam : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![h, 0; -(h ^ 2 * lam / 2), -Real.sin (modeAngle h lam)]

/-- Explicit inverse `!![1/h, 0; -(h lam / 2) / sin θ, -1 / sin θ]`. -/
def conjugatorInv (h lam : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![h⁻¹, 0; -(h * lam / 2) / Real.sin (modeAngle h lam), -(Real.sin (modeAngle h lam))⁻¹]

theorem conjugator_det (h lam : ℝ) :
    (conjugator h lam).det = -(h * Real.sin (modeAngle h lam)) := by
  unfold conjugator
  rw [Matrix.det_fin_two_of]
  ring

theorem conjugator_mul_inv (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : conjugator h lam * conjugatorInv h lam = 1 := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  unfold conjugator conjugatorInv
  rw [Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j
  all_goals simp
  all_goals field_simp
  all_goals ring

theorem conjugatorInv_mul (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : conjugatorInv h lam * conjugator h lam = 1 := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  unfold conjugator conjugatorInv
  rw [Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j
  all_goals simp
  all_goals field_simp
  all_goals ring

theorem conjugatorInv_eq_inv (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : (conjugator h lam)⁻¹ = conjugatorInv h lam :=
  Matrix.inv_eq_right_inv (conjugator_mul_inv h lam hh h0 h4)

/-- The intertwining relation `T Q = Q rot θ`. -/
theorem stepMatrix_mul_conjugator (h lam : ℝ) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) :
    stepMatrix h lam * conjugator h lam = conjugator h lam * rot (modeAngle h lam) := by
  have hc := cos_modeAngle h lam h0.le h4.le
  have hs := sin_sq_modeAngle h lam h0 h4
  unfold stepMatrix conjugator rot
  rw [Matrix.mul_fin_two, Matrix.mul_fin_two]
  ext i j; fin_cases i <;> fin_cases j <;> simp
  · linear_combination (-h) * hc
  · linear_combination (h ^ 2 * lam / 2) * hc + hs
  · linear_combination (Real.sin (modeAngle h lam)) * hc

/-- **(1) Conjugacy to the rotation.**  `T = Q · rot θ · Q⁻¹`. -/
theorem stepMatrix_eq_conj (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) :
    stepMatrix h lam = conjugator h lam * rot (modeAngle h lam) * conjugatorInv h lam := by
  rw [← stepMatrix_mul_conjugator h lam h0 h4, Matrix.mul_assoc,
    conjugator_mul_inv h lam hh h0 h4, Matrix.mul_one]

/-- **(2) The continuous per-mode flow** `rotFlow h lam t = Q · rot (θ t / h) · Q⁻¹`. -/
def rotFlow (h lam t : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  conjugator h lam * rot (modeAngle h lam * t / h) * conjugatorInv h lam

theorem rotFlow_zero (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : rotFlow h lam 0 = 1 := by
  unfold rotFlow
  rw [mul_zero, zero_div, rot_zero, Matrix.mul_one, conjugator_mul_inv h lam hh h0 h4]

theorem rotFlow_add (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (s t : ℝ) :
    rotFlow h lam s * rotFlow h lam t = rotFlow h lam (s + t) := by
  unfold rotFlow
  have hQ := conjugatorInv_mul h lam hh h0 h4
  calc conjugator h lam * rot (modeAngle h lam * s / h) * conjugatorInv h lam *
        (conjugator h lam * rot (modeAngle h lam * t / h) * conjugatorInv h lam)
      = conjugator h lam * (rot (modeAngle h lam * s / h) *
          ((conjugatorInv h lam * conjugator h lam) * rot (modeAngle h lam * t / h))) *
          conjugatorInv h lam := by
        simp only [Matrix.mul_assoc]
    _ = conjugator h lam * rot (modeAngle h lam * (s + t) / h) * conjugatorInv h lam := by
        rw [hQ, Matrix.one_mul, rot_mul]
        congr 3
        ring

theorem rotFlow_step (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) : rotFlow h lam h = stepMatrix h lam := by
  unfold rotFlow
  rw [mul_div_assoc, div_self hh, mul_one, stepMatrix_eq_conj h lam hh h0 h4]

theorem rotFlow_nat_step (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (n : ℕ) : rotFlow h lam (n * h) = stepMatrix h lam ^ n := by
  induction n with
  | zero => rw [Nat.cast_zero, zero_mul, pow_zero, rotFlow_zero h lam hh h0 h4]
  | succ n ih =>
    rw [pow_succ, ← ih, ← rotFlow_step h lam hh h0 h4, rotFlow_add h lam hh h0 h4]
    congr 1
    push_cast
    ring

theorem rotFlow_continuous (h lam : ℝ) : Continuous (rotFlow h lam) := by
  unfold rotFlow
  refine (continuous_const.matrix_mul ?_).matrix_mul continuous_const
  exact rot_continuous.comp (by fun_prop)

/-- The pulled-back Euclidean form `‖Q⁻¹ x‖²`. -/
def pulledForm (h lam : ℝ) (x : Fin 2 → ℝ) : ℝ := sq2 ((conjugatorInv h lam).mulVec x)

theorem pulledForm_rotFlow (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x : Fin 2 → ℝ) :
    pulledForm h lam ((rotFlow h lam t).mulVec x) = pulledForm h lam x := by
  unfold pulledForm rotFlow
  rw [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, ← Matrix.mul_assoc,
    conjugatorInv_mul h lam hh h0 h4, Matrix.one_mul, ← Matrix.mulVec_mulVec, sq2_rot]

/-- The mode form `lam a² + h lam a b + b²` on the state `![a, b]`. -/
def modeForm (h lam : ℝ) (x : Fin 2 → ℝ) : ℝ :=
  lam * x 0 ^ 2 + h * lam * x 0 * x 1 + x 1 ^ 2

/-- **Exact relation** `sin² θ · ‖Q⁻¹ x‖² = lam a² + h lam a b + b²`. -/
theorem sin_sq_mul_pulledForm (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    Real.sin (modeAngle h lam) ^ 2 * pulledForm h lam x = modeForm h lam x := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  have hsq := sin_sq_modeAngle h lam h0 h4
  unfold pulledForm conjugatorInv sq2 modeForm
  rw [mulVec_fin_two]
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  field_simp
  linear_combination (4 * x 0 ^ 2) * hsq

theorem pulledForm_eq_modeForm (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (x : Fin 2 → ℝ) :
    pulledForm h lam x = modeForm h lam x / Real.sin (modeAngle h lam) ^ 2 := by
  have hs := (sin_modeAngle_pos h lam h0 h4).ne'
  rw [← sin_sq_mul_pulledForm h lam hh h0 h4 x]
  field_simp

/-- The flow preserves the mode form. -/
theorem modeForm_rotFlow (h lam : ℝ) (hh : h ≠ 0) (h0 : 0 < h ^ 2 * lam)
    (h4 : h ^ 2 * lam < 4) (t : ℝ) (x : Fin 2 → ℝ) :
    modeForm h lam ((rotFlow h lam t).mulVec x) = modeForm h lam x := by
  rw [← sin_sq_mul_pulledForm h lam hh h0 h4, pulledForm_rotFlow h lam hh h0 h4,
    sin_sq_mul_pulledForm h lam hh h0 h4]


/-! ## 2b. The committed staggered form of a mode history -/

theorem faceInner_curvature_eigen (v : Fin 30 → ℝ) (lam : ℝ)
    (hv : localMaxwellOperator v = lam • v) :
    faceInner (faceCurvature v) (faceCurvature v) = lam * realSeamEnergy v := by
  rw [faceCurvature_codifferential_adjoint]
  show realSeamInner v (localMaxwellOperator v) = _
  rw [hv, seamInner_smul_right, realSeamInner_self_eq_energy]

theorem succ_eq_add_velocity (h : ℝ) (hh : h ≠ 0) (a : ℕ → ℝ) (n : ℕ) :
    a (n + 1) = a n + h * velocity h a n := by
  unfold velocity
  field_simp
  ring

/-- **Energy of a mode history.**  The committed staggered form of
`a n • v` is `(‖v‖² / 2) (lam a² + h lam a b + b²)` at the state
`(a n, b n)`. -/
theorem fieldEnergyScaled_scalarHistory (h : ℝ) (hh : h ≠ 0) (a : ℕ → ℝ)
    (v : Fin 30 → ℝ) (lam : ℝ) (hv : localMaxwellOperator v = lam • v) (n : ℕ) :
    fieldEnergyScaled h (scalarHistory a v) (fun _ ↦ 0) n =
      (realSeamEnergy v / 2) * modeForm h lam ![a n, velocity h a n] := by
  unfold fieldEnergyScaled
  rw [electricField_scalarHistory_velocity, seamEnergy_smul]
  unfold magneticField scalarHistory
  rw [map_smul, map_smul, faceInner_smul_left, faceInner_smul_right,
    faceInner_curvature_eigen v lam hv, succ_eq_add_velocity h hh a n]
  unfold modeForm
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

/-! ## 3. The gradient sector: shear flow on the kernel of `CᵀC` -/

/-- Shear flow `!![1, t; 0, 1]`. -/
def shearFlow (t : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![1, t; 0, 1]

theorem shearFlow_zero : shearFlow 0 = 1 := by
  unfold shearFlow
  ext i j; fin_cases i <;> fin_cases j <;> simp

theorem shearFlow_add (s t : ℝ) : shearFlow s * shearFlow t = shearFlow (s + t) := by
  unfold shearFlow
  rw [Matrix.mul_fin_two]
  ext i j; fin_cases i <;> fin_cases j <;> simp [add_comm]

/-- **(4)** The kernel step is the shear at `t = h`. -/
theorem shearFlow_step (h : ℝ) : shearFlow h = stepMatrix h 0 := by
  unfold shearFlow stepMatrix
  ext i j; fin_cases i <;> fin_cases j <;> simp

theorem shearFlow_nat_step (h : ℝ) (n : ℕ) : shearFlow (n * h) = stepMatrix h 0 ^ n := by
  induction n with
  | zero => rw [Nat.cast_zero, zero_mul, pow_zero, shearFlow_zero]
  | succ n ih =>
    rw [pow_succ, ← ih, ← shearFlow_step h, shearFlow_add]
    congr 1
    push_cast
    ring

theorem shearFlow_continuous : Continuous shearFlow := by
  unfold shearFlow
  refine continuous_matrix fun i j => ?_
  fin_cases i <;> fin_cases j <;> simp <;> fun_prop

theorem modeForm_shearFlow (h t : ℝ) (x : Fin 2 → ℝ) :
    modeForm h 0 ((shearFlow t).mulVec x) = modeForm h 0 x := by
  unfold shearFlow modeForm
  rw [mulVec_fin_two]
  simp

/-! ## 4. The unified per-mode flow -/

/-- Admissible mode data: the kernel, or an eigenvalue inside the open window. -/
def Admissible (h lam : ℝ) : Prop := lam = 0 ∨ (0 < h ^ 2 * lam ∧ h ^ 2 * lam < 4)

open Classical in
/-- The per-mode flow: the shear flow on the kernel, the rotation conjugate
otherwise. -/
def modeFlow (h lam t : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  if lam = 0 then shearFlow t else rotFlow h lam t

theorem modeFlow_of_zero (h t : ℝ) : modeFlow h 0 t = shearFlow t := by
  unfold modeFlow
  rw [if_pos rfl]

theorem modeFlow_of_pos (h lam t : ℝ) (h0 : 0 < h ^ 2 * lam) :
    modeFlow h lam t = rotFlow h lam t := by
  unfold modeFlow
  have hne : lam ≠ 0 := by
    rintro rfl
    simp at h0
  rw [if_neg hne]

theorem modeFlow_zero (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) :
    modeFlow h lam 0 = 1 := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · rw [modeFlow_of_zero, shearFlow_zero]
  · rw [modeFlow_of_pos h lam 0 h0, rotFlow_zero h lam hh h0 h4]

theorem modeFlow_add (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) (s t : ℝ) :
    modeFlow h lam s * modeFlow h lam t = modeFlow h lam (s + t) := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · simp only [modeFlow_of_zero, shearFlow_add]
  · simp only [modeFlow_of_pos h lam _ h0, rotFlow_add h lam hh h0 h4]

theorem modeFlow_step (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) :
    modeFlow h lam h = stepMatrix h lam := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · rw [modeFlow_of_zero, shearFlow_step]
  · rw [modeFlow_of_pos h lam _ h0, rotFlow_step h lam hh h0 h4]

theorem modeFlow_nat_step (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) (n : ℕ) :
    modeFlow h lam (n * h) = stepMatrix h lam ^ n := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · rw [modeFlow_of_zero, shearFlow_nat_step]
  · rw [modeFlow_of_pos h lam _ h0, rotFlow_nat_step h lam hh h0 h4]

theorem modeFlow_continuous (h lam : ℝ) : Continuous (modeFlow h lam) := by
  unfold modeFlow
  by_cases hl : lam = 0
  · simp only [hl, if_true]
    exact shearFlow_continuous
  · simp only [hl, if_false]
    exact rotFlow_continuous h lam

theorem modeForm_modeFlow (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) (t : ℝ)
    (x : Fin 2 → ℝ) : modeForm h lam ((modeFlow h lam t).mulVec x) = modeForm h lam x := by
  rcases hadm with rfl | ⟨h0, h4⟩
  · rw [modeFlow_of_zero, modeForm_shearFlow]
  · rw [modeFlow_of_pos h lam _ h0, modeForm_rotFlow h lam hh h0 h4]

/-- One step of the flow advances the state by the step matrix. -/
theorem modeFlow_succ_step (h lam : ℝ) (hh : h ≠ 0) (hadm : Admissible h lam) (n : ℕ)
    (x : Fin 2 → ℝ) :
    (modeFlow h lam ((n + 1 : ℕ) * h)).mulVec x =
      (stepMatrix h lam).mulVec ((modeFlow h lam (n * h)).mulVec x) := by
  rw [modeFlow_nat_step h lam hh hadm, modeFlow_nat_step h lam hh hadm, Matrix.mulVec_mulVec,
    pow_succ']


/-! ## 5. The assembled flow on a finite family of orthogonal eigenvectors -/

section Assembled

variable {ι : Type}

/-- Componentwise flow on the state `ι → Fin 2 → ℝ`. -/
def assembledFlow (h : ℝ) (lam : ι → ℝ) (t : ℝ) (x : ι → Fin 2 → ℝ) : ι → Fin 2 → ℝ :=
  fun i ↦ (modeFlow h (lam i) t).mulVec (x i)

theorem assembledFlow_zero (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) :
    assembledFlow h lam 0 x = x := by
  funext i
  unfold assembledFlow
  rw [modeFlow_zero h (lam i) hh (hadm i), Matrix.one_mulVec]

theorem assembledFlow_add (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (s t : ℝ) (x : ι → Fin 2 → ℝ) :
    assembledFlow h lam s (assembledFlow h lam t x) = assembledFlow h lam (s + t) x := by
  funext i
  unfold assembledFlow
  rw [Matrix.mulVec_mulVec, modeFlow_add h (lam i) hh (hadm i)]

theorem assembledFlow_add_state (h : ℝ) (lam : ι → ℝ) (t : ℝ) (x y : ι → Fin 2 → ℝ) :
    assembledFlow h lam t (x + y) = assembledFlow h lam t x + assembledFlow h lam t y := by
  funext i
  unfold assembledFlow
  simp only [Pi.add_apply, Matrix.mulVec_add]

theorem assembledFlow_smul_state (h : ℝ) (lam : ι → ℝ) (t c : ℝ) (x : ι → Fin 2 → ℝ) :
    assembledFlow h lam t (c • x) = c • assembledFlow h lam t x := by
  funext i
  unfold assembledFlow
  simp only [Pi.smul_apply, Matrix.mulVec_smul]

theorem assembledFlow_continuous (h : ℝ) (lam : ι → ℝ) (x : ι → Fin 2 → ℝ) :
    Continuous fun t ↦ assembledFlow h lam t x := by
  refine continuous_pi fun i ↦ ?_
  exact (modeFlow_continuous h (lam i)).matrix_mulVec continuous_const

theorem assembledFlow_step (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) (i : ι) :
    assembledFlow h lam h x i = (stepMatrix h (lam i)).mulVec (x i) := by
  unfold assembledFlow
  rw [modeFlow_step h (lam i) hh (hadm i)]

/-- Scalar profile of mode `i` along the flow. -/
def modeProfile (h : ℝ) (lam : ι → ℝ) (x : ι → Fin 2 → ℝ) (i : ι) (n : ℕ) : ℝ :=
  assembledFlow h lam (n * h) x i 0

/-- The velocity of the mode profile is the second state component. -/
theorem velocity_modeProfile (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) (i : ι) (n : ℕ) :
    velocity h (modeProfile h lam x i) n = assembledFlow h lam (n * h) x i 1 := by
  unfold velocity modeProfile assembledFlow
  rw [modeFlow_succ_step h (lam i) hh (hadm i), stepMatrix_mulVec]
  simp only [Fin.isValue, Matrix.cons_val_zero]
  field_simp
  ring

theorem state_modeProfile (h : ℝ) (lam : ι → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) (i : ι) (n : ℕ) :
    ![modeProfile h lam x i n, velocity h (modeProfile h lam x i) n] =
      assembledFlow h lam (n * h) x i := by
  rw [velocity_modeProfile h lam hh hadm]
  ext j; fin_cases j
  · rfl
  · rfl

/-- Each mode profile solves the scalar recurrence. -/
theorem modeProfile_ampere (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i)
    (x : ι → Fin 2 → ℝ) (i : ι) :
    AmpereEvolutionScaled h (scalarHistory (modeProfile h lam x i) (v i)) (fun _ ↦ 0)
      (fun _ ↦ 0) := by
  apply scalarHistory_ampere h hh _ (v i) (lam i) (hv i)
  apply (recurrence_iff_stepMatrix h hh _ (lam i)).mpr
  intro n
  rw [state_modeProfile h lam hh hadm, state_modeProfile h lam hh hadm]
  unfold assembledFlow
  rw [modeFlow_succ_step h (lam i) hh (hadm i)]

variable [Fintype ι]

/-- The seam potential of a state: `∑ i, a i • v i`. -/
def potentialOf (v : ι → Fin 30 → ℝ) (x : ι → Fin 2 → ℝ) : Fin 30 → ℝ :=
  ∑ i, x i 0 • v i

/-- The seam electric field of a state: `-∑ i, b i • v i`. -/
def electricOf (v : ι → Fin 30 → ℝ) (x : ι → Fin 2 → ℝ) : Fin 30 → ℝ :=
  -∑ i, x i 1 • v i

/-- The assembled energy `∑ i, (‖v i‖² / 2) modeForm h (lam i) (x i)`. -/
def assembledEnergy (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x : ι → Fin 2 → ℝ) : ℝ :=
  ∑ i, (realSeamEnergy (v i) / 2) * modeForm h (lam i) (x i)

/-- The history generated by the flow at the step values `n h`. -/
def assembledHistory (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (x : ι → Fin 2 → ℝ) :
    ℕ → Fin 30 → ℝ :=
  fun n ↦ potentialOf v (assembledFlow h lam (n * h) x)

theorem assembledEnergy_flow (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (t : ℝ) (x : ι → Fin 2 → ℝ) :
    assembledEnergy h lam v (assembledFlow h lam t x) = assembledEnergy h lam v x := by
  unfold assembledEnergy assembledFlow
  exact Finset.sum_congr rfl fun i _ ↦ by
    rw [modeForm_modeFlow h (lam i) hh (hadm i)]

theorem assembledHistory_eq_sum (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ)
    (x : ι → Fin 2 → ℝ) (n : ℕ) :
    assembledHistory h lam v x n = ∑ i, scalarHistory (modeProfile h lam x i) (v i) n := rfl

theorem electricFieldScaled_sum (h : ℝ) (A : ι → ℕ → Fin 30 → ℝ) (n : ℕ) :
    electricFieldScaled h (fun n ↦ ∑ i, A i n) (fun _ ↦ 0) n =
      ∑ i, electricFieldScaled h (A i) (fun _ ↦ 0) n := by
  unfold electricFieldScaled
  rw [map_zero]
  funext e
  simp only [Finset.sum_apply, Pi.neg_apply, Pi.smul_apply, Pi.sub_apply, Pi.zero_apply,
    smul_eq_mul, sub_zero]
  rw [← Finset.sum_sub_distrib, Finset.mul_sum, Finset.sum_neg_distrib]

theorem ampereScaled_sum (h : ℝ) (A : ι → ℕ → Fin 30 → ℝ)
    (hA : ∀ i, AmpereEvolutionScaled h (A i) (fun _ ↦ 0) (fun _ ↦ 0)) :
    AmpereEvolutionScaled h (fun n ↦ ∑ i, A i n) (fun _ ↦ 0) (fun _ ↦ 0) := by
  intro n
  rw [electricFieldScaled_sum, electricFieldScaled_sum, ← Finset.sum_sub_distrib]
  have hB : magneticField (fun n ↦ ∑ i, A i n) (n + 1) = ∑ i, magneticField (A i) (n + 1) := by
    unfold magneticField
    exact map_sum faceCurvature _ _
  rw [hB, map_sum, sub_zero, Finset.smul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  have := hA i n
  rwa [sub_zero] at this

/-- **(3) The assembled history solves the committed evolution.** -/
theorem assembledHistory_ampere (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i)
    (x : ι → Fin 2 → ℝ) :
    AmpereEvolutionScaled h (assembledHistory h lam v x) (fun _ ↦ 0) (fun _ ↦ 0) := by
  have : assembledHistory h lam v x = fun n ↦ ∑ i, scalarHistory (modeProfile h lam x i) (v i) n :=
    funext fun n ↦ assembledHistory_eq_sum h lam v x n
  rw [this]
  exact ampereScaled_sum h _ fun i ↦ modeProfile_ampere h lam v hh hadm hv x i

/-- The electric field of the assembled history is the electric readout of the
flowed state. -/
theorem assembledHistory_electricField (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (x : ι → Fin 2 → ℝ) (n : ℕ) :
    electricFieldScaled h (assembledHistory h lam v x) (fun _ ↦ 0) n =
      electricOf v (assembledFlow h lam (n * h) x) := by
  have : assembledHistory h lam v x = fun n ↦ ∑ i, scalarHistory (modeProfile h lam x i) (v i) n :=
    funext fun n ↦ assembledHistory_eq_sum h lam v x n
  rw [this, electricFieldScaled_sum]
  unfold electricOf
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [electricField_scalarHistory_velocity, velocity_modeProfile h lam hh hadm, neg_smul]

theorem realSeamInner_sum_left (f : ι → Fin 30 → ℝ) (g : Fin 30 → ℝ) :
    realSeamInner (∑ i, f i) g = ∑ i, realSeamInner (f i) g := by
  unfold realSeamInner
  simp only [Finset.sum_apply, Finset.sum_mul]
  exact Finset.sum_comm

theorem realSeamInner_sum_right (g : Fin 30 → ℝ) (f : ι → Fin 30 → ℝ) :
    realSeamInner g (∑ i, f i) = ∑ i, realSeamInner g (f i) := by
  unfold realSeamInner
  simp only [Finset.sum_apply, Finset.mul_sum]
  exact Finset.sum_comm

/-- Pairing of two combinations of pairwise orthogonal seam vectors. -/
theorem seamInner_sum_smul (v : ι → Fin 30 → ℝ)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) (α β : ι → ℝ) :
    realSeamInner (∑ i, α i • v i) (∑ j, β j • v j) =
      ∑ i, α i * β i * realSeamEnergy (v i) := by
  rw [realSeamInner_sum_left]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [realSeamInner_sum_right, Finset.sum_eq_single i]
  · rw [seamInner_smul_left, seamInner_smul_right, realSeamInner_self_eq_energy]
    ring
  · intro j _ hji
    rw [seamInner_smul_left, seamInner_smul_right, horth i j (Ne.symm hji), mul_zero, mul_zero]
  · intro hi
    exact absurd (Finset.mem_univ i) hi

/-- **(3) Energy of the assembled history.**  Under pairwise orthogonality
the committed staggered form of the assembled history is the assembled
energy of the initial state at every step. -/
theorem assembledHistory_energy (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) (x : ι → Fin 2 → ℝ) (n : ℕ) :
    fieldEnergyScaled h (assembledHistory h lam v x) (fun _ ↦ 0) n =
      assembledEnergy h lam v x := by
  rw [← assembledEnergy_flow h lam v hh hadm (n * h) x]
  have hsum : assembledHistory h lam v x =
      fun n ↦ ∑ i, scalarHistory (modeProfile h lam x i) (v i) n :=
    funext fun n ↦ assembledHistory_eq_sum h lam v x n
  rw [hsum]
  unfold fieldEnergyScaled
  rw [electricFieldScaled_sum]
  simp only [electricField_scalarHistory_velocity]
  rw [← realSeamInner_self_eq_energy, seamInner_sum_smul v horth]
  unfold magneticField
  rw [faceCurvature_codifferential_adjoint]
  have hL : faceCodifferential (faceCurvature (∑ i, scalarHistory (modeProfile h lam x i) (v i)
      (n + 1))) = ∑ j, (modeProfile h lam x j (n + 1) * lam j) • v j := by
    show localMaxwellOperator _ = _
    rw [map_sum]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    unfold scalarHistory
    rw [map_smul, hv j, smul_smul]
  rw [hL]
  unfold scalarHistory
  rw [seamInner_sum_smul v horth]
  unfold assembledEnergy
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  rw [← state_modeProfile h lam hh hadm x i n, succ_eq_add_velocity h hh (modeProfile h lam x i) n]
  unfold modeForm
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

/-- **(3), (4) Assembled flow on the stated subspace.**  The componentwise flow
on the span of a finite family of pairwise orthogonal admissible eigenvectors
is a one-parameter group of linear maps, continuous in `t`, whose `h`-step
is the committed step and whose `n h` values generate a solution of the
committed zero-current evolution with the committed electric field and
constant committed staggered form. -/
theorem assembled_flow (h : ℝ) (lam : ι → ℝ) (v : ι → Fin 30 → ℝ) (hh : h ≠ 0)
    (hadm : ∀ i, Admissible h (lam i)) (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) :
    (∀ x, assembledFlow h lam 0 x = x) ∧
    (∀ s t x, assembledFlow h lam s (assembledFlow h lam t x) = assembledFlow h lam (s + t) x) ∧
    (∀ t x y, assembledFlow h lam t (x + y) = assembledFlow h lam t x + assembledFlow h lam t y) ∧
    (∀ (t c : ℝ) x, assembledFlow h lam t (c • x) = c • assembledFlow h lam t x) ∧
    (∀ x, Continuous fun t ↦ assembledFlow h lam t x) ∧
    (∀ x i, assembledFlow h lam h x i = (stepMatrix h (lam i)).mulVec (x i)) ∧
    (∀ t x, assembledEnergy h lam v (assembledFlow h lam t x) = assembledEnergy h lam v x) ∧
    (∀ x, AmpereEvolutionScaled h (assembledHistory h lam v x) (fun _ ↦ 0) (fun _ ↦ 0)) ∧
    (∀ x n, electricFieldScaled h (assembledHistory h lam v x) (fun _ ↦ 0) n =
      electricOf v (assembledFlow h lam (n * h) x)) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h lam v x) (fun _ ↦ 0) n =
      assembledEnergy h lam v x) :=
  ⟨assembledFlow_zero h lam hh hadm, assembledFlow_add h lam hh hadm,
    assembledFlow_add_state h lam, assembledFlow_smul_state h lam,
    assembledFlow_continuous h lam, assembledFlow_step h lam hh hadm,
    assembledEnergy_flow h lam v hh hadm, assembledHistory_ampere h lam v hh hadm hv,
    assembledHistory_electricField h lam v hh hadm,
    assembledHistory_energy h lam v hh hadm hv horth⟩

end Assembled

/-! ## 6. Orthogonality of eigenvectors and the five committed modes -/

theorem localMaxwellOperator_symm (u w : Fin 30 → ℝ) :
    realSeamInner (localMaxwellOperator u) w = realSeamInner u (localMaxwellOperator w) := by
  show realSeamInner (faceCodifferential (faceCurvature u)) w =
    realSeamInner u (faceCodifferential (faceCurvature w))
  rw [← faceCurvature_codifferential_adjoint, realSeamInner_comm,
    ← faceCurvature_codifferential_adjoint, faceInner_comm]

/-- Eigenvectors of `CᵀC` with distinct eigenvalues are orthogonal. -/
theorem eigen_orthogonal (u w : Fin 30 → ℝ) (lam mu : ℝ)
    (hu : localMaxwellOperator u = lam • u) (hw : localMaxwellOperator w = mu • w)
    (hne : lam ≠ mu) : realSeamInner u w = 0 := by
  have hs := localMaxwellOperator_symm u w
  rw [hu, hw, seamInner_smul_left, seamInner_smul_right] at hs
  have hz : (lam - mu) * realSeamInner u w = 0 := by linear_combination hs
  rcases mul_eq_zero.mp hz with h1 | h1
  · exact absurd (sub_eq_zero.mp h1) hne
  · exact h1

theorem two_lt_sqrt5 : (2 : ℝ) < Real.sqrt 5 := by
  rw [Real.lt_sqrt (by norm_num)]
  norm_num

theorem gradient_eigen (χ : Fin 12 → ℝ) :
    localMaxwellOperator (realCoboundary χ) = (0 : ℝ) • realCoboundary χ := by
  show faceCodifferential (faceCurvature (realCoboundary χ)) = _
  rw [faceCurvature_coboundary, map_zero, zero_smul]

/-- The four committed eigenvectors and a gradient seam field. -/
def fiveFamily (χ : Fin 12 → ℝ) : Fin 5 → Fin 30 → ℝ :=
  ![twoMode, threeMode, fiveMode, goldenMode, realCoboundary χ]

/-- Their eigenvalues `2, 3, 5, 3 + √5, 0`. -/
def fiveLam : Fin 5 → ℝ := ![2, 3, 5, 3 + Real.sqrt 5, 0]

theorem fiveFamily_eigen (χ : Fin 12 → ℝ) (i : Fin 5) :
    localMaxwellOperator (fiveFamily χ i) = fiveLam i • fiveFamily χ i := by
  fin_cases i
  · exact twoMode_eigen
  · exact threeMode_eigen
  · exact fiveMode_eigen
  · exact goldenMode_eigen
  · exact gradient_eigen χ

theorem fiveLam_distinct (i j : Fin 5) (hij : i ≠ j) : fiveLam i ≠ fiveLam j := by
  have h2 := two_lt_sqrt5
  have h3 := sqrt5_lt_three
  fin_cases i <;> fin_cases j <;> simp [fiveLam] at hij ⊢ <;> intro h <;> linarith

theorem fiveFamily_orth (χ : Fin 12 → ℝ) (i j : Fin 5) (hij : i ≠ j) :
    realSeamInner (fiveFamily χ i) (fiveFamily χ j) = 0 :=
  eigen_orthogonal _ _ _ _ (fiveFamily_eigen χ i) (fiveFamily_eigen χ j) (fiveLam_distinct i j hij)

theorem fiveLam_admissible (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (i : Fin 5) :
    Admissible h (fiveLam i) := by
  have hp : 0 < h ^ 2 := by positivity
  have h2 := two_lt_sqrt5
  fin_cases i
  · right
    show 0 < h ^ 2 * 2 ∧ h ^ 2 * 2 < 4
    constructor <;> nlinarith
  · right
    show 0 < h ^ 2 * 3 ∧ h ^ 2 * 3 < 4
    constructor <;> nlinarith
  · right
    show 0 < h ^ 2 * 5 ∧ h ^ 2 * 5 < 4
    constructor <;> nlinarith
  · right
    show 0 < h ^ 2 * (3 + Real.sqrt 5) ∧ h ^ 2 * (3 + Real.sqrt 5) < 4
    constructor <;> nlinarith
  · left
    rfl

/-- **Five-mode instance.**  On the span of `twoMode, threeMode, fiveMode,
goldenMode` and a gradient `d χ`, at any declared step inside the sharp
window, the assembled flow packet holds. -/
theorem carrier_flow_five_modes (χ : Fin 12 → ℝ) (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    (∀ x, assembledFlow h fiveLam 0 x = x) ∧
    (∀ s t x, assembledFlow h fiveLam s (assembledFlow h fiveLam t x) =
      assembledFlow h fiveLam (s + t) x) ∧
    (∀ t x y, assembledFlow h fiveLam t (x + y) =
      assembledFlow h fiveLam t x + assembledFlow h fiveLam t y) ∧
    (∀ (t c : ℝ) x, assembledFlow h fiveLam t (c • x) = c • assembledFlow h fiveLam t x) ∧
    (∀ x, Continuous fun t ↦ assembledFlow h fiveLam t x) ∧
    (∀ x i, assembledFlow h fiveLam h x i = (stepMatrix h (fiveLam i)).mulVec (x i)) ∧
    (∀ t x, assembledEnergy h fiveLam (fiveFamily χ) (assembledFlow h fiveLam t x) =
      assembledEnergy h fiveLam (fiveFamily χ) x) ∧
    (∀ x, AmpereEvolutionScaled h (assembledHistory h fiveLam (fiveFamily χ) x)
      (fun _ ↦ 0) (fun _ ↦ 0)) ∧
    (∀ x n, electricFieldScaled h (assembledHistory h fiveLam (fiveFamily χ) x) (fun _ ↦ 0) n =
      electricOf (fiveFamily χ) (assembledFlow h fiveLam (n * h) x)) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h fiveLam (fiveFamily χ) x) (fun _ ↦ 0) n =
      assembledEnergy h fiveLam (fiveFamily χ) x) :=
  assembled_flow h fiveLam (fiveFamily χ) hh (fiveLam_admissible h hh h4)
    (fiveFamily_eigen χ) (fiveFamily_orth χ)


/-! ## 7. Stone reading at scope -/

/-- **(5) Stone reading at scope.**  For every finite family of pairwise
orthogonal admissible eigenvectors of `CᵀC` there is a map
`Φ : ℝ → State → State` on the coefficient state space `ι → Fin 2 → ℝ` such
that: every `Φ t` is linear; `Φ 0` is the identity; `Φ s ∘ Φ t = Φ (s + t)`;
every orbit `t ↦ Φ t x` is continuous; the assembled energy is constant along
every orbit; the `h`-value is the committed step on every component; and the
`n h`-values generate a solution of the committed zero-current evolution
whose electric field is the electric readout of the flowed state and whose
committed staggered form is the assembled energy of the initial state.

Not claimed: unitarity with respect to a Hilbert-space inner product on the
state space (none is declared here beyond the committed seam pairing that
enters `assembledEnergy`); a spectral generator on the private algebra; the
identification of `t` with physical time (source clock and duration row);
the extension to the whole curl sector (the missing input is the selection
of nineteen pairwise orthogonal seam vectors from the columns of the six
face projectors `projTwoR, projThreeR, projFiveR, goldenPlusR, goldenMinusR`
through `faceCodifferential`, with `CᵀC` eigen receipts by
`codifferential_eigen`, `normalR_plus`, `normalR_minus`; the theorem applies
to any such family).  The corpus Stone representation of
`Dynamics/StoneConverse.lean`, `EventAlgebra/SchroedingerFrameFlow.lean` and
`EventAlgebra/QuantumAdequacySurface.lean` (`unique_continuous_flow`, the
OL-C2 row surface) consumes a supplied continuous star-automorphism flow on
the private algebra; this theorem supplies a continuous linear coefficient
flow that generates field histories. Because no nonzero hypothesis is imposed
on `v`, the coefficient-to-field map need not be injective; identifying this
with a faithful flow on the actual field span needs nonzero linear independence
or a quotient of redundant coefficient directions. The two constructions are
joined by no theorem here. -/
theorem carrier_flow_stone_reading {ι : Type} [Fintype ι] (h : ℝ) (lam : ι → ℝ)
    (v : ι → Fin 30 → ℝ) (hh : h ≠ 0) (hadm : ∀ i, Admissible h (lam i))
    (hv : ∀ i, localMaxwellOperator (v i) = lam i • v i)
    (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) :
    ∃ Φ : ℝ → (ι → Fin 2 → ℝ) → (ι → Fin 2 → ℝ),
      (∀ t, IsLinearMap ℝ (Φ t)) ∧
      Φ 0 = id ∧
      (∀ s t, Φ s ∘ Φ t = Φ (s + t)) ∧
      (∀ x, Continuous fun t ↦ Φ t x) ∧
      (∀ t x, assembledEnergy h lam v (Φ t x) = assembledEnergy h lam v x) ∧
      (∀ x i, Φ h x i = (stepMatrix h (lam i)).mulVec (x i)) ∧
      (∀ x, AmpereEvolutionScaled h (fun n ↦ potentialOf v (Φ (n * h) x)) (fun _ ↦ 0)
        (fun _ ↦ 0) ∧
        (∀ n, electricFieldScaled h (fun n ↦ potentialOf v (Φ (n * h) x)) (fun _ ↦ 0) n =
          electricOf v (Φ (n * h) x)) ∧
        (∀ n, fieldEnergyScaled h (fun n ↦ potentialOf v (Φ (n * h) x)) (fun _ ↦ 0) n =
          assembledEnergy h lam v x)) := by
  obtain ⟨h0, hadd, hlin, hsmul, hcont, hstep, hen, hamp, hel, hfe⟩ :=
    assembled_flow h lam v hh hadm hv horth
  refine ⟨assembledFlow h lam, fun t ↦ ⟨hlin t, hsmul t⟩, funext h0, fun s t ↦ funext (hadd s t),
    hcont, hen, hstep, fun x ↦ ⟨hamp x, hel x, hfe x⟩⟩

end

end OPH.CarrierEvolutionFlow

#print axioms OPH.CarrierEvolutionFlow.ampere_iff_stepMatrix
#print axioms OPH.CarrierEvolutionFlow.stepMatrix_det
#print axioms OPH.CarrierEvolutionFlow.stepMatrix_trace
#print axioms OPH.CarrierEvolutionFlow.stepMatrix_eq_conj
#print axioms OPH.CarrierEvolutionFlow.conjugatorInv_eq_inv
#print axioms OPH.CarrierEvolutionFlow.rotFlow_add
#print axioms OPH.CarrierEvolutionFlow.rotFlow_nat_step
#print axioms OPH.CarrierEvolutionFlow.rotFlow_continuous
#print axioms OPH.CarrierEvolutionFlow.sin_sq_mul_pulledForm
#print axioms OPH.CarrierEvolutionFlow.fieldEnergyScaled_scalarHistory
#print axioms OPH.CarrierEvolutionFlow.modeForm_rotFlow
#print axioms OPH.CarrierEvolutionFlow.shearFlow_step
#print axioms OPH.CarrierEvolutionFlow.eigen_orthogonal
#print axioms OPH.CarrierEvolutionFlow.assembled_flow
#print axioms OPH.CarrierEvolutionFlow.carrier_flow_five_modes
#print axioms OPH.CarrierEvolutionFlow.carrier_flow_stone_reading
