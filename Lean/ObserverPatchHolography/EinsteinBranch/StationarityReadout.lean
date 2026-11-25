import Mathlib.Data.Matrix.Basis
import Mathlib.Algebra.Algebra.Subalgebra.Basic
import Mathlib.Algebra.Module.LinearMap.Defs
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Util.AssertNoSorry
import ObserverPatchHolography.EinsteinBranch.RegionalMoments

/-! Fixed-volume tests for first variations of independently calibrated instruments.
The vector space consists of admitted source tangents. Physical identification of the
three linear forms and the radius-dependent error estimates belong to the source. -/

namespace OPH.EinsteinBranch.StationarityReadout

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

/-- A volume-changing calibration accounts for every discarded volume component. -/
theorem fixed_volume_iff (vol residual : E →ₗ[ℝ] ℝ) (q : E) (hq : vol q = 1) :
    (∀ x, vol x = 0 → residual x = 0) ↔
      ∀ x, residual x = vol x * residual q := by
  constructor
  · intro h x
    have hv : vol (x - vol x • q) = 0 := by simp [hq]
    have hr := h (x - vol x • q) hv
    simpa [sub_eq_zero, smul_eq_mul] using hr
  · intro h x hx
    rw [h x, hx, zero_mul]

/-- On a complete state-coordinate space, a stationary entropy-plus-area differential
is a linear combination of trace and independently surveyed volume. The trace direction
and the volume calibration may be chosen by subtracting components before application. -/
theorem trace_volume_stationarity (tr vol residual : E →ₗ[ℝ] ℝ) (n q : E)
    (htn : tr n = 1) (hvn : vol n = 0) (htq : tr q = 0) (hvq : vol q = 1) :
    (∀ x, tr x = 0 → vol x = 0 → residual x = 0) ↔
      ∀ x, residual x = tr x * residual n + vol x * residual q := by
  constructor
  · intro h x
    let y := x - tr x • n - vol x • q
    have ht : tr y = 0 := by simp [y, htn, htq]
    have hv : vol y = 0 := by simp [y, hvn, hvq]
    have hr := h y ht hv
    simp [y, smul_eq_mul] at hr
    linarith
  · intro h x ht hv
    simp [h x, ht, hv]

/-- Calibrating on one ordinary, volume-matched transfer fixes the coefficient.
The cross product checks it on all the other volume-matched source changes. -/
theorem ordinary_calibration (vol area entropy : E →ₗ[ℝ] ℝ) (u : E)
    (hv : vol u = 0) (ha : area u ≠ 0) (η : ℝ) :
    (∀ x, vol x = 0 → entropy x + η * area x = 0) ↔
      η = -entropy u / area u ∧
      ∀ x, vol x = 0 → area u * entropy x = entropy u * area x := by
  constructor
  · intro h
    have hu := h u hv
    refine ⟨(eq_div_iff ha).2 (by linarith), ?_⟩
    intro x hx
    have hu' : entropy u = -(η * area u) := by linarith
    have hx' : entropy x = -(η * area x) := by linarith [h x hx]
    rw [hu', hx']
    ring
  · rintro ⟨rfl, h⟩ x hx
    field_simp [ha]
    linarith [h x hx]

/-- A transfer invisible to either ruler keeps its entropy residual for every gain. -/
theorem blind_transfer (area entropy : E →ₗ[ℝ] ℝ) (x : E) (h : area x = 0) (η : ℝ) :
    |entropy x + η * area x| = |entropy x| := by simp [h]

/-- Exact finite stationarity error after freezing the ordinary calibration. -/
theorem calibrated_error (k a ku au : ℝ) (ha : au ≠ 0) :
    |k + (-ku / au) * a| = |au * k - ku * a| / |au| := by
  rw [← abs_div]
  congr 1
  field_simp [ha]
  ring

open scoped BigOperators in
/-- One observer's measured inputs feed the existing Einstein residual bound directly.
Here u is an ordinary calibration, the tested pair (s,a) is held out, and the area
expansion refers to an independent curvature reconstruction on the matched ball. -/
theorem measured_residual {ι : Type*} [Fintype ι] (w t : ι → ℝ)
    (s a ku au G T Λ η ℓ osc gain variation leakage ea : ℝ)
    (hu : au ≠ 0) (hcal : η = -ku / au) (hη : 0 < η) (hℓ : 0 < ℓ)
    (hrep : |s - ∑ i, w i * t i| ≤ leakage)
    (hm : |(∑ i, w i) - 2*Real.pi*Quantitative.coefficient ℓ| ≤ gain)
    (htv : (∑ i, |w i|) ≤ variation) (hosc : ∀ i, |t i - T| ≤ osc) (ho : 0 ≤ osc)
    (ha : |a + Quantitative.coefficient ℓ*(G-Λ)| ≤ ea) :
    |G-Λ-(2*Real.pi/η)*T| ≤
      (leakage + variation*osc + gain*|T| + η*ea + |au*s-ku*a|/|au|) /
        (η*Quantitative.coefficient ℓ) := by
  apply Quantitative.residual_bound s a G T Λ η ℓ _ ea _ hη hℓ
  · exact RegionalMoments.bulk_error w t s _ T osc gain variation leakage hrep hm htv hosc ho
  · exact ha
  · rw [hcal, calibrated_error s a ku au hu]

/-- A conserved operator compatible with a sharp port commutes with its
    terminal return effect. B denotes the backward evolution. -/
theorem conserved_port_effect {R : Type*} [Monoid R] {C P U B : R}
    (port : Commute C P) (flight : Commute C U) (back : Commute C B) :
    Commute C (P * B * P * U * P) :=
  (((port.mul_right back).mul_right port).mul_right flight).mul_right port

/-- A conserved restriction commuting with a tomographically complete set of
    return effects on one parity block is scalar there. -/
theorem scalar_of_complete_returns {F n : Type*} [Field F] [Fintype n] [DecidableEq n]
    (K : Matrix n n F) (s : Set (Matrix n n F))
    (complete : Submodule.span F s = ⊤) (compatible : ∀ E ∈ s, Commute K E) :
    K ∈ Set.range (Matrix.scalar n) := by
  let z := Subalgebra.centralizer F ({K} : Set (Matrix n n F))
  have hspan : Submodule.span F s ≤ z.toSubmodule := by
    apply Submodule.span_le.mpr
    intro E hE
    change E ∈ Subalgebra.centralizer F ({K} : Set (Matrix n n F))
    rw [Subalgebra.mem_centralizer_iff]
    simpa using (compatible E hE).eq
  rw [complete] at hspan
  apply Matrix.mem_range_scalar_of_commute_single
  intro i j _
  have hm : Matrix.single i j (1 : F) ∈ z.toSubmodule := hspan (Submodule.mem_top)
  change Matrix.single i j (1 : F) ∈ Subalgebra.centralizer F ({K} : Set (Matrix n n F)) at hm
  rw [Subalgebra.mem_centralizer_iff] at hm
  exact (hm K (by simp)).symm

/-- On a product of accessible local readouts, agreement of independent estimators
    forces their values to be constant, even for nonlinear estimators. -/
theorem constant_on_joint_survey {X A B Y : Type*}
    (states : Set X) (left : X → A) (right : X → B)
    (u : Set A) (v : Set B) (base : B) (hb : base ∈ v)
    (joint : ∀ a ∈ u, ∀ b ∈ v, ∃ x ∈ states, left x = a ∧ right x = b)
    (f : A → Y) (g : B → Y) (agree : ∀ x ∈ states, f (left x) = g (right x)) :
    ∀ a ∈ u, f a = g base := by
  intro a ha
  obtain ⟨x, hx, hx₁, hx₂⟩ := joint a ha base hb
  simpa only [hx₁, hx₂] using agree x hx

/-- A pair of independently accessible range variations admits no common nonzero
    linear response. Evaluating at (a,0) and (0,b) isolates each response. -/
theorem zero_common_survey_form {X A B : Type*}
    [AddCommGroup X] [Module ℝ X] [AddCommGroup A] [Module ℝ A]
    [AddCommGroup B] [Module ℝ B] (l : X →ₗ[ℝ] A) (r : X →ₗ[ℝ] B)
    (joint : Function.Surjective (fun x => (l x, r x)))
    (f : A →ₗ[ℝ] ℝ) (g : B →ₗ[ℝ] ℝ) (agree : ∀ x, f (l x) = g (r x)) :
    f = 0 ∧ g = 0 := by
  constructor
  · ext a
    obtain ⟨x, hx⟩ := joint (a, 0)
    have hl : l x = a := congrArg Prod.fst hx
    have hr : r x = 0 := congrArg Prod.snd hx
    simpa [hl, hr] using agree x
  · ext b
    obtain ⟨x, hx⟩ := joint (0, b)
    have hl : l x = 0 := congrArg Prod.fst hx
    have hr : r x = b := congrArg Prod.snd hx
    simpa [hl, hr] using (agree x).symm

/-- The finite survey Hessian certificate bounds the geometric Taylor remainder.
Quantum Pinsker converts that bound to a gap relative to the Gibbs state. -/
theorem survey_lagrangian_gap (relative squaredNorm remainder : ℝ)
    (pinsker : squaredNorm / 2 ≤ relative)
    (geometry : remainder ≤ (36 / 100 : ℝ) * squaredNorm) :
    (28 / 100 : ℝ) * relative ≤ relative - remainder := by
  linarith

/-- Add the certified one-copy Lagrangian gaps at a common price. Volume
may pass between surveys. The joint entropy deficit includes the total
correlation, and relative entropy to the product reference has the same
correlation term. The relative-entropy identities are supplied by the state. -/
theorem survey_product_gap {ι : Type*} [Fintype ι]
    (relative objectiveGap volumeOffset : ι → ℝ) (price correlation : ℝ)
    (hcor : 0 ≤ correlation) (hvolume : ∑ i, volumeOffset i = 0)
    (hlocal : ∀ i, (28 / 100 : ℝ) * relative i ≤
      objectiveGap i + price * volumeOffset i) :
    (28 / 100 : ℝ) * (correlation + ∑ i, relative i) ≤ correlation + ∑ i, objectiveGap i := by
  have hsum := Finset.sum_le_sum (s := Finset.univ) fun i _ => hlocal i
  rw [← Finset.mul_sum, Finset.sum_add_distrib, ← Finset.mul_sum, hvolume, mul_zero, add_zero] at hsum
  linarith

/-- Two volume-matched transfers isolate an entropy response invisible to area.
A nonzero determinant also excludes every common linearized area gain. -/
theorem area_volume_blind_tangent (vol area entropy : E →ₗ[ℝ] ℝ) (p q : E)
    (hp : vol p = 0) (hq : vol q = 0) :
    let d := area q • p - area p • q
    vol d = 0 ∧ area d = 0 ∧
      entropy d = area q * entropy p - area p * entropy q := by
  dsimp
  constructor
  · simp [hp, hq]
  constructor
  · simp [smul_eq_mul, mul_comm]
  · simp [smul_eq_mul]

/-- An invertible clock probability map gives a unique linear score for each
calibrated source observable. On a fixed-marginal slice the modular entropy
differential and the independently measured material force are two such observables. -/
theorem unique_clock_score {F : Type*} [AddCommGroup F] [Module ℝ F]
    (read : E ≃ₗ[ℝ] F) (target : E →ₗ[ℝ] ℝ) :
    ∃! score : F →ₗ[ℝ] ℝ, ∀ x, score (read x) = target x := by
  refine ⟨target.comp read.symm.toLinearMap, ?_, ?_⟩
  · intro x
    simp
  · intro other ho
    ext y
    calc
      other y = other (read (read.symm y)) := by rw [read.apply_symm_apply]
      _ = target (read.symm y) := ho _
      _ = _ := rfl

/-- With nonnegative score increments, uniformly earlier and later completion
probabilities give opposite signs to the two changes in mean duration. -/
theorem opposite_duration_responses {ι : Type*} [Fintype ι]
    (a b u v : ι → ℝ) (m : ℝ) (hm : 0 < m)
    (ha : ∀ i, m ≤ a i) (hb : ∀ i, b i ≤ -m)
    (hu : ∀ i, 0 ≤ u i) (hv : ∀ i, 0 ≤ v i)
    (hsu : 0 < ∑ i, u i) (hsv : 0 < ∑ i, v i) :
    -(∑ i, u i * a i) < 0 ∧ 0 < -(∑ i, v i * b i) := by
  have h₁ : (∑ i, u i) * m ≤ ∑ i, u i * a i := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (ha i) (hu i)
  have h₂ : (∑ i, v i * b i) ≤ (∑ i, v i) * (-m) := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hb i) (hv i)
  constructor
  · exact neg_neg_of_pos (lt_of_lt_of_le (mul_pos hsu hm) h₁)
  · exact neg_pos.mpr (lt_of_le_of_lt h₂ (mul_neg_of_pos_of_neg hsv (neg_neg_of_pos hm)))

/-- A regional generalized-entropy error budget on correlated subsystems.
Here corr is the change in total correlation, and extra is the change in
the conditional entropy of all other degrees of freedom in the ball. -/
theorem regional_compensation {n : ℕ} (localGain : Fin n → ℝ)
    (s corr extra area η ε m : ℝ)
    (hs : s = (∑ i, localGain i) - corr + extra)
    (hm : ∀ i, m ≤ localGain i) (hbalance : |s + η * area| ≤ ε) :
    n * m - ε ≤ corr - extra - η * area := by
  have hsum : (∑ _i : Fin n, m) ≤ ∑ i, localGain i :=
    Finset.sum_le_sum fun i _ => hm i
  have htotal : (n : ℝ) * m ≤ ∑ i, localGain i := by
    simpa using hsum
  have hr := (abs_le.mp hbalance).2
  linarith

/-- One common clock acceleration u and one flux conversion b imply the
three independently calibrated odd optical equations have zero determinant. -/
theorem common_clock_flux_determinant (g₀ g₁ g₂ h₀ h₁ h₂ j₀ j₁ j₂ u b : ℝ)
    (e₀ : g₀ + u * h₀ = b * j₀)
    (e₁ : g₁ + u * h₁ = b * j₁)
    (e₂ : g₂ + u * h₂ = b * j₂) :
    (h₀ * j₁ - h₁ * j₀) * g₂ + (h₁ * j₂ - h₂ * j₁) * g₀ +
      (h₂ * j₀ - h₀ * j₂) * g₁ = 0 := by
  have r₀ : g₀ = b * j₀ - u * h₀ := by linarith
  have r₁ : g₁ = b * j₁ - u * h₁ := by linarith
  have r₂ : g₂ = b * j₂ - u * h₂ := by linarith
  rw [r₀, r₁, r₂]
  ring

open scoped BigOperators in
/-- If only observers in a small part of the survey see a signal, a nonzero
squared mean forces a large second moment. Apply to area divided by radius⁴
and a hinge neighborhood whose volume fraction is at most c times radius². -/
theorem localized_center_second_moment {ι : Type*} (s : Finset ι) (w x : ι → ℝ)
    (c r m : ℝ) (hw : ∀ i ∈ s, 0 ≤ w i)
    (mass : ∑ i ∈ s, w i ≤ c * r^2) (positive : 0 < c * r^2)
    (signal : m ≤ (∑ i ∈ s, w i * x i)^2) :
    m / (c * r^2) ≤ ∑ i ∈ s, w i * (x i)^2 := by
  have hnon : ∀ i ∈ s, 0 ≤ w i * (x i)^2 :=
    fun i hi => mul_nonneg (hw i hi) (sq_nonneg _)
  have cs := Finset.sum_sq_le_sum_mul_sum_of_sq_eq_mul
    (r := fun i => w i * x i) s hw hnon (fun i _ => by ring)
  apply (div_le_iff₀ positive).mpr
  calc
    m ≤ _ := signal
    _ ≤ _ := cs
    _ ≤ _ := by simpa [mul_comm] using
      mul_le_mul_of_nonneg_right mass (Finset.sum_nonneg hnon)

/-- On a locally isolated hinge, bounded modular response and fixed-volume
stationarity force the directional change of the cone angle to shrink with
the square of the proper radius. The area lower bound uses the independently
reconstructed cone; the other two estimates belong to the common source. -/
theorem isolated_hinge_angle_budget (ℓ η angle a s S ε : ℝ) (hr : 0 < ℓ)
    (hη : 0 < η) (hcone : (2 / 3 : ℝ) * ℓ ^ 2 * |angle| ≤ |a|)
    (hs : |s| ≤ S * ℓ ^ 4) (hstat : |s + η * a| ≤ ε * ℓ ^ 4) :
    2 * η * |angle| ≤ 3 * (S + ε) * ℓ ^ 2 := by
  have h : η * ((2 / 3 : ℝ) * ℓ ^ 2 * |angle|) ≤ (S + ε) * ℓ ^ 4 := calc
    _ ≤ η * |a| := mul_le_mul_of_nonneg_left hcone hη.le
    _ = |(s + η * a) - s| := by simp [abs_mul, abs_of_pos hη]
    _ ≤ |s + η * a| + |s| := abs_sub _ _
    _ ≤ ε * ℓ ^ 4 + S * ℓ ^ 4 := add_le_add hstat hs
    _ = _ := by ring
  apply le_of_mul_le_mul_left _ (sq_pos_of_pos hr)
  nlinarith only [h]

assert_no_sorry isolated_hinge_angle_budget
#print axioms isolated_hinge_angle_budget

assert_no_sorry localized_center_second_moment
#print axioms localized_center_second_moment

assert_no_sorry common_clock_flux_determinant
#print axioms common_clock_flux_determinant

assert_no_sorry regional_compensation
#print axioms regional_compensation

assert_no_sorry opposite_duration_responses
#print axioms opposite_duration_responses

assert_no_sorry unique_clock_score
#print axioms unique_clock_score

assert_no_sorry area_volume_blind_tangent
#print axioms area_volume_blind_tangent

assert_no_sorry survey_lagrangian_gap
assert_no_sorry survey_product_gap
#print axioms survey_lagrangian_gap
#print axioms survey_product_gap

assert_no_sorry constant_on_joint_survey
assert_no_sorry zero_common_survey_form
#print axioms constant_on_joint_survey
#print axioms zero_common_survey_form

assert_no_sorry conserved_port_effect
assert_no_sorry scalar_of_complete_returns

assert_no_sorry trace_volume_stationarity
assert_no_sorry fixed_volume_iff
assert_no_sorry ordinary_calibration
assert_no_sorry blind_transfer
assert_no_sorry calibrated_error
assert_no_sorry measured_residual

#print axioms conserved_port_effect
#print axioms scalar_of_complete_returns
#print axioms trace_volume_stationarity
#print axioms fixed_volume_iff
#print axioms ordinary_calibration
#print axioms blind_transfer
#print axioms calibrated_error
#print axioms measured_residual

end OPH.EinsteinBranch.StationarityReadout
