import TemporalMaxwellEvolution

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.ScaledMaxwellStability

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution

/-!
# Step-scaled Maxwell evolution on the committed carrier: action origin and stability

WHAT IS PROVED.  Exact finite real linear algebra on the committed
twelve-port, thirty-seam, twenty-oriented-face complex, reusing the exact
committed objects of `DiscreteCoulombGreen`, `PositionSpaceMaxwellAction`,
`LocalFaceMaxwellAction`, and `TemporalMaxwellEvolution`: the port
coboundary `realCoboundary` (d), its transpose `realBoundary` (∂), the
oriented face curvature `faceCurvature` (C), its transpose
`faceCodifferential` (Cᵀ), the local normal operator
`localMaxwellOperator` (CᵀC), the committed local sourced face action
`localSourcedAction`, and the committed unit-step packet
`electricField`, `magneticField`, `AmpereEvolution`.

(A) Scaled staggered packet.  A real step `h` scales the staggered
convention: `electricFieldScaled h A φ n = -(h⁻¹ • (A (n+1) - A n)) - d (φ n)`,
`magneticField` is reused, and `AmpereEvolutionScaled h A φ J` declares
`E (n+1) - E n = h • (Cᵀ (B (n+1)) - J n)`.  The scaled Faraday law
`B (n+1) - B n = -(h • C (E n))` is a theorem of the definitions
(`faraday_law_scaled`); the step `h = 1` recovers the committed packet
exactly (`electricFieldScaled_one`, `ampereEvolutionScaled_one_iff`); the
scaled Gauss-continuity equivalence and its propagation hold
(`gauss_step_iff_scaled`, `gauss_propagation_scaled`); the scaled
time-dependent gauge transformation leaves `E`, `B`, the evolution law and
the staggered form invariant; the scaled wave law
`A (n+2) - 2 A (n+1) + A n + h² • CᵀC (A (n+1)) + h • d (φ (n+1) - φ n) = 0`
holds for zero current (`wave_law_scaled`).

(B) Action origin.  `windowAction h N` is the finite-window sum of the
per-step Lagrangian `(h/2) ‖E n‖² - h • localSourcedAction (J n) (A (n+1))
+ h • ⟨ρ n, φ n⟩`: the declared discrete kinetic term `(h/2) ‖E n‖²` added
to the committed local face action with its committed seam-current
coupling, plus the port-load coupling typed on `Fin 12`.  The exact
expansion `windowAction (A + a) (φ + f) = windowAction A φ + firstVariation
+ quadraticRemainder` holds with an explicit linear functional and an
explicit quadratic remainder (`windowAction_expansion`); for variations
vanishing at the two window endpoints the first variation is the sum over
the interior steps of the Ampere residuals paired against `a`, minus `h`
times the sum over the window of the Gauss residuals paired against `f`
(`firstVariation_interior`).  Stationarity under all
endpoint-vanishing variations of `A` is equivalent to the scaled Ampere
update at every interior step (`action_stationary_A_iff_ampere`), and
stationarity under all variations of `φ` is equivalent to the Gauss
constraint `∂ (E n) = ρ n` at every window step
(`action_stationary_phi_iff_gauss`).  The port load `ρ : Fin 12` and the
seam current `J : Fin 30` stay distinct types throughout.

(C) Stability.  The scaled staggered form
`fieldEnergyScaled h A φ n = (1/2) ‖E n‖² + (1/2) ⟨B n, B (n+1)⟩` obeys the
exact balance `𝓔 (n+1) = 𝓔 n - (h/2) ⟨E n + E (n+1), J n⟩`
(`energy_balance_scaled`) and is conserved for zero current
(`energy_conserved_scaled`).  Under a Courant hypothesis
`∀ v, ‖C v‖² ≤ Λ ‖v‖²` the exact identity
`𝓔 n = (1/2) ‖E n‖² + (1/8) ‖C (A n + A (n+1))‖² - (h²/8) ‖C (E n)‖²`
gives the lower bound
`𝓔 n ≥ (1/2 - h² Λ / 8) ‖E n‖² + (1/8) ‖C (A n + A (n+1))‖²`
(`fieldEnergyScaled_lower_bound`), hence positivity for `h² Λ ≤ 4`
(`fieldEnergyScaled_nonneg`) and, for `h² Λ < 4` and zero current, the
uniform bounds `‖E n‖² ≤ 8 𝓔 0 / (4 - h² Λ)` and
`‖B n‖² ≤ 16 𝓔 0 / (4 - h² Λ)` for every step (`stability_certificate`).
Conversely, for an eigenvector `CᵀC v = λ v`, `v ≠ 0`, with `h² λ > 4`, the
history `A n = r ^ n • v` with `r` the real root of
`r² - (2 - h² λ) r + 1 = 0` below `-1` solves the zero-current evolution
in the gauge `φ = 0`, with `‖E n‖² = (r²)^n ‖E 0‖²`, `‖E 0‖² > 0`, and
`‖E n‖²` unbounded in `n` (`unstable_mode`).

(D) Committed-carrier spectral facts.  Three seams per face and two faces
per seam give the elementary Courant constant
`‖C v‖² ≤ 6 ‖v‖²` (`faceEnergy_curvature_le_six`), so `h² < 2/3` certifies
positivity and boundedness on the committed carrier.  An explicit integer
eigenvector of `CᵀC` with eigenvalue `5` (`fiveMode_eigen`) and an explicit
eigenvector with entries in `ℤ + ℤ √5` and eigenvalue `3 + √5`
(`goldenMode_eigen`) are kernel-checked from the committed incidence
tables; with (C) they give instability of the declared unit step
(`unit_step_instability`) and the necessary conditions `h² ≤ 4/5`
(`instability_above_five`) and `h² ≤ 4/(3 + √5)`
(`instability_above_golden`) for boundedness of every initial datum.  The
sharp constant is reached: the five spectral projectors of the face normal
matrix `N = C Cᵀ` (eigenvalues `0, 2, 3, 5` and the sector
`N² - 6N + 4 = 0`) are explicit rational matrices whose idempotence,
symmetry, eigen-relations and resolution of the identity are kernel
`decide` checks on integer literals; the golden sector obeys
`‖(N - 3) u‖² = 5 ‖u‖²`, so `⟨w, N w⟩ ≤ (3 + √5) ‖w‖²`
(`faceNormal_quadratic_bound`) and, by Cauchy-Schwarz through `Cᵀ`,
`‖C v‖² ≤ (3 + √5) ‖v‖²` (`faceEnergy_curvature_le_golden`).  The
threshold `h² (3 + √5) = 4` is therefore sharp on the committed carrier:
below it every zero-current solution is bounded, above it an unbounded
zero-current solution exists (`courant_threshold_sharp`).

(E) One composed receipt `scaledMaxwellStability_receipt` carries the
scaled identities, the action equivalences, conservation, positivity and
the bounds from the single typed antecedent bundle
`ScaledMaxwellBundle`; `demoScaledBundle` is an explicit nonstatic
inhabitant at `h = 1/2`, `Λ = 6`.

DECLARED INPUTS.  The step `h`, the discrete kinetic term `(h/2) ‖E n‖²`,
the port load `ρ` and the seam current `J` are declared inputs.

PREMISE CONSUMPTION.  PR-66 is relocated: the unit-step Ampere recurrence
of `TemporalMaxwellEvolution` is no separate declaration here; it is the
interior Euler-Lagrange equation of the declared discrete action
`windowAction`, whose only new ingredient relative to the committed static
packet `staticAction` is the kinetic term `(h/2) ‖E n‖²`.  What is
declared is that kinetic term and the step `h`.  PR-53 (physical
propagation, frame, and comparison attachment), PR-54 (source gauge-field,
current, and action attachment) and PR-15 (clock and energy calibration)
stay open and are not consumed.  The equal-weight pairing coincides with
the PR-20 equal seam-counting selection.

WHAT IS NOT PROVED.  The step index `n` is a declared evolution parameter,
not physical time.  The sources `ρ` and `J` are declared inputs, not
source-produced currents.  No photon, propagation speed, Lorentz
covariance, continuum limit, or laboratory readout is claimed.  The
kinetic term is declared, not derived from the three OPH axioms or from
repair dynamics.  At the threshold `h² (3 + √5) = 4` itself the staggered
form is nonnegative and conserved; a uniform bound at exact equality is not
stated.

Axiom audit.  Every proof composes the committed receipts with exact real
linear algebra and kernel `decide` checks on the committed integer
tables; the module adds no project axiom and uses no native decision
procedure.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## Pairing bookkeeping -/

theorem seamInner_neg_right (x y : Fin 30 → ℝ) :
    realSeamInner x (-y) = -realSeamInner x y := by
  unfold realSeamInner
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ by simp

theorem seamInner_smul_right (c : ℝ) (x y : Fin 30 → ℝ) :
    realSeamInner x (c • y) = c * realSeamInner x y := by
  unfold realSeamInner
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun e _ ↦ by simp [mul_left_comm]

theorem seamInner_add_left (x y z : Fin 30 → ℝ) :
    realSeamInner (x + y) z = realSeamInner x z + realSeamInner y z := by
  rw [realSeamInner_comm, realSeamInner_add_right, realSeamInner_comm z x,
    realSeamInner_comm z y]

theorem seamEnergy_add (x y : Fin 30 → ℝ) :
    realSeamEnergy (x + y) =
      realSeamEnergy x + 2 * realSeamInner x y + realSeamEnergy y := by
  unfold realSeamEnergy realSeamInner
  rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ by
    simp only [Pi.add_apply]
    ring

theorem seamEnergy_smul (c : ℝ) (x : Fin 30 → ℝ) :
    realSeamEnergy (c • x) = c ^ 2 * realSeamEnergy x := by
  unfold realSeamEnergy
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun e _ ↦ by
    simp only [Pi.smul_apply, smul_eq_mul]
    ring

theorem seamEnergy_neg (x : Fin 30 → ℝ) :
    realSeamEnergy (-x) = realSeamEnergy x := by
  unfold realSeamEnergy
  exact Finset.sum_congr rfl fun e _ ↦ by simp

theorem seamEnergy_pos_of_ne_zero (x : Fin 30 → ℝ) (hx : x ≠ 0) :
    0 < realSeamEnergy x := by
  rcases (realSeamEnergy_nonneg x).lt_or_eq with hlt | heq
  · exact hlt
  · exact absurd ((realSeamEnergy_eq_zero_iff x).mp heq.symm) hx

theorem portInner_add_right (x y z : Fin 12 → ℝ) :
    realPortInner x (y + z) = realPortInner x y + realPortInner x z := by
  unfold realPortInner
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun p _ ↦ mul_add _ _ _

theorem portInner_sub_left (x y z : Fin 12 → ℝ) :
    realPortInner (x - y) z = realPortInner x z - realPortInner y z := by
  unfold realPortInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun p _ ↦ sub_mul _ _ _

theorem portInner_zero_right (x : Fin 12 → ℝ) : realPortInner x 0 = 0 := by
  unfold realPortInner
  simp

theorem portInner_zero_left (x : Fin 12 → ℝ) : realPortInner 0 x = 0 := by
  unfold realPortInner
  simp

theorem faceInner_smul_left (c : ℝ) (F G : Fin 20 → ℝ) :
    faceInner (c • F) G = c * faceInner F G := by
  unfold faceInner
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun f _ ↦ by simp [mul_assoc]

theorem faceInner_neg_left (F G : Fin 20 → ℝ) :
    faceInner (-F) G = -faceInner F G := by
  unfold faceInner
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl fun f _ ↦ by simp

theorem faceEnergy_smul (c : ℝ) (F : Fin 20 → ℝ) :
    faceEnergy (c • F) = c ^ 2 * faceEnergy F := by
  unfold faceEnergy
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.smul_apply, smul_eq_mul]
    ring

theorem faceEnergy_neg (F : Fin 20 → ℝ) : faceEnergy (-F) = faceEnergy F := by
  unfold faceEnergy
  exact Finset.sum_congr rfl fun f _ ↦ by simp

/-- Polarization for the face pairing:
`⟨F, G⟩ = (1/4) (‖F + G‖² - ‖G - F‖²)`. -/
theorem faceInner_polarization (F G : Fin 20 → ℝ) :
    faceInner F G = (1 / 4) * (faceEnergy (F + G) - faceEnergy (G - F)) := by
  unfold faceInner faceEnergy
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.add_apply, Pi.sub_apply]
    ring

/-- Parallelogram bound for the face energy:
`‖F‖² ≤ (1/2) ‖F + G‖² + (1/2) ‖G - F‖²`. -/
theorem faceEnergy_le_parallelogram (F G : Fin 20 → ℝ) :
    faceEnergy F ≤ (1 / 2) * faceEnergy (F + G) + (1 / 2) * faceEnergy (G - F) := by
  unfold faceEnergy
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_le_sum fun f _ ↦ ?_
  simp only [Pi.add_apply, Pi.sub_apply]
  nlinarith [sq_nonneg (G f)]

/-! ## (A) The scaled staggered packet -/

/-- Scaled electric seam field on the half step between `n` and `n+1`:
`E n = -(h⁻¹ • (A (n+1) - A n)) - d (φ n)`. -/
def electricFieldScaled (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : Fin 30 → ℝ :=
  -(h⁻¹ • (A (n + 1) - A n)) - realCoboundary (φ n)

/-- The scaled Ampere update: `E (n+1) - E n = h • (Cᵀ (B (n+1)) - J n)`. -/
def AmpereEvolutionScaled (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) : Prop :=
  ∀ n : ℕ, electricFieldScaled h A φ (n + 1) - electricFieldScaled h A φ n =
    h • (faceCodifferential (magneticField A (n + 1)) - J n)

/-- The step `h = 1` recovers the committed electric field exactly. -/
theorem electricFieldScaled_one (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : electricFieldScaled 1 A φ n = electricField A φ n := by
  unfold electricFieldScaled electricField
  rw [inv_one, one_smul]

/-- The step `h = 1` recovers the committed Ampere update exactly. -/
theorem ampereEvolutionScaled_one_iff (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) :
    AmpereEvolutionScaled 1 A φ J ↔ AmpereEvolution A φ J := by
  unfold AmpereEvolutionScaled AmpereEvolution
  simp only [electricFieldScaled_one, one_smul]

/-- Linearity of the scaled electric field in the pair `(A, φ)`. -/
theorem electricFieldScaled_add (h : ℝ) (A a : ℕ → Fin 30 → ℝ)
    (φ f : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricFieldScaled h (A + a) (φ + f) n =
      electricFieldScaled h A φ n + electricFieldScaled h a f n := by
  unfold electricFieldScaled
  funext e
  simp only [map_add, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.smul_apply,
    smul_eq_mul]
  ring

/-- The face curvature of the scaled electric field is the backward
magnetic difference over `h`; the scalar term drops by `C ∘ d = 0`. -/
theorem faceCurvature_electricFieldScaled (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    faceCurvature (electricFieldScaled h A φ n) =
      h⁻¹ • (magneticField A n - magneticField A (n + 1)) := by
  unfold electricFieldScaled magneticField
  rw [map_sub, faceCurvature_coboundary, sub_zero, map_neg, map_smul, map_sub,
    ← smul_neg, neg_sub]

/-- Scaled Faraday law as a theorem of the definitions. -/
theorem faraday_law_scaled (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    magneticField A (n + 1) - magneticField A n =
      -(h • faceCurvature (electricFieldScaled h A φ n)) := by
  rw [faceCurvature_electricFieldScaled, smul_smul, mul_inv_cancel₀ hh, one_smul,
    neg_sub]

/-- Scaled continuity, both directions: along the scaled evolution the
Gauss constraint propagates one step exactly when
`ρ (n+1) - ρ n + h • ∂ (J n) = 0`. -/
theorem gauss_step_iff_scaled (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ)
    (hn : realBoundary (electricFieldScaled h A φ n) = ρ n) :
    realBoundary (electricFieldScaled h A φ (n + 1)) = ρ (n + 1) ↔
      ρ (n + 1) - ρ n + h • realBoundary (J n) = 0 := by
  have hkey : realBoundary (electricFieldScaled h A φ (n + 1)) =
      ρ n - h • realBoundary (J n) := by
    have h1 := congrArg realBoundary (hAmp n)
    rw [map_sub, map_smul, map_sub, realBoundary_faceCodifferential, hn, zero_sub,
      smul_neg, sub_eq_iff_eq_add] at h1
    rw [h1]
    abel
  rw [hkey]
  constructor
  · intro heq
    rw [← heq]
    abel
  · intro heq
    calc ρ n - h • realBoundary (J n)
        = ρ (n + 1) - (ρ (n + 1) - ρ n + h • realBoundary (J n)) := by abel
      _ = ρ (n + 1) := by rw [heq, sub_zero]

/-- Under the scaled continuity equation the Gauss constraint propagates
from the initial step to every step. -/
theorem gauss_propagation_scaled (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J)
    (h0 : realBoundary (electricFieldScaled h A φ 0) = ρ 0)
    (hcont : ∀ n : ℕ, ρ (n + 1) - ρ n + h • realBoundary (J n) = 0) :
    ∀ n : ℕ, realBoundary (electricFieldScaled h A φ n) = ρ n := by
  intro n
  induction n with
  | zero => exact h0
  | succ m ih => exact (gauss_step_iff_scaled h A φ J ρ hAmp m ih).mpr (hcont m)

/-- Scaled scalar gauge transformation `φ n ↦ φ n - h⁻¹ • (χ (n+1) - χ n)`,
paired with the committed `gaugeTransformA`. -/
def gaugeTransformPhiScaled (h : ℝ) (φ : ℕ → Fin 12 → ℝ)
    (χ : ℕ → Fin 12 → ℝ) : ℕ → Fin 12 → ℝ :=
  fun n ↦ φ n - h⁻¹ • (χ (n + 1) - χ n)

theorem electricFieldScaled_gauge_invariant (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (χ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricFieldScaled h (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) n
      = electricFieldScaled h A φ n := by
  unfold electricFieldScaled gaugeTransformA gaugeTransformPhiScaled
  funext e
  simp only [map_sub, map_smul, Pi.sub_apply, Pi.neg_apply, Pi.add_apply,
    Pi.smul_apply, smul_eq_mul, realCoboundary_apply]
  ring

theorem ampereEvolutionScaled_gauge_invariant (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (χ : ℕ → Fin 12 → ℝ) :
    AmpereEvolutionScaled h (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) J
      ↔ AmpereEvolutionScaled h A φ J := by
  unfold AmpereEvolutionScaled
  simp only [electricFieldScaled_gauge_invariant, magneticField_gauge_invariant]

/-- Scaled wave law for zero current:
`A (n+2) - 2 A (n+1) + A n + h² • CᵀC (A (n+1)) + h • d (φ (n+1) - φ n) = 0`. -/
theorem wave_law_scaled (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (hAmp : AmpereEvolutionScaled h A φ (fun _ ↦ 0))
    (n : ℕ) :
    A (n + 2) - (2 : ℝ) • A (n + 1) + A n +
      (h ^ 2) • localMaxwellOperator (A (n + 1)) +
      h • realCoboundary (φ (n + 1) - φ n) = 0 := by
  have hE := hAmp n
  have hlm : faceCodifferential (magneticField A (n + 1)) =
      localMaxwellOperator (A (n + 1)) := rfl
  rw [hlm] at hE
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  funext e
  have he := congrFun hE e
  simp only [electricFieldScaled, map_sub, Pi.sub_apply, Pi.neg_apply,
    Pi.add_apply, Pi.smul_apply, smul_eq_mul, realCoboundary_apply,
    Pi.zero_apply, sub_zero] at he ⊢
  linear_combination (-h) * he -
    (A (n + 2) e - 2 * A (n + 1) e + A n e) * hinv

/-! ## (C) The scaled staggered form: balance, conservation -/

/-- Scaled staggered quadratic form `(1/2) ‖E n‖² + (1/2) ⟨B n, B (n+1)⟩`. -/
def fieldEnergyScaled (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : ℝ :=
  (1 / 2) * realSeamEnergy (electricFieldScaled h A φ n) +
    (1 / 2) * faceInner (magneticField A n) (magneticField A (n + 1))

theorem fieldEnergyScaled_one (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : fieldEnergyScaled 1 A φ n = fieldEnergy A φ n := by
  unfold fieldEnergyScaled fieldEnergy
  rw [electricFieldScaled_one]

theorem fieldEnergyScaled_gauge_invariant (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (χ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    fieldEnergyScaled h (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) n
      = fieldEnergyScaled h A φ n := by
  unfold fieldEnergyScaled
  rw [electricFieldScaled_gauge_invariant, magneticField_gauge_invariant,
    magneticField_gauge_invariant]

/-- Exact per-step balance of the scaled staggered form. -/
theorem energy_balance_scaled (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n -
      (h / 2) * realSeamInner
        (electricFieldScaled h A φ n + electricFieldScaled h A φ (n + 1)) (J n) := by
  set E := electricFieldScaled h A φ with hEdef
  set B := magneticField A with hBdef
  have hdiff : realSeamEnergy (E (n + 1)) - realSeamEnergy (E n) =
      realSeamInner (E (n + 1) + E n) (E (n + 1) - E n) :=
    realSeamEnergy_sub_eq_inner _ _
  rw [hAmp n, seamInner_smul_right, realSeamInner_sub_right] at hdiff
  have hadj : realSeamInner (E (n + 1) + E n) (faceCodifferential (B (n + 1))) =
      faceInner (faceCurvature (E (n + 1) + E n)) (B (n + 1)) :=
    (faceCurvature_codifferential_adjoint _ _).symm
  have hCsum : faceCurvature (E (n + 1) + E n) =
      h⁻¹ • (B n - B (n + 1 + 1)) := by
    rw [map_add, hEdef, faceCurvature_electricFieldScaled,
      faceCurvature_electricFieldScaled, ← smul_add, hBdef]
    congr 1
    abel
  rw [hCsum, faceInner_smul_left, faceInner_sub_left] at hadj
  rw [hadj] at hdiff
  have hflip : faceInner (B (n + 1 + 1)) (B (n + 1)) =
      faceInner (B (n + 1)) (B (n + 1 + 1)) := faceInner_comm _ _
  have hEflip : realSeamInner (E (n + 1) + E n) (J n) =
      realSeamInner (E n + E (n + 1)) (J n) := by
    rw [add_comm (E (n + 1)) (E n)]
  rw [hflip, hEflip] at hdiff
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  unfold fieldEnergyScaled
  rw [← hEdef, ← hBdef]
  linear_combination (1 / 2) * hdiff +
    ((1 / 2) * (faceInner (B n) (B (n + 1)) -
      faceInner (B (n + 1)) (B (n + 1 + 1)))) * hinv

/-- Exact conservation of the scaled staggered form for zero current. -/
theorem energy_conserved_scaled (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (hAmp : AmpereEvolutionScaled h A φ (fun _ ↦ 0))
    (n : ℕ) : fieldEnergyScaled h A φ n = fieldEnergyScaled h A φ 0 := by
  induction n with
  | zero => rfl
  | succ m ih =>
    have hb := energy_balance_scaled h hh A φ (fun _ ↦ 0) hAmp m
    rw [realSeamInner_zero_right, mul_zero, sub_zero] at hb
    rw [hb]
    exact ih

/-! ## (B) The discrete action and its Euler-Lagrange equations -/

/-- Per-step Lagrangian: the declared kinetic term `(h/2) ‖E n‖²`, minus
`h` times the committed local sourced face action at the forward step
`A (n+1)` with seam current `J n`, plus `h` times the port-load coupling
`⟨ρ n, φ n⟩`.  The two sources keep their distinct types. -/
def stepLagrangian (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (n : ℕ) : ℝ :=
  (h / 2) * realSeamEnergy (electricFieldScaled h A φ n)
    - h * localSourcedAction (J n) (A (n + 1))
    + h * realPortInner (ρ n) (φ n)

/-- Window action over the steps `0, …, N`. -/
def windowAction (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) : ℝ :=
  ∑ n ∈ Finset.range (N + 1), stepLagrangian h A φ ρ J n

/-- The linear functional of the exact expansion. -/
def firstVariation (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (a : ℕ → Fin 30 → ℝ)
    (f : ℕ → Fin 12 → ℝ) : ℝ :=
  ∑ n ∈ Finset.range (N + 1),
    (h * realSeamInner (electricFieldScaled h A φ n) (electricFieldScaled h a f n)
      - h * realSeamInner (localMaxwellOperator (A (n + 1)) - J n) (a (n + 1))
      + h * realPortInner (ρ n) (f n))

/-- The explicit quadratic remainder of the exact expansion. -/
def quadraticRemainder (h : ℝ) (N : ℕ) (a : ℕ → Fin 30 → ℝ)
    (f : ℕ → Fin 12 → ℝ) : ℝ :=
  ∑ n ∈ Finset.range (N + 1),
    ((h / 2) * realSeamEnergy (electricFieldScaled h a f n)
      - (h / 2) * faceEnergy (faceCurvature (a (n + 1))))

theorem stepLagrangian_expansion (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (a : ℕ → Fin 30 → ℝ) (f : ℕ → Fin 12 → ℝ) (n : ℕ) :
    stepLagrangian h (A + a) (φ + f) ρ J n = stepLagrangian h A φ ρ J n
      + (h * realSeamInner (electricFieldScaled h A φ n) (electricFieldScaled h a f n)
          - h * realSeamInner (localMaxwellOperator (A (n + 1)) - J n) (a (n + 1))
          + h * realPortInner (ρ n) (f n))
      + ((h / 2) * realSeamEnergy (electricFieldScaled h a f n)
          - (h / 2) * faceEnergy (faceCurvature (a (n + 1)))) := by
  unfold stepLagrangian
  rw [electricFieldScaled_add, seamEnergy_add]
  simp only [Pi.add_apply]
  rw [localSourcedAction_expansion, portInner_add_right]
  ring

/-- **Exact expansion of the window action.** -/
theorem windowAction_expansion (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (a : ℕ → Fin 30 → ℝ) (f : ℕ → Fin 12 → ℝ) :
    windowAction h N (A + a) (φ + f) ρ J = windowAction h N A φ ρ J +
      firstVariation h N A φ ρ J a f + quadraticRemainder h N a f := by
  unfold windowAction firstVariation quadraticRemainder
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun n _ ↦ stepLagrangian_expansion h A φ ρ J a f n

/-- The Ampere residual at step `m`: zero exactly when the scaled update
holds across step `m+1`. -/
def ampereResidual (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (m : ℕ) : Fin 30 → ℝ :=
  electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m -
    h • (localMaxwellOperator (A (m + 1)) - J m)

theorem ampereEvolutionScaled_iff_residual (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) :
    AmpereEvolutionScaled h A φ J ↔ ∀ m, ampereResidual h A φ J m = 0 := by
  unfold AmpereEvolutionScaled ampereResidual
  constructor
  · intro hA m
    rw [sub_eq_zero]
    exact hA m
  · intro hR m
    exact sub_eq_zero.mp (hR m)

theorem seamInner_smul_left (c : ℝ) (x y : Fin 30 → ℝ) :
    realSeamInner (c • x) y = c * realSeamInner x y := by
  rw [realSeamInner_comm, seamInner_smul_right, realSeamInner_comm]

/-- Pairing a seam field against the scaled electric field of a variation:
the kinetic term pairs against the forward difference and the port
boundary. -/
theorem kinetic_pairing (h : ℝ) (hh : h ≠ 0) (x : Fin 30 → ℝ)
    (a : ℕ → Fin 30 → ℝ) (f : ℕ → Fin 12 → ℝ) (n : ℕ) :
    h * realSeamInner x (electricFieldScaled h a f n) =
      realSeamInner x (a n) - realSeamInner x (a (n + 1)) -
        h * realPortInner (realBoundary x) (f n) := by
  unfold electricFieldScaled
  rw [realSeamInner_sub_right, seamInner_neg_right, seamInner_smul_right,
    realSeamInner_sub_right, realSeamInner_comm x (realCoboundary (f n)),
    realCoboundary_boundary_adjoint, realPortInner_comm (f n) (realBoundary x)]
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  linear_combination (realSeamInner x (a n) - realSeamInner x (a (n + 1))) * hinv

/-- **Euler-Lagrange form of the first variation.**  For variations of `A`
vanishing at both window endpoints, the first variation is the sum of the
Ampere residuals paired against `a` at the interior steps, minus `h` times
the Gauss residuals paired against `f`. -/
theorem firstVariation_interior (h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (a : ℕ → Fin 30 → ℝ) (f : ℕ → Fin 12 → ℝ)
    (ha0 : a 0 = 0) (haN : a (N + 1) = 0) :
    firstVariation h N A φ ρ J a f =
      (∑ m ∈ Finset.range N,
        realSeamInner (ampereResidual h A φ J m) (a (m + 1))) -
      h * ∑ n ∈ Finset.range (N + 1),
        realPortInner (realBoundary (electricFieldScaled h A φ n) - ρ n) (f n) := by
  unfold firstVariation
  set E := electricFieldScaled h A φ with hE
  set K : ℕ → Fin 30 → ℝ := fun n ↦ localMaxwellOperator (A (n + 1)) - J n with hK
  have hterm : ∀ n, h * realSeamInner (E n) (electricFieldScaled h a f n)
      - h * realSeamInner (localMaxwellOperator (A (n + 1)) - J n) (a (n + 1))
      + h * realPortInner (ρ n) (f n)
      = (realSeamInner (E n) (a n) - realSeamInner (E n + h • K n) (a (n + 1)))
        - h * realPortInner (realBoundary (E n) - ρ n) (f n) := by
    intro n
    rw [kinetic_pairing h hh, seamInner_add_left, seamInner_smul_left,
      portInner_sub_left]
    simp only [hK]
    ring
  rw [Finset.sum_congr rfl fun n _ ↦ hterm n, Finset.sum_sub_distrib,
    Finset.sum_sub_distrib, ← Finset.mul_sum]
  congr 1
  have h1 : (∑ n ∈ Finset.range (N + 1), realSeamInner (E n) (a n)) =
      ∑ m ∈ Finset.range N, realSeamInner (E (m + 1)) (a (m + 1)) := by
    rw [Finset.sum_range_succ', ha0, realSeamInner_zero_right, add_zero]
  have h2 : (∑ n ∈ Finset.range (N + 1),
      realSeamInner (E n + h • K n) (a (n + 1))) =
      ∑ m ∈ Finset.range N, realSeamInner (E m + h • K m) (a (m + 1)) := by
    rw [Finset.sum_range_succ, haN, realSeamInner_zero_right, add_zero]
  rw [h1, h2, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun m _ ↦ ?_
  unfold ampereResidual
  rw [realSeamInner_sub_left, realSeamInner_sub_left, seamInner_add_left]
  simp only [hK, hE]
  ring

/-- **Stationarity in `A` is the scaled Ampere update at the interior
steps.**  Stationarity means: for every variation of `A` vanishing at both
window endpoints, the action moves by exactly the quadratic remainder. -/
theorem action_stationary_A_iff_ampere (h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
      windowAction h N (A + a) φ ρ J =
        windowAction h N A φ ρ J + quadraticRemainder h N a 0) ↔
    (∀ m, m < N → ampereResidual h A φ J m = 0) := by
  have hexp : ∀ a, windowAction h N (A + a) φ ρ J = windowAction h N A φ ρ J +
      firstVariation h N A φ ρ J a 0 + quadraticRemainder h N a 0 := by
    intro a
    have hx := windowAction_expansion h N A φ ρ J a 0
    rwa [add_zero] at hx
  have hfv : ∀ a, a 0 = 0 → a (N + 1) = 0 → firstVariation h N A φ ρ J a 0 =
      ∑ m ∈ Finset.range N, realSeamInner (ampereResidual h A φ J m) (a (m + 1)) := by
    intro a ha0 haN
    rw [firstVariation_interior h hh N A φ ρ J a 0 ha0 haN]
    simp only [Pi.zero_apply, portInner_zero_right, Finset.sum_const_zero, mul_zero,
      sub_zero]
  constructor
  · intro hstat m hm
    let a : ℕ → Fin 30 → ℝ := fun n ↦ if n = m + 1 then ampereResidual h A φ J m else 0
    have ha0 : a 0 = 0 := by simp [a]
    have haN : a (N + 1) = 0 := by
      show (if N + 1 = m + 1 then ampereResidual h A φ J m else 0) = 0
      rw [if_neg (by omega)]
    have h1 := hexp a
    rw [hstat a ha0 haN, hfv a ha0 haN] at h1
    have hsum : (∑ k ∈ Finset.range N,
        realSeamInner (ampereResidual h A φ J k) (a (k + 1))) = 0 := by linarith
    have hk : ∀ k, realSeamInner (ampereResidual h A φ J k) (a (k + 1)) =
        if k = m then realSeamEnergy (ampereResidual h A φ J m) else 0 := by
      intro k
      by_cases hkm : k = m
      · rw [if_pos hkm, hkm]
        show realSeamInner (ampereResidual h A φ J m)
          (if m + 1 = m + 1 then ampereResidual h A φ J m else 0) = _
        rw [if_pos rfl, realSeamInner_self_eq_energy]
      · rw [if_neg hkm]
        have hne : k + 1 ≠ m + 1 := by omega
        show realSeamInner (ampereResidual h A φ J k)
          (if k + 1 = m + 1 then ampereResidual h A φ J m else 0) = 0
        rw [if_neg hne, realSeamInner_zero_right]
    rw [Finset.sum_congr rfl fun k _ ↦ hk k, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hm)] at hsum
    exact (realSeamEnergy_eq_zero_iff _).mp hsum
  · intro hres a ha0 haN
    rw [hexp a, hfv a ha0 haN]
    have hz : (∑ m ∈ Finset.range N,
        realSeamInner (ampereResidual h A φ J m) (a (m + 1))) = 0 := by
      refine Finset.sum_eq_zero fun m hm ↦ ?_
      rw [hres m (Finset.mem_range.mp hm), realSeamInner_zero_left]
    rw [hz, add_zero]

/-- **Stationarity in `φ` is the Gauss constraint at every window step.** -/
theorem action_stationary_phi_iff_gauss (h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    (∀ f : ℕ → Fin 12 → ℝ,
      windowAction h N A (φ + f) ρ J =
        windowAction h N A φ ρ J + quadraticRemainder h N 0 f) ↔
    (∀ n, n < N + 1 → realBoundary (electricFieldScaled h A φ n) = ρ n) := by
  have hexp : ∀ f, windowAction h N A (φ + f) ρ J = windowAction h N A φ ρ J +
      firstVariation h N A φ ρ J 0 f + quadraticRemainder h N 0 f := by
    intro f
    have hx := windowAction_expansion h N A φ ρ J 0 f
    rwa [add_zero] at hx
  have hfv : ∀ f, firstVariation h N A φ ρ J 0 f =
      -(h * ∑ n ∈ Finset.range (N + 1),
        realPortInner (realBoundary (electricFieldScaled h A φ n) - ρ n) (f n)) := by
    intro f
    rw [firstVariation_interior h hh N A φ ρ J 0 f rfl rfl]
    simp only [Pi.zero_apply, realSeamInner_zero_right, Finset.sum_const_zero,
      zero_sub]
  constructor
  · intro hstat n hn
    let f : ℕ → Fin 12 → ℝ := fun k ↦
      if k = n then realBoundary (electricFieldScaled h A φ n) - ρ n else 0
    have h1 := hexp f
    rw [hstat f, hfv f] at h1
    have hsum : h * ∑ k ∈ Finset.range (N + 1),
        realPortInner (realBoundary (electricFieldScaled h A φ k) - ρ k) (f k) = 0 := by
      linarith
    have hsum0 : (∑ k ∈ Finset.range (N + 1),
        realPortInner (realBoundary (electricFieldScaled h A φ k) - ρ k) (f k)) = 0 :=
      (mul_eq_zero.mp hsum).resolve_left hh
    have hk : ∀ k, realPortInner (realBoundary (electricFieldScaled h A φ k) - ρ k) (f k)
        = if k = n then realPortInner
            (realBoundary (electricFieldScaled h A φ n) - ρ n)
            (realBoundary (electricFieldScaled h A φ n) - ρ n) else 0 := by
      intro k
      by_cases hkn : k = n
      · rw [if_pos hkn, hkn]
        show realPortInner _ (if n = n then _ else 0) = _
        rw [if_pos rfl]
      · rw [if_neg hkn]
        show realPortInner _ (if k = n then _ else 0) = 0
        rw [if_neg hkn, portInner_zero_right]
    rw [Finset.sum_congr rfl fun k _ ↦ hk k, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hn)] at hsum0
    exact sub_eq_zero.mp ((realPortInner_self_eq_zero_iff _).mp hsum0)
  · intro hg f
    rw [hexp f, hfv f]
    have hz : (∑ n ∈ Finset.range (N + 1),
        realPortInner (realBoundary (electricFieldScaled h A φ n) - ρ n) (f n)) = 0 := by
      refine Finset.sum_eq_zero fun n hn ↦ ?_
      rw [hg n (Finset.mem_range.mp hn), sub_self, portInner_zero_left]
    rw [hz, mul_zero, neg_zero, add_zero]

/-! ## (C) Positivity, bounds, and the unstable mode -/

/-- Courant hypothesis: `‖C v‖² ≤ Λ ‖v‖²` for every seam field. -/
def CourantBound (Λ : ℝ) : Prop :=
  ∀ v : Fin 30 → ℝ, faceEnergy (faceCurvature v) ≤ Λ * realSeamEnergy v

/-- Exact decomposition of the scaled staggered form. -/
theorem fieldEnergyScaled_eq (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    fieldEnergyScaled h A φ n =
      (1 / 2) * realSeamEnergy (electricFieldScaled h A φ n) +
      (1 / 8) * faceEnergy (faceCurvature (A n + A (n + 1))) -
      (h ^ 2 / 8) * faceEnergy (faceCurvature (electricFieldScaled h A φ n)) := by
  unfold fieldEnergyScaled
  rw [faceInner_polarization, faraday_law_scaled h hh A φ n, faceEnergy_neg,
    faceEnergy_smul]
  have hS : magneticField A n + magneticField A (n + 1) =
      faceCurvature (A n + A (n + 1)) := by
    unfold magneticField
    rw [map_add]
  rw [hS]
  ring

/-- **Lower bound under the Courant hypothesis.** -/
theorem fieldEnergyScaled_lower_bound (h : ℝ) (hh : h ≠ 0) (Λ : ℝ)
    (hΛ : CourantBound Λ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    (1 / 2 - h ^ 2 * Λ / 8) * realSeamEnergy (electricFieldScaled h A φ n) +
      (1 / 8) * faceEnergy (faceCurvature (A n + A (n + 1))) ≤
      fieldEnergyScaled h A φ n := by
  rw [fieldEnergyScaled_eq h hh]
  have hc := hΛ (electricFieldScaled h A φ n)
  have hsq : 0 ≤ h ^ 2 := sq_nonneg h
  have hm := mul_le_mul_of_nonneg_left hc hsq
  linarith

/-- **Positivity** for `h² Λ ≤ 4`. -/
theorem fieldEnergyScaled_nonneg (h : ℝ) (hh : h ≠ 0) (Λ : ℝ)
    (hΛ : CourantBound Λ) (hc : h ^ 2 * Λ ≤ 4) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) : 0 ≤ fieldEnergyScaled h A φ n := by
  have hb := fieldEnergyScaled_lower_bound h hh Λ hΛ A φ n
  have h1 : 0 ≤ (1 / 2 - h ^ 2 * Λ / 8) *
      realSeamEnergy (electricFieldScaled h A φ n) :=
    mul_nonneg (by linarith) (realSeamEnergy_nonneg _)
  have h2 := faceEnergy_nonneg (faceCurvature (A n + A (n + 1)))
  linarith

theorem electricFieldScaled_temporal_gauge (h : ℝ) (A : ℕ → Fin 30 → ℝ) (n : ℕ) :
    electricFieldScaled h A (fun _ ↦ 0) n = -(h⁻¹ • (A (n + 1) - A n)) := by
  unfold electricFieldScaled
  rw [map_zero, sub_zero]

/-- The lower bound in the gauge `φ = 0`, in terms of the seam increment. -/
theorem fieldEnergyScaled_lower_bound_temporal_gauge (h : ℝ) (hh : h ≠ 0) (Λ : ℝ)
    (hΛ : CourantBound Λ) (A : ℕ → Fin 30 → ℝ) (n : ℕ) :
    (1 / (2 * h ^ 2) - Λ / 8) * realSeamEnergy (A (n + 1) - A n) +
      (1 / 8) * faceEnergy (faceCurvature (A n + A (n + 1))) ≤
      fieldEnergyScaled h A (fun _ ↦ 0) n := by
  have hb := fieldEnergyScaled_lower_bound h hh Λ hΛ A (fun _ ↦ 0) n
  rw [electricFieldScaled_temporal_gauge, seamEnergy_neg, seamEnergy_smul] at hb
  have hcoef : (1 / (2 * h ^ 2) - Λ / 8) =
      (1 / 2 - h ^ 2 * Λ / 8) * (h⁻¹) ^ 2 := by
    field_simp
  rw [hcoef]
  linarith

/-- **Stability certificate.**  For `h² Λ < 4` and zero current, the
electric seam energy and the magnetic face energy stay bounded at every
step by explicit multiples of the initial staggered form. -/
theorem stability_certificate (h : ℝ) (hh : h ≠ 0) (Λ : ℝ) (hΛ0 : 0 ≤ Λ)
    (hΛ : CourantBound Λ) (hc : h ^ 2 * Λ < 4) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (hAmp : AmpereEvolutionScaled h A φ (fun _ ↦ 0)) (n : ℕ) :
    realSeamEnergy (electricFieldScaled h A φ n) ≤
        8 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * Λ) ∧
      faceEnergy (magneticField A n) ≤
        16 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * Λ) := by
  have hμ0 : 0 ≤ h ^ 2 * Λ := mul_nonneg (sq_nonneg h) hΛ0
  have hK : 0 < 4 - h ^ 2 * Λ := by linarith
  have hcons := energy_conserved_scaled h hh A φ hAmp n
  have hb := fieldEnergyScaled_lower_bound h hh Λ hΛ A φ n
  rw [hcons] at hb
  have hX0 : 0 ≤ realSeamEnergy (electricFieldScaled h A φ n) := realSeamEnergy_nonneg _
  have hY0 : 0 ≤ faceEnergy (faceCurvature (A n + A (n + 1))) := faceEnergy_nonneg _
  have hE : realSeamEnergy (electricFieldScaled h A φ n) ≤
      8 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * Λ) := by
    rw [le_div_iff₀ hK]
    nlinarith
  refine ⟨hE, ?_⟩
  have hpar := faceEnergy_le_parallelogram (magneticField A n) (magneticField A (n + 1))
  have hS : magneticField A n + magneticField A (n + 1) =
      faceCurvature (A n + A (n + 1)) := by
    unfold magneticField
    rw [map_add]
  rw [hS] at hpar
  have hD : faceEnergy (magneticField A (n + 1) - magneticField A n) ≤
      h ^ 2 * Λ * realSeamEnergy (electricFieldScaled h A φ n) := by
    rw [faraday_law_scaled h hh A φ n, faceEnergy_neg, faceEnergy_smul]
    have hc' := hΛ (electricFieldScaled h A φ n)
    have := mul_le_mul_of_nonneg_left hc' (sq_nonneg h)
    linarith
  rw [le_div_iff₀ hK]
  have hp1 := mul_le_mul_of_nonneg_right hpar hK.le
  have hp2 := mul_le_mul_of_nonneg_right hD hK.le
  have hp3 := mul_nonneg hY0 hμ0
  have hp4 := mul_nonneg hX0 (sq_nonneg (4 - h ^ 2 * Λ))
  nlinarith

/-- The real root of `r² - (2 - μ) r + 1 = 0` below `-1`, for `μ > 4`. -/
def growthRoot (μ : ℝ) : ℝ := ((2 - μ) - Real.sqrt (μ ^ 2 - 4 * μ)) / 2

theorem growthRoot_root (μ : ℝ) (hμ : 4 < μ) :
    growthRoot μ ^ 2 - (2 - μ) * growthRoot μ + 1 = 0 := by
  unfold growthRoot
  have h0 : 0 ≤ μ ^ 2 - 4 * μ := by nlinarith
  have hs := Real.sq_sqrt h0
  linear_combination (1 / 4) * hs

theorem growthRoot_lt_neg_one (μ : ℝ) (hμ : 4 < μ) : growthRoot μ < -1 := by
  unfold growthRoot
  have := Real.sqrt_nonneg (μ ^ 2 - 4 * μ)
  linarith

theorem one_lt_sq_of_lt_neg_one (r : ℝ) (hr : r < -1) : 1 < r ^ 2 := by
  nlinarith [mul_pos_of_neg_of_neg (show r + 1 < 0 by linarith)
    (show r - 1 < 0 by linarith)]

/-- The eigenmode history `A n = r ^ n • v`. -/
def modeHistory (r : ℝ) (v : Fin 30 → ℝ) : ℕ → Fin 30 → ℝ := fun n ↦ r ^ n • v

theorem modeHistory_electricField (h r : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) n =
      (-(h⁻¹ * (r ^ (n + 1) - r ^ n))) • v := by
  unfold electricFieldScaled modeHistory
  rw [map_zero, sub_zero]
  funext e
  simp only [Pi.neg_apply, Pi.smul_apply, Pi.sub_apply, smul_eq_mul]
  ring

theorem modeHistory_ampere (h : ℝ) (hh : h ≠ 0) (v : Fin 30 → ℝ) (lam : ℝ)
    (hv : localMaxwellOperator v = lam • v) (hμ : 4 < h ^ 2 * lam) :
    AmpereEvolutionScaled h (modeHistory (growthRoot (h ^ 2 * lam)) v)
      (fun _ ↦ 0) (fun _ ↦ 0) := by
  intro n
  have hroot := growthRoot_root (h ^ 2 * lam) hμ
  set r := growthRoot (h ^ 2 * lam) with hr
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  have hlm : faceCodifferential (magneticField (modeHistory r v) (n + 1)) =
      r ^ (n + 1) • (lam • v) := by
    show localMaxwellOperator (r ^ (n + 1) • v) = _
    rw [map_smul, hv]
  rw [hlm, modeHistory_electricField, modeHistory_electricField]
  funext e
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul, sub_zero]
  linear_combination (-(h⁻¹ * r ^ n * v e)) * hroot +
    (h * lam * r ^ (n + 1) * v e) * hinv

theorem modeHistory_energy (h r : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) n) =
      (r ^ 2) ^ n *
        realSeamEnergy (electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) 0) := by
  rw [modeHistory_electricField, modeHistory_electricField, seamEnergy_smul,
    seamEnergy_smul]
  ring

theorem modeHistory_energy_zero_pos (h : ℝ) (hh : h ≠ 0) (r : ℝ) (hr : r < -1)
    (v : Fin 30 → ℝ) (hv : v ≠ 0) :
    0 < realSeamEnergy (electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) 0) := by
  rw [modeHistory_electricField, seamEnergy_smul]
  apply mul_pos
  · have hc : -(h⁻¹ * (r ^ (0 + 1) - r ^ 0)) ≠ 0 := by
      simp only [zero_add, pow_one, pow_zero]
      exact neg_ne_zero.mpr (mul_ne_zero (inv_ne_zero hh) (sub_ne_zero.mpr (by linarith)))
    exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hc))
  · exact seamEnergy_pos_of_ne_zero v hv

theorem modeHistory_unbounded (h : ℝ) (hh : h ≠ 0) (r : ℝ) (hr : r < -1)
    (v : Fin 30 → ℝ) (hv : v ≠ 0) (M : ℝ) :
    ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) n) := by
  have hE0 := modeHistory_energy_zero_pos h hh r hr v hv
  have hr2 := one_lt_sq_of_lt_neg_one r hr
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt
    (M / realSeamEnergy (electricFieldScaled h (modeHistory r v) (fun _ ↦ 0) 0)) hr2
  refine ⟨n, ?_⟩
  rw [modeHistory_energy]
  rwa [div_lt_iff₀ hE0] at hn

/-- **The unstable mode.**  For an eigenvector `CᵀC v = λ v`, `v ≠ 0`, with
`h² λ > 4`, the history `r ^ n • v` with `r = growthRoot (h² λ) < -1` solves
the zero-current evolution in the gauge `φ = 0`; its electric seam energy is
`(r²)^n` times a positive initial value and is unbounded in `n`. -/
theorem unstable_mode (h : ℝ) (hh : h ≠ 0) (v : Fin 30 → ℝ) (hv : v ≠ 0) (lam : ℝ)
    (hev : localMaxwellOperator v = lam • v) (hμ : 4 < h ^ 2 * lam) :
    AmpereEvolutionScaled h (modeHistory (growthRoot (h ^ 2 * lam)) v)
        (fun _ ↦ 0) (fun _ ↦ 0) ∧
      growthRoot (h ^ 2 * lam) < -1 ∧
      (∀ n, realSeamEnergy (electricFieldScaled h
          (modeHistory (growthRoot (h ^ 2 * lam)) v) (fun _ ↦ 0) n) =
        (growthRoot (h ^ 2 * lam) ^ 2) ^ n * realSeamEnergy (electricFieldScaled h
          (modeHistory (growthRoot (h ^ 2 * lam)) v) (fun _ ↦ 0) 0)) ∧
      0 < realSeamEnergy (electricFieldScaled h
          (modeHistory (growthRoot (h ^ 2 * lam)) v) (fun _ ↦ 0) 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h
          (modeHistory (growthRoot (h ^ 2 * lam)) v) (fun _ ↦ 0) n) :=
  ⟨modeHistory_ampere h hh v lam hev hμ,
    growthRoot_lt_neg_one _ hμ,
    modeHistory_energy h _ v,
    modeHistory_energy_zero_pos h hh _ (growthRoot_lt_neg_one _ hμ) v hv,
    modeHistory_unbounded h hh _ (growthRoot_lt_neg_one _ hμ) v hv⟩

/-! ## (D) Committed-carrier spectral facts -/

theorem faceIncidenceR_cases (f : Fin 20) (e : Fin 30) :
    faceIncidenceR f e = 0 ∨ faceIncidenceR f e = 1 ∨ faceIncidenceR f e = -1 := by
  unfold faceIncidenceR faceIncidenceZ
  split_ifs <;> simp

theorem faceIncidenceR_cube (f : Fin 20) (e : Fin 30) :
    faceIncidenceR f e ^ 2 * faceIncidenceR f e = faceIncidenceR f e := by
  rcases faceIncidenceR_cases f e with h | h | h <;> rw [h] <;> norm_num

theorem faceIncidenceR_sq_sq (f : Fin 20) (e : Fin 30) :
    (faceIncidenceR f e ^ 2) ^ 2 = faceIncidenceR f e ^ 2 := by
  rcases faceIncidenceR_cases f e with h | h | h <;> rw [h] <;> norm_num

/-- Three seams per face: each face row of the committed incidence has
squared sum three. -/
theorem faceIncidence_sq_row_sum_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * faceIncidenceZ f e) = 3 := by
  decide

theorem faceIncidence_sq_row_sum (f : Fin 20) :
    (∑ e : Fin 30, faceIncidenceR f e ^ 2) = 3 := by
  have h := faceIncidence_sq_row_sum_Z f
  have hc : (∑ e : Fin 30, faceIncidenceR f e ^ 2) =
      ((∑ e : Fin 30, faceIncidenceZ f e * faceIncidenceZ f e : ℤ) : ℝ) := by
    rw [Int.cast_sum]
    exact Finset.sum_congr rfl fun e _ ↦ by
      rw [Int.cast_mul]
      unfold faceIncidenceR
      ring
  rw [hc, h]
  norm_num

/-- Two faces per seam: each seam column of the committed incidence has
squared sum two (the committed diagonal of the local Hessian). -/
theorem faceIncidence_sq_col_sum (e : Fin 30) :
    (∑ f : Fin 20, faceIncidenceR f e ^ 2) = 2 := by
  have h := localKineticZ_diagonal_two e
  unfold localKineticZ at h
  have hc : (∑ f : Fin 20, faceIncidenceR f e ^ 2) =
      ((∑ f : Fin 20, faceIncidenceZ f e * faceIncidenceZ f e : ℤ) : ℝ) := by
    rw [Int.cast_sum]
    exact Finset.sum_congr rfl fun f _ ↦ by
      rw [Int.cast_mul]
      unfold faceIncidenceR
      ring
  rw [hc, h]
  norm_num

/-- **D1.**  Elementary Courant constant on the committed carrier:
`‖C v‖² ≤ 6 ‖v‖²`, from Cauchy-Schwarz on each three-seam face row and
the two faces per seam. -/
theorem faceEnergy_curvature_le_six (v : Fin 30 → ℝ) :
    faceEnergy (faceCurvature v) ≤ 6 * realSeamEnergy v := by
  unfold faceEnergy
  have hface : ∀ f : Fin 20, faceCurvature v f ^ 2 ≤
      3 * ∑ e : Fin 30, faceIncidenceR f e ^ 2 * v e ^ 2 := by
    intro f
    rw [faceCurvature_apply]
    have hrw : (∑ e : Fin 30, faceIncidenceR f e * v e) =
        ∑ e : Fin 30, faceIncidenceR f e ^ 2 * (faceIncidenceR f e * v e) := by
      refine Finset.sum_congr rfl fun e _ ↦ ?_
      rw [← mul_assoc, faceIncidenceR_cube]
    rw [hrw]
    have hcs : (∑ e : Fin 30, faceIncidenceR f e ^ 2 * (faceIncidenceR f e * v e)) ^ 2 ≤
        (∑ e : Fin 30, (faceIncidenceR f e ^ 2) ^ 2) *
          ∑ e : Fin 30, (faceIncidenceR f e * v e) ^ 2 :=
      Finset.sum_mul_sq_le_sq_mul_sq _ _ _
    have h3 : (∑ e : Fin 30, (faceIncidenceR f e ^ 2) ^ 2) = 3 := by
      rw [Finset.sum_congr rfl fun e _ ↦ faceIncidenceR_sq_sq f e]
      exact faceIncidence_sq_row_sum f
    have h4 : (∑ e : Fin 30, (faceIncidenceR f e * v e) ^ 2) =
        ∑ e : Fin 30, faceIncidenceR f e ^ 2 * v e ^ 2 :=
      Finset.sum_congr rfl fun e _ ↦ by ring
    rw [h3, h4] at hcs
    exact hcs
  calc (∑ f : Fin 20, faceCurvature v f ^ 2)
      ≤ ∑ f : Fin 20, 3 * ∑ e : Fin 30, faceIncidenceR f e ^ 2 * v e ^ 2 :=
        Finset.sum_le_sum fun f _ ↦ hface f
    _ = 3 * ∑ e : Fin 30, (∑ f : Fin 20, faceIncidenceR f e ^ 2) * v e ^ 2 := by
        rw [← Finset.mul_sum, Finset.sum_comm]
        congr 1
        exact Finset.sum_congr rfl fun e _ ↦ by rw [Finset.sum_mul]
    _ = 3 * ∑ e : Fin 30, 2 * v e ^ 2 := by
        congr 1
        exact Finset.sum_congr rfl fun e _ ↦ by rw [faceIncidence_sq_col_sum]
    _ = 6 * realSeamEnergy v := by
        unfold realSeamEnergy
        rw [← Finset.mul_sum]
        ring

/-- The committed carrier satisfies the Courant hypothesis with `Λ = 6`. -/
theorem committed_courant : CourantBound 6 := faceEnergy_curvature_le_six

/-- Cast transport: an integer curvature identity lifts to the real face
curvature. -/
theorem faceCurvature_cast (vZ : Fin 30 → ℤ) (wZ : Fin 20 → ℤ)
    (hZ : ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * vZ e) = wZ f) :
    faceCurvature (fun e ↦ (vZ e : ℝ)) = fun f ↦ (wZ f : ℝ) := by
  funext f
  rw [faceCurvature_apply, ← hZ f, Int.cast_sum]
  exact Finset.sum_congr rfl fun e _ ↦ by
    rw [Int.cast_mul]
    rfl

theorem faceCodifferential_cast (wZ : Fin 20 → ℤ) (uZ : Fin 30 → ℤ)
    (hZ : ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * wZ f) = uZ e) :
    faceCodifferential (fun f ↦ (wZ f : ℝ)) = fun e ↦ (uZ e : ℝ) := by
  funext e
  rw [faceCodifferential_apply, ← hZ e, Int.cast_sum]
  exact Finset.sum_congr rfl fun f _ ↦ by
    rw [Int.cast_mul]
    rfl

theorem localMaxwellOperator_cast (vZ uZ : Fin 30 → ℤ) (wZ : Fin 20 → ℤ)
    (h1 : ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * vZ e) = wZ f)
    (h2 : ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * wZ f) = uZ e) :
    localMaxwellOperator (fun e ↦ (vZ e : ℝ)) = fun e ↦ (uZ e : ℝ) := by
  show faceCodifferential (faceCurvature _) = _
  rw [faceCurvature_cast vZ wZ h1, faceCodifferential_cast wZ uZ h2]

/-! ### D2: an integer eigenvector of `CᵀC` with eigenvalue `5` -/

def fiveModeZ : Fin 30 → ℤ :=
  ![0, 0, 0, -1, 1, 1, -1, -2, 2, 0, 1, 0, 0, -1, 0, -2, -1, 2, -2, 0, 1, 1, -2, 0,
    -1, -1, 0, 1, 0, 0]

def fiveFaceZ : Fin 20 → ℤ :=
  ![1, 1, 1, 1, -4, -4, -4, 6, 1, 1, 1, 1, 6, -4, -4, 1, -4, 1, 1, 1]

set_option maxRecDepth 16384 in
theorem fiveMode_curvature_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * fiveModeZ e) = fiveFaceZ f := by
  decide

set_option maxRecDepth 16384 in
theorem fiveMode_codifferential_Z :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * fiveFaceZ f) = 5 * fiveModeZ e := by
  decide

/-- The real cast of the integer mode. -/
def fiveMode : Fin 30 → ℝ := fun e ↦ (fiveModeZ e : ℝ)

/-- **D2.**  `CᵀC fiveMode = 5 • fiveMode`, kernel-checked from the
committed incidence. -/
theorem fiveMode_eigen : localMaxwellOperator fiveMode = (5 : ℝ) • fiveMode := by
  have h := localMaxwellOperator_cast fiveModeZ (fun e ↦ 5 * fiveModeZ e) fiveFaceZ
    fiveMode_curvature_Z fiveMode_codifferential_Z
  rw [show fiveMode = (fun e ↦ (fiveModeZ e : ℝ)) from rfl, h]
  funext e
  simp only [Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

theorem fiveMode_ne_zero : fiveMode ≠ 0 := by
  intro h
  have h3 := congrFun h 3
  have hz : fiveModeZ 3 = -1 := by decide
  have hc : ((fiveModeZ 3 : ℤ) : ℝ) = 0 := h3
  rw [hz] at hc
  norm_num at hc

/-! ### D3: an eigenvector with entries in `ℤ + ℤ √5` and eigenvalue `3 + √5` -/

def goldenAZ : Fin 30 → ℤ :=
  ![0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 1, -1, 0, -1, 1, 0, 1, -1, 0, 0, 0, -1, 1, 0,
    0, 0, 0, 0, 0, 0]

def goldenBZ : Fin 30 → ℤ :=
  ![0, 0, 0, 0, 0, 2, -2, -1, 1, 2, 1, -1, -2, -1, 1, 2, 1, -1, -2, 2, 0, -1, 1, -2,
    0, 2, 0, -2, 0, 0]

def goldenAFaceZ : Fin 20 → ℤ :=
  ![0, 0, 0, 0, 0, -2, -2, 2, -2, 2, -2, 2, -2, 2, 0, 0, 2, 0, 0, 0]

def goldenBFaceZ : Fin 20 → ℤ :=
  ![2, 2, 2, 2, 2, -4, -4, 4, -4, 4, -4, 4, -4, 4, -2, -2, 4, -2, -2, -2]

set_option maxRecDepth 16384 in
theorem goldenA_curvature_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * goldenAZ e) = goldenAFaceZ f := by
  decide

set_option maxRecDepth 16384 in
theorem goldenB_curvature_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * goldenBZ e) = goldenBFaceZ f := by
  decide

set_option maxRecDepth 16384 in
theorem goldenA_codifferential_Z :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * goldenAFaceZ f) =
      3 * goldenAZ e + goldenBZ e := by
  decide

set_option maxRecDepth 16384 in
theorem goldenB_codifferential_Z :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * goldenBFaceZ f) =
      5 * goldenAZ e + 3 * goldenBZ e := by
  decide

def goldenA : Fin 30 → ℝ := fun e ↦ (goldenAZ e : ℝ)
def goldenB : Fin 30 → ℝ := fun e ↦ (goldenBZ e : ℝ)

/-- The golden mode `√5 • goldenA + goldenB`. -/
def goldenMode : Fin 30 → ℝ := Real.sqrt 5 • goldenA + goldenB

theorem goldenA_eigen_pair :
    localMaxwellOperator goldenA = (3 : ℝ) • goldenA + goldenB := by
  have h := localMaxwellOperator_cast goldenAZ (fun e ↦ 3 * goldenAZ e + goldenBZ e)
    goldenAFaceZ goldenA_curvature_Z goldenA_codifferential_Z
  rw [show goldenA = (fun e ↦ (goldenAZ e : ℝ)) from rfl,
    show goldenB = (fun e ↦ (goldenBZ e : ℝ)) from rfl, h]
  funext e
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

theorem goldenB_eigen_pair :
    localMaxwellOperator goldenB = (5 : ℝ) • goldenA + (3 : ℝ) • goldenB := by
  have h := localMaxwellOperator_cast goldenBZ
    (fun e ↦ 5 * goldenAZ e + 3 * goldenBZ e)
    goldenBFaceZ goldenB_curvature_Z goldenB_codifferential_Z
  rw [show goldenA = (fun e ↦ (goldenAZ e : ℝ)) from rfl,
    show goldenB = (fun e ↦ (goldenBZ e : ℝ)) from rfl, h]
  funext e
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

/-- **D3.**  `CᵀC goldenMode = (3 + √5) • goldenMode`, kernel-checked from
the committed incidence and `√5 · √5 = 5`. -/
theorem goldenMode_eigen :
    localMaxwellOperator goldenMode = (3 + Real.sqrt 5) • goldenMode := by
  unfold goldenMode
  rw [map_add, map_smul, goldenA_eigen_pair, goldenB_eigen_pair]
  funext e
  have hs : Real.sqrt 5 * Real.sqrt 5 = 5 := Real.mul_self_sqrt (by norm_num)
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  linear_combination (-(goldenA e)) * hs

theorem goldenMode_ne_zero : goldenMode ≠ 0 := by
  intro h
  have h5 := congrFun h 5
  have ha : goldenAZ 5 = 0 := by decide
  have hb : goldenBZ 5 = 2 := by decide
  have hc : Real.sqrt 5 * ((goldenAZ 5 : ℤ) : ℝ) + ((goldenBZ 5 : ℤ) : ℝ) = 0 := h5
  rw [ha, hb] at hc
  norm_num at hc

/-! ### Necessary conditions from the two committed eigenvectors -/

/-- For `h² · 5 > 4` the committed carrier carries an unbounded zero-current
solution: boundedness of every initial datum forces `h² ≤ 4/5`. -/
theorem instability_above_five (h : ℝ) (hh : h ≠ 0) (hc : 4 < h ^ 2 * 5) :
    ∃ A : ℕ → Fin 30 → ℝ, AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h A (fun _ ↦ 0) n) :=
  ⟨_, (unstable_mode h hh fiveMode fiveMode_ne_zero 5 fiveMode_eigen hc).1,
    (unstable_mode h hh fiveMode fiveMode_ne_zero 5 fiveMode_eigen hc).2.2.2.2⟩

/-- For `h² (3 + √5) > 4` the committed carrier carries an unbounded
zero-current solution: boundedness of every initial datum forces
`h² ≤ 4/(3 + √5)`, the bound `|h| ≤ 2/√(3 + √5)`. -/
theorem instability_above_golden (h : ℝ) (hh : h ≠ 0)
    (hc : 4 < h ^ 2 * (3 + Real.sqrt 5)) :
    ∃ A : ℕ → Fin 30 → ℝ, AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h A (fun _ ↦ 0) n) :=
  ⟨_, (unstable_mode h hh goldenMode goldenMode_ne_zero _ goldenMode_eigen hc).1,
    (unstable_mode h hh goldenMode goldenMode_ne_zero _ goldenMode_eigen hc).2.2.2.2⟩

/-- **Instability of the declared unit step.**  The committed
`AmpereEvolution` of `TemporalMaxwellEvolution` admits a zero-current,
zero-potential solution whose electric seam energy is unbounded. -/
theorem unit_step_instability :
    ∃ A : ℕ → Fin 30 → ℝ, AmpereEvolution A (fun _ ↦ 0) (fun _ ↦ 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricField A (fun _ ↦ 0) n) := by
  obtain ⟨A, hA, hM⟩ := instability_above_five 1 one_ne_zero (by norm_num)
  refine ⟨A, (ampereEvolutionScaled_one_iff A _ _).mp hA, fun M ↦ ?_⟩
  obtain ⟨n, hn⟩ := hM M
  exact ⟨n, by rwa [electricFieldScaled_one] at hn⟩

/-- **Sufficient condition on the committed carrier.**  For `h² < 2/3`
and zero current every solution is bounded by the explicit multiples of
the initial staggered form. -/
theorem committed_carrier_stability (h : ℝ) (hh : h ≠ 0) (hc : h ^ 2 * 6 < 4)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ (fun _ ↦ 0)) (n : ℕ) :
    0 ≤ fieldEnergyScaled h A φ n ∧
      realSeamEnergy (electricFieldScaled h A φ n) ≤
        8 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * 6) ∧
      faceEnergy (magneticField A n) ≤
        16 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * 6) :=
  ⟨fieldEnergyScaled_nonneg h hh 6 committed_courant hc.le A φ n,
    stability_certificate h hh 6 (by norm_num) committed_courant hc A φ hAmp n⟩

/-! ### D4: the sharp Courant constant `3 + √5` by explicit spectral projectors

The face normal matrix `N = C Cᵀ` is a twenty-by-twenty integer matrix.
Its five spectral projectors (eigenvalues `0`, `2`, `3`, `5`, and the
six-dimensional sector `N² - 6N + 4 = 0`, the eigenvalues `3 ± √5` of
multiplicity three each) are explicit rational matrices;
their integer multiples are listed below and every identity used is a
kernel `decide` on the committed tables.  No matrix power is computed in
Lean: each product identity compares one product of two literals against a
third literal. -/

/-- The committed face normal matrix `N = C Cᵀ`: diagonal three, signed dual adjacency. -/
def faceNormalZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![3, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      ![-1, 3, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      ![-1, 0, 3, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      ![0, -1, 0, 3, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      ![0, 0, -1, -1, 3, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
      ![-1, 0, 0, 0, 0, 3, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      ![0, -1, 0, 0, 0, 0, 3, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
      ![0, 0, 0, 0, 0, -1, -1, 3, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
      ![0, 0, -1, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
      ![0, 0, 0, 0, 0, -1, 0, 0, -1, 3, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
      ![0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, 0, -1, 0, 0, 0],
      ![0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 3, 0, 0, 0, 0, 0, -1, 0, 0],
      ![0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, -1, 0, 0, 0],
      ![0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 3, 0, 0, 0, 0, -1, 0],
      ![0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 3, -1, 0, -1, 0, 0],
      ![0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 3, 0, 0, -1, 0],
      ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 3, 0, 0, -1],
      ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 3, 0, -1],
      ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 3, -1],
      ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 3]]

/-- `20` times the spectral projector of `N` for eigenvalue `0`. -/
def projZeroZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

/-- `12` times the spectral projector of `N` for eigenvalue `2`. -/
def projTwoZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![3, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 3],
      ![1, 3, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 3, 1],
      ![1, -1, 3, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 3, -1, 1],
      ![-1, 1, -1, 3, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 3, -1, -1, 1, -1],
      ![-1, -1, 1, 1, 3, -1, -1, 1, -1, -1, -1, -1, 1, -1, 3, 1, -1, 1, -1, -1],
      ![1, -1, -1, -1, -1, 3, -1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 3, -1, -1, 1],
      ![-1, 1, -1, -1, -1, -1, 3, 1, 1, -1, -1, 1, 1, 3, -1, -1, -1, -1, 1, -1],
      ![-1, -1, -1, -1, 1, 1, 1, 3, -1, -1, -1, -1, 3, 1, 1, -1, 1, -1, -1, -1],
      ![-1, -1, 1, -1, -1, -1, 1, -1, 3, 1, 1, 3, -1, 1, -1, -1, -1, 1, -1, -1],
      ![-1, -1, -1, 1, -1, 1, -1, -1, 1, 3, 3, 1, -1, -1, -1, 1, 1, -1, -1, -1],
      ![-1, -1, -1, 1, -1, 1, -1, -1, 1, 3, 3, 1, -1, -1, -1, 1, 1, -1, -1, -1],
      ![-1, -1, 1, -1, -1, -1, 1, -1, 3, 1, 1, 3, -1, 1, -1, -1, -1, 1, -1, -1],
      ![-1, -1, -1, -1, 1, 1, 1, 3, -1, -1, -1, -1, 3, 1, 1, -1, 1, -1, -1, -1],
      ![-1, 1, -1, -1, -1, -1, 3, 1, 1, -1, -1, 1, 1, 3, -1, -1, -1, -1, 1, -1],
      ![-1, -1, 1, 1, 3, -1, -1, 1, -1, -1, -1, -1, 1, -1, 3, 1, -1, 1, -1, -1],
      ![-1, 1, -1, 3, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 3, -1, -1, 1, -1],
      ![1, -1, -1, -1, -1, 3, -1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 3, -1, -1, 1],
      ![1, -1, 3, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 3, -1, 1],
      ![1, 3, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 3, 1],
      ![3, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 3]]

/-- `10` times the spectral projector of `N` for eigenvalue `3`. -/
def projThreeZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![2, 0, 0, -1, -1, 0, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -2],
      ![0, 2, -1, 0, -1, -1, 0, -1, 1, 1, -1, -1, 1, 0, 1, 0, 1, 1, -2, 0],
      ![0, -1, 2, -1, 0, -1, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, 1, -2, 1, 0],
      ![-1, 0, -1, 2, 0, 1, -1, 1, 1, 0, 0, -1, -1, 1, 0, -2, -1, 1, 0, 1],
      ![-1, -1, 0, 0, 2, 1, 1, 0, -1, 1, -1, 1, 0, -1, -2, 0, -1, 0, 1, 1],
      ![0, -1, -1, 1, 1, 2, -1, 0, -1, 0, 0, 1, 0, 1, -1, -1, -2, 1, 1, 0],
      ![-1, 0, 1, -1, 1, -1, 2, 0, 0, 1, -1, 0, 0, -2, -1, 1, 1, -1, 0, 1],
      ![-1, -1, 1, 1, 0, 0, 0, 2, 1, -1, 1, -1, -2, 0, 0, -1, 0, -1, 1, 1],
      ![-1, 1, 0, 1, -1, -1, 0, 1, 2, 0, 0, -2, -1, 0, 1, -1, 1, 0, -1, 1],
      ![-1, 1, -1, 0, 1, 0, 1, -1, 0, 2, -2, 0, 1, -1, -1, 0, 0, 1, -1, 1],
      ![1, -1, 1, 0, -1, 0, -1, 1, 0, -2, 2, 0, -1, 1, 1, 0, 0, -1, 1, -1],
      ![1, -1, 0, -1, 1, 1, 0, -1, -2, 0, 0, 2, 1, 0, -1, 1, -1, 0, 1, -1],
      ![1, 1, -1, -1, 0, 0, 0, -2, -1, 1, -1, 1, 2, 0, 0, 1, 0, 1, -1, -1],
      ![1, 0, -1, 1, -1, 1, -2, 0, 0, -1, 1, 0, 0, 2, 1, -1, -1, 1, 0, -1],
      ![1, 1, 0, 0, -2, -1, -1, 0, 1, -1, 1, -1, 0, 1, 2, 0, 1, 0, -1, -1],
      ![1, 0, 1, -2, 0, -1, 1, -1, -1, 0, 0, 1, 1, -1, 0, 2, 1, -1, 0, -1],
      ![0, 1, 1, -1, -1, -2, 1, 0, 1, 0, 0, -1, 0, -1, 1, 1, 2, -1, -1, 0],
      ![0, 1, -2, 1, 0, 1, -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, -1, 2, -1, 0],
      ![0, -2, 1, 0, 1, 1, 0, 1, -1, -1, 1, 1, -1, 0, -1, 0, -1, -1, 2, 0],
      ![-2, 0, 0, 1, 1, 0, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 2]]

/-- `30` times the spectral projector of `N` for eigenvalue `5`. -/
def projFiveZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![6, -4, -4, 1, 1, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -4, -4, -4, 6],
      ![-4, 6, 1, -4, 1, 1, -4, 1, 1, 1, 1, 1, 1, -4, 1, -4, 1, 1, 6, -4],
      ![-4, 1, 6, 1, -4, 1, 1, 1, -4, 1, 1, -4, 1, 1, -4, 1, 1, 6, 1, -4],
      ![1, -4, 1, 6, -4, 1, 1, 1, 1, -4, -4, 1, 1, 1, -4, 6, 1, 1, -4, 1],
      ![1, 1, -4, -4, 6, 1, 1, -4, 1, 1, 1, 1, -4, 1, 6, -4, 1, -4, 1, 1],
      ![-4, 1, 1, 1, 1, 6, 1, -4, 1, -4, -4, 1, -4, 1, 1, 1, 6, 1, 1, -4],
      ![1, -4, 1, 1, 1, 1, 6, -4, -4, 1, 1, -4, -4, 6, 1, 1, 1, 1, -4, 1],
      ![1, 1, 1, 1, -4, -4, -4, 6, 1, 1, 1, 1, 6, -4, -4, 1, -4, 1, 1, 1],
      ![1, 1, -4, 1, 1, 1, -4, 1, 6, -4, -4, 6, 1, -4, 1, 1, 1, -4, 1, 1],
      ![1, 1, 1, -4, 1, -4, 1, 1, -4, 6, 6, -4, 1, 1, 1, -4, -4, 1, 1, 1],
      ![1, 1, 1, -4, 1, -4, 1, 1, -4, 6, 6, -4, 1, 1, 1, -4, -4, 1, 1, 1],
      ![1, 1, -4, 1, 1, 1, -4, 1, 6, -4, -4, 6, 1, -4, 1, 1, 1, -4, 1, 1],
      ![1, 1, 1, 1, -4, -4, -4, 6, 1, 1, 1, 1, 6, -4, -4, 1, -4, 1, 1, 1],
      ![1, -4, 1, 1, 1, 1, 6, -4, -4, 1, 1, -4, -4, 6, 1, 1, 1, 1, -4, 1],
      ![1, 1, -4, -4, 6, 1, 1, -4, 1, 1, 1, 1, -4, 1, 6, -4, 1, -4, 1, 1],
      ![1, -4, 1, 6, -4, 1, 1, 1, 1, -4, -4, 1, 1, 1, -4, 6, 1, 1, -4, 1],
      ![-4, 1, 1, 1, 1, 6, 1, -4, 1, -4, -4, 1, -4, 1, 1, 1, 6, 1, 1, -4],
      ![-4, 1, 6, 1, -4, 1, 1, 1, -4, 1, 1, -4, 1, 1, -4, 1, 1, 6, 1, -4],
      ![-4, 6, 1, -4, 1, 1, -4, 1, 1, 1, 1, 1, 1, -4, 1, -4, 1, 1, 6, -4],
      ![6, -4, -4, 1, 1, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -4, -4, -4, 6]]

/-- `10` times the spectral projector of `N` onto the sector `N² - 6N + 4 = 0`. -/
def projGoldenZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![3, 0, 0, 1, 1, 0, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -3],
      ![0, 3, 1, 0, 1, 1, 0, 1, -1, -1, 1, 1, -1, 0, -1, 0, -1, -1, -3, 0],
      ![0, 1, 3, 1, 0, 1, -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, -1, -3, -1, 0],
      ![1, 0, 1, 3, 0, -1, 1, -1, -1, 0, 0, 1, 1, -1, 0, -3, 1, -1, 0, -1],
      ![1, 1, 0, 0, 3, -1, -1, 0, 1, -1, 1, -1, 0, 1, -3, 0, 1, 0, -1, -1],
      ![0, 1, 1, -1, -1, 3, 1, 0, 1, 0, 0, -1, 0, -1, 1, 1, -3, -1, -1, 0],
      ![1, 0, -1, 1, -1, 1, 3, 0, 0, -1, 1, 0, 0, -3, 1, -1, -1, 1, 0, -1],
      ![1, 1, -1, -1, 0, 0, 0, 3, -1, 1, -1, 1, -3, 0, 0, 1, 0, 1, -1, -1],
      ![1, -1, 0, -1, 1, 1, 0, -1, 3, 0, 0, -3, 1, 0, -1, 1, -1, 0, 1, -1],
      ![1, -1, 1, 0, -1, 0, -1, 1, 0, 3, -3, 0, -1, 1, 1, 0, 0, -1, 1, -1],
      ![-1, 1, -1, 0, 1, 0, 1, -1, 0, -3, 3, 0, 1, -1, -1, 0, 0, 1, -1, 1],
      ![-1, 1, 0, 1, -1, -1, 0, 1, -3, 0, 0, 3, -1, 0, 1, -1, 1, 0, -1, 1],
      ![-1, -1, 1, 1, 0, 0, 0, -3, 1, -1, 1, -1, 3, 0, 0, -1, 0, -1, 1, 1],
      ![-1, 0, 1, -1, 1, -1, -3, 0, 0, 1, -1, 0, 0, 3, -1, 1, 1, -1, 0, 1],
      ![-1, -1, 0, 0, -3, 1, 1, 0, -1, 1, -1, 1, 0, -1, 3, 0, -1, 0, 1, 1],
      ![-1, 0, -1, -3, 0, 1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 3, -1, 1, 0, 1],
      ![0, -1, -1, 1, 1, -3, -1, 0, -1, 0, 0, 1, 0, 1, -1, -1, 3, 1, 1, 0],
      ![0, -1, -3, -1, 0, -1, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, 1, 3, 1, 0],
      ![0, -3, -1, 0, -1, -1, 0, -1, 1, 1, -1, -1, 1, 0, 1, 0, 1, 1, 3, 0],
      ![-3, 0, 0, -1, -1, 0, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 3]]

/-- The product `N · projGoldenZ`, the image of the golden sector under ten times `N`. -/
def goldenImageZ : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![9, -5, -5, 3, 3, -5, 3, 3, 3, 3, -3, -3, -3, -3, -3, -3, 5, 5, 5, -9],
      ![-5, 9, 3, -5, 3, 3, -5, 3, -3, -3, 3, 3, -3, 5, -3, 5, -3, -3, -9, 5],
      ![-5, 3, 9, 3, -5, 3, -3, -3, -5, 3, -3, 5, 3, 3, 5, -3, -3, -9, -3, 5],
      ![3, -5, 3, 9, -5, -3, 3, -3, -3, 5, -5, 3, 3, -3, 5, -9, 3, -3, 5, -3],
      ![3, 3, -5, -5, 9, -3, -3, 5, 3, -3, 3, -3, -5, 3, -9, 5, 3, 5, -3, -3],
      ![-5, 3, 3, -3, -3, 9, 3, -5, 3, -5, 5, -3, 5, -3, 3, 3, -9, -3, -3, 5],
      ![3, -5, -3, 3, -3, 3, 9, -5, 5, -3, 3, -5, 5, -9, 3, -3, -3, 3, 5, -3],
      ![3, 3, -3, -3, 5, -5, -5, 9, -3, 3, -3, 3, -9, 5, -5, 3, 5, 3, -3, -3],
      ![3, -3, -5, -3, 3, 3, 5, -3, 9, -5, 5, -9, 3, -5, -3, 3, -3, 5, 3, -3],
      ![3, -3, 3, 5, -3, -5, -3, 3, -5, 9, -9, 5, -3, 3, 3, -5, 5, -3, 3, -3],
      ![-3, 3, -3, -5, 3, 5, 3, -3, 5, -9, 9, -5, 3, -3, -3, 5, -5, 3, -3, 3],
      ![-3, 3, 5, 3, -3, -3, -5, 3, -9, 5, -5, 9, -3, 5, 3, -3, 3, -5, -3, 3],
      ![-3, -3, 3, 3, -5, 5, 5, -9, 3, -3, 3, -3, 9, -5, 5, -3, -5, -3, 3, 3],
      ![-3, 5, 3, -3, 3, -3, -9, 5, -5, 3, -3, 5, -5, 9, -3, 3, 3, -3, -5, 3],
      ![-3, -3, 5, 5, -9, 3, 3, -5, -3, 3, -3, 3, 5, -3, 9, -5, -3, -5, 3, 3],
      ![-3, 5, -3, -9, 5, 3, -3, 3, 3, -5, 5, -3, -3, 3, -5, 9, -3, 3, -5, 3],
      ![5, -3, -3, 3, 3, -9, -3, 5, -3, 5, -5, 3, -5, 3, -3, -3, 9, 3, 3, -5],
      ![5, -3, -9, -3, 5, -3, 3, 3, 5, -3, 3, -5, -3, -3, -5, 3, 3, 9, 3, -5],
      ![5, -9, -3, 5, -3, -3, 5, -3, 3, 3, -3, -3, 3, -5, 3, -5, 3, 3, 9, -5],
      ![-9, 5, 5, -3, -3, 5, -3, -3, -3, -3, 3, 3, 3, 3, 3, 3, -5, -5, -5, 9]]

set_option maxRecDepth 16384 in
/-- The face normal table is the committed incidence product `C Cᵀ`. -/
theorem faceNormal_eq_incidence :
    ∀ f g : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * faceIncidenceZ g e) =
      faceNormalZ f g := by
  decide

theorem faceNormalZ_symm : ∀ i j : Fin 20, faceNormalZ i j = faceNormalZ j i := by
  decide

set_option maxRecDepth 16384 in
theorem projZero_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, projZeroZ i k * projZeroZ k j) = 20 * projZeroZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projZero_eigen :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * projZeroZ k j) = 0 * projZeroZ i j := by
  decide

theorem projZero_symm : ∀ i j : Fin 20, projZeroZ i j = projZeroZ j i := by
  decide

set_option maxRecDepth 16384 in
theorem projTwo_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, projTwoZ i k * projTwoZ k j) = 12 * projTwoZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projTwo_eigen :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * projTwoZ k j) = 2 * projTwoZ i j := by
  decide

theorem projTwo_symm : ∀ i j : Fin 20, projTwoZ i j = projTwoZ j i := by
  decide

set_option maxRecDepth 16384 in
theorem projThree_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, projThreeZ i k * projThreeZ k j) = 10 * projThreeZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projThree_eigen :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * projThreeZ k j) = 3 * projThreeZ i j := by
  decide

theorem projThree_symm : ∀ i j : Fin 20, projThreeZ i j = projThreeZ j i := by
  decide

set_option maxRecDepth 16384 in
theorem projFive_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, projFiveZ i k * projFiveZ k j) = 30 * projFiveZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projFive_eigen :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * projFiveZ k j) = 5 * projFiveZ i j := by
  decide

theorem projFive_symm : ∀ i j : Fin 20, projFiveZ i j = projFiveZ j i := by
  decide

set_option maxRecDepth 16384 in
theorem projGolden_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, projGoldenZ i k * projGoldenZ k j) =
      10 * projGoldenZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projGolden_image :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * projGoldenZ k j) = goldenImageZ i j := by
  decide

set_option maxRecDepth 16384 in
theorem projGolden_commute :
    ∀ i j : Fin 20, (∑ k : Fin 20, projGoldenZ i k * faceNormalZ k j) = goldenImageZ i j := by
  decide

set_option maxRecDepth 16384 in
/-- The golden sector is annihilated by `N² - 6N + 4`. -/
theorem projGolden_quadratic :
    ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * goldenImageZ k j) -
      6 * goldenImageZ i j + 4 * projGoldenZ i j = 0 := by
  decide

theorem projGolden_symm : ∀ i j : Fin 20, projGoldenZ i j = projGoldenZ j i := by
  decide

/-- The five projectors resolve the identity:
`(1/20) P₀ + (1/12) P₂ + (1/10) P₃ + (1/30) P₅ + (1/10) P_Q = 1`, cleared to `60`. -/
theorem proj_sum :
    ∀ i j : Fin 20, 3 * projZeroZ i j + 5 * projTwoZ i j + 6 * projThreeZ i j +
      2 * projFiveZ i j + 6 * projGoldenZ i j = if i = j then 60 else 0 := by
  decide

/-! #### Real transport of the integer certificates -/

/-- Real cast of an integer face matrix. -/
def castZ (A : Matrix (Fin 20) (Fin 20) ℤ) : Matrix (Fin 20) (Fin 20) ℝ :=
  Matrix.of fun i j ↦ (A i j : ℝ)

theorem castZ_apply (A : Matrix (Fin 20) (Fin 20) ℤ) (i j : Fin 20) :
    castZ A i j = (A i j : ℝ) := rfl

theorem castZ_mul (A B C : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ)
    (hZ : ∀ i j : Fin 20, (∑ k : Fin 20, A i k * B k j) = d * C i j) :
    castZ A * castZ B = (d : ℝ) • castZ C := by
  ext i j
  rw [Matrix.mul_apply, Matrix.smul_apply, smul_eq_mul, castZ_apply]
  have hc : (∑ k : Fin 20, castZ A i k * castZ B k j) =
      ((∑ k : Fin 20, A i k * B k j : ℤ) : ℝ) := by
    rw [Int.cast_sum]
    exact Finset.sum_congr rfl fun k _ ↦ by
      rw [Int.cast_mul]
      rfl
  rw [hc, hZ i j, Int.cast_mul]

theorem castZ_mul_eq (A B C : Matrix (Fin 20) (Fin 20) ℤ)
    (hZ : ∀ i j : Fin 20, (∑ k : Fin 20, A i k * B k j) = C i j) :
    castZ A * castZ B = castZ C := by
  have h := castZ_mul A B C 1 (fun i j ↦ by rw [hZ i j, one_mul])
  rwa [Int.cast_one, one_smul] at h

theorem castZ_symm (A : Matrix (Fin 20) (Fin 20) ℤ) (hZ : ∀ i j, A i j = A j i) :
    (castZ A)ᵀ = castZ A := by
  ext i j
  rw [Matrix.transpose_apply, castZ_apply, castZ_apply, hZ]

/-- Symmetric matrices move across the face pairing. -/
theorem faceInner_mulVec_symm (M : Matrix (Fin 20) (Fin 20) ℝ) (hM : Mᵀ = M)
    (x y : Fin 20 → ℝ) :
    faceInner (M.mulVec x) y = faceInner x (M.mulVec y) := by
  show (M *ᵥ x) ⬝ᵥ y = x ⬝ᵥ (M *ᵥ y)
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hM]

theorem faceInner_add_right (F G H : Fin 20 → ℝ) :
    faceInner F (G + H) = faceInner F G + faceInner F H := by
  unfold faceInner
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.add_apply]
    ring

theorem faceInner_sub_right (F G H : Fin 20 → ℝ) :
    faceInner F (G - H) = faceInner F G - faceInner F H := by
  unfold faceInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.sub_apply]
    ring

theorem faceInner_smul_right (c : ℝ) (F G : Fin 20 → ℝ) :
    faceInner F (c • G) = c * faceInner F G := by
  rw [faceInner_comm, faceInner_smul_left, faceInner_comm]

theorem faceEnergy_sub (F G : Fin 20 → ℝ) :
    faceEnergy (F - G) = faceEnergy F - 2 * faceInner F G + faceEnergy G := by
  unfold faceEnergy faceInner
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.sub_apply]
    ring

/-- The real face normal matrix. -/
def faceNormalR : Matrix (Fin 20) (Fin 20) ℝ := castZ faceNormalZ

theorem faceNormal_symm : faceNormalRᵀ = faceNormalR :=
  castZ_symm faceNormalZ faceNormalZ_symm

/-- `C Cᵀ` acts on face loads by the face normal matrix. -/
theorem faceNormal_mulVec (w : Fin 20 → ℝ) :
    faceCurvature (faceCodifferential w) = faceNormalR.mulVec w := by
  funext f
  rw [faceCurvature_apply]
  simp only [faceCodifferential_apply]
  have hm : faceNormalR.mulVec w f = ∑ g : Fin 20, (faceNormalZ f g : ℝ) * w g := rfl
  rw [hm]
  calc (∑ e : Fin 30, faceIncidenceR f e * ∑ g : Fin 20, faceIncidenceR g e * w g)
      = ∑ e : Fin 30, ∑ g : Fin 20, faceIncidenceR f e * (faceIncidenceR g e * w g) := by
          refine Finset.sum_congr rfl fun e _ ↦ ?_
          rw [Finset.mul_sum]
    _ = ∑ g : Fin 20, ∑ e : Fin 30, faceIncidenceR f e * (faceIncidenceR g e * w g) :=
          Finset.sum_comm
    _ = ∑ g : Fin 20, (∑ e : Fin 30, faceIncidenceR f e * faceIncidenceR g e) * w g := by
          refine Finset.sum_congr rfl fun g _ ↦ ?_
          rw [Finset.sum_mul]
          exact Finset.sum_congr rfl fun e _ ↦ by ring
    _ = ∑ g : Fin 20, (faceNormalZ f g : ℝ) * w g := by
          refine Finset.sum_congr rfl fun g _ ↦ ?_
          congr 1
          rw [← faceNormal_eq_incidence f g, Int.cast_sum]
          exact Finset.sum_congr rfl fun e _ ↦ by
            rw [Int.cast_mul]
            rfl

/-- A rational projector `d⁻¹ • Q` from its integer multiple. -/
def scaledProj (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ) : Matrix (Fin 20) (Fin 20) ℝ :=
  ((d : ℝ))⁻¹ • castZ Q

theorem scaledProj_symm (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ)
    (hsymm : ∀ i j, Q i j = Q j i) : (scaledProj Q d)ᵀ = scaledProj Q d := by
  unfold scaledProj
  rw [Matrix.transpose_smul, castZ_symm Q hsymm]

theorem scaledProj_idem (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ) (hd : (d : ℝ) ≠ 0)
    (hidem : ∀ i j : Fin 20, (∑ k : Fin 20, Q i k * Q k j) = d * Q i j) :
    scaledProj Q d * scaledProj Q d = scaledProj Q d := by
  unfold scaledProj
  rw [Matrix.smul_mul, Matrix.mul_smul, castZ_mul Q Q Q d hidem, smul_smul, smul_smul]
  congr 1
  field_simp

theorem scaledProj_eigen (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ) (lam : ℤ)
    (heig : ∀ i j : Fin 20, (∑ k : Fin 20, faceNormalZ i k * Q k j) = lam * Q i j) :
    faceNormalR * scaledProj Q d = (lam : ℝ) • scaledProj Q d := by
  unfold scaledProj faceNormalR
  rw [Matrix.mul_smul, castZ_mul faceNormalZ Q Q lam heig, smul_comm]

/-- A symmetric idempotent pairs as an energy. -/
theorem projector_pairing (P : Matrix (Fin 20) (Fin 20) ℝ) (hsymm : Pᵀ = P)
    (hidem : P * P = P) (w : Fin 20 → ℝ) :
    faceInner w (P.mulVec w) = faceEnergy (P.mulVec w) := by
  rw [← faceInner_self_eq_energy, faceInner_mulVec_symm P hsymm, Matrix.mulVec_mulVec,
    hidem]

theorem eigen_pairing (P : Matrix (Fin 20) (Fin 20) ℝ) (lam : ℝ)
    (heig : faceNormalR * P = lam • P) (w : Fin 20 → ℝ) :
    faceInner w (faceNormalR.mulVec (P.mulVec w)) = lam * faceInner w (P.mulVec w) := by
  rw [Matrix.mulVec_mulVec, heig, Matrix.smul_mulVec, faceInner_smul_right]

def projZeroR : Matrix (Fin 20) (Fin 20) ℝ := scaledProj projZeroZ 20
def projTwoR : Matrix (Fin 20) (Fin 20) ℝ := scaledProj projTwoZ 12
def projThreeR : Matrix (Fin 20) (Fin 20) ℝ := scaledProj projThreeZ 10
def projFiveR : Matrix (Fin 20) (Fin 20) ℝ := scaledProj projFiveZ 30
def projGoldenR : Matrix (Fin 20) (Fin 20) ℝ := scaledProj projGoldenZ 10

theorem golden_idem : projGoldenR * projGoldenR = projGoldenR :=
  scaledProj_idem projGoldenZ 10 (by norm_num) projGolden_idem

theorem golden_symm : projGoldenRᵀ = projGoldenR :=
  scaledProj_symm projGoldenZ 10 projGolden_symm

theorem golden_commute : faceNormalR * projGoldenR = projGoldenR * faceNormalR := by
  have h1 : castZ faceNormalZ * castZ projGoldenZ = castZ goldenImageZ :=
    castZ_mul_eq _ _ _ projGolden_image
  have h2 : castZ projGoldenZ * castZ faceNormalZ = castZ goldenImageZ :=
    castZ_mul_eq _ _ _ projGolden_commute
  unfold projGoldenR scaledProj faceNormalR
  rw [Matrix.mul_smul, Matrix.smul_mul, h1, h2]

theorem golden_quadratic_mat :
    faceNormalR * (faceNormalR * projGoldenR) =
      (6 : ℝ) • (faceNormalR * projGoldenR) - (4 : ℝ) • projGoldenR := by
  have hW : castZ faceNormalZ * castZ projGoldenZ = castZ goldenImageZ :=
    castZ_mul_eq _ _ _ projGolden_image
  unfold projGoldenR scaledProj faceNormalR
  rw [Matrix.mul_smul, hW, Matrix.mul_smul, smul_smul, smul_smul]
  ext i j
  simp only [Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul, Matrix.mul_apply, castZ_apply]
  have hq := projGolden_quadratic i j
  have hc : (∑ k : Fin 20, (faceNormalZ i k : ℝ) * (goldenImageZ k j : ℝ)) -
      6 * (goldenImageZ i j : ℝ) + 4 * (projGoldenZ i j : ℝ) = 0 := by
    exact_mod_cast hq
  linear_combination ((10 : ℝ))⁻¹ * hc

theorem golden_quadratic_vec (w : Fin 20 → ℝ) :
    faceNormalR.mulVec (faceNormalR.mulVec (projGoldenR.mulVec w)) =
      (6 : ℝ) • faceNormalR.mulVec (projGoldenR.mulVec w) - (4 : ℝ) • projGoldenR.mulVec w := by
  have h := congrArg (fun M : Matrix (Fin 20) (Fin 20) ℝ ↦ M.mulVec w) golden_quadratic_mat
  simp only [Matrix.sub_mulVec, Matrix.smul_mulVec, ← Matrix.mulVec_mulVec] at h
  exact h

theorem golden_pairing (w : Fin 20 → ℝ) :
    faceInner w (faceNormalR.mulVec (projGoldenR.mulVec w)) =
      faceInner (projGoldenR.mulVec w) (faceNormalR.mulVec (projGoldenR.mulVec w)) := by
  have hNP : faceNormalR * projGoldenR = projGoldenR * (faceNormalR * projGoldenR) := by
    rw [← Matrix.mul_assoc, ← golden_commute, Matrix.mul_assoc, golden_idem]
  have h1 : faceNormalR.mulVec (projGoldenR.mulVec w) =
      projGoldenR.mulVec (faceNormalR.mulVec (projGoldenR.mulVec w)) := by
    rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec, ← hNP]
  conv_lhs => rw [h1]
  rw [← faceInner_mulVec_symm projGoldenR golden_symm]

/-- On the golden sector `⟨u, N u⟩ ≤ (3 + √5) ‖u‖²`, from
`‖N u - 3u‖² = 5 ‖u‖²` and `‖√5 u - (N u - 3u)‖² ≥ 0`. -/
theorem golden_sector_bound (u : Fin 20 → ℝ)
    (hu : faceNormalR.mulVec (faceNormalR.mulVec u) =
      (6 : ℝ) • faceNormalR.mulVec u - (4 : ℝ) • u) :
    faceInner u (faceNormalR.mulVec u) ≤ (3 + Real.sqrt 5) * faceEnergy u := by
  have hs : Real.sqrt 5 * Real.sqrt 5 = 5 := Real.mul_self_sqrt (by norm_num)
  have hNN : faceEnergy (faceNormalR.mulVec u) =
      6 * faceInner u (faceNormalR.mulVec u) - 4 * faceEnergy u := by
    rw [← faceInner_self_eq_energy, faceInner_mulVec_symm faceNormalR faceNormal_symm, hu,
      faceInner_sub_right, faceInner_smul_right, faceInner_smul_right,
      faceInner_self_eq_energy]
  have hdiff : faceEnergy (faceNormalR.mulVec u - (3 : ℝ) • u) = 5 * faceEnergy u := by
    rw [faceEnergy_sub, faceInner_smul_right, faceEnergy_smul, hNN,
      faceInner_comm (faceNormalR.mulVec u) u]
    ring
  have hpos : 0 ≤ faceEnergy (Real.sqrt 5 • u - (faceNormalR.mulVec u - (3 : ℝ) • u)) :=
    faceEnergy_nonneg _
  rw [faceEnergy_sub, faceEnergy_smul, hdiff, faceInner_smul_left, faceInner_sub_right,
    faceInner_smul_right, faceInner_self_eq_energy] at hpos
  have h5 : Real.sqrt 5 ^ 2 = 5 := by rw [sq, hs]
  rw [h5] at hpos
  have hpos2 : 0 < 2 * Real.sqrt 5 := by positivity
  refine le_of_mul_le_mul_left ?_ hpos2
  have key : 2 * Real.sqrt 5 * ((3 + Real.sqrt 5) * faceEnergy u) =
      6 * Real.sqrt 5 * faceEnergy u + 10 * faceEnergy u := by
    linear_combination (2 * faceEnergy u) * hs
  rw [key]
  linarith

theorem proj_sum_R : projZeroR + projTwoR + projThreeR + projFiveR + projGoldenR = 1 := by
  ext i j
  simp only [projZeroR, projTwoR, projThreeR, projFiveR, projGoldenR, scaledProj,
    Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, castZ_apply, Matrix.one_apply]
  have h := proj_sum i j
  by_cases hij : i = j
  · rw [if_pos hij] at h
    rw [if_pos hij]
    have hc : 3 * (projZeroZ i j : ℝ) + 5 * (projTwoZ i j : ℝ) + 6 * (projThreeZ i j : ℝ) +
        2 * (projFiveZ i j : ℝ) + 6 * (projGoldenZ i j : ℝ) = 60 := by exact_mod_cast h
    push_cast
    linarith
  · rw [if_neg hij] at h
    rw [if_neg hij]
    have hc : 3 * (projZeroZ i j : ℝ) + 5 * (projTwoZ i j : ℝ) + 6 * (projThreeZ i j : ℝ) +
        2 * (projFiveZ i j : ℝ) + 6 * (projGoldenZ i j : ℝ) = 0 := by exact_mod_cast h
    push_cast
    linarith

/-- **D4, face side.**  `⟨w, N w⟩ ≤ (3 + √5) ‖w‖²` for every face load. -/
theorem faceNormal_quadratic_bound (w : Fin 20 → ℝ) :
    faceInner w (faceNormalR.mulVec w) ≤ (3 + Real.sqrt 5) * faceEnergy w := by
  have hsq : (2 : ℝ) ≤ Real.sqrt 5 :=
    (Real.le_sqrt (by norm_num) (by norm_num)).mpr (by norm_num)
  have hdecomp : projZeroR.mulVec w + projTwoR.mulVec w + projThreeR.mulVec w +
      projFiveR.mulVec w + projGoldenR.mulVec w = w := by
    have h := congrArg (fun M : Matrix (Fin 20) (Fin 20) ℝ ↦ M.mulVec w) proj_sum_R
    simp only [Matrix.add_mulVec, Matrix.one_mulVec] at h
    exact h
  have p0 := projector_pairing projZeroR (scaledProj_symm _ _ projZero_symm)
    (scaledProj_idem _ _ (by norm_num) projZero_idem) w
  have p2 := projector_pairing projTwoR (scaledProj_symm _ _ projTwo_symm)
    (scaledProj_idem _ _ (by norm_num) projTwo_idem) w
  have p3 := projector_pairing projThreeR (scaledProj_symm _ _ projThree_symm)
    (scaledProj_idem _ _ (by norm_num) projThree_idem) w
  have p5 := projector_pairing projFiveR (scaledProj_symm _ _ projFive_symm)
    (scaledProj_idem _ _ (by norm_num) projFive_idem) w
  have pQ := projector_pairing projGoldenR golden_symm golden_idem w
  have e0 := eigen_pairing projZeroR ((0 : ℤ) : ℝ) (scaledProj_eigen _ _ 0 projZero_eigen) w
  have e2 := eigen_pairing projTwoR ((2 : ℤ) : ℝ) (scaledProj_eigen _ _ 2 projTwo_eigen) w
  have e3 := eigen_pairing projThreeR ((3 : ℤ) : ℝ) (scaledProj_eigen _ _ 3 projThree_eigen) w
  have e5 := eigen_pairing projFiveR ((5 : ℤ) : ℝ) (scaledProj_eigen _ _ 5 projFive_eigen) w
  have eQ := golden_pairing w
  have bQ := golden_sector_bound (projGoldenR.mulVec w) (golden_quadratic_vec w)
  have hE : faceEnergy w = faceEnergy (projZeroR.mulVec w) + faceEnergy (projTwoR.mulVec w) +
      faceEnergy (projThreeR.mulVec w) + faceEnergy (projFiveR.mulVec w) +
      faceEnergy (projGoldenR.mulVec w) := by
    rw [← faceInner_self_eq_energy]
    conv_lhs => rw [← hdecomp]
    rw [faceInner_add_right, faceInner_add_right, faceInner_add_right, faceInner_add_right]
    rw [hdecomp, p0, p2, p3, p5, pQ]
  have hN : faceInner w (faceNormalR.mulVec w) =
      faceInner w (faceNormalR.mulVec (projZeroR.mulVec w)) +
      faceInner w (faceNormalR.mulVec (projTwoR.mulVec w)) +
      faceInner w (faceNormalR.mulVec (projThreeR.mulVec w)) +
      faceInner w (faceNormalR.mulVec (projFiveR.mulVec w)) +
      faceInner w (faceNormalR.mulVec (projGoldenR.mulVec w)) := by
    conv_lhs => rw [← hdecomp]
    rw [Matrix.mulVec_add, Matrix.mulVec_add, Matrix.mulVec_add, Matrix.mulVec_add,
      faceInner_add_right, faceInner_add_right, faceInner_add_right, faceInner_add_right,
      hdecomp]
  rw [hN, e0, e2, e3, e5, eQ, p0, p2, p3, p5, hE]
  push_cast
  have n0 := faceEnergy_nonneg (projZeroR.mulVec w)
  have n2 := faceEnergy_nonneg (projTwoR.mulVec w)
  have n3 := faceEnergy_nonneg (projThreeR.mulVec w)
  have n5 := faceEnergy_nonneg (projFiveR.mulVec w)
  nlinarith [mul_le_mul_of_nonneg_right hsq n0, mul_le_mul_of_nonneg_right hsq n2,
    mul_le_mul_of_nonneg_right hsq n3, mul_le_mul_of_nonneg_right hsq n5]

/-- **D4.**  The sharp Courant constant on the committed carrier:
`‖C v‖² ≤ (3 + √5) ‖v‖²` for every seam field. -/
theorem faceEnergy_curvature_le_golden (v : Fin 30 → ℝ) :
    faceEnergy (faceCurvature v) ≤ (3 + Real.sqrt 5) * realSeamEnergy v := by
  have hX0 : 0 ≤ faceEnergy (faceCurvature v) := faceEnergy_nonneg _
  have hpair : faceEnergy (faceCurvature v) =
      realSeamInner v (faceCodifferential (faceCurvature v)) := by
    rw [← faceInner_self_eq_energy, faceCurvature_codifferential_adjoint]
  have hcs : realSeamInner v (faceCodifferential (faceCurvature v)) ^ 2 ≤
      realSeamEnergy v * realSeamEnergy (faceCodifferential (faceCurvature v)) := by
    unfold realSeamInner realSeamEnergy
    exact Finset.sum_mul_sq_le_sq_mul_sq _ _ _
  have hN : realSeamEnergy (faceCodifferential (faceCurvature v)) =
      faceInner (faceCurvature v) (faceNormalR.mulVec (faceCurvature v)) := by
    rw [← realSeamInner_self_eq_energy, ← faceCurvature_codifferential_adjoint,
      faceNormal_mulVec, faceInner_comm]
  have hbound := faceNormal_quadratic_bound (faceCurvature v)
  rw [← hN] at hbound
  have hv0 := realSeamEnergy_nonneg v
  have hsq : 0 ≤ (3 + Real.sqrt 5) := by positivity
  rcases hX0.lt_or_eq with hpos | hzero
  · have h1 : faceEnergy (faceCurvature v) * faceEnergy (faceCurvature v) ≤
        ((3 + Real.sqrt 5) * realSeamEnergy v) * faceEnergy (faceCurvature v) := by
      calc faceEnergy (faceCurvature v) * faceEnergy (faceCurvature v)
          = realSeamInner v (faceCodifferential (faceCurvature v)) ^ 2 := by
              rw [hpair]
              ring
        _ ≤ realSeamEnergy v * realSeamEnergy (faceCodifferential (faceCurvature v)) := hcs
        _ ≤ realSeamEnergy v * ((3 + Real.sqrt 5) * faceEnergy (faceCurvature v)) :=
              mul_le_mul_of_nonneg_left hbound hv0
        _ = ((3 + Real.sqrt 5) * realSeamEnergy v) * faceEnergy (faceCurvature v) := by
              ring
    exact le_of_mul_le_mul_right h1 hpos
  · rw [← hzero]
    exact mul_nonneg hsq hv0

/-- The committed carrier satisfies the Courant hypothesis with the sharp
constant `Λ = 3 + √5`. -/
theorem committed_courant_sharp : CourantBound (3 + Real.sqrt 5) :=
  faceEnergy_curvature_le_golden

/-- **Sharp threshold on the committed carrier.**  Below
`h² (3 + √5) = 4` every zero-current solution is bounded by explicit
multiples of the initial staggered form; above it an unbounded zero-current
solution exists. -/
theorem courant_threshold_sharp (h : ℝ) (hh : h ≠ 0) :
    (h ^ 2 * (3 + Real.sqrt 5) < 4 →
      ∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
        AmpereEvolutionScaled h A φ (fun _ ↦ 0) → ∀ n,
          0 ≤ fieldEnergyScaled h A φ n ∧
          realSeamEnergy (electricFieldScaled h A φ n) ≤
            8 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * (3 + Real.sqrt 5)) ∧
          faceEnergy (magneticField A n) ≤
            16 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * (3 + Real.sqrt 5)))
    ∧ (4 < h ^ 2 * (3 + Real.sqrt 5) →
      ∃ A : ℕ → Fin 30 → ℝ, AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0) ∧
        ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h A (fun _ ↦ 0) n)) := by
  refine ⟨fun hc A φ hAmp n ↦ ?_, fun hc ↦ instability_above_golden h hh hc⟩
  have hΛ0 : 0 ≤ 3 + Real.sqrt 5 := by positivity
  exact ⟨fieldEnergyScaled_nonneg h hh _ committed_courant_sharp hc.le A φ n,
    stability_certificate h hh _ hΛ0 committed_courant_sharp hc A φ hAmp n⟩

/-! ## (E) The one citable composed receipt -/

/-- The single typed antecedent bundle: a positive step, a Courant
constant with the strict Courant condition, declared histories and
sources, the scaled Ampere evolution, the initial Gauss constraint, and the
scaled continuity equation.  The port load `rho` and the seam current `J`
are separate fields with no conversion map. -/
structure ScaledMaxwellBundle where
  h : ℝ
  Λ : ℝ
  h_pos : 0 < h
  Λ_nonneg : 0 ≤ Λ
  courant : CourantBound Λ
  courant_strict : h ^ 2 * Λ < 4
  A : ℕ → Fin 30 → ℝ
  phi : ℕ → Fin 12 → ℝ
  rho : ℕ → Fin 12 → ℝ
  J : ℕ → Fin 30 → ℝ
  ampere : AmpereEvolutionScaled h A phi J
  gauss_init : realBoundary (electricFieldScaled h A phi 0) = rho 0
  continuity : ∀ n : ℕ, rho (n + 1) - rho n + h • realBoundary (J n) = 0

/-- **The scaled Maxwell stability receipt.**  One typed conjunction from
the single antecedent bundle: Gauss propagation; the step-local continuity
equivalence; the scaled Faraday law; the exact balance of the scaled
staggered form; scaled gauge invariance of `E`, `B`, the evolution law and
the form; for every window the two action equivalences and the derived
stationarity of the bundle's history in `A` and in `φ`; positivity of the
staggered form at every step; and, for zero current, conservation, the
scaled wave law, and the explicit uniform bounds.  The step index is a
declared evolution parameter and the sources are declared; PR-66 is
relocated to the declared action and PR-15, PR-53, PR-54 stay open. -/
theorem scaledMaxwellStability_receipt (S : ScaledMaxwellBundle) :
    (∀ n, realBoundary (electricFieldScaled S.h S.A S.phi n) = S.rho n)
    ∧ (∀ (ρ' : ℕ → Fin 12 → ℝ) (n : ℕ),
        realBoundary (electricFieldScaled S.h S.A S.phi n) = ρ' n →
        (realBoundary (electricFieldScaled S.h S.A S.phi (n + 1)) = ρ' (n + 1) ↔
          ρ' (n + 1) - ρ' n + S.h • realBoundary (S.J n) = 0))
    ∧ (∀ n, magneticField S.A (n + 1) - magneticField S.A n =
        -(S.h • faceCurvature (electricFieldScaled S.h S.A S.phi n)))
    ∧ (∀ n, fieldEnergyScaled S.h S.A S.phi (n + 1) =
        fieldEnergyScaled S.h S.A S.phi n -
          (S.h / 2) * realSeamInner
            (electricFieldScaled S.h S.A S.phi n +
              electricFieldScaled S.h S.A S.phi (n + 1)) (S.J n))
    ∧ (∀ χ : ℕ → Fin 12 → ℝ,
        (∀ n, electricFieldScaled S.h (gaugeTransformA S.A χ)
            (gaugeTransformPhiScaled S.h S.phi χ) n =
          electricFieldScaled S.h S.A S.phi n)
        ∧ (∀ n, magneticField (gaugeTransformA S.A χ) n = magneticField S.A n)
        ∧ AmpereEvolutionScaled S.h (gaugeTransformA S.A χ)
            (gaugeTransformPhiScaled S.h S.phi χ) S.J
        ∧ (∀ n, fieldEnergyScaled S.h (gaugeTransformA S.A χ)
            (gaugeTransformPhiScaled S.h S.phi χ) n =
          fieldEnergyScaled S.h S.A S.phi n))
    ∧ (∀ N : ℕ,
        ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
          windowAction S.h N (S.A + a) S.phi S.rho S.J =
            windowAction S.h N S.A S.phi S.rho S.J + quadraticRemainder S.h N a 0)
          ↔ (∀ m, m < N → ampereResidual S.h S.A S.phi S.J m = 0))
        ∧ ((∀ f : ℕ → Fin 12 → ℝ,
          windowAction S.h N S.A (S.phi + f) S.rho S.J =
            windowAction S.h N S.A S.phi S.rho S.J + quadraticRemainder S.h N 0 f)
          ↔ (∀ n, n < N + 1 →
            realBoundary (electricFieldScaled S.h S.A S.phi n) = S.rho n))
        ∧ (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
          windowAction S.h N (S.A + a) S.phi S.rho S.J =
            windowAction S.h N S.A S.phi S.rho S.J + quadraticRemainder S.h N a 0)
        ∧ (∀ f : ℕ → Fin 12 → ℝ,
          windowAction S.h N S.A (S.phi + f) S.rho S.J =
            windowAction S.h N S.A S.phi S.rho S.J + quadraticRemainder S.h N 0 f))
    ∧ (∀ n, 0 ≤ fieldEnergyScaled S.h S.A S.phi n)
    ∧ (S.J = (fun _ ↦ 0) →
        (∀ n, fieldEnergyScaled S.h S.A S.phi n = fieldEnergyScaled S.h S.A S.phi 0)
        ∧ (∀ n, S.A (n + 2) - (2 : ℝ) • S.A (n + 1) + S.A n +
            (S.h ^ 2) • localMaxwellOperator (S.A (n + 1)) +
            S.h • realCoboundary (S.phi (n + 1) - S.phi n) = 0)
        ∧ (∀ n, realSeamEnergy (electricFieldScaled S.h S.A S.phi n) ≤
              8 * fieldEnergyScaled S.h S.A S.phi 0 / (4 - S.h ^ 2 * S.Λ)
            ∧ faceEnergy (magneticField S.A n) ≤
              16 * fieldEnergyScaled S.h S.A S.phi 0 / (4 - S.h ^ 2 * S.Λ))) := by
  have hh : S.h ≠ 0 := S.h_pos.ne'
  have hgauss := gauss_propagation_scaled S.h S.A S.phi S.J S.rho S.ampere
    S.gauss_init S.continuity
  refine ⟨hgauss, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro ρ' n hn
    exact gauss_step_iff_scaled S.h S.A S.phi S.J ρ' S.ampere n hn
  · intro n
    exact faraday_law_scaled S.h hh S.A S.phi n
  · exact energy_balance_scaled S.h hh S.A S.phi S.J S.ampere
  · intro χ
    exact ⟨electricFieldScaled_gauge_invariant S.h S.A S.phi χ,
      magneticField_gauge_invariant S.A χ,
      (ampereEvolutionScaled_gauge_invariant S.h S.A S.phi S.J χ).mpr S.ampere,
      fieldEnergyScaled_gauge_invariant S.h S.A S.phi χ⟩
  · intro N
    have hres : ∀ m, ampereResidual S.h S.A S.phi S.J m = 0 :=
      (ampereEvolutionScaled_iff_residual S.h S.A S.phi S.J).mp S.ampere
    exact ⟨action_stationary_A_iff_ampere S.h hh N S.A S.phi S.rho S.J,
      action_stationary_phi_iff_gauss S.h hh N S.A S.phi S.rho S.J,
      (action_stationary_A_iff_ampere S.h hh N S.A S.phi S.rho S.J).mpr
        (fun m _ ↦ hres m),
      (action_stationary_phi_iff_gauss S.h hh N S.A S.phi S.rho S.J).mpr
        (fun n _ ↦ hgauss n)⟩
  · intro n
    exact fieldEnergyScaled_nonneg S.h hh S.Λ S.courant S.courant_strict.le S.A S.phi n
  · intro hJ
    have hAmp0 : AmpereEvolutionScaled S.h S.A S.phi (fun _ ↦ 0) := by
      rw [← hJ]
      exact S.ampere
    exact ⟨energy_conserved_scaled S.h hh S.A S.phi hAmp0,
      wave_law_scaled S.h hh S.A S.phi hAmp0,
      fun n ↦ stability_certificate S.h hh S.Λ S.Λ_nonneg S.courant S.courant_strict
        S.A S.phi hAmp0 n⟩

/-! ## An explicit nonstatic inhabitant at `h = 1/2` -/

/-- Zero initial data kicked by the committed face-zero boundary cycle,
evolved by the scaled recursion in the gauge `φ = 0`. -/
def demoScaledA (h : ℝ) : ℕ → Fin 30 → ℝ
  | 0 => 0
  | 1 => demoInitial
  | n + 2 => (2 : ℝ) • demoScaledA h (n + 1) - demoScaledA h n -
      (h ^ 2) • localMaxwellOperator (demoScaledA h (n + 1))

theorem demoScaledA_zero (h : ℝ) : demoScaledA h 0 = 0 := by
  simp only [demoScaledA]

theorem demoScaledA_one (h : ℝ) : demoScaledA h 1 = demoInitial := by
  simp only [demoScaledA]

theorem demoScaled_ampere (h : ℝ) (hh : h ≠ 0) :
    AmpereEvolutionScaled h (demoScaledA h) (fun _ ↦ 0) (fun _ ↦ 0) := by
  intro n
  rw [electricFieldScaled_temporal_gauge, electricFieldScaled_temporal_gauge]
  have hlm : faceCodifferential (magneticField (demoScaledA h) (n + 1)) =
      localMaxwellOperator (demoScaledA h (n + 1)) := rfl
  rw [hlm]
  have hidx : demoScaledA h (n + 1 + 1) = (2 : ℝ) • demoScaledA h (n + 1) -
      demoScaledA h n - (h ^ 2) • localMaxwellOperator (demoScaledA h (n + 1)) := by
    simp only [demoScaledA]
  rw [hidx]
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  funext e
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul, sub_zero]
  linear_combination (h * localMaxwellOperator (demoScaledA h (n + 1)) e) * hinv

/-- The explicit inhabitant of the antecedent bundle at `h = 1/2`, `Λ = 6`:
`h² Λ = 3/2 < 4`. -/
def demoScaledBundle : ScaledMaxwellBundle where
  h := 1 / 2
  Λ := 6
  h_pos := by norm_num
  Λ_nonneg := by norm_num
  courant := committed_courant
  courant_strict := by norm_num
  A := demoScaledA (1 / 2)
  phi := fun _ ↦ 0
  rho := fun _ ↦ 0
  J := fun _ ↦ 0
  ampere := demoScaled_ampere (1 / 2) (by norm_num)
  gauss_init := by
    rw [electricFieldScaled_temporal_gauge, demoScaledA_one, demoScaledA_zero, sub_zero,
      map_neg, map_smul, demoInitial_boundary, smul_zero, neg_zero]
  continuity := by
    intro n
    simp

theorem demoScaledBundle_nonvacuous :
    demoScaledBundle.A 1 ≠ demoScaledBundle.A 0 ∧
      magneticField demoScaledBundle.A 0 = 0 ∧
      magneticField demoScaledBundle.A 1 0 = 3 := by
  have h1 : demoScaledBundle.A 1 = demoA 1 := by
    rw [show demoScaledBundle.A = demoScaledA (1 / 2) from rfl, demoScaledA_one,
      show demoA 1 = demoInitial from by simp only [demoA]]
  have h0 : demoScaledBundle.A 0 = demoA 0 := by
    rw [show demoScaledBundle.A = demoScaledA (1 / 2) from rfl, demoScaledA_zero,
      show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]]
  refine ⟨?_, ?_, ?_⟩
  · rw [h1, h0]
    exact demo_nonstatic
  · show faceCurvature (demoScaledBundle.A 0) = 0
    rw [h0]
    exact demo_magnetic_start
  · show faceCurvature (demoScaledBundle.A 1) 0 = 3
    rw [h1]
    exact demo_magnetic_step_one

end

end OPH.ScaledMaxwellStability

/- Axiom audit: the committed receipts, exact real linear algebra, and
kernel `decide` on the committed integer tables only.  Expected axioms per
line: at most `propext`, `Classical.choice`, `Quot.sound`.  No native
decision procedure is used. -/

#print axioms OPH.ScaledMaxwellStability.electricFieldScaled_one
#print axioms OPH.ScaledMaxwellStability.ampereEvolutionScaled_one_iff
#print axioms OPH.ScaledMaxwellStability.faraday_law_scaled
#print axioms OPH.ScaledMaxwellStability.gauss_step_iff_scaled
#print axioms OPH.ScaledMaxwellStability.gauss_propagation_scaled
#print axioms OPH.ScaledMaxwellStability.electricFieldScaled_gauge_invariant
#print axioms OPH.ScaledMaxwellStability.ampereEvolutionScaled_gauge_invariant
#print axioms OPH.ScaledMaxwellStability.wave_law_scaled
#print axioms OPH.ScaledMaxwellStability.energy_balance_scaled
#print axioms OPH.ScaledMaxwellStability.energy_conserved_scaled
#print axioms OPH.ScaledMaxwellStability.windowAction_expansion
#print axioms OPH.ScaledMaxwellStability.firstVariation_interior
#print axioms OPH.ScaledMaxwellStability.action_stationary_A_iff_ampere
#print axioms OPH.ScaledMaxwellStability.action_stationary_phi_iff_gauss
#print axioms OPH.ScaledMaxwellStability.fieldEnergyScaled_eq
#print axioms OPH.ScaledMaxwellStability.fieldEnergyScaled_lower_bound
#print axioms OPH.ScaledMaxwellStability.fieldEnergyScaled_nonneg
#print axioms OPH.ScaledMaxwellStability.fieldEnergyScaled_lower_bound_temporal_gauge
#print axioms OPH.ScaledMaxwellStability.stability_certificate
#print axioms OPH.ScaledMaxwellStability.unstable_mode
#print axioms OPH.ScaledMaxwellStability.faceEnergy_curvature_le_six
#print axioms OPH.ScaledMaxwellStability.fiveMode_eigen
#print axioms OPH.ScaledMaxwellStability.goldenMode_eigen
#print axioms OPH.ScaledMaxwellStability.instability_above_five
#print axioms OPH.ScaledMaxwellStability.instability_above_golden
#print axioms OPH.ScaledMaxwellStability.unit_step_instability
#print axioms OPH.ScaledMaxwellStability.committed_carrier_stability
#print axioms OPH.ScaledMaxwellStability.faceNormal_eq_incidence
#print axioms OPH.ScaledMaxwellStability.proj_sum
#print axioms OPH.ScaledMaxwellStability.faceNormal_mulVec
#print axioms OPH.ScaledMaxwellStability.golden_sector_bound
#print axioms OPH.ScaledMaxwellStability.faceNormal_quadratic_bound
#print axioms OPH.ScaledMaxwellStability.faceEnergy_curvature_le_golden
#print axioms OPH.ScaledMaxwellStability.committed_courant_sharp
#print axioms OPH.ScaledMaxwellStability.courant_threshold_sharp
#print axioms OPH.ScaledMaxwellStability.scaledMaxwellStability_receipt
#print axioms OPH.ScaledMaxwellStability.demoScaled_ampere
#print axioms OPH.ScaledMaxwellStability.demoScaledBundle_nonvacuous
