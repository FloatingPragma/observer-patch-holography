import ScaledMaxwellStability
import QFT.SourceRecurrenceClock
import Geometry.CommonWorldMaxwellClockJoin

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.CarrierModeOscillators

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.PhysicalCalibrationImport
open OPH.QFT.SourceRecurrenceClock
open OPH.CommonWorldMaxwellClockJoin

/-!
# Periodic mode oscillators of the committed carrier and their step-index
commensurability with the source recurrence clock

STATUS.  Exact finite real linear algebra, elementary trigonometry, and
kernel `decide` checks on the committed integer tables of
`ScaledMaxwellStability`.  Every history below is a solution of the
committed zero-current scaled Ampere evolution `AmpereEvolutionScaled h A 0 0`
in the temporal gauge `φ = 0`.  The step `h` is a declared Courant number.
No register row is discharged.

WHAT IS PROVED.
1. Mode oscillators.  For a real step `h`, an eigenvalue `lam` of the local
   operator `CᵀC` with `0 ≤ h² lam ≤ 4`, the declared mode angle
   `modeAngle h lam = arccos (1 - h² lam / 2)` obeys
   `2 - 2 cos θ = h² lam` (`two_sub_two_cos_modeAngle`).  For an eigenvector
   `CᵀC v = lam • v` the histories `cosHistory θ v n = cos (n θ) • v` and
   `sinHistory θ v n = sin (n θ) • v` solve the zero-current scaled
   evolution in the temporal gauge (`cosHistory_ampere`, `sinHistory_ampere`,
   through the scalar recurrence `scalarHistory_ampere` and the trigonometric
   identities `cos_recurrence`, `sin_recurrence`); the staggered form is
   constant along them (`cosHistory_energy_conserved`,
   `sinHistory_energy_conserved`, by `energy_conserved_scaled`); the history
   and its electric field are bounded by explicit multiples of `‖v‖²`
   (`cosHistory_bounded`, `cosHistory_electric_bounded`, and the sine
   versions).  The committed eigenvectors `fiveMode` (`lam = 5`) and
   `goldenMode` (`lam = 3 + √5`) are instantiated (`fiveMode_oscillator`,
   `goldenMode_oscillator`).
2. Integer eigenvectors for `lam = 2` and `lam = 3`.  `twoMode` and
   `threeMode` are explicit integer seam vectors, the codifferentials of the
   zeroth columns of the committed projectors `projTwoZ` and `projThreeZ`,
   with `CᵀC twoMode = 2 • twoMode` and `CᵀC threeMode = 3 • threeMode`
   kernel-checked (`twoMode_eigen`, `threeMode_eigen`), both nonzero.  The
   transport lemma `codifferential_eigen` records that `Cᵀ` carries every
   face eigenvector of `C Cᵀ` to a seam eigenvector of `CᵀC` with the same
   eigenvalue.
3. Exact integer periods.  `IsPeriod A p` declares the return
   `A (n + p) = A n` of the field configuration for every `n`.  At
   `h² = 1/2` the `lam = 2` cosine and sine histories have period `6`
   (`twoMode_period_six`), at `h² = 2/3` the `lam = 3` histories have period
   `4` (`threeMode_period_four`), at `h² = 1/3` the `lam = 3` histories have
   period `6` (`threeMode_period_six`); each of these steps lies strictly
   inside the sharp window `h² (3 + √5) < 4` (`window_half`,
   `window_two_thirds`, `window_third`, via `sqrt5_lt_three`).
4. Golden structure.  `3 + √5 = 2 φ²` and `3 - √5 = 2 / φ²` for
   `φ = Real.goldenRatio` (`three_add_sqrt5_eq_two_goldenRatio_sq`,
   `three_sub_sqrt5_eq_two_div_goldenRatio_sq`), so the sharp window reads
   `h² < 2 / φ²`, equivalently `h < √2 / φ` for `h > 0`
   (`window_iff_goldenRatio`, `window_iff_lt_sqrt_two_div_goldenRatio`).
   The traces of the five committed integer projectors and of the golden
   image are kernel-checked (`projZero_trace`, `projTwo_trace`,
   `projThree_trace`, `projFive_trace`, `projGolden_trace`,
   `goldenImage_trace`), giving the real traces `1, 5, 4, 4, 6` of the
   rational projectors and trace `18` of `N` on the golden sector
   (`projector_traces`, `golden_sector_normal_trace`).
5. Tick-free frequency ratios.  For any `ClockCalibration` and nonzero
   denominator angle `θ₂`, the ratio of two laboratory frequencies
   `labFrequency θ₁ / labFrequency θ₂` equals `θ₁ / θ₂` and is the same for
   every calibration (`labFrequency_ratio_of_ne_zero`,
   `labFrequency_ratio_tick_free_of_ne_zero`).  The compatibility lemmas
   `labFrequency_ratio` and `labFrequency_ratio_tick_free` retain the
   corresponding totalized real-division identities for all angles; at
   `θ₂ = 0` those identities are only Lean's zero-denominator algebra and are
   not a physical frequency ratio.  The mode angle
   obeys the exact two-sided bound
   `h √lam ≤ modeAngle h lam ≤ h √lam / √(1 - h² lam / 4)`
   (`modeAngle_ge`, `modeAngle_le`), hence the frequency ratio of two modes
   is squeezed between explicit multiples of `√(lam₁ / lam₂)`
   (`modeAngle_ratio_bounds`), and `modeAngle h lam / h` tends to `√lam` as
   `h` tends to `0` from the right (`modeAngle_div_tendsto`).  The golden
   over slowest ratio `√((3 + √5) / (3 - √5))` equals `φ²`
   (`sqrt_golden_ratio_eq_goldenRatio_sq`).
6. Commensurability in one step index.  For any calibration the ratio of the
   laboratory duration of `p` steps to the laboratory mean return time of
   the source recurrence clock is `7155 p / 61511`
   (`period_return_ratio`), the same for every calibration
   (`period_return_ratio_tick_free`), with the values `28620 / 61511` at
   `p = 4` and `42930 / 61511` at `p = 6` (`period_return_ratio_values`).
   In a `MaxwellClockJoinedArchitecture` the join at index `n + p` differs
   from the join at index `n` for every `p > 0` (`join_ne_of_period`), so a
   field period is a return of the configuration in one shared index and
   never a return of the joined record.

ROWS TOUCHED (none discharged).  Source clock and duration row: `h` is a
declared step with no unit and the join reads both islands at one step
index, with no physical duration attached to that index
(`stepDuration_not_forced`, `returnDuration_not_forced`).
Physical spacetime attachment row and light-signal row: the identification
of a mode with a physical oscillation is open.  Coupled-action row: the
kinetic term `(h/2) ‖E‖²` is declared in `ScaledMaxwellStability`.
Laboratory clock and energy calibration import: `ClockCalibration.tau` is
the single declared tick, and only tick-free ratios are stated here.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized source
histories select no velocity curvature or Legendre map, so every Lagrangian
shape is a declared enrichment; the mode oscillators inherit the declared
kinetic term.  Abstract rate non-identifiability: its repair-layer bridge
fails per `RateBridgeObstruction`, which forbids nothing about a
source-hosted process.  The repair-layer theorems `periodic_iterate_is_fixed`
and `acceptedStep_no_cycle` of `SourceRecurrenceClock` state that accepted
repair hosts no periodic orbit; the periodic histories here live in the
field sector on the carrier and are solutions of a declared evolution, so
the two statements have disjoint scope and no contradiction.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`.  Forward differences: the electric seam
field on the half step is `E n = -(h⁻¹ • (A (n+1) - A n)) - d (φ n)` and the
magnetic face field is `B n = C (A n)`; the potential `A` lives on oriented
seams with the committed smaller-to-larger orientation and `φ` on ports.
The temporal gauge `φ = 0` and zero seam current are used throughout.
Angles are in radians per step; `labFrequency θ = θ / (2 π τ)` hertz.

FALSIFIER.  A kernel evaluation of `CᵀC twoMode` or `CᵀC threeMode` that is
not `2 • twoMode` or `3 • threeMode`, a projector trace other than the ones
stated, a cosine history at `h² = 1/2`, `lam = 2` with `A 6 ≠ A 0`, or a
mode angle outside the stated two-sided bound would make the module wrong.

Axiom audit.  The `#print axioms` lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## 1. Mode angle and the scalar recurrence -/

/-- Declared mode angle: `arccos (1 - h² lam / 2)`, the angular advance per
step of a mode with eigenvalue `lam` at declared step `h`. -/
def modeAngle (h lam : ℝ) : ℝ := Real.arccos (1 - h ^ 2 * lam / 2)

theorem cos_modeAngle (h lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4) :
    Real.cos (modeAngle h lam) = 1 - h ^ 2 * lam / 2 := by
  unfold modeAngle
  exact Real.cos_arccos (by linarith) (by linarith)

theorem two_sub_two_cos_modeAngle (h lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam)
    (h4 : h ^ 2 * lam ≤ 4) :
    2 - 2 * Real.cos (modeAngle h lam) = h ^ 2 * lam := by
  rw [cos_modeAngle h lam h0 h4]
  ring

/-- A scalar-profile history `A n = c n • v` on a fixed seam vector. -/
def scalarHistory (c : ℕ → ℝ) (v : Fin 30 → ℝ) : ℕ → Fin 30 → ℝ :=
  fun n ↦ c n • v

theorem scalarHistory_electricField (h : ℝ) (c : ℕ → ℝ) (v : Fin 30 → ℝ)
    (n : ℕ) :
    electricFieldScaled h (scalarHistory c v) (fun _ ↦ 0) n =
      (-(h⁻¹ * (c (n + 1) - c n))) • v := by
  unfold electricFieldScaled scalarHistory
  rw [map_zero, sub_zero]
  funext e
  simp only [Pi.neg_apply, Pi.smul_apply, Pi.sub_apply, smul_eq_mul]
  ring

/-- A scalar profile obeying the three-term recurrence
`c (n+2) - 2 c (n+1) + c n + h² lam c (n+1) = 0` on an eigenvector of
`CᵀC` with eigenvalue `lam` solves the zero-current scaled evolution in the
temporal gauge. -/
theorem scalarHistory_ampere (h : ℝ) (hh : h ≠ 0) (c : ℕ → ℝ)
    (v : Fin 30 → ℝ) (lam : ℝ) (hv : localMaxwellOperator v = lam • v)
    (hc : ∀ n, c (n + 2) - 2 * c (n + 1) + c n + h ^ 2 * lam * c (n + 1) = 0) :
    AmpereEvolutionScaled h (scalarHistory c v) (fun _ ↦ 0) (fun _ ↦ 0) := by
  intro n
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  have hlm : faceCodifferential (magneticField (scalarHistory c v) (n + 1)) =
      c (n + 1) • (lam • v) := by
    show localMaxwellOperator (c (n + 1) • v) = _
    rw [map_smul, hv]
  rw [hlm, scalarHistory_electricField, scalarHistory_electricField]
  funext e
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul, sub_zero]
  have hr := hc n
  linear_combination (-(h⁻¹ * v e)) * hr +
    (h * lam * c (n + 1) * v e) * hinv

theorem scalarHistory_seamEnergy (c : ℕ → ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (scalarHistory c v n) = c n ^ 2 * realSeamEnergy v := by
  unfold scalarHistory
  exact seamEnergy_smul _ _

theorem scalarHistory_electric_seamEnergy (h : ℝ) (c : ℕ → ℝ)
    (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (electricFieldScaled h (scalarHistory c v) (fun _ ↦ 0) n) =
      (h⁻¹ * (c (n + 1) - c n)) ^ 2 * realSeamEnergy v := by
  rw [scalarHistory_electricField, seamEnergy_smul]
  ring

/-! ## 1b. Cosine and sine histories -/

/-- Trigonometric three-term recurrence:
`cos ((n+2) θ) - 2 cos ((n+1) θ) + cos (n θ) + (2 - 2 cos θ) cos ((n+1) θ) = 0`. -/
theorem cos_recurrence (θ : ℝ) (n : ℕ) :
    Real.cos (((n + 2 : ℕ) : ℝ) * θ) - 2 * Real.cos (((n + 1 : ℕ) : ℝ) * θ) +
      Real.cos ((n : ℝ) * θ) +
      (2 - 2 * Real.cos θ) * Real.cos (((n + 1 : ℕ) : ℝ) * θ) = 0 := by
  have h2 : ((n + 2 : ℕ) : ℝ) * θ = ((n + 1 : ℕ) : ℝ) * θ + θ := by push_cast; ring
  have h0 : (n : ℝ) * θ = ((n + 1 : ℕ) : ℝ) * θ - θ := by push_cast; ring
  rw [h2, h0, Real.cos_add, Real.cos_sub]
  ring

theorem sin_recurrence (θ : ℝ) (n : ℕ) :
    Real.sin (((n + 2 : ℕ) : ℝ) * θ) - 2 * Real.sin (((n + 1 : ℕ) : ℝ) * θ) +
      Real.sin ((n : ℝ) * θ) +
      (2 - 2 * Real.cos θ) * Real.sin (((n + 1 : ℕ) : ℝ) * θ) = 0 := by
  have h2 : ((n + 2 : ℕ) : ℝ) * θ = ((n + 1 : ℕ) : ℝ) * θ + θ := by push_cast; ring
  have h0 : (n : ℝ) * θ = ((n + 1 : ℕ) : ℝ) * θ - θ := by push_cast; ring
  rw [h2, h0, Real.sin_add, Real.sin_sub]
  ring

/-- The cosine history `A n = cos (n θ) • v`. -/
def cosHistory (θ : ℝ) (v : Fin 30 → ℝ) : ℕ → Fin 30 → ℝ :=
  scalarHistory (fun n ↦ Real.cos ((n : ℝ) * θ)) v

/-- The sine history `A n = sin (n θ) • v`. -/
def sinHistory (θ : ℝ) (v : Fin 30 → ℝ) : ℕ → Fin 30 → ℝ :=
  scalarHistory (fun n ↦ Real.sin ((n : ℝ) * θ)) v

theorem cosHistory_ampere (h : ℝ) (hh : h ≠ 0) (lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam)
    (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) :
    AmpereEvolutionScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0)
      (fun _ ↦ 0) := by
  apply scalarHistory_ampere h hh _ v lam hv
  intro n
  have hr := cos_recurrence (modeAngle h lam) n
  rw [two_sub_two_cos_modeAngle h lam h0 h4] at hr
  exact hr

theorem sinHistory_ampere (h : ℝ) (hh : h ≠ 0) (lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam)
    (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) :
    AmpereEvolutionScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0)
      (fun _ ↦ 0) := by
  apply scalarHistory_ampere h hh _ v lam hv
  intro n
  have hr := sin_recurrence (modeAngle h lam) n
  rw [two_sub_two_cos_modeAngle h lam h0 h4] at hr
  exact hr

/-- The staggered form is constant along the cosine history
(`energy_conserved_scaled`). -/
theorem cosHistory_energy_conserved (h : ℝ) (hh : h ≠ 0) (lam : ℝ)
    (h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) (n : ℕ) :
    fieldEnergyScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0) n =
      fieldEnergyScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0) 0 :=
  energy_conserved_scaled h hh _ _ (cosHistory_ampere h hh lam h0 h4 v hv) n

theorem sinHistory_energy_conserved (h : ℝ) (hh : h ≠ 0) (lam : ℝ)
    (h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) (n : ℕ) :
    fieldEnergyScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0) n =
      fieldEnergyScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0) 0 :=
  energy_conserved_scaled h hh _ _ (sinHistory_ampere h hh lam h0 h4 v hv) n

theorem sq_le_one_of_abs_le_one {x : ℝ} (hx : |x| ≤ 1) : x ^ 2 ≤ 1 := by
  have h := abs_le.mp hx
  nlinarith [h.1, h.2]

/-- The cosine history is bounded by the seam energy of the mode. -/
theorem cosHistory_bounded (θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (cosHistory θ v n) ≤ realSeamEnergy v := by
  unfold cosHistory
  rw [scalarHistory_seamEnergy]
  have h1 : Real.cos ((n : ℝ) * θ) ^ 2 ≤ 1 :=
    sq_le_one_of_abs_le_one (Real.abs_cos_le_one _)
  have h2 := realSeamEnergy_nonneg v
  nlinarith

theorem sinHistory_bounded (θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (sinHistory θ v n) ≤ realSeamEnergy v := by
  unfold sinHistory
  rw [scalarHistory_seamEnergy]
  have h1 : Real.sin ((n : ℝ) * θ) ^ 2 ≤ 1 :=
    sq_le_one_of_abs_le_one (Real.abs_sin_le_one _)
  have h2 := realSeamEnergy_nonneg v
  nlinarith

/-- The electric seam field of the cosine history is bounded by
`4 h⁻² ‖v‖²`. -/
theorem cosHistory_electric_bounded (h θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (electricFieldScaled h (cosHistory θ v) (fun _ ↦ 0) n) ≤
      4 * (h⁻¹) ^ 2 * realSeamEnergy v := by
  unfold cosHistory
  rw [scalarHistory_electric_seamEnergy]
  have ha := Real.abs_cos_le_one (((n + 1 : ℕ) : ℝ) * θ)
  have hb := Real.abs_cos_le_one ((n : ℝ) * θ)
  have hd : (Real.cos (((n + 1 : ℕ) : ℝ) * θ) - Real.cos ((n : ℝ) * θ)) ^ 2 ≤ 4 := by
    have h1 := abs_le.mp ha
    have h2 := abs_le.mp hb
    nlinarith [h1.1, h1.2, h2.1, h2.2]
  have h2 := realSeamEnergy_nonneg v
  have hsq : (h⁻¹ * (Real.cos (((n + 1 : ℕ) : ℝ) * θ) - Real.cos ((n : ℝ) * θ))) ^ 2 =
      (h⁻¹) ^ 2 * (Real.cos (((n + 1 : ℕ) : ℝ) * θ) - Real.cos ((n : ℝ) * θ)) ^ 2 := by
    ring
  rw [hsq]
  have h3 : 0 ≤ (h⁻¹) ^ 2 := sq_nonneg _
  nlinarith [mul_nonneg h3 h2]

theorem sinHistory_electric_bounded (h θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    realSeamEnergy (electricFieldScaled h (sinHistory θ v) (fun _ ↦ 0) n) ≤
      4 * (h⁻¹) ^ 2 * realSeamEnergy v := by
  unfold sinHistory
  rw [scalarHistory_electric_seamEnergy]
  have ha := Real.abs_sin_le_one (((n + 1 : ℕ) : ℝ) * θ)
  have hb := Real.abs_sin_le_one ((n : ℝ) * θ)
  have hd : (Real.sin (((n + 1 : ℕ) : ℝ) * θ) - Real.sin ((n : ℝ) * θ)) ^ 2 ≤ 4 := by
    have h1 := abs_le.mp ha
    have h2 := abs_le.mp hb
    nlinarith [h1.1, h1.2, h2.1, h2.2]
  have h2 := realSeamEnergy_nonneg v
  have hsq : (h⁻¹ * (Real.sin (((n + 1 : ℕ) : ℝ) * θ) - Real.sin ((n : ℝ) * θ))) ^ 2 =
      (h⁻¹) ^ 2 * (Real.sin (((n + 1 : ℕ) : ℝ) * θ) - Real.sin ((n : ℝ) * θ)) ^ 2 := by
    ring
  rw [hsq]
  have h3 : 0 ≤ (h⁻¹) ^ 2 := sq_nonneg _
  nlinarith [mul_nonneg h3 h2]

/-- The mode oscillator packet: solution, conservation, boundedness, for both
phases. -/
def ModeOscillator (h : ℝ) (lam : ℝ) (v : Fin 30 → ℝ) : Prop :=
  AmpereEvolutionScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0) (fun _ ↦ 0) ∧
  AmpereEvolutionScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0) (fun _ ↦ 0) ∧
  (∀ n, fieldEnergyScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0) n =
    fieldEnergyScaled h (cosHistory (modeAngle h lam) v) (fun _ ↦ 0) 0) ∧
  (∀ n, fieldEnergyScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0) n =
    fieldEnergyScaled h (sinHistory (modeAngle h lam) v) (fun _ ↦ 0) 0) ∧
  (∀ n, realSeamEnergy (cosHistory (modeAngle h lam) v n) ≤ realSeamEnergy v) ∧
  (∀ n, realSeamEnergy (sinHistory (modeAngle h lam) v n) ≤ realSeamEnergy v) ∧
  (∀ n, realSeamEnergy (electricFieldScaled h (cosHistory (modeAngle h lam) v)
    (fun _ ↦ 0) n) ≤ 4 * (h⁻¹) ^ 2 * realSeamEnergy v) ∧
  (∀ n, realSeamEnergy (electricFieldScaled h (sinHistory (modeAngle h lam) v)
    (fun _ ↦ 0) n) ≤ 4 * (h⁻¹) ^ 2 * realSeamEnergy v)

theorem modeOscillator (h : ℝ) (hh : h ≠ 0) (lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam)
    (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) : ModeOscillator h lam v :=
  ⟨cosHistory_ampere h hh lam h0 h4 v hv, sinHistory_ampere h hh lam h0 h4 v hv,
    cosHistory_energy_conserved h hh lam h0 h4 v hv,
    sinHistory_energy_conserved h hh lam h0 h4 v hv,
    cosHistory_bounded _ v, sinHistory_bounded _ v,
    cosHistory_electric_bounded h _ v, sinHistory_electric_bounded h _ v⟩

/-- The committed `lam = 5` eigenvector hosts an oscillator for
`h² 5 ≤ 4`. -/
theorem fiveMode_oscillator (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * 5 ≤ 4) :
    ModeOscillator h 5 fiveMode :=
  modeOscillator h hh 5 (by positivity) h4 fiveMode fiveMode_eigen

/-- The committed golden eigenvector hosts an oscillator for
`h² (3 + √5) ≤ 4`. -/
theorem goldenMode_oscillator (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) ≤ 4) :
    ModeOscillator h (3 + Real.sqrt 5) goldenMode :=
  modeOscillator h hh _ (by positivity) h4 goldenMode goldenMode_eigen

/-! ## 2. Integer eigenvectors of `CᵀC` for `lam = 2` and `lam = 3` -/

/-- `Cᵀ` carries a face eigenvector of `C Cᵀ` to a seam eigenvector of
`CᵀC` with the same eigenvalue. -/
theorem codifferential_eigen (w : Fin 20 → ℝ) (lam : ℝ)
    (hw : faceNormalR.mulVec w = lam • w) :
    localMaxwellOperator (faceCodifferential w) = lam • faceCodifferential w := by
  show faceCodifferential (faceCurvature (faceCodifferential w)) = _
  rw [faceNormal_mulVec, hw, map_smul]

/-- Half the codifferential of the zeroth column of `projTwoZ`. -/
def twoModeZ : Fin 30 → ℤ :=
  ![1, -1, 1, -1, 0, 1, -1, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1,
    -1, -1, 1, 1, -1, 1]

/-- The zeroth column of `projTwoZ`. -/
def twoFaceZ : Fin 20 → ℤ :=
  ![3, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 3]

set_option maxRecDepth 16384 in
theorem twoMode_curvature_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * twoModeZ e) = twoFaceZ f := by
  decide

set_option maxRecDepth 16384 in
theorem twoMode_codifferential_Z :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * twoFaceZ f) = 2 * twoModeZ e := by
  decide

/-- The real cast of the `lam = 2` integer mode. -/
def twoMode : Fin 30 → ℝ := fun e ↦ (twoModeZ e : ℝ)

theorem twoMode_eigen : localMaxwellOperator twoMode = (2 : ℝ) • twoMode := by
  have h := localMaxwellOperator_cast twoModeZ (fun e ↦ 2 * twoModeZ e) twoFaceZ
    twoMode_curvature_Z twoMode_codifferential_Z
  rw [show twoMode = (fun e ↦ (twoModeZ e : ℝ)) from rfl, h]
  funext e
  simp only [Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

theorem twoMode_ne_zero : twoMode ≠ 0 := by
  intro h
  have h0 := congrFun h 0
  have hz : twoModeZ 0 = 1 := by decide
  have hc : ((twoModeZ 0 : ℤ) : ℝ) = 0 := h0
  rw [hz] at hc
  norm_num at hc

/-- The codifferential of the zeroth column of `projThreeZ`. -/
def threeModeZ : Fin 30 → ℤ :=
  ![2, -2, 1, -1, 0, 2, -1, 1, 0, 1, -1, 0, 2, -2, 0, -2, 2, 0, 2, -2, 0, 1, -1,
    -1, 1, 1, -1, -2, 2, -2]

/-- Three times the zeroth column of `projThreeZ`. -/
def threeFaceZ : Fin 20 → ℤ :=
  ![6, 0, 0, -3, -3, 0, -3, -3, -3, -3, 3, 3, 3, 3, 3, 3, 0, 0, 0, -6]

set_option maxRecDepth 16384 in
theorem threeMode_curvature_Z :
    ∀ f : Fin 20, (∑ e : Fin 30, faceIncidenceZ f e * threeModeZ e) = threeFaceZ f := by
  decide

set_option maxRecDepth 16384 in
theorem threeMode_codifferential_Z :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e * threeFaceZ f) = 3 * threeModeZ e := by
  decide

/-- The real cast of the `lam = 3` integer mode. -/
def threeMode : Fin 30 → ℝ := fun e ↦ (threeModeZ e : ℝ)

theorem threeMode_eigen : localMaxwellOperator threeMode = (3 : ℝ) • threeMode := by
  have h := localMaxwellOperator_cast threeModeZ (fun e ↦ 3 * threeModeZ e) threeFaceZ
    threeMode_curvature_Z threeMode_codifferential_Z
  rw [show threeMode = (fun e ↦ (threeModeZ e : ℝ)) from rfl, h]
  funext e
  simp only [Pi.smul_apply, smul_eq_mul]
  push_cast
  ring

theorem threeMode_ne_zero : threeMode ≠ 0 := by
  intro h
  have h0 := congrFun h 0
  have hz : threeModeZ 0 = 2 := by decide
  have hc : ((threeModeZ 0 : ℤ) : ℝ) = 0 := h0
  rw [hz] at hc
  norm_num at hc

theorem twoMode_oscillator (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * 2 ≤ 4) :
    ModeOscillator h 2 twoMode :=
  modeOscillator h hh 2 (by positivity) h4 twoMode twoMode_eigen

theorem threeMode_oscillator (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * 3 ≤ 4) :
    ModeOscillator h 3 threeMode :=
  modeOscillator h hh 3 (by positivity) h4 threeMode threeMode_eigen

/-! ## 3. Exact integer periods -/

/-- Declared notion of period: the field configuration returns after `p`
steps, `A (n + p) = A n` for every `n`.  A period is a count of steps in
the shared index, the quantity the source recurrence clock also counts. -/
def IsPeriod (A : ℕ → Fin 30 → ℝ) (p : ℕ) : Prop := ∀ n, A (n + p) = A n

theorem cosHistory_isPeriod (θ : ℝ) (v : Fin 30 → ℝ) (p : ℕ)
    (hθ : (p : ℝ) * θ = 2 * Real.pi) : IsPeriod (cosHistory θ v) p := by
  intro n
  unfold cosHistory scalarHistory
  congr 1
  push_cast
  rw [add_mul, hθ, Real.cos_add_two_pi]

theorem sinHistory_isPeriod (θ : ℝ) (v : Fin 30 → ℝ) (p : ℕ)
    (hθ : (p : ℝ) * θ = 2 * Real.pi) : IsPeriod (sinHistory θ v) p := by
  intro n
  unfold sinHistory scalarHistory
  congr 1
  push_cast
  rw [add_mul, hθ, Real.sin_add_two_pi]

theorem arccos_one_half : Real.arccos (1 / 2) = Real.pi / 3 := by
  rw [← Real.cos_pi_div_three]
  exact Real.arccos_cos (by positivity) (by linarith [Real.pi_pos])

/-- At `h² = 1/2` the `lam = 2` mode advances by `π / 3` per step. -/
theorem modeAngle_half_two (h : ℝ) (hh : h ^ 2 = 1 / 2) :
    modeAngle h 2 = Real.pi / 3 := by
  unfold modeAngle
  rw [hh]
  norm_num
  exact arccos_one_half

/-- At `h² = 2/3` the `lam = 3` mode advances by `π / 2` per step. -/
theorem modeAngle_two_thirds_three (h : ℝ) (hh : h ^ 2 = 2 / 3) :
    modeAngle h 3 = Real.pi / 2 := by
  unfold modeAngle
  rw [hh]
  norm_num

/-- At `h² = 1/3` the `lam = 3` mode advances by `π / 3` per step. -/
theorem modeAngle_third_three (h : ℝ) (hh : h ^ 2 = 1 / 3) :
    modeAngle h 3 = Real.pi / 3 := by
  unfold modeAngle
  rw [hh]
  norm_num
  exact arccos_one_half

theorem twoMode_period_six (h : ℝ) (hh : h ^ 2 = 1 / 2) :
    IsPeriod (cosHistory (modeAngle h 2) twoMode) 6 ∧
      IsPeriod (sinHistory (modeAngle h 2) twoMode) 6 := by
  rw [modeAngle_half_two h hh]
  exact ⟨cosHistory_isPeriod _ _ 6 (by push_cast; ring),
    sinHistory_isPeriod _ _ 6 (by push_cast; ring)⟩

theorem threeMode_period_four (h : ℝ) (hh : h ^ 2 = 2 / 3) :
    IsPeriod (cosHistory (modeAngle h 3) threeMode) 4 ∧
      IsPeriod (sinHistory (modeAngle h 3) threeMode) 4 := by
  rw [modeAngle_two_thirds_three h hh]
  exact ⟨cosHistory_isPeriod _ _ 4 (by push_cast; ring),
    sinHistory_isPeriod _ _ 4 (by push_cast; ring)⟩

theorem threeMode_period_six (h : ℝ) (hh : h ^ 2 = 1 / 3) :
    IsPeriod (cosHistory (modeAngle h 3) threeMode) 6 ∧
      IsPeriod (sinHistory (modeAngle h 3) threeMode) 6 := by
  rw [modeAngle_third_three h hh]
  exact ⟨cosHistory_isPeriod _ _ 6 (by push_cast; ring),
    sinHistory_isPeriod _ _ 6 (by push_cast; ring)⟩

theorem sqrt5_lt_three : Real.sqrt 5 < 3 :=
  (Real.sqrt_lt' (by norm_num)).mpr (by norm_num)

theorem sqrt5_pos : 0 < Real.sqrt 5 := Real.sqrt_pos.mpr (by norm_num)

/-- `h² = 1/2` lies strictly inside the sharp window. -/
theorem window_half (h : ℝ) (hh : h ^ 2 = 1 / 2) :
    h ^ 2 * (3 + Real.sqrt 5) < 4 := by
  rw [hh]; linarith [sqrt5_lt_three]

theorem window_two_thirds (h : ℝ) (hh : h ^ 2 = 2 / 3) :
    h ^ 2 * (3 + Real.sqrt 5) < 4 := by
  rw [hh]; linarith [sqrt5_lt_three]

theorem window_third (h : ℝ) (hh : h ^ 2 = 1 / 3) :
    h ^ 2 * (3 + Real.sqrt 5) < 4 := by
  rw [hh]; linarith [sqrt5_lt_three]

/-- The three periodic cases are oscillators in the sense of
`ModeOscillator`, with the stated periods. -/
theorem periodic_oscillators (h : ℝ) (hh : h ≠ 0) :
    (h ^ 2 = 1 / 2 → ModeOscillator h 2 twoMode ∧
      IsPeriod (cosHistory (modeAngle h 2) twoMode) 6) ∧
    (h ^ 2 = 2 / 3 → ModeOscillator h 3 threeMode ∧
      IsPeriod (cosHistory (modeAngle h 3) threeMode) 4) ∧
    (h ^ 2 = 1 / 3 → ModeOscillator h 3 threeMode ∧
      IsPeriod (cosHistory (modeAngle h 3) threeMode) 6) :=
  ⟨fun h2 ↦ ⟨twoMode_oscillator h hh (by rw [h2]; norm_num),
      (twoMode_period_six h h2).1⟩,
    fun h2 ↦ ⟨threeMode_oscillator h hh (by rw [h2]; norm_num),
      (threeMode_period_four h h2).1⟩,
    fun h2 ↦ ⟨threeMode_oscillator h hh (by rw [h2]; norm_num),
      (threeMode_period_six h h2).1⟩⟩

/-! ## 4. Golden structure of the spectrum -/

theorem three_add_sqrt5_eq_two_goldenRatio_sq :
    3 + Real.sqrt 5 = 2 * Real.goldenRatio ^ 2 := by
  have hs : Real.sqrt 5 * Real.sqrt 5 = 5 := Real.mul_self_sqrt (by norm_num)
  unfold Real.goldenRatio
  linear_combination (-1 / 2) * hs

theorem goldenRatio_sq_pos : 0 < Real.goldenRatio ^ 2 :=
  pow_pos Real.goldenRatio_pos 2

theorem three_sub_sqrt5_eq_two_div_goldenRatio_sq :
    3 - Real.sqrt 5 = 2 / Real.goldenRatio ^ 2 := by
  have hs : Real.sqrt 5 * Real.sqrt 5 = 5 := Real.mul_self_sqrt (by norm_num)
  rw [eq_div_iff goldenRatio_sq_pos.ne']
  unfold Real.goldenRatio
  linear_combination (1 / 4 - Real.sqrt 5 / 4) * hs

/-- The sharp window `h² (3 + √5) < 4` is `h² < 2 / φ²`. -/
theorem window_iff_goldenRatio (h : ℝ) :
    h ^ 2 * (3 + Real.sqrt 5) < 4 ↔ h ^ 2 < 2 / Real.goldenRatio ^ 2 := by
  rw [← three_sub_sqrt5_eq_two_div_goldenRatio_sq]
  have hs : Real.sqrt 5 * Real.sqrt 5 = 5 := Real.mul_self_sqrt (by norm_num)
  have hpos : 0 < 3 + Real.sqrt 5 := by linarith [sqrt5_pos]
  constructor
  · intro hlt
    nlinarith [sqrt5_lt_three]
  · intro hlt
    nlinarith [sqrt5_lt_three]

/-- For `h > 0` the sharp window is `h < √2 / φ`. -/
theorem window_iff_lt_sqrt_two_div_goldenRatio (h : ℝ) (h0 : 0 < h) :
    h ^ 2 * (3 + Real.sqrt 5) < 4 ↔ h < Real.sqrt 2 / Real.goldenRatio := by
  rw [window_iff_goldenRatio]
  have h2 : 2 / Real.goldenRatio ^ 2 = (Real.sqrt 2 / Real.goldenRatio) ^ 2 := by
    have hs := Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
    set g := Real.goldenRatio
    rw [div_pow (Real.sqrt 2) g 2, hs]
  rw [h2]
  exact pow_lt_pow_iff_left₀ h0.le (by positivity) (by norm_num)

theorem projZero_trace : (∑ i : Fin 20, projZeroZ i i) = 20 := by decide
theorem projTwo_trace : (∑ i : Fin 20, projTwoZ i i) = 60 := by decide
theorem projThree_trace : (∑ i : Fin 20, projThreeZ i i) = 40 := by decide
theorem projFive_trace : (∑ i : Fin 20, projFiveZ i i) = 120 := by decide
theorem projGolden_trace : (∑ i : Fin 20, projGoldenZ i i) = 60 := by decide
theorem goldenImage_trace : (∑ i : Fin 20, goldenImageZ i i) = 180 := by decide

theorem trace_scaledProj (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ) :
    Matrix.trace (scaledProj Q d) = ((d : ℝ))⁻¹ * ((∑ i : Fin 20, Q i i : ℤ) : ℝ) := by
  unfold scaledProj
  rw [Matrix.trace_smul, smul_eq_mul, Int.cast_sum]
  rfl

/-- Real traces of the five rational projectors: `1, 5, 4, 4, 6`.  The
multiplicity reading is an inference outside this theorem: with rank equal
to trace for a symmetric idempotent and the committed resolution of the
identity `proj_sum_R`, these traces are the multiplicities `1, 5, 4, 4` of
the eigenvalues `0, 2, 3, 5` of the face normal operator and the dimension
`6` of the golden sector.  Observation, not a theorem here: `1, 5, 4, 3, 3`
are the dimensions of the irreducible representations of `A5`. -/
theorem projector_traces :
    Matrix.trace projZeroR = 1 ∧ Matrix.trace projTwoR = 5 ∧
      Matrix.trace projThreeR = 4 ∧ Matrix.trace projFiveR = 4 ∧
      Matrix.trace projGoldenR = 6 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · unfold projZeroR; rw [trace_scaledProj, projZero_trace]; norm_num
  · unfold projTwoR; rw [trace_scaledProj, projTwo_trace]; norm_num
  · unfold projThreeR; rw [trace_scaledProj, projThree_trace]; norm_num
  · unfold projFiveR; rw [trace_scaledProj, projFive_trace]; norm_num
  · unfold projGoldenR; rw [trace_scaledProj, projGolden_trace]; norm_num

/-- The trace of `N` on the golden sector is `18 = 3 (3 + √5) + 3 (3 - √5)`;
read with the sector dimension `6` (rank equal to trace), the two golden
eigenvalues each have multiplicity `3`, an inference outside this theorem. -/
theorem golden_sector_normal_trace :
    Matrix.trace (faceNormalR * projGoldenR) = 18 := by
  have h1 : castZ faceNormalZ * castZ projGoldenZ = castZ goldenImageZ :=
    castZ_mul_eq _ _ _ projGolden_image
  unfold projGoldenR scaledProj faceNormalR
  rw [Matrix.mul_smul, h1, Matrix.trace_smul, smul_eq_mul]
  have : Matrix.trace (castZ goldenImageZ) = ((∑ i : Fin 20, goldenImageZ i i : ℤ) : ℝ) := by
    rw [Int.cast_sum]; rfl
  rw [this, goldenImage_trace]
  norm_num

/-! ## 5. Tick-free frequency ratios and the small-step bound -/

/-- **Totalized algebra identity (compatibility lemma).**  Cancelling the
common nonzero tick factor gives the same real-division expression for all
angles.  When `θ₂ = 0`, both sides use Lean's totalized division and this is
not a physically defined frequency ratio.  Use
`labFrequency_ratio_of_ne_zero` for the citable ratio statement. -/
theorem labFrequency_ratio (cal : ClockCalibration) (θ₁ θ₂ : ℝ) :
    cal.labFrequency θ₁ / cal.labFrequency θ₂ = θ₁ / θ₂ := by
  simp only [ClockCalibration.labFrequency_def]
  have hc : 2 * Real.pi * cal.tau ≠ 0 := by
    have := cal.tau_pos; have := Real.pi_pos; positivity
  exact div_div_div_cancel_right₀ hc θ₁ θ₂

/-- A nonzero angular advance gives a nonzero laboratory frequency under
every positive declared tick. -/
theorem labFrequency_ne_zero (cal : ClockCalibration) (θ : ℝ) (hθ : θ ≠ 0) :
    cal.labFrequency θ ≠ 0 := by
  simp only [ClockCalibration.labFrequency_def]
  have hc : 2 * Real.pi * cal.tau ≠ 0 := by
    have := cal.tau_pos; have := Real.pi_pos; positivity
  exact div_ne_zero hθ hc

/-- **Citable laboratory-frequency ratio.**  If the denominator angle is
nonzero, then the denominator laboratory frequency is nonzero and its ratio
with the numerator frequency is the ratio of angular advances per step. -/
theorem labFrequency_ratio_of_ne_zero (cal : ClockCalibration) (θ₁ θ₂ : ℝ)
    (hθ₂ : θ₂ ≠ 0) :
    cal.labFrequency θ₁ / cal.labFrequency θ₂ = θ₁ / θ₂ := by
  have hfreq₂ : cal.labFrequency θ₂ ≠ 0 := labFrequency_ne_zero cal θ₂ hθ₂
  rw [div_eq_iff hfreq₂]
  simp only [ClockCalibration.labFrequency_def]
  have hc : 2 * Real.pi * cal.tau ≠ 0 := by
    have := cal.tau_pos; have := Real.pi_pos; positivity
  field_simp [hθ₂, hc]

/-- **Totalized tick-cancellation identity (compatibility lemma).**  This is
valid as real algebra for every `θ₂`, but at `θ₂ = 0` it is not a physical
frequency ratio.  Use `labFrequency_ratio_tick_free_of_ne_zero` for that
statement. -/
theorem labFrequency_ratio_tick_free (cal cal' : ClockCalibration) (θ₁ θ₂ : ℝ) :
    cal.labFrequency θ₁ / cal.labFrequency θ₂ =
      cal'.labFrequency θ₁ / cal'.labFrequency θ₂ := by
  rw [labFrequency_ratio, labFrequency_ratio]

/-- A laboratory-frequency ratio with nonzero denominator carries no declared
tick. -/
theorem labFrequency_ratio_tick_free_of_ne_zero
    (cal cal' : ClockCalibration) (θ₁ θ₂ : ℝ) (hθ₂ : θ₂ ≠ 0) :
    cal.labFrequency θ₁ / cal.labFrequency θ₂ =
      cal'.labFrequency θ₁ / cal'.labFrequency θ₂ := by
  rw [labFrequency_ratio_of_ne_zero cal θ₁ θ₂ hθ₂,
    labFrequency_ratio_of_ne_zero cal' θ₁ θ₂ hθ₂]

theorem modeAngle_nonneg (h lam : ℝ) : 0 ≤ modeAngle h lam :=
  Real.arccos_nonneg _

theorem modeAngle_le_pi (h lam : ℝ) : modeAngle h lam ≤ Real.pi :=
  Real.arccos_le_pi _

theorem modeAngle_lt_pi (h lam : ℝ) (h4 : h ^ 2 * lam < 4) :
    modeAngle h lam < Real.pi := by
  unfold modeAngle
  exact Real.arccos_lt_pi.mpr (by linarith)

/-- Half-angle sine: `sin (θ / 2) = h √lam / 2`. -/
theorem sin_half_modeAngle (h lam : ℝ) (h0 : 0 ≤ h) (hl : 0 ≤ lam)
    (h4 : h ^ 2 * lam ≤ 4) :
    Real.sin (modeAngle h lam / 2) = h * Real.sqrt lam / 2 := by
  have hx0 : 0 ≤ h ^ 2 * lam := by positivity
  have hcos := cos_modeAngle h lam hx0 h4
  have hθ : modeAngle h lam = 2 * (modeAngle h lam / 2) := by ring
  rw [hθ, Real.cos_two_mul, Real.cos_sq'] at hcos
  have hsq : Real.sin (modeAngle h lam / 2) ^ 2 = (h * Real.sqrt lam / 2) ^ 2 := by
    have hl2 : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hl
    nlinarith [hcos, hl2]
  have hpi : modeAngle h lam ≤ Real.pi := modeAngle_le_pi h lam
  have hnn : 0 ≤ modeAngle h lam := modeAngle_nonneg h lam
  have hs0 : 0 ≤ Real.sin (modeAngle h lam / 2) :=
    Real.sin_nonneg_of_nonneg_of_le_pi (by linarith) (by linarith)
  have hr0 : 0 ≤ h * Real.sqrt lam / 2 := by positivity
  exact (sq_eq_sq₀ hs0 hr0).mp hsq

/-- Half-angle cosine: `cos (θ / 2) = √(1 - h² lam / 4)`. -/
theorem cos_half_modeAngle (h lam : ℝ) (h0 : 0 ≤ h) (hl : 0 ≤ lam)
    (h4 : h ^ 2 * lam ≤ 4) :
    Real.cos (modeAngle h lam / 2) = Real.sqrt (1 - h ^ 2 * lam / 4) := by
  have hsin := sin_half_modeAngle h lam h0 hl h4
  have hpi : modeAngle h lam ≤ Real.pi := modeAngle_le_pi h lam
  have hnn : 0 ≤ modeAngle h lam := modeAngle_nonneg h lam
  have hc0 : 0 ≤ Real.cos (modeAngle h lam / 2) :=
    Real.cos_nonneg_of_neg_pi_div_two_le_of_le (by linarith) (by linarith)
  have hl2 : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hl
  have hsq : Real.cos (modeAngle h lam / 2) ^ 2 = 1 - h ^ 2 * lam / 4 := by
    rw [Real.cos_sq', hsin]
    nlinarith [hl2]
  rw [← hsq, Real.sqrt_sq hc0]

/-- Lower bound: `h √lam ≤ modeAngle h lam`. -/
theorem modeAngle_ge (h lam : ℝ) (h0 : 0 ≤ h) (hl : 0 ≤ lam)
    (h4 : h ^ 2 * lam ≤ 4) :
    h * Real.sqrt lam ≤ modeAngle h lam := by
  have hsin := sin_half_modeAngle h lam h0 hl h4
  have hu0 : 0 ≤ modeAngle h lam / 2 := by linarith [modeAngle_nonneg h lam]
  have hle : Real.sin (modeAngle h lam / 2) ≤ modeAngle h lam / 2 := by
    rcases hu0.lt_or_eq with hlt | heq
    · exact (Real.sin_lt hlt).le
    · rw [← heq, Real.sin_zero]
  linarith

/-- Upper bound: `modeAngle h lam ≤ h √lam / √(1 - h² lam / 4)` strictly
inside the window. -/
theorem modeAngle_le (h lam : ℝ) (h0 : 0 ≤ h) (hl : 0 ≤ lam)
    (h4 : h ^ 2 * lam < 4) :
    modeAngle h lam ≤ h * Real.sqrt lam / Real.sqrt (1 - h ^ 2 * lam / 4) := by
  have hsin := sin_half_modeAngle h lam h0 hl h4.le
  have hcos := cos_half_modeAngle h lam h0 hl h4.le
  have hu0 : 0 ≤ modeAngle h lam / 2 := by linarith [modeAngle_nonneg h lam]
  have hu1 : modeAngle h lam / 2 < Real.pi / 2 := by
    linarith [modeAngle_lt_pi h lam h4]
  have htan := Real.le_tan hu0 hu1
  rw [Real.tan_eq_sin_div_cos, hsin, hcos] at htan
  have hpos : 0 < Real.sqrt (1 - h ^ 2 * lam / 4) := Real.sqrt_pos.mpr (by linarith)
  rw [le_div_iff₀ hpos] at htan ⊢
  linarith

/-- The frequency ratio of two modes at one step is squeezed between
explicit multiples of `√lam₁ / √lam₂`. -/
theorem modeAngle_ratio_bounds (h lam₁ lam₂ : ℝ) (h0 : 0 < h) (hl₁ : 0 ≤ lam₁)
    (hl₂ : 0 < lam₂) (h₁ : h ^ 2 * lam₁ < 4) (h₂ : h ^ 2 * lam₂ < 4) :
    Real.sqrt (1 - h ^ 2 * lam₂ / 4) * (Real.sqrt lam₁ / Real.sqrt lam₂) ≤
        modeAngle h lam₁ / modeAngle h lam₂ ∧
      modeAngle h lam₁ / modeAngle h lam₂ ≤
        (Real.sqrt lam₁ / Real.sqrt lam₂) / Real.sqrt (1 - h ^ 2 * lam₁ / 4) := by
  have hge₁ := modeAngle_ge h lam₁ h0.le hl₁ h₁.le
  have hle₁ := modeAngle_le h lam₁ h0.le hl₁ h₁
  have hge₂ := modeAngle_ge h lam₂ h0.le hl₂.le h₂.le
  have hle₂ := modeAngle_le h lam₂ h0.le hl₂.le h₂
  have hs₂ : 0 < Real.sqrt lam₂ := Real.sqrt_pos.mpr hl₂
  have hs₁ : 0 ≤ Real.sqrt lam₁ := Real.sqrt_nonneg _
  have hθ₂ : 0 < modeAngle h lam₂ := by nlinarith
  have hq₁ : 0 < Real.sqrt (1 - h ^ 2 * lam₁ / 4) := Real.sqrt_pos.mpr (by linarith)
  have hq₂ : 0 < Real.sqrt (1 - h ^ 2 * lam₂ / 4) := Real.sqrt_pos.mpr (by linarith)
  constructor
  · rw [le_div_iff₀ hθ₂]
    have hle₂' := (le_div_iff₀ hq₂).mp hle₂
    calc Real.sqrt (1 - h ^ 2 * lam₂ / 4) * (Real.sqrt lam₁ / Real.sqrt lam₂) *
          modeAngle h lam₂
        ≤ Real.sqrt (1 - h ^ 2 * lam₂ / 4) * (Real.sqrt lam₁ / Real.sqrt lam₂) *
          (h * Real.sqrt lam₂ / Real.sqrt (1 - h ^ 2 * lam₂ / 4)) := by
          apply mul_le_mul_of_nonneg_left hle₂
          positivity
      _ = h * Real.sqrt lam₁ * (Real.sqrt (1 - h ^ 2 * lam₂ / 4) /
          Real.sqrt (1 - h ^ 2 * lam₂ / 4)) * (Real.sqrt lam₂ / Real.sqrt lam₂) := by
          ring
      _ = h * Real.sqrt lam₁ := by
          rw [div_self hq₂.ne', div_self hs₂.ne', mul_one, mul_one]
      _ ≤ modeAngle h lam₁ := hge₁
  · rw [div_le_iff₀ hθ₂]
    calc modeAngle h lam₁
        ≤ h * Real.sqrt lam₁ / Real.sqrt (1 - h ^ 2 * lam₁ / 4) := hle₁
      _ = (Real.sqrt lam₁ / Real.sqrt lam₂) / Real.sqrt (1 - h ^ 2 * lam₁ / 4) *
          (h * Real.sqrt lam₂) := by
          have : (Real.sqrt lam₁ / Real.sqrt lam₂) / Real.sqrt (1 - h ^ 2 * lam₁ / 4) *
              (h * Real.sqrt lam₂) = h * Real.sqrt lam₁ / Real.sqrt (1 - h ^ 2 * lam₁ / 4) *
              (Real.sqrt lam₂ / Real.sqrt lam₂) := by ring
          rw [this, div_self hs₂.ne', mul_one]
      _ ≤ (Real.sqrt lam₁ / Real.sqrt lam₂) / Real.sqrt (1 - h ^ 2 * lam₁ / 4) *
          modeAngle h lam₂ := by
          apply mul_le_mul_of_nonneg_left hge₂
          positivity

/-- The small-step limit `√(lam₁ / lam₂)` of the frequency ratio, from
`modeAngle_div_tendsto`, equals `φ²` for the golden pair:
`√((3 + √5) / (3 - √5)) = φ²`. -/
theorem sqrt_golden_ratio_eq_goldenRatio_sq :
    Real.sqrt ((3 + Real.sqrt 5) / (3 - Real.sqrt 5)) = Real.goldenRatio ^ 2 := by
  rw [three_add_sqrt5_eq_two_goldenRatio_sq, three_sub_sqrt5_eq_two_div_goldenRatio_sq]
  have hg := goldenRatio_sq_pos
  have : 2 * Real.goldenRatio ^ 2 / (2 / Real.goldenRatio ^ 2) =
      (Real.goldenRatio ^ 2) ^ 2 := by
    field_simp
  rw [this, Real.sqrt_sq hg.le]

/-- Small-step limit: `modeAngle h lam / h → √lam` as `h → 0⁺`. -/
theorem modeAngle_div_tendsto (lam : ℝ) (hl : 0 ≤ lam) :
    Filter.Tendsto (fun h : ℝ ↦ modeAngle h lam / h) (nhdsWithin 0 (Set.Ioi 0))
      (nhds (Real.sqrt lam)) := by
  have hδ : (0 : ℝ) < 1 / (lam + 1) := by positivity
  have hmem : Set.Ioo (0 : ℝ) (1 / (lam + 1)) ∈ nhdsWithin (0 : ℝ) (Set.Ioi 0) :=
    Ioo_mem_nhdsGT hδ
  have hwin : ∀ h ∈ Set.Ioo (0 : ℝ) (1 / (lam + 1)), h ^ 2 * lam < 4 := by
    intro h hh
    obtain ⟨h0, h1⟩ := hh
    have h1' : h * (lam + 1) < 1 := by
      rwa [lt_div_iff₀ (by positivity)] at h1
    have hh1 : h < 1 := by nlinarith
    nlinarith [mul_nonneg h0.le hl]
  have hupper : Filter.Tendsto
      (fun h : ℝ ↦ Real.sqrt lam / Real.sqrt (1 - h ^ 2 * lam / 4))
      (nhdsWithin 0 (Set.Ioi 0)) (nhds (Real.sqrt lam)) := by
    have hcont : Filter.Tendsto
        (fun h : ℝ ↦ Real.sqrt lam / Real.sqrt (1 - h ^ 2 * lam / 4)) (nhds 0)
        (nhds (Real.sqrt lam / Real.sqrt (1 - (0 : ℝ) ^ 2 * lam / 4))) := by
      apply Filter.Tendsto.div tendsto_const_nhds
      · exact (Real.continuous_sqrt.comp (by fun_prop)).tendsto 0
      · simp
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, zero_mul,
      zero_div, sub_zero, Real.sqrt_one, div_one] at hcont
    exact hcont.mono_left nhdsWithin_le_nhds
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hupper ?_ ?_
  · filter_upwards [hmem] with h hh
    have h0 : 0 < h := hh.1
    rw [le_div_iff₀ h0]
    have := modeAngle_ge h lam h0.le hl (hwin h hh).le
    linarith
  · filter_upwards [hmem] with h hh
    have h0 : 0 < h := hh.1
    rw [div_le_iff₀ h0]
    have := modeAngle_le h lam h0.le hl (hwin h hh)
    have hq : 0 < Real.sqrt (1 - h ^ 2 * lam / 4) :=
      Real.sqrt_pos.mpr (by linarith [hwin h hh])
    calc modeAngle h lam ≤ h * Real.sqrt lam / Real.sqrt (1 - h ^ 2 * lam / 4) := this
      _ = Real.sqrt lam / Real.sqrt (1 - h ^ 2 * lam / 4) * h := by ring

/-! ## 6. Commensurability with the source recurrence clock in one index -/

/-- For any declared tick, the laboratory duration of `p` steps over the
laboratory mean return time of the source recurrence clock from state `0`
is `7155 p / 61511`.  The ratio depends on the declared step `h` only
through the period `p`; the shared step index is a count, not a physical
duration. -/
theorem period_return_ratio (cal : ClockCalibration) (p : ℕ) :
    cal.labSeconds (p : ℝ) / cal.labSeconds (meanReturn 0) =
      7155 * (p : ℝ) / 61511 := by
  simp only [ClockCalibration.labSeconds_def]
  rw [mul_div_mul_right _ _ cal.tau_pos.ne', meanReturn_values.1]
  rw [div_div_eq_mul_div]
  ring

/-- The period-to-return ratio carries no tick. -/
theorem period_return_ratio_tick_free (cal cal' : ClockCalibration) (p : ℕ) :
    cal.labSeconds (p : ℝ) / cal.labSeconds (meanReturn 0) =
      cal'.labSeconds (p : ℝ) / cal'.labSeconds (meanReturn 0) := by
  rw [period_return_ratio, period_return_ratio]

/-- Exact values at the periods `4` and `6`. -/
theorem period_return_ratio_values (cal : ClockCalibration) :
    cal.labSeconds (4 : ℝ) / cal.labSeconds (meanReturn 0) = 28620 / 61511 ∧
      cal.labSeconds (6 : ℝ) / cal.labSeconds (meanReturn 0) = 42930 / 61511 := by
  constructor
  · have h := period_return_ratio cal 4
    push_cast at h
    rw [h]; norm_num
  · have h := period_return_ratio cal 6
    push_cast at h
    rw [h]; norm_num

/-- In the joined architecture the join at index `n + p` differs from the
join at index `n` for every positive `p`: a field period is a return of the
configuration inside one shared index, never a return of the joined record
(`join_injective`). -/
theorem join_ne_of_period (J : MaxwellClockJoinedArchitecture) (n p : ℕ)
    (hp : 0 < p) : J.join (n + p) ≠ J.join n := by
  intro h
  have := join_injective J h
  omega

/-- Composite receipt for the period-`6` mode at `h² = 1/2`: the field
returns every `6` steps, the `lam = 2` oscillator packet holds, the step is
inside the sharp window, and the step-count ratio to the mean return time is
`42930 / 61511` for every declared tick. -/
theorem twoMode_commensurability (h : ℝ) (hh : h ≠ 0) (h2 : h ^ 2 = 1 / 2)
    (cal : ClockCalibration) :
    IsPeriod (cosHistory (modeAngle h 2) twoMode) 6 ∧
      ModeOscillator h 2 twoMode ∧
      h ^ 2 * (3 + Real.sqrt 5) < 4 ∧
      cal.labSeconds (6 : ℝ) / cal.labSeconds (meanReturn 0) = 42930 / 61511 :=
  ⟨(twoMode_period_six h h2).1, twoMode_oscillator h hh (by rw [h2]; norm_num),
    window_half h h2, (period_return_ratio_values cal).2⟩

/-- Composite receipt for the period-`4` mode at `h² = 2/3`. -/
theorem threeMode_commensurability (h : ℝ) (hh : h ≠ 0) (h2 : h ^ 2 = 2 / 3)
    (cal : ClockCalibration) :
    IsPeriod (cosHistory (modeAngle h 3) threeMode) 4 ∧
      ModeOscillator h 3 threeMode ∧
      h ^ 2 * (3 + Real.sqrt 5) < 4 ∧
      cal.labSeconds (4 : ℝ) / cal.labSeconds (meanReturn 0) = 28620 / 61511 :=
  ⟨(threeMode_period_four h h2).1, threeMode_oscillator h hh (by rw [h2]; norm_num),
    window_two_thirds h h2, (period_return_ratio_values cal).1⟩

end

end OPH.CarrierModeOscillators

#print axioms OPH.CarrierModeOscillators.two_sub_two_cos_modeAngle
#print axioms OPH.CarrierModeOscillators.scalarHistory_ampere
#print axioms OPH.CarrierModeOscillators.cosHistory_ampere
#print axioms OPH.CarrierModeOscillators.sinHistory_ampere
#print axioms OPH.CarrierModeOscillators.cosHistory_energy_conserved
#print axioms OPH.CarrierModeOscillators.cosHistory_electric_bounded
#print axioms OPH.CarrierModeOscillators.modeOscillator
#print axioms OPH.CarrierModeOscillators.fiveMode_oscillator
#print axioms OPH.CarrierModeOscillators.goldenMode_oscillator
#print axioms OPH.CarrierModeOscillators.codifferential_eigen
#print axioms OPH.CarrierModeOscillators.twoMode_eigen
#print axioms OPH.CarrierModeOscillators.threeMode_eigen
#print axioms OPH.CarrierModeOscillators.twoMode_ne_zero
#print axioms OPH.CarrierModeOscillators.threeMode_ne_zero
#print axioms OPH.CarrierModeOscillators.twoMode_period_six
#print axioms OPH.CarrierModeOscillators.threeMode_period_four
#print axioms OPH.CarrierModeOscillators.threeMode_period_six
#print axioms OPH.CarrierModeOscillators.periodic_oscillators
#print axioms OPH.CarrierModeOscillators.window_half
#print axioms OPH.CarrierModeOscillators.three_add_sqrt5_eq_two_goldenRatio_sq
#print axioms OPH.CarrierModeOscillators.three_sub_sqrt5_eq_two_div_goldenRatio_sq
#print axioms OPH.CarrierModeOscillators.window_iff_goldenRatio
#print axioms OPH.CarrierModeOscillators.window_iff_lt_sqrt_two_div_goldenRatio
#print axioms OPH.CarrierModeOscillators.projector_traces
#print axioms OPH.CarrierModeOscillators.golden_sector_normal_trace
#print axioms OPH.CarrierModeOscillators.labFrequency_ne_zero
#print axioms OPH.CarrierModeOscillators.labFrequency_ratio_of_ne_zero
#print axioms OPH.CarrierModeOscillators.labFrequency_ratio_tick_free
#print axioms OPH.CarrierModeOscillators.labFrequency_ratio_tick_free_of_ne_zero
#print axioms OPH.CarrierModeOscillators.modeAngle_ge
#print axioms OPH.CarrierModeOscillators.modeAngle_le
#print axioms OPH.CarrierModeOscillators.modeAngle_ratio_bounds
#print axioms OPH.CarrierModeOscillators.sqrt_golden_ratio_eq_goldenRatio_sq
#print axioms OPH.CarrierModeOscillators.modeAngle_div_tendsto
#print axioms OPH.CarrierModeOscillators.period_return_ratio
#print axioms OPH.CarrierModeOscillators.period_return_ratio_tick_free
#print axioms OPH.CarrierModeOscillators.period_return_ratio_values
#print axioms OPH.CarrierModeOscillators.join_ne_of_period
#print axioms OPH.CarrierModeOscillators.twoMode_commensurability
#print axioms OPH.CarrierModeOscillators.threeMode_commensurability
