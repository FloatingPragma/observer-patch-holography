import Geometry.ProperTimeInternalAction
import Geometry.WorldlineHopTransport
import Geometry.CommonWorldMaxwellClockJoin
import QFT.SourceRecurrenceClock

set_option autoImplicit false

open scoped BigOperators

/-!
# The source clock rate along seam-step worldlines and refinement invariance
(issues 736, 739)

STATUS.  Comparison module on the source clock and duration row.  Two
declared accrual rules for an internal process hosted on a seam-step
worldline are placed side by side: index accrual, one internal step per
shared step index at the declared unit `τ` per index (declared here as
the reading of the abstract step index of `QFT/SourceRecurrenceClock.lean`
as the shared step index of the join in
`Geometry/CommonWorldMaxwellClockJoin.lean`; the join's own `stepClock`
advances per index at a declared `stepDuration`, and that duration is not
forced, `joined_stepDuration_not_forced`), and proper-length accrual, the
declared proper-time principle of `Geometry/ProperTimeInternalAction.lean`.  The module proves the exact
difference of the two rules on every window, the exact behaviour of each
under the declared midpoint refinement, the exact dilation factor of a
uniformly moving seam-step worldline, and the exact mean return time of
the recurrence clock under each rule.  Both rules are declared; nothing
below selects one of them from the source.  The time-dilation reading of
section (3) is a declared identification of proper length with elapsed
internal time, and the factor is normalised by the proper length of the
resting worldline over the same index count (division by `τ` relative to
proper length per index step).

WHAT IS PROVED.

(1) Two accrual rules and their exact difference.  `indexAccrual τ E M x
= internalAction 0 τ E M x = E * (M * τ)`, index accrual at the declared
unit `τ` per index (the member `a = 0, b = τ`), and `lengthAccrual E M x =
internalAction 1 0 E M x = E * properLength M x` (the member `a = 1,
b = 0`).  The comparison is at the same unit: `E M τ` against
`E properLength`.  On the generated path
of a seam-step worldline at unit `τ ≥ 0` the proper length over `M`
steps is `c * √(τ² - 4) + r * τ` with `c` crossings and `r` rests
(`properLength_generated_counts`), and the two accruals differ by exactly
`E * c * (τ - √(τ² - 4))` (`accrual_difference`): rests contribute `E τ`
under both rules, crossings differ.  For `E ≠ 0` and `2 < τ` the difference vanishes exactly when the
window has no crossing (`accrual_difference_eq_zero_iff`,
`accruals_agree_iff_resting`), and it is nonzero on every window with a
crossing (`accrual_difference_ne_zero`).

(2) Refinement.  The midpoint refinement of a generated path has forward
increments of Lorentz square `(τ² - stepNormSq) / 4`, that is
`(τ² - 4) / 4` at a crossing and `τ² / 4` at a rest
(`refine_generated_lorentzQ`), and spatial increments of norm squared
`1` at a crossing (`refine_generated_spatialNormSq`).  No seam step has
spatial norm squared `1` (`no_seam_step_norm_one`), so the refinement of
a worldline with a crossing is the generated path of no seam-step
worldline at any unit (`refine_not_generated`): the declared refinement
acts on the Lorentz-module image and leaves the port class.  Index
accrual is the member `a = 0, b = τ` of the declared family and fails
refinement invariance at every nonzero `E` and nonzero `τ`
(`indexAccrual_not_refinementInvariant`, through
`refinementInvariant_iff_b_zero`; at `τ = 1` this is
`additive_term_not_refinementInvariant` of the imported module); on the
image it doubles under one refinement (`indexAccrual_refine`).  Proper-length accrual is the member
`a = 1, b = 0`, passes (`lengthAccrual_refinementInvariant`), and is
fixed by one refinement on the image (`lengthAccrual_refine`).  The
composite statement `index_accrual_compatible_iff_resting`: for `E ≠ 0`
and `2 < τ`, index accrual at unit `τ` (a non-invariant member) takes the
value of the refinement-invariant member on a window exactly when the
window has no crossing.  The negative `rules_both_declared`: index accrual
at a fixed unit is the same function of the window on every worldline
(its definition ignores the path), proper-length accrual separates the
one-step crossing worldline from a resting worldline, so the two rules
are distinct functionals on the class.  That neither rule is selected by
the source chain or the join is a scope statement of this header, proved
nowhere below.

(3) The dilation factor.  For the uniform worldline `uniformWorldline e
k` (one crossing of seam `e` every `k + 1` steps, alternating forward and
backward, `k` rests between) the crossing count over `p` periods is `p`
(`crossingCount_uniform`) and the proper length over `p` periods is
`p * (k * τ + √(τ² - 4))` (`properLength_uniform`).  With
`dilationFactor τ k = (k * τ + √(τ² - 4)) / ((k + 1) * τ)` the proper
length over `p` periods is `dilationFactor τ k` times the proper length
of a resting worldline over the same index count, `(k + 1) * p * τ`
(`properLength_uniform_dilation`).  The factor is strictly below `1` for
`2 < τ` (`dilationFactor_lt_one`), strictly positive
(`dilationFactor_pos`), tends to `1` as `k` grows
(`dilationFactor_tendsto_one`), and equals `√5 / 3` at `τ = 3, k = 0`
(`dilationFactor_three_zero`).

(4) Mean return under the two rules, under declared readings.  Under the
declared reading of one internal unit as one shared step (index rule) the
mean return time of the committed stable class is `61511/7155` shared
steps, a definition that ignores the worldline (`indexRule_meanReturn`).  Under proper-length accrual the same
number counts units of proper length in the rest step `τ`, and on the
uniform worldline the corresponding index count is
`(61511/7155) * ((k + 1) * τ) / (k * τ + √(τ² - 4))`
(`lengthRule_return_indexSteps`, a real number, not a step count), the
inverse dilation factor times the index-rule value (`lengthRule_return_eq_inv_dilation`); the proper
length accrued over that many index steps, read on the uniform worldline
through `properLength_uniform_dilation`, is the index-rule value times
`τ` (`lengthRule_return_consistent`).  Both readings are declared; no theorem
relates `returnCount` or `meanReturn` to a worldline, a window, or
`properLength`, and the content of this section is `meanReturn_values`
together with the algebra of `dilationFactor`.

PRIOR WORK.  `Geometry/ProperTimeInternalAction.lean` proves
`refinementInvariant_iff_b_zero`, `additive_term_not_refinementInvariant`
(`¬ RefinementInvariant 0 1 E`, reused below for the unit-`1` case),
`properLength_term_refinementInvariant`, `properLength_refine`, and the
slope selection; its `restPath` is the unit-`1` rest worldline.
`Geometry/WorldlineHopTransport.lean` proves `clockAction_generated` (the
quadratic functional, `τ² - 4` per crossing) and
`seam_step_timelike_iff`; the proper-length functional of a generated
path appears here for the first time.  `QFT/SourceRecurrenceClock.lean`
proves `meanReturn_values` on the abstract step index and
`returnCount_block_additive`; the crossing count below is the same
distinct-event count on a step sequence.
`Geometry/CommonWorldMaxwellClockJoin.lean` proves `joined_clock_advance`
(the `stepClock` phase advances by `mass * (stepTime δ m - stepTime δ n)`
at a declared `stepDuration δ`) and `joined_stepDuration_not_forced` (two
joined inhabitants over one bridged base with different `δ`): the join's
own clock is index accrual at a declared unit per index, and the join
contains no source chain.  The per-index advance of the source chain
along a worldline is declared in this module, proved nowhere.  No prior module states a dilation factor for a seam-step
worldline or compares the two accrual rules.

ROWS TOUCHED.  The source clock and duration row (both accrual rules are
declared, the unit `τ` is declared, the tick of one chain step stays a
declared calibration); the coupled-action row (the internal-action family
is a declared enrichment); the physical spacetime attachment row (the
Lorentz module is the declared image of the port class; no port is
attached to a point); the light-signal row (no signal propagation is
attached to the threshold `τ = 2`); the laboratory clock and energy
calibration import (no unit, calibration, or readout is attached to `τ`
or `E`); the gravitation-route energy identification (no identification
of `E` with an energy is made).  The module discharges none of these
rows (`dischargedRows_empty`).

NEGATIVES CITED.  The Legendre non-identifiability at scope
(`legendre_scope_cited`, re-cited from
`Geometry/InternalEnergyInertia.lean`): realized histories select no
velocity curvature or Legendre map, so every accrual shape is a declared
enrichment.  The rate non-identifiability at scope
(`rate_scope_cited`, re-cited from
`ObserverPatchHolography/RateNonidentifiability.lean`): the declared
reading of an abstract transition system fixes no first-locking count;
it forbids nothing about a process hosted on the source chain, and it is
cited here for scope only.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - spatialNormSq v.2`; forward differences; the
generated path advances `τ` per step in the scalar coordinate and one
signed seam vector of norm squared `4` per crossing.  `properLength` is
the sum of real square roots of the Lorentz squares of the forward
increments (clamped to zero below the threshold).  A window of `M` steps
is the index set `0, …, M - 1`.  `crossingCount σ M` counts the steps
`k < M` with `σ k ≠ rest`.  The uniform worldline crosses at every index
divisible by `k + 1`, forward on even crossings and backward on odd ones,
so it stays admissible on the two endpoints of one seam.

FALSIFIER.  The module is wrong if some crossing step of a generated
path has proper length other than `√(τ² - 4)`, if the refined path of a
worldline with a crossing is the generated path of some seam-step
worldline, if the uniform worldline is inadmissible at some step, if the
dilation factor at `τ = 3, k = 0` differs from `√5 / 3`, if
`RefinementInvariant 0 τ E` holds at some `τ ≠ 0`, `E ≠ 0`, or if the mean
return of the committed stable class differs from `61511/7155`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

namespace OPH.SourceClockRateAlongWorldlines

open OPH.C1Lorentz (Spatial Herm2 spatialNormSq lorentzQ)
open OPH.ProperTimeInternalAction
open OPH.WorldlineHopTransport
open OPH.InternalEnergyInertia (OpenRow)
open OPH.SeamCurrentCarrierQuotient (seamLeft seamRight)
open OPH.QFT.SourceRecurrenceClock (meanReturn meanReturn_values)

noncomputable section

/-! ## (1) Two declared accrual rules on a seam-step worldline -/

/-- Index accrual (declared): one internal step per shared step index at
the declared unit `τ` per index.  It is the declared reading of the
abstract chain step index of `SourceRecurrenceClock` as the shared step
index of the join; the join's own `stepClock` advances per index at a
declared `stepDuration` (`joined_clock_advance`,
`joined_stepDuration_not_forced`).  It is the member `a = 0, b = τ` of the
declared family, and its value ignores the path. -/
def indexAccrual (τ E : ℝ) (M : ℕ) (x : ℕ → Herm2) : ℝ := internalAction 0 τ E M x

/-- Proper-length accrual (declared): the internal process advances per
unit proper length, the proper-time principle.  It is the member
`a = 1, b = 0` of the declared family. -/
def lengthAccrual (E : ℝ) (M : ℕ) (x : ℕ → Herm2) : ℝ := internalAction 1 0 E M x

theorem indexAccrual_eq (τ E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    indexAccrual τ E M x = E * (M * τ) := by
  unfold indexAccrual internalAction
  ring

theorem lengthAccrual_eq (E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    lengthAccrual E M x = E * properLength M x := by
  unfold lengthAccrual internalAction
  ring

theorem lengthAccrual_eq_properTime (E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    lengthAccrual E M x = properTimeInternalAction E M x :=
  (properTimeInternalAction_eq E M x).symm

/-- The number of crossing steps among the first `M` steps of a step
sequence: a distinct-event count in the sense of
`SourceRecurrenceClock.returnCount`. -/
def crossingCount (σ : ℕ → SeamStep) : ℕ → ℕ
  | 0 => 0
  | n + 1 => crossingCount σ n + if σ n ≠ .rest then 1 else 0

/-- The number of rest steps among the first `M` steps. -/
def restCount (σ : ℕ → SeamStep) : ℕ → ℕ
  | 0 => 0
  | n + 1 => restCount σ n + if σ n = .rest then 1 else 0

theorem crossingCount_add_restCount (σ : ℕ → SeamStep) (M : ℕ) :
    crossingCount σ M + restCount σ M = M := by
  induction M with
  | zero => rfl
  | succ M ih =>
    simp only [crossingCount, restCount]
    by_cases h : σ M = .rest <;>
      simp only [h, ne_eq, not_true_eq_false, not_false_eq_true, if_true, if_false] <;> omega

theorem crossingCount_block_additive (σ : ℕ → SeamStep) (n m : ℕ) :
    crossingCount σ (n + m) = crossingCount σ n + crossingCount (fun j => σ (n + j)) m := by
  induction m with
  | zero => simp [crossingCount]
  | succ m ih =>
    rw [← Nat.add_assoc]
    simp only [crossingCount]
    rw [ih]
    ring

theorem crossingCount_eq_zero_iff (σ : ℕ → SeamStep) (M : ℕ) :
    crossingCount σ M = 0 ↔ ∀ k, k < M → σ k = .rest := by
  induction M with
  | zero => simp [crossingCount]
  | succ M ih =>
    simp only [crossingCount]
    constructor
    · intro h k hk
      have h1 : crossingCount σ M = 0 := by omega
      have h2 : (if σ M ≠ .rest then 1 else 0) = 0 := by omega
      rcases Nat.lt_succ_iff_lt_or_eq.mp hk with hk | hk
      · exact ih.mp h1 k hk
      · subst hk
        by_contra hne
        simp [hne] at h2
    · intro h
      rw [ih.mpr (fun k hk => h k (Nat.lt_succ_of_lt hk))]
      simp [h M (Nat.lt_succ_self M)]

/-- The proper length of one step of a generated path: `√(τ² - stepNormSq)`. -/
theorem sqrt_lorentzQ_generated_step (τ : ℝ) (w : SeamStepWorldline) (k : ℕ) :
    Real.sqrt (lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k)) =
      Real.sqrt (τ ^ 2 - stepNormSq (w.steps k)) := by
  rw [lorentzQ_generated_step]

/-- **Proper length of a generated path by counts.**  At a nonnegative
unit, the proper length over `M` steps is the crossing count times
`√(τ² - 4)` plus the rest count times `τ`. -/
theorem properLength_generated_counts (τ : ℝ) (hτ : 0 ≤ τ) (w : SeamStepWorldline)
    (M : ℕ) :
    properLength M (generatedPath τ w) =
      (crossingCount w.steps M : ℝ) * Real.sqrt (τ ^ 2 - 4) +
        (restCount w.steps M : ℝ) * τ := by
  induction M with
  | zero => simp [properLength, crossingCount, restCount]
  | succ M ih =>
    unfold properLength at ih ⊢
    rw [Finset.sum_range_succ, ih, sqrt_lorentzQ_generated_step]
    simp only [crossingCount, restCount]
    by_cases h : w.steps M = .rest
    · rw [h]
      simp only [stepNormSq, sub_zero, Real.sqrt_sq hτ, ne_eq, not_true_eq_false, if_false,
        if_true]
      push_cast
      ring
    · rw [stepNormSq_of_ne_rest _ h]
      simp only [ne_eq, h, not_false_eq_true, if_true, if_false]
      push_cast
      ring

/-- **The exact difference of the two rules.**  Index accrual minus
proper-length accrual on a window with `c` crossings is
`E * c * (τ - √(τ² - 4))`: rests contribute `E τ` under both rules,
crossings contribute `E τ` under index accrual and `E √(τ² - 4)` under
proper-length accrual.  Index accrual is at the declared unit `τ` per
index, so the comparison is between `E * M * τ` and `E * properLength`. -/
theorem accrual_difference (τ : ℝ) (hτ : 0 ≤ τ) (E : ℝ) (w : SeamStepWorldline) (M : ℕ) :
    indexAccrual τ E M (generatedPath τ w) - lengthAccrual E M (generatedPath τ w) =
      E * ((crossingCount w.steps M : ℝ) * (τ - Real.sqrt (τ ^ 2 - 4))) := by
  rw [indexAccrual_eq, lengthAccrual_eq, properLength_generated_counts τ hτ]
  have hM : (M : ℝ) = (crossingCount w.steps M : ℝ) + (restCount w.steps M : ℝ) := by
    have h := crossingCount_add_restCount w.steps M
    exact_mod_cast h.symm
  rw [hM]
  ring

theorem sqrt_sub_four_lt (τ : ℝ) (hτ : 2 < τ) : Real.sqrt (τ ^ 2 - 4) < τ := by
  rw [Real.sqrt_lt' (by linarith)]
  linarith

theorem sqrt_sub_four_pos (τ : ℝ) (hτ : 2 < τ) : 0 < Real.sqrt (τ ^ 2 - 4) := by
  rw [Real.sqrt_pos]
  nlinarith

/-- **Nonzero on every window with a crossing.**  For `E ≠ 0` and every
timelike unit `2 < τ`, the difference is nonzero as soon as the window
contains one crossing. -/
theorem accrual_difference_ne_zero (τ : ℝ) (hτ : 2 < τ) (E : ℝ) (hE : E ≠ 0)
    (w : SeamStepWorldline) (M : ℕ) (hc : crossingCount w.steps M ≠ 0) :
    indexAccrual τ E M (generatedPath τ w) - lengthAccrual E M (generatedPath τ w) ≠ 0 := by
  rw [accrual_difference τ (by linarith) E w M]
  have h1 : (0 : ℝ) < τ - Real.sqrt (τ ^ 2 - 4) := by linarith [sqrt_sub_four_lt τ hτ]
  have h2 : (0 : ℝ) < (crossingCount w.steps M : ℝ) := by
    exact_mod_cast Nat.pos_of_ne_zero hc
  exact mul_ne_zero hE (mul_ne_zero h2.ne' h1.ne')

/-- **Zero exactly on resting windows.** -/
theorem accrual_difference_eq_zero_iff (τ : ℝ) (hτ : 2 < τ) (E : ℝ) (hE : E ≠ 0)
    (w : SeamStepWorldline) (M : ℕ) :
    indexAccrual τ E M (generatedPath τ w) - lengthAccrual E M (generatedPath τ w) = 0 ↔
      crossingCount w.steps M = 0 := by
  constructor
  · intro h
    by_contra hc
    exact accrual_difference_ne_zero τ hτ E hE w M hc h
  · intro hc
    rw [accrual_difference τ (by linarith) E w M, hc]
    simp

/-- The two rules agree on a window exactly when every step of the window
is a rest. -/
theorem accruals_agree_iff_resting (τ : ℝ) (hτ : 2 < τ) (E : ℝ) (hE : E ≠ 0)
    (w : SeamStepWorldline) (M : ℕ) :
    indexAccrual τ E M (generatedPath τ w) = lengthAccrual E M (generatedPath τ w) ↔
      ∀ k, k < M → w.steps k = .rest := by
  rw [← sub_eq_zero, accrual_difference_eq_zero_iff τ hτ E hE, crossingCount_eq_zero_iff]

/-! ## (2) Refinement: action on the image, departure from the port class -/

/-- Every forward increment of the refined generated path is one half of
the parent increment at index `n / 2`. -/
theorem refine_generated_increment (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    refine (generatedPath τ w) (n + 1) - refine (generatedPath τ w) n =
      (1 / 2 : ℝ) • (generatedPath τ w (n / 2 + 1) - generatedPath τ w (n / 2)) := by
  rcases Nat.even_or_odd' n with ⟨k, rfl | rfl⟩
  · have hk : 2 * k / 2 = k := by omega
    rw [hk]
    exact refine_increment_first (generatedPath τ w) k
  · have hk : (2 * k + 1) / 2 = k := by omega
    rw [hk]
    exact refine_increment_second (generatedPath τ w) k

/-- **Lorentz squares of the refined path.**  Each half increment carries
one quarter of the parent Lorentz square: `(τ² - 4) / 4` at a crossing,
`τ² / 4` at a rest. -/
theorem refine_generated_lorentzQ (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    lorentzQ (refine (generatedPath τ w) (n + 1) - refine (generatedPath τ w) n) =
      (τ ^ 2 - stepNormSq (w.steps (n / 2))) / 4 := by
  rw [refine_generated_increment, OPH.C1Lorentz.lorentzQ_smul, lorentzQ_generated_step]
  ring

/-- The spatial half increment of the refined path has norm squared one
quarter of the step norm squared: `1` at a crossing, `0` at a rest. -/
theorem refine_generated_spatialNormSq (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    spatialNormSq (refine (generatedPath τ w) (n + 1) - refine (generatedPath τ w) n).2 =
      stepNormSq (w.steps (n / 2)) / 4 := by
  rw [refine_generated_increment, Prod.smul_snd, Prod.snd_sub, generatedPath_increment,
    OPH.C1Lorentz.spatialNormSq_smul, spatialNormSq_stepVector]
  ring

/-- **No seam step has spatial norm squared one.**  A step is a rest
(norm squared `0`) or a crossing (norm squared `4`). -/
theorem no_seam_step_norm_one (s : SeamStep) : stepNormSq s ≠ 1 := by
  cases s with
  | rest => simp [stepNormSq]
  | forward _ => simp [stepNormSq, seamNormSq]
  | backward _ => simp [stepNormSq, seamNormSq]

/-- **Half seams are not seams.**  The midpoint refinement of the generated
path of a worldline with a crossing at step `k` is the generated path of
no seam-step worldline at any unit: the refined spatial increment at index
`2k` has norm squared `1`, and no seam step does.  The declared refinement
acts on the Lorentz-module image and leaves the port class. -/
theorem refine_not_generated (τ : ℝ) (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest) (τ' : ℝ) (w' : SeamStepWorldline) :
    generatedPath τ' w' ≠ refine (generatedPath τ w) := by
  intro heq
  have h1 : spatialNormSq (generatedPath τ' w' (2 * k + 1) - generatedPath τ' w' (2 * k)).2 =
      stepNormSq (w'.steps (2 * k)) := by
    rw [Prod.snd_sub, generatedPath_increment, spatialNormSq_stepVector]
  have h2 : spatialNormSq (refine (generatedPath τ w) (2 * k + 1) -
      refine (generatedPath τ w) (2 * k)).2 = 1 := by
    rw [refine_generated_spatialNormSq]
    have hk2 : 2 * k / 2 = k := by omega
    rw [hk2, stepNormSq_of_ne_rest _ hk]
    norm_num
  rw [heq, h2] at h1
  exact no_seam_step_norm_one _ h1.symm

/-- **Index accrual fails refinement invariance.**  It is the member
`a = 0, b = τ` of the declared family, and `refinementInvariant_iff_b_zero`
gives `b = 0` for every invariant member at nonzero `E`. -/
theorem indexAccrual_not_refinementInvariant (τ E : ℝ) (hτ : τ ≠ 0) (hE : E ≠ 0) :
    ¬ RefinementInvariant 0 τ E := by
  rw [refinementInvariant_iff_b_zero 0 τ E hE]
  exact hτ

/-- The unit-`1` case is `additive_term_not_refinementInvariant` of the
imported module, reused here. -/
theorem indexAccrual_not_refinementInvariant_one (E : ℝ) (hE : E ≠ 0) :
    ¬ RefinementInvariant 0 1 E :=
  additive_term_not_refinementInvariant E hE

/-- **Proper-length accrual passes refinement invariance.**  It is the
member `a = 1, b = 0`. -/
theorem lengthAccrual_refinementInvariant (E : ℝ) : RefinementInvariant 1 0 E :=
  properLength_term_refinementInvariant E

/-- On the image, one refinement doubles the index accrual. -/
theorem indexAccrual_refine (τ E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    indexAccrual τ E (2 * M) (refine x) = 2 * indexAccrual τ E M x := by
  rw [indexAccrual_eq, indexAccrual_eq]
  push_cast
  ring

/-- On the image, one refinement fixes the proper-length accrual. -/
theorem lengthAccrual_refine (E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    lengthAccrual E (2 * M) (refine x) = lengthAccrual E M x := by
  rw [lengthAccrual_eq, lengthAccrual_eq, properLength_refine]

/-- **Compatibility of index accrual with refinement invariance.**  Index
accrual at unit `τ` fails invariance and proper-length accrual passes it;
on the generated path of a seam-step worldline at a timelike unit, index
accrual (a non-invariant member) takes the value of the invariant member
exactly on windows with no crossing.  Compatibility on a moving worldline
requires the internal process to advance per unit proper length, the
declared proper-time principle. -/
theorem index_accrual_compatible_iff_resting (τ : ℝ) (hτ : 2 < τ) (E : ℝ) (hE : E ≠ 0)
    (w : SeamStepWorldline) (M : ℕ) :
    ¬ RefinementInvariant 0 τ E ∧ RefinementInvariant 1 0 E ∧
      (indexAccrual τ E M (generatedPath τ w) = lengthAccrual E M (generatedPath τ w) ↔
        ∀ k, k < M → w.steps k = .rest) :=
  ⟨indexAccrual_not_refinementInvariant τ E (by linarith) hE, lengthAccrual_refinementInvariant E,
    accruals_agree_iff_resting τ hτ E hE w M⟩

/-- The resting seam-step worldline at a port: every step is a rest. -/
def restWorldline (u : Fin 12) : SeamStepWorldline where
  start := u
  steps := fun _ => .rest
  adm := fun _ => trivial

theorem restWorldline_crossingCount (u : Fin 12) (M : ℕ) :
    crossingCount (restWorldline u).steps M = 0 :=
  (crossingCount_eq_zero_iff _ M).mpr fun _ _ => rfl

theorem restWorldline_restCount (u : Fin 12) (M : ℕ) :
    restCount (restWorldline u).steps M = M := by
  have h := crossingCount_add_restCount (restWorldline u).steps M
  rw [restWorldline_crossingCount] at h
  omega

/-- The proper length of the resting worldline over `M` steps is `M * τ`:
the index count times the declared unit. -/
theorem properLength_restWorldline (τ : ℝ) (hτ : 0 ≤ τ) (u : Fin 12) (M : ℕ) :
    properLength M (generatedPath τ (restWorldline u)) = M * τ := by
  rw [properLength_generated_counts τ hτ, restWorldline_crossingCount,
    restWorldline_restCount]
  simp

theorem crossingWorldline_crossingCount_one (e : Fin 30) :
    crossingCount (crossingWorldline e).steps 1 = 1 := by
  simp [crossingCount, crossingWorldline]

/-- **The two rules are distinct functionals on the class.**  Proved: index
accrual at a fixed unit is one and the same function of the window on
every worldline (its definition ignores the path, so this is `rfl` after
unfolding), while proper-length accrual separates the one-step crossing
worldline from a resting worldline at every timelike unit and nonzero
`E`.  Scope (prose, proved nowhere): both rules are declared, the
selection between them is the declared principle of
`ProperTimeInternalAction`, and this module attributes no selection to
the source chain or the join. -/
theorem rules_both_declared (E : ℝ) (hE : E ≠ 0) :
    (∀ (τ : ℝ) (w w' : SeamStepWorldline) (M : ℕ),
        indexAccrual τ E M (generatedPath τ w) = indexAccrual τ E M (generatedPath τ w')) ∧
      (∀ (τ : ℝ), 2 < τ → ∀ (e : Fin 30) (u : Fin 12),
        lengthAccrual E 1 (generatedPath τ (crossingWorldline e)) ≠
          lengthAccrual E 1 (generatedPath τ (restWorldline u))) := by
  refine ⟨fun τ w w' M => by rw [indexAccrual_eq, indexAccrual_eq], ?_⟩
  intro τ hτ e u heq
  rw [lengthAccrual_eq, lengthAccrual_eq, properLength_restWorldline τ (by linarith),
    properLength_generated_counts τ (by linarith), crossingWorldline_crossingCount_one] at heq
  have hr : restCount (crossingWorldline e).steps 1 = 0 := by
    have h := crossingCount_add_restCount (crossingWorldline e).steps 1
    rw [crossingWorldline_crossingCount_one] at h
    omega
  rw [hr] at heq
  push_cast at heq
  have h1 := sqrt_sub_four_lt τ hτ
  have h2 : E * (1 * Real.sqrt (τ ^ 2 - 4) + 0 * τ) - E * (1 * τ) = 0 := by rw [heq]; ring
  have h3 : E * (Real.sqrt (τ ^ 2 - 4) - τ) = 0 := by linarith [h2]
  rcases mul_eq_zero.mp h3 with h | h
  · exact hE h
  · linarith

/-! ## (3) The uniform worldline and its dilation factor -/

/-- The number of indices `j < n` divisible by `k + 1`. -/
def crossingsBefore (k : ℕ) : ℕ → ℕ
  | 0 => 0
  | n + 1 => crossingsBefore k n + if n % (k + 1) = 0 then 1 else 0

/-- The declared uniform step sequence on seam `e` with `k` rests between
crossings: at every index divisible by `k + 1` a crossing, forward on
even crossings and backward on odd ones, a rest otherwise. -/
def uniformSteps (e : Fin 30) (k : ℕ) (n : ℕ) : SeamStep :=
  if n % (k + 1) = 0 then
    (if crossingsBefore k n % 2 = 0 then .forward e else .backward e)
  else .rest

/-- The port the uniform sequence occupies at index `n`: the left endpoint
after an even number of crossings, the right endpoint after an odd one. -/
def uniformPort (e : Fin 30) (k : ℕ) (n : ℕ) : Fin 12 :=
  if crossingsBefore k n % 2 = 0 then seamLeft e else seamRight e

theorem crossingsBefore_succ (k n : ℕ) :
    crossingsBefore k (n + 1) = crossingsBefore k n + if n % (k + 1) = 0 then 1 else 0 := rfl

theorem uniform_portSeq (e : Fin 30) (k : ℕ) (n : ℕ) :
    portSeq (seamLeft e) (uniformSteps e k) n = uniformPort e k n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show stepTarget (portSeq (seamLeft e) (uniformSteps e k) n) (uniformSteps e k n) =
      uniformPort e k (n + 1)
    rw [ih]
    unfold uniformSteps uniformPort
    rw [crossingsBefore_succ]
    by_cases h1 : n % (k + 1) = 0
    · by_cases h2 : crossingsBefore k n % 2 = 0
      · have h3 : (crossingsBefore k n + 1) % 2 ≠ 0 := by omega
        simp [h1, h2, h3, stepTarget]
      · have h3 : (crossingsBefore k n + 1) % 2 = 0 := by omega
        simp [h1, h2, h3, stepTarget]
    · simp [h1, stepTarget]

/-- The uniform worldline: start at the left endpoint of `e`, uniform
steps.  Admissibility follows from the port identity. -/
def uniformWorldline (e : Fin 30) (k : ℕ) : SeamStepWorldline where
  start := seamLeft e
  steps := uniformSteps e k
  adm := by
    intro n
    rw [uniform_portSeq]
    unfold StepAdmissible uniformSteps uniformPort
    by_cases h1 : n % (k + 1) = 0
    · by_cases h2 : crossingsBefore k n % 2 = 0
      · simp [h1, h2]
      · simp [h1, h2]
    · simp [h1]

theorem uniformSteps_ne_rest_iff (e : Fin 30) (k : ℕ) (n : ℕ) :
    uniformSteps e k n ≠ .rest ↔ n % (k + 1) = 0 := by
  unfold uniformSteps
  by_cases h1 : n % (k + 1) = 0
  · by_cases h2 : crossingsBefore k n % 2 = 0 <;> simp [h1, h2]
  · simp [h1]

/-- Within one period starting at a multiple of `k + 1`, the crossing
count over the first `j + 1 ≤ k + 1` steps is one. -/
theorem crossingCount_uniform_block (e : Fin 30) (k p : ℕ) :
    ∀ j, j ≤ k →
      crossingCount (fun i => uniformSteps e k ((k + 1) * p + i)) (j + 1) = 1 := by
  intro j
  induction j with
  | zero =>
    intro _
    simp only [crossingCount, ne_eq, uniformSteps_ne_rest_iff, Nat.add_zero,
      Nat.mul_mod_right, if_true]
  | succ j ih =>
    intro hj
    simp only [crossingCount] at ih ⊢
    rw [ih (by omega)]
    have h : ¬ ((k + 1) * p + (j + 1)) % (k + 1) = 0 := by
      rw [Nat.mul_add_mod, Nat.mod_eq_of_lt (by omega)]
      omega
    simp only [ne_eq, uniformSteps_ne_rest_iff, h, if_false]

/-- **Crossing count of the uniform worldline**: `p` crossings over `p`
periods of `k + 1` steps. -/
theorem crossingCount_uniform (e : Fin 30) (k : ℕ) (p : ℕ) :
    crossingCount (uniformWorldline e k).steps ((k + 1) * p) = p := by
  induction p with
  | zero => simp [crossingCount]
  | succ p ih =>
    rw [Nat.mul_succ, crossingCount_block_additive, ih]
    show p + crossingCount (fun i => uniformSteps e k ((k + 1) * p + i)) (k + 1) = p + 1
    rw [crossingCount_uniform_block e k p k le_rfl]

theorem restCount_uniform (e : Fin 30) (k : ℕ) (p : ℕ) :
    restCount (uniformWorldline e k).steps ((k + 1) * p) = k * p := by
  have h := crossingCount_add_restCount (uniformWorldline e k).steps ((k + 1) * p)
  rw [crossingCount_uniform] at h
  have h2 : (k + 1) * p = k * p + p := by ring
  rw [h2] at h ⊢
  omega

/-- **Proper length of the uniform worldline** over `p` periods:
`p * (k * τ + √(τ² - 4))`. -/
theorem properLength_uniform (τ : ℝ) (hτ : 0 ≤ τ) (e : Fin 30) (k p : ℕ) :
    properLength ((k + 1) * p) (generatedPath τ (uniformWorldline e k)) =
      (p : ℝ) * ((k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4)) := by
  rw [properLength_generated_counts τ hτ, crossingCount_uniform, restCount_uniform]
  push_cast
  ring

/-- The dilation factor of the uniform worldline: the proper length per
period over the proper length of a resting worldline per period, expressed
through the seam norm squared `4` and the declared unit `τ`. -/
def dilationFactor (τ : ℝ) (k : ℕ) : ℝ :=
  ((k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4)) / (((k : ℝ) + 1) * τ)

/-- **Exact dilation.**  Over `p` periods the proper length of the uniform
worldline is the dilation factor times the proper length `(k + 1) p τ` of
the resting worldline over the same index count. -/
theorem properLength_uniform_dilation (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k p : ℕ)
    (u : Fin 12) :
    properLength ((k + 1) * p) (generatedPath τ (uniformWorldline e k)) =
      dilationFactor τ k * properLength ((k + 1) * p) (generatedPath τ (restWorldline u)) := by
  rw [properLength_uniform τ hτ.le, properLength_restWorldline τ hτ.le]
  unfold dilationFactor
  have hk : ((k : ℝ) + 1) * τ ≠ 0 := by positivity
  field_simp
  push_cast
  ring

/-- The proper length per index step on the uniform worldline is
`dilationFactor τ k * τ`. -/
theorem properLength_uniform_per_index (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k p : ℕ)
    (hp : p ≠ 0) :
    properLength ((k + 1) * p) (generatedPath τ (uniformWorldline e k)) /
        (((k + 1) * p : ℕ) : ℝ) = dilationFactor τ k * τ := by
  rw [properLength_uniform τ hτ.le]
  unfold dilationFactor
  have hk : ((k : ℝ) + 1) * τ ≠ 0 := by positivity
  have hp' : (p : ℝ) ≠ 0 := by exact_mod_cast hp
  have hkp : (((k + 1) * p : ℕ) : ℝ) ≠ 0 := by positivity
  push_cast at hkp ⊢
  field_simp

/-- **Strictly below one** at every timelike unit. -/
theorem dilationFactor_lt_one (τ : ℝ) (hτ : 2 < τ) (k : ℕ) : dilationFactor τ k < 1 := by
  unfold dilationFactor
  rw [div_lt_one (by positivity)]
  linarith [sqrt_sub_four_lt τ hτ]

theorem dilationFactor_pos (τ : ℝ) (hτ : 2 < τ) (k : ℕ) : 0 < dilationFactor τ k := by
  unfold dilationFactor
  have := sqrt_sub_four_pos τ hτ
  positivity

/-- The dilation factor as one minus a per-period deficit. -/
theorem dilationFactor_eq (τ : ℝ) (hτ : 0 < τ) (k : ℕ) :
    dilationFactor τ k = 1 - (τ - Real.sqrt (τ ^ 2 - 4)) / τ * (1 / ((k : ℝ) + 1)) := by
  unfold dilationFactor
  have hk : ((k : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

/-- **The factor tends to one** as the number of rests per crossing grows. -/
theorem dilationFactor_tendsto_one (τ : ℝ) (hτ : 0 < τ) :
    Filter.Tendsto (fun k : ℕ => dilationFactor τ k) Filter.atTop (nhds 1) := by
  have hfun : (fun k : ℕ => dilationFactor τ k) =
      fun k : ℕ => 1 - (τ - Real.sqrt (τ ^ 2 - 4)) / τ * (1 / ((k : ℝ) + 1)) :=
    funext fun k => dilationFactor_eq τ hτ k
  rw [hfun]
  have h := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul
    ((τ - Real.sqrt (τ ^ 2 - 4)) / τ)
  have h2 := tendsto_const_nhds (x := (1 : ℝ)) |>.sub h
  simpa using h2

/-- **At `τ = 3, k = 0`** the factor is `√5 / 3`. -/
theorem dilationFactor_three_zero : dilationFactor 3 0 = Real.sqrt 5 / 3 := by
  unfold dilationFactor
  norm_num

/-! ## (4) Mean return of the recurrence clock under the two rules

Row statement.  Everything below is on the source clock and duration row:
the mean return of `SourceRecurrenceClock` lives on the abstract step
index, and each rule declares how that index is read along a worldline.
The index rule reads it as the shared step index of the join; the
proper-length rule reads one unit of the internal index as one rest step
`τ` of proper length.  No duration is attached to either reading. -/

/-- The index-rule mean return along a worldline (declared reading): the
mean return of the committed stable class on the abstract step index, read
as the shared step index.  The definition ignores its worldline argument;
no theorem relates `meanReturn` to a worldline. -/
def indexRuleReturn (_w : SeamStepWorldline) : ℝ := meanReturn 0

/-- **Index rule, declared reading.**  The value is `61511/7155`
(`meanReturn_values`), and "on every worldline" states that the definition
ignores the worldline (second conjunct is `rfl`). -/
theorem indexRule_meanReturn :
    (∀ w : SeamStepWorldline, indexRuleReturn w = 61511 / 7155) ∧
      ∀ w w' : SeamStepWorldline, indexRuleReturn w = indexRuleReturn w' :=
  ⟨fun _ => meanReturn_values.1, fun _ _ => rfl⟩

/-- The proper-length-rule mean return on the uniform worldline, counted in
shared index steps (declared reading): the internal process advances one
unit per rest step `τ` of proper length, the uniform worldline accrues
`dilationFactor τ k * τ` of proper length per index step
(`properLength_uniform_per_index`), so `meanReturn 0` internal units are
read as `meanReturn 0 / dilationFactor τ k` index steps.  The value is a
real number, not a step count; the definition encodes the reading. -/
def lengthRuleReturnIndexSteps (τ : ℝ) (k : ℕ) : ℝ := meanReturn 0 / dilationFactor τ k

/-- **Proper-length rule on the uniform worldline, declared reading.**  The
index count is `(61511/7155) * ((k + 1) τ) / (k τ + √(τ² - 4))`, the
consequence of the declared reading, of `meanReturn_values`, and of the
algebra of `dilationFactor`. -/
theorem lengthRule_return_indexSteps (τ : ℝ) (hτ : 2 < τ) (k : ℕ) :
    lengthRuleReturnIndexSteps τ k =
      (61511 / 7155 : ℝ) * (((k : ℝ) + 1) * τ) / ((k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4)) := by
  unfold lengthRuleReturnIndexSteps dilationFactor
  rw [meanReturn_values.1]
  have h1 : ((k : ℝ) + 1) * τ ≠ 0 := by
    have : (0 : ℝ) < τ := by linarith
    positivity
  have h2 : (k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4) ≠ 0 := by
    have := sqrt_sub_four_pos τ hτ
    have : (0 : ℝ) ≤ (k : ℝ) * τ := by
      have : (0 : ℝ) < τ := by linarith
      positivity
    linarith
  field_simp

/-- The proper-length-rule value is the inverse dilation factor times the
index-rule value. -/
theorem lengthRule_return_eq_inv_dilation (τ : ℝ) (k : ℕ) (w : SeamStepWorldline) :
    lengthRuleReturnIndexSteps τ k = (dilationFactor τ k)⁻¹ * indexRuleReturn w := by
  unfold lengthRuleReturnIndexSteps indexRuleReturn
  rw [div_eq_inv_mul]

/-- **More shared steps on a moving worldline.**  At every timelike unit
and every `k`, the proper-length rule returns after strictly more shared
index steps than the index rule. -/
theorem lengthRule_return_gt_index (τ : ℝ) (hτ : 2 < τ) (k : ℕ) (w : SeamStepWorldline) :
    indexRuleReturn w < lengthRuleReturnIndexSteps τ k := by
  unfold lengthRuleReturnIndexSteps indexRuleReturn
  rw [lt_div_iff₀ (dilationFactor_pos τ hτ k)]
  have h0 := OPH.QFT.SourceRecurrenceClock.meanReturn_pos 0
  have h1 := dilationFactor_lt_one τ hτ k
  nlinarith

/-- **Consistency with the proper length of the uniform worldline.**  The
proper length accrued per index step, `dilationFactor τ k * τ`
(`properLength_uniform_per_index`), times the proper-length-rule index
count, is `meanReturn 0 * τ`: the index-rule value read in rest steps of
proper length. -/
theorem lengthRule_return_consistent (τ : ℝ) (hτ : 2 < τ) (k : ℕ) :
    dilationFactor τ k * τ * lengthRuleReturnIndexSteps τ k = meanReturn 0 * τ := by
  unfold lengthRuleReturnIndexSteps
  have h := (dilationFactor_pos τ hτ k).ne'
  field_simp

/-- On the resting worldline the two rules coincide: the proper length per
index step is `τ`, so `meanReturn 0` rest steps of proper length are
`meanReturn 0` shared steps.  This is the `k`-free reading of
`accruals_agree_iff_resting` at the recurrence clock. -/
theorem rules_agree_on_rest (τ : ℝ) (hτ : 0 ≤ τ) (u : Fin 12) (M : ℕ) :
    properLength M (generatedPath τ (restWorldline u)) = (M : ℝ) * τ :=
  properLength_restWorldline τ hτ u M

/-! ## Negatives cited and rows touched -/

/-- The Legendre non-identifiability at its scope, re-cited from
`Geometry/InternalEnergyInertia.lean`: realized histories select no
velocity curvature or Legendre map, so both accrual rules are declared
enrichments and the selection between them is by a declared principle. -/
theorem legendre_scope_cited :
    (¬ ∃ vel : ℝ → ℝ → ℝ, OPH.Variational.SolvesMomentum
        OPH.Variational.chainLogLagrangian vel) ∧
      OPH.Variational.chainCurvedLagrangian 1 ≠ OPH.Variational.chainCurvedLagrangian 2 :=
  OPH.InternalEnergyInertia.legendre_nonidentifiability_cited

/-- The rate non-identifiability at its scope, re-cited from
`ObserverPatchHolography/RateNonidentifiability.lean`: the declared
reading of an abstract transition system (locked set, reachability,
basin) returns no first-locking count.  It concerns that reading only and
forbids nothing about a process hosted on the source chain; the index
rule above is a declaration on the join and uses no functional of that
reading. -/
theorem rate_scope_cited {S : Type} [DecidableEq S] (T : S → S) (x : S) (n : ℕ)
    (hlock : OPH.RateNonidentifiability.FirstLock T n x) (hpos : 0 < n)
    (rate : OPH.RateNonidentifiability.DeclaredReading S → ℕ)
    (hrate : ∀ s : ℕ, OPH.RateNonidentifiability.FirstLock
      (OPH.RateNonidentifiability.stutterStep T s)
      (rate (OPH.RateNonidentifiability.projReading T s)) (x, Fin.last s)) :
    False :=
  OPH.RateNonidentifiability.no_rate_functional T x n hlock hpos rate hrate

/-- The rows this module touches, as register labels: the source clock
and duration row (both rules and the unit `τ` declared), the coupled
action (the internal-action family is a declared enrichment), and the
physical spacetime attachment row (the Lorentz module is the declared
image of the port class).  The light-signal row, the laboratory clock and
energy calibration import, and the gravitation-route energy identification
are named in the header and carry no `OpenRow` label.  A label is not a
discharge. -/
def touchedRows : List OpenRow :=
  [OpenRow.sourceClock, OpenRow.coupledAction, OpenRow.spacetimeAttachment]

/-- The rows this module discharges: none. -/
def dischargedRows : List OpenRow := []

theorem dischargedRows_empty : dischargedRows = [] := rfl

end

end OPH.SourceClockRateAlongWorldlines

/- Axiom audit: expected at most `propext`, `Classical.choice`, `Quot.sound`
per line; no native decision procedure. -/

#print axioms OPH.SourceClockRateAlongWorldlines.properLength_generated_counts
#print axioms OPH.SourceClockRateAlongWorldlines.accrual_difference
#print axioms OPH.SourceClockRateAlongWorldlines.accrual_difference_ne_zero
#print axioms OPH.SourceClockRateAlongWorldlines.accrual_difference_eq_zero_iff
#print axioms OPH.SourceClockRateAlongWorldlines.accruals_agree_iff_resting
#print axioms OPH.SourceClockRateAlongWorldlines.refine_generated_lorentzQ
#print axioms OPH.SourceClockRateAlongWorldlines.refine_generated_spatialNormSq
#print axioms OPH.SourceClockRateAlongWorldlines.no_seam_step_norm_one
#print axioms OPH.SourceClockRateAlongWorldlines.refine_not_generated
#print axioms OPH.SourceClockRateAlongWorldlines.indexAccrual_not_refinementInvariant
#print axioms OPH.SourceClockRateAlongWorldlines.indexAccrual_not_refinementInvariant_one
#print axioms OPH.SourceClockRateAlongWorldlines.lengthAccrual_refinementInvariant
#print axioms OPH.SourceClockRateAlongWorldlines.indexAccrual_refine
#print axioms OPH.SourceClockRateAlongWorldlines.lengthAccrual_refine
#print axioms OPH.SourceClockRateAlongWorldlines.index_accrual_compatible_iff_resting
#print axioms OPH.SourceClockRateAlongWorldlines.rules_both_declared
#print axioms OPH.SourceClockRateAlongWorldlines.uniform_portSeq
#print axioms OPH.SourceClockRateAlongWorldlines.crossingCount_uniform
#print axioms OPH.SourceClockRateAlongWorldlines.properLength_uniform
#print axioms OPH.SourceClockRateAlongWorldlines.properLength_uniform_dilation
#print axioms OPH.SourceClockRateAlongWorldlines.properLength_uniform_per_index
#print axioms OPH.SourceClockRateAlongWorldlines.dilationFactor_lt_one
#print axioms OPH.SourceClockRateAlongWorldlines.dilationFactor_pos
#print axioms OPH.SourceClockRateAlongWorldlines.dilationFactor_tendsto_one
#print axioms OPH.SourceClockRateAlongWorldlines.dilationFactor_three_zero
#print axioms OPH.SourceClockRateAlongWorldlines.indexRule_meanReturn
#print axioms OPH.SourceClockRateAlongWorldlines.lengthRule_return_indexSteps
#print axioms OPH.SourceClockRateAlongWorldlines.lengthRule_return_eq_inv_dilation
#print axioms OPH.SourceClockRateAlongWorldlines.lengthRule_return_gt_index
#print axioms OPH.SourceClockRateAlongWorldlines.lengthRule_return_consistent
#print axioms OPH.SourceClockRateAlongWorldlines.legendre_scope_cited
#print axioms OPH.SourceClockRateAlongWorldlines.rate_scope_cited
#print axioms OPH.SourceClockRateAlongWorldlines.dischargedRows_empty
