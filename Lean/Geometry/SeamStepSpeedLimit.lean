import CarrierModeOscillators
import Geometry.WorldlineHopTransport

set_option autoImplicit false

open scoped BigOperators

namespace OPH.SeamStepSpeedLimit

open OPH.ScaledMaxwellStability
open OPH.DiscreteCoulombGreen OPH.TemporalMaxwellEvolution
open OPH.CarrierModeOscillators
open OPH.WorldlineHopTransport
open OPH.SeamCurrentCarrierQuotient (seamLeft seamRight)
open OPH.CommonWorldInstrumentJoin
open OPH.CommonWorldMaxwellClockJoin
open OPH.C1Lorentz (Spatial Herm2 spatialNormSq lorentzQ)

/-!
# The seam-step speed limit from the field stability window

STATUS.  Candidate module on the light-signal, coupled-action, source-clock
and spacetime-attachment rows.  The committed field evolution has a sharp
stability window in the step `h` (`courant_threshold_sharp`, two-sided:
below the threshold `h² (3 + √5) < 4`, i.e. `h² < 2/φ²`, every
zero-current solution is bounded for every datum; above it an unbounded
zero-current solution exists; the boundary `h² (3 + √5) = 4` is covered
by neither direction), and the seam-step worldlines of the coupled-action row have a
sharp timelike threshold in the declared time unit `τ` per step
(`seam_step_timelike_iff`: a crossing step is timelike iff `τ² > 4`).  The
Maxwell-clock join reads the field step and the record step at one index
but carries the step `h` and the duration `stepDuration` as two separate
declared fields.  This module declares the identification `τ = h` of the
worldline unit with the field step (`UnitIdentifiedWithStep`), keeps the
alternative `τ ≠ h` explicit, and proves what follows on each side: under
the identification the two thresholds are incompatible with the exact gap
`2φ² = 3 + √5`, rest-diluted worldlines need at least two rests per
crossing in every stable evolution, the block speed of a timelike block is
below `1` with supremum `1` over the window, and without the identification
the two thresholds are independent.  Every statement is an inequality or
identity in `ℝ`; nothing here selects the identification.

WHAT IS PROVED.

1. Declared identification and exact incompatibility.  Under
   `UnitIdentifiedWithStep τ h` (the declared equation `τ = h`) a crossing
   step is timelike iff `4 < h²` (`crossing_timelike_iff_identified`); the
   window `h² (3 + √5) < 4` gives `h² < 1` (`window_sq_lt_one`); no `h`
   lies in the window with `4 < h²` (`window_timelike_incompatible`); in
   every stable evolution every crossing step of every seam-step worldline
   at the identified unit is spacelike, and the evolution is bounded for
   every datum (`stable_evolution_crossing_spacelike`); the ratio of the
   two thresholds is `4 / (2/φ²) = 2φ² = 3 + √5` (`threshold_gap`).  For a
   joined architecture whose duration equals its field step the same
   conclusion holds with no window hypothesis, because every joined
   architecture's step satisfies the sharp certificate
   (`joined_step_certified`; `joined_identified_crossing_spacelike`).
2. Rest-diluted worldlines.  A block of `k` rests followed by one crossing
   has time `(k+1) τ` and spatial displacement of norm squared `4`, Lorentz
   square `((k+1) τ)² - 4` (`block_lorentzQ`); it is timelike iff
   `2 < (k+1) τ` (`block_timelike_iff`).  Under the identification inside
   the window, a timelike block has `√2 φ < k + 1` (`block_rests_gt`),
   hence `2 ≤ k` (`block_two_rests`); `k = 2` is timelike iff `2/3 < h`
   (`two_rest_block_timelike_iff`), and `2/3 < √2/φ` (`two_thirds_lt_edge`),
   so the two-rest block is timelike for every `h` in `(2/3, √2/φ)`.  The
   two-rest worldline is exhibited (`twoRestWorldline`, alternating a
   forward and a backward crossing of seam `0` after two rests) and its
   blocks are timelike for such `h` (`twoRestWorldline_block_timelike`).
3. Speed reading.  The block speed is `2 / ((k+1) τ)` (`blockSpeed`); for
   `0 < τ` it is below `1` iff the block is timelike
   (`blockSpeed_lt_one_iff`) and equals `1` iff the block is null
   (`blockSpeed_eq_one_iff`).  Every timelike block of a stable evolution
   has speed below `1` (`stable_block_speed_lt_one`); for every `ε > 0`
   there are `h` in the window and `k` with a timelike block of speed above
   `1 - ε` (`block_speed_sup_one`); no timelike block has speed `1`
   (`timelike_block_speed_ne_one`).  So `1` is the supremum of transported
   timelike block speeds over the window and is never attained.
4. Without the identification.  For `0 < τ` a crossing step is timelike
   iff `2 < τ` (`crossing_timelike_iff_pos`); at `τ = 3`, `h = 1/2` the
   crossing is timelike and the evolution is bounded for every datum
   (`thresholds_independent`); for every `h` in the window every `τ` with
   `2 < τ` gives timelike crossings and every `τ` with `τ < 2` gives
   spacelike crossings (`window_does_not_constrain_unit`); the committed
   joined witness carries `stepDuration = 1` and `h = 1/2`, so it is an
   inhabitant outside the identification
   (`committedJoinedWitness_not_identified`); at the free unit `τ = 2`
   every crossing is null (`duration_two_null_crossing`), and the corpus
   inhabitant `doubleStepJoinedWitness` has that duration
   (`doubleStepJoinedWitness_null_crossing`).

PRIOR WORK.  `courant_threshold_sharp` and `instability_above_golden`
(`Lean/Screen/ScaledMaxwellStability.lean`) are the window; `window_iff_goldenRatio`,
`window_iff_lt_sqrt_two_div_goldenRatio`,
`three_add_sqrt5_eq_two_goldenRatio_sq` and
`three_sub_sqrt5_eq_two_div_goldenRatio_sq`
(`Lean/Screen/CarrierModeOscillators.lean`) are the golden-ratio forms;
`seam_step_timelike_iff`, `seam_step_timelike_of_two_lt`,
`lorentzQ_generated_step`, `clockAction_generated` and `two_units_two_paths`
(`Lean/Geometry/WorldlineHopTransport.lean`) are the per-step threshold and
the non-forcing of the unit; `scaledBundle_courant_ge_golden`
(`Lean/Geometry/CommonWorldInstrumentJoin.lean`) places every scaled bundle
inside the sharp window; `MaxwellClockJoinedArchitecture`,
`joined_step_certified` (the step of every joined architecture satisfies
`h² (3 + √5) < 4`, used directly here), `joined_stepDuration_not_forced`,
`committedJoinedWitness` and `doubleStepJoinedWitness`
(`Lean/Geometry/CommonWorldMaxwellClockJoin.lean`) carry the duration as a
separate field.  `sqrt_five_lt_three` restates
`OPH.OrientedFaceBracketSelector.sqrt_five_lt_three`
(`Lean/Screen/OrientedFaceBracketSelector.lean`), whose module is not
imported here.  The per-step threshold restated
here at `τ = h` (`crossing_timelike_iff_identified`) is a labelled
restatement of `seam_step_timelike_iff` for the identified unit.

ROWS TOUCHED.  The source clock and duration row (the identification
`τ = h` is declared here; the duration of the join and the unit of the
worldline are two declared reals, and no theorem selects their equality);
the light-signal row (the field step `h` is the step of the committed
evolution; no signal propagation is attached to the threshold `τ = 2` or
to the speed `1`); the coupled-action row (the seam-step class, the block
speed and the rest dilution are declared readings of the transported
worldline); the physical spacetime attachment row (seam units and step
units are carrier coordinates; no attachment to a spacetime point or to a
metre or second is supplied); the laboratory clock and energy calibration
import (no calibration of `h`, `τ`, or the speed is attached); the
gravitation-route energy identification (no energy is identified).  The
module discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized histories
select no velocity curvature or Legendre map, so the clock action, the
block speed and the identification are declared enrichments; cited at
scope only.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`; the generated path has time coordinate
`τ n` at step `n`; a seam vector has norm squared `4` and length `2`, one
seam unit of length `2`; the block speed is spatial length over time,
`2 / ((k+1) τ)`; `φ = Real.goldenRatio`; the window edge is
`√2/φ = √(3 - √5) ≈ 0.874`; `√2 φ = √(3 + √5) ≈ 2.288`.

FALSIFIER.  The module is wrong if some `h` satisfies both
`h² (3 + √5) < 4` and `4 < h²`, if the ratio of the thresholds differs
from `3 + √5`, if the two-rest worldline has an inadmissible step, or if
a timelike block of a stable evolution at the identified unit has speed
`≥ 1`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

noncomputable section

/-! ## (1) The declared identification and the exact incompatibility -/

/-- The sharp stability window of the committed field evolution in the
step `h`: `h² (3 + √5) < 4`.  The proposition is the antecedent of
`courant_threshold_sharp`. -/
def StableWindow (h : ℝ) : Prop := h ^ 2 * (3 + Real.sqrt 5) < 4

/-- **Declared identification.**  The worldline's time unit per step is
the field step: `τ = h`.  This is a declared modelling choice on the
source clock and duration row; the alternative `τ ≠ h` is kept explicit in
section (4). -/
def UnitIdentifiedWithStep (τ h : ℝ) : Prop := τ = h

theorem unitIdentifiedWithStep_refl (h : ℝ) : UnitIdentifiedWithStep h h := rfl

/-- A crossing step at the identified unit is timelike iff `4 < h²`
(labelled restatement of `seam_step_timelike_iff` at `τ = h`). -/
theorem crossing_timelike_iff_identified (τ h : ℝ) (hid : UnitIdentifiedWithStep τ h)
    (w : SeamStepWorldline) (k : ℕ) (hk : w.steps k ≠ .rest) :
    0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) ↔ 4 < h ^ 2 := by
  rw [seam_step_timelike_iff τ w k hk, evalPhi_threshold, hid]
  norm_num

theorem two_lt_sqrt_five : 2 < Real.sqrt 5 := by
  rw [show (2 : ℝ) = Real.sqrt 4 by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

/-- `√5 < 3`; restates `OPH.OrientedFaceBracketSelector.sqrt_five_lt_three`
from a module not imported here. -/
theorem sqrt_five_lt_three : Real.sqrt 5 < 3 := by
  rw [show (3 : ℝ) = Real.sqrt 9 by
    rw [show (9 : ℝ) = 3 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

/-- The window gives `h² < 1`: the field step is shorter than half the
null crossing time `2` of one seam. -/
theorem window_sq_lt_one (h : ℝ) (hw : StableWindow h) : h ^ 2 < 1 := by
  unfold StableWindow at hw
  have h5 := two_lt_sqrt_five
  nlinarith [sq_nonneg h]

/-- The window in golden-ratio form, `h² < 2/φ² < 1`. -/
theorem window_sq_lt_two_div_goldenRatio_sq_lt_one (h : ℝ) (hw : StableWindow h) :
    h ^ 2 < 2 / Real.goldenRatio ^ 2 ∧ 2 / Real.goldenRatio ^ 2 < 1 := by
  refine ⟨(window_iff_goldenRatio h).mp hw, ?_⟩
  rw [← three_sub_sqrt5_eq_two_div_goldenRatio_sq]
  linarith [two_lt_sqrt_five]

/-- **Exact incompatibility.**  No step `h` lies in the stability window
and has a timelike per-step crossing at the identified unit. -/
theorem window_timelike_incompatible (h : ℝ) : ¬ (StableWindow h ∧ 4 < h ^ 2) := by
  rintro ⟨hw, h4⟩
  linarith [window_sq_lt_one h hw]

/-- In the window every per-step crossing at the identified unit is
spacelike. -/
theorem window_crossing_spacelike (h : ℝ) (hw : StableWindow h)
    (w : SeamStepWorldline) (k : ℕ) (hk : w.steps k ≠ .rest) :
    lorentzQ (generatedPath h w (k + 1) - generatedPath h w k) < 0 := by
  rw [lorentzQ_generated_step, stepNormSq_of_ne_rest _ hk]
  linarith [window_sq_lt_one h hw]

/-- **Stable evolution admits no per-step crossing.**  For `h ≠ 0` in the
window: every zero-current solution is bounded for every datum
(`courant_threshold_sharp`), and every crossing step of every seam-step
worldline at the identified unit is spacelike.  Conversely, if some
crossing step at the identified unit is timelike then the window fails
and an unbounded zero-current solution exists. -/
theorem stable_evolution_crossing_spacelike (h : ℝ) (hh : h ≠ 0) (hw : StableWindow h) :
    (∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
        AmpereEvolutionScaled h A φ (fun _ ↦ 0) → ∀ n,
          realSeamEnergy (electricFieldScaled h A φ n) ≤
            8 * fieldEnergyScaled h A φ 0 / (4 - h ^ 2 * (3 + Real.sqrt 5)))
    ∧ (∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
        lorentzQ (generatedPath h w (k + 1) - generatedPath h w k) < 0) :=
  ⟨fun A φ hA n ↦ ((courant_threshold_sharp h hh).1 hw A φ hA n).2.1,
    fun w k hk ↦ window_crossing_spacelike h hw w k hk⟩

/-- A timelike per-step crossing at the identified unit puts `h` above
the window, where an unbounded zero-current solution exists. -/
theorem timelike_crossing_unstable (h : ℝ) (hh : h ≠ 0) (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest)
    (ht : 0 < lorentzQ (generatedPath h w (k + 1) - generatedPath h w k)) :
    ∃ A : ℕ → Fin 30 → ℝ, AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricFieldScaled h A (fun _ ↦ 0) n) := by
  have h4 : 4 < h ^ 2 :=
    (crossing_timelike_iff_identified h h rfl w k hk).mp ht
  refine (courant_threshold_sharp h hh).2 ?_
  have h5 := two_lt_sqrt_five
  nlinarith

/-- **The exact gap.**  The ratio of the timelike threshold `4` for `h²`
to the stability threshold `2/φ²` is `2φ² = 3 + √5`. -/
theorem threshold_gap :
    4 / (2 / Real.goldenRatio ^ 2) = 2 * Real.goldenRatio ^ 2 ∧
      2 * Real.goldenRatio ^ 2 = 3 + Real.sqrt 5 := by
  refine ⟨?_, three_add_sqrt5_eq_two_goldenRatio_sq.symm⟩
  have hg := goldenRatio_sq_pos
  field_simp
  ring

/-- **Joined architecture at the identified duration.**  If a joined
architecture's declared duration equals its field step, every crossing
step of every seam-step worldline at that duration is spacelike, with no
window hypothesis: the step of every joined architecture satisfies the
sharp certificate (`joined_step_certified`). -/
theorem joined_identified_crossing_spacelike (J : MaxwellClockJoinedArchitecture)
    (hid : UnitIdentifiedWithStep J.stepDuration J.scaled.h)
    (w : SeamStepWorldline) (k : ℕ) (hk : w.steps k ≠ .rest) :
    lorentzQ (generatedPath J.stepDuration w (k + 1) - generatedPath J.stepDuration w k) < 0 := by
  have hw : StableWindow J.scaled.h := joined_step_certified J
  rw [hid]
  exact window_crossing_spacelike _ hw w k hk

/-! ## (2) Rest-diluted worldlines -/

/-- A rest block of a seam-step worldline: `k` rest steps from index `n`,
then one crossing at index `n + k`; `k + 1` steps in all. -/
def RestBlock (w : SeamStepWorldline) (n k : ℕ) : Prop :=
  (∀ i, i < k → w.steps (n + i) = .rest) ∧ w.steps (n + k) ≠ .rest

theorem spatialAt_rest_run (w : SeamStepWorldline) (n : ℕ) :
    ∀ i, (∀ j, j < i → w.steps (n + j) = .rest) → spatialAt w (n + i) = spatialAt w n := by
  intro i
  induction i with
  | zero => intro _; rfl
  | succ i ih =>
    intro hr
    rw [← add_assoc, spatialAt_succ, ih (fun j hj ↦ hr j (Nat.lt_succ_of_lt hj)),
      hr i (Nat.lt_succ_self i), stepVector_rest, add_zero]

/-- The spatial displacement across a rest block is the step vector of the
crossing: one seam, norm squared `4`. -/
theorem block_displacement (w : SeamStepWorldline) (n k : ℕ) (hb : RestBlock w n k) :
    spatialAt w (n + k + 1) - spatialAt w n = stepVector (w.steps (n + k)) := by
  rw [spatialAt_succ, spatialAt_rest_run w n k hb.1, add_sub_cancel_left]

/-- The time coordinate advances by `(k+1) τ` across a block. -/
theorem block_time (τ : ℝ) (w : SeamStepWorldline) (n k : ℕ) :
    (generatedPath τ w (n + k + 1)).1 - (generatedPath τ w n).1 = (k + 1 : ℝ) * τ := by
  rw [generatedPath_fst, generatedPath_fst]
  push_cast
  ring

/-- **Block Lorentz square.**  Across a rest block the generated path has
time `(k+1) τ`, spatial norm squared `4`, and Lorentz square
`((k+1) τ)² - 4`. -/
theorem block_lorentzQ (τ : ℝ) (w : SeamStepWorldline) (n k : ℕ) (hb : RestBlock w n k) :
    lorentzQ (generatedPath τ w (n + k + 1) - generatedPath τ w n) =
      ((k + 1 : ℝ) * τ) ^ 2 - 4 := by
  unfold lorentzQ
  rw [Prod.fst_sub, Prod.snd_sub, generatedPath_fst, generatedPath_fst,
    generatedPath_snd, generatedPath_snd, block_displacement w n k hb,
    spatialNormSq_stepVector, stepNormSq_of_ne_rest _ hb.2]
  push_cast
  ring

/-- A block is timelike iff `2 < (k+1) τ`, for a positive unit. -/
theorem block_timelike_iff (τ : ℝ) (hτ : 0 < τ) (w : SeamStepWorldline) (n k : ℕ)
    (hb : RestBlock w n k) :
    0 < lorentzQ (generatedPath τ w (n + k + 1) - generatedPath τ w n) ↔
      2 < (k + 1 : ℝ) * τ := by
  rw [block_lorentzQ τ w n k hb]
  have hx : 0 < (k + 1 : ℝ) * τ := by positivity
  constructor
  · intro h
    nlinarith
  · intro h
    nlinarith

theorem sqrt_two_mul_goldenRatio_sq :
    (Real.sqrt 2 * Real.goldenRatio) ^ 2 = 3 + Real.sqrt 5 := by
  rw [mul_pow, Real.sq_sqrt (by norm_num), three_add_sqrt5_eq_two_goldenRatio_sq]

theorem sqrt_two_div_goldenRatio_sq :
    (Real.sqrt 2 / Real.goldenRatio) ^ 2 = 2 / Real.goldenRatio ^ 2 := by
  rw [div_pow, Real.sq_sqrt (by norm_num)]

/-- `2 < √2 φ` (numerically `2.288`). -/
theorem two_lt_sqrt_two_mul_goldenRatio : 2 < Real.sqrt 2 * Real.goldenRatio := by
  have hpos : 0 < Real.sqrt 2 * Real.goldenRatio := by positivity
  have hsq := sqrt_two_mul_goldenRatio_sq
  nlinarith [two_lt_sqrt_five]

/-- The window edge times `√2 φ` is `2`. -/
theorem edge_mul_sqrt_two_goldenRatio :
    (Real.sqrt 2 / Real.goldenRatio) * (Real.sqrt 2 * Real.goldenRatio) = 2 := by
  have hg := Real.goldenRatio_pos.ne'
  have h2 := Real.mul_self_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
  field_simp
  linarith

/-- **Rest count bound.**  Under the identification inside the window, a
timelike block has `√2 φ < k + 1`. -/
theorem block_rests_gt (h : ℝ) (h0 : 0 < h) (hw : StableWindow h)
    (w : SeamStepWorldline) (n k : ℕ) (hb : RestBlock w n k)
    (ht : 0 < lorentzQ (generatedPath h w (n + k + 1) - generatedPath h w n)) :
    Real.sqrt 2 * Real.goldenRatio < (k + 1 : ℝ) := by
  have h2 : 2 < (k + 1 : ℝ) * h := (block_timelike_iff h h0 w n k hb).mp ht
  have he : h < Real.sqrt 2 / Real.goldenRatio :=
    (window_iff_lt_sqrt_two_div_goldenRatio h h0).mp hw
  have hk1 : (0 : ℝ) < k + 1 := by positivity
  have hg : 0 < Real.sqrt 2 * Real.goldenRatio := by positivity
  have h3 : 2 < (k + 1 : ℝ) * (Real.sqrt 2 / Real.goldenRatio) :=
    lt_trans h2 (mul_lt_mul_of_pos_left he hk1)
  have h4 := mul_lt_mul_of_pos_right h3 hg
  rw [mul_assoc, edge_mul_sqrt_two_goldenRatio] at h4
  linarith

/-- **At least two rests per crossing** at a step in the stability window,
at the identified unit: a timelike block has `2 ≤ k`. -/
theorem block_two_rests (h : ℝ) (h0 : 0 < h) (hw : StableWindow h)
    (w : SeamStepWorldline) (n k : ℕ) (hb : RestBlock w n k)
    (ht : 0 < lorentzQ (generatedPath h w (n + k + 1) - generatedPath h w n)) :
    2 ≤ k := by
  have h1 := block_rests_gt h h0 hw w n k hb ht
  have h2 := two_lt_sqrt_two_mul_goldenRatio
  have h3 : (1 : ℝ) < k := by linarith
  have h4 : (1 : ℕ) < k := by exact_mod_cast h3
  omega

/-- The two-rest block is timelike iff `2/3 < h`. -/
theorem two_rest_block_timelike_iff (h : ℝ) (h0 : 0 < h) (w : SeamStepWorldline) (n : ℕ)
    (hb : RestBlock w n 2) :
    0 < lorentzQ (generatedPath h w (n + 2 + 1) - generatedPath h w n) ↔ 2 / 3 < h := by
  rw [block_timelike_iff h h0 w n 2 hb]
  norm_num
  constructor <;> intro hh <;> linarith

/-- `2/3 < √2/φ`: the two-rest threshold lies inside the window. -/
theorem two_thirds_lt_edge : (2 : ℝ) / 3 < Real.sqrt 2 / Real.goldenRatio := by
  refine lt_of_pow_lt_pow_left₀ 2 (by positivity) ?_
  rw [sqrt_two_div_goldenRatio_sq, ← three_sub_sqrt5_eq_two_div_goldenRatio_sq]
  have h5 : Real.sqrt 5 < 23 / 9 := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  norm_num
  linarith

/-- The two-rest block is timelike for every `h` in `(2/3, √2/φ)`, an
interval inside the window. -/
theorem two_rest_admissible_interval (h : ℝ) (hl : 2 / 3 < h)
    (hu : h < Real.sqrt 2 / Real.goldenRatio) (w : SeamStepWorldline) (n : ℕ)
    (hb : RestBlock w n 2) :
    StableWindow h ∧ 0 < lorentzQ (generatedPath h w (n + 2 + 1) - generatedPath h w n) := by
  have h0 : 0 < h := by linarith
  exact ⟨(window_iff_lt_sqrt_two_div_goldenRatio h h0).mpr hu,
    (two_rest_block_timelike_iff h h0 w n hb).mpr hl⟩

/-! ### The two-rest worldline -/

/-- Steps of the two-rest worldline: two rests, then a crossing of seam `0`,
forward on even blocks and backward on odd blocks. -/
def twoRestSteps (n : ℕ) : SeamStep :=
  if n % 3 = 2 then (if (n / 3) % 2 = 0 then .forward 0 else .backward 0) else .rest

/-- The port of the two-rest worldline at step `n`. -/
def twoRestPort (n : ℕ) : Fin 12 :=
  if (n / 3) % 2 = 0 then seamLeft 0 else seamRight 0

theorem twoRest_portSeq : ∀ n, portSeq (seamLeft 0) twoRestSteps n = twoRestPort n := by
  intro n
  induction n with
  | zero => rfl
  | succ n ih =>
    show stepTarget (portSeq (seamLeft 0) twoRestSteps n) (twoRestSteps n) = twoRestPort (n + 1)
    rw [ih]
    unfold twoRestSteps twoRestPort
    by_cases h3 : n % 3 = 2
    · have hd : (n + 1) / 3 = n / 3 + 1 := by omega
      rw [hd]
      by_cases hp : n / 3 % 2 = 0
      · have hp' : (n / 3 + 1) % 2 ≠ 0 := by omega
        simp only [if_pos h3, if_pos hp, if_neg hp']
        rfl
      · have hp' : (n / 3 + 1) % 2 = 0 := by omega
        simp only [if_pos h3, if_neg hp, if_pos hp']
        rfl
    · have hd : (n + 1) / 3 = n / 3 := by omega
      rw [hd]
      simp only [if_neg h3]
      rfl

/-- **The two-rest worldline**: start at `seamLeft 0`, two rests, cross
seam `0` forward, two rests, cross it backward, and so on. -/
def twoRestWorldline : SeamStepWorldline where
  start := seamLeft 0
  steps := twoRestSteps
  adm := by
    intro n
    rw [twoRest_portSeq]
    unfold twoRestSteps twoRestPort
    by_cases h3 : n % 3 = 2
    · by_cases hp : n / 3 % 2 = 0
      · simp only [if_pos h3, if_pos hp]
        exact rfl
      · simp only [if_pos h3, if_neg hp]
        exact rfl
    · simp only [if_neg h3]
      exact trivial

/-- Every triple `3j, 3j+1, 3j+2` of the two-rest worldline is a rest block
with `k = 2`. -/
theorem twoRestWorldline_block (j : ℕ) : RestBlock twoRestWorldline (3 * j) 2 := by
  refine ⟨fun i hi ↦ ?_, ?_⟩
  · show twoRestSteps (3 * j + i) = .rest
    unfold twoRestSteps
    have : (3 * j + i) % 3 ≠ 2 := by omega
    simp only [if_neg this]
  · show twoRestSteps (3 * j + 2) ≠ .rest
    unfold twoRestSteps
    have h3 : (3 * j + 2) % 3 = 2 := by omega
    simp only [if_pos h3]
    by_cases hp : (3 * j + 2) / 3 % 2 = 0
    · simp only [if_pos hp]
      intro h
      cases h
    · simp only [if_neg hp]
      intro h
      cases h

/-- **Two-rest worldline blocks are timelike** at every `h` with `2/3 < h`;
for `h` below the window edge the evolution is stable as well. -/
theorem twoRestWorldline_block_timelike (h : ℝ) (hl : 2 / 3 < h) (j : ℕ) :
    0 < lorentzQ (generatedPath h twoRestWorldline (3 * j + 2 + 1) -
      generatedPath h twoRestWorldline (3 * j)) :=
  (two_rest_block_timelike_iff h (by linarith) twoRestWorldline (3 * j)
    (twoRestWorldline_block j)).mpr hl

/-! ## (3) The speed reading -/

/-- The block speed: one seam over the block time `(k+1) τ`, in the length
unit of the Lorentz module (one seam = length `2`) per time unit. -/
def blockSpeed (τ : ℝ) (k : ℕ) : ℝ := 2 / ((k + 1 : ℝ) * τ)

theorem blockSpeed_pos (τ : ℝ) (hτ : 0 < τ) (k : ℕ) : 0 < blockSpeed τ k := by
  unfold blockSpeed
  positivity

theorem blockSpeed_lt_one_iff (τ : ℝ) (hτ : 0 < τ) (k : ℕ) :
    blockSpeed τ k < 1 ↔ 2 < (k + 1 : ℝ) * τ := by
  unfold blockSpeed
  rw [div_lt_one (by positivity)]

theorem blockSpeed_eq_one_iff (τ : ℝ) (hτ : 0 < τ) (k : ℕ) :
    blockSpeed τ k = 1 ↔ (k + 1 : ℝ) * τ = 2 := by
  unfold blockSpeed
  rw [div_eq_one_iff_eq (by positivity)]
  exact eq_comm

/-- **Speed and causal character.**  A block is timelike iff its speed is
below `1`, and null iff its speed is `1`. -/
theorem block_timelike_iff_speed_lt_one (τ : ℝ) (hτ : 0 < τ) (w : SeamStepWorldline)
    (n k : ℕ) (hb : RestBlock w n k) :
    (0 < lorentzQ (generatedPath τ w (n + k + 1) - generatedPath τ w n) ↔
      blockSpeed τ k < 1) ∧
    (lorentzQ (generatedPath τ w (n + k + 1) - generatedPath τ w n) = 0 ↔
      blockSpeed τ k = 1) := by
  refine ⟨by rw [block_timelike_iff τ hτ w n k hb, blockSpeed_lt_one_iff τ hτ k], ?_⟩
  rw [block_lorentzQ τ w n k hb, blockSpeed_eq_one_iff τ hτ k]
  have hx : 0 < (k + 1 : ℝ) * τ := by positivity
  constructor
  · intro h
    nlinarith
  · intro h
    rw [h]
    norm_num

/-- Every timelike block at a step in the stability window, at the
identified unit, has speed below `1` and at least two rests. -/
theorem stable_block_speed_lt_one (h : ℝ) (h0 : 0 < h) (hw : StableWindow h)
    (w : SeamStepWorldline) (n k : ℕ) (hb : RestBlock w n k)
    (ht : 0 < lorentzQ (generatedPath h w (n + k + 1) - generatedPath h w n)) :
    blockSpeed h k < 1 ∧ 2 ≤ k :=
  ⟨(block_timelike_iff_speed_lt_one h h0 w n k hb).1.mp ht,
    block_two_rests h h0 hw w n k hb ht⟩

/-- No timelike block has speed `1`. -/
theorem timelike_block_speed_ne_one (τ : ℝ) (hτ : 0 < τ) (w : SeamStepWorldline)
    (n k : ℕ) (hb : RestBlock w n k)
    (ht : 0 < lorentzQ (generatedPath τ w (n + k + 1) - generatedPath τ w n)) :
    blockSpeed τ k ≠ 1 :=
  ne_of_lt ((block_timelike_iff_speed_lt_one τ hτ w n k hb).1.mp ht)

/-- **The supremum is `1`.**  For every `ε > 0` there are a step `h` in
the window and a rest count `k` whose block is timelike with speed above
`1 - ε`; the witness uses `k = 8` rests and the step
`h = 2 / (9 (1 - ε/2))` for `ε < 1`, and `k = 2`, `h = 7/10` for `1 ≤ ε`. -/
theorem block_speed_sup_one (ε : ℝ) (hε : 0 < ε) :
    ∃ (h : ℝ) (k : ℕ), StableWindow h ∧ 0 < h ∧ 2 < (k + 1 : ℝ) * h ∧
      1 - ε < blockSpeed h k := by
  by_cases h1 : 1 ≤ ε
  · refine ⟨7 / 10, 2, ?_, by norm_num, by norm_num, ?_⟩
    · unfold StableWindow
      nlinarith [sqrt_five_lt_three]
    · have := blockSpeed_pos (7 / 10) (by norm_num) 2
      linarith
  · push Not at h1
    set d : ℝ := 1 - ε / 2 with hd
    have hd0 : 1 / 2 < d := by rw [hd]; linarith
    have hd1 : d < 1 := by rw [hd]; linarith
    have hdpos : 0 < d := by linarith
    refine ⟨2 / (9 * d), 8, ?_, by positivity, ?_, ?_⟩
    · unfold StableWindow
      have hlt : 2 / (9 * d) < 4 / 9 := by
        rw [div_lt_iff₀ (by positivity)]
        nlinarith
      have hpos : 0 < 2 / (9 * d) := by positivity
      have hsq : (2 / (9 * d)) ^ 2 < (4 / 9) ^ 2 := by
        exact pow_lt_pow_left₀ hlt hpos.le (by norm_num)
      nlinarith [sqrt_five_lt_three]
    · have : ((8 : ℕ) + 1 : ℝ) * (2 / (9 * d)) = 2 / d := by
        norm_num
        field_simp
      rw [this, lt_div_iff₀ hdpos]
      linarith
    · unfold blockSpeed
      have : ((8 : ℕ) + 1 : ℝ) * (2 / (9 * d)) = 2 / d := by
        norm_num
        field_simp
      have h2 : (2 : ℝ) / (2 / d) = d := by
        field_simp
      rw [this, h2, hd]
      linarith

/-! ## (4) Without the identification -/

/-- For a positive free unit a crossing step is timelike iff `2 < τ`; the
statement contains no field step. -/
theorem crossing_timelike_iff_pos (τ : ℝ) (hτ : 0 < τ) (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest) :
    0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) ↔ 2 < τ := by
  rw [seam_step_timelike_iff τ w k hk, evalPhi_threshold]
  constructor
  · intro h
    nlinarith
  · intro h
    nlinarith

/-- **The two thresholds are independent.**  At `τ = 3` and `h = 1/2` the
unit is not identified with the step, every crossing step is timelike, the
step lies in the window, and every zero-current solution is bounded for
every datum. -/
theorem thresholds_independent :
    ¬ UnitIdentifiedWithStep 3 (1 / 2) ∧
    (∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
      0 < lorentzQ (generatedPath 3 w (k + 1) - generatedPath 3 w k)) ∧
    StableWindow (1 / 2) ∧
    (∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
        AmpereEvolutionScaled (1 / 2) A φ (fun _ ↦ 0) → ∀ n,
          realSeamEnergy (electricFieldScaled (1 / 2) A φ n) ≤
            8 * fieldEnergyScaled (1 / 2) A φ 0 / (4 - (1 / 2) ^ 2 * (3 + Real.sqrt 5))) := by
  have hw : StableWindow (1 / 2) := by
    unfold StableWindow
    nlinarith [sqrt_five_lt_three]
  refine ⟨?_, fun w k hk ↦ seam_step_timelike_of_two_lt 3 (by norm_num) w k hk, hw,
    fun A φ hA n ↦ ((courant_threshold_sharp (1 / 2) (by norm_num)).1 hw A φ hA n).2.1⟩
  unfold UnitIdentifiedWithStep
  norm_num

/-- **The window constrains nothing about the free unit.**  For every `h`
in the window, every unit above `2` gives timelike crossings and every
nonnegative unit below `2` gives spacelike crossings, on every seam-step
worldline. -/
theorem window_does_not_constrain_unit (h : ℝ) (_hw : StableWindow h) :
    (∀ τ : ℝ, 2 < τ → ∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
      0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k)) ∧
    (∀ τ : ℝ, 0 ≤ τ → τ < 2 → ∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
      lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) < 0) := by
  refine ⟨fun τ hτ w k hk ↦ seam_step_timelike_of_two_lt τ hτ w k hk, fun τ h0 h2 w k hk ↦ ?_⟩
  rw [lorentzQ_generated_step, stepNormSq_of_ne_rest _ hk]
  nlinarith

/-- **The incompatibility is a consequence of the identification alone.**
Under the identification every crossing in the window is spacelike; without
it there are a unit and a step in the window with every crossing timelike. -/
theorem incompatibility_from_identification :
    (∀ τ h : ℝ, UnitIdentifiedWithStep τ h → StableWindow h →
      ∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
        lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) < 0) ∧
    (∃ τ h : ℝ, ¬ UnitIdentifiedWithStep τ h ∧ StableWindow h ∧
      ∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
        0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k)) := by
  refine ⟨fun τ h hid hw w k hk ↦ ?_, ⟨3, 1 / 2, thresholds_independent.1,
    thresholds_independent.2.2.1, thresholds_independent.2.1⟩⟩
  rw [hid]
  exact window_crossing_spacelike h hw w k hk

/-- The committed joined witness carries `stepDuration = 1` and field step
`h = 1/2`: an inhabitant of the joined architecture outside the
identification.  Its duration `1` lies below the threshold `2`, so its
per-step crossings are spacelike as well; the witness shows only that the
join carries the two values as separate fields. -/
theorem committedJoinedWitness_not_identified :
    ¬ UnitIdentifiedWithStep committedJoinedWitness.stepDuration
      committedJoinedWitness.scaled.h := by
  show ¬ ((1 : ℝ) = 1 / 2)
  norm_num

/-- At the free unit `τ = 2` every crossing step is null: `2² - 4 = 0`. -/
theorem duration_two_null_crossing (w : SeamStepWorldline) (k : ℕ) (hk : w.steps k ≠ .rest) :
    lorentzQ (generatedPath 2 w (k + 1) - generatedPath 2 w k) = 0 := by
  rw [lorentzQ_generated_step, stepNormSq_of_ne_rest _ hk]
  norm_num

/-- The corpus inhabitant `doubleStepJoinedWitness` has `stepDuration = 2`,
so at its duration every crossing step is null; its duration is distinct
from its field step `1/2`, so the null crossings sit outside the
identification. -/
theorem doubleStepJoinedWitness_null_crossing (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest) :
    lorentzQ (generatedPath doubleStepJoinedWitness.stepDuration w (k + 1) -
      generatedPath doubleStepJoinedWitness.stepDuration w k) = 0 ∧
    ¬ UnitIdentifiedWithStep doubleStepJoinedWitness.stepDuration
      doubleStepJoinedWitness.scaled.h := by
  refine ⟨duration_two_null_crossing w k hk, ?_⟩
  show ¬ ((2 : ℝ) = 1 / 2)
  norm_num

/-! ## The one citable composed receipt -/

/-- **Seam-step speed-limit receipt.**  Under the declared identification
`τ = h`: the window and a timelike per-step crossing are incompatible; the
thresholds differ by the factor `3 + √5`; every timelike block of a stable
evolution has at least two rests and speed below `1`; the two-rest worldline
has timelike blocks for `2/3 < h`; the supremum of timelike block speeds
over the window is `1`, never attained.  Without the identification the
thresholds are independent. -/
theorem seamStepSpeedLimit_receipt :
    (∀ h : ℝ, ¬ (StableWindow h ∧ 4 < h ^ 2)) ∧
    (4 / (2 / Real.goldenRatio ^ 2) = 3 + Real.sqrt 5) ∧
    (∀ h : ℝ, 0 < h → StableWindow h → ∀ (w : SeamStepWorldline) (n k : ℕ), RestBlock w n k →
      0 < lorentzQ (generatedPath h w (n + k + 1) - generatedPath h w n) →
        2 ≤ k ∧ blockSpeed h k < 1) ∧
    (∀ h : ℝ, 2 / 3 < h → ∀ j : ℕ,
      0 < lorentzQ (generatedPath h twoRestWorldline (3 * j + 2 + 1) -
        generatedPath h twoRestWorldline (3 * j))) ∧
    (∀ ε : ℝ, 0 < ε → ∃ (h : ℝ) (k : ℕ), StableWindow h ∧ 0 < h ∧
      2 < (k + 1 : ℝ) * h ∧ 1 - ε < blockSpeed h k) ∧
    (∃ τ h : ℝ, ¬ UnitIdentifiedWithStep τ h ∧ StableWindow h ∧
      ∀ (w : SeamStepWorldline) (k : ℕ), w.steps k ≠ .rest →
        0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k)) :=
  ⟨window_timelike_incompatible,
    by rw [threshold_gap.1, threshold_gap.2],
    fun h h0 hw w n k hb ht ↦ ⟨block_two_rests h h0 hw w n k hb ht,
      (stable_block_speed_lt_one h h0 hw w n k hb ht).1⟩,
    twoRestWorldline_block_timelike,
    block_speed_sup_one,
    incompatibility_from_identification.2⟩

end

#print axioms OPH.SeamStepSpeedLimit.crossing_timelike_iff_identified
#print axioms OPH.SeamStepSpeedLimit.window_sq_lt_one
#print axioms OPH.SeamStepSpeedLimit.window_timelike_incompatible
#print axioms OPH.SeamStepSpeedLimit.stable_evolution_crossing_spacelike
#print axioms OPH.SeamStepSpeedLimit.timelike_crossing_unstable
#print axioms OPH.SeamStepSpeedLimit.threshold_gap
#print axioms OPH.SeamStepSpeedLimit.joined_identified_crossing_spacelike
#print axioms OPH.SeamStepSpeedLimit.block_lorentzQ
#print axioms OPH.SeamStepSpeedLimit.block_timelike_iff
#print axioms OPH.SeamStepSpeedLimit.block_rests_gt
#print axioms OPH.SeamStepSpeedLimit.block_two_rests
#print axioms OPH.SeamStepSpeedLimit.two_rest_block_timelike_iff
#print axioms OPH.SeamStepSpeedLimit.two_thirds_lt_edge
#print axioms OPH.SeamStepSpeedLimit.two_rest_admissible_interval
#print axioms OPH.SeamStepSpeedLimit.twoRestWorldline_block
#print axioms OPH.SeamStepSpeedLimit.twoRestWorldline_block_timelike
#print axioms OPH.SeamStepSpeedLimit.block_timelike_iff_speed_lt_one
#print axioms OPH.SeamStepSpeedLimit.stable_block_speed_lt_one
#print axioms OPH.SeamStepSpeedLimit.timelike_block_speed_ne_one
#print axioms OPH.SeamStepSpeedLimit.block_speed_sup_one
#print axioms OPH.SeamStepSpeedLimit.crossing_timelike_iff_pos
#print axioms OPH.SeamStepSpeedLimit.thresholds_independent
#print axioms OPH.SeamStepSpeedLimit.window_does_not_constrain_unit
#print axioms OPH.SeamStepSpeedLimit.incompatibility_from_identification
#print axioms OPH.SeamStepSpeedLimit.committedJoinedWitness_not_identified
#print axioms OPH.SeamStepSpeedLimit.doubleStepJoinedWitness_null_crossing
#print axioms OPH.SeamStepSpeedLimit.seamStepSpeedLimit_receipt

end OPH.SeamStepSpeedLimit
