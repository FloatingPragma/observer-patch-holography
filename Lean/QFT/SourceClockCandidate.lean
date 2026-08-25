import PhysicalCalibrationImport
import QFT.SourceHistoryThreeSlotLocalGNS

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A source-counted overlapping-window action statistic

WHAT IS PROVED.  The legacy-named structure `SourceClockCandidate`
packages a histogram and counting observables from the committed
source-history packet, together with exact rational averages.  It does
not construct an ordered clock process.

(1) The statistic.  `tickCount` is the per-history multiplicity
histogram `sourceWindowCount`, whose total is `1754`; it is not a
time-ordered tick stream.  For `repairClockCandidate`, `eventCount` is
the sum of the two domain-wall indicators in each retained length-three
window.  The count-weighted total `197 = 94 + 103` is therefore an
overlapping-window, two-edge action-incidence count: `94` counts the
first edge of a window and `103` the second.  Because consecutive
windows overlap, `197` is not a count of distinct underlying source
transitions or repair events.  The module proves the corresponding
means `197/1754`, `94/1754`, and `103/1754`, their exact split, and their
agreement with the committed diagonal-state expectations.

(2) Conditional unit conversion.  `windowIdentification` is an
unconstrained proposition and is never proved.  More strongly,
`toClockCalibration` deliberately ignores its candidate argument and
uses only an externally declared positive real `tau_w`.  The quantities
named `labRate` are consequently just the algebraic quotients of the
dimensionless incidence averages by that declaration.  The frequency
and cancellation theorems are exact identities conditional on the same
declaration; they do not identify a physical duration or validate an
event-rate interpretation.

(3) Delimitation.
`conversion_dictionary_depends_only_on_declared_duration` says only that
the generated conversion dictionary is determined by the declared
duration; it does not reduce the empirical calibration burden using the
histogram.  `wallOnlyClockCandidate` is an alternative observable on the
same histogram, not a second observed subprocess.  Nothing here gives
window order, successor data, cadence, periodicity, a physical seconds
value, a source-selected observable, or a clock.  The retained legacy
names record the intended calibration proposal, not a theorem that the
proposal is physically realized.
-/

namespace OPH.QFT.SourceClock

noncomputable section

open OPH.PhysicalCalibrationImport
open OPH.InformationProjection
open OPH.QFT.SourceHistoryGNSDynamics
open OPH.QFT.SourceHistoryThreeSlotLocalGNS

/-! ## (1) The candidate: a source-counted window statistic -/

/-- A legacy-named package for a committed window histogram, a chosen
counting observable, and its exact empirical averages.  The histogram
contains no ordering or successor relation and hence is not itself a
tick process.  The single open identification field
`windowIdentification` is carried as an unproved proposition. -/
structure SourceClockCandidate where
  /-- The per-history multiplicity histogram of retained windows. -/
  tickCount : Fin 8 → ℕ
  /-- A chosen per-window counting observable. -/
  eventCount : Fin 8 → ℕ
  /-- The `(0,1)`-wall component of the event observable. -/
  wall01Count : Fin 8 → ℕ
  /-- The `(1,2)`-wall component of the event observable. -/
  wall12Count : Fin 8 → ℕ
  /-- The exact source-counted mean of the chosen observable per window. -/
  ratePerWindow : ℚ
  /-- The exact mean first-edge incidence per window. -/
  wallRate01 : ℚ
  /-- The exact mean second-edge incidence per window. -/
  wallRate12 : ℚ
  /-- Source-fixed: the histogram is the committed retained-window count. -/
  tick_committed : tickCount = sourceWindowCount
  /-- Source-fixed: the event observable splits exactly into the two
  wall components. -/
  event_split : ∀ g : Fin 8, eventCount g = wall01Count g + wall12Count g
  /-- Source-fixed: this field is the exact mean of the counting observable
  under the committed empirical law. -/
  rate_committed : ratePerWindow = ∑ g, sourceTauEmpQ g * (eventCount g : ℚ)
  /-- Source-fixed: the `(0,1)`-wall rate is the exact mean of its
  component. -/
  wall01_committed : wallRate01 = ∑ g, sourceTauEmpQ g * (wall01Count g : ℚ)
  /-- Source-fixed: the `(1,2)`-wall rate is the exact mean of its
  component. -/
  wall12_committed : wallRate12 = ∑ g, sourceTauEmpQ g * (wall12Count g : ℚ)
  /-- Source-fixed positivity of the mean statistic. -/
  rate_pos : 0 < ratePerWindow
  /-- An unproved identification schema.  No field relates this
  proposition to a duration, and `toClockCalibration` below does not
  consume it. -/
  windowIdentification : Prop

/-- **Wall split of the mean statistics.**  For every candidate the two
component means sum exactly to the total mean; this is an identity of the
source-fixed fields, not an extra assumption. -/
theorem SourceClockCandidate.wall_split (c : SourceClockCandidate) :
    c.wallRate01 + c.wallRate12 = c.ratePerWindow := by
  rw [c.wall01_committed, c.wall12_committed, c.rate_committed,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun g _ => ?_
  rw [c.event_split g]
  push_cast
  ring

/-- **What the source fixes.**  Every candidate carries the committed
window histogram; its statistic is the exact empirical mean of its
chosen counting observable, with the wall split and positivity shown.
This receipt contains neither window order nor a clock identification. -/
theorem SourceClockCandidate.source_fixed_receipt (c : SourceClockCandidate) :
    c.tickCount = sourceWindowCount
      ∧ c.ratePerWindow = ∑ g, sourceTauEmpQ g * (c.eventCount g : ℚ)
      ∧ c.wallRate01 + c.wallRate12 = c.ratePerWindow
      ∧ 0 < c.ratePerWindow :=
  ⟨c.tick_committed, c.rate_committed, c.wall_split, c.rate_pos⟩

/-- The legacy-named repair statistic on the committed window histogram.
It counts the two domain-wall incidences inside each overlapping window,
with exact mean `197/1754` and position split `94/1754 + 103/1754`.
Thus `197` is an overlapping two-edge incidence total, not the number of
distinct underlying repair events.  The identification field is only a
stand-in schema saying that some positive real exists. -/
def repairClockCandidate : SourceClockCandidate where
  tickCount := sourceWindowCount
  eventCount := sourceAction
  wall01Count := bond01Energy
  wall12Count := bond12Energy
  ratePerWindow := 197 / 1754
  wallRate01 := 94 / 1754
  wallRate12 := 103 / 1754
  tick_committed := rfl
  event_split := by decide
  rate_committed := sourceTauEmpQ_meanAction.symm
  wall01_committed := sourceTauEmpQ_bond01_mean.symm
  wall12_committed := sourceTauEmpQ_bond12_mean.symm
  rate_pos := by norm_num
  windowIdentification := ∃ tau_w : ℝ, 0 < tau_w

/-- The repair candidate uses the committed window histogram. -/
@[simp] theorem repairClockCandidate_tickCount :
    repairClockCandidate.tickCount = sourceWindowCount := rfl

/-- The repair candidate's observable is the two-edge window action. -/
@[simp] theorem repairClockCandidate_eventCount :
    repairClockCandidate.eventCount = sourceAction := rfl

/-- The exact action-incidence mean `197/1754` per overlapping window. -/
@[simp] theorem repairClockCandidate_ratePerWindow :
    repairClockCandidate.ratePerWindow = 197 / 1754 := rfl

/-- The exact first-edge incidence mean `94/1754` per window. -/
@[simp] theorem repairClockCandidate_wallRate01 :
    repairClockCandidate.wallRate01 = 94 / 1754 := rfl

/-- The exact second-edge incidence mean `103/1754` per window. -/
@[simp] theorem repairClockCandidate_wallRate12 :
    repairClockCandidate.wallRate12 = 103 / 1754 := rfl

/-- The committed histogram totals `1754` retained overlapping windows. -/
theorem repairClockCandidate_tick_total :
    ∑ g, repairClockCandidate.tickCount g = 1754 := by
  show ∑ g, sourceWindowCount g = 1754
  rw [Fin.sum_univ_eight]
  exact sourceWindowCount_total

/-- The count-weighted action-incidence total is `197`.  It sums both
edges of overlapping windows and is not an underlying event count. -/
theorem repairClockCandidate_event_total :
    ∑ g, repairClockCandidate.tickCount g * repairClockCandidate.eventCount g
      = 197 := by
  show ∑ g, sourceWindowCount g * sourceAction g = 197
  rw [Fin.sum_univ_eight]
  exact sourceWindowCount_action_total

/-- The wall decomposition of the action mean, as the structure identity
and as the exact rational identity `94/1754 + 103/1754 = 197/1754`. -/
theorem repairClockCandidate_wall_decomposition :
    repairClockCandidate.wallRate01 + repairClockCandidate.wallRate12
        = repairClockCandidate.ratePerWindow
      ∧ (94 / 1754 : ℚ) + 103 / 1754 = 197 / 1754 :=
  ⟨repairClockCandidate.wall_split, by norm_num⟩

/-- **The statistics are committed source expectations.**  The action
incidence mean is the expectation of the repair Hamiltonian in the
committed source state, and its position components are the expectations
of the two Ising domain-wall bonds. -/
theorem repairClockCandidate_rates_are_source_expectations :
    (repairClockCandidate.ratePerWindow : ℂ)
        = EventAlgebra.expectation sourceHistoryDensity
            sourceHistoryHamiltonian
      ∧ (repairClockCandidate.wallRate01 : ℂ)
        = EventAlgebra.expectation sourceHistoryDensity bond01
      ∧ (repairClockCandidate.wallRate12 : ℂ)
        = EventAlgebra.expectation sourceHistoryDensity bond12 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [sourceHistoryDensity_meanEnergy, repairClockCandidate_ratePerWindow]
    norm_num
  · rw [sourceHistoryDensity_bond01_expectation,
      repairClockCandidate_wallRate01]
    norm_num
  · rw [sourceHistoryDensity_bond12_expectation,
      repairClockCandidate_wallRate12]
    norm_num

/-- **The identification field is unconstrained.**  For every
proposition there is a candidate carrying it in the identification field
while agreeing with the repair candidate on the histogram and the mean
statistic: the structure fixes these data and never asserts, selects, or
discharges the identification. -/
theorem identification_field_free (p : Prop) :
    ∃ c : SourceClockCandidate,
      c.windowIdentification = p
        ∧ c.tickCount = repairClockCandidate.tickCount
        ∧ c.ratePerWindow = repairClockCandidate.ratePerWindow :=
  ⟨{ repairClockCandidate with windowIdentification := p }, rfl, rfl, rfl⟩

/-! ## (2) Composition with the calibration import -/

/-- A declared positive real makes a `ClockCalibration`.  The candidate
argument is deliberately unused, and neither its histogram nor its
`windowIdentification` field is consulted.  This is an external unit
declaration, not a source-derived clock calibration. -/
def SourceClockCandidate.toClockCalibration (_c : SourceClockCandidate)
    (tau_w : ℝ) (h : 0 < tau_w) : ClockCalibration where
  tau := tau_w
  tau_pos := h

/-- The composed calibration stores the externally declared `tau_w`. -/
@[simp] theorem toClockCalibration_tau (c : SourceClockCandidate)
    (tau_w : ℝ) (h : 0 < tau_w) :
    (c.toClockCalibration tau_w h).tau = tau_w := rfl

/-- The conditional quotient `ratePerWindow / tau` under a declared
calibration.  Calling it a laboratory rate does not establish that the
histogram is an ordered process or that `tau` is a window duration. -/
def SourceClockCandidate.labRate (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.ratePerWindow : ℝ) / cal.tau

/-- The first-edge mean divided by the declared calibration parameter. -/
def SourceClockCandidate.labWallRate01 (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.wallRate01 : ℝ) / cal.tau

/-- The second-edge mean divided by the declared calibration parameter. -/
def SourceClockCandidate.labWallRate12 (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.wallRate12 : ℝ) / cal.tau

/-- Under any declared calibration the action-incidence mean gives the
exact conditional quotient `(197/1754) / tau`. -/
theorem repairClockCandidate_labRate (cal : ClockCalibration) :
    repairClockCandidate.labRate cal = 197 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labRate, repairClockCandidate_ratePerWindow]
  norm_num

/-- Through `toClockCalibration`, a declared positive `tau_w` gives the
exact conditional quotient `(197/1754) / tau_w`; the candidate is not
used to identify `tau_w`. -/
theorem repairClockCandidate_labRate_declared (tau_w : ℝ) (h : 0 < tau_w) :
    repairClockCandidate.labRate
        (repairClockCandidate.toClockCalibration tau_w h)
      = 197 / 1754 / tau_w :=
  repairClockCandidate_labRate _

/-- The first-edge incidence mean gives `(94/1754) / tau` conditionally. -/
theorem repairClockCandidate_labWallRate01 (cal : ClockCalibration) :
    repairClockCandidate.labWallRate01 cal = 94 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labWallRate01,
    repairClockCandidate_wallRate01]
  norm_num

/-- The second-edge incidence mean gives `(103/1754) / tau` conditionally. -/
theorem repairClockCandidate_labWallRate12 (cal : ClockCalibration) :
    repairClockCandidate.labWallRate12 cal = 103 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labWallRate12,
    repairClockCandidate_wallRate12]
  norm_num

/-- **Exact quotient sum identity.**  For every candidate and declared
calibration the two component quotients sum to the total quotient. -/
theorem SourceClockCandidate.labWall_split (c : SourceClockCandidate)
    (cal : ClockCalibration) :
    c.labWallRate01 cal + c.labWallRate12 cal = c.labRate cal := by
  simp only [SourceClockCandidate.labWallRate01,
    SourceClockCandidate.labWallRate12, SourceClockCandidate.labRate]
  rw [← add_div, ← Rat.cast_add, c.wall_split]

/-- **Algebraic composition with the frequency dictionary.**  The
conditional quotient equals the import's frequency conversion of
`2 * pi * ratePerWindow`.  This identity supplies no process or physical
frequency interpretation. -/
theorem SourceClockCandidate.labRate_eq_labFrequency
    (c : SourceClockCandidate) (cal : ClockCalibration) :
    c.labRate cal
      = cal.labFrequency (2 * Real.pi * (c.ratePerWindow : ℝ)) := by
  have h2pi : (2 * Real.pi : ℝ) ≠ 0 :=
    mul_ne_zero two_ne_zero Real.pi_ne_zero
  simp only [SourceClockCandidate.labRate,
    ClockCalibration.labFrequency_def]
  rw [mul_div_mul_left _ _ h2pi]

/-- **Algebraic cancellation.**  The conditional quotient times the
declared duration of `n` units is `n * ratePerWindow`.  This is exact
real algebra, not an independently calibrated event count. -/
theorem SourceClockCandidate.labRate_mul_labSeconds
    (c : SourceClockCandidate) (cal : ClockCalibration) (n : ℝ) :
    c.labRate cal * cal.labSeconds n = n * (c.ratePerWindow : ℝ) := by
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  simp only [SourceClockCandidate.labRate, ClockCalibration.labSeconds_def]
  calc (c.ratePerWindow : ℝ) / cal.tau * (n * cal.tau)
      = n * ((c.ratePerWindow : ℝ) * (cal.tau / cal.tau)) := by ring
    _ = n * (c.ratePerWindow : ℝ) := by rw [div_self hτ, mul_one]

/-! ## (3) Dependence on the single declaration -/

/-- Two clock-calibration records with the same declared `tau` are equal:
the dictionary carries no data beyond that parameter. -/
theorem calibration_eq_of_tau_eq {cal cal' : ClockCalibration}
    (h : cal.tau = cal'.tau) : cal = cal' := by
  obtain ⟨t, ht⟩ := cal
  obtain ⟨t', ht'⟩ := cal'
  have h' : t = t' := h
  subst h'
  rfl

/-- **Dependence on the declaration alone.**  Two calls agreeing on the
declared duration agree everywhere, even for different candidates,
because `toClockCalibration` ignores the candidate.  This is not an
empirical reduction supplied by the source histogram. -/
theorem reduction_to_single_identification (c c' : SourceClockCandidate)
    (t t' : ℝ) (h : 0 < t) (h' : 0 < t') (he : t = t') :
    c.toClockCalibration t h = c'.toClockCalibration t' h'
      ∧ (∀ x : ℝ, (c.toClockCalibration t h).labSeconds x
          = (c'.toClockCalibration t' h').labSeconds x)
      ∧ (∀ m : ℝ, (c.toClockCalibration t h).labFrequency m
          = (c'.toClockCalibration t' h').labFrequency m) := by
  have htau : (c.toClockCalibration t h).tau
      = (c'.toClockCalibration t' h').tau := by
    simp only [toClockCalibration_tau]
    exact he
  have hcal := calibration_eq_of_tau_eq htau
  exact ⟨hcal, fun x => by rw [hcal], fun m => by rw [hcal]⟩

/-- Accurate public alias: the generated conversion dictionary depends only
on the declared duration.  The legacy theorem name above is retained for API
stability and does not express a reduction of the empirical calibration
premise. -/
theorem conversion_dictionary_depends_only_on_declared_duration
    (c c' : SourceClockCandidate) (t t' : ℝ) (h : 0 < t) (h' : 0 < t')
    (he : t = t') :
    c.toClockCalibration t h = c'.toClockCalibration t' h'
      ∧ (∀ x : ℝ, (c.toClockCalibration t h).labSeconds x
          = (c'.toClockCalibration t' h').labSeconds x)
      ∧ (∀ m : ℝ, (c.toClockCalibration t h).labFrequency m
          = (c'.toClockCalibration t' h').labFrequency m) :=
  reduction_to_single_identification c c' t t' h h' he

/-- **Non-forcing.**  The identification is an import, never derived:
for any positive-mean candidate, two calibrations with distinct declared
ticks give distinct quotients.  No internal statement fixes a physical
rate or makes the histogram into a process. -/
theorem SourceClockCandidate.labRate_not_forced (c : SourceClockCandidate)
    (cal cal' : ClockCalibration) (hne : cal.tau ≠ cal'.tau) :
    c.labRate cal ≠ c.labRate cal' := by
  simp only [SourceClockCandidate.labRate]
  have hrpos : (0 : ℝ) < (c.ratePerWindow : ℝ) := by
    exact_mod_cast c.rate_pos
  intro heq
  apply hne
  have h2 := (div_eq_div_iff cal.tau_pos.ne' cal'.tau_pos.ne').mp heq
  exact (mul_left_cancel₀ hrpos.ne' h2).symm

/-! ## (4) Delimitation receipts -/

/-- **Window statistic, not a rate per second.**  The empirical mean is
declaration-invariant, while its quotient changes with the declared
duration.  The packet fixes an overlapping-window statistic; no theorem
produces a duration or an ordered process. -/
theorem rate_is_per_window_not_per_second (c : SourceClockCandidate)
    (cal cal' : ClockCalibration) (hne : cal.tau ≠ cal'.tau) :
    c.ratePerWindow = c.ratePerWindow
      ∧ c.labRate cal ≠ c.labRate cal' :=
  ⟨rfl, c.labRate_not_forced cal cal' hne⟩

/-- An alternative observable on the same committed histogram, counting
only first-edge wall incidences, with mean `94/1754`.  This is not a
separately observed subprocess. -/
def wallOnlyClockCandidate : SourceClockCandidate where
  tickCount := sourceWindowCount
  eventCount := bond01Energy
  wall01Count := bond01Energy
  wall12Count := fun _ => 0
  ratePerWindow := 94 / 1754
  wallRate01 := 94 / 1754
  wallRate12 := 0
  tick_committed := rfl
  event_split := fun g => (Nat.add_zero (bond01Energy g)).symm
  rate_committed := sourceTauEmpQ_bond01_mean.symm
  wall01_committed := sourceTauEmpQ_bond01_mean.symm
  wall12_committed := by simp
  rate_pos := by norm_num
  windowIdentification := ∃ tau_w : ℝ, 0 < tau_w

/-- The wall-only observable uses the same committed histogram. -/
@[simp] theorem wallOnlyClockCandidate_tickCount :
    wallOnlyClockCandidate.tickCount = sourceWindowCount := rfl

/-- The counting observable of the wall-only candidate is the
`(0,1)`-wall count. -/
@[simp] theorem wallOnlyClockCandidate_eventCount :
    wallOnlyClockCandidate.eventCount = bond01Energy := rfl

/-- The exact wall-only incidence mean `94/1754` per window. -/
@[simp] theorem wallOnlyClockCandidate_ratePerWindow :
    wallOnlyClockCandidate.ratePerWindow = 94 / 1754 := rfl

/-- The first-edge incidence mean of the wall-only observable. -/
@[simp] theorem wallOnlyClockCandidate_wallRate01 :
    wallOnlyClockCandidate.wallRate01 = 94 / 1754 := rfl

/-- The second-edge incidence mean of the wall-only observable vanishes. -/
@[simp] theorem wallOnlyClockCandidate_wallRate12 :
    wallOnlyClockCandidate.wallRate12 = 0 := rfl

/-- The wall-only mean gives the conditional quotient `(94/1754) / tau`. -/
theorem wallOnlyClockCandidate_labRate (cal : ClockCalibration) :
    wallOnlyClockCandidate.labRate cal = 94 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labRate,
    wallOnlyClockCandidate_ratePerWindow]
  norm_num

/-- **Observable selection is a declared choice.**  Two observables on
the same committed histogram have different exact means, `94/1754` and
`197/1754`.  The packet does not select one as a physical clock signal. -/
theorem candidate_selection_declared :
    repairClockCandidate.tickCount = wallOnlyClockCandidate.tickCount
      ∧ repairClockCandidate.ratePerWindow
          ≠ wallOnlyClockCandidate.ratePerWindow := by
  refine ⟨rfl, ?_⟩
  rw [repairClockCandidate_ratePerWindow, wallOnlyClockCandidate_ratePerWindow]
  norm_num

/-- Under one shared declaration the two candidates attach to distinct
laboratory rates: the declared candidate selection carries through to
the composed dictionary. -/
theorem candidate_choice_changes_lab_rate (cal : ClockCalibration) :
    repairClockCandidate.labRate cal ≠ wallOnlyClockCandidate.labRate cal := by
  rw [repairClockCandidate_labRate, wallOnlyClockCandidate_labRate]
  intro h
  rw [div_eq_div_iff cal.tau_pos.ne' cal.tau_pos.ne'] at h
  have h2 := mul_right_cancel₀ cal.tau_pos.ne' h
  norm_num at h2

/-! ## Bundled receipt -/

/-- **Legacy source-clock-candidate receipt.**  The bundled statement
contains a `1754`-window histogram, `197` overlapping two-edge action
incidences (not `197` distinct repair events), the exact mean and wall
split, conditional quotients under an arbitrary declared calibration,
and an alternative observable on the same histogram.  It does not add
order, cadence, a duration identification, or a realized clock process. -/
theorem sourceClockCandidate_receipt :
    (∑ g, repairClockCandidate.tickCount g = 1754)
      ∧ (∑ g, repairClockCandidate.tickCount g
            * repairClockCandidate.eventCount g = 197)
      ∧ repairClockCandidate.ratePerWindow = 197 / 1754
      ∧ repairClockCandidate.wallRate01 + repairClockCandidate.wallRate12
          = repairClockCandidate.ratePerWindow
      ∧ (∀ cal : ClockCalibration,
          repairClockCandidate.labRate cal = 197 / 1754 / cal.tau)
      ∧ (∀ cal : ClockCalibration,
          repairClockCandidate.labWallRate01 cal
              + repairClockCandidate.labWallRate12 cal
            = repairClockCandidate.labRate cal)
      ∧ (∀ (c c' : SourceClockCandidate) (t t' : ℝ)
          (h : 0 < t) (h' : 0 < t'), t = t' →
            c.toClockCalibration t h = c'.toClockCalibration t' h')
      ∧ (∀ (c : SourceClockCandidate) (cal cal' : ClockCalibration),
          cal.tau ≠ cal'.tau → c.labRate cal ≠ c.labRate cal')
      ∧ repairClockCandidate.tickCount = wallOnlyClockCandidate.tickCount
      ∧ repairClockCandidate.ratePerWindow
          ≠ wallOnlyClockCandidate.ratePerWindow :=
  ⟨repairClockCandidate_tick_total, repairClockCandidate_event_total,
    repairClockCandidate_ratePerWindow, repairClockCandidate.wall_split,
    repairClockCandidate_labRate,
    fun cal => repairClockCandidate.labWall_split cal,
    fun c c' t t' h h' he =>
      (conversion_dictionary_depends_only_on_declared_duration
        c c' t t' h h' he).1,
    fun c cal cal' hne => c.labRate_not_forced cal cal' hne,
    candidate_selection_declared.1, candidate_selection_declared.2⟩

/-! ## Axiom audit -/

#print axioms SourceClockCandidate.wall_split
#print axioms SourceClockCandidate.source_fixed_receipt
#print axioms repairClockCandidate_tickCount
#print axioms repairClockCandidate_eventCount
#print axioms repairClockCandidate_ratePerWindow
#print axioms repairClockCandidate_wallRate01
#print axioms repairClockCandidate_wallRate12
#print axioms repairClockCandidate_tick_total
#print axioms repairClockCandidate_event_total
#print axioms repairClockCandidate_wall_decomposition
#print axioms repairClockCandidate_rates_are_source_expectations
#print axioms identification_field_free
#print axioms toClockCalibration_tau
#print axioms repairClockCandidate_labRate
#print axioms repairClockCandidate_labRate_declared
#print axioms repairClockCandidate_labWallRate01
#print axioms repairClockCandidate_labWallRate12
#print axioms SourceClockCandidate.labWall_split
#print axioms SourceClockCandidate.labRate_eq_labFrequency
#print axioms SourceClockCandidate.labRate_mul_labSeconds
#print axioms calibration_eq_of_tau_eq
#print axioms reduction_to_single_identification
#print axioms SourceClockCandidate.labRate_not_forced
#print axioms rate_is_per_window_not_per_second
#print axioms wallOnlyClockCandidate_tickCount
#print axioms wallOnlyClockCandidate_eventCount
#print axioms wallOnlyClockCandidate_ratePerWindow
#print axioms wallOnlyClockCandidate_wallRate01
#print axioms wallOnlyClockCandidate_wallRate12
#print axioms wallOnlyClockCandidate_labRate
#print axioms candidate_selection_declared
#print axioms candidate_choice_changes_lab_rate
#print axioms sourceClockCandidate_receipt

end

end OPH.QFT.SourceClock
