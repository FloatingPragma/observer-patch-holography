import PhysicalCalibrationImport
import QFT.SourceHistoryThreeSlotLocalGNS

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A source-derived clock candidate for the calibration import

WHAT IS PROVED.  The calibration import of
`Thermodynamics/PhysicalCalibrationImport.lean` reduces PR-15 to one
declared tick; its open item is a source-produced process able to serve
as the tick candidate.  This module types such a candidate from the
committed source-history packet of
`InformationProjection/SourceHistoryPacket.lean` and proves its
composition and reduction receipts against the import dictionary.

(1) The candidate.  `SourceClockCandidate` is a source-produced counting
process on the committed packet: the retained-window sequence is the
tick stream (`tick_committed` forces the per-path window multiplicities
to the committed counts `sourceWindowCount`, with total `1754` windows,
`repairClockCandidate_tick_total`), the per-window counting observable
splits into two candidate-supplied wall components whose values the
repair candidate pins to the committed domain-wall counts
(`event_split`), and the exact per-window rates are the means of those
observables under the committed empirical law (`rate_committed`,
`wall01_committed`, `wall12_committed`), positive by `rate_pos`.  The
wall split of the rates is an identity of the structure (`wall_split`).
The inhabitant `repairClockCandidate` counts the repair moves: its rate
invariants are the exact source-counted mean action `197/1754` per
window with the wall decomposition `94/1754 + 103/1754`, and these
equal the committed source expectations of the repair Hamiltonian and
its two Ising bonds
(`repairClockCandidate_rates_are_source_expectations`).  What the source
does not fix is the physical duration of one window, typed as the single
open identification field `windowIdentification`, a proposition the
module never asserts (`identification_field_free`).

(2) Composition with the import.  A declared window duration `tau_w`
makes the candidate a `ClockCalibration` (`toClockCalibration`), and
every dimensionful statement is phrased against a declared calibration
whose tick is read as the window duration.  Under it the source repair
rate attaches to the exact laboratory rate `(197/1754) / tau_w` events
per second (`repairClockCandidate_labRate`,
`repairClockCandidate_labRate_declared`), the wall rates to
`(94/1754) / tau_w` and `(103/1754) / tau_w`
(`repairClockCandidate_labWallRate01`,
`repairClockCandidate_labWallRate12`) with their exact sum identity
(`labWall_split`); the counting rate composes with the import's
frequency dictionary through the exact `2 * pi` conversion
(`labRate_eq_labFrequency`), and rate times duration cancels exactly to
the dimensionless expected event count (`labRate_mul_labSeconds`).
Every such statement is conditional on the declaration and
exact-symbolic; no numerical evaluation enters.

(3) Reduction.  The composed dictionary is fully determined by the
single identification: two candidate calibrations agreeing on the
declared duration agree everywhere
(`reduction_to_single_identification`), and the identification is an
import, never derived: for any candidate, two calibrations with
distinct ticks give distinct laboratory rates (`labRate_not_forced`).
The candidate improves on the bare declared tick exactly here: with the
bare import both the tick process and its duration are declared, while
the candidate's tick stream and rate invariants are source-produced
committed data (`source_fixed_receipt`), so the open empirical content
shrinks from an arbitrary declared tick to one duration of one
committed process.  The shrinkage is a change of bookkeeping, not of
physics: one real number is open in both readings.

(4) Delimitation.  The packet fixes rates per window, not per second:
the per-window rate is declaration-invariant while the laboratory rate
changes with the declared duration, and no theorem produces the
duration (`rate_is_per_window_not_per_second`).  The candidate is one
committed process among possibly many: the wall-only subprocess is a
second candidate on the same committed tick stream with the different
exact rate invariant `94/1754` (`wallOnlyClockCandidate`,
`candidate_selection_declared`), and the two candidates attach to
distinct laboratory rates under one shared declaration
(`candidate_choice_changes_lab_rate`), so candidate selection is a
declared choice.

WHAT IS NOT PROVED HERE.  PR-15 keeps its import disposition: no
physical seconds value enters, the identification of one retained
window with a laboratory duration is the open empirical content of the
import, and no theorem selects, supplies, or discharges it.  Uniformity
of the window duration across windows is part of the same open
identification: the composition declares one duration for every window,
and the source proves no periodicity in laboratory time.  The candidate
is a proposal pending registration, not a discharged premise, a frozen
record, or a scored comparison.  Candidate selection among
source-produced counting processes is a declared choice, witnessed
inside the module.  The stand-in propositions placed in the
identification field of the exhibited candidates are schemas, not
physical identifications.  Nothing here calibrates any prediction
surface.
-/

namespace OPH.QFT.SourceClock

noncomputable section

open OPH.PhysicalCalibrationImport
open OPH.InformationProjection
open OPH.QFT.SourceHistoryGNSDynamics
open OPH.QFT.SourceHistoryThreeSlotLocalGNS

/-! ## (1) The candidate: a source-produced counting process -/

/-- A source-produced counting process from the committed history
packet, offered as a tick candidate for the calibration import.  The
tick stream is the committed retained-window sequence
(`tick_committed`), the counting observable splits into the two Ising
domain-wall components (`event_split`), and the per-window rates are the
exact means of those observables under the committed empirical law
(`rate_committed`, `wall01_committed`, `wall12_committed`), with the
rate positive (`rate_pos`).  These fields are what the source genuinely
fixes.  What the source does not fix is the physical duration of one
window, typed as the single open identification field
`windowIdentification`. -/
structure SourceClockCandidate where
  /-- The tick stream: the per-path window multiplicities of the
  retained-window sequence of the committed packet. -/
  tickCount : Fin 8 → ℕ
  /-- The per-window event observable counted by the candidate. -/
  eventCount : Fin 8 → ℕ
  /-- The `(0,1)`-wall component of the event observable. -/
  wall01Count : Fin 8 → ℕ
  /-- The `(1,2)`-wall component of the event observable. -/
  wall12Count : Fin 8 → ℕ
  /-- The exact source-counted event rate per window. -/
  ratePerWindow : ℚ
  /-- The exact `(0,1)`-wall rate per window. -/
  wallRate01 : ℚ
  /-- The exact `(1,2)`-wall rate per window. -/
  wallRate12 : ℚ
  /-- Source-fixed: the tick stream is the committed retained-window
  count. -/
  tick_committed : tickCount = sourceWindowCount
  /-- Source-fixed: the event observable splits exactly into the two
  wall components. -/
  event_split : ∀ g : Fin 8, eventCount g = wall01Count g + wall12Count g
  /-- Source-fixed: the rate is the exact mean of the event observable
  under the committed empirical law. -/
  rate_committed : ratePerWindow = ∑ g, sourceTauEmpQ g * (eventCount g : ℚ)
  /-- Source-fixed: the `(0,1)`-wall rate is the exact mean of its
  component. -/
  wall01_committed : wallRate01 = ∑ g, sourceTauEmpQ g * (wall01Count g : ℚ)
  /-- Source-fixed: the `(1,2)`-wall rate is the exact mean of its
  component. -/
  wall12_committed : wallRate12 = ∑ g, sourceTauEmpQ g * (wall12Count g : ℚ)
  /-- Source-fixed positivity of the rate. -/
  rate_pos : 0 < ratePerWindow
  /-- The single open identification field: the declarer's statement
  assigning one retained window a physical duration in SI seconds.  The
  structure carries the statement without asserting it: no field demands
  a proof, and no theorem in the corpus supplies, selects, or discharges
  it.  A dimensionful reading enters only through a declared
  `ClockCalibration`. -/
  windowIdentification : Prop

/-- **Wall split of the rate invariants.**  For every candidate the two
wall rates sum exactly to the event rate: the split is an identity of
the source-fixed fields, not an extra assumption. -/
theorem SourceClockCandidate.wall_split (c : SourceClockCandidate) :
    c.wallRate01 + c.wallRate12 = c.ratePerWindow := by
  rw [c.wall01_committed, c.wall12_committed, c.rate_committed,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun g _ => ?_
  rw [c.event_split g]
  push_cast
  ring

/-- **What the source fixes.**  Every candidate carries the committed
tick stream, its rate is the exact mean of its counting observable under
the committed empirical law, the wall rates sum to it, and it is
positive.  The one field the source does not fix is the identification
field. -/
theorem SourceClockCandidate.source_fixed_receipt (c : SourceClockCandidate) :
    c.tickCount = sourceWindowCount
      ∧ c.ratePerWindow = ∑ g, sourceTauEmpQ g * (c.eventCount g : ℚ)
      ∧ c.wallRate01 + c.wallRate12 = c.ratePerWindow
      ∧ 0 < c.ratePerWindow :=
  ⟨c.tick_committed, c.rate_committed, c.wall_split, c.rate_pos⟩

/-- The repair clock candidate: the committed retained-window sequence
as the tick stream, counting the repair moves of each window.  Its rate
invariants are the exact source-counted mean action `197/1754` per
window and the wall decomposition `94/1754 + 103/1754`.  The
identification field carries a stand-in schema, the bare statement that
some positive duration could be declared; a physical declaration would
replace it with the laboratory identification of one window, which no
theorem supplies. -/
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

/-- The tick stream of the repair candidate is the committed count. -/
@[simp] theorem repairClockCandidate_tickCount :
    repairClockCandidate.tickCount = sourceWindowCount := rfl

/-- The counting observable of the repair candidate is the committed
repair-move count. -/
@[simp] theorem repairClockCandidate_eventCount :
    repairClockCandidate.eventCount = sourceAction := rfl

/-- The exact repair rate invariant `197/1754` per window. -/
@[simp] theorem repairClockCandidate_ratePerWindow :
    repairClockCandidate.ratePerWindow = 197 / 1754 := rfl

/-- The exact `(0,1)`-wall rate invariant `94/1754` per window. -/
@[simp] theorem repairClockCandidate_wallRate01 :
    repairClockCandidate.wallRate01 = 94 / 1754 := rfl

/-- The exact `(1,2)`-wall rate invariant `103/1754` per window. -/
@[simp] theorem repairClockCandidate_wallRate12 :
    repairClockCandidate.wallRate12 = 103 / 1754 := rfl

/-- The tick stream totals the committed `1754` retained windows. -/
theorem repairClockCandidate_tick_total :
    ∑ g, repairClockCandidate.tickCount g = 1754 := by
  show ∑ g, sourceWindowCount g = 1754
  rw [Fin.sum_univ_eight]
  exact sourceWindowCount_total

/-- The count-weighted event total is the committed `197`: the rate
invariant is the exact ratio of two committed integers. -/
theorem repairClockCandidate_event_total :
    ∑ g, repairClockCandidate.tickCount g * repairClockCandidate.eventCount g
      = 197 := by
  show ∑ g, sourceWindowCount g * sourceAction g = 197
  rw [Fin.sum_univ_eight]
  exact sourceWindowCount_action_total

/-- The wall decomposition of the repair rate, as the structure identity
and as the exact rational identity `94/1754 + 103/1754 = 197/1754`. -/
theorem repairClockCandidate_wall_decomposition :
    repairClockCandidate.wallRate01 + repairClockCandidate.wallRate12
        = repairClockCandidate.ratePerWindow
      ∧ (94 / 1754 : ℚ) + 103 / 1754 = 197 / 1754 :=
  ⟨repairClockCandidate.wall_split, by norm_num⟩

/-- **The rate invariants are committed source expectations.**  The
repair rate per window is the expectation of the repair Hamiltonian in
the committed source state, and the wall rates are the expectations of
the two Ising domain-wall bonds. -/
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
while agreeing with the repair candidate on the tick stream and the rate
invariant: the structure fixes rates and never asserts, selects, or
discharges the identification. -/
theorem identification_field_free (p : Prop) :
    ∃ c : SourceClockCandidate,
      c.windowIdentification = p
        ∧ c.tickCount = repairClockCandidate.tickCount
        ∧ c.ratePerWindow = repairClockCandidate.ratePerWindow :=
  ⟨{ repairClockCandidate with windowIdentification := p }, rfl, rfl, rfl⟩

/-! ## (2) Composition with the calibration import -/

/-- A declared window duration makes the candidate a `ClockCalibration`:
the composed dictionary carries the declared duration as its tick.  The
candidate argument is deliberately unused in the tick field: the
duration is entirely the declaration's contribution, while the candidate
contributes the committed tick stream and its exact per-window rate
invariants, none of which is a duration. -/
def SourceClockCandidate.toClockCalibration (_c : SourceClockCandidate)
    (tau_w : ℝ) (h : 0 < tau_w) : ClockCalibration where
  tau := tau_w
  tau_pos := h

/-- The composed calibration ticks at the declared duration. -/
@[simp] theorem toClockCalibration_tau (c : SourceClockCandidate)
    (tau_w : ℝ) (h : 0 < tau_w) :
    (c.toClockCalibration tau_w h).tau = tau_w := rfl

/-- The laboratory event rate of a candidate under a declared
calibration whose tick is read as the window duration:
`ratePerWindow / tau` events per second, exact-symbolic and conditional
on the declaration. -/
def SourceClockCandidate.labRate (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.ratePerWindow : ℝ) / cal.tau

/-- The laboratory `(0,1)`-wall rate under a declared calibration. -/
def SourceClockCandidate.labWallRate01 (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.wallRate01 : ℝ) / cal.tau

/-- The laboratory `(1,2)`-wall rate under a declared calibration. -/
def SourceClockCandidate.labWallRate12 (c : SourceClockCandidate)
    (cal : ClockCalibration) : ℝ :=
  (c.wallRate12 : ℝ) / cal.tau

/-- **Laboratory attachment of the repair rate.**  Under any declared
calibration the repair rate attaches to the exact laboratory rate
`(197/1754) / tau` events per second. -/
theorem repairClockCandidate_labRate (cal : ClockCalibration) :
    repairClockCandidate.labRate cal = 197 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labRate, repairClockCandidate_ratePerWindow]
  norm_num

/-- The laboratory attachment through the composed calibration: under a
declared window duration `tau_w` the repair rate is exactly
`(197/1754) / tau_w` events per second. -/
theorem repairClockCandidate_labRate_declared (tau_w : ℝ) (h : 0 < tau_w) :
    repairClockCandidate.labRate
        (repairClockCandidate.toClockCalibration tau_w h)
      = 197 / 1754 / tau_w :=
  repairClockCandidate_labRate _

/-- The `(0,1)`-wall rate attaches to `(94/1754) / tau` events per
second. -/
theorem repairClockCandidate_labWallRate01 (cal : ClockCalibration) :
    repairClockCandidate.labWallRate01 cal = 94 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labWallRate01,
    repairClockCandidate_wallRate01]
  norm_num

/-- The `(1,2)`-wall rate attaches to `(103/1754) / tau` events per
second. -/
theorem repairClockCandidate_labWallRate12 (cal : ClockCalibration) :
    repairClockCandidate.labWallRate12 cal = 103 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labWallRate12,
    repairClockCandidate_wallRate12]
  norm_num

/-- **Exact laboratory sum identity.**  For every candidate and every
declared calibration the two laboratory wall rates sum exactly to the
laboratory event rate. -/
theorem SourceClockCandidate.labWall_split (c : SourceClockCandidate)
    (cal : ClockCalibration) :
    c.labWallRate01 cal + c.labWallRate12 cal = c.labRate cal := by
  simp only [SourceClockCandidate.labWallRate01,
    SourceClockCandidate.labWallRate12, SourceClockCandidate.labRate]
  rw [← add_div, ← Rat.cast_add, c.wall_split]

/-- **Composition with the frequency dictionary.**  The laboratory
counting rate is the import's laboratory frequency of the internal
angular rate `2 * pi * ratePerWindow`: the `2 * pi` factor records that
`labFrequency` converts angular rates while the candidate counts
events. -/
theorem SourceClockCandidate.labRate_eq_labFrequency
    (c : SourceClockCandidate) (cal : ClockCalibration) :
    c.labRate cal
      = cal.labFrequency (2 * Real.pi * (c.ratePerWindow : ℝ)) := by
  have h2pi : (2 * Real.pi : ℝ) ≠ 0 :=
    mul_ne_zero two_ne_zero Real.pi_ne_zero
  simp only [SourceClockCandidate.labRate,
    ClockCalibration.labFrequency_def]
  rw [mul_div_mul_left _ _ h2pi]

/-- **Dimensionful cancellation.**  The laboratory rate times the
laboratory duration of `n` windows is the exact expected event count
`n * ratePerWindow`, independent of the declared duration: events per
second times seconds cancels as exact real algebra. -/
theorem SourceClockCandidate.labRate_mul_labSeconds
    (c : SourceClockCandidate) (cal : ClockCalibration) (n : ℝ) :
    c.labRate cal * cal.labSeconds n = n * (c.ratePerWindow : ℝ) := by
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  simp only [SourceClockCandidate.labRate, ClockCalibration.labSeconds_def]
  calc (c.ratePerWindow : ℝ) / cal.tau * (n * cal.tau)
      = n * ((c.ratePerWindow : ℝ) * (cal.tau / cal.tau)) := by ring
    _ = n * (c.ratePerWindow : ℝ) := by rw [div_self hτ, mul_one]

/-! ## (3) Reduction to the single identification -/

/-- Two clock calibrations with the same tick are equal: the calibration
dictionary carries no data beyond its tick. -/
theorem calibration_eq_of_tau_eq {cal cal' : ClockCalibration}
    (h : cal.tau = cal'.tau) : cal = cal' := by
  obtain ⟨t, ht⟩ := cal
  obtain ⟨t', ht'⟩ := cal'
  have h' : t = t' := h
  subst h'
  rfl

/-- **Reduction to the single identification.**  Two candidate
calibrations agreeing on the declared duration agree everywhere: the
composed calibrations are equal, and every derived conversion, duration
and frequency alike, coincides.  The composed dictionary is fully
determined by the one declared duration. -/
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

/-- **Non-forcing.**  The identification is an import, never derived:
for any candidate, two calibrations with distinct ticks give distinct
laboratory rates.  No internal statement fixes the laboratory rate of
the committed process without a declared duration. -/
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

/-- **Rates per window, not per second.**  The candidate's per-window
rate invariant is declaration-invariant, while the composed laboratory
rate changes with the declared duration: under two distinct ticks the
per-window rate agrees and the laboratory rates differ.  The packet
fixes rates per window; no theorem produces the duration. -/
theorem rate_is_per_window_not_per_second (c : SourceClockCandidate)
    (cal cal' : ClockCalibration) (hne : cal.tau ≠ cal'.tau) :
    c.ratePerWindow = c.ratePerWindow
      ∧ c.labRate cal ≠ c.labRate cal' :=
  ⟨rfl, c.labRate_not_forced cal cal' hne⟩

/-- The wall-only subprocess as a second candidate: the same committed
tick stream, counting only the `(0,1)`-wall events, with the different
exact rate invariant `94/1754` per window.  The identification field
carries the same stand-in schema as the repair candidate. -/
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

/-- The tick stream of the wall-only candidate is the same committed
count. -/
@[simp] theorem wallOnlyClockCandidate_tickCount :
    wallOnlyClockCandidate.tickCount = sourceWindowCount := rfl

/-- The counting observable of the wall-only candidate is the
`(0,1)`-wall count. -/
@[simp] theorem wallOnlyClockCandidate_eventCount :
    wallOnlyClockCandidate.eventCount = bond01Energy := rfl

/-- The exact wall-only rate invariant `94/1754` per window. -/
@[simp] theorem wallOnlyClockCandidate_ratePerWindow :
    wallOnlyClockCandidate.ratePerWindow = 94 / 1754 := rfl

/-- The `(0,1)`-wall rate of the wall-only candidate. -/
@[simp] theorem wallOnlyClockCandidate_wallRate01 :
    wallOnlyClockCandidate.wallRate01 = 94 / 1754 := rfl

/-- The `(1,2)`-wall rate of the wall-only candidate vanishes. -/
@[simp] theorem wallOnlyClockCandidate_wallRate12 :
    wallOnlyClockCandidate.wallRate12 = 0 := rfl

/-- The wall-only rate attaches to `(94/1754) / tau` events per
second. -/
theorem wallOnlyClockCandidate_labRate (cal : ClockCalibration) :
    wallOnlyClockCandidate.labRate cal = 94 / 1754 / cal.tau := by
  simp only [SourceClockCandidate.labRate,
    wallOnlyClockCandidate_ratePerWindow]
  norm_num

/-- **Candidate selection is a declared choice.**  The wall-only
subprocess is a second source-produced candidate on the same committed
tick stream with a different exact rate invariant: `94/1754` against
`197/1754`.  The packet does not select the counting observable. -/
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

/-- **Source clock candidate receipt.**  The bundled statement: the
repair candidate carries the committed tick stream with total `1754`
windows and `197` counted repair events, the exact rate invariant
`197/1754` per window with its wall split, the laboratory attachments
`(197/1754) / tau` with the exact wall sum identity under every
declared calibration, the reduction of the composed dictionary to the
single declared duration, the non-forcing of that duration, and the
wall-only subprocess as a second candidate with a different rate
invariant on the same committed tick stream. -/
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
      (reduction_to_single_identification c c' t t' h h' he).1,
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
