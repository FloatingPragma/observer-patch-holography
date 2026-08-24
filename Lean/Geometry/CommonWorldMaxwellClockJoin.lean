import Geometry.CommonWorldInstrumentJoin
import CertifiedScaledStepInstrument

set_option autoImplicit false

open scoped BigOperators Matrix

/-!
# CW1 Maxwell-clock step join on the instrumented record (issue #740)

WHAT IS PROVED.  The instrumented common-world record
`InstrumentedCommonWorldArchitecture` of
`Geometry/CommonWorldInstrumentJoin.lean` carries its certified Maxwell
bundle and its clock worldline on one record, linked by the shared step
index type `ℕ` and by no clause: the boundary paragraph of that module
records that no clause relates a Maxwell quantity at step `n` to the clock
or the worldline event at step `n`.  This module extends that record by one
join field and proves four clauses of every inhabitant.

The join field.  `MaxwellClockJoinedArchitecture` extends the instrumented
record by a single function `join : ℕ → ℝ × Herm2` whose first component is
the scaled staggered form `fieldEnergyScaled scaled.h scaled.A scaled.phi n`
of the record's certified bundle and whose second component is the clock
worldline event `frameWorldline frame (stepTime stepDuration n)` of the
record's own frame at the record's declared step time.  One function of one
index produces both components, so the two islands are read at the same
index by construction rather than by convention.

Clause 1, same-index binding (`CarriesJoinedStepIndex`).  The two
projection identities hold (`join_energy`, `join_stepEvent`), the join is
the pair of the two committed readings (`join_eq_pair`), and the join is
injective in the step index (`join_injective`), so the index is recovered
from the joined object.  The binding is load-bearing: at every index the
joined event differs from the event of the successor index
(`join_event_ne_succ`), so no inhabitant's join is the successor-shifted
reading (`shifted_join_excluded`).

Clause 2, causal-energy compatibility (`CarriesJoinedCausalEnergy`).  Along
the join, successive events are causally ordered in the record's cone order
and the ordering extends to every pair of indices `n ≤ m`
(`joined_causal_mono`, from the landed `stepEvent_causal_mono`); the joined
value is nonnegative at every index (`joined_energy_nonneg`); the record's
step satisfies the sharp certificate `h² (3 + √5) < 4`
(`joined_step_certified`, from the landed eigenvalue bound
`scaledBundle_courant_ge_golden`), so the record's step is the step of a
`CertifiedStepInstrument` of `Screen/CertifiedScaledStepInstrument.lean`
(`certifiedStepOfJoin`); and for zero seam current the electric seam energy
and the magnetic face energy at index `n` are bounded by `8` and `16` times
the joined value at index `0` over the record's margin `4 - h² Λ`
(`joined_causal_energy_bound`), which states the landed uniform bound with
the joined object on the right.

Clause 3, conservation transport (`CarriesJoinedConservation`).  The clock
phase between two step indices advances by the declared rotation of the
elapsed step time (`joined_clock_advance`), and for zero seam current the
joined value at the later index of every causally ordered pair equals the
joined value at the earlier index (`joined_conservation`).  This is the
first record-level clause relating the screen island and the kinematics
island of the record: one conserved screen quantity transported along the
record's own clock.  On the committed inhabitant the two sides are not both
constant: the joined value is constant at every index while the clock does
not repeat with period one step
(`committedJoinedWitness_energy_constant_clock_moving`, through
`unit_phase_not_multiple_of_two_pi` and the landed
`stepClock_periodic_iff`).

Clause 4, non-supply (`JoinSuppliesNoCarrierMap`).  The join's event
component is a function of the record's frame and step duration alone and
its value component is a function of the record's certified bundle alone:
two inhabitants agreeing on frame and step duration have the same event
component at every index, and two inhabitants agreeing on the certified
bundle have the same value component at every index
(`joinSuppliesNoCarrierMap`).  No field of the join takes a port, a seam,
or a face to a point, an interval, or a cone of the Lorentz module, and the
clause is the form in which the audit checks that against the structure.

INHABITANTS.  `joinOf` produces a joined record over any instrumented
record.  `committedJoinedWitness` is `joinOf instrumentedCommittedWitness`
and restricts to the committed instrumented inhabitant definitionally
(`committedJoinedWitness_extends`).  `fourFifthsJoinedWitness` is the same
record with the certified bundle `fourFifthsBundle` at the certified step
`4/5` and the sharp carrier constant `3 + √5`; its step is the step of
`fourFifthsInstrument` (`fourFifthsJoined_certifiedStep`) and its joined
value carries the rational certified bounds `‖E n‖² ≤ (25/2) (join 0).1`
and `‖B n‖² ≤ 25 (join 0).1` (`fourFifthsJoined_rational_bound`).

DECLARED DATA.  Every field of the instrumented record stays declared as
recorded there, the step duration included.  The join field is constrained
by the two projection identities and adds no further datum:
`joined_stepDuration_not_forced` exhibits two joined inhabitants over one
bridged base with one certified bundle and one value component whose step
durations differ and whose joined events at index one differ, so the join
does not fix the clock calibration (PR-15 open).

FALSIFIER.  The module fails if some inhabitant's join misses a projection
identity, if some join is noninjective in the index, if some inhabitant
admits the successor-shifted event reading, if some pair of joined events
of ordered indices is causally unordered, if some record's step misses the
sharp certificate, if the zero-current bound or the conservation transport
fails at some index, if the two non-supply agreements fail, or if the
extended record type is uninhabited over the committed inhabitant.

WHAT IS NOT PROVED HERE.  The step dictionary stays a declared calibration:
`stepDuration` is a positive real with no unit, the join identifies no
physical time and no physical space, and no clock is calibrated (PR-15
open).  After this module the common-world joins that stay missing are
exactly four: no common action across carriers, that is no one variational
principle carrying the screen action of the scaled Maxwell island and the
kinematic flow of the Lorentz module on one carrier (PR-54 open); no port,
seam, or face of the screen carrier mapped to a point, interval, or cone of
the Lorentz module (PR-53 open), which clause 4 states as a non-supply of
the join itself; no matter dynamics beyond the committed structure
premises; and no observer readout, the instrument carrying channel clauses
and a repeatability proxy with no public outcome, readback, or provenance
(PR-64 readback and PR-65 open).  The joined value is the staggered form of
the record's declared bundle and is not an energy in any unit; the
conservation clause is conditional on the record's zero seam current.
Nothing here closes issue #740, and no observation-ledger row is promoted.

Axiom audit.  Every proof composes committed receipts with exact
mathematics; the module adds no project axiom and uses no native decision
procedure.  The guard lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

namespace OPH.CommonWorldMaxwellClockJoin

open OPH.C1Lorentz OPH.CommonWorld OPH.CommonWorldIslandBridge
open OPH.CausalComposition
open OPH.TemporalMaxwellEvolution OPH.ScaledMaxwellStability
open OPH.LocalFaceMaxwellAction OPH.DiscreteCoulombGreen
open OPH.CommonWorldInstrumentJoin OPH.CertifiedScaledStepInstrument

noncomputable section

/-! ## The joined step record -/

/-- The Maxwell-clock joined architecture: the instrumented common-world
record extended by one join field.  The field is a single function of a
single step index producing the pair of the scaled staggered-form value of
the record's certified bundle and the worldline event of the record's clock
at the record's declared step time, so both islands are read at the same
index by construction.  Every field of the instrumented record stays a
declared hypothesis as recorded there. -/
structure MaxwellClockJoinedArchitecture extends
    InstrumentedCommonWorldArchitecture where
  /-- The joined step: one pair per step index. -/
  join : ℕ → ℝ × Herm2
  /-- The value component reads the record's certified bundle at the index. -/
  join_energy : ∀ n : ℕ,
    (join n).1 = fieldEnergyScaled scaled.h scaled.A scaled.phi n
  /-- The event component reads the record's clock worldline at the step
  time of the same index. -/
  join_event : ∀ n : ℕ,
    (join n).2 = frameWorldline frame (stepTime stepDuration n)

/-- The joined record over any instrumented record, with the two committed
readings as the join. -/
def joinOf (W : InstrumentedCommonWorldArchitecture) :
    MaxwellClockJoinedArchitecture where
  toInstrumentedCommonWorldArchitecture := W
  join := fun n =>
    (fieldEnergyScaled W.scaled.h W.scaled.A W.scaled.phi n, stepEvent W n)
  join_energy := fun _ => rfl
  join_event := fun _ => rfl

/-- The joined record over `W` restricts to `W` definitionally. -/
theorem joinOf_extends (W : InstrumentedCommonWorldArchitecture) :
    (joinOf W).toInstrumentedCommonWorldArchitecture = W := rfl

/-! ## Clause 1: same-index binding -/

/-- The event component of the join at index `n` is the record's step event
at index `n`. -/
theorem join_stepEvent (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    (J.join n).2 = stepEvent J.toInstrumentedCommonWorldArchitecture n :=
  J.join_event n

/-- The join at index `n` is the pair of the two committed readings at that
one index. -/
theorem join_eq_pair (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    J.join n =
      (fieldEnergyScaled J.scaled.h J.scaled.A J.scaled.phi n,
        stepEvent J.toInstrumentedCommonWorldArchitecture n) := by
  rw [← J.join_energy n, ← join_stepEvent J n]

/-- The join is injective in the step index: the index is recovered from
the joined object through its event component. -/
theorem join_injective (J : MaxwellClockJoinedArchitecture) :
    Function.Injective J.join := by
  intro n m h
  have h2 : (J.join n).2 = (J.join m).2 := by rw [h]
  rw [join_stepEvent, join_stepEvent] at h2
  exact stepEvent_injective J.toInstrumentedCommonWorldArchitecture h2

/-- The join at index `n` carries the event of index `n` and not the event
of the successor index. -/
theorem join_event_ne_succ (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    (J.join n).2 ≠ stepEvent J.toInstrumentedCommonWorldArchitecture (n + 1) := by
  rw [join_stepEvent]
  intro h
  have h1 := stepEvent_injective J.toInstrumentedCommonWorldArchitecture h
  omega

/-- The same-index binding is load-bearing: no inhabitant's join reads the
event of the successor index into the pair of index `n`. -/
theorem shifted_join_excluded (J : MaxwellClockJoinedArchitecture)
    (h : ∀ n : ℕ,
      (J.join n).2 = stepEvent J.toInstrumentedCommonWorldArchitecture (n + 1)) :
    False :=
  join_event_ne_succ J 0 (h 0)

/-- The same-index clause of the composed receipt. -/
def CarriesJoinedStepIndex (J : MaxwellClockJoinedArchitecture) : Prop :=
  (∀ n : ℕ, J.join n =
    (fieldEnergyScaled J.scaled.h J.scaled.A J.scaled.phi n,
      stepEvent J.toInstrumentedCommonWorldArchitecture n)) ∧
  (∀ n : ℕ, (J.join n).1 =
    fieldEnergyScaled J.scaled.h J.scaled.A J.scaled.phi n) ∧
  (∀ n : ℕ, (J.join n).2 = stepEvent J.toInstrumentedCommonWorldArchitecture n) ∧
  Function.Injective J.join

theorem carriesJoinedStepIndex (J : MaxwellClockJoinedArchitecture) :
    CarriesJoinedStepIndex J :=
  ⟨join_eq_pair J, J.join_energy, join_stepEvent J, join_injective J⟩

/-! ## Clause 2: causal-energy compatibility along the join -/

/-- Joined events of ordered indices are causally ordered in the record's
cone order. -/
theorem joined_causal_mono (J : MaxwellClockJoinedArchitecture) {n m : ℕ}
    (h : n ≤ m) : causalLE (J.join n).2 (J.join m).2 := by
  rw [join_stepEvent, join_stepEvent]
  exact stepEvent_causal_mono J.toInstrumentedCommonWorldArchitecture h

/-- The joined value is nonnegative at every index. -/
theorem joined_energy_nonneg (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    0 ≤ (J.join n).1 := by
  rw [J.join_energy]
  exact fieldEnergyScaled_nonneg J.scaled.h J.scaled.h_pos.ne' J.scaled.Λ
    J.scaled.courant J.scaled.courant_strict.le J.scaled.A J.scaled.phi n

/-- The record's step satisfies the sharp certificate `h² (3 + √5) < 4`:
the record's Courant constant dominates the sharp eigenvalue. -/
theorem joined_step_certified (J : MaxwellClockJoinedArchitecture) :
    J.scaled.h ^ 2 * (3 + Real.sqrt 5) < 4 := by
  have hgold := scaledBundle_courant_ge_golden J.scaled
  have hstrict := J.scaled.courant_strict
  have hsq : 0 ≤ J.scaled.h ^ 2 := sq_nonneg _
  nlinarith

/-- The record's step is the step of a certified step instrument of
`Screen/CertifiedScaledStepInstrument.lean`. -/
def certifiedStepOfJoin (J : MaxwellClockJoinedArchitecture) :
    CertifiedStepInstrument :=
  instrumentOfStrictStep J.scaled.h J.scaled.h_pos (joined_step_certified J)

theorem certifiedStepOfJoin_step (J : MaxwellClockJoinedArchitecture) :
    (certifiedStepOfJoin J).step = J.scaled.h := rfl

/-- The zero-current scaled evolution law read off the record: the
certified bundle and the record's unit-step bundle carry one seam
current. -/
theorem joined_scaled_ampere_zero (J : MaxwellClockJoinedArchitecture)
    (hJ : J.maxwell.J = (fun _ ↦ 0)) :
    AmpereEvolutionScaled J.scaled.h J.scaled.A J.scaled.phi (fun _ ↦ 0) := by
  have h := J.scaled.ampere
  rw [J.scaled_current, hJ] at h
  exact h

/-- **Causal-energy compatibility.**  For zero seam current, at every index
the successive joined events are causally ordered while the electric seam
energy and the magnetic face energy stay below the certified multiples of
the joined value at index `0`. -/
theorem joined_causal_energy_bound (J : MaxwellClockJoinedArchitecture)
    (hJ : J.maxwell.J = (fun _ ↦ 0)) (n : ℕ) :
    causalLE (J.join n).2 (J.join (n + 1)).2 ∧
      realSeamEnergy (electricFieldScaled J.scaled.h J.scaled.A J.scaled.phi n) ≤
        8 * (J.join 0).1 / (4 - J.scaled.h ^ 2 * J.scaled.Λ) ∧
      faceEnergy (magneticField J.scaled.A n) ≤
        16 * (J.join 0).1 / (4 - J.scaled.h ^ 2 * J.scaled.Λ) := by
  obtain ⟨hE, hB⟩ := stability_certificate J.scaled.h J.scaled.h_pos.ne' J.scaled.Λ
    J.scaled.Λ_nonneg J.scaled.courant J.scaled.courant_strict J.scaled.A
    J.scaled.phi (joined_scaled_ampere_zero J hJ) n
  rw [J.join_energy 0]
  exact ⟨joined_causal_mono J (Nat.le_succ n), hE, hB⟩

/-- The causal-energy clause of the composed receipt. -/
def CarriesJoinedCausalEnergy (J : MaxwellClockJoinedArchitecture) : Prop :=
  (∀ n m : ℕ, n ≤ m → causalLE (J.join n).2 (J.join m).2) ∧
  (∀ n : ℕ, 0 ≤ (J.join n).1) ∧
  J.scaled.h ^ 2 * (3 + Real.sqrt 5) < 4 ∧
  (certifiedStepOfJoin J).step = J.scaled.h ∧
  (J.maxwell.J = (fun _ ↦ 0) → ∀ n : ℕ,
    causalLE (J.join n).2 (J.join (n + 1)).2 ∧
      realSeamEnergy (electricFieldScaled J.scaled.h J.scaled.A J.scaled.phi n) ≤
        8 * (J.join 0).1 / (4 - J.scaled.h ^ 2 * J.scaled.Λ) ∧
      faceEnergy (magneticField J.scaled.A n) ≤
        16 * (J.join 0).1 / (4 - J.scaled.h ^ 2 * J.scaled.Λ))

theorem carriesJoinedCausalEnergy (J : MaxwellClockJoinedArchitecture) :
    CarriesJoinedCausalEnergy J :=
  ⟨fun _ _ h => joined_causal_mono J h, joined_energy_nonneg J,
    joined_step_certified J, certifiedStepOfJoin_step J,
    fun hJ n => joined_causal_energy_bound J hJ n⟩

/-! ## Clause 3: conservation transport along the clock -/

/-- The clock phase between two step indices advances by the declared
rotation of the elapsed step time. -/
theorem joined_clock_advance (J : MaxwellClockJoinedArchitecture) (n m : ℕ) :
    stepClock J.toInstrumentedCommonWorldArchitecture m =
      stepClock J.toInstrumentedCommonWorldArchitecture n *
        Complex.exp (-Complex.I *
          ((J.mass *
            (stepTime J.stepDuration m - stepTime J.stepDuration n) : ℝ) : ℂ)) := by
  rw [stepClock_eq_exp, stepClock_eq_exp, ← Complex.exp_add]
  congr 1
  push_cast
  ring

/-- **Conservation transport.**  For zero seam current, along every
causally ordered pair of joined events the joined value at the later index
equals the joined value at the earlier index. -/
theorem joined_conservation (J : MaxwellClockJoinedArchitecture)
    (hJ : J.maxwell.J = (fun _ ↦ 0)) {n m : ℕ} (h : n ≤ m) :
    causalLE (J.join n).2 (J.join m).2 ∧ (J.join m).1 = (J.join n).1 := by
  refine ⟨joined_causal_mono J h, ?_⟩
  rw [J.join_energy m, J.join_energy n]
  rw [energy_conserved_scaled J.scaled.h J.scaled.h_pos.ne' J.scaled.A J.scaled.phi
      (joined_scaled_ampere_zero J hJ) m,
    energy_conserved_scaled J.scaled.h J.scaled.h_pos.ne' J.scaled.A J.scaled.phi
      (joined_scaled_ampere_zero J hJ) n]

/-- The conservation-transport clause of the composed receipt: the clock
phase advances by the declared rotation between any two indices, and for
zero seam current the joined value is unchanged along every causally
ordered pair. -/
def CarriesJoinedConservation (J : MaxwellClockJoinedArchitecture) : Prop :=
  (∀ n m : ℕ, stepClock J.toInstrumentedCommonWorldArchitecture m =
    stepClock J.toInstrumentedCommonWorldArchitecture n *
      Complex.exp (-Complex.I *
        ((J.mass *
          (stepTime J.stepDuration m - stepTime J.stepDuration n) : ℝ) : ℂ))) ∧
  (J.maxwell.J = (fun _ ↦ 0) → ∀ n m : ℕ, n ≤ m →
    causalLE (J.join n).2 (J.join m).2 ∧ (J.join m).1 = (J.join n).1)

theorem carriesJoinedConservation (J : MaxwellClockJoinedArchitecture) :
    CarriesJoinedConservation J :=
  ⟨joined_clock_advance J, fun hJ _ _ h => joined_conservation J hJ h⟩

/-! ## Clause 4: the join supplies no carrier map -/

/-- The non-supply clause: the join's event component is a function of the
record's frame and step duration alone, and its value component is a
function of the record's certified bundle alone.  No field of the join
takes a port, a seam, or a face to a point, an interval, or a cone of the
Lorentz module; the two agreement statements are the form in which that is
checked against the structure. -/
def JoinSuppliesNoCarrierMap (J : MaxwellClockJoinedArchitecture) : Prop :=
  (∀ n : ℕ, (J.join n).2 = frameWorldline J.frame (stepTime J.stepDuration n)) ∧
  (∀ K : MaxwellClockJoinedArchitecture, K.frame = J.frame →
    K.stepDuration = J.stepDuration → ∀ n : ℕ, (K.join n).2 = (J.join n).2) ∧
  (∀ K : MaxwellClockJoinedArchitecture, K.scaled = J.scaled →
    ∀ n : ℕ, (K.join n).1 = (J.join n).1)

theorem joinSuppliesNoCarrierMap (J : MaxwellClockJoinedArchitecture) :
    JoinSuppliesNoCarrierMap J := by
  refine ⟨J.join_event, fun K hf hδ n => ?_, fun K hs n => ?_⟩
  · rw [K.join_event, J.join_event, hf, hδ]
  · rw [K.join_energy, J.join_energy, hs]

/-! ## The composed receipt -/

/-- **The Maxwell-clock joined receipt (issue #740, one added join; CW1
stays open).**  Every joined record carries the four join clauses of the
instrumented record and the four clauses of the step join: same-index
binding, causal-energy compatibility, conservation transport along the
clock, and the non-supply of a carrier map.  The ten base clauses of the
bridged receipt transport unchanged through
`instrumentedCommonWorld_receipt` at
`toInstrumentedCommonWorldArchitecture`. -/
theorem maxwellClockJoined_receipt (J : MaxwellClockJoinedArchitecture) :
    CarriesCertifiedEvolution J.toInstrumentedCommonWorldArchitecture ∧
      CarriesInstrumentJoin J.toInstrumentedCommonWorldArchitecture ∧
      CarriesStepDictionary J.toInstrumentedCommonWorldArchitecture ∧
      CarriesFirstLawJoin J.toInstrumentedCommonWorldArchitecture ∧
      CarriesJoinedStepIndex J ∧
      CarriesJoinedCausalEnergy J ∧
      CarriesJoinedConservation J ∧
      JoinSuppliesNoCarrierMap J :=
  ⟨carriesCertifiedEvolution _, carriesInstrumentJoin _, carriesStepDictionary _,
    carriesFirstLawJoin _, carriesJoinedStepIndex J, carriesJoinedCausalEnergy J,
    carriesJoinedConservation J, joinSuppliesNoCarrierMap J⟩

/-! ## The committed inhabitant -/

/-- The committed joined inhabitant: the committed instrumented inhabitant
with the join of its own two readings. -/
def committedJoinedWitness : MaxwellClockJoinedArchitecture :=
  joinOf instrumentedCommittedWitness

/-- The committed joined inhabitant restricts to the committed instrumented
inhabitant definitionally. -/
theorem committedJoinedWitness_extends :
    committedJoinedWitness.toInstrumentedCommonWorldArchitecture =
      instrumentedCommittedWitness := rfl

/-- The extended record is inhabited. -/
theorem maxwellClockJoined_inhabited :
    Nonempty MaxwellClockJoinedArchitecture := ⟨committedJoinedWitness⟩

/-- The composed receipt holds of the committed joined inhabitant. -/
theorem committedJoinedWitness_receipt :
    CarriesCertifiedEvolution
        committedJoinedWitness.toInstrumentedCommonWorldArchitecture ∧
      CarriesInstrumentJoin
        committedJoinedWitness.toInstrumentedCommonWorldArchitecture ∧
      CarriesStepDictionary
        committedJoinedWitness.toInstrumentedCommonWorldArchitecture ∧
      CarriesFirstLawJoin
        committedJoinedWitness.toInstrumentedCommonWorldArchitecture ∧
      CarriesJoinedStepIndex committedJoinedWitness ∧
      CarriesJoinedCausalEnergy committedJoinedWitness ∧
      CarriesJoinedConservation committedJoinedWitness ∧
      JoinSuppliesNoCarrierMap committedJoinedWitness :=
  maxwellClockJoined_receipt committedJoinedWitness

/-- The unit phase is no integer multiple of `2π`: from `π > 3`, a nonzero
multiple has absolute value at least `6` and the zero multiple is `0`. -/
theorem unit_phase_not_multiple_of_two_pi :
    ¬ ∃ m : ℤ, (1 : ℝ) = m * (2 * Real.pi) := by
  rintro ⟨m, hm⟩
  have hpi : (3 : ℝ) < Real.pi := Real.pi_gt_three
  rcases le_or_gt (m : ℝ) 0 with h | h
  · nlinarith
  · have hz : (0 : ℤ) < m := by exact_mod_cast h
    have h1 : (1 : ℤ) ≤ m := by omega
    have h1' : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast h1
    nlinarith

/-- The committed inhabitant's clock does not repeat with period one step:
its mass and step duration give the phase `1`, which is no integer multiple
of `2π`. -/
theorem committedJoinedWitness_clock_not_unit_periodic :
    ¬ ∀ n : ℕ,
      stepClock committedJoinedWitness.toInstrumentedCommonWorldArchitecture
          (n + 1) =
        stepClock committedJoinedWitness.toInstrumentedCommonWorldArchitecture n := by
  intro hper
  obtain ⟨m, hm⟩ :=
    (stepClock_periodic_iff committedJoinedWitness.toInstrumentedCommonWorldArchitecture
      1).mp hper
  have hphase :
      committedJoinedWitness.mass *
          stepTime committedJoinedWitness.stepDuration 1 = 1 := by
    show (1 : ℝ) * stepTime (1 : ℝ) 1 = 1
    simp [stepTime]
  rw [hphase] at hm
  exact unit_phase_not_multiple_of_two_pi ⟨m, hm⟩

/-- **The two islands on one record.**  On the committed inhabitant the
joined value is constant at every index while the clock does not repeat
with period one step: one conserved screen quantity transported along a
nonrepeating clock. -/
theorem committedJoinedWitness_energy_constant_clock_moving :
    (∀ n : ℕ, (committedJoinedWitness.join n).1 =
      (committedJoinedWitness.join 0).1) ∧
    ¬ ∀ n : ℕ,
      stepClock committedJoinedWitness.toInstrumentedCommonWorldArchitecture
          (n + 1) =
        stepClock committedJoinedWitness.toInstrumentedCommonWorldArchitecture n :=
  ⟨fun n => (joined_conservation committedJoinedWitness rfl (Nat.zero_le n)).2,
    committedJoinedWitness_clock_not_unit_periodic⟩

/-! ## The certified-step inhabitant -/

/-- The committed instrumented record with the certified bundle at the step
`4/5` and the sharp carrier constant `3 + √5`.  The bridged base and every
other field are those of the committed instrumented inhabitant. -/
def fourFifthsWorld : InstrumentedCommonWorldArchitecture :=
  { instrumentedCommittedWitness with
    scaled := fourFifthsBundle
    scaled_initial_seam := by
      constructor
      · show demoScaledA (4 / 5) 0 = demoA 0
        rw [demoScaledA_zero]
        exact (show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]).symm
      · show demoScaledA (4 / 5) 1 = demoA 1
        rw [demoScaledA_one]
        exact (show demoA 1 = demoInitial from by simp only [demoA]).symm
    scaled_current := rfl }

/-- The certified-step joined inhabitant. -/
def fourFifthsJoinedWitness : MaxwellClockJoinedArchitecture :=
  joinOf fourFifthsWorld

/-- The certified-step inhabitant carries the same bridged base as the
committed inhabitant. -/
theorem fourFifthsJoinedWitness_base :
    fourFifthsJoinedWitness.toBridgedCommonWorldArchitecture =
      committedJoinedWitness.toBridgedCommonWorldArchitecture := rfl

/-- The certified-step inhabitant's step is the step of the committed
certified instrument `fourFifthsInstrument`. -/
theorem fourFifthsJoined_certifiedStep :
    (certifiedStepOfJoin fourFifthsJoinedWitness).step =
      fourFifthsInstrument.step := rfl

/-- The rational certified bounds on the certified-step inhabitant, stated
with the joined value at index `0` on the right: `‖E n‖² ≤ (25/2) (join 0).1`
and `‖B n‖² ≤ 25 (join 0).1`. -/
theorem fourFifthsJoined_rational_bound (n : ℕ) :
    realSeamEnergy (electricFieldScaled fourFifthsJoinedWitness.scaled.h
        fourFifthsJoinedWitness.scaled.A fourFifthsJoinedWitness.scaled.phi n) ≤
      (25 / 2) * (fourFifthsJoinedWitness.join 0).1 ∧
    faceEnergy (magneticField fourFifthsJoinedWitness.scaled.A n) ≤
      25 * (fourFifthsJoinedWitness.join 0).1 := by
  rw [fourFifthsJoinedWitness.join_energy 0]
  exact certifiedStep_energy_bounds_rational (demoScaledA (4 / 5)) (fun _ ↦ 0)
    (demoScaled_ampere (4 / 5) (by norm_num)) n

/-! ## The join does not fix the clock calibration -/

/-- Two joined records with one frame and distinct step durations carry
distinct joined events at index one. -/
theorem joined_event_ne_of_stepDuration_ne (J K : MaxwellClockJoinedArchitecture)
    (hf : K.frame = J.frame) (hδ : K.stepDuration ≠ J.stepDuration) :
    (K.join 1).2 ≠ (J.join 1).2 := by
  rw [K.join_event, J.join_event, hf]
  intro h
  have h1 : stepTime K.stepDuration 1 = stepTime J.stepDuration 1 :=
    frameWorldline_injective J.frame h
  rw [stepTime_one, stepTime_one] at h1
  exact hδ h1

/-- The committed joined inhabitant at step duration `2`. -/
def doubleStepJoinedWitness : MaxwellClockJoinedArchitecture :=
  joinOf { instrumentedCommittedWitness with
    stepDuration := 2
    stepDuration_pos := two_pos }

/-- **The join does not fix the clock calibration.**  Two joined
inhabitants share one bridged base, one certified bundle, and one value
component of the join, while their step durations differ and their joined
events at index one differ.  PR-15 is open: no committed theorem selects
the step duration. -/
theorem joined_stepDuration_not_forced :
    ∃ J₁ J₂ : MaxwellClockJoinedArchitecture,
      J₁.toBridgedCommonWorldArchitecture = J₂.toBridgedCommonWorldArchitecture ∧
        J₁.scaled = J₂.scaled ∧
        (∀ n : ℕ, (J₁.join n).1 = (J₂.join n).1) ∧
        J₁.stepDuration ≠ J₂.stepDuration ∧
        (J₁.join 1).2 ≠ (J₂.join 1).2 := by
  refine ⟨doubleStepJoinedWitness, committedJoinedWitness, rfl, rfl,
    fun _ => rfl, ?_, ?_⟩
  · show (2 : ℝ) ≠ 1
    norm_num
  · refine joined_event_ne_of_stepDuration_ne committedJoinedWitness
      doubleStepJoinedWitness rfl ?_
    show (2 : ℝ) ≠ 1
    norm_num

end

end OPH.CommonWorldMaxwellClockJoin

/- Axiom audit: committed receipts and exact mathematics only.  Expected
axioms per line: at most `propext`, `Classical.choice`, `Quot.sound`.  No
native decision procedure is used. -/

#print axioms OPH.CommonWorldMaxwellClockJoin.joinOf_extends
#print axioms OPH.CommonWorldMaxwellClockJoin.join_stepEvent
#print axioms OPH.CommonWorldMaxwellClockJoin.join_eq_pair
#print axioms OPH.CommonWorldMaxwellClockJoin.join_injective
#print axioms OPH.CommonWorldMaxwellClockJoin.join_event_ne_succ
#print axioms OPH.CommonWorldMaxwellClockJoin.shifted_join_excluded
#print axioms OPH.CommonWorldMaxwellClockJoin.carriesJoinedStepIndex
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_causal_mono
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_energy_nonneg
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_step_certified
#print axioms OPH.CommonWorldMaxwellClockJoin.certifiedStepOfJoin_step
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_scaled_ampere_zero
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_causal_energy_bound
#print axioms OPH.CommonWorldMaxwellClockJoin.carriesJoinedCausalEnergy
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_clock_advance
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_conservation
#print axioms OPH.CommonWorldMaxwellClockJoin.carriesJoinedConservation
#print axioms OPH.CommonWorldMaxwellClockJoin.joinSuppliesNoCarrierMap
#print axioms OPH.CommonWorldMaxwellClockJoin.maxwellClockJoined_receipt
#print axioms OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness_extends
#print axioms OPH.CommonWorldMaxwellClockJoin.maxwellClockJoined_inhabited
#print axioms OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness_receipt
#print axioms OPH.CommonWorldMaxwellClockJoin.unit_phase_not_multiple_of_two_pi
#print axioms OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness_clock_not_unit_periodic
#print axioms OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness_energy_constant_clock_moving
#print axioms OPH.CommonWorldMaxwellClockJoin.fourFifthsJoinedWitness_base
#print axioms OPH.CommonWorldMaxwellClockJoin.fourFifthsJoined_certifiedStep
#print axioms OPH.CommonWorldMaxwellClockJoin.fourFifthsJoined_rational_bound
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_event_ne_of_stepDuration_ne
#print axioms OPH.CommonWorldMaxwellClockJoin.joined_stepDuration_not_forced
