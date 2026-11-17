import Geometry.SourceDerivedSpacetimeCarrier
import ObserverPatchHolography.Provenance.HistoryCausalInvariance

/-!
# Executed native repairs in the authenticated source-carrier interface

Four actual endpoint repairs on seams (0,1),(2,3),(0,1),(0,1), for every
binary coin sequence, instantiate the existing executed semantic-log
constructor. Payloads are computed from the preceding register state by
floor/ceiling averaging. Both endpoints are version-stamped even when the
current is zero: every retained instrument invocation is a commit in this
declared complete execution log. This is instrumentation of classical
branches, not a claim that an intrinsic observer reads all source values.

The source midpoint chart supplies spatial readback; canonical ancestry
height supplies a dimensionless ordinal. The resulting forward placement
uses the registered rank-three quotient. No laboratory clock, influence
interpretation, continuous Lorentz symmetry or continuum is inferred.
-/

namespace OPH.NativeRepairSemanticPlacement

open OPH.Provenance OPH.SourceDerivedSpacetime
  OPH.PrimitivePortFrameQuotient

abbrev Port := Fin 12
abbrev Event := Fin 4
abbrev Coins := Event → Bool

/-- The declared binary total-two preparation. -/
def start : VersionedState Port ℕ Event :=
  ⟨fun p ↦ if p = 0 ∨ p = 2 then 1 else 0, fun _ ↦ none⟩

/-- Fixed external schedule; no state-dependent scheduling reads are hidden. -/
def endpoints (e : Event) : Port × Port :=
  if e = 1 then (2, 3) else (0, 1)

/-- The native low/high endpoint payload at one classical branch. -/
def repaired (q : VersionedState Port ℕ Event) (e : Event) (coin : Bool) : Port → ℕ :=
  let (a, b) := endpoints e
  let s := q.value a + q.value b
  fun p ↦ if p = a then (if coin then (s+1)/2 else s/2)
    else if p = b then (if coin then s/2 else (s+1)/2) else q.value p

/-- Every accepted invocation reads and commits both endpoint registers. -/
def spec (q : VersionedState Port ℕ Event) (e : Event) (coin : Bool) :
    CommitSpec Port ℕ Event where
  eventId := e
  readSet := {(endpoints e).1, (endpoints e).2}
  writeSet := {(endpoints e).1, (endpoints e).2}
  causalSupp := {(endpoints e).1, (endpoints e).2}
  supp_subset_read := fun _ h ↦ h
  writeValue := repaired q e coin

def s0 (coins : Coins) := spec start 0 (coins 0)
def q1 (coins : Coins) := (s0 coins).apply start
def s1 (coins : Coins) := spec (q1 coins) 1 (coins 1)
def q2 (coins : Coins) := (s1 coins).apply (q1 coins)
def s2 (coins : Coins) := spec (q2 coins) 2 (coins 2)
def q3 (coins : Coins) := (s2 coins).apply (q2 coins)
def s3 (coins : Coins) := spec (q3 coins) 3 (coins 3)

def specs (coins : Coins) : List (CommitSpec Port ℕ Event) :=
  [s0 coins, s1 coins, s2 coins, s3 coins]

def coinsOfIndex (i : Fin 16) : Coins :=
  fun e ↦ i.val / 2^(3-e.val) % 2 == 1

def current (q : VersionedState Port ℕ Event) (e : Event) (coin : Bool) : ℤ :=
  (q.value (endpoints e).1 : ℤ) - (repaired q e coin (endpoints e).1 : ℤ)

/-- The admitted endpoint values determine the written payload and visible
current. No outside load is a hidden read dependency of this rule. -/
theorem endpoint_support_sound (q r : VersionedState Port ℕ Event)
    (e : Event) (coin : Bool)
    (ha : q.value (endpoints e).1 = r.value (endpoints e).1)
    (hb : q.value (endpoints e).2 = r.value (endpoints e).2) :
    (repaired q e coin (endpoints e).1, repaired q e coin (endpoints e).2,
        current q e coin) =
      (repaired r e coin (endpoints e).1, repaired r e coin (endpoints e).2,
        current r e coin) := by
  simp [repaired, current, ha, hb]

def branchRow (i : Fin 16) : List ℕ × List ℤ × List ℕ :=
  let c := coinsOfIndex i
  (List.ofFn (fun e ↦ if c e then 1 else 0),
    [current start 0 (c 0), current (q1 c) 1 (c 1),
     current (q2 c) 2 (c 2), current (q3 c) 3 (c 3)],
    List.ofFn (execState start (specs c)).value)

/-- Literal cross-language boundary: Python compares every row to its native
execution, and the following theorem proves every row from these Lean specs.
-/
def branchTable : List (List ℕ × List ℤ × List ℕ) := [
    ([0, 0, 0, 0], [1, 1, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 0, 0, 1], [1, 1, 0, -1], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 0, 1, 0], [1, 1, -1, 1], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 0, 1, 1], [1, 1, -1, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 1, 0, 1], [1, 0, 0, -1], [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 1, 1, 0], [1, 0, -1, 1], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([0, 1, 1, 1], [1, 0, -1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 0, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 0, 0, 1], [0, 1, 1, -1], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 0, 1, 0], [0, 1, 0, 1], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 0, 1, 1], [0, 1, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 1, 0, 0], [0, 0, 1, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 1, 0, 1], [0, 0, 1, -1], [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 1, 1, 0], [0, 0, 0, 1], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ([1, 1, 1, 1], [0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])]

set_option maxRecDepth 10000 in
theorem branch_table_exact : (List.finRange 16).map branchRow = branchTable := by decide

theorem every_coin_sequence_represented (coins : Coins) :
    ∃ i : Fin 16, coinsOfIndex i = coins := by
  revert coins
  decide

theorem fresh (coins : Coins) : FreshStart start (specs coins) := by
  intro s hs p
  simp [start]

theorem identifiers (coins : Coins) :
    (specs coins).map CommitSpec.eventId = [0, 1, 2, 3] := rfl

theorem ids_nodup (coins : Coins) :
    ((specs coins).map CommitSpec.eventId).Nodup := by
  rw [identifiers]
  decide

theorem complete (coins : Coins) : ExecutedCarrierComplete (specs coins) := by
  intro e
  rw [identifiers]
  fin_cases e <;> simp

/-- This log is obtained from the executed native payloads, not postulated. -/
noncomputable def nativeLog (coins : Coins) : SemanticEventLog Port ℕ Event :=
  semanticEventLogOfExecution start (specs coins)
    (fresh coins) (ids_nodup coins) (complete coins)

theorem parent_iff (coins : Coins) (e f : Event) :
    (nativeLog coins).ParentEdge e f ↔ (e = 0 ∧ f = 2) ∨ (e = 2 ∧ f = 3) := by
  rw [nativeLog, semanticEventLogOfExecution_parentEdge_iff,
    ← execParents_iff_authenticated (fresh coins) (ids_nodup coins), execParents_iff]
  fin_cases e <;> fin_cases f <;>
    simp [specs, s0, s1, s2, s3, q1, q2, q3, spec, endpoints, start,
      execChildCites_cons, execChildCites_nil, SpecCites, CommitSpec.apply]

/-- Zero physical current does not erase a retained, version-stamped
instrument invocation. These provenance edges do not assert a signal. -/
theorem zero_current_commits_have_parents :
    (branchRow 15).2.1 = [0, 0, 0, 0] ∧
      (nativeLog (coinsOfIndex 15)).ParentEdge 0 2 ∧
      (nativeLog (coinsOfIndex 15)).ParentEdge 2 3 := by
  refine ⟨by decide, ?_, ?_⟩ <;> simp [parent_iff]

theorem root_heights (coins : Coins) :
    (nativeLog coins).sourceHeight 0 = 0 ∧ (nativeLog coins).sourceHeight 1 = 0 := by
  constructor <;> apply (SemanticEventLog.sourceHeight_eq_zero_iff _ _).2 <;>
    intro p <;> simp [parent_iff]

theorem remaining_heights (coins : Coins) :
    (nativeLog coins).sourceHeight 2 = 1 ∧ (nativeLog coins).sourceHeight 3 = 2 := by
  have h2 : (nativeLog coins).sourceHeight 2 = 1 := by
    rw [SemanticEventLog.sourceHeight_eq]
    simp [Fin.univ_succ, parent_iff, (root_heights coins).1]
  refine ⟨h2, ?_⟩
  rw [SemanticEventLog.sourceHeight_eq]
  simp [Fin.univ_succ, parent_iff, h2]

/-- The midpoint of source ports 0 and 1 in the existing Cartesian chart.
The common source radius is sqrt(2+phi); ports 2 and 3 have its negative.
-/
noncomputable def midpoint : FrameQuotient :=
  quotientEquivVec3.symm ![0, Real.goldenRatio / Real.sqrt (2+Real.goldenRatio), 0]

noncomputable def spatial (e : Event) : FrameQuotient :=
  if e = 1 then -midpoint else midpoint

/-- The positive source midpoint norm; its square is (1+1/sqrt(5))/2. -/
noncomputable def timeScale : ℝ := Real.goldenRatio / Real.sqrt (2+Real.goldenRatio)

theorem timeScale_pos : 0 < timeScale := by
  unfold timeScale
  have hp : 0 < Real.goldenRatio := by
    unfold Real.goldenRatio
    positivity
  exact div_pos hp (Real.sqrt_pos.2 (by linarith))

theorem midpoint_norm : quotientGram midpoint midpoint = timeScale^2 := by
  simp [midpoint, quotientGram, timeScale, OPH.PrimitivePortTranslationBridge.dot,
    Fin.sum_univ_succ, pow_two]

/-- Native authenticated parent edges have equal seam midpoints in this
four-event witness. The all-seam speed bound is checked separately by the
executable finite Gram certificate. -/
noncomputable def placement (coins : Coins) : RankSpatialCausalPlacement (nativeLog coins) where
  timeScale := timeScale
  timeScale_pos := timeScale_pos
  spatial := spatial
  edge_speed_sq := by
    intro e f h
    have same : spatial f = spatial e := by
      rcases (parent_iff coins e f).mp h with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> simp [spatial]
    rw [same, sub_self]
    simpa [quotientGram, OPH.PrimitivePortTranslationBridge.dot] using
      sq_nonneg (timeScale * (((nativeLog coins).sourceHeight f : ℝ) -
        ((nativeLog coins).sourceHeight e : ℝ)))

/-- The existing forward theorem applies to an actual native executed log. -/
theorem native_generatedBeforeEq_sourceCausalLE (coins : Coins) {e f : Event}
    (h : (nativeLog coins).GeneratedBeforeEq e f) :
    sourceCausalLE ((placement coins).eventPoint e) ((placement coins).eventPoint f) :=
  (placement coins).generatedBeforeEq_sourceCausalLE h

/-- The native four-event history retains a future-null comparison which is
not authenticated ancestry. Forward placement is not order reconstruction. -/
theorem native_cone_converse_fails (coins : Coins) :
    sourceCausalLE ((placement coins).eventPoint 1) ((placement coins).eventPoint 3) ∧
      ¬ (nativeLog coins).GeneratedBeforeEq 1 3 := by
  constructor
  · constructor
    · simp [sourceLorentzQ, sourceLorentzB, RankSpatialCausalPlacement.eventPoint,
        placement, spatial, (root_heights coins).2, (remaining_heights coins).2,
        quotientGram, midpoint, timeScale, OPH.PrimitivePortTranslationBridge.dot,
        Fin.sum_univ_succ] <;> ring_nf
      exact le_refl _
    · simp [RankSpatialCausalPlacement.eventPoint, placement,
        (root_heights coins).2, (remaining_heights coins).2]
      exact le_of_lt timeScale_pos
  · rintro (h | h)
    · exact (by decide : (1 : Event) ≠ 3) h
    · have no_from_one : ∀ {e f : Event}, (nativeLog coins).GeneratedBefore e f → e ≠ 1 := by
        intro e f hg
        induction hg with
        | single hedge =>
            rcases (parent_iff coins _ _).mp hedge with ⟨rfl, _⟩ | ⟨rfl, _⟩ <;> decide
        | tail _ _ ih => exact ih
      exact no_from_one h rfl

#print axioms parent_iff
#print axioms branch_table_exact
#print axioms native_generatedBeforeEq_sourceCausalLE
#print axioms native_cone_converse_fails

end OPH.NativeRepairSemanticPlacement
