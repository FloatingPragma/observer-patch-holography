import ObservableNormalForms.MechanismVariants

/-!
# Conditional A2/A5 mechanism-variant example

This finite model instantiates Will's shared A2/A5 mechanism-design fixture
with the OPH observable-normal-form and repair calculus.  The grace-period
and externally-owned-account interpretations are explicit, replaceable
adapters rather than claims about every Harberger implementation.  The model
keeps broad and flash targets, raw and quotiented authority, outcome and
trace-sensitive observations, and strategic predicates separate.
-/

namespace ObservableNormalForms.Examples.AmdA2A5Conditional

open Relation
open MechanismVariants

inductive Actor where
  | incumbent
  | buyer
  deriving DecidableEq, Repr

inductive CallerKind where
  | eoa
  | contract
  deriving DecidableEq, Repr

inductive Variant where
  | baseline
  | a2
  | a5
  deriving DecidableEq, Repr

structure BaseState where
  owner : Actor
  deriving DecidableEq, Repr

inductive A2State where
  | active (owner : Actor)
  | pending (owner acquirer : Actor)
  deriving DecidableEq, Repr

def AmdState : Variant → Type
  | .baseline => BaseState
  | .a2 => A2State
  | .a5 => BaseState

inductive BaselineStep : BaseState → BaseState → Prop where
  | acquire (caller : CallerKind) :
      BaselineStep ⟨.incumbent⟩ ⟨.buyer⟩

inductive A2Step : A2State → A2State → Prop where
  | initiate : A2Step (.active .incumbent) (.pending .incumbent .buyer)
  | complete : A2Step (.pending .incumbent .buyer) (.active .buyer)
  | block : A2Step (.pending .incumbent .buyer) (.active .incumbent)

inductive A5Step : BaseState → BaseState → Prop where
  | acquireEoa : A5Step ⟨.incumbent⟩ ⟨.buyer⟩

def amdStep : (i : Variant) → AmdState i → AmdState i → Prop
  | .baseline => BaselineStep
  | .a2 => A2Step
  | .a5 => A5Step

def amdInitial : (i : Variant) → AmdState i → Prop
  | .baseline, s => s = ⟨.incumbent⟩
  | .a2, s => s = .active .incumbent
  | .a5, s => s = ⟨.incumbent⟩

inductive BaselineTrace where
  | stay
  | acquireEoa
  | acquireContract
  deriving DecidableEq, Repr

inductive A2Trace where
  | stay
  | ownerDeclines (caller : CallerKind)
  | ownerBlocks (caller : CallerKind)
  deriving DecidableEq, Repr

inductive A5Trace where
  | stay
  | acquireEoa
  deriving DecidableEq, Repr

def AmdTrace : Variant → Type
  | .baseline => BaselineTrace
  | .a2 => A2Trace
  | .a5 => A5Trace

def traceStart : (i : Variant) → AmdTrace i → AmdState i
  | .baseline, _ => ⟨.incumbent⟩
  | .a2, _ => .active .incumbent
  | .a5, _ => ⟨.incumbent⟩

def traceFinish : (i : Variant) → AmdTrace i → AmdState i
  | .baseline, .stay => ⟨.incumbent⟩
  | .baseline, .acquireEoa => ⟨.buyer⟩
  | .baseline, .acquireContract => ⟨.buyer⟩
  | .a2, .stay => .active .incumbent
  | .a2, .ownerDeclines _ => .active .buyer
  | .a2, .ownerBlocks _ => .active .incumbent
  | .a5, .stay => ⟨.incumbent⟩
  | .a5, .acquireEoa => ⟨.buyer⟩

theorem trace_starts_initial :
    ∀ (i : Variant) (t : AmdTrace i), amdInitial i (traceStart i t) := by
  intro i t
  cases i <;> cases t <;> rfl

theorem trace_reachable :
    ∀ (i : Variant) (t : AmdTrace i),
      ReflTransGen (amdStep i) (traceStart i t) (traceFinish i t) := by
  intro i t
  cases i with
  | baseline =>
      cases t with
      | stay => exact .refl
      | acquireEoa => exact .tail .refl (BaselineStep.acquire .eoa)
      | acquireContract => exact .tail .refl (BaselineStep.acquire .contract)
  | a2 =>
      cases t with
      | stay => exact .refl
      | ownerDeclines _ =>
          exact .tail (.tail .refl A2Step.initiate) A2Step.complete
      | ownerBlocks _ =>
          exact .tail (.tail .refl A2Step.initiate) A2Step.block
  | a5 =>
      cases t with
      | stay => exact .refl
      | acquireEoa => exact .tail .refl A5Step.acquireEoa

def dynamics : Dynamics Variant where
  State := AmdState
  Step := amdStep
  Initial := amdInitial
  Trace := AmdTrace
  start := traceStart
  finish := traceFinish
  startsInitial := trace_starts_initial
  traceReachable := trace_reachable

/-- Broad reading: any transfer without the modeled pending recourse is bad. -/
def BadBroad : TargetBad dynamics
  | .baseline, .stay => False
  | .baseline, .acquireEoa => True
  | .baseline, .acquireContract => True
  | .a2, .stay => False
  | .a2, .ownerDeclines _ => False
  | .a2, .ownerBlocks _ => False
  | .a5, .stay => False
  | .a5, .acquireEoa => True

/-- Narrow reading: only an atomic contract-origin acquisition is bad. -/
def BadFlash : TargetBad dynamics
  | .baseline, .acquireContract => True
  | _, _ => False

def AllValid : ValidTrace dynamics := fun _ _ => True

def AllCoreOK : CoreOK dynamics := fun _ _ => True

def BroadAdmissible (i : Variant) (t : AmdTrace i) : Prop :=
  Admissible dynamics BadBroad AllValid AllCoreOK i t

def FlashAdmissible (i : Variant) (t : AmdTrace i) : Prop :=
  Admissible dynamics BadFlash AllValid AllCoreOK i t

def finalOwner : (i : Variant) → AmdState i → Actor
  | .baseline, s => s.owner
  | .a2, .active owner => owner
  | .a2, .pending owner _ => owner
  | .a5, s => s.owner

/-- Outcome profile forgets caller kind, pending state, and transition order. -/
def ObsOutcome : Observation dynamics Actor :=
  fun i t => finalOwner i (traceFinish i t)

inductive Route where
  | stayed
  | immediateEoa
  | immediateContract
  | pendingCompleteEoa
  | pendingCompleteContract
  | pendingBlockEoa
  | pendingBlockContract
  deriving DecidableEq, Repr

structure FullObservation where
  outcome : Actor
  route : Route
  deriving DecidableEq, Repr

/-- Full profile retains the modeled acquisition route. -/
def ObsFull : Observation dynamics FullObservation
  | .baseline, .stay => ⟨.incumbent, .stayed⟩
  | .baseline, .acquireEoa => ⟨.buyer, .immediateEoa⟩
  | .baseline, .acquireContract => ⟨.buyer, .immediateContract⟩
  | .a2, .stay => ⟨.incumbent, .stayed⟩
  | .a2, .ownerDeclines .eoa => ⟨.buyer, .pendingCompleteEoa⟩
  | .a2, .ownerDeclines .contract => ⟨.buyer, .pendingCompleteContract⟩
  | .a2, .ownerBlocks .eoa => ⟨.incumbent, .pendingBlockEoa⟩
  | .a2, .ownerBlocks .contract => ⟨.incumbent, .pendingBlockContract⟩
  | .a5, .stay => ⟨.incumbent, .stayed⟩
  | .a5, .acquireEoa => ⟨.buyer, .immediateEoa⟩

def stayObservation : FullObservation := ⟨.incumbent, .stayed⟩
def eoaObservation : FullObservation := ⟨.buyer, .immediateEoa⟩
def contractObservation : FullObservation := ⟨.buyer, .immediateContract⟩

def ProtectedOutcome : Set Actor := Set.univ

/-- Baseline stay, EOA acquisition, and contract-wallet acquisition are fixed
before comparing the variants. -/
def ProtectedFull : Set FullObservation :=
  {stayObservation, eoaObservation, contractObservation}

/-- Hybrid profile forgets the pending delay but retains caller kind. -/
structure HybridObservation where
  outcome : Actor
  caller : Option CallerKind
  deriving DecidableEq, Repr

def ObsHybrid : Observation dynamics HybridObservation
  | .baseline, .stay => ⟨.incumbent, none⟩
  | .baseline, .acquireEoa => ⟨.buyer, some .eoa⟩
  | .baseline, .acquireContract => ⟨.buyer, some .contract⟩
  | .a2, .stay => ⟨.incumbent, none⟩
  | .a2, .ownerDeclines caller => ⟨.buyer, some caller⟩
  | .a2, .ownerBlocks caller => ⟨.incumbent, some caller⟩
  | .a5, .stay => ⟨.incumbent, none⟩
  | .a5, .acquireEoa => ⟨.buyer, some .eoa⟩

def stayHybrid : HybridObservation := ⟨.incumbent, none⟩
def eoaHybrid : HybridObservation := ⟨.buyer, some .eoa⟩
def contractHybrid : HybridObservation := ⟨.buyer, some .contract⟩

def ProtectedHybrid : Set HybridObservation :=
  {stayHybrid, eoaHybrid, contractHybrid}

def broadOutcomeCut (i : Variant) : Set Actor :=
  BehaviorCut dynamics ProtectedOutcome BroadAdmissible ObsOutcome i

def broadFullCut (i : Variant) : Set FullObservation :=
  BehaviorCut dynamics ProtectedFull BroadAdmissible ObsFull i

def flashFullCut (i : Variant) : Set FullObservation :=
  BehaviorCut dynamics ProtectedFull FlashAdmissible ObsFull i

def flashHybridCut (i : Variant) : Set HybridObservation :=
  BehaviorCut dynamics ProtectedHybrid FlashAdmissible ObsHybrid i

theorem a2_ownerDeclines_premises_activate :
    BroadAdmissible .a2 (.ownerDeclines .eoa) ∧
      Reachable dynamics .a2
        (traceStart .a2 (.ownerDeclines .eoa))
        (traceFinish .a2 (.ownerDeclines .eoa)) := by
  constructor
  · simp [BroadAdmissible, Admissible, AllValid, AllCoreOK, BadBroad]
  · exact dynamics.traceReachable .a2 (.ownerDeclines .eoa)

theorem a2_broad_target_eliminated :
    ∀ t : AmdTrace .a2, ¬ BadBroad .a2 t := by
  intro t
  cases t <;> simp [BadBroad]

theorem a5_flash_target_eliminated :
    ∀ t : AmdTrace .a5, ¬ BadFlash .a5 t := by
  intro t
  cases t <;> simp [BadFlash]

theorem a5_broad_target_failure :
    ∃ t : AmdTrace .a5, BadBroad .a5 t := by
  exact ⟨.acquireEoa, by simp [BadBroad]⟩

theorem baseline_broad_and_flash_targets_activate :
    BadBroad .baseline .acquireContract ∧
      BadFlash .baseline .acquireContract := by
  simp [BadBroad, BadFlash]

theorem each_variant_has_initial_state :
    ∀ i : Variant, ∃ s : AmdState i, amdInitial i s := by
  intro i
  cases i with
  | baseline => exact ⟨⟨.incumbent⟩, rfl⟩
  | a2 => exact ⟨.active .incumbent, rfl⟩
  | a5 => exact ⟨⟨.incumbent⟩, rfl⟩

theorem each_variant_has_nontrivial_step :
    ∀ i : Variant, ∃ x y, amdStep i x y := by
  intro i
  cases i with
  | baseline => exact ⟨⟨.incumbent⟩, ⟨.buyer⟩, BaselineStep.acquire .eoa⟩
  | a2 => exact ⟨.active .incumbent, .pending .incumbent .buyer, A2Step.initiate⟩
  | a5 => exact ⟨⟨.incumbent⟩, ⟨.buyer⟩, A5Step.acquireEoa⟩

theorem a5_eoa_premises_activate :
    FlashAdmissible .a5 .acquireEoa ∧
      Reachable dynamics .a5
        (traceStart .a5 .acquireEoa) (traceFinish .a5 .acquireEoa) := by
  constructor
  · simp [FlashAdmissible, Admissible, AllValid, AllCoreOK, BadFlash]
  · exact dynamics.traceReachable .a5 .acquireEoa

theorem a2_broad_outcome_preserving :
    MechanismVariants.ObservationPreserving dynamics ProtectedOutcome BroadAdmissible
      ObsOutcome .a2 := by
  ext o
  simp only [Set.mem_empty_iff_false, iff_false]
  intro ho
  apply ho.2
  cases o with
  | incumbent =>
      exact ⟨.stay, by simp [BroadAdmissible, Admissible, AllValid,
        AllCoreOK, BadBroad, ObsOutcome, finalOwner, traceFinish]⟩
  | buyer =>
      exact ⟨.ownerDeclines .eoa, by simp [BroadAdmissible, Admissible, AllValid,
        AllCoreOK, BadBroad, ObsOutcome, finalOwner, traceFinish]⟩

theorem a2_buyer_not_cut : Actor.buyer ∉ broadOutcomeCut .a2 := by
  intro h
  exact h.2 ⟨.ownerDeclines .eoa, by simp [BroadAdmissible, Admissible,
    AllValid, AllCoreOK, BadBroad, ObsOutcome, finalOwner, traceFinish]⟩

theorem a5_buyer_cut : Actor.buyer ∈ broadOutcomeCut .a5 := by
  constructor
  · exact Set.mem_univ _
  · rintro ⟨t, ht, hobs⟩
    cases t <;> simp [BroadAdmissible, Admissible, AllValid, AllCoreOK,
      BadBroad, ObsOutcome, finalOwner, traceFinish] at ht hobs

theorem a2_pairwise_lt_a5_outcome_broad :
    BehaviorLT broadOutcomeCut .a2 .a5 := by
  apply behaviorLT_of_subset_of_witness (o := Actor.buyer)
  · intro o ho
    have hempty : broadOutcomeCut .a2 = ∅ := by
      simpa [broadOutcomeCut, MechanismVariants.ObservationPreserving] using
        a2_broad_outcome_preserving
    rw [hempty] at ho
    simp only [Set.mem_empty_iff_false] at ho
  · exact a2_buyer_not_cut
  · exact a5_buyer_cut

theorem a2_immediate_eoa_cut_full :
    eoaObservation ∈ broadFullCut .a2 := by
  constructor
  · simp [ProtectedFull, eoaObservation]
  · rintro ⟨t, _, hobs⟩
    cases t with
    | stay => simp [ObsFull, eoaObservation] at hobs
    | ownerDeclines caller =>
        cases caller <;> simp [ObsFull, eoaObservation] at hobs
    | ownerBlocks caller =>
        cases caller <;> simp [ObsFull, eoaObservation] at hobs

theorem a5_immediate_eoa_cut_broad :
    eoaObservation ∈ broadFullCut .a5 := by
  constructor
  · simp [ProtectedFull, eoaObservation]
  · rintro ⟨t, ht, hobs⟩
    cases t <;> simp [BroadAdmissible, Admissible, AllValid, AllCoreOK,
      BadBroad, ObsFull, eoaObservation] at ht hobs

theorem a5_contract_wallet_cut_flash :
    contractObservation ∈ flashFullCut .a5 := by
  constructor
  · simp [ProtectedFull, contractObservation]
  · rintro ⟨t, _, hobs⟩
    cases t <;> simp [ObsFull, contractObservation] at hobs

theorem a2_flash_hybrid_preserving :
    MechanismVariants.ObservationPreserving dynamics ProtectedHybrid
      FlashAdmissible ObsHybrid .a2 := by
  ext o
  simp only [Set.mem_empty_iff_false, iff_false]
  intro ho
  have hcases :
      o = stayHybrid ∨ o = eoaHybrid ∨ o = contractHybrid := by
    simpa [ProtectedHybrid] using ho.1
  apply ho.2
  rcases hcases with hstay | heoa | hcontract
  · subst o
    exact ⟨.stay, by simp [FlashAdmissible, Admissible, AllValid,
      AllCoreOK, BadFlash, ObsHybrid, stayHybrid]⟩
  · subst o
    exact ⟨.ownerDeclines .eoa, by simp [FlashAdmissible, Admissible,
      AllValid, AllCoreOK, BadFlash, ObsHybrid, eoaHybrid]⟩
  · subst o
    exact ⟨.ownerDeclines .contract, by simp [FlashAdmissible, Admissible,
      AllValid, AllCoreOK, BadFlash, ObsHybrid, contractHybrid]⟩

theorem a5_contract_hybrid_cut_flash :
    contractHybrid ∈ flashHybridCut .a5 := by
  constructor
  · simp [ProtectedHybrid, contractHybrid]
  · rintro ⟨t, _, hobs⟩
    cases t <;> simp [ObsHybrid, contractHybrid] at hobs

theorem a2_pairwise_lt_a5_hybrid_flash :
    BehaviorLT flashHybridCut .a2 .a5 := by
  apply behaviorLT_of_subset_of_witness (o := contractHybrid)
  · intro o ho
    have hempty : flashHybridCut .a2 = ∅ := by
      simpa [flashHybridCut, MechanismVariants.ObservationPreserving] using
        a2_flash_hybrid_preserving
    rw [hempty] at ho
    simp only [Set.mem_empty_iff_false] at ho
  · intro hcut
    exact hcut.2 ⟨.ownerDeclines .contract, by
      simp [FlashAdmissible, Admissible, AllValid, AllCoreOK, BadFlash,
        ObsHybrid, contractHybrid]⟩
  · exact a5_contract_hybrid_cut_flash

theorem outcome_equal_immediate_pending :
    ObsOutcome .baseline .acquireEoa =
      ObsOutcome .a2 (.ownerDeclines .eoa) := rfl

theorem full_distinguishes_immediate_pending :
    ObsFull .baseline .acquireEoa ≠
      ObsFull .a2 (.ownerDeclines .eoa) := by
  decide

theorem hybrid_forgets_pending_but_retains_caller :
    ObsHybrid .baseline .acquireEoa =
        ObsHybrid .a2 (.ownerDeclines .eoa) ∧
      ObsHybrid .baseline .acquireContract ≠
        ObsHybrid .baseline .acquireEoa := by
  exact ⟨rfl, by decide⟩

theorem observation_profile_changes_behavior_judgment :
    Actor.buyer ∉ broadOutcomeCut .a2 ∧
      Actor.buyer ∈ broadOutcomeCut .a5 ∧
      eoaObservation ∈ broadFullCut .a2 ∧
      eoaObservation ∈ broadFullCut .a5 :=
  ⟨a2_buyer_not_cut, a5_buyer_cut,
    a2_immediate_eoa_cut_full, a5_immediate_eoa_cut_broad⟩

inductive RawCapability where
  | acquire
  | reassess
  | blockAcquisition
  deriving DecidableEq, Repr

def rawCaps : RawCaps Variant RawCapability
  | .baseline => {.acquire, .reassess}
  | .a2 => {.acquire, .reassess, .blockAcquisition}
  | .a5 => {.acquire, .reassess}

def RawActionView : AuthorityView RawCapability RawCapability := id

inductive AuthorityClass where
  | acquisition
  | ownerControl
  deriving DecidableEq, Repr

def AuthorityClassView : AuthorityView RawCapability AuthorityClass
  | .acquire => .acquisition
  | .reassess => .ownerControl
  | .blockAcquisition => .ownerControl

theorem a2_zero_new_raw_capability_fails :
    ¬ ZeroNewAuthority rawCaps RawActionView .baseline .a2 := by
  intro h
  have hnew : RawCapability.blockAcquisition ∈
      ViewedAuthority rawCaps RawActionView .a2 := by
    exact ⟨.blockAcquisition, by simp [rawCaps], rfl⟩
  have hold := h hnew
  rcases hold with ⟨cap, hcap, hview⟩
  have hcapEq : cap = RawCapability.blockAcquisition := by
    simpa [RawActionView] using hview
  subst cap
  simp [rawCaps] at hcap

theorem a2_zero_new_authority_class :
    ZeroNewAuthority rawCaps AuthorityClassView .baseline .a2 := by
  rintro _ ⟨cap, hcap, rfl⟩
  cases cap <;>
    simp [rawCaps, AuthorityClassView, ViewedAuthority] at hcap ⊢

theorem a5_zero_new_authority_class :
    ZeroNewAuthority rawCaps AuthorityClassView .baseline .a5 := by
  rintro _ ⟨cap, hcap, rfl⟩
  cases cap <;>
    simp [rawCaps, AuthorityClassView, ViewedAuthority] at hcap ⊢

/-- A conservative strategy adapter may reject delayed acquisition even when
the outcome profile loses no protected behavior. -/
def DelayedStrategy : StrategicOK Variant
  | .baseline => True
  | .a2 => False
  | .a5 => True

theorem observation_preservation_does_not_supply_strategy :
    MechanismVariants.ObservationPreserving dynamics ProtectedOutcome BroadAdmissible
        ObsOutcome .a2 ∧
      ¬ DelayedStrategy .a2 := by
  exact ⟨a2_broad_outcome_preserving, by simp [DelayedStrategy]⟩

inductive Strategy where
  | acquire
  | abstain
  deriving DecidableEq, Repr

def baselineUtility : Strategy → Int
  | .acquire => 1
  | .abstain => 0

/-- Synthetic delay/escrow cost used only to refute an automatic implication
from outcome preservation to strategic preservation. -/
def delayedA2Utility : Strategy → Int
  | .acquire => -1
  | .abstain => 0

theorem strategic_choice_changes_under_delay_cost :
    baselineUtility .abstain < baselineUtility .acquire ∧
      delayedA2Utility .acquire < delayedA2Utility .abstain := by
  decide

inductive Support where
  | full
  | missing
  deriving DecidableEq, Repr

def declaredSupports : Supports Variant Support := fun _ => Set.univ

abbrev WriteSpace : Variant → Support → Type := fun _ _ => Unit
abbrev CollarSpace : Variant → Support → Type := fun _ _ => Bool

def supportRelation : Relation WriteSpace CollarSpace
  | _, .full => Set.univ
  | _, .missing => {z | z.2 = false}

/-- The protected outcome representative carried by the fixed collar. -/
def outcomeCollar : Actor → Bool
  | .incumbent => false
  | .buyer => true

theorem outcomeCollar_injective : Function.Injective outcomeCollar := by
  intro a b
  cases a <;> cases b <;> simp [outcomeCollar]

def a2FullEncode
    (t : {t : AmdTrace .a2 // FlashAdmissible .a2 t}) : Unit × Bool :=
  ((), outcomeCollar (ObsOutcome .a2 t.1))

theorem a2FullEncode_outcomeCollar
    (t : {t : AmdTrace .a2 // FlashAdmissible .a2 t}) :
    (a2FullEncode t).2 = outcomeCollar (ObsOutcome .a2 t.1) :=
  rfl

def a2FullRelationEncoding :
    RelationEncoding dynamics FlashAdmissible WriteSpace CollarSpace
      supportRelation .a2 .full where
  encode := a2FullEncode
  exactRange := by
    ext z
    constructor
    · intro _
      rcases z with ⟨u, b⟩
      cases u
      cases b with
      | false =>
          exact ⟨⟨.stay, by simp [FlashAdmissible, Admissible, AllValid,
            AllCoreOK, BadFlash]⟩, rfl⟩
      | true =>
          exact ⟨⟨.ownerDeclines .eoa, by simp [FlashAdmissible, Admissible,
            AllValid, AllCoreOK, BadFlash]⟩, rfl⟩
    · rintro ⟨t, rfl⟩
      simp [supportRelation]

theorem a2FullRelationEncoding_exactRange :
    supportRelation .a2 .full = Set.range a2FullRelationEncoding.encode :=
  a2FullRelationEncoding.exactRange

def supportCost : SupportCost Variant Support Nat
  | _, .full => 1
  | _, .missing => 2

theorem declared_support_is_active : Support.full ∈ declaredSupports .a2 := by
  simp [declaredSupports]

theorem declared_missing_support_is_active :
    Support.missing ∈ declaredSupports .a2 := by
  simp [declaredSupports]

def comparisonFamily : Set Variant := {.a2, .a5}

/-- Broad-target/outcome-observation comparison provenance. -/
def broadOutcomePolicy :
    ComparisonPolicy dynamics Actor RawCapability AuthorityClass Support Nat where
  target := BadBroad
  valid := AllValid
  core := AllCoreOK
  observe := ObsOutcome
  protectedSet := ProtectedOutcome
  rawCaps := rawCaps
  authorityView := AuthorityClassView
  authorityBaseline := .baseline
  supports := declaredSupports
  selectedSupport := fun _ => .full
  supportCost := supportCost
  family := comparisonFamily

/-- Flash-target/hybrid-observation comparison provenance. -/
def flashHybridPolicy :
    ComparisonPolicy dynamics HybridObservation RawCapability AuthorityClass Support Nat where
  target := BadFlash
  valid := AllValid
  core := AllCoreOK
  observe := ObsHybrid
  protectedSet := ProtectedHybrid
  rawCaps := rawCaps
  authorityView := AuthorityClassView
  authorityBaseline := .baseline
  supports := declaredSupports
  selectedSupport := fun _ => .full
  supportCost := supportCost
  family := comparisonFamily

theorem broadOutcomePolicy_a2_eligible : broadOutcomePolicy.Eligible .a2 := by
  refine ⟨?_, ?_, ?_⟩
  · simp [broadOutcomePolicy, comparisonFamily]
  · simp [broadOutcomePolicy, declaredSupports]
  · simpa [broadOutcomePolicy] using a2_zero_new_authority_class

theorem broadOutcomePolicy_a5_eligible : broadOutcomePolicy.Eligible .a5 := by
  refine ⟨?_, ?_, ?_⟩
  · simp [broadOutcomePolicy, comparisonFamily]
  · simp [broadOutcomePolicy, declaredSupports]
  · simpa [broadOutcomePolicy] using a5_zero_new_authority_class

theorem flashHybridPolicy_a2_eligible : flashHybridPolicy.Eligible .a2 := by
  refine ⟨?_, ?_, ?_⟩
  · simp [flashHybridPolicy, comparisonFamily]
  · simp [flashHybridPolicy, declaredSupports]
  · simpa [flashHybridPolicy] using a2_zero_new_authority_class

theorem flashHybridPolicy_a5_eligible : flashHybridPolicy.Eligible .a5 := by
  refine ⟨?_, ?_, ?_⟩
  · simp [flashHybridPolicy, comparisonFamily]
  · simp [flashHybridPolicy, declaredSupports]
  · simpa [flashHybridPolicy] using a5_zero_new_authority_class

theorem broadOutcomePolicy_a2_lt_a5 :
    broadOutcomePolicy.PairwiseLT .a2 .a5 := by
  refine ⟨broadOutcomePolicy_a2_eligible, broadOutcomePolicy_a5_eligible, ?_⟩
  simpa [ComparisonPolicy.cut, ComparisonPolicy.admissible, broadOutcomePolicy,
    broadOutcomeCut, BroadAdmissible] using a2_pairwise_lt_a5_outcome_broad

theorem flashHybridPolicy_a2_lt_a5 :
    flashHybridPolicy.PairwiseLT .a2 .a5 := by
  refine ⟨flashHybridPolicy_a2_eligible, flashHybridPolicy_a5_eligible, ?_⟩
  simpa [ComparisonPolicy.cut, ComparisonPolicy.admissible, flashHybridPolicy,
    flashHybridCut, FlashAdmissible] using a2_pairwise_lt_a5_hybrid_flash

theorem support_costs_are_nontrivial :
    supportCost .a2 .full < supportCost .a2 .missing := by
  decide

theorem support_relations_are_distinct :
    ((), true) ∈ supportRelation .a2 .full ∧
      ((), true) ∉ supportRelation .a2 .missing := by
  simp [supportRelation]

def fullStrongRepair : StrongRepair (supportRelation .a2 .full) where
  toFun := fun z => ⟨z, by simp [supportRelation]⟩
  preservesCollar := fun _ => rfl
  fixes := by
    intro _ _
    rfl

theorem fixed_relation_positive_repairability :
    Nonempty (StrongRepair (supportRelation .a2 .full)) :=
  ⟨fullStrongRepair⟩

theorem fixed_relation_full_collar_covered :
    AdmissibleCollar (supportRelation .a2 .full) Set.univ := by
  intro d _
  exact ⟨(), by simp [supportRelation]⟩

/-- Full fixed-relation repair is equivalent to representatives for every
protected outcome, through the declared exact trace encoding. -/
theorem a2_full_repair_iff_protected_outcome_representatives :
    Nonempty (StrongRepair (supportRelation .a2 .full)) ↔
      ∀ o, o ∈ ProtectedOutcome →
        ∃ t : {t : AmdTrace .a2 // FlashAdmissible .a2 t},
          ObsOutcome .a2 t.1 = o := by
  constructor
  · intro hrepair o _
    have hsurj :=
      (strongRepair_exists_iff_encodedCollarSurjective
        a2FullRelationEncoding).mp hrepair
    rcases hsurj (outcomeCollar o) with ⟨t, ht⟩
    refine ⟨t, outcomeCollar_injective ?_⟩
    exact (a2FullEncode_outcomeCollar t).symm.trans ht
  · intro hrepresentatives
    apply (strongRepair_exists_iff_encodedCollarSurjective
      a2FullRelationEncoding).mpr
    intro b
    cases b with
    | false =>
        rcases hrepresentatives .incumbent (by simp [ProtectedOutcome]) with ⟨t, ht⟩
        refine ⟨t, ?_⟩
        change (a2FullEncode t).2 = false
        rw [a2FullEncode_outcomeCollar, ht]
        rfl
    | true =>
        rcases hrepresentatives .buyer (by simp [ProtectedOutcome]) with ⟨t, ht⟩
        refine ⟨t, ?_⟩
        change (a2FullEncode t).2 = true
        rw [a2FullEncode_outcomeCollar, ht]
        rfl

theorem missing_collar_no_repair :
    IsEmpty (StrongRepair (supportRelation .a2 .missing)) := by
  apply no_strongRepair_of_missing_admissibleCollar (d := true)
  simp [supportRelation]

theorem protected_sets_are_nonempty :
    ProtectedOutcome.Nonempty ∧ ProtectedFull.Nonempty ∧
      ProtectedHybrid.Nonempty := by
  refine ⟨⟨.incumbent, Set.mem_univ _⟩, ?_, ?_⟩
  · exact ⟨stayObservation, by simp [ProtectedFull]⟩
  · exact ⟨stayHybrid, by simp [ProtectedHybrid]⟩

theorem variants_are_distinct : Variant.a2 ≠ Variant.a5 := by
  decide

/-- One finite bundle showing that initiality, nontrivial dynamics, valid
admissible traces, bad-target witnesses, protected observations, and support
premises are simultaneously satisfiable. -/
structure PremiseActivation where
  initial : ∀ i : Variant, ∃ s : AmdState i, amdInitial i s
  step : ∀ i : Variant, ∃ x y, amdStep i x y
  a2Admissible : BroadAdmissible .a2 (.ownerDeclines .eoa)
  a5Admissible : FlashAdmissible .a5 .acquireEoa
  baselineBad : BadBroad .baseline .acquireContract ∧
    BadFlash .baseline .acquireContract
  a5BroadBad : BadBroad .a5 .acquireEoa
  protectedNonempty : ProtectedOutcome.Nonempty ∧ ProtectedFull.Nonempty ∧
    ProtectedHybrid.Nonempty
  supportActive : Support.full ∈ declaredSupports .a2
  collarCovered : AdmissibleCollar (supportRelation .a2 .full) Set.univ

def finiteControlPremises : PremiseActivation where
  initial := each_variant_has_initial_state
  step := each_variant_has_nontrivial_step
  a2Admissible := a2_ownerDeclines_premises_activate.1
  a5Admissible := a5_eoa_premises_activate.1
  baselineBad := baseline_broad_and_flash_targets_activate
  a5BroadBad := by simp [BadBroad]
  protectedNonempty := protected_sets_are_nonempty
  supportActive := declared_support_is_active
  collarCovered := fixed_relation_full_collar_covered

theorem finite_control_premises_satisfiable : Nonempty PremiseActivation :=
  ⟨finiteControlPremises⟩

inductive FamilyVariant where
  | a2
  | a5
  | future
  deriving DecidableEq, Repr

inductive BehaviorAtom where
  | commonLoss
  | extraLoss
  deriving DecidableEq, Repr

def familyCut : FamilyVariant → Set BehaviorAtom
  | .a2 => {.commonLoss}
  | .a5 => {.commonLoss, .extraLoss}
  | .future => ∅

def declaredFamily : Set FamilyVariant := Set.univ

theorem a2_pairwise_lt_a5_in_declared_family :
    BehaviorLT familyCut .a2 .a5 := by
  constructor
  · intro o ho
    cases o <;> simp [familyCut] at ho ⊢
  · intro hreverse
    have hExtra : BehaviorAtom.extraLoss ∈ familyCut .a5 := by
      simp [familyCut]
    have := hreverse hExtra
    simp [familyCut] at this

theorem pairwise_does_not_establish_family_minimum :
    ¬ IsBehaviorMinimum declaredFamily familyCut .a2 := by
  intro hminimum
  have hsub := hminimum.2 .future (by simp [declaredFamily])
  have hCommon : BehaviorAtom.commonLoss ∈ familyCut .a2 := by
    simp [familyCut]
  have := hsub hCommon
  simp [familyCut] at this

/-- One Lean record whose value carries all declared control outcomes
simultaneously.  Its comparison fields deliberately retain two separate typed
comparison policies; the record does not assert one common target, authority,
observation, support, or family policy. -/
structure FullControlActivation where
  broadA2TargetEliminated : ∀ t : AmdTrace .a2, ¬ BadBroad .a2 t
  flashA5TargetEliminated : ∀ t : AmdTrace .a5, ¬ BadFlash .a5 t
  broadA5TargetWitness : ∃ t : AmdTrace .a5, BadBroad .a5 t
  baselineTargetWitnesses : BadBroad .baseline .acquireContract ∧
    BadFlash .baseline .acquireContract
  broadA2FixtureActivation : BroadAdmissible .a2 (.ownerDeclines .eoa) ∧
    ObsOutcome .a2 (.ownerDeclines .eoa) = .buyer ∧
    ¬ BadBroad .a2 (.ownerDeclines .eoa)
  flashA5FixtureActivation : FlashAdmissible .a5 .acquireEoa ∧
    ObsOutcome .a5 .acquireEoa = .buyer ∧
    BadBroad .a5 .acquireEoa
  rawCapabilityInclusionFails : ¬ ZeroNewAuthority rawCaps RawActionView .baseline .a2
  a2AuthorityClassInclusion : ZeroNewAuthority rawCaps AuthorityClassView .baseline .a2
  a5AuthorityClassInclusion : ZeroNewAuthority rawCaps AuthorityClassView .baseline .a5
  protectedSetsNonempty : ProtectedOutcome.Nonempty ∧ ProtectedFull.Nonempty ∧
    ProtectedHybrid.Nonempty
  outcomeObservationEqual : ObsOutcome .baseline .acquireEoa =
    ObsOutcome .a2 (.ownerDeclines .eoa)
  fullObservationDistinguishes : ObsFull .baseline .acquireEoa ≠
    ObsFull .a2 (.ownerDeclines .eoa)
  hybridObservationControl : ObsHybrid .baseline .acquireEoa =
      ObsHybrid .a2 (.ownerDeclines .eoa) ∧
    ObsHybrid .baseline .acquireContract ≠ ObsHybrid .baseline .acquireEoa
  observationProfileCutJudgment : Actor.buyer ∉ broadOutcomeCut .a2 ∧
    Actor.buyer ∈ broadOutcomeCut .a5 ∧
    eoaObservation ∈ broadFullCut .a2 ∧
    eoaObservation ∈ broadFullCut .a5
  broadOutcomePairwiseOrder : broadOutcomePolicy.PairwiseLT .a2 .a5
  flashHybridPairwiseOrder : flashHybridPolicy.PairwiseLT .a2 .a5
  declaredFamilyPairwiseOrder : BehaviorLT familyCut .a2 .a5
  notDeclaredFamilyMinimum : ¬ IsBehaviorMinimum declaredFamily familyCut .a2
  variantIndicesDistinct : Variant.a2 ≠ Variant.a5
  observationalDoesNotImplyStrategic :
    MechanismVariants.ObservationPreserving dynamics ProtectedOutcome BroadAdmissible
        ObsOutcome .a2 ∧
      ¬ DelayedStrategy .a2
  utilityRankingChangesWithDelay : baselineUtility .abstain < baselineUtility .acquire ∧
    delayedA2Utility .acquire < delayedA2Utility .abstain
  eachVariantInitialState : ∀ i : Variant, ∃ s : AmdState i, amdInitial i s
  eachVariantNontrivialStep : ∀ i : Variant, ∃ x y, amdStep i x y
  eachTraceStartsInitial : ∀ (i : Variant) (t : AmdTrace i),
    amdInitial i (traceStart i t)
  eachTraceReachable : ∀ (i : Variant) (t : AmdTrace i),
    ReflTransGen (amdStep i) (traceStart i t) (traceFinish i t)
  fullSupportDeclared : Support.full ∈ declaredSupports .a2
  missingSupportDeclared : Support.missing ∈ declaredSupports .a2
  fullSupportExactTraceRange : supportRelation .a2 .full =
    Set.range a2FullRelationEncoding.encode
  fullRepairRepresentativeEquivalence : Nonempty (StrongRepair (supportRelation .a2 .full)) ↔
    ∀ o, o ∈ ProtectedOutcome →
      ∃ t : {t : AmdTrace .a2 // FlashAdmissible .a2 t},
        ObsOutcome .a2 t.1 = o
  fullSupportRepairExists : Nonempty (StrongRepair (supportRelation .a2 .full))
  fullSupportCollarCoverage : AdmissibleCollar (supportRelation .a2 .full) Set.univ
  missingSupportNoRepair : IsEmpty (StrongRepair (supportRelation .a2 .missing))
  supportCostsDistinct : supportCost .a2 .full < supportCost .a2 .missing
  supportRelationsDistinct : ((), true) ∈ supportRelation .a2 .full ∧
    ((), true) ∉ supportRelation .a2 .missing
  fullWriteSpaceNonempty : Nonempty (WriteSpace .a2 .full)
  fullCollarSpaceNonempty : Nonempty (CollarSpace .a2 .full)

def fullControlActivation : FullControlActivation where
  broadA2TargetEliminated := a2_broad_target_eliminated
  flashA5TargetEliminated := a5_flash_target_eliminated
  broadA5TargetWitness := a5_broad_target_failure
  baselineTargetWitnesses := baseline_broad_and_flash_targets_activate
  broadA2FixtureActivation := by
    exact ⟨a2_ownerDeclines_premises_activate.1, rfl, by simp [BadBroad]⟩
  flashA5FixtureActivation := by
    exact ⟨a5_eoa_premises_activate.1, rfl, by simp [BadBroad]⟩
  rawCapabilityInclusionFails := a2_zero_new_raw_capability_fails
  a2AuthorityClassInclusion := a2_zero_new_authority_class
  a5AuthorityClassInclusion := a5_zero_new_authority_class
  protectedSetsNonempty := protected_sets_are_nonempty
  outcomeObservationEqual := outcome_equal_immediate_pending
  fullObservationDistinguishes := full_distinguishes_immediate_pending
  hybridObservationControl := hybrid_forgets_pending_but_retains_caller
  observationProfileCutJudgment := observation_profile_changes_behavior_judgment
  broadOutcomePairwiseOrder := broadOutcomePolicy_a2_lt_a5
  flashHybridPairwiseOrder := flashHybridPolicy_a2_lt_a5
  declaredFamilyPairwiseOrder := a2_pairwise_lt_a5_in_declared_family
  notDeclaredFamilyMinimum := pairwise_does_not_establish_family_minimum
  variantIndicesDistinct := variants_are_distinct
  observationalDoesNotImplyStrategic := observation_preservation_does_not_supply_strategy
  utilityRankingChangesWithDelay := strategic_choice_changes_under_delay_cost
  eachVariantInitialState := each_variant_has_initial_state
  eachVariantNontrivialStep := each_variant_has_nontrivial_step
  eachTraceStartsInitial := trace_starts_initial
  eachTraceReachable := trace_reachable
  fullSupportDeclared := declared_support_is_active
  missingSupportDeclared := declared_missing_support_is_active
  fullSupportExactTraceRange := a2FullRelationEncoding_exactRange
  fullRepairRepresentativeEquivalence := a2_full_repair_iff_protected_outcome_representatives
  fullSupportRepairExists := fixed_relation_positive_repairability
  fullSupportCollarCoverage := fixed_relation_full_collar_covered
  missingSupportNoRepair := missing_collar_no_repair
  supportCostsDistinct := support_costs_are_nontrivial
  supportRelationsDistinct := support_relations_are_distinct
  fullWriteSpaceNonempty := inferInstance
  fullCollarSpaceNonempty := inferInstance

theorem full_control_activation_satisfiable : Nonempty FullControlActivation :=
  ⟨fullControlActivation⟩

/-! The designated example build doubles as its theorem-level axiom audit. -/

#print axioms a2_ownerDeclines_premises_activate
#print axioms a5_eoa_premises_activate
#print axioms baseline_broad_and_flash_targets_activate
#print axioms a2_broad_target_eliminated
#print axioms a5_flash_target_eliminated
#print axioms a5_broad_target_failure
#print axioms observation_profile_changes_behavior_judgment
#print axioms a2_pairwise_lt_a5_hybrid_flash
#print axioms a2_zero_new_raw_capability_fails
#print axioms a2_zero_new_authority_class
#print axioms observation_preservation_does_not_supply_strategy
#print axioms strategic_choice_changes_under_delay_cost
#print axioms a2FullRelationEncoding
#print axioms fixed_relation_positive_repairability
#print axioms missing_collar_no_repair
#print axioms finite_control_premises_satisfiable
#print axioms a2_pairwise_lt_a5_in_declared_family
#print axioms pairwise_does_not_establish_family_minimum

end ObservableNormalForms.Examples.AmdA2A5Conditional
