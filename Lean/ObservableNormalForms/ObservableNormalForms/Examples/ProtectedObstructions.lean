import ObservableNormalForms.ProtectedObstructions
import ObservableNormalForms.ObserverConfluence
import ObservableNormalForms.Examples.AmdA2A5Conditional

/-!
# Protected-obstruction finite fixtures

Canonical two-file example source: named transport, strict reversal, M0,
TwoBit, C4, and proper-target C6 activation fixtures.
-/

namespace ObservableNormalForms.ProtectedObstructions.NamedFixtures

open scoped BigOperators
open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel

abbrev Two := Fin 2
abbrev Three := Fin 3

/-- Deterministic swap on two states. -/
noncomputable def swap2 : FiniteMarkovKernel Two where
  probability x := fun y => ![![0, 1], ![1, 0]] x y
  probability_nonneg x y := by
    fin_cases x <;> fin_cases y <;> norm_num
  probability_sum_one x := by
    fin_cases x <;> norm_num [Fin.sum_univ_succ]

/-- Deterministic identity kernel on two states. -/
noncomputable def stutter2 : FiniteMarkovKernel Two where
  probability x y := if y = x then 1 else 0
  probability_nonneg x y := by split <;> norm_num
  probability_sum_one x := by
    fin_cases x <;> norm_num [Fin.sum_univ_succ]

/-- Proper-target branching kernel: source `0` branches equally to absorbing
states `1` and `2`; there is no source self-loop. -/
noncomputable def branch3 : FiniteMarkovKernel Three where
  probability x := fun y =>
    ![![0, 1 / 2, 1 / 2],
      ![0, 1, 0],
      ![0, 0, 1]] x y
  probability_nonneg x y := by
    fin_cases x <;> fin_cases y <;> norm_num
  probability_sum_one x := by
    fin_cases x <;> norm_num [Fin.sum_univ_succ]

/-- Quotient chain for the map `0 |-> 0`, `{1,2} |-> 1`. -/
noncomputable def quotient2 : FiniteMarkovKernel Two where
  probability x := fun y =>
    ![![0, 1],
      ![0, 1]] x y
  probability_nonneg x y := by
    fin_cases x <;> fin_cases y <;> norm_num
  probability_sum_one x := by
    fin_cases x <;> norm_num [Fin.sum_univ_succ]

def targetLast2 (x : Two) : Prop := x = 1
def targetOne3 (x : Three) : Prop := x = 1
def targetEndpoints3 (x : Three) : Prop := x = 1 ∨ x = 2

instance : DecidablePred targetLast2 := fun x => inferInstanceAs (Decidable (x = 1))
instance : DecidablePred targetOne3 := fun x => inferInstanceAs (Decidable (x = 1))
instance : DecidablePred targetEndpoints3 := fun x =>
  inferInstanceAs (Decidable (x = 1 ∨ x = 2))

/-! ## tau=0 initial-target fixture -/

theorem tau0_firstHitAtMass :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass swap2 targetLast2 0 (1 : Two) (1 : Two) = 1 := by
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero]
  simp [targetLast2]

theorem tau0_hitBy :
    ObservableNormalForms.ProtectedObstructions.hitBy swap2 targetLast2 0 (1 : Two) = 1 := by
  rw [ObservableNormalForms.ProtectedObstructions.hitBy_zero]
  simp [targetLast2]

/-! ## Unreachable target under an identity kernel -/

theorem unreachable_hitBy_zero (N : ℕ) :
    ObservableNormalForms.ProtectedObstructions.hitBy stutter2 targetLast2 N (0 : Two) = 0 := by
  apply ObservableNormalForms.ProtectedObstructions.hitBy_eq_zero_on_closed stutter2 targetLast2 ({0} : Set Two)
  · intro u hu
    simp only [Set.mem_singleton_iff] at hu
    subst u
    simp [targetLast2]
  · intro u hu v huv
    simp only [Set.mem_singleton_iff] at hu ⊢
    subst u
    by_contra hv
    simp [ObservableNormalForms.ProtectedObstructions.SupportStep, stutter2, hv] at huv
  · simp

theorem unreachable_hitProbability_zero :
    ObservableNormalForms.ProtectedObstructions.hitProbability stutter2 targetLast2 (0 : Two) = 0 := by
  apply ObservableNormalForms.ProtectedObstructions.hitProbability_eq_zero_on_closed
    stutter2 targetLast2 ({0} : Set Two)
  · intro u hu
    simp only [Set.mem_singleton_iff] at hu
    subst u
    simp [targetLast2]
  · intro u hu v huv
    simp only [Set.mem_singleton_iff] at hu ⊢
    subst u
    by_contra hv
    simp [ObservableNormalForms.ProtectedObstructions.SupportStep, stutter2, hv] at huv
  · simp

/-! ## M1: positive and almost-sure proper-target one-step hit -/

theorem m1_firstHitAtMass_one :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass swap2 targetLast2 1 (0 : Two) (1 : Two) = 1 := by
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_succ, if_neg]
  · norm_num [swap2, ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, targetLast2,
      Fin.sum_univ_succ]
  · simp [targetLast2]

theorem m1_hitBy_one :
    ObservableNormalForms.ProtectedObstructions.hitBy swap2 targetLast2 1 (0 : Two) = 1 := by
  rw [ObservableNormalForms.ProtectedObstructions.hitBy_succ_recurrence, if_neg]
  · norm_num [swap2, ObservableNormalForms.ProtectedObstructions.hitBy_zero, targetLast2, Fin.sum_univ_succ]
  · simp [targetLast2]

theorem m1_hitProbability_one :
    ObservableNormalForms.ProtectedObstructions.hitProbability swap2 targetLast2 (0 : Two) = 1 := by
  apply le_antisymm (ObservableNormalForms.ProtectedObstructions.hitProbability_le_one swap2 targetLast2 0)
  exact (show (1 : ℝ) ≤ ObservableNormalForms.ProtectedObstructions.hitProbability swap2 targetLast2 0 by
    simpa [m1_hitBy_one] using
      ObservableNormalForms.ProtectedObstructions.hitBy_le_hitProbability swap2 targetLast2 1 (0 : Two))

/-! ## M2: positive but not almost-sure, with a proper closed trap -/

theorem m2_target_mass_half :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass branch3 targetOne3 1 (0 : Three) (1 : Three) = 1 / 2 := by
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_succ, if_neg]
  · norm_num [branch3, ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, targetOne3,
      Fin.sum_univ_succ]
  · simp [targetOne3]

def m2Trap : ObservableNormalForms.ProtectedObstructions.ReachableClosedTrap branch3 targetOne3 (0 : Three) where
  carrier := {2}
  nonempty := ⟨2, by simp⟩
  target_free := by
    intro u hu
    simp only [Set.mem_singleton_iff] at hu
    subst u
    simp [targetOne3]
  support_closed := by
    intro u hu v huv
    simp only [Set.mem_singleton_iff] at hu ⊢
    subst u
    fin_cases v
    · exfalso
      change 0 < (0 : ℝ) at huv
      norm_num at huv
    · exfalso
      change 0 < (0 : ℝ) at huv
      norm_num at huv
    · rfl
  reachable := by
    refine ⟨2, by simp, ?_⟩
    exact Relation.ReflTransGen.single ⟨by simp [targetOne3], by
      change 0 < (1 / 2 : ℝ)
      norm_num⟩

theorem m2_positive :
    0 < ObservableNormalForms.ProtectedObstructions.endpointHitProbability branch3 targetOne3 (0 : Three) (1 : Three) := by
  apply (lt_ciSup_iff
    (ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass_bddAbove branch3 targetOne3 0 1)).2
  refine ⟨1, ?_⟩
  norm_num [ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass, Finset.sum_range_succ,
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, m2_target_mass_half, targetOne3]

theorem m2_not_almostSure :
    ObservableNormalForms.ProtectedObstructions.hitProbability branch3 targetOne3 (0 : Three) ≠ 1 :=
  ObservableNormalForms.ProtectedObstructions.reachableClosedTrap_imp_hitProbability_ne_one
    branch3 targetOne3 0 m2Trap

/-! ## M3: a.s. proper-target hit with two positive identity classes -/

theorem m3_endpoint_one_half :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass branch3 targetEndpoints3 1
      (0 : Three) (1 : Three) = 1 / 2 := by
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_succ, if_neg]
  · norm_num [branch3, ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, targetEndpoints3,
      Fin.sum_univ_succ] <;> decide
  · simp [targetEndpoints3]

theorem m3_endpoint_two_half :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass branch3 targetEndpoints3 1
      (0 : Three) (2 : Three) = 1 / 2 := by
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_succ, if_neg]
  · norm_num [branch3, ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, targetEndpoints3,
      Fin.sum_univ_succ] <;> decide
  · simp [targetEndpoints3]

theorem m3_hitBy_one :
    ObservableNormalForms.ProtectedObstructions.hitBy branch3 targetEndpoints3 1 (0 : Three) = 1 := by
  rw [ObservableNormalForms.ProtectedObstructions.hitBy_succ_recurrence, if_neg]
  · norm_num [branch3, ObservableNormalForms.ProtectedObstructions.hitBy_zero, targetEndpoints3,
      Fin.sum_univ_succ]
  · simp [targetEndpoints3]

theorem m3_hitProbability_one :
    ObservableNormalForms.ProtectedObstructions.hitProbability branch3 targetEndpoints3 (0 : Three) = 1 := by
  apply le_antisymm (ObservableNormalForms.ProtectedObstructions.hitProbability_le_one branch3 targetEndpoints3 0)
  exact (show (1 : ℝ) ≤ ObservableNormalForms.ProtectedObstructions.hitProbability branch3 targetEndpoints3 0 by
    simpa [m3_hitBy_one] using
      ObservableNormalForms.ProtectedObstructions.hitBy_le_hitProbability branch3 targetEndpoints3 1 (0 : Three))

theorem m3_two_positive_distinct_endpoints :
    (0 < ObservableNormalForms.ProtectedObstructions.endpointHitProbability branch3 targetEndpoints3 0 1) ∧
    (0 < ObservableNormalForms.ProtectedObstructions.endpointHitProbability branch3 targetEndpoints3 0 2) ∧
    (1 : Three) ≠ 2 := by
  constructor
  · apply (lt_ciSup_iff
      (ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass_bddAbove branch3 targetEndpoints3 0 1)).2
    refine ⟨1, ?_⟩
    norm_num [ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass, Finset.sum_range_succ,
      ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, m3_endpoint_one_half, targetEndpoints3]
  · constructor
    · apply (lt_ciSup_iff
        (ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass_bddAbove branch3 targetEndpoints3 0 2)).2
      refine ⟨1, ?_⟩
      have h02 : (0 : Three) ≠ 2 := by decide
      norm_num [ObservableNormalForms.ProtectedObstructions.cumulativeEndpointMass, Finset.sum_range_succ,
        ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, m3_endpoint_two_half, targetEndpoints3,
        h02]
    · decide

/-! ## Nonidentity `3 -> 2` quotient morphism fixture -/

def fixtureStateMap : Three → Two
  | 0 => 0
  | 1 => 1
  | 2 => 1

theorem fixtureStateMap_surjective : Function.Surjective fixtureStateMap := by
  intro b
  fin_cases b
  · exact ⟨0, rfl⟩
  · exact ⟨1, rfl⟩

theorem fixtureStateMap_noninjective : ¬ Function.Injective fixtureStateMap := by
  intro h
  have heq : (1 : Three) = 2 := h (show fixtureStateMap 1 = fixtureStateMap 2 by rfl)
  exact (by decide : (1 : Three) ≠ 2) heq

theorem fixture_target_reflect (a : Three) :
    targetEndpoints3 a ↔ targetLast2 (fixtureStateMap a) := by
  fin_cases a <;> simp [targetEndpoints3, targetLast2, fixtureStateMap]

theorem fixture_target_coverage (b : Two) (hb : targetLast2 b) :
    ∃ a : Three, targetEndpoints3 a ∧ fixtureStateMap a = b := by
  fin_cases b
  · simp [targetLast2] at hb
  · exact ⟨1, by simp [targetEndpoints3], rfl⟩

/-- Exact expectation-form lumping, matching the sealed structural transport
interface and containing no endpoint-pushforward premise. -/
theorem fixture_kernel_lumping (a : Three) (g : Two → ℝ) :
    (∑ a' : Three, branch3.probability a a' * g (fixtureStateMap a')) =
      ∑ b' : Two, quotient2.probability (fixtureStateMap a) b' * g b' := by
  fin_cases a <;>
    simp [branch3, quotient2, fixtureStateMap, Fin.sum_univ_succ] <;> ring

theorem fixture_fiber_nonempty (b : Two) :
    ∃ a : Three, fixtureStateMap a = b :=
  fixtureStateMap_surjective b

/-! The remaining frozen structural surfaces, kept explicit rather than
smuggled into a semantic conclusion field. -/

def initialA (a : Three) : Prop := a = 0
def initialB (b : Two) : Prop := b = 0
def protectedA (_ : Unit) : Prop := True
def protectedB (_ : Unit) : Prop := True
def observationA (_ : Three) : Unit := ()
def observationB (_ : Two) : Unit := ()
def silentA (a a' : Three) : Prop := fixtureStateMap a = fixtureStateMap a'

theorem fixture_protected_reflect (u : Unit) :
    protectedA u ↔ protectedB u := by simp [protectedA, protectedB]

theorem fixture_protected_coverage (u : Unit) (hu : protectedB u) :
    ∃ u' : Unit, protectedA u' := by
  exact ⟨(), by simp [protectedA]⟩

theorem fixture_initial_reflect (a : Three) :
    initialA a ↔ initialB (fixtureStateMap a) := by
  fin_cases a <;> simp [initialA, initialB, fixtureStateMap]

theorem fixture_initial_coverage (b : Two) (hb : initialB b) :
    ∃ a : Three, initialA a ∧ fixtureStateMap a = b := by
  fin_cases b
  · exact ⟨0, rfl, rfl⟩
  · simp [initialB] at hb

theorem fixture_observation_commute (a : Three) :
    observationB (fixtureStateMap a) = observationA a := by rfl

theorem fixture_quotient_exact (a a' : Three) :
    silentA a a' ↔ fixtureStateMap a = fixtureStateMap a' := by
  rfl

theorem fixture_endpoint_silent : silentA 1 2 := by rfl

theorem fixture_source_endpoint_push :
    ObservableNormalForms.ProtectedObstructions.firstHitAtMass quotient2 targetLast2 1 0 1 =
      ObservableNormalForms.ProtectedObstructions.firstHitAtMass branch3 targetEndpoints3 1 0 1 +
      ObservableNormalForms.ProtectedObstructions.firstHitAtMass branch3 targetEndpoints3 1 0 2 := by
  rw [m3_endpoint_one_half, m3_endpoint_two_half]
  rw [ObservableNormalForms.ProtectedObstructions.firstHitAtMass_succ, if_neg]
  · norm_num [quotient2, ObservableNormalForms.ProtectedObstructions.firstHitAtMass_zero, targetLast2,
      Fin.sum_univ_succ]
  · simp [targetLast2]

end ObservableNormalForms.ProtectedObstructions.NamedFixtures

/-! ## Required nonidentity three-to-two-state fixture -/

namespace ObservableNormalForms.ProtectedObstructions.NonidentityExactness.Fixture

inductive SourceState where | s0 | s1 | t
  deriving DecidableEq, Fintype

inductive TargetState where | s | t
  deriving DecidableEq, Fintype

def stateMap : SourceState → TargetState
  | .s0 | .s1 => .s
  | .t => .t

def sourceKernel : FiniteMarkovKernel SourceState where
  probability _ y := if y = .t then 1 else 0
  probability_nonneg := by intro x y; split <;> norm_num
  probability_sum_one := by intro x; simp

def targetKernel : FiniteMarkovKernel TargetState where
  probability _ y := if y = .t then 1 else 0
  probability_nonneg := by intro x y; split <;> norm_num
  probability_sum_one := by intro x; simp

def sourceSilent : Setoid SourceState where
  r x y := stateMap x = stateMap y
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h₁ h₂ => h₁.trans h₂⟩

def targetSilent : Setoid TargetState where
  r x y := x = y
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h₁ h₂ => h₁.trans h₂⟩

def sourceModel : FixedBehaviorModel Unit SourceState where
  protectedSet := {()}
  initial _ := {.s0, .s1}
  initial_nonempty := by intro b hb; exact ⟨.s0, by simp⟩
  target _ q := q == .t
  kernel := sourceKernel
  silent := sourceSilent

def targetModel : FixedBehaviorModel Unit TargetState where
  protectedSet := {()}
  initial _ := {.s}
  initial_nonempty := by intro b hb; exact ⟨.s, by simp⟩
  target _ q := q == .t
  kernel := targetKernel
  silent := targetSilent

def morphism : StructuralMorphism sourceModel targetModel where
  stateMap := stateMap
  behaviorMap := id
  protected_reflect := by intro b; exact Iff.rfl
  protected_cover := by
    intro d hd
    refine ⟨(), ?_, Subsingleton.elim _ _⟩
    change () ∈ ({()} : Finset Unit)
    simp
  initial_reflect := by intro b x hb; cases x <;> simp [sourceModel, targetModel, stateMap]
  initial_cover := by
    intro b y hb hy
    cases y with
    | s => exact ⟨.s0, by simp [sourceModel], rfl⟩
    | t => simp [targetModel] at hy
  kernel_lumping := by
    intro x φ
    cases x <;> simp [sourceModel, targetModel, sourceKernel, targetKernel,
      stateMap, FiniteMarkovKernel.apply]
  target_reflect := by intro b x hb; cases x <;> simp [sourceModel, targetModel, stateMap,
    FixedBehaviorModel.IsTarget]
  target_cover := by
    intro b y hb hy
    cases y with
    | s => simp [targetModel, FixedBehaviorModel.IsTarget] at hy
    | t => exact ⟨.t, by simp [sourceModel, FixedBehaviorModel.IsTarget], rfl⟩
  quotient_exact := by intro c d; rfl

theorem stateMap_surjective : Function.Surjective stateMap := by
  intro y
  cases y with
  | s => exact ⟨.s0, rfl⟩
  | t => exact ⟨.t, rfl⟩

theorem stateMap_not_injective : ¬ Function.Injective stateMap := by
  intro h
  have : SourceState.s0 = SourceState.s1 := h rfl
  cases this

theorem nonidentity_layer0_exact :
    sourceModel.InL0 () ↔ targetModel.InL0 () := morphism.layer0_exact

theorem nonidentity_layer1_exact :
    sourceModel.InL1 () ↔ targetModel.InL1 () := morphism.layer1_exact

theorem nonidentity_layer2_exact :
    sourceModel.InL2 () ↔ targetModel.InL2 () := morphism.layer2_exact

theorem nonidentity_layer3_exact :
    sourceModel.InL3 () ↔ targetModel.InL3 () := morphism.layer3_exact

theorem nonidentity_fiber_exact :
    sourceModel.FiberNonempty () ↔ targetModel.FiberNonempty () :=
  morphism.fiber_nonempty_exact_derived (by simp [sourceModel])

end ObservableNormalForms.ProtectedObstructions.NonidentityExactness.Fixture

/-! # Common-carrier strict-reversal activation -/

namespace ObservableNormalForms.ProtectedObstructions.PublicAdapter.GenericCompletion.StrictReversal

open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel
open ObservableNormalForms.MechanismVariants
open ObservableNormalForms.ProtectedObstructions
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

noncomputable section

abbrev State := Bool × Bool

private def eqSetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

def observe (q : State) : Bool := q.1

def consistentA (q : State) : Prop := q.2 = true
def consistentB (q : State) : Prop := q.1 = true ∧ q.2 = true

def rewriteA (q r : State) : Prop := ¬ consistentA q ∧ r = q
def rewriteB (q r : State) : Prop := ¬ consistentB q ∧ r = q

def identityKernel : FiniteMarkovKernel State where
  probability q r := if r = q then 1 else 0
  probability_nonneg := by intros; split <;> norm_num
  probability_sum_one := by intro q; simp

def destinationB (q : State) : State :=
  if q = (true, false) then (true, true) else q

def progressKernel : FiniteMarkovKernel State where
  probability q r := if r = destinationB q then 1 else 0
  probability_nonneg := by intros; split <;> norm_num
  probability_sum_one := by intro q; simp

def coreA : CoreModel Bool State where
  protectedSet := Finset.univ
  initial b := {(b, false)}
  initial_nonempty := by intro b hb; exact ⟨(b, false), by simp⟩
  target b q := q.2 && decide (q.1 = b)
  kernel := identityKernel
  silent := eqSetoid State

def coreB : CoreModel Bool State where
  protectedSet := Finset.univ
  initial b := {(b, false)}
  initial_nonempty := by intro b hb; exact ⟨(b, false), by simp⟩
  target b q := q.1 && q.2 && decide (q.1 = b)
  kernel := progressKernel
  silent := eqSetoid State

def profileA : Profile Bool State where
  core := coreA
  protected_nonempty := ⟨false, by simp [coreA]⟩
  observe := observe
  consistent := consistentA
  target_iff := by
    intro b q
    rcases q with ⟨x, d⟩
    cases x <;> cases d <;> cases b <;>
      simp [coreA, FixedBehaviorModel.IsTarget, consistentA, observe]
  initial_observe := by
    intro b q hb hq
    have : q = (b, false) := by simpa [coreA] using hq
    subst q
    rfl
  rewrite := rewriteA
  rewrite_observe := by
    intro q r h
    simpa [rewriteA, observe, h.2]
  normal_iff_consistent := by
    rintro ⟨b, d⟩
    cases b <;> cases d <;>
      simp [ObservableNormalForms.IsNormalForm, rewriteA, consistentA]

def profileB : Profile Bool State where
  core := coreB
  protected_nonempty := ⟨false, by simp [coreB]⟩
  observe := observe
  consistent := consistentB
  target_iff := by
    intro b q
    rcases q with ⟨x, d⟩
    cases x <;> cases d <;> cases b <;>
      simp [coreB, FixedBehaviorModel.IsTarget, consistentB, observe]
  initial_observe := by
    intro b q hb hq
    have : q = (b, false) := by simpa [coreB] using hq
    subst q
    rfl
  rewrite := rewriteB
  rewrite_observe := by
    intro q r h
    simpa [rewriteB, observe, h.2]
  normal_iff_consistent := by
    rintro ⟨b, d⟩
    cases b <;> cases d <;>
      simp [ObservableNormalForms.IsNormalForm, rewriteB, consistentB]

theorem profileA_source_hit_zero (b : Bool) :
    profileA.core.hit b (b, false) = 0 := by
  apply hitProbability_eq_zero_on_closed profileA.core.kernel
    (profileA.core.IsTarget b) {q | q.2 = false}
  · intro q hq ht
    have hpair : q.2 = true ∧ q.1 = b := by
      simpa [profileA, coreA, FixedBehaviorModel.IsTarget] using ht
    have ht' : q.2 = true := hpair.1
    cases q.2 <;> simp_all
  · intro q hq r hr
    have hre : r = q := by
      by_contra hne
      simp [Profile.SupportStep, profileA, coreA, identityKernel, hne] at hr
    simpa [hre] using hq
  · simp

theorem profileA_cut0 : profileA.cut0 = ∅ := by
  ext b
  cases b <;> simp [Profile.cut0, Profile.layer0, profileA, coreA,
    FixedBehaviorModel.InL0, FixedBehaviorModel.FiberNonempty,
    FixedBehaviorModel.IsTarget, consistentA, observe]

theorem profileA_cut1 : profileA.cut1 = Set.univ := by
  apply Set.eq_univ_of_forall
  intro b
  refine ⟨by simp [profileA, coreA], ?_⟩
  rintro ⟨h0, q, hq, hpos⟩
  have hq' : q = (b, false) := by simpa [profileA, coreA] using hq
  subst q
  rw [profileA_source_hit_zero] at hpos
  exact (lt_irrefl 0 hpos)

theorem profileB_cut0 : profileB.cut0 = {false} := by
  ext b
  cases b <;> simp [Profile.cut0, Profile.layer0, profileB, coreB,
    FixedBehaviorModel.InL0, FixedBehaviorModel.FiberNonempty,
    FixedBehaviorModel.IsTarget, consistentB, observe]

private def hitPath : Fin 1 → State := fun _ => (true, true)

theorem profileB_true_positive : profileB.core.Positive true := by
  refine ⟨(true, false), by simp [profileB, coreB], ?_⟩
  apply (profileB.hit_pos_iff_endpoint true (true, false)).mpr
  refine ⟨(true, true), ?_⟩
  apply (profileB.endpoint_pos_iff_path (true, false) true (true, true)).mpr
  exact ⟨1, hitPath, by
      simp [FirstHitAt, hitPath, profileB, coreB, FixedBehaviorModel.IsTarget,
        consistentB, observe], by simp [pathEndpoint, hitPath], by
      simp [SupportedFinPath, hitPath, profileB, coreB, progressKernel,
        destinationB]⟩

theorem profileB_cut1 : profileB.cut1 = {false} := by
  ext b
  cases b
  · simp [Profile.cut1, Profile.layer1, profileB, coreB,
      FixedBehaviorModel.InL1, FixedBehaviorModel.InL0,
      FixedBehaviorModel.FiberNonempty, FixedBehaviorModel.IsTarget]
  · simp only [Set.mem_singleton_iff, Bool.true_eq_false]
    change true ∈ profileB.cut1 ↔ False
    constructor
    · intro h
      have hL0 : profileB.core.InL0 true := by
        refine ⟨by simp [profileB, coreB], ?_⟩
        exact ⟨(true, true), by
          simp [profileB, coreB, FixedBehaviorModel.IsTarget]⟩
      exact h.2 ⟨hL0, profileB_true_positive⟩
    · exact False.elim

def dataA : CommonProfileData Bool State Finset.univ where
  initial := profileA.core.initial
  initial_nonempty := by
    intro b hb
    exact profileA.core.initial_nonempty b (by simp [profileA, coreA])
  target := profileA.core.target
  kernel := profileA.core.kernel
  silent := profileA.core.silent
  observe := profileA.observe
  consistent := profileA.consistent
  target_iff := by
    intro b q
    simpa [FixedBehaviorModel.IsTarget] using profileA.target_iff b q
  initial_observe := by
    intro b q hb hq
    exact profileA.initial_observe (by simp [profileA, coreA]) hq
  rewrite := profileA.rewrite
  rewrite_observe := profileA.rewrite_observe
  normal_iff_consistent := profileA.normal_iff_consistent

def dataB : CommonProfileData Bool State Finset.univ where
  initial := profileB.core.initial
  initial_nonempty := by
    intro b hb
    exact profileB.core.initial_nonempty b (by simp [profileB, coreB])
  target := profileB.core.target
  kernel := profileB.core.kernel
  silent := profileB.core.silent
  observe := profileB.observe
  consistent := profileB.consistent
  target_iff := by
    intro b q
    simpa [FixedBehaviorModel.IsTarget] using profileB.target_iff b q
  initial_observe := by
    intro b q hb hq
    exact profileB.initial_observe (by simp [profileB, coreB]) hq
  rewrite := profileB.rewrite
  rewrite_observe := profileB.rewrite_observe
  normal_iff_consistent := profileB.normal_iff_consistent

def family : VariantFamily Bool Bool State where
  protectedCarrier := Finset.univ
  protected_nonempty := ⟨false, by simp⟩
  data i := if i then dataB else dataA

theorem row0_false : family.cutAt 0 false = ∅ := by
  rw [VariantFamily.cutAt_zero]
  simpa [family, VariantFamily.profile, CommonProfileData.toProfile, dataA]
    using profileA_cut0

theorem row0_true : family.cutAt 0 true = {false} := by
  rw [VariantFamily.cutAt_zero]
  simpa [family, VariantFamily.profile, CommonProfileData.toProfile, dataB]
    using profileB_cut0

theorem row1_false : family.cutAt 1 false = Set.univ := by
  rw [VariantFamily.cutAt_one]
  simpa [family, VariantFamily.profile, CommonProfileData.toProfile, dataA]
    using profileA_cut1

theorem row1_true : family.cutAt 1 true = {false} := by
  rw [VariantFamily.cutAt_one]
  simpa [family, VariantFamily.profile, CommonProfileData.toProfile, dataB]
    using profileB_cut1

/-- The actual T7 row 0 strictly prefers `false` to `true`. -/
theorem strict_forward : BehaviorLT (family.cutAt 0) false true := by
  unfold BehaviorLT BehaviorLE
  constructor
  · rw [row0_false]
    exact Set.empty_subset _
  · intro h
    have hf : false ∈ family.cutAt 0 true := by rw [row0_true]; simp
    have := h hf
    rw [row0_false] at this
    exact this

/-- The actual T7 row 1 reverses that strict comparison. -/
theorem strict_reverse : BehaviorLT (family.cutAt 1) true false := by
  unfold BehaviorLT BehaviorLE
  constructor
  · rw [row1_false]
    exact Set.subset_univ _
  · intro h
    have ht : true ∈ family.cutAt 1 false := by rw [row1_false]; simp
    have := h ht
    rw [row1_true] at this
    simpa using this

theorem not_profileLE_false_true : ¬ family.LE false true := by
  intro h
  exact strict_reverse.2 (h 1)

theorem not_profileLE_true_false : ¬ family.LE true false := by
  intro h
  exact strict_forward.2 (h 0)

end

end ObservableNormalForms.ProtectedObstructions.PublicAdapter.GenericCompletion.StrictReversal

/-! # Native M0, TwoBit, and C4 activations -/

namespace ObservableNormalForms.ProtectedObstructions.PublicAdapter.NativeTwoBitC4

open Relation
open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

noncomputable section

private def m0EqSetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

namespace M0

def kernel : FiniteMarkovKernel Bool where
  probability q r := if r = q then 1 else 0
  probability_nonneg := by
    intro q r
    split <;> norm_num
  probability_sum_one := by
    intro q
    cases q <;> norm_num

def core : CoreModel Bool Bool where
  protectedSet := {true}
  initial b := {b}
  initial_nonempty := by
    intro b hb
    exact ⟨b, by simp⟩
  target b q := decide (q = false ∧ q = b)
  kernel := kernel
  silent := m0EqSetoid Bool

def observe (q : Bool) : Bool := q
def consistent (q : Bool) : Prop := q = false
def rewrite (q r : Bool) : Prop := q = true ∧ r = true

def profile : Profile Bool Bool where
  core := core
  protected_nonempty := ⟨true, by simp [core]⟩
  observe := observe
  consistent := consistent
  target_iff := by
    intro b q
    simp [core, FixedBehaviorModel.IsTarget, consistent, observe, and_comm]
  initial_observe := by
    intro b q hb hq
    simpa [core, observe] using hq
  rewrite := rewrite
  rewrite_observe := by
    intro q r hqr
    rcases hqr with ⟨rfl, rfl⟩
    rfl
  normal_iff_consistent := by
    intro q
    cases q <;> simp [ObservableNormalForms.IsNormalForm, rewrite, consistent]

theorem protected_fiber_empty : profile.fiber true = ∅ := by
  ext q
  cases q <;>
    simp [Profile.fiber, profile, core, FixedBehaviorModel.IsTarget]

theorem protected_no_consistent_endpoint :
    ¬ ∃ q, profile.consistent q ∧ profile.observe q = true :=
  (profile.t1_fiber_empty_iff true).mp protected_fiber_empty

theorem protected_in_delta0 : true ∈ profile.delta0 := by
  refine ⟨by simp [core, profile], ?_⟩
  intro h0
  rcases h0.2 with ⟨q, hq⟩
  have hmem : q ∈ profile.fiber true := by
    simpa [Profile.fiber] using hq
  rw [protected_fiber_empty] at hmem
  exact Finset.notMem_empty q hmem

end M0

end

open Relation
open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

noncomputable section

private def twoBitEqSetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

namespace TwoBit

def kernel : FiniteMarkovKernel TwoBitRepair.State where
  probability q r := if r = (q.1, false) then 1 else 0
  probability_nonneg := by
    intro q r
    split <;> norm_num
  probability_sum_one := by
    intro q
    rcases q with ⟨b, d⟩
    cases b <;> cases d <;>
      norm_num [Fintype.sum_prod_type]

theorem kernel_apply (V : TwoBitRepair.State → ℝ) (q : TwoBitRepair.State) :
    kernel.apply V q = V (q.1, false) := by
  rcases q with ⟨b, d⟩
  cases b <;> cases d <;>
    simp [FiniteMarkovKernel.apply, kernel]

def FineTargetPred (b : Bool) (q : TwoBitRepair.State) : Prop :=
  q.2 = false ∧ q.1 = b

instance fineTargetPredDecidable (b : Bool) : DecidablePred (FineTargetPred b) :=
  fun q => inferInstanceAs (Decidable (q.2 = false ∧ q.1 = b))

def fineTarget (b : Bool) (q : TwoBitRepair.State) : Bool :=
  decide (FineTargetPred b q)

def fineCore : CoreModel Bool TwoBitRepair.State where
  protectedSet := Finset.univ
  initial b := {(b, true)}
  initial_nonempty := by
    intro b hb
    exact ⟨(b, true), by simp⟩
  target := fineTarget
  kernel := kernel
  silent := twoBitEqSetoid TwoBitRepair.State

def fineProfile : Profile Bool TwoBitRepair.State where
  core := fineCore
  protected_nonempty := ⟨false, by simp [fineCore]⟩
  observe := TwoBitRepair.observe
  consistent := fun q => q ∈ TwoBitRepair.consistent
  target_iff := by
    intro b q
    simp [fineCore, fineTarget, FineTargetPred, FixedBehaviorModel.IsTarget,
      TwoBitRepair.consistent, TwoBitRepair.observe, and_comm]
  initial_observe := by
    intro b q hb hq
    have hq' : q = (b, true) := by simpa [fineCore] using hq
    subst q
    rfl
  rewrite := TwoBitRepair.step
  rewrite_observe := TwoBitRepair.step_observationPreserving
  normal_iff_consistent := TwoBitRepair.step_completeFor_consistent

def fineScheduler (b : Bool) : fineProfile.SchedulerLaws b where
  supportSound := by
    intro q r hreach hsupport
    rcases q with ⟨bq, dq⟩
    rcases r with ⟨br, dr⟩
    cases bq <;> cases dq <;> cases br <;> cases dr <;>
      simp [Profile.SupportStep, fineProfile, fineCore, kernel,
        TwoBitRepair.step] at hsupport ⊢
  rewriteComplete := by
    intro q r hreach hrewrite
    rcases q with ⟨bq, dq⟩
    rcases r with ⟨br, dr⟩
    cases bq <;> cases dq <;> cases br <;> cases dr <;>
      simp [Profile.SupportStep, fineProfile, fineCore, kernel,
        TwoBitRepair.step] at hrewrite ⊢

theorem fine_first_hit_path (b : Bool) :
    fineProfile.PositiveFirstHitPath (b, true) b (b, false) := by
  refine ⟨1, fun _ => (b, false), ?_, ?_, ?_⟩
  · simp [FirstHitAt, fineProfile, fineCore, fineTarget, FineTargetPred,
      FixedBehaviorModel.IsTarget]
  · rfl
  · simp [SupportedFinPath, fineProfile, fineCore, kernel]

theorem fine_endpoint_positive (b : Bool) :
    0 < fineProfile.core.endpointMass b (b, true) (b, false) :=
  (fineProfile.endpoint_pos_iff_path (b, true) b (b, false)).mpr
    (fine_first_hit_path b)

theorem fine_hit_one (b : Bool) :
    fineProfile.core.hit b (b, true) = 1 := by
  change fineCore.hit b (b, true) = 1
  unfold FixedBehaviorModel.hit
  rw [hitProbability_bellman]
  have hs : ¬ fineCore.IsTarget b (b, true) := by
    simp [fineCore, fineTarget, FineTargetPred, FixedBehaviorModel.IsTarget]
  rw [if_neg hs]
  change kernel.apply
    (fun z => hitProbability kernel (fineCore.IsTarget b) z) (b, true) = 1
  rw [kernel_apply]
  rw [hitProbability_bellman]
  have ht : fineCore.IsTarget b (b, false) := by
    simp [fineCore, fineTarget, FineTargetPred, FixedBehaviorModel.IsTarget]
  rw [if_pos ht]

theorem fine_in_layer2 (b : Bool) : b ∈ fineProfile.layer2 := by
  change fineCore.InL2 b
  refine ⟨⟨⟨by simp [fineCore], ?_⟩, ?_⟩, ?_⟩
  · exact ⟨(b, false), by
      simp [fineCore, fineTarget, FineTargetPred, FixedBehaviorModel.IsTarget]⟩
  · exact ⟨(b, true), by simp [fineCore], by
      have hh : fineCore.hit b (b, true) = 1 := by
        simpa only [fineProfile] using fine_hit_one b
      rw [hh]
      norm_num⟩
  · intro q hq
    have hq' : q = (b, true) := by simpa [fineCore] using hq
    subst q
    exact fine_hit_one b

theorem fine_layer2_eq_layer3 : fineProfile.layer2 = fineProfile.layer3 := by
  apply fineProfile.t5_observableDetermination_collapse
  intro b c d hc hd
  have hc' : c.2 = false ∧ c.1 = b := by
    simpa [fineProfile, fineCore, fineTarget, FineTargetPred,
      FixedBehaviorModel.IsTarget] using hc
  have hd' : d.2 = false ∧ d.1 = b := by
    simpa [fineProfile, fineCore, fineTarget, FineTargetPred,
      FixedBehaviorModel.IsTarget] using hd
  change c = d
  apply Prod.ext
  · exact hc'.2.trans hd'.2.symm
  · exact hc'.1.trans hd'.1.symm

theorem fine_survivor (b : Bool) : b ∈ fineProfile.layer3 := by
  rw [← fine_layer2_eq_layer3]
  exact fine_in_layer2 b

theorem fine_scheduler_activation (b : Bool) :
    Scheduler.PathClosure fineProfile.rewrite (b, true) (b, false) := by
  have hinit : (b, true) ∈ fineProfile.core.initial b := by
    change (b, true) ∈ fineCore.initial b
    simp [fineCore]
  apply fineProfile.t6_support_to_rewrite (fineScheduler b) hinit
  exact .step
    (by simp [fineProfile, fineCore, fineTarget, FineTargetPred,
      FixedBehaviorModel.IsTarget])
    (by simp [Profile.SupportStep, fineProfile, fineCore, kernel])
    (.done (by
      simp [fineProfile, fineCore, fineTarget, FineTargetPred,
        FixedBehaviorModel.IsTarget]))

def CoarseTargetPred (q : TwoBitRepair.State) : Prop := q.2 = false

instance coarseTargetPredDecidable : DecidablePred CoarseTargetPred :=
  fun q => inferInstanceAs (Decidable (q.2 = false))

def coarseTarget (q : TwoBitRepair.State) : Bool := decide (CoarseTargetPred q)

def coarseCore : CoreModel Unit TwoBitRepair.State where
  protectedSet := {()}
  initial _ := {(false, true), (true, true)}
  initial_nonempty := by
    intro b hb
    exact ⟨(false, true), by simp⟩
  target _ := coarseTarget
  kernel := kernel
  silent := twoBitEqSetoid TwoBitRepair.State

def coarseProfile : Profile Unit TwoBitRepair.State where
  core := coarseCore
  protected_nonempty := ⟨(), by simp [coarseCore]⟩
  observe := TwoBitRepair.coarseObserve
  consistent := fun q => q ∈ TwoBitRepair.consistent
  target_iff := by
    intro b q
    cases b
    simp [coarseCore, coarseTarget, CoarseTargetPred, FixedBehaviorModel.IsTarget,
      TwoBitRepair.consistent, TwoBitRepair.coarseObserve]
  initial_observe := by
    intro b q hb hq
    cases b
    rfl
  rewrite := TwoBitRepair.step
  rewrite_observe := by
    intro q r hqr
    rfl
  normal_iff_consistent := TwoBitRepair.step_completeFor_consistent

def coarseScheduler : coarseProfile.SchedulerLaws () where
  supportSound := by
    intro q r hreach hsupport
    rcases q with ⟨bq, dq⟩
    rcases r with ⟨br, dr⟩
    cases bq <;> cases dq <;> cases br <;> cases dr <;>
      simp [Profile.SupportStep, coarseProfile, coarseCore, kernel,
        TwoBitRepair.step] at hsupport ⊢
  rewriteComplete := by
    intro q r hreach hrewrite
    rcases q with ⟨bq, dq⟩
    rcases r with ⟨br, dr⟩
    cases bq <;> cases dq <;> cases br <;> cases dr <;>
      simp [Profile.SupportStep, coarseProfile, coarseCore, kernel,
        TwoBitRepair.step] at hrewrite ⊢

theorem coarse_first_hit_path (b : Bool) :
    coarseProfile.PositiveFirstHitPath (b, true) () (b, false) := by
  refine ⟨1, fun _ => (b, false), ?_, ?_, ?_⟩
  · simp [FirstHitAt, coarseProfile, coarseCore, coarseTarget, CoarseTargetPred,
      FixedBehaviorModel.IsTarget]
  · rfl
  · simp [SupportedFinPath, coarseProfile, coarseCore, kernel]

theorem coarse_hit_one (b : Bool) :
    coarseProfile.core.hit () (b, true) = 1 := by
  change coarseCore.hit () (b, true) = 1
  unfold FixedBehaviorModel.hit
  rw [hitProbability_bellman]
  have hs : ¬ coarseCore.IsTarget () (b, true) := by
    simp [coarseCore, coarseTarget, CoarseTargetPred, FixedBehaviorModel.IsTarget]
  rw [if_neg hs]
  change kernel.apply
    (fun z => hitProbability kernel (coarseCore.IsTarget ()) z) (b, true) = 1
  rw [kernel_apply]
  rw [hitProbability_bellman]
  have ht : coarseCore.IsTarget () (b, false) := by
    simp [coarseCore, coarseTarget, CoarseTargetPred, FixedBehaviorModel.IsTarget]
  rw [if_pos ht]

theorem coarse_in_layer2 : () ∈ coarseProfile.layer2 := by
  change coarseCore.InL2 ()
  refine ⟨⟨⟨by simp [coarseCore], ?_⟩, ?_⟩, ?_⟩
  · exact ⟨(false, false), by
      simp [coarseCore, coarseTarget, CoarseTargetPred,
        FixedBehaviorModel.IsTarget]⟩
  · exact ⟨(false, true), by simp [coarseCore], by
      have hh : coarseCore.hit () (false, true) = 1 := by
        simpa only [coarseProfile] using coarse_hit_one false
      rw [hh]
      norm_num⟩
  · intro q hq
    simp [coarseCore] at hq
    rcases hq with rfl | rfl <;> exact coarse_hit_one _

theorem coarse_not_in_layer3 : () ∉ coarseProfile.layer3 := by
  intro h3
  change coarseCore.InL3 () at h3
  have hall := h3.2.2
  have hleft : 0 < coarseProfile.core.endpointMass ()
      (false, true) (false, false) :=
    (coarseProfile.endpoint_pos_iff_path (false, true) () (false, false)).mpr
      (coarse_first_hit_path false)
  have hright : 0 < coarseProfile.core.endpointMass ()
      (true, true) (true, false) :=
    (coarseProfile.endpoint_pos_iff_path (true, true) () (true, false)).mpr
      (coarse_first_hit_path true)
  have heq := hall (false, true) (by simp [coarseCore]) (false, false) hleft
    (true, true) (by simp [coarseCore]) (true, false) hright
  exact Bool.false_ne_true (congrArg Prod.fst heq)

theorem coarse_in_delta3 : () ∈ coarseProfile.delta3 :=
  ⟨coarse_in_layer2, coarse_not_in_layer3⟩

theorem native_fine_vs_coarse_boundary :
    ObserverEndpointUniqueModulo TwoBitRepair.step TwoBitRepair.observe Eq ∧
      ¬ ObserverEndpointUniqueModulo TwoBitRepair.step
        TwoBitRepair.coarseObserve Eq :=
  ⟨TwoBitRepair.observerEndpointUnique,
    TwoBitRepair.coarse_observerEndpointUnique_fails⟩

end TwoBit

end

open ObservableNormalForms
open ObservableNormalForms.MechanismVariants
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

noncomputable section

private def c4EqSetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

namespace C4

open ObservableNormalForms.Examples.AmdA2A5Conditional

local instance supportFintype : Fintype Support where
  elems := {.full, .missing}
  complete s := by
    cases s <;> simp

def unitKernel : FiniteMarkovKernel Unit where
  probability _ _ := 1
  probability_nonneg := by
    intro q r
    norm_num
  probability_sum_one := by
    intro q
    simp

/-- The actual collar-fiber predicate; this is the only C4 target source. -/
def FiberPred (R : Support → Set (Unit × Bool))
    (sd : Support × Bool) : Prop :=
  ∃ w : Unit, (w, sd.2) ∈ R sd.1

@[reducible] noncomputable def fiberPredDecidable
    (R : Support → Set (Unit × Bool)) : DecidablePred (FiberPred R) := by
  classical
  exact fun sd => inferInstance

/-- A classical Bool reflection of `FiberPred`. The explicit local classical
decision is confined to converting an arbitrary Set predicate to the Bool
field required by `FixedBehaviorModel`. -/
noncomputable def fiberTarget
    (R : Support → Set (Unit × Bool)) (sd : Support × Bool) : Bool := by
  exact @decide (FiberPred R sd) (fiberPredDecidable R sd)

noncomputable def coreOf
    (R : Support → Set (Unit × Bool)) : CoreModel (Support × Bool) Unit where
  protectedSet := Finset.univ
  initial _ := {()}
  initial_nonempty := by
    intro b hb
    exact ⟨(), by simp⟩
  target sd _ := fiberTarget R sd
  kernel := unitKernel
  silent := c4EqSetoid Unit

def nativeRelation (s : Support) : Set (Unit × Bool) :=
  supportRelation .a2 s

abbrev nativeCore : CoreModel (Support × Bool) Unit := coreOf nativeRelation

def AllL0 (R : Support → Set (Unit × Bool)) (s : Support) : Prop :=
  ∀ d : Bool, (coreOf R).InL0 (s, d)

def StaticK0 (R : Support → Set (Unit × Bool)) (s : Support) : Prop :=
  ∃ d : Bool, ¬ (coreOf R).InL0 (s, d)

theorem target_iff_fiberPred
    (R : Support → Set (Unit × Bool)) (s : Support) (d : Bool) :
    (coreOf R).IsTarget (s, d) () ↔ FiberPred R (s, d) := by
  classical
  simp [coreOf, fiberTarget, FixedBehaviorModel.IsTarget]

theorem inL0_iff_native_fiber
    (R : Support → Set (Unit × Bool)) (s : Support) (d : Bool) :
    (coreOf R).InL0 (s, d) ↔ ∃ w : Unit, (w, d) ∈ R s := by
  constructor
  · rintro ⟨hp, q, hq⟩
    cases q
    exact (target_iff_fiberPred R s d).mp hq
  · intro hfiber
    refine ⟨by simp [coreOf], (), ?_⟩
    exact (target_iff_fiberPred R s d).mpr hfiber

theorem allL0_iff_admissibleCollar
    (R : Support → Set (Unit × Bool)) (s : Support) :
    AllL0 R s ↔ AdmissibleCollar (R s) Set.univ := by
  constructor
  · intro hall d hd
    exact (inL0_iff_native_fiber R s d).mp (hall d)
  · intro hcoverage d
    apply (inL0_iff_native_fiber R s d).mpr
    exact hcoverage d (Set.mem_univ d)

theorem allL0_iff_strongRepair
    (R : Support → Set (Unit × Bool)) (s : Support) :
    AllL0 R s ↔ Nonempty (StrongRepair (R s)) := by
  exact (allL0_iff_admissibleCollar R s).trans
    (strongRepair_exists_iff_full_admissibleCollar (R s)).symm

theorem staticK0_iff_not_strongRepair
    (R : Support → Set (Unit × Bool)) (s : Support) :
    StaticK0 R s ↔ ¬ Nonempty (StrongRepair (R s)) := by
  constructor
  · rintro ⟨d, hd⟩ hrepair
    exact hd ((allL0_iff_strongRepair R s).mpr hrepair d)
  · intro hno
    by_contra hcut
    apply hno
    apply (allL0_iff_strongRepair R s).mp
    intro d
    by_contra hd
    exact hcut ⟨d, hd⟩

theorem staticK0_iff_isEmpty
    (R : Support → Set (Unit × Bool)) (s : Support) :
    StaticK0 R s ↔ IsEmpty (StrongRepair (R s)) := by
  constructor
  · intro hk
    constructor
    intro repair
    exact (staticK0_iff_not_strongRepair R s).mp hk ⟨repair⟩
  · intro hempty
    apply (staticK0_iff_not_strongRepair R s).mpr
    rintro ⟨repair⟩
    exact hempty.false repair

theorem full_L0_iff_strongRepair :
    AllL0 nativeRelation .full ↔
      Nonempty (StrongRepair (supportRelation .a2 .full)) :=
  allL0_iff_strongRepair nativeRelation .full

theorem full_native_nonvacuity : AllL0 nativeRelation .full :=
  full_L0_iff_strongRepair.mpr fixed_relation_positive_repairability

theorem missing_native_hole :
    ¬ ∃ w : Unit, (w, true) ∈ nativeRelation .missing := by
  simp [nativeRelation, supportRelation]

theorem missing_in_staticK0 : StaticK0 nativeRelation .missing := by
  refine ⟨true, ?_⟩
  intro hL0
  exact missing_native_hole
    ((inL0_iff_native_fiber nativeRelation .missing true).mp hL0)

theorem missing_K0_iff_noStrongRepair :
    StaticK0 nativeRelation .missing ↔
      IsEmpty (StrongRepair (supportRelation .a2 .missing)) :=
  staticK0_iff_isEmpty nativeRelation .missing

def filledRelation (_ : Support) : Set (Unit × Bool) := Set.univ

theorem filled_relation_allL0 (s : Support) : AllL0 filledRelation s := by
  intro d
  apply (inL0_iff_native_fiber filledRelation s d).mpr
  exact ⟨(), Set.mem_univ _⟩

theorem filled_relation_not_staticK0 (s : Support) :
    ¬ StaticK0 filledRelation s := by
  rintro ⟨d, hd⟩
  exact hd (filled_relation_allL0 s d)

end C4

end

end ObservableNormalForms.ProtectedObstructions.PublicAdapter.NativeTwoBitC4

/-! # Proper-target C6 M2/M3 activations -/

namespace ObservableNormalForms.ProtectedObstructions.PublicAdapter.NativeC6

open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel
open ObservableNormalForms.MechanismVariants
open ObservableNormalForms.ProtectedObstructions
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness
open ObservableNormalForms.ProtectedObstructions.NamedFixtures
open ObservableNormalForms.ProtectedObstructions.PublicAdapter

noncomputable section

private def equalitySetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

def targetOne3Bool (q : Three) : Bool := q == 1
def targetEndpoints3Bool (q : Three) : Bool := q == 1 || q == 2

def m2Core : FixedBehaviorModel Unit Three where
  protectedSet := {()}
  initial _ := {0}
  initial_nonempty := by intro b hb; exact ⟨0, by simp⟩
  target _ := targetOne3Bool
  kernel := branch3
  silent := equalitySetoid Three

def m3Core : FixedBehaviorModel Unit Three where
  protectedSet := {()}
  initial _ := {0}
  initial_nonempty := by intro b hb; exact ⟨0, by simp⟩
  target _ := targetEndpoints3Bool
  kernel := branch3
  silent := equalitySetoid Three

def m2Rewrite (x y : Three) : Prop :=
  (x = 0 ∧ (y = 1 ∨ y = 2)) ∨ (x = 2 ∧ y = 2)

def m3Rewrite (x y : Three) : Prop :=
  x = 0 ∧ (y = 1 ∨ y = 2)

theorem m2_isTarget_iff (q : Three) :
    m2Core.IsTarget () q ↔ targetOne3 q := by
  fin_cases q <;> decide

theorem m3_isTarget_iff (q : Three) :
    m3Core.IsTarget () q ↔ targetEndpoints3 q := by
  fin_cases q <;> decide

theorem m2_normal_iff : ∀ q : Three,
    IsNormalForm m2Rewrite q ↔ targetOne3 q := by
  intro q
  fin_cases q
  · constructor
    · intro hn; exact False.elim (hn 1 (by simp [m2Rewrite]))
    · intro h; simp [targetOne3] at h
  · constructor
    · intro _; rfl
    · intro _ y hy; simp [m2Rewrite] at hy
  · constructor
    · intro hn; exact False.elim (hn 2 (by simp [m2Rewrite]))
    · intro h; simp [targetOne3] at h

theorem m3_normal_iff : ∀ q : Three,
    IsNormalForm m3Rewrite q ↔ targetEndpoints3 q := by
  intro q
  fin_cases q
  · constructor
    · intro hn; exact False.elim (hn 1 (by simp [m3Rewrite]))
    · intro h; simp [targetEndpoints3] at h
  · constructor
    · intro _; simp [targetEndpoints3]
    · intro _ y hy; simp [m3Rewrite] at hy
  · constructor
    · intro _; simp [targetEndpoints3]
    · intro _ y hy; simp [m3Rewrite] at hy

theorem m2_rewrite_observation :
    ObservationPreserving m2Rewrite (fun _ : Three => ()) := by
  intro x y hxy
  rfl

theorem m3_rewrite_observation :
    ObservationPreserving m3Rewrite (fun _ : Three => ()) := by
  intro x y hxy
  rfl

def m2Profile : Profile Unit Three where
  core := m2Core
  protected_nonempty := ⟨(), by simp [m2Core]⟩
  observe := fun _ => ()
  consistent := targetOne3
  target_iff := by
    intro b q
    cases b
    rw [m2_isTarget_iff]
    simp
  initial_observe := by intros; rfl
  rewrite := m2Rewrite
  rewrite_observe := m2_rewrite_observation
  normal_iff_consistent := m2_normal_iff

def m3Profile : Profile Unit Three where
  core := m3Core
  protected_nonempty := ⟨(), by simp [m3Core]⟩
  observe := fun _ => ()
  consistent := targetEndpoints3
  target_iff := by
    intro b q
    cases b
    rw [m3_isTarget_iff]
    simp
  initial_observe := by intros; rfl
  rewrite := m3Rewrite
  rewrite_observe := m3_rewrite_observation
  normal_iff_consistent := m3_normal_iff

theorem m2_support_sound (x y : Three) :
    0 < branch3.probability x y → m2Rewrite x y ∨ x = y := by
  fin_cases x <;> fin_cases y <;> norm_num [branch3, m2Rewrite] <;> decide

theorem m2_rewrite_complete (x y : Three) :
    m2Rewrite x y → 0 < branch3.probability x y := by
  fin_cases x <;> fin_cases y <;> norm_num [branch3, m2Rewrite] <;> decide

theorem m3_support_sound (x y : Three) :
    0 < branch3.probability x y → m3Rewrite x y ∨ x = y := by
  fin_cases x <;> fin_cases y <;> norm_num [branch3, m3Rewrite] <;> decide

theorem m3_rewrite_complete (x y : Three) :
    m3Rewrite x y → 0 < branch3.probability x y := by
  fin_cases x <;> fin_cases y <;> norm_num [branch3, m3Rewrite] <;> decide

theorem source_self_loop_zero : branch3.probability (0 : Three) 0 = 0 := by
  norm_num [branch3]

theorem source_equality_does_not_imply_support :
    (0 : Three) = 0 ∧ ¬ 0 < branch3.probability (0 : Three) 0 := by
  exact ⟨rfl, by norm_num [branch3]⟩

def m2SchedulerLaws : m2Profile.SchedulerLaws () where
  supportSound := by
    intro x y reached hxy
    change 0 < branch3.probability x y at hxy
    change m2Rewrite x y ∨ x = y
    exact m2_support_sound x y hxy
  rewriteComplete := by
    intro x y reached hxy
    change m2Rewrite x y at hxy
    change 0 < branch3.probability x y
    exact m2_rewrite_complete x y hxy

def m3SchedulerLaws : m3Profile.SchedulerLaws () where
  supportSound := by
    intro x y reached hxy
    change 0 < branch3.probability x y at hxy
    change m3Rewrite x y ∨ x = y
    exact m3_support_sound x y hxy
  rewriteComplete := by
    intro x y reached hxy
    change m3Rewrite x y at hxy
    change 0 < branch3.probability x y
    exact m3_rewrite_complete x y hxy

def m2FirstHit :
    Scheduler.FirstHitPath m2Profile.SupportStep
      (m2Core.IsTarget ()) (0 : Three) (1 : Three) := by
  apply Scheduler.FirstHitPath.step
  · intro h
    have : targetOne3 (0 : Three) := (m2_isTarget_iff 0).mp h
    exact (by decide : ¬ targetOne3 (0 : Three)) this
  · change 0 < branch3.probability 0 1
    exact m2_rewrite_complete 0 1 (by simp [m2Rewrite])
  · exact Scheduler.FirstHitPath.done ((m2_isTarget_iff 1).mpr (by decide))

def m3FirstHitRight :
    Scheduler.FirstHitPath m3Profile.SupportStep
      (m3Core.IsTarget ()) (0 : Three) (2 : Three) := by
  apply Scheduler.FirstHitPath.step
  · intro h
    have : targetEndpoints3 (0 : Three) := (m3_isTarget_iff 0).mp h
    exact (by decide : ¬ targetEndpoints3 (0 : Three)) this
  · change 0 < branch3.probability 0 2
    exact m3_rewrite_complete 0 2 (by simp [m3Rewrite])
  · exact Scheduler.FirstHitPath.done ((m3_isTarget_iff 2).mpr (by decide))

theorem m2_t6_support_to_rewrite :
    Scheduler.PathClosure m2Rewrite (0 : Three) (1 : Three) := by
  exact m2Profile.t6_support_to_rewrite m2SchedulerLaws
    (by change (0 : Three) ∈ ({0} : Finset Three); simp) m2FirstHit

theorem m2_t6_rewrite_to_support :
    Scheduler.FirstHitPath m2Profile.SupportStep
      (m2Core.IsTarget ()) (0 : Three) (1 : Three) := by
  apply m2Profile.t6_rewrite_to_support m2SchedulerLaws
  · change (0 : Three) ∈ ({0} : Finset Three); simp
  · exact (m2_isTarget_iff 1).mpr (by decide)
  · exact Scheduler.PathClosure.head
      (show m2Rewrite (0 : Three) 1 from Or.inl ⟨rfl, Or.inl rfl⟩) (.refl _)

theorem m3_t6_both_directions :
    Scheduler.PathClosure m3Rewrite (0 : Three) (2 : Three) ∧
    Scheduler.FirstHitPath m3Profile.SupportStep
      (m3Core.IsTarget ()) (0 : Three) (2 : Three) := by
  constructor
  · exact m3Profile.t6_support_to_rewrite m3SchedulerLaws
      (by change (0 : Three) ∈ ({0} : Finset Three); simp) m3FirstHitRight
  · apply m3Profile.t6_rewrite_to_support m3SchedulerLaws
    · change (0 : Three) ∈ ({0} : Finset Three); simp
    · exact (m3_isTarget_iff 2).mpr (by decide)
    · exact Scheduler.PathClosure.head
        (show m3Rewrite (0 : Three) 2 from ⟨rfl, Or.inr rfl⟩) (.refl _)

theorem m2_firstHit_half :
    firstHitAtMass branch3 (m2Core.IsTarget ()) 1
      (0 : Three) (1 : Three) = 1 / 2 := by
  rw [firstHitAtMass_succ, if_neg]
  · norm_num [branch3, firstHitAtMass_zero, m2Core,
      targetOne3Bool, FixedBehaviorModel.IsTarget, Fin.sum_univ_succ]
  · decide

def m2TrapCore : ReachableClosedTrap branch3 (m2Core.IsTarget ()) (0 : Three) where
  carrier := {2}
  nonempty := ⟨2, by simp⟩
  target_free := by
    intro u hu
    simp only [Set.mem_singleton_iff] at hu
    subst u
    decide
  support_closed := by
    intro u hu v huv
    simp only [Set.mem_singleton_iff] at hu ⊢
    subst u
    fin_cases v
    · exfalso; change 0 < (0 : ℝ) at huv; norm_num at huv
    · exfalso; change 0 < (0 : ℝ) at huv; norm_num at huv
    · rfl
  reachable := by
    refine ⟨2, by simp, ?_⟩
    exact Relation.ReflTransGen.single ⟨by decide, by
      change 0 < (1 / 2 : ℝ); norm_num⟩

theorem m3_firstHit_one_half :
    firstHitAtMass branch3 (m3Core.IsTarget ()) 1
      (0 : Three) (1 : Three) = 1 / 2 := by
  rw [firstHitAtMass_succ, if_neg]
  · norm_num [branch3, firstHitAtMass_zero, m3Core,
      targetEndpoints3Bool, FixedBehaviorModel.IsTarget,
      Fin.sum_univ_succ] <;> decide
  · decide

theorem m3_firstHit_two_half :
    firstHitAtMass branch3 (m3Core.IsTarget ()) 1
      (0 : Three) (2 : Three) = 1 / 2 := by
  rw [firstHitAtMass_succ, if_neg]
  · norm_num [branch3, firstHitAtMass_zero, m3Core,
      targetEndpoints3Bool, FixedBehaviorModel.IsTarget,
      Fin.sum_univ_succ] <;> decide
  · decide

theorem m3_hitBy_one :
    hitBy branch3 (m3Core.IsTarget ()) 1 (0 : Three) = 1 := by
  rw [hitBy_succ_recurrence, if_neg]
  · norm_num [branch3, hitBy_zero, m3Core, targetEndpoints3Bool,
      FixedBehaviorModel.IsTarget, Fin.sum_univ_succ]
  · decide

private def pathToOne : Fin 1 → Three := fun _ => 1
private def pathToTwo : Fin 1 → Three := fun _ => 2

theorem m2_endpoint_positive :
    0 < m2Core.endpointMass () (0 : Three) (1 : Three) := by
  unfold FixedBehaviorModel.endpointMass
  apply (endpointHitProbability_pos_iff branch3 (m2Core.IsTarget ()) 0 1).2
  refine ⟨1, pathToOne, ?_, ?_, ?_⟩
  · simp [FirstHitAt, pathToOne, m2Core, targetOne3Bool,
      FixedBehaviorModel.IsTarget]
  · simp [pathToOne, pathEndpoint]
  · norm_num [SupportedFinPath, pathToOne, branch3]

theorem m2_hit_not_one : m2Core.hit () (0 : Three) ≠ 1 := by
  unfold FixedBehaviorModel.hit
  exact reachableClosedTrap_imp_hitProbability_ne_one
    branch3 (m2Core.IsTarget ()) 0 m2TrapCore

theorem m3_endpoint_one_positive :
    0 < m3Core.endpointMass () (0 : Three) (1 : Three) := by
  unfold FixedBehaviorModel.endpointMass
  apply (endpointHitProbability_pos_iff branch3 (m3Core.IsTarget ()) 0 1).2
  refine ⟨1, pathToOne, ?_, ?_, ?_⟩
  · simp [FirstHitAt, pathToOne, m3Core, targetEndpoints3Bool,
      FixedBehaviorModel.IsTarget]
  · simp [pathToOne, pathEndpoint]
  · norm_num [SupportedFinPath, pathToOne, branch3]

theorem m3_endpoint_two_positive :
    0 < m3Core.endpointMass () (0 : Three) (2 : Three) := by
  unfold FixedBehaviorModel.endpointMass
  apply (endpointHitProbability_pos_iff branch3 (m3Core.IsTarget ()) 0 2).2
  refine ⟨1, pathToTwo, ?_, ?_, ?_⟩
  · simp [FirstHitAt, pathToTwo, m3Core, targetEndpoints3Bool,
      FixedBehaviorModel.IsTarget]
  · simp [pathToTwo, pathEndpoint]
  · change 0 < (1 / 2 : ℝ) ∧ True
    norm_num

theorem m3_hit_one : m3Core.hit () (0 : Three) = 1 := by
  unfold FixedBehaviorModel.hit
  apply le_antisymm (hitProbability_le_one branch3 (m3Core.IsTarget ()) 0)
  exact (show (1 : ℝ) ≤ hitProbability branch3 (m3Core.IsTarget ()) 0 by
    simpa [m3_hitBy_one] using
      hitBy_le_hitProbability branch3 (m3Core.IsTarget ()) 1 (0 : Three))

theorem m2_in_layer1 : () ∈ m2Profile.layer1 := by
  refine ⟨⟨by change () ∈ ({()} : Finset Unit); simp, ?_⟩, ?_⟩
  · exact ⟨1, (m2_isTarget_iff 1).mpr (by decide)⟩
  · refine ⟨0, ?_, ?_⟩
    · change (0 : Three) ∈ ({0} : Finset Three); simp
    · exact (m2Profile.hit_pos_iff_endpoint () 0).2 ⟨1, m2_endpoint_positive⟩

theorem m2_not_in_layer2 : () ∉ m2Profile.layer2 := by
  intro h2
  apply m2_hit_not_one
  apply h2.2 0
  change (0 : Three) ∈ ({0} : Finset Three)
  simp

theorem m2_delta2 : () ∈ m2Profile.delta2 := ⟨m2_in_layer1, m2_not_in_layer2⟩

theorem m3_in_layer2 : () ∈ m3Profile.layer2 := by
  refine ⟨?_, ?_⟩
  · refine ⟨⟨by change () ∈ ({()} : Finset Unit); simp, ?_⟩, ?_⟩
    · exact ⟨1, (m3_isTarget_iff 1).mpr (by decide)⟩
    · refine ⟨0, ?_, ?_⟩
      · change (0 : Three) ∈ ({0} : Finset Three); simp
      · exact (m3Profile.hit_pos_iff_endpoint () 0).2
          ⟨1, m3_endpoint_one_positive⟩
  · intro q hq
    have hq0 : q = 0 := by
      change q ∈ ({0} : Finset Three) at hq
      simpa using hq
    subst q
    exact m3_hit_one

theorem m3_not_in_layer3 : () ∉ m3Profile.layer3 := by
  intro h3
  have hsilent := h3.2.2
      (0 : Three) (by change (0 : Three) ∈ ({0} : Finset Three); simp)
      (1 : Three) m3_endpoint_one_positive
      (0 : Three) (by change (0 : Three) ∈ ({0} : Finset Three); simp)
      (2 : Three) m3_endpoint_two_positive
  change (1 : Three) = 2 at hsilent
  exact (by decide : (1 : Three) ≠ 2) hsilent

theorem m3_delta3 : () ∈ m3Profile.delta3 := ⟨m3_in_layer2, m3_not_in_layer3⟩

theorem m2_t7_level1_live :
    () ∈ RealizedBehavior m2Profile.supportedDynamics
      m2Profile.Admissible1 m2Profile.traceObserve () := by
  rw [m2Profile.t7_realizedBehavior1]
  exact m2_in_layer1

theorem m2_t7_level2_rejects :
    () ∉ RealizedBehavior m2Profile.supportedDynamics
      m2Profile.Admissible2 m2Profile.traceObserve () := by
  rw [m2Profile.t7_realizedBehavior2]
  exact m2_not_in_layer2

theorem m3_t7_level2_live :
    () ∈ RealizedBehavior m3Profile.supportedDynamics
      m3Profile.Admissible2 m3Profile.traceObserve () := by
  rw [m3Profile.t7_realizedBehavior2]
  exact m3_in_layer2

theorem m3_t7_level3_rejects :
    () ∉ RealizedBehavior m3Profile.supportedDynamics
      m3Profile.Admissible3 m3Profile.traceObserve () := by
  rw [m3Profile.t7_realizedBehavior3]
  exact m3_not_in_layer3

end

end ObservableNormalForms.ProtectedObstructions.PublicAdapter.NativeC6
