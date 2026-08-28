import ObserverPatchHolography.Execution.CumulativeAttemptCapacity
import ObserverPatchHolography.Execution.AdaptiveRunCounterexamples

/-!
# Sharp controls for cumulative attempt capacity

The committed `TwoCell` source supplies a nonzero bounded-waste inhabitant and
an eventually normalizing scheduler with arbitrarily delayed first repair.
The independent-defect family below is the sharpness carrier: each site owns
one self-loop mismatch, so a genuine local repair can remove exactly one
mismatch and cannot shortcut another site.
-/

namespace OPH.Execution.CumulativeCapacityExamples

open OPH
open OPH.Execution
open OPH.Locality

noncomputable section

/-! ## Committed TwoCell controls -/

/-- One quiescent remote attempt alternates with one probe attempt. -/
def alternatingProbe :
    AdaptiveScheduler OPH.Locality.TwoCell.twoCellCarrier :=
  fun n _ => if n % 2 = 0 then true else false

theorem alternatingProbe_boundedWaste_one :
    BoundedWaste OPH.Locality.TwoCell.twoCellCarrier 1 alternatingProbe := by
  intro start x henabled
  have hprobe : localRepair OPH.Locality.TwoCell.twoCellCarrier false x ≠ x := by
    obtain ⟨i, hi⟩ := henabled
    cases i with
    | false => exact hi
    | true =>
        exact False.elim
          (hi (OPH.Locality.TwoCell.localRepair_remote_quiescent x))
  have hmodlt : start % 2 < 2 := Nat.mod_lt _ (by omega)
  by_cases hone : start % 2 = 1
  · refine ⟨0, by omega, ?_⟩
    change localRepair OPH.Locality.TwoCell.twoCellCarrier
      (alternatingProbe start x) x ≠ x
    simpa [alternatingProbe, hone] using hprobe
  · have hzero : start % 2 = 0 := by omega
    have hnext : (start + 1) % 2 = 1 := by omega
    refine ⟨1, by omega, ?_⟩
    let sigmaTail : AdaptiveScheduler OPH.Locality.TwoCell.twoCellCarrier :=
      fun k => alternatingProbe (start + k)
    have hrunOne : adaptiveRun 1 sigmaTail x = x := by
      rw [adaptiveRun_last_step]
      change localRepair OPH.Locality.TwoCell.twoCellCarrier
        (alternatingProbe start x) x = x
      simpa [alternatingProbe, hzero] using
        OPH.Locality.TwoCell.localRepair_remote_quiescent x
    change adaptiveRun 2 sigmaTail x ≠ adaptiveRun 1 sigmaTail x
    rw [adaptiveRun_last_step, hrunOne]
    change localRepair OPH.Locality.TwoCell.twoCellCarrier
      (alternatingProbe (start + 1) x) x ≠ x
    simpa [alternatingProbe, hnext] using hprobe

/-- Stutter for `delay` attempts, then use the genuine probe repair. -/
def delayThenProbe (delay : Nat) :
    AdaptiveScheduler OPH.Locality.TwoCell.twoCellCarrier :=
  fun n _ => if n < delay then true else false

theorem twoCellFlagUp_initialRank_eq_one :
    mismatchCount OPH.Locality.TwoCell.flagUp = 1 := by
  have hpos : 0 < mismatchCount OPH.Locality.TwoCell.flagUp := by
    have hlt := mismatchCount_localRepair_lt
      OPH.Locality.TwoCell.twoCellCarrier false OPH.Locality.TwoCell.flagUp
      OPH.Execution.TwoCell.flagUp_probe_fires
    omega
  have hcard : (OPH.brokenSet OPH.Locality.TwoCell.flagUp).card ≤ 1 := by
    calc
      (OPH.brokenSet OPH.Locality.TwoCell.flagUp).card ≤
          (Finset.univ : Finset Unit).card :=
        Finset.card_le_card (Finset.subset_univ _)
      _ = 1 := by simp
  change 0 < (OPH.brokenSet OPH.Locality.TwoCell.flagUp).card at hpos
  change (OPH.brokenSet OPH.Locality.TwoCell.flagUp).card = 1
  omega

theorem alwaysProbe_normal_at_one :
    NormalForm OPH.Locality.TwoCell.twoCellCarrier
      (adaptiveRun 1 OPH.Execution.TwoCell.alwaysProbe
        OPH.Locality.TwoCell.flagUp) := by
  obtain ⟨N, hN, hnormal, _⟩ :=
    OPH.Execution.TwoCell.alwaysProbe_eventually_normal
      OPH.Locality.TwoCell.flagUp
  rw [twoCellFlagUp_initialRank_eq_one] at hN
  have hNpos : 0 < N := by
    by_contra hnot
    have hzero : N = 0 := Nat.eq_zero_of_not_pos hnot
    subst N
    exact OPH.Execution.TwoCell.flagUp_not_normal hnormal
  have hNone : N = 1 := by omega
  simpa [hNone] using hnormal

theorem delayThenProbe_stutters_through (delay : Nat) :
    ∀ n, n ≤ delay →
      adaptiveRun n (delayThenProbe delay) OPH.Locality.TwoCell.flagUp =
        OPH.Locality.TwoCell.flagUp := by
  intro n hn
  induction n with
  | zero => rfl
  | succ n ih =>
      have hnle : n ≤ delay := Nat.le_trans (Nat.le_succ n) hn
      have hnlt : n < delay := Nat.lt_of_succ_le hn
      rw [adaptiveRun_last_step, ih hnle]
      change localRepair OPH.Locality.TwoCell.twoCellCarrier
        (delayThenProbe delay n OPH.Locality.TwoCell.flagUp)
        OPH.Locality.TwoCell.flagUp = OPH.Locality.TwoCell.flagUp
      simpa [delayThenProbe, hnlt] using
        OPH.Locality.TwoCell.localRepair_remote_quiescent
          OPH.Locality.TwoCell.flagUp

theorem delayThenProbe_normal_at_threshold (delay : Nat) :
    NormalForm OPH.Locality.TwoCell.twoCellCarrier
      (adaptiveRun (delay + 1) (delayThenProbe delay)
        OPH.Locality.TwoCell.flagUp) := by
  rw [adaptiveRun_last_step,
    delayThenProbe_stutters_through delay delay (le_refl _)]
  simpa [delayThenProbe, OPH.Execution.TwoCell.alwaysProbe,
    adaptiveRun_last_step] using alwaysProbe_normal_at_one

theorem delayed_normalizing_attempt_no_go (delay : Nat) :
    mismatchCount OPH.Locality.TwoCell.flagUp = 1 ∧
    (∀ n, n ≤ delay →
      ¬ NormalForm OPH.Locality.TwoCell.twoCellCarrier
        (adaptiveRun n (delayThenProbe delay)
          OPH.Locality.TwoCell.flagUp)) ∧
    NormalForm OPH.Locality.TwoCell.twoCellCarrier
      (adaptiveRun (delay + 1) (delayThenProbe delay)
        OPH.Locality.TwoCell.flagUp) := by
  refine ⟨twoCellFlagUp_initialRank_eq_one, ?_,
    delayThenProbe_normal_at_threshold delay⟩
  intro n hn hnormal
  rw [delayThenProbe_stutters_through delay n hn] at hnormal
  exact OPH.Execution.TwoCell.flagUp_not_normal hnormal

/-! ## Parametric independent-defect sharpness carrier -/

def independentDefectCarrier (width : Nat) : OPHCarrier where
  Patch := Fin width
  State := fun _ => Bool
  Edge := Fin width
  src := id
  tgt := id
  Iface := fun _ => Bool
  projSrc := fun _ s => s
  projTgt := fun _ _ => true
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases h : a = b
    · rw [if_pos h]
      exact ⟨fun _ => h, fun _ => rfl⟩
    · rw [if_neg h]
      exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 h⟩

def allBroken (width : Nat) : Records (independentDefectCarrier width) :=
  fun _ => false

def allRepaired (width : Nat) : Records (independentDefectCarrier width) :=
  fun _ => true

def prefixRepaired (width repaired : Nat) :
    Records (independentDefectCarrier width) :=
  fun i => decide (i.val < repaired)

def sequentialScheduler (width : Nat) [NeZero width] :
    AdaptiveScheduler (independentDefectCarrier width) :=
  fun n _ => ⟨n % width, Nat.mod_lt _ (NeZero.pos width)⟩

instance independentEdgeDecidableEq (width : Nat) :
    DecidableEq (independentDefectCarrier width).Edge :=
  inferInstanceAs (DecidableEq (Fin width))

theorem independent_solvesAt_iff (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width))
    (s : Bool) :
    SolvesAt (C := independentDefectCarrier width) i x s ↔ s = true := by
  constructor
  · intro h
    have hi := h i (Or.inl rfl)
    simpa [edgeConsistentAt, independentDefectCarrier] using hi
  · rintro rfl e hinc
    change e = i ∨ e = i at hinc
    rcases hinc with rfl | rfl <;>
      simp [edgeConsistentAt, independentDefectCarrier]

theorem independent_trigger_iff (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width)) :
    LocalTrigger (C := independentDefectCarrier width) i x ↔ x i = false := by
  constructor
  · rintro ⟨e, hinc, hbroken⟩
    change e = i ∨ e = i at hinc
    rcases hinc with rfl | rfl <;>
      simpa [edgeConsistentAt, independentDefectCarrier] using hbroken
  · intro hfalse
    refine ⟨i, Or.inl rfl, ?_⟩
    simp [edgeConsistentAt, independentDefectCarrier, hfalse]

theorem independent_edgeConsistent_iff (width : Nat)
    (e : Fin width) (x : Records (independentDefectCarrier width)) :
    edgeConsistentAt e x ↔ x e = true := by
  rfl

theorem independent_mem_brokenSet_iff (width : Nat)
    (e : Fin width) (x : Records (independentDefectCarrier width)) :
    e ∈ brokenSet x ↔ x e = false := by
  constructor
  · intro hmem
    cases hxe : x e with
    | false => rfl
    | true =>
        exfalso
        exact (mem_brokenSet_iff_not_consistent
          (C := independentDefectCarrier width)).mp hmem
          ((independent_edgeConsistent_iff width e x).mpr hxe)
  · intro hfalse
    apply (mem_brokenSet_iff_not_consistent
      (C := independentDefectCarrier width)).mpr
    intro hconsistent
    have htrue := (independent_edgeConsistent_iff width e x).mp hconsistent
    simp [hfalse] at htrue

theorem independent_localRepair_eq_update_true (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width)) :
    localRepair (independentDefectCarrier width) i x =
      Function.update x i true := by
  cases hxi : x i with
  | false =>
      have hfire :
          LocalTrigger (C := independentDefectCarrier width) i x ∧
            LocallySolvable (C := independentDefectCarrier width) i x :=
        ⟨(independent_trigger_iff width i x).mpr hxi,
          ⟨true, (independent_solvesAt_iff width i x true).mpr rfl⟩⟩
      rw [localRepair_of_fire (independentDefectCarrier width) i x hfire]
      have hchosen : Classical.choose hfire.2 = true :=
        (independent_solvesAt_iff width i x _).mp
          (Classical.choose_spec hfire.2)
      rw [hchosen]
  | true =>
      have hquiet :
          ¬ (LocalTrigger (C := independentDefectCarrier width) i x ∧
            LocallySolvable (C := independentDefectCarrier width) i x) := by
        rintro ⟨htrigger, _⟩
        have := (independent_trigger_iff width i x).mp htrigger
        simp [hxi] at this
      rw [localRepair_of_quiescent (independentDefectCarrier width) i x hquiet]
      funext j
      by_cases hji : j = i
      · subst j
        simp [hxi]
      · simp [Function.update_of_ne hji]

theorem independent_localRepair_ne_iff (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width)) :
    localRepair (independentDefectCarrier width) i x ≠ x ↔ x i = false := by
  rw [independent_localRepair_eq_update_true]
  constructor
  · intro hne
    cases hxi : x i with
    | false => rfl
    | true =>
        exfalso
        apply hne
        funext j
        by_cases hji : j = i
        · subst j
          simp [hxi]
        · simp [Function.update_of_ne hji]
  · intro hfalse heq
    have hi := congrFun heq i
    simp [hfalse] at hi

theorem independent_brokenSet_after_repair (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width)) :
    brokenSet (localRepair (independentDefectCarrier width) i x) =
      (brokenSet x).erase i := by
  classical
  rw [independent_localRepair_eq_update_true]
  ext e
  constructor
  · intro hmem
    have hvalue := (independent_mem_brokenSet_iff width e _).mp hmem
    have hei : e ≠ i := by
      intro heq
      subst e
      simp at hvalue
    have hsame : Function.update x i true (show Fin width from e) = x e :=
      Function.update_of_ne hei true x
    apply Finset.mem_erase.mpr
    refine ⟨hei, (independent_mem_brokenSet_iff width e x).mpr ?_⟩
    change Function.update x i true (show Fin width from e) = false at hvalue
    rw [hsame] at hvalue
    exact hvalue
  · intro hmem
    obtain ⟨hei, hbroken⟩ := Finset.mem_erase.mp hmem
    have hvalue := (independent_mem_brokenSet_iff width e x).mp hbroken
    have hsame : Function.update x i true (show Fin width from e) = x e :=
      Function.update_of_ne hei true x
    apply (independent_mem_brokenSet_iff width e _).mpr
    change Function.update x i true (show Fin width from e) = false
    rw [hsame]
    exact hvalue

theorem independent_localRepair_exact_unit (width : Nat)
    (i : Fin width) (x : Records (independentDefectCarrier width))
    (hfire : localRepair (independentDefectCarrier width) i x ≠ x) :
    mismatchCount (localRepair (independentDefectCarrier width) i x) + 1 =
      mismatchCount x := by
  classical
  have hfalse := (independent_localRepair_ne_iff width i x).mp hfire
  have himem : i ∈ brokenSet x := by
    exact (independent_mem_brokenSet_iff width i x).mpr hfalse
  change (brokenSet (localRepair (independentDefectCarrier width) i x)).card + 1 =
    (brokenSet x).card
  rw [independent_brokenSet_after_repair width i x]
  exact Finset.card_erase_add_one himem

theorem independent_allBroken_initialRank (width : Nat) :
    mismatchCount (allBroken width) = width := by
  classical
  change (brokenSet (allBroken width)).card = width
  have hset : brokenSet (allBroken width) = (Finset.univ : Finset (Fin width)) := by
    apply Finset.eq_univ_iff_forall.mpr
    intro e
    exact (independent_mem_brokenSet_iff width e _).mpr rfl
  rw [hset]
  change (Finset.univ : Finset (Fin width)).card = width
  simp

theorem independent_normal_rank_zero (width : Nat)
    (x : Records (independentDefectCarrier width))
    (hnormal : NormalForm (independentDefectCarrier width) x) :
    mismatchCount x = 0 := by
  classical
  have hall : ∀ i, x i = true := by
    intro i
    cases hxi : x i with
    | false =>
        exfalso
        exact hnormal _ ⟨i, rfl,
          (independent_localRepair_ne_iff width i x).mpr hxi⟩
    | true => rfl
  change (brokenSet x).card = 0
  have hset : brokenSet x = ∅ := by
    ext e
    simp only [Finset.notMem_empty, iff_false]
    intro hmem
    have hfalse := (independent_mem_brokenSet_iff width e x).mp hmem
    simp [hall e] at hfalse
  simp [hset]

theorem independent_rank_zero_normal (width : Nat)
    (x : Records (independentDefectCarrier width))
    (hrank : mismatchCount x = 0) :
    NormalForm (independentDefectCarrier width) x := by
  intro y hstep
  obtain ⟨i, rfl, hfire⟩ := hstep
  have hfalse := (independent_localRepair_ne_iff width i x).mp hfire
  have himem := (independent_mem_brokenSet_iff width i x).mpr hfalse
  change (brokenSet x).card = 0 at hrank
  have hempty : brokenSet x = ∅ := Finset.card_eq_zero.mp hrank
  rw [hempty] at himem
  exact Finset.notMem_empty i himem

def firstBrokenScheduler (width : Nat) [NeZero width] :
    AdaptiveScheduler (independentDefectCarrier width) := by
  classical
  exact fun _ x =>
    if h : ∃ i, x i = false then Classical.choose h
    else ⟨0, NeZero.pos width⟩

theorem firstBrokenScheduler_selects_false (width : Nat) [NeZero width]
    (n : Nat) (x : Records (independentDefectCarrier width))
    (h : ∃ i, x i = false) :
    x (firstBrokenScheduler width n x) = false := by
  classical
  have hselect : firstBrokenScheduler width n x = Classical.choose h := by
    unfold firstBrokenScheduler
    dsimp
    split
    · congr
    · contradiction
  rw [hselect]
  exact Classical.choose_spec h

theorem firstBrokenScheduler_workConserving (width : Nat) [NeZero width] :
    WorkConserving (independentDefectCarrier width)
      (firstBrokenScheduler width) := by
  intro n x henabled
  have hfalse : ∃ i, x i = false := by
    obtain ⟨i, hi⟩ := henabled
    exact ⟨i, (independent_localRepair_ne_iff width i x).mp hi⟩
  apply (independent_localRepair_ne_iff width _ x).mpr
  exact firstBrokenScheduler_selects_false width n x hfalse

theorem independent_adaptive_change_exact_unit (width : Nat)
    (sigma : AdaptiveScheduler (independentDefectCarrier width))
    (x : Records (independentDefectCarrier width)) (n : Nat)
    (hchange : adaptiveRun (n + 1) sigma x ≠ adaptiveRun n sigma x) :
    mismatchCount (adaptiveRun (n + 1) sigma x) + 1 =
      mismatchCount (adaptiveRun n sigma x) := by
  rw [adaptiveRun_last_step] at hchange ⊢
  exact independent_localRepair_exact_unit width _ _ hchange

theorem firstBroken_run_rank (width : Nat) [NeZero width] :
    ∀ n, n ≤ width →
      mismatchCount
          (adaptiveRun n (firstBrokenScheduler width) (allBroken width)) =
        width - n := by
  intro n hn
  induction n with
  | zero => simpa using independent_allBroken_initialRank width
  | succ n ih =>
      have hnle : n ≤ width := Nat.le_trans (Nat.le_succ n) hn
      have hnlt : n < width := Nat.lt_of_succ_le hn
      have hprev := ih hnle
      have hprevPos :
          0 < mismatchCount
            (adaptiveRun n (firstBrokenScheduler width) (allBroken width)) := by
        rw [hprev]
        omega
      have hnotnormal :
          ¬ NormalForm (independentDefectCarrier width)
            (adaptiveRun n (firstBrokenScheduler width) (allBroken width)) := by
        intro hnormal
        have hzero := independent_normal_rank_zero width _ hnormal
        omega
      have henabled := (not_normal_iff_exists_firing
        (independentDefectCarrier width) _).mp hnotnormal
      have hfire := firstBrokenScheduler_workConserving width n _ henabled
      have hchange :
          adaptiveRun (n + 1) (firstBrokenScheduler width) (allBroken width) ≠
            adaptiveRun n (firstBrokenScheduler width) (allBroken width) := by
        rw [adaptiveRun_last_step]
        exact hfire
      have hexact := independent_adaptive_change_exact_unit width
        (firstBrokenScheduler width) (allBroken width) n hchange
      omega

theorem firstBroken_normal_at_width (width : Nat) [NeZero width] :
    NormalForm (independentDefectCarrier width)
      (adaptiveRun width (firstBrokenScheduler width) (allBroken width)) := by
  apply independent_rank_zero_normal width
  simpa using firstBroken_run_rank width width (le_refl _)

theorem firstBroken_attempt_threshold_iff (width : Nat) [NeZero width]
    (budget : Nat) :
    ReachesNormalWithinAttemptBudget (independentDefectCarrier width)
      (firstBrokenScheduler width) (allBroken width) budget ↔
      width ≤ budget := by
  constructor
  · rintro ⟨n, hn, hnormal⟩
    have hcost : n ≤ budget := by simpa [cumulativeAttemptCost] using hn
    by_cases hnwidth : n ≤ width
    · have hrank := firstBroken_run_rank width n hnwidth
      have hzero := independent_normal_rank_zero width _ hnormal
      omega
    · omega
  · intro hbudget
    exact ⟨width, by simpa [cumulativeAttemptCost] using hbudget,
      firstBroken_normal_at_width width⟩

/-- The existing TwoCell source is a structurally nonisomorphic second sharp
instance: it has two patches rather than the one-patch width-one independent
carrier. -/
theorem twoCell_attempt_threshold_iff (budget : Nat) :
    ReachesNormalWithinAttemptBudget OPH.Locality.TwoCell.twoCellCarrier
      OPH.Execution.TwoCell.alwaysProbe OPH.Locality.TwoCell.flagUp budget ↔
      1 ≤ budget := by
  constructor
  · rintro ⟨n, hn, hnormal⟩
    have hnpos : 0 < n := by
      by_contra hnot
      have hzero : n = 0 := Nat.eq_zero_of_not_pos hnot
      subst n
      exact OPH.Execution.TwoCell.flagUp_not_normal hnormal
    simpa [cumulativeAttemptCost] using le_trans hnpos hn
  · intro hbudget
    refine ⟨1, by simpa [cumulativeAttemptCost] using hbudget,
      alwaysProbe_normal_at_one⟩

theorem sharpInstances_patchCardinality_differs :
    Fintype.card OPH.Locality.TwoCell.twoCellCarrier.Patch = 2 ∧
    Fintype.card (independentDefectCarrier 1).Patch = 1 := by
  decide

#print axioms alternatingProbe_boundedWaste_one
#print axioms delayed_normalizing_attempt_no_go
#print axioms firstBrokenScheduler_workConserving
#print axioms firstBroken_attempt_threshold_iff
#print axioms twoCell_attempt_threshold_iff
#print axioms sharpInstances_patchCardinality_differs

end

end OPH.Execution.CumulativeCapacityExamples
