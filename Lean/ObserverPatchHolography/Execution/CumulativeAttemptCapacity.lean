import ObserverPatchHolography.Execution.AdaptiveRunStratification
import ObserverPatchHolography.Execution.RankedAttemptCapacity

/-!
# Cumulative attempt capacity for canonical adaptive repair

`adaptiveRun` invokes one site choice at every natural-number index, even when
the chosen `localRepair` stutters. This module charges one unit for every such
attempt and separates that charge from the number of genuine record changes.

The main positive condition is `BoundedWaste C q sigma`: from every scheduler
index and every reducible record, at least one genuine repair occurs among the
next `q + 1` attempts. Under this condition a normal form is reached by attempt
`(q + 1) * mismatchCount initial`. A separate budget premise turns that bound
into an actual capacity guarantee.

The module does not select a scheduler, infer fairness, assign physical time or
rate, or give attempt capacity an energy, bandwidth, fee, or hardware meaning.
-/

namespace OPH.Execution

open OPH
open OPH.Locality

noncomputable section

/-- Every adaptive scheduler invocation costs one attempt unit, including an
equality stutter. -/
def cumulativeAttemptCost (attempts : Nat) : Nat := attempts

/-- Number of genuine record changes in the first `attempts` scheduler
invocations. -/
def cumulativeGenuineChangeCost (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) : Nat := by
  classical
  exact (Finset.range attempts).filter (fun k =>
    adaptiveRun (k + 1) sigma x ≠ adaptiveRun k sigma x) |>.card

/-- Number of equality stutters in the first `attempts` invocations. -/
def cumulativeStutterCost (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) : Nat :=
  attempts - cumulativeGenuineChangeCost C sigma x attempts

theorem cumulativeGenuineChangeCost_le_attempts (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) :
    cumulativeGenuineChangeCost C sigma x attempts ≤ attempts := by
  classical
  unfold cumulativeGenuineChangeCost
  calc
    ((Finset.range attempts).filter fun k =>
        adaptiveRun (k + 1) sigma x ≠ adaptiveRun k sigma x).card ≤
        (Finset.range attempts).card := Finset.card_filter_le _ _
    _ = attempts := Finset.card_range attempts

theorem cumulativeAttemptCost_eq_genuine_add_stutter (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) :
    cumulativeAttemptCost attempts =
      cumulativeGenuineChangeCost C sigma x attempts +
        cumulativeStutterCost C sigma x attempts := by
  have hle := cumulativeGenuineChangeCost_le_attempts C sigma x attempts
  simp only [cumulativeAttemptCost, cumulativeStutterCost]
  omega

/-- The number of genuine changes already paid for, plus the mismatch rank
remaining after those attempts, never exceeds the initial mismatch rank. -/
theorem cumulativeGenuineChangeCost_add_rank_le_initial (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) :
    cumulativeGenuineChangeCost C sigma x attempts +
        mismatchCount (adaptiveRun attempts sigma x) ≤ mismatchCount x := by
  classical
  induction attempts with
  | zero =>
      simp [cumulativeGenuineChangeCost]
  | succ n ih =>
      by_cases hchange :
          adaptiveRun (n + 1) sigma x ≠ adaptiveRun n sigma x
      · have hstrict := adaptiveRun_change_strict_rank C n sigma x hchange
        have hcost :
            cumulativeGenuineChangeCost C sigma x (n + 1) =
              cumulativeGenuineChangeCost C sigma x n + 1 := by
          simp [cumulativeGenuineChangeCost, Finset.range_add_one,
            Finset.filter_insert, Finset.card_insert_of_notMem,
            hchange]
        rw [hcost]
        omega
      · have hsame :
            adaptiveRun (n + 1) sigma x = adaptiveRun n sigma x :=
          not_ne_iff.mp hchange
        have hcost :
            cumulativeGenuineChangeCost C sigma x (n + 1) =
              cumulativeGenuineChangeCost C sigma x n := by
          simp [cumulativeGenuineChangeCost, Finset.range_add_one,
            Finset.filter_insert, hchange]
        rw [hcost, hsame]
        exact ih

/-- The initial mismatch rank directly bounds the cumulative number of
genuine record changes, independently of equality stutters. -/
theorem cumulativeGenuineChangeCost_le_initialMismatch (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (attempts : Nat) :
    cumulativeGenuineChangeCost C sigma x attempts ≤ mismatchCount x := by
  have hbound :=
    cumulativeGenuineChangeCost_add_rank_le_initial C sigma x attempts
  omega

/-- A normal form is reached without exceeding the declared cumulative
attempt budget. -/
def ReachesNormalWithinAttemptBudget (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) (budget : Nat) : Prop :=
  ∃ n, cumulativeAttemptCost n ≤ budget ∧
    NormalForm C (adaptiveRun n sigma x)

/-- From every scheduler index and every reducible record, at least one
genuine change occurs within the next `waste + 1` attempts. This is a
quantitative scheduling premise. It is not scheduler selection or a rate. -/
def BoundedWaste (C : OPHCarrier) (waste : Nat)
    (sigma : AdaptiveScheduler C) : Prop :=
  ∀ start x, (∃ i, localRepair C i x ≠ x) →
    ∃ offset, offset ≤ waste ∧
      adaptiveRun (offset + 1) (fun k => sigma (start + k)) x ≠
        adaptiveRun offset (fun k => sigma (start + k)) x

def adaptiveRankedAttemptSystem (C : OPHCarrier) :
    OPH.RankedAttempt.System (Records C) (Site C) where
  step := localRepair C
  rank := mismatchCount
  change_rank_lt := by
    intro i x hchange
    exact mismatchCount_localRepair_lt C i x hchange

@[simp] theorem rankedAttemptRun_eq_adaptiveRun
    (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    OPH.RankedAttempt.run (adaptiveRankedAttemptSystem C) n sigma x =
      adaptiveRun n sigma x := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [OPH.RankedAttempt.run_last_step, adaptiveRun_last_step, ih]
      rfl

theorem rankedQuiescent_iff_normalForm
    (C : OPHCarrier) (x : Records C) :
    OPH.RankedAttempt.Quiescent (adaptiveRankedAttemptSystem C) x ↔
      NormalForm C x := by
  constructor
  · intro hquiet
    by_contra hnotNormal
    obtain ⟨i, hchange⟩ :=
      (not_normal_iff_exists_firing C x).mp hnotNormal
    exact hchange (hquiet i)
  · intro hnormal i
    by_contra hchange
    exact (not_normal_iff_exists_firing C x).mpr ⟨i, hchange⟩ hnormal

theorem boundedWaste_iff_rankedBoundedWaste
    (C : OPHCarrier) (waste : Nat)
    (sigma : AdaptiveScheduler C) :
    BoundedWaste C waste sigma ↔
      OPH.RankedAttempt.BoundedWaste
        (adaptiveRankedAttemptSystem C) waste sigma := by
  constructor
  · intro h start x hnotQuiet
    have hnotNormal : ¬ NormalForm C x := by
      intro hnormal
      exact hnotQuiet ((rankedQuiescent_iff_normalForm C x).2 hnormal)
    obtain ⟨offset, hoffset, hchange⟩ :=
      h start x ((not_normal_iff_exists_firing C x).mp hnotNormal)
    refine ⟨offset, hoffset, ?_⟩
    simpa only [rankedAttemptRun_eq_adaptiveRun] using hchange
  · intro h start x henabled
    have hnotNormal : ¬ NormalForm C x :=
      (not_normal_iff_exists_firing C x).mpr henabled
    have hnotQuiet :
        ¬ OPH.RankedAttempt.Quiescent (adaptiveRankedAttemptSystem C) x := by
      intro hquiet
      exact hnotNormal ((rankedQuiescent_iff_normalForm C x).1 hquiet)
    obtain ⟨offset, hoffset, hchange⟩ := h start x hnotQuiet
    refine ⟨offset, hoffset, ?_⟩
    simpa only [rankedAttemptRun_eq_adaptiveRun] using hchange

theorem workConserving_boundedWaste_zero (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (hwork : WorkConserving C sigma) :
    BoundedWaste C 0 sigma := by
  intro start x henabled
  refine ⟨0, le_refl 0, ?_⟩
  change localRepair C (sigma start x) x ≠ x
  exact hwork start x henabled

theorem boundedWaste_exists_normal_by_rank (C : OPHCarrier) (waste : Nat) :
    ∀ (x : Records C) (sigma : AdaptiveScheduler C),
      BoundedWaste C waste sigma →
      ∃ N, N ≤ (waste + 1) * mismatchCount x ∧
        NormalForm C (adaptiveRun N sigma x) := by
  intro x sigma hbounded
  have hranked :=
    (boundedWaste_iff_rankedBoundedWaste C waste sigma).mp hbounded
  obtain ⟨N, hN, hquiet⟩ :=
    OPH.RankedAttempt.boundedWaste_exists_quiescent_by_rank
      (adaptiveRankedAttemptSystem C) waste x sigma hranked
  refine ⟨N, hN, ?_⟩
  exact (rankedQuiescent_iff_normalForm C _).1
    (by simpa only [rankedAttemptRun_eq_adaptiveRun] using hquiet)

theorem boundedWaste_eventually_normal (C : OPHCarrier) (waste : Nat)
    (x : Records C) (sigma : AdaptiveScheduler C)
    (hbounded : BoundedWaste C waste sigma) :
    ∃ N, N ≤ (waste + 1) * mismatchCount x ∧
      NormalForm C (adaptiveRun N sigma x) ∧
      ∀ n, N ≤ n → adaptiveRun n sigma x = adaptiveRun N sigma x := by
  obtain ⟨N, hN, hnormal⟩ :=
    boundedWaste_exists_normal_by_rank C waste x sigma hbounded
  exact ⟨N, hN, hnormal,
    normal_at_eventually_constant C N sigma x hnormal⟩

theorem boundedWaste_reaches_within_attempt_budget (C : OPHCarrier)
    (waste budget : Nat) (x : Records C) (sigma : AdaptiveScheduler C)
    (hbounded : BoundedWaste C waste sigma)
    (hbudget : (waste + 1) * mismatchCount x ≤ budget) :
    ReachesNormalWithinAttemptBudget C sigma x budget := by
  obtain ⟨N, hN, hnormal⟩ :=
    boundedWaste_exists_normal_by_rank C waste x sigma hbounded
  exact ⟨N, le_trans hN hbudget, hnormal⟩

#print axioms cumulativeAttemptCost_eq_genuine_add_stutter
#print axioms cumulativeGenuineChangeCost_add_rank_le_initial
#print axioms cumulativeGenuineChangeCost_le_initialMismatch
#print axioms workConserving_boundedWaste_zero
#print axioms boundedWaste_eventually_normal
#print axioms boundedWaste_reaches_within_attempt_budget

end

end OPH.Execution
