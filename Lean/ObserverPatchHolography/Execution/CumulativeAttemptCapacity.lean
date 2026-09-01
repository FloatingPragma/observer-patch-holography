import ObserverPatchHolography.Execution.AdaptiveRunStratification

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
  generalize hd : mismatchCount x = d
  induction d using Nat.strong_induction_on generalizing x sigma with
  | h d ih =>
      by_cases hnormal : NormalForm C x
      · exact ⟨0, Nat.zero_le _, by simpa using hnormal⟩
      · have henabled := (not_normal_iff_exists_firing C x).mp hnormal
        obtain ⟨offset, hoffset, hchange⟩ := hbounded 0 x henabled
        simp at hchange
        let y : Records C := adaptiveRun (offset + 1) sigma x
        let sigmaTail : AdaptiveScheduler C :=
          fun k => sigma ((offset + 1) + k)
        have hylt : mismatchCount y < d := by
          have hstrict :=
            adaptiveRun_change_strict_rank C offset sigma x hchange
          have hprefix := adaptiveRun_rank_le_initial C offset sigma x
          dsimp [y]
          omega
        have hboundedTail : BoundedWaste C waste sigmaTail := by
          intro start z hfire
          obtain ⟨later, hlater, hchanges⟩ :=
            hbounded ((offset + 1) + start) z hfire
          refine ⟨later, hlater, ?_⟩
          simpa [sigmaTail, Nat.add_assoc] using hchanges
        obtain ⟨N, hN, hnormalN⟩ :=
          ih (mismatchCount y) hylt y sigmaTail hboundedTail rfl
        refine ⟨(offset + 1) + N, ?_, ?_⟩
        · have hrank : mismatchCount y + 1 ≤ d := Nat.succ_le_iff.mpr hylt
          have hmul := Nat.mul_le_mul_left (waste + 1) hrank
          rw [Nat.mul_add, Nat.mul_one] at hmul
          omega
        · rw [adaptiveRun_add]
          simpa [y, sigmaTail] using hnormalN

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
