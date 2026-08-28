import ObserverPatchHolography.Execution.AdaptiveRunStratification

/-!
# Nonvacuity and premise boundaries for adaptive repair

The work-conserving branch is inhabited by a source-level scheduler on the
committed TwoCell carrier. The `remoteReader` scheduler is the exact
negative control: it keeps attempting a quiescent remote site while the
probe remains enabled, so its run is constant but never normal.
-/

namespace OPH.Execution.TwoCell

open OPH
open OPH.Execution
open OPH.Locality
open OPH.Locality.TwoCell

noncomputable section

def alwaysProbe : AdaptiveScheduler twoCellCarrier :=
  fun _ _ => false

theorem alwaysProbe_workConserving :
    WorkConserving twoCellCarrier alwaysProbe := by
  intro n x henabled
  obtain ⟨i, hi⟩ := henabled
  cases i
  · exact hi
  · exact False.elim (hi (localRepair_remote_quiescent x))

theorem alwaysProbe_eventually_normal (x : Records twoCellCarrier) :
    ∃ N, N ≤ mismatchCount x ∧
      NormalForm twoCellCarrier (adaptiveRun N alwaysProbe x) ∧
      ∀ n, N ≤ n →
        adaptiveRun n alwaysProbe x = adaptiveRun N alwaysProbe x :=
  workConserving_eventually_normal twoCellCarrier alwaysProbe
    alwaysProbe_workConserving x

theorem alwaysProbe_pathwiseWeakFair (x : Records twoCellCarrier) :
    PathwiseWeakFair twoCellCarrier alwaysProbe x :=
  workConserving_pathwiseWeakFair twoCellCarrier alwaysProbe
    alwaysProbe_workConserving x

theorem remoteReader_flagUp_stutters (n : Nat) :
    adaptiveRun n remoteReader flagUp = flagUp := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [adaptiveRun_succ]
      have hq : localRepair twoCellCarrier
          (remoteReader 0 flagUp) flagUp = flagUp :=
        localRepair_remote_quiescent _
      rw [hq]
      change adaptiveRun n remoteReader flagUp = flagUp
      exact ih

theorem flagUp_probe_fires :
    localRepair twoCellCarrier false flagUp ≠ flagUp := by
  intro h
  have hf := congrFun h false
  rw [localRepair_probe_broken flagUp rfl] at hf
  exact Bool.false_ne_true (by simpa [flagUp] using hf.symm)

theorem flagUp_not_normal :
    ¬ NormalForm twoCellCarrier flagUp := by
  intro hnormal
  exact hnormal _ ⟨false, rfl, flagUp_probe_fires⟩

theorem remoteReader_no_eventual_normalForm :
    ¬ ∃ N, NormalForm twoCellCarrier
      (adaptiveRun N remoteReader flagUp) := by
  rintro ⟨N, hnormal⟩
  rw [remoteReader_flagUp_stutters] at hnormal
  exact flagUp_not_normal hnormal

theorem remoteReader_not_pathwiseWeakFair :
    ¬ PathwiseWeakFair twoCellCarrier remoteReader flagUp := by
  intro hfair
  have hcontinuous : ∀ n, 0 ≤ n →
      localRepair twoCellCarrier false
          (adaptiveRun n remoteReader flagUp) ≠
        adaptiveRun n remoteReader flagUp := by
    intro n _
    rw [remoteReader_flagUp_stutters]
    exact flagUp_probe_fires
  obtain ⟨m, _, hm⟩ := hfair false 0 hcontinuous
  rw [remoteReader_flagUp_stutters] at hm
  simp [remoteReader, flagUp] at hm

#print axioms alwaysProbe_eventually_normal
#print axioms alwaysProbe_pathwiseWeakFair
#print axioms remoteReader_no_eventual_normalForm
#print axioms remoteReader_not_pathwiseWeakFair

end

end OPH.Execution.TwoCell
