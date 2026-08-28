import ObserverPatchHolography.Locality.AdaptiveScheduler

namespace OPH.Execution

open OPH
open OPH.Locality
open Relation

noncomputable section

theorem adaptiveRun_last_step {C : OPHCarrier} (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    adaptiveRun (n + 1) sigma x =
      localRepair C (sigma n (adaptiveRun n sigma x))
        (adaptiveRun n sigma x) := by
  induction n generalizing sigma x with
  | zero => rfl
  | succ n ih =>
      rw [adaptiveRun_succ, ih]
      rfl

theorem adaptiveRun_add (C : OPHCarrier) (m n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    adaptiveRun (m + n) sigma x =
      adaptiveRun n (fun k => sigma (m + k))
        (adaptiveRun m sigma x) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [Nat.add_succ, adaptiveRun_last_step,
        adaptiveRun_last_step, ih]

theorem localRepair_mismatchCount_le (C : OPHCarrier)
    (i : Site C) (x : Records C) :
    mismatchCount (localRepair C i x) ≤ mismatchCount x := by
  by_cases hfire : localRepair C i x ≠ x
  · exact Nat.le_of_lt (mismatchCount_localRepair_lt C i x hfire)
  · rw [not_ne_iff.mp hfire]

theorem adaptiveRun_step_rank_le (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    mismatchCount (adaptiveRun (n + 1) sigma x) ≤
      mismatchCount (adaptiveRun n sigma x) := by
  rw [adaptiveRun_last_step]
  exact localRepair_mismatchCount_le C _ _

theorem adaptiveRun_change_strict_rank (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C)
    (hchange : adaptiveRun (n + 1) sigma x ≠ adaptiveRun n sigma x) :
    mismatchCount (adaptiveRun (n + 1) sigma x) <
      mismatchCount (adaptiveRun n sigma x) := by
  rw [adaptiveRun_last_step] at hchange ⊢
  exact mismatchCount_localRepair_lt C _ _ hchange

theorem adaptiveRun_rank_le_initial (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    mismatchCount (adaptiveRun n sigma x) ≤ mismatchCount x := by
  induction n with
  | zero => rfl
  | succ n ih =>
      exact le_trans (adaptiveRun_step_rank_le C n sigma x) ih

theorem adaptiveRun_reachable (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C) :
    ReflTransGen (acceptedStep C) x (adaptiveRun n sigma x) := by
  induction n with
  | zero => exact ReflTransGen.refl
  | succ n ih =>
      rw [adaptiveRun_last_step]
      by_cases hfire :
          localRepair C (sigma n (adaptiveRun n sigma x))
            (adaptiveRun n sigma x) ≠ adaptiveRun n sigma x
      · exact ReflTransGen.tail ih
          ⟨sigma n (adaptiveRun n sigma x), rfl, hfire⟩
      · rw [not_ne_iff.mp hfire]
        exact ih

theorem not_normal_iff_exists_firing (C : OPHCarrier)
    (x : Records C) :
    ¬ NormalForm C x ↔ ∃ i, localRepair C i x ≠ x := by
  constructor
  · intro hnormal
    by_contra hnone
    apply hnormal
    intro y hstep
    obtain ⟨i, _, hfire⟩ := hstep
    exact hnone ⟨i, hfire⟩
  · rintro ⟨i, hfire⟩ hnormal
    exact hnormal (localRepair C i x) ⟨i, rfl, hfire⟩

def WorkConserving (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) : Prop :=
  ∀ n x, (∃ i, localRepair C i x ≠ x) →
    localRepair C (sigma n x) x ≠ x

theorem workConserving_descends_if_not_normal (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (hwork : WorkConserving C sigma)
    (n : Nat) (x : Records C)
    (hnormal : ¬ NormalForm C (adaptiveRun n sigma x)) :
    mismatchCount (adaptiveRun (n + 1) sigma x) <
      mismatchCount (adaptiveRun n sigma x) := by
  apply adaptiveRun_change_strict_rank
  rw [adaptiveRun_last_step]
  exact hwork n _ ((not_normal_iff_exists_firing C _).mp hnormal)

theorem workConserving_exists_normal_by_rank (C : OPHCarrier) :
    ∀ (x : Records C) (sigma : AdaptiveScheduler C),
      WorkConserving C sigma →
      ∃ N, N ≤ mismatchCount x ∧
        NormalForm C (adaptiveRun N sigma x) := by
  intro x sigma hwork
  generalize hd : mismatchCount x = d
  induction d using Nat.strong_induction_on generalizing x sigma with
  | h d ih =>
      by_cases hnormal : NormalForm C x
      · exact ⟨0, Nat.zero_le _, by simpa using hnormal⟩
      · have henabled := (not_normal_iff_exists_firing C x).mp hnormal
        have hfire : localRepair C (sigma 0 x) x ≠ x :=
          hwork 0 x henabled
        let x1 : Records C := localRepair C (sigma 0 x) x
        let sigma1 : AdaptiveScheduler C := fun k => sigma (k + 1)
        have hlt : mismatchCount x1 < mismatchCount x :=
          mismatchCount_localRepair_lt C (sigma 0 x) x hfire
        have hltD : mismatchCount x1 < d := hd ▸ hlt
        have hwork1 : WorkConserving C sigma1 := by
          intro n y hy
          exact hwork (n + 1) y hy
        obtain ⟨N, hN, hNnormal⟩ :=
          ih (mismatchCount x1) hltD x1 sigma1 hwork1 rfl
        refine ⟨N + 1, ?_, ?_⟩
        · omega
        · simpa [x1, sigma1, adaptiveRun_succ] using hNnormal

theorem adaptiveRun_of_normal (C : OPHCarrier) (n : Nat)
    (sigma : AdaptiveScheduler C) {x : Records C}
    (hnormal : NormalForm C x) :
    adaptiveRun n sigma x = x := by
  induction n generalizing sigma with
  | zero => rfl
  | succ n ih =>
      rw [adaptiveRun_succ]
      have hquiet : localRepair C (sigma 0 x) x = x := by
        by_contra hfire
        exact hnormal (localRepair C (sigma 0 x) x)
          ⟨sigma 0 x, rfl, hfire⟩
      rw [hquiet]
      exact ih (fun k => sigma (k + 1))

theorem normal_at_eventually_constant (C : OPHCarrier) (N : Nat)
    (sigma : AdaptiveScheduler C) (x : Records C)
    (hnormal : NormalForm C (adaptiveRun N sigma x)) :
    ∀ n, N ≤ n →
      adaptiveRun n sigma x = adaptiveRun N sigma x := by
  intro n hn
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
  rw [adaptiveRun_add]
  exact adaptiveRun_of_normal C k (fun j => sigma (N + j)) hnormal

theorem workConserving_eventually_normal (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (hwork : WorkConserving C sigma)
    (x : Records C) :
    ∃ N, N ≤ mismatchCount x ∧
      NormalForm C (adaptiveRun N sigma x) ∧
      ∀ n, N ≤ n →
        adaptiveRun n sigma x = adaptiveRun N sigma x := by
  obtain ⟨N, hN, hnormal⟩ :=
    workConserving_exists_normal_by_rank C x sigma hwork
  exact ⟨N, hN, hnormal,
    normal_at_eventually_constant C N sigma x hnormal⟩

theorem adaptiveRun_eventually_constant (C : OPHCarrier) :
    ∀ (x : Records C) (sigma : AdaptiveScheduler C),
      ∃ N, ∀ n, N ≤ n →
        adaptiveRun n sigma x = adaptiveRun N sigma x := by
  intro x sigma
  generalize hd : mismatchCount x = d
  induction d using Nat.strong_induction_on generalizing x sigma with
  | h d ih =>
      by_cases hchange :
          ∃ k, adaptiveRun (k + 1) sigma x ≠ adaptiveRun k sigma x
      · obtain ⟨k, hk⟩ := hchange
        let y : Records C := adaptiveRun (k + 1) sigma x
        let sigmaTail : AdaptiveScheduler C :=
          fun j => sigma ((k + 1) + j)
        have hylt : mismatchCount y < d := by
          have hstrict := adaptiveRun_change_strict_rank C k sigma x hk
          have hprefix := adaptiveRun_rank_le_initial C k sigma x
          dsimp [y]
          omega
        obtain ⟨N, hN⟩ :=
          ih (mismatchCount y) hylt y sigmaTail rfl
        refine ⟨(k + 1) + N, ?_⟩
        intro n hn
        obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le hn
        calc
          adaptiveRun (((k + 1) + N) + r) sigma x =
              adaptiveRun ((k + 1) + (N + r)) sigma x := by
                congr 1 <;> omega
          _ = adaptiveRun (N + r) sigmaTail y := by
                simpa [y, sigmaTail] using
                  adaptiveRun_add C (k + 1) (N + r) sigma x
          _ = adaptiveRun N sigmaTail y :=
                hN (N + r) (Nat.le_add_right N r)
          _ = adaptiveRun ((k + 1) + N) sigma x := by
                simpa [y, sigmaTail] using
                  (adaptiveRun_add C (k + 1) N sigma x).symm
      · have hsame : ∀ k,
            adaptiveRun (k + 1) sigma x = adaptiveRun k sigma x := by
          intro k
          exact not_ne_iff.mp (not_exists.mp hchange k)
        refine ⟨0, ?_⟩
        intro n hn
        induction n with
        | zero => rfl
        | succ n ihN =>
            rw [hsame n, ihN (Nat.zero_le n)]

def PathwiseWeakFair (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C) : Prop :=
  ∀ i N,
    (∀ n, N ≤ n →
      localRepair C i (adaptiveRun n sigma x) ≠
        adaptiveRun n sigma x) →
    ∃ m, N ≤ m ∧ sigma m (adaptiveRun m sigma x) = i

theorem workConserving_pathwiseWeakFair (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (hwork : WorkConserving C sigma)
    (x : Records C) : PathwiseWeakFair C sigma x := by
  intro i N hcontinuous
  obtain ⟨K, hKrank, hnormal, hconstant⟩ :=
    workConserving_eventually_normal C sigma hwork x
  let m := max N K
  have hmN : N ≤ m := le_max_left _ _
  have hmK : K ≤ m := le_max_right _ _
  have hfire := hcontinuous m hmN
  have hquiet : localRepair C i (adaptiveRun K sigma x) =
      adaptiveRun K sigma x := by
    by_contra hchange
    exact hnormal _ ⟨i, rfl, hchange⟩
  rw [hconstant m hmK] at hfire
  exact (hfire hquiet).elim

theorem pathwiseWeakFair_eventual_normalForm (C : OPHCarrier)
    (sigma : AdaptiveScheduler C) (x : Records C)
    (hfair : PathwiseWeakFair C sigma x) :
    ∃ N, NormalForm C (adaptiveRun N sigma x) ∧
      ∀ n, N ≤ n →
        adaptiveRun n sigma x = adaptiveRun N sigma x := by
  obtain ⟨N, hconst⟩ := adaptiveRun_eventually_constant C x sigma
  refine ⟨N, ?_, hconst⟩
  by_contra hnormal
  obtain ⟨i, hi⟩ :=
    (not_normal_iff_exists_firing C _).mp hnormal
  have hcontinuous : ∀ n, N ≤ n →
      localRepair C i (adaptiveRun n sigma x) ≠
        adaptiveRun n sigma x := by
    intro n hn
    rw [hconst n hn]
    exact hi
  obtain ⟨m, hm, hselected⟩ := hfair i N hcontinuous
  have hchange :
      adaptiveRun (m + 1) sigma x ≠ adaptiveRun m sigma x := by
    rw [adaptiveRun_last_step, hselected]
    exact hcontinuous m hm
  apply hchange
  rw [hconst m hm,
    hconst (m + 1) (Nat.le_trans hm (Nat.le_succ m))]

#print axioms adaptiveRun_eventually_constant
#print axioms workConserving_eventually_normal
#print axioms workConserving_pathwiseWeakFair
#print axioms pathwiseWeakFair_eventual_normalForm

end

end OPH.Execution
