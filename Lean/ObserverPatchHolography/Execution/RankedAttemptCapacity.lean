import Mathlib

/-!
# Ranked attempt capacity

This module isolates the standard bounded-waste induction shared by adaptive
repair and fixed compiled-node attempts.  Domain-specific ranks, normal forms,
fairness, compiler structure, and sharpness remain in their owning modules.
-/

universe u v

namespace OPH.RankedAttempt

structure System (X : Type u) (A : Type v) where
  step : A → X → X
  rank : X → Nat
  change_rank_lt :
    ∀ (a : A) (x : X), step a x ≠ x → rank (step a x) < rank x

abbrev Scheduler (X : Type u) (A : Type v) :=
  Nat → X → A

def run {X : Type u} {A : Type v} (S : System X A) :
    Nat → Scheduler X A → X → X
  | 0, _sigma, x => x
  | n + 1, sigma, x =>
      let current := run S n sigma x
      S.step (sigma n current) current

theorem run_last_step {X : Type u} {A : Type v}
    (S : System X A) (n : Nat) (sigma : Scheduler X A) (x : X) :
    run S (n + 1) sigma x =
      S.step (sigma n (run S n sigma x)) (run S n sigma x) := rfl

theorem run_add {X : Type u} {A : Type v}
    (S : System X A) (m n : Nat) (sigma : Scheduler X A) (x : X) :
    run S (m + n) sigma x =
      run S n (fun k y => sigma (m + k) y) (run S m sigma x) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [Nat.add_succ, run_last_step, run_last_step, ih]

theorem run_change_rank_lt {X : Type u} {A : Type v}
    (S : System X A) (n : Nat) (sigma : Scheduler X A) (x : X)
    (hchange : run S (n + 1) sigma x ≠ run S n sigma x) :
    S.rank (run S (n + 1) sigma x) < S.rank (run S n sigma x) := by
  rw [run_last_step] at hchange ⊢
  exact S.change_rank_lt _ _ hchange

theorem run_step_rank_le {X : Type u} {A : Type v}
    (S : System X A) (n : Nat) (sigma : Scheduler X A) (x : X) :
    S.rank (run S (n + 1) sigma x) ≤ S.rank (run S n sigma x) := by
  by_cases hchange : run S (n + 1) sigma x ≠ run S n sigma x
  · exact Nat.le_of_lt (run_change_rank_lt S n sigma x hchange)
  · rw [not_ne_iff.mp hchange]

theorem run_rank_le_initial {X : Type u} {A : Type v}
    (S : System X A) (n : Nat) (sigma : Scheduler X A) (x : X) :
    S.rank (run S n sigma x) ≤ S.rank x := by
  induction n with
  | zero => rfl
  | succ n ih => exact le_trans (run_step_rank_le S n sigma x) ih

def Quiescent {X : Type u} {A : Type v}
    (S : System X A) (x : X) : Prop :=
  ∀ a, S.step a x = x

def BoundedWaste {X : Type u} {A : Type v}
    (S : System X A) (waste : Nat)
    (sigma : Scheduler X A) : Prop :=
  ∀ start x, ¬ Quiescent S x →
    ∃ offset, offset ≤ waste ∧
      run S (offset + 1) (fun k y => sigma (start + k) y) x ≠
        run S offset (fun k y => sigma (start + k) y) x

theorem run_eq_of_quiescent {X : Type u} {A : Type v}
    (S : System X A) (x : X) (sigma : Scheduler X A)
    (hquiet : Quiescent S x) :
    ∀ n, run S n sigma x = x := by
  intro n
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [run_last_step, ih]
      exact hquiet _

theorem boundedWaste_exists_quiescent_by_rank
    {X : Type u} {A : Type v}
    (S : System X A) (waste : Nat) :
    ∀ (x : X) (sigma : Scheduler X A),
      BoundedWaste S waste sigma →
      ∃ N, N ≤ (waste + 1) * S.rank x ∧
        Quiescent S (run S N sigma x) := by
  intro x sigma hbounded
  generalize hd : S.rank x = d
  induction d using Nat.strong_induction_on generalizing x sigma with
  | h d ih =>
      by_cases hquiet : Quiescent S x
      · exact ⟨0, Nat.zero_le _, by simpa using hquiet⟩
      · obtain ⟨offset, hoffset, hchange⟩ := hbounded 0 x hquiet
        simp at hchange
        let y : X := run S (offset + 1) sigma x
        let sigmaTail : Scheduler X A :=
          fun k z => sigma ((offset + 1) + k) z
        have hylt : S.rank y < d := by
          have hstrict := run_change_rank_lt S offset sigma x hchange
          have hprefix := run_rank_le_initial S offset sigma x
          dsimp [y]
          omega
        have hboundedTail : BoundedWaste S waste sigmaTail := by
          intro start z hnotQuiet
          obtain ⟨later, hlater, hchanges⟩ :=
            hbounded ((offset + 1) + start) z hnotQuiet
          refine ⟨later, hlater, ?_⟩
          simpa [sigmaTail, Nat.add_assoc] using hchanges
        obtain ⟨N, hN, hquietN⟩ :=
          ih (S.rank y) hylt y sigmaTail hboundedTail rfl
        refine ⟨(offset + 1) + N, ?_, ?_⟩
        · have hrank : S.rank y + 1 ≤ d := Nat.succ_le_iff.mpr hylt
          have hmul := Nat.mul_le_mul_left (waste + 1) hrank
          rw [Nat.mul_add, Nat.mul_one] at hmul
          omega
        · rw [run_add]
          simpa [y, sigmaTail] using hquietN

theorem boundedWaste_eventually_quiescent
    {X : Type u} {A : Type v}
    (S : System X A) (waste : Nat)
    (x : X) (sigma : Scheduler X A)
    (hbounded : BoundedWaste S waste sigma) :
    ∃ N,
      N ≤ (waste + 1) * S.rank x ∧
      Quiescent S (run S N sigma x) ∧
      ∀ n, N ≤ n → run S n sigma x = run S N sigma x := by
  obtain ⟨N, hN, hquiet⟩ :=
    boundedWaste_exists_quiescent_by_rank S waste x sigma hbounded
  refine ⟨N, hN, hquiet, ?_⟩
  intro n hn
  obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le hn
  rw [run_add]
  exact run_eq_of_quiescent S (run S N sigma x)
    (fun k y => sigma (N + k) y) hquiet r

#print axioms boundedWaste_exists_quiescent_by_rank
#print axioms boundedWaste_eventually_quiescent

end OPH.RankedAttempt
