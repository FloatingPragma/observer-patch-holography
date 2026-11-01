import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Group.Commute.Defs

/-!
# Finite unital instruments preserve commuting reference weights

These algebraic lemmas retain the full finite Kraus sum. The concrete
number-conserving Weyl/contact and native repair premises, source bounds,
faithfulness and feedback counterexample are checked in
`model/continuing_reference`; they are not encoded as a Lean matrix here.
No positivity, fixed-point classification or physical clock theorem is
claimed by the generic finite-sum lemmas alone.
-/
namespace OPH.QFT.ContinuingReference

variable {A I : Type*} [Semiring A] [Fintype I]

def instrument (left right : I → A) (x : A) : A :=
  ∑ i, left i * x * right i

theorem commuting_reference_fixed (left right : I → A) (x : A)
    (commuting : ∀ i, left i * x = x * left i)
    (unital : ∑ i, left i * right i = 1) :
    instrument left right x = x := by
  calc
    instrument left right x = ∑ i, x * (left i * right i) := by
      unfold instrument
      apply Finset.sum_congr rfl
      intro i _
      rw [commuting i, mul_assoc]
    _ = x * (∑ i, left i * right i) := (Finset.mul_sum _ _ _).symm
    _ = x := by rw [unital, mul_one]

theorem commuting_power_fixed (left right : I → A) (x : A)
    (commuting : ∀ i, Commute (left i) x)
    (unital : ∑ i, left i * right i = 1) (n : ℕ) :
    instrument left right (x^n) = x^n := by
  apply commuting_reference_fixed left right (x^n) _ unital
  intro i
  exact (commuting i).pow_right n

def history : List (A → A) → A → A
  | [], x => x
  | f :: fs, x => history fs (f x)

omit [Semiring A] in
theorem fixed_history (fs : List (A → A)) (x : A)
    (fixed : ∀ f ∈ fs, f x = x) : history fs x = x := by
  induction fs with
  | nil => rfl
  | cons f fs ih =>
    simp only [history]
    rw [fixed f (List.mem_cons_self)]
    apply ih
    intro g hg
    exact fixed g (List.mem_cons_of_mem f hg)

#print axioms commuting_reference_fixed
#print axioms commuting_power_fixed
#print axioms fixed_history

end OPH.QFT.ContinuingReference
