import Mathlib.Algebra.BigOperators.Group.Finset.Basic

/-!
The declared two-register classical controller has a record-conditioned
inverse. Incoming rows select a single instrument, so a common fixed weight
remains fixed. Concrete native completeness, positivity, source preparation
and the eight-step quantum minorization are supplied by the exact Python
attachment, not by these generic finite-sum lemmas.
-/
namespace OPH.QFT.AdaptiveReference

def advance (m : Bool × Bool) (r : Bool) : Bool × Bool :=
  (Bool.xor m.2 r, m.1)

def retreat (m : Bool × Bool) (r : Bool) : Bool × Bool :=
  (m.2, Bool.xor m.1 r)

theorem retreat_advance (m : Bool × Bool) (r : Bool) :
    retreat (advance m r) r = m := by
  rcases m with ⟨c, d⟩
  cases c <;> cases d <;> cases r <;> rfl

theorem advance_retreat (m : Bool × Bool) (r : Bool) :
    advance (retreat m r) r = m := by
  rcases m with ⟨c, d⟩
  cases c <;> cases d <;> cases r <;> rfl

theorem incoming_selector (m : Bool × Bool) (r : Bool) :
    (retreat m r).1 = m.2 := rfl

variable {A R : Type*} [AddCommMonoid A] [Fintype R]

def incoming (instrument : Bool → R → A → A) (record : R → Bool)
    (weights : Bool × Bool → A) (m : Bool × Bool) : A :=
  ∑ r, instrument (retreat m (record r)).1 r (weights (retreat m (record r)))

theorem uniform_weight_fixed (instrument : Bool → R → A → A)
    (record : R → Bool) (x : A)
    (fixed : ∀ c, ∑ r, instrument c r x = x) (m : Bool × Bool) :
    incoming instrument record (fun _ => x) m = x := by
  exact fixed m.2

#print axioms retreat_advance
#print axioms advance_retreat
#print axioms incoming_selector
#print axioms uniform_weight_fixed

end OPH.QFT.AdaptiveReference
