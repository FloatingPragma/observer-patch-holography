import ObservableNormalForms.ObserverConfluence

/-! The quotient observed by a public priority sweep with total load two.
The native column and flight certificate supplying the map to the three
classes is checked in check.py. -/

namespace TwoLoadMark

open ObservableNormalForms Relation

inductive State
  | blankBinary | blankDouble | markedBinary
  deriving DecidableEq

open State

def step (x y : State) : Prop := x = blankDouble ∧ y = markedBinary
def consistent : Set State := {x | x ≠ blankDouble}
def boundary (x : State) : Bool := x != blankBinary
def rank : State → Nat
  | blankDouble => 1
  | _ => 0

theorem preserves : ObservationPreserving step boundary := by
  rintro x y ⟨rfl, rfl⟩
  rfl

theorem complete : CompleteFor step consistent := by
  intro x
  cases x <;> simp [IsNormalForm, step, consistent]

theorem descent {x y : State} (h : step x y) : rank y < rank x := by
  rcases h with ⟨rfl, rfl⟩
  decide

theorem diamond {x y z : State} (h : step x y) (k : step x z) :
    ∃ w, ReflTransGen step y w ∧ ReflTransGen step z w := by
  have eq : y = z := h.2.trans k.2.symm
  subst z
  exact ⟨y, .refl, .refl⟩

theorem identifies : BoundaryIdentifiesModulo consistent boundary Eq := by
  intro x y hx hy h
  cases x <;> cases y <;> simp_all [consistent, boundary]

theorem normalizes (x : State) : WeaklyNormalizing step x := by
  cases x with
  | blankBinary => exact ⟨_, .refl, (complete _).mpr (by simp [consistent])⟩
  | markedBinary => exact ⟨_, .refl, (complete _).mpr (by simp [consistent])⟩
  | blankDouble =>
    exact ⟨markedBinary, .single ⟨rfl, rfl⟩, (complete _).mpr (by simp [consistent])⟩

theorem uniqueEndpoint : ObserverEndpointUniqueModulo step boundary Eq :=
  (boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo preserves complete).mp identifies

end TwoLoadMark
