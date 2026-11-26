import Mathlib.Data.List.Basic
import Mathlib.Logic.Relation

namespace OPH.SeamFlags

/-- The first load is supplied separately. Each subsequent load gives one seam read. -/
def flags (first : Bool) : List Bool → List Bool
  | [] => []
  | b :: rest => (first.xor b) :: flags b rest

/-- Recover the other loads along an observed path once a root load is specified. -/
def recover (first : Bool) : List Bool → List Bool
  | [] => []
  | b :: rest => (first.xor b) :: recover (first.xor b) rest

def arrival (first : Bool) : List Bool → Bool
  | [] => first
  | b :: rest => arrival (first.xor b) rest

theorem flags_recover (first : Bool) (bits : List Bool) :
    flags first (recover first bits) = bits := by
  induction bits generalizing first with
  | nil => rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, recover, ih]

theorem recover_flags (first : Bool) (tail : List Bool) :
    recover first (flags first tail) = tail := by
  induction tail generalizing first with
  | nil => rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, recover, ih]

/-- Replaying the complete measured path returns to its first load on a closed loop. -/
theorem arrival_close (first last : Bool) (tail : List Bool) :
    arrival first (flags first (tail ++ [last])) = last := by
  induction tail generalizing first with
  | nil => cases first <;> cases last <;> rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, arrival, ih]

/-- A connected receiver network determines a binary source word up to a common flip. -/
theorem connected_source_fiber {V : Type*} (edge : V → V → Prop) (root : V)
    (connected : ∀ v, Relation.ReflTransGen edge root v) (f g : V → Bool) :
    (∀ a b, edge a b → ((f a).xor (f b) = (g a).xor (g b))) ↔
      (f = g ∨ f = fun v => !(g v)) := by
  constructor
  · intro same
    have step : ∀ a b, edge a b → ((f a).xor (g a) = (f b).xor (g b)) := by
      intro a b hab
      have h : ∀ x y z w : Bool, x.xor y = z.xor w → x.xor z = y.xor w := by decide
      exact h _ _ _ _ (same a b hab)
    have invariant : ∀ v, (f root).xor (g root) = (f v).xor (g v) := by
      intro v
      induction connected v with
      | refl => rfl
      | tail _ hab ih => exact ih.trans (step _ _ hab)
    cases h : (f root).xor (g root) with
    | false =>
      left
      funext v
      have hv := (invariant v).symm.trans h
      exact (show ∀ x y : Bool, x.xor y = false → x = y by decide) _ _ hv
    | true =>
      right
      funext v
      have hv := (invariant v).symm.trans h
      exact (show ∀ x y : Bool, x.xor y = true → x = !y by decide) _ _ hv
  · rintro (rfl | rfl) a b _
    · rfl
    · exact (show ∀ x y : Bool, (!x).xor (!y) = x.xor y by decide) _ _

end OPH.SeamFlags
