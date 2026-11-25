import Mathlib.Algebra.Group.Hom.Basic

/-!
# Recombination of reflected local exposure histories

The guarded native face has a reflection interchanging its two side ports. Each
reflected pulse intertwines with the direct pulse under that reflection. For
an odd orbital source, the two final arms therefore differ by minus reflection
for a list of any length. Their unnormalized bright and dark amplitudes select
opposite mirror parities. Finite rates of following the split energy line, its
geometric sign, and physical rotations require the source and spectrum.
-/

namespace OPH.SourceLaw.OrbitalMirrorCycle

variable {V : Type*} [AddCommGroup V]

/-- A pulse stores the actual operation of each arm and their reflection relation. -/
structure Pulse (mirror : V →+ V) where
  direct : V →+ V
  reflected : V →+ V
  intertwines : ∀ x, reflected (mirror x) = mirror (direct x)

/-- The first pulse in the list acts first on both arms. -/
def flight {mirror : V →+ V} (arm : Pulse mirror → V →+ V)
    (history : List (Pulse mirror)) (x : V) : V :=
  history.foldl (fun state pulse => arm pulse state) x

theorem reflected_history {mirror : V →+ V} (history : List (Pulse mirror)) (x : V) :
    flight Pulse.reflected history (-(mirror x)) =
      -(mirror (flight Pulse.direct history x)) := by
  induction history generalizing x with
  | nil => rfl
  | cons p ps ih =>
    change flight Pulse.reflected ps (p.reflected (-(mirror x))) =
      -(mirror (flight Pulse.direct ps (p.direct x)))
    rw [map_neg, p.intertwines, ih]

/-- Twice each output amplitude of a balanced recombiner for the odd source. -/
theorem output_parities {mirror : V →+ V} (history : List (Pulse mirror)) (x : V)
    (odd : mirror x = -x) :
    let l := flight Pulse.direct history x
    let r := flight Pulse.reflected history x
    l + r = l - mirror l ∧ l - r = l + mirror l := by
  have hr := reflected_history history x
  rw [odd, neg_neg] at hr
  dsimp only
  rw [hr, sub_neg_eq_add, sub_eq_add_neg]
  exact ⟨rfl, rfl⟩

end OPH.SourceLaw.OrbitalMirrorCycle
