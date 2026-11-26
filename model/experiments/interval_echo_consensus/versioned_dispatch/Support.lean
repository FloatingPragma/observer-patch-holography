import ObserverPatchHolography.Provenance.SemanticEventProvenance
import Mathlib.Tactic

namespace OPH.VersionedDispatch

/-- Agreement on the named registers forces the measured result at this launch. -/
def Sufficient {I V O : Type} [DecidableEq I] (f : (I → V) → O)
    (x : I → V) (s : Finset I) := ∀ z, (∀ i ∈ s, z i = x i) → f z = f x

/-- A sufficient certificate whose every individual input is pivotal is the unique least one. -/
theorem least_iff {I V O : Type} [DecidableEq I] {f : (I → V) → O}
    {x : I → V} {s : Finset I} (hs : Sufficient f x s)
    (piv : ∀ i ∈ s, ∃ z, (∀ j, j ≠ i → z j = x j) ∧ f z ≠ f x)
    (t : Finset I) : Sufficient f x t ↔ s ⊆ t := by
  constructor
  · intro ht i hi
    by_contra h
    obtain ⟨z, hz, diff⟩ := piv i hi
    exact diff (ht z (fun j hj => hz j (fun eq => h (eq ▸ hj))))
  · intro incl z hz
    exact hs z (fun i hi => hz i (incl hi))

abbrev Input := Fin 6 → Bool

def delivery (x : Input) : Option (Bool × Bool × Bool) :=
  if x 0 then some (x 1, x (if x 1 then 3 else 2), x 4) else none

def cert (x : Input) : Finset (Fin 6) :=
  if x 0 then {0, 1, (if x 1 then 3 else 2), 4} else {0}

theorem sufficient (x : Input) : Sufficient delivery x (cert x) := by
  intro z h
  cases hg : x 0 <;> cases hb : x 1 <;> simp [cert, hg, hb] at h <;>
    simp [delivery, *]

theorem pivotal (x : Input) (i : Fin 6) (hi : i ∈ cert x) :
    ∃ z, (∀ j, j ≠ i → z j = x j) ∧ delivery z ≠ delivery x := by
  refine ⟨Function.update x i (!(x i)), fun j hj => Function.update_of_ne hj _ _, ?_⟩
  fin_cases i <;> cases hg : x 0 <;> cases hb : x 1 <;>
    simp_all [cert, delivery, Function.update]

theorem exact_support (x : Input) (t : Finset (Fin 6)) :
    Sufficient delivery x t ↔ cert x ⊆ t :=
  least_iff (sufficient x) (pivotal x) t

/-- The continuation consumes the composite dispatch value as one register. -/
def response : Option (Bool × Bool × Bool) → Option Bool
  | none => none
  | some (_, payload, onward) => some (xor payload onward)

theorem continuation_pivotal (r : Option (Bool × Bool × Bool)) :
    ∃ s, response s ≠ response r := by
  cases r with
  | none => exact ⟨some (false, false, false), by decide⟩
  | some v => exact ⟨none, by simp [response]⟩

end OPH.VersionedDispatch
