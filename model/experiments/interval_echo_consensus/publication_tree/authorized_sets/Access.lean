import Mathlib.Data.Set.Basic
import Mathlib.Tactic.ByContra

namespace OPH.TwelveReaderCheckpoint

variable {X P V : Type*}

/-- The authorized group collects the named readers' records for the same source. -/
def observe (record : X → P → V) (group : Set P) (x : X) : group → V :=
  fun i => record x i.val

/-- Exact recoverability on an enrolled code depends on the reads of that group. -/
theorem reads_separate_iff (record : X → P → V) (group : Set P) :
    Function.Injective (observe record group) ↔
      ∀ x y : X, x ≠ y → ∃ i ∈ group, record x i ≠ record y i := by
  classical
  constructor
  · intro inj x y hxy
    by_contra! eq
    exact hxy (inj (funext fun i => eq i.val i.property))
  · intro sep x y eq
    by_contra hxy
    obtain ⟨i, hi, h⟩ := sep x y hxy
    exact h (congrFun eq ⟨i, hi⟩)

/-- When every pair differs at two addresses, losing one named address is correctable. -/
theorem omit_one (record : X → P → V)
    (separate : ∀ x y : X, x ≠ y → ∃ i j, i ≠ j ∧ record x i ≠ record y i ∧
                                                            record x j ≠ record y j)
    (lost : P) : Function.Injective (observe record ({lost}ᶜ)) := by
  apply (reads_separate_iff ..).2
  intro x y h
  obtain ⟨i, j, hij, hi, hj⟩ := separate x y h
  by_cases hil : i = lost
  · subst lost
    exact ⟨j, fun eq => hij (show j = i from eq).symm, hj⟩
  · exact ⟨i, fun eq => hil (show i = lost from eq), hi⟩

/-- A source-to-record map tested on each setting intertwines arbitrary setting words. -/
theorem read_word {S A B : Type*} (run : S → A → A) (action : S → B → B)
    (read : A → B) (step : ∀ s a, read (run s a) = action s (read a))
    (word : List S) (a : A) :
    read (word.foldl (fun x s => run s x) a) =
      word.foldl (fun x s => action s x) (read a) := by
  induction word generalizing a with
  | nil => rfl
  | cons s rest ih => simpa [List.foldl, step] using ih (run s a)

end OPH.TwelveReaderCheckpoint
