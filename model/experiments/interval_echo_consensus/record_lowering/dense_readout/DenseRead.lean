import Mathlib.Data.Fintype.Card
import Mathlib.Tactic

/-! Named bits, immutable old records and local lossless decoding at a fixed calendar. -/

namespace OPH.DenseRead

/-- An intervention is a local exchange of the balanced live rails. -/
def edit {A : Type*} [DecidableEq A] (bits : A → Bool) (i : A) (v : Bool) :
    A → Bool := Function.update bits i v

/-- The full prior cache is a separate register. -/
def boundary {A : Type*} [DecidableEq A] (x : (A → Bool) × (A → Bool))
    (i : A) (v : Bool) := (edit x.1 i v, x.2)

theorem boundary_local {A : Type*} [DecidableEq A]
    (x : (A → Bool) × (A → Bool)) (i : A) (v : Bool) :
    (boundary x i v).1 i = v ∧ (boundary x i v).2 = x.2 ∧
    ∀ j ≠ i, (boundary x i v).1 j = x.1 j := by
  refine ⟨by simp [boundary, edit], rfl, ?_⟩
  intro j hj
  simp [boundary, edit, Function.update, hj]

/-- Arbitrary fixed timing and control metadata can be included in each tape. -/
theorem local_injective {A Meta : Type*} (m : Meta) :
    Function.Injective (fun x : A → Bool => (m, x)) := by
  intro a b h
  exact congrArg Prod.snd h

theorem stale_collision {A Meta : Type*} [DecidableEq A] (i : A)
    (x : (A → Bool) × (A → Bool)) (m : Meta) :
    boundary x i false ≠ boundary x i true ∧
    (m, (boundary x i false).2) = (m, (boundary x i true).2) := by
  refine ⟨?_, rfl⟩
  intro heq
  have h := congrArg (fun r : (A → Bool) × (A → Bool) => r.1 i) heq
  simp [boundary, edit] at h

/-- Any lossless end record, including all control and timing, needs this capacity. -/
theorem transcript_card {k : ℕ} {Tape : Type*} [Fintype Tape]
    (emit : (Fin k → Bool) → Tape) (read : Tape → (Fin k → Bool))
    (recover : Function.LeftInverse read emit) : 2 ^ k ≤ Fintype.card Tape := by
  simpa using Fintype.card_le_of_injective emit recover.injective

/-- A farther in-menu arrival and an earlier out-of-menu arrival exclude any deadline. -/
theorem deadline_inversion (inside outside : ℕ) (early : outside ≤ inside) :
    ¬ ∃ deadline : ℕ, inside ≤ deadline ∧ deadline < outside := by omega

end OPH.DenseRead
