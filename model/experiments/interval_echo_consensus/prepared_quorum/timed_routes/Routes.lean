import Mathlib

namespace OPH.QuorumRoutes

/-- Contracting the last edge of a five-cell ring leaves a cycle of four
validator regions. Every surviving pair has a route avoiding the silent one. -/
def adjacent (a b : Fin 4) : Prop := (a.val + 1) % 4 = b.val ∨ (b.val + 1) % 4 = a.val

instance (a b : Fin 4) : Decidable (adjacent a b) := by unfold adjacent; infer_instance

theorem every_triple : ∀ (f a b : Fin 4), a ≠ f → b ≠ f →
    a = b ∨ adjacent a b ∨ ∃ c : Fin 4, c ≠ f ∧ adjacent a c ∧ adjacent c b := by
  decide

/-- A finite alphabet of complete deadline-bounded receiver transcripts cannot
recover arbitrarily many view tags. This also permits variable-length frames. -/
theorem unbounded_view {Transcript : Type*} [Fintype Transcript]
    (send : ℕ → Transcript) (decode : Transcript → ℕ)
    (faithful : ∀ v, decode (send v) = v) : False := by
  let f : Fin (Fintype.card Transcript + 1) → Transcript := fun v => send v.val
  have inj : Function.Injective f := by
    intro a b h
    apply Fin.ext
    exact (faithful a.val).symm.trans ((congrArg decode h).trans (faithful b.val))
  have bound := Fintype.card_le_of_injective f inj
  rw [Fintype.card_fin] at bound
  omega

end OPH.QuorumRoutes
