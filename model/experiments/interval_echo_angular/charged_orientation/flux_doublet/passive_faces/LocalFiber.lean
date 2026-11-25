/- The exact native face-square projector and commuting-seam claims are checked on
the integer matrices in measure.py. These lemmas isolate the argument that extends
one repeatable face to every local fiber of a connected shell. -/
namespace ObserverPatchHolography.PassiveFaces

universe u

/-- A projector whose image of a sector covers that sector fixes it pointwise. -/
theorem fixed_of_onto {V : Type u} (P : V → V) (F : V → Prop)
    (projector : ∀ v, P (P v) = P v)
    (onto : ∀ v, F v → ∃ w, F w ∧ P w = v) :
    ∀ v, F v → P v = v := by
  intro v hv
  obtain ⟨w, _, hw⟩ := onto v hv
  rw [← hw]
  exact projector w

/-- Propagate the head constraint across a surjective local hopping block. -/
theorem transfer_fixed {V : Type u} (P T : V → V) (A B : V → Prop)
    (commute : ∀ v, P (T v) = T (P v))
    (fixed : ∀ v, A v → P v = v)
    (onto : ∀ v, B v → ∃ w, A w ∧ T w = v) :
    ∀ v, B v → P v = v := by
  intro v hv
  obtain ⟨w, hw, ht⟩ := onto v hv
  rw [← ht, commute, fixed w hw]

end ObserverPatchHolography.PassiveFaces
