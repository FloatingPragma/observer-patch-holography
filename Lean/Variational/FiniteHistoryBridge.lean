import InformationProjection.PathGibbs
import Variational.DiscreteNoether

namespace OPH.Variational

/-!
# The finite-history / real-variation bridge obstruction

The finite Gibbs package and the real discrete Euler--Lagrange package do not
compose merely because both use finite-length histories.  The Gibbs theorem
requires a finite state type, while the differential variational theorem uses
real-valued records and quantifies over every real single-site replacement.

This module makes the missing bridge precise.  For a fixed real-valued path
and site, the map sending a replacement value to the updated path is
injective.  Consequently no finite family of real-valued paths can be closed
under all of the single-site variations required by the current
Euler--Lagrange hypothesis.  A future composition theorem therefore needs a
genuine enrichment or limiting construction (for example, an embedding plus
a proved density/transfer theorem), not a silent identification of the two
path spaces.

This is an obstruction theorem about the present interfaces.  It does not
rule out a valid finite-to-continuous bridge, and it does not construct one.
-/

variable {N : ℕ}

/-- Updating a fixed path at a fixed site remembers the replacement value. -/
theorem updateAt_injective (γ : Fin (N + 1) → ℝ) (i : Fin (N + 1)) :
    Function.Injective (fun x : ℝ => Function.update γ i x) := by
  intro x y hxy
  have hi := congrFun hxy i
  simpa using hi

/-- **Finite-history variation obstruction.** Given any finite family of
real-valued paths, a fixed path, and a site, some real single-site variation
of that path lies outside the family.  Thus a finite Gibbs support cannot by
itself supply the universal real variations assumed by the differential
Euler--Lagrange theorem. -/
theorem exists_real_site_variation_outside
    (Γ : Finset (Fin (N + 1) → ℝ))
    (γ : Fin (N + 1) → ℝ) (i : Fin (N + 1)) :
    ∃ x : ℝ, Function.update γ i x ∉ Γ := by
  by_contra h
  simp only [not_exists, not_not] at h
  have hsubset : Set.range (fun x : ℝ => Function.update γ i x)
      ⊆ (Γ : Set (Fin (N + 1) → ℝ)) := by
    rintro _ ⟨x, rfl⟩
    exact h x
  exact (Set.infinite_range_of_injective (updateAt_injective γ i)).not_finite
    (Γ.finite_toSet.subset hsubset)

/-- Equivalent closure form: no finite real-path family contains every
single-site variation of even one fixed path. -/
theorem not_all_real_site_variations_mem
    (Γ : Finset (Fin (N + 1) → ℝ))
    (γ : Fin (N + 1) → ℝ) (i : Fin (N + 1)) :
    ¬ ∀ x : ℝ, Function.update γ i x ∈ Γ := by
  intro h
  obtain ⟨x, hx⟩ := exists_real_site_variation_outside Γ γ i
  exact hx (h x)

end OPH.Variational

#print axioms OPH.Variational.updateAt_injective
#print axioms OPH.Variational.exists_real_site_variation_outside
#print axioms OPH.Variational.not_all_real_site_variations_mem
