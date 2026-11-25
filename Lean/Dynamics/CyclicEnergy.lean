import Mathlib.Algebra.Module.LinearMap.Polynomial

/-! # Symmetries of a source with a cyclic energy operator

The native charged path supplies the cyclicity hypothesis by the nonzero
diagonal of its integer endpoint Krylov matrix. The conclusion concerns
maps preserving the energy operator of that same reachable sector.
-/

namespace OPH.Dynamics.CyclicEnergy

variable {K V : Type*} [Field K] [AddCommGroup V] [Module K V]

/-- A finite linear combination of powers of the energy operator. -/
def combination (H : Module.End K V) (n : ℕ) (c : ℕ → K) : Module.End K V :=
  ∑ i ∈ Finset.range n, c i • H ^ i

theorem commute_combination {H A : Module.End K V} (ha : Commute A H)
    (n : ℕ) (c : ℕ → K) : Commute A (combination H n c) := by
  apply Commute.sum_right
  intro i _
  exact (ha.pow_right i).smul_right (c i)

/-- The endpoint fixes every energy symmetry on every vector in its cyclic span. -/
theorem cyclic_centralizer_commutes (H A B : Module.End K V) (v : V)
    (cyclic : ∀ x : V, ∃ (n : ℕ) (c : ℕ → K), combination H n c v = x)
    (ha : Commute A H) (hb : Commute B H) : Commute A B := by
  obtain ⟨n, c, hv⟩ := cyclic (A v)
  have hp : Commute (combination H n c) H :=
    (commute_combination (Commute.refl H) n c).symm
  have heq : A = combination H n c := by
    apply LinearMap.ext
    intro x
    obtain ⟨m, d, hx⟩ := cyclic x
    rw [← hx]
    change (A * combination H m d) v = (combination H n c * combination H m d) v
    rw [(commute_combination ha m d).eq, (commute_combination hp m d).eq]
    change combination H m d (A v) = combination H m d (combination H n c v)
    rw [hv]
  rw [heq]
  exact (commute_combination hb n c).symm

end OPH.Dynamics.CyclicEnergy
