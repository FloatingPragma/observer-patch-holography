import Mathlib.LinearAlgebra.Basis.Basic

/-!
In a basis of separately recorded one-dimensional lines, two coherent
linear actions that leave each line invariant have the same order amplitude.
-/

namespace EnergyOrder

theorem commute_of_fixed_lines
    {ι R V : Type*} [CommSemiring R] [AddCommMonoid V] [Module R V]
    (b : Module.Basis ι R V) (a c : V →ₗ[R] V)
    (u v : ι → R)
    (ha : ∀ i, a (b i) = u i • b i)
    (hc : ∀ i, c (b i) = v i • b i) :
    a.comp c = c.comp a := by
  apply b.ext
  intro i
  simp [LinearMap.comp_apply, ha, hc, smul_smul, mul_comm]

theorem commute_of_simple_energy
    {ι R V : Type*} [CommSemiring R] [AddCommMonoid V] [Module R V]
    (b : Module.Basis ι R V) (h a c : V →ₗ[R] V) (energy : ι → R)
    (hb : ∀ i, h (b i) = energy i • b i)
    (simple : ∀ i x, h x = energy i • x → ∃ k : R, x = k • b i)
    (ha : h.comp a = a.comp h) (hc : h.comp c = c.comp h) :
    a.comp c = c.comp a := by
  have preserve (f : V →ₗ[R] V) (hf : h.comp f = f.comp h) (i : ι) :
      ∃ k : R, f (b i) = k • b i := by
    apply simple i
    calc
      h (f (b i)) = f (h (b i)) := DFunLike.congr_fun hf (b i)
      _ = energy i • f (b i) := by rw [hb, map_smul]
  choose u hu using preserve a ha
  choose v hv using preserve c hc
  exact commute_of_fixed_lines b a c u v hu hv

end EnergyOrder
