import Mathlib

namespace OPH.A5SixAxes

/-! # The six-axis A5 action and the dimension-six branch module skeleton

Issue #604.  The compact-Lie trichotomy excludes the centre case with
semisimple part of dimension six by a module argument: the two `su(2)`
ideals would have to carry `1 + 5`, and the five-dimensional summand is
irreducible, so no decomposition into two three-dimensional invariant
summands exists.  This file supplies the concrete six-point `A5` action and
the module-theoretic half of that argument.

CONTENT.

* `t`, `s`, `L60`: the icosahedral six-axis action realized as
  `PSL(2, F5)` on the projective line `{0, 1, 2, 3, 4, ∞}` (index `5` is
  the point at infinity), with `t : z ↦ z + 1` and `s : z ↦ -1/z`.  `L60`
  lists all sixty group elements as explicit permutations of `Fin 6`.
* Kernel-`decide` facts: `L60` has sixty distinct members, contains the
  identity and both generators, is closed under inverse and product, and
  acts 2-transitively (`two_transitive`).
* `V5`: the sum-zero hyperplane of the permutation module `Fin 6 → ℚ`,
  invariant under the action (`V5_invariant`), of dimension five
  (`finrank_V5`).
* `no_three_dim_invariant`: GIVEN that `V5` is irreducible under the
  listed action, no invariant subspace of the permutation module has
  dimension three.  The proof is the dimension count
  `finrank (W ⊔ V5) + finrank (W ⊓ V5) = 3 + 5` with `finrank (W ⊔ V5) ≤ 6`,
  so `W ⊓ V5` is a nonzero invariant subspace of `V5`, hence all of `V5`,
  hence `finrank W ≥ 5`.  This excludes the two-summand decomposition
  demanded by the dimension-six branch, since a summand of a
  three-plus-three splitting is a three-dimensional invariant subspace.

BOUNDARY.  The irreducibility of `V5` under the listed action enters as an
explicit hypothesis (`hirr`).  Its proof (exclusion of invariant lines and
planes; the plane case is the order-five crystallographic restriction over
`ℚ`) is the open remainder of issue #604, recorded there.  Nothing in this
file asserts the hypothesis.  No physical content is at stake in this
receipt lane. -/

/-- Build a permutation of `Fin 6` from its value vector and inverse
vector; the two inverse laws are discharged by `decide` at each literal. -/
def perm (f g : Fin 6 → Fin 6)
    (h₁ : Function.LeftInverse g f := by decide)
    (h₂ : Function.RightInverse g f := by decide) : Equiv.Perm (Fin 6) :=
  ⟨f, g, h₁, h₂⟩

/-- Generator `t : z ↦ z + 1` on `P¹(F₅)`, the five-cycle fixing `∞`. -/
def t : Equiv.Perm (Fin 6) := perm ![1, 2, 3, 4, 0, 5] ![4, 0, 1, 2, 3, 5]

/-- Generator `s : z ↦ -1/z` on `P¹(F₅)`, the involution `(0 ∞)(1 4)`. -/
def s : Equiv.Perm (Fin 6) := perm ![5, 4, 2, 3, 1, 0] ![5, 4, 2, 3, 1, 0]

/-- All sixty elements of the generated copy of `A5 ≅ PSL(2, F5)` acting on
the six axes. -/
def L60 : List (Equiv.Perm (Fin 6)) := [
    perm ![0, 1, 2, 3, 4, 5] ![0, 1, 2, 3, 4, 5],
    perm ![0, 1, 4, 5, 2, 3] ![0, 1, 4, 5, 2, 3],
    perm ![0, 2, 1, 3, 5, 4] ![0, 2, 1, 3, 5, 4],
    perm ![0, 2, 5, 4, 1, 3] ![0, 4, 1, 5, 3, 2],
    perm ![0, 3, 4, 2, 5, 1] ![0, 5, 3, 1, 2, 4],
    perm ![0, 3, 5, 1, 4, 2] ![0, 3, 5, 1, 4, 2],
    perm ![0, 4, 1, 5, 3, 2] ![0, 2, 5, 4, 1, 3],
    perm ![0, 4, 3, 2, 1, 5] ![0, 4, 3, 2, 1, 5],
    perm ![0, 5, 2, 4, 3, 1] ![0, 5, 2, 4, 3, 1],
    perm ![0, 5, 3, 1, 2, 4] ![0, 3, 4, 2, 5, 1],
    perm ![1, 0, 2, 5, 4, 3] ![1, 0, 2, 5, 4, 3],
    perm ![1, 0, 4, 3, 2, 5] ![1, 0, 4, 3, 2, 5],
    perm ![1, 2, 0, 5, 3, 4] ![2, 0, 1, 4, 5, 3],
    perm ![1, 2, 3, 4, 0, 5] ![4, 0, 1, 2, 3, 5],
    perm ![1, 3, 2, 4, 5, 0] ![5, 0, 2, 1, 3, 4],
    perm ![1, 3, 5, 0, 2, 4] ![3, 0, 4, 1, 5, 2],
    perm ![1, 4, 0, 3, 5, 2] ![2, 0, 5, 3, 1, 4],
    perm ![1, 4, 5, 2, 0, 3] ![4, 0, 3, 5, 1, 2],
    perm ![1, 5, 3, 0, 4, 2] ![3, 0, 5, 2, 4, 1],
    perm ![1, 5, 4, 2, 3, 0] ![5, 0, 3, 4, 2, 1],
    perm ![2, 0, 1, 4, 5, 3] ![1, 2, 0, 5, 3, 4],
    perm ![2, 0, 5, 3, 1, 4] ![1, 4, 0, 3, 5, 2],
    perm ![2, 1, 0, 4, 3, 5] ![2, 1, 0, 4, 3, 5],
    perm ![2, 1, 3, 5, 0, 4] ![4, 1, 0, 2, 5, 3],
    perm ![2, 3, 1, 5, 4, 0] ![5, 2, 0, 1, 4, 3],
    perm ![2, 3, 4, 0, 1, 5] ![3, 4, 0, 1, 2, 5],
    perm ![2, 4, 3, 0, 5, 1] ![3, 5, 0, 2, 1, 4],
    perm ![2, 4, 5, 1, 3, 0] ![5, 3, 0, 4, 1, 2],
    perm ![2, 5, 0, 3, 4, 1] ![2, 5, 0, 3, 4, 1],
    perm ![2, 5, 4, 1, 0, 3] ![4, 3, 0, 5, 2, 1],
    perm ![3, 0, 4, 1, 5, 2] ![1, 3, 5, 0, 2, 4],
    perm ![3, 0, 5, 2, 4, 1] ![1, 5, 3, 0, 4, 2],
    perm ![3, 1, 2, 0, 5, 4] ![3, 1, 2, 0, 5, 4],
    perm ![3, 1, 5, 4, 2, 0] ![5, 1, 4, 0, 3, 2],
    perm ![3, 2, 1, 0, 4, 5] ![3, 2, 1, 0, 4, 5],
    perm ![3, 2, 4, 5, 1, 0] ![5, 4, 1, 0, 2, 3],
    perm ![3, 4, 0, 1, 2, 5] ![2, 3, 4, 0, 1, 5],
    perm ![3, 4, 2, 5, 0, 1] ![4, 5, 2, 0, 1, 3],
    perm ![3, 5, 0, 2, 1, 4] ![2, 4, 3, 0, 5, 1],
    perm ![3, 5, 1, 4, 0, 2] ![4, 2, 5, 0, 3, 1],
    perm ![4, 0, 1, 2, 3, 5] ![1, 2, 3, 4, 0, 5],
    perm ![4, 0, 3, 5, 1, 2] ![1, 4, 5, 2, 0, 3],
    perm ![4, 1, 0, 2, 5, 3] ![2, 1, 3, 5, 0, 4],
    perm ![4, 1, 5, 3, 0, 2] ![4, 1, 5, 3, 0, 2],
    perm ![4, 2, 3, 1, 5, 0] ![5, 3, 1, 2, 0, 4],
    perm ![4, 2, 5, 0, 3, 1] ![3, 5, 1, 4, 0, 2],
    perm ![4, 3, 0, 5, 2, 1] ![2, 5, 4, 1, 0, 3],
    perm ![4, 3, 2, 1, 0, 5] ![4, 3, 2, 1, 0, 5],
    perm ![4, 5, 1, 3, 2, 0] ![5, 2, 4, 3, 0, 1],
    perm ![4, 5, 2, 0, 1, 3] ![3, 4, 2, 5, 0, 1],
    perm ![5, 0, 2, 1, 3, 4] ![1, 3, 2, 4, 5, 0],
    perm ![5, 0, 3, 4, 2, 1] ![1, 5, 4, 2, 3, 0],
    perm ![5, 1, 3, 2, 4, 0] ![5, 1, 3, 2, 4, 0],
    perm ![5, 1, 4, 0, 3, 2] ![3, 1, 5, 4, 2, 0],
    perm ![5, 2, 0, 1, 4, 3] ![2, 3, 1, 5, 4, 0],
    perm ![5, 2, 4, 3, 0, 1] ![4, 5, 1, 3, 2, 0],
    perm ![5, 3, 0, 4, 1, 2] ![2, 4, 5, 1, 3, 0],
    perm ![5, 3, 1, 2, 0, 4] ![4, 2, 3, 1, 5, 0],
    perm ![5, 4, 1, 0, 2, 3] ![3, 2, 4, 5, 1, 0],
    perm ![5, 4, 2, 3, 1, 0] ![5, 4, 2, 3, 1, 0],  ]

set_option maxHeartbeats 2000000 in
theorem length_L60 : L60.length = 60 ∧ L60.Nodup := by decide

theorem one_mem : (1 : Equiv.Perm (Fin 6)) ∈ L60 := by decide

theorem t_mem : t ∈ L60 := by decide

theorem s_mem : s ∈ L60 := by decide

set_option maxHeartbeats 2000000 in
set_option maxRecDepth 16384 in
theorem inv_closed : ∀ g ∈ L60, g⁻¹ ∈ L60 := by decide

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
theorem mul_closed : ∀ g ∈ L60, ∀ h ∈ L60, g * h ∈ L60 := by decide

set_option maxHeartbeats 4000000 in
/-- The listed action is 2-transitive on the six axes. -/
theorem two_transitive :
    ∀ i j k l : Fin 6, i ≠ j → k ≠ l → ∃ g ∈ L60, g i = k ∧ g j = l := by
  decide

/-! ## The permutation module and the dimension count -/

/-- Coordinate-sum functional of the permutation module. -/
def sumF : (Fin 6 → ℚ) →ₗ[ℚ] ℚ where
  toFun v := ∑ i, v i
  map_add' v w := by simp [Finset.sum_add_distrib]
  map_smul' c v := by simp [Finset.mul_sum]

/-- The sum-zero hyperplane. -/
def V5 : Submodule ℚ (Fin 6 → ℚ) := LinearMap.ker sumF

/-- A permutation acts linearly on the module by right composition with its
inverse. -/
def pAct (g : Equiv.Perm (Fin 6)) : (Fin 6 → ℚ) →ₗ[ℚ] (Fin 6 → ℚ) where
  toFun v i := v (g⁻¹ i)
  map_add' v w := rfl
  map_smul' c v := rfl

/-- Invariance under every listed group element. -/
def InvariantUnder (W : Submodule ℚ (Fin 6 → ℚ)) : Prop :=
  ∀ g ∈ L60, ∀ w ∈ W, pAct g w ∈ W

theorem V5_invariant : InvariantUnder V5 := by
  intro g _ w hw
  rw [V5, LinearMap.mem_ker] at hw ⊢
  show ∑ i, w (g⁻¹ i) = 0
  rw [Equiv.sum_comp g⁻¹ w]
  exact hw

theorem finrank_V5 : Module.finrank ℚ V5 = 5 := by
  have hsurj : Function.Surjective sumF := by
    intro q
    refine ⟨fun i => if i = 0 then q else 0, ?_⟩
    simp [sumF, Finset.sum_ite_eq']
  have hrank := LinearMap.finrank_range_add_finrank_ker sumF
  rw [LinearMap.range_eq_top.mpr hsurj] at hrank
  simp only [finrank_top, Module.finrank_self, Module.finrank_pi,
    Fintype.card_fin] at hrank
  show Module.finrank ℚ (LinearMap.ker sumF) = 5
  omega

/-- The dimension-six branch module skeleton: if `V5` is irreducible under
the listed action, the permutation module has no three-dimensional
invariant subspace; in particular `1 + 5` admits no decomposition into two
three-dimensional invariant summands. -/
theorem no_three_dim_invariant
    (hirr : ∀ U : Submodule ℚ (Fin 6 → ℚ),
      U ≤ V5 → InvariantUnder U → U = ⊥ ∨ U = V5)
    (W : Submodule ℚ (Fin 6 → ℚ)) (hW : InvariantUnder W)
    (h3 : Module.finrank ℚ W = 3) : False := by
  have hinf_inv : InvariantUnder (W ⊓ V5) := by
    intro g hg w hw
    exact ⟨hW g hg w hw.1, V5_invariant g hg w hw.2⟩
  have hsum := Submodule.finrank_sup_add_finrank_inf_eq W V5
  have hsup : Module.finrank ℚ ↥(W ⊔ V5) ≤ 6 := by
    have h := Submodule.finrank_le (W ⊔ V5)
    simpa only [Module.finrank_pi, Fintype.card_fin] using h
  have hV5 := finrank_V5
  have h2 : 2 ≤ Module.finrank ℚ ↥(W ⊓ V5) := by omega
  rcases hirr (W ⊓ V5) inf_le_right hinf_inv with hbot | heq
  · rw [hbot] at h2
    simp at h2
  · have hle : V5 ≤ W := by
      rw [← heq]
      exact inf_le_left
    have hmono : Module.finrank ℚ V5 ≤ Module.finrank ℚ W :=
      Submodule.finrank_mono hle
    omega

end OPH.A5SixAxes

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5SixAxes.length_L60
#print axioms OPH.A5SixAxes.two_transitive
#print axioms OPH.A5SixAxes.mul_closed
#print axioms OPH.A5SixAxes.inv_closed
#print axioms OPH.A5SixAxes.V5_invariant
#print axioms OPH.A5SixAxes.finrank_V5
#print axioms OPH.A5SixAxes.no_three_dim_invariant
