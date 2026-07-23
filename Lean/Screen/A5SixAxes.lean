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
* `V5_irreducible`: the irreducibility hypothesis is discharged.  The
  engine is the kernel-checked sharp fiber count `count_fibers` of the
  listed action (stabilizer size ten, sharp 2-transitive fiber size two);
  averaging any nonzero sum-zero vector over the stabilizer cosets of a
  coordinate where it does not vanish produces every difference vector
  `e k - e l` inside an invariant subspace, and these span `V5`.
* `no_three_dim_invariant_unconditional`, `no_three_plus_three_split`:
  the dimension-six branch exclusion with the hypothesis discharged.

BOUNDARY.  The conditional form `no_three_dim_invariant` is retained so the
module argument and the irreducibility input remain separately auditable.
No physical content is at stake in this receipt lane. -/

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
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

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

/-! ## Irreducibility of `V5` over `ℚ`

The declared hypothesis of `no_three_dim_invariant` is discharged here.
The engine is the sharp fiber count of the listed action at the base point
`0`: the sixty elements distribute over the joint conditions `g 0 = k`,
`g j = i` with multiplicities forced by sharp 2-transitivity — ten on the
compatible diagonal (`i = k`, `j = 0`), zero on the two mixed cases, two
otherwise.  A nonzero vector of an invariant subspace is first transported
by a listed element so that it does not vanish at `0` (`exists_to_zero`);
averaging it over a stabilizer coset then yields every difference vector
`e k - e l` inside the subspace, and these span `V5`.  Fixing the base
point and using only forward applications keeps the kernel obligations
small. -/

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Sharp fiber count of the listed action at base point `0`,
kernel-checked over all `6³` index choices. -/
theorem count_fibers :
    ∀ k i j : Fin 6,
      L60.countP (fun g => decide (g 0 = k ∧ g j = i)) =
        if i = k then (if j = 0 then 10 else 0)
        else (if j = 0 then 0 else 2) := by
  decide

/-- Every point is carried to the base point `0` by a listed element. -/
theorem exists_to_zero : ∀ a : Fin 6, ∃ g ∈ L60, g a = 0 := by decide

/-- Stabilizer-coset average of `u` with target `k`, written as a single
list sum with indicator terms so that membership in an invariant submodule
is termwise. -/
def avg (u : Fin 6 → ℚ) (k : Fin 6) : Fin 6 → ℚ :=
  (L60.map (fun g => if g 0 = k then pAct g u else 0)).sum

theorem list_sum_mem (U : Submodule ℚ (Fin 6 → ℚ)) :
    ∀ l : List (Fin 6 → ℚ), (∀ x ∈ l, x ∈ U) → l.sum ∈ U := by
  intro l
  induction l with
  | nil => intro _; rw [List.sum_nil]; exact U.zero_mem
  | cons x t ih =>
      intro h
      rw [List.sum_cons]
      exact U.add_mem (h x (by simp)) (ih fun y hy => h y (by simp [hy]))

theorem avg_mem (U : Submodule ℚ (Fin 6 → ℚ)) (hinv : InvariantUnder U)
    {u : Fin 6 → ℚ} (hu : u ∈ U) (k : Fin 6) : avg u k ∈ U := by
  refine list_sum_mem U _ ?_
  intro x hx
  obtain ⟨g, hg, rfl⟩ := List.mem_map.mp hx
  split_ifs
  · exact hinv g hg u hu
  · exact U.zero_mem

/-- Pointwise value of an indicator list sum, expressed through the fiber
counts of the underlying list. -/
theorem sum_indicator_apply (S : List (Equiv.Perm (Fin 6))) (u : Fin 6 → ℚ)
    (k i : Fin 6) :
    (S.map (fun g => if g 0 = k then pAct g u else 0)).sum i =
      ∑ j : Fin 6,
        (S.countP (fun g => decide (g 0 = k ∧ g j = i)) : ℚ) * u j := by
  induction S with
  | nil => simp
  | cons g S ih =>
      rw [List.map_cons, List.sum_cons, Pi.add_apply, ih]
      simp only [List.countP_cons]
      push_cast
      rw [Finset.sum_congr rfl
        (fun j _ => add_mul (S.countP (fun g => decide (g 0 = k ∧ g j = i)) : ℚ)
          _ (u j)), Finset.sum_add_distrib]
      have hterm : ((if g 0 = k then pAct g u else 0) : Fin 6 → ℚ) i =
          ∑ j : Fin 6,
            (if decide (g 0 = k ∧ g j = i) = true then (1 : ℚ) else 0) * u j := by
        simp only [decide_eq_true_eq]
        by_cases h : g 0 = k
        · have hcond : ∀ j : Fin 6,
              (if g 0 = k ∧ g j = i then (1 : ℚ) else 0) * u j =
                if g⁻¹ i = j then u j else 0 := by
            intro j
            by_cases hj : g⁻¹ i = j
            · subst hj
              rw [if_pos ⟨h, Equiv.apply_symm_apply g i⟩, if_pos rfl, one_mul]
            · have hne : g j ≠ i := fun he =>
                hj (by rw [← he]; exact Equiv.symm_apply_apply g j)
              rw [if_neg fun hc => hne hc.2, if_neg hj, zero_mul]
          rw [Finset.sum_congr rfl fun j _ => hcond j, Finset.sum_ite_eq]
          simp [h]
          rfl
        · simp [h]
      rw [hterm]
      ring

/-- Value of the stabilizer-coset average on a sum-zero vector. -/
theorem avg_apply {u : Fin 6 → ℚ} (hsum : ∑ j : Fin 6, u j = 0) (k i : Fin 6) :
    avg u k i = if i = k then 10 * u 0 else -2 * u 0 := by
  rw [avg, sum_indicator_apply]
  simp only [count_fibers k i]
  by_cases h : i = k
  · simp only [h, if_true]
    have : ∀ j : Fin 6,
        (((if j = 0 then 10 else 0 : ℕ)) : ℚ) * u j =
          if j = 0 then 10 * u j else 0 := by
      intro j; by_cases hj : j = 0 <;> simp [hj]
    rw [Finset.sum_congr rfl fun j _ => this j, Finset.sum_ite_eq']
    simp
  · simp only [h, if_false]
    have : ∀ j : Fin 6,
        (((if j = 0 then 0 else 2 : ℕ)) : ℚ) * u j =
          2 * u j - (if j = 0 then 2 * u j else 0) := by
      intro j; by_cases hj : j = 0 <;> simp [hj]
    rw [Finset.sum_congr rfl fun j _ => this j, Finset.sum_sub_distrib,
      ← Finset.mul_sum, hsum, Finset.sum_ite_eq']
    simp

/-- The explicit difference vector `e k - e l`. -/
def dvec (k l : Fin 6) : Fin 6 → ℚ :=
  fun i => (if i = k then 1 else 0) - (if i = l then 1 else 0)

theorem avg_sub_avg {u : Fin 6 → ℚ} (hsum : ∑ j : Fin 6, u j = 0)
    (k l : Fin 6) :
    avg u k - avg u l = (12 * u 0) • dvec k l := by
  funext i
  rw [Pi.sub_apply, avg_apply hsum, avg_apply hsum, Pi.smul_apply, dvec,
    smul_eq_mul]
  split_ifs <;> ring

theorem dvec_mem (U : Submodule ℚ (Fin 6 → ℚ)) (hU : U ≤ V5)
    (hinv : InvariantUnder U) {u : Fin 6 → ℚ} (hu : u ∈ U)
    (h0 : u 0 ≠ 0) (k l : Fin 6) : dvec k l ∈ U := by
  have hsum : ∑ j : Fin 6, u j = 0 := by
    have h := hU hu
    rw [V5, LinearMap.mem_ker] at h
    exact h
  have h12 : (12 : ℚ) * u 0 ≠ 0 := mul_ne_zero (by norm_num) h0
  have hmem : avg u k - avg u l ∈ U :=
    U.sub_mem (avg_mem U hinv hu k) (avg_mem U hinv hu l)
  rw [avg_sub_avg hsum] at hmem
  have hsc := U.smul_mem (12 * u 0)⁻¹ hmem
  rwa [smul_smul, inv_mul_cancel₀ h12, one_smul] at hsc

/-- `V5` is irreducible over `ℚ` under the listed action: the declared
hypothesis of `no_three_dim_invariant`, discharged. -/
theorem V5_irreducible :
    ∀ U : Submodule ℚ (Fin 6 → ℚ),
      U ≤ V5 → InvariantUnder U → U = ⊥ ∨ U = V5 := by
  intro U hU hinv
  rcases eq_or_ne U ⊥ with hbot | hbot
  · exact Or.inl hbot
  · refine Or.inr (le_antisymm hU ?_)
    obtain ⟨w, hw, hw0⟩ := (Submodule.ne_bot_iff U).mp hbot
    obtain ⟨a, hwa⟩ : ∃ a, w a ≠ 0 := Function.ne_iff.mp hw0
    obtain ⟨g, hg, hga⟩ := exists_to_zero a
    set u : Fin 6 → ℚ := pAct g w with hu_def
    have hu : u ∈ U := hinv g hg w hw
    have ha : u 0 ≠ 0 := by
      have hgi : g⁻¹ 0 = a := by
        rw [← hga]
        exact Equiv.symm_apply_apply g a
      show w (g⁻¹ 0) ≠ 0
      rw [hgi]
      exact hwa
    intro v hv
    have hvs : ∑ j : Fin 6, v j = 0 := by
      rw [V5, LinearMap.mem_ker] at hv
      exact hv
    have hrepr : v = ∑ j : Fin 6, v j • dvec j 0 := by
      funext i
      rw [Finset.sum_apply]
      simp only [Pi.smul_apply, dvec, smul_eq_mul]
      have expand : ∀ j : Fin 6,
          v j * ((if i = j then (1 : ℚ) else 0) - (if i = 0 then 1 else 0)) =
            (if i = j then v j else 0) - (if i = 0 then v j else 0) := by
        intro j; split_ifs <;> ring
      rw [Finset.sum_congr rfl fun j _ => expand j, Finset.sum_sub_distrib,
        Finset.sum_ite_eq]
      by_cases h0 : i = 0 <;> simp [h0, hvs]
    rw [hrepr]
    exact Submodule.sum_mem U fun j _ =>
      U.smul_mem (v j) (dvec_mem U hU hinv hu ha j 0)

/-- Unconditional dimension-six branch exclusion: the permutation module of
the six-axis action has no three-dimensional invariant subspace. -/
theorem no_three_dim_invariant_unconditional
    (W : Submodule ℚ (Fin 6 → ℚ)) (hW : InvariantUnder W)
    (h3 : Module.finrank ℚ W = 3) : False :=
  no_three_dim_invariant V5_irreducible W hW h3

/-- Paper-facing form of the dimension-six branch: `1 ⊕ 5` admits no
decomposition into two three-dimensional invariant summands. -/
theorem no_three_plus_three_split :
    ¬ ∃ W₁ W₂ : Submodule ℚ (Fin 6 → ℚ),
        InvariantUnder W₁ ∧ InvariantUnder W₂ ∧
        Module.finrank ℚ W₁ = 3 ∧ Module.finrank ℚ W₂ = 3 ∧
        W₁ ⊓ W₂ = ⊥ ∧ W₁ ⊔ W₂ = ⊤ := by
  rintro ⟨W₁, _W₂, h₁, _, d₁, _⟩
  exact no_three_dim_invariant_unconditional W₁ h₁ d₁

end OPH.A5SixAxes

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5SixAxes.length_L60
#print axioms OPH.A5SixAxes.two_transitive
#print axioms OPH.A5SixAxes.mul_closed
#print axioms OPH.A5SixAxes.inv_closed
#print axioms OPH.A5SixAxes.V5_invariant
#print axioms OPH.A5SixAxes.finrank_V5
#print axioms OPH.A5SixAxes.no_three_dim_invariant
#print axioms OPH.A5SixAxes.count_fibers
#print axioms OPH.A5SixAxes.V5_irreducible
#print axioms OPH.A5SixAxes.no_three_dim_invariant_unconditional
#print axioms OPH.A5SixAxes.no_three_plus_three_split
