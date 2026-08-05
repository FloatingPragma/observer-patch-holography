import Dynamics.PublicMarkov
import Mathlib.Topology.LocallyConstant.Basic

/-!
# Automorphisms of a finite public record algebra

This file closes the algebraic part of the public reversibility boundary.
Every complex star-algebra automorphism of a finite function algebra is
pullback by one permutation of the labels.  The proof is elementary and
finite: automorphisms preserve the orthogonal idempotent record basis, and at
each output label exactly one basis idempotent has value one.

The parameter of a family of automorphisms is not identified with physical
time here.  Pointwise operator continuity is stated explicitly and is used to
turn the classified labels into locally constant maps; it supplies no clock,
source law, or physical evolution.
-/

namespace OPH.Dynamics

open scoped BigOperators

universe u

variable {ι : Type u} [Fintype ι] [DecidableEq ι]

private theorem complex_idempotent_eq_zero_or_one {z : ℂ}
    (hz : z * z = z) : z = 0 ∨ z = 1 := by
  have hfactor : z * (z - 1) = 0 := by
    calc
      z * (z - 1) = z * z - z := by ring
      _ = 0 := by rw [hz]; simp
  rcases mul_eq_zero.mp hfactor with hzero | hone
  · exact Or.inl hzero
  · exact Or.inr (sub_eq_zero.mp hone)

omit [Fintype ι] in
private theorem recordBasis_mul_self (j : ι) :
    recordBasis j * recordBasis j = recordBasis j := by
  funext x
  by_cases h : x = j <;> simp [recordBasis, h]

omit [Fintype ι] in
private theorem recordBasis_mul_eq_zero {j k : ι} (hjk : j ≠ k) :
    recordBasis j * recordBasis k = 0 := by
  funext x
  by_cases hj : x = j
  · subst x
    simp [recordBasis, hjk]
  · simp [recordBasis, hj]

omit [Fintype ι] in
/-- Images of record-basis idempotents remain pointwise zero-one valued. -/
theorem starAutomorphism_recordBasis_eq_zero_or_one
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x y : ι) :
    F (recordBasis y) x = 0 ∨ F (recordBasis y) x = 1 := by
  apply complex_idempotent_eq_zero_or_one
  have hmap :
      F (recordBasis y) * F (recordBasis y) = F (recordBasis y) := by
    rw [← map_mul, recordBasis_mul_self]
  exact congrFun hmap x

/-- At every output label, exactly one input record idempotent is selected by
the automorphism. -/
theorem existsUnique_recordBasis_image_eq_one
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x : ι) :
    ∃! y : ι, F (recordBasis y) x = 1 := by
  classical
  have hsum : ∑ y : ι, F (recordBasis y) x = 1 := by
    calc
      (∑ y : ι, F (recordBasis y) x) =
          (∑ y : ι, F (recordBasis y)) x := by
            rw [Finset.sum_apply]
      _ = F (∑ y : ι, recordBasis y) x := by rw [map_sum]
      _ = F 1 x := by rw [sum_recordBasis]
      _ = 1 := by simp
  have hexists : ∃ y : ι, F (recordBasis y) x = 1 := by
    by_contra hnone
    have hzero : ∀ y : ι, F (recordBasis y) x = 0 := by
      intro y
      rcases starAutomorphism_recordBasis_eq_zero_or_one F x y with hy | hy
      · exact hy
      · exact False.elim (hnone ⟨y, hy⟩)
    have : (0 : ℂ) = 1 := by simpa [hzero] using hsum
    exact zero_ne_one this
  obtain ⟨y, hy⟩ := hexists
  refine ⟨y, hy, ?_⟩
  intro z hz
  by_contra hyz
  have horth : F (recordBasis y) * F (recordBasis z) = 0 := by
    rw [← map_mul,
      recordBasis_mul_eq_zero (j := y) (k := z) (Ne.symm hyz), map_zero]
  have hpoint := congrFun horth x
  simp only [Pi.mul_apply, Pi.zero_apply] at hpoint
  rw [hy, hz] at hpoint
  simpa using hpoint

/-- The input label selected at an output label by a public star
automorphism. -/
noncomputable def publicLabelMap
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x : ι) : ι :=
  Classical.choose (existsUnique_recordBasis_image_eq_one F x)

@[simp]
theorem publicLabelMap_spec
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x : ι) :
    F (recordBasis (publicLabelMap F x)) x = 1 :=
  (Classical.choose_spec (existsUnique_recordBasis_image_eq_one F x)).1

theorem eq_publicLabelMap_of_recordBasis_image_eq_one
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x y : ι)
    (hy : F (recordBasis y) x = 1) :
    y = publicLabelMap F x :=
  (Classical.choose_spec (existsUnique_recordBasis_image_eq_one F x)).2 y hy

/-- The image of a record basis function is the indicator of the selected
label. -/
theorem starAutomorphism_recordBasis_apply
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (x y : ι) :
    F (recordBasis y) x = if y = publicLabelMap F x then 1 else 0 := by
  by_cases h : y = publicLabelMap F x
  · subst y
    simp
  · rw [if_neg h]
    rcases starAutomorphism_recordBasis_eq_zero_or_one F x y with hy | hy
    · exact hy
    · exact False.elim (h (eq_publicLabelMap_of_recordBasis_image_eq_one F x y hy))

/-- Every public star automorphism acts by pullback along its selected-label
map. -/
theorem starAutomorphism_apply_eq
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) (f : ι → ℂ) (x : ι) :
    F f x = f (publicLabelMap F x) := by
  classical
  calc
    F f x = F (∑ y : ι, f y • recordBasis y) x := by
      rw [sum_smul_recordBasis]
    _ = (∑ y : ι, F (f y • recordBasis y)) x := by rw [map_sum]
    _ = ∑ y : ι, f y * F (recordBasis y) x := by
      simp [Finset.sum_apply]
    _ = f (publicLabelMap F x) := by
      simp [starAutomorphism_recordBasis_apply]

/-- The selected-label map is injective.  Surjectivity of the algebra
automorphism lets arbitrary record basis functions separate output labels. -/
theorem publicLabelMap_injective
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) :
    Function.Injective (publicLabelMap F) := by
  classical
  intro x z hxz
  by_contra hxne
  obtain ⟨f, hf⟩ := F.surjective (recordBasis x)
  have hvalues : recordBasis x x = recordBasis x z := by
    calc
      recordBasis x x = F f x := (congrFun hf x).symm
      _ = f (publicLabelMap F x) := starAutomorphism_apply_eq F f x
      _ = f (publicLabelMap F z) := by rw [hxz]
      _ = F f z := (starAutomorphism_apply_eq F f z).symm
      _ = recordBasis x z := congrFun hf z
  simpa [recordBasis, hxne, Ne.symm hxne] using hvalues

/-- The permutation of public labels classified by a star automorphism. -/
noncomputable def publicLabelPerm
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) : Equiv.Perm ι :=
  Equiv.ofBijective (publicLabelMap F)
    ⟨publicLabelMap_injective F,
      Finite.surjective_of_injective (publicLabelMap_injective F)⟩

/-- Classification theorem: every star automorphism of a finite public
function algebra is pullback by one label permutation. -/
theorem publicStarAutomorphism_is_labelPermutation
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ)) :
    ∃ σ : Equiv.Perm ι, ∀ f x, F f x = f (σ x) := by
  refine ⟨publicLabelPerm F, ?_⟩
  intro f x
  exact starAutomorphism_apply_eq F f x

omit [Fintype ι] in
/-- The classifying label permutation is unique. -/
theorem publicStarAutomorphism_labelPermutation_unique
    (F : (ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ))
    {σ τ : Equiv.Perm ι}
    (hσ : ∀ f x, F f x = f (σ x))
    (hτ : ∀ f x, F f x = f (τ x)) :
    σ = τ := by
  ext x
  by_contra hne
  have h := (hσ (recordBasis (σ x)) x).symm.trans
    (hτ (recordBasis (σ x)) x)
  have : (1 : ℂ) = 0 := by
    simpa [recordBasis, hne, Ne.symm hne] using h
  exact one_ne_zero this

/-! ## Arbitrary continuous public automorphism flows -/

/-- A pointwise-continuous real-parameter group of public star
automorphisms.  Continuity is required only after applying the automorphism to
one function and evaluating at one label. -/
structure ContinuousPublicStarFlow (ι : Type u) [Fintype ι]
    [DecidableEq ι] where
  toAut : ℝ → ((ι → ℂ) ≃⋆ₐ[ℂ] (ι → ℂ))
  map_zero : toAut 0 = StarAlgEquiv.refl
  map_add : ∀ s t f, toAut (s + t) f = toAut s (toAut t f)
  continuous_apply : ∀ f x, Continuous fun t => toAut t f x

namespace ContinuousPublicStarFlow

/-- Pointwise operator continuity forces the permutation selected at a fixed
output label to be locally constant in the real parameter. -/
theorem publicLabelMap_isLocallyConstant
    (A : ContinuousPublicStarFlow ι) (x : ι) :
    IsLocallyConstant fun t : ℝ => publicLabelMap (A.toAut t) x := by
  classical
  refine IsLocallyConstant.iff_isOpen_fiber.mpr ?_
  intro y
  let g : ℝ → ℂ := fun t => A.toAut t (recordBasis y) x
  have hg : Continuous g := A.continuous_apply (recordBasis y) x
  have hopen : IsOpen (g ⁻¹' Metric.ball (1 : ℂ) (1 / 2 : ℝ)) :=
    hg.isOpen_preimage _ Metric.isOpen_ball
  have hfiber :
      (fun t : ℝ => publicLabelMap (A.toAut t) x) ⁻¹' ({y} : Set ι) =
        g ⁻¹' Metric.ball (1 : ℂ) (1 / 2 : ℝ) := by
    ext t
    simp only [Set.mem_preimage, Set.mem_singleton_iff, Metric.mem_ball]
    constructor
    · intro ht
      have hone : A.toAut t (recordBasis y) x = 1 := by
        rw [← ht]
        exact publicLabelMap_spec (A.toAut t) x
      have hone' : g t = 1 := by simpa [g] using hone
      change dist (g t) 1 < 1 / 2
      rw [hone']
      norm_num
    · intro ht
      change dist (A.toAut t (recordBasis y) x) 1 < 1 / 2 at ht
      rcases starAutomorphism_recordBasis_eq_zero_or_one
          (A.toAut t) x y with hzero | hone
      · rw [hzero] at ht
        norm_num [Complex.dist_eq] at ht
      · exact (eq_publicLabelMap_of_recordBasis_image_eq_one
          (A.toAut t) x y hone).symm
  rw [hfiber]
  exact hopen

/-- The label selected by a pointwise-continuous public automorphism group is
constant along the real parameter and equals the original output label. -/
theorem publicLabelMap_eq
    (A : ContinuousPublicStarFlow ι) (t : ℝ) (x : ι) :
    publicLabelMap (A.toAut t) x = x := by
  have hconst := (A.publicLabelMap_isLocallyConstant x)
    |>.apply_eq_of_preconnectedSpace t 0
  have hzero : publicLabelMap (A.toAut 0) x = x := by
    have hxone : A.toAut 0 (recordBasis x) x = 1 := by
      rw [A.map_zero]
      simp [recordBasis]
    exact (eq_publicLabelMap_of_recordBasis_image_eq_one
      (A.toAut 0) x x hxone).symm
  exact hconst.trans hzero

/-- Every pointwise-continuous real-parameter group of finite public
star-algebra automorphisms is the identity. -/
theorem toAut_eq_refl
    (A : ContinuousPublicStarFlow ι) (t : ℝ) :
    A.toAut t = StarAlgEquiv.refl := by
  ext f x
  rw [starAutomorphism_apply_eq, A.publicLabelMap_eq]
  rfl

end ContinuousPublicStarFlow

-- Axiom audit: declarations may use only the standard Mathlib basis.
#print axioms OPH.Dynamics.publicStarAutomorphism_is_labelPermutation
#print axioms OPH.Dynamics.publicStarAutomorphism_labelPermutation_unique
#print axioms OPH.Dynamics.ContinuousPublicStarFlow.toAut_eq_refl

end OPH.Dynamics
