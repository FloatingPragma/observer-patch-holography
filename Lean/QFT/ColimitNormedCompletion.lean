import QFT.LocallyCovariantLimit
import Mathlib.Analysis.CStarAlgebra.Matrix
import Mathlib.Analysis.CStarAlgebra.Spectrum
import Mathlib.Analysis.CStarAlgebra.Classes
import Mathlib.Analysis.Normed.Algebra.Spectrum
import Mathlib.Analysis.Normed.Module.Completion
import Mathlib.Topology.Algebra.UniformRing
import Mathlib.Topology.Algebra.StarSubalgebra
import Mathlib.LinearAlgebra.Eigenspace.Minpoly

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The inductive-limit C*-seminorm and the normed star completion of the tower colimit

This module sits on top of the committed filtered-colimit layer
`QFT.LocallyCovariantLimit` and adds the analytic layer that the committed
module's claim boundary excludes: a limit seminorm on `TowerColimit T`, its
kernel, the normed star structure it induces, and the completion.  The
module imports the committed layer and Mathlib and nothing else of the
project.

## What is proved

* The stage norms are the L2 operator norms of the finite matrix stages.
  Every connecting map of the tower is a unital star algebra homomorphism
  between finite matrix C*-algebras, hence norm nonincreasing
  (`algebraRefine_norm_le`).
* `colimitNorm` is the inductive-limit seminorm on the committed colimit:
  the value at a class is the infimum over all declared refinements of the
  stage norm of the refined representative (`stageValues`, `germNorm`).
  Representative independence is `germNorm_congr`, proved from
  contractivity of the connecting maps.
* The value set of one representative is finite
  (`stageValues_finite`): the norm of a star homomorphic image of `X` is the
  square root of a spectral value of `star X * X`, spectra shrink under
  algebra homomorphisms, and the spectrum of a finite matrix is finite.
  The defining infimum is therefore attained at some declared refinement
  (`exists_stage_norm_eq`), and every class has a stable representative at
  which no further refinement changes the norm (`exists_stableRep`).
* Seminorm laws: `colimitNorm_zero`, `colimitNorm_add_le`,
  `colimitNorm_smul`, `colimitNorm_neg`, submultiplicativity
  `colimitNorm_mul_le`, star invariance `colimitNorm_star`, and the
  C*-identity `colimitNorm_star_mul_self`.
* The exact kernel theorem `colimitNorm_mk_eq_zero_iff`: a stage element has
  limit seminorm zero exactly when some declared refinement sends it to
  zero.  Composed with the committed receipt `colimitMk_eq_zero_iff` this
  reads `colimitNorm_mk_eq_zero_iff_mk_eq_zero` and
  `colimitNorm_eq_zero_iff`: the seminorm kernel is the committed
  refinement-null kernel, and that kernel is the zero subgroup, so the
  seminorm is a norm on the committed colimit.  The norm topology of the
  colimit is T2 (`colimit_t2Space`).
* Normed structure on the colimit: `NormedRing`, `NormedAlgebra ℂ`,
  `NormedStarGroup`, `CStarRing` instances, with the norm definitionally
  equal to `colimitNorm` (`norm_def`).
* The completion `UniformSpace.Completion (TowerColimit T)` carries a star
  ring structure extended by uniform continuity from the isometric star of
  the colimit (`completionStarRing`, `completion_star_coe`), a continuous
  star, the star module law over `ℂ`, the normed algebra law, the isometry
  of the star (`completion_norm_star`), and the C*-identity transported
  from the dense image (`completion_norm_star_mul_self`).  The bundle
  `CStarAlgebra (Completion (TowerColimit T))` is registered
  (`completionCStarAlgebra`).
* The canonical map `completionMkHom` is a star algebra homomorphism that
  is isometric, injective, and has dense range; `stageToCompletion` is the
  composite cocone map from each stage, with
  `stageToCompletion_refine` and `norm_stageToCompletion_le`.
* Transport: `completionLocalAlgebra` is the closure in the completion of
  the image of a committed limit local algebra.  Isotony
  (`completionLocalAlgebra_mono`), refinement independence
  (`completionLocalAlgebra_refine`), locality on disjoint regions
  (`completionLocalAlgebra_commute`), and the committed order-theoretic
  Cauchy-embedding equality
  (`completionLocalAlgebra_eq_of_isCauchyEmbedding`) hold for the
  closures.  The committed statement is an equality of subalgebras, and
  the closure of the image is a function of that equality.
* The hypotheses of the four transported statements are inhabited on the
  committed support-graded net.
  `supportGradedNet_completionLocalAlgebra_top_eq_full` instantiates the
  Cauchy-embedding equality at the committed Cauchy witness `topToFull`,
  `supportGradedNet_completionLocalAlgebra_bottom_le_top` instantiates
  isotony at `bottomToTop`,
  `supportGradedNet_completionLocalAlgebra_commute` instantiates locality
  at the committed disjoint pair, and
  `supportGradedNet_completionMkHom_e0_mem` together with
  `supportGradedNet_norm_completionMkHom_e0` places an element of norm one
  in the transported local algebra of the top region.
* Nonvacuity: on the committed witness tower the limit seminorm of a class
  equals the stage norm of its representative
  (`witnessTower_colimitNorm_mk`), the separating projector `e0` has limit
  seminorm one (`witnessTower_colimitNorm_e0`), its image in the
  completion has norm one (`witnessTower_norm_stageToCompletion_e0`), and
  the unit has norm one (`witnessTower_norm_one`).
* Controls: on the committed collapsing tower the coarse unit is a nonzero
  stage element of stage norm one whose limit seminorm is zero
  (`collapsingTower_colimitNorm_one`), so the kernel clause has content.
  The star clause of the connecting maps is load-bearing for the
  contractivity step that makes the seminorm well defined: a unital
  algebra homomorphism of the two-by-two stage that conjugates by a
  diagonal scaling fails contractivity (`scalingConj_not_contractive`),
  and therefore fails star preservation
  (`scalingConj_not_star_preserving`).

## Declared inputs

One input is declared: the choice of the L2 operator norm on each matrix
stage, read from Mathlib's scoped `Matrix.Norms.L2Operator` instances, and
the `CStarAlgebra` bundle assembled from those instances
(`matrixStageCStarAlgebra`).  Mathlib scopes its four matrix norms so that
no global norm is fixed for `Matrix`, and this module keeps that
discipline: `matrixStageCStarAlgebra` is a `def` whose instance attribute
sits in the scope `OPH.QFT.StageNorms`, which the module opens for its own
elaboration.  An importer that wants the stage bundle writes
`open scoped OPH.QFT.StageNorms`.

The five completion star declarations `completionStarRing`,
`completionContinuousStar`, `completionNormedStarGroup`,
`completionCStarRing`, and `completionStarModule` range over every `A`
carrying `[NormedRing A] [StarRing A] [NormedStarGroup A]`, so their
instance attributes sit in the scope `OPH.QFT.CompletionStar` for the same
reason.  The instances this module exports globally are attached to its
own objects: `instTowerColimitNorm`, `instTowerColimitNormedRing`,
`instTowerColimitNormedSpace`, `instTowerColimitNormedAlgebra`,
`instTowerColimitNormedStarGroup`, and `instTowerColimitCStarRing` on
`TowerColimit T`, together with `completionNormedAlgebra` and
`completionCStarAlgebra` on `Completion (TowerColimit T)`.

Mathlib registers no star structure on `UniformSpace.Completion` for any
base, so the star of the completion is constructed here.  Mathlib
registers the normed algebra law on completions for a commutative base
only, so that law is constructed here as well.  The stage algebras, the
connecting star algebra homomorphisms, the colimit, the germ relation, the
limit local algebras, the Cauchy embedding class, and the committed
witnesses are read from the committed layer without change.  No
C*-structure is assumed as a typeclass input: the stages are finite
complex matrix algebras and the connecting maps are star algebra
homomorphisms in the committed interface.

## Premise consumption per register row

The module consumes no physical premise row of
`tracking/premise_register.json`.  It uses the algebraic and
order-theoretic content of the committed colimit together with Mathlib's
C*-algebra library.  Two rows bound what the module does not reach.
PR-58 records that no registered premise supplies a time-indexed or
Cauchy-region interface with a net-compatible evolution; the statement
transported here is the committed order-theoretic equality over
`IsCauchyEmbedding` and nothing more.  PR-52 owns the physical spacetime
and causal reading, which this module does not use.  The measurement-facing
rows PR-02, PR-43, and PR-44 are untouched.

## What is NOT proved

No Lorentzian, spacetime, causal, or continuum reading is attached to any
object here.  The completion is a C*-algebra obtained from a filtered
colimit of finite complex matrix algebras; no field, field operator,
vacuum state, Wightman or Haag-Kastler axiom, or physical Hilbert space
representation is constructed.  The local algebras of the completion are
closures of images of the committed limit local algebras; no claim is made
that they exhaust the completion or satisfy Haag duality, additivity, or
nuclearity.  No quotient object is constructed anywhere in the module:
the seminorm kernel is the zero subgroup, so the colimit carries the norm
directly.  No CP or CPTP structure is inferred.  The C*-identity is
proved for the constructed norm and its extension to the completion; no
uniqueness of that norm among C*-norms on the colimit is claimed.  The
seminorm laws are proved through the attained infimum, which uses the
finite dimension of the committed stages; the module states nothing about
stage systems with infinite-dimensional stages, and the committed stage
type contains none.  The non-star control is a homomorphism of one stage
algebra and is not a member of any `ConsensusTower`, whose interface types
every connecting map as a star algebra homomorphism; the control shows
which clause the contractivity lemma consumes.
-/

namespace OPH.QFT

open OPH.Tower
open CategoryTheory
open UniformSpace
open scoped Matrix.Norms.L2Operator

universe u

/- Instance search on the colimit and its completion passes through the
full Mathlib C*-algebra hierarchy; the committed layer sets the same budget
for the same reason. -/
set_option synthInstance.maxHeartbeats 400000

/-! ## Declared input: the C*-structure of the matrix stages -/

/-- **Declared input.**  The C*-algebra bundle on each finite matrix stage,
assembled from Mathlib's scoped L2 operator norm instances.  No field is
supplied by hand: every parent structure is the scoped instance. -/
@[reducible] noncomputable def matrixStageCStarAlgebra (n : ℕ) :
    CStarAlgebra (Matrix (Fin n) (Fin n) ℂ) where

/- Mathlib scopes its matrix norms so that no global norm is fixed for
`Matrix`.  The stage bundle follows that discipline: the instance
attribute sits in the scope `OPH.QFT.StageNorms`, which the line below
opens for the rest of this module. -/
scoped[OPH.QFT.StageNorms] attribute [instance] matrixStageCStarAlgebra

open scoped OPH.QFT.StageNorms

section StageCStar

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- Every connecting map of the tower is norm nonincreasing: it is a star
algebra homomorphism between C*-algebras. -/
theorem algebraRefine_norm_le {r s : ι} (hrs : r ≤ s)
    (X : ConsensusTower.PrivateAlgebra T r) :
    ‖T.algebraRefine hrs X‖ ≤ ‖X‖ :=
  NonUnitalStarAlgHom.norm_apply_le _ _

end StageCStar

/-! ## The inductive-limit seminorm -/

section LimitSeminorm

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- The stage value set of a representative: the stage norms of all its
declared refinements. -/
def stageValues {r : ι} (X : ConsensusTower.PrivateAlgebra T r) : Set ℝ :=
  {a | ∃ (t : ι) (h : r ≤ t), a = ‖T.algebraRefine h X‖}

/-- The stage norm at any declared refinement is a stage value. -/
theorem norm_refine_mem_stageValues {r t : ι} (h : r ≤ t)
    (X : ConsensusTower.PrivateAlgebra T r) :
    ‖T.algebraRefine h X‖ ∈ stageValues (T := T) X :=
  ⟨t, h, rfl⟩

/-- The stage norm of the representative itself is a stage value. -/
theorem norm_mem_stageValues {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    ‖X‖ ∈ stageValues (T := T) X :=
  ⟨r, le_rfl, by rw [T.algebra_refine_refl]⟩

/-- The stage value set is nonempty. -/
theorem stageValues_nonempty {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    (stageValues (T := T) X).Nonempty :=
  ⟨_, norm_mem_stageValues X⟩

/-- Every stage value is nonnegative. -/
theorem stageValues_nonneg {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    ∀ a ∈ stageValues (T := T) X, 0 ≤ a := by
  rintro a ⟨t, h, rfl⟩
  exact norm_nonneg _

/-- The stage value set is bounded below. -/
theorem stageValues_bddBelow {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    BddBelow (stageValues (T := T) X) :=
  ⟨0, stageValues_nonneg X⟩

/-- Refining the representative shrinks its value set. -/
theorem stageValues_refine_subset {r s : ι} (hrs : r ≤ s)
    (X : ConsensusTower.PrivateAlgebra T r) :
    stageValues (T := T) (T.algebraRefine hrs X) ⊆ stageValues X := by
  rintro a ⟨t, hst, rfl⟩
  exact ⟨t, hrs.trans hst, by rw [T.algebra_refine_trans]⟩

/-- The infimum of the value set is refinement independent.  One inequality
is the inclusion above; the other uses a common refinement and
contractivity of the connecting maps. -/
theorem sInf_stageValues_refine {r s : ι} (hrs : r ≤ s)
    (X : ConsensusTower.PrivateAlgebra T r) :
    sInf (stageValues (T := T) (T.algebraRefine hrs X)) = sInf (stageValues X) := by
  apply le_antisymm
  · refine le_csInf (stageValues_nonempty X) ?_
    rintro a ⟨u, hru, rfl⟩
    obtain ⟨v, huv, hsv⟩ := T.commonRefinement u s
    calc sInf (stageValues (T.algebraRefine hrs X))
        ≤ ‖T.algebraRefine hsv (T.algebraRefine hrs X)‖ :=
          csInf_le (stageValues_bddBelow _) (norm_refine_mem_stageValues hsv _)
      _ = ‖T.algebraRefine huv (T.algebraRefine hru X)‖ := by
          rw [T.algebra_refine_trans, T.algebra_refine_trans]
      _ ≤ ‖T.algebraRefine hru X‖ := algebraRefine_norm_le huv _
  · exact csInf_le_csInf (stageValues_bddBelow X) (stageValues_nonempty _)
      (stageValues_refine_subset hrs X)

/-- The germ-level seminorm: the infimum of the stage value set. -/
noncomputable def germNorm (p : TowerGerm T) : ℝ :=
  sInf (stageValues p.2)

/-- Representative independence of the germ seminorm. -/
theorem germNorm_congr {p q : TowerGerm T} (h : germRel T p q) :
    germNorm p = germNorm q := by
  obtain ⟨t, hp, hq, e⟩ := h
  show sInf (stageValues p.2) = sInf (stageValues q.2)
  rw [← sInf_stageValues_refine hp p.2, ← sInf_stageValues_refine hq q.2, e]

/-- **The inductive-limit seminorm** on the committed filtered-colimit
algebra. -/
noncomputable def colimitNorm (x : TowerColimit T) : ℝ :=
  Quotient.liftOn x germNorm fun _ _ h => germNorm_congr h

/-- The limit seminorm of a stage class is the infimum of the stage value
set of its representative. -/
theorem colimitNorm_mk {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    colimitNorm (colimitMk (T := T) X) = sInf (stageValues X) :=
  rfl

/-- The limit seminorm of a stage class is at most the stage norm at any
declared refinement of its representative. -/
theorem colimitNorm_mk_le_refine {r t : ι} (h : r ≤ t)
    (X : ConsensusTower.PrivateAlgebra T r) :
    colimitNorm (colimitMk (T := T) X) ≤ ‖T.algebraRefine h X‖ :=
  csInf_le (stageValues_bddBelow X) (norm_refine_mem_stageValues h X)

/-- The limit seminorm of a stage class is at most the stage norm of its
representative. -/
theorem colimitNorm_mk_le {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    colimitNorm (colimitMk (T := T) X) ≤ ‖X‖ :=
  csInf_le (stageValues_bddBelow X) (norm_mem_stageValues X)

/-- The limit seminorm is nonnegative. -/
theorem colimitNorm_nonneg (x : TowerColimit T) : 0 ≤ colimitNorm x := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  exact Real.sInf_nonneg (stageValues_nonneg X)

/-- A representative whose stage norm is constant along every declared
refinement realizes that constant as the limit seminorm of its class. -/
theorem colimitNorm_mk_eq_of_forall {r : ι} (X : ConsensusTower.PrivateAlgebra T r)
    {a : ℝ} (h : ∀ (t : ι) (ht : r ≤ t), ‖T.algebraRefine ht X‖ = a) :
    colimitNorm (colimitMk (T := T) X) = a := by
  have hX : ‖X‖ = a := by
    have := h r le_rfl
    rwa [T.algebra_refine_refl] at this
  apply le_antisymm
  · rw [← hX]
    exact colimitNorm_mk_le X
  · refine le_csInf (stageValues_nonempty X) ?_
    rintro b ⟨t, ht, rfl⟩
    rw [h t ht]

end LimitSeminorm

/-! ## Finiteness of the value set and attainment of the infimum -/

section Attainment

/-- The norm of the image of `X` under a star algebra homomorphism of
finite matrix algebras is zero or the square root of a spectral value of
`star X * X`: the C*-identity turns the squared norm into the norm of a
selfadjoint element, which equals its spectral radius, which is attained
on the spectrum, and spectra shrink under algebra homomorphisms. -/
theorem norm_starAlgHom_mem_spectralSet {n m : ℕ}
    (φ : Matrix (Fin n) (Fin n) ℂ →⋆ₐ[ℂ] Matrix (Fin m) (Fin m) ℂ)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    ‖φ X‖ ∈ insert (0 : ℝ)
      ((fun z : ℂ => Real.sqrt ‖z‖) '' (spectrum ℂ (star X * X))) := by
  rcases Nat.eq_zero_or_pos m with hm | hm
  · subst hm
    left
    have h0 : φ X = 0 := Subsingleton.elim (α := Matrix (Fin 0) (Fin 0) ℂ) _ _
    rw [h0, norm_zero]
  · haveI : Nonempty (Fin m) := ⟨⟨0, hm⟩⟩
    haveI : Nontrivial (Matrix (Fin m) (Fin m) ℂ) := inferInstance
    right
    set a : Matrix (Fin m) (Fin m) ℂ := star (φ X) * φ X with ha
    have hsa : IsSelfAdjoint a := IsSelfAdjoint.star_mul_self _
    obtain ⟨z, hz, hzn⟩ :=
      spectrum.exists_nnnorm_eq_spectralRadius_of_nonempty (spectrum.nonempty a)
    have hrad : spectralRadius ℂ a = ‖a‖₊ := hsa.spectralRadius_eq_nnnorm
    have hnn : (‖z‖₊ : ENNReal) = (‖a‖₊ : ENNReal) := by rw [hzn, hrad]
    have hz' : ‖z‖ = ‖a‖ := by
      have h := congrArg ENNReal.toNNReal hnn
      simpa using congrArg NNReal.toReal h
    have hsub : spectrum ℂ a ⊆ spectrum ℂ (star X * X) := by
      have hcalc : a = φ (star X * X) := by rw [ha, map_mul, map_star]
      rw [hcalc]
      exact AlgHom.spectrum_apply_subset φ _
    refine ⟨z, hsub hz, ?_⟩
    have hna : ‖a‖ = ‖φ X‖ * ‖φ X‖ := CStarRing.norm_star_mul_self
    show Real.sqrt ‖z‖ = ‖φ X‖
    rw [hz', hna, Real.sqrt_mul_self (norm_nonneg _)]

/-- The spectral value set of a finite matrix is finite. -/
theorem spectralSet_finite {n : ℕ} (X : Matrix (Fin n) (Fin n) ℂ) :
    (insert (0 : ℝ)
      ((fun z : ℂ => Real.sqrt ‖z‖) '' (spectrum ℂ (star X * X)))).Finite :=
  ((Matrix.finite_spectrum (star X * X)).image _).insert 0

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- **Finiteness of the value set.**  Every stage value of `X` lies in the
finite spectral set of `star X * X`. -/
theorem stageValues_finite {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    (stageValues (T := T) X).Finite := by
  refine (spectralSet_finite X).subset ?_
  rintro a ⟨t, ht, rfl⟩
  exact norm_starAlgHom_mem_spectralSet (T.algebraRefine ht) X

/-- **Attainment.**  The defining infimum of the limit seminorm is the stage
norm at some declared refinement. -/
theorem exists_stage_norm_eq {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    ∃ (t : ι) (h : r ≤ t),
      ‖T.algebraRefine h X‖ = colimitNorm (colimitMk (T := T) X) := by
  obtain ⟨t, h, e⟩ := (stageValues_nonempty X).csInf_mem (stageValues_finite X)
  exact ⟨t, h, e.symm⟩

/-- A stable representative: no declared refinement changes its stage
norm. -/
def IsStableRep {r : ι} (X : ConsensusTower.PrivateAlgebra T r) : Prop :=
  ∀ (t : ι) (h : r ≤ t), ‖T.algebraRefine h X‖ = ‖X‖

/-- A stable representative realizes its stage norm as the limit seminorm of
its class. -/
theorem IsStableRep.colimitNorm_mk {r : ι} {X : ConsensusTower.PrivateAlgebra T r}
    (hX : IsStableRep (T := T) X) :
    colimitNorm (colimitMk (T := T) X) = ‖X‖ :=
  colimitNorm_mk_eq_of_forall X hX

/-- Stability of a representative passes to its declared refinements. -/
theorem IsStableRep.refine {r s : ι} {X : ConsensusTower.PrivateAlgebra T r}
    (hX : IsStableRep (T := T) X) (hrs : r ≤ s) :
    IsStableRep (T := T) (T.algebraRefine hrs X) := by
  intro t hst
  rw [T.algebra_refine_trans, hX, hX]

/-- Every class has a stable representative above any given regulator. -/
theorem exists_stableRep (x : TowerColimit T) (r₀ : ι) :
    ∃ (t : ι) (_ : r₀ ≤ t) (X : ConsensusTower.PrivateAlgebra T t),
      x = colimitMk X ∧ IsStableRep X := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  obtain ⟨s, hrs, hs⟩ := exists_stage_norm_eq X
  obtain ⟨t, hst, hr₀t⟩ := T.commonRefinement s r₀
  have hstable : IsStableRep (T := T) (T.algebraRefine hrs X) := by
    intro u hsu
    apply le_antisymm (algebraRefine_norm_le hsu _)
    rw [hs, T.algebra_refine_trans]
    exact colimitNorm_mk_le_refine (hrs.trans hsu) X
  refine ⟨t, hr₀t, T.algebraRefine hst (T.algebraRefine hrs X), ?_, hstable.refine hst⟩
  rw [colimitMk_refine, colimitMk_refine]

/-- Two classes have stable representatives at one common regulator. -/
theorem exists_stableRep₂ (x y : TowerColimit T) :
    ∃ (t : ι) (X Y : ConsensusTower.PrivateAlgebra T t),
      x = colimitMk X ∧ y = colimitMk Y ∧ IsStableRep X ∧ IsStableRep Y := by
  obtain ⟨r, -, X, rfl, hX⟩ := exists_stableRep x (towerBase T)
  obtain ⟨t, hrt, Y, rfl, hY⟩ := exists_stableRep y r
  exact ⟨t, T.algebraRefine hrt X, Y, (colimitMk_refine hrt X).symm, rfl,
    hX.refine hrt, hY⟩

end Attainment

/-! ## Seminorm laws and the C*-identity -/

section SeminormLaws

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- The limit seminorm of zero is zero. -/
theorem colimitNorm_zero : colimitNorm (0 : TowerColimit T) = 0 := by
  rw [← colimitMk_zero (towerBase T)]
  apply le_antisymm _ (colimitNorm_nonneg _)
  calc colimitNorm (colimitMk (0 : ConsensusTower.PrivateAlgebra T (towerBase T)))
      ≤ ‖(0 : ConsensusTower.PrivateAlgebra T (towerBase T))‖ := colimitNorm_mk_le _
    _ = 0 := norm_zero

/-- The triangle inequality for the limit seminorm. -/
theorem colimitNorm_add_le (x y : TowerColimit T) :
    colimitNorm (x + y) ≤ colimitNorm x + colimitNorm y := by
  obtain ⟨t, X, Y, rfl, rfl, hX, hY⟩ := exists_stableRep₂ x y
  rw [colimitMk_add, hX.colimitNorm_mk, hY.colimitNorm_mk]
  exact (colimitNorm_mk_le _).trans (norm_add_le X Y)

/-- Submultiplicativity. -/
theorem colimitNorm_mul_le (x y : TowerColimit T) :
    colimitNorm (x * y) ≤ colimitNorm x * colimitNorm y := by
  obtain ⟨t, X, Y, rfl, rfl, hX, hY⟩ := exists_stableRep₂ x y
  rw [colimitMk_mul, hX.colimitNorm_mk, hY.colimitNorm_mk]
  exact (colimitNorm_mk_le _).trans (norm_mul_le X Y)

/-- Absolute homogeneity of the limit seminorm over `ℂ`. -/
theorem colimitNorm_smul (c : ℂ) (x : TowerColimit T) :
    colimitNorm (c • x) = ‖c‖ * colimitNorm x := by
  obtain ⟨t, -, X, rfl, hX⟩ := exists_stableRep x (towerBase T)
  rw [colimitMk_smul, hX.colimitNorm_mk]
  apply colimitNorm_mk_eq_of_forall
  intro u hu
  rw [map_smul, norm_smul, hX u hu]

/-- The limit seminorm is invariant under negation. -/
theorem colimitNorm_neg (x : TowerColimit T) : colimitNorm (-x) = colimitNorm x := by
  obtain ⟨t, -, X, rfl, hX⟩ := exists_stableRep x (towerBase T)
  rw [colimitMk_neg, hX.colimitNorm_mk]
  apply colimitNorm_mk_eq_of_forall
  intro u hu
  rw [map_neg, norm_neg, hX u hu]

/-- Star invariance. -/
theorem colimitNorm_star (x : TowerColimit T) : colimitNorm (star x) = colimitNorm x := by
  obtain ⟨t, -, X, rfl, hX⟩ := exists_stableRep x (towerBase T)
  rw [colimitMk_star, hX.colimitNorm_mk]
  apply colimitNorm_mk_eq_of_forall
  intro u hu
  rw [map_star, norm_star, hX u hu]

/-- **The C*-identity of the limit seminorm.** -/
theorem colimitNorm_star_mul_self (x : TowerColimit T) :
    colimitNorm (star x * x) = colimitNorm x * colimitNorm x := by
  obtain ⟨t, -, X, rfl, hX⟩ := exists_stableRep x (towerBase T)
  rw [colimitMk_star, colimitMk_mul, hX.colimitNorm_mk]
  apply colimitNorm_mk_eq_of_forall
  intro u hu
  rw [map_mul, map_star, CStarRing.norm_star_mul_self, hX u hu]

end SeminormLaws

/-! ## The exact kernel theorem -/

section Kernel

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- **The exact kernel theorem.**  A stage element has limit seminorm zero
exactly when some declared refinement sends it to zero. -/
theorem colimitNorm_mk_eq_zero_iff {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    colimitNorm (colimitMk (T := T) X) = 0 ↔
      ∃ (t : ι) (h : r ≤ t), T.algebraRefine h X = 0 := by
  constructor
  · intro h0
    obtain ⟨t, h, e⟩ := exists_stage_norm_eq X
    exact ⟨t, h, norm_eq_zero.mp (e.trans h0)⟩
  · rintro ⟨t, h, e⟩
    apply le_antisymm _ (colimitNorm_nonneg _)
    calc colimitNorm (colimitMk (T := T) X)
        ≤ ‖T.algebraRefine h X‖ := colimitNorm_mk_le_refine h X
      _ = 0 := by rw [e, norm_zero]

/-- The seminorm kernel is the committed refinement-null kernel. -/
theorem colimitNorm_mk_eq_zero_iff_mk_eq_zero {r : ι}
    (X : ConsensusTower.PrivateAlgebra T r) :
    colimitNorm (colimitMk (T := T) X) = 0 ↔ colimitMk X = (0 : TowerColimit T) := by
  rw [colimitNorm_mk_eq_zero_iff, colimitMk_eq_zero_iff]

/-- **Definiteness.**  The limit seminorm is a norm on the committed
colimit: the seminorm kernel is the zero subgroup. -/
theorem colimitNorm_eq_zero_iff (x : TowerColimit T) : colimitNorm x = 0 ↔ x = 0 := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  exact colimitNorm_mk_eq_zero_iff_mk_eq_zero X

end Kernel

/-! ## Normed star structure on the colimit -/

section NormedColimit

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- The norm of the colimit is the limit seminorm. -/
noncomputable instance instTowerColimitNorm : Norm (TowerColimit T) :=
  ⟨colimitNorm⟩

/-- The norm of the colimit unfolds to `colimitNorm`. -/
theorem norm_def (x : TowerColimit T) : ‖x‖ = colimitNorm x :=
  rfl

/-- The four core laws that make the limit seminorm a norm on the colimit:
nonnegativity, absolute homogeneity, the triangle inequality, and
definiteness.  `NormedSpace.Core` is a `Prop`, so this is a proof and
carries no data. -/
theorem colimitNormCore (T : ConsensusTower ι) : NormedSpace.Core ℂ (TowerColimit T) where
  norm_nonneg := colimitNorm_nonneg
  norm_smul := colimitNorm_smul
  norm_triangle := colimitNorm_add_le
  norm_eq_zero_iff := colimitNorm_eq_zero_iff

/-- The normed ring structure of the colimit.  The ring part is the
committed instance; the metric is the one induced by the limit norm. -/
noncomputable instance instTowerColimitNormedRing : NormedRing (TowerColimit T) :=
  { NormedAddCommGroup.ofCore (colimitNormCore T),
      (inferInstance : Ring (TowerColimit T)) with
    norm_mul_le := colimitNorm_mul_le }

/-- The colimit is a normed space over `ℂ` for the limit norm; the scalar
law is `colimitNorm_smul`. -/
noncomputable instance instTowerColimitNormedSpace : NormedSpace ℂ (TowerColimit T) where
  norm_smul_le c x := le_of_eq (colimitNorm_smul c x)

/-- The colimit is a normed `ℂ`-algebra for the limit norm. -/
noncomputable instance instTowerColimitNormedAlgebra : NormedAlgebra ℂ (TowerColimit T) where
  norm_smul_le := norm_smul_le

/-- The star of the colimit is isometric for the limit norm; the law is
`colimitNorm_star`. -/
instance instTowerColimitNormedStarGroup : NormedStarGroup (TowerColimit T) :=
  ⟨fun x => le_of_eq (colimitNorm_star x)⟩

/-- The colimit with the limit norm satisfies the C*-identity. -/
instance instTowerColimitCStarRing : CStarRing (TowerColimit T) :=
  ⟨fun x => le_of_eq (colimitNorm_star_mul_self x).symm⟩

/-- The norm topology of the colimit is T2.  The content sits in
`colimitNorm_eq_zero_iff`; this statement reads the separation off the
metric induced by the limit norm. -/
theorem colimit_t2Space : T2Space (TowerColimit T) :=
  inferInstance

/-- The colimit inclusion of a stage element is norm nonincreasing. -/
theorem norm_colimitMk_le {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    ‖colimitMk (T := T) X‖ ≤ ‖X‖ :=
  colimitNorm_mk_le X

end NormedColimit

/-! ## Star structure on the completion of a seminormed star ring

Mathlib registers no star structure on `UniformSpace.Completion` for any
base.  The star of a seminormed star ring with isometric star extends by
uniform continuity, and the ring laws transport from the dense image.
The five declarations of this section range over every normed star ring,
so each instance attribute sits in the scope `OPH.QFT.CompletionStar` and
the global instance graph is unchanged. -/

section CompletionStar

variable {A : Type*} [NormedRing A] [StarRing A] [NormedStarGroup A]

/-- The star of a normed star group is an isometry. -/
theorem star_isometry_of_normedStarGroup : Isometry (star : A → A) :=
  AddMonoidHomClass.isometry_of_norm (starAddEquiv : A ≃+ A) fun a => norm_star a

/-- The uniformly continuous extension of the star acts on the dense image by
the star of the base. -/
theorem completion_map_star_coe (a : A) :
    Completion.map (star : A → A) (a : Completion A) = ((star a : A) : Completion A) :=
  Completion.map_coe (star_isometry_of_normedStarGroup (A := A)).uniformContinuous a

omit [NormedStarGroup A] in
/-- The uniformly continuous extension of the star is continuous. -/
theorem completion_map_star_continuous :
    Continuous (Completion.map (star : A → A) : Completion A → Completion A) :=
  Completion.continuous_map

/-- The star ring structure on the completion, extended from the dense
image. -/
@[reducible] noncomputable def completionStarRing : StarRing (Completion A) where
  star := Completion.map (star : A → A)
  star_involutive x := by
    induction x using Completion.induction_on with
    | hp =>
      exact isClosed_eq
        (completion_map_star_continuous.comp completion_map_star_continuous) continuous_id
    | ih a =>
      rw [completion_map_star_coe (A := A), completion_map_star_coe (A := A), star_star]
  star_mul x y := by
    induction x, y using Completion.induction_on₂ with
    | hp =>
      exact isClosed_eq (completion_map_star_continuous.comp continuous_mul)
        (continuous_mul.comp ((completion_map_star_continuous.comp continuous_snd).prodMk
          (completion_map_star_continuous.comp continuous_fst)))
    | ih a b =>
      rw [← Completion.coe_mul, completion_map_star_coe (A := A),
        completion_map_star_coe (A := A), completion_map_star_coe (A := A),
        ← Completion.coe_mul, star_mul]
  star_add x y := by
    induction x, y using Completion.induction_on₂ with
    | hp =>
      exact isClosed_eq (completion_map_star_continuous.comp continuous_add)
        (continuous_add.comp ((completion_map_star_continuous.comp continuous_fst).prodMk
          (completion_map_star_continuous.comp continuous_snd)))
    | ih a b =>
      rw [← Completion.coe_add, completion_map_star_coe (A := A),
        completion_map_star_coe (A := A), completion_map_star_coe (A := A),
        ← Completion.coe_add, star_add]

scoped[OPH.QFT.CompletionStar] attribute [instance] completionStarRing

open scoped OPH.QFT.CompletionStar

/-- The star of the completion acts on the dense image by the star of the
base. -/
theorem completion_star_coe (a : A) :
    star (a : Completion A) = ((star a : A) : Completion A) :=
  completion_map_star_coe a

/-- The star of the completion is continuous. -/
@[reducible] def completionContinuousStar : ContinuousStar (Completion A) :=
  ⟨completion_map_star_continuous⟩

/-- The star of the completion is isometric. -/
theorem completion_norm_star (x : Completion A) : ‖star x‖ = ‖x‖ := by
  induction x using Completion.induction_on with
  | hp =>
    exact isClosed_eq (continuous_norm.comp completion_map_star_continuous) continuous_norm
  | ih a => rw [completion_star_coe (A := A), Completion.norm_coe, Completion.norm_coe, norm_star]

/-- The star of the completion is isometric for the extended norm. -/
@[reducible] def completionNormedStarGroup : NormedStarGroup (Completion A) :=
  ⟨fun x => le_of_eq (completion_norm_star x)⟩

/-- **The C*-identity transports to the completion** from the dense image,
by continuity of the norm, the star, and the multiplication. -/
theorem completion_norm_star_mul_self [CStarRing A] (x : Completion A) :
    ‖star x * x‖ = ‖x‖ * ‖x‖ := by
  induction x using Completion.induction_on with
  | hp =>
    exact isClosed_eq
      (continuous_norm.comp (continuous_mul.comp
        (completion_map_star_continuous.prodMk continuous_id)))
      (continuous_norm.mul continuous_norm)
  | ih a =>
    rw [completion_star_coe (A := A), ← Completion.coe_mul, Completion.norm_coe,
      Completion.norm_coe, CStarRing.norm_star_mul_self]

/-- The completion of a C*-ring satisfies the C*-identity. -/
@[reducible] def completionCStarRing [CStarRing A] : CStarRing (Completion A) :=
  ⟨fun x => le_of_eq (completion_norm_star_mul_self x).symm⟩

/-- The star module law transports to the completion. -/
@[reducible] def completionStarModule [NormedSpace ℂ A] [StarModule ℂ A] :
    StarModule ℂ (Completion A) :=
  ⟨fun c x => by
    induction x using Completion.induction_on with
    | hp =>
      exact isClosed_eq (completion_map_star_continuous.comp (continuous_const_smul c))
        ((continuous_const_smul (star c)).comp completion_map_star_continuous)
    | ih a =>
      rw [← Completion.coe_smul, completion_star_coe (A := A), completion_star_coe (A := A),
        ← Completion.coe_smul, star_smul]⟩

scoped[OPH.QFT.CompletionStar]
  attribute [instance] completionContinuousStar completionNormedStarGroup
    completionCStarRing completionStarModule

end CompletionStar

/- The scope carries the completion star instances through the remaining
sections of this module. -/
open scoped OPH.QFT.CompletionStar

/-! ## The completion of the colimit as a C*-algebra -/

section ColimitCompletion

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- Mathlib registers the normed algebra law on completions only for a
commutative base; the noncommutative case is the same one-line proof. -/
noncomputable instance completionNormedAlgebra :
    NormedAlgebra ℂ (Completion (TowerColimit T)) where
  norm_smul_le := norm_smul_le

/-- **The completion of the colimit is a C*-algebra.** -/
noncomputable instance completionCStarAlgebra :
    CStarAlgebra (Completion (TowerColimit T)) where

/-- The canonical map from the colimit into its completion, as a star
algebra homomorphism. -/
noncomputable def completionMkHom (T : ConsensusTower ι) :
    TowerColimit T →⋆ₐ[ℂ] Completion (TowerColimit T) where
  toFun := ((↑) : TowerColimit T → Completion (TowerColimit T))
  map_one' := Completion.coe_one _
  map_mul' := Completion.coe_mul
  map_zero' := Completion.coe_zero
  map_add' := Completion.coe_add
  commutes' _ := rfl
  map_star' a := (completion_star_coe a).symm

/-- The canonical map acts by the coercion into the completion. -/
theorem completionMkHom_apply (x : TowerColimit T) :
    completionMkHom T x = (x : Completion (TowerColimit T)) :=
  rfl

/-- The canonical map into the completion is isometric. -/
theorem completionMkHom_isometry : Isometry (completionMkHom T) :=
  Completion.coe_isometry

/-- The canonical map into the completion is injective. -/
theorem completionMkHom_injective : Function.Injective (completionMkHom T) :=
  completionMkHom_isometry.injective

/-- The canonical map into the completion has dense range. -/
theorem completionMkHom_denseRange : DenseRange (completionMkHom T) :=
  Completion.denseRange_coe

/-- The canonical map into the completion preserves the norm. -/
theorem norm_completionMkHom (x : TowerColimit T) : ‖completionMkHom T x‖ = ‖x‖ :=
  Completion.norm_coe x

/-- The cocone map from a stage into the completion. -/
noncomputable def stageToCompletion (T : ConsensusTower ι) (r : ι) :
    ConsensusTower.PrivateAlgebra T r →⋆ₐ[ℂ] Completion (TowerColimit T) :=
  (completionMkHom T).comp (colimitMkHom T r)

/-- The stage cocone map acts by the colimit class of its argument. -/
theorem stageToCompletion_apply (r : ι) (X : ConsensusTower.PrivateAlgebra T r) :
    stageToCompletion T r X = (colimitMk X : Completion (TowerColimit T)) :=
  rfl

/-- Refinement independence of the stage cocone maps into the completion. -/
theorem stageToCompletion_refine {r s : ι} (hrs : r ≤ s)
    (X : ConsensusTower.PrivateAlgebra T r) :
    stageToCompletion T s (T.algebraRefine hrs X) = stageToCompletion T r X := by
  rw [stageToCompletion_apply, stageToCompletion_apply, colimitMk_refine]

/-- The completion norm of a stage image is the limit seminorm of its colimit
class. -/
theorem norm_stageToCompletion (r : ι) (X : ConsensusTower.PrivateAlgebra T r) :
    ‖stageToCompletion T r X‖ = colimitNorm (colimitMk X) :=
  Completion.norm_coe _

/-- The stage cocone map into the completion is norm nonincreasing. -/
theorem norm_stageToCompletion_le (r : ι) (X : ConsensusTower.PrivateAlgebra T r) :
    ‖stageToCompletion T r X‖ ≤ ‖X‖ := by
  rw [norm_stageToCompletion]
  exact colimitNorm_mk_le X

/-- The image of the colimit in the completion is the union of the stage
images. -/
theorem range_completionMkHom_eq_iUnion :
    Set.range (completionMkHom T) = ⋃ r : ι, Set.range (stageToCompletion T r) := by
  ext z
  simp only [Set.mem_range, Set.mem_iUnion]
  constructor
  · rintro ⟨x, rfl⟩
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    exact ⟨r, X, rfl⟩
  · rintro ⟨r, X, rfl⟩
    exact ⟨colimitMk X, rfl⟩

/-- The stage images are jointly dense in the completion. -/
theorem dense_iUnion_range_stageToCompletion :
    Dense (⋃ r : ι, Set.range (stageToCompletion T r)) := by
  rw [← range_completionMkHom_eq_iUnion]
  exact completionMkHom_denseRange

end ColimitCompletion

/-! ## Transport of isotony, locality, refinement independence, time slice -/

section Transport

/-- Commutation passes to closures: two nested closed-set arguments. -/
theorem commute_of_mem_closure {B : Type*} [TopologicalSpace B] [T2Space B] [Mul B]
    [ContinuousMul B]
    {s t : Set B} (h : ∀ a ∈ s, ∀ b ∈ t, Commute a b)
    {x y : B} (hx : x ∈ closure s) (hy : y ∈ closure t) : Commute x y := by
  have h1 : ∀ b ∈ t, Commute x b := by
    intro b hb
    have hcl : IsClosed {a : B | a * b = b * a} :=
      isClosed_eq (continuous_mul_const b) (continuous_const_mul b)
    exact closure_minimal (fun a ha => h a ha b hb) hcl hx
  have hcl : IsClosed {b : B | x * b = b * x} :=
    isClosed_eq (continuous_const_mul x) (continuous_mul_const x)
  exact closure_minimal (fun b hb => h1 b hb) hcl hy

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}
variable {N : FiniteCausalObserverNet T}

/-- The local algebra of an admissible event region in the completion: the
closure of the image of the committed limit local algebra. -/
noncomputable def completionLocalAlgebra (N : FiniteCausalObserverNet T)
    (U : EventRegion N) : StarSubalgebra ℂ (Completion (TowerColimit T)) :=
  ((limitLocalAlgebra N U).map (completionMkHom T)).topologicalClosure

/-- Limit-level regional observables land in the transported local algebra of
their event region. -/
theorem completionMkHom_mem_completionLocalAlgebra (U : EventRegion N)
    {x : TowerColimit T} (hx : x ∈ limitLocalAlgebra N U) :
    completionMkHom T x ∈ completionLocalAlgebra N U :=
  StarSubalgebra.le_topologicalClosure _ (StarSubalgebra.mem_map.mpr ⟨x, hx, rfl⟩)

/-- **Isotony in the completion.** -/
theorem completionLocalAlgebra_mono {U V : EventRegion N} (f : U ⟶ V) :
    completionLocalAlgebra N U ≤ completionLocalAlgebra N V :=
  StarSubalgebra.topologicalClosure_mono
    (StarSubalgebra.map_mono (limitLocalAlgebra_mono f))

/-- **Refinement independence in the completion.** -/
theorem completionLocalAlgebra_refine {r s : ι} (hrs : r ≤ s) (U : N.Region r) :
    completionLocalAlgebra N ⟨s, N.regionRefine hrs U⟩ =
      completionLocalAlgebra N ⟨r, U⟩ := by
  unfold completionLocalAlgebra
  rw [limitLocalAlgebra_refine]

/-- **Locality in the completion.**  Closures of images of algebras of
declared-disjoint regions commute elementwise. -/
theorem completionLocalAlgebra_commute {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V) {x y : Completion (TowerColimit T)}
    (hx : x ∈ completionLocalAlgebra N ⟨r, U⟩)
    (hy : y ∈ completionLocalAlgebra N ⟨r, V⟩) :
    Commute x y := by
  have hx' : x ∈ closure
      (((limitLocalAlgebra N ⟨r, U⟩).map (completionMkHom T) :
        StarSubalgebra ℂ (Completion (TowerColimit T))) : Set (Completion (TowerColimit T))) := hx
  have hy' : y ∈ closure
      (((limitLocalAlgebra N ⟨r, V⟩).map (completionMkHom T) :
        StarSubalgebra ℂ (Completion (TowerColimit T))) : Set (Completion (TowerColimit T))) := hy
  refine commute_of_mem_closure ?_ hx' hy'
  intro a ha b hb
  obtain ⟨a₀, ha₀, rfl⟩ := StarSubalgebra.mem_map.mp ha
  obtain ⟨b₀, hb₀, rfl⟩ := StarSubalgebra.mem_map.mp hb
  exact (limitLocalAlgebra_commute hUV ha₀ hb₀).map (completionMkHom T)

/-- **The committed Cauchy-embedding equality transports to the
closures.**  The committed order-theoretic statement is an equality of
limit local algebras; closure of the image is a function of that
equality. -/
theorem completionLocalAlgebra_eq_of_isCauchyEmbedding {U V : EventRegion N}
    (f : U ⟶ V) (hf : IsCauchyEmbedding f) :
    completionLocalAlgebra N U = completionLocalAlgebra N V := by
  unfold completionLocalAlgebra
  rw [limitLocalAlgebra_eq_of_isCauchyEmbedding f hf]

end Transport

/-! ## Nonvacuity on the committed witness tower -/

section Witness

/-- On the witness tower every connecting map is the identity. -/
theorem witnessTower_algebraRefine (h : () ≤ ())
    (X : ConsensusTower.PrivateAlgebra witnessTower ()) :
    witnessTower.algebraRefine h X = X :=
  witnessTower.algebra_refine_refl () X

/-- Every representative of the witness tower is stable. -/
theorem witnessTower_isStableRep (X : ConsensusTower.PrivateAlgebra witnessTower ()) :
    IsStableRep (T := witnessTower) X := by
  intro t h
  rw [witnessTower_algebraRefine]

/-- On the witness tower the limit seminorm of a class is the stage norm of
its representative. -/
theorem witnessTower_colimitNorm_mk (X : ConsensusTower.PrivateAlgebra witnessTower ()) :
    colimitNorm (colimitMk (T := witnessTower) X) = ‖X‖ :=
  (witnessTower_isStableRep X).colimitNorm_mk

/-- The separating projector is a projection. -/
theorem e0_star_mul_self : star e0 * e0 = e0 := by
  unfold e0
  rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  by_cases hi : i = 0 <;> simp [hi]

/-- The separating projector is nonzero. -/
theorem e0_ne_zero : e0 ≠ 0 := by
  intro h
  have := congrFun (congrFun h 0) 0
  simp [e0] at this

/-- A nonzero projection has operator norm one. -/
theorem norm_e0 : ‖e0‖ = 1 := by
  have h1 : ‖e0‖ * ‖e0‖ = ‖e0‖ := by
    rw [← CStarRing.norm_star_mul_self (x := e0), e0_star_mul_self]
  have h0 : ‖e0‖ ≠ 0 := norm_ne_zero_iff.mpr e0_ne_zero
  exact mul_left_cancel₀ h0 (h1.trans (mul_one _).symm)

/-- **Nonvacuity.**  The class of the separating projector has limit
seminorm one on the witness tower. -/
theorem witnessTower_colimitNorm_e0 :
    colimitNorm (colimitMk (T := witnessTower) (r := ()) e0) = 1 := by
  rw [witnessTower_colimitNorm_mk]
  exact norm_e0

/-- The limit seminorm of the separating projector is positive. -/
theorem witnessTower_colimitNorm_e0_pos :
    0 < colimitNorm (colimitMk (T := witnessTower) (r := ()) e0) := by
  rw [witnessTower_colimitNorm_e0]
  exact zero_lt_one

/-- The unit of the witness colimit has norm one. -/
theorem witnessTower_norm_one : ‖(1 : TowerColimit witnessTower)‖ = 1 := by
  rw [← colimitMk_one (), norm_def, witnessTower_colimitNorm_mk]
  exact (norm_one : ‖(1 : Matrix (Fin 2) (Fin 2) ℂ)‖ = 1)

/-- **Nonvacuity at the completion level.**  The image of the separating
projector in the completion of the witness colimit has norm one. -/
theorem witnessTower_norm_stageToCompletion_e0 :
    ‖stageToCompletion witnessTower () e0‖ = 1 := by
  rw [norm_stageToCompletion, witnessTower_colimitNorm_e0]

end Witness

/-! ## Nonvacuity of the transported hypotheses

The four transported statements of the section `Transport` quantify over
nets, regions, disjoint pairs, and Cauchy embeddings.  The committed
support-graded net supplies a witness for each quantifier, and the
transported local algebra of its top region carries an element of norm
one. -/

section TransportWitness

/-- **Nonvacuity of the transported Cauchy-embedding equality.**  The
committed singleton-to-full embedding instantiates
`completionLocalAlgebra_eq_of_isCauchyEmbedding`. -/
theorem supportGradedNet_completionLocalAlgebra_top_eq_full :
    completionLocalAlgebra supportGradedNet topRegion =
      completionLocalAlgebra supportGradedNet fullRegion :=
  completionLocalAlgebra_eq_of_isCauchyEmbedding topToFull topToFull_isCauchyEmbedding

/-- **Nonvacuity of the transported isotony.**  The committed
bottom-to-top embedding instantiates `completionLocalAlgebra_mono`. -/
theorem supportGradedNet_completionLocalAlgebra_bottom_le_top :
    completionLocalAlgebra supportGradedNet bottomRegion ≤
      completionLocalAlgebra supportGradedNet topRegion :=
  completionLocalAlgebra_mono bottomToTop

/-- **Nonvacuity of the transported locality.**  The committed disjoint
pair of region labels instantiates `completionLocalAlgebra_commute`. -/
theorem supportGradedNet_completionLocalAlgebra_commute
    {x y : Completion (TowerColimit witnessTower)}
    (hx : x ∈ completionLocalAlgebra supportGradedNet ⟨(), ({0} : Finset (Fin 2))⟩)
    (hy : y ∈ completionLocalAlgebra supportGradedNet ⟨(), ({1} : Finset (Fin 2))⟩) :
    Commute x y :=
  completionLocalAlgebra_commute supportGradedNet_has_disjoint_pair hx hy

/-- The class of the separating projector lies in the transported local
algebra of the top region. -/
theorem supportGradedNet_completionMkHom_e0_mem :
    completionMkHom witnessTower (colimitMk (T := witnessTower) (r := ()) e0) ∈
      completionLocalAlgebra supportGradedNet topRegion :=
  completionMkHom_mem_completionLocalAlgebra topRegion
    (colimitMk_mem_limitLocalAlgebra topRegion ⟨e0, e0_mem_singleton⟩)

/-- **The transported local algebra of the top region is nontrivial.**  It
carries an element of norm one. -/
theorem supportGradedNet_norm_completionMkHom_e0 :
    ‖completionMkHom witnessTower (colimitMk (T := witnessTower) (r := ()) e0)‖ = 1 := by
  rw [norm_completionMkHom, norm_def, witnessTower_colimitNorm_e0]

end TransportWitness

/-! ## Controls -/

section Controls

/-- **Control 1: the kernel clause has content.**  On the collapsing tower
the coarse unit is a nonzero stage element of stage norm one whose limit
seminorm is zero. -/
theorem collapsingTower_colimitNorm_one :
    colimitNorm (colimitMk (1 : ConsensusTower.PrivateAlgebra collapsingTower false)) = 0 ∧
      (1 : ConsensusTower.PrivateAlgebra collapsingTower false) ≠ 0 ∧
      ‖(1 : ConsensusTower.PrivateAlgebra collapsingTower false)‖ = 1 :=
  ⟨(colimitNorm_mk_eq_zero_iff_mk_eq_zero _).mpr collapsingTower_one_refinementNull.1,
    collapsingTower_one_refinementNull.2,
    (norm_one : ‖(1 : Matrix (Fin 1) (Fin 1) ℂ)‖ = 1)⟩

/-- The diagonal scaling unit of the two-by-two stage. -/
noncomputable def scalingUnit : (Matrix (Fin 2) (Fin 2) ℂ)ˣ where
  val := Matrix.diagonal fun i => if i = 0 then (2 : ℂ) else 1
  inv := Matrix.diagonal fun i => if i = 0 then (2 : ℂ)⁻¹ else 1
  val_inv := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    funext i
    by_cases hi : i = 0 <;> simp [hi]
  inv_val := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    funext i
    by_cases hi : i = 0 <;> simp [hi]

/-- Conjugation by a unit as a unital algebra homomorphism. -/
noncomputable def conjAlgHom (D : (Matrix (Fin 2) (Fin 2) ℂ)ˣ) :
    Matrix (Fin 2) (Fin 2) ℂ →ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ where
  toFun A := (D : Matrix (Fin 2) (Fin 2) ℂ) * A * ((D⁻¹ : (Matrix (Fin 2) (Fin 2) ℂ)ˣ) : Matrix (Fin 2) (Fin 2) ℂ)
  map_one' := by rw [mul_one, Units.mul_inv]
  map_mul' A B := by
    simp only [mul_assoc]
    rw [Units.inv_mul_cancel_left]
  map_zero' := by rw [mul_zero, zero_mul]
  map_add' A B := by rw [mul_add, add_mul]
  commutes' c := by
    rw [Algebra.algebraMap_eq_smul_one, mul_smul_comm, smul_mul_assoc, mul_one, Units.mul_inv]

/-- Conjugation acts by the two-sided product with the unit and its inverse. -/
theorem conjAlgHom_apply (D : (Matrix (Fin 2) (Fin 2) ℂ)ˣ) (A : Matrix (Fin 2) (Fin 2) ℂ) :
    conjAlgHom D A =
      (D : Matrix (Fin 2) (Fin 2) ℂ) * A * ((D⁻¹ : (Matrix (Fin 2) (Fin 2) ℂ)ˣ) : Matrix (Fin 2) (Fin 2) ℂ) :=
  rfl

/-- The non-star control: conjugation by the diagonal scaling. -/
noncomputable def scalingConj : Matrix (Fin 2) (Fin 2) ℂ →ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ :=
  conjAlgHom scalingUnit

/-- The control doubles the off-diagonal matrix unit. -/
theorem scalingConj_single :
    scalingConj (Matrix.single (0 : Fin 2) (1 : Fin 2) (1 : ℂ)) = (2 : ℂ) • Matrix.single (0 : Fin 2) (1 : Fin 2) (1 : ℂ) := by
  have hval : (scalingUnit : Matrix (Fin 2) (Fin 2) ℂ) =
      Matrix.diagonal fun i => if i = 0 then (2 : ℂ) else 1 := rfl
  have hinv : ((scalingUnit⁻¹ : (Matrix (Fin 2) (Fin 2) ℂ)ˣ) : Matrix (Fin 2) (Fin 2) ℂ) =
      Matrix.diagonal fun i => if i = 0 then (2 : ℂ)⁻¹ else 1 := rfl
  rw [scalingConj, conjAlgHom_apply, hval, hinv]
  ext i j
  rw [Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.smul_apply, smul_eq_mul]
  fin_cases i <;> fin_cases j <;> simp [Matrix.single]

/-- The off-diagonal matrix unit is nonzero. -/
theorem single_ne_zero' : Matrix.single (0 : Fin 2) (1 : Fin 2) (1 : ℂ) ≠ (0 : Matrix (Fin 2) (Fin 2) ℂ) := by
  intro h
  have := congrFun (congrFun h 0) 1
  simp at this

/-- **Control 2: the star clause is load-bearing for contractivity.**  A
unital algebra homomorphism of the two-by-two stage that does not preserve
the star increases the norm of a matrix unit. -/
theorem scalingConj_not_contractive :
    ∃ X : Matrix (Fin 2) (Fin 2) ℂ, ‖X‖ < ‖scalingConj X‖ := by
  refine ⟨Matrix.single (0 : Fin 2) (1 : Fin 2) (1 : ℂ), ?_⟩
  rw [scalingConj_single, norm_smul, RCLike.norm_two]
  have hpos : 0 < ‖Matrix.single (0 : Fin 2) (1 : Fin 2) (1 : ℂ)‖ := norm_pos_iff.mpr single_ne_zero'
  linarith

/-- The control is not star preserving.  A star-preserving map with the
same action bundles as a star algebra homomorphism, the contractivity
lemma consumed by `algebraRefine_norm_le` applies to that bundle, and the
contractivity failure excludes it. -/
theorem scalingConj_not_star_preserving :
    ¬ ∀ X : Matrix (Fin 2) (Fin 2) ℂ, scalingConj (star X) = star (scalingConj X) := by
  intro hstar
  obtain ⟨X, hX⟩ := scalingConj_not_contractive
  let φ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ :=
    { scalingConj with map_star' := hstar }
  have hle : ‖φ X‖ ≤ ‖X‖ := NonUnitalStarAlgHom.norm_apply_le φ X
  have hφ : φ X = scalingConj X := rfl
  rw [hφ] at hle
  exact absurd hX (not_lt.mpr hle)

end Controls

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.matrixStageCStarAlgebra
#print axioms OPH.QFT.algebraRefine_norm_le
#print axioms OPH.QFT.germNorm_congr
#print axioms OPH.QFT.colimitNorm
#print axioms OPH.QFT.colimitNorm_mk
#print axioms OPH.QFT.norm_starAlgHom_mem_spectralSet
#print axioms OPH.QFT.spectralSet_finite
#print axioms OPH.QFT.stageValues_finite
#print axioms OPH.QFT.exists_stage_norm_eq
#print axioms OPH.QFT.exists_stableRep
#print axioms OPH.QFT.exists_stableRep₂
#print axioms OPH.QFT.colimitNorm_zero
#print axioms OPH.QFT.colimitNorm_add_le
#print axioms OPH.QFT.colimitNorm_mul_le
#print axioms OPH.QFT.colimitNorm_smul
#print axioms OPH.QFT.colimitNorm_neg
#print axioms OPH.QFT.colimitNorm_star
#print axioms OPH.QFT.colimitNorm_star_mul_self
#print axioms OPH.QFT.colimitNorm_mk_eq_zero_iff
#print axioms OPH.QFT.colimitNorm_mk_eq_zero_iff_mk_eq_zero
#print axioms OPH.QFT.colimitNorm_eq_zero_iff
#print axioms OPH.QFT.instTowerColimitNorm
#print axioms OPH.QFT.norm_def
#print axioms OPH.QFT.colimitNormCore
#print axioms OPH.QFT.instTowerColimitNormedRing
#print axioms OPH.QFT.instTowerColimitNormedSpace
#print axioms OPH.QFT.instTowerColimitNormedAlgebra
#print axioms OPH.QFT.instTowerColimitNormedStarGroup
#print axioms OPH.QFT.instTowerColimitCStarRing
#print axioms OPH.QFT.colimit_t2Space
#print axioms OPH.QFT.completionStarRing
#print axioms OPH.QFT.completion_star_coe
#print axioms OPH.QFT.completionContinuousStar
#print axioms OPH.QFT.completion_norm_star
#print axioms OPH.QFT.completionNormedStarGroup
#print axioms OPH.QFT.completion_norm_star_mul_self
#print axioms OPH.QFT.completionCStarRing
#print axioms OPH.QFT.completionStarModule
#print axioms OPH.QFT.completionNormedAlgebra
#print axioms OPH.QFT.completionCStarAlgebra
#print axioms OPH.QFT.completionMkHom
#print axioms OPH.QFT.completionMkHom_isometry
#print axioms OPH.QFT.completionMkHom_injective
#print axioms OPH.QFT.completionMkHom_denseRange
#print axioms OPH.QFT.norm_completionMkHom
#print axioms OPH.QFT.stageToCompletion
#print axioms OPH.QFT.stageToCompletion_refine
#print axioms OPH.QFT.norm_stageToCompletion
#print axioms OPH.QFT.norm_stageToCompletion_le
#print axioms OPH.QFT.range_completionMkHom_eq_iUnion
#print axioms OPH.QFT.dense_iUnion_range_stageToCompletion
#print axioms OPH.QFT.commute_of_mem_closure
#print axioms OPH.QFT.completionLocalAlgebra
#print axioms OPH.QFT.completionLocalAlgebra_mono
#print axioms OPH.QFT.completionLocalAlgebra_refine
#print axioms OPH.QFT.completionLocalAlgebra_commute
#print axioms OPH.QFT.completionLocalAlgebra_eq_of_isCauchyEmbedding
#print axioms OPH.QFT.witnessTower_colimitNorm_mk
#print axioms OPH.QFT.witnessTower_colimitNorm_e0
#print axioms OPH.QFT.witnessTower_colimitNorm_e0_pos
#print axioms OPH.QFT.witnessTower_norm_stageToCompletion_e0
#print axioms OPH.QFT.witnessTower_norm_one
#print axioms OPH.QFT.supportGradedNet_completionLocalAlgebra_top_eq_full
#print axioms OPH.QFT.supportGradedNet_completionLocalAlgebra_bottom_le_top
#print axioms OPH.QFT.supportGradedNet_completionLocalAlgebra_commute
#print axioms OPH.QFT.supportGradedNet_completionMkHom_e0_mem
#print axioms OPH.QFT.supportGradedNet_norm_completionMkHom_e0
#print axioms OPH.QFT.collapsingTower_colimitNorm_one
#print axioms OPH.QFT.scalingConj_not_contractive
#print axioms OPH.QFT.scalingConj_not_star_preserving
