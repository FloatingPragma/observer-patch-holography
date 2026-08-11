import QFT.JointSlotFactorisation
import QFT.NoncommutativeWitness

/-!
# Regional nets with conditional-expectation restrictions

The committed `FiniteCausalObserverNet` interface carries contravariant
star-homomorphic restriction maps.  Above a scalar bottom region that
interface is uninhabitable for full matrix factors: the repository's own
obstruction (`no_scalar_hom_of_matrix`) excludes every unital algebra
homomorphism from a matrix algebra of dimension at least two to the
scalars.  The E1 continuation therefore names a principled redesign:
replace the impossible scalar-character restrictions by state-compatible
conditional expectations.  This module defines that interface.

A `CPRegionalNet` carries a finite region poset with declared overlaps
and disjointness, region-indexed star subalgebras with isotony and
locality, and, in place of restriction homomorphisms, one conditional
expectation per region: a linear map into the regional algebra that
fixes it pointwise, preserves positive semidefiniteness and the trace,
and satisfies the tower law along region inclusion.  Idempotence,
unitality, and the retraction of isotony are theorems of the interface
rather than extra fields.  A `coverage` field states the genuine joint
law the earlier `FiniteCover` type lacked: the regional algebras of a
declared generating family generate the top algebra.

The module also proves the generic two-slot coverage theorem (the left
and right slot embeddings of a product index generate the full matrix
algebra) and packages the obstruction as the redesign rationale: a
regional algebra carrying an embedded full n-by-n matrix factor with
2 <= n admits no unital algebra homomorphism to the scalars, so the
star-homomorphic restriction interface cannot be inhabited over such
regions with a scalar bottom.

The interface fields demand positivity and trace preservation; the
name's complete-positivity reading is certified at the witness level by
explicit matrix-unit Kraus identities, and no interface field states
amplified positivity.  The strength of the coverage law lives in the
declared generating family: a family containing a region with the full
algebra discharges it trivially, and the two-slot witness declares the
proper pair of slot regions instead.

**Boundary.**  This is a finite algebraic interface with a declared
region system; no tower stage, spacetime, causal relation, channel
semantics, instrument, clock, or physical claim is attached, and no
source object is constructed here.  Issue #692 remains open for a
source-attached regional assembly, correlation or descent receipt, and
nonconstant source realization.  Closed issue #712 consumes the older
finite-net interface and does not supply this migration.  Instrument,
continuum, and clock scopes remain at #693, #700, and #703.
-/

namespace OPH.QFT

open Matrix
open Kronecker
open scoped ComplexOrder

noncomputable section

variable {γ : Type*} [Fintype γ] [DecidableEq γ]

/-- A finite regional net whose restriction data are conditional
expectations rather than star homomorphisms. -/
structure CPRegionalNet (γ : Type*) [Fintype γ] [DecidableEq γ] where
  /-- Finite observer-accessible regions. -/
  Region : Type
  regionFintype : Fintype Region
  regionNonempty : Nonempty Region
  /-- Region inclusion with its poset laws. -/
  regionLE : Region → Region → Prop
  regionLE_refl : ∀ U, regionLE U U
  regionLE_trans : ∀ {U V W}, regionLE U V → regionLE V W → regionLE U W
  regionLE_antisymm : ∀ {U V}, regionLE U V → regionLE V U → U = V
  /-- A declared overlap region with the meet universal property. -/
  overlap : Region → Region → Region
  overlap_le_left : ∀ U V, regionLE (overlap U V) U
  overlap_le_right : ∀ U V, regionLE (overlap U V) V
  le_overlap : ∀ {W U V}, regionLE W U → regionLE W V →
    regionLE W (overlap U V)
  /-- Declared disjointness; no geometric reading is built in. -/
  disjoint : Region → Region → Prop
  disjoint_symm : ∀ {U V}, disjoint U V → disjoint V U
  disjoint_irrefl : ∀ U, ¬ disjoint U U
  /-- Region-indexed star subalgebras of one finite matrix algebra. -/
  localAlgebra : Region → StarSubalgebra ℂ (Matrix γ γ ℂ)
  /-- Isotony. -/
  isotony : ∀ {U V}, regionLE U V → localAlgebra U ≤ localAlgebra V
  /-- Algebraic locality on the declared disjointness relation. -/
  locality : ∀ {U V}, disjoint U V →
    ∀ X ∈ localAlgebra U, ∀ Y ∈ localAlgebra V, Commute X Y
  /-- **The redesigned restriction data**: one conditional expectation
  per region. -/
  expect : Region → Matrix γ γ ℂ →ₗ[ℂ] Matrix γ γ ℂ
  expect_mem : ∀ U X, expect U X ∈ localAlgebra U
  expect_fixes : ∀ U, ∀ X ∈ localAlgebra U, expect U X = X
  expect_posSemidef : ∀ U {X : Matrix γ γ ℂ},
    X.PosSemidef → (expect U X).PosSemidef
  expect_trace : ∀ U X, (expect U X).trace = X.trace
  /-- The tower law along region inclusion. -/
  expect_tower : ∀ {U V}, regionLE U V → ∀ X,
    expect U (expect V X) = expect U X
  /-- The genuine joint coverage law: a declared generating family of
  regions generates the top algebra. -/
  generating : Finset Region
  coverage : StarAlgebra.adjoin ℂ
    (⋃ U ∈ generating, (localAlgebra U : Set (Matrix γ γ ℂ))) = ⊤

attribute [instance] CPRegionalNet.regionFintype CPRegionalNet.regionNonempty

namespace CPRegionalNet

variable {N : CPRegionalNet γ}

/-- Conditional expectations are idempotent. -/
theorem expect_idem (U : N.Region) (X : Matrix γ γ ℂ) :
    N.expect U (N.expect U X) = N.expect U X :=
  N.expect_fixes U _ (N.expect_mem U X)

/-- Conditional expectations are unital. -/
theorem expect_one (U : N.Region) : N.expect U 1 = 1 :=
  N.expect_fixes U 1 (one_mem _)

/-- The expectation of a smaller region retracts the isotony inclusion:
observables of the smaller region are untouched. -/
theorem expect_retracts {U V : N.Region} (_hUV : N.regionLE U V)
    (X : Matrix γ γ ℂ) (hX : X ∈ N.localAlgebra U) :
    N.expect U X = X :=
  N.expect_fixes U X hX

/-- Restriction along inclusion, as the composite map into the smaller
algebra.  This replaces the star-homomorphic `restrict` field. -/
def restrictCP (U : N.Region) :
    Matrix γ γ ℂ →ₗ[ℂ] Matrix γ γ ℂ :=
  N.expect U

/-- The tower law in composition form. -/
theorem restrictCP_comp {U V : N.Region} (hUV : N.regionLE U V)
    (X : Matrix γ γ ℂ) :
    N.restrictCP U (N.restrictCP V X) = N.restrictCP U X :=
  N.expect_tower hUV X

end CPRegionalNet

/-! ## The generic two-slot coverage theorem -/

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]

omit [Fintype α] [Fintype β] in
/-- Kronecker product of two unit matrices is a joint unit matrix. -/
theorem single_kronecker_single (i k : α) (j l : β) :
    (Matrix.single i k (1 : ℂ)) ⊗ₖ (Matrix.single j l (1 : ℂ)) =
      Matrix.single ((i, j) : α × β) ((k, l) : α × β) 1 := by
  ext z z'
  simp only [Matrix.kroneckerMap_apply, Matrix.single, Matrix.of_apply]
  by_cases h1 : i = z.1 <;> by_cases h2 : k = z'.1 <;>
    by_cases h3 : j = z.2 <;> by_cases h4 : l = z'.2 <;>
      simp [h1, h2, h3, h4, Prod.ext_iff]

/-- **Joint coverage.**  The two slot embeddings of a product index
generate the full matrix algebra: every joint matrix unit is the product
of a left-slot unit and a right-slot unit. -/
theorem slot_ranges_generate_top :
    StarAlgebra.adjoin ℂ
        (((slotLeft β (α := α)).range : Set (Matrix (α × β) (α × β) ℂ)) ∪
          ((slotRight α (β := β)).range : Set (Matrix (α × β) (α × β) ℂ))) =
      ⊤ := by
  rw [StarSubalgebra.eq_top_iff]
  intro X
  rw [Matrix.matrix_eq_sum_single X]
  refine sum_mem fun p _ => sum_mem fun q _ => ?_
  have hsingle : Matrix.single p q (X p q) =
      X p q • Matrix.single p q (1 : ℂ) := by
    rw [Matrix.smul_single, smul_eq_mul, mul_one]
  rw [hsingle]
  refine SMulMemClass.smul_mem _ ?_
  have hfactor : Matrix.single p q (1 : ℂ) =
      (Matrix.single p.1 q.1 (1 : ℂ)) ⊗ₖ (Matrix.single p.2 q.2 (1 : ℂ)) := by
    rw [single_kronecker_single]
  have hleft : (Matrix.single p.1 q.1 (1 : ℂ)) ⊗ₖ (1 : Matrix β β ℂ) ∈
      StarAlgebra.adjoin ℂ
        (((slotLeft β (α := α)).range : Set (Matrix (α × β) (α × β) ℂ)) ∪
          ((slotRight α (β := β)).range : Set (Matrix (α × β) (α × β) ℂ))) :=
    StarAlgebra.subset_adjoin ℂ _ (Or.inl ⟨Matrix.single p.1 q.1 1, rfl⟩)
  have hright : (1 : Matrix α α ℂ) ⊗ₖ (Matrix.single p.2 q.2 (1 : ℂ)) ∈
      StarAlgebra.adjoin ℂ
        (((slotLeft β (α := α)).range : Set (Matrix (α × β) (α × β) ℂ)) ∪
          ((slotRight α (β := β)).range : Set (Matrix (α × β) (α × β) ℂ))) :=
    StarAlgebra.subset_adjoin ℂ _ (Or.inr ⟨Matrix.single p.2 q.2 1, rfl⟩)
  have hmul := mul_mem hleft hright
  rw [← Matrix.mul_kronecker_mul, mul_one, one_mul] at hmul
  rw [hfactor]
  exact hmul

/-! ## The redesign rationale -/

/-- The slot embedding is injective. -/
theorem slotLeft_injective [Nonempty β] :
    Function.Injective (slotLeft β (α := α)) := by
  intro A B hAB
  ext i j
  obtain ⟨b⟩ := ‹Nonempty β›
  have h := congrFun (congrFun hAB (i, b)) (j, b)
  simpa [slotLeft, Matrix.kroneckerMap_apply, Matrix.one_apply] using h

/-- **Redesign rationale.**  A regional algebra carrying an embedded
full n-by-n matrix factor with `2 <= n` admits no unital algebra
homomorphism to the scalars, so the star-homomorphic restriction
interface cannot be inhabited above a scalar bottom region over such
regional algebras.  The conditional-expectation interface carries the
same regional data without that obstruction. -/
theorem no_scalar_restriction_of_matrix_factor {n : ℕ} (hn : 2 ≤ n)
    (A : Type*) [Ring A] [Algebra ℂ A]
    (embed : Matrix (Fin n) (Fin n) ℂ →ₐ[ℂ] A)
    (φ : A →ₐ[ℂ] ℂ) : False :=
  no_scalar_hom_of_matrix n hn (φ.comp embed)

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.CPRegionalNet.expect_idem
#print axioms OPH.QFT.CPRegionalNet.expect_one
#print axioms OPH.QFT.CPRegionalNet.restrictCP_comp
#print axioms OPH.QFT.CPRegionalNet.expect_retracts
#print axioms OPH.QFT.slotLeft_injective
#print axioms OPH.QFT.single_kronecker_single
#print axioms OPH.QFT.slot_ranges_generate_top
#print axioms OPH.QFT.no_scalar_restriction_of_matrix_factor
