import QFT.TwoSlotCPNetWitness
import Tower.ConsensusTower

/-!
# The two-slot diamond anchored over a consensus-tower stage

The two-slot witness lives on the matrices indexed by the joint slot
space `Fin 13 × Fin 14`; the A3 tower interface indexes its private
algebra by `Fin (dim)`.  This module closes that gap.  A basis
equivalence `Fin 13 × Fin 14 ≃ Fin 182` induces a star-algebra
equivalence of the two ambient algebras, conditional-expectation nets
transport along any such equivalence (a generic functor with all
receipts carried), and one constant consensus-tower stage is built at
dimension 182 whose commutative public layer is spanned by the
reindexed lifted checkpoint partition of observer 86 and whose selected
state is the uniform density matrix.  The transported diamond is a
conditional-expectation net whose ambient algebra is definitionally the
tower stage's private algebra, every member of the tower's public
partition lies in the anchored left region, the anchored partition
members are proper, and the anchored left region is a proper
subalgebra.

**Boundary.**  The tower stage is the constant adaptor: one regulator,
one observer, the declared partition, the uniform state, and the zero
generator.  It carries no nonconstant refinement, no source-produced
tower dynamics, and no physical claim.  The slot assembly and label
conventions stay declared postprocessors over post-hoc transcriptions,
ineligible as validation.  The source attachment of the assembly and
any nonconstant tower realization stay open under issue #692.
-/

namespace OPH.QFT

open Matrix
open Kronecker
open OPH.Tower
open scoped ComplexOrder

noncomputable section

/-! ## The basis equivalence and the ambient star equivalence -/

/-- The basis equivalence of the joint slot space with `Fin 182`. -/
def pairFin : PairIndex ≃ Fin 182 :=
  finProdFinEquiv.trans (finCongr (by norm_num))

/-- Reindexing commutes with the conjugate transpose. -/
theorem reindex_conjTranspose {m n : Type*} [Fintype m] [Fintype n]
    [DecidableEq m] [DecidableEq n] (e : m ≃ n) (M : Matrix m m ℂ) :
    (Matrix.reindex e e M)ᴴ = Matrix.reindex e e Mᴴ := by
  ext i j
  simp [Matrix.reindex_apply, Matrix.conjTranspose_apply]

/-- The ambient star-algebra equivalence induced by the basis
equivalence. -/
def ambientEquiv :
    Matrix PairIndex PairIndex ℂ ≃⋆ₐ[ℂ] Matrix (Fin 182) (Fin 182) ℂ :=
  StarAlgEquiv.ofAlgEquiv (Matrix.reindexAlgEquiv ℂ ℂ pairFin)
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose pairFin M)

theorem ambientEquiv_apply (M : Matrix PairIndex PairIndex ℂ) :
    ambientEquiv M = Matrix.reindex pairFin pairFin M := rfl

/-- Reindexing along an equivalence preserves the trace. -/
theorem trace_reindex_equiv {m n : Type*} [Fintype m] [Fintype n]
    [DecidableEq m] [DecidableEq n] (e : m ≃ n) (M : Matrix m m ℂ) :
    (Matrix.reindex e e M).trace = M.trace := by
  unfold Matrix.trace
  rw [← Equiv.sum_comp e]
  simp [Matrix.reindex_apply, Matrix.diag]

/-- Reindexing along an equivalence preserves positive
semidefiniteness. -/
theorem posSemidef_reindex_equiv {m n : Type*} [Fintype m] [Fintype n]
    [DecidableEq m] [DecidableEq n] (e : m ≃ n)
    {M : Matrix m m ℂ} (hM : M.PosSemidef) :
    (Matrix.reindex e e M).PosSemidef := by
  have h : Matrix.reindex e e M = M.submatrix e.symm e.symm := rfl
  rw [h]
  exact hM.submatrix _

/-! ## Generic transport of conditional-expectation nets -/

variable {γ δ : Type*} [Fintype γ] [DecidableEq γ] [Fintype δ] [DecidableEq δ]

/-- A conditional-expectation net transports along a star-algebra
equivalence of ambient matrix algebras, provided the equivalence
preserves positive semidefiniteness in both directions and the trace in
the forward direction; the reverse trace direction is derived.  All
interface receipts are carried. -/
def CPRegionalNet.transport (N : CPRegionalNet γ)
    (Φ : Matrix γ γ ℂ ≃⋆ₐ[ℂ] Matrix δ δ ℂ)
    (hpos : ∀ {X : Matrix γ γ ℂ}, X.PosSemidef → (Φ X).PosSemidef)
    (hpos' : ∀ {Y : Matrix δ δ ℂ}, Y.PosSemidef → (Φ.symm Y).PosSemidef)
    (htrace : ∀ X : Matrix γ γ ℂ, (Φ X).trace = X.trace) :
    CPRegionalNet δ where
  Region := N.Region
  regionFintype := N.regionFintype
  regionNonempty := N.regionNonempty
  regionLE := N.regionLE
  regionLE_refl := N.regionLE_refl
  regionLE_trans := N.regionLE_trans
  regionLE_antisymm := N.regionLE_antisymm
  overlap := N.overlap
  overlap_le_left := N.overlap_le_left
  overlap_le_right := N.overlap_le_right
  le_overlap := N.le_overlap
  disjoint := N.disjoint
  disjoint_symm := N.disjoint_symm
  disjoint_irrefl := N.disjoint_irrefl
  localAlgebra U := StarSubalgebra.map (Φ : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix δ δ ℂ) (N.localAlgebra U)
  isotony hUV := StarSubalgebra.map_mono (N.isotony hUV)
  locality hUV X hX Y hY := by
    obtain ⟨X0, hX0, rfl⟩ := hX
    obtain ⟨Y0, hY0, rfl⟩ := hY
    have h := N.locality hUV X0 hX0 Y0 hY0
    unfold Commute SemiconjBy at h ⊢
    rw [← map_mul, ← map_mul, h]
  expect U :=
    { toFun := fun X => Φ (N.expect U (Φ.symm X))
      map_add' := fun X Y => by rw [map_add, map_add, map_add]
      map_smul' := fun c X => by
        rw [map_smul, map_smul, map_smul, RingHom.id_apply] }
  expect_mem U X := ⟨N.expect U (Φ.symm X), N.expect_mem U _, rfl⟩
  expect_fixes U X hX := by
    obtain ⟨X0, hX0, rfl⟩ := hX
    show Φ (N.expect U (Φ.symm (Φ X0))) = Φ X0
    rw [StarAlgEquiv.symm_apply_apply, N.expect_fixes U X0 hX0]
  expect_posSemidef U {X} hX := by
    show (Φ (N.expect U (Φ.symm X))).PosSemidef
    exact hpos (N.expect_posSemidef U (hpos' hX))
  expect_trace U X := by
    show (Φ (N.expect U (Φ.symm X))).trace = X.trace
    rw [htrace, N.expect_trace]
    have h := htrace (Φ.symm X)
    rw [StarAlgEquiv.apply_symm_apply] at h
    exact h.symm
  expect_tower hUV X := by
    show Φ (N.expect _ (Φ.symm (Φ (N.expect _ (Φ.symm X))))) =
      Φ (N.expect _ (Φ.symm X))
    rw [StarAlgEquiv.symm_apply_apply, N.expect_tower hUV]
  generating := N.generating
  coverage := by
    have himg : (⋃ U ∈ N.generating,
        ((StarSubalgebra.map (Φ : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix δ δ ℂ)
          (N.localAlgebra U)) : Set (Matrix δ δ ℂ))) =
        (Φ : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix δ δ ℂ) '' (⋃ U ∈ N.generating,
          (N.localAlgebra U : Set (Matrix γ γ ℂ))) := by
      rw [Set.image_iUnion₂]
      rfl
    rw [himg, ← StarAlgHom.map_adjoin, N.coverage]
    rw [StarSubalgebra.eq_top_iff]
    intro Y
    exact ⟨Φ.symm Y, trivial, StarAlgEquiv.apply_symm_apply Φ Y⟩

/-! ## The checkpoint partition at dimension 182 -/

/-- The checkpoint projectors are pairwise orthogonal. -/
theorem chk86_orthogonal {v w : Fin 4} (hvw : v ≠ w) :
    chk86 v * chk86 w = 0 := by
  unfold chk86
  rw [Matrix.diagonal_mul_diagonal]
  have hzero : (fun s => (if obs86.tbl s 1 = (v : ℕ) then (1 : ℂ) else 0) *
      (if obs86.tbl s 1 = (w : ℕ) then 1 else 0)) = fun _ => 0 := by
    funext s
    by_cases hv : obs86.tbl s 1 = (v : ℕ)
    · have hw : ¬ obs86.tbl s 1 = (w : ℕ) := by
        intro hw
        exact hvw (Fin.ext (by omega))
      rw [if_pos hv, if_neg hw, mul_zero]
    · rw [if_neg hv, zero_mul]
  rw [hzero, Matrix.diagonal_zero]

/-- The reindexed lifted checkpoint partition of observer 86: four
orthogonal projections summing to the identity at dimension 182. -/
def anchoredCheckpointPartition :
    EventAlgebra.ProjectivePartition 182 4 where
  proj v := ambientEquiv (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ))
  isEvent v := by
    constructor
    · show (ambientEquiv _)ᴴ = _
      rw [← Matrix.star_eq_conjTranspose, ← map_star]
      congr 1
      rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_kronecker,
        Matrix.conjTranspose_one, chk86_conjTranspose]
    · rw [← map_mul]
      congr 1
      rw [← Matrix.mul_kronecker_mul, one_mul, chk86_idem]
  orthogonal v w hvw := by
    rw [← map_mul, ← Matrix.mul_kronecker_mul, one_mul,
      chk86_orthogonal hvw, Matrix.zero_kronecker, map_zero]
  complete := by
    rw [← map_sum]
    have hsum : (∑ v : Fin 4,
        chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ)) =
        (∑ v : Fin 4, chk86 v) ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ) := by
      ext z z'
      rw [Matrix.sum_apply]
      simp [Matrix.kroneckerMap_apply, Matrix.sum_apply, Finset.sum_mul]
    rw [hsum, chk86_sum, Matrix.one_kronecker_one, map_one]

/-- The uniform state at dimension 182. -/
def anchoredState : EventAlgebra.StateMatrix 182 :=
  ⟨(182 : ℂ)⁻¹ • 1, by
    constructor
    · refine Matrix.PosSemidef.one.smul ?_
      rw [show ((182 : ℂ))⁻¹ = ((182 : ℝ)⁻¹ : ℂ) by push_cast; ring]
      exact_mod_cast inv_nonneg.mpr (by norm_num : (0:ℝ) ≤ 182)
    · rw [Matrix.trace_smul, Matrix.trace_one]
      norm_num⟩

/-- The anchored tower stage: the constant consensus tower over the
checkpoint partition and the uniform state. -/
def anchoredTower : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower anchoredCheckpointPartition
    anchoredState

/-- The tower stage's private algebra is the anchor ambient. -/
theorem anchoredTower_privateAlgebra :
    ConsensusTower.PrivateAlgebra anchoredTower () =
      Matrix (Fin 182) (Fin 182) ℂ := rfl

/-! ## The anchored diamond -/

/-- The two-slot diamond, transported onto the tower stage's private
algebra. -/
def anchoredNet : CPRegionalNet (Fin 182) :=
  twoSlotNet.transport ambientEquiv
    (fun hX => posSemidef_reindex_equiv pairFin hX)
    (fun hY => by
      have h := posSemidef_reindex_equiv pairFin.symm hY
      exact h)
    (fun X => trace_reindex_equiv pairFin X)

/-- Every member of the tower's public checkpoint partition lies in the
anchored left region. -/
theorem partition_members_in_left_region (v : Fin 4) :
    (anchoredCheckpointPartition.proj v) ∈
      anchoredNet.localAlgebra
        (show anchoredNet.Region from TwoSlotRegion.left) := by
  exact ⟨chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ), ⟨chk86 v, rfl⟩, rfl⟩

/-- The anchored partition members are proper: none vanishes and none
is the identity, so the tower's record layer carries four genuine
labels. -/
theorem anchoredCheckpointPartition_proper (v : Fin 4) :
    anchoredCheckpointPartition.proj v ≠ 0 ∧
      anchoredCheckpointPartition.proj v ≠ 1 := by
  obtain ⟨⟨s, hs⟩, ⟨s', hs'⟩⟩ := chk86_proper v
  constructor
  · intro h
    have hentry := congrFun (congrFun h (pairFin (s, 0))) (pairFin (s, 0))
    simp only [anchoredCheckpointPartition, ambientEquiv_apply,
      Matrix.reindex_apply, Matrix.submatrix_apply,
      Equiv.symm_apply_apply, Matrix.kroneckerMap_apply,
      Matrix.zero_apply] at hentry
    rw [Matrix.one_apply_eq] at hentry
    unfold chk86 at hentry
    rw [Matrix.diagonal_apply_eq, if_pos hs, one_mul] at hentry
    exact one_ne_zero hentry
  · intro h
    have hentry := congrFun (congrFun h (pairFin (s', 0))) (pairFin (s', 0))
    simp only [anchoredCheckpointPartition, ambientEquiv_apply,
      Matrix.reindex_apply, Matrix.submatrix_apply,
      Equiv.symm_apply_apply, Matrix.kroneckerMap_apply] at hentry
    rw [Matrix.one_apply_eq, Matrix.one_apply_eq] at hentry
    unfold chk86 at hentry
    rw [Matrix.diagonal_apply_eq, if_neg hs', zero_mul] at hentry
    exact zero_ne_one hentry

/-- The anchored left region is a proper subalgebra of the tower
stage's private algebra: the right-slot step operator lies outside
it. -/
theorem anchoredNet_left_ne_top :
    anchoredNet.localAlgebra
        (show anchoredNet.Region from TwoSlotRegion.left) ≠ ⊤ := by
  intro htop
  have hmem : ambientEquiv ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
      Matrix.single (0 : Fin 14) (1 : Fin 14) (1 : ℂ)) ∈
      anchoredNet.localAlgebra
        (show anchoredNet.Region from TwoSlotRegion.left) := by
    rw [htop]
    trivial
  obtain ⟨X, hX, hEq⟩ := hmem
  obtain ⟨A, rfl⟩ := hX
  have hAB : A ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ) =
      (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
        Matrix.single (0 : Fin 14) (1 : Fin 14) (1 : ℂ) :=
    ambientEquiv.injective hEq
  have hentry := congrFun (congrFun hAB ((0 : Fin 13), (0 : Fin 14)))
    ((0 : Fin 13), (1 : Fin 14))
  simp only [Matrix.kroneckerMap_apply] at hentry
  rw [Matrix.one_apply_eq, one_mul, Matrix.single_apply_same] at hentry
  have hoff : (1 : Matrix (Fin 14) (Fin 14) ℂ) 0 1 = 0 :=
    Matrix.one_apply_ne (by decide)
  rw [hoff, mul_zero] at hentry
  exact zero_ne_one hentry

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.CPRegionalNet.transport
#print axioms OPH.QFT.chk86_orthogonal
#print axioms OPH.QFT.anchoredCheckpointPartition
#print axioms OPH.QFT.anchoredState
#print axioms OPH.QFT.anchoredTower
#print axioms OPH.QFT.anchoredNet
#print axioms OPH.QFT.partition_members_in_left_region
#print axioms OPH.QFT.trace_reindex_equiv
#print axioms OPH.QFT.posSemidef_reindex_equiv
#print axioms OPH.QFT.anchoredCheckpointPartition_proper
#print axioms OPH.QFT.anchoredNet_left_ne_top
