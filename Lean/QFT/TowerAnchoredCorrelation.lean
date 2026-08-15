import QFT.TowerAnchoredDiamond
import QFT.SourceCorrelationCapstone

/-!
# The 86/247 correlation net anchored over a consensus-tower stage

`QFT.TowerAnchoredDiamond` transports the 86/88 diamond onto one
constant consensus-tower stage at dimension 182;
`QFT.SourceCorrelationCapstone` records (lines on the committed gap)
that no theorem transports the 86/247 diamond to a tower or relates the
two carriers.  This module removes that asymmetry with the same
construction and proves the first committed theorem relating the two
carriers.

**What is proved (finite matrix algebra on committed data only).**

* The basis equivalence `Fin 13 × Fin 13 ≃ Fin 169` induces an ambient
  star-algebra equivalence, and the committed 86/247 diamond
  `twoSlotNet247` transports along it through the committed generic
  `CPRegionalNet.transport`, exactly as the 86/88 diamond does at
  dimension 182.  One constant consensus-tower stage is built at
  dimension 169 whose commutative public layer is spanned by the
  reindexed lifted checkpoint partition of observer 86 — the **same**
  committed `chk86` family the 182-stage consumes — and whose selected
  state is the uniform density matrix.  The template's nontriviality
  receipts are carried: every partition member lies in the anchored
  left region, every partition member is proper, and the anchored left
  region is a proper subalgebra.  Both committed carriers now carry one
  common tower-anchored construction (`bothCarriers_common_anchoring`).
* **Hinge agreement.**  Observer 86 is the left observer of both
  committed pairs.  The 86/88 witness's left slot data are read from
  the committed `obs86` payload (its joint path's left components are
  definitionally `obs86.path`), and the 86/247 joint path's left
  components are kernel-checked equal to `obs86.path`
  (`jointPath247_marginals`).  So both sides marginalize the same
  committed `obs86` payload, and the exact agreement holds: the two
  committed joint paths have equal left label components, the induced
  left occupation count vectors are equal, the induced diagonal hinge
  states are equal and equal to the committed `marginal86State`, and
  the left partial trace of the committed counted correlation state is
  exactly the hinge state of the 86/88 side
  (`towerCarriers_share_obs86_marginal`).

**What is NOT proved.**  The two carriers are not identified: the
ambient dimensions 169 and 182 differ, and no map between the two
anchored nets or towers is constructed.  The hinge agreement is an
exact equality of left-marginal occupation data of the shared observer,
not a carrier or net identification.  Both tower stages are the
constant adaptor: one regulator, one observer, the declared partition,
the uniform state, and the zero generator; tower anchoring is
structural bookkeeping on committed finite data, with no nonconstant
refinement, no source-produced tower dynamics, and no physical claim.
The slot assembly and label conventions stay declared postprocessors
over post-hoc transcriptions, ineligible as validation (premise rows
PR-32, PR-33, PR-34 as bundled in the adequacy surface).  Physical
time, causality, and spacetime attachment stay open under PR-52; the
time-slice interface stays open under PR-58; the nonconstant source
realization is scoped to issue #728.
-/

set_option autoImplicit false
set_option maxRecDepth 65536

namespace OPH.QFT

open Matrix
open Kronecker
open OPH.Tower
open scoped ComplexOrder

noncomputable section

/-! ## The basis equivalence and the ambient star equivalence at 169 -/

/-- The basis equivalence of the 86/247 joint slot space with
`Fin 169`. -/
def pairFin247 : PairIndex247 ≃ Fin 169 :=
  finProdFinEquiv.trans (finCongr (by norm_num))

/-- The ambient star-algebra equivalence induced by the basis
equivalence, the 86/247 mirror of `ambientEquiv`. -/
def ambientEquiv247 :
    Matrix PairIndex247 PairIndex247 ℂ ≃⋆ₐ[ℂ] Matrix (Fin 169) (Fin 169) ℂ :=
  StarAlgEquiv.ofAlgEquiv (Matrix.reindexAlgEquiv ℂ ℂ pairFin247)
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose pairFin247 M)

theorem ambientEquiv247_apply (M : Matrix PairIndex247 PairIndex247 ℂ) :
    ambientEquiv247 M = Matrix.reindex pairFin247 pairFin247 M := rfl

/-! ## The checkpoint partition at dimension 169 -/

/-- The reindexed lifted checkpoint partition of observer 86 on the
86/247 carrier: four orthogonal projections summing to the identity at
dimension 169, built from the same committed `chk86` family as the
182-stage partition. -/
def anchoredCheckpointPartition247 :
    EventAlgebra.ProjectivePartition 169 4 where
  proj v := ambientEquiv247 (chk86 v ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ))
  isEvent v := by
    constructor
    · show (ambientEquiv247 _)ᴴ = _
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
        chk86 v ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) =
        (∑ v : Fin 4, chk86 v) ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ) := by
      ext z z'
      rw [Matrix.sum_apply]
      simp [Matrix.kroneckerMap_apply, Matrix.sum_apply, Finset.sum_mul]
    rw [hsum, chk86_sum, Matrix.one_kronecker_one, map_one]

/-- The uniform state at dimension 169. -/
def anchoredState247 : EventAlgebra.StateMatrix 169 :=
  ⟨(169 : ℂ)⁻¹ • 1, by
    constructor
    · refine Matrix.PosSemidef.one.smul ?_
      rw [show ((169 : ℂ))⁻¹ = ((169 : ℝ)⁻¹ : ℂ) by push_cast; ring]
      exact_mod_cast inv_nonneg.mpr (by norm_num : (0:ℝ) ≤ 169)
    · rw [Matrix.trace_smul, Matrix.trace_one]
      norm_num⟩

/-- The anchored tower stage at 169: the constant consensus tower over
the checkpoint partition and the uniform state. -/
def anchoredTower247 : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower anchoredCheckpointPartition247
    anchoredState247

/-- The 169-stage's private algebra is the anchor ambient. -/
theorem anchoredTower247_privateAlgebra :
    ConsensusTower.PrivateAlgebra anchoredTower247 () =
      Matrix (Fin 169) (Fin 169) ℂ := rfl

/-! ## The anchored 86/247 diamond -/

/-- The committed 86/247 diamond, transported onto the 169-stage's
private algebra through the committed generic
`CPRegionalNet.transport`. -/
def anchoredNet247 : CPRegionalNet (Fin 169) :=
  twoSlotNet247.transport ambientEquiv247
    (fun hX => posSemidef_reindex_equiv pairFin247 hX)
    (fun hY => by
      have h := posSemidef_reindex_equiv pairFin247.symm hY
      exact h)
    (fun X => trace_reindex_equiv pairFin247 X)

/-- Every member of the 169-stage's public checkpoint partition lies in
the anchored left region. -/
theorem partition_members_in_left_region247 (v : Fin 4) :
    (anchoredCheckpointPartition247.proj v) ∈
      anchoredNet247.localAlgebra
        (show anchoredNet247.Region from TwoSlotRegion.left) := by
  exact ⟨chk86 v ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ), ⟨chk86 v, rfl⟩, rfl⟩

/-- The anchored partition members are proper: none vanishes and none
is the identity, so the 169-stage's record layer carries four genuine
labels. -/
theorem anchoredCheckpointPartition247_proper (v : Fin 4) :
    anchoredCheckpointPartition247.proj v ≠ 0 ∧
      anchoredCheckpointPartition247.proj v ≠ 1 := by
  obtain ⟨⟨s, hs⟩, ⟨s', hs'⟩⟩ := chk86_proper v
  constructor
  · intro h
    have hentry := congrFun (congrFun h (pairFin247 (s, 0))) (pairFin247 (s, 0))
    simp only [anchoredCheckpointPartition247, ambientEquiv247_apply,
      Matrix.reindex_apply, Matrix.submatrix_apply,
      Equiv.symm_apply_apply, Matrix.kroneckerMap_apply,
      Matrix.zero_apply] at hentry
    rw [Matrix.one_apply_eq] at hentry
    unfold chk86 at hentry
    rw [Matrix.diagonal_apply_eq, if_pos hs, one_mul] at hentry
    exact one_ne_zero hentry
  · intro h
    have hentry := congrFun (congrFun h (pairFin247 (s', 0))) (pairFin247 (s', 0))
    simp only [anchoredCheckpointPartition247, ambientEquiv247_apply,
      Matrix.reindex_apply, Matrix.submatrix_apply,
      Equiv.symm_apply_apply, Matrix.kroneckerMap_apply] at hentry
    rw [Matrix.one_apply_eq, Matrix.one_apply_eq] at hentry
    unfold chk86 at hentry
    rw [Matrix.diagonal_apply_eq, if_neg hs', zero_mul] at hentry
    exact zero_ne_one hentry

/-- The anchored left region is a proper subalgebra of the 169-stage's
private algebra: the right-slot step operator lies outside it. -/
theorem anchoredNet247_left_ne_top :
    anchoredNet247.localAlgebra
        (show anchoredNet247.Region from TwoSlotRegion.left) ≠ ⊤ := by
  intro htop
  have hmem : ambientEquiv247 ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
      Matrix.single (0 : Fin 13) (1 : Fin 13) (1 : ℂ)) ∈
      anchoredNet247.localAlgebra
        (show anchoredNet247.Region from TwoSlotRegion.left) := by
    rw [htop]
    trivial
  obtain ⟨X, hX, hEq⟩ := hmem
  obtain ⟨A, rfl⟩ := hX
  have hAB : A ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ) =
      (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
        Matrix.single (0 : Fin 13) (1 : Fin 13) (1 : ℂ) :=
    ambientEquiv247.injective hEq
  have hentry := congrFun (congrFun hAB ((0 : Fin 13), (0 : Fin 13)))
    ((0 : Fin 13), (1 : Fin 13))
  simp only [Matrix.kroneckerMap_apply] at hentry
  rw [Matrix.one_apply_eq, one_mul, Matrix.single_apply_same] at hentry
  have hoff : (1 : Matrix (Fin 13) (Fin 13) ℂ) 0 1 = 0 :=
    Matrix.one_apply_ne (by decide)
  rw [hoff, mul_zero] at hentry
  exact zero_ne_one hentry

/-! ## Both carriers in one anchoring form -/

/-- Both anchored partitions consume the same committed `chk86`
checkpoint family of observer 86: each is that family, lifted into the
respective left slot and reindexed by the respective basis
equivalence. -/
theorem anchored_partitions_share_chk86 :
    (∀ v : Fin 4, anchoredCheckpointPartition.proj v =
        ambientEquiv (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ))) ∧
    (∀ v : Fin 4, anchoredCheckpointPartition247.proj v =
        ambientEquiv247 (chk86 v ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ))) :=
  ⟨fun _ => rfl, fun _ => rfl⟩

/-- **One common tower-anchored construction on both committed
carriers.**  The 86/88 diamond and the 86/247 diamond are each the
committed generic `CPRegionalNet.transport` of the respective committed
two-slot net along the respective ambient basis star-equivalence; each
anchor ambient is definitionally the private algebra of one constant
consensus-tower stage (dimensions 182 and 169); on both sides every
member of the stage's public checkpoint partition lies in the anchored
left region and the anchored left region is proper.  This removes the
asymmetry recorded in `QFT.SourceCorrelationCapstone`: the constant-
tower transport now exists for both committed carriers in the same
form.  The two carriers are **not** identified — the dimensions 169
and 182 differ and no map between the anchored nets is constructed. -/
theorem bothCarriers_common_anchoring :
    (anchoredNet = twoSlotNet.transport ambientEquiv
        (fun hX => posSemidef_reindex_equiv pairFin hX)
        (fun hY => posSemidef_reindex_equiv pairFin.symm hY)
        (fun X => trace_reindex_equiv pairFin X)) ∧
    (anchoredNet247 = twoSlotNet247.transport ambientEquiv247
        (fun hX => posSemidef_reindex_equiv pairFin247 hX)
        (fun hY => posSemidef_reindex_equiv pairFin247.symm hY)
        (fun X => trace_reindex_equiv pairFin247 X)) ∧
    (ConsensusTower.PrivateAlgebra anchoredTower () =
        Matrix (Fin 182) (Fin 182) ℂ) ∧
    (ConsensusTower.PrivateAlgebra anchoredTower247 () =
        Matrix (Fin 169) (Fin 169) ℂ) ∧
    (∀ v : Fin 4, anchoredCheckpointPartition.proj v ∈
        anchoredNet.localAlgebra
          (show anchoredNet.Region from TwoSlotRegion.left)) ∧
    (∀ v : Fin 4, anchoredCheckpointPartition247.proj v ∈
        anchoredNet247.localAlgebra
          (show anchoredNet247.Region from TwoSlotRegion.left)) ∧
    (anchoredNet.localAlgebra
        (show anchoredNet.Region from TwoSlotRegion.left) ≠ ⊤) ∧
    (anchoredNet247.localAlgebra
        (show anchoredNet247.Region from TwoSlotRegion.left) ≠ ⊤) :=
  ⟨rfl, rfl, anchoredTower_privateAlgebra, anchoredTower247_privateAlgebra,
    partition_members_in_left_region, partition_members_in_left_region247,
    anchoredNet_left_ne_top, anchoredNet247_left_ne_top⟩

/-! ## The shared-observer hinge agreement -/

/-- The left occupation count read from the committed 86/88 joint
path. -/
def leftHingeCount88 (i : Fin 13) : ℕ :=
  (Finset.univ.filter fun t : Fin 32 => (jointPath t).1 = i).card

/-- The left occupation count read from the committed 86/247 joint
path. -/
def leftHingeCount247 (i : Fin 13) : ℕ :=
  (Finset.univ.filter fun t : Fin 32 => (jointPath247 t).1 = i).card

/-- The diagonal hinge state induced by the 86/88 left counts. -/
def leftHingeState88 : Matrix (Fin 13) (Fin 13) ℂ :=
  Matrix.diagonal fun i => (leftHingeCount88 i : ℂ) / 32

/-- The diagonal hinge state induced by the 86/247 left counts. -/
def leftHingeState247 : Matrix (Fin 13) (Fin 13) ℂ :=
  Matrix.diagonal fun i => (leftHingeCount247 i : ℂ) / 32

/-- The two committed joint paths carry equal left label components:
the 86/88 side definitionally, the 86/247 side by the committed kernel
check `jointPath247_marginals`. -/
theorem jointPaths_share_left_components (t : Fin 32) :
    (jointPath t).1 = (jointPath247 t).1 :=
  ((jointPath247_marginals t).1).symm

/-- The 86/88 left hinge count is the committed `count86`: the 86/88
joint path's left components are definitionally `obs86.path`. -/
theorem leftHingeCount88_eq_count86 (i : Fin 13) :
    leftHingeCount88 i = count86 i := rfl

/-- The 86/247 left hinge count is the committed `count86`. -/
theorem leftHingeCount247_eq_count86 (i : Fin 13) :
    leftHingeCount247 i = count86 i := by
  unfold leftHingeCount247 count86
  apply congrArg Finset.card
  ext t
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  rw [(jointPath247_marginals t).1]

/-- The 86/88 hinge state is the committed left marginal state. -/
theorem leftHingeState88_eq_marginal86State :
    leftHingeState88 = marginal86State := rfl

/-- The two hinge states are equal. -/
theorem leftHingeState247_eq_leftHingeState88 :
    leftHingeState247 = leftHingeState88 := by
  unfold leftHingeState247 leftHingeState88
  have h : (fun i : Fin 13 => (leftHingeCount247 i : ℂ) / 32) =
      fun i : Fin 13 => (leftHingeCount88 i : ℂ) / 32 := by
    funext i
    rw [leftHingeCount247_eq_count86 i, leftHingeCount88_eq_count86 i]
  rw [h]

/-- **The shared-observer hinge agreement.**  Observer 86 is the left
observer of both committed pairs, and both committed joint paths
marginalize the same committed `obs86` payload.  Exactly: the two joint
paths have equal left label components at every step; the induced left
occupation count vectors are equal; the induced diagonal hinge states
are equal and equal the committed `marginal86State`; and the left
partial trace of the committed counted correlation state of the 86/247
pair is exactly the hinge state read from the 86/88 side.  This is the
first committed theorem relating the two carriers.  It is an exact
left-marginal agreement of the shared observer's occupation data, not
an identification of the carriers, nets, or towers. -/
theorem towerCarriers_share_obs86_marginal :
    (∀ t : Fin 32, (jointPath t).1 = (jointPath247 t).1) ∧
    (∀ i : Fin 13, leftHingeCount88 i = leftHingeCount247 i) ∧
    leftHingeState88 = leftHingeState247 ∧
    leftHingeState88 = marginal86State ∧
    ptraceSnd correlationState = leftHingeState88 := by
  refine ⟨jointPaths_share_left_components, ?_,
    leftHingeState247_eq_leftHingeState88.symm,
    leftHingeState88_eq_marginal86State, ?_⟩
  · intro i
    rw [leftHingeCount88_eq_count86 i, leftHingeCount247_eq_count86 i]
  · rw [correlationState_marginal_left]
    exact leftHingeState88_eq_marginal86State.symm

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
-- Expected axioms per declaration: propext, Classical.choice, Quot.sound.
#print axioms OPH.QFT.pairFin247
#print axioms OPH.QFT.ambientEquiv247
#print axioms OPH.QFT.ambientEquiv247_apply
#print axioms OPH.QFT.anchoredCheckpointPartition247
#print axioms OPH.QFT.anchoredState247
#print axioms OPH.QFT.anchoredTower247
#print axioms OPH.QFT.anchoredTower247_privateAlgebra
#print axioms OPH.QFT.anchoredNet247
#print axioms OPH.QFT.partition_members_in_left_region247
#print axioms OPH.QFT.anchoredCheckpointPartition247_proper
#print axioms OPH.QFT.anchoredNet247_left_ne_top
#print axioms OPH.QFT.anchored_partitions_share_chk86
#print axioms OPH.QFT.bothCarriers_common_anchoring
#print axioms OPH.QFT.leftHingeCount88
#print axioms OPH.QFT.leftHingeCount247
#print axioms OPH.QFT.leftHingeState88
#print axioms OPH.QFT.leftHingeState247
#print axioms OPH.QFT.jointPaths_share_left_components
#print axioms OPH.QFT.leftHingeCount88_eq_count86
#print axioms OPH.QFT.leftHingeCount247_eq_count86
#print axioms OPH.QFT.leftHingeState88_eq_marginal86State
#print axioms OPH.QFT.leftHingeState247_eq_leftHingeState88
#print axioms OPH.QFT.towerCarriers_share_obs86_marginal
