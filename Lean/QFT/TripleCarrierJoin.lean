import QFT.PathTimeSliceInterface
import QFT.TowerAnchoredCorrelation

/-!
# A common carrier for the committed 86/88 and 86/247 paths

The committed quantum surface previously related the `86/88` carrier
`Fin 13 × Fin 14` and the `86/247` carrier `Fin 13 × Fin 13` only through
their common observer-86 marginal.  This module constructs an actual common
finite carrier from the same source rows:

```
((86 label, 88 label), 247 label)
```

at each of the same 32 steps.  Its two coordinate projections are genuine
surjective, noninjective maps, not dimension-changing casts.
They recover both committed pair paths pointwise.  The source-counted diagonal
state on the triple carrier has exact partial traces equal to the counted
states on both pair carriers; on the 86/247 side this is the already committed
`correlationState`.

The observer-86 checkpoint data are compatible with the same join: after
normalized partial trace over the unused coordinate, each triple checkpoint
member recovers the corresponding un-reindexed checkpoint member on both
pair carriers, and hence the corresponding members of
`anchoredCheckpointPartition` and `anchoredCheckpointPartition247` after the
committed basis reindexings.

This does **not** construct a morphism of either anchored regional net or
constant consensus tower.  In particular, partial trace is a matrix-level
channel, not a unital star-algebra homomorphism, and the existing tower states
are independently declared uniform states.  Net/tower morphisms therefore
remain an explicit residual.  The two marginal maps are typed below as
linear, positive, and trace preserving.  No complete-positivity interface is
instantiated here, so this module must not be cited as a typed CPTP receipt.
The shared source-step alignment selects one concrete triple coupling; the
two pair marginals do not in general reconstruct an arbitrary triple matrix,
and no uniqueness of positive state couplings is claimed.

Finally, on each of the 31 actual adjacent source transitions, the product of
the three path-read transpositions transports the triple path and its diagonal
slice projector from `t.castSucc` to `t.succ`.  Its two marginal maps
intertwine this evolution with the corresponding pair evolutions.  The domain
is deliberately `Fin 31`: no final-to-first cyclic edge is invented.  This is
an evolved-data finite compatibility theorem.  It is not physical time, a
Cauchy surface, a source-selected dynamics, a continuum net, or a discharge
of PR-52/PR-58.  The slot reading and the common-step assembly remain the
declared finite postprocessor already present in the imported source modules.
-/

set_option autoImplicit false
set_option maxRecDepth 65536

namespace OPH.QFT

open Matrix
open Kronecker
open scoped ComplexOrder

noncomputable section

/-! ## The carrier and its two nontrivial quotient maps -/

/-- The common three-observer carrier, associated so that ordinary
`ptraceSnd` drops observer 247 and recovers the 86/88 carrier. -/
abbrev TripleCarrier := PairIndex × Fin 13

/-- Forget observer 247. -/
def tripleToPair88 (q : TripleCarrier) : PairIndex := q.1

/-- Forget observer 88 while retaining the common observer 86. -/
def tripleToPair247 (q : TripleCarrier) : PairIndex247 := (q.1.1, q.2)

/-- Every 86/88 label pair has a lift to the triple carrier. -/
theorem tripleToPair88_surjective : Function.Surjective tripleToPair88 := by
  intro p
  exact ⟨(p, 0), rfl⟩

/-- Every 86/247 label pair has a lift to the triple carrier. -/
theorem tripleToPair247_surjective : Function.Surjective tripleToPair247 := by
  intro p
  exact ⟨((p.1, 0), p.2), rfl⟩

/-- The 86/88 projection genuinely forgets information. -/
theorem tripleToPair88_not_injective : ¬ Function.Injective tripleToPair88 := by
  intro h
  have heq : (((0, 0), 0) : TripleCarrier) = (((0, 0), 1) : TripleCarrier) :=
    h rfl
  exact (by decide : (((0, 0), 0) : TripleCarrier) ≠ (((0, 0), 1) : TripleCarrier)) heq

/-- The 86/247 projection genuinely forgets information. -/
theorem tripleToPair247_not_injective : ¬ Function.Injective tripleToPair247 := by
  intro h
  have heq : (((0, 0), 0) : TripleCarrier) = (((0, 1), 0) : TripleCarrier) :=
    h rfl
  exact (by decide : (((0, 0), 0) : TripleCarrier) ≠ (((0, 1), 0) : TripleCarrier)) heq

/-! ## The same-step triple path -/

/-- The three committed observer labels at each of the same 32 source steps. -/
def triplePath (t : Fin 32) : TripleCarrier :=
  (jointPath t, obs247.path t)

/-- Projection to the 86/88 pair recovers its committed path pointwise. -/
theorem triplePath_projects_88 (t : Fin 32) :
    tripleToPair88 (triplePath t) = jointPath t := rfl

/-- Projection to the 86/247 pair recovers its committed path pointwise. -/
theorem triplePath_projects_247 (t : Fin 32) :
    tripleToPair247 (triplePath t) = jointPath247 t := by
  apply Prod.ext
  · exact (jointPath247_marginals t).1.symm
  · exact (jointPath247_marginals t).2.symm

/-! ## Source-counted laws and exact state marginals -/

/-- Occupation count on the common triple carrier. -/
def tripleCount (q : TripleCarrier) : ℕ :=
  ∑ t : Fin 32, if triplePath t = q then 1 else 0

/-- Occupation count on the committed 86/88 path. -/
def jointCount88 (p : PairIndex) : ℕ :=
  ∑ t : Fin 32, if jointPath t = p then 1 else 0

/-- The existing filter-card definition of the 86/247 table is the
equivalent indicator sum. -/
theorem jointCount_eq_indicator (p : PairIndex247) :
    jointCount p = ∑ t : Fin 32, if jointPath247 t = p then 1 else 0 := by
  unfold jointCount
  rw [Finset.card_filter]

/-- The triple occupation table contains exactly the 32 committed steps. -/
theorem tripleCount_total : ∑ q : TripleCarrier, tripleCount q = 32 := by
  unfold tripleCount
  rw [Finset.sum_comm]
  simp

/-- Marginalizing observer 247 recovers the exact 86/88 occupation table. -/
theorem tripleCount_marginal_88 (p : PairIndex) :
    ∑ k : Fin 13, tripleCount (p, k) = jointCount88 p := by
  unfold tripleCount jointCount88
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro t _
  by_cases h : jointPath t = p
  · have ht : triplePath t = (p, obs247.path t) := by
      rw [triplePath, h]
    simp [ht, h]
  · have ht : ∀ k : Fin 13, triplePath t ≠ (p, k) := by
      intro k hk
      exact h (congrArg Prod.fst hk)
    simp [h, ht]

/-- Marginalizing observer 88 recovers the exact committed 86/247 table. -/
theorem tripleCount_marginal_247 (p : PairIndex247) :
    ∑ j : Fin 14, tripleCount ((p.1, j), p.2) = jointCount p := by
  rw [jointCount_eq_indicator]
  unfold tripleCount
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro t _
  by_cases h : jointPath247 t = p
  · have htriple : triplePath t = ((p.1, obs88.path t), p.2) := by
      apply Prod.ext
      · apply Prod.ext
        · exact (jointPath247_marginals t).1.symm.trans (congrArg Prod.fst h)
        · rfl
      · exact (jointPath247_marginals t).2.symm.trans (congrArg Prod.snd h)
    simp only [htriple]
    simp [h]
  · have ht : ∀ j : Fin 14, triplePath t ≠ ((p.1, j), p.2) := by
      intro j hj
      apply h
      apply Prod.ext
      · exact (jointPath247_marginals t).1.trans
          (congrArg (fun q => q.1.1) hj)
      · exact (jointPath247_marginals t).2.trans
          (congrArg (fun q => q.2) hj)
    simp [h, ht]

/-- The source-counted diagonal state on the common carrier. -/
def tripleCorrelationState : Matrix TripleCarrier TripleCarrier ℂ :=
  Matrix.diagonal fun q => (tripleCount q : ℂ) / 32

/-- The source-counted diagonal state on the 86/88 carrier. -/
def correlationState88 : Matrix PairIndex PairIndex ℂ :=
  Matrix.diagonal fun p => (jointCount88 p : ℂ) / 32

/-- Partial trace over the middle (observer-88) coordinate. -/
def ptraceMiddle
    (M : Matrix TripleCarrier TripleCarrier ℂ) :
    Matrix PairIndex247 PairIndex247 ℂ :=
  fun p p' => ∑ j : Fin 14, M ((p.1, j), p.2) ((p'.1, j), p'.2)

/-- Reassociation/permutation of the trace indices used by
`ptraceMiddle`. -/
def middleTraceIndexEquiv : PairIndex247 × Fin 14 ≃ TripleCarrier where
  toFun q := ((q.1.1, q.2), q.1.2)
  invFun q := ((q.1.1, q.2), q.1.2)
  left_inv q := by cases q with | mk p j => cases p; rfl
  right_inv q := by cases q with | mk p k => cases p; rfl

/-- The middle partial trace preserves the ordinary matrix trace. -/
theorem ptraceMiddle_trace (M : Matrix TripleCarrier TripleCarrier ℂ) :
    (ptraceMiddle M).trace = M.trace := by
  change (∑ p : PairIndex247, ∑ j : Fin 14,
      M (middleTraceIndexEquiv (p, j)) (middleTraceIndexEquiv (p, j))) =
    ∑ q : TripleCarrier, M q q
  calc
    _ = ∑ q : PairIndex247 × Fin 14,
        M (middleTraceIndexEquiv q) (middleTraceIndexEquiv q) :=
      (Fintype.sum_prod_type
        (fun q : PairIndex247 × Fin 14 =>
          M (middleTraceIndexEquiv q) (middleTraceIndexEquiv q))).symm
    _ = _ := Equiv.sum_comp middleTraceIndexEquiv (fun q => M q q)

/-- Middle partial trace of a diagonal triple matrix is the diagonal of
the observer-88 fibre sums. -/
theorem ptraceMiddle_diagonal (f : TripleCarrier → ℂ) :
    ptraceMiddle (Matrix.diagonal f) =
      Matrix.diagonal fun p : PairIndex247 => ∑ j : Fin 14, f ((p.1, j), p.2) := by
  ext p p'
  unfold ptraceMiddle
  by_cases h : p = p'
  · subst p'
    simp
  · have hpoint : ∀ j : Fin 14,
        ((p.1, j), p.2) ≠ ((p'.1, j), p'.2) := by
      intro j hj
      apply h
      exact Prod.ext (congrArg (fun q => q.1.1) hj)
        (congrArg (fun q => q.2) hj)
    rw [Matrix.diagonal_apply_ne _ h]
    simp [Matrix.diagonal_apply_ne, hpoint]

/-- The 86/88 marginal as an explicit linear map. -/
def tripleMarginal88 :
    Matrix TripleCarrier TripleCarrier ℂ →ₗ[ℂ] Matrix PairIndex PairIndex ℂ where
  toFun := ptraceSnd
  map_add' M N := by
    ext p p'
    simp [ptraceSnd, Finset.sum_add_distrib]
  map_smul' c M := by
    ext p p'
    simp [ptraceSnd, Finset.mul_sum]

/-- The 86/247 marginal as an explicit linear map. -/
def tripleMarginal247 :
    Matrix TripleCarrier TripleCarrier ℂ →ₗ[ℂ]
      Matrix PairIndex247 PairIndex247 ℂ where
  toFun := ptraceMiddle
  map_add' M N := by
    ext p p'
    simp [ptraceMiddle, Finset.sum_add_distrib]
  map_smul' c M := by
    ext p p'
    simp [ptraceMiddle, Finset.mul_sum]

theorem tripleMarginal88_apply (M : Matrix TripleCarrier TripleCarrier ℂ) :
    tripleMarginal88 M = ptraceSnd M := rfl

theorem tripleMarginal247_apply (M : Matrix TripleCarrier TripleCarrier ℂ) :
    tripleMarginal247 M = ptraceMiddle M := rfl

/-- The 86/88 marginal preserves positive semidefiniteness.  The proof
exhibits it as a sum of principal compressions. -/
theorem tripleMarginal88_posSemidef
    {M : Matrix TripleCarrier TripleCarrier ℂ} (hM : M.PosSemidef) :
    (tripleMarginal88 M).PosSemidef := by
  have heq : tripleMarginal88 M =
      ∑ k : Fin 13, M.submatrix (fun p : PairIndex => (p, k))
        (fun p : PairIndex => (p, k)) := by
    ext p p'
    rw [Matrix.sum_apply]
    rfl
  rw [heq]
  exact Matrix.posSemidef_sum _ fun k _ => hM.submatrix _

/-- The 86/247 marginal preserves positive semidefiniteness, again as a
sum of principal compressions. -/
theorem tripleMarginal247_posSemidef
    {M : Matrix TripleCarrier TripleCarrier ℂ} (hM : M.PosSemidef) :
    (tripleMarginal247 M).PosSemidef := by
  have heq : tripleMarginal247 M =
      ∑ j : Fin 14, M.submatrix
        (fun p : PairIndex247 => ((p.1, j), p.2))
        (fun p : PairIndex247 => ((p.1, j), p.2)) := by
    ext p p'
    rw [Matrix.sum_apply]
    rfl
  rw [heq]
  exact Matrix.posSemidef_sum _ fun j _ => hM.submatrix _

/-- Both explicit marginal maps preserve trace. -/
theorem tripleMarginals_trace_preserving
    (M : Matrix TripleCarrier TripleCarrier ℂ) :
    (tripleMarginal88 M).trace = M.trace ∧
      (tripleMarginal247 M).trace = M.trace :=
  ⟨ptraceSnd_trace M, ptraceMiddle_trace M⟩

/-- A nonzero off-diagonal triple matrix erased by both pair marginals. -/
def tripleMarginalKernelWitness : Matrix TripleCarrier TripleCarrier ℂ :=
  Matrix.single (((0, 0), 0) : TripleCarrier)
    (((1, 1), 1) : TripleCarrier) 1

theorem tripleMarginalKernelWitness_ne_zero :
    tripleMarginalKernelWitness ≠ 0 := by
  intro h
  have hentry := congrFun
    (congrFun h (((0, 0), 0) : TripleCarrier)) (((1, 1), 1) : TripleCarrier)
  simp [tripleMarginalKernelWitness] at hentry

theorem tripleMarginal88_kernelWitness :
    tripleMarginal88 tripleMarginalKernelWitness = 0 := by
  rw [tripleMarginal88_apply]
  ext p p'
  change (∑ k : Fin 13,
    tripleMarginalKernelWitness (p, k) (p', k)) = 0
  apply Finset.sum_eq_zero
  intro k _
  by_cases hk : k = 0
  · subst k
    simp [tripleMarginalKernelWitness, Matrix.single]
  · have hrow : (((0, 0), 0) : TripleCarrier) ≠ (p, k) := by
      intro h
      exact hk (congrArg Prod.snd h).symm
    simp [tripleMarginalKernelWitness, Matrix.single, hrow]

theorem tripleMarginal247_kernelWitness :
    tripleMarginal247 tripleMarginalKernelWitness = 0 := by
  rw [tripleMarginal247_apply]
  ext p p'
  change (∑ j : Fin 14,
    tripleMarginalKernelWitness ((p.1, j), p.2) ((p'.1, j), p'.2)) = 0
  apply Finset.sum_eq_zero
  intro j _
  by_cases hj : j = 0
  · subst j
    simp [tripleMarginalKernelWitness, Matrix.single]
  · have hrow : (((0, 0), 0) : TripleCarrier) ≠ ((p.1, j), p.2) := by
      intro h
      exact hj (congrArg (fun z => z.1.2) h).symm
    simp [tripleMarginalKernelWitness, Matrix.single, hrow]

/-- The two pair marginals do not jointly reconstruct arbitrary triple
matrices.  This is a delimitation of the join, not a construction of two
distinct positive state couplings. -/
theorem tripleMarginals_not_jointly_injective :
    ¬ Function.Injective (fun M : Matrix TripleCarrier TripleCarrier ℂ =>
      (tripleMarginal88 M, tripleMarginal247 M)) := by
  intro h
  have hpair :
      (tripleMarginal88 tripleMarginalKernelWitness,
        tripleMarginal247 tripleMarginalKernelWitness) =
      (tripleMarginal88 0, tripleMarginal247 0) := by
    rw [tripleMarginal88_kernelWitness, tripleMarginal247_kernelWitness]
    simp
  exact tripleMarginalKernelWitness_ne_zero (h hpair)

/-- The triple state is positive semidefinite. -/
theorem tripleCorrelationState_posSemidef : tripleCorrelationState.PosSemidef := by
  refine Matrix.posSemidef_diagonal_iff.mpr fun q => ?_
  have h : (0 : ℝ) ≤ (tripleCount q : ℝ) / 32 := by positivity
  rw [show ((tripleCount q : ℂ) / 32) =
      (((tripleCount q : ℝ) / 32 : ℝ) : ℂ) by push_cast; ring]
  exact_mod_cast h

/-- The triple state has trace one. -/
theorem tripleCorrelationState_trace : tripleCorrelationState.trace = 1 := by
  unfold tripleCorrelationState
  rw [Matrix.trace_diagonal, ← Finset.sum_div]
  rw [show (∑ q : TripleCarrier, (tripleCount q : ℂ)) =
      ((∑ q : TripleCarrier, tripleCount q : ℕ) : ℂ) by push_cast; rfl]
  rw [tripleCount_total]
  norm_num

/-- Dropping observer 247 recovers the exact counted 86/88 state. -/
theorem tripleCorrelationState_marginal_88 :
    ptraceSnd tripleCorrelationState = correlationState88 := by
  unfold tripleCorrelationState correlationState88
  rw [ptraceSnd_diagonal]
  congr 1
  funext p
  rw [← Finset.sum_div]
  congr 1
  exact_mod_cast tripleCount_marginal_88 p

/-- Dropping observer 88 recovers the already committed counted 86/247
correlation state exactly. -/
theorem tripleCorrelationState_marginal_247 :
    ptraceMiddle tripleCorrelationState = correlationState := by
  unfold tripleCorrelationState correlationState
  rw [ptraceMiddle_diagonal]
  congr 1
  funext p
  rw [← Finset.sum_div]
  congr 1
  exact_mod_cast tripleCount_marginal_247 p

/-- The recovered 86/88 counted matrix is a positive semidefinite state. -/
theorem correlationState88_posSemidef : correlationState88.PosSemidef := by
  rw [← tripleCorrelationState_marginal_88]
  exact tripleMarginal88_posSemidef tripleCorrelationState_posSemidef

/-- The recovered 86/88 counted state has trace one. -/
theorem correlationState88_trace : correlationState88.trace = 1 := by
  rw [← tripleCorrelationState_marginal_88, ptraceSnd_trace,
    tripleCorrelationState_trace]

/-! ## Checkpoint compatibility (not a net/tower morphism) -/

/-- The observer-86 checkpoint member lifted to both auxiliary slots. -/
def tripleCheckpoint (v : Fin 4) : Matrix TripleCarrier TripleCarrier ℂ :=
  (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ)) ⊗ₖ
    (1 : Matrix (Fin 13) (Fin 13) ℂ)

/-- Partial trace over the middle factor of a three-factor Kronecker
product.  This is a matrix-level marginal identity, not an algebra
homomorphism. -/
theorem ptraceMiddle_triple_kronecker
    (A : Matrix (Fin 13) (Fin 13) ℂ)
    (B : Matrix (Fin 14) (Fin 14) ℂ)
    (C : Matrix (Fin 13) (Fin 13) ℂ) :
    ptraceMiddle ((A ⊗ₖ B) ⊗ₖ C) = B.trace • (A ⊗ₖ C) := by
  ext p p'
  simp only [ptraceMiddle, Matrix.kroneckerMap_apply, Matrix.smul_apply,
    smul_eq_mul]
  unfold Matrix.trace
  simp only [Matrix.diag]
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- Normalized trace over observer 247 recovers the 86/88 checkpoint
member before the committed basis reindexing. -/
theorem tripleCheckpoint_marginal_88 (v : Fin 4) :
    (13 : ℂ)⁻¹ • ptraceSnd (tripleCheckpoint v) =
      chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ) := by
  unfold tripleCheckpoint
  rw [ptraceSnd_slotLeft, smul_smul]
  norm_num

/-- Normalized trace over observer 88 recovers the 86/247 checkpoint
member before the committed basis reindexing. -/
theorem tripleCheckpoint_marginal_247 (v : Fin 4) :
    (14 : ℂ)⁻¹ • ptraceMiddle (tripleCheckpoint v) =
      chk86 v ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ) := by
  unfold tripleCheckpoint
  rw [ptraceMiddle_triple_kronecker, Matrix.trace_one, smul_smul]
  norm_num

/-- The 86/88 anchored checkpoint is the reindexing of the first exact
triple marginal. -/
theorem tripleCheckpoint_recovers_anchored_88 (v : Fin 4) :
    ambientEquiv ((13 : ℂ)⁻¹ • ptraceSnd (tripleCheckpoint v)) =
      anchoredCheckpointPartition.proj v := by
  rw [tripleCheckpoint_marginal_88]
  rfl

/-- The 86/247 anchored checkpoint is the reindexing of the second exact
triple marginal. -/
theorem tripleCheckpoint_recovers_anchored_247 (v : Fin 4) :
    ambientEquiv247 ((14 : ℂ)⁻¹ • ptraceMiddle (tripleCheckpoint v)) =
      anchoredCheckpointPartition247.proj v := by
  rw [tripleCheckpoint_marginal_247]
  rfl

/-! ## One-step evolved-data compatibility on the common carrier -/

/-- Every one of the 31 adjacent row pairs is present in each observer's
committed transition-count table. -/
theorem source_adjacent_transition_counts_nonzero :
    (∀ t : Fin 31,
      obs86.cnt (obs86.path t.castSucc) (obs86.path t.succ) ≠ 0) ∧
    (∀ t : Fin 31,
      obs88.cnt (obs88.path t.castSucc) (obs88.path t.succ) ≠ 0) ∧
    (∀ t : Fin 31,
      obs247.cnt (obs247.path t.castSucc) (obs247.path t.succ) ≠ 0) := by
  decide

/-- In contrast, the final-to-initial pair is absent from all three
committed transition-count tables.  This is why the evolution below has
domain `Fin 31` rather than a cyclic `Fin 32`. -/
theorem source_terminal_to_initial_counts_zero :
    obs86.cnt (obs86.path 31) (obs86.path 0) = 0 ∧
    obs88.cnt (obs88.path 31) (obs88.path 0) = 0 ∧
    obs247.cnt (obs247.path 31) (obs247.path 0) = 0 := by
  decide

/-- Observer 88's path-read transposition at one of the 31 actual
adjacent source transitions. -/
def walkSwap88 (t : Fin 31) : Fin 14 ≃ Fin 14 :=
  Equiv.swap (obs88.path t.castSucc) (obs88.path (t.castSucc + 1))

/-- The corresponding pair evolution on the committed 86/88 carrier. -/
def walkStepEquiv88 (t : Fin 31) : PairIndex ≃ PairIndex :=
  (walkSwapLeft t.castSucc).prodCongr (walkSwap88 t)

/-- The product of all three path-read transpositions. -/
def tripleStepEquiv (t : Fin 31) : TripleCarrier ≃ TripleCarrier :=
  (walkStepEquiv88 t).prodCongr (walkSwapRight t.castSucc)

/-- The common-carrier one-step star automorphism. -/
def tripleStepEvolve (t : Fin 31) :
    Matrix TripleCarrier TripleCarrier ℂ ≃⋆ₐ[ℂ]
      Matrix TripleCarrier TripleCarrier ℂ :=
  StarAlgEquiv.ofAlgEquiv (Matrix.reindexAlgEquiv ℂ ℂ (tripleStepEquiv t))
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose (tripleStepEquiv t) M)

/-- The pair-86/88 one-step star automorphism. -/
def pair88StepEvolve (t : Fin 31) :
    Matrix PairIndex PairIndex ℂ ≃⋆ₐ[ℂ] Matrix PairIndex PairIndex ℂ :=
  StarAlgEquiv.ofAlgEquiv (Matrix.reindexAlgEquiv ℂ ℂ (walkStepEquiv88 t))
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose (walkStepEquiv88 t) M)

/-- Casting one of the 31 transition indices into the 32 source rows and
adding one reaches its non-wrapping successor. -/
theorem castSucc_add_one_eq_succ (t : Fin 31) :
    t.castSucc + 1 = t.succ := by
  apply Fin.ext
  simp [Fin.add_def]
  omega

/-- The triple evolution carries the common path across one actual
adjacent source transition. -/
theorem tripleStepEquiv_carries_path (t : Fin 31) :
    tripleStepEquiv t (triplePath t.castSucc) = triplePath t.succ := by
  unfold tripleStepEquiv walkStepEquiv88 triplePath jointPath walkSwap88
    walkSwapLeft walkSwapRight
  simp only [Equiv.prodCongr_apply]
  change
    ((Equiv.swap (obs86.path t.castSucc) (obs86.path (t.castSucc + 1))
          (obs86.path t.castSucc),
        Equiv.swap (obs88.path t.castSucc) (obs88.path (t.castSucc + 1))
          (obs88.path t.castSucc)),
      Equiv.swap (obs247.path t.castSucc) (obs247.path (t.castSucc + 1))
        (obs247.path t.castSucc)) =
    ((obs86.path t.succ, obs88.path t.succ), obs247.path t.succ)
  rw [castSucc_add_one_eq_succ]
  rw [Equiv.swap_apply_left, Equiv.swap_apply_left, Equiv.swap_apply_left]

/-- The diagonal slice projector on the common path. -/
def tripleWalkSliceProjector (t : Fin 32) :
    Matrix TripleCarrier TripleCarrier ℂ :=
  Matrix.single (triplePath t) (triplePath t) 1

/-- The corresponding diagonal slice projector on the 86/88 path. -/
def pair88WalkSliceProjector (t : Fin 32) : Matrix PairIndex PairIndex ℂ :=
  Matrix.single (jointPath t) (jointPath t) 1

/-- Partial trace of one diagonal matrix unit over the last coordinate. -/
theorem ptraceSnd_single_diagonal (q : TripleCarrier) :
    ptraceSnd (Matrix.single q q (1 : ℂ)) =
      Matrix.single q.1 q.1 1 := by
  ext p p'
  unfold ptraceSnd
  rw [Finset.sum_eq_single q.2]
  · rcases q with ⟨q1, q2⟩
    simp [Matrix.single]
  · intro k _ hk
    have hq : q ≠ (p, k) := by
      intro h
      exact hk (congrArg Prod.snd h).symm
    simp [Matrix.single, hq]
  · intro h
    exact absurd (Finset.mem_univ q.2) h

/-- Partial trace of one diagonal matrix unit over the middle coordinate. -/
theorem ptraceMiddle_single_diagonal (q : TripleCarrier) :
    ptraceMiddle (Matrix.single q q (1 : ℂ)) =
      Matrix.single (tripleToPair247 q) (tripleToPair247 q) 1 := by
  ext p p'
  unfold ptraceMiddle tripleToPair247
  rw [Finset.sum_eq_single q.1.2]
  · rcases q with ⟨⟨qi, qj⟩, qk⟩
    rcases p with ⟨pi, pk⟩
    rcases p' with ⟨pi', pk'⟩
    simp [Matrix.single]
  · intro j _ hj
    have hq : q ≠ ((p.1, j), p.2) := by
      intro h
      exact hj (congrArg (fun z => z.1.2) h).symm
    simp [Matrix.single, hq]
  · intro h
    exact absurd (Finset.mem_univ q.1.2) h

/-- The 86/88 marginal of the common slice is exactly the pair slice. -/
theorem tripleWalkSliceProjector_marginal_88 (t : Fin 32) :
    tripleMarginal88 (tripleWalkSliceProjector t) =
      pair88WalkSliceProjector t := by
  unfold tripleMarginal88 tripleWalkSliceProjector pair88WalkSliceProjector
  change ptraceSnd (Matrix.single (triplePath t) (triplePath t) 1) = _
  rw [ptraceSnd_single_diagonal]
  rfl

/-- The 86/247 marginal of the common slice is exactly the already
committed pair slice. -/
theorem tripleWalkSliceProjector_marginal_247 (t : Fin 32) :
    tripleMarginal247 (tripleWalkSliceProjector t) =
      walkSliceProjector t := by
  unfold tripleMarginal247 tripleWalkSliceProjector walkSliceProjector
  change ptraceMiddle (Matrix.single (triplePath t) (triplePath t) 1) = _
  rw [ptraceMiddle_single_diagonal, triplePath_projects_247]

/-- The common evolution transports the common slice projector using
evolved, rather than raw, data. -/
theorem tripleStepEvolve_transports_slice (t : Fin 31) :
    tripleStepEvolve t (tripleWalkSliceProjector t.castSucc) =
      tripleWalkSliceProjector t.succ := by
  unfold tripleStepEvolve tripleWalkSliceProjector
  simp only [StarAlgEquiv.ofAlgEquiv_apply, Matrix.reindexAlgEquiv_apply]
  rw [reindex_single, tripleStepEquiv_carries_path]

/-- Partial trace to the 86/88 carrier intertwines the common evolution
with the pair's path-read evolution. -/
theorem ptraceSnd_tripleStepEvolve (t : Fin 31)
    (M : Matrix TripleCarrier TripleCarrier ℂ) :
    ptraceSnd (tripleStepEvolve t M) = pair88StepEvolve t (ptraceSnd M) := by
  ext p p'
  unfold ptraceSnd tripleStepEvolve pair88StepEvolve
  simp only [StarAlgEquiv.ofAlgEquiv_apply, Matrix.reindexAlgEquiv_apply,
    Matrix.reindex_apply, Matrix.submatrix_apply, tripleStepEquiv]
  rw [← Equiv.sum_comp (walkSwapRight t.castSucc)]
  apply Finset.sum_congr rfl
  intro i _
  simp [Equiv.prodCongr_symm]

/-- Partial trace to the 86/247 carrier intertwines the common evolution
with the already committed 86/247 path-read evolution. -/
theorem ptraceMiddle_tripleStepEvolve (t : Fin 31)
    (M : Matrix TripleCarrier TripleCarrier ℂ) :
    ptraceMiddle (tripleStepEvolve t M) =
      stepEvolve t.castSucc (ptraceMiddle M) := by
  ext p p'
  unfold ptraceMiddle tripleStepEvolve stepEvolve
  simp only [StarAlgEquiv.ofAlgEquiv_apply, Matrix.reindexAlgEquiv_apply,
    Matrix.reindex_apply, Matrix.submatrix_apply, tripleStepEquiv,
    walkStepEquiv88, walkStepEquiv]
  rw [← Equiv.sum_comp (walkSwap88 t)]
  apply Finset.sum_congr rfl
  intro j _
  simp [Equiv.prodCongr_symm]

/-- The common-carrier **path/state/checkpoint** join in one receipt:
both path projections, both state marginals, both checkpoint marginals,
slice transport, and evolution/marginal intertwining.  Every statement
is finite and exact.  No anchored-net or tower-state morphism is included. -/
theorem tripleCarrier_join_receipt :
    (Function.Surjective tripleToPair88 ∧
      ¬ Function.Injective tripleToPair88) ∧
    (Function.Surjective tripleToPair247 ∧
      ¬ Function.Injective tripleToPair247) ∧
    (∀ t : Fin 32, tripleToPair88 (triplePath t) = jointPath t) ∧
    (∀ t : Fin 32, tripleToPair247 (triplePath t) = jointPath247 t) ∧
    ptraceSnd tripleCorrelationState = correlationState88 ∧
    ptraceMiddle tripleCorrelationState = correlationState ∧
    (∀ v : Fin 4,
      ambientEquiv ((13 : ℂ)⁻¹ • ptraceSnd (tripleCheckpoint v)) =
        anchoredCheckpointPartition.proj v) ∧
    (∀ v : Fin 4,
      ambientEquiv247 ((14 : ℂ)⁻¹ • ptraceMiddle (tripleCheckpoint v)) =
        anchoredCheckpointPartition247.proj v) ∧
    (∀ t : Fin 31,
      tripleStepEvolve t (tripleWalkSliceProjector t.castSucc) =
        tripleWalkSliceProjector t.succ) ∧
    (∀ t : Fin 32,
      tripleMarginal88 (tripleWalkSliceProjector t) =
        pair88WalkSliceProjector t) ∧
    (∀ t : Fin 32,
      tripleMarginal247 (tripleWalkSliceProjector t) =
        walkSliceProjector t) ∧
    (∀ (t : Fin 31) (M : Matrix TripleCarrier TripleCarrier ℂ),
      ptraceSnd (tripleStepEvolve t M) = pair88StepEvolve t (ptraceSnd M)) ∧
    (∀ (t : Fin 31) (M : Matrix TripleCarrier TripleCarrier ℂ),
      ptraceMiddle (tripleStepEvolve t M) =
        stepEvolve t.castSucc (ptraceMiddle M)) :=
  ⟨⟨tripleToPair88_surjective, tripleToPair88_not_injective⟩,
    ⟨tripleToPair247_surjective, tripleToPair247_not_injective⟩,
    triplePath_projects_88, triplePath_projects_247,
    tripleCorrelationState_marginal_88, tripleCorrelationState_marginal_247,
    tripleCheckpoint_recovers_anchored_88,
    tripleCheckpoint_recovers_anchored_247,
    tripleStepEvolve_transports_slice,
    tripleWalkSliceProjector_marginal_88,
    tripleWalkSliceProjector_marginal_247,
    ptraceSnd_tripleStepEvolve,
    ptraceMiddle_tripleStepEvolve⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
-- Expected axioms per declaration: propext, Classical.choice, Quot.sound.
#print axioms OPH.QFT.tripleToPair88_surjective
#print axioms OPH.QFT.tripleToPair247_surjective
#print axioms OPH.QFT.tripleToPair88_not_injective
#print axioms OPH.QFT.tripleToPair247_not_injective
#print axioms OPH.QFT.triplePath_projects_88
#print axioms OPH.QFT.triplePath_projects_247
#print axioms OPH.QFT.tripleCount_marginal_88
#print axioms OPH.QFT.tripleCount_marginal_247
#print axioms OPH.QFT.ptraceMiddle_trace
#print axioms OPH.QFT.tripleMarginal88
#print axioms OPH.QFT.tripleMarginal247
#print axioms OPH.QFT.tripleMarginal88_posSemidef
#print axioms OPH.QFT.tripleMarginal247_posSemidef
#print axioms OPH.QFT.tripleMarginals_trace_preserving
#print axioms OPH.QFT.tripleMarginalKernelWitness_ne_zero
#print axioms OPH.QFT.tripleMarginal88_kernelWitness
#print axioms OPH.QFT.tripleMarginal247_kernelWitness
#print axioms OPH.QFT.tripleMarginals_not_jointly_injective
#print axioms OPH.QFT.tripleCorrelationState_posSemidef
#print axioms OPH.QFT.tripleCorrelationState_trace
#print axioms OPH.QFT.tripleCorrelationState_marginal_88
#print axioms OPH.QFT.tripleCorrelationState_marginal_247
#print axioms OPH.QFT.correlationState88_posSemidef
#print axioms OPH.QFT.correlationState88_trace
#print axioms OPH.QFT.tripleCheckpoint_marginal_88
#print axioms OPH.QFT.tripleCheckpoint_marginal_247
#print axioms OPH.QFT.tripleCheckpoint_recovers_anchored_88
#print axioms OPH.QFT.tripleCheckpoint_recovers_anchored_247
#print axioms OPH.QFT.source_adjacent_transition_counts_nonzero
#print axioms OPH.QFT.source_terminal_to_initial_counts_zero
#print axioms OPH.QFT.tripleStepEquiv_carries_path
#print axioms OPH.QFT.tripleStepEvolve_transports_slice
#print axioms OPH.QFT.tripleWalkSliceProjector_marginal_88
#print axioms OPH.QFT.tripleWalkSliceProjector_marginal_247
#print axioms OPH.QFT.ptraceSnd_tripleStepEvolve
#print axioms OPH.QFT.ptraceMiddle_tripleStepEvolve
#print axioms OPH.QFT.tripleCarrier_join_receipt
