import QFT.ConvexAffinityBridge

/-!
# A source-produced convex preparation family from the B12 run

The `ConvexAffinityBridge` module derives outcome affinity for every
preparation family with a declared convex mixing law; the open B13 core
is source production of such a family.  This module mirrors, as exact
literals, a four-fibre subfamily of the convex family carried by the B12
preregistered source run, together with the boundary theorem naming what
this realization cannot supply.

## Source provenance

Run `b12_prereg_16k_20260806` of the sim repository (seed `20260806`,
run git commit `b39b78faf894894ebe573571e0902ccfaaeac32a`) commits the
conditional-resampling realization receipt
`conditional_resampling_realization_receipt.json` (sha256
`d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749`),
whose pinned reference is the realized joint frequency table over 32
record classes and 8 realized companion classes on 16384 patches, built
from `freezeout_fields.npz` (sha256
`b962c5b80205a17d5d6bc023f5f5d487bc23c1327793978e7d1bc69494fac49e`).
The extraction script `scripts/extract_b13_convex_family.py` re-derives
the table with the receipt's own binning producer and emits
`docs/B13_CONVEX_FAMILY_PAYLOAD.json`: the 32 exact fibre-conditional
laws, the realized record marginal, and declared record-rebalancing
mixtures, each verified in exact rational arithmetic to equal the
companion marginal of the correspondingly reweighted table.  The fibre
laws are realized states on the companion alphabet, and convex mixtures
of fibre laws weighted by any record distribution are exactly the
companion marginals of record-rebalanced versions of the run's own
table, so the mixing operation is realized by the run's record
statistics rather than declared by hand.

The four laws mirrored here are the fibre-conditional laws of record
classes 0, 1, 2, 3 (the payload's declared selection rule: the first
four record classes in label order), each a vector of integer counts
over the fibre mass 512.  The realized companion alphabet has the eight
labels `0, 1, 3, 6, 7, 10, 13, 15`, mirrored as `Fin 8` in label order.
Each fibre carries mass 512 of the restricted total 2048, so the
realized record marginal restricted to the chosen fibres is exactly
`(1/4, 1/4, 1/4, 1/4)`.

## Content

* the four fibre laws as exact rational literals, with probability
  receipts and exact pairwise distinctness;
* their embedding as diagonal density matrices on the companion
  alphabet, certified by `EventAlgebra.IsState`;
* `sourceConvexFamily`, a `ConvexPreparationFamily` over the simplex of
  record weights on the four fibres, whose mixing operation is
  record-marginal rebalancing and whose mixing law is proved exactly
  (`sourceMixture_eq_diagonal`: a mixture of diagonals is the diagonal
  of the mixture);
* the affinity receipt `sourceConvexFamily_outcome_affine`, consumed
  from the committed `ConvexPreparationFamily.outcome_affine`;
* the earned-family theorem `sourceConvexFamily_earned_from_run`: the
  simplex vertices realize the run's fibre laws, the payload's two
  committed simplex points (the restricted realized record marginal and
  the declared rebalanced point) are family states whose diagonals are
  the payload's exact mixed laws, and the four fibre states are pairwise
  distinct;
* the boundary theorems: every family outcome factors through the
  diagonal compression of the effect
  (`sourceConvexFamily_outcome_eq_diagonalPart`), the realized effect
  algebra is commutative (`sourceRealizedEffects_commute`), and no
  noncommuting context pair exists in this realization
  (`source_no_noncommuting_context_pair`).

## Claim boundary

The literals mirror realized frequencies of one committed run; the
family is a mathematical mirror of source data, and its mixing operation
is the run's record-marginal rebalancing.  What the full Born bridge
additionally needs is interlocking contexts: shared effects across
incompatible measurements.  The boundary theorems prove that this
diagonal realization cannot supply them, because every effect it can
present commutes with every other and every outcome reads only the
diagonal of the probed effect.  No physical instrument, measurement
attachment, or Born-rule closure is claimed.
-/

namespace OPH.QFT

open EventAlgebra
open Matrix
open scoped ComplexOrder

noncomputable section

/-! ## The mirrored fibre laws -/

/-- The fibre-conditional laws of record classes 0, 1, 2, 3 of the B12
run, as exact rationals over the fibre mass 512.  Companion coordinates
are indexed in realized-label order `0, 1, 3, 6, 7, 10, 13, 15`.  The
integer numerators are the realized cell counts of the pinned joint
frequency table. -/
def sourceFibreLaw : Fin 4 → Fin 8 → ℚ :=
  ![![17/512, 40/512, 49/512, 92/512, 17/512, 138/512, 99/512, 60/512],
    ![12/512, 34/512, 82/512, 102/512, 21/512, 105/512, 102/512, 54/512],
    ![19/512, 25/512, 71/512, 87/512, 22/512, 139/512, 90/512, 59/512],
    ![17/512, 37/512, 83/512, 89/512, 12/512, 109/512, 98/512, 67/512]]

/-- Every mirrored law entry is nonnegative. -/
theorem sourceFibreLaw_nonneg : ∀ (i : Fin 4) (j : Fin 8),
    0 ≤ sourceFibreLaw i j := by
  intro i j
  fin_cases i <;> fin_cases j <;> norm_num [sourceFibreLaw]

/-- Every mirrored law is normalized: its entries sum to one. -/
theorem sourceFibreLaw_sum : ∀ i : Fin 4, ∑ j, sourceFibreLaw i j = 1 := by
  intro i
  fin_cases i <;> norm_num [sourceFibreLaw, Fin.sum_univ_succ]

/-- The four mirrored laws are pairwise distinct, exactly.  Companion
coordinate 1 separates every pair: its four values `40/512`, `34/512`,
`25/512`, `37/512` are pairwise distinct. -/
theorem sourceFibreLaw_pairwise_ne : ∀ i j : Fin 4, i ≠ j →
    sourceFibreLaw i ≠ sourceFibreLaw j := by
  intro i j hij h
  have h1 := congrFun h 1
  fin_cases i <;> fin_cases j <;>
    first
      | exact absurd rfl hij
      | norm_num [sourceFibreLaw] at h1

/-! ## Diagonal embedding as certified states -/

/-- A diagonal matrix with nonnegative real entries of unit sum is a
certified state: positivity is entrywise and the trace is the entry
sum. -/
theorem isState_diagonal_prob (v : Fin 8 → ℝ) (h0 : ∀ j, 0 ≤ v j)
    (h1 : ∑ j, v j = 1) :
    IsState (Matrix.diagonal fun j => ((v j : ℝ) : ℂ)) := by
  constructor
  · refine Matrix.PosSemidef.diagonal ?_
    intro j
    exact Complex.zero_le_real.mpr (h0 j)
  · rw [Matrix.trace_diagonal, ← Complex.ofReal_sum, h1, Complex.ofReal_one]

/-- The diagonal density matrix of fibre `i`: the exact law embedded on
the diagonal of the companion alphabet. -/
def sourceFibreState (i : Fin 4) : Matrix (Fin 8) (Fin 8) ℂ :=
  Matrix.diagonal fun j => ((sourceFibreLaw i j : ℝ) : ℂ)

/-- Every fibre state is a certified state. -/
theorem sourceFibreState_isState (i : Fin 4) : IsState (sourceFibreState i) :=
  isState_diagonal_prob _ (fun j => by exact_mod_cast sourceFibreLaw_nonneg i j)
    (by exact_mod_cast sourceFibreLaw_sum i)

/-- The four fibre states are pairwise distinct: equality of the
diagonal matrices forces equality of the laws entry by entry. -/
theorem sourceFibreState_pairwise_ne : ∀ i j : Fin 4, i ≠ j →
    sourceFibreState i ≠ sourceFibreState j := by
  intro i j hij h
  refine sourceFibreLaw_pairwise_ne i j hij ?_
  funext k
  have hk := congrFun (congrFun h k) k
  simp only [sourceFibreState, Matrix.diagonal_apply_eq] at hk
  exact_mod_cast hk

/-! ## The record-rebalancing simplex and the source family -/

/-- The simplex of record weights over the four chosen fibres: the
parameter space of record-marginal rebalancing. -/
def SourceRecordSimplex : Type :=
  {w : Fin 4 → ℝ // (∀ i, 0 ≤ w i) ∧ ∑ i, w i = 1}

/-- The mixing law in closed form: a convex mixture of the diagonal
fibre states is the diagonal matrix of the mixed law. -/
theorem sourceMixture_eq_diagonal (w : Fin 4 → ℝ) :
    ∑ i, (w i : ℂ) • sourceFibreState i
      = Matrix.diagonal
          fun j => ((∑ i, w i * (sourceFibreLaw i j : ℝ) : ℝ) : ℂ) := by
  ext j k
  by_cases hjk : j = k
  · subst hjk
    simp only [Matrix.sum_apply, Matrix.smul_apply, sourceFibreState,
      Matrix.diagonal_apply_eq, smul_eq_mul]
    push_cast
    rfl
  · simp only [Matrix.sum_apply, Matrix.smul_apply, sourceFibreState,
      Matrix.diagonal_apply_ne _ hjk, smul_zero, Finset.sum_const_zero]

/-- The source convex family: states are record-weight mixtures of the
four diagonal fibre states, and the mixing operation is record-marginal
rebalancing on the weight simplex.  Both the states and the mixing
operation are mirrors of realized run data; the mixing law is proved
from the diagonal closed form. -/
def sourceConvexFamily : ConvexPreparationFamily 8 SourceRecordSimplex where
  mix l hl w w' :=
    ⟨fun i => l * w.1 i + (1 - l) * w'.1 i, by
      constructor
      · intro i
        have h1l : (0 : ℝ) ≤ 1 - l := by linarith [hl.2]
        exact add_nonneg (mul_nonneg hl.1 (w.2.1 i))
          (mul_nonneg h1l (w'.2.1 i))
      · rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
          w.2.2, w'.2.2]
        ring⟩
  state w := ∑ i, (w.1 i : ℂ) • sourceFibreState i
  state_isState w := by
    rw [sourceMixture_eq_diagonal w.1]
    refine isState_diagonal_prob _ (fun j => ?_) ?_
    · exact Finset.sum_nonneg fun i _ =>
        mul_nonneg (w.2.1 i) (by exact_mod_cast sourceFibreLaw_nonneg i j)
    · rw [Finset.sum_comm]
      have hrow : ∀ i : Fin 4, ∑ j, w.1 i * (sourceFibreLaw i j : ℝ)
          = w.1 i := by
        intro i
        rw [← Finset.mul_sum,
          show (∑ j, (sourceFibreLaw i j : ℝ)) = 1 by
            exact_mod_cast sourceFibreLaw_sum i,
          mul_one]
      rw [Finset.sum_congr rfl fun i _ => hrow i]
      exact w.2.2
  state_mix l hl w w' := by
    simp only [Finset.smul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    push_cast
    module

/-! ## Realization receipts: vertices and committed simplex points -/

/-- The vertex of the simplex concentrated on fibre `i`. -/
def sourceVertex (i : Fin 4) : SourceRecordSimplex :=
  ⟨fun k => if k = i then 1 else 0, by
    refine ⟨fun k => ?_, by simp⟩
    by_cases hk : k = i <;> simp [hk]⟩

/-- The family state at vertex `i` is exactly the run's fibre state. -/
theorem sourceConvexFamily_state_vertex (i : Fin 4) :
    sourceConvexFamily.state (sourceVertex i) = sourceFibreState i := by
  simp [sourceConvexFamily, sourceVertex, ite_smul]

/-- The simplex point of an exact rational weight vector. -/
def sourceRatPoint (v : Fin 4 → ℚ) (h0 : ∀ i, 0 ≤ v i)
    (h1 : ∑ i, v i = 1) : SourceRecordSimplex :=
  ⟨fun i => ((v i : ℚ) : ℝ), by
    refine ⟨fun i => ?_, ?_⟩
    · show (0 : ℝ) ≤ ((v i : ℚ) : ℝ)
      exact_mod_cast h0 i
    · show (∑ i, ((v i : ℚ) : ℝ)) = 1
      exact_mod_cast h1⟩

/-- The family state at a rational simplex point is the diagonal matrix
of the exact rational mixed law. -/
theorem sourceConvexFamily_state_ratPoint (v : Fin 4 → ℚ)
    (h0 : ∀ i, 0 ≤ v i) (h1 : ∑ i, v i = 1) :
    sourceConvexFamily.state (sourceRatPoint v h0 h1)
      = Matrix.diagonal
          fun j => (((∑ i, v i * sourceFibreLaw i j : ℚ) : ℝ) : ℂ) := by
  show (∑ i, (((v i : ℚ) : ℝ) : ℂ) • sourceFibreState i) = _
  rw [sourceMixture_eq_diagonal]
  congr 1
  funext j
  push_cast
  rfl

/-- The restricted realized record marginal: each chosen fibre carries
mass 512 of the restricted total 2048, so the renormalized weights are
exactly `(1/4, 1/4, 1/4, 1/4)`. -/
def sourceUniformWeights : Fin 4 → ℚ := ![1/4, 1/4, 1/4, 1/4]

/-- The declared rebalanced simplex point of the payload. -/
def sourceRebalancedWeights : Fin 4 → ℚ := ![1/2, 1/4, 1/8, 1/8]

theorem sourceUniformWeights_nonneg : ∀ i, 0 ≤ sourceUniformWeights i := by
  intro i
  fin_cases i <;> norm_num [sourceUniformWeights]

theorem sourceUniformWeights_sum : ∑ i, sourceUniformWeights i = 1 := by
  norm_num [sourceUniformWeights, Fin.sum_univ_succ]

theorem sourceRebalancedWeights_nonneg :
    ∀ i, 0 ≤ sourceRebalancedWeights i := by
  intro i
  fin_cases i <;> norm_num [sourceRebalancedWeights]

theorem sourceRebalancedWeights_sum :
    ∑ i, sourceRebalancedWeights i = 1 := by
  norm_num [sourceRebalancedWeights, Fin.sum_univ_succ]

/-- The payload's exact mixed companion law at the restricted realized
record marginal. -/
def sourceUniformMixedLaw : Fin 8 → ℚ :=
  ![65/2048, 17/256, 285/2048, 185/1024, 9/256, 491/2048, 389/2048,
    15/128]

/-- The payload's exact mixed companion law at the declared rebalanced
point. -/
def sourceRebalancedMixedLaw : Fin 8 → ℚ :=
  ![1/32, 145/2048, 257/2048, 187/1024, 9/256, 505/2048, 197/1024,
    237/2048]

/-- Exact identity between the convex mixture at the restricted realized
marginal and the payload's mixed law. -/
theorem sourceUniformMixedLaw_exact :
    ∀ j, ∑ i, sourceUniformWeights i * sourceFibreLaw i j
      = sourceUniformMixedLaw j := by
  intro j
  fin_cases j <;>
    norm_num [sourceUniformWeights, sourceFibreLaw, sourceUniformMixedLaw,
      Fin.sum_univ_succ]

/-- Exact identity between the convex mixture at the declared rebalanced
point and the payload's mixed law. -/
theorem sourceRebalancedMixedLaw_exact :
    ∀ j, ∑ i, sourceRebalancedWeights i * sourceFibreLaw i j
      = sourceRebalancedMixedLaw j := by
  intro j
  fin_cases j <;>
    norm_num [sourceRebalancedWeights, sourceFibreLaw,
      sourceRebalancedMixedLaw, Fin.sum_univ_succ]

/-- The restricted realized record marginal as a simplex point. -/
def sourceRestrictedMarginalPoint : SourceRecordSimplex :=
  sourceRatPoint sourceUniformWeights sourceUniformWeights_nonneg
    sourceUniformWeights_sum

/-- The declared rebalanced simplex point. -/
def sourceRebalancedPoint : SourceRecordSimplex :=
  sourceRatPoint sourceRebalancedWeights sourceRebalancedWeights_nonneg
    sourceRebalancedWeights_sum

/-- The family state at the restricted realized record marginal is the
diagonal of the payload's exact mixed law. -/
theorem sourceConvexFamily_state_restrictedMarginal :
    sourceConvexFamily.state sourceRestrictedMarginalPoint
      = Matrix.diagonal
          fun j => ((sourceUniformMixedLaw j : ℝ) : ℂ) := by
  unfold sourceRestrictedMarginalPoint
  rw [sourceConvexFamily_state_ratPoint]
  congr 1
  funext j
  rw [sourceUniformMixedLaw_exact j]

/-- The family state at the declared rebalanced point is the diagonal of
the payload's exact mixed law. -/
theorem sourceConvexFamily_state_rebalanced :
    sourceConvexFamily.state sourceRebalancedPoint
      = Matrix.diagonal
          fun j => ((sourceRebalancedMixedLaw j : ℝ) : ℂ) := by
  unfold sourceRebalancedPoint
  rw [sourceConvexFamily_state_ratPoint]
  congr 1
  funext j
  rw [sourceRebalancedMixedLaw_exact j]

/-! ## The affinity receipt -/

/-- The affinity receipt for the source family, consumed from the
committed bridge theorem: the outcome functional of every effect matrix
is affine in the record-rebalancing mixing parameter. -/
theorem sourceConvexFamily_outcome_affine (E : Matrix (Fin 8) (Fin 8) ℂ)
    (l : ℝ) (hl : l ∈ Set.Icc (0 : ℝ) 1) (p q : SourceRecordSimplex) :
    bornWeight (sourceConvexFamily.state (sourceConvexFamily.mix l hl p q)) E
      = (l : ℂ) * bornWeight (sourceConvexFamily.state p) E
        + ((1 - l : ℝ) : ℂ) * bornWeight (sourceConvexFamily.state q) E :=
  sourceConvexFamily.outcome_affine E l hl p q

/-- The earned-family theorem.  The family is realized from the run's
fibre laws with the mixing operation realized by record-marginal
rebalancing: the simplex vertices are the four fibre states extracted
from the pinned joint table, the restricted realized record marginal and
the declared rebalanced point are family states whose diagonals are the
payload's exact mixed laws, and the four fibre states are pairwise
distinct.  Provenance: run `b12_prereg_16k_20260806`, receipt sha256
`d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749`,
payload `docs/B13_CONVEX_FAMILY_PAYLOAD.json` of the sim repository. -/
theorem sourceConvexFamily_earned_from_run :
    (∀ i, sourceConvexFamily.state (sourceVertex i) = sourceFibreState i)
    ∧ sourceConvexFamily.state sourceRestrictedMarginalPoint
        = Matrix.diagonal (fun j => ((sourceUniformMixedLaw j : ℝ) : ℂ))
    ∧ sourceConvexFamily.state sourceRebalancedPoint
        = Matrix.diagonal (fun j => ((sourceRebalancedMixedLaw j : ℝ) : ℂ))
    ∧ (∀ i j : Fin 4, i ≠ j → sourceFibreState i ≠ sourceFibreState j) :=
  ⟨sourceConvexFamily_state_vertex,
    sourceConvexFamily_state_restrictedMarginal,
    sourceConvexFamily_state_rebalanced,
    sourceFibreState_pairwise_ne⟩

/-! ## The commutative-context boundary -/

/-- The effect algebra this realization can present: the diagonal
matrices on the companion alphabet. -/
def sourceRealizedEffectAlgebra : Set (Matrix (Fin 8) (Fin 8) ℂ) :=
  Set.range Matrix.diagonal

/-- The diagonal compression of a matrix. -/
def sourceDiagonalPart (E : Matrix (Fin 8) (Fin 8) ℂ) :
    Matrix (Fin 8) (Fin 8) ℂ :=
  Matrix.diagonal fun j => E j j

/-- The diagonal compression lands in the realized effect algebra. -/
theorem sourceDiagonalPart_mem (E : Matrix (Fin 8) (Fin 8) ℂ) :
    sourceDiagonalPart E ∈ sourceRealizedEffectAlgebra :=
  ⟨_, rfl⟩

/-- Every family outcome factors through the diagonal compression: the
family's states are diagonal, so the trace pairing reads only the
diagonal of the probed effect.  The family therefore cannot distinguish
an effect from its diagonal part. -/
theorem sourceConvexFamily_outcome_eq_diagonalPart
    (p : SourceRecordSimplex) (E : Matrix (Fin 8) (Fin 8) ℂ) :
    bornWeight (sourceConvexFamily.state p) E
      = bornWeight (sourceConvexFamily.state p) (sourceDiagonalPart E) := by
  have hp : sourceConvexFamily.state p
      = Matrix.diagonal
          fun j => ((∑ i, p.1 i * (sourceFibreLaw i j : ℝ) : ℝ) : ℂ) :=
    sourceMixture_eq_diagonal p.1
  rw [hp]
  simp [bornWeight, Matrix.trace, Matrix.diag, Matrix.diagonal_mul,
    sourceDiagonalPart]

/-- The realized effect algebra is commutative: any two diagonal effects
commute. -/
theorem sourceRealizedEffects_commute :
    ∀ E ∈ sourceRealizedEffectAlgebra, ∀ F ∈ sourceRealizedEffectAlgebra,
      E * F = F * E := by
  rintro E ⟨e, rfl⟩ F ⟨f, rfl⟩
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  exact congrArg _ (funext fun j => mul_comm (e j) (f j))

/-- The family's states commute pairwise: all of them are diagonal. -/
theorem sourceConvexFamily_states_commute (p q : SourceRecordSimplex) :
    sourceConvexFamily.state p * sourceConvexFamily.state q
      = sourceConvexFamily.state q * sourceConvexFamily.state p := by
  show (∑ i, (p.1 i : ℂ) • sourceFibreState i)
      * (∑ i, (q.1 i : ℂ) • sourceFibreState i)
    = (∑ i, (q.1 i : ℂ) • sourceFibreState i)
      * (∑ i, (p.1 i : ℂ) • sourceFibreState i)
  rw [sourceMixture_eq_diagonal p.1, sourceMixture_eq_diagonal q.1,
    Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  exact congrArg _ (funext fun j => mul_comm _ _)

/-- The boundary theorem: no noncommuting context pair exists in this
realization.  Interlocking contexts require shared effects across
incompatible measurements; the realized effect algebra is commutative,
so this diagonal family supplies none.  This names the exact gap between
the earned convex-mixing production and the full B13 Born bridge. -/
theorem source_no_noncommuting_context_pair :
    ¬ ∃ E F : Matrix (Fin 8) (Fin 8) ℂ,
        E ∈ sourceRealizedEffectAlgebra ∧ F ∈ sourceRealizedEffectAlgebra
          ∧ E * F ≠ F * E := by
  rintro ⟨E, F, hE, hF, hne⟩
  exact hne (sourceRealizedEffects_commute E hE F hF)

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.
#print axioms sourceFibreLaw_nonneg
#print axioms sourceFibreLaw_sum
#print axioms sourceFibreLaw_pairwise_ne
#print axioms isState_diagonal_prob
#print axioms sourceFibreState_isState
#print axioms sourceFibreState_pairwise_ne
#print axioms sourceMixture_eq_diagonal
#print axioms sourceConvexFamily
#print axioms sourceConvexFamily_state_vertex
#print axioms sourceConvexFamily_state_ratPoint
#print axioms sourceUniformMixedLaw_exact
#print axioms sourceRebalancedMixedLaw_exact
#print axioms sourceConvexFamily_state_restrictedMarginal
#print axioms sourceConvexFamily_state_rebalanced
#print axioms sourceConvexFamily_outcome_affine
#print axioms sourceConvexFamily_earned_from_run
#print axioms sourceDiagonalPart_mem
#print axioms sourceConvexFamily_outcome_eq_diagonalPart
#print axioms sourceRealizedEffects_commute
#print axioms sourceConvexFamily_states_commute
#print axioms source_no_noncommuting_context_pair

end OPH.QFT
