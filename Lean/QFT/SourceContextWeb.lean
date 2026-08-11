import Mathlib

/-!
# Source context web: a declared algebraic adapter over B12 gauge labels

This module consumes the context-web payload
`docs/BORN_CONTEXT_WEB_PAYLOAD.json` of the `oph-physics-sim`
repository (schema `oph.sim.born_context_web_payload.v1`, produced by
`scripts/extract_born_context_web.py`) and mirrors, as exact literals
over `ℝ` with `Real.sqrt 3`, a source-attached algebraic projector web on
the designated split fibre of observer 64: the declared two-dimensional
irreducible S3 representation images of the run's realized edge gauge
elements, the conjugated binary contexts they algebraically generate over the diagonal
record/companion projectors of the fibre block, the exact
noncommutation certificates, and the outcome-frequency boundary that
separates the one context with realized run statistics from the
contexts without any.

## Provenance (verbatim from the payload's provenance block)

* `run_id`: `b12_prereg_16k_20260806`, `seed`: `20260806`,
  `run_git_commit`: `b39b78faf894894ebe573571e0902ccfaaeac32a`
* `receipt_sha256`:
  `d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749`
* `freezeout_sha256`:
  `b962c5b80205a17d5d6bc023f5f5d487bc23c1327793978e7d1bc69494fac49e`
* `gauge_state_sha256`:
  `5f2cd276b84b1917c4d321a532db6efd12442931fbc80fd764cbcf1918dadca2`
* `e1_payload_sha256`:
  `005223dc4fa7442bf10e6e8b446a1e5ebba08e41796d6ac89c45251de090bf45`
* extraction: `record_field = record_signature`,
  `companion_field = cumulative_repair_load`, gauge element tables from
  `oph_fpe.finite_groups`, incidence rule: an edge row counts once when
  at least one endpoint lies in the node set.

## The realized data

The run's final edge gauge state assigns one S3 element to each of the
97404 undirected edges.  On the 71 edges incident to observer 64's
truncated support the realized element counts in the index order of
`oph_fpe.finite_groups.S3_ELEMENTS` are `16, 12, 12, 11, 10, 10`; on
the 21 edges incident to the designated fibre nodes 64 and 43 they are
`8, 4, 3, 2, 3, 1`.  Every S3 element is realized on the support and on
the fibre nodes, the two permutation-order-3 elements included.  The
committed joint frequency table realizes the fibre block's two diagonal
outcomes `(record 18, companion 13)` and `(record 18, companion 3)`
with counts 111 and 68 out of the fibre mass 512 on 16384 patches.

## The irreducible representation

The exact two-dimensional irreducible representation fixes the
generator assignment `rho((1,2,0)) = rotation by 120 degrees` and
`rho((0,2,1)) = diag(1, -1)` under the composition convention
`(g*h)[i] = g[h[i]]` of `oph_fpe.finite_groups.S3_MUL`.  The extraction
script verifies the homomorphism property on all 36 ordered pairs, the
orthogonality of every image, and the trace-character match in exact
fraction arithmetic over Q(sqrt 3); this module mirrors the six image
matrices with entries of the form `a + b * sqrt3` for rational `a, b`
and proves orthogonality, the order-3 identities, and the projector
computations exactly from `sqrt3 ^ 2 = 3`.

## Receipts

* `gaugeIrrep_unitary` and the legacy-named `realized_images_unitary`:
  every declared representation image whose label occurs in the source
  census is exactly orthogonal, so unitary as a real matrix.
* `rotationPlus_order_three` with `s3Perm_threeCycles_order_three`: the
  permutation-order-3 elements land on matrices of multiplicative
  order 3, the 120-degree rotations.
* `conjProjector_three` / `conjProjector_four` and companions: the
  exact conjugated projector candidates for every occurring nonidentity label.
* `commutator_three` / `commutator_four` and
  `conjProjector_three_noncommuting` /
  `conjProjector_four_noncommuting`: the exact commutators of the
  order-3 conjugated projectors against the record projector, with the
  nonvanishing certificates.
* `realizedWeb_contains_noncommuting_context` (legacy API name): an
  occurring label has a declared conjugated context that fails to
  commute with the diagonal context.
* `conjProjector_one_eq_record`, `conjProjector_two_eq_three`,
  `conjProjector_four_eq_five`: the exact coincidence structure of the
  web; conjugation by the diagonal reflection reproduces the diagonal
  context, and the five nonidentity elements produce exactly two
  contexts beyond it.
* `conjContext_complete` and `conjProjector_idem`: each conjugated pair
  is a genuine binary context of idempotent effects summing to one.
* `realized_statistics_only_diagonal`,
  `additivityReceipt_requires_rotated_mass`, and
  `additivityReceipt_rotated_underdetermined`: the run's realized
  outcome statistics inhabit only the diagonal context;
  `AdditivityReceiptData` states the exact data requirement of the
  noncontextual additivity receipt over the web, every such datum must
  carry positive outcome mass on every rotated context, and the run
  data leaves those fields undetermined.

## Claim boundary

This module attaches a declared algebraic projector web to run data: the
source realizes gauge labels and counts, while the two-dimensional irrep,
projectors, and conjugation rule are mathematical adapters.  They give exact
representation images and genuinely noncommuting candidate contexts over the
declared fibre block, not source-produced instruments.  The noncontextual
additivity receipt over the web is the remaining Born obligation, and
it needs realized per-context outcome frequencies on the rotated
contexts.  Producing those statistics is a named simulator capability
gap: the simulator has no producer that realizes a gauge-conjugated
binary readout on the committed fibre states.  The gap must never be
filled by computing Born predictions from the diagonal state; the
payload records the same boundary in its machine-readable
`outcome_frequency_boundary` field.  No probability claim beyond the native
diagonal frequencies, no rotated measurement attachment, and no Born-rule
closure is asserted.
-/

namespace OPH.QFT

open Matrix

noncomputable section

/-! ## The exact scalar `sqrt 3` -/

/-- The square root of three, the one irrational scalar of the web. -/
def sqrt3 : ℝ := Real.sqrt 3

theorem sqrt3_sq : sqrt3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)

theorem sqrt3_cube : sqrt3 ^ 3 = 3 * sqrt3 := by
  have h : sqrt3 ^ 3 = sqrt3 ^ 2 * sqrt3 := by ring
  rw [h, sqrt3_sq]

theorem sqrt3_pos : 0 < sqrt3 := Real.sqrt_pos.mpr (by norm_num)

theorem sqrt3_ne_zero : sqrt3 ≠ 0 := ne_of_gt sqrt3_pos

/-! ## The permutation literals and the mirrored counts -/

/-- The six S3 elements in the index order of
`oph_fpe.finite_groups.S3_ELEMENTS`: the tuple at index `g` sends
`i` to `s3Perm g i`. -/
def s3Perm : Fin 6 → Fin 3 → Fin 3 :=
  ![![0, 1, 2], ![0, 2, 1], ![1, 0, 2], ![1, 2, 0], ![2, 0, 1], ![2, 1, 0]]

/-- Both permutation-order-3 elements cube to the identity. -/
theorem s3Perm_threeCycles_order_three :
    (∀ i, s3Perm 3 (s3Perm 3 (s3Perm 3 i)) = i)
      ∧ (∀ i, s3Perm 4 (s3Perm 4 (s3Perm 4 i)) = i) := by
  constructor <;> decide

/-- Realized element counts on the 71 edges incident to observer 64's
truncated support, payload field
`edge_gauge_extraction.per_observer[2].element_counts`. -/
def s64SupportGaugeCounts : Fin 6 → ℕ := ![16, 12, 12, 11, 10, 10]

/-- Realized element counts on the 21 edges incident to the designated
fibre nodes 64 and 43, payload field
`edge_gauge_extraction.per_split_fibre[0].element_counts`. -/
def designatedFibreGaugeCounts : Fin 6 → ℕ := ![8, 4, 3, 2, 3, 1]

theorem s64SupportGaugeCounts_total : ∑ g, s64SupportGaugeCounts g = 71 := by
  decide

theorem designatedFibreGaugeCounts_total :
    ∑ g, designatedFibreGaugeCounts g = 21 := by
  decide

/-- Every S3 element is realized on observer 64's support edges. -/
theorem gaugeElements_all_realized : ∀ g, s64SupportGaugeCounts g ≠ 0 := by
  decide

/-- Every S3 element is realized on the designated fibre edges. -/
theorem gaugeElements_all_realized_on_fibre :
    ∀ g, designatedFibreGaugeCounts g ≠ 0 := by
  decide

/-! ## The irreducible representation images -/

/-- Image of `(0,2,1)`, the transposition fixing the record axis:
`diag(1, -1)`. -/
def reflectionOne : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, -1]

/-- Image of `(1,0,2)`: the reflection across the 60-degree axis. -/
def reflectionTwo : Matrix (Fin 2) (Fin 2) ℝ :=
  !![-(1/2), sqrt3/2; sqrt3/2, 1/2]

/-- Image of `(1,2,0)`: the rotation by 120 degrees. -/
def rotationPlus : Matrix (Fin 2) (Fin 2) ℝ :=
  !![-(1/2), -(sqrt3/2); sqrt3/2, -(1/2)]

/-- Image of `(2,0,1)`: the rotation by 240 degrees. -/
def rotationMinus : Matrix (Fin 2) (Fin 2) ℝ :=
  !![-(1/2), sqrt3/2; -(sqrt3/2), -(1/2)]

/-- Image of `(2,1,0)`: the reflection across the 150-degree axis. -/
def reflectionFive : Matrix (Fin 2) (Fin 2) ℝ :=
  !![-(1/2), -(sqrt3/2); -(sqrt3/2), 1/2]

/-- The declared two-dimensional irreducible representation on the
source label indices, mirroring the payload's `irrep.elements` matrices. -/
def gaugeIrrep : Fin 6 → Matrix (Fin 2) (Fin 2) ℝ :=
  ![1, reflectionOne, reflectionTwo, rotationPlus, rotationMinus,
    reflectionFive]

theorem reflectionOne_orthogonal :
    reflectionOne * reflectionOneᵀ = 1 ∧ reflectionOneᵀ * reflectionOne = 1 := by
  constructor <;>
    · ext i j
      fin_cases i <;> fin_cases j <;>
        simp [reflectionOne, Matrix.mul_apply, Fin.sum_univ_two]

theorem reflectionTwo_orthogonal :
    reflectionTwo * reflectionTwoᵀ = 1 ∧ reflectionTwoᵀ * reflectionTwo = 1 := by
  have h := sqrt3_sq
  constructor <;>
    · ext i j
      fin_cases i <;> fin_cases j <;>
        simp [reflectionTwo, Matrix.mul_apply, Fin.sum_univ_two] <;>
        first | ring1 | linear_combination (1/4 : ℝ) * h

theorem rotationPlus_orthogonal :
    rotationPlus * rotationPlusᵀ = 1 ∧ rotationPlusᵀ * rotationPlus = 1 := by
  have h := sqrt3_sq
  constructor <;>
    · ext i j
      fin_cases i <;> fin_cases j <;>
        simp [rotationPlus, Matrix.mul_apply, Fin.sum_univ_two] <;>
        first | ring1 | linear_combination (1/4 : ℝ) * h

theorem rotationMinus_orthogonal :
    rotationMinus * rotationMinusᵀ = 1 ∧ rotationMinusᵀ * rotationMinus = 1 := by
  have h := sqrt3_sq
  constructor <;>
    · ext i j
      fin_cases i <;> fin_cases j <;>
        simp [rotationMinus, Matrix.mul_apply, Fin.sum_univ_two] <;>
        first | ring1 | linear_combination (1/4 : ℝ) * h

theorem reflectionFive_orthogonal :
    reflectionFive * reflectionFiveᵀ = 1 ∧ reflectionFiveᵀ * reflectionFive = 1 := by
  have h := sqrt3_sq
  constructor <;>
    · ext i j
      fin_cases i <;> fin_cases j <;>
        simp [reflectionFive, Matrix.mul_apply, Fin.sum_univ_two] <;>
        first | ring1 | linear_combination (1/4 : ℝ) * h

/-- Every image of the representation is exactly orthogonal. -/
theorem gaugeIrrep_unitary (g : Fin 6) :
    gaugeIrrep g * (gaugeIrrep g)ᵀ = 1 ∧ (gaugeIrrep g)ᵀ * gaugeIrrep g = 1 := by
  fin_cases g
  · exact ⟨by simp [gaugeIrrep], by simp [gaugeIrrep]⟩
  · exact reflectionOne_orthogonal
  · exact reflectionTwo_orthogonal
  · exact rotationPlus_orthogonal
  · exact rotationMinus_orthogonal
  · exact reflectionFive_orthogonal

/-- Every declared representation image is unitary; the premise records
that its label occurs on observer 64's support edges. -/
theorem realized_images_unitary :
    ∀ g : Fin 6, s64SupportGaugeCounts g ≠ 0 →
      gaugeIrrep g * (gaugeIrrep g)ᵀ = 1
        ∧ (gaugeIrrep g)ᵀ * gaugeIrrep g = 1 :=
  fun g _ => gaugeIrrep_unitary g

/-- Squaring the 120-degree rotation gives the 240-degree rotation. -/
theorem rotationPlus_sq : rotationPlus * rotationPlus = rotationMinus := by
  have h := sqrt3_sq
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rotationPlus, rotationMinus, Matrix.mul_apply,
      Fin.sum_univ_two] <;>
    first | ring1 | linear_combination (-(1/4 : ℝ)) * h

theorem rotationPlus_transpose : rotationPlusᵀ = rotationMinus := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rotationPlus, rotationMinus, Matrix.transpose_apply]

/-- The image of a permutation-order-3 element has multiplicative
order 3: it is a genuine 120-degree rotation. -/
theorem rotationPlus_order_three :
    rotationPlus * rotationPlus * rotationPlus = 1 := by
  rw [rotationPlus_sq, ← rotationPlus_transpose]
  exact (gaugeIrrep_unitary 3).2

/-! ## The diagonal context of the designated fibre block -/

/-- Declared projector onto the first block basis state, labelled by the
realized class pair `(record 18, companion 13)` at node 64. -/
def recordProjector : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, 0]

/-- Declared projector onto the second block basis state, labelled by the
realized class pair `(record 18, companion 3)` at node 43. -/
def companionProjector : Matrix (Fin 2) (Fin 2) ℝ := !![0, 0; 0, 1]

theorem record_add_companion : recordProjector + companionProjector = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [recordProjector, companionProjector]

theorem recordProjector_idem :
    recordProjector * recordProjector = recordProjector := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [recordProjector, Matrix.mul_apply, Fin.sum_univ_two]

/-! ## The conjugated contexts -/

/-- The record-side projector of the context conjugated by the image of
gauge element `g`. -/
def conjProjector (g : Fin 6) : Matrix (Fin 2) (Fin 2) ℝ :=
  gaugeIrrep g * recordProjector * (gaugeIrrep g)ᵀ

/-- The companion-side projector of the conjugated context. -/
def conjCompanionProjector (g : Fin 6) : Matrix (Fin 2) (Fin 2) ℝ :=
  gaugeIrrep g * companionProjector * (gaugeIrrep g)ᵀ

/-- Each conjugated pair is a complete binary context: the two
projectors sum to one. -/
theorem conjContext_complete (g : Fin 6) :
    conjProjector g + conjCompanionProjector g = 1 := by
  have h := (gaugeIrrep_unitary g).1
  simp only [conjProjector, conjCompanionProjector]
  rw [← add_mul, ← mul_add, record_add_companion, mul_one]
  exact h

/-- Each conjugated record projector is idempotent. -/
theorem conjProjector_idem (g : Fin 6) :
    conjProjector g * conjProjector g = conjProjector g := by
  have h := (gaugeIrrep_unitary g).2
  simp only [conjProjector]
  calc gaugeIrrep g * recordProjector * (gaugeIrrep g)ᵀ
        * (gaugeIrrep g * recordProjector * (gaugeIrrep g)ᵀ)
      = gaugeIrrep g * recordProjector * ((gaugeIrrep g)ᵀ * gaugeIrrep g)
        * recordProjector * (gaugeIrrep g)ᵀ := by
        simp only [mul_assoc]
    _ = gaugeIrrep g * recordProjector * (gaugeIrrep g)ᵀ := by
        rw [h, mul_one,
          mul_assoc (gaugeIrrep g) recordProjector recordProjector,
          recordProjector_idem]

/-- Conjugation by the diagonal reflection reproduces the diagonal
record projector exactly: this declared algebraic context coincides with the
diagonal one and carries no separate effect data. -/
theorem conjProjector_one_eq_record : conjProjector 1 = recordProjector := by
  show reflectionOne * recordProjector * reflectionOneᵀ = recordProjector
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [reflectionOne, recordProjector, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.transpose_apply]

theorem conjProjector_two :
    conjProjector 2 = !![1/4, -(sqrt3/4); -(sqrt3/4), 3/4] := by
  have h := sqrt3_sq
  show reflectionTwo * recordProjector * reflectionTwoᵀ = _
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [reflectionTwo, recordProjector, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.transpose_apply] <;>
    first | ring1 | linear_combination (1/4 : ℝ) * h

theorem conjProjector_three :
    conjProjector 3 = !![1/4, -(sqrt3/4); -(sqrt3/4), 3/4] := by
  have h := sqrt3_sq
  show rotationPlus * recordProjector * rotationPlusᵀ = _
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rotationPlus, recordProjector, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.transpose_apply] <;>
    first | ring1 | linear_combination (1/4 : ℝ) * h

theorem conjProjector_four :
    conjProjector 4 = !![1/4, sqrt3/4; sqrt3/4, 3/4] := by
  have h := sqrt3_sq
  show rotationMinus * recordProjector * rotationMinusᵀ = _
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rotationMinus, recordProjector, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.transpose_apply] <;>
    first | ring1 | linear_combination (1/4 : ℝ) * h

theorem conjProjector_five :
    conjProjector 5 = !![1/4, sqrt3/4; sqrt3/4, 3/4] := by
  have h := sqrt3_sq
  show reflectionFive * recordProjector * reflectionFiveᵀ = _
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [reflectionFive, recordProjector, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.transpose_apply] <;>
    first | ring1 | linear_combination (1/4 : ℝ) * h

/-- The web's coincidence structure: elements 2 and 3 generate one
rotated context. -/
theorem conjProjector_two_eq_three : conjProjector 2 = conjProjector 3 := by
  rw [conjProjector_two, conjProjector_three]

/-- The web's coincidence structure: elements 4 and 5 generate the
other rotated context. -/
theorem conjProjector_four_eq_five : conjProjector 4 = conjProjector 5 := by
  rw [conjProjector_four, conjProjector_five]

/-! ## The noncommutation certificates -/

/-- Exact commutator of the order-3 conjugated projector against the
record projector, mirroring the payload's
`commutator_with_record_projector` for element 3. -/
theorem commutator_three :
    conjProjector 3 * recordProjector - recordProjector * conjProjector 3
      = !![0, sqrt3/4; -(sqrt3/4), 0] := by
  rw [conjProjector_three]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [recordProjector, Matrix.sub_apply]

/-- Exact commutator for the other order-3 element. -/
theorem commutator_four :
    conjProjector 4 * recordProjector - recordProjector * conjProjector 4
      = !![0, -(sqrt3/4); sqrt3/4, 0] := by
  rw [conjProjector_four]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [recordProjector, Matrix.sub_apply]

/-- The conjugated context of the first order-3 element fails to commute
with the diagonal context: the commutator entry `sqrt3/4` is nonzero. -/
theorem conjProjector_three_noncommuting :
    conjProjector 3 * recordProjector ≠ recordProjector * conjProjector 3 := by
  intro hc
  have h01 : (conjProjector 3 * recordProjector
      - recordProjector * conjProjector 3) 0 1 = 0 := by
    rw [hc, sub_self, Matrix.zero_apply]
  rw [commutator_three] at h01
  have hval : sqrt3 / 4 = 0 := by simpa using h01
  have := sqrt3_pos
  linarith

/-- The conjugated context of the second order-3 element fails to
commute with the diagonal context. -/
theorem conjProjector_four_noncommuting :
    conjProjector 4 * recordProjector ≠ recordProjector * conjProjector 4 := by
  intro hc
  have h01 : (conjProjector 4 * recordProjector
      - recordProjector * conjProjector 4) 0 1 = 0 := by
    rw [hc, sub_self, Matrix.zero_apply]
  rw [commutator_four] at h01
  have hval : -(sqrt3 / 4) = 0 := by simpa using h01
  have := sqrt3_pos
  linarith

/-- The declared source-attached web contains genuinely noncommuting
projector candidates: an S3 label occurring on the designated fibre
edges indexes a conjugation that fails to commute with the record projector. -/
theorem realizedWeb_contains_noncommuting_context :
    ∃ g : Fin 6, designatedFibreGaugeCounts g ≠ 0
      ∧ conjProjector g * recordProjector
        ≠ recordProjector * conjProjector g :=
  ⟨3, by decide, conjProjector_three_noncommuting⟩

/-! ## The outcome-frequency boundary -/

/-- The contexts of the declared source-attached algebraic web: the
diagonal record/companion context and one conjugated context per gauge
label index. -/
inductive WebContext : Type where
  | diagonal : WebContext
  | conjugated : Fin 6 → WebContext
  deriving DecidableEq

/-- The realized outcome counts of the run, per context.  The diagonal
context carries the committed joint-table cells `(record 18,
companion 13) = 111` and `(record 18, companion 3) = 68`; no conjugated
context carries realized statistics. -/
def realizedOutcomeCounts : WebContext → Option (ℕ × ℕ)
  | WebContext.diagonal => some (111, 68)
  | WebContext.conjugated _ => none

/-- The two diagonal outcome counts, payload field
`context_web.diagonal_context.realized_outcome_counts`. -/
def diagonalOutcomeCounts : Fin 2 → ℕ := ![111, 68]

/-- Fibre mass of record class 18 in the committed joint table. -/
def designatedFibreMass : ℕ := 512

theorem diagonalOutcomeCounts_mass :
    diagonalOutcomeCounts 0 + diagonalOutcomeCounts 1 = 179 := rfl

theorem diagonal_counts_within_fibre :
    diagonalOutcomeCounts 0 + diagonalOutcomeCounts 1
      ≤ designatedFibreMass := by
  decide

theorem realizedOutcomeCounts_diagonal :
    realizedOutcomeCounts WebContext.diagonal
      = some (diagonalOutcomeCounts 0, diagonalOutcomeCounts 1) := rfl

/-- The boundary theorem: the run's realized outcome statistics inhabit
only the diagonal context. -/
theorem realized_statistics_only_diagonal (c : WebContext) :
    (realizedOutcomeCounts c).isSome ↔ c = WebContext.diagonal := by
  cases c <;> simp [realizedOutcomeCounts]

/-- The exact data requirement of the noncontextual additivity receipt
over the declared algebraic web: per-context outcome counts with positive mass on
every context, the rotated contexts included, with the diagonal fields
pinned to the run's committed values.  The run inhabits the diagonal
fields alone (`realized_statistics_only_diagonal`); filling the rotated
fields is a named simulator capability gap, and computing them as Born
predictions from the diagonal state is prohibited. -/
structure AdditivityReceiptData where
  counts : WebContext → ℕ × ℕ
  mass : WebContext → ℕ
  mass_pos : ∀ c, 0 < mass c
  counts_sum : ∀ c, (counts c).1 + (counts c).2 = mass c
  diagonal_counts_run : counts WebContext.diagonal = (111, 68)
  diagonal_mass_run : mass WebContext.diagonal = 179

/-- Every additivity-receipt datum must carry realized outcome mass on
every rotated context: the requirement the run does not meet. -/
theorem additivityReceipt_requires_rotated_mass
    (d : AdditivityReceiptData) (g : Fin 6) :
    0 < d.mass (WebContext.conjugated g) :=
  d.mass_pos _

/-- The run data leave the rotated fields of the additivity receipt
undetermined: two receipt data agree with every realized statistic and
disagree on a rotated context.  The receipt is an open obligation, not
a consequence of this module's literals. -/
theorem additivityReceipt_rotated_underdetermined :
    ∃ d₁ d₂ : AdditivityReceiptData,
      d₁.counts (WebContext.conjugated 3)
        ≠ d₂.counts (WebContext.conjugated 3) := by
  refine ⟨⟨fun c => if c = WebContext.diagonal then (111, 68) else (1, 0),
      fun c => if c = WebContext.diagonal then 179 else 1,
      ?_, ?_, by simp, by simp⟩,
    ⟨fun c => if c = WebContext.diagonal then (111, 68) else (0, 1),
      fun c => if c = WebContext.diagonal then 179 else 1,
      ?_, ?_, by simp, by simp⟩, ?_⟩
  · intro c
    by_cases hc : c = WebContext.diagonal <;> simp [hc]
  · intro c
    by_cases hc : c = WebContext.diagonal <;> simp [hc]
  · intro c
    by_cases hc : c = WebContext.diagonal <;> simp [hc]
  · intro c
    by_cases hc : c = WebContext.diagonal <;> simp [hc]
  · simp

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.
#print axioms s3Perm_threeCycles_order_three
#print axioms gaugeElements_all_realized
#print axioms gaugeElements_all_realized_on_fibre
#print axioms gaugeIrrep_unitary
#print axioms realized_images_unitary
#print axioms rotationPlus_order_three
#print axioms conjContext_complete
#print axioms conjProjector_idem
#print axioms conjProjector_one_eq_record
#print axioms conjProjector_two_eq_three
#print axioms conjProjector_four_eq_five
#print axioms commutator_three
#print axioms commutator_four
#print axioms conjProjector_three_noncommuting
#print axioms conjProjector_four_noncommuting
#print axioms realizedWeb_contains_noncommuting_context
#print axioms realized_statistics_only_diagonal
#print axioms additivityReceipt_requires_rotated_mass
#print axioms additivityReceipt_rotated_underdetermined

end OPH.QFT
