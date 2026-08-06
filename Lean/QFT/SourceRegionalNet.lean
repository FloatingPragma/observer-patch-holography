import QFT.SimEarnedWitness
import QFT.CoverageReceipt

/-!
# Source regional net: the enriched earned net from the B12 split fibres

This module consumes the fibre payload
`docs/E1_REGIONAL_PAYLOAD.json` of the `oph-physics-sim` repository
(schema `oph.sim.e1_regional_payload.v1`, produced by
`scripts/extract_e1_regional_payload.py`) and builds, over the earned
tower `simEarnedTower` of `QFT/SimEarnedWitness.lean`, a second finite
causal observer net whose regional algebra at the split-fibre support is
genuinely noncommutative.  Every off-diagonal direction of that algebra
is indexed by a realized split fibre of the committed run; no synthetic
matrix unit enters.

## Provenance (verbatim from the payload's provenance block)

* `run_id`: `b12_prereg_16k_20260806`, `seed`: `20260806`,
  `run_git_commit`: `b39b78faf894894ebe573571e0902ccfaaeac32a`
* `receipt_sha256`:
  `d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749`
* `freezeout_sha256`:
  `b962c5b80205a17d5d6bc023f5f5d487bc23c1327793978e7d1bc69494fac49e`
* `observer_views_sha256`:
  `3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca`
* `witness_payload_sha256`:
  `e880f8d75f1dada71a05e47ded177897c93f3279147a165020fab127cf89352f`
* extraction: `record_field = record_signature`,
  `companion_field = cumulative_repair_load`, `support_truncation = 8`,
  `split_fibre_rule = a record class of a truncated support is a split
  fibre when it carries at least two distinct companion classes there`,
  `designated_block_rule = first split fibre in observer file order and
  truncated support position order`.

## The realized fibre structure at truncation 8

The payload reports, per observer, which record classes carry several
companion classes on the truncated support, with exact counts:

* observer 36: no split fibre (7 record classes, all single-companion;
  record class 20 appears twice with the one companion class 6);
* observer 60: no split fibre (8 pairwise distinct record classes);
* observer 64: two split fibres.  Record class 18 is realized at
  support positions 0 and 3 (nodes 64 and 43) with companion classes 13
  and 3; record class 28 is realized at support positions 4 and 7
  (nodes 51 and 77) with companion classes 13 and 6.  The payload's
  designated block is the fibre of record class 18, the first split
  fibre in observer file order and support position order;
* observer 92: no split fibre (8 pairwise distinct record classes).

In the `Fin 33` basis layout of `QFT/SimEarnedWitness.lean` the fibre
of record class 18 occupies the basis indices `{17, 20}` and the fibre
of record class 28 occupies `{21, 24}`; the receipts
`splitFibre18_basis` and `splitFibre28_basis` pin this translation.

## The enrichment rule and the collapse receipts

The declared enrichment rule is uniform in the payload literals: two
basis indices of one observer support are related exactly when they lie
in one realized split fibre of that observer (`splitPair`), and the
regional algebra of the support is the block algebra of that relation:
matrices supported on related pairs whose diagonal is constant across
the ambient indices.  For observer 64 the rule yields two free
two-by-two blocks, one per realized split fibre, and the noncommutative
receipts below are earned from the run.  For observers 36, 60, and 92
the rule yields the relation `i = j`, and the collapse receipts
(`enrichedRegion_collapse_o36` and companions) state that the enriched
algebra equals the diagonal regional algebra of the sim-earned net:
the truncated payload data cannot support a noncommutative receipt
there, and the collapse is proved instead of being repaired by
synthetic data.

The net `sourceRegionalNet` decorates `simEarnedTower` with these
regional algebras: scalars at the bottom, the unchanged diagonal spans
at the three collapsed supports, the enriched block algebra at `s64`,
and the block-diagonal enriched algebra at the top.  Restriction,
isotony, locality on declared-disjoint supports, and identity repair
are supplied exactly as for `simEarnedNet`.

## Receipts

* `sourceRegionalS64_noncommutative` and
  `sourceRegionalNet_noncommutative_iff`: the regional algebra of the
  `s64` support is noncommutative, and it is the one support region of
  the net with that property.
* `m2Block18_mem` / `m2Block28_mem` / `fibreUnit_not_mem_s64`: the full
  matrix-unit block of each realized split fibre lies in the enriched
  algebra, and every matrix unit off the two fibres is excluded, so the
  noncommutative sector is exactly the realized one.
* `diagonal_lt_sourceRegionalS64`: the enrichment is strict over the
  diagonal regional algebra of `simEarnedNet`.
* `sourceRegionalNet_disjoint_supports` and
  `sourceRegionalNet_separation`: region separation; every pair of
  distinct supports is declared disjoint and the enriched algebra
  commutes elementwise with every remote regional algebra.
* `sourceRegionalCover_coverageLaw` and
  `sourceRegionalCover_familyJoin`: through the machinery of
  `QFT/CoverageReceipt.lean`, the four-support family under the top
  region satisfies the coverage law, and its family join is computed to
  equal the top regional algebra.  The negative control
  `sourceRegionalDiagonalCover_not_coverageLaw` shows that the family
  without the `s64` support fails the law: the missing region carries
  the entire noncommutative sector.

## Claim boundary

Realized finite data from one committed run under a declared binning
and the witness payload's observer selection; the fibre tables are
exact counts on this run's freezeout fields, the split-fibre and
designation rules are declared extraction parameters with no canonical
status, and no physical claim and no probability claim beyond realized
frequencies are asserted.  The run statistic that bounds the
noncommutative sector is the payload field
`truncated_split_fibre_count` at `support_truncation = 8`: the sector
grows exactly with the observers and record classes for which
`splitRecord` becomes true (`enrichedRel_of_no_split_fibre`).  The
payload's `full_support_fibre_statistics` report 28, 23, 26, and 24
multi-companion record classes on the four full 96-node supports, so a
deeper truncation of the same committed run is the named route to a
payload in which every observer carries a split fibre; until such a
payload is extracted and mirrored, the three collapsed supports stay
diagonal here.
-/

namespace OPH.QFT

open Matrix
open OPH.Tower

set_option maxRecDepth 32768

/-! ## The split-fibre rule as a decidable predicate on the literals -/

/-- The payload's `split_fibre_rule` on the truncated support of one
observer: record class `r` is a split fibre when two support positions
carry `r` with distinct companion classes. -/
def splitRecord (o : PayloadObserver) (r : Fin 32) : Bool :=
  decide (∃ j k : Fin 8, PayloadObserver.recordClassAt o j = r ∧
    PayloadObserver.recordClassAt o k = r ∧
    PayloadObserver.companionClassAt o j ≠ PayloadObserver.companionClassAt o k)

/-- Mirror of the payload field `truncated_split_fibre_record_classes`
for observer 64: exactly the record classes 18 and 28 are split. -/
theorem splitRecord_o64_char : ∀ r : Fin 32,
    splitRecord .o64 r = (decide (r = 18) || decide (r = 28)) := by decide

/-- Observer 36 has no split fibre at truncation 8. -/
theorem splitRecord_o36_char : ∀ r : Fin 32, splitRecord .o36 r = false := by
  decide

/-- Observer 60 has no split fibre at truncation 8. -/
theorem splitRecord_o60_char : ∀ r : Fin 32, splitRecord .o60 r = false := by
  decide

/-- Observer 92 has no split fibre at truncation 8. -/
theorem splitRecord_o92_char : ∀ r : Fin 32, splitRecord .o92 r = false := by
  decide

/-- Mirror of the payload's `enlargement_statistic.values`: the
truncated split-fibre counts of the four observers are 0, 0, 2, 0. -/
theorem truncatedSplitFibreCounts :
    (Finset.univ.filter fun r : Fin 32 => splitRecord .o36 r = true).card = 0 ∧
    (Finset.univ.filter fun r : Fin 32 => splitRecord .o60 r = true).card = 0 ∧
    (Finset.univ.filter fun r : Fin 32 => splitRecord .o64 r = true).card = 2 ∧
    (Finset.univ.filter fun r : Fin 32 => splitRecord .o92 r = true).card = 0 := by
  decide

/-- Mirror of the payload's first split fibre: record class 18 is
realized exactly at support positions 0 and 3 of observer 64, with the
distinct companion classes 13 and 3.  This fibre is the payload's
designated block under the declared first-fibre rule. -/
theorem splitFibre18_receipt :
    PayloadObserver.recordClassAt .o64 0 = 18 ∧
    PayloadObserver.recordClassAt .o64 3 = 18 ∧
    PayloadObserver.companionClassAt .o64 0 = 13 ∧
    PayloadObserver.companionClassAt .o64 3 = 3 ∧
    (∀ j : Fin 8, PayloadObserver.recordClassAt .o64 j = 18 →
      j = 0 ∨ j = 3) := by
  decide

/-- Mirror of the payload's second split fibre: record class 28 is
realized exactly at support positions 4 and 7 of observer 64, with the
distinct companion classes 13 and 6. -/
theorem splitFibre28_receipt :
    PayloadObserver.recordClassAt .o64 4 = 28 ∧
    PayloadObserver.recordClassAt .o64 7 = 28 ∧
    PayloadObserver.companionClassAt .o64 4 = 13 ∧
    PayloadObserver.companionClassAt .o64 7 = 6 ∧
    (∀ j : Fin 8, PayloadObserver.recordClassAt .o64 j = 28 →
      j = 4 ∨ j = 7) := by
  decide

/-- The designated-block rule receipt: position 0 is the first support
position of observer 64 whose record class is split. -/
theorem designatedBlock_receipt :
    ((List.finRange 8).filter fun j =>
      splitRecord .o64 (PayloadObserver.recordClassAt .o64 j)).head? =
        some 0 := by
  decide

/-- The fibre of record class 18 occupies basis indices 17 and 20. -/
theorem splitFibre18_basis :
    nodeIndex .o64 0 = 17 ∧ nodeIndex .o64 3 = 20 := by decide

/-- The fibre of record class 28 occupies basis indices 21 and 24. -/
theorem splitFibre28_basis :
    nodeIndex .o64 4 = 21 ∧ nodeIndex .o64 7 = 24 := by decide

/-! ## The enriched block relation -/

/-- Two basis indices lie in one realized split fibre of observer `o`:
both are support nodes of `o`, they carry one record class, and that
class satisfies the payload's split rule. -/
def splitPair (o : PayloadObserver) (i j : Fin 33) : Bool :=
  decide (nodeOwner i = some o) && decide (nodeOwner j = some o) &&
    decide (nodeRecordClass i = nodeRecordClass j) &&
    splitRecord o (nodeRecordClass i)

/-- The enriched block relation of one observer support: equality
together with the realized split pairs. -/
def enrichedRel (o : PayloadObserver) (i j : Fin 33) : Bool :=
  decide (i = j) || splitPair o i j

/-- The ambient marker of one observer support: indices owned by any
other observer or by no observer. -/
def enrichedAmb (o : PayloadObserver) (i : Fin 33) : Bool :=
  decide (nodeOwner i ≠ some o)

/-- The enriched block relation of the top region: equality together
with the realized split pairs of every payload observer. -/
def enrichedRelTop (i j : Fin 33) : Bool :=
  decide (i = j) || (splitPair .o36 i j || splitPair .o60 i j ||
    splitPair .o64 i j || splitPair .o92 i j)

/-- The block representative map realized by the payload: the two split
fibres `{17, 20}` and `{21, 24}` collapse to their least members and
every other index represents itself. -/
def splitRep (i : Fin 33) : Fin 33 :=
  if i = 20 then 17 else if i = 24 then 21 else i

/-- Characterization receipt: the enriched relation of observer 64 is
the kernel of the representative map. -/
theorem enrichedRel_o64_char : ∀ i j,
    enrichedRel .o64 i j = decide (splitRep i = splitRep j) := by decide

/-- Collapse characterization for observer 36. -/
theorem enrichedRel_o36_char : ∀ i j,
    enrichedRel .o36 i j = decide (i = j) := by decide

/-- Collapse characterization for observer 60. -/
theorem enrichedRel_o60_char : ∀ i j,
    enrichedRel .o60 i j = decide (i = j) := by decide

/-- Collapse characterization for observer 92. -/
theorem enrichedRel_o92_char : ∀ i j,
    enrichedRel .o92 i j = decide (i = j) := by decide

/-- The top relation is the kernel of the same representative map: only
observer 64 contributes split pairs. -/
theorem enrichedRelTop_char : ∀ i j,
    enrichedRelTop i j = decide (splitRep i = splitRep j) := by decide

/-- The top relation coincides with the enriched relation of observer
64. -/
theorem enrichedRelTop_eq_o64 : ∀ i j,
    enrichedRelTop i j = enrichedRel .o64 i j := fun i j =>
  (enrichedRelTop_char i j).trans (enrichedRel_o64_char i j).symm

/-- Ambient indices of the `s64` support are relation singletons. -/
theorem enrichedAmb_o64_singleton : ∀ i j, enrichedAmb .o64 i = true →
    enrichedRel .o64 i j = true → i = j := by decide

/-- Index 0 is a singleton of the top relation. -/
theorem enrichedRelTop_zero_singleton : ∀ k,
    enrichedRelTop 0 k = true → k = 0 := by decide

/-- Index 0 is a singleton of the `s64` relation. -/
theorem enrichedRel_o64_zero_singleton : ∀ k,
    enrichedRel .o64 0 k = true → k = 0 := by decide

/-- Off-diagonal related pairs of the top relation are owned by
observer 64 on both sides. -/
theorem enrichedRelTop_offdiag_owner : ∀ i j,
    enrichedRelTop i j = true → i ≠ j →
      nodeOwner i = some .o64 ∧ nodeOwner j = some .o64 := by decide

/-- Off-diagonal related pairs of the `s64` relation are owned by
observer 64 on both sides. -/
theorem enrichedRel_o64_offdiag_owner : ∀ i j,
    enrichedRel .o64 i j = true → i ≠ j →
      nodeOwner i = some .o64 ∧ nodeOwner j = some .o64 := by decide

/-- The `s64` relation stays inside the support of observer 64. -/
theorem enrichedRel_o64_owner_closed : ∀ i j,
    enrichedRel .o64 i j = true → nodeOwner i = some .o64 →
      nodeOwner j = some .o64 := by decide

/-- Indices outside the support of observer 64 are relation
singletons. -/
theorem enrichedRel_o64_not_owner_singleton : ∀ i j,
    ¬ nodeOwner i = some .o64 → enrichedRel .o64 i j = true → i = j := by
  decide

/-- The collapse image of the `s36` support consists of singletons of
the top relation. -/
theorem collapse_singleton_o36 : ∀ i k,
    enrichedRelTop (collapseMap .o36 i) k = true →
      k = collapseMap .o36 i := by decide

/-- The collapse image of the `s60` support consists of singletons of
the top relation. -/
theorem collapse_singleton_o60 : ∀ i k,
    enrichedRelTop (collapseMap .o60 i) k = true →
      k = collapseMap .o60 i := by decide

/-- The collapse image of the `s92` support consists of singletons of
the top relation. -/
theorem collapse_singleton_o92 : ∀ i k,
    enrichedRelTop (collapseMap .o92 i) k = true →
      k = collapseMap .o92 i := by decide

/-- The ambient basis index is the one index owned by no observer. -/
theorem nodeOwner_none_iff : ∀ i : Fin 33, nodeOwner i = none ↔ i = 0 := by
  decide

theorem collapseMap_of_owner {o : PayloadObserver} {i : Fin 33}
    (h : nodeOwner i = some o) : collapseMap o i = i := by
  simp [collapseMap, h]

theorem collapseMap_of_not_owner {o : PayloadObserver} {i : Fin 33}
    (h : ¬ nodeOwner i = some o) : collapseMap o i = 0 := by
  simp [collapseMap, h]

/-- Boundary receipt naming the run statistic that carries the
noncommutative sector: on any observer with no split record class at
the declared truncation, the enrichment rule collapses to equality, so
its block algebra is the diagonal one.  A future payload enlarges the
noncommutative sector exactly where it makes `splitRecord` true. -/
theorem enrichedRel_of_no_split_fibre (o : PayloadObserver)
    (h : ∀ r : Fin 32, splitRecord o r = false) :
    ∀ i j, enrichedRel o i j = decide (i = j) := by
  intro i j
  unfold enrichedRel splitPair
  rw [h (nodeRecordClass i)]
  simp

private theorem bool_eq_false {b : Bool} (h : ¬ b = true) : b = false := by
  cases b
  · rfl
  · exact absurd rfl h

private theorem ne_true_of_eq_false {b : Bool} (h : b = false) :
    ¬ b = true := fun h2 => by
  rw [h] at h2
  exact Bool.noConfusion h2

/-! ## Block relations and their subalgebras -/

/-- A Boolean block relation on the basis: an equivalence given by a
reflexive, symmetric, transitive Boolean relation, together with a
marker of ambient indices that must be relation singletons.  The
subalgebra of a block relation consists of the matrices supported on
related pairs whose diagonal is constant across the ambient marker. -/
structure BlockRelation where
  rel : Fin 33 → Fin 33 → Bool
  amb : Fin 33 → Bool
  rel_refl : ∀ i, rel i i = true
  rel_symm : ∀ i j, rel i j = rel j i
  rel_trans : ∀ i j k, rel i j = true → rel j k = true → rel i k = true
  amb_singleton : ∀ i j, amb i = true → rel i j = true → i = j

namespace BlockRelation

/-- Build a block relation from a representative-map characterization:
the equivalence laws follow from the kernel form. -/
def ofRep (rel : Fin 33 → Fin 33 → Bool) (amb : Fin 33 → Bool)
    (f : Fin 33 → Fin 33) (hchar : ∀ i j, rel i j = decide (f i = f j))
    (hamb : ∀ i j, amb i = true → rel i j = true → i = j) :
    BlockRelation where
  rel := rel
  amb := amb
  rel_refl i := by rw [hchar]; exact decide_eq_true rfl
  rel_symm i j := by
    rw [hchar, hchar]
    exact decide_eq_decide.mpr ⟨Eq.symm, Eq.symm⟩
  rel_trans i j k h1 h2 := by
    rw [hchar] at h1 h2 ⊢
    exact decide_eq_true ((of_decide_eq_true h1).trans (of_decide_eq_true h2))
  amb_singleton := hamb

/-- The block subalgebra of a block relation: matrices supported on
related pairs, constant on the diagonal across ambient indices. -/
def subalgebra (B : BlockRelation) :
    StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ) where
  carrier := {X | (∀ i j, B.rel i j = false → X i j = 0) ∧
    (∀ i j, B.amb i = true → B.amb j = true → X i i = X j j)}
  zero_mem' := ⟨fun _ _ _ => rfl, fun _ _ _ _ => rfl⟩
  one_mem' := by
    refine ⟨fun i j h => ?_, fun i j _ _ => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((B.rel_refl i).symm.trans h)
      exact Matrix.one_apply_ne hij
    · rw [Matrix.one_apply_eq, Matrix.one_apply_eq]
  add_mem' := by
    rintro X Y ⟨hX0, hXs⟩ ⟨hY0, hYs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.add_apply, hX0 i j h, hY0 i j h, add_zero]
    · rw [Matrix.add_apply, Matrix.add_apply, hXs i j hi hj, hYs i j hi hj]
  mul_mem' := by
    rintro X Y ⟨hX0, hXs⟩ ⟨hY0, hYs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.mul_apply]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hik : B.rel i k = true
      · have hkj : B.rel k j = false := by
          cases hkj' : B.rel k j
          · rfl
          · exact absurd (B.rel_trans i k j hik hkj')
              (fun htrue => Bool.noConfusion (h.symm.trans htrue))
        rw [hY0 k j hkj, mul_zero]
      · rw [hX0 i k (bool_eq_false hik), zero_mul]
    · rw [Matrix.mul_apply, Matrix.mul_apply]
      rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
          (fun k _ hk => by
            have h0 : B.rel i k = false := by
              cases hik : B.rel i k
              · rfl
              · exact absurd (B.amb_singleton i k hi hik).symm hk
            rw [hX0 i k h0, zero_mul]),
        Finset.sum_eq_single_of_mem j (Finset.mem_univ j)
          (fun k _ hk => by
            have h0 : B.rel j k = false := by
              cases hjk : B.rel j k
              · rfl
              · exact absurd (B.amb_singleton j k hj hjk).symm hk
            rw [hX0 j k h0, zero_mul]),
        hXs i j hi hj, hYs i j hi hj]
  algebraMap_mem' := by
    intro z
    refine ⟨fun i j h => ?_, fun i j _ _ => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((B.rel_refl i).symm.trans h)
      rw [← diagonal_const_eq_algebraMap, Matrix.diagonal_apply_ne _ hij]
    · rw [← diagonal_const_eq_algebraMap, Matrix.diagonal_apply_eq,
        Matrix.diagonal_apply_eq]
  star_mem' := by
    rintro X ⟨hX0, hXs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.star_apply, hX0 j i ((B.rel_symm j i).trans h), star_zero]
    · rw [Matrix.star_apply, Matrix.star_apply, hXs i j hi hj]

theorem mem_subalgebra_iff (B : BlockRelation)
    {X : Matrix (Fin 33) (Fin 33) ℂ} :
    X ∈ B.subalgebra ↔
      (∀ i j, B.rel i j = false → X i j = 0) ∧
        (∀ i j, B.amb i = true → B.amb j = true → X i i = X j j) :=
  Iff.rfl

end BlockRelation

/-- The block relation of the enriched `s64` support region. -/
def o64BlockRelation : BlockRelation :=
  .ofRep (enrichedRel .o64) (enrichedAmb .o64) splitRep
    enrichedRel_o64_char enrichedAmb_o64_singleton

/-- The block relation of the enriched top region: no ambient scalar
constraint. -/
def topBlockRelation : BlockRelation :=
  .ofRep enrichedRelTop (fun _ => false) splitRep enrichedRelTop_char
    (fun _ _ h _ => Bool.noConfusion h)

/-- Reflexivity restated on the payload-derived relation. -/
theorem enrichedRel_o64_refl : ∀ i, enrichedRel .o64 i i = true :=
  o64BlockRelation.rel_refl

/-- Symmetry restated on the payload-derived relation. -/
theorem enrichedRel_o64_symm : ∀ i j,
    enrichedRel .o64 i j = enrichedRel .o64 j i :=
  o64BlockRelation.rel_symm

/-- The block relation of the collapsed `s36` support. -/
def o36BlockRelation : BlockRelation :=
  .ofRep (enrichedRel .o36) (enrichedAmb .o36) id enrichedRel_o36_char
    (fun i j _ h =>
      of_decide_eq_true ((enrichedRel_o36_char i j).symm.trans h))

/-- The block relation of the collapsed `s60` support. -/
def o60BlockRelation : BlockRelation :=
  .ofRep (enrichedRel .o60) (enrichedAmb .o60) id enrichedRel_o60_char
    (fun i j _ h =>
      of_decide_eq_true ((enrichedRel_o60_char i j).symm.trans h))

/-- The block relation of the collapsed `s92` support. -/
def o92BlockRelation : BlockRelation :=
  .ofRep (enrichedRel .o92) (enrichedAmb .o92) id enrichedRel_o92_char
    (fun i j _ h =>
      of_decide_eq_true ((enrichedRel_o92_char i j).symm.trans h))

/-- The enriched regional algebra of the split-fibre support. -/
def sourceRegionalS64 : StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ) :=
  o64BlockRelation.subalgebra

/-- The enriched top algebra: block-diagonal matrices with two free
two-by-two blocks on the realized split fibres and free diagonal
elsewhere. -/
def sourceRegionalTop : StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ) :=
  topBlockRelation.subalgebra

/-! ## The collapse receipts

On each observer support with no realized split fibre, the block
algebra of the enrichment rule equals the diagonal regional algebra of
`simEarnedNet`.  The truncated payload data cannot support a
noncommutative receipt at these supports, so the collapse is proved. -/

private theorem block_collapse_of_diag {o : PayloadObserver}
    (B : BlockRelation) (hrel : B.rel = enrichedRel o)
    (hamb : B.amb = enrichedAmb o)
    (hchar : ∀ i j, enrichedRel o i j = decide (i = j)) :
    B.subalgebra = diagonalClassSubalgebra (fReg o) := by
  apply le_antisymm
  · rintro X ⟨h0, hs⟩
    refine ⟨fun x => x.elim (X 0 0) fun k => X k k, ?_⟩
    ext i j
    by_cases hij : i = j
    · subst hij
      rw [Matrix.diagonal_apply_eq]
      by_cases ho : nodeOwner i = some o
      · simp [fReg, ho]
      · have hi : B.amb i = true := by
          rw [hamb]; exact decide_eq_true ho
        have h00 : B.amb 0 = true := by
          rw [hamb]
          refine decide_eq_true ?_
          rw [nodeOwner_zero]
          exact fun h => nomatch h
        have hnone : fReg o i = none := by simp [fReg, ho]
        rw [hnone]
        exact hs 0 i h00 hi
    · rw [Matrix.diagonal_apply_ne _ hij]
      refine (h0 i j ?_).symm
      rw [hrel, hchar]
      exact decide_eq_false hij
  · rintro X ⟨c, rfl⟩
    refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((B.rel_refl i).symm.trans h0)
      exact Matrix.diagonal_apply_ne _ hij
    · rw [hamb] at hi hj
      have hoi : ¬ nodeOwner i = some o := of_decide_eq_true hi
      have hoj : ¬ nodeOwner j = some o := of_decide_eq_true hj
      rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]
      simp [fReg, hoi, hoj]

/-- Collapse receipt: the enriched algebra of the `s36` support is the
diagonal regional algebra of `simEarnedNet`. -/
theorem enrichedRegion_collapse_o36 :
    o36BlockRelation.subalgebra = diagonalClassSubalgebra (fReg .o36) :=
  block_collapse_of_diag o36BlockRelation rfl rfl enrichedRel_o36_char

/-- Collapse receipt for the `s60` support. -/
theorem enrichedRegion_collapse_o60 :
    o60BlockRelation.subalgebra = diagonalClassSubalgebra (fReg .o60) :=
  block_collapse_of_diag o60BlockRelation rfl rfl enrichedRel_o60_char

/-- Collapse receipt for the `s92` support. -/
theorem enrichedRegion_collapse_o92 :
    o92BlockRelation.subalgebra = diagonalClassSubalgebra (fReg .o92) :=
  block_collapse_of_diag o92BlockRelation rfl rfl enrichedRel_o92_char

/-! ## The noncommutative receipts on the realized split fibres -/

/-- The matrix unit at one basis-index pair. -/
def fibreUnit (a b : Fin 33) : Matrix (Fin 33) (Fin 33) ℂ :=
  fun i j => if i = a ∧ j = b then 1 else 0

/-- Every matrix unit of a realized split pair lies in the enriched
regional algebra. -/
theorem fibreUnit_mem_s64 {a b : Fin 33} (h : splitPair .o64 a b = true) :
    fibreUnit a b ∈ sourceRegionalS64 := by
  have hrel : enrichedRel .o64 a b = true := by
    unfold enrichedRel
    rw [h, Bool.or_true]
  have ha : nodeOwner a = some .o64 := by
    unfold splitPair at h
    simp only [Bool.and_eq_true] at h
    exact of_decide_eq_true h.1.1.1
  have hb : nodeOwner b = some .o64 := by
    unfold splitPair at h
    simp only [Bool.and_eq_true] at h
    exact of_decide_eq_true h.1.1.2
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · by_cases hij : i = a ∧ j = b
    · obtain ⟨rfl, rfl⟩ := hij
      exact absurd hrel fun htrue => Bool.noConfusion (h0.symm.trans htrue)
    · simp [fibreUnit, hij]
  · have hia : i ≠ a := fun he => (of_decide_eq_true hi) (he ▸ ha)
    have hja : j ≠ a := fun he => (of_decide_eq_true hj) (he ▸ ha)
    simp only [fibreUnit]
    rw [if_neg fun hc => hia hc.1, if_neg fun hc => hja hc.1]

/-- Exclusion receipt: every matrix unit off the realized relation is
outside the enriched regional algebra, so the noncommutative sector is
exactly the two realized fibres. -/
theorem fibreUnit_not_mem_s64 {a b : Fin 33}
    (h : enrichedRel .o64 a b = false) :
    fibreUnit a b ∉ sourceRegionalS64 := by
  rintro ⟨h0, _⟩
  have := h0 a b h
  simp [fibreUnit] at this

/-- The full matrix-unit block of the designated fibre (record class
18, basis indices 17 and 20) lies in the enriched algebra. -/
theorem m2Block18_mem : ∀ a b : Fin 33, (a = 17 ∨ a = 20) →
    (b = 17 ∨ b = 20) → fibreUnit a b ∈ sourceRegionalS64 := by
  rintro a b (rfl | rfl) (rfl | rfl) <;>
    exact fibreUnit_mem_s64 (by decide)

/-- The full matrix-unit block of the second realized fibre (record
class 28, basis indices 21 and 24) lies in the enriched algebra. -/
theorem m2Block28_mem : ∀ a b : Fin 33, (a = 21 ∨ a = 24) →
    (b = 21 ∨ b = 24) → fibreUnit a b ∈ sourceRegionalS64 := by
  rintro a b (rfl | rfl) (rfl | rfl) <;>
    exact fibreUnit_mem_s64 (by decide)

/-- The enriched regional algebra of the split-fibre support is
noncommutative: the raising and lowering units of the designated fibre
fail to commute. -/
theorem sourceRegionalS64_noncommutative :
    ∃ X Y : sourceRegionalS64,
      (X : Matrix (Fin 33) (Fin 33) ℂ) * Y ≠
        (Y : Matrix (Fin 33) (Fin 33) ℂ) * X := by
  refine ⟨⟨fibreUnit 17 20, m2Block18_mem 17 20 (Or.inl rfl) (Or.inr rfl)⟩,
    ⟨fibreUnit 20 17, m2Block18_mem 20 17 (Or.inr rfl) (Or.inl rfl)⟩, ?_⟩
  intro h
  have hL : (fibreUnit 17 20 * fibreUnit 20 17) 17 17 = 1 := by
    rw [Matrix.mul_apply]
    rw [Finset.sum_eq_single_of_mem (20 : Fin 33) (Finset.mem_univ _)
      (fun k _ hk => by simp [fibreUnit, hk])]
    simp [fibreUnit]
  have hR : (fibreUnit 20 17 * fibreUnit 17 20) 17 17 = 0 := by
    rw [Matrix.mul_apply]
    refine Finset.sum_eq_zero fun k _ => ?_
    simp [fibreUnit]
  have h17 := congrFun (congrFun h 17) 17
  rw [hL, hR] at h17
  exact one_ne_zero h17

/-- The enrichment is strict: the diagonal regional algebra of
`simEarnedNet` at the `s64` support is a proper subalgebra of the
enriched one. -/
theorem diagonal_lt_sourceRegionalS64 :
    diagonalClassSubalgebra (fReg .o64) < sourceRegionalS64 := by
  rw [SetLike.lt_iff_le_and_exists]
  constructor
  · rintro X ⟨c, rfl⟩
    refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((o64BlockRelation.rel_refl i).symm.trans h0)
      exact Matrix.diagonal_apply_ne _ hij
    · rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]
      have hoi : ¬ nodeOwner i = some .o64 := of_decide_eq_true hi
      have hoj : ¬ nodeOwner j = some .o64 := of_decide_eq_true hj
      simp [fReg, hoi, hoj]
  · refine ⟨fibreUnit 17 20,
      m2Block18_mem 17 20 (Or.inl rfl) (Or.inr rfl), ?_⟩
    rintro ⟨c, hc⟩
    have h1720 := congrFun (congrFun hc 17) 20
    rw [Matrix.diagonal_apply_ne _ (by decide : (17 : Fin 33) ≠ 20)] at h1720
    simp [fibreUnit] at h1720

/-! ## The enriched regional assignment, isotony, and locality -/

/-- The enriched regional algebra assignment: scalars at the bottom,
the unchanged diagonal spans at the three collapsed supports, the
enriched block algebra at the split-fibre support, and the enriched
block-diagonal algebra at the top. -/
def sourceRegionalAlgebra :
    SimRegion → StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ)
  | .bot => diagonalClassSubalgebra fBot
  | .s36 => diagonalClassSubalgebra (fReg .o36)
  | .s60 => diagonalClassSubalgebra (fReg .o60)
  | .s64 => sourceRegionalS64
  | .s92 => diagonalClassSubalgebra (fReg .o92)
  | .top => sourceRegionalTop

/-- Every diagonal matrix lies in the enriched top algebra. -/
theorem diagonal_mem_sourceRegionalTop (d : Fin 33 → ℂ) :
    Matrix.diagonal d ∈ sourceRegionalTop := by
  refine ⟨fun i j h0 => ?_, fun i j hi _ => Bool.noConfusion hi⟩
  have hij : i ≠ j := fun he => by
    subst he
    exact Bool.noConfusion ((topBlockRelation.rel_refl i).symm.trans h0)
  exact Matrix.diagonal_apply_ne _ hij

/-- The scalar bottom algebra lies in every block subalgebra. -/
theorem scalars_le_block (B : BlockRelation) :
    diagonalClassSubalgebra fBot ≤ B.subalgebra := by
  rintro X ⟨c, rfl⟩
  refine ⟨fun i j h0 => ?_, fun i j _ _ => ?_⟩
  · have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion ((B.rel_refl i).symm.trans h0)
    exact Matrix.diagonal_apply_ne _ hij
  · rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]

/-- Every diagonal support span lies in the enriched top algebra. -/
theorem span_fReg_le_top (o : PayloadObserver) :
    diagonalClassSubalgebra (fReg o) ≤ sourceRegionalTop := by
  rintro X ⟨c, rfl⟩
  exact diagonal_mem_sourceRegionalTop _

/-- The enriched support algebra lies in the enriched top algebra. -/
theorem s64_le_top : sourceRegionalS64 ≤ sourceRegionalTop := by
  rintro X ⟨h0, _⟩
  refine ⟨fun i j h => h0 i j ((enrichedRelTop_eq_o64 i j).symm.trans h),
    fun i j hi _ => Bool.noConfusion hi⟩

/-- Isotony of the enriched regional assignment. -/
theorem sourceRegionalAlgebra_isotony {U V : SimRegion}
    (h : SimRegion.le U V) :
    sourceRegionalAlgebra U ≤ sourceRegionalAlgebra V := by
  cases U <;> cases V <;>
    first
      | exact le_rfl
      | exact absurd h (by decide)
      | exact diagonalClassSubalgebra_comp_le (fun _ => ()) _
      | exact scalars_le_block _
      | exact span_fReg_le_top _
      | exact s64_le_top

/-- Two members of diagonal spans commute. -/
theorem span_commute_span {α β : Type} {f : Fin 33 → α} {g : Fin 33 → β}
    {X Y : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ diagonalClassSubalgebra f)
    (hY : Y ∈ diagonalClassSubalgebra g) : Commute X Y := by
  obtain ⟨c, rfl⟩ := hX
  obtain ⟨d, rfl⟩ := hY
  show Matrix.diagonal _ * Matrix.diagonal _ =
    Matrix.diagonal _ * Matrix.diagonal _
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm _ _

/-- A member of a remote diagonal support span commutes with every
member of the enriched algebra: the remote span is constant across the
`s64` support, where the enriched off-diagonal entries live. -/
theorem span_commute_block (o : PayloadObserver) (ho : o ≠ .o64)
    {X Y : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ diagonalClassSubalgebra (fReg o))
    (hY : Y ∈ sourceRegionalS64) : Commute X Y := by
  obtain ⟨c, rfl⟩ := hX
  obtain ⟨h0, _⟩ := hY
  show Matrix.diagonal _ * Y = Y * Matrix.diagonal _
  ext i j
  rw [Matrix.diagonal_mul, Matrix.mul_diagonal]
  by_cases hrel : enrichedRel .o64 i j = true
  · by_cases hij : i = j
    · subst hij
      exact mul_comm _ _
    · obtain ⟨hoi, hoj⟩ := enrichedRel_o64_offdiag_owner i j hrel hij
      have hne : some PayloadObserver.o64 ≠ some o :=
        fun hcontra => ho (Option.some.inj hcontra).symm
      have hfi : fReg o i = none := by
        simp [fReg, hoi, hne]
      have hfj : fReg o j = none := by
        simp [fReg, hoj, hne]
      rw [hfi, hfj, mul_comm]
  · rw [h0 i j (bool_eq_false hrel), mul_zero, zero_mul]

/-- Elementwise commutation across every declared-disjoint region pair
of the enriched assignment. -/
theorem sourceRegional_locality {U V : SimRegion} (h : SimRegion.disj U V)
    (X : sourceRegionalAlgebra U) (Y : sourceRegionalAlgebra V) :
    Commute (X : Matrix (Fin 33) (Fin 33) ℂ)
      (Y : Matrix (Fin 33) (Fin 33) ℂ) := by
  cases U <;> cases V <;>
    first
      | exact absurd rfl h.1
      | exact absurd rfl h.2.1
      | exact absurd rfl h.2.2.1
      | exact absurd rfl h.2.2.2.1
      | exact absurd rfl h.2.2.2.2
      | exact span_commute_span X.2 Y.2
      | exact span_commute_block _ (by decide) X.2 Y.2
      | exact (span_commute_block _ (by decide) Y.2 X.2).symm

/-! ## Restriction maps

Restrictions onto diagonal targets pull the diagonal back along a
basis self-map whose image consists of relation singletons; the
restriction onto the enriched support keeps the block entries and
collapses the ambient diagonal to the ambient scalar. -/

/-- The character-style restriction from a block subalgebra onto a
diagonal span, along a basis self-map with singleton image. -/
noncomputable def BlockRelation.charHom {β : Type} (B : BlockRelation)
    (g : Fin 33 → β) (σ : Fin 33 → Fin 33)
    (hσ : ∀ i k, B.rel (σ i) k = true → k = σ i)
    (hmem : ∀ X ∈ B.subalgebra, diagPull σ X ∈ diagonalClassSubalgebra g) :
    B.subalgebra →⋆ₐ[ℂ] diagonalClassSubalgebra g where
  toFun X := ⟨diagPull σ X.1, hmem X.1 X.2⟩
  map_one' := Subtype.ext (diagPull_one σ)
  map_mul' X Y := Subtype.ext (by
    show diagPull σ (X.1 * Y.1) = diagPull σ X.1 * diagPull σ Y.1
    unfold diagPull
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    rw [Matrix.mul_apply]
    refine Finset.sum_eq_single_of_mem (σ i) (Finset.mem_univ _)
      fun k _ hk => ?_
    have hrel : B.rel (σ i) k = false := by
      cases hbk : B.rel (σ i) k
      · rfl
      · exact absurd (hσ i k hbk) hk
    rw [X.2.1 (σ i) k hrel, zero_mul])
  map_zero' := Subtype.ext (diagPull_zero σ)
  map_add' X Y := Subtype.ext (diagPull_add σ X.1 Y.1)
  map_star' X := Subtype.ext (diagPull_star σ X.1)
  commutes' z := Subtype.ext (diagPull_algebraMap σ z)

/-- The constant collapse lands in the scalar bottom span. -/
theorem constPull_mem_bot (X : Matrix (Fin 33) (Fin 33) ℂ) :
    diagPull (fun _ => 0) X ∈ diagonalClassSubalgebra fBot :=
  ⟨fun _ => X 0 0, rfl⟩

/-- The support collapse lands in the diagonal support span. -/
theorem topPull_mem_span (o : PayloadObserver)
    (X : Matrix (Fin 33) (Fin 33) ℂ) :
    diagPull (collapseMap o) X ∈ diagonalClassSubalgebra (fReg o) := by
  refine ⟨fun x => x.elim (X 0 0) fun k => X k k, ?_⟩
  show Matrix.diagonal _ = Matrix.diagonal _
  congr 1
  funext i
  by_cases h : nodeOwner i = some o
  · simp [fReg, h, collapseMap]
  · simp [fReg, h, collapseMap]

/-- The character onto the bottom scalars from a block subalgebra. -/
noncomputable def sourceCharToBot (B : BlockRelation)
    (h0 : ∀ k, B.rel 0 k = true → k = 0) :
    B.subalgebra →⋆ₐ[ℂ] diagonalClassSubalgebra fBot :=
  B.charHom fBot (fun _ => 0) (fun _ k => h0 k)
    (fun X _ => constPull_mem_bot X)

/-- The restriction from the enriched top onto one collapsed diagonal
support. -/
noncomputable def sourceTopToSupp (o : PayloadObserver)
    (hσ : ∀ i k, enrichedRelTop (collapseMap o i) k = true →
      k = collapseMap o i) :
    sourceRegionalTop →⋆ₐ[ℂ] diagonalClassSubalgebra (fReg o) :=
  topBlockRelation.charHom (fReg o) (collapseMap o) hσ
    (fun X _ => topPull_mem_span o X)

/-- The matrix value of the enriched restriction: block entries are
kept on the `s64` support, the ambient diagonal collapses to the
ambient scalar, and everything else vanishes. -/
def blockPull64 (X : Matrix (Fin 33) (Fin 33) ℂ) :
    Matrix (Fin 33) (Fin 33) ℂ :=
  fun i j => if enrichedRel .o64 i j = true then
    X (collapseMap .o64 i) (collapseMap .o64 j) else 0

theorem blockPull64_mem (X : Matrix (Fin 33) (Fin 33) ℂ) :
    blockPull64 X ∈ sourceRegionalS64 := by
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · have h0' : enrichedRel .o64 i j = false := h0
    simp only [blockPull64]
    rw [if_neg (ne_true_of_eq_false h0')]
  · have hoi : ¬ nodeOwner i = some .o64 := of_decide_eq_true hi
    have hoj : ¬ nodeOwner j = some .o64 := of_decide_eq_true hj
    simp only [blockPull64]
    rw [if_pos (enrichedRel_o64_refl i), if_pos (enrichedRel_o64_refl j),
      collapseMap_of_not_owner hoi, collapseMap_of_not_owner hoj]

/-- The enriched restriction fixes every member of the enriched
support algebra. -/
theorem blockPull64_fixes {X : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ sourceRegionalS64) : blockPull64 X = X := by
  ext i j
  simp only [blockPull64]
  by_cases hrel : enrichedRel .o64 i j = true
  · rw [if_pos hrel]
    by_cases hij : i = j
    · subst hij
      by_cases ho : nodeOwner i = some .o64
      · rw [collapseMap_of_owner ho]
      · rw [collapseMap_of_not_owner ho]
        exact hX.2 0 i (by decide) (decide_eq_true ho)
    · obtain ⟨hoi, hoj⟩ := enrichedRel_o64_offdiag_owner i j hrel hij
      rw [collapseMap_of_owner hoi, collapseMap_of_owner hoj]
  · rw [if_neg hrel]
    exact (hX.1 i j (bool_eq_false hrel)).symm

/-- The restriction from the enriched top onto the enriched support:
a star-algebra homomorphism that keeps the two realized blocks. -/
noncomputable def enrichedPullback :
    sourceRegionalTop →⋆ₐ[ℂ] sourceRegionalS64 where
  toFun X := ⟨blockPull64 X.1, blockPull64_mem X.1⟩
  map_one' := Subtype.ext (by
    show blockPull64 1 = 1
    ext i j
    simp only [blockPull64]
    by_cases hrel : enrichedRel .o64 i j = true
    · rw [if_pos hrel]
      by_cases hij : i = j
      · subst hij
        rw [Matrix.one_apply_eq, Matrix.one_apply_eq]
      · obtain ⟨hoi, hoj⟩ := enrichedRel_o64_offdiag_owner i j hrel hij
        rw [collapseMap_of_owner hoi, collapseMap_of_owner hoj]
    · have hij : i ≠ j := fun he => by
        subst he
        exact hrel (o64BlockRelation.rel_refl i)
      rw [if_neg hrel, Matrix.one_apply_ne hij])
  map_mul' X Y := Subtype.ext (by
    show blockPull64 (X.1 * Y.1) = blockPull64 X.1 * blockPull64 Y.1
    have hXtop : ∀ i j, enrichedRel .o64 i j = false → X.1 i j = 0 :=
      fun i j h => X.2.1 i j ((enrichedRelTop_eq_o64 i j).trans h)
    have hYtop : ∀ i j, enrichedRel .o64 i j = false → Y.1 i j = 0 :=
      fun i j h => Y.2.1 i j ((enrichedRelTop_eq_o64 i j).trans h)
    ext i j
    by_cases hrel : enrichedRel .o64 i j = true
    · by_cases hoi : nodeOwner i = some .o64
      · have hoj : nodeOwner j = some .o64 := by
          by_cases hij : i = j
          · subst hij; exact hoi
          · exact (enrichedRel_o64_offdiag_owner i j hrel hij).2
        have hL : blockPull64 (X.1 * Y.1) i j = ∑ k, X.1 i k * Y.1 k j := by
          simp only [blockPull64]
          rw [if_pos hrel, collapseMap_of_owner hoi,
            collapseMap_of_owner hoj, Matrix.mul_apply]
        have hR : ∀ k, blockPull64 X.1 i k * blockPull64 Y.1 k j =
            X.1 i k * Y.1 k j := by
          intro k
          by_cases hik : enrichedRel .o64 i k = true
          · have hok : nodeOwner k = some .o64 :=
              enrichedRel_o64_owner_closed i k hik hoi
            by_cases hkj : enrichedRel .o64 k j = true
            · simp only [blockPull64]
              rw [if_pos hik, if_pos hkj, collapseMap_of_owner hoi,
                collapseMap_of_owner hok, collapseMap_of_owner hoj]
            · simp only [blockPull64]
              rw [if_neg hkj, mul_zero,
                hYtop k j (bool_eq_false hkj), mul_zero]
          · simp only [blockPull64]
            rw [if_neg hik, zero_mul,
              hXtop i k (bool_eq_false hik), zero_mul]
        rw [hL, Matrix.mul_apply]
        exact Finset.sum_congr rfl fun k _ => (hR k).symm
      · have hij : i = j := enrichedRel_o64_not_owner_singleton i j hoi hrel
        subst hij
        have hσ : collapseMap .o64 i = 0 := collapseMap_of_not_owner hoi
        have hL : blockPull64 (X.1 * Y.1) i i = X.1 0 0 * Y.1 0 0 := by
          simp only [blockPull64]
          rw [if_pos hrel, hσ, Matrix.mul_apply]
          refine Finset.sum_eq_single_of_mem 0 (Finset.mem_univ _)
            fun m _ hm => ?_
          have h0m : enrichedRel .o64 0 m = false :=
            bool_eq_false fun htrue =>
              hm (enrichedRel_o64_zero_singleton m htrue)
          rw [hXtop 0 m h0m, zero_mul]
        have hR : (blockPull64 X.1 * blockPull64 Y.1) i i =
            X.1 0 0 * Y.1 0 0 := by
          rw [Matrix.mul_apply]
          rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ _)
            (fun k _ hk => by
              have hik : enrichedRel .o64 i k = false :=
                bool_eq_false fun htrue =>
                  hk (enrichedRel_o64_not_owner_singleton i k hoi htrue).symm
              simp only [blockPull64]
              rw [if_neg (ne_true_of_eq_false hik), zero_mul])]
          simp only [blockPull64]
          rw [if_pos hrel, if_pos hrel, hσ]
        exact hL.trans hR.symm
    · have hf : enrichedRel .o64 i j = false := bool_eq_false hrel
      have hL : blockPull64 (X.1 * Y.1) i j = 0 := by
        simp only [blockPull64]
        rw [if_neg hrel]
      have hR : (blockPull64 X.1 * blockPull64 Y.1) i j = 0 := by
        rw [Matrix.mul_apply]
        refine Finset.sum_eq_zero fun k _ => ?_
        by_cases hik : enrichedRel .o64 i k = true
        · by_cases hkj : enrichedRel .o64 k j = true
          · exact absurd (o64BlockRelation.rel_trans i k j hik hkj) hrel
          · simp only [blockPull64]
            rw [if_neg hkj, mul_zero]
        · simp only [blockPull64]
          rw [if_neg hik, zero_mul]
      exact hL.trans hR.symm)
  map_zero' := Subtype.ext (by
    show blockPull64 0 = 0
    ext i j
    simp [blockPull64])
  map_add' X Y := Subtype.ext (by
    show blockPull64 (X.1 + Y.1) = blockPull64 X.1 + blockPull64 Y.1
    ext i j
    simp only [blockPull64, Matrix.add_apply]
    by_cases hrel : enrichedRel .o64 i j = true
    · rw [if_pos hrel, if_pos hrel, if_pos hrel]
    · rw [if_neg hrel, if_neg hrel, if_neg hrel, add_zero])
  map_star' X := Subtype.ext (by
    show blockPull64 (star X.1) = star (blockPull64 X.1)
    ext i j
    have hsym : enrichedRel .o64 i j = enrichedRel .o64 j i :=
      enrichedRel_o64_symm i j
    simp only [blockPull64, Matrix.star_apply]
    by_cases hrel : enrichedRel .o64 i j = true
    · rw [if_pos hrel, if_pos (hsym.symm.trans hrel)]
    · rw [if_neg hrel, if_neg fun h => hrel (hsym.trans h), star_zero])
  commutes' z := Subtype.ext (by
    show blockPull64 (algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) z) =
      algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) z
    rw [← diagonal_const_eq_algebraMap]
    ext i j
    simp only [blockPull64]
    by_cases hrel : enrichedRel .o64 i j = true
    · rw [if_pos hrel]
      by_cases hij : i = j
      · subst hij
        rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]
      · obtain ⟨hoi, hoj⟩ := enrichedRel_o64_offdiag_owner i j hrel hij
        rw [collapseMap_of_owner hoi, collapseMap_of_owner hoj]
    · have hij : i ≠ j := fun he => by
        subst he
        exact hrel (o64BlockRelation.rel_refl i)
      rw [if_neg hrel, Matrix.diagonal_apply_ne _ hij])

/-- The supplied restriction system of the enriched net: identities on
equal regions, characters onto the bottom, support collapses from the
top, and the block restriction onto the enriched support. -/
noncomputable def sourceRegionalRestrict : ∀ U V : SimRegion,
    SimRegion.le V U →
      (sourceRegionalAlgebra U →⋆ₐ[ℂ] sourceRegionalAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .s36, .s36, _ => StarAlgHom.id ℂ _
  | .s60, .s60, _ => StarAlgHom.id ℂ _
  | .s64, .s64, _ => StarAlgHom.id ℂ _
  | .s92, .s92, _ => StarAlgHom.id ℂ _
  | .top, .top, _ => StarAlgHom.id ℂ _
  | .s36, .bot, _ => simCharToBot (fReg .o36)
  | .s60, .bot, _ => simCharToBot (fReg .o60)
  | .s64, .bot, _ =>
      sourceCharToBot o64BlockRelation enrichedRel_o64_zero_singleton
  | .s92, .bot, _ => simCharToBot (fReg .o92)
  | .top, .bot, _ =>
      sourceCharToBot topBlockRelation enrichedRelTop_zero_singleton
  | .top, .s36, _ => sourceTopToSupp .o36 collapse_singleton_o36
  | .top, .s60, _ => sourceTopToSupp .o60 collapse_singleton_o60
  | .top, .s64, _ => enrichedPullback
  | .top, .s92, _ => sourceTopToSupp .o92 collapse_singleton_o92
  | .bot, .s36, h => nomatch h
  | .bot, .s60, h => nomatch h
  | .bot, .s64, h => nomatch h
  | .bot, .s92, h => nomatch h
  | .bot, .top, h => nomatch h
  | .s36, .s60, h => nomatch h
  | .s36, .s64, h => nomatch h
  | .s36, .s92, h => nomatch h
  | .s36, .top, h => nomatch h
  | .s60, .s36, h => nomatch h
  | .s60, .s64, h => nomatch h
  | .s60, .s92, h => nomatch h
  | .s60, .top, h => nomatch h
  | .s64, .s36, h => nomatch h
  | .s64, .s60, h => nomatch h
  | .s64, .s92, h => nomatch h
  | .s64, .top, h => nomatch h
  | .s92, .s36, h => nomatch h
  | .s92, .s60, h => nomatch h
  | .s92, .s64, h => nomatch h
  | .s92, .top, h => nomatch h

/-! ## The enriched earned net -/

/-- The source regional net: the finite causal observer net over the
earned tower whose `s64` regional algebra carries the two realized
matrix blocks of the split fibres, while every other region keeps a
commutative algebra. -/
noncomputable def sourceRegionalNet : FiniteCausalObserverNet simEarnedTower where
  Region := fun _ => SimRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨.bot⟩
  regionLE := fun _ => SimRegion.le
  regionLE_refl := fun _ U => SimRegion.le_refl U
  regionLE_trans := fun _ {U V W} h1 h2 => SimRegion.le_trans U V W h1 h2
  regionLE_antisymm := fun _ {U V} h1 h2 => SimRegion.le_antisymm U V h1 h2
  overlap := fun _ => SimRegion.meet
  overlap_le_left := fun _ U V => SimRegion.meet_le_left U V
  overlap_le_right := fun _ U V => SimRegion.meet_le_right U V
  le_overlap := fun _ {W U V} h1 h2 => SimRegion.le_meet W U V h1 h2
  disjoint := fun _ => SimRegion.disj
  disjoint_symm := fun _ {U V} h =>
    ⟨h.1.symm, h.2.2.2.1, h.2.2.2.2, h.2.1, h.2.2.1⟩
  disjoint_irrefl := fun _ U h => h.1 rfl
  localAlgebra := fun _ => sourceRegionalAlgebra
  isotony := fun _ {U V} h => sourceRegionalAlgebra_isotony h
  locality := fun _ {U V} h X Y => sourceRegional_locality h X Y
  restrict := fun _ {U V} h => sourceRegionalRestrict U V h
  restrict_refl := by
    intro r U X
    cases U <;> rfl
  restrict_trans := by
    intro r U V W hVU hWV X
    cases U <;> cases V <;> cases W <;>
      first
        | rfl
        | exact absurd hVU (by decide)
        | exact absurd hWV (by decide)
  restrict_inclusion := by
    intro r U V hUV X
    cases U <;> cases V <;>
      first
        | rfl
        | exact absurd hUV (by decide)
        | exact Subtype.ext (blockPull64_fixes X.2)
        | (obtain ⟨val, hval⟩ := X
           obtain ⟨c, rfl⟩ := mem_diagonalClassSubalgebra_iff.mp hval
           first
             | rfl
             | (apply Subtype.ext
                exact congrArg Matrix.diagonal (funext fun i =>
                  (Matrix.diagonal_apply_eq _ _).trans
                    (congrArg c (fReg_collapseMap _ i)))))
  regionRefine := fun _ U => U
  region_refine_refl := by intros; rfl
  region_refine_trans := by intros; rfl
  region_refine_mono := by intros; assumption
  overlap_natural := by intros; rfl
  disjoint_natural := by intros; assumption
  localAlgebra_natural := by intros; assumption
  repair := fun _ _ => LinearMap.id
  repair_idempotent := by intros; rfl
  repair_fixes_region := by intros; rfl
  repair_fixes_disjoint := by intros; rfl
  repair_natural := by intros; rfl

/-! ## Region separation receipts -/

/-- Every pair of distinct observer supports is declared disjoint in
the enriched net; the literal receipt
`supportNodes_pairwise_disjoint` backs the declaration. -/
theorem sourceRegionalNet_disjoint_supports :
    ∀ o o' : PayloadObserver, o ≠ o' →
      sourceRegionalNet.disjoint () (SimRegion.regionOf o)
        (SimRegion.regionOf o') := by
  intro o o' h
  cases o <;> cases o' <;>
    first
      | exact absurd rfl h
      | exact ⟨by decide, by decide, by decide, by decide, by decide⟩

/-- Separation of the enriched region: every regional observable of a
remote support commutes elementwise with every observable of the
noncommutative `s64` algebra. -/
theorem sourceRegionalNet_separation (o : PayloadObserver)
    (ho : o ≠ .o64)
    (X : sourceRegionalNet.localAlgebra () (SimRegion.regionOf o))
    (Y : sourceRegionalNet.localAlgebra () (SimRegion.regionOf .o64)) :
    Commute X.1 Y.1 :=
  sourceRegionalNet.locality ()
    (sourceRegionalNet_disjoint_supports o .o64 ho) X Y

/-- The noncommutative sector of the enriched net is exactly the
split-fibre support: a support region carries a noncommuting pair of
regional observables precisely when it is the region of observer 64. -/
theorem sourceRegionalNet_noncommutative_iff (o : PayloadObserver) :
    (∃ X Y : sourceRegionalNet.localAlgebra () (SimRegion.regionOf o),
      X.1 * Y.1 ≠ Y.1 * X.1) ↔ o = .o64 := by
  cases o
  case o64 =>
    exact ⟨fun _ => rfl, fun _ => sourceRegionalS64_noncommutative⟩
  all_goals
    refine ⟨?_, fun h => absurd h (by decide)⟩
    rintro ⟨X, Y, hne⟩
    exact (hne (congrArg Subtype.val
      (diagonalClassSubalgebra_mul_comm _ X Y))).elim

/-! ## The coverage computation for the support family under the top

Through the machinery of `QFT/CoverageReceipt.lean`: the four-support
family satisfies the coverage law for the top region, its family join
is computed to equal the enriched top algebra, and the family without
the `s64` support fails the law because the missing region carries the
whole noncommutative sector. -/

/-- The diagonal projector onto one observer support. -/
def supportProjector (o : PayloadObserver) :
    Matrix (Fin 33) (Fin 33) ℂ :=
  Matrix.diagonal fun i => if nodeOwner i = some o then 1 else 0

theorem supportProjector_mem_span (o : PayloadObserver) :
    supportProjector o ∈ diagonalClassSubalgebra (fReg o) := by
  refine ⟨fun x => x.elim 0 fun _ => 1, ?_⟩
  show Matrix.diagonal _ = Matrix.diagonal _
  congr 1
  funext i
  by_cases h : nodeOwner i = some o
  · simp [fReg, h]
  · simp [fReg, h]

theorem supportProjector_mem_s64 :
    supportProjector .o64 ∈ sourceRegionalS64 := by
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion ((enrichedRel_o64_refl i).symm.trans h0)
    exact Matrix.diagonal_apply_ne _ hij
  · simp only [supportProjector, Matrix.diagonal_apply_eq]
    rw [if_neg (of_decide_eq_true hi), if_neg (of_decide_eq_true hj)]

/-- The diagonal unit at the ambient basis index. -/
def ambientUnit : Matrix (Fin 33) (Fin 33) ℂ :=
  Matrix.diagonal fun i => if i = 0 then 1 else 0

/-- The ambient unit is the identity minus the four support
projectors: the supports partition the basis away from the ambient
index. -/
theorem ambientUnit_eq :
    ambientUnit = 1 - (supportProjector .o36 + supportProjector .o60 +
      supportProjector .o64 + supportProjector .o92) := by
  unfold ambientUnit supportProjector
  rw [← Matrix.diagonal_one, Matrix.diagonal_add, Matrix.diagonal_add,
    Matrix.diagonal_add, Matrix.diagonal_sub]
  congr 1
  funext i
  rcases h : nodeOwner i with _ | o
  · have h0 : i = 0 := (nodeOwner_none_iff i).mp h
    subst h0
    simp
  · have h0 : i ≠ 0 := by
      intro he
      rw [(nodeOwner_none_iff i).mpr he] at h
      simp at h
    cases o <;> simp [h0]

/-- The truncation of a matrix onto one observer support: entries on
related pairs led by the support, zero elsewhere. -/
def supportTruncation (o : PayloadObserver)
    (X : Matrix (Fin 33) (Fin 33) ℂ) : Matrix (Fin 33) (Fin 33) ℂ :=
  fun i j => if nodeOwner i = some o ∧ enrichedRelTop i j = true then
    X i j else 0

theorem supportTruncation_mem_span (o : PayloadObserver)
    (ho : o ≠ .o64) (X : Matrix (Fin 33) (Fin 33) ℂ) :
    supportTruncation o X ∈ diagonalClassSubalgebra (fReg o) := by
  refine ⟨fun x => x.elim 0 fun k => X k k, ?_⟩
  ext i j
  by_cases hij : i = j
  · subst hij
    rw [Matrix.diagonal_apply_eq]
    by_cases h : nodeOwner i = some o
    · simp only [supportTruncation]
      rw [if_pos ⟨h, topBlockRelation.rel_refl i⟩]
      simp [fReg, h]
    · simp only [supportTruncation]
      rw [if_neg fun hc => h hc.1]
      simp [fReg, h]
  · rw [Matrix.diagonal_apply_ne _ hij]
    symm
    simp only [supportTruncation]
    rw [if_neg]
    rintro ⟨h1, h2⟩
    obtain ⟨ho64, _⟩ := enrichedRelTop_offdiag_owner i j h2 hij
    rw [h1] at ho64
    exact ho (Option.some.inj ho64)

theorem supportTruncation_mem_s64 (X : Matrix (Fin 33) (Fin 33) ℂ) :
    supportTruncation .o64 X ∈ sourceRegionalS64 := by
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · have h0' : enrichedRel .o64 i j = false := h0
    simp only [supportTruncation]
    rw [if_neg]
    rintro ⟨_, h2⟩
    rw [enrichedRelTop_eq_o64 i j] at h2
    exact Bool.noConfusion (h0'.symm.trans h2)
  · simp only [supportTruncation]
    rw [if_neg fun hc => (of_decide_eq_true hi) hc.1,
      if_neg fun hc => (of_decide_eq_true hj) hc.1]

/-- Exact decomposition of an enriched top observable into the four
support truncations and the ambient scalar term. -/
theorem top_decomposition {X : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ sourceRegionalTop) :
    X = supportTruncation .o36 X + supportTruncation .o60 X +
      supportTruncation .o64 X + supportTruncation .o92 X +
        X 0 0 • ambientUnit := by
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, supportTruncation,
    ambientUnit]
  by_cases hij : i = j
  · subst hij
    rw [Matrix.diagonal_apply_eq]
    rcases h : nodeOwner i with _ | o
    · have h0 : i = 0 := (nodeOwner_none_iff i).mp h
      subst h0
      simp
    · have h0 : i ≠ 0 := by
        intro he
        rw [(nodeOwner_none_iff i).mpr he] at h
        simp at h
      have hrefl : enrichedRelTop i i = true := topBlockRelation.rel_refl i
      cases o <;> simp [h0, hrefl]
  · rw [Matrix.diagonal_apply_ne _ hij]
    by_cases hrel : enrichedRelTop i j = true
    · obtain ⟨hoi, _⟩ := enrichedRelTop_offdiag_owner i j hrel hij
      simp [hoi, hrel]
    · have h0 : X i j = 0 := hX.1 i j (bool_eq_false hrel)
      simp [hrel, h0]

/-- The four-support family under the top region. -/
noncomputable def sourceRegionalCover :
    sourceRegionalNet.FiniteCover () SimRegion.top where
  regions := ({.s36, .s60, .s64, .s92} : Finset SimRegion)
  nonempty := by
    refine ⟨SimRegion.s36, ?_⟩
    change SimRegion.s36 ∈ ({.s36, .s60, .s64, .s92} : Finset SimRegion)
    simp
  subregion := fun U _ => by cases U <;> trivial

/-- Every enriched top observable is generated by the four support
algebras: the support truncations lie in the member algebras, and the
ambient unit is the identity minus the four support projectors. -/
theorem top_mem_familyJoin {X : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ sourceRegionalTop) :
    X ∈ sourceRegionalNet.familyJoin sourceRegionalCover := by
  have h36 : SimRegion.s36 ∈ sourceRegionalCover.regions := by
    show SimRegion.s36 ∈ ({.s36, .s60, .s64, .s92} : Finset SimRegion)
    simp
  have h60 : SimRegion.s60 ∈ sourceRegionalCover.regions := by
    show SimRegion.s60 ∈ ({.s36, .s60, .s64, .s92} : Finset SimRegion)
    simp
  have h64 : SimRegion.s64 ∈ sourceRegionalCover.regions := by
    show SimRegion.s64 ∈ ({.s36, .s60, .s64, .s92} : Finset SimRegion)
    simp
  have h92 : SimRegion.s92 ∈ sourceRegionalCover.regions := by
    show SimRegion.s92 ∈ ({.s36, .s60, .s64, .s92} : Finset SimRegion)
    simp
  have hT36 := sourceRegionalNet.le_familyJoin sourceRegionalCover h36
    (supportTruncation_mem_span .o36 (by decide) X)
  have hT60 := sourceRegionalNet.le_familyJoin sourceRegionalCover h60
    (supportTruncation_mem_span .o60 (by decide) X)
  have hT64 := sourceRegionalNet.le_familyJoin sourceRegionalCover h64
    (supportTruncation_mem_s64 X)
  have hT92 := sourceRegionalNet.le_familyJoin sourceRegionalCover h92
    (supportTruncation_mem_span .o92 (by decide) X)
  have hP36 := sourceRegionalNet.le_familyJoin sourceRegionalCover h36
    (supportProjector_mem_span .o36)
  have hP60 := sourceRegionalNet.le_familyJoin sourceRegionalCover h60
    (supportProjector_mem_span .o60)
  have hP64 := sourceRegionalNet.le_familyJoin sourceRegionalCover h64
    supportProjector_mem_s64
  have hP92 := sourceRegionalNet.le_familyJoin sourceRegionalCover h92
    (supportProjector_mem_span .o92)
  have hAmb : ambientUnit ∈
      sourceRegionalNet.familyJoin sourceRegionalCover := by
    rw [ambientUnit_eq]
    exact sub_mem (one_mem _)
      (add_mem (add_mem (add_mem hP36 hP60) hP64) hP92)
  rw [top_decomposition hX]
  refine add_mem (add_mem (add_mem (add_mem hT36 hT60) hT64) hT92) ?_
  have hsmul : X 0 0 • ambientUnit =
      algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) (X 0 0) * ambientUnit :=
    Algebra.smul_def _ _
  rw [hsmul]
  exact mul_mem
    ((sourceRegionalNet.familyJoin sourceRegionalCover).algebraMap_mem
      (X 0 0)) hAmb

/-- **The coverage law, computed** through the machinery of
`QFT/CoverageReceipt.lean`. -/
theorem sourceRegionalCover_coverageLaw :
    sourceRegionalNet.CoverageLaw sourceRegionalCover :=
  ⟨fun _X hX => top_mem_familyJoin hX⟩

/-- The family join of the four supports is computed: it equals the
enriched top algebra. -/
theorem sourceRegionalCover_familyJoin :
    sourceRegionalNet.familyJoin sourceRegionalCover =
      sourceRegionalNet.localAlgebra () SimRegion.top :=
  sourceRegionalCover_coverageLaw.join_eq

/-- The three-support family without the enriched region. -/
noncomputable def sourceRegionalDiagonalCover :
    sourceRegionalNet.FiniteCover () SimRegion.top where
  regions := ({.s36, .s60, .s92} : Finset SimRegion)
  nonempty := by
    refine ⟨SimRegion.s36, ?_⟩
    change SimRegion.s36 ∈ ({.s36, .s60, .s92} : Finset SimRegion)
    simp
  subregion := fun U _ => by cases U <;> trivial

/-- The raising unit of the designated fibre lies in the enriched top
algebra. -/
theorem fibreUnit_17_20_mem_top : fibreUnit 17 20 ∈ sourceRegionalTop := by
  refine ⟨fun i j h0 => ?_, fun i j hi _ => Bool.noConfusion hi⟩
  by_cases hij : i = 17 ∧ j = 20
  · obtain ⟨rfl, rfl⟩ := hij
    have : enrichedRelTop 17 20 = true := by decide
    exact absurd this fun htrue => Bool.noConfusion (h0.symm.trans htrue)
  · simp [fibreUnit, hij]

/-- **Negative coverage control.**  The family without the `s64`
support fails the coverage law: its join is diagonal, and the top
algebra carries the raising unit of the designated fibre.  The dropped
region carries the entire noncommutative sector of the net. -/
theorem sourceRegionalDiagonalCover_not_coverageLaw :
    ¬ sourceRegionalNet.CoverageLaw sourceRegionalDiagonalCover := by
  intro hcov
  have hjoin : sourceRegionalNet.familyJoin sourceRegionalDiagonalCover ≤
      diagonalClassSubalgebra fTop := by
    refine sourceRegionalNet.familyJoin_le fun U hU => ?_
    have hU' : U ∈ ({.s36, .s60, .s92} : Finset SimRegion) := hU
    cases U
    case bot => exact absurd hU' (by decide)
    case s36 => exact diagonalClassSubalgebra_comp_le (fReg .o36) fTop
    case s60 => exact diagonalClassSubalgebra_comp_le (fReg .o60) fTop
    case s64 => exact absurd hU' (by decide)
    case s92 => exact diagonalClassSubalgebra_comp_le (fReg .o92) fTop
    case top => exact absurd hU' (by decide)
  have hmem : fibreUnit 17 20 ∈ diagonalClassSubalgebra fTop :=
    hjoin (hcov.covered_le_join fibreUnit_17_20_mem_top)
  obtain ⟨c, hc⟩ := hmem
  have h1720 := congrFun (congrFun hc 17) 20
  rw [Matrix.diagonal_apply_ne _ (by decide : (17 : Fin 33) ≠ 20)] at h1720
  simp [fibreUnit] at h1720

/-! ## Enlargement boundary on the net -/

/-- The regional algebra of each collapsed support is the block
algebra of the uniform enrichment rule: the net hides nothing, and the
collapse receipts identify these algebras with the diagonal spans. -/
theorem sourceRegionalNet_localAlgebra_collapse :
    sourceRegionalNet.localAlgebra () (SimRegion.regionOf .o36) =
        o36BlockRelation.subalgebra ∧
      sourceRegionalNet.localAlgebra () (SimRegion.regionOf .o60) =
        o60BlockRelation.subalgebra ∧
      sourceRegionalNet.localAlgebra () (SimRegion.regionOf .o92) =
        o92BlockRelation.subalgebra :=
  ⟨enrichedRegion_collapse_o36.symm, enrichedRegion_collapse_o60.symm,
    enrichedRegion_collapse_o92.symm⟩

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.splitRecord_o64_char
#print axioms OPH.QFT.truncatedSplitFibreCounts
#print axioms OPH.QFT.splitFibre18_receipt
#print axioms OPH.QFT.splitFibre28_receipt
#print axioms OPH.QFT.designatedBlock_receipt
#print axioms OPH.QFT.splitFibre18_basis
#print axioms OPH.QFT.splitFibre28_basis
#print axioms OPH.QFT.enrichedRegion_collapse_o36
#print axioms OPH.QFT.enrichedRegion_collapse_o60
#print axioms OPH.QFT.enrichedRegion_collapse_o92
#print axioms OPH.QFT.m2Block18_mem
#print axioms OPH.QFT.m2Block28_mem
#print axioms OPH.QFT.fibreUnit_not_mem_s64
#print axioms OPH.QFT.sourceRegionalS64_noncommutative
#print axioms OPH.QFT.diagonal_lt_sourceRegionalS64
#print axioms OPH.QFT.sourceRegionalNet_disjoint_supports
#print axioms OPH.QFT.sourceRegionalNet_separation
#print axioms OPH.QFT.sourceRegionalNet_noncommutative_iff
#print axioms OPH.QFT.top_mem_familyJoin
#print axioms OPH.QFT.sourceRegionalCover_coverageLaw
#print axioms OPH.QFT.sourceRegionalCover_familyJoin
#print axioms OPH.QFT.sourceRegionalDiagonalCover_not_coverageLaw
#print axioms OPH.QFT.enrichedRel_of_no_split_fibre
#print axioms OPH.QFT.sourceRegionalNet_localAlgebra_collapse
