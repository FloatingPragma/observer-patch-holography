import QFT.ObserverAccessCut

/-!
# Sim-earned access-cut witness from the B12 preregistered run

This module mirrors, as Lean literals, the exact finite payload
`docs/SIM_EARNED_WITNESS_PAYLOAD.json` of the `oph-physics-sim`
repository, and builds from those literals an earned instance of the E6
access-cut interface.  Interpretation of every literal is governed by the
payload's sibling `docs/SIM_EARNED_WITNESS_NOTES.md`.

## Provenance (verbatim from the payload's provenance block)

* `schema`: `oph.sim.lean_witness_payload.v1`
* `run_id`: `b12_prereg_16k_20260806`
* `run_git_commit`: `b39b78faf894894ebe573571e0902ccfaaeac32a`
* `seed`: `20260806`
* `receipt_path`:
  `runs/b12_prereg_16k_20260806/conditional_resampling_realization_receipt.json`
* `receipt_sha256`:
  `d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749`
* `receipt_schema`: `oph.sim.conditional_resampling_realization.v1`
* `freezeout_path`: `runs/b12_prereg_16k_20260806/freezeout_fields.npz`
* `freezeout_sha256`:
  `b962c5b80205a17d5d6bc023f5f5d487bc23c1327793978e7d1bc69494fac49e`
* `observer_views_path`: `runs/b12_prereg_16k_20260806/observer_views.jsonl`
* `observer_views_sha256`:
  `3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca`
* extraction: `script = scripts/extract_lean_witness_payload.py`,
  `binning_producer =
  oph_fpe.dynamics.conditional_resampling.realization_inputs_from_freezeout`,
  `record_field = record_signature`,
  `companion_field = cumulative_repair_load`,
  `max_record_classes = 32`, `max_companion_classes = 16`,
  `observer_selection_rule = first 4 patch_observer rows of
  observer_views.jsonl in file order`, `support_truncation = 8`.

## Declared mapping from payload fields to Lean fields

* Each of the four `patch_observer` rows becomes one constructor of
  `PayloadObserver`; `observer_id` is the value of `payloadObserverId`.
* `support_nodes_first8` becomes `supportNodes` (a `List ℕ` literal),
  `support_nodes_first8_record_classes` becomes `recordClassAt`
  (values in `Fin 32`), and `support_nodes_first8_companion_classes`
  becomes `companionClassAt` (values in `Fin 16`).
* `realized_record_class_counts` becomes `realizedRecordClassCounts`,
  `support_record_class_counts` becomes `supportRecordClassCounts`, and
  `selected_support_joint_table.cells` becomes `jointTable`, all as exact
  integer list literals.
* The finite basis `Fin 33` lays out one ambient index `0` for the screen
  outside the four truncated supports and one index per truncated support
  node (`nodeIndex`); the ambient index is declared padding required by
  the restriction system and carries no payload datum.
* The tower `simEarnedTower` has `Record () = Fin 32`, the payload's
  record-class type, and `Observer () = PayloadObserver`.  The public
  algebra of an observer is the diagonal span of its record classes on
  its truncated support; the accessible algebra of `simEarnedCut` is the
  finer diagonal span of its joint (record, companion) classes; the
  regional algebra of its support region is the full diagonal span on the
  support.  The region type `SimRegion` has one region per observer
  support plus a top; the meet closure of that specification forces one
  bottom region, which carries the scalar algebra.
* The tower state is the declared pure state on the ambient index and the
  repair maps are identities: the payload pins no state and no repair
  action, so these fields are declared data, not mirrored data.

## Named gaps

Where the truncated payload data cannot support a receipt, the collapse
is stated as a theorem instead of being repaired by synthetic data:

* `strictnessGap_public_eq_accessible_o36` (and `_o60`, `_o92`): on these
  truncated supports every record class carries a single companion class,
  so the public and accessible algebras coincide and the strict
  interposition receipt is unavailable.  Observer 64 is the one observer
  whose truncated support separates a record class by companion class,
  and `simEarnedCut_public_lt_accessible_o64` records that strictness.
* `regionalGap_accessible_eq_regional_o60` (and `_o64`, `_o92`): on these
  truncated supports the joint classes are pairwise distinct, so the
  accessible algebra equals the full regional diagonal algebra.  Observer
  36 is the one observer with a repeated joint class, and
  `simEarnedCut_accessible_lt_regional_o36` records that strictness.
* The joint table aggregates the full 96-node supports (245 union
  patches), while the per-observer literals are truncated to 8 nodes, so
  no marginal-agreement receipt between `jointTable` and the truncated
  node classes is stated; `jointTable_covers_truncated_pairs` states the
  weaker containment that every truncated (record, companion) pair
  occupies a cell of the table.

## Claim boundary

The literals are earned finite data from the preregistered bounded run
identified above and carry no physical claim.  Verbatim from the
payload: "Realized finite data from one committed run under a declared
binning and a declared observer selection rule.  The integers are counts
and class labels on this run's freezeout fields; no physical claim, no
probability claim beyond realized frequencies, and no canonical status
for the selection are asserted."  Mirroring the integers in Lean carries
no claim that the mirrored structure is forced by the run; the tower,
net, and cut built here are one declared reading of the literals through
the committed E1/E6 interfaces.
-/

namespace OPH.QFT

open OPH.Tower
open EventAlgebra
open Matrix

/- The joint-table literal is a 145-cell list; kernel evaluation of the
`decide` receipts over it needs a larger recursion budget than the
default. -/
set_option maxRecDepth 8192

/-! ## The four payload observers as literals -/

/-- The four selected `patch_observer` rows of the payload, named by their
`observer_id` values. -/
inductive PayloadObserver : Type
  | o36
  | o60
  | o64
  | o92
  deriving DecidableEq, Fintype

namespace PayloadObserver

/-- Payload field `observer_id`. -/
def payloadObserverId : PayloadObserver → ℕ
  | .o36 => 36
  | .o60 => 60
  | .o64 => 64
  | .o92 => 92

/-- Payload field `support_patch_count`: the full support size before the
declared truncation to 8 nodes. -/
def supportPatchCount : PayloadObserver → ℕ := fun _ => 96

/-- Payload field `support_nodes_first8`: the truncated support-node list
of each observer. -/
def supportNodes : PayloadObserver → List ℕ
  | .o36 => [36, 10, 15, 23, 28, 31, 41, 44]
  | .o60 => [60, 26, 34, 39, 47, 52, 68, 73]
  | .o64 => [64, 30, 38, 43, 51, 56, 72, 77]
  | .o92 => [92, 37, 50, 58, 71, 79, 84, 100]

/-- Payload field `support_nodes_first8_record_classes`: the record class
of each truncated support node, as a `Fin 32` value. -/
def recordClassAt : PayloadObserver → Fin 8 → Fin 32
  | .o36 => ![2, 23, 20, 20, 16, 22, 19, 7]
  | .o60 => ![29, 20, 23, 1, 4, 11, 21, 30]
  | .o64 => ![18, 2, 8, 18, 28, 24, 15, 28]
  | .o92 => ![21, 20, 1, 30, 18, 7, 6, 27]

/-- Payload field `support_nodes_first8_companion_classes`: the companion
class of each truncated support node, as a `Fin 16` value. -/
def companionClassAt : PayloadObserver → Fin 8 → Fin 16
  | .o36 => ![3, 3, 6, 6, 1, 3, 10, 7]
  | .o60 => ![3, 13, 10, 6, 10, 7, 10, 7]
  | .o64 => ![13, 13, 15, 3, 13, 13, 3, 6]
  | .o92 => ![1, 13, 10, 0, 15, 6, 10, 3]

/-- Payload field `support_record_class_counts`: the record-class counts
over the full 96-node support of each observer. -/
def supportRecordClassCounts : PayloadObserver → List (ℕ × ℕ)
  | .o36 => [(0, 2), (1, 5), (2, 4), (3, 3), (4, 2), (5, 2), (6, 2),
      (7, 5), (8, 2), (9, 1), (10, 3), (11, 3), (12, 3), (13, 2), (14, 3),
      (15, 1), (16, 3), (17, 3), (18, 3), (19, 4), (20, 5), (21, 3),
      (22, 2), (23, 5), (24, 2), (25, 3), (26, 4), (28, 4), (29, 7),
      (30, 4), (31, 1)]
  | .o60 => [(0, 2), (1, 9), (2, 3), (3, 1), (4, 3), (5, 1), (6, 4),
      (7, 1), (8, 2), (9, 2), (10, 1), (11, 5), (12, 2), (13, 1), (14, 4),
      (15, 2), (16, 3), (17, 1), (18, 4), (19, 3), (20, 7), (21, 4),
      (22, 6), (23, 3), (24, 2), (25, 1), (26, 4), (27, 1), (28, 2),
      (29, 3), (30, 7), (31, 2)]
  | .o64 => [(0, 2), (1, 5), (2, 4), (3, 4), (4, 3), (5, 2), (6, 4),
      (7, 1), (8, 3), (9, 2), (10, 3), (11, 2), (12, 5), (13, 1), (14, 3),
      (15, 2), (16, 4), (17, 4), (18, 5), (19, 5), (20, 4), (21, 3),
      (22, 1), (23, 3), (24, 5), (25, 3), (26, 1), (28, 3), (29, 5),
      (30, 2), (31, 2)]
  | .o92 => [(0, 1), (1, 4), (2, 3), (3, 3), (4, 2), (5, 3), (6, 4),
      (7, 7), (8, 2), (9, 3), (10, 3), (11, 3), (12, 2), (13, 2), (14, 1),
      (15, 4), (17, 2), (18, 4), (19, 1), (20, 5), (21, 6), (22, 3),
      (23, 5), (24, 2), (25, 3), (26, 3), (27, 3), (28, 3), (29, 2),
      (30, 5), (31, 2)]

end PayloadObserver

/-- Payload field `realized_record_class_counts`: the realized
record-class counts over the whole 16384-patch screen. -/
def realizedRecordClassCounts : List (ℕ × ℕ) :=
  [(0, 512), (1, 512), (2, 512), (3, 512), (4, 512), (5, 512), (6, 512),
   (7, 512), (8, 512), (9, 512), (10, 512), (11, 512), (12, 512),
   (13, 512), (14, 512), (15, 512), (16, 512), (17, 512), (18, 512),
   (19, 512), (20, 512), (21, 512), (22, 512), (23, 512), (24, 512),
   (25, 512), (26, 512), (27, 512), (28, 512), (29, 512), (30, 512),
   (31, 512)]

/-- Payload field `companion_class_labels_realized`: the eight companion
labels realized on the screen. -/
def companionLabelsRealized : List ℕ := [0, 1, 3, 6, 7, 10, 13, 15]

/-- Payload field `selected_support_joint_table.cells`: the exact joint
(record class, companion class, count) table over the union of the four
full supports, 145 occupied cells of total mass 245. -/
def jointTable : List (ℕ × ℕ × ℕ) :=
  [(0, 0, 1), (0, 10, 3), (1, 3, 4), (1, 6, 3), (1, 7, 3), (1, 10, 1),
   (1, 13, 2), (1, 15, 1), (2, 3, 2), (2, 6, 1), (2, 13, 4), (3, 3, 2),
   (3, 7, 1), (3, 10, 1), (3, 13, 1), (3, 15, 1), (4, 10, 2), (4, 13, 2),
   (4, 15, 1), (5, 6, 2), (5, 7, 1), (5, 15, 3), (6, 1, 1), (6, 3, 1),
   (6, 6, 1), (6, 10, 2), (6, 13, 2), (6, 15, 1), (7, 3, 2), (7, 6, 2),
   (7, 7, 1), (7, 13, 3), (7, 15, 4), (8, 10, 2), (8, 13, 2), (8, 15, 1),
   (9, 3, 1), (9, 7, 1), (9, 13, 2), (9, 15, 1), (10, 3, 1), (10, 10, 3),
   (10, 13, 1), (10, 15, 2), (11, 1, 1), (11, 3, 1), (11, 6, 1),
   (11, 7, 1), (11, 10, 4), (12, 3, 1), (12, 10, 6), (12, 13, 1),
   (12, 15, 1), (13, 6, 1), (13, 10, 1), (14, 3, 4), (14, 10, 1),
   (14, 13, 1), (14, 15, 1), (15, 1, 1), (15, 3, 2), (15, 7, 1),
   (15, 10, 1), (15, 13, 1), (15, 15, 1), (16, 1, 1), (16, 3, 1),
   (16, 13, 3), (16, 15, 1), (17, 3, 1), (17, 7, 2), (17, 13, 2),
   (17, 15, 1), (18, 3, 2), (18, 6, 1), (18, 7, 1), (18, 10, 3),
   (18, 13, 1), (18, 15, 2), (19, 1, 1), (19, 6, 4), (19, 10, 2),
   (19, 13, 2), (20, 0, 1), (20, 6, 3), (20, 7, 1), (20, 13, 4),
   (20, 15, 4), (21, 1, 1), (21, 3, 3), (21, 6, 1), (21, 10, 2),
   (21, 13, 2), (21, 15, 1), (22, 3, 2), (22, 6, 2), (22, 7, 1),
   (22, 13, 1), (22, 15, 2), (23, 3, 4), (23, 6, 2), (23, 10, 2),
   (23, 13, 1), (23, 15, 2), (24, 3, 1), (24, 6, 1), (24, 7, 1),
   (24, 10, 1), (24, 13, 1), (24, 15, 2), (25, 3, 1), (25, 7, 2),
   (25, 10, 1), (25, 13, 1), (25, 15, 2), (26, 1, 1), (26, 3, 1),
   (26, 6, 2), (26, 7, 2), (26, 13, 2), (27, 3, 1), (27, 6, 2),
   (28, 3, 2), (28, 6, 1), (28, 10, 3), (28, 13, 3), (28, 15, 1),
   (29, 1, 1), (29, 3, 1), (29, 6, 1), (29, 7, 1), (29, 10, 4),
   (29, 13, 2), (29, 15, 1), (30, 0, 1), (30, 3, 2), (30, 6, 1),
   (30, 7, 1), (30, 10, 1), (30, 15, 2), (31, 3, 1), (31, 6, 1),
   (31, 7, 2), (31, 13, 1), (31, 15, 1)]

/-! ## Literal receipts

Each statement below is decided by kernel computation on the literals
alone; the payload's `class_structure` block is the pinned reference the
extraction script checked these totals against. -/

/-- The screen realizes all 32 record classes, in label order. -/
theorem realizedRecordClassCounts_labels :
    realizedRecordClassCounts.map Prod.fst = List.range 32 := by decide

/-- Every realized record class has exactly 512 patches. -/
theorem realizedRecordClassCounts_uniform :
    ∀ p ∈ realizedRecordClassCounts, p.2 = 512 := by decide

/-- The realized record-class counts sum to the total patch count 16384. -/
theorem realizedRecordClassCounts_total :
    (realizedRecordClassCounts.map Prod.snd).sum = 16384 := by decide

/-- Each truncated support list has exactly 8 nodes. -/
theorem supportNodes_length :
    ∀ o : PayloadObserver, (PayloadObserver.supportNodes o).length = 8 := by
  decide

/-- The truncation receipt: each 8-node literal list is a proper
truncation of the observer's 96-patch support. -/
theorem supportNodes_truncated :
    ∀ o : PayloadObserver,
      (PayloadObserver.supportNodes o).length <
        PayloadObserver.supportPatchCount o := by
  decide

/-- The four truncated support lists are pairwise disjoint: no payload
node occurs in two of them.  This literal fact backs the declared
disjointness relation of the net below. -/
theorem supportNodes_pairwise_disjoint :
    ∀ o o' : PayloadObserver, o ≠ o' →
      ∀ x ∈ PayloadObserver.supportNodes o,
        x ∉ PayloadObserver.supportNodes o' := by
  decide

/-- The full-support record-class counts of each observer sum to its
96-patch support size. -/
theorem supportRecordClassCounts_total :
    ∀ o : PayloadObserver,
      ((PayloadObserver.supportRecordClassCounts o).map Prod.snd).sum =
        PayloadObserver.supportPatchCount o := by
  decide

/-- Every companion class realized on a truncated support node is one of
the eight companion labels realized on the screen. -/
theorem companionClasses_realized :
    ∀ (o : PayloadObserver) (j : Fin 8),
      (PayloadObserver.companionClassAt o j : ℕ) ∈
        companionLabelsRealized := by
  decide

/-- The joint table has 145 occupied cells. -/
theorem jointTable_cells : jointTable.length = 145 := by decide

/-- The joint table's total mass is the 245-patch support union. -/
theorem jointTable_mass :
    (jointTable.map (fun c => c.2.2)).sum = 245 := by decide

/-- Every companion label occupied in the joint table is realized on the
screen. -/
theorem jointTable_companions_realized :
    ∀ c ∈ jointTable, c.2.1 ∈ companionLabelsRealized := by decide

/-- Every truncated (record class, companion class) pair of every
observer occupies a cell of the joint table.  This is the containment
direction available at truncation 8; marginal agreement is a named gap
recorded in the module docstring. -/
theorem jointTable_covers_truncated_pairs :
    ∀ (o : PayloadObserver) (j : Fin 8), ∃ c ∈ jointTable,
      c.1 = (PayloadObserver.recordClassAt o j : ℕ) ∧
        c.2.1 = (PayloadObserver.companionClassAt o j : ℕ) := by
  decide

/-- Nonconstancy of the record classes, from the literal counts: two
distinct record classes both have nonzero realized count. -/
theorem realized_record_classes_nonconstant :
    ∃ p ∈ realizedRecordClassCounts, ∃ q ∈ realizedRecordClassCounts,
      p.1 ≠ q.1 ∧ 0 < p.2 ∧ 0 < q.2 := by
  exact ⟨(0, 512), by decide, (1, 512), by decide, by decide⟩

/-- Nonconstancy of the record classes on every truncated support: each
observer reads at least two distinct record classes. -/
theorem truncated_record_classes_nonconstant :
    ∀ o : PayloadObserver, ∃ j : Fin 8,
      PayloadObserver.recordClassAt o j ≠
        PayloadObserver.recordClassAt o 0 := by
  decide

/-! ## The finite basis layout

`Fin 33` carries one ambient index `0` for the screen outside the four
truncated supports and one index per truncated support node, in observer
order.  The ambient index is declared padding: the restriction system of
the net needs one basis point owned by no observer. -/

/-- The block position of each observer in the basis layout. -/
def obsIdx : PayloadObserver → Fin 4
  | .o36 => 0
  | .o60 => 1
  | .o64 => 2
  | .o92 => 3

/-- The basis index of the `j`-th truncated support node of observer
`o`: index `0` is ambient and the four 8-node blocks follow. -/
def nodeIndex (o : PayloadObserver) (j : Fin 8) : Fin 33 :=
  ⟨1 + 8 * (obsIdx o).val + j.val, by
    have h1 := (obsIdx o).isLt
    have h2 := j.isLt
    omega⟩

/-- The owning observer of a basis index; the ambient index has none. -/
def nodeOwner : Fin 33 → Option PayloadObserver := fun i =>
  if i.val = 0 then none
  else if i.val ≤ 8 then some .o36
  else if i.val ≤ 16 then some .o60
  else if i.val ≤ 24 then some .o64
  else some .o92

/-- The in-block position of a basis index. -/
def nodePos : Fin 33 → Fin 8 := fun i =>
  ⟨(i.val - 1) % 8, Nat.mod_lt _ (by norm_num)⟩

/-- The record class read at a basis index; the ambient index carries the
junk value `0`, which no construction below reads. -/
def nodeRecordClass : Fin 33 → Fin 32 := fun i =>
  match nodeOwner i with
  | none => 0
  | some o => PayloadObserver.recordClassAt o (nodePos i)

/-- The companion class read at a basis index; the ambient value is junk
in the same sense as `nodeRecordClass`. -/
def nodeCompanionClass : Fin 33 → Fin 16 := fun i =>
  match nodeOwner i with
  | none => 0
  | some o => PayloadObserver.companionClassAt o (nodePos i)

/-- The layout is coherent: `nodeIndex` lands in the block of its
observer at its position. -/
theorem nodeIndex_coherent : ∀ (o : PayloadObserver) (j : Fin 8),
    nodeOwner (nodeIndex o j) = some o ∧ nodePos (nodeIndex o j) = j := by
  decide

/-- The basis reads back exactly the payload record-class literals. -/
theorem nodeRecordClass_nodeIndex : ∀ (o : PayloadObserver) (j : Fin 8),
    nodeRecordClass (nodeIndex o j) = PayloadObserver.recordClassAt o j := by
  decide

/-- The basis reads back exactly the payload companion-class literals. -/
theorem nodeCompanionClass_nodeIndex : ∀ (o : PayloadObserver) (j : Fin 8),
    nodeCompanionClass (nodeIndex o j) =
      PayloadObserver.companionClassAt o j := by
  decide

/-! ## Class functions

Each algebra below is the diagonal span induced by one class function on
the basis.  The regional function remembers the node, the accessible
function remembers the joint (record, companion) class, and the public
function remembers the record class alone; off the observer's support all
three collapse to a single class. -/

/-- The regional class function of an observer: the node itself on the
observer's support, one collapsed class elsewhere. -/
def fReg (o : PayloadObserver) : Fin 33 → Option (Fin 33) := fun i =>
  if nodeOwner i = some o then some i else none

/-- The accessible class function of an observer: the joint
(record class, companion class) pair on the support. -/
def fAcc (o : PayloadObserver) : Fin 33 → Option (Fin 32 × Fin 16) :=
  fun i => (fReg o i).map fun j => (nodeRecordClass j, nodeCompanionClass j)

/-- The public class function of an observer: the record class alone. -/
def fPub (o : PayloadObserver) : Fin 33 → Option (Fin 32) := fun i =>
  (fAcc o i).map Prod.fst

/-- The class function of the top region: every node is its own class. -/
def fTop : Fin 33 → Fin 33 := fun i => i

/-- The class function of the bottom region: one class. -/
def fBot : Fin 33 → Unit := fun _ => ()

/-- The collapse of the basis onto an observer's support: fixed on the
support, the ambient index elsewhere. -/
def collapseMap (o : PayloadObserver) : Fin 33 → Fin 33 := fun i =>
  if nodeOwner i = some o then i else 0

theorem nodeOwner_zero : nodeOwner 0 = none := rfl

theorem collapseMap_zero (o : PayloadObserver) : collapseMap o 0 = 0 := by
  simp [collapseMap, nodeOwner_zero]

/-- The regional class function is blind to the collapse. -/
theorem fReg_collapseMap (o : PayloadObserver) (i : Fin 33) :
    fReg o (collapseMap o i) = fReg o i := by
  by_cases h : nodeOwner i = some o
  · simp [collapseMap, h]
  · simp [collapseMap, h, fReg, nodeOwner_zero]

/-! ## Diagonal spans induced by class functions -/

/-- The diagonal of a constant vector is the image of the scalar under
the algebra map. -/
theorem diagonal_const_eq_algebraMap (z : ℂ) :
    Matrix.diagonal (fun _ : Fin 33 => z) =
      algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) z := by
  rw [Matrix.algebraMap_eq_diagonal]
  rfl

/-- The diagonal span induced by a class function `f`: all matrices
`diagonal (c ∘ f)` for class-indexed coefficient vectors `c`.  Two basis
indices are inseparable by this algebra exactly when `f` identifies
them. -/
def diagonalClassSubalgebra {α : Type} (f : Fin 33 → α) :
    StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ) where
  carrier := {X | ∃ c : α → ℂ, Matrix.diagonal (fun i => c (f i)) = X}
  zero_mem' := ⟨fun _ => 0, by simp⟩
  one_mem' := ⟨fun _ => 1, by simp⟩
  add_mem' := by
    rintro X Y ⟨c, rfl⟩ ⟨d, rfl⟩
    exact ⟨fun a => c a + d a, by simp [Matrix.diagonal_add]⟩
  mul_mem' := by
    rintro X Y ⟨c, rfl⟩ ⟨d, rfl⟩
    exact ⟨fun a => c a * d a, by simp [Matrix.diagonal_mul_diagonal]⟩
  algebraMap_mem' := fun z => ⟨fun _ => z, diagonal_const_eq_algebraMap z⟩
  star_mem' := by
    rintro X ⟨c, rfl⟩
    exact ⟨fun a => star (c a), by
      simp [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]⟩

theorem mem_diagonalClassSubalgebra_iff {α : Type} {f : Fin 33 → α}
    {X : Matrix (Fin 33) (Fin 33) ℂ} :
    X ∈ diagonalClassSubalgebra f ↔
      ∃ c : α → ℂ, Matrix.diagonal (fun i => c (f i)) = X :=
  Iff.rfl

/-- Coarsening: post-composing the class function shrinks the span. -/
theorem diagonalClassSubalgebra_comp_le {α β : Type} (h : α → β)
    (f : Fin 33 → α) :
    diagonalClassSubalgebra (fun i => h (f i)) ≤
      diagonalClassSubalgebra f := by
  rintro X ⟨c, rfl⟩
  exact ⟨fun a => c (h a), rfl⟩

/-- Every diagonal span is commutative. -/
theorem diagonalClassSubalgebra_mul_comm {α : Type} (f : Fin 33 → α)
    (X Y : diagonalClassSubalgebra f) : X * Y = Y * X := by
  apply Subtype.ext
  obtain ⟨c, hc⟩ := X.2
  obtain ⟨d, hd⟩ := Y.2
  show (X : Matrix (Fin 33) (Fin 33) ℂ) * Y = (Y : Matrix (Fin 33) (Fin 33) ℂ) * X
  rw [← hc, ← hd, Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm _ _

/-- The raising matrix unit on the first two basis indices. -/
def simRaise : Matrix (Fin 33) (Fin 33) ℂ :=
  fun k l => if k = 0 ∧ l = 1 then 1 else 0

/-- No diagonal span contains the raising unit. -/
theorem simRaise_not_mem {α : Type} (f : Fin 33 → α) :
    simRaise ∉ diagonalClassSubalgebra f := by
  rintro ⟨c, hc⟩
  have h01 := congrFun (congrFun hc 0) 1
  rw [Matrix.diagonal_apply_ne _ (by decide : (0 : Fin 33) ≠ 1)] at h01
  simp [simRaise] at h01

/-! ## The region lattice

One region per observer support plus a top; the meet closure of that
specification forces one bottom region, which represents the undeclared
pairwise intersections and carries the scalar algebra. -/

/-- Region labels: the bottom region, the four truncated supports, and
the top region. -/
inductive SimRegion : Type
  | bot
  | s36
  | s60
  | s64
  | s92
  | top
  deriving DecidableEq, Fintype

namespace SimRegion

/-- Region inclusion: bottom below everything, everything below top,
supports pairwise incomparable. -/
def le : SimRegion → SimRegion → Prop
  | _, .top => True
  | .bot, _ => True
  | .s36, .s36 => True
  | .s60, .s60 => True
  | .s64, .s64 => True
  | .s92, .s92 => True
  | _, _ => False

instance : (U V : SimRegion) → Decidable (le U V) := fun U V => by
  cases U <;> cases V <;>
    first
      | exact (inferInstance : Decidable True)
      | exact (inferInstance : Decidable False)

theorem le_refl : ∀ U, le U U := by decide

theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide

theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide

/-- The declared overlap: the meet in the region lattice. -/
def meet (U V : SimRegion) : SimRegion :=
  if le U V then U else if le V U then V else .bot

theorem meet_le_left : ∀ U V, le (meet U V) U := by decide

theorem meet_le_right : ∀ U V, le (meet U V) V := by decide

theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

/-- The support region of each payload observer. -/
def regionOf : PayloadObserver → SimRegion
  | .o36 => .s36
  | .o60 => .s60
  | .o64 => .s64
  | .o92 => .s92

/-- Declared disjointness: two distinct support regions.  The literal
receipt `supportNodes_pairwise_disjoint` backs this declaration. -/
def disj (U V : SimRegion) : Prop :=
  U ≠ V ∧ U ≠ .bot ∧ U ≠ .top ∧ V ≠ .bot ∧ V ≠ .top

/-- The top region is outside the cone below every support region, so
the containment laws stated at support regions are falsifiable. -/
theorem top_not_le_regionOf : ∀ o, ¬ le .top (regionOf o) := by decide

end SimRegion

/-! ## Regional algebras and the restriction system -/

/-- The regional algebra assignment: scalars at the bottom, the full
support diagonal at each support region, the full diagonal at the top. -/
def simRegionAlgebra : SimRegion → StarSubalgebra ℂ (Matrix (Fin 33) (Fin 33) ℂ)
  | .bot => diagonalClassSubalgebra fBot
  | .s36 => diagonalClassSubalgebra (fReg .o36)
  | .s60 => diagonalClassSubalgebra (fReg .o60)
  | .s64 => diagonalClassSubalgebra (fReg .o64)
  | .s92 => diagonalClassSubalgebra (fReg .o92)
  | .top => diagonalClassSubalgebra fTop

theorem simRegionAlgebra_regionOf (o : PayloadObserver) :
    simRegionAlgebra (SimRegion.regionOf o) =
      diagonalClassSubalgebra (fReg o) := by
  cases o <;> rfl

/-- Every regional member is a diagonal matrix. -/
theorem simRegionAlgebra_mem_diagonal {U : SimRegion}
    {X : Matrix (Fin 33) (Fin 33) ℂ} (hX : X ∈ simRegionAlgebra U) :
    ∃ d : Fin 33 → ℂ, Matrix.diagonal d = X := by
  cases U <;> · obtain ⟨c, hc⟩ := hX; exact ⟨_, hc⟩

/-- Isotony of the regional assignment along region inclusion. -/
theorem simRegionAlgebra_isotony {U V : SimRegion} (hUV : SimRegion.le U V) :
    simRegionAlgebra U ≤ simRegionAlgebra V := by
  cases U <;> cases V <;>
    first
      | exact le_rfl
      | exact absurd hUV (by decide)
      | exact diagonalClassSubalgebra_comp_le (fun _ => ()) _
      | exact diagonalClassSubalgebra_comp_le (fReg _) fTop

/-- The matrix value shared by all restriction maps: pull the diagonal
back along a basis self-map. -/
def diagPull (σ : Fin 33 → Fin 33) (X : Matrix (Fin 33) (Fin 33) ℂ) :
    Matrix (Fin 33) (Fin 33) ℂ :=
  Matrix.diagonal fun i => X (σ i) (σ i)

theorem diagPull_one (σ : Fin 33 → Fin 33) :
    diagPull σ (1 : Matrix (Fin 33) (Fin 33) ℂ) = 1 := by
  simp [diagPull]

theorem diagPull_add (σ : Fin 33 → Fin 33)
    (X Y : Matrix (Fin 33) (Fin 33) ℂ) :
    diagPull σ (X + Y) = diagPull σ X + diagPull σ Y := by
  simp [diagPull, Matrix.diagonal_add]

theorem diagPull_zero (σ : Fin 33 → Fin 33) :
    diagPull σ (0 : Matrix (Fin 33) (Fin 33) ℂ) = 0 := by
  simp [diagPull]

theorem diagPull_star (σ : Fin 33 → Fin 33)
    (X : Matrix (Fin 33) (Fin 33) ℂ) :
    diagPull σ (star X) = star (diagPull σ X) := by
  simp [diagPull, Matrix.star_eq_conjTranspose,
    Matrix.diagonal_conjTranspose]

theorem diagPull_algebraMap (σ : Fin 33 → Fin 33) (z : ℂ) :
    diagPull σ (algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) z) =
      algebraMap ℂ (Matrix (Fin 33) (Fin 33) ℂ) z := by
  rw [← diagonal_const_eq_algebraMap]
  simp [diagPull]

/-- Multiplicativity of the pullback on members of one diagonal span. -/
theorem diagPull_mul {α : Type} {f : Fin 33 → α} (σ : Fin 33 → Fin 33)
    {X Y : Matrix (Fin 33) (Fin 33) ℂ}
    (hX : X ∈ diagonalClassSubalgebra f)
    (hY : Y ∈ diagonalClassSubalgebra f) :
    diagPull σ (X * Y) = diagPull σ X * diagPull σ Y := by
  obtain ⟨c, rfl⟩ := hX
  obtain ⟨d, rfl⟩ := hY
  simp [diagPull, Matrix.diagonal_mul_diagonal]

/-- The restriction star-homomorphism between two diagonal spans induced
by a basis self-map, given that the pullback lands in the target. -/
def simRestrictHom {α β : Type} (f : Fin 33 → α) (g : Fin 33 → β)
    (σ : Fin 33 → Fin 33)
    (hmem : ∀ X ∈ diagonalClassSubalgebra f,
      diagPull σ X ∈ diagonalClassSubalgebra g) :
    diagonalClassSubalgebra f →⋆ₐ[ℂ] diagonalClassSubalgebra g where
  toFun X := ⟨diagPull σ X.1, hmem X.1 X.2⟩
  map_one' := Subtype.ext (by simp [diagPull_one])
  map_mul' X Y := Subtype.ext (by
    simp only [MulMemClass.coe_mul]
    exact diagPull_mul σ X.2 Y.2)
  map_zero' := Subtype.ext (by simp [diagPull_zero])
  map_add' X Y := Subtype.ext (by
    simp only [AddMemClass.coe_add]
    exact diagPull_add σ X.1 Y.1)
  commutes' z := Subtype.ext (diagPull_algebraMap σ z)
  map_star' X := Subtype.ext (diagPull_star σ X.1)

@[simp]
theorem simRestrictHom_apply {α β : Type} (f : Fin 33 → α) (g : Fin 33 → β)
    (σ : Fin 33 → Fin 33) (hmem) (X : diagonalClassSubalgebra f) :
    (simRestrictHom f g σ hmem X : Matrix (Fin 33) (Fin 33) ℂ) =
      diagPull σ X.1 :=
  rfl

/-- The character onto the bottom scalars: evaluation at the ambient
index. -/
def simCharToBot {α : Type} (f : Fin 33 → α) :
    diagonalClassSubalgebra f →⋆ₐ[ℂ] diagonalClassSubalgebra fBot :=
  simRestrictHom f fBot (fun _ => 0)
    (fun X _ => ⟨fun _ => X 0 0, rfl⟩)

/-- The pullback from the top diagonal onto a support region: fixed on
the support, ambient evaluation elsewhere. -/
def simPullbackToSupp (o : PayloadObserver) :
    diagonalClassSubalgebra fTop →⋆ₐ[ℂ] diagonalClassSubalgebra (fReg o) :=
  simRestrictHom fTop (fReg o) (collapseMap o)
    (fun X _ => ⟨fun x => x.elim (X 0 0) (fun j => X j j), by
      show Matrix.diagonal _ = Matrix.diagonal _
      congr 1
      funext i
      by_cases h : nodeOwner i = some o
      · simp [fReg, collapseMap, h]
      · simp [fReg, collapseMap, h]⟩)

theorem simCharToBot_apply {α : Type} (f : Fin 33 → α)
    (X : diagonalClassSubalgebra f) :
    (simCharToBot f X : Matrix (Fin 33) (Fin 33) ℂ) =
      Matrix.diagonal (fun _ => (X : Matrix (Fin 33) (Fin 33) ℂ) 0 0) :=
  rfl

theorem simPullbackToSupp_apply (o : PayloadObserver)
    (X : diagonalClassSubalgebra fTop) :
    (simPullbackToSupp o X : Matrix (Fin 33) (Fin 33) ℂ) =
      Matrix.diagonal (fun i => (X : Matrix (Fin 33) (Fin 33) ℂ)
        (collapseMap o i) (collapseMap o i)) :=
  rfl

/-- The supplied restriction system: identities on equal regions, the
bottom character on drops to the bottom, and the support pullbacks on
drops from the top.  The remaining constructor pairs are excluded by the
region order. -/
noncomputable def simRestrict : ∀ U V : SimRegion, SimRegion.le V U →
    (simRegionAlgebra U →⋆ₐ[ℂ] simRegionAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .s36, .s36, _ => StarAlgHom.id ℂ _
  | .s60, .s60, _ => StarAlgHom.id ℂ _
  | .s64, .s64, _ => StarAlgHom.id ℂ _
  | .s92, .s92, _ => StarAlgHom.id ℂ _
  | .top, .top, _ => StarAlgHom.id ℂ _
  | .s36, .bot, _ => simCharToBot (fReg .o36)
  | .s60, .bot, _ => simCharToBot (fReg .o60)
  | .s64, .bot, _ => simCharToBot (fReg .o64)
  | .s92, .bot, _ => simCharToBot (fReg .o92)
  | .top, .bot, _ => simCharToBot fTop
  | .top, .s36, _ => simPullbackToSupp .o36
  | .top, .s60, _ => simPullbackToSupp .o60
  | .top, .s64, _ => simPullbackToSupp .o64
  | .top, .s92, _ => simPullbackToSupp .o92
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

/-! ## The earned tower, net, and access cut -/

/-- The selected tower state: the declared pure state on the ambient
index.  The payload pins no state; this field is declared data. -/
noncomputable def simEarnedState : StateMatrix 33 :=
  ⟨Matrix.diagonal (Pi.single 0 1),
    ((diagonalPartition 33).isEvent 0).posSemidef, by
      rw [Matrix.trace_diagonal]
      simp [Pi.single_apply]⟩

/-- The earned consensus tower: `Record () = Fin 32` is the payload's
record-class type, `Observer () = PayloadObserver` carries the four
payload observers, and the public algebra of each observer is the
diagonal span of its record classes on its truncated support. -/
noncomputable def simEarnedTower : ConsensusTower Unit where
  indexNonempty := ⟨()⟩
  commonRefinement := fun _ _ => ⟨(), le_rfl, le_rfl⟩
  dim := fun _ => 33
  Observer := fun _ => PayloadObserver
  observerFinite := fun _ => inferInstance
  Record := fun _ => Fin 32
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ _ => OPH.TimeOrderLedger.ObserverRecordOrder.discrete (Fin 32)
  publicAlgebra := fun _ o => diagonalClassSubalgebra (fPub o)
  public_mul_comm := fun _ o X Y => diagonalClassSubalgebra_mul_comm (fPub o) X Y
  recordElement := fun _ o c =>
    ⟨Matrix.diagonal fun i => if fPub o i = some c then 1 else 0,
      fun x => if x = some c then 1 else 0, rfl⟩
  state := fun _ _ => simEarnedState.1
  state_isState := fun _ _ => simEarnedState.2
  generator := fun _ _ => 0
  observerRefine := fun _ o => o
  recordRefine := fun _ x => x
  algebraRefine := fun _ => StarAlgHom.id ℂ _
  observer_refine_refl := by intros; rfl
  observer_refine_trans := by intros; rfl
  record_refine_refl := by intros; rfl
  record_refine_trans := by intros; rfl
  algebra_refine_refl := by intros; rfl
  algebra_refine_trans := by intros; rfl
  record_order_natural := by
    simp [OPH.TimeOrderLedger.ObserverRecordOrder.discrete]
  public_natural := by
    intro r s hrs o X hX
    exact hX
  record_element_natural := by intros; rfl
  state_natural := by intros; rfl
  generator_natural := by intros; simp

/-- The record element of a class label is its class projector. -/
theorem simEarnedTower_recordElement (o : PayloadObserver) (c : Fin 32) :
    (simEarnedTower.recordElement () o c).1 =
      Matrix.diagonal (fun i => if fPub o i = some c then 1 else 0) :=
  rfl

/-- The tower's record type has exactly the payload's 32 record
classes. -/
theorem simEarnedTower_record_card :
    Fintype.card (simEarnedTower.Record ()) = 32 :=
  rfl

/-- The tower carries exactly the four payload observers. -/
theorem simEarnedTower_observer_card :
    Fintype.card (simEarnedTower.Observer ()) = 4 :=
  rfl

/-- Elementwise commutation of regional members across any two regions:
every regional member is diagonal. -/
theorem simRegionAlgebra_commute {U V : SimRegion}
    (X : simRegionAlgebra U) (Y : simRegionAlgebra V) :
    Commute (X : Matrix (Fin 33) (Fin 33) ℂ)
      (Y : Matrix (Fin 33) (Fin 33) ℂ) := by
  obtain ⟨dX, hX⟩ := simRegionAlgebra_mem_diagonal X.2
  obtain ⟨dY, hY⟩ := simRegionAlgebra_mem_diagonal Y.2
  show (X : Matrix (Fin 33) (Fin 33) ℂ) * Y =
    (Y : Matrix (Fin 33) (Fin 33) ℂ) * X
  rw [← hX, ← hY, Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm _ _

/-- The earned finite causal observer net over the earned tower: the
region lattice `SimRegion`, the diagonal regional algebras, the
character-and-pullback restriction system, distinct-support
disjointness, and identity repair. -/
noncomputable def simEarnedNet : FiniteCausalObserverNet simEarnedTower where
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
  localAlgebra := fun _ => simRegionAlgebra
  isotony := fun _ {U V} h => simRegionAlgebra_isotony h
  locality := fun _ {U V} _ X Y => simRegionAlgebra_commute X Y
  restrict := fun _ {U V} h => simRestrict U V h
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

/-- The earned access cut: observer regions are the payload support
regions and accessible algebras are the joint-class diagonal spans. -/
noncomputable def simEarnedCut :
    ObserverAccessCut simEarnedTower simEarnedNet () where
  observerRegion := SimRegion.regionOf
  accessibleAlgebra := fun o => diagonalClassSubalgebra (fAcc o)
  accessible_le_region := fun o => by
    cases o <;>
      exact diagonalClassSubalgebra_comp_le
        (Option.map fun j => (nodeRecordClass j, nodeCompanionClass j))
        (fReg _)
  public_le_accessible := fun o =>
    diagonalClassSubalgebra_comp_le (Option.map Prod.fst) (fAcc o)
  restrict_preserves := by
    intro o U V hU hVU X hX
    have hU' : SimRegion.le U (SimRegion.regionOf o) := hU
    have hVU' : SimRegion.le V U := hVU
    cases U <;> cases V <;>
      first
        | exact hX
        | exact absurd hVU' (by decide)
        | exact ⟨fun _ => X.1 (0 : Fin 33) (0 : Fin 33), rfl⟩
        | (cases o <;> exact absurd hU' (by decide))

/-! ## The receipt battery

The receipts mirror the battery of the synthetic access-cut witnesses:
a declared disjoint region pair, the two interposition strictness
receipts where the payload data provides them, the repair-stability
receipt, and the accessible-glue receipt.  Where the truncated data
collapses an interposition, the collapse is proved as a named gap
equality instead of being repaired by synthetic data. -/

/-- The net keeps a declared pair of disjoint observer support
regions. -/
theorem simEarnedNet_has_disjoint_pair :
    simEarnedNet.disjoint () (SimRegion.regionOf PayloadObserver.o36)
      (SimRegion.regionOf PayloadObserver.o60) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-- The separating element for observer 64: the projector of the joint
class (18, 13), which splits record class 18 by companion class. -/
def o64Separator : Matrix (Fin 33) (Fin 33) ℂ :=
  Matrix.diagonal fun i => if fAcc .o64 i = some (18, 13) then 1 else 0

theorem o64Separator_mem_accessible :
    o64Separator ∈ diagonalClassSubalgebra (fAcc .o64) :=
  ⟨fun x => if x = some (18, 13) then 1 else 0, rfl⟩

theorem o64Separator_not_mem_public :
    o64Separator ∉ diagonalClassSubalgebra (fPub .o64) := by
  rintro ⟨c, hc⟩
  have h17 := congrFun (congrFun hc 17) 17
  have h20 := congrFun (congrFun hc 20) 20
  simp only [o64Separator, Matrix.diagonal_apply_eq] at h17 h20
  have e1 : fPub PayloadObserver.o64 17 = some 18 := by decide
  have e2 : fAcc PayloadObserver.o64 17 = some (18, 13) := by decide
  have e3 : fPub PayloadObserver.o64 20 = some 18 := by decide
  have e4 : fAcc PayloadObserver.o64 20 = some (18, 3) := by decide
  rw [e1, e2, if_pos rfl] at h17
  rw [e3, e4, if_neg (by decide)] at h20
  rw [h17] at h20
  exact one_ne_zero h20

/-- Strictness, earned from the payload: observer 64's truncated support
realizes record class 18 with two distinct companion classes, so its
public record algebra sits strictly inside its accessible algebra. -/
theorem simEarnedCut_public_lt_accessible_o64 :
    simEarnedTower.publicAlgebra () PayloadObserver.o64 <
      simEarnedCut.accessibleAlgebra PayloadObserver.o64 := by
  rw [SetLike.lt_iff_le_and_exists]
  exact ⟨simEarnedCut.public_le_accessible _, o64Separator,
    o64Separator_mem_accessible, o64Separator_not_mem_public⟩

/-- The companion class attached to each record class realized on a
truncated support, read off the first matching position; unrealized
classes get the junk value `0`. -/
def companionOfRecord (o : PayloadObserver) (r : Fin 32) : Fin 16 :=
  (((List.finRange 8).find? fun j =>
    PayloadObserver.recordClassAt o j == r).map
      (PayloadObserver.companionClassAt o)).getD 0

theorem fAcc_factor_o36 : ∀ i, fAcc PayloadObserver.o36 i =
    Option.map (fun r => (r, companionOfRecord PayloadObserver.o36 r))
      (fPub PayloadObserver.o36 i) := by decide

theorem fAcc_factor_o60 : ∀ i, fAcc PayloadObserver.o60 i =
    Option.map (fun r => (r, companionOfRecord PayloadObserver.o60 r))
      (fPub PayloadObserver.o60 i) := by decide

theorem fAcc_factor_o92 : ∀ i, fAcc PayloadObserver.o92 i =
    Option.map (fun r => (r, companionOfRecord PayloadObserver.o92 r))
      (fPub PayloadObserver.o92 i) := by decide

/-- Named gap: observer 36's truncated support realizes each record
class with a single companion class, so the public and accessible
algebras coincide and no strict interposition receipt exists for it. -/
theorem strictnessGap_public_eq_accessible_o36 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o36 =
      simEarnedTower.publicAlgebra () PayloadObserver.o36 := by
  refine le_antisymm ?_ (simEarnedCut.public_le_accessible _)
  show diagonalClassSubalgebra (fAcc PayloadObserver.o36) ≤
    diagonalClassSubalgebra (fPub PayloadObserver.o36)
  rw [funext fAcc_factor_o36]
  exact diagonalClassSubalgebra_comp_le
    (Option.map fun r => (r, companionOfRecord PayloadObserver.o36 r))
    (fPub PayloadObserver.o36)

/-- Named gap: the record classes on observer 60's truncated support are
pairwise distinct, so its public and accessible algebras coincide. -/
theorem strictnessGap_public_eq_accessible_o60 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o60 =
      simEarnedTower.publicAlgebra () PayloadObserver.o60 := by
  refine le_antisymm ?_ (simEarnedCut.public_le_accessible _)
  show diagonalClassSubalgebra (fAcc PayloadObserver.o60) ≤
    diagonalClassSubalgebra (fPub PayloadObserver.o60)
  rw [funext fAcc_factor_o60]
  exact diagonalClassSubalgebra_comp_le
    (Option.map fun r => (r, companionOfRecord PayloadObserver.o60 r))
    (fPub PayloadObserver.o60)

/-- Named gap: the record classes on observer 92's truncated support are
pairwise distinct, so its public and accessible algebras coincide. -/
theorem strictnessGap_public_eq_accessible_o92 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o92 =
      simEarnedTower.publicAlgebra () PayloadObserver.o92 := by
  refine le_antisymm ?_ (simEarnedCut.public_le_accessible _)
  show diagonalClassSubalgebra (fAcc PayloadObserver.o92) ≤
    diagonalClassSubalgebra (fPub PayloadObserver.o92)
  rw [funext fAcc_factor_o92]
  exact diagonalClassSubalgebra_comp_le
    (Option.map fun r => (r, companionOfRecord PayloadObserver.o92 r))
    (fPub PayloadObserver.o92)

/-- The separating element for observer 36's regional side: the node
projector of basis index 3, payload node 15, one of the two nodes
sharing the joint class (20, 6). -/
def o36Separator : Matrix (Fin 33) (Fin 33) ℂ :=
  Matrix.diagonal fun i => if fReg .o36 i = some 3 then 1 else 0

theorem o36Separator_mem_regional :
    o36Separator ∈ diagonalClassSubalgebra (fReg .o36) :=
  ⟨fun x => if x = some 3 then 1 else 0, rfl⟩

theorem o36Separator_not_mem_accessible :
    o36Separator ∉ diagonalClassSubalgebra (fAcc .o36) := by
  rintro ⟨c, hc⟩
  have h3 := congrFun (congrFun hc 3) 3
  have h4 := congrFun (congrFun hc 4) 4
  simp only [o36Separator, Matrix.diagonal_apply_eq] at h3 h4
  have e1 : fAcc PayloadObserver.o36 3 = some (20, 6) := by decide
  have e2 : fReg PayloadObserver.o36 3 = some 3 := by decide
  have e3 : fAcc PayloadObserver.o36 4 = some (20, 6) := by decide
  have e4 : fReg PayloadObserver.o36 4 = some 4 := by decide
  rw [e1, e2, if_pos rfl] at h3
  rw [e3, e4, if_neg (by decide)] at h4
  rw [h3] at h4
  exact one_ne_zero h4

/-- Strictness, earned from the payload: observer 36's truncated support
carries two nodes with the same joint class (20, 6), so its accessible
algebra sits strictly inside the regional algebra of its support. -/
theorem simEarnedCut_accessible_lt_regional_o36 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o36 <
      simEarnedNet.localAlgebra ()
        (SimRegion.regionOf PayloadObserver.o36) := by
  rw [SetLike.lt_iff_le_and_exists]
  exact ⟨simEarnedCut.accessible_le_region _, o36Separator,
    o36Separator_mem_regional, o36Separator_not_mem_accessible⟩

/-- The node recovered from a joint class on a support: the first basis
index realizing the class. -/
def nodeOfJointClass (o : PayloadObserver) (p : Fin 32 × Fin 16) :
    Option (Fin 33) :=
  (List.finRange 33).find? fun i => fAcc o i == some p

theorem fReg_factor_o60 : ∀ i, fReg PayloadObserver.o60 i =
    Option.bind (fAcc PayloadObserver.o60 i)
      (nodeOfJointClass PayloadObserver.o60) := by decide

theorem fReg_factor_o64 : ∀ i, fReg PayloadObserver.o64 i =
    Option.bind (fAcc PayloadObserver.o64 i)
      (nodeOfJointClass PayloadObserver.o64) := by decide

theorem fReg_factor_o92 : ∀ i, fReg PayloadObserver.o92 i =
    Option.bind (fAcc PayloadObserver.o92 i)
      (nodeOfJointClass PayloadObserver.o92) := by decide

/-- Named gap: the joint classes on observer 60's truncated support are
pairwise distinct, so its accessible algebra equals the full regional
algebra of its support and no strict regional receipt exists. -/
theorem regionalGap_accessible_eq_regional_o60 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o60 =
      simEarnedNet.localAlgebra ()
        (SimRegion.regionOf PayloadObserver.o60) := by
  refine le_antisymm (simEarnedCut.accessible_le_region _) ?_
  show diagonalClassSubalgebra (fReg PayloadObserver.o60) ≤
    diagonalClassSubalgebra (fAcc PayloadObserver.o60)
  rw [funext fReg_factor_o60]
  exact diagonalClassSubalgebra_comp_le
    (fun x => Option.bind x (nodeOfJointClass PayloadObserver.o60))
    (fAcc PayloadObserver.o60)

/-- Named gap: the joint classes on observer 64's truncated support are
pairwise distinct, so its accessible algebra equals its full regional
algebra. -/
theorem regionalGap_accessible_eq_regional_o64 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o64 =
      simEarnedNet.localAlgebra ()
        (SimRegion.regionOf PayloadObserver.o64) := by
  refine le_antisymm (simEarnedCut.accessible_le_region _) ?_
  show diagonalClassSubalgebra (fReg PayloadObserver.o64) ≤
    diagonalClassSubalgebra (fAcc PayloadObserver.o64)
  rw [funext fReg_factor_o64]
  exact diagonalClassSubalgebra_comp_le
    (fun x => Option.bind x (nodeOfJointClass PayloadObserver.o64))
    (fAcc PayloadObserver.o64)

/-- Named gap: the joint classes on observer 92's truncated support are
pairwise distinct, so its accessible algebra equals its full regional
algebra. -/
theorem regionalGap_accessible_eq_regional_o92 :
    simEarnedCut.accessibleAlgebra PayloadObserver.o92 =
      simEarnedNet.localAlgebra ()
        (SimRegion.regionOf PayloadObserver.o92) := by
  refine le_antisymm (simEarnedCut.accessible_le_region _) ?_
  show diagonalClassSubalgebra (fReg PayloadObserver.o92) ≤
    diagonalClassSubalgebra (fAcc PayloadObserver.o92)
  rw [funext fReg_factor_o92]
  exact diagonalClassSubalgebra_comp_le
    (fun x => Option.bind x (nodeOfJointClass PayloadObserver.o92))
    (fAcc PayloadObserver.o92)

/-- Upper strictness: every accessible algebra is a proper subalgebra of
the private matrix algebra; the raising unit separates them. -/
theorem simEarnedCut_accessible_lt_top (o : PayloadObserver) :
    simEarnedCut.accessibleAlgebra o <
      (⊤ : StarSubalgebra ℂ
        (ConsensusTower.PrivateAlgebra simEarnedTower ())) := by
  rw [SetLike.lt_iff_le_and_exists]
  exact ⟨le_top, simRaise, StarSubalgebra.mem_top, simRaise_not_mem (fAcc o)⟩

/-- The cut satisfies the repair-stability receipt: the net's repair
maps are identities. -/
theorem simEarnedCut_repairStable : RepairStableAccessCut simEarnedCut := by
  constructor
  intro o U hU X hX
  exact hX

/-- The cut satisfies the accessible-glue receipt over the singleton
descent packet of the earned net. -/
theorem simEarnedCut_accessibleGlueClosure (o : PayloadObserver) :
    simEarnedCut.AccessibleGlueClosure o simEarnedNet.singletonDescent :=
  simEarnedCut.accessibleGlueClosure_singleton o

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.realizedRecordClassCounts_total
#print axioms OPH.QFT.supportNodes_pairwise_disjoint
#print axioms OPH.QFT.jointTable_mass
#print axioms OPH.QFT.jointTable_covers_truncated_pairs
#print axioms OPH.QFT.truncated_record_classes_nonconstant
#print axioms OPH.QFT.nodeRecordClass_nodeIndex
#print axioms OPH.QFT.nodeCompanionClass_nodeIndex
#print axioms OPH.QFT.simEarnedNet_has_disjoint_pair
#print axioms OPH.QFT.simEarnedCut_public_lt_accessible_o64
#print axioms OPH.QFT.strictnessGap_public_eq_accessible_o36
#print axioms OPH.QFT.strictnessGap_public_eq_accessible_o60
#print axioms OPH.QFT.strictnessGap_public_eq_accessible_o92
#print axioms OPH.QFT.simEarnedCut_accessible_lt_regional_o36
#print axioms OPH.QFT.regionalGap_accessible_eq_regional_o60
#print axioms OPH.QFT.regionalGap_accessible_eq_regional_o64
#print axioms OPH.QFT.regionalGap_accessible_eq_regional_o92
#print axioms OPH.QFT.simEarnedCut_accessible_lt_top
#print axioms OPH.QFT.simEarnedCut_repairStable
#print axioms OPH.QFT.simEarnedCut_accessibleGlueClosure
