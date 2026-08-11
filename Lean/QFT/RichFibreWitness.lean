import QFT.FiniteCausalObserverNet

/-!
# Rich-fibre witness: literals of the preregistered run `e1_prereg2_64k_20260810`

Exact mirror of the committed extraction payload of `oph-physics-sim`
(schema `oph.sim.e1_rich_fibre_payload.v2`, seed 20260811).  The observers were selected by the
preregistered greedy-disjoint rule (support node identities only), the
windows are the declared first-20 truncations, and every receipt below is
decided by kernel computation on the literals alone.  The module asserts
nothing beyond the payload at the source-window level: no synthetic source
window node, class, or fibre enters.  The downstream adapter does add one
auxiliary non-source base-point basis index.

Provenance (verbatim from the payload's provenance block):

* `run_id`: `e1_prereg2_64k_20260810`, `seed`: `20260811`,
  `run_git_commit`: `3f44ca6167feb8fc73e7b20042dc6b1464e59e23`
* `receipt_sha256`: `6a1de14805d8d784a519fc7e5f5673f8dbf827a575c8da59aba5e8b672496926`
* `freezeout_sha256`: `cd87da6ccd60651790d6e488105ab387fa553cd79fdc43149d03031944d074ab`
* `observer_views_sha256`: `07ab3e555fdd8e0d50091794e3e1b20c4e98c96e3282fbee68c0c3b413fc0619`
-/

namespace OPH.QFT

set_option maxRecDepth 16384

/-! ## The selected observers as literals -/

/-- The greedy-disjoint `patch_observer` rows, named by `observer_id`. -/
inductive RichObserver : Type
  | o86
  | o88
  | o247
  | o384
  deriving DecidableEq, Fintype

namespace RichObserver

/-- Payload field `observer_id`. -/
def observerId : RichObserver → ℕ
  | .o86 => 86
  | .o88 => 88
  | .o247 => 247
  | .o384 => 384

/-- Payload field `file_order_scan_rank`. -/
def scanRank : RichObserver → ℕ
  | .o86 => 1
  | .o88 => 2
  | .o247 => 4
  | .o384 => 8

/-- Payload field `truncated_support_nodes`. -/
def supportNodes : RichObserver → List ℕ
  | .o86 => [86, 44, 52, 65, 73, 78, 94, 99, 107, 120, 128, 141, 154, 10, 18, 23, 31, 36, 57, 26]
  | .o88 => [88, 46, 54, 67, 75, 80, 96, 101, 109, 122, 130, 143, 156, 12, 17, 20, 25, 33, 38, 59]
  | .o247 => [247, 158, 179, 192, 213, 226, 234, 260, 268, 281, 302, 315, 336, 90, 103, 116, 124, 137, 145, 171]
  | .o384 => [384, 274, 295, 316, 329, 350, 363, 405, 418, 439, 473, 494, 528, 185, 206, 219, 240, 253, 261, 308]

/-- Payload field `truncated_record_classes`. -/
def recordClassAt : RichObserver → Fin 20 → Fin 32
  | .o86 => ![19, 26, 27, 4, 6, 24, 27, 25, 1, 16, 24, 28, 18, 11, 23, 16, 13, 25, 5, 24]
  | .o88 => ![16, 16, 19, 5, 22, 0, 29, 8, 3, 2, 28, 1, 9, 11, 5, 10, 25, 8, 17, 0]
  | .o247 => ![9, 26, 15, 4, 28, 11, 21, 5, 13, 31, 27, 5, 15, 0, 17, 28, 26, 22, 19, 15]
  | .o384 => ![10, 7, 3, 27, 20, 20, 19, 24, 13, 29, 9, 20, 11, 15, 22, 11, 29, 21, 15, 25]

/-- Payload field `truncated_companion_classes`. -/
def companionClassAt : RichObserver → Fin 20 → Fin 16
  | .o86 => ![15, 15, 4, 8, 4, 8, 8, 15, 1, 8, 3, 8, 15, 12, 4, 8, 4, 3, 8, 3]
  | .o88 => ![12, 12, 8, 15, 15, 12, 3, 15, 15, 8, 12, 15, 12, 1, 12, 1, 15, 1, 8, 8]
  | .o247 => ![15, 12, 8, 1, 3, 1, 8, 15, 3, 1, 15, 12, 15, 8, 12, 15, 3, 12, 12, 15]
  | .o384 => ![12, 3, 8, 8, 1, 15, 15, 15, 15, 8, 0, 8, 12, 8, 12, 15, 12, 1, 12, 8]

/-- Payload field `truncated_split_fibres`. -/
def splitFibreClasses : RichObserver → List ℕ
  | .o86 => [27, 24, 25]
  | .o88 => [5, 0, 8]
  | .o247 => [26, 15, 28, 5]
  | .o384 => [20, 29, 11, 15]

/-- Payload field `designated_block_record_class`. -/
def designatedBlockClass : RichObserver → ℕ
  | .o86 => 27
  | .o88 => 5
  | .o247 => 26
  | .o384 => 20

end RichObserver

open RichObserver

/-! ## Literal receipts, kernel-decided -/

/-- Each truncated window has exactly 20 nodes. -/
theorem richSupport_length :
    ∀ o : RichObserver, (supportNodes o).length = 20 := by decide

/-- The preregistered disjointness gate, restated on the literals: the
windows are pairwise disjoint. -/
theorem richSupport_pairwise_disjoint :
    ∀ o o' : RichObserver, o ≠ o' →
      ∀ n ∈ supportNodes o, n ∉ supportNodes o' := by decide
/-- The split-fibre census: a listed split class occupies at least two
window positions with at least two distinct companion classes. -/
theorem richSplit_census : ∀ o : RichObserver, ∀ c ∈ splitFibreClasses o,
    2 ≤ ((List.finRange 20).filter
      (fun j => (recordClassAt o j : ℕ) = c)).length ∧
    2 ≤ (((List.finRange 20).filter
      (fun j => (recordClassAt o j : ℕ) = c)).map
        (fun j => (companionClassAt o j : ℕ))).dedup.length := by decide

/-- Gate G1 on the literals: every observer realizes at least one split
fibre. -/
theorem richSplit_g1 :
    ∀ o : RichObserver, 1 ≤ (splitFibreClasses o).length := by decide

/-- The exact split-fibre count of each observer. -/
theorem richSplit_counts :
    (splitFibreClasses RichObserver.o86).length = 3 ∧ (splitFibreClasses RichObserver.o88).length = 3 ∧ (splitFibreClasses RichObserver.o247).length = 4 ∧ (splitFibreClasses RichObserver.o384).length = 4 := by decide

/-! ## Realized split-fibre window positions -/

/-- Observer 86, split record class 27: window positions [2, 6] with companion classes [4, 8]. -/
theorem split_o86_c27 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o86 j : ℕ) = 27)).map (fun j => (j : ℕ)) = [2, 6] := by decide

/-- Observer 86, split record class 24: window positions [5, 10, 19] with companion classes [8, 3, 3]. -/
theorem split_o86_c24 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o86 j : ℕ) = 24)).map (fun j => (j : ℕ)) = [5, 10, 19] := by decide

/-- Observer 86, split record class 25: window positions [7, 17] with companion classes [15, 3]. -/
theorem split_o86_c25 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o86 j : ℕ) = 25)).map (fun j => (j : ℕ)) = [7, 17] := by decide

/-- Observer 88, split record class 5: window positions [3, 14] with companion classes [15, 12]. -/
theorem split_o88_c5 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o88 j : ℕ) = 5)).map (fun j => (j : ℕ)) = [3, 14] := by decide

/-- Observer 88, split record class 0: window positions [5, 19] with companion classes [12, 8]. -/
theorem split_o88_c0 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o88 j : ℕ) = 0)).map (fun j => (j : ℕ)) = [5, 19] := by decide

/-- Observer 88, split record class 8: window positions [7, 17] with companion classes [15, 1]. -/
theorem split_o88_c8 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o88 j : ℕ) = 8)).map (fun j => (j : ℕ)) = [7, 17] := by decide

/-- Observer 247, split record class 26: window positions [1, 16] with companion classes [12, 3]. -/
theorem split_o247_c26 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o247 j : ℕ) = 26)).map (fun j => (j : ℕ)) = [1, 16] := by decide

/-- Observer 247, split record class 15: window positions [2, 12, 19] with companion classes [8, 15, 15]. -/
theorem split_o247_c15 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o247 j : ℕ) = 15)).map (fun j => (j : ℕ)) = [2, 12, 19] := by decide

/-- Observer 247, split record class 28: window positions [4, 15] with companion classes [3, 15]. -/
theorem split_o247_c28 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o247 j : ℕ) = 28)).map (fun j => (j : ℕ)) = [4, 15] := by decide

/-- Observer 247, split record class 5: window positions [7, 11] with companion classes [15, 12]. -/
theorem split_o247_c5 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o247 j : ℕ) = 5)).map (fun j => (j : ℕ)) = [7, 11] := by decide

/-- Observer 384, split record class 20: window positions [4, 5, 11] with companion classes [1, 15, 8]. -/
theorem split_o384_c20 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o384 j : ℕ) = 20)).map (fun j => (j : ℕ)) = [4, 5, 11] := by decide

/-- Observer 384, split record class 29: window positions [9, 16] with companion classes [8, 12]. -/
theorem split_o384_c29 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o384 j : ℕ) = 29)).map (fun j => (j : ℕ)) = [9, 16] := by decide

/-- Observer 384, split record class 11: window positions [12, 15] with companion classes [12, 15]. -/
theorem split_o384_c11 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o384 j : ℕ) = 11)).map (fun j => (j : ℕ)) = [12, 15] := by decide

/-- Observer 384, split record class 15: window positions [13, 18] with companion classes [8, 12]. -/
theorem split_o384_c15 :
    ((List.finRange 20).filter (fun j => (recordClassAt RichObserver.o384 j : ℕ) = 15)).map (fun j => (j : ℕ)) = [13, 18] := by decide

/-! ## Ambient indexing: `Fin 81` = 4 windows of 20 plus a base point -/

def richObsIdx : RichObserver → Fin 4
  | .o86 => 0
  | .o88 => 1
  | .o247 => 2
  | .o384 => 3

/-- Window position `j` of observer `o` occupies ambient index
`richObsIdx o * 20 + j + 1`; index 0 is the base point. -/
def richNodeIndex (o : RichObserver) (j : Fin 20) : Fin 81 :=
  ⟨(richObsIdx o : ℕ) * 20 + (j : ℕ) + 1, by omega⟩

theorem richNodeIndex_injective :
    ∀ (o o' : RichObserver) (j j' : Fin 20),
      richNodeIndex o j = richNodeIndex o' j' → o = o' ∧ j = j' := by decide

def richNodeOwner : Fin 81 → Option RichObserver := fun i =>
  if i = 0 then none
  else some (match ((i : ℕ) - 1) / 20 with
    | 0 => RichObserver.o86
    | 1 => RichObserver.o88
    | 2 => RichObserver.o247
    | 3 => RichObserver.o384
    | _ => RichObserver.o384)

def richNodePos : Fin 81 → Fin 20 := fun i =>
  ⟨((i : ℕ) - 1) % 20, by omega⟩

theorem richNodeOwner_index : ∀ (o : RichObserver) (j : Fin 20),
    richNodeOwner (richNodeIndex o j) = some o := by decide

theorem richNodePos_index : ∀ (o : RichObserver) (j : Fin 20),
    richNodePos (richNodeIndex o j) = j := by decide

/-- Ambient record class: the literal class at the owning window
position, and the zero default at the base point. -/
def richNodeRecordClass : Fin 81 → Fin 32 := fun i =>
  match richNodeOwner i with
  | none => 0
  | some o => recordClassAt o (richNodePos i)

def richNodeCompanionClass : Fin 81 → Fin 16 := fun i =>
  match richNodeOwner i with
  | none => 0
  | some o => companionClassAt o (richNodePos i)

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.richSupport_length
#print axioms OPH.QFT.richSupport_pairwise_disjoint
#print axioms OPH.QFT.richSplit_census
#print axioms OPH.QFT.richSplit_counts
#print axioms OPH.QFT.richNodeIndex_injective
