import ObserverPatchHolography.Provenance.HistoryCausalInvariance

/-! Sharp rail readings on a fixed list of native appointments. The driver
supplies the owner and the realized exchange bit. A repeated rail reading is
a fresh public version. The compiled log uses longest authenticated ancestry. -/
namespace OPH.Provenance.GuardedRail

variable {k n : ℕ}

def initial (rail : Fin k → Bool) : VersionedState (Fin k) Bool (Fin n) :=
  ⟨rail, fun _ => none⟩

def readRail (q : VersionedState (Fin k) Bool (Fin n))
    (e : Fin n) (a : Fin k) (flip : Bool) : CommitSpec (Fin k) Bool (Fin n) where
  eventId := e
  readSet := {a}
  writeSet := {a}
  causalSupp := {a}
  supp_subset_read := fun _ h => h
  writeValue := fun _ => xor (q.value a) flip

def compile (owner : Fin n → Fin k) (flip : Fin n → Bool)
    (q : VersionedState (Fin k) Bool (Fin n)) :
    List (Fin n) → List (CommitSpec (Fin k) Bool (Fin n))
  | [] => []
  | e :: es =>
    let s := readRail q e (owner e) (flip e)
    s :: compile owner flip (s.apply q) es

theorem identifiers (owner : Fin n → Fin k) (flip : Fin n → Bool)
    (q : VersionedState (Fin k) Bool (Fin n)) (es : List (Fin n)) :
    (compile owner flip q es).map CommitSpec.eventId = es := by
  induction es generalizing q with
  | nil => rfl
  | cons e es ih => simp [compile, readRail, ih]

noncomputable def railLog (rail : Fin k → Bool) (owner : Fin n → Fin k)
    (flip : Fin n → Bool) : SemanticEventLog (Fin k) Bool (Fin n) :=
  let q := initial (n := n) rail
  let specs := compile owner flip q (List.finRange n)
  semanticEventLogOfExecutionWithSourceHeight q specs
    (by intro s _ a; simp [q, initial])
    (by rw [identifiers]; exact List.nodup_finRange n)
    (by intro e; rw [identifiers]; simp)

theorem authenticated_rank (rail : Fin k → Bool) (owner : Fin n → Fin k)
    (flip : Fin n → Bool) {e f : Fin n} (h : (railLog rail owner flip).ParentEdge e f) :
    (railLog rail owner flip).rank e < (railLog rail owner flip).rank f :=
  (railLog rail owner flip).rank_lt_of_parent h

/-- The conditional sharp transfer acts on a one-site constant and sign by
one and one half. The product word receives the product of its sign factors. -/
theorem sharpTransfer (c s : ℚ) :
    ((3/4:ℚ)*(c+s) + (1/4:ℚ)*(c-s) = c + s/2) ∧
    ((1/4:ℚ)*(c+s) + (3/4:ℚ)*(c-s) = c - s/2) := by
  constructor <;> ring

theorem retainedWord {α β : Type*} [DecidableEq β] (child : α ↪ β)
    (word : Finset α) :
    (1/2:ℚ) ^ (word.map child).card = (1/2:ℚ) ^ word.card := by simp

theorem fourChildParity : (1/2:ℚ)^4 ≠ 1/2 := by norm_num

end OPH.Provenance.GuardedRail
