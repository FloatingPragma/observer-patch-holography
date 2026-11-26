import ObserverPatchHolography.Execution.CumulativeAttemptCapacity

/-! # Measured ballots for canonical single-site Boolean repair

Each edge compares its two public bits for equality. Its endpoints can read
the comparison by an echo on extremal ternary rails. A ballot lights when
the site has a neighbor and every incident comparison reports a difference.
A scanner consults a covering list of stored ballots and chooses its first
positive read. The allowance counts all scheduler invocations.
-/

namespace OPH.Execution.MeasuredBallots

open OPH OPH.Locality

noncomputable section

variable {n m : ℕ}

/-- A finite overlap graph with one declared orientation per comparison. -/
structure Layout (n m : ℕ) where
  left : Fin m → Fin n
  right : Fin m → Fin n
  apart : ∀ e, left e ≠ right e

def Layout.carrier (g : Layout n m) : OPHCarrier where
  Patch := Fin n
  State := fun _ => Bool
  Edge := Fin m
  src := g.left
  tgt := g.right
  Iface := fun _ => Bool
  projSrc := fun _ b => b
  projTgt := fun _ b => b
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := by simp
  dist_eq_zero := by intros; split_ifs <;> simp_all

def Incident (g : Layout n m) (i : Fin n) (e : Fin m) : Prop :=
  g.left e = i ∨ g.right e = i

def Ballot (g : Layout n m) (x : Fin n → Bool) (i : Fin n) : Prop :=
  (∃ e, Incident g i e) ∧
    ∀ e, Incident g i e → x (g.left e) ≠ x (g.right e)

private theorem new_bit {a b : Bool} (h : a ≠ b) : (!a) = b := by
  cases a <;> cases b <;> simp_all

/-- A flip interchanges "equal" and "different" on each incident edge. -/
theorem flip_agrees (g : Layout n m) (x : Fin n → Bool) (i : Fin n) (e : Fin m)
    (he : Incident g i e) :
    edgeConsistentAt (C := g.carrier) e (Function.update x i (!(x i))) ↔
      x (g.left e) ≠ x (g.right e) := by
  change Function.update x i (!(x i)) (g.left e) =
    Function.update x i (!(x i)) (g.right e) ↔ _
  rcases he with h | h
  · have hn : g.right e ≠ i := h ▸ Ne.symm (g.apart e)
    rw [h, Function.update_self, Function.update_of_ne hn]
    cases x i <;> cases x (g.right e) <;> decide
  · have hn : g.left e ≠ i := h ▸ g.apart e
    rw [h, Function.update_self, Function.update_of_ne hn]
    cases x i <;> cases x (g.left e) <;> decide

/-- The measured Boolean ballot is exactly the finite-stack firing test. -/
theorem ballot_iff_fires (g : Layout n m) (x : Fin n → Bool) (i : Fin n) :
    Ballot g x i ↔ localRepair g.carrier i x ≠ x := by
  rw [localRepair_ne_iff]
  constructor
  · rintro ⟨⟨e, he⟩, hall⟩
    exact ⟨⟨e, he, hall e he⟩,
      ⟨!(x i), fun f hf => (flip_agrees g x i f hf).2 (hall f hf)⟩⟩
  · rintro ⟨⟨e, he, hbad⟩, ⟨v, hv⟩⟩
    have different : x i ≠ v := by
      intro h
      have hfixed : Function.update x i v = x := Function.update_eq_self_iff.mpr h.symm
      exact hbad ((congrArg (fun y : Records g.carrier => edgeConsistentAt e y) hfixed).mp
        (hv e he))
    have hflip := new_bit different
    refine ⟨⟨e, he⟩, fun f hf => (flip_agrees g x i f hf).1 ?_⟩
    rw [hflip]
    exact hv f hf

/-- Any solving replacement is unique at an enabled Boolean site. -/
theorem canonical_flip (g : Layout n m) (x : Fin n → Bool) (i : Fin n)
    (h : Ballot g x i) :
    localRepair g.carrier i x = Function.update x i (!(x i)) := by
  have hf := (localRepair_ne_iff g.carrier i x).1 ((ballot_iff_fires g x i).1 h)
  rw [localRepair_of_fire g.carrier i x hf]
  have hd : x i ≠ Classical.choose hf.2 := by
    intro heq
    have hsame : Function.update x i (Classical.choose hf.2) = x :=
      Function.update_eq_self_iff.mpr heq.symm
    exact (ballot_iff_fires g x i).1 h ((localRepair_of_fire g.carrier i x hf).trans hsame)
  exact congrArg (fun v : Bool => Function.update x i v) (new_bit hd).symm

/-- Two owners share an observed comparison, including the owner with itself when incident. -/
def Touches (g : Layout n m) (i j : Fin n) : Prop :=
  ∃ e, Incident g i e ∧ Incident g j e

private theorem off_edge (g : Layout n m) (x : Fin n → Bool) (i : Fin n) (e : Fin m)
    (h : ¬ Incident g i e) :
    (Function.update x i (!(x i)) (g.left e) ≠ Function.update x i (!(x i)) (g.right e)) ↔
      x (g.left e) ≠ x (g.right e) := by
  rw [Function.update_of_ne (fun he => h (Or.inl he)),
      Function.update_of_ne (fun he => h (Or.inr he))]

/-- The initial survey remains accurate by withdrawing precisely the lit touching ballots. -/
theorem ballot_after_flip (g : Layout n m) (x : Fin n → Bool) (i j : Fin n)
    (hi : Ballot g x i) :
    Ballot g (Function.update x i (!(x i))) j ↔ Ballot g x j ∧ ¬ Touches g i j := by
  constructor
  · intro hj
    have apart : ¬ Touches g i j := by
      rintro ⟨e, he, hj'⟩
      exact hj.2 e hj' ((flip_agrees g x i e he).2 (hi.2 e he))
    refine ⟨⟨hj.1, fun e he => ?_⟩, apart⟩
    exact (off_edge g x i e (fun hi' => apart ⟨e, hi', he⟩)).1 (hj.2 e he)
  · rintro ⟨hj, apart⟩
    refine ⟨hj.1, fun e he => ?_⟩
    exact (off_edge g x i e (fun hi' => apart ⟨e, hi', he⟩)).2 (hj.2 e he)

/-- The fallback also pays one invocation when every ballot is dark. -/
def scan {α : Type} (lit : α → Bool) (fallback : α) : List α → α
  | [] => fallback
  | i :: rest => if lit i then i else scan lit fallback rest

theorem scan_lit {α : Type} (lit : α → Bool) (fallback : α) (l : List α)
    (h : ∃ i ∈ l, lit i = true) : lit (scan lit fallback l) = true := by
  induction l with
  | nil => simp at h
  | cons a rest ih =>
    by_cases ha : lit a = true
    · simp [scan, ha]
    · rw [scan, if_neg ha]
      obtain ⟨i, hmem, hi⟩ := h
      rcases List.mem_cons.mp hmem with rfl | hr
      · exact (ha hi).elim
      · exact ih ⟨i, hr, hi⟩

/-- A covering scan of accurate ballots implements a work-conserving canonical scheduler. -/
theorem scan_workConserving (g : Layout n m) (fallback : ℕ → Fin n)
    (visits : ℕ → List (Fin n)) (covers : ∀ t i, i ∈ visits t)
    (reader : (Fin n → Bool) → Fin n → Bool)
    (accurate : ∀ x i, reader x i = true ↔ Ballot g x i) :
    WorkConserving g.carrier (fun t x => scan (reader x) (fallback t) (visits t)) := by
  rintro t x ⟨i, hi⟩
  apply (ballot_iff_fires g x _).1
  apply (accurate x _).1
  exact scan_lit _ _ _ ⟨i, covers t i, (accurate x i).2 ((ballot_iff_fires g x i).2 hi)⟩

/-- One charged invocation per initially broken comparison suffices for a normal record. -/
theorem measured_service (g : Layout n m) (fallback : ℕ → Fin n)
    (visits : ℕ → List (Fin n)) (covers : ∀ t i, i ∈ visits t)
    (reader : (Fin n → Bool) → Fin n → Bool)
    (accurate : ∀ x i, reader x i = true ↔ Ballot g x i)
    (x : Records g.carrier) :
    let sigma : AdaptiveScheduler g.carrier := fun t y => scan (reader y) (fallback t) (visits t)
    PathwiseWeakFair g.carrier sigma x ∧ BoundedWaste g.carrier 0 sigma ∧
      ReachesNormalWithinAttemptBudget g.carrier sigma x (mismatchCount x) := by
  dsimp only
  have hw := scan_workConserving g fallback visits covers reader accurate
  have hz := workConserving_boundedWaste_zero _ _ hw
  exact ⟨workConserving_pathwiseWeakFair _ _ hw x, hz,
    boundedWaste_reaches_within_attempt_budget _ 0 _ x _ hz (by simp)⟩

end

end OPH.Execution.MeasuredBallots
#print axioms OPH.Execution.MeasuredBallots.canonical_flip
#print axioms OPH.Execution.MeasuredBallots.scan_workConserving
#print axioms OPH.Execution.MeasuredBallots.ballot_after_flip
#print axioms OPH.Execution.MeasuredBallots.measured_service
