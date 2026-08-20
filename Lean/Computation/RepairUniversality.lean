import Mathlib

/-!
# Input-specialized repair satisfiability and width-five permutation programs

Finite observer federations over Boolean registers: each observer holds a
finite patch, an acceptance verdict computed from that patch alone, and a
repair move that rewrites inside the patch.  Main statements:

* `repair_universality` (a historical theorem name) states that for every
  Boolean function `f` and every input `x` there is an input-specialized
  federation whose declared sweep reaches a satisfying state with output
  `f x`;
* `consensus_agree` / `compile_output_unique` — consensus states of a
  compiled federation agree on every program register, so the reached
  consensus output is the unique one;
* `bar_computes` / `bar_length` state that every formula compiles to a
  width-five permutation word program in the ambient group `S₅`, using
  conjugates of a fixed five-cycle, of length at most `4 ^ depth`, whose
  product reads out the formula value.  The file does not type the program
  in `A₅` or prove an even-permutation membership invariant;
* `word_universality` / `echosahedral_universality` — the packaged
  universality statements for both engines;
* kernel `decide` receipts for the SHA-256 primitive functions `Ch` and
  `Maj`, a full-adder sum, and end-to-end federation and word-program runs.

The first construction is input-specialized straight-line constraint
satisfaction, not one fixed repair computer reading its input from an
initial state.  Its `RepairStep` relation records only a failing observer and
an outside-patch frame condition; it permits stuttering and does not require
the selected observer to accept after the step.  The file does not formalize
the canonical OPH `acceptedStep`, the r2020 OMEGA hardware contract, analog
noise margins, convergence rates on unpinned instances, uniform circuit
families, or any speedup claim.
-/

namespace OPH.RepairUniversality

/-! ## States, observers, accepted repair -/

/-- Register states: Boolean records indexed by natural-number registers. -/
abbrev State := ℕ → Bool

/-- One observer: a finite patch, an acceptance verdict reading only the
patch, and the record of that locality. -/
structure Obs where
  patch : List ℕ
  ok : State → Bool
  ok_local : ∀ s t : State, (∀ i ∈ patch, s i = t i) → ok s = ok t

/-- A federation is a finite list of observers. -/
abbrev Federation := List Obs

/-- Consensus: every observer accepts. -/
def Consensus (F : Federation) (s : State) : Prop := ∀ o ∈ F, o.ok s = true

/-- Weak patch step used by this file: some observer fails before the step
and registers outside its patch are unchanged.  No change or post-step
acceptance is required.  This relation is strictly weaker than the canonical
OPH accepted-repair relation. -/
def RepairStep (F : Federation) (s s' : State) : Prop :=
  ∃ o ∈ F, o.ok s = false ∧ ∀ i, i ∉ o.patch → s' i = s i

/-- Explicit boundary control: the weak relation admits a stuttering step at
every state with a failing observer. -/
theorem repairStep_allows_stuttering (F : Federation) (s : State) (o : Obs)
    (ho : o ∈ F) (hfail : o.ok s = false) : RepairStep F s s :=
  ⟨o, ho, hfail, fun _ _ ↦ rfl⟩

/-- Finitely many accepted repair steps. -/
inductive Reaches (F : Federation) : State → State → Prop
  | refl (s : State) : Reaches F s s
  | step {s t u : State} : RepairStep F s t → Reaches F t u → Reaches F s u

/-! ## Straight-line constraint nodes -/

/-- A straight-line constraint node: register `reg` is forced to the value
`val`, computed from the dependency registers alone. -/
structure Node where
  reg : ℕ
  deps : List ℕ
  val : State → Bool
  val_local : ∀ s t : State, (∀ i ∈ deps, s i = t i) → val s = val t

/-- The observer of a node: patch `reg :: deps`, verdict `s reg == val s`. -/
def Node.obs (n : Node) : Obs where
  patch := n.reg :: n.deps
  ok s := s n.reg == n.val s
  ok_local := by
    intro s t h
    have hr : s n.reg = t n.reg := h n.reg (List.mem_cons_self)
    have hv : n.val s = n.val t :=
      n.val_local s t fun i hi => h i (List.mem_cons_of_mem _ hi)
    rw [hr, hv]

/-- Well-formed node lists: registers strictly increase along the list and
every dependency points strictly below its node's register. -/
def NodesWF (L : List Node) : Prop :=
  L.Pairwise (fun a b => a.reg < b.reg) ∧ ∀ n ∈ L, ∀ i ∈ n.deps, i < n.reg

/-- Overwrite one register. -/
def write (s : State) (r : ℕ) (b : Bool) : State := fun i => if i = r then b else s i

/-- The sweep schedule: visit the nodes in order; each node rewrites its
register to the forced value. -/
def sweepFrom (s : State) : List Node → State
  | [] => s
  | n :: L => sweepFrom (write s n.reg (n.val s)) L

theorem sweepFrom_append (s : State) (A B : List Node) :
    sweepFrom s (A ++ B) = sweepFrom (sweepFrom s A) B := by
  induction A generalizing s with
  | nil => rfl
  | cons n L ih =>
      simp only [List.cons_append, sweepFrom]
      exact ih _

theorem sweepFrom_frame (s : State) (L : List Node) (i : ℕ)
    (h : ∀ n ∈ L, n.reg ≠ i) : sweepFrom s L i = s i := by
  induction L generalizing s with
  | nil => rfl
  | cons n L ih =>
      have hi : n.reg ≠ i := h n (List.mem_cons_self)
      rw [sweepFrom, ih _ fun m hm => h m (List.mem_cons_of_mem n hm)]
      simp only [write]
      exact if_neg fun heq => hi heq.symm

theorem sweepFrom_forced (s : State) (L : List Node) (hWF : NodesWF L) :
    ∀ n ∈ L, sweepFrom s L n.reg = n.val (sweepFrom s L) := by
  induction L generalizing s with
  | nil => intro n hn; cases hn
  | cons m L ih =>
      obtain ⟨hpw, hdep⟩ := hWF
      rw [List.pairwise_cons] at hpw
      intro n hn
      rcases List.mem_cons.mp hn with rfl | hn
      · rw [sweepFrom]
        have hfr : sweepFrom (write s n.reg (n.val s)) L n.reg
            = write s n.reg (n.val s) n.reg :=
          sweepFrom_frame _ _ _ fun p hp => ne_of_gt (hpw.1 p hp)
        have hval : n.val (sweepFrom (write s n.reg (n.val s)) L) = n.val s := by
          apply n.val_local
          intro i hi
          have hilt : i < n.reg := hdep n (List.mem_cons_self) i hi
          rw [sweepFrom_frame _ _ _ fun p hp => ne_of_gt (lt_trans hilt (hpw.1 p hp))]
          simp only [write]
          exact if_neg (Nat.ne_of_lt hilt)
        rw [hfr, hval]
        simp [write]
      · rw [sweepFrom]
        exact ih _ ⟨hpw.2, fun p hp => hdep p (List.mem_cons_of_mem m hp)⟩ n hn

/-- The sweep result is a consensus state of the node federation. -/
theorem sweep_consensus (s : State) (L : List Node) (hWF : NodesWF L) :
    Consensus (L.map Node.obs) (sweepFrom s L) := by
  intro o ho
  rcases List.mem_map.mp ho with ⟨n, hn, rfl⟩
  show (sweepFrom s L n.reg == n.val (sweepFrom s L)) = true
  rw [sweepFrom_forced s L hWF n hn]
  exact beq_self_eq_true _

/-- The sweep schedule is a chain of accepted repair steps: every rewrite is
performed by a currently failing observer inside its own patch. -/
theorem sweep_reaches (F : Federation) (s : State) (L : List Node)
    (hsub : ∀ n ∈ L, n.obs ∈ F) : Reaches F s (sweepFrom s L) := by
  induction L generalizing s with
  | nil => exact Reaches.refl s
  | cons n L ih =>
      rw [sweepFrom]
      by_cases h : n.obs.ok s = true
      · have h' : s n.reg = n.val s := by
          simpa [Node.obs, beq_iff_eq] using h
        have hs : write s n.reg (n.val s) = s := by
          funext i
          simp only [write]
          split
          · next hi => rw [hi]; exact h'.symm
          · rfl
        rw [hs]
        exact ih _ fun m hm => hsub m (List.mem_cons_of_mem n hm)
      · have hfalse : n.obs.ok s = false := by simpa using h
        have hstep : RepairStep F s (write s n.reg (n.val s)) := by
          refine ⟨n.obs, hsub n (List.mem_cons_self), hfalse, ?_⟩
          intro i hi
          simp only [write]
          refine if_neg fun heq => hi ?_
          rw [heq]
          exact List.mem_cons_self
        exact Reaches.step hstep (ih _ fun m hm => hsub m (List.mem_cons_of_mem n hm))

/-- Consensus states satisfy every node constraint. -/
theorem consensus_forced {L : List Node} {s : State}
    (hs : Consensus (L.map Node.obs) s) : ∀ n ∈ L, s n.reg = n.val s := by
  intro n hn
  have h := hs n.obs (List.mem_map.mpr ⟨n, hn, rfl⟩)
  simpa [Node.obs, beq_iff_eq] using h

/-- Two consensus states of a covered well-formed program agree on every
program register: the consensus is unique where it matters. -/
theorem consensus_agree (L : List Node) (hWF : NodesWF L)
    (hcov : ∀ n ∈ L, ∀ i ∈ n.deps, ∃ m ∈ L, m.reg = i)
    (s t : State) (hs : Consensus (L.map Node.obs) s)
    (ht : Consensus (L.map Node.obs) t) :
    ∀ n ∈ L, s n.reg = t n.reg := by
  have key : ∀ r : ℕ, ∀ n ∈ L, n.reg = r → s n.reg = t n.reg := by
    intro r
    induction r using Nat.strong_induction_on with
    | _ r ih =>
      intro n hn hr
      rw [consensus_forced hs n hn, consensus_forced ht n hn]
      apply n.val_local
      intro i hi
      obtain ⟨m, hm, hmr⟩ := hcov n hn i hi
      have hlt : i < n.reg := hWF.2 n hn i hi
      have h2 := ih i (hr ▸ hlt) m hm hmr
      rwa [hmr] at h2
  exact fun n hn => key n.reg n hn rfl

/-! ## NAND formulas -/

/-- NAND formulas with variables and pinned constants. -/
inductive Formula (k : ℕ) where
  | var (i : Fin k)
  | const (b : Bool)
  | nand (φ ψ : Formula k)

namespace Formula

variable {k : ℕ}

def evalF : Formula k → (Fin k → Bool) → Bool
  | var i, x => x i
  | const b, _ => b
  | nand φ ψ, x => !(evalF φ x && evalF ψ x)

def depth : Formula k → ℕ
  | var _ => 0
  | const _ => 0
  | nand φ ψ => max (depth φ) (depth ψ) + 1

def mapVar {m : ℕ} (h : Fin k → Fin m) : Formula k → Formula m
  | var i => var (h i)
  | const b => const b
  | nand φ ψ => nand (mapVar h φ) (mapVar h ψ)

theorem evalF_mapVar {m : ℕ} (h : Fin k → Fin m) (φ : Formula k) (x : Fin m → Bool) :
    evalF (mapVar h φ) x = evalF φ (x ∘ h) := by
  induction φ with
  | var i => rfl
  | const b => rfl
  | nand φ ψ ihφ ihψ => simp only [mapVar, evalF, ihφ, ihψ]

/-- Shannon expansion: every Boolean function on finitely many bits is the
value function of a NAND formula. -/
theorem exists_formula : ∀ {k : ℕ} (f : (Fin k → Bool) → Bool),
    ∃ φ : Formula k, ∀ x, evalF φ x = f x := by
  intro k
  induction k with
  | zero =>
      intro f
      refine ⟨const (f fun i => i.elim0), fun x => ?_⟩
      have hx : x = fun i => i.elim0 := funext fun i => i.elim0
      simp only [evalF]
      rw [hx]
  | succ k ihk =>
      intro f
      obtain ⟨φ₀, h₀⟩ := ihk fun y => f (Fin.cons false y)
      obtain ⟨φ₁, h₁⟩ := ihk fun y => f (Fin.cons true y)
      refine ⟨nand (nand (var 0) (mapVar Fin.succ φ₁))
        (nand (nand (var 0) (var 0)) (mapVar Fin.succ φ₀)), fun x => ?_⟩
      have e₀ : evalF (mapVar Fin.succ φ₀) x = f (Fin.cons false (Fin.tail x)) := by
        rw [evalF_mapVar]
        exact h₀ (Fin.tail x)
      have e₁ : evalF (mapVar Fin.succ φ₁) x = f (Fin.cons true (Fin.tail x)) := by
        rw [evalF_mapVar]
        exact h₁ (Fin.tail x)
      have hx : f (Fin.cons (x 0) (Fin.tail x)) = f x := by
        rw [Fin.cons_self_tail]
      cases hx0 : x 0 with
      | false =>
          rw [← hx, hx0]
          simp [evalF, e₀, hx0]
      | true =>
          rw [← hx, hx0]
          simp [evalF, e₁, hx0]

end Formula

/-! ## Compiler: formula → straight-line federation -/

/-- Compiler output: emitted nodes, output register, next free register. -/
structure Out where
  nodes : List Node
  out : ℕ
  next : ℕ

variable {k : ℕ}

/-- Compile the interior of a formula, allocating registers from `c`.
Input registers `0, …, k-1` are supplied by `inputs`. -/
def go : Formula k → ℕ → Out
  | .var i, c =>
      ⟨[⟨c, [i.val], fun s => s i.val, fun _ _ h => h i.val (by simp)⟩], c, c + 1⟩
  | .const b, c => ⟨[⟨c, [], fun _ => b, fun _ _ _ => rfl⟩], c, c + 1⟩
  | .nand φ ψ, c =>
      ⟨(go φ c).nodes ++ ((go ψ (go φ c).next).nodes ++
          [⟨(go ψ (go φ c).next).next, [(go φ c).out, (go ψ (go φ c).next).out],
            fun s => !(s (go φ c).out && s (go ψ (go φ c).next).out),
            fun s t h => by
              dsimp only
              rw [h (go φ c).out (by simp), h (go ψ (go φ c).next).out (by simp)]⟩]),
        (go ψ (go φ c).next).next, (go ψ (go φ c).next).next + 1⟩

/-- The facts the compiler guarantees, relative to allocation start `c`. -/
structure GoFacts (k c : ℕ) (φ : Formula k) (x : Fin k → Bool) (o : Out) : Prop where
  c_le_out : c ≤ o.out
  out_lt_next : o.out < o.next
  range : ∀ n ∈ o.nodes, c ≤ n.reg ∧ n.reg < o.next
  pairwise : o.nodes.Pairwise fun a b => a.reg < b.reg
  deps_lt : ∀ n ∈ o.nodes, ∀ i ∈ n.deps, i < n.reg
  out_mem : o.out ∈ o.nodes.map Node.reg
  deps_cov : ∀ n ∈ o.nodes, ∀ i ∈ n.deps, i < k ∨ i ∈ o.nodes.map Node.reg
  sweep_out : ∀ s : State, (∀ j : Fin k, s j.val = x j) →
    sweepFrom s o.nodes o.out = Formula.evalF φ x

theorem go_facts (φ : Formula k) : ∀ (c : ℕ), k ≤ c → ∀ x : Fin k → Bool,
    GoFacts k c φ x (go φ c) := by
  induction φ with
  | var i =>
      intro c hkc x
      refine ⟨le_refl c, Nat.lt_succ_self c, ?_, ?_, ?_, ?_, ?_, ?_⟩
      · intro n hn
        rcases List.mem_singleton.mp hn with rfl
        exact ⟨le_refl c, Nat.lt_succ_self c⟩
      · simp [go]
      · intro n hn j hj
        rcases List.mem_singleton.mp hn with rfl
        rcases List.mem_singleton.mp hj with rfl
        exact lt_of_lt_of_le i.isLt hkc
      · simp [go]
      · intro n hn j hj
        rcases List.mem_singleton.mp hn with rfl
        rcases List.mem_singleton.mp hj with rfl
        exact Or.inl i.isLt
      · intro s hs
        simp [go, sweepFrom, write, Formula.evalF]
        exact hs i
  | const b =>
      intro c hkc x
      refine ⟨le_refl c, Nat.lt_succ_self c, ?_, ?_, ?_, ?_, ?_, ?_⟩
      · intro n hn
        rcases List.mem_singleton.mp hn with rfl
        exact ⟨le_refl c, Nat.lt_succ_self c⟩
      · simp [go]
      · intro n hn j hj
        rcases List.mem_singleton.mp hn with rfl
        cases hj
      · simp [go]
      · intro n hn j hj
        rcases List.mem_singleton.mp hn with rfl
        cases hj
      · intro s _
        simp [go, sweepFrom, write, Formula.evalF]
  | nand φ ψ ihφ ihψ =>
      intro c hkc x
      have Fφ := ihφ c hkc x
      have hkc' : k ≤ (go φ c).next :=
        le_trans hkc (le_trans Fφ.c_le_out (le_of_lt Fφ.out_lt_next))
      have Fψ := ihψ (go φ c).next hkc' x
      have hφψ : (go φ c).next ≤ (go ψ (go φ c).next).next :=
        le_trans Fψ.c_le_out (le_of_lt Fψ.out_lt_next)
      have houtφ : (go φ c).out < (go ψ (go φ c).next).next :=
        lt_of_lt_of_le Fφ.out_lt_next hφψ
      have hcnext : c ≤ (go ψ (go φ c).next).next :=
        le_trans (le_trans Fφ.c_le_out (le_of_lt Fφ.out_lt_next)) hφψ
      refine ⟨hcnext, Nat.lt_succ_self _, ?_, ?_, ?_, ?_, ?_, ?_⟩
      · intro n hn
        simp only [go, List.mem_append, List.mem_singleton] at hn
        rcases hn with hn | hn | rfl
        · rcases Fφ.range n hn with ⟨h1, h2⟩
          exact ⟨h1, lt_trans (lt_of_lt_of_le h2 hφψ) (Nat.lt_succ_self _)⟩
        · rcases Fψ.range n hn with ⟨h1, h2⟩
          exact ⟨le_trans (le_trans Fφ.c_le_out (le_of_lt Fφ.out_lt_next)) h1,
            lt_trans h2 (Nat.lt_succ_self _)⟩
        · exact ⟨hcnext, Nat.lt_succ_self _⟩
      · simp only [go]
        rw [List.pairwise_append]
        refine ⟨Fφ.pairwise, ?_, ?_⟩
        · rw [List.pairwise_append]
          refine ⟨Fψ.pairwise, by simp, ?_⟩
          intro a ha b hb
          rcases List.mem_singleton.mp hb with rfl
          exact (Fψ.range a ha).2
        · intro a ha b hb
          have haφ : a.reg < (go φ c).next := (Fφ.range a ha).2
          simp only [List.mem_append, List.mem_singleton] at hb
          rcases hb with hb | rfl
          · exact lt_of_lt_of_le haφ (Fψ.range b hb).1
          · exact lt_of_lt_of_le haφ hφψ
      · intro n hn j hj
        simp only [go, List.mem_append, List.mem_singleton] at hn
        rcases hn with hn | hn | rfl
        · exact Fφ.deps_lt n hn j hj
        · exact Fψ.deps_lt n hn j hj
        · simp only [List.mem_cons, List.not_mem_nil, or_false] at hj
          rcases hj with rfl | rfl
          · exact houtφ
          · exact Fψ.out_lt_next
      · simp only [go, List.map_append, List.mem_append]
        right; right
        simp
      · intro n hn j hj
        simp only [go, List.mem_append, List.mem_singleton] at hn
        rcases hn with hn | hn | rfl
        · rcases Fφ.deps_cov n hn j hj with h | h
          · exact Or.inl h
          · refine Or.inr ?_
            simp only [go, List.map_append, List.mem_append]
            exact Or.inl h
        · rcases Fψ.deps_cov n hn j hj with h | h
          · exact Or.inl h
          · refine Or.inr ?_
            simp only [go, List.map_append, List.mem_append]
            exact Or.inr (Or.inl h)
        · simp only [List.mem_cons, List.not_mem_nil, or_false] at hj
          rcases hj with rfl | rfl
          · refine Or.inr ?_
            simp only [go, List.map_append, List.mem_append]
            exact Or.inl Fφ.out_mem
          · refine Or.inr ?_
            simp only [go, List.map_append, List.mem_append]
            exact Or.inr (Or.inl Fψ.out_mem)
      · intro s hs
        simp only [go, sweepFrom_append]
        have hs₁ : ∀ j : Fin k, sweepFrom s (go φ c).nodes j.val = x j := by
          intro j
          rw [sweepFrom_frame _ _ _ fun p hp =>
            ne_of_gt (lt_of_lt_of_le j.isLt (le_trans hkc (Fφ.range p hp).1))]
          exact hs j
        have hφout : sweepFrom (sweepFrom s (go φ c).nodes)
            (go ψ (go φ c).next).nodes (go φ c).out = Formula.evalF φ x := by
          rw [sweepFrom_frame _ _ _ fun p hp =>
            ne_of_gt (lt_of_lt_of_le Fφ.out_lt_next (Fψ.range p hp).1)]
          exact Fφ.sweep_out s hs
        have hψout : sweepFrom (sweepFrom s (go φ c).nodes)
            (go ψ (go φ c).next).nodes (go ψ (go φ c).next).out = Formula.evalF ψ x :=
          Fψ.sweep_out _ hs₁
        simp only [sweepFrom, write, if_true]
        rw [hφout, hψout]
        rfl

/-! ## Input pinning and full compilation -/

/-- Input nodes for registers `0, …, m-1`, pinning register `j` to `x j`. -/
def buildInputs (x : Fin k → Bool) : (m : ℕ) → m ≤ k → List Node
  | 0, _ => []
  | m + 1, h =>
      buildInputs x m (Nat.le_of_succ_le h) ++
        [⟨m, [], fun _ => x ⟨m, h⟩, fun _ _ _ => rfl⟩]

theorem buildInputs_reg_lt (x : Fin k → Bool) :
    ∀ (m : ℕ) (h : m ≤ k), ∀ n ∈ buildInputs x m h, n.reg < m := by
  intro m
  induction m with
  | zero => intro h n hn; cases hn
  | succ m ih =>
      intro h n hn
      simp only [buildInputs, List.mem_append, List.mem_singleton] at hn
      rcases hn with hn | rfl
      · exact lt_trans (ih _ n hn) (Nat.lt_succ_self m)
      · exact Nat.lt_succ_self m

theorem buildInputs_deps (x : Fin k → Bool) :
    ∀ (m : ℕ) (h : m ≤ k), ∀ n ∈ buildInputs x m h, n.deps = [] := by
  intro m
  induction m with
  | zero => intro h n hn; cases hn
  | succ m ih =>
      intro h n hn
      simp only [buildInputs, List.mem_append, List.mem_singleton] at hn
      rcases hn with hn | rfl
      · exact ih _ n hn
      · rfl

theorem buildInputs_pairwise (x : Fin k → Bool) :
    ∀ (m : ℕ) (h : m ≤ k), (buildInputs x m h).Pairwise fun a b => a.reg < b.reg := by
  intro m
  induction m with
  | zero => intro h; exact List.Pairwise.nil
  | succ m ih =>
      intro h
      rw [buildInputs, List.pairwise_append]
      refine ⟨ih _, by simp, ?_⟩
      intro a ha b hb
      rcases List.mem_singleton.mp hb with rfl
      exact buildInputs_reg_lt x m _ a ha

theorem buildInputs_covers (x : Fin k → Bool) :
    ∀ (m : ℕ) (h : m ≤ k), ∀ i, i < m → ∃ n ∈ buildInputs x m h, n.reg = i := by
  intro m
  induction m with
  | zero => intro h i hi; exact absurd hi (Nat.not_lt_zero i)
  | succ m ih =>
      intro h i hi
      rcases Nat.lt_succ_iff_lt_or_eq.mp hi with hi | rfl
      · obtain ⟨n, hn, hnr⟩ := ih _ i hi
        exact ⟨n, List.mem_append.mpr (Or.inl hn), hnr⟩
      · refine ⟨_, List.mem_append.mpr (Or.inr (List.mem_singleton.mpr rfl)), rfl⟩

theorem buildInputs_sweep (x : Fin k → Bool) (s : State) :
    ∀ (m : ℕ) (h : m ≤ k) (j : Fin k), j.val < m →
      sweepFrom s (buildInputs x m h) j.val = x j := by
  intro m
  induction m with
  | zero => intro h j hj; exact absurd hj (Nat.not_lt_zero _)
  | succ m ih =>
      intro h j hj
      rw [buildInputs, sweepFrom_append]
      rcases Nat.lt_succ_iff_lt_or_eq.mp hj with hj | hj
      · simp only [sweepFrom, write]
        rw [if_neg (Nat.ne_of_lt hj)]
        exact ih _ j hj
      · simp only [sweepFrom, write]
        rw [if_pos hj]
        exact congrArg x (Fin.ext hj.symm)

/-- The input federation: registers `0, …, k-1` pinned to `x`. -/
def inputs (x : Fin k → Bool) : List Node := buildInputs x k le_rfl

/-- The compiled program: pinned inputs followed by the formula circuit. -/
def compile (φ : Formula k) (x : Fin k → Bool) : List Node :=
  inputs x ++ (go φ k).nodes

/-- The output register of the compiled program. -/
def outReg (φ : Formula k) : ℕ := (go φ k).out

/-- The compiled observer federation. -/
def federation (φ : Formula k) (x : Fin k → Bool) : Federation :=
  (compile φ x).map Node.obs

theorem compile_wf (φ : Formula k) (x : Fin k → Bool) : NodesWF (compile φ x) := by
  have F := go_facts φ k le_rfl x
  constructor
  · rw [compile, List.pairwise_append]
    refine ⟨buildInputs_pairwise x k le_rfl, F.pairwise, ?_⟩
    intro a ha b hb
    exact lt_of_lt_of_le (buildInputs_reg_lt x k le_rfl a ha) (F.range b hb).1
  · intro n hn i hi
    simp only [compile] at hn
    rcases List.mem_append.mp hn with hn | hn
    · rw [buildInputs_deps x k le_rfl n hn] at hi
      cases hi
    · exact F.deps_lt n hn i hi

theorem inputs_sweep (x : Fin k → Bool) (s : State) (j : Fin k) :
    sweepFrom s (inputs x) j.val = x j :=
  buildInputs_sweep x s k le_rfl j j.isLt

theorem compile_sweep_out (φ : Formula k) (x : Fin k → Bool) (s₀ : State) :
    sweepFrom s₀ (compile φ x) (outReg φ) = Formula.evalF φ x := by
  rw [compile, sweepFrom_append]
  exact (go_facts φ k le_rfl x).sweep_out _ fun j => inputs_sweep x s₀ j

theorem compile_consensus (φ : Formula k) (x : Fin k → Bool) (s₀ : State) :
    Consensus (federation φ x) (sweepFrom s₀ (compile φ x)) :=
  sweep_consensus s₀ _ (compile_wf φ x)

theorem compile_reaches (φ : Formula k) (x : Fin k → Bool) (s₀ : State) :
    Reaches (federation φ x) s₀ (sweepFrom s₀ (compile φ x)) :=
  sweep_reaches _ s₀ _ fun n hn => List.mem_map.mpr ⟨n, hn, rfl⟩

theorem compile_cov (φ : Formula k) (x : Fin k → Bool) :
    ∀ n ∈ compile φ x, ∀ i ∈ n.deps, ∃ m ∈ compile φ x, m.reg = i := by
  intro n hn i hi
  have F := go_facts φ k le_rfl x
  simp only [compile] at hn ⊢
  rcases List.mem_append.mp hn with hn | hn
  · rw [buildInputs_deps x k le_rfl n hn] at hi
    cases hi
  · rcases F.deps_cov n hn i hi with h | h
    · obtain ⟨m, hm, hmr⟩ := buildInputs_covers x k le_rfl i h
      exact ⟨m, List.mem_append.mpr (Or.inl hm), hmr⟩
    · obtain ⟨m, hm, hmr⟩ := List.mem_map.mp h
      exact ⟨m, List.mem_append.mpr (Or.inr hm), hmr⟩

/-- Consensus output uniqueness for compiled federations. -/
theorem compile_output_unique (φ : Formula k) (x : Fin k → Bool) (s t : State)
    (hs : Consensus (federation φ x) s) (ht : Consensus (federation φ x) t) :
    s (outReg φ) = t (outReg φ) := by
  have F := go_facts φ k le_rfl x
  obtain ⟨m, hm, hmr⟩ := List.mem_map.mp F.out_mem
  have h := consensus_agree (compile φ x) (compile_wf φ x) (compile_cov φ x) s t hs ht
    m (List.mem_append.mpr (Or.inr hm))
  rw [outReg, ← hmr]
  exact h

/-- Input-specialized satisfiability theorem (historical name).  For every
Boolean function and input, the construction hardwires that input into a
different federation and its weak-step sweep reaches a satisfying state with
the required output.  This does not construct one fixed federation that reads
immutable input registers from the initial state, and it does not use the
canonical OPH accepted-repair relation. -/
theorem repair_universality (k : ℕ) (f : (Fin k → Bool) → Bool) :
    ∃ φ : Formula k, ∀ (x : Fin k → Bool) (s₀ : State),
      Reaches (federation φ x) s₀ (sweepFrom s₀ (compile φ x)) ∧
      Consensus (federation φ x) (sweepFrom s₀ (compile φ x)) ∧
      sweepFrom s₀ (compile φ x) (outReg φ) = f x := by
  obtain ⟨φ, hφ⟩ := Formula.exists_formula f
  refine ⟨φ, fun x s₀ => ⟨compile_reaches φ x s₀, compile_consensus φ x s₀, ?_⟩⟩
  rw [compile_sweep_out, hφ]

/-! ## Width-five word programs in the ambient symmetric group -/

/-- Permutations of five letters; the ambient group of the word programs. -/
abbrev P5 := Equiv.Perm (Fin 5)

/-- The reference five-cycle `(0 1 2 3 4)`. -/
def σ₀ : P5 :=
  ⟨![1, 2, 3, 4, 0], ![4, 0, 1, 2, 3], by intro i; fin_cases i <;> rfl,
    by intro i; fin_cases i <;> rfl⟩

/-- First conjugator of the commutator witness. -/
def τa : P5 :=
  ⟨![0, 1, 3, 4, 2], ![0, 1, 4, 2, 3], by intro i; fin_cases i <;> rfl,
    by intro i; fin_cases i <;> rfl⟩

/-- Second conjugator of the commutator witness. -/
def τb : P5 :=
  ⟨![0, 3, 2, 4, 1], ![0, 4, 2, 1, 3], by intro i; fin_cases i <;> rfl,
    by intro i; fin_cases i <;> rfl⟩

/-- The kernel-checked commutator witness: the two conjugates of the
reference five-cycle have commutator `σ₀⁻¹`.  This single identity carries
the whole Barrington induction. -/
theorem comm_witness :
    (τa * σ₀ * τa⁻¹) * (τb * σ₀ * τb⁻¹) * (τa * σ₀ * τa⁻¹)⁻¹ * (τb * σ₀ * τb⁻¹)⁻¹ = σ₀⁻¹ :=
  Equiv.ext fun i => by fin_cases i <;> rfl

/-- One program instruction: read one input bit, emit one of two group
elements. -/
structure Instr (k : ℕ) where
  idx : Fin k
  g0 : P5
  g1 : P5

/-- The group element an instruction contributes on input `x`. -/
def Instr.act (ins : Instr k) (x : Fin k → Bool) : P5 :=
  if x ins.idx then ins.g1 else ins.g0

/-- Run a program: multiply the contributed elements in order. -/
def runP : List (Instr k) → (Fin k → Bool) → P5
  | [], _ => 1
  | ins :: P, x => ins.act x * runP P x

/-- A program computes `f` with target `σ` when its product is `σ` on
accepting inputs and the identity otherwise. -/
def Computes (P : List (Instr k)) (f : (Fin k → Bool) → Bool) (σ : P5) : Prop :=
  ∀ x, runP P x = if f x then σ else 1

theorem runP_append (A B : List (Instr k)) (x : Fin k → Bool) :
    runP (A ++ B) x = runP A x * runP B x := by
  induction A with
  | nil => simp [runP]
  | cons ins A ih => simp [runP, ih, mul_assoc]

/-- Reverse the program and invert every entry: the run inverts. -/
def invP (P : List (Instr k)) : List (Instr k) :=
  (P.map fun ins => ⟨ins.idx, ins.g0⁻¹, ins.g1⁻¹⟩).reverse

theorem invP_length (P : List (Instr k)) : (invP P).length = P.length := by
  simp [invP]

theorem runP_invP (P : List (Instr k)) (x : Fin k → Bool) :
    runP (invP P) x = (runP P x)⁻¹ := by
  induction P with
  | nil => simp [invP, runP]
  | cons ins P ih =>
      have hrw : invP (ins :: P) = invP P ++ [⟨ins.idx, ins.g0⁻¹, ins.g1⁻¹⟩] := by
        simp [invP]
      rw [hrw, runP_append, ih]
      simp only [runP, Instr.act, mul_inv_rev]
      cases hx : x ins.idx <;> simp

/-- Fuse a constant right factor into the final instruction; the length is
unchanged. -/
def fuse (g : P5) : List (Instr k) → List (Instr k)
  | [] => []
  | [ins] => [⟨ins.idx, ins.g0 * g, ins.g1 * g⟩]
  | ins :: P => ins :: fuse g P

theorem fuse_length (g : P5) (P : List (Instr k)) : (fuse g P).length = P.length := by
  induction P with
  | nil => rfl
  | cons ins P ih =>
      cases P with
      | nil => rfl
      | cons ins' P' =>
          show (ins :: fuse g (ins' :: P')).length = _
          simp only [List.length_cons] at ih ⊢
          omega

theorem runP_fuse (g : P5) (P : List (Instr k)) (hP : P ≠ []) (x : Fin k → Bool) :
    runP (fuse g P) x = runP P x * g := by
  induction P with
  | nil => cases hP rfl
  | cons ins P ih =>
      cases P with
      | nil =>
          simp only [fuse, runP, Instr.act]
          cases hx : x ins.idx <;> simp
      | cons ins' P' =>
          show runP (ins :: fuse g (ins' :: P')) x = _
          simp only [runP, ih (by simp), mul_assoc]

section Bar

variable [NeZero k]

/-- Barrington compilation: `bar φ τ` computes `φ` with target
`τ * σ₀ * τ⁻¹`.  The NAND case runs the two branch programs at shifted
conjugators, undoes them, and fuses the target into the last instruction. -/
def bar : Formula k → P5 → List (Instr k)
  | .var i, τ => [⟨i, 1, τ * σ₀ * τ⁻¹⟩]
  | .const b, τ =>
      [⟨0, if b then τ * σ₀ * τ⁻¹ else 1, if b then τ * σ₀ * τ⁻¹ else 1⟩]
  | .nand φ ψ, τ =>
      fuse (τ * σ₀ * τ⁻¹)
        (bar φ (τ * τa) ++ (bar ψ (τ * τb) ++
          (invP (bar φ (τ * τa)) ++ invP (bar ψ (τ * τb)))))

theorem bar_length_pos (φ : Formula k) (τ : P5) : 0 < (bar φ τ).length := by
  induction φ generalizing τ with
  | var i => simp [bar]
  | const b => simp [bar]
  | nand φ ψ ihφ ihψ =>
      simp only [bar, fuse_length, List.length_append]
      have := ihφ (τ * τa)
      omega

private theorem conj_key (τ : P5) :
    ((τ * τa) * σ₀ * (τ * τa)⁻¹) * (((τ * τb) * σ₀ * (τ * τb)⁻¹) *
      (((τ * τa) * σ₀ * (τ * τa)⁻¹)⁻¹ * ((τ * τb) * σ₀ * (τ * τb)⁻¹)⁻¹)) *
      (τ * σ₀ * τ⁻¹) = 1 := by
  have expand : ((τ * τa) * σ₀ * (τ * τa)⁻¹) * (((τ * τb) * σ₀ * (τ * τb)⁻¹) *
      (((τ * τa) * σ₀ * (τ * τa)⁻¹)⁻¹ * ((τ * τb) * σ₀ * (τ * τb)⁻¹)⁻¹)) *
      (τ * σ₀ * τ⁻¹)
      = τ * (((τa * σ₀ * τa⁻¹) * (τb * σ₀ * τb⁻¹) * (τa * σ₀ * τa⁻¹)⁻¹ *
          (τb * σ₀ * τb⁻¹)⁻¹) * σ₀) * τ⁻¹ := by
    group
  rw [expand, comm_witness]
  group

theorem bar_computes (φ : Formula k) (τ : P5) :
    Computes (bar φ τ) (fun x => Formula.evalF φ x) (τ * σ₀ * τ⁻¹) := by
  induction φ generalizing τ with
  | var i =>
      intro x
      cases hx : x i <;> simp [bar, runP, Instr.act, Formula.evalF, hx]
  | const b =>
      intro x
      cases b <;> simp [bar, runP, Instr.act, Formula.evalF]
  | nand φ ψ ihφ ihψ =>
      intro x
      have hlen : 0 < (bar φ (τ * τa) ++ (bar ψ (τ * τb) ++
          (invP (bar φ (τ * τa)) ++ invP (bar ψ (τ * τb))))).length := by
        rw [List.length_append]
        exact Nat.lt_of_lt_of_le (bar_length_pos φ (τ * τa)) (Nat.le_add_right _ _)
      have hne : (bar φ (τ * τa) ++ (bar ψ (τ * τb) ++
          (invP (bar φ (τ * τa)) ++ invP (bar ψ (τ * τb))))) ≠ [] := by
        intro hnil
        rw [hnil] at hlen
        simp at hlen
      simp only [bar]
      rw [runP_fuse _ _ hne, runP_append, runP_append,
        runP_append, runP_invP, runP_invP, ihφ (τ * τa) x, ihψ (τ * τb) x]
      simp only [Formula.evalF]
      rcases Bool.eq_false_or_eq_true (Formula.evalF φ x) with h1 | h1 <;>
        rcases Bool.eq_false_or_eq_true (Formula.evalF ψ x) with h2 | h2 <;>
        simp only [h1, h2, Bool.true_and, Bool.false_and, Bool.not_true,
          Bool.not_false, reduceIte, if_neg Bool.false_ne_true, inv_one]
      · exact conj_key τ
      · group
      · group
      · group

/-- Barrington length bound: `4 ^ depth`. -/
theorem bar_length (φ : Formula k) (τ : P5) :
    (bar φ τ).length ≤ 4 ^ Formula.depth φ := by
  induction φ generalizing τ with
  | var i => simp [bar, Formula.depth]
  | const b => simp [bar, Formula.depth]
  | nand φ ψ ihφ ihψ =>
      have h1 := ihφ (τ * τa)
      have h2 := ihψ (τ * τb)
      have h3 : (4 : ℕ) ^ Formula.depth φ ≤ 4 ^ max (Formula.depth φ) (Formula.depth ψ) :=
        Nat.pow_le_pow_right (by norm_num) (le_max_left _ _)
      have h4 : (4 : ℕ) ^ Formula.depth ψ ≤ 4 ^ max (Formula.depth φ) (Formula.depth ψ) :=
        Nat.pow_le_pow_right (by norm_num) (le_max_right _ _)
      have h5 := invP_length (bar φ (τ * τa))
      have h6 := invP_length (bar ψ (τ * τb))
      simp only [bar, fuse_length, List.length_append, Formula.depth, pow_succ]
      omega

/-- Word universality: every Boolean function is computed by a width-five
permutation word program with target the reference five-cycle. -/
theorem word_universality (k : ℕ) [NeZero k] (f : (Fin k → Bool) → Bool) :
    ∃ P : List (Instr k), Computes P f σ₀ := by
  obtain ⟨φ, hφ⟩ := Formula.exists_formula f
  refine ⟨bar φ 1, fun x => ?_⟩
  have h := bar_computes φ 1 x
  simpa [hφ x] using h

/-- Echosahedral universality: every finite Boolean function is computed
both by the accepted-repair dynamics of an observer federation and by a
width-five permutation word program over conjugates of the reference
five-cycle. -/
theorem echosahedral_universality (k : ℕ) [NeZero k] (f : (Fin k → Bool) → Bool) :
    (∃ φ : Formula k, ∀ (x : Fin k → Bool) (s₀ : State),
        Reaches (federation φ x) s₀ (sweepFrom s₀ (compile φ x)) ∧
        Consensus (federation φ x) (sweepFrom s₀ (compile φ x)) ∧
        sweepFrom s₀ (compile φ x) (outReg φ) = f x) ∧
      ∃ P : List (Instr k), Computes P f σ₀ :=
  ⟨repair_universality k f, word_universality k f⟩

end Bar

/-! ## Derived connectives and SHA-256 primitive receipts -/

def notF (φ : Formula k) : Formula k := .nand φ φ
def andF (φ ψ : Formula k) : Formula k := notF (.nand φ ψ)
def orF (φ ψ : Formula k) : Formula k := .nand (notF φ) (notF ψ)
def xorF (φ ψ : Formula k) : Formula k := andF (orF φ ψ) (.nand φ ψ)

/-- SHA-256 choice function as a NAND formula. -/
def ChF : Formula 3 := xorF (andF (.var 0) (.var 1)) (andF (notF (.var 0)) (.var 2))

/-- SHA-256 majority function as a NAND formula. -/
def MajF : Formula 3 :=
  xorF (andF (.var 0) (.var 1)) (xorF (andF (.var 0) (.var 2)) (andF (.var 1) (.var 2)))

/-- Full-adder sum as a NAND formula. -/
def AdderSumF : Formula 3 := xorF (.var 0) (xorF (.var 1) (.var 2))

theorem ChF_spec : ∀ x : Fin 3 → Bool,
    Formula.evalF ChF x = ((x 0 && x 1) ^^ (!(x 0) && x 2)) := by decide

theorem MajF_spec : ∀ x : Fin 3 → Bool,
    Formula.evalF MajF x = ((x 0 && x 1) ^^ ((x 0 && x 2) ^^ (x 1 && x 2))) := by decide

theorem AdderSumF_spec : ∀ x : Fin 3 → Bool,
    Formula.evalF AdderSumF x = (x 0 ^^ (x 1 ^^ x 2)) := by decide

/-- The federation compiled from `Ch` repairs, from any initial state, to
the SHA-256 choice value. -/
theorem ch_repair_receipt (x : Fin 3 → Bool) (s₀ : State) :
    sweepFrom s₀ (compile ChF x) (outReg ChF) = ((x 0 && x 1) ^^ (!(x 0) && x 2)) := by
  rw [compile_sweep_out, ChF_spec]

/-- The federation compiled from `Maj` repairs, from any initial state, to
the SHA-256 majority value. -/
theorem maj_repair_receipt (x : Fin 3 → Bool) (s₀ : State) :
    sweepFrom s₀ (compile MajF x) (outReg MajF)
      = ((x 0 && x 1) ^^ ((x 0 && x 2) ^^ (x 1 && x 2))) := by
  rw [compile_sweep_out, MajF_spec]

/-- The federation compiled from the full-adder sum repairs, from any
initial state, to the sum bit. -/
theorem adder_repair_receipt (x : Fin 3 → Bool) (s₀ : State) :
    sweepFrom s₀ (compile AdderSumF x) (outReg AdderSumF) = (x 0 ^^ (x 1 ^^ x 2)) := by
  rw [compile_sweep_out, AdderSumF_spec]

/-- Kernel receipt: end-to-end sweep of the compiled two-bit NAND
federation from the all-false state. -/
theorem nand_kernel_receipt : ∀ x : Fin 2 → Bool,
    sweepFrom (fun _ => false) (compile (.nand (.var 0) (.var 1)) x)
      (outReg (k := 2) (.nand (.var 0) (.var 1))) = !(x 0 && x 1) := by decide

/-- Kernel receipt: end-to-end word-program run for two-bit NAND,
pointwise on the five letters. -/
theorem nand_word_kernel_receipt : ∀ x : Fin 2 → Bool, ∀ i : Fin 5,
    runP (bar (k := 2) (.nand (.var 0) (.var 1)) 1) x i
      = (if !(x 0 && x 1) then σ₀ else 1) i := by decide

/-! ## Axiom audit -/

#print axioms repair_universality
#print axioms repairStep_allows_stuttering
#print axioms compile_output_unique
#print axioms bar_computes
#print axioms bar_length
#print axioms word_universality
#print axioms echosahedral_universality
#print axioms comm_witness
#print axioms nand_kernel_receipt
#print axioms nand_word_kernel_receipt

end OPH.RepairUniversality
