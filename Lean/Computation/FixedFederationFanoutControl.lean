import Computation.FixedFederationExecutionExamples

/-!
# Fanout control: the single-consumer hypothesis is load-bearing

`FixedFederationComplexity` proves a triangular accepted-step bound for the
formula compiler output, and that proof consumes
`AtMostOneDownstreamConsumer`: every generated register is read by at most one
later node.  This module supplies the typed negative control for that
hypothesis.

`fanoutChain n` is a well-formed straight-line node list in which register `k`
is read by every later node.  Each node is affine in its dependencies, so
flipping any single upstream register flips the forced value of every node
above it.  A bottom-up recursive schedule then produces a genuine canonical
accepted path of length `2 ^ n - 1` from `allFalse` to a consensus state.

The generic bound `canonicalAcceptedSteps_lt_pow` says every well-formed node
list admits fewer than `2 ^ L.length` canonical accepted steps.  The chain
family attains `2 ^ n - 1`, so the generic exponential bound is sharp and not
merely an artifact of the `defectRank` encoding.  Since `fanoutChain 3`
violates `AtMostOneDownstreamConsumer`, the quadratic bound of
`fixedProgram_acceptedSteps_quadratic` cannot be extended to general
well-formed federations by dropping that hypothesis.

Nonclaims.  Exponential worst-case repair length in asynchronous Boolean
networks is not new mathematics; Gadouleau and Richard (2018,
arXiv:1804.01931) establish fixing-length results of order `n * 2 ^ n` for
acyclic asynchronous graph families under a different metric and a different
step relation.  `fanoutChain` is not compiler output, carries no input, and
has no physical interpretation.  Its only role is to witness that the
single-consumer hypothesis is load-bearing for the quadratic bound.  Nothing
here asserts a lower bound for `fixedProgram`, a speedup claim, a hardware
contract, or any statement about scheduler attempts as opposed to genuine
canonical accepted steps.
-/

namespace OPH.RepairUniversality.FixedFederation

open OPH.RepairUniversality

noncomputable section

/-! ## Parity over the registers strictly below an index -/

/-- Exclusive-or of the register values strictly below `m`. -/
def parityBelow : Nat → State → Bool
  | 0, _ => false
  | m + 1, s => xor (s m) (parityBelow m s)

theorem parityBelow_congr (m : Nat) (s t : State)
    (h : ∀ i, i < m → s i = t i) :
    parityBelow m s = parityBelow m t := by
  induction m with
  | zero => rfl
  | succ m ih =>
      have hm : s m = t m := h m (Nat.lt_succ_self m)
      have hlow : parityBelow m s = parityBelow m t :=
        ih fun i hi => h i (by omega)
      simp only [parityBelow, hm, hlow]

/-- Parity below `m` is affine: flipping exactly one register strictly below
`m` flips the parity. -/
theorem parityBelow_flip (s t : State) (j : Nat)
    (hj : t j = !s j) (hother : ∀ i, i ≠ j → t i = s i) :
    ∀ m, j < m → parityBelow m t = !(parityBelow m s) := by
  intro m
  induction m with
  | zero =>
      intro h
      exact absurd h (by omega)
  | succ m ih =>
      intro hjm
      by_cases hcase : j = m
      · have hlow : parityBelow m t = parityBelow m s :=
          parityBelow_congr m t s fun i hi => hother i (by omega)
        have htop : t m = !s m := by
          rw [← hcase]; exact hj
        simp only [parityBelow, htop, hlow]
        cases hs : s m <;> cases hp : parityBelow m s <;> simp
      · have hlt : j < m := by omega
        have htop : t m = s m := hother m (by omega)
        have hrec : parityBelow m t = !(parityBelow m s) := ih hlt
        simp only [parityBelow, htop, hrec]
        cases hs : s m <;> cases hp : parityBelow m s <;> simp

/-! ## The fanout chain family -/

/-- Node `m` writes register `m` and reads every register strictly below it.
Its forced value is the negated parity of those registers, so the node is
affine in each dependency. -/
def fanoutChainNode (m : Nat) : Node where
  reg := m
  deps := List.range m
  val := fun s => !(parityBelow m s)
  val_local := by
    intro s t h
    have hp : parityBelow m s = parityBelow m t :=
      parityBelow_congr m s t fun i hi => h i (List.mem_range.mpr hi)
    rw [hp]

/-- The `n`-node fanout chain. -/
def fanoutChain (n : Nat) : List Node :=
  (List.range n).map fanoutChainNode

theorem fanoutChain_length (n : Nat) : (fanoutChain n).length = n := by
  simp [fanoutChain]

theorem fanoutChainNode_mem {m n : Nat} (h : m < n) :
    fanoutChainNode m ∈ fanoutChain n :=
  List.mem_map.mpr ⟨m, List.mem_range.mpr h, rfl⟩

theorem fanoutChain_mem_iff {x : Node} {n : Nat} :
    x ∈ fanoutChain n ↔ ∃ m, m < n ∧ x = fanoutChainNode m := by
  constructor
  · intro hx
    rcases List.mem_map.mp hx with ⟨m, hm, rfl⟩
    exact ⟨m, List.mem_range.mp hm, rfl⟩
  · rintro ⟨m, hm, rfl⟩
    exact fanoutChainNode_mem hm

theorem range_pairwise_lt (n : Nat) :
    (List.range n).Pairwise (fun a b => a < b) := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [List.range_succ]
      refine List.pairwise_append.mpr ⟨ih, by simp, ?_⟩
      intro a ha b hb
      simp only [List.mem_singleton] at hb
      subst hb
      exact List.mem_range.mp ha

theorem fanoutChain_wf (n : Nat) : NodesWF (fanoutChain n) := by
  refine ⟨?_, ?_⟩
  · show ((List.range n).map fanoutChainNode).Pairwise
      (fun a b => a.reg < b.reg)
    exact List.pairwise_map.mpr (range_pairwise_lt n)
  · intro x hx i hi
    rcases fanoutChain_mem_iff.mp hx with ⟨m, _hm, rfl⟩
    show i < m
    exact List.mem_range.mp hi

theorem fanoutChainNode_deps_lt (m : Nat) :
    ∀ i ∈ (fanoutChainNode m).deps, i < (fanoutChainNode m).reg := by
  intro i hi
  show i < m
  exact List.mem_range.mp hi

theorem fanoutChainNode_patch_le (m : Nat) :
    ∀ i ∈ (fanoutChainNode m).obs.patch, i ≤ m := by
  intro i hi
  have hi' : i ∈ m :: List.range m := hi
  rcases List.mem_cons.mp hi' with rfl | hmem
  · exact le_rfl
  · exact Nat.le_of_lt (List.mem_range.mp hmem)

/-- The verdict of node `m` reads only the registers up to `m`. -/
theorem fanoutChainNode_ok_congr (m : Nat) (s t : State)
    (h : ∀ r, r ≤ m → s r = t r) :
    (fanoutChainNode m).obs.ok s = (fanoutChainNode m).obs.ok t :=
  (fanoutChainNode m).obs.ok_local s t fun i hi =>
    h i (fanoutChainNode_patch_le m i hi)

theorem parityBelow_allFalse (m : Nat) : parityBelow m allFalse = false := by
  induction m with
  | zero => rfl
  | succ m ih => simp [parityBelow, allFalse, ih]

/-- Every node of the chain fails at the all-false state. -/
theorem fanoutChainNode_fails_allFalse (m : Nat) :
    (fanoutChainNode m).obs.ok allFalse = false := by
  show (allFalse m == !(parityBelow m allFalse)) = false
  rw [parityBelow_allFalse]
  simp [allFalse]

/-- Repairing a lower node breaks every accepting node above it: this is the
fanout that the compiler's single-consumer structure rules out. -/
theorem fanoutChainNode_fails_after_lower_repair
    {j m : Nat} (hjm : j < m) (s : State)
    (hj : (fanoutChainNode j).obs.ok s = false)
    (hm : (fanoutChainNode m).obs.ok s = true) :
    (fanoutChainNode m).obs.ok (repairNode (fanoutChainNode j) s) = false := by
  have hflip : repairNode (fanoutChainNode j) s j = !s j :=
    repairNode_output_flip_of_fails (fanoutChainNode j) s hj
  have hother : ∀ i, i ≠ j → repairNode (fanoutChainNode j) s i = s i := by
    intro i hi
    simp [repairNode, write, fanoutChainNode, hi]
  have hpar : parityBelow m (repairNode (fanoutChainNode j) s)
      = !(parityBelow m s) :=
    parityBelow_flip s (repairNode (fanoutChainNode j) s) j hflip hother m hjm
  have hreg : repairNode (fanoutChainNode j) s m = s m :=
    hother m (by omega)
  have hmval : s m = !(parityBelow m s) := by
    have : (s m == !(parityBelow m s)) = true := hm
    simpa [beq_iff_eq] using this
  show (repairNode (fanoutChainNode j) s m ==
    !(parityBelow m (repairNode (fanoutChainNode j) s))) = false
  rw [hreg, hpar, hmval]
  cases hp : parityBelow m s <;> simp

/-! ## The bottom-up doubling schedule -/

/-- `settleTail j c` drives nodes `j, ..., j + c - 1` from all-failing to all-
accepting.  Repairing node `j` flips register `j` and therefore breaks every
node above it again, so the inner block appears twice. -/
def settleTail (j : Nat) : Nat → List Node
  | 0 => []
  | c + 1 =>
      settleTail (j + 1) c ++ (fanoutChainNode j :: settleTail (j + 1) c)

theorem settleTail_length (j c : Nat) :
    (settleTail j c).length + 1 = 2 ^ c := by
  induction c generalizing j with
  | zero => rfl
  | succ c ih =>
      have h1 := ih (j + 1)
      show (settleTail (j + 1) c ++
        (fanoutChainNode j :: settleTail (j + 1) c)).length + 1 = 2 ^ (c + 1)
      rw [List.length_append, List.length_cons, pow_succ]
      omega

/-- The schedule invariant: from a state where nodes `j, ..., j + c - 1` all
fail, `settleTail j c` is a genuine accepted schedule that leaves all of them
accepting and never writes a register below `j`. -/
theorem settleTail_accepted (N : Nat) :
    ∀ (c j : Nat), j + c ≤ N → ∀ s : State,
      (∀ m, j ≤ m → m < j + c → (fanoutChainNode m).obs.ok s = false) →
      AcceptedNodeSchedule (fanoutChain N) s (settleTail j c) ∧
      (∀ m, j ≤ m → m < j + c →
        (fanoutChainNode m).obs.ok (runNodeSchedule s (settleTail j c)) = true) ∧
      (∀ r, r < j → runNodeSchedule s (settleTail j c) r = s r) := by
  intro c
  induction c with
  | zero =>
      intro j _hb s _hf
      refine ⟨?_, ?_, ?_⟩
      · show AcceptedNodeSchedule (fanoutChain N) s []
        trivial
      · intro m h1 h2
        omega
      · intro r _
        rfl
  | succ c ih =>
      intro j hb s hf
      have hbInner : j + 1 + c ≤ N := by omega
      have hfInner : ∀ m, j + 1 ≤ m → m < j + 1 + c →
          (fanoutChainNode m).obs.ok s = false :=
        fun m h1 h2 => hf m (by omega) (by omega)
      obtain ⟨hacc1, hok1, hfr1⟩ := ih (j + 1) hbInner s hfInner
      have hjfail :
          (fanoutChainNode j).obs.ok
            (runNodeSchedule s (settleTail (j + 1) c)) = false := by
        have hcongr := fanoutChainNode_ok_congr j
          (runNodeSchedule s (settleTail (j + 1) c)) s
          (fun r hr => hfr1 r (by omega))
        rw [hcongr]
        exact hf j le_rfl (by omega)
      have hfInner2 : ∀ m, j + 1 ≤ m → m < j + 1 + c →
          (fanoutChainNode m).obs.ok
            (repairNode (fanoutChainNode j)
              (runNodeSchedule s (settleTail (j + 1) c))) = false :=
        fun m h1 h2 =>
          fanoutChainNode_fails_after_lower_repair (by omega) _
            hjfail (hok1 m h1 h2)
      obtain ⟨hacc2, hok2, hfr2⟩ :=
        ih (j + 1) hbInner
          (repairNode (fanoutChainNode j)
            (runNodeSchedule s (settleTail (j + 1) c))) hfInner2
      have hrun :
          runNodeSchedule s (settleTail j (c + 1)) =
            runNodeSchedule
              (repairNode (fanoutChainNode j)
                (runNodeSchedule s (settleTail (j + 1) c)))
              (settleTail (j + 1) c) := by
        show runNodeSchedule s
            (settleTail (j + 1) c ++
              (fanoutChainNode j :: settleTail (j + 1) c)) = _
        rw [runNodeSchedule_append]
        rfl
      refine ⟨?_, ?_, ?_⟩
      · show AcceptedNodeSchedule (fanoutChain N) s
          (settleTail (j + 1) c ++ (fanoutChainNode j :: settleTail (j + 1) c))
        refine acceptedNodeSchedule_append hacc1 ?_
        exact ⟨fanoutChainNode_mem (by omega), hjfail, hacc2⟩
      · intro m h1 h2
        rw [hrun]
        by_cases hjm : m = j
        · rw [hjm]
          have hcongr := fanoutChainNode_ok_congr j
            (runNodeSchedule
              (repairNode (fanoutChainNode j)
                (runNodeSchedule s (settleTail (j + 1) c)))
              (settleTail (j + 1) c))
            (repairNode (fanoutChainNode j)
              (runNodeSchedule s (settleTail (j + 1) c)))
            (fun r hr => hfr2 r (by omega))
          rw [hcongr]
          exact repairNode_accepts (fanoutChainNode j) _
            (fanoutChainNode_deps_lt j)
        · exact hok2 m (by omega) (by omega)
      · intro r hr
        have hne : r ≠ j := Nat.ne_of_lt hr
        rw [hrun, hfr2 r (by omega)]
        have hskip :
            repairNode (fanoutChainNode j)
                (runNodeSchedule s (settleTail (j + 1) c)) r =
              runNodeSchedule s (settleTail (j + 1) c) r := by
          simp [repairNode, write, fanoutChainNode, hne]
        rw [hskip]
        exact hfr1 r (by omega)

theorem fanoutChain_settle_valid (n : Nat) :
    AcceptedNodeSchedule (fanoutChain n) allFalse (settleTail 0 n) ∧
    Consensus ((fanoutChain n).map Node.obs)
      (runNodeSchedule allFalse (settleTail 0 n)) := by
  obtain ⟨hacc, hok, _hfr⟩ :=
    settleTail_accepted n n 0 (by omega) allFalse
      fun m _h1 _h2 => fanoutChainNode_fails_allFalse m
  refine ⟨hacc, ?_⟩
  intro o ho
  rcases List.mem_map.mp ho with ⟨x, hx, rfl⟩
  rcases fanoutChain_mem_iff.mp hx with ⟨m, hm, rfl⟩
  exact hok m (Nat.zero_le m) (by omega)

/-! ## The two certificates -/

/-- The chain family realizes an accepted path of length exactly `2 ^ n - 1`
that ends in consensus. -/
theorem fanoutChain_exponential_lower (n : Nat) :
    ∃ t : State,
      CanonicalAcceptedSteps (fanoutChain n) (2 ^ n - 1) allFalse t ∧
      Consensus ((fanoutChain n).map Node.obs) t := by
  obtain ⟨hacc, hcons⟩ := fanoutChain_settle_valid n
  have hsteps := acceptedNodeSchedule_steps hacc
  have hlen : (settleTail 0 n).length = 2 ^ n - 1 := by
    have h := settleTail_length 0 n
    omega
  rw [hlen] at hsteps
  exact ⟨_, hsteps, hcons⟩

/-- The generic exponential budget: `defectRank` pays for every canonical
accepted step of any well-formed node list. -/
theorem canonicalAcceptedSteps_defectRank_budget
    (L : List Node) (hWF : NodesWF L) {m : Nat} {s t : State}
    (h : CanonicalAcceptedSteps L m s t) :
    m + defectRank L t ≤ defectRank L s := by
  induction h with
  | refl s => simp
  | tail hpath hstep ih =>
      have hlt := canonicalAcceptedStep_defectRank_lt L hWF hstep
      omega

/-- Every well-formed node list admits fewer than `2 ^ L.length` canonical
accepted steps.  No single-consumer hypothesis is used. -/
theorem canonicalAcceptedSteps_lt_pow
    (L : List Node) (hWF : NodesWF L) {m : Nat} {s t : State}
    (h : CanonicalAcceptedSteps L m s t) :
    m < 2 ^ L.length := by
  have hbudget := canonicalAcceptedSteps_defectRank_budget L hWF h
  have hcap := defectRank_lt_pow L s
  omega

/-- The chain family is outside the hypothesis of the quadratic bound:
register `0` is read by two later nodes already at size three. -/
theorem fanoutChain_not_single_consumer :
    ¬ AtMostOneDownstreamConsumer (fanoutChain 3) := by
  intro h
  have hcount := h (fanoutChainNode 0) (fanoutChainNode_mem (by omega))
  have htwo :
      (dependencyOccurrences (fanoutChain 3)).count (fanoutChainNode 0).reg
        = 2 := by decide
  omega

/-- Positive control.  At two nodes the same construction still satisfies
`AtMostOneDownstreamConsumer`, and its accepted path length `2 ^ 2 - 1 = 3`
coincides with the triangular budget `triangle 2 = 3`.  The negative control
below is therefore not an artifact of the construction: the separation begins
exactly at the first size where the hypothesis fails. -/
theorem fanoutChain_two_single_consumer :
    AtMostOneDownstreamConsumer (fanoutChain 2) ∧
    triangle (fanoutChain 2).length = 3 := by
  refine ⟨?_, ?_⟩
  · show ∀ n ∈ fanoutChain 2,
      (dependencyOccurrences (fanoutChain 2)).count n.reg ≤ 1
    decide
  · rw [fanoutChain_length]
    decide

/-- At size three the fanout family already realizes more canonical accepted
steps than the triangular budget that `fixedProgram_acceptedSteps_triangle`
grants to single-consumer programs of the same node count.  This is the exact
sense in which `AtMostOneDownstreamConsumer` is load-bearing: dropping it
makes the triangular bound false, not merely unproved. -/
theorem fanoutChain_exceeds_triangle :
    ∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 3) 7 allFalse t ∧
      Consensus ((fanoutChain 3).map Node.obs) t ∧
      triangle (fanoutChain 3).length = 6 ∧
      ¬ (7 ≤ triangle (fanoutChain 3).length) := by
  obtain ⟨t, hsteps, hcons⟩ := fanoutChain_exponential_lower 3
  refine ⟨t, by simpa using hsteps, hcons, ?_, ?_⟩
  · rw [fanoutChain_length]
    decide
  · rw [fanoutChain_length]
    decide

/-- The generic exponential bound is sharp: the chain family attains the
largest step count the bound allows. -/
theorem fanoutChain_sharp_exponential (n : Nat) :
    (∀ (L : List Node), NodesWF L → ∀ {m : Nat} {s t : State},
      CanonicalAcceptedSteps L m s t → m < 2 ^ L.length) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain n) (2 ^ n - 1) allFalse t ∧
      Consensus ((fanoutChain n).map Node.obs) t ∧
      NodesWF (fanoutChain n) ∧
      (fanoutChain n).length = n) := by
  refine ⟨fun L hWF _m _s _t h => canonicalAcceptedSteps_lt_pow L hWF h, ?_⟩
  obtain ⟨t, hsteps, hcons⟩ := fanoutChain_exponential_lower n
  exact ⟨t, hsteps, hcons, fanoutChain_wf n, fanoutChain_length n⟩

/-! ## Kernel-checked finite receipts -/

theorem fanoutChain_defectRank_allFalse_two :
    defectRank (fanoutChain 2) allFalse = 3 := by decide

theorem fanoutChain_defectRank_allFalse_three :
    defectRank (fanoutChain 3) allFalse = 7 := by decide

theorem fanoutChain_defectRank_allFalse_four :
    defectRank (fanoutChain 4) allFalse = 15 := by decide

theorem settleTail_length_two : (settleTail 0 2).length = 3 := by decide

theorem settleTail_length_three : (settleTail 0 3).length = 7 := by decide

theorem settleTail_length_four : (settleTail 0 4).length = 15 := by decide

example :
    ∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 2) 3 allFalse t ∧
      Consensus ((fanoutChain 2).map Node.obs) t := by
  simpa using fanoutChain_exponential_lower 2

example :
    ∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 3) 7 allFalse t ∧
      Consensus ((fanoutChain 3).map Node.obs) t := by
  simpa using fanoutChain_exponential_lower 3

example :
    ∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 4) 15 allFalse t ∧
      Consensus ((fanoutChain 4).map Node.obs) t := by
  simpa using fanoutChain_exponential_lower 4

#print axioms fanoutChain_wf
#print axioms fanoutChain_exponential_lower
#print axioms canonicalAcceptedSteps_lt_pow
#print axioms fanoutChain_not_single_consumer
#print axioms fanoutChain_two_single_consumer
#print axioms fanoutChain_exceeds_triangle
#print axioms fanoutChain_sharp_exponential
#print axioms settleTail_accepted
#print axioms settleTail_length
#print axioms parityBelow_flip
#print axioms fanoutChainNode_fails_after_lower_repair

end

end OPH.RepairUniversality.FixedFederation
