import Computation.FixedFederationComplexity

/-!
# Fixed-federation execution controls

The fixtures in this module use actual `go` compiler output.  They separate
tail-relative fairness from recurrence and expose the accepted-step lower
mechanism without assigning a physical duration to an attempt.
-/

namespace OPH.RepairUniversality.FixedFederation

open OPH.RepairUniversality

noncomputable section

/-! ## Recurrence is strictly stronger than tail-relative fairness -/

def recurrenceGapFormula : Formula 0 :=
  .nand (.const false) (.const false)

def recurrenceGapNode0 : Node where
  reg := 0
  deps := []
  val := fun _ => false
  val_local := by intro _ _ _; rfl

def recurrenceGapNode1 : Node where
  reg := 1
  deps := []
  val := fun _ => false
  val_local := by intro _ _ _; rfl

theorem recurrenceGapNode0_mem :
    recurrenceGapNode0 ∈ fixedProgram recurrenceGapFormula := by
  simp [recurrenceGapFormula, recurrenceGapNode0, fixedProgram, go]

theorem recurrenceGapNode1_mem :
    recurrenceGapNode1 ∈ fixedProgram recurrenceGapFormula := by
  simp [recurrenceGapFormula, recurrenceGapNode1, fixedProgram, go]

def recurrenceGapState : State :=
  sweepFrom allFalse (fixedProgram recurrenceGapFormula)

def recurrenceGapScheduler :
    NodeScheduler (fixedProgram recurrenceGapFormula) :=
  fun _ _ => ⟨recurrenceGapNode0, recurrenceGapNode0_mem⟩

theorem recurrenceGap_consensus :
    Consensus (fixedFederation recurrenceGapFormula) recurrenceGapState :=
  fixedSweep_consensus recurrenceGapFormula allFalse

theorem recurrenceGap_run_constant (n : Nat) :
    attemptRun (fixedProgram recurrenceGapFormula) n
      recurrenceGapScheduler recurrenceGapState = recurrenceGapState := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [attemptRun_last_step, ih]
      change repairNode recurrenceGapNode0 recurrenceGapState = recurrenceGapState
      apply repairNode_eq_self_of_accepts
      exact recurrenceGap_consensus recurrenceGapNode0.obs
        (List.mem_map.mpr ⟨recurrenceGapNode0,
          recurrenceGapNode0_mem, rfl⟩)

theorem recurrenceGap_pathwiseWeakFair :
    NodePathwiseWeakFair (fixedProgram recurrenceGapFormula)
      recurrenceGapScheduler recurrenceGapState := by
  apply stableConsensusTail_nodePathwiseWeakFair
  exact ⟨0, recurrenceGap_consensus, fun n _ => recurrenceGap_run_constant n⟩

theorem recurrenceGap_not_siteRecurrent :
    ¬ NodeSiteRecurrent (fixedProgram recurrenceGapFormula)
      recurrenceGapScheduler recurrenceGapState := by
  intro hsite
  let member : {n : Node // n ∈ fixedProgram recurrenceGapFormula} :=
    ⟨recurrenceGapNode1, recurrenceGapNode1_mem⟩
  obtain ⟨m, _hm, hreg⟩ := hsite member 0
  simp [recurrenceGapScheduler, member,
    recurrenceGapNode0, recurrenceGapNode1] at hreg

theorem recurrenceGap_not_memberRecurrent :
    ¬ NodeMemberRecurrent (fixedProgram recurrenceGapFormula)
      recurrenceGapScheduler recurrenceGapState := by
  intro hmember
  exact recurrenceGap_not_siteRecurrent
    (nodeMemberRecurrent_nodeSiteRecurrent hmember)

theorem recurrence_strictly_stronger_than_tail_fairness :
    ∃ (k : Nat) (phi : Formula k) (sigma : NodeScheduler (fixedProgram phi))
      (s : State),
      NodePathwiseWeakFair (fixedProgram phi) sigma s ∧
      ¬ NodeMemberRecurrent (fixedProgram phi) sigma s ∧
      ¬ NodeSiteRecurrent (fixedProgram phi) sigma s := by
  exact ⟨0, recurrenceGapFormula, recurrenceGapScheduler, recurrenceGapState,
    recurrenceGap_pathwiseWeakFair,
    recurrenceGap_not_memberRecurrent,
    recurrenceGap_not_siteRecurrent⟩

/-! ## Weak fairness has no finite attempt horizon -/

def delayedRoundRobin (delay : Nat) :
    NodeScheduler (fixedProgram recurrenceGapFormula) :=
  fun n s =>
    if n < delay then
      ⟨recurrenceGapNode0, recurrenceGapNode0_mem⟩
    else
      fixedRoundRobinScheduler recurrenceGapFormula n s

theorem delayedRoundRobin_memberRecurrent (delay : Nat) (s : State) :
    NodeMemberRecurrent (fixedProgram recurrenceGapFormula)
      (delayedRoundRobin delay) s := by
  intro member N
  obtain ⟨m, hm, hselected⟩ :=
    fixedRoundRobin_memberRecurrent recurrenceGapFormula s member
      (max N delay)
  have hNm : N ≤ m := le_trans (le_max_left N delay) hm
  have hdelay : delay ≤ m := le_trans (le_max_right N delay) hm
  refine ⟨m, hNm, ?_⟩
  simp only [delayedRoundRobin, if_neg (not_lt_of_ge hdelay)]
  simpa [fixedRoundRobinScheduler, roundRobinScheduler] using hselected

theorem delayedRoundRobin_pathwiseWeakFair (delay : Nat) (s : State) :
    NodePathwiseWeakFair (fixedProgram recurrenceGapFormula)
      (delayedRoundRobin delay) s :=
  nodeMemberRecurrent_nodePathwiseWeakFair
    (delayedRoundRobin_memberRecurrent delay s)

theorem recurrenceGapNode0_accepts_allFalse :
    recurrenceGapNode0.obs.ok allFalse = true := by
  rfl

theorem delayedRoundRobin_run_before (delay n : Nat) (hn : n ≤ delay) :
    attemptRun (fixedProgram recurrenceGapFormula) n
      (delayedRoundRobin delay) allFalse = allFalse := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [attemptRun_last_step, ih (Nat.le_trans (Nat.le_succ n) hn)]
      have hndelay : n < delay := by omega
      simp only [delayedRoundRobin, if_pos hndelay]
      change repairNode recurrenceGapNode0 allFalse = allFalse
      apply repairNode_eq_self_of_accepts
      exact recurrenceGapNode0_accepts_allFalse

def recurrenceGapOutNode : Node where
  reg := 2
  deps := [0, 1]
  val := fun s => !(s 0 && s 1)
  val_local := by
    intro s t h
    rw [h 0 (by simp), h 1 (by simp)]

theorem recurrenceGap_allFalse_not_consensus :
    ¬ Consensus (fixedFederation recurrenceGapFormula) allFalse := by
  intro hcons
  have hout : recurrenceGapOutNode ∈ fixedProgram recurrenceGapFormula := by
    simp [recurrenceGapOutNode, recurrenceGapFormula, fixedProgram, go]
  have h := hcons recurrenceGapOutNode.obs
    (List.mem_map.mpr ⟨recurrenceGapOutNode, hout, rfl⟩)
  simp [recurrenceGapOutNode, Node.obs, allFalse] at h

/-- Every finite attempt budget can be exceeded by a finite stuttering prefix
of an actual fixed-program scheduler that remains pathwise weak-fair. -/
theorem weakFair_no_uniform_attempt_bound (budget : Nat) :
    ∃ sigma : NodeScheduler (fixedProgram recurrenceGapFormula),
      NodePathwiseWeakFair (fixedProgram recurrenceGapFormula) sigma allFalse ∧
      ∀ n, n ≤ budget →
        ¬ Consensus (fixedFederation recurrenceGapFormula)
          (attemptRun (fixedProgram recurrenceGapFormula) n sigma allFalse) := by
  let sigma := delayedRoundRobin (budget + 1)
  refine ⟨sigma, delayedRoundRobin_pathwiseWeakFair (budget + 1) allFalse, ?_⟩
  intro n hn
  rw [show attemptRun (fixedProgram recurrenceGapFormula) n sigma allFalse =
      allFalse by
    exact delayedRoundRobin_run_before (budget + 1) n (by omega)]
  exact recurrenceGap_allFalse_not_consensus

/-! ## Typed premise controls -/

/-- Starting the allocator below the input range can collide with an input
register and falsifies the exact compiler use profile. -/
def allocatorCollisionFormula : Formula 1 :=
  .var ⟨0, by omega⟩

theorem allocatorSeparation_needed :
    ¬ GoUseProfile allocatorCollisionFormula 0 := by
  intro h
  have hout := h.output_unused
  simp [allocatorCollisionFormula, go, dependencyOccurrences] at hout

/-- Without a bounded-waste premise, even a weak-fair actual scheduler need
not reach consensus inside the triangular attempt horizon. -/
theorem boundedWaste_premise_needed :
    ∃ sigma : NodeScheduler (fixedProgram recurrenceGapFormula),
      NodePathwiseWeakFair (fixedProgram recurrenceGapFormula) sigma allFalse ∧
      ¬ ∃ N,
        N ≤ (0 + 1) * triangle (fixedProgram recurrenceGapFormula).length ∧
        Consensus (fixedFederation recurrenceGapFormula)
          (attemptRun (fixedProgram recurrenceGapFormula) N sigma allFalse) ∧
        ∀ n, N ≤ n →
          attemptRun (fixedProgram recurrenceGapFormula) n sigma allFalse =
            attemptRun (fixedProgram recurrenceGapFormula) N sigma allFalse := by
  let budget := triangle (fixedProgram recurrenceGapFormula).length
  obtain ⟨sigma, hfair, hdelay⟩ :=
    weakFair_no_uniform_attempt_bound budget
  refine ⟨sigma, hfair, ?_⟩
  intro hbound
  obtain ⟨N, hN, hcons, _⟩ := hbound
  apply hdelay N
  · simpa [budget] using hN
  · exact hcons

#print axioms recurrence_strictly_stronger_than_tail_fairness
#print axioms weakFair_no_uniform_attempt_bound
#print axioms allocatorSeparation_needed
#print axioms boundedWaste_premise_needed

/-! ## Actual nested-negation compiler family -/

def combConstNode (r : Nat) (b : Bool) : Node where
  reg := r
  deps := []
  val := fun _ => b
  val_local := by intro _ _ _; rfl

def combNandNode (r left right : Nat) : Node where
  reg := r
  deps := [left, right]
  val := fun s => !(s left && s right)
  val_local := by
    intro s t h
    rw [h left (by simp), h right (by simp)]

def negComb : Nat → Formula 0
  | 0 => .const false
  | d + 1 => .nand (negComb d) (.const true)

def combNodes : Nat → List Node
  | 0 => [combConstNode 0 false]
  | d + 1 =>
      combNodes d ++
        [combConstNode (2 * d + 1) true,
          combNandNode (2 * d + 2) (2 * d) (2 * d + 1)]

theorem go_negComb (d : Nat) :
    go (negComb d) 0 =
      ⟨combNodes d, 2 * d, 2 * d + 1⟩ := by
  induction d with
  | zero => rfl
  | succ d ih =>
      simp only [negComb, go, ih, combNodes]
      congr 1

theorem fixedProgram_negComb (d : Nat) :
    fixedProgram (negComb d) = combNodes d := by
  simpa [fixedProgram] using congrArg Out.nodes (go_negComb d)

theorem fixedOutReg_negComb (d : Nat) :
    fixedOutReg (negComb d) = 2 * d := by
  simpa [fixedOutReg] using congrArg Out.out (go_negComb d)

theorem combNodes_length (d : Nat) :
    (combNodes d).length = 2 * d + 1 := by
  induction d with
  | zero => rfl
  | succ d ih =>
      simp [combNodes, ih]
      omega

theorem fixedProgram_negComb_length (d : Nat) :
    (fixedProgram (negComb d)).length = 2 * d + 1 := by
  rw [fixedProgram_negComb, combNodes_length]

def combBlocks : Nat → List (List Nat)
  | 0 => [[0]]
  | d + 1 =>
      [[2 * d + 1, 2 * d + 2]] ++
        (combBlocks d).map fun block => block ++ [2 * d + 2]

def combSchedule (d : Nat) : List Nat :=
  (combBlocks d).flatten

def lowerSteps : Nat → Nat
  | 0 => 1
  | d + 1 => lowerSteps d + d + 3

theorem lowerSteps_double (d : Nat) :
    2 * lowerSteps d = d * d + 5 * d + 2 := by
  induction d with
  | zero => rfl
  | succ d ih =>
      calc
        2 * lowerSteps (d + 1) =
            2 * (lowerSteps d + d + 3) := by rfl
        _ = 2 * lowerSteps d + 2 * d + 6 := by omega
        _ = (d * d + 5 * d + 2) + 2 * d + 6 := by rw [ih]
        _ = (d + 1) * (d + 1) + 5 * (d + 1) + 2 := by ring

theorem combBlocks_count (d : Nat) :
    (combBlocks d).length = d + 1 := by
  induction d with
  | zero => rfl
  | succ d ih =>
      simp [combBlocks, ih]

theorem flatten_append_singleton_length
    {α : Type} (blocks : List (List α)) (x : α) :
    (blocks.map (fun block => block ++ [x])).flatten.length =
      blocks.flatten.length + blocks.length := by
  induction blocks with
  | nil => simp
  | cons block blocks ih =>
      simp [ih]
      omega

theorem combSchedule_length (d : Nat) :
    (combSchedule d).length = lowerSteps d := by
  induction d with
  | zero => rfl
  | succ d ih =>
      change (combBlocks d).flatten.length = lowerSteps d at ih
      simp only [combSchedule, combBlocks, List.flatten_append,
        List.flatten_cons, List.flatten_nil, List.length_append,
        List.length_cons, List.length_nil]
      rw [flatten_append_singleton_length, ih, combBlocks_count]
      simp [lowerSteps]
      omega

theorem negComb_size_lowerSteps_identity (d : Nat) :
    let n := (fixedProgram (negComb d)).length
    let m := lowerSteps d
    8 * m + 1 = n * n + 8 * n := by
  simp only [fixedProgram_negComb_length]
  have hs := lowerSteps_double d
  calc
    8 * lowerSteps d + 1 = 4 * (2 * lowerSteps d) + 1 := by ring
    _ = 4 * (d * d + 5 * d + 2) + 1 := by rw [hs]
    _ = (2 * d + 1) * (2 * d + 1) + 8 * (2 * d + 1) := by ring

/-! ## Accepted trace for the lower family -/

def combInitial : State
  | 0 => true
  | n + 1 => !combInitial n

def combLeaf : Node := combConstNode 0 false

def combTrue (q : Nat) : Node :=
  combConstNode (2 * q + 1) true

def combGate (q : Nat) : Node :=
  combNandNode (2 * q + 2) (2 * q) (2 * q + 1)

def gateTail (q : Nat) : Nat → List Node
  | 0 => []
  | n + 1 => combGate q :: gateTail (q + 1) n

def combBlockNodes (d : Nat) : Nat → List Node
  | 0 => combLeaf :: gateTail 0 d
  | q + 1 => combTrue q :: gateTail q (d - q)

def combNodeSchedule (d : Nat) : List Node :=
  ((List.range (d + 1)).reverse.map (combBlockNodes d)).flatten

def combScheduleFrom (d : Nat) : Nat → List Node
  | 0 => combBlockNodes d 0
  | i + 1 => combBlockNodes d (i + 1) ++ combScheduleFrom d i

def combAcceptedSchedule (d : Nat) : List Node :=
  combScheduleFrom d d

def runNodeSchedule : State → List Node → State
  | s, [] => s
  | s, n :: ns => runNodeSchedule (repairNode n s) ns

def AcceptedNodeSchedule (L : List Node) : State → List Node → Prop
  | _, [] => True
  | s, n :: ns =>
      n ∈ L ∧ n.obs.ok s = false ∧
        AcceptedNodeSchedule L (repairNode n s) ns

theorem acceptedNodeSchedule_append {L : List Node}
    {s : State} {A B : List Node}
    (hA : AcceptedNodeSchedule L s A)
    (hB : AcceptedNodeSchedule L (runNodeSchedule s A) B) :
    AcceptedNodeSchedule L s (A ++ B) := by
  induction A generalizing s with
  | nil => simpa [runNodeSchedule, AcceptedNodeSchedule] using hB
  | cons n A ih =>
      rcases hA with ⟨hn, hfail, hA⟩
      exact ⟨hn, hfail, ih hA hB⟩

theorem runNodeSchedule_append (s : State) (A B : List Node) :
    runNodeSchedule s (A ++ B) =
      runNodeSchedule (runNodeSchedule s A) B := by
  induction A generalizing s with
  | nil => rfl
  | cons n A ih =>
      simp only [List.cons_append, runNodeSchedule]
      exact ih _

theorem canonicalAcceptedSteps_trans {L : List Node}
    {m n : Nat} {s t u : State}
    (h₁ : CanonicalAcceptedSteps L m s t)
    (h₂ : CanonicalAcceptedSteps L n t u) :
    CanonicalAcceptedSteps L (m + n) s u := by
  induction h₂ with
  | refl => simpa using h₁
  | tail hpath hstep ih =>
      simpa [Nat.add_assoc] using
        CanonicalAcceptedSteps.tail (ih h₁) hstep

theorem acceptedNodeSchedule_steps {L : List Node}
    {s : State} {ns : List Node}
    (h : AcceptedNodeSchedule L s ns) :
    CanonicalAcceptedSteps L ns.length s (runNodeSchedule s ns) := by
  induction ns generalizing s with
  | nil => exact CanonicalAcceptedSteps.refl s
  | cons n ns ih =>
      rcases h with ⟨hn, hfail, htail⟩
      have hone : CanonicalAcceptedSteps L 1 s (repairNode n s) := by
        simpa using CanonicalAcceptedSteps.tail
          (CanonicalAcceptedSteps.refl s) ⟨n, hn, hfail, rfl⟩
      have hrest := ih htail
      simpa [runNodeSchedule, Nat.add_comm] using
        canonicalAcceptedSteps_trans hone hrest

def CombReady (d i : Nat) (s : State) : Prop :=
  i ≤ d ∧
  s 0 = true ∧
  (∀ q, q < i →
    s (2 * q + 1) = false ∧ s (2 * q + 2) = true) ∧
  (∀ q, i ≤ q → q < d →
    s (2 * q + 1) = true ∧ (combGate q).obs.ok s = true)

def CombDone (d : Nat) (s : State) : Prop :=
  s 0 = false ∧
  ∀ q, q < d →
    s (2 * q + 1) = true ∧ (combGate q).obs.ok s = true

theorem combInitial_pair (q : Nat) :
    combInitial (2 * q) = true ∧
      combInitial (2 * q + 1) = false := by
  induction q with
  | zero => decide
  | succ q ih =>
      constructor <;>
        simp only [combInitial] at ih ⊢ <;>
        simp_all

theorem combInitial_ready (d : Nat) :
    CombReady d d combInitial := by
  refine ⟨le_rfl, (combInitial_pair 0).1, ?_, ?_⟩
  · intro q hq
    exact ⟨(combInitial_pair q).2,
      by simpa [show 2 * q + 2 = 2 * (q + 1) by omega]
        using (combInitial_pair (q + 1)).1⟩
  · intro q hd hq
    omega

theorem combLeaf_mem (d : Nat) : combLeaf ∈ combNodes d := by
  induction d with
  | zero => simp [combLeaf, combNodes]
  | succ d ih => simp [combNodes, ih]

theorem combTrue_mem {q d : Nat} (h : q < d) :
    combTrue q ∈ combNodes d := by
  induction d with
  | zero => omega
  | succ d ih =>
      by_cases hqd : q = d
      · subst q
        simp [combNodes, combTrue]
      · have hlt : q < d := by omega
        simp [combNodes, ih hlt]

theorem combGate_mem {q d : Nat} (h : q < d) :
    combGate q ∈ combNodes d := by
  induction d with
  | zero => omega
  | succ d ih =>
      by_cases hqd : q = d
      · subst q
        simp [combNodes, combGate]
      · have hlt : q < d := by omega
        simp [combNodes, ih hlt]

theorem combGate_fails_after_left_flip
    (q : Nat) (base current : State)
    (hrBase : base (2 * q + 1) = true)
    (hrCurrent : current (2 * q + 1) = true)
    (haccept : (combGate q).obs.ok base = true)
    (hleft : current (2 * q) = !base (2 * q))
    (hout : current (2 * q + 2) = base (2 * q + 2)) :
    (combGate q).obs.ok current = false ∧
    repairNode (combGate q) current (2 * q + 2) =
      !base (2 * q + 2) := by
  cases h : base (2 * q) <;>
    simp_all [combGate, combNandNode, Node.obs, repairNode, write]

theorem combGate_deps_lt (q : Nat) :
    ∀ i ∈ (combGate q).deps, i < (combGate q).reg := by
  intro i hi
  simp [combGate, combNandNode] at hi ⊢
  rcases hi with rfl | rfl <;> omega

theorem gateTail_accepted
    {d q n : Nat} (hbound : q + n ≤ d)
    (base current : State)
    (hleft : current (2 * q) = !base (2 * q))
    (habove : ∀ r, 2 * q < r → current r = base r)
    (hbase : ∀ j, j < n →
      base (2 * (q + j) + 1) = true ∧
      (combGate (q + j)).obs.ok base = true) :
    AcceptedNodeSchedule (combNodes d) current (gateTail q n) ∧
    let final := runNodeSchedule current (gateTail q n)
    (∀ j, j < n →
      final (2 * (q + j) + 1) = true ∧
      (combGate (q + j)).obs.ok final = true) ∧
    (∀ r, r < 2 * q + 2 → final r = current r) := by
  induction n generalizing q current with
  | zero =>
      simp [gateTail, runNodeSchedule, AcceptedNodeSchedule]
  | succ n ih =>
      have hqd : q < d := by omega
      have hsideBase := (hbase 0 (Nat.zero_lt_succ n)).1
      have hacceptBase := (hbase 0 (Nat.zero_lt_succ n)).2
      have hsideCurrent : current (2 * q + 1) = true := by
        rw [habove (2 * q + 1) (by omega)]
        simpa using hsideBase
      have houtCurrent :
          current (2 * q + 2) = base (2 * q + 2) :=
        habove (2 * q + 2) (by omega)
      have hflip := combGate_fails_after_left_flip q base current
        hsideBase hsideCurrent hacceptBase hleft houtCurrent
      let next := repairNode (combGate q) current
      have hnextLeft : next (2 * (q + 1)) = !base (2 * (q + 1)) := by
        simpa [next, show 2 * (q + 1) = 2 * q + 2 by omega] using hflip.2
      have hnextAbove : ∀ r, 2 * (q + 1) < r → next r = base r := by
        intro r hr
        rw [show next r = current r by
          simp [next, repairNode, write, combGate, combNandNode]
          omega]
        exact habove r (by omega)
      have hbaseTail : ∀ j, j < n →
          base (2 * ((q + 1) + j) + 1) = true ∧
          (combGate ((q + 1) + j)).obs.ok base = true := by
        intro j hj
        simpa [Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
          hbase (j + 1) (by omega)
      obtain ⟨htailAccepted, htailProcessed, htailFrame⟩ :=
        ih (q := q + 1) (current := next) (by omega)
          hnextLeft hnextAbove hbaseTail
      refine ⟨⟨combGate_mem hqd, hflip.1, htailAccepted⟩, ?_, ?_⟩
      · intro j hj
        rcases j with _ | j
        · have hsideNext : next (2 * q + 1) = true := by
            simp [next, repairNode, write, combGate, combNandNode,
              hsideCurrent]
          have hsideFinal := htailFrame (2 * q + 1) (by omega)
          have hacceptNext : (combGate q).obs.ok next = true := by
            exact repairNode_accepts (combGate q) current (combGate_deps_lt q)
          have hobsEq :
              (combGate q).obs.ok
                  (runNodeSchedule next (gateTail (q + 1) n)) =
                (combGate q).obs.ok next := by
            apply (combGate q).obs.ok_local
            intro i hi
            apply htailFrame i
            simp [combGate, combNandNode, Node.obs] at hi
            rcases hi with rfl | rfl | rfl <;> omega
          exact ⟨hsideFinal.trans hsideNext, hobsEq.trans hacceptNext⟩
        · simpa [Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
            htailProcessed j (by omega)
      · intro r hr
        calc
          runNodeSchedule current (gateTail q (n + 1)) r =
              runNodeSchedule next (gateTail (q + 1) n) r := rfl
          _ = next r := htailFrame r (by omega)
          _ = current r := by
            simp [next, repairNode, write, combGate, combNandNode]
            omega

theorem repairNode_output_flip_of_fails
    (n : Node) (s : State) (hfail : n.obs.ok s = false) :
    repairNode n s n.reg = !s n.reg := by
  cases hreg : s n.reg <;> cases hval : n.val s <;>
    simp_all [Node.obs, repairNode, write]

theorem combGate_fails_after_right_flip
    (q : Nat) (base current : State)
    (hleftBase : base (2 * q) = true)
    (hleft : current (2 * q) = base (2 * q))
    (hright : current (2 * q + 1) = !base (2 * q + 1))
    (hout : current (2 * q + 2) = base (2 * q + 2))
    (haccept : (combGate q).obs.ok base = true) :
    (combGate q).obs.ok current = false := by
  cases hl : base (2 * q) <;> cases hr : base (2 * q + 1) <;>
    simp_all [combGate, combNandNode, Node.obs]

theorem combBlock_transition
    {d i : Nat} {s : State}
    (hready : CombReady d i s) :
    let block := combBlockNodes d i
    let t := runNodeSchedule s block
    AcceptedNodeSchedule (combNodes d) s block ∧
      match i with
      | 0 => CombDone d t
      | q + 1 => CombReady d q t := by
  rcases hready with ⟨hid, hzero, hunprocessed, hprocessed⟩
  cases i with
  | zero =>
      let afterLeaf := repairNode combLeaf s
      have hleafFail : combLeaf.obs.ok s = false := by
        simp [combLeaf, combConstNode, Node.obs, hzero]
      have hleft : afterLeaf 0 = !s 0 := by
        simpa [afterLeaf, combLeaf] using
          repairNode_output_flip_of_fails combLeaf s hleafFail
      have habove : ∀ r, 0 < r → afterLeaf r = s r := by
        intro r hr
        simp [afterLeaf, combLeaf, combConstNode, repairNode, write]
        omega
      have hbase : ∀ j, j < d →
          s (2 * (0 + j) + 1) = true ∧
          (combGate (0 + j)).obs.ok s = true := by
        intro j hj
        simpa using hprocessed j (Nat.zero_le j) hj
      obtain ⟨htail, htailProcessed, htailFrame⟩ :=
        gateTail_accepted (d := d) (q := 0) (n := d) (by omega)
          s afterLeaf hleft habove hbase
      refine ⟨⟨combLeaf_mem d, hleafFail, htail⟩, ?_⟩
      refine ⟨?_, ?_⟩
      · have hframeZero := htailFrame 0 (by omega)
        calc
          runNodeSchedule s (combBlockNodes d 0) 0 =
              runNodeSchedule afterLeaf (gateTail 0 d) 0 := rfl
          _ = afterLeaf 0 := hframeZero
          _ = false := by simpa [hzero] using hleft
      · intro q hq
        simpa using htailProcessed q hq
  | succ q =>
      have hqd : q < d := by omega
      have hqValues := hunprocessed q (Nat.lt_succ_self q)
      let afterSide := repairNode (combTrue q) s
      have hsideFail : (combTrue q).obs.ok s = false := by
        simp [combTrue, combConstNode, Node.obs, hqValues.1]
      have hsideTrue : afterSide (2 * q + 1) = true := by
        simp [afterSide, combTrue, combConstNode, repairNode, write]
      have hleftSame : afterSide (2 * q) = s (2 * q) := by
        simp [afterSide, combTrue, combConstNode, repairNode, write]
      have houtSame : afterSide (2 * q + 2) = s (2 * q + 2) := by
        simp [afterSide, combTrue, combConstNode, repairNode, write]
      have hleftBase : s (2 * q) = true := by
        cases q with
        | zero => simpa using hzero
        | succ p =>
            simpa [show 2 * (p + 1) = 2 * p + 2 by omega] using
              (hunprocessed p (by omega)).2
      have hgateAcceptS : (combGate q).obs.ok s = true := by
        simp [combGate, combNandNode, Node.obs, hqValues.1, hqValues.2]
      have hrightFlip : afterSide (2 * q + 1) = !s (2 * q + 1) := by
        simp [hsideTrue, hqValues.1]
      have hgateFail : (combGate q).obs.ok afterSide = false :=
        combGate_fails_after_right_flip q s afterSide
          hleftBase hleftSame hrightFlip houtSame hgateAcceptS
      let afterGate := repairNode (combGate q) afterSide
      have hgateLeftFlip :
          afterGate (2 * (q + 1)) = !afterSide (2 * (q + 1)) := by
        simpa [afterGate, show 2 * (q + 1) = 2 * q + 2 by omega] using
          repairNode_output_flip_of_fails (combGate q) afterSide hgateFail
      have hgateAbove : ∀ r, 2 * (q + 1) < r →
          afterGate r = afterSide r := by
        intro r hr
        simp [afterGate, combGate, combNandNode, repairNode, write]
        omega
      have houterBase : ∀ j, j < d - (q + 1) →
          afterSide (2 * ((q + 1) + j) + 1) = true ∧
          (combGate ((q + 1) + j)).obs.ok afterSide = true := by
        intro j hj
        have hidx : q + 1 + j < d := by omega
        have houter := hprocessed (q + 1 + j) (by omega) hidx
        have hsideEq :
            afterSide (2 * (q + 1 + j) + 1) =
              s (2 * (q + 1 + j) + 1) := by
          simp [afterSide, combTrue, combConstNode, repairNode, write]
          omega
        have hobsEq :
            (combGate (q + 1 + j)).obs.ok afterSide =
              (combGate (q + 1 + j)).obs.ok s := by
          apply (combGate (q + 1 + j)).obs.ok_local
          intro r hr
          simp [combGate, combNandNode, Node.obs] at hr
          rcases hr with rfl | rfl | rfl <;>
            simp [afterSide, combTrue, combConstNode, repairNode, write] <;>
            omega
        exact ⟨hsideEq.trans houter.1, hobsEq.trans houter.2⟩
      obtain ⟨htail, htailProcessed, htailFrame⟩ :=
        gateTail_accepted (d := d) (q := q + 1)
          (n := d - (q + 1)) (by omega)
          afterSide afterGate hgateLeftFlip hgateAbove houterBase
      have hlen : d - q = (d - (q + 1)) + 1 := by omega
      have hrest :
          AcceptedNodeSchedule (combNodes d) afterSide (gateTail q (d - q)) := by
        rw [hlen, gateTail]
        exact ⟨combGate_mem hqd, hgateFail, htail⟩
      have haccepted :
          AcceptedNodeSchedule (combNodes d) s (combBlockNodes d (q + 1)) :=
        ⟨combTrue_mem hqd, hsideFail, hrest⟩
      have hblockRun :
          runNodeSchedule s (combBlockNodes d (q + 1)) =
            runNodeSchedule afterGate (gateTail (q + 1) (d - (q + 1))) := by
        rw [show combBlockNodes d (q + 1) =
            combTrue q :: combGate q :: gateTail (q + 1) (d - (q + 1)) by
          simp [combBlockNodes, hlen, gateTail]]
        rfl
      refine ⟨haccepted, ?_⟩
      rw [hblockRun]
      refine ⟨by omega, ?_, ?_, ?_⟩
      · have hzeroSide : afterSide 0 = s 0 := by
          simp [afterSide, combTrue, combConstNode, repairNode, write]
        have hzeroGate : afterGate 0 = afterSide 0 := by
          simp [afterGate, combGate, combNandNode, repairNode, write]
        have hzeroTail := htailFrame 0 (by omega)
        exact hzeroTail.trans (hzeroGate.trans (hzeroSide.trans hzero))
      · intro r hr
        have hold := hunprocessed r (by omega)
        have hsideFrame : afterSide (2 * r + 1) = s (2 * r + 1) := by
          simp [afterSide, combTrue, combConstNode, repairNode, write]
          omega
        have hgateFrame : afterGate (2 * r + 1) = afterSide (2 * r + 1) := by
          simp [afterGate, combGate, combNandNode, repairNode, write]
          omega
        have htailSide := htailFrame (2 * r + 1) (by omega)
        have houtSide : afterSide (2 * r + 2) = s (2 * r + 2) := by
          simp [afterSide, combTrue, combConstNode, repairNode, write]
          omega
        have houtGate : afterGate (2 * r + 2) = afterSide (2 * r + 2) := by
          simp [afterGate, combGate, combNandNode, repairNode, write]
          omega
        have htailOut := htailFrame (2 * r + 2) (by omega)
        exact ⟨htailSide.trans (hgateFrame.trans (hsideFrame.trans hold.1)),
          htailOut.trans (houtGate.trans (houtSide.trans hold.2))⟩
      · intro r hqr hrd
        by_cases hrq : r = q
        · subst r
          have hsideTail := htailFrame (2 * q + 1) (by omega)
          have hgateAccept : (combGate q).obs.ok afterGate = true :=
            repairNode_accepts (combGate q) afterSide (combGate_deps_lt q)
          have hobsEq :
              (combGate q).obs.ok
                  (runNodeSchedule afterGate
                    (gateTail (q + 1) (d - (q + 1)))) =
                (combGate q).obs.ok afterGate := by
            apply (combGate q).obs.ok_local
            intro x hx
            apply htailFrame x
            simp [combGate, combNandNode, Node.obs] at hx
            rcases hx with rfl | rfl | rfl <;> omega
          have hsideAfterGate : afterGate (2 * q + 1) = true := by
            simp [afterGate, combGate, combNandNode, repairNode, write,
              hsideTrue]
          exact ⟨hsideTail.trans hsideAfterGate, hobsEq.trans hgateAccept⟩
        · have hr : q + 1 ≤ r := by omega
          obtain ⟨j, rfl⟩ := Nat.exists_eq_add_of_le hr
          simpa [Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
            htailProcessed j (by omega)

theorem combScheduleFrom_accepted
    {d i : Nat} {s : State} (hready : CombReady d i s) :
    AcceptedNodeSchedule (combNodes d) s (combScheduleFrom d i) ∧
    CombDone d (runNodeSchedule s (combScheduleFrom d i)) := by
  induction i generalizing s with
  | zero =>
      simpa [combScheduleFrom] using combBlock_transition hready
  | succ i ih =>
      have hblock := combBlock_transition hready
      have hreadyNext :
          CombReady d i (runNodeSchedule s (combBlockNodes d (i + 1))) := by
        simpa using hblock.2
      obtain ⟨hrest, hdone⟩ := ih hreadyNext
      refine ⟨acceptedNodeSchedule_append hblock.1 hrest, ?_⟩
      change CombDone d
        (runNodeSchedule s
          (combBlockNodes d (i + 1) ++ combScheduleFrom d i))
      rw [runNodeSchedule_append]
      exact hdone

theorem combAcceptedSchedule_valid (d : Nat) :
    AcceptedNodeSchedule (combNodes d) combInitial (combAcceptedSchedule d) ∧
    CombDone d (runNodeSchedule combInitial (combAcceptedSchedule d)) := by
  exact combScheduleFrom_accepted (combInitial_ready d)

theorem combDone_node_accepts {d : Nat} {s : State}
    (hdone : CombDone d s) :
    ∀ n ∈ combNodes d, n.obs.ok s = true := by
  induction d with
  | zero =>
      intro n hn
      simp [combNodes] at hn
      subst n
      simp [combConstNode, Node.obs, hdone.1]
  | succ d ih =>
      intro n hn
      simp only [combNodes, List.mem_append, List.mem_cons] at hn
      rcases hn with hn | rfl | rfl | hn
      · apply ih
        · exact ⟨hdone.1, fun q hq => hdone.2 q (by omega)⟩
        · exact hn
      · have hside := (hdone.2 d (by omega)).1
        simp [combConstNode, Node.obs, hside]
      · exact (hdone.2 d (by omega)).2
      · cases hn

theorem combDone_consensus {d : Nat} {s : State}
    (hdone : CombDone d s) :
    Consensus ((combNodes d).map Node.obs) s := by
  intro o ho
  rcases List.mem_map.mp ho with ⟨n, hn, rfl⟩
  exact combDone_node_accepts hdone n hn

theorem gateTail_length (q n : Nat) :
    (gateTail q n).length = n := by
  induction n generalizing q with
  | zero => rfl
  | succ n ih =>
      simp [gateTail, ih]

def combStepCount (d : Nat) : Nat → Nat
  | 0 => d + 1
  | i + 1 => d - i + 1 + combStepCount d i

theorem combBlockNodes_length_zero (d : Nat) :
    (combBlockNodes d 0).length = d + 1 := by
  simp [combBlockNodes, gateTail_length]

theorem combBlockNodes_length_succ (d i : Nat) :
    (combBlockNodes d (i + 1)).length = d - i + 1 := by
  simp [combBlockNodes, gateTail_length]

theorem combScheduleFrom_length (d i : Nat) :
    (combScheduleFrom d i).length = combStepCount d i := by
  induction i with
  | zero => exact combBlockNodes_length_zero d
  | succ i ih =>
      simp [combScheduleFrom, combStepCount,
        combBlockNodes_length_succ, ih]

theorem combStepCount_shift (d i : Nat) (hi : i ≤ d) :
    combStepCount (d + 1) i = combStepCount d i + (i + 1) := by
  induction i with
  | zero => simp [combStepCount]
  | succ i ih =>
      simp only [combStepCount]
      rw [ih (by omega)]
      omega

theorem combStepCount_diagonal (d : Nat) :
    combStepCount d d = lowerSteps d := by
  induction d with
  | zero => rfl
  | succ d ih =>
      simp only [combStepCount, lowerSteps]
      rw [combStepCount_shift d d le_rfl, ih]
      omega

theorem combAcceptedSchedule_length (d : Nat) :
    (combAcceptedSchedule d).length = lowerSteps d := by
  rw [combAcceptedSchedule, combScheduleFrom_length,
    combStepCount_diagonal]

/-- The recursive nested-NAND family has an actual accepted path whose exact
length is quadratic in the emitted program size. -/
theorem negComb_quadratic_lower (d : Nat) :
    ∃ t : State,
      CanonicalAcceptedSteps
        (fixedProgram (negComb d)) (lowerSteps d) combInitial t ∧
      Consensus (fixedFederation (negComb d)) t := by
  let t := runNodeSchedule combInitial (combAcceptedSchedule d)
  have hvalid := combAcceptedSchedule_valid d
  have hsteps := acceptedNodeSchedule_steps hvalid.1
  rw [combAcceptedSchedule_length] at hsteps
  refine ⟨t, ?_, ?_⟩
  · simpa [t, fixedProgram_negComb] using hsteps
  · simpa [t, fixedFederation, fixedProgram_negComb] using
      combDone_consensus hvalid.2

/-- Upper and lower certificates use the same metric, the number of nodes in
the actual `fixedProgram` compiler output. -/
theorem fixedProgram_sharp_quadratic_certificates :
    (∀ {k : Nat} (phi : Formula k) {m : Nat} {s t : State},
      CanonicalAcceptedSteps (fixedProgram phi) m s t →
      m ≤ (fixedProgram phi).length *
        ((fixedProgram phi).length + 1) / 2) ∧
    (∀ d : Nat, ∃ t : State,
      CanonicalAcceptedSteps
        (fixedProgram (negComb d)) (lowerSteps d) combInitial t ∧
      Consensus (fixedFederation (negComb d)) t ∧
      8 * lowerSteps d + 1 =
        (fixedProgram (negComb d)).length *
          (fixedProgram (negComb d)).length +
        8 * (fixedProgram (negComb d)).length) := by
  constructor
  · intro k phi m s t h
    exact fixedProgram_acceptedSteps_quadratic phi h
  · intro d
    obtain ⟨t, hsteps, hcons⟩ := negComb_quadratic_lower d
    exact ⟨t, hsteps, hcons, negComb_size_lowerSteps_identity d⟩

#print axioms negComb_quadratic_lower
#print axioms fixedProgram_sharp_quadratic_certificates

#print axioms combScheduleFrom_accepted
#print axioms combDone_consensus

#print axioms combBlock_transition

#print axioms combGate_fails_after_left_flip
#print axioms gateTail_accepted

#print axioms go_negComb
#print axioms fixedProgram_negComb_length
#print axioms combSchedule_length
#print axioms negComb_size_lowerSteps_identity

end


end OPH.RepairUniversality.FixedFederation
