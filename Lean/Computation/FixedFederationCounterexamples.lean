import Computation.FixedFederationProgress

/-!
# Fixed-federation negative controls

These controls keep the historical weak attempt relation separate from the
canonical accepted subrelation.  The first fixture is a one-member fair
stuttering execution.  The second fixture shows why ordinary defect count is
not the well-founded measure for a straight-line dependency federation.
-/

namespace OPH.RepairUniversality.FixedFederation

noncomputable section

/-- A scheduled witness for the historical weak semantics. -/
def WeakNodeAttempt (L : List Node) (selected : Node)
    (s t : State) : Prop :=
  selected ∈ L ∧ selected.obs.ok s = false ∧
    ∀ i, i ∉ selected.obs.patch → t i = s i

def WeakAttemptRun (L : List Node) (sigma : Nat → Node)
    (rho : Nat → State) : Prop :=
  ∀ n, WeakNodeAttempt L (sigma n) (rho n) (rho (n + 1))

/-- The node-specialized weak attempt is an actual step of the historical
`RepairStep` relation. -/
theorem weakNodeAttempt_to_repairStep {L : List Node} {selected : Node}
    {s t : State} (h : WeakNodeAttempt L selected s t) :
    RepairStep (L.map Node.obs) s t := by
  obtain ⟨hmember, hfail, hframe⟩ := h
  exact ⟨selected.obs, List.mem_map.mpr ⟨selected, hmember, rfl⟩,
    hfail, hframe⟩

theorem weakAttemptRun_isRepairStepRun {L : List Node}
    {sigma : Nat → Node} {rho : Nat → State}
    (h : WeakAttemptRun L sigma rho) :
    ∀ n, RepairStep (L.map Node.obs) (rho n) (rho (n + 1)) := by
  intro n
  exact weakNodeAttempt_to_repairStep (h n)

def MemberFair (L : List Node) (sigma : Nat → Node) : Prop :=
  ∀ n ∈ L, ∀ N, ∃ m, N ≤ m ∧ sigma m = n

def SiteFair (L : List Node) (sigma : Nat → Node) : Prop :=
  ∀ n ∈ L, ∀ N, ∃ m, N ≤ m ∧ (sigma m).reg = n.reg

def ContinuouslyEnabledFair (L : List Node) (sigma : Nat → Node)
    (rho : Nat → State) : Prop :=
  ∀ n ∈ L, ∀ N,
    (∀ m, N ≤ m → n.obs.ok (rho m) = false) →
    ∃ m, N ≤ m ∧ sigma m = n

def forceTrueNode : Node where
  reg := 0
  deps := []
  val := fun _ => true
  val_local := by intro _ _ _; rfl

def allFalse : State := fun _ => false

def weakStutterSchedule : Nat → Node := fun _ => forceTrueNode

def weakStutterRun : Nat → State := fun _ => allFalse

theorem weakStutter_isRun :
    WeakAttemptRun [forceTrueNode] weakStutterSchedule weakStutterRun := by
  intro n
  refine ⟨by simp [weakStutterSchedule], ?_, ?_⟩
  · simp [weakStutterSchedule, forceTrueNode, Node.obs,
      weakStutterRun, allFalse]
  · intro i _hi
    rfl

theorem weakStutter_isRepairStepRun :
    ∀ n, RepairStep [forceTrueNode.obs]
      (weakStutterRun n) (weakStutterRun (n + 1)) := by
  simpa using weakAttemptRun_isRepairStepRun weakStutter_isRun

theorem weakStutter_memberFair :
    MemberFair [forceTrueNode] weakStutterSchedule := by
  intro n hn N
  simp only [List.mem_singleton] at hn
  subst n
  exact ⟨N, le_rfl, rfl⟩

theorem weakStutter_siteFair :
    SiteFair [forceTrueNode] weakStutterSchedule := by
  intro n hn N
  simp only [List.mem_singleton] at hn
  subst n
  exact ⟨N, le_rfl, rfl⟩

theorem weakStutter_continuouslyEnabledFair :
    ContinuouslyEnabledFair [forceTrueNode]
      weakStutterSchedule weakStutterRun := by
  intro n hn N _henabled
  simp only [List.mem_singleton] at hn
  subst n
  exact ⟨N, le_rfl, rfl⟩

theorem weakStutter_never_consensus :
    ∀ n, ¬ Consensus [forceTrueNode.obs] (weakStutterRun n) := by
  intro n hcons
  have h := hcons forceTrueNode.obs (by simp)
  simp [forceTrueNode, Node.obs, weakStutterRun, allFalse] at h

/-- Even the conjunction of member recurrence, site recurrence, and
continuously-enabled fairness cannot force progress for the historical weak
step semantics, because the scheduled step itself may stutter. -/
theorem weak_fair_stuttering_no_go :
    ∃ (L : List Node) (sigma : Nat → Node) (rho : Nat → State),
      WeakAttemptRun L sigma rho ∧
      (∀ n, RepairStep (L.map Node.obs) (rho n) (rho (n + 1))) ∧
      MemberFair L sigma ∧
      SiteFair L sigma ∧
      ContinuouslyEnabledFair L sigma rho ∧
      (∀ n, ¬ Consensus (L.map Node.obs) (rho n)) := by
  exact ⟨[forceTrueNode], weakStutterSchedule, weakStutterRun,
    weakStutter_isRun, weakStutter_isRepairStepRun,
    weakStutter_memberFair, weakStutter_siteFair,
    weakStutter_continuouslyEnabledFair, weakStutter_never_consensus⟩

/-- Historical input pinning really produces different compiled devices. This
typed control is the contrast class excluded by `fixedProgram`. -/
def inputPinnedControlProgram (b : Bool) : List Node :=
  compile (Formula.const false) (fun _ : Fin 1 => b)

theorem inputPinnedControlProgram_depends_on_input :
    inputPinnedControlProgram false ≠ inputPinnedControlProgram true := by
  intro h
  have hhead := congrArg
    (fun L : List Node => L.head?.map fun n => n.val allFalse) h
  simp [inputPinnedControlProgram, compile, inputs, buildInputs] at hhead

def singletonCanonicalScheduler : NodeScheduler [forceTrueNode] :=
  fun _ _ => ⟨forceTrueNode, by simp⟩

def singletonStart : State := allFalse

theorem singletonStart_not_consensus :
    ¬ Consensus [forceTrueNode.obs] singletonStart := by
  intro hcons
  have h := hcons forceTrueNode.obs (by simp)
  simp [singletonStart, allFalse, forceTrueNode, Node.obs] at h

theorem singletonCanonical_pathwiseWeakFair :
    NodePathwiseWeakFair [forceTrueNode]
      singletonCanonicalScheduler singletonStart := by
  intro member N _henabled
  have hmember : member.1 = forceTrueNode := by
    simpa using member.2
  refine ⟨N, le_rfl, ?_⟩
  simp [singletonCanonicalScheduler, hmember]

theorem singletonCanonical_positive :
    ∃ N,
      Consensus [forceTrueNode.obs]
        (attemptRun [forceTrueNode] N singletonCanonicalScheduler singletonStart) ∧
      ∀ n, N ≤ n →
        attemptRun [forceTrueNode] n singletonCanonicalScheduler singletonStart =
          attemptRun [forceTrueNode] N singletonCanonicalScheduler singletonStart := by
  exact nodePathwiseWeakFair_eventually_consensus [forceTrueNode]
    (by simp [NodesWF, forceTrueNode]) singletonCanonicalScheduler singletonStart
    singletonCanonical_pathwiseWeakFair

/-- Count currently failing node observers.  This is deliberately not the
well-founded rank for dependency-ordered canonical repair. -/
def federationDefectCount (L : List Node) (s : State) : Nat :=
  (L.filter fun n => n.obs.ok s = false).length

def copyRootOne : Node where
  reg := 1
  deps := [0]
  val := fun s => s 0
  val_local := by
    intro s t h
    exact h 0 (by simp)

def copyRootTwo : Node where
  reg := 2
  deps := [0]
  val := fun s => s 0
  val_local := by
    intro s t h
    exact h 0 (by simp)

def fanoutNodes : List Node :=
  [forceTrueNode, copyRootOne, copyRootTwo]

def fanoutAfterRoot : State :=
  repairNode forceTrueNode allFalse

theorem fanoutNodes_wf : NodesWF fanoutNodes := by
  constructor
  · simp [fanoutNodes, forceTrueNode, copyRootOne, copyRootTwo]
  · intro n hn i hi
    simp [fanoutNodes] at hn
    rcases hn with rfl | rfl | rfl
    · simp [forceTrueNode] at hi
    · have hieq : i = 0 := by simpa [copyRootOne] using hi
      simp [copyRootOne, hieq]
    · have hieq : i = 0 := by simpa [copyRootTwo] using hi
      simp [copyRootTwo, hieq]

theorem fanout_root_step :
    CanonicalAcceptedStep fanoutNodes allFalse fanoutAfterRoot := by
  refine ⟨forceTrueNode, by simp [fanoutNodes], ?_, rfl⟩
  simp [forceTrueNode, Node.obs, allFalse]

/-- Repairing the one upstream defect can expose two downstream defects.
Hence raw failing-member count is not the global descent measure. -/
theorem fanout_defect_count_increases :
    federationDefectCount fanoutNodes allFalse = 1 ∧
    federationDefectCount fanoutNodes fanoutAfterRoot = 2 := by
  decide

def fanoutRootScheduler : NodeScheduler fanoutNodes :=
  fun _ _ => ⟨forceTrueNode, by simp [fanoutNodes]⟩

/-- A strictly weaker fairness notion that only protects nodes failing from
time zero.  It does not protect a node that becomes enabled later. -/
def NodeGloballyEnabledFair (L : List Node) (sigma : NodeScheduler L)
    (s : State) : Prop :=
  ∀ member : {n : Node // n ∈ L},
    (∀ m, member.1.obs.ok (attemptRun L m sigma s) = false) →
    ∀ N, ∃ m, N ≤ m ∧
      (sigma m (attemptRun L m sigma s)).1 = member.1

theorem fanoutRoot_run_succ (n : Nat) :
    attemptRun fanoutNodes (n + 1) fanoutRootScheduler allFalse =
      fanoutAfterRoot := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [attemptRun_last_step, ih]
      change repairNode forceTrueNode fanoutAfterRoot = fanoutAfterRoot
      apply repairNode_eq_self_of_accepts
      simp [fanoutAfterRoot, forceTrueNode, Node.obs, repairNode, write]

theorem fanoutAfterRoot_not_consensus :
    ¬ Consensus (fanoutNodes.map Node.obs) fanoutAfterRoot := by
  intro hcons
  have h := hcons copyRootOne.obs (by simp [fanoutNodes])
  simp [fanoutAfterRoot, forceTrueNode, copyRootOne,
    Node.obs, repairNode, write, allFalse] at h

theorem fanoutRoot_never_consensus_after_start (n : Nat) :
    ¬ Consensus (fanoutNodes.map Node.obs)
      (attemptRun fanoutNodes (n + 1) fanoutRootScheduler allFalse) := by
  rw [fanoutRoot_run_succ]
  exact fanoutAfterRoot_not_consensus

theorem fanoutRoot_globallyEnabledFair :
    NodeGloballyEnabledFair fanoutNodes fanoutRootScheduler allFalse := by
  intro member hglobal N
  have hcases := member.2
  simp [fanoutNodes] at hcases
  rcases hcases with hroot | hleft | hright
  · have hfalse : False := by
      have hfail := hglobal 1
      simp [hroot, fanoutRoot_run_succ, fanoutAfterRoot, forceTrueNode,
        Node.obs, repairNode, write] at hfail
    exact hfalse.elim
  · have hfalse : False := by
      have hfail := hglobal 0
      simp [hleft, attemptRun, copyRootOne, Node.obs, allFalse] at hfail
    exact hfalse.elim
  · have hfalse : False := by
      have hfail := hglobal 0
      simp [hright, attemptRun, copyRootTwo, Node.obs, allFalse] at hfail
    exact hfalse.elim

/-- Fairness restricted to nodes enabled from the initial instant is too weak:
the downstream nodes become failing after the root repair and are starved. -/
theorem globallyEnabledFair_no_go :
    ∃ (L : List Node) (sigma : NodeScheduler L) (s : State),
      NodesWF L ∧ NodeGloballyEnabledFair L sigma s ∧
      ∀ n, ¬ Consensus (L.map Node.obs)
        (attemptRun L (n + 1) sigma s) := by
  exact ⟨fanoutNodes, fanoutRootScheduler, allFalse, fanoutNodes_wf,
    fanoutRoot_globallyEnabledFair, fanoutRoot_never_consensus_after_start⟩

/-- A still weaker scheduler contract: a node that fails on a tail need only
have been selected once before that tail began. -/
def SelectedBeforeFailingTailFair (L : List Node) (sigma : NodeScheduler L)
    (s : State) : Prop :=
  ∀ member : {n : Node // n ∈ L}, ∀ N,
    (∀ n, N ≤ n →
      member.1.obs.ok (attemptRun L n sigma s) = false) →
    ∃ m, m < N ∧
      (sigma m (attemptRun L m sigma s)).1 = member.1

def preTailNodes : List Node := [forceTrueNode, copyRootOne]

def preTailScheduler : NodeScheduler preTailNodes
  | 0, _ => ⟨copyRootOne, by simp [preTailNodes]⟩
  | _ + 1, _ => ⟨forceTrueNode, by simp [preTailNodes]⟩

theorem preTail_run_one :
    attemptRun preTailNodes 1 preTailScheduler allFalse = allFalse := by
  rw [attemptRun_last_step]
  change repairNode copyRootOne allFalse = allFalse
  apply repairNode_eq_self_of_accepts
  simp [copyRootOne, Node.obs, allFalse]

theorem preTail_run_add_two (n : Nat) :
    attemptRun preTailNodes (n + 2) preTailScheduler allFalse =
      fanoutAfterRoot := by
  induction n with
  | zero =>
      rw [show 0 + 2 = 1 + 1 by omega, attemptRun_last_step,
        preTail_run_one]
      rfl
  | succ n ih =>
      rw [show (n + 1) + 2 = (n + 2) + 1 by omega,
        attemptRun_last_step, ih]
      change repairNode forceTrueNode fanoutAfterRoot = fanoutAfterRoot
      apply repairNode_eq_self_of_accepts
      simp [fanoutAfterRoot, forceTrueNode, Node.obs, repairNode, write]

theorem preTailScheduler_selectedBeforeFailingTailFair :
    SelectedBeforeFailingTailFair preTailNodes preTailScheduler allFalse := by
  intro member N htail
  have hcases := member.2
  simp [preTailNodes] at hcases
  rcases hcases with hroot | hcopy
  · have hfalse : False := by
      have hfail := htail (N + 2) (by omega)
      rw [preTail_run_add_two] at hfail
      simp [hroot, fanoutAfterRoot, forceTrueNode, Node.obs,
        repairNode, write] at hfail
    exact hfalse.elim
  · have htwo : 2 ≤ N := by
      by_contra hnot
      have hsmall : N = 0 ∨ N = 1 := by omega
      rcases hsmall with rfl | rfl
      · have hfail := htail 0 le_rfl
        simp [hcopy, attemptRun, copyRootOne, Node.obs, allFalse] at hfail
      · have hfail := htail 1 le_rfl
        rw [preTail_run_one] at hfail
        simp [hcopy, copyRootOne, Node.obs, allFalse] at hfail
    refine ⟨0, by omega, ?_⟩
    simp [preTailScheduler, hcopy]

theorem preTailScheduler_never_consensus_after_activation (n : Nat) :
    ¬ Consensus (preTailNodes.map Node.obs)
      (attemptRun preTailNodes (n + 2) preTailScheduler allFalse) := by
  rw [preTail_run_add_two]
  intro hcons
  have h := hcons copyRootOne.obs (by simp [preTailNodes])
  simp [fanoutAfterRoot, forceTrueNode, copyRootOne,
    Node.obs, repairNode, write, allFalse] at h

/-- One selection before activation cannot replace selection on the failing
tail: after the root repair, the downstream node is starved forever. -/
theorem selectedBeforeFailingTailFair_no_go :
    ∃ (L : List Node) (sigma : NodeScheduler L) (s : State),
      NodesWF L ∧ SelectedBeforeFailingTailFair L sigma s ∧
      ∀ n, ¬ Consensus (L.map Node.obs)
        (attemptRun L (n + 2) sigma s) := by
  exact ⟨preTailNodes, preTailScheduler, allFalse,
    (by simp [NodesWF, preTailNodes, forceTrueNode, copyRootOne]),
    preTailScheduler_selectedBeforeFailingTailFair,
    preTailScheduler_never_consensus_after_activation⟩

theorem fanoutRoot_not_pathwiseWeakFair :
    ¬ NodePathwiseWeakFair fanoutNodes fanoutRootScheduler allFalse := by
  intro hfair
  let member : {n : Node // n ∈ fanoutNodes} :=
    ⟨copyRootOne, by simp [fanoutNodes]⟩
  have hcontinuous : ∀ n, 1 ≤ n →
      member.1.obs.ok
        (attemptRun fanoutNodes n fanoutRootScheduler allFalse) = false := by
    intro n hn
    obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le hn
    rw [Nat.add_comm]
    simp only [member, fanoutRoot_run_succ]
    simp [fanoutAfterRoot, forceTrueNode, copyRootOne,
      Node.obs, repairNode, write, allFalse]
  obtain ⟨m, _hm, hselected⟩ := hfair member 1 hcontinuous
  have hreg := congrArg Node.reg hselected
  simp [fanoutRootScheduler, member, forceTrueNode, copyRootOne] at hreg

#print axioms weak_fair_stuttering_no_go
#print axioms weakAttemptRun_isRepairStepRun
#print axioms inputPinnedControlProgram_depends_on_input
#print axioms fanout_defect_count_increases
#print axioms singletonCanonical_positive
#print axioms singletonStart_not_consensus
#print axioms fanoutRoot_not_pathwiseWeakFair
#print axioms globallyEnabledFair_no_go
#print axioms selectedBeforeFailingTailFair_no_go

end

end OPH.RepairUniversality.FixedFederation
