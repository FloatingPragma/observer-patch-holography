import Computation.RepairUniversality
import ObservableNormalForms.ObserverConfluence

/-!
# Fixed-federation canonical repair progress

The historical `RepairStep` relation records only a failing observer and an
outside-patch frame condition.  This module leaves that weak relation intact
and introduces the deterministic accepted subrelation induced by a compiled
node.  Formula nodes, unlike the historical input-pinning nodes, form one
federation independent of the input carried by the initial register state.
-/

namespace OPH.RepairUniversality.FixedFederation

open Relation
open ObservableNormalForms

noncomputable section

/-- The deterministic repair attached to one straight-line constraint node. -/
def repairNode (n : Node) (s : State) : State :=
  write s n.reg (n.val s)

/-- A genuine canonical node repair.  The selected node must fail and the
successor is its declared deterministic write. -/
def CanonicalAcceptedStep (L : List Node) (s t : State) : Prop :=
  ∃ n ∈ L, n.obs.ok s = false ∧ t = repairNode n s

/-- Every canonical node repair is a step of the historical weak relation. -/
theorem canonicalAcceptedStep_to_repairStep {L : List Node} {s t : State}
    (h : CanonicalAcceptedStep L s t) :
    RepairStep (L.map Node.obs) s t := by
  obtain ⟨n, hn, hfail, rfl⟩ := h
  refine ⟨n.obs, List.mem_map.mpr ⟨n, hn, rfl⟩, hfail, ?_⟩
  intro i hi
  simp only [repairNode, write]
  exact if_neg fun hir => hi (by rw [hir]; exact List.mem_cons_self)

/-- Canonical normality is exactly consensus for the declared node list. -/
theorem isNormalForm_canonicalAcceptedStep_iff_consensus
    (L : List Node) (s : State) :
    IsNormalForm (CanonicalAcceptedStep L) s ↔
      Consensus (L.map Node.obs) s := by
  constructor
  · intro hnormal o ho
    rcases List.mem_map.mp ho with ⟨n, hn, rfl⟩
    by_cases hok : n.obs.ok s = true
    · exact hok
    · have hfail : n.obs.ok s = false := by simpa using hok
      exact False.elim (hnormal (repairNode n s) ⟨n, hn, hfail, rfl⟩)
  · intro hcons t hstep
    obtain ⟨n, hn, hfail, _⟩ := hstep
    have hok : n.obs.ok s = true :=
      hcons n.obs (List.mem_map.mpr ⟨n, hn, rfl⟩)
    rw [hok] at hfail
    contradiction

/-- A sweep over a sublist is reachable in the canonical relation of the
ambient node list. -/
theorem sweepFrom_reachable_in (F : List Node) (s : State) :
    ∀ L : List Node, (∀ n ∈ L, n ∈ F) →
      ReflTransGen (CanonicalAcceptedStep F) s (sweepFrom s L) := by
  intro L hsub
  induction L generalizing s with
  | nil => exact ReflTransGen.refl
  | cons n L ih =>
      by_cases hok : n.obs.ok s = true
      · have hforced : s n.reg = n.val s := by
          simpa [Node.obs, beq_iff_eq] using hok
        have hwrite : write s n.reg (n.val s) = s := by
          funext i
          simp only [write]
          split
          · next hi => rw [hi]; exact hforced.symm
          · rfl
        simpa only [sweepFrom, hwrite] using
          ih (s := s) (fun (m : Node) hm =>
            hsub m (List.mem_cons_of_mem n hm))
      · have hfail : n.obs.ok s = false := by simpa using hok
        have hstep : CanonicalAcceptedStep F s (repairNode n s) :=
          ⟨n, hsub n List.mem_cons_self, hfail, rfl⟩
        have htail := ih (s := repairNode n s)
          (fun m hm => hsub m (List.mem_cons_of_mem n hm))
        simpa only [sweepFrom, repairNode] using ReflTransGen.head hstep htail

theorem sweepFrom_reachable (L : List Node) (s : State) :
    ReflTransGen (CanonicalAcceptedStep L) s (sweepFrom s L) :=
  sweepFrom_reachable_in L s L (fun _ h => h)

/-- A well-formed fixed node federation is weakly normalizing under the
canonical accepted subrelation. -/
theorem canonicalAcceptedStep_weaklyNormalizing
    (L : List Node) (hWF : NodesWF L) (s : State) :
    WeaklyNormalizing (CanonicalAcceptedStep L) s := by
  refine ⟨sweepFrom s L, sweepFrom_reachable L s, ?_⟩
  exact (isNormalForm_canonicalAcceptedStep_iff_consensus L _).2
    (sweep_consensus s L hWF)

/-- The lexicographic defect word encoded as a natural number.  Earlier
nodes are more significant than all later nodes together. -/
def defectRank : List Node → State → Nat
  | [], _ => 0
  | n :: L, s =>
      (if n.obs.ok s = false then 2 ^ L.length else 0) + defectRank L s

theorem defectRank_lt_pow (L : List Node) (s : State) :
    defectRank L s < 2 ^ L.length := by
  induction L with
  | nil => simp [defectRank]
  | cons n L ih =>
      have hpow : 0 < 2 ^ L.length := pow_pos (by omega) _
      by_cases hfail : n.obs.ok s = false
      · simp [defectRank, hfail, pow_succ]
        omega
      · simp [defectRank, hfail, pow_succ]
        omega

theorem repairNode_accepts (n : Node) (s : State)
    (hdeps : ∀ i ∈ n.deps, i < n.reg) :
    n.obs.ok (repairNode n s) = true := by
  have hval : n.val (repairNode n s) = n.val s := by
    apply n.val_local
    intro i hi
    simp [repairNode, write, Nat.ne_of_lt (hdeps i hi)]
  change ((repairNode n s) n.reg == n.val (repairNode n s)) = true
  rw [hval]
  simp [repairNode, write]

theorem earlier_observer_unchanged (a b : Node) (s : State)
    (hab : a.reg < b.reg)
    (hdeps : ∀ i ∈ a.deps, i < a.reg) :
    a.obs.ok (repairNode b s) = a.obs.ok s := by
  have hreg : repairNode b s a.reg = s a.reg := by
    simp [repairNode, write, Nat.ne_of_lt hab]
  have hval : a.val (repairNode b s) = a.val s := by
    apply a.val_local
    intro i hi
    simp [repairNode, write,
      Nat.ne_of_lt (lt_trans (hdeps i hi) hab)]
  simp [Node.obs, hreg, hval]

/-- Every genuine canonical node repair strictly decreases the source-derived
dependency rank, even when it exposes arbitrarily many later defects. -/
theorem canonicalAcceptedStep_defectRank_lt
    (L : List Node) (hWF : NodesWF L) {s t : State}
    (hstep : CanonicalAcceptedStep L s t) :
    defectRank L t < defectRank L s := by
  induction L with
  | nil =>
      obtain ⟨n, hn, _hfail, _⟩ := hstep
      cases hn
  | cons a L ih =>
      obtain ⟨n, hn, hfail, rfl⟩ := hstep
      rcases hWF with ⟨hpw, hdepsAll⟩
      rw [List.pairwise_cons] at hpw
      rcases List.mem_cons.mp hn with hna | hn
      · subst n
        have hpost := repairNode_accepts a s
          (hdepsAll a List.mem_cons_self)
        have htail := defectRank_lt_pow L (repairNode a s)
        simp [defectRank, hfail, hpost]
        omega
      · have hhead := earlier_observer_unchanged a n s
          (hpw.1 n hn) (hdepsAll a List.mem_cons_self)
        have htailWF : NodesWF L := by
          exact ⟨hpw.2,
            fun p hp => hdepsAll p (List.mem_cons_of_mem a hp)⟩
        have htailStep : CanonicalAcceptedStep L s (repairNode n s) :=
          ⟨n, hn, hfail, rfl⟩
        have htail := ih htailWF htailStep
        simp only [defectRank, hhead]
        exact Nat.add_lt_add_left htail _

/-- Accepted-step termination needs no fairness. -/
theorem canonicalAcceptedStep_wellFounded
    (L : List Node) (hWF : NodesWF L) :
    WellFounded (fun t s => CanonicalAcceptedStep L s t) := by
  exact Subrelation.wf
    (fun {_ _} h => canonicalAcceptedStep_defectRank_lt L hWF h)
    (InvImage.wf _ wellFounded_lt)

/-- A scheduler for attempts on one fixed node federation. -/
def NodeScheduler (L : List Node) :=
  Nat → State → {n : Node // n ∈ L}

/-- Attempts may stutter when the selected member already accepts. -/
def attemptRun (L : List Node) : Nat → NodeScheduler L → State → State
  | 0, _sigma, s => s
  | n + 1, sigma, s =>
      let current := attemptRun L n sigma s
      repairNode (sigma n current).1 current

theorem attemptRun_last_step (L : List Node) (n : Nat)
    (sigma : NodeScheduler L) (s : State) :
    attemptRun L (n + 1) sigma s =
      repairNode (sigma n (attemptRun L n sigma s)).1
        (attemptRun L n sigma s) := rfl

theorem attemptRun_add (L : List Node) (m n : Nat)
    (sigma : NodeScheduler L) (s : State) :
    attemptRun L (m + n) sigma s =
      attemptRun L n (fun k x => sigma (m + k) x)
        (attemptRun L m sigma s) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [Nat.add_succ, attemptRun_last_step,
        attemptRun_last_step, ih]

theorem repairNode_eq_self_of_accepts (n : Node) (s : State)
    (hok : n.obs.ok s = true) : repairNode n s = s := by
  have hforced : s n.reg = n.val s := by
    simpa [Node.obs, beq_iff_eq] using hok
  funext i
  simp only [repairNode, write]
  split
  · next hi => rw [hi]; exact hforced.symm
  · rfl

/-- Once every declared node accepts, arbitrary further attempts stutter. -/
theorem attemptRun_eq_self_of_consensus
    (L : List Node) (sigma : NodeScheduler L) (s : State)
    (hcons : Consensus (L.map Node.obs) s) :
    ∀ steps, attemptRun L steps sigma s = s := by
  intro steps
  induction steps with
  | zero => rfl
  | succ n ih =>
      rw [attemptRun_last_step, ih]
      apply repairNode_eq_self_of_accepts
      exact hcons (sigma n s).1.obs
        (List.mem_map.mpr ⟨(sigma n s).1, (sigma n s).2, rfl⟩)

theorem repairNode_ne_self_of_fails (n : Node) (s : State)
    (hfail : n.obs.ok s = false) : repairNode n s ≠ s := by
  intro heq
  have hreg := congrFun heq n.reg
  have hforced : s n.reg = n.val s := by
    simpa [repairNode, write] using hreg.symm
  have hok : n.obs.ok s = true := by
    simp [Node.obs, hforced]
  rw [hok] at hfail
  contradiction

theorem attemptRun_step_rank_le (L : List Node) (hWF : NodesWF L)
    (n : Nat) (sigma : NodeScheduler L) (s : State) :
    defectRank L (attemptRun L (n + 1) sigma s) ≤
      defectRank L (attemptRun L n sigma s) := by
  rw [attemptRun_last_step]
  let selected := (sigma n (attemptRun L n sigma s)).1
  by_cases hchange :
      repairNode selected (attemptRun L n sigma s) = attemptRun L n sigma s
  · rw [hchange]
  · have hfail : selected.obs.ok (attemptRun L n sigma s) = false := by
      cases hok : selected.obs.ok (attemptRun L n sigma s) with
      | false => rfl
      | true =>
          exact False.elim (hchange
            (repairNode_eq_self_of_accepts selected _ hok))
    exact Nat.le_of_lt (canonicalAcceptedStep_defectRank_lt L hWF
      ⟨selected, (sigma n (attemptRun L n sigma s)).2,
        hfail, rfl⟩)

theorem attemptRun_change_rank_lt (L : List Node) (hWF : NodesWF L)
    (n : Nat) (sigma : NodeScheduler L) (s : State)
    (hchange : attemptRun L (n + 1) sigma s ≠
      attemptRun L n sigma s) :
    defectRank L (attemptRun L (n + 1) sigma s) <
      defectRank L (attemptRun L n sigma s) := by
  rw [attemptRun_last_step] at hchange ⊢
  let selected := (sigma n (attemptRun L n sigma s)).1
  have hfail : selected.obs.ok (attemptRun L n sigma s) = false := by
    cases hok : selected.obs.ok (attemptRun L n sigma s) with
    | false => rfl
    | true =>
        exact False.elim (hchange
          (repairNode_eq_self_of_accepts selected _ hok))
  exact canonicalAcceptedStep_defectRank_lt L hWF
    ⟨selected, (sigma n (attemptRun L n sigma s)).2,
      hfail, rfl⟩

theorem attemptRun_rank_le_initial (L : List Node) (hWF : NodesWF L)
    (n : Nat) (sigma : NodeScheduler L) (s : State) :
    defectRank L (attemptRun L n sigma s) ≤ defectRank L s := by
  induction n with
  | zero => rfl
  | succ n ih =>
      exact le_trans (attemptRun_step_rank_le L hWF n sigma s) ih

theorem attemptRun_reachable (L : List Node) (n : Nat)
    (sigma : NodeScheduler L) (s : State) :
    ReflTransGen (CanonicalAcceptedStep L) s (attemptRun L n sigma s) := by
  induction n with
  | zero => exact ReflTransGen.refl
  | succ n ih =>
      rw [attemptRun_last_step]
      let selected := (sigma n (attemptRun L n sigma s)).1
      by_cases hchange :
          repairNode selected (attemptRun L n sigma s) = attemptRun L n sigma s
      · change ReflTransGen (CanonicalAcceptedStep L) s
          (repairNode selected (attemptRun L n sigma s))
        rw [hchange]
        exact ih
      · have hfail : selected.obs.ok (attemptRun L n sigma s) = false := by
          cases hok : selected.obs.ok (attemptRun L n sigma s) with
          | false => rfl
          | true =>
              exact False.elim (hchange
                (repairNode_eq_self_of_accepts selected _ hok))
        exact ReflTransGen.tail ih
          ⟨selected, (sigma n (attemptRun L n sigma s)).2,
            hfail, rfl⟩

theorem attemptRun_eventually_constant
    (L : List Node) (hWF : NodesWF L) :
    ∀ (s : State) (sigma : NodeScheduler L),
      ∃ N, ∀ n, N ≤ n →
        attemptRun L n sigma s = attemptRun L N sigma s := by
  intro s sigma
  generalize hd : defectRank L s = d
  induction d using Nat.strong_induction_on generalizing s sigma with
  | h d ih =>
      by_cases hchange :
          ∃ k, attemptRun L (k + 1) sigma s ≠ attemptRun L k sigma s
      · obtain ⟨k, hk⟩ := hchange
        let y : State := attemptRun L (k + 1) sigma s
        let sigmaTail : NodeScheduler L :=
          fun j x => sigma ((k + 1) + j) x
        have hylt : defectRank L y < d := by
          have hstrict := attemptRun_change_rank_lt L hWF k sigma s hk
          have hprefix := attemptRun_rank_le_initial L hWF k sigma s
          dsimp [y]
          omega
        obtain ⟨N, hN⟩ :=
          ih (defectRank L y) hylt y sigmaTail rfl
        refine ⟨(k + 1) + N, ?_⟩
        intro n hn
        obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le hn
        calc
          attemptRun L (((k + 1) + N) + r) sigma s =
              attemptRun L ((k + 1) + (N + r)) sigma s := by
                congr 1
                omega
          _ = attemptRun L (N + r) sigmaTail y := by
                simpa [y, sigmaTail] using
                  attemptRun_add L (k + 1) (N + r) sigma s
          _ = attemptRun L N sigmaTail y :=
                hN (N + r) (Nat.le_add_right N r)
          _ = attemptRun L ((k + 1) + N) sigma s := by
                simpa [y, sigmaTail] using
                  (attemptRun_add L (k + 1) N sigma s).symm
      · have hsame : ∀ k,
            attemptRun L (k + 1) sigma s = attemptRun L k sigma s := by
          intro k
          exact not_ne_iff.mp (not_exists.mp hchange k)
        refine ⟨0, ?_⟩
        intro n _hn
        induction n with
        | zero => rfl
        | succ n ihN => rw [hsame n, ihN (Nat.zero_le n)]

/-- Path-relative weak fairness: a node that remains failing on an entire
tail is eventually selected on that tail. -/
def NodePathwiseWeakFair (L : List Node) (sigma : NodeScheduler L)
    (s : State) : Prop :=
  ∀ member : {n : Node // n ∈ L}, ∀ N,
    (∀ n, N ≤ n →
      member.1.obs.ok (attemptRun L n sigma s) = false) →
    ∃ m, N ≤ m ∧
      (sigma m (attemptRun L m sigma s)).1 = member.1

theorem nodePathwiseWeakFair_eventually_consensus
    (L : List Node) (hWF : NodesWF L)
    (sigma : NodeScheduler L) (s : State)
    (hfair : NodePathwiseWeakFair L sigma s) :
    ∃ N, Consensus (L.map Node.obs) (attemptRun L N sigma s) ∧
      ∀ n, N ≤ n →
        attemptRun L n sigma s = attemptRun L N sigma s := by
  obtain ⟨N, hconstant⟩ := attemptRun_eventually_constant L hWF s sigma
  refine ⟨N, ?_, hconstant⟩
  by_contra hcons
  have hnotnormal :
      ¬ IsNormalForm (CanonicalAcceptedStep L) (attemptRun L N sigma s) := by
    intro hnormal
    exact hcons ((isNormalForm_canonicalAcceptedStep_iff_consensus L _).1 hnormal)
  simp only [IsNormalForm] at hnotnormal
  push Not at hnotnormal
  obtain ⟨t, ht⟩ := hnotnormal
  obtain ⟨member, hmember, hfail, _⟩ := ht
  let selected : {n : Node // n ∈ L} := ⟨member, hmember⟩
  have hcontinuous : ∀ n, N ≤ n →
      selected.1.obs.ok (attemptRun L n sigma s) = false := by
    intro n hn
    rw [hconstant n hn]
    exact hfail
  obtain ⟨m, hm, hselected⟩ := hfair selected N hcontinuous
  have hchange : attemptRun L (m + 1) sigma s ≠
      attemptRun L m sigma s := by
    rw [attemptRun_last_step, hselected]
    exact repairNode_ne_self_of_fails member _ (hcontinuous m hm)
  apply hchange
  rw [hconstant m hm,
    hconstant (m + 1) (Nat.le_trans hm (Nat.le_succ m))]

/-- Formula nodes only.  Input registers are read from the initial state, so
this list and its observer federation do not depend on an input value. -/
def fixedProgram {k : Nat} (phi : Formula k) : List Node :=
  (go phi k).nodes

def fixedFederation {k : Nat} (phi : Formula k) : Federation :=
  (fixedProgram phi).map Node.obs

def fixedOutReg {k : Nat} (phi : Formula k) : Nat :=
  (go phi k).out

def InputObservation (k : Nat) (s : State) : Fin k → Bool :=
  fun j => s j.val

def CarriesInput {k : Nat} (x : Fin k → Bool) (s : State) : Prop :=
  InputObservation k s = x

theorem fixedProgram_wf {k : Nat} (phi : Formula k) :
    NodesWF (fixedProgram phi) := by
  have F := go_facts phi k le_rfl (fun _ => false)
  exact ⟨F.pairwise, F.deps_lt⟩

theorem fixedProgram_inputSeparated {k : Nat} (phi : Formula k) :
    ∀ n ∈ fixedProgram phi, k ≤ n.reg := by
  intro n hn
  exact (go_facts phi k le_rfl (fun _ => false)).range n hn |>.1

theorem fixedSweep_carriesInput {k : Nat} (phi : Formula k)
    (x : Fin k → Bool) (s : State) (hinput : CarriesInput x s) :
    CarriesInput x (sweepFrom s (fixedProgram phi)) := by
  funext j
  change sweepFrom s (fixedProgram phi) j.val = x j
  rw [sweepFrom_frame]
  · exact congrFun hinput j
  · intro n hn
    exact ne_of_gt (lt_of_lt_of_le j.isLt (fixedProgram_inputSeparated phi n hn))

theorem fixedSweep_consensus {k : Nat} (phi : Formula k) (s : State) :
    Consensus (fixedFederation phi) (sweepFrom s (fixedProgram phi)) :=
  sweep_consensus s (fixedProgram phi) (fixedProgram_wf phi)

theorem fixedSweep_output {k : Nat} (phi : Formula k)
    (x : Fin k → Bool) (s : State) (hinput : CarriesInput x s) :
    sweepFrom s (fixedProgram phi) (fixedOutReg phi) =
      Formula.evalF phi x := by
  exact (go_facts phi k le_rfl x).sweep_out s (fun j => congrFun hinput j)

theorem fixedSweep_reachable {k : Nat} (phi : Formula k) (s : State) :
    ReflTransGen (CanonicalAcceptedStep (fixedProgram phi)) s
      (sweepFrom s (fixedProgram phi)) :=
  sweepFrom_reachable (fixedProgram phi) s

/-- Two consensus states of the same fixed program and with the same input
observation agree on the output register. -/
theorem fixedConsensus_output_unique {k : Nat} (phi : Formula k)
    (s t : State)
    (hs : Consensus (fixedFederation phi) s)
    (ht : Consensus (fixedFederation phi) t)
    (hinput : InputObservation k s = InputObservation k t) :
    s (fixedOutReg phi) = t (fixedOutReg phi) := by
  have F := go_facts phi k le_rfl (InputObservation k s)
  have key : ∀ r : Nat, ∀ n ∈ fixedProgram phi, n.reg = r →
      s n.reg = t n.reg := by
    intro r
    induction r using Nat.strong_induction_on with
    | _ r ih =>
      intro n hn hr
      rw [consensus_forced hs n hn, consensus_forced ht n hn]
      apply n.val_local
      intro i hi
      rcases F.deps_cov n hn i hi with hExternal | hInternal
      · exact congrFun hinput ⟨i, hExternal⟩
      · obtain ⟨m, hm, hmr⟩ := List.mem_map.mp hInternal
        have hlt : i < n.reg := F.deps_lt n hn i hi
        have hsame := ih i (hr ▸ hlt) m hm hmr
        rwa [hmr] at hsame
  obtain ⟨n, hn, hreg⟩ := List.mem_map.mp F.out_mem
  simpa [fixedOutReg, hreg] using key (go phi k).out n hn hreg

/-- The protected input observation identifies all fixed-program consensus
states modulo equality of their output register. -/
theorem fixedBoundaryIdentifiesOutput {k : Nat} (phi : Formula k) :
    BoundaryIdentifiesModulo
      {s | Consensus (fixedFederation phi) s}
      (InputObservation k)
      (fun s t => s (fixedOutReg phi) = t (fixedOutReg phi)) := by
  intro s t hs ht hinput
  exact fixedConsensus_output_unique phi s t hs ht hinput

/-- The canonical fixed-program relation preserves every input register. -/
theorem fixedCanonicalStep_observationPreserving {k : Nat}
    (phi : Formula k) :
    ObservationPreserving
      (CanonicalAcceptedStep (fixedProgram phi))
      (InputObservation k) := by
  intro s t hstep
  obtain ⟨n, hn, _hfail, rfl⟩ := hstep
  funext j
  simp only [InputObservation, repairNode, write]
  exact (if_neg (Nat.ne_of_lt
    (lt_of_lt_of_le j.isLt (fixedProgram_inputSeparated phi n hn)))).symm

theorem fixedRepairNode_preservesInputObservation {k : Nat}
    (phi : Formula k) (n : Node) (hn : n ∈ fixedProgram phi)
    (s : State) :
    InputObservation k (repairNode n s) = InputObservation k s := by
  funext j
  simp only [InputObservation, repairNode, write]
  exact if_neg (Nat.ne_of_lt
    (lt_of_lt_of_le j.isLt (fixedProgram_inputSeparated phi n hn)))

theorem fixedAttemptRun_preservesInputObservation {k : Nat}
    (phi : Formula k) (steps : Nat)
    (sigma : NodeScheduler (fixedProgram phi)) (s : State) :
    InputObservation k (attemptRun (fixedProgram phi) steps sigma s) =
      InputObservation k s := by
  induction steps with
  | zero => rfl
  | succ n ih =>
      rw [attemptRun_last_step]
      exact (fixedRepairNode_preservesInputObservation phi
        (sigma n (attemptRun (fixedProgram phi) n sigma s)).1
        (sigma n (attemptRun (fixedProgram phi) n sigma s)).2 _).trans ih

theorem fixedCanonicalStep_completeForConsensus {k : Nat}
    (phi : Formula k) :
    CompleteFor
      (CanonicalAcceptedStep (fixedProgram phi))
      {s | Consensus (fixedFederation phi) s} :=
  isNormalForm_canonicalAcceptedStep_iff_consensus (fixedProgram phi)

/-- Same-input normal endpoints of the fixed federation have the same output.
This is the direct observation-determined-normal-forms consumer. -/
theorem fixedObserverEndpointUniqueOutput {k : Nat} (phi : Formula k) :
    ObserverEndpointUniqueModulo
      (CanonicalAcceptedStep (fixedProgram phi))
      (InputObservation k)
      (fun s t => s (fixedOutReg phi) = t (fixedOutReg phi)) :=
  (boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    (fixedCanonicalStep_observationPreserving phi)
    (fixedCanonicalStep_completeForConsensus phi)).mp
      (fixedBoundaryIdentifiesOutput phi)

/-- One fixed federation computes the same Boolean function under every
pathwise weak-fair attempt scheduler.  Fairness removes attempt stutter;
accepted-step termination itself follows from the dependency rank. -/
theorem fixed_federation_fair_universality (k : Nat)
    (f : (Fin k → Bool) → Bool) :
    ∃ phi : Formula k, ∀ (x : Fin k → Bool) (s : State)
      (sigma : NodeScheduler (fixedProgram phi)),
      CarriesInput x s →
      NodePathwiseWeakFair (fixedProgram phi) sigma s →
      ∃ N,
        Consensus (fixedFederation phi)
          (attemptRun (fixedProgram phi) N sigma s) ∧
        (∀ n, N ≤ n →
          attemptRun (fixedProgram phi) n sigma s =
            attemptRun (fixedProgram phi) N sigma s) ∧
        CarriesInput x (attemptRun (fixedProgram phi) N sigma s) ∧
        attemptRun (fixedProgram phi) N sigma s (fixedOutReg phi) = f x := by
  obtain ⟨phi, hphi⟩ := Formula.exists_formula f
  refine ⟨phi, ?_⟩
  intro x s sigma hinput hfair
  obtain ⟨N, hcons, hconstant⟩ :=
    nodePathwiseWeakFair_eventually_consensus
      (fixedProgram phi) (fixedProgram_wf phi) sigma s hfair
  have hrunInput :
      CarriesInput x (attemptRun (fixedProgram phi) N sigma s) := by
    exact (fixedAttemptRun_preservesInputObservation phi N sigma s).trans hinput
  have hsweepInput := fixedSweep_carriesInput phi x s hinput
  have houtput := fixedConsensus_output_unique phi
    (attemptRun (fixedProgram phi) N sigma s)
    (sweepFrom s (fixedProgram phi)) hcons (fixedSweep_consensus phi s)
    (hrunInput.trans hsweepInput.symm)
  refine ⟨N, hcons, hconstant, hrunInput, ?_⟩
  exact houtput.trans ((fixedSweep_output phi x s hinput).trans (hphi x))

/-- Every Boolean function has one fixed federation whose canonical repair
computes it for all inputs carried by initial states. -/
theorem fixed_federation_universality (k : Nat)
    (f : (Fin k → Bool) → Bool) :
    ∃ phi : Formula k, ∀ (x : Fin k → Bool) (s : State),
      CarriesInput x s →
      ReflTransGen (CanonicalAcceptedStep (fixedProgram phi)) s
        (sweepFrom s (fixedProgram phi)) ∧
      Consensus (fixedFederation phi) (sweepFrom s (fixedProgram phi)) ∧
      CarriesInput x (sweepFrom s (fixedProgram phi)) ∧
      sweepFrom s (fixedProgram phi) (fixedOutReg phi) = f x := by
  obtain ⟨phi, hphi⟩ := Formula.exists_formula f
  refine ⟨phi, ?_⟩
  intro x s hinput
  exact ⟨fixedSweep_reachable phi s, fixedSweep_consensus phi s,
    fixedSweep_carriesInput phi x s hinput,
    (fixedSweep_output phi x s hinput).trans (hphi x)⟩

#print axioms fixed_federation_universality
#print axioms fixed_federation_fair_universality
#print axioms canonicalAcceptedStep_wellFounded
#print axioms fixedObserverEndpointUniqueOutput

end

end OPH.RepairUniversality.FixedFederation
