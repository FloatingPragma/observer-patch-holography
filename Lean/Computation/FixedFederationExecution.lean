import Computation.FixedFederationCounterexamples

/-!
# Fixed-federation execution classification

This module classifies the qualitative fairness and genuine accepted-step
cost of the input-independent formula federation.  Its scheduler is a
mathematical function.  No runtime, clock, rate, arbitration, or hardware
realization is constructed here.
-/

namespace OPH.RepairUniversality.FixedFederation

open OPH.RepairUniversality

noncomputable section

/-! ## Exact tail-fairness boundary -/

/-- Stable consensus makes the forever-failing antecedent of pathwise weak
fairness false.  This reverse direction is a structural, vacuous implication;
it does not establish recurrent selection. -/
theorem stableConsensusTail_nodePathwiseWeakFair
    (L : List Node) (sigma : NodeScheduler L) (s : State)
    (hstable :
      ∃ K,
        Consensus (L.map Node.obs) (attemptRun L K sigma s) ∧
        ∀ n, K ≤ n →
          attemptRun L n sigma s = attemptRun L K sigma s) :
    NodePathwiseWeakFair L sigma s := by
  rcases hstable with ⟨K, hcons, hconst⟩
  intro member N hfail
  have ht :
      member.1.obs.ok
        (attemptRun L (max N K) sigma s) = true := by
    rw [hconst (max N K) (le_max_right N K)]
    exact hcons member.1.obs
      (List.mem_map.mpr ⟨member.1, member.2, rfl⟩)
  have hf := hfail (max N K) (le_max_left N K)
  rw [ht] at hf
  contradiction

/-- On a well-formed fixed federation, tail-relative weak fairness is exactly
the condition that upgrades eventual attempt stabilization to consensus.  The
reverse implication remains vacuous after consensus. -/
theorem nodePathwiseWeakFair_iff_eventuallyStableConsensus
    (L : List Node) (hWF : NodesWF L)
    (sigma : NodeScheduler L) (s : State) :
    NodePathwiseWeakFair L sigma s ↔
      ∃ K,
        Consensus (L.map Node.obs) (attemptRun L K sigma s) ∧
        ∀ n, K ≤ n →
          attemptRun L n sigma s = attemptRun L K sigma s := by
  constructor
  · exact nodePathwiseWeakFair_eventually_consensus L hWF sigma s
  · exact stableConsensusTail_nodePathwiseWeakFair L sigma s

/-! ## Recurrence at member and register-site granularity -/

def NodeMemberRecurrent (L : List Node) (sigma : NodeScheduler L)
    (s : State) : Prop :=
  ∀ member : {n : Node // n ∈ L}, ∀ N,
    ∃ m, N ≤ m ∧
      (sigma m (attemptRun L m sigma s)).1 = member.1

def NodeSiteRecurrent (L : List Node) (sigma : NodeScheduler L)
    (s : State) : Prop :=
  ∀ member : {n : Node // n ∈ L}, ∀ N,
    ∃ m, N ≤ m ∧
      (sigma m (attemptRun L m sigma s)).1.reg = member.1.reg

def RegistersUnique (L : List Node) : Prop :=
  ∀ ⦃a : Node⦄, a ∈ L → ∀ ⦃b : Node⦄, b ∈ L →
    a.reg = b.reg → a = b

theorem nodeMemberRecurrent_nodeSiteRecurrent {L : List Node}
    {sigma : NodeScheduler L} {s : State}
    (h : NodeMemberRecurrent L sigma s) :
    NodeSiteRecurrent L sigma s := by
  intro member N
  obtain ⟨m, hm, hselected⟩ := h member N
  exact ⟨m, hm, congrArg Node.reg hselected⟩

theorem nodeMemberRecurrent_nodePathwiseWeakFair {L : List Node}
    {sigma : NodeScheduler L} {s : State}
    (h : NodeMemberRecurrent L sigma s) :
    NodePathwiseWeakFair L sigma s := by
  intro member N _hfail
  exact h member N

theorem nodesWF_registersUnique {L : List Node} (hWF : NodesWF L) :
    RegistersUnique L := by
  induction L with
  | nil =>
      intro a ha
      cases ha
  | cons head tail ih =>
      rw [NodesWF, List.pairwise_cons] at hWF
      rcases hWF with ⟨⟨hhead, htail⟩, hdeps⟩
      have ihUnique : RegistersUnique tail :=
        ih ⟨htail, fun n hn => hdeps n (List.mem_cons_of_mem head hn)⟩
      intro a ha b hb hreg
      rcases List.mem_cons.mp ha with rfl | ha
      · rcases List.mem_cons.mp hb with rfl | hb
        · rfl
        · have hlt := hhead b hb
          omega
      · have hheadA : head.reg < a.reg := hhead a ha
        rcases List.mem_cons.mp hb with rfl | hb
        · exact False.elim (by omega)
        · exact ihUnique ha hb hreg

theorem nodeSiteRecurrent_nodeMemberRecurrent_of_registersUnique
    {L : List Node} {sigma : NodeScheduler L} {s : State}
    (hunique : RegistersUnique L)
    (h : NodeSiteRecurrent L sigma s) :
    NodeMemberRecurrent L sigma s := by
  intro member N
  obtain ⟨m, hm, hreg⟩ := h member N
  refine ⟨m, hm, ?_⟩
  exact hunique (sigma m (attemptRun L m sigma s)).2 member.2 hreg

theorem nodeMemberRecurrent_iff_nodeSiteRecurrent
    (L : List Node) (hWF : NodesWF L)
    (sigma : NodeScheduler L) (s : State) :
    NodeMemberRecurrent L sigma s ↔ NodeSiteRecurrent L sigma s := by
  constructor
  · exact nodeMemberRecurrent_nodeSiteRecurrent
  · exact nodeSiteRecurrent_nodeMemberRecurrent_of_registersUnique
      (nodesWF_registersUnique hWF)

/-! ## Explicit state-blind fair scheduler -/

theorem fixedProgram_nonempty {k : Nat} (phi : Formula k) :
    fixedProgram phi ≠ [] := by
  have hout := (go_facts phi k le_rfl (fun _ => false)).out_mem
  rcases List.mem_map.mp hout with ⟨n, hn, _⟩
  exact List.ne_nil_of_mem hn

def roundRobinScheduler (L : List Node) (hL : L ≠ []) :
    NodeScheduler L :=
  fun step _ =>
    let i : Fin L.length :=
      ⟨step % L.length,
        Nat.mod_lt step (List.length_pos_iff.mpr hL)⟩
    ⟨L.get i, List.get_mem L i⟩

/-- Before the first wraparound, round robin executes exactly the source-order
prefix sweep. -/
theorem roundRobinScheduler_prefix_eq_sweepFrom
    (L : List Node) (hL : L ≠ []) (s : State) (steps : Nat)
    (hsteps : steps ≤ L.length) :
    attemptRun L steps (roundRobinScheduler L hL) s =
      sweepFrom s (L.take steps) := by
  induction steps with
  | zero => rfl
  | succ n ih =>
      have hnlt : n < L.length := Nat.lt_of_succ_le hsteps
      rw [attemptRun_last_step, ih (Nat.le_of_lt hnlt)]
      change repairNode
        (L.get ⟨n % L.length,
          Nat.mod_lt n (List.length_pos_iff.mpr hL)⟩)
        (sweepFrom s (L.take n)) = _
      have hindex :
          (⟨n % L.length,
            Nat.mod_lt n (List.length_pos_iff.mpr hL)⟩ : Fin L.length) =
            ⟨n, hnlt⟩ := by
        apply Fin.ext
        exact Nat.mod_eq_of_lt hnlt
      rw [hindex, List.get_eq_getElem,
        List.take_succ_eq_append_getElem hnlt, sweepFrom_append]
      rfl

/-- One complete round-robin cycle is definitionally the existing
source-order sweep. -/
theorem roundRobinScheduler_cycle_eq_sweepFrom
    (L : List Node) (hL : L ≠ []) (s : State) :
    attemptRun L L.length (roundRobinScheduler L hL) s =
      sweepFrom s L := by
  simpa using
    roundRobinScheduler_prefix_eq_sweepFrom L hL s L.length le_rfl

theorem roundRobinScheduler_memberRecurrent
    (L : List Node) (hL : L ≠ []) (s : State) :
    NodeMemberRecurrent L (roundRobinScheduler L hL) s := by
  intro member N
  obtain ⟨i, hi⟩ := List.mem_iff_get.mp member.2
  let m := N * L.length + i.val
  have hlen : 1 ≤ L.length := List.length_pos_iff.mpr hL
  have hNm : N ≤ m := by
    calc
      N = N * 1 := by simp
      _ ≤ N * L.length := Nat.mul_le_mul_left N hlen
      _ ≤ N * L.length + i.val := Nat.le_add_right _ _
  refine ⟨m, hNm, ?_⟩
  change L.get
      ⟨m % L.length,
        Nat.mod_lt m (List.length_pos_iff.mpr hL)⟩ = member.1
  have hindex :
      (⟨m % L.length,
        Nat.mod_lt m (List.length_pos_iff.mpr hL)⟩ : Fin L.length) = i := by
    apply Fin.ext
    simp [m, Nat.add_mod, Nat.mod_eq_of_lt i.isLt]
  rw [hindex, hi]

/-- From every starting index, round robin selects any declared member within
one finite cycle. -/
theorem roundRobinScheduler_selects_within
    (L : List Node) (hL : L ≠ [])
    (member : {n : Node // n ∈ L}) (start : Nat) (s : State) :
    ∃ offset, offset < L.length ∧
      (roundRobinScheduler L hL (start + offset) s).1 = member.1 := by
  obtain ⟨i, hi⟩ := List.mem_iff_get.mp member.2
  have hlen : 0 < L.length := List.length_pos_iff.mpr hL
  let a := start % L.length
  let offset := if a ≤ i.val then i.val - a else L.length - a + i.val
  have ha : a < L.length := Nat.mod_lt start hlen
  have hoffset : offset < L.length := by
    dsimp [offset]
    split <;> omega
  have hmod : (start + offset) % L.length = i.val := by
    rw [Nat.add_mod]
    rw [Nat.mod_eq_of_lt hoffset]
    dsimp [offset, a]
    split
    · rename_i h
      have hadd : start % L.length + (i.val - start % L.length) = i.val :=
        Nat.add_sub_of_le h
      rw [hadd, Nat.mod_eq_of_lt i.isLt]
    · rename_i h
      have hle : start % L.length ≤ L.length := Nat.le_of_lt ha
      have hadd : start % L.length + (L.length - start % L.length) = L.length :=
        Nat.add_sub_of_le hle
      rw [← Nat.add_assoc, hadd, Nat.add_mod_left, Nat.mod_eq_of_lt i.isLt]
  refine ⟨offset, hoffset, ?_⟩
  change L.get
      ⟨(start + offset) % L.length,
        Nat.mod_lt (start + offset) hlen⟩ = member.1
  have hindex :
      (⟨(start + offset) % L.length,
        Nat.mod_lt (start + offset) hlen⟩ : Fin L.length) = i := by
    apply Fin.ext
    exact hmod
  rw [hindex, hi]

def fixedRoundRobinScheduler {k : Nat} (phi : Formula k) :
    NodeScheduler (fixedProgram phi) :=
  roundRobinScheduler (fixedProgram phi) (fixedProgram_nonempty phi)

theorem fixedRoundRobin_cycle_eq_fixedSweep
    {k : Nat} (phi : Formula k) (s : State) :
    attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s =
      sweepFrom s (fixedProgram phi) :=
  roundRobinScheduler_cycle_eq_sweepFrom
    (fixedProgram phi) (fixedProgram_nonempty phi) s

/-- Source-order round robin reaches consensus after at most one cycle, hence
within a number of attempts linear in the emitted node count. -/
theorem fixedRoundRobin_consensus_after_one_cycle
    {k : Nat} (phi : Formula k) (s : State) :
    Consensus (fixedFederation phi)
      (attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s) := by
  rw [fixedRoundRobin_cycle_eq_fixedSweep]
  exact fixedSweep_consensus phi s

theorem fixedRoundRobin_output_after_one_cycle
    {k : Nat} (phi : Formula k) (x : Fin k → Bool) (s : State)
    (hinput : CarriesInput x s) :
    attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s (fixedOutReg phi) =
      Formula.evalF phi x := by
  rw [fixedRoundRobin_cycle_eq_fixedSweep]
  exact fixedSweep_output phi x s hinput

theorem fixedRoundRobin_stable_after_one_cycle
    {k : Nat} (phi : Formula k) (s : State) :
    ∀ n, (fixedProgram phi).length ≤ n →
      attemptRun (fixedProgram phi) n (fixedRoundRobinScheduler phi) s =
        attemptRun (fixedProgram phi) (fixedProgram phi).length
          (fixedRoundRobinScheduler phi) s := by
  intro n hn
  obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le hn
  rw [attemptRun_add]
  exact attemptRun_eq_self_of_consensus
    (fixedProgram phi)
    (fun k x => fixedRoundRobinScheduler phi ((fixedProgram phi).length + k) x)
    (attemptRun (fixedProgram phi) (fixedProgram phi).length
      (fixedRoundRobinScheduler phi) s)
    (fixedRoundRobin_consensus_after_one_cycle phi s) m

theorem fixedRoundRobin_memberRecurrent {k : Nat} (phi : Formula k)
    (s : State) :
    NodeMemberRecurrent (fixedProgram phi)
      (fixedRoundRobinScheduler phi) s :=
  roundRobinScheduler_memberRecurrent
    (fixedProgram phi) (fixedProgram_nonempty phi) s

theorem fixedRoundRobin_siteRecurrent {k : Nat} (phi : Formula k)
    (s : State) :
    NodeSiteRecurrent (fixedProgram phi)
      (fixedRoundRobinScheduler phi) s :=
  nodeMemberRecurrent_nodeSiteRecurrent
    (fixedRoundRobin_memberRecurrent phi s)

theorem fixedRoundRobin_pathwiseWeakFair {k : Nat} (phi : Formula k)
    (s : State) :
    NodePathwiseWeakFair (fixedProgram phi)
      (fixedRoundRobinScheduler phi) s :=
  nodeMemberRecurrent_nodePathwiseWeakFair
    (fixedRoundRobin_memberRecurrent phi s)

/-- At least one state-blind mathematical scheduler is recurrent and computes
the fixed formula correctly.  This is not a physical realization theorem. -/
theorem exists_mathematical_fair_scheduler {k : Nat}
    (phi : Formula k) (x : Fin k → Bool) (s : State)
    (hinput : CarriesInput x s) :
    ∃ sigma : NodeScheduler (fixedProgram phi),
      NodeMemberRecurrent (fixedProgram phi) sigma s ∧
      NodeSiteRecurrent (fixedProgram phi) sigma s ∧
      NodePathwiseWeakFair (fixedProgram phi) sigma s ∧
      ∃ N,
        Consensus (fixedFederation phi)
          (attemptRun (fixedProgram phi) N sigma s) ∧
        (∀ n, N ≤ n →
          attemptRun (fixedProgram phi) n sigma s =
            attemptRun (fixedProgram phi) N sigma s) ∧
        CarriesInput x (attemptRun (fixedProgram phi) N sigma s) ∧
        attemptRun (fixedProgram phi) N sigma s (fixedOutReg phi) =
          Formula.evalF phi x := by
  let sigma := fixedRoundRobinScheduler phi
  have hmember : NodeMemberRecurrent (fixedProgram phi) sigma s :=
    fixedRoundRobin_memberRecurrent phi s
  have hsite : NodeSiteRecurrent (fixedProgram phi) sigma s :=
    nodeMemberRecurrent_nodeSiteRecurrent hmember
  have hfair : NodePathwiseWeakFair (fixedProgram phi) sigma s :=
    nodeMemberRecurrent_nodePathwiseWeakFair hmember
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
  refine ⟨sigma, hmember, hsite, hfair, N, hcons, hconstant,
    hrunInput, ?_⟩
  exact houtput.trans (fixedSweep_output phi x s hinput)

#print axioms stableConsensusTail_nodePathwiseWeakFair
#print axioms nodePathwiseWeakFair_iff_eventuallyStableConsensus
#print axioms nodeMemberRecurrent_iff_nodeSiteRecurrent
#print axioms fixedRoundRobin_memberRecurrent
#print axioms roundRobinScheduler_cycle_eq_sweepFrom
#print axioms fixedRoundRobin_consensus_after_one_cycle
#print axioms exists_mathematical_fair_scheduler

end


end OPH.RepairUniversality.FixedFederation
