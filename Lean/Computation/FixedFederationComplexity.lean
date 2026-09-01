import Computation.FixedFederationExecution
import ObserverPatchHolography.Execution.RankedAttemptCapacity

/-!
# Fixed-program accepted-step complexity

The formula compiler emits a tree-shaped straight-line program.  Generated
registers have at most one downstream consumer, so linear weights replace the
generic exponential defect rank.  The resulting bound counts genuine
canonical accepted steps, not scheduler attempts.
-/

namespace OPH.RepairUniversality.FixedFederation

open OPH.RepairUniversality

noncomputable section

def dependencyOccurrences (L : List Node) : List Nat :=
  L.flatMap fun n => n.deps

@[simp] theorem dependencyOccurrences_nil :
    dependencyOccurrences [] = [] := rfl

@[simp] theorem dependencyOccurrences_cons (n : Node) (L : List Node) :
    dependencyOccurrences (n :: L) =
      n.deps ++ dependencyOccurrences L := rfl

@[simp] theorem dependencyOccurrences_append (A B : List Node) :
    dependencyOccurrences (A ++ B) =
      dependencyOccurrences A ++ dependencyOccurrences B := by
  simp [dependencyOccurrences]

theorem go_dependency_scope {k : Nat}
    (phi : Formula k) (c : Nat) (hc : k ≤ c)
    {r : Nat}
    (hr : r ∈ dependencyOccurrences (go phi c).nodes) :
    r < k ∨ (c ≤ r ∧ r < (go phi c).next) := by
  rcases List.mem_flatMap.mp hr with ⟨n, hn, hrdep⟩
  have F := go_facts phi c hc (fun _ => false)
  rcases F.deps_cov n hn r hrdep with hexternal | hinternal
  · exact Or.inl hexternal
  · obtain ⟨m, hm, hmr⟩ := List.mem_map.mp hinternal
    right
    simpa [hmr] using F.range m hm

theorem go_dependency_count_eq_zero_of_lt_start {k : Nat}
    (phi : Formula k) (c : Nat) (hc : k ≤ c)
    (r : Nat) (hkr : k ≤ r) (hrc : r < c) :
    (dependencyOccurrences (go phi c).nodes).count r = 0 := by
  apply List.count_eq_zero.mpr
  intro hr
  rcases go_dependency_scope phi c hc hr with hext | hint
  · omega
  · omega

theorem go_dependency_count_eq_zero_of_next_le {k : Nat}
    (phi : Formula k) (c : Nat) (hc : k ≤ c)
    (r : Nat) (hr : (go phi c).next ≤ r) :
    (dependencyOccurrences (go phi c).nodes).count r = 0 := by
  apply List.count_eq_zero.mpr
  intro hmem
  have F := go_facts phi c hc (fun _ => false)
  have hcnext : c < (go phi c).next :=
    lt_of_le_of_lt F.c_le_out F.out_lt_next
  rcases go_dependency_scope phi c hc hmem with hext | hint
  · omega
  · omega

structure GoUseProfile {k : Nat} (phi : Formula k) (c : Nat) : Prop where
  output_unused :
    (dependencyOccurrences (go phi c).nodes).count (go phi c).out = 0
  nonoutput_used_once :
    ∀ r ∈ (go phi c).nodes.map Node.reg,
      r ≠ (go phi c).out →
      (dependencyOccurrences (go phi c).nodes).count r = 1

theorem go_use_profile {k : Nat}
    (phi : Formula k) (c : Nat) (hc : k ≤ c) :
    GoUseProfile phi c := by
  induction phi generalizing c with
  | var i =>
      constructor
      · have hic : i.val ≠ c :=
          Nat.ne_of_lt (lt_of_lt_of_le i.isLt hc)
        simp [go, dependencyOccurrences, hic]
      · intro r hr hne
        simp [go] at hr
        subst r
        exact False.elim (hne rfl)
  | const b =>
      constructor
      · simp [go, dependencyOccurrences]
      · intro r hr hne
        simp [go] at hr
        subst r
        exact False.elim (hne rfl)
  | nand phi psi ihPhi ihPsi =>
      let A := go phi c
      let B := go psi A.next
      have Fphi := go_facts phi c hc (fun _ => false)
      have hcB : k ≤ A.next := by
        dsimp [A]
        exact le_trans hc
          (le_trans Fphi.c_le_out (Nat.le_of_lt Fphi.out_lt_next))
      have Fpsi := go_facts psi A.next hcB (fun _ => false)
      have PA := ihPhi c hc
      have PB := ihPsi A.next hcB
      have hAout_lt_Bout : A.out < B.out := by
        dsimp [A, B]
        exact lt_of_lt_of_le Fphi.out_lt_next Fpsi.c_le_out
      have hBout_lt_parent : B.out < B.next := Fpsi.out_lt_next
      have hAout_lt_parent : A.out < B.next :=
        lt_trans hAout_lt_Bout hBout_lt_parent
      have hAparent :
          (dependencyOccurrences A.nodes).count B.next = 0 := by
        dsimp [A, B]
        exact go_dependency_count_eq_zero_of_next_le phi c hc _
          (le_trans Fpsi.c_le_out (Nat.le_of_lt Fpsi.out_lt_next))
      have hBparent :
          (dependencyOccurrences B.nodes).count B.next = 0 := by
        dsimp [A, B]
        exact go_dependency_count_eq_zero_of_next_le
          psi (go phi c).next hcB _ le_rfl
      constructor
      · simp only [go, dependencyOccurrences_append,
          dependencyOccurrences_cons, dependencyOccurrences_nil,
          List.append_nil, List.count_append]
        simp [A, B, hAparent, hBparent,
          ne_of_lt hAout_lt_parent,
          ne_of_lt hBout_lt_parent]
      · intro r hrmem hrne
        simp only [go, List.map_append, List.mem_append] at hrmem
        rcases hrmem with hrA | hrB | hrParent
        · obtain ⟨a, ha, har⟩ := List.mem_map.mp hrA
          subst r
          have hcr : c ≤ a.reg := by
            exact (Fphi.range a ha).1
          have hrA_next : a.reg < A.next := by
            exact (Fphi.range a ha).2
          have hBzero :
              (dependencyOccurrences B.nodes).count a.reg = 0 := by
            dsimp [B]
            exact go_dependency_count_eq_zero_of_lt_start
              psi A.next hcB a.reg (le_trans hc hcr) hrA_next
          have hr_ne_Bout : a.reg ≠ B.out :=
            ne_of_lt (lt_of_lt_of_le hrA_next Fpsi.c_le_out)
          by_cases hrAout : a.reg = A.out
          · have hAzeroA :
                (dependencyOccurrences (go phi c).nodes).count a.reg = 0 := by
              simpa [hrAout] using PA.output_unused
            have hBzeroA :
                (dependencyOccurrences (go psi (go phi c).next).nodes).count
                  a.reg = 0 := by
              exact hBzero
            simp only [go, dependencyOccurrences_append,
              dependencyOccurrences_cons, dependencyOccurrences_nil,
              List.append_nil, List.count_append]
            rw [hAzeroA, hBzeroA]
            have heqA : (go phi c).out = a.reg := by
              simpa [A] using hrAout.symm
            have hneB : (go psi (go phi c).next).out ≠ a.reg :=
              Ne.symm (by simpa [A, B] using hr_ne_Bout)
            simp [heqA, hneB]
          · have hAone := PA.nonoutput_used_once a.reg hrA hrAout
            simp only [go, dependencyOccurrences_append,
              dependencyOccurrences_cons, dependencyOccurrences_nil,
              List.append_nil, List.count_append]
            rw [hAone, hBzero]
            have hneA : a.reg ≠ (go phi c).out := by
              simpa [A] using hrAout
            have hneB : a.reg ≠ (go psi (go phi c).next).out := by
              simpa [A, B] using hr_ne_Bout
            simp [Ne.symm hneA, Ne.symm hneB]
        · obtain ⟨b, hb, hbr⟩ := List.mem_map.mp hrB
          subst r
          have hA_next_r : A.next ≤ b.reg := by
            exact (Fpsi.range b hb).1
          have hAzero :
              (dependencyOccurrences A.nodes).count b.reg = 0 := by
            dsimp [A]
            exact go_dependency_count_eq_zero_of_next_le
              phi c hc b.reg hA_next_r
          have hr_ne_Aout : b.reg ≠ A.out :=
            ne_of_gt (lt_of_lt_of_le Fphi.out_lt_next hA_next_r)
          by_cases hrBout : b.reg = B.out
          · have hBzeroB :
                (dependencyOccurrences (go psi (go phi c).next).nodes).count
                  b.reg = 0 := by
              simpa [hrBout] using PB.output_unused
            have hAzeroB :
                (dependencyOccurrences (go phi c).nodes).count b.reg = 0 := by
              exact hAzero
            simp only [go, dependencyOccurrences_append,
              dependencyOccurrences_cons, dependencyOccurrences_nil,
              List.append_nil, List.count_append]
            rw [hAzeroB, hBzeroB]
            have heqB : (go psi (go phi c).next).out = b.reg := by
              simpa [A, B] using hrBout.symm
            have hneA : (go phi c).out ≠ b.reg :=
              Ne.symm (by simpa [A] using hr_ne_Aout)
            simp [heqB, hneA]
          · have hBone := PB.nonoutput_used_once b.reg hrB hrBout
            simp only [go, dependencyOccurrences_append,
              dependencyOccurrences_cons, dependencyOccurrences_nil,
              List.append_nil, List.count_append]
            rw [hAzero, hBone]
            have hneA : b.reg ≠ (go phi c).out := by
              simpa [A] using hr_ne_Aout
            have hneB : b.reg ≠ (go psi (go phi c).next).out := by
              simpa [A, B] using hrBout
            simp [Ne.symm hneA, Ne.symm hneB]
        · apply False.elim
          apply hrne
          simpa [go] using hrParent

def AtMostOneDownstreamConsumer (L : List Node) : Prop :=
  ∀ n ∈ L, (dependencyOccurrences L).count n.reg ≤ 1

theorem fixedProgram_atMostOneDownstreamConsumer {k : Nat}
    (phi : Formula k) :
    AtMostOneDownstreamConsumer (fixedProgram phi) := by
  intro n hn
  have P := go_use_profile phi k le_rfl
  have hreg : n.reg ∈ (go phi k).nodes.map Node.reg :=
    List.mem_map.mpr ⟨n, hn, rfl⟩
  change (dependencyOccurrences (go phi k).nodes).count n.reg ≤ 1
  by_cases hout : n.reg = (go phi k).out
  · rw [hout, P.output_unused]
    omega
  · rw [P.nonoutput_used_once n.reg hreg hout]

theorem observer_ok_write_eq_of_not_used
    (a : Node) (s : State) (r : Nat) (b : Bool)
    (hreg : a.reg ≠ r) (hdep : r ∉ a.deps) :
    a.obs.ok (write s r b) = a.obs.ok s := by
  have hr : write s r b a.reg = s a.reg := by
    simp [write, hreg]
  have hv : a.val (write s r b) = a.val s := by
    apply a.val_local
    intro i hi
    simp only [write]
    rw [if_neg]
    intro hir
    apply hdep
    exact hir ▸ hi
  simp [Node.obs, hr, hv]

def linearDefectRank : List Node → State → Nat
  | [], _ => 0
  | n :: L, s =>
      (if n.obs.ok s = false then L.length + 1 else 0) +
        linearDefectRank L s

theorem linearDefectRank_write_eq_of_not_used
    (L : List Node) (s : State) (r : Nat) (b : Bool)
    (hregs : ∀ n ∈ L, n.reg ≠ r)
    (hnot : r ∉ dependencyOccurrences L) :
    linearDefectRank L (write s r b) = linearDefectRank L s := by
  induction L with
  | nil => rfl
  | cons a L ih =>
      have hhead : r ∉ a.deps := by
        intro hr
        apply hnot
        simp only [dependencyOccurrences_cons, List.mem_append]
        exact Or.inl hr
      have htail : r ∉ dependencyOccurrences L := by
        intro hr
        apply hnot
        simp only [dependencyOccurrences_cons, List.mem_append]
        exact Or.inr hr
      have hobs := observer_ok_write_eq_of_not_used
        a s r b (hregs a List.mem_cons_self) hhead
      have hregsTail : ∀ n ∈ L, n.reg ≠ r :=
        fun n hn => hregs n (List.mem_cons_of_mem a hn)
      simp [linearDefectRank, hobs, ih hregsTail htail]

theorem linearDefectRank_write_le_add_length
    (L : List Node) (s : State) (r : Nat) (b : Bool)
    (hregs : ∀ n ∈ L, n.reg ≠ r)
    (hcount : (dependencyOccurrences L).count r ≤ 1) :
    linearDefectRank L (write s r b) ≤
      linearDefectRank L s + L.length := by
  induction L with
  | nil => simp [linearDefectRank]
  | cons a L ih =>
      have hregA : a.reg ≠ r := hregs a List.mem_cons_self
      have hregsTail : ∀ n ∈ L, n.reg ≠ r :=
        fun n hn => hregs n (List.mem_cons_of_mem a hn)
      have hsplit :
          a.deps.count r + (dependencyOccurrences L).count r ≤ 1 := by
        simpa [dependencyOccurrences, List.count_append] using hcount
      by_cases hused : r ∈ a.deps
      · have hpos : 0 < a.deps.count r := List.count_pos_iff.mpr hused
        have htailZero : (dependencyOccurrences L).count r = 0 := by
          omega
        have htailNot : r ∉ dependencyOccurrences L :=
          List.count_eq_zero.mp htailZero
        have htailEq := linearDefectRank_write_eq_of_not_used
          L s r b hregsTail htailNot
        simp only [linearDefectRank, htailEq, List.length_cons]
        split <;> split <;> omega
      · have hobs := observer_ok_write_eq_of_not_used
          a s r b hregA hused
        have htailCount : (dependencyOccurrences L).count r ≤ 1 := by
          omega
        have htailLe := ih hregsTail htailCount
        simp only [linearDefectRank, hobs, List.length_cons]
        omega

theorem atMostOneDownstreamConsumer_tail
    (a : Node) (L : List Node)
    (h : AtMostOneDownstreamConsumer (a :: L)) :
    AtMostOneDownstreamConsumer L := by
  intro n hn
  have hfull := h n (List.mem_cons_of_mem a hn)
  have hsplit :
      a.deps.count n.reg +
        (dependencyOccurrences L).count n.reg ≤ 1 := by
    simpa [dependencyOccurrences, List.count_append] using hfull
  omega

theorem canonicalAcceptedStep_linearDefectRank_lt
    (L : List Node) (hWF : NodesWF L)
    (hsingle : AtMostOneDownstreamConsumer L)
    {s t : State} (hstep : CanonicalAcceptedStep L s t) :
    linearDefectRank L t < linearDefectRank L s := by
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
        have hregsTail : ∀ n ∈ L, n.reg ≠ a.reg := by
          intro n hn
          exact ne_of_gt (hpw.1 n hn)
        have hfull := hsingle a List.mem_cons_self
        have htailCount :
            (dependencyOccurrences L).count a.reg ≤ 1 := by
          have hsplit :
              a.deps.count a.reg +
                (dependencyOccurrences L).count a.reg ≤ 1 := by
            simpa [dependencyOccurrences, List.count_append] using hfull
          omega
        have htailLe := linearDefectRank_write_le_add_length
          L s a.reg (a.val s) hregsTail htailCount
        simp only [repairNode] at hpost htailLe ⊢
        have htailStrict :
            linearDefectRank L (write s a.reg (a.val s)) <
              L.length + 1 + linearDefectRank L s := by
          omega
        simpa [linearDefectRank, hfail, hpost] using htailStrict
      · have hhead := earlier_observer_unchanged a n s
          (hpw.1 n hn) (hdepsAll a List.mem_cons_self)
        have htailWF : NodesWF L :=
          ⟨hpw.2, fun p hp =>
            hdepsAll p (List.mem_cons_of_mem a hp)⟩
        have htailSingle :=
          atMostOneDownstreamConsumer_tail a L hsingle
        have htailStep :
            CanonicalAcceptedStep L s (repairNode n s) :=
          ⟨n, hn, hfail, rfl⟩
        have htail := ih htailWF htailSingle htailStep
        simp only [linearDefectRank, hhead]
        exact Nat.add_lt_add_left htail _

inductive CanonicalAcceptedSteps (L : List Node) :
    Nat → State → State → Prop
  | refl (s : State) : CanonicalAcceptedSteps L 0 s s
  | tail {m : Nat} {s t u : State} :
      CanonicalAcceptedSteps L m s t →
      CanonicalAcceptedStep L t u →
      CanonicalAcceptedSteps L (m + 1) s u

theorem canonicalAcceptedSteps_rank_budget
    (L : List Node) (hWF : NodesWF L)
    (hsingle : AtMostOneDownstreamConsumer L)
    {m : Nat} {s t : State}
    (h : CanonicalAcceptedSteps L m s t) :
    m + linearDefectRank L t ≤ linearDefectRank L s := by
  induction h with
  | refl s => simp
  | tail hpath hstep ih =>
      have hlt := canonicalAcceptedStep_linearDefectRank_lt
        L hWF hsingle hstep
      omega

def triangle : Nat → Nat
  | 0 => 0
  | n + 1 => triangle n + (n + 1)

theorem linearDefectRank_le_triangle (L : List Node) (s : State) :
    linearDefectRank L s ≤ triangle L.length := by
  induction L with
  | nil => simp [linearDefectRank, triangle]
  | cons a L ih =>
      by_cases hfail : a.obs.ok s = false
      · simp [linearDefectRank, triangle, hfail]
        omega
      · simp [linearDefectRank, triangle, hfail]
        omega

theorem triangle_double (n : Nat) :
    2 * triangle n = n * (n + 1) := by
  induction n with
  | zero => simp [triangle]
  | succ n ih =>
      calc
        2 * triangle (n + 1) =
            2 * (triangle n + (n + 1)) := by rfl
        _ = 2 * triangle n + 2 * (n + 1) := by omega
        _ = n * (n + 1) + 2 * (n + 1) := by rw [ih]
        _ = (n + 1) * ((n + 1) + 1) := by ring

theorem triangle_eq_closed (n : Nat) :
    triangle n = n * (n + 1) / 2 := by
  have h := congrArg (fun z : Nat => z / 2) (triangle_double n)
  simpa [Nat.mul_comm] using h

theorem fixedProgram_acceptedSteps_triangle
    {k : Nat} (phi : Formula k)
    {m : Nat} {s t : State}
    (h : CanonicalAcceptedSteps (fixedProgram phi) m s t) :
    m ≤ triangle (fixedProgram phi).length := by
  have hbudget := canonicalAcceptedSteps_rank_budget
    (fixedProgram phi) (fixedProgram_wf phi)
    (fixedProgram_atMostOneDownstreamConsumer phi) h
  have hcap := linearDefectRank_le_triangle (fixedProgram phi) s
  omega

theorem fixedProgram_acceptedSteps_quadratic
    {k : Nat} (phi : Formula k)
    {m : Nat} {s t : State}
    (h : CanonicalAcceptedSteps (fixedProgram phi) m s t) :
    m ≤ (fixedProgram phi).length *
      ((fixedProgram phi).length + 1) / 2 := by
  rw [← triangle_eq_closed]
  exact fixedProgram_acceptedSteps_triangle phi h

/-! ## Conditional attempt bounds through the shared ranked provider -/

def fixedProgramPotential {k : Nat} (phi : Formula k) (s : State) : Nat :=
  linearDefectRank (fixedProgram phi) s

def fixedProgramRankedAttemptSystem {k : Nat} (phi : Formula k) :
    OPH.RankedAttempt.System
      State {n : Node // n ∈ fixedProgram phi} where
  step := fun member s => repairNode member.1 s
  rank := fixedProgramPotential phi
  change_rank_lt := by
    intro member s hchange
    have hfail : member.1.obs.ok s = false := by
      cases hok : member.1.obs.ok s with
      | false => rfl
      | true =>
          exact False.elim (hchange
            (repairNode_eq_self_of_accepts member.1 s hok))
    exact canonicalAcceptedStep_linearDefectRank_lt
      (fixedProgram phi) (fixedProgram_wf phi)
      (fixedProgram_atMostOneDownstreamConsumer phi)
      ⟨member.1, member.2, hfail, rfl⟩

@[simp] theorem rankedAttemptRun_eq_fixedAttemptRun
    {k : Nat} (phi : Formula k) (n : Nat)
    (sigma : NodeScheduler (fixedProgram phi)) (s : State) :
    OPH.RankedAttempt.run
        (fixedProgramRankedAttemptSystem phi) n sigma s =
      attemptRun (fixedProgram phi) n sigma s := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [OPH.RankedAttempt.run_last_step, attemptRun_last_step, ih]
      rfl

theorem fixedProgram_rankedQuiescent_iff_consensus
    {k : Nat} (phi : Formula k) (s : State) :
    OPH.RankedAttempt.Quiescent
        (fixedProgramRankedAttemptSystem phi) s ↔
      Consensus (fixedFederation phi) s := by
  constructor
  · intro hquiet o ho
    rcases List.mem_map.mp ho with ⟨n, hn, rfl⟩
    by_contra hnot
    have hfail : n.obs.ok s = false := by simpa using hnot
    let member : {n : Node // n ∈ fixedProgram phi} := ⟨n, hn⟩
    exact repairNode_ne_self_of_fails n s hfail (hquiet member)
  · intro hcons member
    apply repairNode_eq_self_of_accepts
    exact hcons member.1.obs
      (List.mem_map.mpr ⟨member.1, member.2, rfl⟩)

def NodeBoundedWaste (L : List Node) (waste : Nat)
    (sigma : NodeScheduler L) : Prop :=
  ∀ start s, ¬ Consensus (L.map Node.obs) s →
    ∃ offset, offset ≤ waste ∧
      attemptRun L (offset + 1) (fun k x => sigma (start + k) x) s ≠
        attemptRun L offset (fun k x => sigma (start + k) x) s

theorem fixedProgram_nodeBoundedWaste_iff_ranked
    {k : Nat} (phi : Formula k) (waste : Nat)
    (sigma : NodeScheduler (fixedProgram phi)) :
    NodeBoundedWaste (fixedProgram phi) waste sigma ↔
      OPH.RankedAttempt.BoundedWaste
        (fixedProgramRankedAttemptSystem phi) waste sigma := by
  constructor
  · intro h start s hnotQuiet
    have hnotConsensus : ¬ Consensus (fixedFederation phi) s := by
      intro hcons
      exact hnotQuiet ((fixedProgram_rankedQuiescent_iff_consensus phi s).2 hcons)
    obtain ⟨offset, hoffset, hchange⟩ := h start s hnotConsensus
    refine ⟨offset, hoffset, ?_⟩
    simpa only [rankedAttemptRun_eq_fixedAttemptRun] using hchange
  · intro h start s hnotConsensus
    have hnotQuiet :
        ¬ OPH.RankedAttempt.Quiescent
          (fixedProgramRankedAttemptSystem phi) s := by
      intro hquiet
      exact hnotConsensus
        ((fixedProgram_rankedQuiescent_iff_consensus phi s).1 hquiet)
    obtain ⟨offset, hoffset, hchange⟩ := h start s hnotQuiet
    refine ⟨offset, hoffset, ?_⟩
    simpa only [rankedAttemptRun_eq_fixedAttemptRun] using hchange

theorem fixedProgram_nodeBoundedWaste_eventually_consensus
    {k : Nat} (phi : Formula k) (waste : Nat)
    (s : State) (sigma : NodeScheduler (fixedProgram phi))
    (hbounded : NodeBoundedWaste (fixedProgram phi) waste sigma) :
    ∃ N,
      N ≤ (waste + 1) * fixedProgramPotential phi s ∧
      Consensus (fixedFederation phi)
        (attemptRun (fixedProgram phi) N sigma s) ∧
      ∀ n, N ≤ n →
        attemptRun (fixedProgram phi) n sigma s =
          attemptRun (fixedProgram phi) N sigma s := by
  have hranked :=
    (fixedProgram_nodeBoundedWaste_iff_ranked phi waste sigma).mp hbounded
  obtain ⟨N, hN, hquiet, hconstant⟩ :=
    OPH.RankedAttempt.boundedWaste_eventually_quiescent
      (fixedProgramRankedAttemptSystem phi) waste s sigma hranked
  refine ⟨N, hN, ?_, ?_⟩
  · exact (fixedProgram_rankedQuiescent_iff_consensus phi _).1
      (by simpa only [rankedAttemptRun_eq_fixedAttemptRun] using hquiet)
  · intro n hn
    simpa only [rankedAttemptRun_eq_fixedAttemptRun] using hconstant n hn

theorem fixedProgram_nodeBoundedWaste_triangle
    {k : Nat} (phi : Formula k) (waste : Nat)
    (s : State) (sigma : NodeScheduler (fixedProgram phi))
    (hbounded : NodeBoundedWaste (fixedProgram phi) waste sigma) :
    ∃ N,
      N ≤ (waste + 1) * triangle (fixedProgram phi).length ∧
      Consensus (fixedFederation phi)
        (attemptRun (fixedProgram phi) N sigma s) ∧
      ∀ n, N ≤ n →
        attemptRun (fixedProgram phi) n sigma s =
          attemptRun (fixedProgram phi) N sigma s := by
  obtain ⟨N, hN, hcons, hconstant⟩ :=
    fixedProgram_nodeBoundedWaste_eventually_consensus
      phi waste s sigma hbounded
  have hpotential := linearDefectRank_le_triangle (fixedProgram phi) s
  refine ⟨N, ?_, hcons, hconstant⟩
  exact le_trans hN (Nat.mul_le_mul_left (waste + 1) hpotential)

#print axioms fixedProgram_atMostOneDownstreamConsumer
#print axioms canonicalAcceptedStep_linearDefectRank_lt
#print axioms fixedProgram_acceptedSteps_quadratic
#print axioms fixedProgram_nodeBoundedWaste_triangle

end

end OPH.RepairUniversality.FixedFederation
