import Computation.FixedFederationExecutionExamples
import ObserverPatchHolography.Execution.RankedAttemptCapacity

/-!
# Fixed-federation execution endpoints

The Tower boundary consumes the exact fairness classification, explicit
mathematical scheduler, sharp accepted-step infrastructure, and the shared
bounded-waste provider.  Attempt bounds remain conditional and carry no
physical unit.
-/

namespace OPH.Tower.FixedFederationExecutionEndpoint

open OPH.RepairUniversality
open OPH.RepairUniversality.FixedFederation

noncomputable section

/-- The explicit state-blind scheduler reaches a stable, correct endpoint.
This is an existence result for a mathematical function only. -/
theorem roundRobin_reaches_correct_output {k : Nat}
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
          Formula.evalF phi x :=
  exists_mathematical_fair_scheduler phi x s hinput

/-- The emission-order round-robin scheduler reaches stable consensus and the
exact formula output after one emitted-node cycle.  This is a linear attempt
bound for this explicit mathematical scheduler, not a physical runtime. -/
theorem roundRobin_reaches_correct_output_within_linear {k : Nat}
    (phi : Formula k) (x : Fin k → Bool) (s : State)
    (hinput : CarriesInput x s) :
    Consensus (fixedFederation phi)
      (attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s) ∧
    (∀ n, (fixedProgram phi).length ≤ n →
      attemptRun (fixedProgram phi) n (fixedRoundRobinScheduler phi) s =
        attemptRun (fixedProgram phi) (fixedProgram phi).length
          (fixedRoundRobinScheduler phi) s) ∧
    CarriesInput x
      (attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s) ∧
    attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s (fixedOutReg phi) =
      Formula.evalF phi x := by
  refine ⟨fixedRoundRobin_consensus_after_one_cycle phi s,
    fixedRoundRobin_stable_after_one_cycle phi s, ?_,
    fixedRoundRobin_output_after_one_cycle phi x s hinput⟩
  exact (fixedAttemptRun_preservesInputObservation phi
    (fixedProgram phi).length (fixedRoundRobinScheduler phi) s).trans hinput

/-- The actual compiler accepted-step upper is consumed at the Tower boundary.
This wrapper receives no separate novelty credit. -/
theorem canonicalAcceptedSteps_within_quadratic {k : Nat}
    (phi : Formula k) {m : Nat} {s t : State}
    (hsteps : CanonicalAcceptedSteps (fixedProgram phi) m s t) :
    m ≤ (fixedProgram phi).length *
      ((fixedProgram phi).length + 1) / 2 :=
  fixedProgram_acceptedSteps_quadratic phi hsteps

/-- A supplied fixed-node bounded-waste premise yields a triangular-size
attempt horizon and a stable same-input unique output. -/
theorem boundedWaste_reaches_unique_output {k : Nat}
    (phi : Formula k) (x : Fin k → Bool) (s : State)
    (sigma : NodeScheduler (fixedProgram phi)) (waste : Nat)
    (hinput : CarriesInput x s)
    (hbounded : NodeBoundedWaste (fixedProgram phi) waste sigma) :
    ∃ N,
      N ≤ (waste + 1) * triangle (fixedProgram phi).length ∧
      Consensus (fixedFederation phi)
        (attemptRun (fixedProgram phi) N sigma s) ∧
      (∀ n, N ≤ n →
        attemptRun (fixedProgram phi) n sigma s =
          attemptRun (fixedProgram phi) N sigma s) ∧
      CarriesInput x (attemptRun (fixedProgram phi) N sigma s) ∧
      ∀ t : State,
        Consensus (fixedFederation phi) t →
        CarriesInput x t →
        attemptRun (fixedProgram phi) N sigma s (fixedOutReg phi) =
          t (fixedOutReg phi) := by
  have hranked :=
    (fixedProgram_nodeBoundedWaste_iff_ranked phi waste sigma).mp hbounded
  obtain ⟨N, hN, hquiet, hconstantRanked⟩ :=
    OPH.RankedAttempt.boundedWaste_eventually_quiescent
      (fixedProgramRankedAttemptSystem phi) waste s sigma hranked
  have hcons : Consensus (fixedFederation phi)
      (attemptRun (fixedProgram phi) N sigma s) :=
    (fixedProgram_rankedQuiescent_iff_consensus phi _).1
      (by simpa only [rankedAttemptRun_eq_fixedAttemptRun] using hquiet)
  have hconstant : ∀ n, N ≤ n →
      attemptRun (fixedProgram phi) n sigma s =
        attemptRun (fixedProgram phi) N sigma s := by
    intro n hn
    simpa only [rankedAttemptRun_eq_fixedAttemptRun] using
      hconstantRanked n hn
  have hpotential := linearDefectRank_le_triangle (fixedProgram phi) s
  have hNtriangle :
      N ≤ (waste + 1) * triangle (fixedProgram phi).length :=
    le_trans hN (Nat.mul_le_mul_left (waste + 1) hpotential)
  have hrunInput :
      CarriesInput x (attemptRun (fixedProgram phi) N sigma s) :=
    (fixedAttemptRun_preservesInputObservation phi N sigma s).trans hinput
  refine ⟨N, hNtriangle, hcons, hconstant, hrunInput, ?_⟩
  intro t ht htInput
  exact fixedConsensus_output_unique phi
    (attemptRun (fixedProgram phi) N sigma s) t hcons ht
    (hrunInput.trans htInput.symm)

/-- The explicit round-robin scheduler discharges the fixed-node bounded-waste
premise, yielding a finite mathematical attempt horizon and unique output. -/
theorem roundRobin_reaches_bounded_unique_output {k : Nat}
    (phi : Formula k) (x : Fin k → Bool) (s : State)
    (hinput : CarriesInput x s) :
    ∃ N,
      N ≤ (((fixedProgram phi).length - 1) + 1) *
        triangle (fixedProgram phi).length ∧
      Consensus (fixedFederation phi)
        (attemptRun (fixedProgram phi) N (fixedRoundRobinScheduler phi) s) ∧
      (∀ n, N ≤ n →
        attemptRun (fixedProgram phi) n (fixedRoundRobinScheduler phi) s =
          attemptRun (fixedProgram phi) N (fixedRoundRobinScheduler phi) s) ∧
      CarriesInput x
        (attemptRun (fixedProgram phi) N (fixedRoundRobinScheduler phi) s) ∧
      ∀ t : State,
        Consensus (fixedFederation phi) t →
        CarriesInput x t →
        attemptRun (fixedProgram phi) N (fixedRoundRobinScheduler phi) s
            (fixedOutReg phi) = t (fixedOutReg phi) :=
  boundedWaste_reaches_unique_output phi x s
    (fixedRoundRobinScheduler phi) ((fixedProgram phi).length - 1)
    hinput (fixedRoundRobin_nodeBoundedWaste phi)

#print axioms roundRobin_reaches_correct_output
#print axioms roundRobin_reaches_correct_output_within_linear
#print axioms canonicalAcceptedSteps_within_quadratic
#print axioms boundedWaste_reaches_unique_output
#print axioms roundRobin_reaches_bounded_unique_output

end

end OPH.Tower.FixedFederationExecutionEndpoint
