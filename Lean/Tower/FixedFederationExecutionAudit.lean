import Tower.FixedFederationExecutionEndpoint
import Tower.CumulativeCapacityEndpoint
import Tower.FixedFederationFanoutEndpoint

/-!
# Typed audit for fixed-federation execution classification

These examples pin the exact public theorem types consumed by the Tower.  The
module is compiled by the normal Lean build and by the Lean CI kernel-check
step, so theorem-name text alone cannot satisfy the release check.
-/

namespace OPH.Tower.FixedFederationExecutionAudit

open OPH.RepairUniversality
open OPH.RepairUniversality.FixedFederation
open OPH.Tower.FixedFederationExecutionEndpoint
open OPH.Tower.FixedFederationFanoutEndpoint

noncomputable section

example (L : List Node) (hL : L ≠ []) (s : State) :
    attemptRun L L.length (roundRobinScheduler L hL) s =
      sweepFrom s L :=
  roundRobinScheduler_cycle_eq_sweepFrom L hL s

example {k : Nat} (phi : Formula k) (s : State) :
    Consensus (fixedFederation phi)
      (attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s) :=
  fixedRoundRobin_consensus_after_one_cycle phi s

example {k : Nat} (phi : Formula k) (x : Fin k → Bool) (s : State)
    (hinput : CarriesInput x s) :
    attemptRun (fixedProgram phi) (fixedProgram phi).length
        (fixedRoundRobinScheduler phi) s (fixedOutReg phi) =
      Formula.evalF phi x :=
  fixedRoundRobin_output_after_one_cycle phi x s hinput

example {k : Nat} (phi : Formula k) (x : Fin k → Bool) (s : State)
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
      Formula.evalF phi x :=
  roundRobin_reaches_correct_output_within_linear phi x s hinput

example {k : Nat} (phi : Formula k) {m : Nat} {s t : State}
    (hsteps : CanonicalAcceptedSteps (fixedProgram phi) m s t) :
    m ≤ (fixedProgram phi).length *
      ((fixedProgram phi).length + 1) / 2 :=
  fixedProgram_acceptedSteps_quadratic phi hsteps

example : ¬ GoUseProfile allocatorCollisionFormula 0 :=
  allocatorSeparation_needed

example :
    ∃ sigma : NodeScheduler (fixedProgram recurrenceGapFormula),
      NodePathwiseWeakFair (fixedProgram recurrenceGapFormula) sigma allFalse ∧
      ¬ ∃ N,
        N ≤ (0 + 1) * triangle (fixedProgram recurrenceGapFormula).length ∧
        Consensus (fixedFederation recurrenceGapFormula)
          (attemptRun (fixedProgram recurrenceGapFormula) N sigma allFalse) ∧
        ∀ n, N ≤ n →
          attemptRun (fixedProgram recurrenceGapFormula) n sigma allFalse =
            attemptRun (fixedProgram recurrenceGapFormula) N sigma allFalse :=
  boundedWaste_premise_needed

example (n : Nat) :
    (∀ (L : List Node), NodesWF L → ∀ {m : Nat} {s t : State},
      CanonicalAcceptedSteps L m s t → m < 2 ^ L.length) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain n) (2 ^ n - 1) allFalse t ∧
      Consensus ((fanoutChain n).map Node.obs) t ∧
      NodesWF (fanoutChain n) ∧
      (fanoutChain n).length = n) :=
  fanoutChain_sharp_exponential n

example :
    (∀ {k : Nat} (phi : Formula k) {m : Nat} {s t : State},
      CanonicalAcceptedSteps (fixedProgram phi) m s t →
      m ≤ triangle (fixedProgram phi).length) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 3) 7 allFalse t ∧
      Consensus ((fanoutChain 3).map Node.obs) t ∧
      triangle (fanoutChain 3).length = 6 ∧
      ¬ (7 ≤ triangle (fanoutChain 3).length)) ∧
    NodesWF (fanoutChain 3) ∧
    ¬ AtMostOneDownstreamConsumer (fanoutChain 3) :=
  triangular_bound_needs_single_consumer

end

end OPH.Tower.FixedFederationExecutionAudit
