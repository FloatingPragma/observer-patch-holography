import Computation.FixedFederationProgress

/-!
# Fixed computation-federation endpoints

This Tower consumer packages the fixed-program progress theorem as a stable
same-input output endpoint.  It does not construct the finite `PublicWorld`
presentation used by `Tower.FixedPointEndpoint`; the raw state remains
`Nat -> Bool`, and unused registers are intentionally outside the endpoint
equivalence.
-/

namespace OPH.Tower.FixedFederationEndpoint

open OPH.RepairUniversality
open OPH.RepairUniversality.FixedFederation

noncomputable section

/-- Every pathwise weak-fair attempt scheduler reaches a stable consensus
whose output agrees with every other consensus carrying the same input. -/
theorem fairAttempt_reaches_unique_output {k : Nat} (phi : Formula k)
    (x : Fin k → Bool) (s : State)
    (sigma : NodeScheduler (fixedProgram phi))
    (hinput : CarriesInput x s)
    (hfair : NodePathwiseWeakFair (fixedProgram phi) sigma s) :
    ∃ N,
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
  obtain ⟨N, hcons, hconstant⟩ :=
    nodePathwiseWeakFair_eventually_consensus
      (fixedProgram phi) (fixedProgram_wf phi) sigma s hfair
  have hrunInput :
      CarriesInput x (attemptRun (fixedProgram phi) N sigma s) := by
    exact (fixedAttemptRun_preservesInputObservation phi N sigma s).trans hinput
  refine ⟨N, hcons, hconstant, hrunInput, ?_⟩
  intro t ht htInput
  exact fixedConsensus_output_unique phi
    (attemptRun (fixedProgram phi) N sigma s) t hcons ht
    (hrunInput.trans htInput.symm)

/-- The exact observation-determined-normal-forms consumer exposed at the
Tower boundary. -/
theorem sameInput_normalEndpoints_unique {k : Nat} (phi : Formula k) :
    ObservableNormalForms.ObserverEndpointUniqueModulo
      (CanonicalAcceptedStep (fixedProgram phi))
      (InputObservation k)
      (fun s t => s (fixedOutReg phi) = t (fixedOutReg phi)) :=
  fixedObserverEndpointUniqueOutput phi

#print axioms fairAttempt_reaches_unique_output
#print axioms sameInput_normalEndpoints_unique

end

end OPH.Tower.FixedFederationEndpoint
