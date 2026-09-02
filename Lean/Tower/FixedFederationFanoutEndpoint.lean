import Computation.FixedFederationFanoutControl

/-!
# Tower endpoint: the compiler separation in one metric

This consumer states the upper and lower certificates of the fixed-federation
execution classification in a single metric, the number of nodes in the node
list.

* generic well-formed federations admit fewer than `2 ^ n` canonical accepted
  steps, and the `fanoutChain` family attains `2 ^ n - 1`, so that bound is
  sharp;
* the formula compiler output admits at most `n * (n + 1) / 2` canonical
  accepted steps;
* the separating family fails `AtMostOneDownstreamConsumer`, which is exactly
  the hypothesis the quadratic bound consumes.

The endpoint therefore records that the quadratic classification is a
consequence of the compiler's tree shape and not of well-formedness alone.
It makes no claim about scheduler attempts, about a lower bound for
`fixedProgram`, or about any physical realization.
-/

namespace OPH.Tower.FixedFederationFanoutEndpoint

open OPH.RepairUniversality
open OPH.RepairUniversality.FixedFederation

noncomputable section

/-- The full separation, stated in the node-count metric. -/
theorem fanout_vs_compiler_sharp_separation :
    (∀ (L : List Node), NodesWF L → ∀ {m : Nat} {s t : State},
      CanonicalAcceptedSteps L m s t → m < 2 ^ L.length) ∧
    (∀ n : Nat, ∃ t : State,
      NodesWF (fanoutChain n) ∧
      (fanoutChain n).length = n ∧
      CanonicalAcceptedSteps (fanoutChain n) (2 ^ n - 1) allFalse t ∧
      Consensus ((fanoutChain n).map Node.obs) t) ∧
    (∀ {k : Nat} (phi : Formula k) {m : Nat} {s t : State},
      CanonicalAcceptedSteps (fixedProgram phi) m s t →
      m ≤ (fixedProgram phi).length *
        ((fixedProgram phi).length + 1) / 2) ∧
    ¬ AtMostOneDownstreamConsumer (fanoutChain 3) := by
  refine ⟨fun L hWF _m _s _t h => canonicalAcceptedSteps_lt_pow L hWF h,
    ?_, ?_, fanoutChain_not_single_consumer⟩
  · intro n
    obtain ⟨t, hsteps, hcons⟩ := fanoutChain_exponential_lower n
    exact ⟨t, fanoutChain_wf n, fanoutChain_length n, hsteps, hcons⟩
  · intro k phi m s t h
    exact fixedProgram_acceptedSteps_quadratic phi h

/-- The hypothesis is load-bearing in the strongest available sense: at three
nodes the fanout family realizes seven canonical accepted steps while the
triangular budget for a three-node single-consumer program is six. -/
theorem triangular_bound_needs_single_consumer :
    (∀ {k : Nat} (phi : Formula k) {m : Nat} {s t : State},
      CanonicalAcceptedSteps (fixedProgram phi) m s t →
      m ≤ triangle (fixedProgram phi).length) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 3) 7 allFalse t ∧
      Consensus ((fanoutChain 3).map Node.obs) t ∧
      triangle (fanoutChain 3).length = 6 ∧
      ¬ (7 ≤ triangle (fanoutChain 3).length)) ∧
    NodesWF (fanoutChain 3) ∧
    ¬ AtMostOneDownstreamConsumer (fanoutChain 3) := by
  refine ⟨?_, fanoutChain_exceeds_triangle, fanoutChain_wf 3,
    fanoutChain_not_single_consumer⟩
  intro k phi m s t h
  exact fixedProgram_acceptedSteps_triangle phi h

/-- The same statement specialized to the exact sizes that are also checked
by kernel evaluation in the computation module. -/
theorem fanoutChain_small_witnesses :
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 2) 3 allFalse t ∧
      Consensus ((fanoutChain 2).map Node.obs) t) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 3) 7 allFalse t ∧
      Consensus ((fanoutChain 3).map Node.obs) t) ∧
    (∃ t : State,
      CanonicalAcceptedSteps (fanoutChain 4) 15 allFalse t ∧
      Consensus ((fanoutChain 4).map Node.obs) t) := by
  refine ⟨?_, ?_, ?_⟩
  · simpa using fanoutChain_exponential_lower 2
  · simpa using fanoutChain_exponential_lower 3
  · simpa using fanoutChain_exponential_lower 4

#print axioms fanout_vs_compiler_sharp_separation
#print axioms triangular_bound_needs_single_consumer
#print axioms fanoutChain_small_witnesses

end

end OPH.Tower.FixedFederationFanoutEndpoint
