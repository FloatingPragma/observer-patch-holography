# Issue 551: Complete-Class Capacity Non-Identifiability

## Result

The complete declared fixed-packet structure lifts across the
capacity-indexed generation-register family from one target-clean source
rule, and the completed class does not entail a unique slack zero.

The lift transports every structural component of the issue #548 packet per
rung: terminal-fiber completeness with the structural one-fault gate and
rung faults, observer and interface atom maps, exact equalizer public
sections, endogenous reachability histories, the frozen publicness family,
continuation-manifest closure, joint kernels with local-marginal
consistency, A2 meaning-map naturality including the extension square, A3
feasible-set state determinacy, capacity-extension embeddings with no new
confusability, fixed-rung refinement stability, and exact fiber-product
sewing.

Two admissibility readings of the continuation manifest carry the verdict:

| Reading | Admissible completions | Slack zero set | Consequence |
|---|---|---|---|
| source-closed | compositions of the forty declared reversible generators | every positive rung, for every admissible completion | slack vanishes identically; no unique zero exists |
| widened | completion kernels beyond the generator vocabulary, when every transported control passes | identity: all rungs; erasure: rung 1; two-class cap: rungs 1 and 2 | surviving zero sets are inequivalent; no unique zero is entailed |

Two candidate directions fail named transported controls, with exact
witnesses in the receipt: hidden spectator multiplicity above one fails A3
state determinacy (4,096 raw families per public class at multiplicity two),
and parity oscillation fails the A2 extension square and the
no-new-confusability control at the odd-to-even rung step.

Under both readings the completed declared source class does not select a
unique capacity. A rung selector is an additional named source law, which is
a new physical premise rather than an open calculation.

## Evidence

- `complete_packet_capacity_lift.py` builds the per-rung structures, checks
  every transported control, and emits the receipt and certificate under
  `runtime/`.
- `verify_complete_packet_lift_independent.py` recomputes the capacities,
  zero sets, survivor logic, and verdict implication without importing the
  producer.
- `capacity_indexed_source_family.py` remains the bounded-layer producer;
  the complete lift cross-checks every shared capacity value against it.
- `direct_n_closure_verdict.py` consumes both layers and emits the issue
  #505 verdict `LOCKED_NONIDENTIFIABILITY_COMPLETED_CAPACITY_SOURCE_CLASS`.
- `Lean/ObserverPatchHolography/CapacityNonidentifiability.lean` proves the
  arithmetic skeleton: `sourceClosed_no_unique_positive_fixed_rung`,
  `oscillation_fixed_iff_odd`, `oscillation_identity_differentFixedSets`,
  and `completeClass_doesNotEntailUniqueZero`, in addition to the bounded
  all-rung theorems.
- `test_complete_packet_capacity_lift.py` and
  `test_direct_n_closure_verdict.py` exercise the controls, the adversarial
  mutations, byte stability, and the independent verifier.

## Scientific boundary

No numeric \(N\) is emitted, no cosmological comparison is permitted, and
the strange-loop identity is not rejected. The frozen family is the
generation-register extension of the issue #548 packet; non-entailment
within one admissible declared family implies non-entailment for every
wider same-antecedent class. The common \(k=1\) screen packet is not a
cosmic selection, and promoting any capacity to a universe-level carrier
requires a separately named source law together with the physical
attachments owned by the gravitational chain.
