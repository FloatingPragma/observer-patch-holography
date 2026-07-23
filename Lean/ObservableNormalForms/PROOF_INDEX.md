# Observation-Determined Normal Forms proof index

This index follows the current manuscript by theorem name and LaTeX label;
numeric theorem counters are intentionally omitted because editorial insertions
can renumber shared theorem environments.

Status legend:

- ✅ theorem statement and proof formalized, with no `sorry`
- 🟡 load-bearing mathematical core formalized; manuscript-specific packaging
  remains outside Lean
- ⬜ not formalized

## Exact observable fibers

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Universal property, (a) ↔ (b), `thm:universal` | `observableDetermined_iff_unique_partialNormalizer` | ✅ | Unique proof-carrying `Option` normalizer, including obstruction, fixed-point, soundness, and boundary-extensionality clauses. |
| Universal property, (a) ↔ (c) | — | ⬜ | Image-subtype bijection is standard but not encoded. |
| Reachable normal forms live in the observable fiber, `prop:reachable-fiber` | `reachable_normalForm_mem_fiber` | ✅ | Reflexive-transitive closure, one-step observation preservation, and exact normal-form completeness. |
| Empty/singleton alternative, `cor:empty-singleton` | `empty_fiber_no_reachable_normalForm`; `singleton_fiber_forces_normalForm`; `singleton_fiber_weakNormalization_confluent` | ✅ | All three clauses; weak normalization is assumed only on the singleton fiber. |
| Audited terminal alternative, `thm:audited-terminal` | `existsUnique_auditedTerminal` | ✅ | Unique terminal, preserved observation, and consistency iff the input fiber is realizable, under exactly the current theorem hypotheses. |
| Observation-leak separation example | `ObservationLeakCounterexample.singleton_fiber_insufficient_without_observation_preservation` | ✅ | Complete three-state system with a singleton consistent fiber whose rewrite escapes to another observation class. |

## Observation-relative endpoint uniqueness

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Cross-source uniqueness modulo an equivalence, `thm:cross-source-modulo` | `boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo` | ✅ | Exact equivalence between identification of consistent states modulo an arbitrary relation and cross-source normal-endpoint uniqueness, assuming observation preservation and exact normal-form completeness. |
| Weak-normalization endpoint existence | `exists_equivalent_observer_endpoints` | ✅ | Every equally observed source pair has at least one pair of normal endpoints, which are equivalent; universal comparison is supplied by the preceding iff theorem. |
| Fine/coarse two-bit separation | `TwoBitRepair.observerEndpointUnique`; `TwoBitRepair.exists_equal_endpoints`; `TwoBitRepair.coarse_confluent`; `TwoBitRepair.coarse_boundary_does_not_identify`; `TwoBitRepair.coarse_observerEndpointUnique_fails` | ✅ | Complete one-step repair system with a protected bit: the relation remains same-source confluent, while discarding the protected bit fails both boundary identification and cross-source endpoint uniqueness. |

## Observational stability

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Minimal intrinsic moduli, `prop:minimal-moduli` | — | ⬜ | Finite maxima, minimality, zero values, and finite separation radius remain to be encoded. |
| Heterogeneous two-output estimate, `thm:master-bound` | `heterogeneous_two_output_estimate` | ✅ | Full triangle proof from attained error-bound witnesses, inverse-observation bound, monotone modulus, and Lipschitz bound. Finite intrinsic moduli instantiate these certificates. |
| Symmetric estimate, `eq:symmetric-bound` | `symmetric_two_output_estimate` | ✅ | Equal-residual specialization. |
| Rate transfer, `cor:rate-transfer` | — | ⬜ | Rate substitution not encoded. |
| Sharpness, `thm:sharpness` | — | ⬜ | Exact zero-residual modulus and sharp three-point example not encoded. |
| Approximate schedule independence, `cor:schedule` | `approximate_schedule_independence` | ✅ | Endpoint theorem; makes no confluence or common-path assumption. |
| High-probability schedules, `cor:probabilistic` | — | ⬜ | Union-bound and diameter-expectation wrapper not encoded. |
| Finite Markov drift receipt, `prop:markov-receipt` | `FiniteMarkovKernel.finite_markov_drift_iteration`; `finite_markov_tail_bound` | 🟡 | Exact finite Markov-operator drift iteration and one-time finite-distribution Markov inequality. The kernel-law/event wrapper and positive-support observation-preservation clause are not packaged as one theorem. |
| Compact-space extension / uniform family, `prop:compact-extension`, `thm:uniform-family` | — | ⬜ | Not encoded. |
| Exact product calculus, `thm:product-calculus` | — | ⬜ | Current manuscript correctly requires a nonempty finite index family. |
| Sensor enrichment, `cor:sensor-enrichment` | `inverse_bound_of_sensor_enrichment` | 🟡 | Certificate-level monotonicity proved; equality of intrinsic finite maxima is not encoded. |

## Naturality and refinement towers

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| One-step approximate naturality, `thm:one-step-naturality` | `one_step_approximate_naturality` | ✅ | Full semantic-defect proof with normalizer consistency and observation contracts. |
| Exact naturality, `cor:exact-naturality` | `exact_naturality_from_uniqueness` | ✅ | No metric assumptions needed. |
| Consistency-model perturbation, `cor:model-perturbation` | — | ⬜ | Hausdorff perturbation wrapper not encoded. |
| Telescoping refinement bound, `thm:tower` | `telescoping_refinement_bound` | 🟡 | Arbitrary-depth metric telescope. A dependent family of stage types and restriction compositions is not packaged. |
| Summable/modulus tower corollaries, `cor:summable`, `cor:modulus-tower` | — | ⬜ | Not encoded. |
| Receipt-to-exact comparison, `lem:receipt-to-exact` | `heterogeneous_two_output_estimate` | 🟡 | Follows by the manuscript specialization; a dedicated receipt structure/wrapper is not encoded. |
| Same-level implementation agreement, `prop:same-level-implementations` | `same_level_implementation_agreement` | ✅ | Full metric/Lipschitz core with no refinement-tail term. |
| Pathwise telescope to an anchor, `lem:pathwise-anchor` | `telescoping_refinement_bound` | 🟡 | The arbitrary chain sum is proved; path-dependent `α_j` and dependent restrictions are not packaged. |
| Fine-to-coarse certificate, `cor:fine-to-coarse-solver` | `projective_implementation_bound_from_tower_receipt` | 🟡 | Solver receipt plus an exact-tower receipt under a Lipschitz restriction. |
| Anchored cross-level comparison, `thm:anchored-cross-level` | `anchored_cross_level_metric_core` | 🟡 | Complete five-segment metric core. The stage-indexed restriction, modulus, and pathwise-`A` wrapper remains outside Lean. |
| Nested compatible levels, `cor:nested-cross-level` | `nested_compatible_levels_metric_core` | 🟡 | Complete three-segment metric core after compatibility removes the first path and anchor-mismatch terms. |
| Cofinal common limit, `cor:cofinal-projective` | — | ⬜ | Cauchy/completeness, cofinal-subsequence, and inverse-system compatibility proof not encoded. |
| Inverse-limit normalizer, `thm:inverse-limit` | — | ⬜ | Dependent inverse-limit construction not encoded. |
| Earlier-draft two-tail comparison | `two_implementations_bound_from_tower_receipts` | 🟡 | Supporting precursor only; it is not identified as a current paper theorem. |

## Repair and selection

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Collar-section criterion, `thm:collar-section` | `strongRepair_exists_iff_collarProjection_surjective` | ✅ | Uses the current nonempty-write hypothesis and strengthens the finite result classically by dropping finiteness. |
| No-repair certificate | `no_strongRepair_of_missing_collar` | ✅ | Nonempty write space. |
| Robust no-repair margin, `prop:robust-no-repair` | `robust_no_repair_margin`; `repairMargin_pos_of_compact` | ✅ | Uses the current nonempty-relation hypothesis; abstract collar map and dominating metric. |
| Empty-domain/empty-relation audits | `empty_write_space_counterexample`; `empty_relation_repairMargin_zero` | ✅ | Machine-checks why the nonemptiness qualifications are load-bearing. |
| Equivariant section / stabilizer obstruction, `thm:equivariant-section`, `thm:stabilizer` | — | ⬜ | Not encoded. |

## Finite conditional resampling

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Fiber averaging is conditional expectation, `thm:fiber-conditional-expectation` | `FiniteWeightedObservation.transition_nonneg`; `FiniteWeightedObservation.transition_sum_one`; `FiniteWeightedObservation.resample_eq_fiber_average`; `FiniteWeightedObservation.observationMeasurable_iff_exists_factor`; `FiniteWeightedObservation.resample_eq_self_iff_observationMeasurable`; `FiniteWeightedObservation.resample_idempotent`; `FiniteWeightedObservation.resample_weighted_self_adjoint`; `FiniteWeightedObservation.resample_weighted_energy_identity`; `FiniteWeightedObservation.resample_weighted_energy_le`; `FiniteWeightedObservation.kernel_eq_conditionalResamplingKernel_iff_recognition` | 🟡 | The finite algebraic characterization is formalized for strictly positive finite state weights: stochasticity, exact weighted fiber formula and fixed space, projector properties, Pythagorean identity, weighted-`L2` contraction, and the exact R1 fiber-support/R2 equal-row/R3 weighted-detailed-balance matrix-recognition equivalence. Equality with Mathlib's measure-theoretic `condexp` operator onto `σ(B)` is not wrapped. |

## Ranked functional systems and examples

| Paper result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Synchronous settling by dependency depth, `cor:functional-synchronous` | `RankedSynchronousSystem.synchronousEvolve_agrees_through_rank`; `synchronous_depth_settling` | ✅ | Heterogeneous site values and extensional strict-rank causality; one additional rank settles per round. |
| Generated-extension uniqueness | `RankedSynchronousSystem.generatedExtension_unique` | ✅ | Same boundary plus common finite rank bound. |
| Width-three Rule 90 | `Rule90.kernel_exact`; `image_exact`; `read01_injective_on_image`; `read02_not_injective_on_image`; `no_total_reverse_repair` | ✅ | Self-contained exact kernel/image/readout/reverse-repair statements; imports no other project library. |

Explicitly not formalized:

- asynchronous ranked-functional update counts and the path-gain estimate;
- linear quotient/rank/singular-value results;
- succinct complexity classifications.
