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
| Minimal intrinsic moduli, `prop:minimal-moduli` | `residualModulus_mono`; `residualModulus_zero`; `residualModulus_le_of_errorBoundWitness`; `errorBoundWitness_residualModulus`; `inverseModulus_mono`; `inverseModulus_le_of_inverseObservationBound`; `inverseObservationBound_inverseModulus`; `inverseModulus_zero_iff_injOn`; `exists_separation_radius` | ✅ | Intrinsic finite maxima defined (`residualModulus`, `inverseModulus`, `settledModulus`); every extremum's admissible set carries a discharged nonemptiness witness; the assumed Stability certificates are instantiated by `intrinsic_two_output_estimate`. |
| Heterogeneous two-output estimate, `thm:master-bound` | `heterogeneous_two_output_estimate` | ✅ | Full triangle proof from attained error-bound witnesses, inverse-observation bound, monotone modulus, and Lipschitz bound. Finite intrinsic moduli instantiate these certificates. |
| Symmetric estimate, `eq:symmetric-bound` | `symmetric_two_output_estimate` | ✅ | Equal-residual specialization. |
| Rate transfer, `cor:rate-transfer` | `rate_transfer` | ✅ | Power-rate substitution (real exponents via `rpow`) into the intrinsic master estimate. |
| Sharpness, `thm:sharpness` | `settledModulus_zero_residual`; `settledModulus_le_symmetric_bound`; `SharpnessWitness.symmetric_bound_attained`; `SharpnessWitness.coefficient_two_sharp`; `NegativeControl.bound_not_attained` | ✅ | Exact zero-residual equality `Θ(0,r) = ωB(r)`; three-point witness with computed values `ηΦ(1)=1`, `ωB(0)=0`, `Θ(1,0)=2` defeating every coefficient `a < 2`; two-point negative control where the bound is strictly slack. |
| Approximate schedule independence, `cor:schedule` | `approximate_schedule_independence` | ✅ | Endpoint theorem; makes no confluence or common-path assumption. |
| High-probability schedules, `cor:probabilistic` | — | ⬜ | Union-bound and diameter-expectation wrapper not encoded. |
| Finite Markov drift receipt, `prop:markov-receipt` | `FiniteMarkovKernel.finite_markov_endpoint_receipt`; `FiniteMarkovKernel.endpointDistribution_tail_bound`; `FiniteMarkovKernel.finiteExpectation_endpointDistribution`; `FiniteMarkovKernel.endpointDistribution_eq_pathWeight_sum`; `FiniteMarkovKernel.observation_eq_of_endpointDistribution_pos` | ✅ | Complete package: explicit finite path law (product weights on `Fin n → S`), endpoint distribution proved as its endpoint marginal, expectation identity to the drift iterate, one-time tail and settling bounds, and positive-support (rewrite-or-stutter) observation preservation along supported paths. One-time endpoint bound only, matching the proposition's limitation clause; support-checker rejection tests in `SupportAudit`. |
| Stationary and occupation receipts, `cor:markov-occupation` | `FiniteMarkovKernel.stationary_expectation_bound`; `FiniteMarkovKernel.stationary_tail_bound`; `FiniteMarkovKernel.occupation_tail_bound`; `FiniteMarkovKernel.occupation_tail_bound_average`; `PersistentNoiseControl.noisy_confinement_vanishes` | ✅ | Stationary expectation/tail bounds and summed/averaged occupation bounds. No infinite-horizon confinement claim: the persistent-noise control proves a drift kernel with `ξ > 0` whose one-time excursion probability stays `1/2` at every time while window confinement `(1/2)^N → 0` (`confinedMass_eq_pathWeight_sum` grounds the window event in the same path law). |
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
| Canonical adaptive repair stratification, `thm:oph-adaptive-repair-stratification` | `OPH.Execution.adaptiveRun_eventually_constant`; `OPH.Execution.pathwiseWeakFair_eventual_normalForm`; `OPH.Execution.workConserving_eventually_normal` | ✅ | Every adaptive attempt stream is eventually constant; pathwise weak fairness supplies normality; work conservation gives a normal horizon bounded by the initial mismatch count. No stochastic kernel, rate, physical scheduler, or refinement result follows. |
| Cumulative attempt-capacity classification, `thm:oph-cumulative-attempt-capacity` | `OPH.Execution.cumulativeGenuineChangeCost_le_initialMismatch`; `boundedWaste_eventually_normal`; `CumulativeCapacityExamples.delayed_normalizing_attempt_no_go`; `delayThenProbe_pathwiseWeakFair`; `firstBroken_attempt_threshold_iff`; `twoCell_attempt_threshold_iff` | ✅ | Every scheduler invocation costs one unit, including stutter. The initial mismatch rank directly bounds genuine changes, but no single mismatch-only finite attempt horizon uniformly bounds all normalizing schedulers, including a pathwise weak-fair delayed family. Bounded waste gives the `(q + 1)`-scaled horizon; two nonisomorphic canonical finite sources make the work-conserving threshold sharp. No physical time, rate, or resource interpretation follows. |
| Adaptive canonical public endpoint | `OPH.Tower.AdaptiveEndpoint.adaptive_endpoint_mem_fixedPointObject`; `OPH.Tower.AdaptiveEndpoint.adaptive_public_endpoint_exists_unique` | ✅ | Completeness makes the weak-fair stable normal form consistent; confluence identifies it with the canonical `Repair` endpoint and therefore the same public fixed-point object. |
| Capacity-bounded canonical public endpoint | `OPH.Tower.CumulativeCapacityEndpoint.boundedWaste_endpoint_within_attempt_budget`; `boundedWaste_public_endpoint_exists_unique` | ✅ | Bounded waste and sufficient attempt budget produce the existing consistent raw and public fixed-point endpoint by an index no larger than the budget; the scheduler and budget sources remain declared premises. |

## Conditional mechanism-variant comparison

The generic declarations are part of the library root. The finite A2/A5
module is a designated build target and a conditional reference model for
Will's shared A2/A5 mechanism-design fixture using the OPH calculus. Its domain
mappings are explicit, replaceable assumptions. It is imported by
`AxiomAudit.lean` so that the complete theorem set receives dependency output.

| Result | Lean declaration | Status | Scope |
|---|---|---:|---|
| Variant-indexed trace projection | `MechanismVariants.trace_to_reflTransGen`; `stateObservation_eq_of_trace` | ✅ | Data-carrying traces start in an initial state and project to the declared `ReflTransGen`; state-observation equality reuses the existing step-preservation lemma. Path-sensitive observation still needs its own trace adapter. |
| Encoded collar and fixed-relation repair | `MechanismVariants.admissibleCollar_iff_encodedCollarSurjective`; `strongRepair_exists_iff_encodedCollarSurjective`; `strongRepair_exists_iff_full_admissibleCollar`; `no_strongRepair_of_missing_admissibleCollar` | ✅ | `RelationEncoding.exactRange` now composes the actual admissible-trace encoding with full collar coverage and `StrongRepair`. The negative theorem remains a static fixed-`R` hole. No same-initial executable repair, all-path completeness, or full `Repair` width follows. |
| Typed comparison provenance and family boundary | `MechanismVariants.ComparisonPolicy`; `behaviorLT_of_subset_of_witness`; `IsBehaviorMinimum` | ✅ | A policy fixes target, observation, protected set, authority view, selected support, cost source, and family. Cost is provenance only and is absent from `cut` and `PairwiseLT`. A family minimum separately requires comparison with every declared member. |
| Conditional A2 broad and A5 flash targets | `AmdA2A5Conditional.a2_broad_target_eliminated`; `a5_flash_target_eliminated`; `a5_broad_target_failure` | ✅ | Finite branch controls only. A5 blocks the modeled contract-origin flash target but the modeled EOA broad-target trace remains. |
| Raw and quotiented authority | `AmdA2A5Conditional.a2_zero_new_raw_capability_fails`; `a2_zero_new_authority_class`; `a5_zero_new_authority_class` | ✅ | The A2 model adds raw `blockAcquisition`; one explicit quotient maps it to an existing owner-control class. Both policy candidates pass only under that declared quotient, which is a replaceable assumption rather than confirmed authority semantics. |
| Outcome, full, and hybrid observation | `AmdA2A5Conditional.observation_profile_changes_behavior_judgment`; `a2_pairwise_lt_a5_hybrid_flash`; `hybrid_forgets_pending_but_retains_caller` | ✅ | Outcome forgets caller and delay; full observation retains both; hybrid forgets pending delay but retains caller kind. The pairwise result is profile- and target-relative. |
| Strategy separation | `AmdA2A5Conditional.observation_preservation_does_not_supply_strategy`; `strategic_choice_changes_under_delay_cost` | ✅ | A synthetic two-strategy delay-cost control changes the preferred action while protected outcome behavior is preserved. No automatic observation-to-strategy bridge is claimed. |
| Separate typed pairwise policies | `AmdA2A5Conditional.broadOutcomePolicy_a2_lt_a5`; `flashHybridPolicy_a2_lt_a5` | ✅ | Broad/outcome and flash/hybrid use distinct `ComparisonPolicy` values. Each theorem carries family membership, selected-support membership, and quotient-authority eligibility. Their cost field records provenance and is not a load-bearing order premise. |
| Encoded repair and protected representatives | `AmdA2A5Conditional.a2FullRelationEncoding_exactRange`; `a2FullEncode_outcomeCollar`; `a2_full_repair_iff_protected_outcome_representatives` | ✅ | The full-support relation is the exact range of admissible A2 traces, the Bool collar is injectively tied to `ObsOutcome`, and fixed-relation `StrongRepair` exists iff every protected outcome has an admissible representative. |
| Nonvacuity, missing support, and simultaneous controls | `AmdA2A5Conditional.full_control_activation_satisfiable`; `fixed_relation_positive_repairability`; `missing_collar_no_repair`; `support_relations_are_distinct` | ✅ | One `FullControlActivation` value carries all declared outcomes, including two separate policies. Missing support is an active, distinct static relation with a collar hole; no fake exact-range encoding or complete support-family claim is made. |
| Pairwise is not minimum | `AmdA2A5Conditional.a2_pairwise_lt_a5_in_declared_family`; `pairwise_does_not_establish_family_minimum` | ✅ | A synthetic third variant has a smaller cut, so A2 is strictly below A5 pairwise but is not minimum of that declared family. No unique or global minimum is claimed. |

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

## Protected-obstruction mechanism transport

Lean is canonical for the stochastic foundation and transport results in this
section.
A separately maintained Isabelle/HOL companion is available at the
[pinned public reference](https://github.com/Oraclizer/formal-verification/blob/d61f4dcf8b3f55a50ddeb2494049af322c5b7ec1/Protected_Behavior_Obstructions/README.md).
It is explicitly `PARTIAL / NO SAME`: it checks related set/profile consequences
under stated locale assumptions and is not an independent reconstruction of
the Lean kernel, `pathWeight`, first-hit law, scheduler, behavior-cut adapter,
or quantitative pushforward.

Finite first-hit expansion, positive-path reachability, the finite
closed-class characterization of almost-sure hitting, least nonnegative
Bellman fixed points, and strong lumpability are standard finite Markov-chain
facts. The OPH-specific result is the protected-observation and active-source
profile that separates four failure modes, recovers the existing behavior-cut
interface of `MechanismVariants.lean` through explicit supported traces, and
transports the complete source/target/quotient packet without hiding those
assumptions. In that recovery, the L2 and L3 admissibility profiles carry
source-quantified conditions beyond the trace witness; almost-sure hitting has
no single-trace witness, so those coordinates enter the profile as side
conditions permitted by the interface.

| Result | Lean declaration | Lean status | Isabelle boundary |
|---|---|---:|---|
| First-hit foundation | `ProtectedObstructions.FirstHitAt`; `firstHitAtMass`; `cumulativeEndpointMass`; `hitBy`; `hitProbability`; `endpointHitProbability`; `hitProbability_bellman` | ✅ | Lean-only canonical construction from `FiniteMarkovKernel.pathWeight` and a bounded monotone limit. |
| T0 partition | `ProtectedObstructions.PublicAdapter.Profile.t0_partition`; `t0_pairwise_disjoint`; `GenericCompletion.Profile.cut0_subset_cut1`; `cut1_subset_cut2`; `cut2_subset_cut3`; `delta1_eq_cut1_diff_cut0`; `delta2_eq_cut2_diff_cut1`; `delta3_eq_cut3_diff_cut2` | ✅ | `PARTIAL`: Isabelle proves nested-set partition algebra over supplied profile functions; Lean composes it with the kernel-derived layer definitions. |
| T1 empty fiber | `ProtectedObstructions.PublicAdapter.Profile.t1_fiber_empty_iff`; `t1_no_normal_endpoint_of_empty_fiber` | ✅ | `PARTIAL`: Isabelle proves only static target-set emptiness; no stochastic content is claimed for this row. |
| T2 positive first hit | `ProtectedObstructions.endpointHitProbability_pos_iff`; `PublicAdapter.Profile.t2_positive_reach_iff` | ✅ | `PARTIAL`: Isabelle theorem `protected_profile.T2_positive_endpoint` consumes `hit_pos_iff_endpoint`; it does not derive positivity from a kernel or `pathWeight`. |
| T3 almost-sure first hit | `ProtectedObstructions.hitProbability_eq_one_iff_no_reachableClosedTrap`; `PublicAdapter.Profile.t3_almostSure_iff_noClosedTrap` | ✅ | `PARTIAL`: Isabelle theorem `protected_profile.T3_all_sources_no_closed_trap` consumes `hit_one_iff_no_closed_trap`; it has no pre-hit reach graph or finite-kernel proof. |
| Least nonnegative fixed point and trap control | `ProtectedObstructions.hitProbability_least_nonnegative_fixedPoint`; `nonunique_nonnegativeBellmanFixedPoint_of_reachableClosedTrap` | ✅ | Lean-only. The least-fixed-point theorem is derived from the path law; no general uniqueness is claimed in a closed non-target trap. |
| T4 endpoint selection | `ProtectedObstructions.PublicAdapter.Profile.t4_selection_iff_endpointUnionCert` | ✅ | `PARTIAL`: Isabelle unfolds a certificate over an opaque supplied `first_hit` function. |
| T5 observable collapse | `ProtectedObstructions.PublicAdapter.Profile.t5_observableDetermination_collapse` | ✅ | `PARTIAL`: Isabelle proves a profile consequence under locale laws and an explicit determinacy premise. |
| T6 scheduler correspondence | `ProtectedObstructions.PublicAdapter.Profile.t6_support_to_rewrite`; `t6_rewrite_to_support`; `target_is_normal` | ✅ | Lean-only. Soundness, completeness, and target normality are source proofs; T6 has no repeated conclusion fields. |
| T7 behavior-cut recovery | `ProtectedObstructions.PublicAdapter.Profile.t7_behaviorCut0`; `t7_behaviorCut1`; `t7_behaviorCut2`; `t7_behaviorCut3`; `GenericCompletion.VariantFamily.cutAt_zero`; `cutAt_one`; `cutAt_two`; `cutAt_three` | ✅ | Lean-only OPH `RealizedBehavior`/`BehaviorCut` adapter with one structural protected carrier; the L2 and L3 admissibility profiles carry source-quantified side conditions beyond the trace witness. |
| T8 product preorder | `ProtectedObstructions.PublicAdapter.t8_coordinate_strict_iff`; `GenericCompletion.VariantFamily.le_refl`; `le_trans` | ✅ | `PARTIAL`: Isabelle proves four-coordinate preorder facts while `profile_correspondence` assumes exact cut membership; Lean proves the strict-coordinate characterization. |
| Structural morphism exactness | `ProtectedObstructions.NonidentityExactness.StructuralMorphism.fiber_nonempty_exact_derived`; `hitProbability_exact_derived`; `endpointHitProbability_pushforward_derived`; `positive_exact_derived`; `almostSure_exact_derived`; `endpointUnique_exact_derived`; `layer0_exact`; `layer1_exact`; `layer2_exact`; `layer3_exact` | ✅ | `PARTIAL`: Isabelle proves only protected/initial/target/silent set-level identity and composition. It has no kernel lumping or law pushforward. |
| Full-support quantitative transport | `ProtectedObstructions.PublicAdapter.GenericCompletion.InitialLaw.pooledHit_pos_iff`; `pooledEndpoint_pos_iff`; `pooledEndpointUnique_iff`; `pooledHit_eq_one_iff`; `Quantitative.pooledHit_exact_derived`; `pooledEndpoint_pushforward_derived` | ✅ | Lean-only for positive stochastic transport and pooled endpoint pushforward. |
| Nonidentity and native fixtures | `ProtectedObstructions.NonidentityExactness.Fixture.nonidentity_layer3_exact`; `PublicAdapter.GenericCompletion.StrictReversal.strict_forward`; `strict_reverse`; `NativeTwoBitC4.M0.protected_in_delta0`; `TwoBit.fine_survivor`; `coarse_in_delta3`; `C4.full_L0_iff_strongRepair`; `NativeC6.m2_delta2`; `m3_delta3` | ✅ | The built Isabelle M2/M3 values remain `PARTIAL` direct finite functions rather than kernel-derived probabilities. |

### Native activation and counterexample route

| Reader-facing phenomenon | Key declarations | Demonstrates | Does not demonstrate |
|---|---|---|---|
| Fine/coarse TwoBit separation | `NativeTwoBitC4.TwoBit.fine_survivor`; `coarse_in_delta3`; `TwoBitRepair.coarse_confluent` | Observable determination reaches the final layer for the fine observation, while coarse confluence alone leaves endpoint-class ambiguity. | Confluence does not imply protected selection under every observation. |
| Full versus missing fixed relation | `NativeTwoBitC4.C4.full_L0_iff_strongRepair`; `missing_in_staticK0`; `filled_relation_not_staticK0` | The first cut matches the static collar-coverage obstruction for the concrete relation. | No liveness, stochastic-selection, or strategic-completeness claim. |
| Closed trap versus ambiguous endpoint | `NativeC6.m2_delta2`; `m3_delta3`; `m2_t6_support_to_rewrite`; `m3_t6_both_directions` | Proper-target fixtures separate positive reach, almost-sure reach, and quotient-endpoint selection while activating both scheduler directions. | No rate, expected-time, mixing, or infinite-state conclusion. |
| Cross-layer strict reversal | `GenericCompletion.StrictReversal.strict_forward`; `strict_reverse`; `not_profileLE_false_true`; `not_profileLE_true_false` | Different cut coordinates can reverse pairwise preference, so the comparison is a product preorder. | No scalar score, total order, or global optimum. |
| Nonidentity exact transport | `NonidentityExactness.Fixture.nonidentity_layer3_exact`; `Fixture.nonidentity_fiber_exact`; `NamedFixtures.fixture_source_endpoint_push` | The complete finite morphism packet preserves the protected profile across a nonidentity state map. | No arbitrary implementation or product-refinement invariance. |

The finite explorer independently checks the minimal fixture counts and the
source-bound negative controls:

```sh
python tools/verify_protected_obstruction_models.py
```

Explicitly not formalized:

- asynchronous ranked-functional update counts and the path-gain estimate;
- linear quotient/rank/singular-value results;
- succinct complexity classifications.
