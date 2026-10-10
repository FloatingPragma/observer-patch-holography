# OPH model search: what is constructed, what is missing

Status review: **10 October 2026**. Upstream reference:
`bdafda7de78b2a97f04883b062a7facbaabd84f8`.

This document records the physical-model research thread, including its source
reread and corrections. It is a construction and dependency ledger, not a new
axiom, theorem, prediction registration, or declaration that OPH describes
nature. The [axiom reference](AXIOM_REFERENCE.md),
[claim registry](../claims/claim_registry.yaml),
[physical-identification registry](../claims/physical_identification_registry.json)
and their linked proofs remain authoritative. Local research packets retain
their original source revisions; their results are not silently rebased.

**Present assessment:** OPH has concrete twelve-port carriers, repair
implementations, observer/record constructions, geometry and clock results,
and several field-theory constructions. A single physically identified
realization joining all the required components has not been established by
this campaign. The generic models in this campaign undercredited those existing results and
drifted away from the actual carrier. The next construction starts with the
documented OPH source.

## 1. Working rule: the simplest compatible construction

**Assume the simplest construction compatible with the actual OPH carrier and
all retained evidence as the working hypothesis. Start with the existing constructions. Add only what a specific failed interface demonstrably requires.**

This is a research preference, not a theorem that nature must select the
shortest formula. “Simple” does not mean deleting load-bearing hypotheses or
replacing the carrier with an easier unrelated graph.

1. **Use the smallest adequate existing source.** Preserve the oriented
   twelve-port incidence, actual local repair, records, neutral continuation,
   observer interfaces and declared refinement. Use a larger federation only
   when the question needs it.
2. **Try a readout or composition of existing objects first.** Before adding
   metric variables, memory, noise, kinetic terms or counterterms, test whether
   an existing response, record or operation supplies the missing object.
3. **Expose every choice.** List state variables, law, preparation, scheduler,
   measure, quotient, units, parameters and measurement maps. State which are
   axiom inputs, derived, branch hypotheses, calibrated or target-informed.
4. **Require a reason for every addition.** Name the obstruction, show why the
   unchanged construction fails to meet the required contract, and explain
   which minimal addition addresses it. Preserve existing working results.
5. **Test one common construction.** Matter, geometry, clock and entropy must
   be read from compatible states and operations with explicit maps. Separate
   successful models cannot be joined just by sharing notation or symmetry.
6. **Keep the comparison capable of failure.** Attempt a countermodel within
   the admitted assumptions; compare against an established-physics baseline
   and simpler controls. Do not retune after a failed comparison and call the
   original prediction successful.
7. **Use bounded decisions.** Each attempt has one missing connection, one
   decisive calculation, a computation cap and a stopping condition. The list
   of missing items below is not authorization for an unlimited prerequisite
   programme.

Repair can depend on the instantaneous state, local records and changing context.
That does not require an externally changing law. If a proposed law learns or
adapts, its adaptation state and update must also belong to the model.

## 2. What a complete candidate must specify

The microphysics paper separates these objects in its typed simulator:

| Object | Required meaning |
| --- | --- |
| Carrier and observer | Accessible algebras/state, oriented ports, overlap restrictions, readout, records, response and checkpoints. A carrier is not automatically an operational observer. |
| Proposal and accepted repair | Allowed local proposals; validated accepted transactions; protected quantities; actual descent score; completeness and confluence domain. |
| Normalizer | The quotient-visible map taking admissible data to the declared consistent normal form. |
| Equilibrium dynamics | State/ensemble, measure, transition law or generator, rates and invariants. This need not be the strict-descent relation. |
| Quantum dynamics | Noncommutative algebra/state, representation, positive transfer or self-adjoint Hamiltonian, instruments and relevant limits. |
| History and clock | Authenticated event dependence, event correspondence, duration readout and calibration. A repair attempt counter is not automatically physical time. |
| Geometry, fields and instruments | Operational geometry, matter/current/stress readouts, physical region/area maps and laboratory measurement interpretation. |
| Refinement | Maps between resolutions and quantitative compatibility of all objects used in the claim. A stated fixed-cutoff macroscopic regime is an alternative to an unjustified continuum claim. |

See the [microphysics paper](../paper/screen_microphysics_and_observer_synchronization.tex),
the [consensus paper](../paper/reality_as_consensus_protocol.tex), and their
[shared construction](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex).
OPH is invariant under presentation changes preserving the complete visible
contract. It is not neutral under arbitrary changes of carrier or law.

## 3. Constructed items and their exact boundaries

“Constructed” below means a specified mathematical object, theorem under its
stated premises, or implemented finite realization. It does not mean a
laboratory identification. These rows summarize source results; the document
does not independently reprove every theorem or rebuild Lean.

### Carrier, repair, observers and records

| ID | Constructed | Premises and boundary |
| --- | --- | --- |
| C01 | **Oriented icosahedral carrier.** Current A1 includes twelve primitive ports, thirty edges and twenty oriented faces, with typed seam/triple-overlap and support/refinement data. | The incidence is an axiom input, not a new prediction of the model search. Euclidean coordinates, scales and particle labels are not specified by the combinatorial packet. [Axiom reference](AXIOM_REFERENCE.md). |
| C02 | **Carrier geometry.** Inverse pairing, proper A5 symmetry, exact Gram relations `G²=4G`, `tr G=12`, rank three; a response-selected Euclidean completion of conservative integer records. | Counting selection uses its stated integer-load, total-charge and cost realization. The response completion takes the normalized slow-response limit; finite response kernels have ranks eleven/six on the indicated sectors. Physical position and common field/refinement attachment remain distinct. [Microphysics](../paper/screen_microphysics_and_observer_synchronization.tex), [flagship](../flagship/from_observer_consensus_to_standard_physics.tex), [selector implementation](../code/a5_closure/echosahedral_selector_certificate.py). |
| C03 | **Explicit recovery and accepted consensus.** Markov-splice/Petz-type candidates, finite decoding, protected transactional validation, strict descent and quotient normalizers. A complete rooted-tree realization is executable. | Recovery reference/channel/domain, decoder, protected data, completeness and local diamonds must hold for the chosen realization. Atomic commits or arbitrary local mismatch reduction alone do not imply confluence. [Consensus](../paper/reality_as_consensus_protocol.tex), [finite implementation](../code/consensus/README.md), [Lean primitives](../Lean/ObserverPatchHolography/Primitives.lean). |
| C04 | **Native scalar repair implementation.** Real seam means, conservative integer unit transfers, and nearest-integer pair completion are specified; large twelve-port federation runs retain independent checks. | These are distinguishable laws within a declared scalar realization and gluing. The scalar engine alone is not the full observable algebra/current/quantum model. The large runs establish their retained settlement results, not physical identification. [Federation evidence](../evidence/exact_federation_L6_canonical_20260909/README.md), [tower evidence](../evidence/exact_federation_tower_20260924/README.md). |
| C05 | **Activity after settlement.** Fair nearest-integer continuation has a unique uniform fixed-occupancy equilibrium, exact density covariance, occupation-history limit and local attempt-record corrections on a fixed connected graph. | Rates, sector, clock convention, record gain/memory and preparation are declared. The stochastic limit theorem and finite Lean witnesses have different scopes. A joint physical refinement limit is not supplied by fixed-graph stationarity. [Native record theorem](../paper/tex_fragments/NATIVE_REPAIR_RECORD_SPECTRUM.tex), [verification interface](../code/observer_dynamics/README.md). |
| C06 | **Distinct operational observers and continuation.** Finite witnesses have unequal local records agreeing on proper overlaps, record-conditioned control and checkpoint behavior; a nonzero private generator can fix public records. | Agreement equates shared restrictions, not all private states. The cited proper-meet witness is at one regulator; common cross-observer refinement and higher overlaps need additional evidence. Gauge-hidden labels alone are not individuality. Stable marginal distributions are not stable stored records: the current law-space criterion compares the two-time joint record, under a continuation map. [Consensus](../paper/reality_as_consensus_protocol.tex), [operational witness](../Lean/QFT/OperationalOverlapEvidence.lean). |
| C07 | **Authenticated history.** Writer/read-value certificates generate informational precedence and canonical longest-parent-chain height. An abstract finite-poset compiler also exists. | The producer must supply genuine, complete semantic dependencies. Compiler expressivity is not selection of a physical history, and identifying informational ancestry with physical causality requires an additional map. [Consensus](../paper/reality_as_consensus_protocol.tex). |
| C08 | **A finite carrier-to-support construction and refinement theory.** One artifact has twelve carrier charts, thirty seam algebras, twenty nontrivial triple restrictions, observer controls, confluent seam repair and an oriented sphere support limit. Repair-morphism and inverse-limit theorems specify normalizer compatibility. | This bridge is not absent. Its identification with an arbitrary native federation, modular state tower and physical field family needs explicit maps. Each proposed refinement must satisfy the theorem's move-word/coherence conditions or controlled vanishing defects. [Shared construction](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex), [consensus](../paper/reality_as_consensus_protocol.tex). |

### Geometry, quantum fields and gravity

| ID | Constructed | Premises and boundary |
| --- | --- | --- |
| C09 | **Finite Lorentz carrier and flat causal/count limit.** The source rank-three space plus a real axis gives Lorentz inertia `(1,3)`. Prepared conservative source-record populations with specified neighbour reads have a proved flat causal-order and count-volume limit. | Ambient dimension is constructed; arbitrary event-poset dimension is not thereby measured. Population, read law, duration and density choices are explicit. Curved/common interacting physical refinement remains separate. [Spacetime paper](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex), [flagship](../flagship/from_observer_consensus_to_standard_physics.tex). |
| C10 | **Clock constructions.** Fourth-root interval-count ratios recover proper-duration ratios in the flat source family. Separate operational matter-clock constructions and a regular-path Jacobi action clock exist. | These are real constructions, not a missing clock theory. Their agreement with each other and with a running laboratory clock is an additional physical identification. Ordinal height is explicitly not proper time. [Flagship](../flagship/from_observer_consensus_to_standard_physics.tex), [operational clocks](research/MASSIVE_OPERATIONAL_CLOCKS.md). |
| C11 | **Finite quantum and entropy identities.** Born–Lüders identities, declared two-wing bounds, central records, edge-center entropy splitting, modular first laws and conditional repair-channel results. | The algebra/state, representation and instrument matter. An effect table alone does not select the Lüders instrument. Conditioning requires positive event weight and preserves intra-block state; only the appropriate partition-averaged conditional state becomes the normalized block projector. A finite identity does not construct the common quantum scaling limit. [Microphysics](../paper/screen_microphysics_and_observer_synchronization.tex), [quantum repair](../code/quantum_information/ALGEBRAS_AND_REPAIR.md). |
| C12 | **Scalar waves, detectors and feedback.** Source-address/metric scalar actions support classical and Fock waves, controlled detector limits and authenticated intervention/feedback examples. | Action, preparation, boundary, controls, model time and quantization are supplied where stated. These consumers must be joined to the native repair source and physical instruments. [Flagship](../flagship/from_observer_consensus_to_standard_physics.tex), [spacetime paper](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex). |
| C13 | **Seam fields and conditional Maxwell continuum.** Finite charge/holonomy/kinetic constructions and a controlled Maxwell limit exist for a declared reversible curl-pair evolution and spatial domain. | The scalar Markov repair mean does not itself select reversible Maxwell dynamics. Currents, operational space/time, physical units and common refinement need a source attachment. [Microphysics](../paper/screen_microphysics_and_observer_synchronization.tex), [Maxwell module](../Lean/Screen/SeamMaxwellContinuum.lean). |
| C14 | **Interacting field construction.** A declared charged-scalar/Maxwell Whitney action has controlled classical sectors and a gauge-reduced, essentially self-adjoint finite-mesh quantum Hamiltonian with neutral unitary evolution. | Full classical backreaction estimates, interacting quantum continuum, selection/attachment of quantization and common physical source are separate. These finite constructions predate the local tensor-model search. [Flagship](../flagship/from_observer_consensus_to_standard_physics.tex), [observer paper](../paper/observers_are_all_you_need.tex). |
| C15 | **Gauge type and conditional matter structure.** Complete reversible response and endogenous overlap transport force the abstract local `u(1) ⊕ su(2) ⊕ su(3)` type. Declared matter tensors support anomaly/descent and common Z6-kernel results. | The complete response/transport contract is essential. Matrix-current realization, matter action, physical global form, family/seam attachment, laboratory current, scalar multiplicity and poles are not obtained by matching the Lie-algebra name. [Microphysics](../paper/screen_microphysics_and_observer_synchronization.tex), [SM correspondence](registers/SM_LAGRANGIAN_CORRESPONDENCE.md). |
| C16 | **Conditional Einstein reconstruction.** Exact null tomography, edge-entropy and small-ball coefficient results compose with directional balance and independently conserved stress to give Einstein form with constant Λ. A separate source-carrier inverse-square shell law exists. | One family must realize physical geometry, normalized modular charges, stress/Ward identity, area normalization, generalized-entropy stationarity, controlled shrinking-region remainders, universal coupling and scale. Shell flux assumptions and their join to Einstein gravity remain explicit. [Flagship](../flagship/from_observer_consensus_to_standard_physics.tex), [gravity premise ladder](registers/GRAVITY_PREMISE_LADDER.md). |

### Parameters, cosmology and empirical evidence

| ID | Constructed | Premises and boundary |
| --- | --- | --- |
| C17 | **P coordinates and conditional closure calculations.** Cell-area/capacity conventions, `ΔP=P−φ`, the declared relation `α_out=ΔP/√π`, and certified roots for specified endpoint maps. | Distinguish source, inner-coupling and measured-endpoint maps. Measured-α-conditioned P is not an independent α prediction. The registry's `p_detuning_same_quantity_readback` remains undischarged; no derived identification with a repair coin, temperature or physical rate exists here. [P derivation](../code/P_derivation/README.md), [physical identifications](../claims/physical_identification_registry.json). |
| C18 | **Conditional cosmological constructions.** Flat and specified expanding source-record/count representations, conditional primordial and radial reconstruction results, and dark-source scaling relations. | A supplied expansion profile is not derived cosmological dynamics. Finite expanding diamonds do not generally give exact proper-time ratios; redshift/count identities have specific epoch/diamond conditions. Primordial source, physical stress, radial lift, transfer/recombination, abundance and lensing attachments require their stated additional inputs. [Cosmology](../cosmology/README.md), [FLRW theorem](../paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex). |
| C19 | **Reproducible conditional comparisons and hardware experiments.** The repository contains numerical postdictions, frozen-register machinery, apparatus/controller evidence and bounded independent verification. | Numerical reproduction, programmed algorithm behavior, engineering resonance and a universal physical-constant measurement are different claims. Alex's P reference identified in this thread was a bench-test plan; no independent P-measurement receipt was located in the material inspected. This does not establish that no later experiment exists. Apply the public evidence rule to any promoted hardware claim. [Hardware policy](policies/HARDWARE_EVIDENCE_BUNDLE_H.md), [postdictions](POSTDICTION_LEDGER.md), [frozen register](../claims/frozen_prediction_register.json). |

## 4. The actual simple rule we should start from

Keep the native twelve-port seam graph and its scalar loads `x_i`. These
operations exist and must not be conflated:

| Operation | Update on a seam `(i,j)` | Exact property |
| --- | --- | --- |
| Real mean | `(a,b) → ((a+b)/2,(a+b)/2)` | Preserves total; lowers `V=Σx_i²` by `(a−b)²/2`. Uniform sampling of the thirty internal seams has mean operator `I−L/60`. |
| Integer unit transfer | `(a,b) → (a−1,b+1)` when `a−b≥2` | Preserves integer total; lowers V by `2(a−b−1)`. |
| Integer nearest agreement | Replace by `floor((a+b)/2), ceil((a+b)/2)`; declared orientation assigns the ceiling | Composes unit descents for `|a−b|≥2`; for difference one, swaps or waits; for zero, waits. A fair tie supplies the documented stationary continuation. |

For nearest agreement the potential drop is
`(d²−(d mod 2))/2`, with `d=|a−b|`. The full seam-disagreement sum
`Φ=Σ_seams(x_i−x_j)²` is **not** this Lyapunov function and can increase
while V decreases. Temporary local disagreement therefore does not by itself
require a new fundamental noise mechanism.

In a balanced fixed-occupancy sector on a connected graph, write
`n_i∈{0,1}`, `Σn_i=r`, `z=n−(r/p)1`, with `0<r<p`. Independent Poisson
edge attempts with positive declared rates and swap probability one half give

```text
Π = I − 11ᵀ/p
κ = r(p−r)/(p(p−1))
Cov(z(t),z(0)) = κ exp(−|t|L/2) Π.
```

This stationary activity is proved for the stated law. The terminal
component-multiset quotient forgets position, whereas local load histories
retain it. A nonconstant `z_i` cannot factor through a quotient constant on
its whole occupancy sector. A physical use of those histories must specify
the richer record/observer quotient and maps. That is an interface to
construct, not a reason to discard the stationary theorem.

The flagship records a second, distinct quotient boundary: the nonlinear
integer repair does not descend through the signed response quotient even
though its conditional mean does. The working state and ordered history
therefore require their stated additional data. A response quotient alone
cannot be substituted for the whole microscopic dynamics.

Reference implementations and controls are retained in the
[archived native source](../evidence/observer_dynamics_20260925/oph-physics-sim/)
and the [observer-dynamics interface](../code/observer_dynamics/README.md).
The scalar engine does not consume P as a jump probability. Its use as the
baseline is a disclosed branch choice, not a proof that A1–A3 uniquely select
every detail of that algorithm; the
[canonical repair-law RFC](CANONICAL_REPAIR_LAW_RFC.md) states additional
selection premises.

For the affirmative hardware/controller evidence in C19, the
[IBM quantum archive](../code/ibm_quantum_cloud/README.md) supplies repository
receipts and its declared instrument/null scope. It does not measure P.

## 5. Missing connections, with completion criteria

These are dependencies for the claims that use them, not tasks to execute
indiscriminately. A bounded physical proposal may address a subset; a claim
to a common physical model must supply every connection it consumes.

| ID | Missing or incomplete item | Concrete completion criterion |
| --- | --- | --- |
| M01 | **One source-bound candidate specification** | Pin one state space, carrier/federation, algebras, preparations, proposals, commits, equilibrium/quantum law, parameters, records and refinement. Give an explicit restriction or intertwining map for every imported component. |
| M02 | **The selected source's observer/quotient contract** | Show that actual public readouts, protected boundary, observer restrictions, histories and checkpoints factor through the declared quotient; test proper overlaps and interventions. Instantiate the acceptance, completeness and confluence assumptions actually used. Do not substitute a coarse multiset for a history-visible state. |
| M03 | **Operational geometry from that evolving source** | Construct an observable metric/collar/volume or probe-response readout and its state dependence, quotient invariance and refinement. The native load repair keeps its supplied graph/mean operator fixed; renaming load noise “curvature” or setting volume to an exponential of load does not establish this. [Existing obstruction](../paper/tex_fragments/NATIVE_GEOMETRIC_SOURCE_IDENTIFIABILITY.tex). |
| M04 | **A common physical clock and causal interpretation** | Identify which events and influences are physical; prove correspondence among authenticated records, geometric count clock, field/action time and running instruments. Declare the rate/scale calibration and compare moving as well as resting clocks. Accepted-step counts alone do not discharge this. |
| M05 | **Source-bound fields, currents and backreaction** | Realize the field variables and complete response/current contract on the same source; derive or explicitly supply its action. Test stable propagation, conserved currents/Gauss relations, geometry-to-matter and matter-to-geometry response when claimed. Shared scheduling or equal parameter names is insufficient. |
| M06 | **Common quantum state and controlled regime** | Supply the algebra, state, physical sectors and positive/unitary dynamics together. For a continuum claim, control their common refinement, relevant correlators and renormalization; distinguish classical sectors, fixed-background free limits and the full interacting limit. |
| M07 | **Physical entropy/area stationarity** | On the same regional state, define physical entropy, stress, area and fixed-volume variations; establish the area normalization and stationarity for an adequate family of independent perturbations, with controlled remainders. MaxEnt over a declared ensemble or stationarity of a Markov distribution alone is insufficient. |
| M08 | **Dynamical gravity on that family** | Supply the shared modular/geometric/stress inputs and all G1–G6 conditions in the flagship, or an explicitly different justified gravity route. Test conservation and gravitational constraints in the claimed regime. Merely introducing lapse/shear coefficients, or obtaining a common wave speed, does not complete gravity. |
| M09 | **Physical parameters and P readback** | Compute the claimed quantity from the chosen law using the same cell, metric, units and endpoint map; distinguish derived values from calibration. Prove any P-to-detuning/amplitude/rate identification. Do not set exploration probability `η=P−φ` by naming convention. |
| M10 | **Matter spectrum and physical interactions** | Supply the actual matter realization, chirality/family/scalar content and stable poles relevant to the proposed test. Running claims need charged response, thresholds and matching. A declared charge table or an abstract gauge Lie algebra does not by itself deliver these. |
| M11 | **Cosmological or dark-sector attachment, if used** | Derive or explicitly specify common expansion/source dynamics, stress, physical scales and observational transfer. Provide abundance, pressure/slip and lensing for dark-source claims. Preserve distinctions between conditional relations, fitted comparisons and predictions. |
| M12 | **One unavoidable observable and eligible test** | Enumerate all freedom allowed in the fixed candidate, prove a measurable relation/bound that survives it, attempt an admissible countermodel and obtain an independent check. Specify calibration, target exposure, baseline, justified uncertainty, sensitivity and rejection rule; freeze before eligible held-out/future comparison. |

For M07–M08, entropy stationarity is a substantive hypothesis about physical
regional variations. The Einstein implication is present. The work
is to realize its premises together. No claim here establishes that OPH must
add elementary gravitons, a particular kinetic term, or any specific new
microscopic geometry variable.

For noisy consensus in M02, the corrected theorem bounds
block-indexed expected distance under its contraction certificate. The
stated calendar-time bound additionally needs deterministic block endpoints;
random stopping blocks need separate occupation-time control. Block type is chosen at block start; its complete implemented conditional
law, including any within-block adaptation, must equal the certified kernel. A favorable
marginal kernel or stationary histogram is not sufficient.

## 6. Research approaches and minimal repair results

The research has tested several approaches in parallel. The following
summarizes analytical and numerical findings from the unpublished campaign;
it does not add a public verification bundle or an empirical confirmation.

### Approaches tried

| Approach | Result and limitation |
| --- | --- |
| Supplied quantum, Hamiltonian and covariant field models | Binding, energy, propagation and matter–geometry calculations provide consistency benchmarks. An imported Einstein–scalar–Maxwell action does not derive that action from patch repair. |
| Adaptive equilibrium repair and induced interactions | Specified local laws retain charged matter and induce gauge interactions in controlled settings. Naive geometry extensions can pin volume, and some global constraints fail quantum-transfer positivity. Results from different laws cannot be combined without a common construction. |
| Coupled geometry and quantum dynamics | Positive finite Hamiltonians with reciprocal matter/geometry response and a common matter/Maxwell principal metric were constructed. A frozen-background free-scalar continuum and finite joint entropy response were established. Geometry kinetics, couplings and physical time remain supplied. |
| Physical and entropy attachments | Tested direct attachments failed atomic-clock or weak-field comparisons. A one-sided area proxy failed a common stationarity coefficient across two sourced ground states at fixed mean volume to first order. These failures constrain those specified attachments; they do not exclude every physical area map or refute OPH. |
| Minimal local record and phase rules | Exact finite models isolate agreement, frustration, exploration, activity and equilibrium. They clarify which mechanisms suffice in their stated domains, but have no demonstrated map preserving the full twelve-port OPH contract. |

The positive tensor-capable Hamiltonian is the strongest benchmark in this
campaign for quantum positivity and coupled entropy/geometry diagnostics.
The documented native twelve-port construction remains the source-aligned
starting point. An icosahedral mean-repair calculation did reproduce its
scalar operator `I−L/60`; the added inertial laws and general Hamiltonian
models were separate supplied constructions.

### Results for minimal equilibrium detuning and repair laws

Here “detuning” must distinguish three quantities: a physical departure from
a preferred relation, disagreement between copies of a shared record, and the
probability of making an exploratory update. The tests do not identify these
quantities with each other or with OPH P.

| Question | What was shown | Boundary |
| --- | --- | --- |
| Can agreement preserve nontrivial structure? | Endpoint records can agree while phase relations retain nonzero loop holonomy and unequal residuals. | Agreement concerns shared data. Frustration can be static; it does not by itself imply motion. |
| Can purely greedy repair fail? | In the restricted single-patch Z3 cycle rule, random greedy ties leave single-error traps. A temporary second error permits an escape path; an explicit construction extends to finite odd cycles. | This is an obstruction of that move menu. Native OPH uses different admissible repairs and scores; some local mismatches can increase even during valid native descent. |
| Can arbitrarily rare departures restore convergence? | Allowing every increment on an inconsistent adjacent pair gives almost-sure agreement on finite odd Z3 cycles. Mixing this rule into greedy repair with any positive exploration probability gives the same almost-sure agreement, while the mean wait for the first uphill event from a single-error trap grows inversely with that probability. | Even cycles have an additional conserved alternating sum. Rare updates do not guarantee efficient repair, small amplitudes or a uniform bound on total disagreement. No identification of the exploration probability with `P−φ` was derived. |
| Can motion survive agreement? | A compact single-patch rule combines record repair with residual swaps; the tested finite joint model has only agreed recurrent states, with continuing residual motion and multiple stationary sectors. | Its component marginals are autonomous. Sharing one instruction does not establish reciprocal coupling, a unique equilibrium or physical time. Native OPH's neutral swaps and stationary record correlations provide the source-specific counterpart in section 4. |
| Does a small thermal departure solve both preparation and equilibrium? | A common heat-bath weight has an exact equilibrium, but its conditional choices factor. Finite temperature recreates record errors. On the triangle, the cold limiting measure favors agreement while the strictly cold dynamics can trap a single record error. | The long-time and zero-temperature limits need not commute. A desirable equilibrium does not prove that a chosen repair process reaches it. |
| Does residual transport produce physical propagation? | A one-defect sector has an exact local hopping law and a conditional diffusive continuum under a declared event-rate scaling. Reversible swaps also give motion without a thermal bath. | Diffusion is not a relativistic wave equation. Conservation, stationary activity, positive Euclidean transfer and unitary quantum dynamics require distinct checks. |

These results support a focused possibility: local agreement repair can
coexist with persistent distinctions and activity, and occasional exploratory
departures can overcome particular finite algorithmic traps. They do not
establish that fundamental OPH repair requires such departures, select their
physical size or rate, or prove a model of the Universe. The construction rule
is to test the simplest source-compatible version of each mechanism and
preserve the exact distinction between a proposal, an accepted repair and
equilibrium evolution.

### Dongyang Stephen Chen's contributions

Both **Genesis / The Ladder of Symmetry** and **The Balance Bridge** were
supplied by Dongyang Stephen Chen. Genesis motivated minimal local relational
records and exact reconciliation; the Balance Bridge motivated loop-sensitive
tests, symmetry comparisons and checking dynamics against null controls.
The particular local rules and proofs in this campaign were added constructions,
not claimed outputs of unseen Genesis/CycleHologram or Prism code.

The supplied Genesis moment formula has a useful bounded-support result:
with distinct nonzero field nodes and nonzero weights, `2k+1` moments separate
support at most k from support exactly k+1. Universal rejection of arbitrary
overflow does not follow; finite-field collisions were retained. Likewise,
shared algebraic form, a unit mode, phase holonomy or a reported power law does
not by itself select physical dynamics or establish criticality. The supplied
Prism measurements and software-performance claims were not independently
reproduced. These inputs remain useful heuristics and component ideas.

## 7. Concrete next decision and issue completion

The next source-alignment check is **M01–M04 on a bounded native instance**:
export the actual state, transitions, protected data, record restrictions and
readout maps; keep strict descents, swaps and waits distinct; determine exactly
which proposed physical readout survives the declared quotient. Then test its
connection to the existing source geometry/count clock and field action. A
constant fixed-geometry readout is a useful negative control, not a reason
to introduce an arbitrary metric field immediately.

Before a new numerical run, record:

```text
SOURCE:       pinned carrier, state, law, records and refinement
EXISTING:     construction being reused, with proof/evidence references
MISSING:      one named connection from M01–M12
MINIMAL MOVE: readout/composition first; otherwise one justified addition
PRESERVATION: totals, descent/neutral classes, observer/readout identities
TEST:         explicit success condition, counterexample and baseline
INPUTS:       supplied, derived, calibrated and target-exposed quantities
CAP:          bounded computation and decision/stop condition
CLAIM:        exact finite / conditional continuum / physical proposal / tested
```

This check does not replace the actual acceptance criteria of the two physics
issues, whose bodies were read for this review:

- **[#1025: choose one test](https://github.com/FloatingPragma/observer-patch-holography/issues/1025):**
  audit at most three existing physical routes (matter/current relation,
  shared-calibration redshift/deflection, or attached dispersion) and select
  at most one law, observable and bounded calculation. End with a feasible
  specified candidate or **no testable physical model established**. Its
  recorded review date is 9 October; this ledger does not claim that decision
  has been completed.
- **[#1026: derive and test the selected model](https://github.com/FloatingPragma/observer-patch-holography/issues/1026):**
  proceed only after the positive selection required by #1025. Fix the model,
  derive an unavoidable measurable relation, attempt an admissible
  countermodel, independently check it, and state a feasible calibrated
  rejection test. The issue limits the attempt budget and calls for final
  review by 18 October. A testable proposal awaiting data is unconfirmed;
  a justified failed comparison rejects the specified candidate; no viable
  candidate is a legitimate negative outcome.
- **[#1033: audit existing evidence](https://github.com/FloatingPragma/observer-patch-holography/issues/1033):**
  repair reproduced defects with independent controls and downstream impact
  review. Audit completion does not discharge #1025 or #1026.

All three issues were open when read on 10 October. Existing cosmology
comparisons with seen data remain postdictions, fits or diagnostic baselines
according to their ledgers; unarmed or conditionally registered tests retain
their eligibility requirements. Mathematical construction, a testable
physical proposal, and empirical support are separate deliverables.

In particular, FZ-13/FZ-15 retain their conditional registration and anchoring
requirements, and FZ-14 has no evaluated event likelihood in the reviewed
register. These are not completed prospective confirmations. The
[cosmology postdiction ledger](../code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md)
classifies its seen-data comparisons separately.
