# Observer-Patch Holography: Lean 4 Formalisation

Lean 4 / Mathlib formalisation effort for *Paradise as Fixed-Point Consensus*
(B. Müller, 2026; source in `paper/paradise_as_fixed_point_consensus.tex`),
with **Proposition 4.2** as the primary target.

## Scope

This project is an active Lean 4 / Mathlib formalisation and proof-audit
surface for the OPH consensus layer. Contents:

- Lake project with pinned `leanprover/lean4:v4.29.1` and Mathlib `v4.29.1`.
- An abstract-rewriting skeleton (Newman's lemma, normal-form uniqueness,
  descent termination, fixed-point zero-potential corollary) in
  `ObserverPatchHolography/AbstractRewriting.lean`.
- A concrete OPH carrier layer in `ObserverPatchHolography/Primitives.lean`
  for records, patch/interface data, observable overlap maps, the mismatch
  potential `Φ`, gauge equivalence, edge-consistency, and non-vacuity witnesses.
- Machine-checked proof-audit material for the consensus/reconstruction layer,
  including boundary-fiber observer uniqueness, commutation-based confluence,
  concrete countermodels separating confluence from observer-facing uniqueness,
  and axiom audits for the discharged reconstruction statements.
- A sorry-free bridge-boundary module
  `ObserverPatchHolography/BridgeBoundaries.lean`: disjoint writes that
  separately preserve a nonlinear protected observable can fail under
  composition, and the even/odd three-bit laws are distinct despite identical
  one- and two-bit marginals. These theorems keep the local-diamond and
  union-collar coherence receipts explicit in the papers. See
  `BRIDGE_BOUNDARY_INDEX.md`.
- Constructed, admission-free asynchronous repair machinery:
  `localRepair`, `Repair`, and `repair_respects_gauge`, with descent,
  termination, reachability, and normal-form receipts. Confluence and
  completeness remain explicit carrier properties and are not universal.
- A sorry-free scalar completed-repair theorem in
  `ObserverPatchHolography/ScalarSeamRepair.lean`. Seam support, endpoint
  agreement, and endpoint-sum preservation force pair averaging. Uniform
  averaging over a finite seam family is $I-L/(2|E|)$, hence $I-L/60$ for
  thirty seams. This proves the finite linear branch and does not identify it
  with the universal or physical repair law.
- A sorry-free finite repair-word information-projection theorem in
  `ObserverPatchHolography/RepairWordSchedule.lean`. On a complete,
  quotient-deduplicated finite event alphabet, KL minimization against source
  counting measure over the full word simplex forces the uniform product law
  at every length. The module proves exact declared-subset masses, and hence
  the masses of any externally supplied orbit partition, plus prefix
  consistency without a geometric transitivity premise. Completeness of the
  alphabet, its orbit classification, source counting, and temporal freedom
  remain proposed A1-R/A2-R inputs rather than conclusions of the current
  basis.
- A sorry-free V2 A1 ontology ledger in `Time/TimeOrderLedger.lean`.
  Universe closure, repair execution order, observer record order, modular
  parameter, worldline realization, clock readout, proper-time interval, and
  the optional global-time branch have separate types, and all 56 ordered
  transitive coercions between distinct layers fail. Explicit interpretations
  use `NamedRealizationMap` without an implicit function coercion.
  The module also proves positive affine regraduation preserves strict clock
  monotonicity and that record order does not fix a clock origin. It supplies
  no worldline, physical calibration, global time function, or modular-time
  identification. See `A1_TIME_ORDER_LEDGER.md`.
- A sorry-free bounded observer-time package in
  `Time/ObserverHistory.lean`, `Time/ClockReadout.lean`,
  `Time/WorldlineRealization.lean`, `Time/ProperTimeCalibration.lean`, and
  `Time/ClockComparison.lean`. Bare order-compatible scalar readouts admit arbitrary strictly
  increasing regradings, and a cubic three-tick control is provably
  non-affine. Under a supplied affine unit-timelike event law the clock
  increment along that same history is additive and squares to the
  chart-invariant Lorentz quadratic interval. Two events determine a unique
  positive-affine interpolation of four supplied readings; a distinct third
  point gives a nondegenerate no-new-fit-parameter check. A held-out reading
  needs separate predesignation and custody. A nontrivial affine algebraic
  control passes the check and the cubic control fails it. A separate
  three-record, two-chart witness inhabits the geometric interface.
  No source history, physical instrument, three-clock network, refinement
  limit, or SI unit is constructed. See `D1_OBSERVER_TIME_CALIBRATION.md`.
- A sorry-free B1 public-record package in
  `EventAlgebra/PublicRecordAlgebra.lean` and
  `EventAlgebra/NoBroadcastingAdapter.lean`. The span of a projective
  partition is a commutative star subalgebra exactly equivalent to complex
  functions on its nonzero labels. A common linear isometry that sharply
  copies two distinct states from one normalized blank forces them to be
  orthogonal. The mixed-state no-broadcasting implication is an explicit adapter
  premise, not a theorem of the package. See `B1_PUBLIC_RECORD_ALGEBRA.md`.
- A sorry-free A3 root interface in `Tower/ConsensusTower.lean`. One directed
  tower packages finite observer and record fibers, A1 record orders, private
  matrix algebras, commutative public subalgebras, selected states and
  generators, and explicit functorial refinement laws. State compatibility is
  the contravariant trace-pairing restriction equation. A constant adaptor
  packages an existing projective partition and certified state with discrete
  order and zero generator. It constructs no nonconstant source tower,
  public-world endpoint, causal net, geometry, clock, or continuum limit. See
  `A3_CONSENSUS_TOWER.md`.
- A sorry-free B2 finite publicization package in
  `EventAlgebra/PartitionAverageCP.lean`,
  `EventAlgebra/TwoScalePublicRepair.lean`,
  `Thermodynamics/PoissonizedRepair.lean`,
  `Thermodynamics/PoissonizedRepairOperatorExp.lean`, and
  `Dynamics/ConditionalExpectationGenerator.lean`. It proves normalized
  averaging Kraus data, trace preservation, exact residual relaxation and
  semigroup laws, a generator flow equation, the literal bounded-operator
  exponential formula, the projector-rate matrix dissipator identity, and the
  stable commutant. The formal CP/CPTP predicate, multi-collar no-cancellation
  proof, source-derived rate, clock, and physical channel are not supplied. See
  `B2_PUBLICIZATION_DYNAMICS.md`.
- A sorry-free bounded B3 public/private dynamics package in
  `Dynamics/PublicMarkov.lean` and `Dynamics/PrivateInner.lean`. Positive
  unital active-record maps are exactly row-stochastic kernels, continuous
  label-permutation flows are trivial, every automorphism of one finite full
  private matrix block is unitarily inner, and a fixed self-adjoint Hamiltonian
  gives a unitary real-parameter von Neumann flow. Classification of arbitrary
  public star automorphisms as label permutations, arbitrary central-block
  decomposition, and the converse continuous-group-to-Hamiltonian theorem are
  open in issue #679. See `B3_PUBLIC_PRIVATE_DYNAMICS.md`.
- A sorry-free integer completion bridge in
  `ObserverPatchHolography/DirectedSeamRepair.lean`. It proves total-sector
  preservation, the exact parity residual, equality of the two-direction
  rational mean with pair averaging, and strict quadratic descent for a
  nontrivial incoming mismatch. Its companion
  `ObserverPatchHolography/DirectedSeamRepairProgress.lean` proves the
  concrete carrier path bound, the exact minimum-shell characterization for
  arbitrary signed fixed-total loads, local nonincrease, neutral transport,
  and strict descent along one-, two-, and three-seam progress words.
  Schedule selection and the probabilistic hitting statement remain separate.
- A sorry-free **#304 boundary-fiber carrier witness** in
  `ObserverPatchHolography/Rule90.lean` (PR #385): the linear Rule 90 carrier
  discharges the `Hfib` binder of `boundary_fiber_observer_unique` on a proper
  information-set boundary, with a bad-boundary counterexample, a non-trivial
  gauge, and a local-repair no-go (`H1`–`H3` route only). A carrier-level
  witness; it does **not** advance the Prop 4.2 target. See `PROOF_INDEX.md`.
- A sorry-free **Part-A coupling-algebra layer** (13 lemmas, standard axioms
  only): `ObserverPatchHolography/BridgeEquivalence.lean` (bridge
  count/tick equivalence, 5 lemmas),
  `ObserverPatchHolography/CapacityFixedPoint.lean` (capacity
  fixed-point uniqueness schema, 4 lemmas), and
  `ObserverPatchHolography/SeedPi.lean` (CAP-P seed statement,
  4 lemmas). These formalise the ALGEBRAIC layer of the coupling theorem and
  carry no physical-derivation content; the physical identities I1/I2 are
  outside the formalised set. Numeric interval enclosures stay in the Python
  certificates (`code/capacity_readback/`, `code/P_derivation/`); no
  floating-point numerics enter Lean. See `PROOF_INDEX.md`.
- A sorry-free **#578 corrected Einstein-branch kernel** under
  `ObserverPatchHolography/EinsteinBranch/`: bare finite-tower
  non-entailment, typed common-domain and boundary-fibre composition,
  explicit nine-direction null tomography, finite bulk/edge entropy and
  MaxEnt algebra, exact small-ball coefficient arithmetic, timelike/null
  tensor upgrades, Ward/Bianchi constancy, and strict manifest deletion
  logic. Continuum, asymptotic, UC, VR, scale, and tower-nonemptiness inputs
  stay explicit in the main theorem type. See `EINSTEIN_BRANCH_INDEX.md`.
- A standalone, application-neutral proof package in
  `ObservableNormalForms/`.  Its generic endpoint theorem is connected
  to the concrete local-repair interface by
  `ObserverPatchHolography/Bridges/ObservableNormalForms.lean`.  This
  bridge characterizes the #304 application premise for arbitrary carriers;
  the premise itself is discharged on the declared domain by the
  `BoundaryFiber` module below.
- A sorry-free **#304 application theorem**
  (`ObserverPatchHolography/BoundaryFiber.lean`): the verified
  rooted-tree packet-net domain (*Reality* Definition `def:tree-packet-domain`)
  as an `OPHCarrier` family (`TreePacketNet`), the declared physical
  boundary/sector map `B_OPH` as the concrete root-packet readback
  (`TreePacketNet.BOPH`), and the gate theorem
  `BOPH_injective_modulo_gauge`: consistent records with equal `B_OPH`
  readback are `gaugeEquiv`, for every net in the class. The generic
  hypotheses (`ObservationPreserving`, `CompleteFor`) are discharged for the
  domain's own tree repair, giving the endpoint payoff
  `BOPH_observerEndpointUnique` with no confluence input. Witnesses: the
  paper's exported four-vertex `ℤ₃`/`ℤ₂` instance with a proper gauge
  quotient, and a deficient-readback failure. Outside the declared class,
  boundary identifiability stays a named per-net premise (countermodels:
  `demoCarrier_Hfib_fails`, `rule90_Hfib_bad_fails`). See
  `BOUNDARY_FIBER_APPLICATION.md` for the hypothesis-by-hypothesis bridge.
- A sorry-free finite **Physical A5-to-SM non-identifiability boundary** in
  `Screen/PhysicalA5ForcingNoGo.lean`: Euler total charge does not reconstruct
  a unique defect profile, the twelve-port source reduct does not reconstruct
  a unique current completion, and the compact-current reduct does not exclude
  a sterile matter completion. These are finite no-go theorems about the
  exposed reduct, not a claim that a richer operational producer packet cannot
  select the physical completion.
- A sorry-free conditional **B14 oriented-face discriminator** in
  `Screen/OrientedFaceBracketSelector.lean`: the twenty declared oriented
  faces give exactly `B_face = 60 R13`, and an exact coefficient witnesses
  Jacobi failure. On 792 upper-triangular structure-constant coordinates, an
  independently replayed exact primal-dual certificate finds the same compact
  family `G` nearest under total-absolute, Hilbert--Schmidt-squared, and
  worst-coordinate edits. Lean checks the serialized radical values and all
  three strict orderings in `three_norm_unique_nearest_G`. The `F` total-
  absolute value is an unattained infimum. Every norm and minimum-repair rule
  is a basis-dependent supplied premise, and the comparison covers only the
  classified compact locus. Issues #705 and #697 track source bracket
  selection and current realization. See `B14_ORIENTED_FACE_SELECTOR.md`.
- A sorry-free **B14 invariant-metric phase diagram** in
  `Screen/OrientedFaceInvariantMetric.lean` over the independently replayed
  `invariant_metric_phase` certificate: the commutant of the port action is
  four dimensional and spanned by the symmetric spectral projectors, so the
  positive sector-scale cone is the complete family of invariant carrier
  inner products; the induced bracket metric is channel-diagonal, and the
  certified closed forms give quantified real theorems: `P` is excluded for
  every invariant metric, every sector-balanced metric and every metric with
  `beta/delta` in `[1/50, 6]` selects `G` uniquely, `F` wins in a nonempty
  region across the exact three-scale tie surface with witness `(8,1,1)`, and Galois
  conjugation with the sector swap exchanges the two mirror families. The
  nearest-point rule and carrier-induced class remain declared; nothing here
  is source selection. See `B14_ORIENTED_FACE_SELECTOR.md`.
- A sorry-free abstract **A2 holonomy bridge** in
  `Screen/A2HolonomyBridge.lean`: the projective implementer group, the
  supplied port-response flow, and its algebraically generated subgroup are
  separate objects. `CarrierHolonomyFullEndogenous` supplies one closed path
  for every proper carrier action and factors that path's projective
  implementer through the response-generated subgroup and the kernel of its
  response action. Response naturality is evaluated on the same path. Other
  closed paths need not be endogenous. The kernel is the formal
  centralizer/spectator factor and is definitionally silent on the response
  codomain. A nontrivial carrier response therefore requires a nontrivial
  response-generated action. An injective linear response on the twelve-port
  space has a kernel-checked twelve-dimensional image.
  The assembled theorem reaches centre dimension one and factor dimensions
  `{3,8}` only through explicit adapter premises for the centre bound and the
  centreless four-factor fixed-space description. Construction of the actual
  projective unitary group and response flow, source extraction, the two
  Lie-theory adapters, reductive decomposition, and compact-simple
  classification remain explicit inputs.
- A sorry-free inverse-continuum boundary in
  `Screen/BipoSHInverseBoundary.lean`. Uniform finite repair obeys the exact
  operational identity $L=2|E|(I-R)$. Explicit two-coordinate and copy-space
  counterexamples show why finite positivity and scalar transfer cancellation
  do not supply a continuum inverse or remove radial multiplicity. Uniform
  coercivity, quotient control, and a physical intervention/readout remain
  separate premises.
- A sorry-free frame-quotient theorem in
  `Screen/BipoSHFrameInvariant.lean`. A coherent orthogonal or unitary frame
  isometry preserves the frozen unsquared rank-six norm statistic and its
  squared form. Nonzero multiplicity-one band amplitudes cancel in the
  squared statistic. The module does not construct or physicalize that
  isometric transfer, and it does not exclude copy-space mixing.
- A sorry-free conditional determinant kernel in
  `Screen/VolumeReadoutBridge.lean`. It proves collar shift cancellation,
  conformal determinant scaling, the positive volume-ratio implication
  $q=\zeta$, common positive-density-factor cancellation, calibrated
  low-mode projection, and a bare-collar counterexample. The quotient-visible
  uniform-density cut, positive spatial-metric attachment, physical rechart
  law, stress/adiabaticity, and cosmological transfer are open problems.
- A sorry-free conditional primitive-port translation bridge in
  `Screen/PrimitivePortTranslationBridge.lean`. On a declared continuous
  three-dimensional scalar field, it diagonalizes a paired twelve-port shift
  operator on plane waves and proves positivity, port-relabeling covariance,
  and exact passive-frame covariance. A supplied exact directional moment
  packet yields the frozen fourth- and sixth-order coefficients. Derivation of
  the frame and moment packet, canonical source selection, the positive length
  scale, physical time evolution, sector identification, frame and boost
  physics, nuisance isolation, and exclusivity remain outside the theorem.
  The abstract shifts do not define a locally finite spatial lattice.
- A sorry-free primitive-port scale boundary in
  `Screen/PrimitivePortScaleBoundary.lean`. It proves exact rescaling laws for
  the frozen correction coefficients, the fixed scale-free ratios, the
  conditional dimensionless metric relation, its positive-square-root form,
  and an algebraic counterfamily for the metric coefficient. The result does
  not construct source dynamics at arbitrary scale, determine the physical
  length, or select `kappa = 1`.
- A sorry-free primitive-port metric attachment boundary in
  `Screen/PrimitivePortMetricAttachment.lean`. A declared equal-area shell
  attachment gives the exact metric coefficient `orbitCard/(4*pi)`. The
  twelve-port, twenty-face, and thirty-edge orbit readings therefore give
  `3/pi`, `5/pi`, and `15/(2*pi)`, which are pairwise distinct. On the port
  attachment the module proves `a/ell = sqrt(3*P/pi)`, propagates any certified
  `P` interval to an exact squared-hop interval, and substitutes the metric
  relation into the frozen packet to give
  `C4 = -3*P*ell^2/(20*pi)`,
  `B0 = 3*P^2*ell^4/(280*pi^2)`, and
  `B6 = 2*P^2*ell^4/(875*pi^2)`. It does not identify a screen cell with any
  incidence orbit, identify the shell radius with the translation hop, or
  construct the required quotient-natural metric readout.
- A sorry-free primitive-port operator selection boundary in
  `Screen/PrimitivePortOperatorSelectionBoundary.lean`. It proves that A3's
  uniform state weights close the equal-weight operator only after an explicit
  state-to-first-hop attachment. On a complete direct-port action, the
  A1-derived antipodal frame implies reciprocal pairing exactly. The module
  also constructs a normalized nonnegative two-range family with the same
  port covariance and quadratic coefficient but a variable quartic
  coefficient. Port-state symmetry does not by itself select kinetic
  coefficients or first-hop exclusivity.
- A sorry-free complete-family cancellation theorem in
  `Screen/KineticFamilyCancellation.lean`. For the conditionally declared
  matter-trace kinetic ray `(10/3, 2, 2)`, each complete family's imported
  one-loop beta contribution is collinear with that ray, so `nG` cancels
  exactly from `det(x,k,b(nG,nH))`. The theorem derives the general
  scalar-count plane and its `nH = 1` integer specialization. It does not
  select the physical kinetic form or particle census, and it does not supply
  thresholds, higher loops, or a measured coupling attachment.
- A sorry-free strict W/Z pole-consumer scale quotient in
  `Screen/ElectroweakPoleScaleQuotient.lean`. It proves that the explicit
  common pole scale cancels from the ratio of the two factored strict complex
  poles at fixed normalized corrections. This covers passive unit rescaling;
  active vacuum-scale changes can move threshold ratios inside the normalized
  corrections. Its exact unconstrained factored-coordinate counterfamily
  displays the algebraic freedom; it is not an accepted consumer packet.
  Source-selected coupling and self-energy inputs remain necessary for a
  numerical ratio.
- A sorry-free conditional photon-lepton threshold packet in
  `Screen/SeamCurrentPhotonLeptonThreshold.lean`. It proves the exact global
  cosine-symbol upper bound in the explicit Euclidean carrier metric, the
  leading fixed-share photon-lepton coefficient combination, and the rank-one
  equal-share map with its two-dimensional fiber. Under separately stated
  physical-photon, additive-conservation, and Lorentz-invariant charged-lepton
  premises, photon decay into an electron-positron pair is kinematically excluded and the FZ incoming-energy
  domain is contained in the Lorentz-invariant photon domain. On the
  Lorentz-invariant charged-lepton branch, equal sharing uniquely maximizes
  the leading head-on, collinear residual. General independent-lepton and
  full anisotropic optimization are open. No interaction,
  opacity, source, shower, or detector statement follows.
- A sorry-free **finite event-algebra journal-neutral core** (eleven modules
  under `EventAlgebra/`, compiled by lake target `EventAlgebra`, 152 audited
  declarations, standard axioms
  only): events as Hermitian idempotents, states as positive trace-one
  matrices, Born weights (reality, nonnegativity, normalisation,
  additivity, complement bound, monotonicity), Lüders conditioning
  (state preservation, repeatability, idempotence, compatibility,
  classical restriction, fixed-point characterisation), the conditional
  expectation onto a commutative center (projector laws, state
  preservation, trace selfadjointness, Pythagoras, contractivity,
  uniqueness, compatibility with conditioning on central events), the
  expectation functional, supplied-state Robertson inequality, complete
  supplied-partition operational quotient, the exact positive uniform
  independent-sign random-unitary representation of pinching, its global-
  sign negative control, the totalized-log support countermodel, and the
  Tsirelson bound `‖S‖ ≤ 2√2` in abstract unital C*-rings with a matrix
  instantiation. The support countermodel rejects the naive totalized raw
  trace formula, not support-aware relative entropy. **This core is
  OPH-vocabulary-free by design** (namespace `EventAlgebra`, Mathlib-only
  imports, no repository vocabulary). The larger `EventAlgebra` umbrella also
  imports OPH-facing adapter modules which are outside this neutral core. The
  core is the journal-neutral surface for submission; every lemma is tagged
  **algebra-only** or **consumes a tracial state** in its doc comment. Inventory in `PROOF_INDEX.md`
  ("Finite event algebras"); Mathlib friction log in
  `EventAlgebra/MATHLIB_NOTES.md`. Not a Prop 4.2 / Def 4.1 item.
- An admission-free **B4 fixed-word locality package**:
  `ObserverPatchHolography/Locality/DependencyCone.lean` checks an
  `n`-move closed-neighborhood dependency upper bound against the concrete
  `localRepair`; `ObserverPatchHolography/Locality/NoSignalling.lean` proves
  generic finite marginal and partial-trace identities on a supplied
  bipartite split. B4's sole live gate is #692 (E1), which owns only finite
  coverage and OPH region factorization. A separate conditional E2 helper,
  `ObserverPatchHolography/Locality/AdaptiveScheduler.lean`, proves an
  adaptive cone/no-influence/refinement-naturality packet only for supplied
  `ConsultsOnly` scheduler and consultation-region data. It produces no source
  scheduler, state/channel semantics, physical distance, or clock. Those
  production and identification obligations continue under #693 (E2), the operational clock under
  #703, and continuum causal/time-slice and physical spacelike
  attachment under #700 (E3); these are downstream rather than B4 claim
  gates. See
  `B4_LOCALITY_BOUNDARY.md`.
- An admission-free **B5 finite conservation package**:
  `Screen/RegionalContinuity.lean` proves the exact regional and global
  balance laws with the seam boundary treated as net inflow;
  `Screen/DiscreteGauss.lean` proves neutral-load solvability and identifies
  the full solution fibre with an affine translate of the
  nineteen-dimensional cycle kernel; `Dynamics/ProtectedCharge.lean` gives
  the reusable dual-channel fixed-observable criterion. These are finite
  incidence and linear statements. Physical stress-energy realization and
  the continuum Ward identity are explicit inputs. See
  `B5_WARD_BRIDGE.md`.
- An admission-free but **incomplete B7 history/variation package**:
  `InformationProjection/PathGibbs.lean` proves conditional finite path-space
  information-projection, modal ordering, and full strict-nonminimal
  inverse-noise concentration. `Variational/DiscreteEulerLagrange.lean` and
  `Variational/DiscreteNoether.lean` prove separate scalar real-path local
  identities and a chain-wide constant scalar current when the premises hold
  at every interior record. `Variational/FiniteHistoryBridge.lean` proves that
  no finite real-path family is closed under all real single-site variations;
  `Variational/FiniteRealTransfer.lean` gives the receipt-gated composition
  that survives. `InformationProjection/LogTransitionAction.lean` derives the
  exact finite log-transition action, and `Variational/LegendreBridge.lean`
  proves its conditional mechanics interfaces. Finally,
  `Variational/RealizedHistoryLegendreNoGo.lean` proves that the same binary
  source law admits distinct strictly convex real enrichments and distinct
  Hamiltonians while the canonical bilinear extension has no global momentum
  solver. The source does not select the real curvature.
  `InformationProjection/ReferenceNormalForm.lean` fixes the reference gauge
  within a declared normal form: independent target-scrambling invariance at
  fixed source, row-constant weight, and constant step
  action each characterize the uniform kernel among row-stochastic kernels,
  so every invariant Markov reference has path law `stepUniformRef` exactly,
  with biased target-weight and simultaneous-relabeling controls; the latter
  proves ordinary simultaneous relabeling does not force uniformity. The
  stronger invariance principle is a
  normal-form convention, not a source product. Issue #683 tracks
  source selection of the reference principle, physical action, clock,
  amplitude, fields, continuum,
  and observable current. See `B7_HISTORY_BRIDGE.md`.
- A sorry-free **E1 rich-fibre regional net** in
  `QFT/RichFibreWitness.lean` and `QFT/RichFibreRegionalNet.lean`: the
  literal mirror of the second preregistered 64k run's payload (four
  greedy-disjoint observers, twenty-node windows, split-fibre counts
  3, 3, 4, 4, kernel-decided disjointness and census receipts) and the
  finite causal observer net over it, with a genuinely noncommutative
  block algebra and embedded two-by-two designated-fibre matrix corner at every
  observer region, elementwise commutation across distinct windows,
  character-and-block restrictions, computed four-window coverage of the
  top algebra, and the drop-one coverage negative control. The enrichment
  rule, region lattice, restriction system, state, and identity repair are
  declared readings of the custody-pinned payload. No regional tensor-factor
  receipt is proved, so E1 and the B4 attachment are unresolved; downstream
  CP/CPTP, scheduler, continuum, and clock scopes stay with E2, E3, and E5.
  The adapter preflight in `QFT/SourceOperatorGeneration.lean` and
  `QFT/JointSlotFactorisation.lean` adds a bounded post-hoc packet over
  the same retained bundle: the star algebra generated by each rich
  observer's counted transition operator and diagonal field projectors
  equals the full matrix algebra on its realized alphabet, a
  diagonal-only control proves the transition operator load-bearing, and
  on the declared joint slot space of one observer pair the lifted
  generated algebras are exactly the two tensor factors, the complete
  checkpoint-projector family pinches without changing the other slot's
  marginal, and the conditional expectation onto a slot carries Kraus,
  positivity, trace, fixed-slot, scalarisation, and idempotence
  receipts. The extraction is post-hoc and ineligible as validation, and
  the slot assembly is a declared postprocessor; issue #692 keeps the
  source-attached net contract.
  `QFT/CPRestrictionNet.lean` and `QFT/TwoSlotCPNetWitness.lean` execute
  the authorized restriction redesign: conditional expectations replace
  the star-homomorphic restrictions that a kernel-checked obstruction
  excludes above matrix factors, a joint-coverage law is part of the
  interface, and the designated pair inhabits the resulting diamond with
  Kraus-certified expectations, a proper local pinch channel fixing the
  disjoint slot, and coverage by the two slot regions. Tower anchoring
  stays with issue #692.
  See `E1_FINITE_CAUSAL_OBSERVER_NET.md`.
- An admission-free bounded **B13 effect and phase-lift boundary**:
  `EventAlgebra/FiniteEffectClosureBoundary.lean` proves that continuity and
  normalized antipodal binary contexts do not force Born affinity.
  `QFT/SourcePhaseLiftBridge.lean` proves that subtracting the normalized
  complex commutator of the source-attached algebraic pair from `I/2` gives exactly the
  missing positive Pauli-Y projector and completes fixed-trace two-by-two
  tomography, while even a
  generous closure under real coarse graining and real Kraus pullbacks remains
  Pauli-Y blind and cannot produce that effect. This is an algebraic target,
  not a source operation or outcome receipt.
  `QFT/ConjugationGauge.lean` proves that simultaneous entrywise conjugation
  of state and effect preserves real Born weights, while state-only blindness
  requires a conjugation-fixed effect; it fixes the complexified real web
  pointwise, and exchanges the two candidate completions, so the
  two-candidate ambiguity is one conjugation orbit and every
  conjugation-fixed frame conflates each state with its conjugate.
  `Thermodynamics/RepairCurrentOrientation.lean` kernel-decides a post-hoc
  repair-load raw-count diagnostic on retained data from a run with an
  unrelated B12 contract: the designated pair and cycle, their maximality,
  lexicographically least tie breaks, exact row-normalized selected-cycle
  check, and the reversal-odd
  orientation bit. `QFT/SourceOrientedCompletion.lean` transports the bit
  onto the completion torsor under a declared typed convention with a
  kernel-decided applicability condition; either torsor element completes
  state tomography, and the finite Busch--Gleason composition pins the
  represented state by three oriented weights. The transport is a declared
  convention. Neither the statistic nor designation rule was preregistered;
  the y-magnitude readout, operational additivity, and fresh prospective
  validation stay with issue #702. See
  `B13_EFFECT_CLOSURE_BOUNDARY.md`.
- An admission-free **B8 finite transport package**:
  `Thermodynamics/GreenKubo.lean` proves the reversible Dirichlet identity,
  Onsager symmetry, finite-matrix positive semidefiniteness, and an exact
  finite correlation sum with its Poisson remainder. It also proves that,
  under a strictly positive reference weight, a full-fibre
  conditional-resampling projector has either zero positive-lag
  memory on fibre-centered currents or nonstabilizing partial sums.
  `Thermodynamics/GraphDiffusion.lean` proves closed-graph summation by parts,
  typed Fick and Fourier laws, canonical amount and energy updates, exact
  source balance, and source-free conservation. These are finite conditional
  transport identities; source dynamics, physical distance, clock,
  constitutive calibration, stability, and a hydrodynamic limit are separate
  outstanding attachments.
  See `B8_TRANSPORT_KERNEL.md`.

The open foundational endpoint is:

> **Proposition 4.2 (Fixed-point reading of reality).** Define the public
> world as the quotient-normal form `World = NF(x) / ∼_gauge`, where `NF(x)`
> is the terminal state reached by accepted repair and `∼_gauge` identifies
> hidden local presentations with the same declared observable overlap data
> (Definition 4.1). Then `World ∈ Fix(Repair)` and `Repair(World) = World`.
> When OPH confluence and completeness conditions hold, this terminal public
> state is independent of update schedule on the physical quotient.

`Primitives.lean` supplies `OPHCarrier.Patch`, `Records`, `Obs`, `Site`, the
local accepted-step relation, a choice-canonical terminating `Repair`, the
exact mismatch potential `Φ`, gauge equivalence, its equivalence proof, and
repair congruence. It also proves strict descent, reachability, normal-form
production, and termination. `Confluence` and `Completeness` are typed
properties rather than universal conclusions.

A theorem-grade statement matching Proposition 4.2 requires the
paper-level public `World` quotient of terminal normal forms by gauge, its
fixed-point theorem, and schedule independence on the physical quotient under
the explicit confluence and completeness premises. The abstract-rewriting
result must be transferred to the structured OPH relation without replacing
those premises by a deterministic schedule.

The full quotient-normal-form theorem is work in progress as a single
`World` statement. The abstract-rewriting module is the generic skeleton, while
`Primitives` discharges concrete carrier and reconstruction subclaims and keeps
the remaining asynchronous repair obligations explicit so they cannot be
silently elided. See `PROOF_INDEX.md` for the proof-to-paper map and completion
tracker.

## Building

    cd Lean
    lake exe cache get        # fetch pre-built Mathlib oleans
    lake build                # build the proof libraries and Screen modules
                              # (ObservableNormalForms, ObserverPatchHolography,
                              #  EventAlgebra, OPHThermodynamics, OPHScreen,
                              #  OPHConstruction)

The `Main` console entry point is optional and not part of the proof receipt;
build it separately with `lake build oph:exe` if needed.

Lean CI runs on pull requests and pushes that touch `Lean/**`, and can also be
started manually. It builds from a cleaned project target, rejects admissions,
global axioms, and `native_decide` in the Einstein-branch sources. It checks
the audited Einstein theorem subset and locks the repository-wide inventory
of 23 compiler-trusted `native_decide` proofs.

## Provenance

- PR #299 (closed 2026-05-18 unmerged) shipped the abstract-rewriting
  skeleton as a claimed Proposition 4.2 formalisation. Audit verdict: the
  proofs are sorry-free but generic; they do not reach OPH-specific
  structure.
- This scaffold ports those proofs into a properly-built Lake project,
  applies accurate labels, and lays out the gap to be closed.
- Jonathan Hill contributed the substantive Lean formalisation and proof-audit
  work that closed concrete carrier/reconstruction subclaims, added non-vacuity
  witnesses, separated confluence from observer-facing uniqueness, and exposed
  the remaining asynchronous repair obligations explicitly.
- Coordination: "OPH LEAN Proofs" working group (Bernhard Mueller, Ben
  Cassie, Dula, Jonathan Hill). Cross-audit between auditors is required before
  PRs are merged.

## License

This formalisation surface is part of the OPH public repository and is
licensed under [Apache-2.0](../LICENSE). The main [LICENSE](../../LICENSE)
gives the repository-wide licensing map.
