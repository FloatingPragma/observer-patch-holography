# OPH model construction from local repair

The [executable-model hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md) is the
common physical target. The [completion DAG](../docs/COMPLETION_PROGRAM.md)
orders its construction gates, acceptance evidence and scoped stopping rules.
One declared repair/memory process must supply the clocks, particle/energy
interpretation and gravitational premises together; separate successful
benchmarks do not discharge that joint requirement.

The [matter and common-field construction](observer_particles/DERIVATION.md)
connects the native response matrices to a conditional fifteen-state internal
matter representation, one-Higgs interaction tensors and compatible scalar,
CAR and vector field operators. Its complete controls retain unselected
family coefficients, extra modes and geometry-response requirements. It
instantiates restricted finite hypotheses of the theorem library; it does
not derive the measured particle masses. The tested mass prescriptions fail
the charged-lepton hierarchy: the declared linear Weyl–Higgs candidate gives
coefficient ratios `1:1.1903:1.4982`. Under a direct pole-ratio identification,
the broad reference windows require middle/light in `[206,207]` and
heavy/light in `[3470,3485]`; the complex affine class then bounds heavy/light
by `543.55`. These are conditional diagnostics of unmatched coefficients;
a physical comparison requires joint kinetic/mass extraction and matching
([mass-response report](../docs/research/OPH_PARTICLE_MASS_RESPONSE.md)).
The [field-response guide](../docs/research/OPH_FIELD_RECONSTRUCTION.md)
distinguishes record reconstruction, supplied-action spectroscopy and the
physical identification of particles, including residues, currents and
regulator multiplicities.

This is the construction guide for the [model implementation](README.md).
The [repair-law guide](CANDIDATES.md) explains which parts of the repair law
are forced and which are free; the generated [premise map](PREMISES.md)
connects the constructions to the theorem library.

The premise map inventories every canonical claim row and tracked Lean source
module. Typed connections distinguish carrier, law, preparation, records,
controls, clock and regime; changing one requires an explicit mapping.
This checks coverage and source custody, not mathematical entailment. A
theorem without a checked model instantiation remains an open connection.
The [physical audit](PHYSICAL_AUDIT.md) gives the current completion status,
including the actual rank-three occupation-field action, spherical-screen
Green bounds, a common positive continuum generator and compatible
continuous field readouts, endpoint causal
placement, law-uniform finite preparation, complete-history uncertainty and
resource-clock calibration.

**Build the simplest shared repair process that satisfies the required
constraints.**

## What is being built

OPH proposes that the physical world is what overlapping observers agree on
once they have repaired their disagreements. The construction makes that
proposal executable. It fixes a finite carrier, a local repair law, observers
that read records and act on them, and the services those observers use. It
then runs the specification and derives what its observers can measure. The
aim is a specification whose observers' joint statistics admit standard
physics as an effective description.

Mathematically, an instance is a tuple

```text
(carrier H, accessible algebra A, repair instruments {Φ_{e,r}}, scheduler σ,
 observer instruments and feedback maps, preparation domain, declared services).
```

Each seam `e` carries an instrument whose outcomes `r` are committed records.
An observer is a finite tape of such records together with a controller that
chooses the next operation from the tape alone. Everything an observer can
report is a probability computed from these maps.

The working method has four steps. Build the simplest instance that meets the
stated requirements. Derive its consequences exactly. Keep every supplied
input visible. Find the observations that distinguish it from the
alternatives the requirements leave open. A finite implementation is complete
when its declared observation, record and feedback contract runs and passes
its checks. Proving that self-consistency uniquely selects that
implementation is a separate question and is not a prerequisite for
constructing it. A complete physical model additionally needs compatible
physical identifications and tests on the same process.

## Completion target and acceptance criteria

The concrete hypothesis prioritizes localized native defects for stable
massive matter and a derived maintenance-to-rest-energy map. Unstable and
massless excitations need their own same-process dynamics. The full theory
with a massless photon has no strictly positive global vacuum mass gap;
massive-sector thresholds and pure Yang--Mills gap statements have their
own domains. Modular parameters require reference states and physical clock
calibration. Coarse-scale gravity requires the same construction to realize
G1--G6, including geometric generalized-entropy stationarity and conserved
matter stress. These are acceptance obligations, not premises
discharged by running the current C/E code. See the
[common hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md).

**Build an executable realization of the OPH observer algebra whose
compatible inter-observer measurements reproduce the observed Universe.**
Their joint observations must admit a nontrivial 3+1-dimensional QFT as a
controlled effective description, including the observed matter and
interactions, and classical physics and dynamical gravity under appropriate
coarse-graining. Four-dimensional QFT and 3+1-dimensional spacetime refer
to the same dimension count of that description.

**No fundamental QFT or spacetime recovery is required.** The target concerns
what observers can measure and jointly interpret. Fields, particles and
spacetime geometry can be effective descriptions of those observations,
without being primitive objects or an underlying observer-independent
substrate. Implementation data structures do not change that ontology.
Success requires agreement of the accessible joint and conditional
statistics and intervention responses in each description's validity
regime, with controlled errors and corrections. It does not require every
effective description to remain valid at arbitrarily fine scales.

Meaning assignment connects actual preparations, instruments, interventions
and records to shared physical quantities. Consensus and repair must preserve
compatible joint and conditional statistics on overlapping observer domains.
This is the proposed operational sense in which an objective reality emerges.
It does not follow merely from observers agreeing on an interpretation.
See the [interpretation contract](../docs/INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity).

| Completion obligation | Acceptance evidence |
| --- | --- |
| Executable observer algebra | Explicit state, operation, preparation, measurement, retained-record and feedback maps; their algebraic identities and observer-access limits; accounted memory and service assumptions; actual execution of the claimed protocols. A finite two-event implementation discharges its own contract, not every observer or federation interface. |
| Shared operational meaning | Interpretation maps specified independently of the validation outcomes, compatible across observers and refinement, and preserving the declared intervention statistics and quantum correlations. Shared labels alone are insufficient; incompatible unperformed measurements need no common table of classical outcomes. |
| Effective quantum-field description | A common refining family of observer statistics and interpretation maps whose errors relative to a nontrivial 3+1-dimensional QFT tend to zero in the declared continuum domain, including the observed matter and interactions. The comparison concerns measurable correlations and intervention responses, not recovery of fundamental field objects. Where an effective vacuum-sector mass gap is defined, it is nonnegative and may be zero. |
| Effective classical and gravitational descriptions | The same observer construction yields observations that admit compatible classical, geometric, stress and clock descriptions in their stated regimes, reproducing the gravitational predictions with controlled corrections. A fundamental metric or manifold is not a prerequisite. Separate calculations on a supplied geometry need a map to these observer statistics before discharging this obligation. |
| Physical predictions | Explicit parameter/input ancestry, calibration and uncertainty; predictions and rejection rules fixed before independent comparison data; agreement with the observed statistics, including alternatives and distinguishing tests. Reusing a calibration datum is not an independent prediction. |

Exact equality of limiting observable laws is the mathematical target.
Finite simulations and experiments are assessed with declared approximation
and measurement errors. They do not prove arbitrary empirical precision or
force a particular random outcome in an independently repeated experiment.
Agreement must be constrained by measurements; a consistent world model can
describe different physics from nature.

### The observable correspondence to establish

Fix an observer region, preparation family, accessible instruments, feedback
rules and finite resource domain before comparing descriptions. A physical
dictionary maps these operations and their complete outputs to the proposed
effective field description. It must include local clock comparisons, signal
records and field readouts together. Different observers' dictionaries must
agree on their shared records and allowed joint interventions. Choosing
seconds and metres sets units; it does not determine these response laws.

For a protocol `p` in that fixed menu, let `rho_native(p)` and `rho_eff(p)`
be its complete record/quantum outputs, including failures, unresolved waits
and any admitted untouched reference. Map both outputs to one declared
comparison algebra with decoders `D_native` and `D_eff`. The finite target is

```text
sup_p (1/2) ||D_native(rho_native(p)) - D_eff(rho_eff(p))||_1 <= epsilon.
```

The supremum also covers the related initial states and references declared
by the contract. For classical records this is total variation of the full
joint law. Matching individual expectations, two-point functions or clock
means is weaker. Conditioning on a rare successful outcome needs its own
probability-dependent error bound; it cannot erase failed branches from the
comparison. Decoders must retain the observations being claimed: mapping
both outputs to a constant is not an observation-faithful correspondence.

The finite source results provide parts of this contract. A physical QFT
attachment must additionally supply an operational domain with three spatial
and one temporal dimension, the appropriate covariance/locality and state
conditions, and measured interaction responses. Its approximation errors
must be controlled across the claimed domains or refinement. A finite
positive transfer and its spectral logarithm alone do not provide those
hypotheses. The [premise map](PREMISES.md) records which theorem interfaces
are instantiated and which require extra inputs.

The empirical step uses the *same* dictionary and calibrated parameters on
independent observation protocols. Existing public data can supply that test;
a physical emulator or new laboratory is unnecessary. An assumed
experimental outcome is a conditional forecast and cannot validate the
dictionary. The [public-data audit](../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
keeps the ACT comparison separate from native-observer derivation.

The [finite common local-clock observer](native_modular_reference/README.md#a-common-observer-from-two-local-stopping-clocks)
prepares one faithful reference using two incident counters and a reusable
native coin auxiliary. Its bounded selectors and waiting deadlines retain
all failures. One source-compiled modular action and two continuing native
current instruments operate on that prepared mixed carrier, with complete
quantum/reference error and resource accounting. This derives preparation,
scheduler randomness and the readout instrument on the declared finite
domain. Spectator-uniform completed cap operations remove exterior-state
selection; coherent proper-code services remain architecture inputs. The ideal
reference's modular action preserves all
core-containing cap algebras; finite preparation and source words carry the
stated errors. On the binary total-two shell, requiring every bare two-port
algebra to be preserved instead forces a product reference even among
initially noncommuting states, and
connected symmetric equilibrium selects its tracial member. This exact
boundary prevents treating unrestricted regional compatibility as automatic.

The native spherical construction supplies one such mathematical connection:
an ancestral-cell Poincare bound makes the Green subsequence converge in
operator norm to an injective compact `G`. Its inverse `A` determines the
limiting spatial covariance, native heat propagator and admitted coherent
wave amplitudes. Its conservative regular strongly local form has domain
`H1` and spectral dimension two on the spherical observation domain.
An [explicit positive conductance graph](physical_space_time/OCCUPATION_FIELD.md#a-local-conductance-graph-with-the-same-continuum-green-limit)
has the same Green limits, with finite operator error `182*2^-level`.
[Actual density paths](observer_time_fields/REFINEMENT.md#continuous-paths-of-the-actual-native-density-and-green-field)
converge after interpolation in `C([0,T],H^-3(S2))`, with covariance
`exp(-|tau-u|A)` at `tau=attempts/(1260*16^level)`; the declared Green
filter gives `H^-2` paths. This simultaneous refinement/count limit
retains uniform fixed-number preparation, IID seams and full readback.
Normalized [incident counts and inverse clocks](relative_observer_geometry/README.md#incident-clocks-on-the-simultaneous-density-continuum-scale)
sample these same paths, including [predictable port choices](relative_observer_geometry/README.md#predictable-port-access-on-the-density-continuum-scale) with degree-weighted
ticks. This path result needs no clock–field independence; finite-vector
shared-seam jitter retains its separate fixed-event-set contract.
Whole-sequence uniqueness and physical clock calibration are unproved.
With the constant mode restored, a heat off-diagonal estimate gives finite
wave propagation in the supplied flat comparison distance. An
[active incident clock](relative_observer_geometry/README.md#active-coherent-waves-in-a-native-incident-clock)
drives one coherent step per tick on a protected register and gives
`cos((2d/sqrt(630))*tau*sqrt(A))` on bounded spectral bands. That step is a
compiled per-mode unitary on supplied codes, with angles compiled from the
spectrum of `A` and one step per native incident. Passive repair supplies no
wave sector: [ballistic rescaling freezes the native density paths](observer_time_fields/REFINEMENT.md#continuous-paths-of-the-actual-native-density-and-green-field)
and no shared clock conversion turns the density log into a wave dispersion. Its exact
two-attempt Ramsey probability is `1297/1320`; carrier pauses, coherent
services and a vanishing accumulated implementation error are explicit.
The [source-coded quantum fluctuation protocol](observer_time_fields/REFINEMENT.md#faithful-state-clock-compatibility)
retains ordered Weyl phases and kick responses. Faithful mode preparation
matches the active clock precisely when `log((1+r)/(1-r))=b*Omega` with
common `b` and matching axis. The native spectral bound makes the
thermal Fock state of `N` supplied qubit codes normal; its modular
conversion `s=-(2d/sqrt(630))*tau/b` on the encoded mean-zero algebra is
the finite KMS identity of that chosen Gibbs state. It does not supply a spatial-patch boost,
selected temperature or laboratory scale. Mode encodings, protected controls
and the separate zero-mode choice retain their contracts.
Controlled repair gives derivative `2(I-exp(-tau*A))f`;
at fixed scaled selection-slot time its [disjoint-support response](observer_time_fields/REFINEMENT.md#heat-response-and-the-boundary-of-a-finite-signal-cone)
is nonzero with an exponentially small short-time bound. An exact finite
signal cone is incompatible with that same support/control/time dictionary;
its physical locality is unproved. The [boundary-action test](physical_space_time/OCCUPATION_FIELD.md#a-boundary-trace-test-for-a-radial-extension)
excludes a boundedly invertible `L2` energy identification with a smooth compact three-dimensional
Dirichlet bulk's static boundary trace (with smooth boundary). Internal rank-three components
vanish in bounded-energy families under this scaling. The
[normalized-lag invariant](physical_space_time/OCCUPATION_FIELD.md#a-complete-linear-readout-preserves-the-decay-spectrum)
also excludes complete linear identification with a compact three-dimensional
free scalar's Euclidean law under one clock scale. Both tests detect that the
native lag law is diffusive and second order: native rates grow like the
eigenvalues of `A`, Euclidean rates of a relativistic free scalar like their
square roots, so the same argument excludes such a scalar in every dimension
two or higher. Neither is a test of three-dimensionality. These restricted
tests leave independently encoded radial fields outside their scope.
Positive subordination of native lag records gives `exp(-r*sqrt(A))`.
A [declared spectral filter of `A` applied to whitened consecutive long-record blocks](physical_space_time/OCCUPATION_FIELD.md#independent-native-records-give-a-radial-product-field)
synthesizes a Dirichlet product Gaussian field with covariance
`(-d_r^2+Abar)^-1`. Consecutive blocks of one history feed consecutive
radial slices; any independent Gaussian blocks of the same covariance give
the same law, so the radial direction carries no native information beyond
the spectrum of `A`, and the same filter attaches an interval to any
two-dimensional generator.
The [thermal join](physical_space_time/OCCUPATION_FIELD.md#thermal-product-covariance-and-the-shared-wave-operator)
applies a declared thermal filter, distinguishing the quantum covariance
from the classical inverse, and controls finite native counts,
radial grids and the full-port/Schur correspondence. The [two-axis clock](observer_time_fields/REFINEMENT.md#one-product-field-state-and-active-clock)
combines native angular phases and declared radial pulses, giving the same
product generator, thermal CCR state, response and modular conversion with
bounded-band errors. The chosen hyperbolic wave action on this product operator has a
[measurable Lorentz-signature metric density](relative_observer_geometry/README.md#lorentzian-density-of-the-product-wave):
the signature is the choice of action, and the measurable spatial tensor
is the derived part.
Complete native total-one calibration records also determine conditional
coordinate increments. Spatial refinement at fixed density time, followed by
short-time polarization and chart localization, recovers the angular tensor.
The one-attempt-first limit and a finite cell-reset kernel fail explicit controls.
For the declared product wave, intrinsic energy distance equals the earliest
possible nonzero operator response between positive-measure regions, as for
any strongly local regular Dirichlet form. This
uses integrated short-time heat asymptotics and heat-wave transmutation;
a chosen detector can remain dark. It requires neither pointwise smooth
coefficients nor an identification with native repair speed.
The same metric is smooth if the full action density is smooth. Uniformly
elliptic smooth approximants converge in form and propagator senses without
establishing smoothness or curvature convergence of the source limit.
Distinct radial graphs give identical complete radial-average angular laws
and different radial responses: those observations cannot select an interval
or its transport. The interval, preparation, spectral access and controls
remain supplied. The restored angular constant has positive product frequency.
An [executed edge filter](predictions/README.md#preparing-the-finite-edge-reference-from-source-pulses) prepares approximants to `exp(-tL)/Z`
from supplied `I12/12`, retaining failures and bounding free-entropy error.
Its [capped subordination](predictions/README.md#positive-subordination-on-the-prepared-register)
approximates square-root Gibbs weights with complete failures and explicit
state/modular-gap bounds. Exact uncapped sampling has infinite mean cost.
The [same-register continuation](predictions/README.md#prepared-reference-clock) applies a
[source-compiled clock and pair detector](predictions/prepared_clock_source.py)
to that actual finite preparation. For its ideal successful polynomial `M`,
`H=-log(M/m0)` fixes the clock; native incidence and first-current feedback
drive the complete two-attempt experiment. Sampler and herald failures,
padding outcomes and quantum references remain retained. The modular identity
`sigma_s=alpha_(-s)` holds by this definition of `H` from the prepared
density; it concerns the ideal successful twelve-dimensional work marginal
after forgetting the request and tracing the reference. The finite clock is
not `sqrt(L)`; raw source checks are numerical guards, not interval or
physical apparatus certificates.
All eleven constant-to-mode fluctuation pairs inherit their covariance,
commutator and clock frequencies from this one density. Three representative
bands have complete source pointer readouts of an ordered Weyl product;
independent successful copies give a finite-copy error bound. The joint
copy/band estimate scales as `Z/N` for fixed positive gap and probe bounds;
`N` proportional to dimension can fail. Postselection error is amplified by
the copy count, so fixed source precision cannot justify this limit. Scaling
coordinate and momentum by inverse square root and square root frequency
respectively gives the thermal covariance and unit-normalized sine response.
The positive-log bound then connects the same clock to square-root dynamics.
For `D_l=4^l L_l`, the finite native spectrum gives positive eigenvalues
`lambda_n>=n/81940`, a uniform prepared partition bound
`Z<=(1+327760/b²)/m0`, and an explicit positive gap. Native realizability
requires `p=630*h*16^l<=1`. A joint schedule `k=2l+c`, `R=2^k`,
`h=b²*4^-k` meets it, with complete finite-word encoder, filter and pointer
costs. The retained zero-request branch has mass at least `2^(-2R)`;
precision scaled to this floor controls the actual marginal's logarithm,
including imperfect initial purification and every failure. Quantified
arbitrary-precision load/export services remain premises. This is a
constructive resource theorem, not an execution of the growing product
Fock preparation; the positive radial product additionally needs an anchor
sector and controlled spectral tails. Geometric entropy stationarity, physical gauge collars and
full-group sewing remain separate obligations.
The full native-return reference approaches normalized trace
at fixed modular parameter, while its density-wave sector has vanishing
reference weight. The
[common-generator proof](observer_time_fields/REFINEMENT.md#one-continuum-operator-for-fields-native-decay-and-coherent-waves)
and [controlled response](predictions/README.md#a-controlled-continuum-mean-with-an-executed-event-clock)
state their observation access, parameter conventions and limit orders.

The [quantum incident-clock construction](observer_correspondence/QUANTUM_CLOCK.md)
complements the [exact full-state representation](quantum_observer_correspondence/README.md)
with an ordinary CP local-output comparison on the native carrier.
A positive mixture of actual remote and incident repairs
produces the effective instrument; the proof includes adaptive local
interventions and untouched references. Coherent two-tick executions retain
all 81 outputs, error budgets, service costs and reference correlations.
The correspondence uses explicit preparation encoding, trivial on the
retained inputs, and does not identify global source states. This completes
the stated finite quantum observation contract while leaving the physical
QFT and calibration obligations above in force.

Hardware, protocol and implementation choices are open design variables.
No uniqueness theorem for the chosen rule is a prerequisite. For the current
P/N construction, closure-derived P and N are accepted inputs under their
declared contracts; that assumption does not discharge the other physical
interfaces. Supplied background geometry, action, state, interpretation or
clock data must stay visible until the same construction accounts for them.

These obligations can advance in parallel. Passing CI verifies the retained
software contracts; it does not mark the physical target complete. The
current native observer, P detector and source-compiled two-site scalar
apparatus are completed finite components under their stated premises.
The controlled effective interacting-field and gravitational descriptions
of the observer statistics, and independent empirical completion, are open.

## Why this specification

Every choice below is either forced by a stated requirement, chosen as the
simplest member of an explicitly classified family, or declared as an input
with its alternatives kept. The [repair-law guide](CANDIDATES.md) gives the
law-level arguments in full.

| Choice | What the model uses | Reason | What stays open |
| --- | --- | --- | --- |
| Carrier | Twelve ports and thirty seams of the source icosahedron; integer loads with total two | Declared by the architecture premise; total two is the smallest sector where a single load waits or moves and two loads collide | Selection of the architecture from an unrestricted menu; extended carriers and refinement |
| Classical seam rule | Endpoints become floor and ceiling of their mean; a fair coin breaks ties | Forced by conservation, one-step balance, a tie-break that reads no load, and endpoint symmetry | Nothing on its domain; the four requirements are premises |
| Record law | C executes; E is the comparison law | C records the classical rule's own random variables, the coin and the current, and nothing more | Which coherence `g` nature keeps; the admissible laws form one convex family |
| Scheduler | Uniform independent seam attempts as baseline; fixed nonuniform IID weights can be recovered from exact pair-clock laws | Clock ratios and coincidences determine weights; a periodic control shows they do not establish independence | Temporal independence and identification with physical clock observations |
| Observer | A finite record tape and a controller that reads committed records and chooses the next seam | The smallest object that records, predicts and acts on its own records | Larger observers, record delivery and refinement |
| Services | Coherent preparation, readout, mixers, record storage and compiled source words | Declared and costed in every error budget | Physical realization and authenticated service errors |
| Closure inputs | P and N as derived by the closure equations | Accepted under their declared contracts | Physical endpoint corrections |
| Field benchmarks | Supplied scalar, gauge and matter actions on a common source | Conditional comparisons with explicit errors | Derivation from native repair |

**The record law in one paragraph.** Requiring seam locality, the exact
classical law, local number-phase covariance and endpoint symmetry leaves a
convex family of quantum record laws. Each member is a table of five
coherence multipliers, one positive matrix per transported current. Four
members stand out: C0 records only the current, C adds the coin, E adds the
seam parity and T adds the seam total. Exactly E and T are idempotent, and
every law that keeps an empty seam coherent with a full seam completes to E
under repetition. Under the same current-record interface all members share
every classical observation. One
coherent interference experiment reads the multiplier `g` directly, with
predictions `3/8 ± Re(g)/2` at phases `0` and `π`. The model executes C as
its convention among the classified laws: its records are the classical
rule's own random variables, the coin and the current, and nothing more. The
[record-law classification](record_laws/README.md) proves these statements
on the 78-state carrier and checks them against the executed C and E
instruments.

**The observer in one paragraph.** An observer here is the smallest system
that keeps records and uses them. The [native observer](observer/README.md)
commits the coin and current of one repair, forecasts the next current from
that record and chooses its next seam from the current alone. Its controller
sees no hidden load. Under the same current-record interface every admissible
record law records the current, so the same controller runs under each of
them. Larger observers compose such
tapes, federations join elementary carriers, and detectors attach to the
record through declared couplings.

## Current executable result

The [primitive observer runtime](primitive_observer/README.md) executes
complete instruments on mixed states, committed records, feedback and
untouched references. Source-derived native dynamics and supplied physical
benchmark programs use this same execution interface. Their preparations,
clocks, interpretation maps and error budgets remain explicit.

### Local construction and its extension

The [coherent screen observer](coherent_screen_observer/README.md) executes
one quantum extension on the matched level-three and level-five screens,
evaluating its channel on their common 222-port ball. The input is a logical qubit on two root ports, with any untouched
reference, and a separate Bernoulli environment. Its two native incident
ticks, seam-category records, local phase feedback and complete reduced
quantum output compose without a state reset. A bounded operator dual
evaluates the common 222-port channel; the existing complete-output locality
theorem controls its comparison with both full screens. This discharges a
finite coherent screen comparison under that preparation contract, while
retaining the all-level and physical interpretation obligations.

The local proof obligations start on one twelve-port carrier: specify its
accessible algebra, derive repair and record laws, compare its local clocks,
and reconstruct spatial relations from those observations. Observer regions
may overlap inside the same carrier. A federation is useful for tests involving
more carriers, larger distances or refinement; none is required merely to
define an observer or compare two local record views.

An extended spacetime claim requires a compatible family of operational
domains and an effective regime with controlled errors. The existing
spherical federation is a candidate for that extension. Its geometry and
refinement must earn their physical interpretation through observations;
they are not prerequisites imposed on every finite construction. A single
carrier can establish a finite contract without establishing that extension.

Action representations also have different domains. The endpoint Kraus
instrument determines an exact reversible representation of its finite
histories. Choosing a unitary completion and a Hermitian logarithm is an
encoding choice. It supplies no prediction for an extra intermediate-time
measurement until such access and evolution are included in the observer
contract. A valid construction need not be uniquely selected; its advertised
observations must follow from its declared operations.

The binary current calibration and quantum clocks share the same 66-state
sector. Native propagation, occupation covariance and the slow spatial Gram
can therefore be compared without using a doubled-load preparation. A
separate current-fluctuation action is derived in one port's tick count, with
an exact Palm-conditioned triangle variance `3m/5`. Actual host-counter
records test a software clock correspondence with held-out blocks and a
workload control. The [joint clock bridge](native_clock/CLOCK_BRIDGE.md)
derives incident-hit/current-count laws and finite conversion errors from any
initial configuration, including zero hits and a separate feedback-count law.
Its reanalysis preserves the host calibration and workload rejection. These
results retain their total-two preparation domain; neither a traceable lab
clock nor the full-screen field protocol is calibrated by that comparison.

### Native action, federations and observer-relative clocks

The [quantum federation construction](federated_quantum_sector/README.md)
joins identical elementary carriers using the declared source wiring. Native
Kraus maps preserve signed-difference operator blocks, including individual
outcomes and reference correlations. Its complete cross-carrier experiments
retain coherent interventions, all records and conditional quantum states.
The [local action](native_action/README.md) is fixed by these instruments:
`H_native=sum_e(I-Phi_e)`. Its positive local terms reduce to exclusion plus
coherence killing. The original discrete transfer is `I-H_native/M`; the
optional per-edge Poisson extension is `exp(-s H_native)`. This is a local
Liouville action, not a selected relativistic Hamiltonian or physical clock.

The same construction controls truncation of distant ports for complete
local record, quantum and reference outputs. For an observed region A,
maximum degree Delta and retained radius R, the Poisson-parameter error is
at most `min(1, |A|(Delta s)^(R+1)/(R+1)!)`; the fixed-count bound is given
separately. A conditional infinite-federation process follows on a fixed
bounded-degree graph with compatible local states and operations. The
[declared screen refinements](federated_observer/REFINEMENT.md) require
compatible access and dynamics; graph rank or finite positivity alone does
not establish a relativistic spacetime limit.

Actual source neighborhoods are identified at levels 1/3, 2/4 and 3/5,
with preserved local port labels and explicit complete-output comparison
bounds. This supplies finite cross-size observation maps rather than assuming
that parent/child coarse averages intertwine dynamics. Matched initial
restrictions are necessary; a uniform global total-two preparation dilutes
locally as the screen grows. An all-level routing argument and a controlled
nontrivial physical scaling regime are open.

The [incident-clock attachment](native_action/INCIDENT_REFINEMENT.md) compares
these actual screens after the same number of local repair opportunities,
without supplying a duration per operation. Its stopped-word identity and
source boundary counts bound the complete quantum and record output, including
local feedback and untouched references. Four incident ticks give a pair
allowance below `6.656e-9`. A compatible density-`1/3` preparation supplies a
certified nonzero two-time field covariance, while full-screen executions
retain the actual seam scheduler, current records, feedback repair costs,
timeouts and serialized continuation. This closes the finite clock/refinement
connection; the all-level physical interpretation obligations remain separate.

The [finite-density field attachment](observer_time_fields/REFINEMENT.md)
uses that same preparation, root clock, current support and native feedback
to certify a complete four-bit field law, including every finite-cap timeout
prefix. The positive incident transfer preserves this law with the Bernoulli
stationary boundary or with each screen's own canonical fixed-number
restriction. The latter removes a product-preparation premise: native
connected exchange derives the unique diagonal equilibrium conditional on
the conserved total. Its exact degree-two exclusion dual retains conservation
correlations without enumerating the full binary space. At density `1/3`,
feedback reverses the connected four-field statistic on both screens, while
their complete field laws have a resolved finite-size difference. Two field
counters and controller memory survive serialized source continuation;
arithmetic, cutoff, spatial and preparation differences are accounted
separately. Selecting the conserved density and realizing equilibrium
preparation remain inputs, with no large-screen mixing-time claim.

The [temporal-record representation](relative_observer_geometry/README.md#temporal-hilbert-space-reconstructed-from-field-records)
reconstructs a seven-dimensional cyclic positive transfer from exact field
correlations in the elementary total-two, two-port-clock domain. Exact
flatness proves its recurrence under the established positive-transfer
class; ordinary samples do not certify that rank. On the distinct full-screen
product preparation and twelve-port root clock, bounded correlations certify
a two-vector temporal subspace and positive transfer compression. Its Ritz
values are not closed modes or a physical energy/clock identification.

The [local event clock](relational_event_clock/README.md) counts actual
occupation changes at one observer port. Its stopped instrument retains
remote repairs during the wait, finite timeouts and conditional quantum
outputs. The [joint clock process](joint_clock_process/README.md) compares
two such clocks and an observer-chosen field algebra on the same continuing
cross-carrier history. Exact elimination of the hidden factor gives a
record-conditioned quantum memory recurrence, valid with arbitrary untouched
references. Completed events in this specified experiment admit a four-mark
measure-and-prepare description; unresolved waits retain quantum information.
Instrument-valued clock conversion also preserves the count of B events
between A events. Direct native elimination and elimination through the
joint event instrument agree, including B-triggered field feedback. The
passive two-A-tick renewal cycle has B-count mean two and variance 16;
equating mean rates would discard observable clock fluctuations.
Neither reduction identifies native attempts or event counts with modular
flow. A separate [active event clock](transfer_modular_probe/README.md)
triggers one compiled modular step per incident event, using the generator
of the native three-state transfer. Its controller pauses the native carrier
during coherent probe service. Closed native histories derive a faithful
reference for that generator; passive repair itself has a strict modular
automorphism obstruction. The [recovery contract](OBSERVER_RELATIVE_RECOVERY.md)
states these access and scheduling conditions.

The [native field identities](common_source/NATIVE_FIELDS.md) derive exact
continuity, current noise and mean diffusion under the uniform scheduler.
[Quantum transfer](native_quantum_history/README.md),
[record-clock laws](native_clock/README.md),
[joint current/phase readouts](native_relational_readout/README.md),
[replicated fluctuations](native_fluctuations/README.md) and
[one-history circulation](native_circulation/README.md) provide additional
controlled finite or fluctuation descriptions of that source. The recovered
rank-three frame maps seam currents to an axial martingale with stationary
covariance `(4N/99) I3`. From every binary initial state, its rescaled path
converges to three-dimensional Brownian motion with increment-action
coefficient `99/8`. Its norm has the derived radial generator
`(2/99)(d_r^2+2/r d_r)`; the joint polar operator also has
`(2/99)r^-2 Delta_S2`. This is the radius of a retained current counter,
with time parameter rescaled native attempt count. It neither supplies
spatial translations nor selects the separately encoded radial field.
Uniform scheduling, current access and growing record capacity remain
explicit; stationary preparation is needed only for the exact finite covariance. Replication is not spatial refinement, and Gaussian limits
do not select a fundamental field action.

### Continuing measurement predictions

The [compiled overlapping-clock history](compiled_clock_history/README.md)
replaces the coherent correspondence's ideal preparation, intervention and
pointer with emitted source controls. The same continuing native C process
retains both clocks, feedback, timeouts and all code-complement aborts.
Its 512-address execution uses exact rational arithmetic on the returned
binary source arrays and certifies all quantum/reference outputs within
`1e-11` of the established ideal correspondence. Maximum apparatus costs
are 51 pulses, or 54 with the compiled Bell preparation, with two code
load/export services per pulse and at most five native repair attempts.
This supplies a finite compiler-to-observer connection, with physical code
services and the separate original scalar benchmark explicit.

[Native equilibration](native_equilibration/README.md) replaces an exact
stationary preparation by a finite sequence of discarded repairs within a
declared fixed-number sector, uniformly over untouched references. The
[sequential prediction protocol](closed_prediction_protocol/README.md)
applies that preparation before each trial on the same continuing carrier.
It bounds the dependent complete transcript without an IID or exact reset
premise. The locked design uses 379 discarded attempts before each of 3,927
trials and three unread attempts in each useful phase probe. Its exact-source
wrong-candidate probability is at most 0.75%; a separate, unauthenticated
physical service allowance would suffice for a 1% target. Failed heralds
remain trials. The producer executes two consecutive trials per candidate
and reference controls, not the entire proposed experimental campaign.

The [repair-generated scheduler](repair_scheduler/README.md) removes the
separate random-seam service and exact initial preparation from a complete
version of this experiment. Seventeen deterministic sweeps first absorb
arbitrary total-two collision inputs and prepare the binary uniform state
within `259/10^9`, uniformly over untouched references. A separate native
carrier's coarse coin instrument supplies uniform seam choices without
resetting that carrier or assuming its initial independence. Bounded
rejection retains explicit timeout outcomes. Two continuing trials per law
execute all 144 public records; the fixed 3,927-trial schedule has
wrong-candidate bound below `0.006061`, with both carriers' costs counted.
The native fair instrument, auxiliary access, coherent mixer, record handling
and physical identification remain inputs. The deterministic scan changes
preparation only; it does not replace the uniform-scheduler clock laws.

The [native configuration controller](native_clock/configuration.py) derives
conclusive basis measurement and basis preparation from signed repair
currents. Its deterministic adaptive seam choices use the committed tape;
zero-current prefixes never license a guess. Finite retries retain a timeout
outcome and give a uniform error bound, including arbitrary references.
Preparation discards its private tape; retaining that tape would preserve
conditional reference information and change the comparison instrument.
This active protocol does not replace the passive IID clock laws above.

The [labelled field preparation](observer_time_fields/README.md#native-preparation-of-the-labelled-field-state)
composes this reset with the native auxiliary selector. It supplies the
uniform total-two state together with its two initial public field bits,
within a complete reference-safe error. The same carriers then execute four
incident-clock field reads and committed-bit feedback, with all preparation,
selection and waiting failures included. This implements the initial
state/record premise of the elementary positive-transfer theorem without a
stationary-state oracle or independent blank. The screen and continuum
preparations retain their separate contracts.

The [passive-age readout](native_clock/README.md#passive-age-from-two-current-derived-endpoints)
uses two current-derived configurations around an unread interval. An exact
native eigenfield gives score mean `(9/10)^n`, uniformly over the starting
state. Complete continuing trials invert that mean to a finite set of attempt
counts, with martingale confidence and both carriers' resource costs.
This service needs no coherent clock or stationary preparation. It measures
a fixed unknown interval length across trials; it does not supply its stopping
gate, an autonomous oscillator or a physical duration. Both services use
declared joint record access and finite controller memory.

The [native auxiliary pointer](native_interface/pointer.py) uses a prepared
orientation bit and the existing source controls to measure a work bit.
Tracing the auxiliary gives the full work-bit Lüders instrument, rather
than just its probabilities, on native and padding addresses. In contrast,
coarsening a complete configuration measurement dephases more than a bit
measurement. The [full Palm clock](native_modular_reference/README.md)
composes native preparation, its source word and this pointer with complete
timeouts and a common budget. Coherent code routing and controls remain
architecture inputs: committed repair from a joint diagonal accessible input
with fresh independent diagonal ancillas cannot create their missing
coherence. Source-array errors and physical service errors
remain separate. The [finite-law discriminator](predictions/README.md) and its
[encoded realization](observer/source_discrimination/README.md) supply
additional complete protocols under their own preparation and service
contracts. Under the same current-record interface, mean loads and native
classical histories assign the same probabilities under every admissible
record law; the discriminator reads the coherence `g`, which separates C
from E. No public laboratory dataset is
identified with these native-law experiments.

The [complete native state readout](native_interface/README.md#complete-state-identification-from-native-records)
uses those current records and compiled quadratures to identify all 4,356
Hermitian coordinates of a 66-state preparation. Its 4,291 settings retain
all timeouts and export aborts. Reconstruction and modular-action error
bounds make the finite reference operationally identifiable. Bounded native
preparation before every trial supplies a fixed target on a continuing carrier
without independent-copy or successful-trial conditioning assumptions. Coherent
services and physical instrument accuracy remain inputs.

The [common entropy envelope](native_equilibration/common_variation.py)
follows from the mixed-complement preparation and one additional core repair.
Two conserved occupation moments determine a unique coherent entropy maximum
and an explicit smooth family at faithful interior moments. All 1,024 core-containing cap decompositions
share that envelope; 18,601 native work attempts give a reference-uniform
finite preparation bound. This is a specified preparation protocol. Counted
edge multiplicities have no derived gravitational-area identification.

The [complete native history](native_equilibration/README.md#coherent-history-stationarity-and-recovery)
also derives stationarity for coherent state variations through adaptive
records and disturbed endpoints. Its exact Petz response distinguishes
recoverable coherences from directions erased by the readout. Subtracting
the actual expected coin cost and state-independent flag normalization gives
a nonpositive quantum entropy Hessian at fixed trace and both occupation
moments. For the declared census, the equality algebra is
`M28 ⊕ 4M8 ⊕ 6C`; complete C timeouts leave 558 matrix units output-flat.
All 1,024 core-containing cap completions are checked. These are initial-state
variations on a fixed native carrier; physical geometry and stress remain
unidentified.

The [native-screen slab vacuum](native_equilibration/vacuum_geometry.py)
supplies a same-action joint state/geometry control. At fixed proper slab
volume, its bare regional entropy and two-face area have derivatives
-1.1004136186 and -4/3. Their signs agree throughout the nontrivial uniform
shape family, so no positive bare-area term makes this finite slab stationary.
Port-local canonical stress changes also fail to determine its regional
entropy slope. This does not test a renormalized geodesic small ball or
construct physical stress, curvature or one isotropic refinement cutoff.

The [proper-ball modular calculation](native_equilibration/ball_modular.py)
instead uses the full three-dimensional Gaussian vacuum and an equal-spacing
Cartesian cutoff in a declared proper frame. It computes regional modular
charges, entropy first-law secants and nested-region relative-entropy
monotonicity. The fixed conformal stress comparison has finite discrepancies
of 2.124% and 5.716% for two momentum profiles at the retained fine cut;
all grid/box controls and a finer-cut precision refusal remain visible.
An independent space-form calculation derives the fixed-volume area
coefficient `-4*pi*R^4/15` per `G_tau,tau`. A common native continuum limit,
renormalized area density and joint state/geometry stationarity are not
supplied by those finite comparisons.

The [same-action fermion vacuum](common_source/fermion_ball.py) extends the
regional comparison to the calibrated M15 action: 21 Dirac fields (six quark flavors in three colors and three charged leptons)
and three massless Weyl fields share one antiperiodic spatial cutoff and
ball. An exact region-local particle-hole map preserves both CAR relations
and charge. Regional modular response and nested relative entropy follow
from the restricted global vacuum projector. A 32-CAR source preparation
executes 13,616 native pulses; the larger grids are covariance calculations.
At the retained fine cut, massless and running-top comparisons with the
conformal stress weight differ by -0.1953% and +22.83%, respectively.
The massive comparison is not a conformal modular identity. Gauge-region
constraints, physical vacuum selection and joint stationarity remain open.

The [native geometry/vacuum repair](native_equilibration/vacuum_repair.py)
combines a uniform native geometry bit with conditional low-energy field
preparation on a fixed-volume slab family. Exact Rayleigh enclosures certify
the executed band preparations; a separate small joint reset retains every
input and native record. Fresh blanks and entropy export are required. The
stationary switching cycle consumes strictly positive work, including a
separate bound on the literal source preparations. A common repaired state
therefore does not supply autonomous cooling, reciprocal geometric response
or the required renormalized entropy stationarity.

The [reciprocal thermal construction](native_equilibration/reciprocal_vacuum.py)
replaces uniform geometric resampling with fixed-gap thermal collisions and
derives a joint Gibbs fixed point. It uses the same small scalar geometry
pair and a fixed-charge sector of the correlation/exchange field. Matter
free energy determines nonuniform geometry probabilities. Native proposals,
energy readout and selected complete collisions are source-executed, with
rejections and bath records retained. Fresh thermal inputs and the full
reservoir-bank controller remain supplied services. This finite thermal
stationarity is a probability variation, not a metric or entropy--area
variation, and does not certify autonomous controller work.

The [zero-temperature wedge construction](native_equilibration/wedge_modular.py)
uses the free native-coefficient staggered action and its filled vacuum.
An explicit transverse reduction preserves the half-space and both Dirac
tastes. The known exact lattice modular theorem maps to an independently
defined lapse/boost weight, with a directed continuum-normalization bound.
A Legendre-polynomial estimate supplies local covariance control as the
preparation volume grows. An eight-CAR source preparation separately has
an exact binary energy-gap fidelity certificate and full-Fock entropy
first-law check. Transverse factorization gives a bare planar entropy
density on a specified geometric area. Larger grids are ideal vacuum
calculations; physical service accuracy and vacuum preparation remain
requirements. The free half-space attachment does not establish a
renormalized area coefficient or joint small-ball metric stationarity.

The [bare-cap clock probes](native_modular_reference/cap_local_clock.py)
operate on the same locally prepared mixed reference. Uniform rotations over
all spectator configurations remove exterior-state selection from completed
operations. Both isolated current means read the same two modular gaps; their
ideal complete-trial expectation has a bounded inverse on `s` in `[1/2,1]`.
The retained compiled continuation is certified at `s=1`.
Continuing probes retain disturbance and feedback. Intermediate source-pulse
locality, physical service fidelity and conversion to seconds remain separate.

The [native reference readout](native_modular_reference/README.md#learn-the-common-reference-from-one-native-current-setting)
identifies this repair-derived entropy family without supplied reference
weights or coherent tomography. A finite census derives occupation readout
from the existing coin/current records: on eleven connected seams, a change
of coin exposes the original occupation contrasts after accounting for
recorded transfers. The fixed total then determines the initial configuration.
Its exact charge effects, state-independent timeout and terminal carrier are
derived from both C/E instruments; an occupation projector is not a new
primitive. This readout is destructive and does not identify coherences.

Complete census trials on one continuing carrier supply martingale confidence
for the three reference weights. Positive intervals admit a source-compiled
common modular clock and two continuing cap reads. The fixed protocol bounds
the **successful-preparation** history error, separately from herald failures
and statistical failure. It halts explicitly when the confidence box or source
certificate fails. The retained accelerated kernel demonstration uses `2^20`
preparations; every microscopic preparation attempt is not executed. These
finite effects and faithful-reference connections instantiate the specified
observer interfaces. Coherent enactment, physical service accuracy, laboratory
seconds and a common physical field dictionary remain separate obligations.

The [signed phase readout](native_modular_reference/README.md#infer-an-oriented-clock-parameter-from-complete-native-records)
connects that common clock to finite observer records. A second native core
read removes the unwanted frequency from a cap current; cosine/sine analyzers
then give an oriented interval without fitting the unknown reference
amplitude. Complete fixed-length batches carry the same native state forward,
retain failures and propagate preparation, source and sampling error through
an exact inverse. Both caps use the same inferred generator. The executable
certificate includes two repair-law alternatives, raw full-source histories,
an all-column dyadic sampler audit and complete attempt/service budgets.
The admitted parameter window resolves phase aliases; individual records are
not monotone. The compiler enacts a programmed modular parameter, so
this closes a readout/calibration interface without deriving autonomous
repair time or a laboratory duration.

The [public Bell-record comparison](predictions/BELL_RECORDS.md) tests the
scope of a classical interpretation. The diagonal native branch obeys
`S<=2` under independent settings and local factorized responses. Replaying
the deposited Hensen records gives 196 wins in 245 eligible trials and
`S=2.42250`; the published fixed-count calculation has a tail bound near
`0.03917`, conditional on its setting calibration and experimental access
contract. This retrospective null rejection excludes that extension, not
the full quantum C/E construction, and supplies no lab-clock identification
or a prospective validation of OPH.

The [measurement bounds](measurement_bounds/README.md),
[empirical diagnostics](empirical_diagnostics/README.md) and
[action-independent readout identities](readout_controls/README.md)
constrain the relation between a declared apparatus and data. They do not
supply an OPH-specific empirical confirmation. Missing raw counts are not
inferred, and calibration is not refitted from held-out outcomes.

### Supplied interacting work and quantum continuation

The [common interacting source](common_source/README.md) carries electric
flux, a charged CAR pair and two neutral scalar sites on one work register.
Its matter readout controls a feedback operation; the scalar readout controls
a native repair. The complete history retains 1,056 native/quantum records
and their terminal states. Nonzero feedback contrasts have finite and
[infinite-occupation certificates](common_source/INFINITE_HISTORY.md).
The action, initial state, coherent services and Hamiltonian model time are
declared inputs; compiling them does not identify native repair with that
Hamiltonian. The [theorem attachments](common_source/DERIVATION.md#compatible-theorem-use-and-open-physical-connections)
state the action and physical identification hypotheses.

The [mixed-state continuation](common_source/OPEN_HISTORY.md) implements
an explicitly supplied unobserved environment. Discarded phase draws occur
at declared boundaries; their labels cannot enter feedback or postselection.
It preserves visible histories, disturbed states, reference correlations and
checkpoint continuation. Environment variance, scheduling and physical
service errors remain inputs. The [typed protocol](common_source/COMPOSITION.md) and its
[runtime](common_source/protocol.py) execute smooth reads and feedback without
resetting the work register.
[Quantum native continuation](common_source/QUANTUM_NATIVE.md) and
[encoded recovery](common_source/ENCODED_MEMORY.md) attach full native
instruments and active logical memory under explicit code-service contracts.

The [combined recorded observer](common_source/ENCODED_MEMORY.md#optimal-native-carrier-code)
sharpens that memory interface: the full 78-state total-two carrier has maximum exact
record-correctable subspace dimension six across all thirty native seams.
The bound covers arbitrary complex codes and record-dependent CPTP recovery;
a six-matching source code attains it. Six corrected carriers then support a
finite fermion history with a declared three-cell transport word, two number
reads, record feedback, reference entanglement and all selector timeouts.
The same history supplies charged response, local continuity, cycle phases
and an information ledger: native correction preserves quantum information,
while the first field read loses exactly `(15/16) log(2)`. The separate
unencoded entropy identity extends the original native population reference
to coherent endpoint/reference outputs, including unresolved C timeouts.
The [combined receipt](common_source/record_observer_results.json) retains
these domains and complete complex field outputs together.

This completes the specified finite repair-and-transport construction, with
its code services, chosen logical ring and preparation exposed. Native records
alone leave the transport undetermined. The result supplies no physical clock
calibration, four-dimensional spacetime, Standard Model spectrum, gravitational
entropy identification or independent empirical confirmation; those remain
the completion criteria above.

The [protected Weyl–Higgs continuation](common_source/chiral_candidate/README.md#protected-propagation-at-nonzero-momentum)
uses the same code on 630 occupations of one fixed seven-momentum field.
All internal, Spin and family components remain in 96 independently checked
Fock factors. Two-duration number-read protocols and quantum-clock records
diagnose nonzero-momentum propagation along all three coordinate axes, with
feedback, reference outputs, selector timeouts and a separate full-service
budget. The spectral box, trial couplings and coherent services are supplied;
the [measured mass hierarchy fails](../docs/research/OPH_PARTICLE_MASS_RESPONSE.md).
Physical frequencies, momentum, interacting gauge/Higgs dynamics and a common
matter-plus-gravity state remain separate obligations.
The [second-shell test](common_source/chiral_candidate/dispersion_identifiability.py)
exhibits a local two-hop completion with exactly the same restricted
seven-mode histories and low-momentum derivative, yet a different
coordinate-momentum 2/3 joint detector law. Its additional lattice modes and
source resources remain explicit. The existing finite modes therefore do
not select a unique spatial law.

On the same retained matter representation, the [finite overlap
operator](common_source/overlap_regulator.py) checks an actual four-direction
Euclidean Wilson kernel, modified chiral projector and nonabelian gauge
covariance. A declared admissible-link domain gives a volume-independent
ideal locality tail, and complete native free-kernel queries retain both
quantum outcomes. The global fermion measure, interacting continuum and
Lorentzian many-body attachment are unestablished. The
[passive schedule classification](common_source/full_repair_response/stationary_selection.py)
also proves that, at the existing faithful biased reference, stationarity
excludes all eight core-to-outer seams and keeps the tested static Higgs
family seed rank at most two; active correction and different mass routes
are outside that result.

The [fixed active-cycle calculation](common_source/parameter_selection.py)
derives the full logical transfer $I_{36}$ and zero discrete response seed for the
existing recorded code and recovery. All one-port loads and native edge
costs compress to scalars. The original faithful reference is not stationary
under this active cycle on the same carrier. A two-hop weighted correlation
has nonzero occupation differences but leaks out of the code and is not an
existing coherent term in the fixed repair. This fixed cycle alone supplies
neither its coherent coupling nor a physical field dictionary. Specified
phase-control extensions have their own error and scheduling contracts.

The [specified coherent continuation](common_source/repair_continuation.py)
derives a controlled effective Hamiltonian for a two-hop phase followed by
the actual compiled correction. Directed rational contraction certificates
cover its full reachable matching space. The finite-proposal error includes
arbitrary references; numerical and physical per-step errors accumulate
separately. Its fixed-length outer schedule continues after a retained idle
result. A globally terminal timeout has a distinct stopped limit. No one
endpoint channel is treated as selecting an unwrapped Hamiltonian.

The [correlation/exchange field](common_source/repair_field_attachment.py)
uses a finite six-cycle average to preserve the code and positive exchange
between complete carriers to derive hopping with a density interaction.
Independent occupation modes, charge continuity, source readout and the
multiple-particle sectors use the same one-dimensional action. Its native
gaps have ascending ratios $1:3:12$. Twelve million simulated detector
records from a reduced encoded processor enclose those gaps without a mass
packet. Local native attachment controls are executed separately; a full
many-carrier service history is not claimed. The control architecture,
sector initialization and clock interpretation remain inputs, and the
field is not a physical Standard Model or a chiral continuum.

The [persistent geometric object](common_source/native_gauge_object.py)
starts with two incidence-defined, A5-covariant matching codes. Each carries
`1 + 5` and corrects every recorded native C/E branch. Two protected quintets
contain unique geometric `3 + 3′` with exact native intertwiners.
[Source preparation and transport](common_source/native_code_transport.py)
execute the opposite-code pair with padded registers and complete outputs.
Active recovery, coherent services and selection of that geometric input
domain remain explicit. No topological charge or fermion statistics follows
from the code or its antisymmetric color tensor.

The [same-object response](common_source/native_persistent_response.py)
derives port-exchange interactions on all 10,626 native load-four states and
solves leakage constraints before computing the current algebra. The tested
opposite-code response gives `su(3) + su(3) + u(1)` of dimension 17; the
complete 12-current premise of the finite Standard Model theorem is unmet.
The common phase is silent on the geometric sector's matrix algebra, while
coherent references in the retained complement require a larger interface.
Finite pulse approximation, continuous response accessibility and physical
gauge-current selection have distinct contracts. The tested four-carrier
interaction has a singlet block mixing both triplet references and leaks
from the geometric product before additional protection.

The [response comparison](common_source/native_response_comparison.py)
sharpens this mismatch on the same geometric frame: the exchange and
twelve-current laws intersect in eleven dimensions; their union generates
eighteen. Five exchange directions are observable quadrupoles on `3′`,
while the required relative central phase is absent. The common exchange
center cannot replace it. An exact Lie quotient cannot turn a simple
`su(3)` factor into `so(3)`. The
[quantum-channel reduction](common_source/native_response_reduction.py)
also shows that a unital completely positive map on the right qutrit
fixing all three spin-one rotation generators must fix its entire matrix
algebra. On the full six-state space this also requires preservation of the
right block projector. In the covariant
multipole family, erasing the quintet requires reducing the rotation
multiplier to at most `2/3`; this bound is attainable. These are precise
reduction obstructions, not a no-go for approximate emergent dynamics.

A [shared-reference construction](common_source/native_response_reference.py)
supplies a different operational contract. Add the native class-sum
relative-phase control and prepare a maximally entangled `3′` reference
on two geometric objects. The maximal collective current algebra preserving
that reference ray has dimension thirteen; quotienting its common phase
leaves exactly the existing twelve-current law. Without the relative-phase
control it leaves only eleven observable currents. The reference and its
preservation rule are supplied apparatus choices: the tested four-carrier
interaction selects neither this pure reference nor that restriction.
This constructs a measurable compatibility condition, not physical
inaccessibility of the excluded controls or a derivation of the whole
Standard Model.

The [biased-reference construction](common_source/native_response_biased.py)
then uses the actual four-carrier interaction and native class sum to form
`H_lambda = -B + lambda(C tensor I + I tensor C)`. For every rational
`lambda >= 1`, its complete 36-state geometric compression has a unique
mixed Bell ground with gap at least `8/75`. At `lambda = 16` the first
triplet weight is about `2.4113e-8`; it is positive at every finite bias.
An optional geometric penalty on the full 1,296-state color space has
explicit energy, leakage, ground-distance and uniform-time stability
bounds. At penalty `2^44` the ground-distance bound is below `4.47e-10`.
This is a separate exact Hamiltonian bound: that protective service and
the full Hamiltonian's finite pulse implementation are not executed.

The [source cooling stage](common_source/native_response_cooling.py)
prepares the projected ground on a fresh bath and uses an actual native
C record of probability `1/12` to trigger a SWAP. It retains all counter
records and reference/bath outputs; the old quantum information remains
in the bath. The unconditional reference marginal approaches the prepared
state with remaining weight `(11/12)^n`, also for entangled inputs.
The fresh baths and feedback coupling are supplied architecture, and
cooling prepares the projected reference rather than the exact protected
full-color ground. Reference monitoring or cooling alone leaves every
data current available.

The [disturbance controller](common_source/native_response_dynamics.py)
adds the missing coupling to those data currents. It computes the actual
reference-test Gram `G_p` in the complete eighteen-current frame and
applies `x <- (I - eta G_p)x` to every coefficient. No twelve-current menu
or rank cutoff enters this law. The full Gram spectrum is
`0` (seven), `4p/3` (five), `4(1-p)/3` (five), and `8p(1-p)/3` (one).
At fixed finite bias the exact zero-cost algebra has six observable
currents, and the finite implemented filter remains invertible. The
controlled joint limit `p -> 0`, `m -> infinity`, `mp -> 0` instead
converges to a bracket-closed thirteen-dimensional range, whose common-phase
quotient is the twelve-current law. This filters Hamiltonian coefficients;
it is not the prohibited quantum channel fixing rotations while erasing
their quadrupoles. The result requires bounded drives mediated by the
declared controller, not arbitrary rescaling or bypass of its attenuation.
With unrestricted input amplitudes the finite filter reaches all
eighteen currents. The class-sum phase is also an added work actuator;
using it only as a reference bias leaves eleven observable directions in
the filtered limit of the original seventeen-current algebra.

The [finite source history](common_source/native_response_controller.py)
uses the same native C trigger for bath replacement and controller updates.
It retains all 512 attempts, caps the update count at sixteen, and executes
all eighteen primitive data controls for every possible terminal count,
with subsequent geometric encoding and complete native C/E repair.
At bias sixteen and step `1/2`, the count timeout has exact binomial-tail
probability about `4.28715e-7`; it is retained rather than postselected.
For a unit-norm primitive at drive parameter `1/8`, the complete ideal
history is within `9.87e-7` in channel norm of the limiting control, with
matching bath/clock outputs and arbitrary work references. Numerical
source errors accumulate separately through all 512 cycles; physical
controller, storage and actuator errors remain additional services.
The Gram is computed from the source model, not estimated from detector
samples, and the native-event clock is not calibrated laboratory time.
This is an implemented finite architecture approaching the twelve-current
band, not autonomous selection by unchanged C/E repair or a proof of the
exact finite Standard Model theorem's complete-response premise.

The [specified current law](common_source/native_sm_current.py) chooses a
different coherent interaction on the same persistent object. Native port
coordinates generate a closed `u(3) + so(3)` law with all sixty proper
actions implemented internally. Its center is a relative phase between
the triplets, observable in the retained cross-triplet pointer records.
The source execution keeps the same C/E instruments and all raw outputs.
It establishes completeness over its declared current-control menu; the
universal compiler and detector services remain apparatus, so completeness
over all physical interventions and gauge redundancy do not follow.
The [matter construction](common_source/native_matter_attachment.py)
applies a specified weak spin lift and exterior-sector selection to these
same currents. It computes the fifteen internal matter weights and three
Higgs lines and intertwines them with the declared Weyl–Higgs field.
Physical spin/statistics, the family register and numerical Yukawa values
are distinct inputs.

The [free cubic extension](common_source/native_staggered_field.py) derives
oriented staggered hopping from an exact inverse-exchange echo and explicit
parity controls. Its Fourier reduction gives two Dirac copies per native
sector and a fixed-band continuum response with model-coordinate speed
$2J$. All 256 Fock columns of a small open cube, its half-filled source
preparation and a continuing four-outcome instrument are executed.
The larger periodic Fourier identity and growing service schedule have
separate algebraic and resource scopes. The free gaps retain ratios $1:3:12$.

A [composed two-carrier history](common_source/native_field_history.py)
instead keeps the actual nine-bit carriers through two exchange intervals,
detectors and complete C/E recovery. All four input columns and 319 terminal
records include raw aborts and padding; explicit tensor-compression errors
join the input/reference bound. These number reads test exchange and
backaction, with no sensitivity to a uniform onsite term at fixed charge.
The twelve-million-record mass readout remains a separate reduced-processor
experiment; physical proper-code services remain admitted.

The [interacting alternating path](common_source/native_charged_spectrum.py)
retains the exchange density term. Its exact two-cell charged poles are
$\sqrt{J^2+m^2}-J$, with nonzero CAR residues. Fixed-charge exclusion and a
path Poisson equation give the finite-volume second-order half-filled gap
$m^2/[2J(L-1)]$, with an explicit nonuniform remainder domain. These poles
follow from the declared action. Their ratios fail the lepton comparison
and no relativistic bulk mass is established. This
action differs from the echoed free cubic field.

The separate [calibrated matter action](common_source/calibrated_matter.py)
checks a common-scale empirical Standard Model parameter point on the same
internal tensors and actual current matrices. The complex CKM current and
charged mass coordinates pass exact checks. Its optional Weinberg sector
discloses the supplied neutrino splittings, mixing and unmeasured choices.
The [external matching replay](common_source/calibrated_matching.py) checks
calibration with SMDR; its output masses are not native observer records,
and its unfitted W value misses the selected one-standard-uncertainty
interval. The original native candidate contains none of these empirical
inputs. Source selection of masses, interacting evolution and a physical
measurement dictionary remain separate obligations.

The [moving calibrated field](common_source/calibrated_transport.py) uses
1,170 CAR modes on 13 momenta. Fixed preparation and three fixed Pauli
analyzers remove mass-dependent probe controls; adaptive interval rotation
uses earlier records classically. Rest and two momentum shells share the
same clock and supplied action. Complete simulated records, including
timeouts, give conditional energy and invariant-mass intervals; source
frames, spectator costs and physical-service admissions remain explicit.

The [calibrated observer](common_source/calibrated_spectrum_observer.py)
executes all four input columns of the 42 commuting rest-pair source factors,
with six neutral zero modes retained, and samples complete field/clock records.
It retains 84 million simulated trials across 42 settings, including selector
timeouts. The blind decoder reads only counts and declared error/clock contracts;
all nine supplied running gaps lie inside its intervals. Actual CKM frame words
and a complete 64-column weak down-sector comparison support the basis mapping;
the full 90-mode frame error remains an explicit uniform admission. The
conditional sampling bound is about 1.17e-8, with physical service accuracy
assumed separately. The large repeated-service budget is charged, not executed.
This is finite calibrated source/readout evidence; it supplies neither measured
pole precision nor parameter selection, interacting scattering or physical
protection.

The [native-clocked common observer](clocked_common_observer/README.md)
uses an actual current event to trigger that work's read-and-feedback suffix.
Work evolves during the wait. Its 128-slot execution retains 4,097 stopped
and timeout outputs, with complete finite-target error at most 0.000019,
including the specified binary arithmetic. All-input channel and tail bounds
control the finite composition. Work duration per attempt, the supplied
action and physical calibration remain separate inputs.

The [causality identities](common_source/CAUSALITY.md) control committed
records and permitted intervention order. Error bounds apply to their stated
instruments and preparation domains.
The [complete-output comparison](common_source/READOUT.md#complete-comparisons-with-the-actual-finite-pointer)
shows why agreement of one sign probability need not imply a close quantum
output: the sharp detector has a long occupation tail. The [occupation refinement](common_source/OCCUPATION_REFINEMENT.md) and
[smooth-detector refinement](common_source/SMOOTH_REFINEMENT.md) distinguish
their pointer domains. The
[smooth-history bound](common_source/SMOOTH_HISTORY.md) has a different
instrument and a specified preparation contract. The
[direct numerical certificate](common_source/OUTPUT_CERTIFICATE.md) and
[complete observer certificate](primitive_observer/COMPLETE_CERTIFICATE.md)
account for their actual executed arrays; neither is a physical service
certificate. [Observer readouts](common_source/OBSERVER_READOUTS.md) and
[thermal-clock comparison](common_source/THERMAL_CLOCK.md) retain the
declared reference state, access and clock assumptions.

### Effective field benchmarks and empirical readout

The [two-site source compiler](observer/source_scalar/README.md) executes
preparation, evolution and both detector outcomes from source words, with
complete occupation and splitting bounds. Its separate action and apparatus
do not inherit another spatial family’s continuum result.

The [shared scalar observations](scalar_observations/README.md) compare
complete field detector, preparation and feedback laws across a fixed-band
spatial family. Their [occupation proof](scalar_observations/OCCUPATION_DERIVATION.md)
and [receipt](scalar_observations/occupation_results.json) control projection
onto finite source codes. The exact-gradient endpoint is
`GalerkinParameters(grid=None)`: increasing quadrature resolution does not
add modes. Extending the field needs additional modes, matched observation
maps and controlled omitted-mode effects.

The [interacting band extension](interacting_band_extension/README.md)
performs one such extension, adding the cosine and sine at momentum
`(3,0,0)`. The common quartic action fixes their new couplings. Averaging
their specified oscillator blank induces a low-mode mass-squared shift
`g/4`; the comparison includes that shift instead of equating bare and
matched parameters. A block-Duhamel bound controls complete low-mode
observations through measurements and feedback, with the added modes
retained between reads and traced only at the comparison output. Its
finite occupation domain and selected anisotropic band remain explicit.
This does not remove all omitted modes or prove relativistic locality.
The [same-product quartic construction](interacting_band_extension/README.md#interacting-readouts-on-the-same-radial-product)
uses angular-constant profiles `sqrt(2) sin(r)` and `sqrt(2) sin(2r)` on the
declared interval `(0,pi)`. Their projected interaction is
`W=3(q1^4+q2^4)/2+6q1²q2²`, with full oscillator polynomials compressed
after multiplication. At occupation cutoff `4x4`, native source pulses
prepare the interacting Gibbs state and execute complete two-read feedback
histories jointly with its reference. The ideal state determines its modular
flow; raw preparation and propagation residuals remain separate. The
coupling and model-time intervals are declared. A native-incident protocol
advances this same Hamiltonian on zero-inclusive ticks, retaining every
timeout and disturbed work/reference output. Its refinement controls the
conversion to model time under the declared paused-service architecture.
The two-mode positive form has a faithful trace-class Gibbs state and
convergent complete adaptive histories as occupation grows. Canonical
reference embeddings and a common occupation/clock/precision schedule
are necessary. The d=4 to d=8 complete-output distance is about 0.56045;
the small executed code supplies no accurate unrestricted-output claim.
This projection is not an invariant spatial band or a proof of an
interacting ultraviolet continuum.

The [quantitative thermal certificate](interacting_band_extension/README.md#quantitative-preparation-and-complete-history-certificate)
supplies an explicit occupation budget for this same interacting
observer. Exact finite lower forms and a summable harmonic complement bound
the infinite Gibbs preparation. Directed full-generator path residuals
retain the reference, both reads and record feedback. The declared cold
`beta=4,d=8` history has complete ideal error below `0.0025`; its executed
source word has error below `0.005` under the stated source numerical
guards. The original hot `beta=1/4,d=4` benchmark retains its limitations.
The same fixed controller executes the cold preparation without repeating
it, retains all 3,600 native tapes and 4,900 terminal outcomes, and is
checked on all 64 work-input columns. One occupation bound covers every
completed zero-, one- and two-read timeout prefix. The separate four-level
central-weight calibration keeps its own state family, reference scope
and synthetic count examples; its coefficients are not transferred to the
cold protocol.

These premises compose on the same two-mode code: positivity gives a
faithful trace-class Gibbs state; that state supplies the full-algebra
modular flow; complete instruments preserve all records and references;
the occupation residual and incident-clock schedule supply separately
charged refinement errors. Alternative repair instruments, algebras or
preparations require their own attachment. This is not simultaneous
satisfaction of every finite theorem's alternative hypotheses.

The [coherence budget](native_interface/README.md#quantitative-preparation-floors-for-the-native-only-grammar)
also limits what repair can supply. Committed C/E repair with joint diagonal
accessible inputs and fresh independent diagonal ancillas cannot prepare
the declared coherent reference to arbitrary accuracy, even with adaptive
heralding. The successful-work lower bound survives as a complete-output
lower bound after charging all failure mass. Coherent source services
therefore remain a necessary resource for these particular targets; an
uncommitted coherent record dilation is a different operational grammar.

<a id="integrated-effective-description"></a>
### Integrated effective-description verdict

The native angular form, its independent radial-product construction and
its active observer protocol admit one consistent scalar effective theory.
With `B=-d_r^2+Abar`, its free thermal covariance is
`coth(beta sqrt(B)/2)/(2 sqrt(B))`, its retarded response is
`sin(t sqrt(B))/sqrt(B)`, and by the finite KMS identity the same Gibbs
state has modular parameter `s=-t/beta`. Only `Abar` is derived from native
repair; the radial interval, the qubit codes carrying the field, the compiled
per-mode Hamiltonians and the temperature are supplied. The chosen
hyperbolic product wave action has a Lorentz-signature coefficient density
on one time and three spatial coordinates. On a smooth
representative of the complete density this is the ordinary static scalar
field theory with the declared radial boundary condition. The full
coefficient density, including the measure, must be checked. The actual
radial-projection measure is continuous but not C1 across original face
edges in the fixed sphere atlas. Since `Ctt=-mu`, that exact dictionary
does not satisfy the smooth-density premise. Smooth comparison forms
approximate the measurable wave theory; a different smooth dictionary
would require its own justification. The reflecting radial boundary is
part of this theory. No boundary-free Hadamard or renormalized-stress
claim is inferred.

The actual product nevertheless admits smooth effective descriptions.
Angular smoothing of its positive Lipschitz measure and measurable elliptic
tensor, with the local unitary readout map `U_epsilon=sqrt(rho_epsilon/rho)`,
gives norm-resolvent convergence on the common Dirichlet cylinder. Thermal
field covariances and sine responses converge in operator norm; momentum
covariances converge on the Dirichlet form domain. For any fixed finite
Weyl/Ramsey/feedback menu over bounded wave times and any positive error
tolerance, one common smooth static `1+3` metric reproduces all its joint
probabilities and finite detector/reference outputs within that tolerance.
Source and incident-clock errors enter the same triangle budget. This is
an existence result for a fixed sequence of smooth descriptions, with no
data fitting, computable smoothing rate or exact smooth-native claim. Both
sides of the approximation are free fields assembled on the same supplied
codes; they differ in the smoothness of the angular coefficients.
The [proof and readout domain](relative_observer_geometry/README.md#lorentzian-density-of-the-product-wave)
retain the radial boundary, state and active-control assumptions.

Clock conversions retain their protocol: the free product controller uses
`t=(2d/sqrt(630))*tau`, while the interacting controller's word durations
`1/R_l` give `t_model=tau` in its limit. The common action is stated in
Hamiltonian time, with the corresponding coordinate pullback for each
controller. An unchanged physical conversion across the protocols is not
proved by counting their incidents.

The finite interacting protocol has a
[single fixed controller instrument](interacting_band_extension/README.md#one-fixed-instrument-controls-the-clock-and-field)
on its native timer, control registers and quantum work. It reproduces all
3,600 native tapes and 4,900 joint native/field outputs on the full sixteen-dimensional
valid work code with classical timer states, including untouched references.
A separate complete export control retains off-code failures. Its ready/busy
states enforce the carrier pause; timeouts and unfinished executions remain
explicit outputs. At most 775 controller instructions execute the finite
program after preparation; the source-pulse total includes preparation. This realizes its sequence and feedback inside a fixed
instrument; the controller architecture, initialization and protected source
services remain declared, and instruction count supplies no physical duration.
Committed-prefix intertwiners connect this same product program to the finite
persistence-observer interface: its actual first-write state, classical
prefix algebra and normalized block readback supply the observer fields,
with the declared evidence labels and buffer refinements. This does not
identify its registers with the separate native checkpoint carrier or
extend persistence to the creation of a record.

The [actual pointer and reusable source-control menu](common_source/DERIVATION.md#product-readout-derived-access)
generate `B(C8) direct_sum B(C8)`, with first-mode parity coherence
inaccessible. This derives algebraic work access and its effect-to-record
map; stable tomography and the reusable service menu remain separate.
The [complete recorded histories](interacting_band_extension/README.md#one-reference-clock-for-complete-intervened-histories)
use the fixed Gibbs modular parameter `s=-t/beta` of the compiled
Hamiltonian, by the finite KMS identity, with compiled errors retained. The same state fixes detector/Q-kick cross-noise and response;
the finite-kick check uses the existing fixed-interval source histories.
Central parity weights can change observations without changing this
restricted modular flow. With the Hamiltonian and inverse temperature fixed,
[complete native records identify that weight](common_source/DERIVATION.md#central-weight-from-complete-records)
within the sector-Gibbs family, with exact probability enclosures and explicit
sampling and service errors. The inferred interval
[predicts a separate elapsed-time readout](interacting_band_extension/README.md#a-held-out-record-law-after-central-weight-calibration)
without a second fitted parameter. The executed control uses the same Gibbs
preparation; illustrative counts are synthetic, and this short protocol has
weak state sensitivity. Neither the inference nor halt-state stationarity
selects the preparation or identifies a physical clock.

The same local quartic action has the two-mode Galerkin Hamiltonian above.
Its occupation limit is an interacting quantum theory whose finite
measurement histories converge, and the active native clock has a
compatible diagonal limit. This is a finite-spatial-band effective
description. It does not establish a continuum interacting ultraviolet
limit, a conformal theory, the Standard Model or Einstein dynamics.

Selection is a different mathematical question with an explicit obstruction:
three connected radial extensions have identical
[complete adaptive angular histories](relative_observer_geometry/README.md#radial-nonidentifiability-under-complete-histories)
on the constant-radial input code, retaining native marks, timeouts and
untouched references. Angular-only controls preserve this code; nonconstant
radial preparations and full-product Gibbs states are outside the equality.
A supplied radial-contrast probe distinguishes the extensions at finite time.
The native observer-centered shell quotient gives an operational radius
within the original carrier, and the axial current counter has its own
derived radial diffusion limit. Neither selects the separate radial-product
field or its boundary conditions. The shell quotient cannot be tensored back
with all angular modes and counted as an independently derived dimension.
The active wave protocol also differs from passive diffusive repair dynamics. Thus the
present native observations admit this conditional effective construction
without determining it uniquely. Physical calibration may fix units and
parameters, after which unused measurements must test the fixed dictionary;
synthetic parameter recovery alone is not such a test.

The [relativistic scalar regulator](relativistic_scalar/README.md) executes
an eight-site three-dimensional quartic action with an unrestricted-
oscillator comparison. Its ideal record TV error is below 0.000928 at the
stated finite parameters. The supplied positive mass gives a finite-box
gap bound, and the free single-particle lattice dispersion converges to the
relativistic dispersion on a fixed momentum band (analytic, infinite lattice,
no interaction; not a refinement of the executed eight-site cell).
An interacting relativistic limit requires additional spatial/energy control
and parameter matching. The [spatial interaction model](spatial_interaction/README.md)
has controlled growing-volume observation laws and a nonzero interaction
witness for a supplied repulsive Bose action. Growing a box at fixed spacing
does not remove its spatial regulator.

The [observer reconstruction](observer_reconstruction/README.md) recovers
free propagation and finite geometry from executed detector laws and a
specified reference-clock readout, with withheld predictions and controlled
occupation errors. Its quadratic action and mode apparatus are supplied.
A cube adjacency graph does not by itself determine physical dimension;
the [identifiability checks](identifiability_checks/README.md) include a
planar control and alternative actions invisible to a restricted apparatus.

The [published readout comparison](scalar_observations/experimental_readout/README.md)
calibrates a Gaussian response relation on zero-duration trapped-ion data;
100 withheld off-axis complex points, or 200 reported real and imaginary spin
means, test that axis-factorization relation with spin-mean RMSE 0.07639. The
[executable detector](scalar_observations/experimental_readout/EXECUTABLE_READOUT.md)
contracts source preparation, rotation and detector words under the declared
laboratory model; on the same 200 means its primary spin-mean RMSE is
0.06647, against a calibration residual of 0.07397 and 0.06440 for the best
unfitted nominal coordinate convention. Published coordinate conflicts, postselection and
normalization limitations remain explicit. This component postdiction
establishes neither the seven-mode dynamics nor a physical clock.
The [shared-preparation test](scalar_observations/experimental_readout/shared_calibration.py)
propagates exact outer calibration intervals to all 200 reserved means and
separately checks one axis-selected state against every calibration and
validation box. Both published coordinate conventions admit such a state,
with the same parameters used by the retained native source execution.
The conditional simultaneous allowance assumes 250 independent accepted
shots per component, the declared reported-mean rounding bound and zero
additional physical bias; the aggregate deposit
cannot authenticate these assumptions. This is a retrospective component
comparison, with no resolved coordinate convention or all-scale verdict.

The [occupation-control construction](occupation_control/README.md) bounds
the continuous polynomial-path defect, occupation boundary, endpoint jumps,
thermal initialization, displacement and pointer for the supplied pulse
apparatus. It includes arithmetic under its stated IEEE-754 contract, not
the separate native compiler or laboratory error. At 64 levels the case-a
estimator radius is below 9.91e-8 at 0.4 ms, 0.0163 at 2 ms, 0.211 at 5 ms
and 1.68 at 10 ms; the 32-level conditional bound is unavailable at the last
two durations. Agreement between cutoffs is not an error certificate.
The independently controlled unrestricted-versus-exact-finite radius can
compose with a valid native execution bound.

The [identification design](identification_design/README.md) computes the
response rank and nullspace for a declared correction family. Four-level
access cannot detect its `binom(N,4)` term; coherent levels four and five
expose those directions, apart from the global-energy gauge. A 41-pulse
source preparation and phase read distinguish two supplied actions, while
occupation commitment before the read erases the distinction. The efficient
sequence is to identify blind directions, choose a separating intervention,
then compile its locked measurement. A component sampling bound does not
transfer to a noncommuting field model or physical apparatus.

### Source execution and numerical contracts

Exact finite-code identities give controlled occupation parity in two native
pulses and dimension-64 number drift in seven. They do not import unrestricted
oscillator commutators at a cutoff. Parameter-keyed compilation caches share
source words, never observed trajectories, fitted parameters or prediction
outcomes. Canonical verifiers execute every declared calculation and
fail on any violated contract.

Pulse matrices use the fixed-guard polar representation of the source
library. Raw, incremental and corrected-target checks bound the representation
change; no state or measurement outcome chooses that correction. Density
matrices are not renormalized to repair execution. These numerical controls
are distinct from physical gate accuracy. Frozen evidence retains its
original bytes and scope; live receipts reproduce from canonical producers.

## Native repair and its observer contract

Use the twelve-port source, with 78 nonnegative integer configurations of
total load two. Its 66 binary configurations form an invariant shell. For
an oriented seam `(a,b)` and fair placement coin `c`, native C uses

```text
d = x_a − x_b,
j = [d − (2c−1)(|d| mod 2)]/2,
y = x + j(unit_b − unit_a).
```

The load is conserved. Its quantum current instrument has Kraus operators
`L_(c,j)=(1/sqrt(2)) Σ_(x:j(x,c)=j) |y(x,c)><x|`.
The record is `(c,j)`. C is one member of the admissible family in the
[record-law classification](record_laws/README.md). Adding a seam-parity readout gives E; a distinct coherent current-only
instrument gives C0; recording both endpoint loads gives T. Classically
discarding a recorded coin leaves C unchanged and does not give C0. All of them have identical classical transitions.
The classical rule fixes the diagonal of every member and none of its
coherences.

The [native observer](observer/README.md) first applies C on `(0,1)`. Its committed first record alone
selects `(0,5)` when `j1>0`, otherwise `(0,1)`, and supplies the preissued
forecast of a positive second current. No source loads or future outcomes
are accessible to that controller. Full simulator checkpoints and external
occupation audits are distinct from the observer's record access.

The [derivation](observer/DERIVATION.md) gives one fixed CPTP map on the
phase/tape/source carrier and shows classical-prefix preservation. Genuine
feedback can reduce coherence between different committed record values while
preserving their classical probabilities. Halting is an identity operation
on the terminal sector, not a claim of unlimited continuing memory.

Exact controls separate record information from knowledge of the chosen
action. Native-record prediction improves on the action-known prior and a
within-action shuffled record. The forecast is the enumerated law's own
conditional probability and the controller a design choice, so both
comparisons are consistency checks of the construction. Native feedback gives final port-0 occupation
`35/264`, versus `1/6` when the selector is disconnected or independently
shuffled with the same action frequencies. Occupation here is an external
audit; the native current remains the implemented observation.

## P-dependent response and finite capacity

For this design, **P and N are accepted as derived from OPH's self-reference
equations**. Apparatus, protocol and simulation hardware remain construction
choices. The [detector implementation](observer/p_response/README.md)
evaluates the existing forward-closure P input and source-only
structured-running map `α(P)`, then sets `g(P)=sqrt(4π α(P))` in natural units.
The observer does not fit P or N to its desired output.
Its chosen charge assignment is `Q_v=g(P)N_v`, so native current j transfers
charge `g(P)j` while conserving the total.

A prepared path qubit uses the first current to orient a calibrated loop
exposure χ. A Y measurement writes `y=±1` with

```text
Pr(y | j1,P,χ) = [1 + y sin(g(P)j1χ)]/2,
E[g(P)j2] = (5g(P)/264) sin(g(P)χ).
```

Only y chooses the second seam. At χ=1 the mean native second current is
approximately `0.00564883566843`. Disconnecting the detector or independently
shuffling y removes the mean response while preserving its half/half outcome
and seam frequencies. Weak exposure gives
`E[g j2]=(5π/66) α(P)χ + O(χ³)`: P affects an actual response magnitude.
The loop geometry can separately use `A_loop=K_loop P ell_star²`.

The chosen charge, prepared probe, calibrated flux, readout and controller
are specified apparatus. This is a software execution of that conditional
design, not a report of fabricated electromagnetic hardware. The inherited
source-only α branch retains its documented physical endpoint status;
executing the adapter does not turn its output into an independent precision
measurement. The [endpoint check](predictions/endpoint_audit.py) retains
the approximately 300 ppm inverse-alpha deficit of that incomplete
Thomson-endpoint output, whose hadronic part is the declared quark-screening
ansatz `S = 1 - N_c alpha_3/pi`; a zero-correction precision endpoint is
therefore unsupported. Replacing the ansatz with the measured hadronic
payload gives inverse alpha `136.3828` on `[136.3670, 136.3985]`, outside
the measured value; with source anchor and lepton transport fixed, the
bookkeeping attributes the remaining gap to the frozen anchor. Neither result
rejects other OPH constructions or fits P to the target. P does not alone
determine a whole observer's Hilbert dimension,
local exposure, temperature or event rate.

The [capacity calculation](observer/CAPACITY.md) reuses the existing
correctable-public-record theorem. Under the binary preparation, 11 terminal
native words and 16 reachable prefixes allow a four-bit classical prefix
encoding. On the full 78-state source those counts are 17 and 24. These are
reachable record counts, not the dimension of the whole machine. The
P-detector's richer audit has 24 terminal words and 37 prefixes; its actuator
receives a single bit. Program, source, probe and service capacity are
additional.

For a declared capacity carrier of dimension D, `N=ln D` bounds a correctable
classical record code by `ln M0≤N`. The cosmological relation
`N=3π/(Λ ell_star²)` needs the stated physical capacity/area identification.
A finite controller does not establish global saturation or cosmic closure.
See the existing [capacity implementation](../code/capacity_readback/README.md)
and [N-closure theorem](../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex).

## Geometry, response and local event clocks

The native current process supplies a finite port geometry and dimensionless
clock relation under its declared uniform binary reference and IID seam
schedule. These are readouts of the same recorded events, with a finite port
chart and attempt parameter whose physical identification requires an
independent correspondence.
The [repair-law guide](CANDIDATES.md) retains the exact baseline relations.

N2 lets the native current drive the existing local six-level response;
N3 is its simpler flat alternative. They preserve the full native load
instrument, so internal preparations cannot feed back on that marginal.
N4 instead conjugates binary C/E by a coherent load–response exchange.
It preserves locality and total occupation and supplies finite reciprocal
feedback, with an added coupling θ and specified internal preparation.
N4's record label j is not a bare-load current.

The [central-address extension](../code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair)
places the full quantum algebra in twelve identical blocks. Refreshing the
central address and repairing a uniformly selected incident seam reproduces
the entire original channel after decoding. With its declared collective
response, chart implementation, scalar seams and all-level refinement, this
gives a complete abstract A1–A3 source on the binary shell. Reversible
response words and chart paths have control sequences fixed independently
of measured records; feedback uses the declared repair/reset templates.
Allowing records to select reversible words would add 144 independent
address-controlled directions and change this source grammar. The minimum
twelve-block dimension assumes a unital embedding of the full working
matrix algebra and exactly twelve primitive central ports; it is not a
universal physical encoding bound. It does not extend N4 to all 78 loads.

Native and N4 versions share the full tracial quantum reference, collective
response and normalized occupation geometry while their finite kinetic
invariants differ. These requirements therefore allow both implementations;
uniqueness is unnecessary for constructing either. A central control address
is not an observed occupied port. For the prescribed occupied-port effects,
every instrument erases coherences between disjoint load supports, with
recovery error at least `1/2` on an explicit orthogonal pair. Endpoint
comparison consumes a separate record service; scalar seams carry no remote
record. The explicit occupied-port endpoint instrument recovers the
certified global-attempt statistic under the tracial preparation; adding
intermediate occupation reads changes that experiment. Physical realization and service costs remain additional contracts.

The [native-response calculation](../code/native_response/README.md) uses
one common tracial state and fixed bare occupations to recover the same
port Gram geometry and a certified two-attempt memory excess for N4.
For one local occupation under global IID attempts, the corrected statistic
is `b_local=c2−c1²−(1−c1)²/10`. The subtraction removes ordinary spatial-mode
mixing and gives a zero native null. All outcomes are retained and there is
no intermediate measurement or reset.

The [local-tick extension](local_tick/README.md) counts every attempted seam
incident to the chosen port, including unchanged loads. All intervening
remote attempts remain in the evolution. Under uniform IID seam selection
the incident probability is `p=1/6`; the first-hit Schrödinger channel is
`K=B(I−A)^−1` for remote/incident channel sums A and B. It requires a different
correction:

```text
r = 1−c1,
b_tick = c2−c1²−r²(209r+180)/(1331r²+3300r+1800).
```

Native repair has `c1=5/11`, `c2=93/407`, `b_tick=0`. Reusing the existing
edge-Gram certificate `G≥mI` proves `b_tick≥3m/50` for N4. Certified lower
bounds are `0.00001119970056` at θ=1/16 and `0.00022670580642` at θ=1/8.
These are lower bounds, not numerical evaluations of the full N4 correlation.
C and E share c1 and these bounds; the local-tick result also gives
`c2_C≥c2_E`, without asserting strictness.

Preparation, occupation readout, incident-event recognition and calibration
remain specified services. The current-law detector does not implement this
N4 occupation apparatus. Local attempt count, model time, host runtime,
modular parameter and calibrated observer duration are distinct. The
[working terms](../docs/INTERPRETATION_AND_TERMS.md) impose no preferred
universal physical clock.

<a id="retained-scalar-and-entropy-constructions"></a>

## Scalar and entropy constructions

### Executed scalar detector and native continuation

The [scalar adapter](observer/scalar_response/README.md) prepares the original
finite-action vacuum with signed momentum mean `4j f`, where the first native
current is `j∈{−1,0,1}` and f is the inherited left-slab profile. The supplied
positive scalar Hamiltonian evolves this field. A binary Weyl instrument
reads the inherited detector profile at model time `t∈[19/20,1]`, retains
its conditional output state and commits `y=±1`. Only y chooses the second
native seam. The first current and the detector bit remain readable after
feedback. Conditional states use finite superpositions of Weyl-displaced
coherent states, with their coefficients and Gram overlaps retained; a full
Fock-space vector need not be allocated. The instrument acts on the field;
it is not a draw from a stored theoretical receipt or execution of a compiled
native quantum-control word.

Write `r_q(t)=Pr(y=+1|j=1)−1/2` at refinement q. Symmetric signed preparation
gives the native continuation statistic

```text
E[j2] = 5 r_q(t)/132.
```

The comparison uses the existing golden refinements `q=13,34,89,233`, the same
preparation and detector profiles, and the same model-time window. If the
parent scalar bound is `|r_q−r_cont|≤ε_q`, the complete labelled outcome-history
law differs by total variation at most `5ε_q/33`, and the second-current mean
differs by at most `5ε_q/132`. These transfer bounds concern this specified
observable and protocol; they do not assert convergence of the full Fock
state or an interacting theory. The first two grids retain the parent's
trivial error bound; q=89 and q=233 give nontrivial bounds. The analytic parent
establishes the limit independently of these four numerical executions.

The [local evolution](observer/scalar_response/local_evolution.py) executes
the canonical kick–drift–kick recurrence with the local tensor stencil at
q=13 and q=34. It advances the mean and detector pullback columns; the
original-vacuum covariance then supplies the squeezed detector variance.
This local receipt performs one field read and two native repairs.
The q=13 step-halving control compares at the same actual model time. The
separate [receipt](observer/scalar_response/local_results.json) retains
checkpoint summaries and sampled native records. It is not an authenticated
log of every field update or a compiled native quantum-control word.

The existing time-refinement certificate adds a temporal detector error to
the spatial one on all four grids, including the change in the original
vacuum under split evolution. The same `5/33` and `5/132` factors transfer
that combined bound to complete histories and the mean current. A nearest
grid-time comparison also includes the parent's timing bound. Only q=13
and q=34 execute the local loop here; their inherited spatial bound is
trivial. The tighter q=89 and q=233 spatial/temporal bounds remain analytic
comparisons for this extension. Numerical agreement with the spectral
oracle is a software check, not a certified floating-point error bound.

The [prediction protocol](observer/scalar_response/PREDICTION.md) fixes the
feedback statistic with disconnected, zero-preparation and sign-reversal
controls. Independent mathematical reconstruction and seeded execution test
the software prediction. A physical comparison additionally needs a realized
field, apparatus, calibration and independently acquired records. The scalar
action, quantization, vacuum, preparation, Born readout and model clock remain
explicit inputs; the finite native observer does not select them.

### Field and record continuation

The [continued scalar observer](observer/scalar_response/CONTINUATION.md)
retains one field through two noncommuting reads at model times `19/20` and
`1`, with three repairs of the same native loads. All 2,112 labelled histories
are kept. Its 209 reachable record prefixes have an explicit eight-bit
encoding; source, field and service memory are separate. Four coherent Weyl
terms retain the terminal quantum state, while every one of the 12,487,168
q=233 modes contributes to their Gram kernel.

This extension uses continuous full-mode spectral evolution. The local
split-update engine above and this two-read conditional-state continuation
are separate implementations with separate error bounds.
The two-site source compiler executes a separate single-read apparatus;
compiling this two-read history requires its own implementation and bounds.

The directed rational bound puts the q=233 complete-history TV error
below 0.02571. It includes the zero-current histories, whose two-bit
correlation is not controlled by the single-read bias bound. It also
certifies a positive second-read disturbance above 0.000022281162 relative
to omitting the first read on the positive-current preparation. A commuting
Gaussian control retains ordinary temporal correlation but removes that
disturbance. Independent state-vector, conditional-field donor, reset and
record-deletion controls distinguish these mechanisms. The existing supplied
action, vacuum, instrument and time contract is preserved.

### Scalar action and entropy premises

The separate scalar benchmark factors weighted agreement into its supplied
stiffness. If repair capacity also measures inertia, linear dynamics, energy
conservation, time reversal and the declared velocity reading give
`M f_ddot=−(K_D+m²M)f` and the canonical Hamiltonian. These hypotheses determine
the equation within that class. The `I+εA` kinetic counterfamily has the same
normalized classical configuration Gibbs marginal with different frequencies;
the full phase-space reference and quantum vacuum change. Three-mode curvature
and a distance-two response distinguish it after common clock/mass/speed
refits. Nonzero deformation violates the same-capacity inertia premise;
configuration equilibrium alone does not determine inertia.

Local coordinates `χ=(m q+i p,Dq)` represent the original canonical fields
without independent edge oscillators. The existing free-scalar continuum,
local force and retained-pointer instruments remain conditional results;
interactions and physical matter identification are further work. See the
[source instruments](../code/source_scalar_instruments/README.md) and
[finite implementation](../code/source_scalar_finite_instrument/README.md).

The retained q=3 eight-mode preparation/history uses the same 8,192-level
onsite code, with combined mathematical trace-distance error below `0.001237`
before physical service errors. Its ramp error is below `0.000433450` at
model duration 512; `2^24` local time steps add below `0.000695044`.
Existing source controls give a finite synthesis recipe, and exact commuting
compilation removes an unnecessary inner approximation loop. These are
construction/error results; no enormous native pulse word or 106-qubit
field/pointer state was executed. The [premise map](PREMISES.md#re-10)
retains the detailed preparation, history and resource assumptions.

Regional maximum-entropy completion, nested inclusion towers, protected
records and charge-resolved collars remain available from the
[finite algebra results](../code/quantum_information/ALGEBRAS_AND_REPAIR.md).
The [native mixing theorem](observer/QUANTUM_MIXING.md) makes the
unread binary fixed-number state converge to its tracial entropy maximum under
the specified independent schedule. At `I/d`, entropy has zero first
variation for every trace-zero state perturbation on that fixed algebra.
The [joint counting identity](../code/geometry/ENTROPY_FIRST_LAW_AUDIT.md#joint-counting-stationarity-and-its-geometric-boundary)
`S_count(p,d)=log(M)-D(p||d/M)`, with `M=sum d`, also makes joint state/count
variations stationary at saturation when `M` stays fixed, in the stated
positive-real relaxation. The native sector counts satisfy that finite
identity. The [full return reference](../code/geometry/ENTROPY_FIRST_LAW_AUDIT.md#native-reference-entropy-and-region-boundaries)
has entropy deficit at most `2025/(64V)` for `V` ports uniformly over caps; its regional
mutual information grows at most logarithmically in `V`, excluding a
nonzero coefficient times crossing-seam count on a refining base-face
patch. This is a two-dimensional regional boundary, not whole-screen area.
Fixed total count and `log(d)` have no established physical volume/area
map. A charge center supplies no physical gauge-edge factor or stress
coupling. Nonnested completions need not commute. The
[gravity premise ladder](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs) retains
the physical assumptions needed for the existing gravitational conclusions.

## Connecting finite objects to the physical theorem routes

The same word can name different objects in these constructions. Use the
actual carrier, operation and accessible record when applying an existing
theorem. The [generated map](PREMISES.md) links their executable sources and
proofs; the [observation ledger](../docs/registers/OBSERVATION_LEDGER_V3.md)
records the broader physical targets without promoting their status from a
finite software run.

| Constructed object | Result that can be reused | Required connection for a broader claim |
| --- | --- | --- |
| Common flux-dressed CAR/scalar register and two recorded instruments | Source execution, operator Gauss, committed records, finite temporal and occupation bounds, exact sharp-readout contrasts, state-specific complete-output cutoff enclosures and unavoidable occupation-code output errors | Spatial refinement, vacuum/coupling identification, numerical/physical source services and responsive stress/geometry on this same process; post-sharp clock continuation needs its own error analysis |
| Observer representations and appended clocks on the finite sharp history | GNS and local-reference modular readouts, fixed full-work Gibbs/modular correspondence, source-compiled clock probes and CP regional restrictions | Declared observer access, factorization and references; reference selection by repair, physical clock calibration, growing spatial family and occupation limits of the appended clocks |
| Separate smooth Ramsey detector on the common action | Quantitative complete-output occupation bound, finite-energy backreaction and a separate prepared four-read history certificate | Different instrument and bit dictionary; the uniform moment bound is vacuous at small executed cutoffs, and the sharp-history clock does not transfer |
| C source, phase and finite record tape | All-input CPTP law, exact forecast and feedback controls; every finite `OperationalObserver` clause on the actual postcommit carrier, with classical access and declared blank-buffer refinement | Record creation is outside the persistence generator. Typed prediction retains the same evidence; next-current forecasting is separate. A compiled Lean inhabitant, cross-observer overlaps and physical identification are not supplied |
| P-dependent path-qubit detector | Complementary binary effects, branch charge covariance and bit-controlled native feedback under accepted closure inputs | Physical realization of the declared apparatus; a same-quantity screen/Thomson readback and precision endpoint claim keep their existing contracts |
| Native-current scalar preparation, binary field instrument and native continuation | Executed conditional instrument, feedback mean `5r_q/132`, inherited history and mean refinement bounds | Source identification of the scalar action, vacuum, preparation/readout and clock; a physically calibrated independent detector comparison |
| Reachable native record words on the declared carrier | Existing correctable-record evaluator and explicit classical prefix encoding | A compatible physical capacity cut and global closure family before identifying this count with cosmological N |
| Twelve central addresses with binary native/N4 quantum blocks | Full-channel decoding, complete abstract scalar-seam grammar, shared full reference and normalized geometry with distinct finite kinetics | Physical realization, retained-record delivery and compatible apparatus; the address is not an occupation measurement |
| N4 tracial bare occupation at incident attempts | Corrected first-hit null and inherited strict finite memory bounds | N4 occupation preparation, detector and incident-event recognition; the current-law detector does not measure this occupation |
| Supplied canonical scalar action, vacuum and Weyl observables | Existing finite instrument and controlled free-field continuum theorems | Explicit maps from the chosen native implementation with controlled state, observable and service errors along the same refinement family |
| Native C/E endpoint maps and their federation extensions | Local Liouville action, exact instrument blocks and complete quantum observation locality | Compatible screen refinement, physical clock/field correspondence and controlled relativistic scaling |
| Joint local clock and field instruments | Exact record-conditioned quantum memory; completed-event work/reference measure-and-prepare closure | Declared support, schedule, preparation and joint apparatus; unresolved waits must remain in the quantum history |
| Retained matter representation on a finite Euclidean overlap lattice | Gauge-covariant modified-chiral kinetic projector, one free zero species and volume-independent ideal local approximation on a small-link domain | Global chiral measure, quantum gauge/Higgs action, Lorentzian continuation, full Fock execution and source-selected masses |
| Native-angular Gaussian slab vacuum with variable spatial metric | Same-action bare regional entropy/area obstruction and a port-local stress degeneracy control | Renormalized small-ball entropy, physical stress/area map, isotropic ultraviolet refinement and responsive Einstein geometry |
| Coherent overlapping observers on the native total-two carrier | Exact complete-output correspondence for a transported matrix cone, including arbitrary untouched references and number-preserving union-supported quantum interventions | Declared coherent preparation and apparatus; the inverse encoding is not positive on arbitrary effective matrices, and the correspondence supplies no relativistic field dynamics |
| Interacting scalar projections | A seven-to-nine-mode matching bound; a same-product two-mode Gibbs theory with convergent unrestricted-occupation histories; a cold eight-level source/controller history with a complete occupation budget and a native-incident clock schedule | Fixed spatial band, radial design, coupling/temperature, protected services and physical units; the original hot four-level code retains its substantial thermal allowance, and neither construction supplies an ultraviolet completion |
| Continuing native preparation and likelihood test | Uniform finite preparation and dependent-transcript decision bounds | Physical apparatus mapping and authenticated conditional service errors; the finite checks are not an executed laboratory campaign |
| Declared gauge, edge, stress and geometry structures | Existing conditional gauge/gravity implications | Joint physical carrier, interactions, conserved stress and compatible geometric variations on that family |

No observation rung or physical premise is discharged merely by putting these
sources in the same repository. Conversely, the remaining physical joins do
not reopen the completed finite observer contract. Model maintenance belongs
here; changes to the broader papers follow their separate review policy.

<a id="retained-audit-controls"></a>

## Counterexamples and identifiability

The [fixed audit package](../code/programme_obstruction_audit/README.md)
retains exact author and independent controls for the following boundaries.
Run `python3 code/programme_obstruction_audit/run_all.py` from the repository
root with assertions enabled. Its seven scripts validate the frozen inputs
and write fresh receipts outside the package. These controls reject the
named identifications or overly broad inferences; none rejects A1–A3 or
requires a unique designed implementation.

| Audit issue | Exact retained boundary | Construction consequence |
| --- | --- | --- |
| [#1087: complete carrier](https://github.com/FloatingPragma/observer-patch-holography/issues/1087) | `Z(M66)=Z(M78)=C I` cannot supply twelve nonzero primitive central ports. Binary occupations are noncentral, overlap and sum to `2I`; total-two occupations also admit eigenvalue two. | Keep the native matrix algebra as a sector. The twelve-block central-address extension supplies the missing central algebra and real port-controlled operations; occupation measurement remains a separate contract. |
| [#1088: marked triplets](https://github.com/FloatingPragma/observer-patch-holography/issues/1088) | The slow Gram satisfies `G²=4G`, `tr G=12`, `(5I−A)G=(5−sqrt(5))G`. One fixed order-five action has slow-frame trace `(1−sqrt(5))/2` and lower response-triplet trace `(1+sqrt(5))/2`, excluding an invertible same-action intertwiner. | Keep spatial and internal response domains distinct. Whole-law Galois substitution or an explicit action twist needs its own attachment; incidence recharting alone does not identify the marked triplets. |
| [#1089: feedback class](https://github.com/FloatingPragma/observer-patch-holography/issues/1089) | At fixed load `I_odd/2`, the rational number-preserving exchange gives endpoint probabilities `769/1250` and `481/1250`, difference `144/625`. Four Kraus operators certify CPTP; its private commutant has dimension one, excluding every single-product exponential representation, including alternate logarithm branches. | Preserve the all-input reduced-load and single-product no-feedback theorems within their stated classes. The escape is a finite channel control; importing geometry, current or clock semantics requires a common observable map. |
| [#1090: configuration equilibrium](https://github.com/FloatingPragma/observer-patch-holography/issues/1090) | For `A=I+L_path3`, the positive `ε=±1/16` controls retain the normalized configuration Gibbs marginal but change the full reference and vacuum. For stiffness response `Y=γ(A+εA²)`, `Y13=γε` and the third-mode fit residual is `6γε`; actual acceleration is `−Yq`. | The same-capacity inertia theorem retains its premise. Separate this configuration-only control from the native/N4 counterpair, which fixes the full quantum reference and normalized geometry but changes finite kinetics. Neither supplies two physical realizations. |

The historical #1088 control uses `closed_paths[1]` of its frozen response
fixture. The [marked-triplet certificate](../code/a5_closure/marked_triplet_attachment_certificate.py)
transports a different frozen order-five element and has slow trace
`(1+sqrt(5))/2` and lower trace `(1−sqrt(5))/2`. Both exclude the same-action
identification; character signs must be compared at the same marking.

The [independent joint reconstruction](../code/programme_obstruction_audit/joint_independent.py)
checks the marked path and occupation identities; the
[independent channel reconstruction](../code/programme_obstruction_audit/channel_independent.py)
checks all 64 operator matrix units and solves the private commutant.
The [kinetic control](../code/programme_obstruction_audit/kinetic_control.py)
retains all eight prescribed rational cases, including rejection of singular
`ε=−1/4`. The source extension does not erase these negative controls.

## Remaining construction work

The missing joins concern one compatible family of observer observations,
not an inventory of additional finite theorems:

- Establish compatible federation neighborhoods, states, accessible algebras
  and instruments across the intended scaling family. The conditional
  quasilocal quantum construction has these hypotheses; the existing screen
  refinement does not discharge them.
- Connect the derived native local action and event-clock histories to an
  effective interacting field/spacetime description through complete joint
  observation and intervention laws. Preserve quantum memory where a reduced
  field is not Markovian. Modular parameters and event records need a checked
  correspondence on that process, not a preferred universal clock.
- Identify the physical sector, apparatus and parameters with independently
  calibrated measurements. Specify data ancestry, permitted fits and units,
  then lock predictions and rejection rules before comparison. Source control
  words implement a supplied action without selecting its physical meaning.
- Realize the one-event interference apparatus that reads the record law's
  coherence `g`, authenticate its full service budget and collect complete
  records. Finite native preparation and a dependence bound replace exact
  reset and IID premises; physical controls and a laboratory mapping are
  open.
- Connect interacting matter, conserved stress and responsive geometry to the
  same observations and observer-relative clock comparisons. Fixed-volume
  occupation control and growing volume at fixed spacing address different
  limits; neither alone removes a spatial regulator.
- Account for long-running record storage, coherent services and fault
  budgets while preserving finite protocol behavior and references. Checked
  simulator restoration is not a physical reset or autonomous memory service.
  N4 additionally requires its stated occupation preparation/readout and
  event-recognition apparatus.

Useful tests target actual alternatives. The
[binary coherent continuation](observer/COHERENT_CONTINUATION.md) implements
a reversible 3,960-arc lift and a compressed-amplitude recurrence. Different
continuations or cross-record phases give distinct Ramsey predictions while
agreeing on the committed native instrument. Its eleven-mode density sector
has an exact 720-arc representation and a
[114,910-pulse source word](observer/COHERENT_CONTINUATION.md#source-compilation-of-the-finite-density-walk).
All reversible basis inputs and clean-scratch address outputs are checked;
the full-register replay uses raw pulse-executed routines and retains
padding. The complete numerical operator guard is `s/(1-s)`,
`s=829/5000000000`; it is not an interval or physical certificate.
The full 3,960-arc compiler, walk-controlled Ramsey source word and growing
refinement compiler are unproved. The polarized fluctuation-pointer execution
has a separate, completed finite scope. A selecting restriction must constrain
the additional coherent interface freedom. The
[identifiability checks](identifiability_checks/README.md) likewise separate
apparatus restrictions, unseen action directions and geometry ambiguities.
The stronger effective-Hamiltonian comparator remains relevant to empirical
protocol design; missing timing, normalization or raw counts cannot be
inferred by fitting away residuals.

The [mathematical reuse guide](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md)
requires compatibility defined independently of the update and a mapping of
the actual domain, records and scheduler to the nearest theorem. Complete
history, instantaneous fixed state, stationary distribution and continuing
operation are distinct claims. A self-consistent history can contain changing
records and interventions.

## Source selection: the three primitive obligations

The [source-selection guide](SOURCE_SELECTION.md) organizes source selection
around carrier architecture, observer interface and implementation, and the
complete repair rule. For the repair rule the selection problem is sharp:
the classical seam rule is forced by four stated requirements, the record law
ranges over an exactly classified convex family, and one interference
experiment reads one of its coherence multipliers directly. The
[repair-law guide](CANDIDATES.md) states the principle that singles out each
distinguished member. Architecture, observer implementation and repair law
may supply the instance inputs. Action, state, scheduling and control consequences require
explicit derivations from those inputs or a separate declared status; renaming a
supplied input does not discharge it. The existing finite field programs
retain their comparison scope.

Selection is up to the observations admitted by the interface. Scalar phases
on committed Kraus branches cancel through every finite adaptive recorded
history. They become observable when uncommitted record interference is
admitted. The existing source words execute a fixed read with `3/4` versus
`1/4` plus probability for two phase implementations of the same committed
C/E instrument; the commit-first control gives `1/2` for both. The distinct
record dictionary is explicit. Coherent continuation alternatives also
survive all proper carrier rotations, real orthogonality, native support and
a common algebraic time reversal. A source-selection criterion must exclude an
actual alternative under the relevant interface, rather than infer a
preferred wave from these shared structural properties.

## Work ownership and claim accounting

`model/` is the live home of observer and repair-law construction. Its
implementation, tests, candidate comparison and this guide evolve together.
GitHub issues own tasks; the [claim registry](../claims/claim_registry.yaml)
owns established claim scope; [PREMISES.md](PREMISES.md) is generated
navigation, not a proof that every theorem antecedent has been discharged.

The [completion program #1134](https://github.com/FloatingPragma/observer-patch-holography/issues/1134)
and its [dependency DAG](../docs/COMPLETION_PROGRAM.md) own the linked
construction gates, including theorem/model integration and active-library
triage. Earlier coordination issues #1025, #1026, #1076 and #1077 are
superseded; their evidence and unfinished benchmark contracts are retained.
[Correctness #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
and the [broad obstruction search #1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
remain ongoing. Some research branches can proceed in parallel; acceptance
of a joined result requires its mapped predecessors.
Finite results retain their stated validity unless evidence invalidates
them; larger claims add declared requirements. A negative result rejects its stated family.

## Research programme and paper-release hold

Publish tested model code and reproducible results in this directory. The
[paper set](../docs/PAPER_SET.md) records the model milestone around which
the prose papers are organized and the role of each paper. **A paper
release follows the owner's approval of a recorded model milestone.**
Audits, necessary correctness fixes, errata, preview builds and
reproducibility checks continue between releases. Frozen evidence and
deposited releases retain their original bytes. Passing CI or publishing the
model code does not itself release a paper.
