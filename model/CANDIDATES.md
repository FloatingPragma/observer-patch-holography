# Repair laws: what is fixed, what is free and how to tell

**Find the simplest shared repair rule that satisfies the known constraints.**

An OPH world is a network of overlapping patches that keep repairing their
disagreements. The [runnable model](README.md) needs one concrete repair law
to execute. This page sets out what such a law has to specify, which parts
the stated requirements force, which parts they leave free, which members of
the free family the model executes, and which observations can tell the
members apart. The [construction guide](CONSTRUCTION.md) describes the
executable model built on the chosen law, and the [premise map](PREMISES.md)
connects both to the theorem library.

The [current physical audit](PHYSICAL_AUDIT.md) connects the single-carrier
geometry, action and clock results. None needs a federation for its finite
derivation. The [physical identification obligations](SOURCE_SELECTION.md#current-bridge-to-physical-observations)
separate these executed laws from an observer-relative QFT or spacetime
description and from the [actual public-data comparisons](../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md).

## What a repair law specifies

A repair law answers four questions about one event on one seam:

1. **Where does the load go?** This is the classical seam rule.
2. **What does the event write down?** This is the record law, a quantum
   instrument whose outcomes are the records.
3. **What happens to coherent records that nobody reads?** This is the
   continuation and record-access contract.
4. **Does anything besides the load respond?** This is an optional private
   response factor at each port.

Mathematically, write `H` for the span of the load configurations `|x>`. A
law assigns to each seam `e` an instrument `{Φ_{e,r}}`, completely positive
maps indexed by the record `r`, whose sum `Φ_e` is a channel. Its classical
rule is the diagonal action `Φ_e(|x><x|)`. Its record law is the whole
instrument. A scheduler chooses the seams. Continuation and access fix which
further operations may act on unread records, and response extensions enlarge
`H` by private factors.

| Layer | What it fixes | Status | Section |
| --- | --- | --- | --- |
| Carrier | Twelve ports, thirty seams, the total-two load sector | Declared by the architecture premise | [Carrier](#the-carrier) |
| Classical seam rule | Where the load goes | Characterized by four declared requirements; on total two they fix only the tie-break coin | [Classical rule](#the-classical-seam-rule) |
| Record law | Which coherences survive an event | A convex family with five coherence multipliers under declared local number-phase and endpoint-swap covariance; C executed as the convention, E as the comparison | [Record laws](#record-laws) |
| Continuation and access | What unread records do next | Free beyond the committed instrument; exact countercontrols | [Continuation](#coherent-continuation-and-record-access) |
| Response | Private factors driven by the current | N2, N3 and N4 compared | [Response](#response-extensions-n2-n3-and-n4) |
| Field benchmark | A supplied scalar action | Conditional benchmark S1 | [S1](#s1-the-canonical-scalar-benchmark) |

Each layer follows the
[mathematical reuse contract](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md):
define the agreement target independently of the update, state the domain,
scheduler and fault assumptions, and map each theorem to the actual objects.
A heat-flow theorem does not turn arbitrary quantum instruments or
conservative evolution into that heat flow. The
[source-selection contract](SOURCE_SELECTION.md) separates supplied instance
inputs, derived consequences and selection up to admitted observations.

## The carrier

Use the twelve ports and thirty seams of the source icosahedron. The tested
load space is the 78 nonnegative integer configurations with total load two.
Its 66 binary configurations, with at most one load per port, form an
invariant shell. Total two is the smallest sector in which a single load can
wait or move and two loads can collide, so the binary shell carries motion
and the collision states test the integer rule. The architecture, the sector
and the initial disturbance are supplied. The integer rule leaves an exactly
balanced configuration fixed; it does not explain the origin of a departure
from balance.

## The classical seam rule

For an oriented seam `e=(a,b)` and a fair coin `c`:

```text
d = x_a − x_b,
j = [d − (2c−1)(|d| mod 2)]/2,
y = x + j(unit_b − unit_a).
```

The two endpoints become the floor and the ceiling of their mean. When their
total is odd the coin decides which endpoint receives the extra unit. All
other loads are unchanged. The operation conserves the total and never
increases `sum_p x_p²`. A difference of one permits a wait or a swap, so the
binary shell keeps moving rather than settling pathwise. Uniform independent
seam attempts are the baseline schedule; prescribed seams are diagnostic
controls.

**Proposition 1 (the classical rule is characterized).** Let a stochastic seam rule
(i) conserve the seam total, (ii) end at a configuration minimizing
`x_a² + x_b²`, (iii) break any remaining tie with a random choice that does
not depend on the loads, and (iv) be symmetric under exchange of the two
endpoints. Then it is the rule above with a fair coin.

*Proof.* By (i) and (ii) the outputs are the floor and ceiling of the mean.
An even total has one output. An odd total has two, and by (iii) endpoint
`a` receives the ceiling with a constant probability `q`. Exchange of the
endpoints maps `q` to `1−q`, so (iv) gives `q=1/2`. ∎

Each requirement does work. Without (ii), partial moves and idle rules are
allowed. Without (iii), the rule "the endpoint that has more keeps the extra
unit" satisfies the others; it freezes every binary configuration,
so the shell carries no motion and the records carry no geometry. Requirement
(iii) is the basis rule applied to the tie: nothing recovered prefers one
endpoint, so the choice is maximally random. The four requirements are
premises of the construction. On the executed total-two carrier the seam
difference is at most two, so every move-toward-balance rule coincides with
floor and ceiling there and the proposition fixes only the tie-break coin.

On classical diagonal states of the binary shell the rule is the lazy
symmetric simple exclusion process with two particles; the twelve collision
states of the full carrier are transient. Its generator per seam attempt is
`(1/60) Σ_e(S_e−I)`. The exclusion-process corollary of
[Caputo, Liggett and Richthammer, J. Amer. Math. Soc. 23 (2010) 831-851, §4.1.1, equation (4.1)](https://arxiv.org/abs/0906.1238)
gives the exact classical gap `(5−sqrt(5))/60 ≈ 0.046066` on this carrier,
equal to the one-particle gap. Under uniform binary equilibrium and uniform
seam attempts the native records recover the carrier Laplacian, its
rank-three slow Gram and a dimensionless rate relation. Those readouts do not
identify physical space or proper time. See the
[native update](CONSTRUCTION.md#native-repair-and-its-observer-contract) and
[geometry and clock readouts](CONSTRUCTION.md#geometry-response-and-local-event-clocks).

## Record laws

The classical rule fixes where load goes. A quantum event can also keep or
destroy the coherence between configurations that move the same way, and
which coherences survive is decided by what the event writes down. Recording
more destroys more. The [record-law classification](record_laws/README.md)
works this out exactly on the 78-state carrier.

### The admissible family

Require of a one-seam channel `Φ`: seam locality, the exact classical law,
local number-phase covariance (`Φ` commutes with `exp(i Σ_p t_p N_p)` for
every real `t`) and endpoint-swap covariance. Number-phase covariance is an
additional condition beyond classical conservation. It guarantees an
unread-channel representation with definite-current Kraus maps. Both
covariances are declared premises of the construction: number-phase
covariance is a local U(1) phase symmetry of the repair channel at every
port, a local charge-superselection rule, and the classification below, the
C/E contrasts and every discriminating prediction rest on it. Requiring
the actually committed records to resolve those currents is a further
instrument contract, checked by the
[record-basis control](record_laws/INSTRUMENT_CONTRACT.md). Then

```text
Φ(|x><x'|) = Σ_j G_j(κ(x),κ(x')) |x + jD><x' + jD|,
```

where `κ(x)=(x_a,x_b)` is the seam key, `D` moves one load from `a` to `b`,
and each `G_j` is a positive semidefinite matrix whose diagonal is the
classical probability of current `j`. Every such table defines an admissible
law. On the total-two carrier five multipliers remain:

| Multiplier | Coherence it keeps | Bound |
| --- | --- | --- |
| `k` | empty seam with full seam | `\|k\| ≤ 1` |
| `g` | empty seam with a waiting single load | `\|g\| ≤ 1/√2` |
| `h` | full seam with a waiting single load | `\|h\| ≤ 1/√2`; `h=g` when `k=1` |
| `m` | the two waiting orientations | real, `\|m\| ≤ 1/2` |
| `n` | a moving single load with a moving collision | `\|n\| ≤ 1/√2` |

A law that records the endpoint's pre-repair load, `(c,x_a)`, keeps the
classical rule exactly and fails covariance: one of its branches carries two
different currents. The classification retains it as a control.

### Four distinguished laws

| Law | Each event records | Values per seam | `k` | `g` | Idempotent | Principle that singles it out |
| --- | --- | ---: | ---: | ---: | --- | --- |
| C0 | the current | 3 | 1 | `1/√2` | no | Record only what conservation requires |
| C | coin and current | 6 | 1 | `1/2` | no | The tie-break is itself a recorded event; the coarsest such law |
| E | coin and pre-repair contrast | 10 | 1 | 0 | yes | One repair completes the seam; the coarsest idempotent law |
| T | coin and both endpoint loads | 12 | 0 | 0 | yes | The seam occupation is recorded in full |

Each principle in the last column singles out one law; none is derived
from the basis, and the construction executes C as its convention.

They form a ladder. Each step records one more classical quantity through a
diagonal readout of the repaired seam:

```text
C = C0 followed by a readout of the coin,
E = C followed by a readout of the seam parity,
T = C followed by a readout of the seam total.
```

The [encoded source realization](observer/source_discrimination/README.md)
implements the second step as one three-pulse CNOT after the shared repair
circuit. Forgetting a label afterwards does not undo its dephasing, so C is
not a classical coarse-graining of E.

Repeating the repair of one seam multiplies each coherence by its multiplier
again: after `r` repairs the table is `(k^r, g^r, h^r, m^r, n conj(h)^(r−1))`.
Hence exactly E and T are idempotent, every law with `k=1` completes to E,
and every law with `|k|<1` completes to T. For C this is

```text
E_e² = E_e,   E_e C_e = C_e E_e = E_e,
C_e^n = E_e + 2^(1−n)(C_e−E_e),   n ≥ 1.
```

On the binary shell E is a trace-preserving conditional expectation. The
[finite entropy-completion theorem](../code/quantum_information/ALGEBRAS_AND_REPAIR.md)
requires a specified retained algebra and a faithful modular-compatible
reference, and its [binary native attachment](primitive_observer/DERIVATION_COMPLETION.md)
identifies E as the resulting tracial expectation. Completion therefore
explains why E is the endpoint of every coherent-seam law; it does not choose
the finite-event law.

### What the executed laws share and where they differ

The model executes C as its convention among the classified laws and E as the comparison law. Older
documents and retained receipts call C the N1 current law and E the N1-E
contrast law. C0 and T are classified with exact one-seam predictions; every
admissible table has the same network gap `(5−√5)/60`, and C0's slowest
binary coherence factor is `43/60 + 2√2/15 ≈ 0.905` per attempt.

| Quantity | C | E |
| --- | --- | --- |
| Classical load laws, histories and geometry | identical | identical |
| Full channel gap on the shell and all 78 states ([quantum mixing](observer/QUANTUM_MIXING.md)) | `(5−√5)/60` | `(5−√5)/60` |
| Slowest binary coherence, factor per attempt | `5/6` | `7/10` |
| Forgotten-record channel entanglement-breaking on the binary shell | by 134 attempts | by 134 attempts |
| Identify the 78-state input from complete records ([record recovery](observer/RECORD_RECOVERY.md)) | no bounded record-only policy; an almost-surely terminating one | eight prescribed attempts, seven on the shell |
| Joint clock field-plus probability after three attempts ([joint clocks](joint_clock_process/README.md)) | `35/54` | `14/27` |
| Amplitude of the modular-phase response ([modular dynamics](relative_modular_dynamics/README.md)) | `1/140` | `1/220` |
| One-event interference at phases `0` and `π` | `5/8`, `1/8` | `3/8`, `3/8` |

On the full 78-state carrier both channels retain provable reference
entanglement at every finite attempt count for a collision/binary input, with
a signal tending to zero, so replacing the carrier by its attracting shell
changes an exact finite-time quantum property. Both laws admit an exactly
correctable encoded qubit after one event. Every surviving matrix unit keeps
its signed load difference, which reduces each binary block to exclusion on
the remaining vertices plus a constant loss.

### What observations decide

**Common classical load/current records do not select C/E.** Every admissible
channel has the same diagonal action. For C and E the [observability results](identifiability_checks/README.md)
prove more: for every load basis state the joint law of coin, current and
post-repair load agrees once E's extra label is summed out, and the
equality continues through finite adaptive classical protocols and through
quantum work controlled by common records. Any admissible classical record
dataset has likelihood ratio `L_C/L_E=1` wherever defined.

**One coherent event reads a multiplier.** Prepare
`(|X> + e^{iθ}|Z>)/√2`, with loads at ports 2 and 3 in `X` and at ports 1
and 2 in `Z`, repair seam `(0,1)` once and measure the projector onto
`(|X>+|Z>)/√2`. Summed over records,

```text
Pr(plus) = 3/8 + Re(g e^{−iθ})/2.
```

| Law | `θ=0` | `θ=π` | `θ=π/2` |
| --- | ---: | ---: | ---: |
| C0 | `3/8+√2/4 ≈ 0.729` | `3/8−√2/4 ≈ 0.021` | `3/8` |
| C | `5/8` | `1/8` | `3/8` |
| E and T | `3/8` | `3/8` | `3/8` |

The settings `0`, `π` and `π/2` read `Re g` and `Im g`, so the experiment
measures the coherence of whichever law nature uses. Replacing the empty
seam's partner by a full seam gives `Pr(plus)=1/2+Re(k e^{−iθ})/2`, which
separates T from E. The [fixed prediction protocol](predictions/README.md)
locks the first test with 8,192 attempted trials per setting, every joint
record bin and an apparatus error budget of `1/64`. The
[continuing-carrier protocol](closed_prediction_protocol/README.md) runs it
on one continuing carrier with 379 discarded preparation attempts before
each of 3,927 trials and an exact-source wrong-choice bound of `3/400`.
The [law-family warm-up bound](record_laws/README.md#preparation-by-finite-native-warm-up)
also applies to every admitted current law, including complex phases: 221
binary or 318 full-carrier attempts give stationary-state error at most
`1/1000`, including an untouched reference, after discarding the warm-up tape.

Passive mixing also forbids a nontrivial modular automorphism on an invariant
observer algebra. The [active event clock](transfer_modular_probe/README.md)
adds a fixed coherent step per incident event and a declared pause schedule;
closed histories supply its faithful three-mode reference. This construction
does not select a quantum record law or identify an atomic clock.

**No retained dataset reads `g`.** The public oscillator and Bell records
have no identified native apparatus map, so they supply no likelihood for any
record law; the [empirical status](../docs/EMPIRICAL_STATUS.md) gives the
details. Choosing among the laws therefore takes either a principle from the
table above or a measurement of `g`.

### Complete current-instrument identification

The [nine-quadrature identification](repair_identification/README.md) extends
the C/E readout to all five multipliers, including their imaginary parts and
the collision coherence n. It identifies the entire classified coarse
current instrument on 78 states and bounds its complete adaptive-history
error by seven times the largest joint-probability error per repair.

A compiled control with `(k,g,h,m,n)=(1,0,0,0,1/4)` agrees with E on every
binary-shell input operator, so binary-sector protocols cannot exclude it.
The collision probe separates it from E with probabilities 1/2 versus 3/8.
Its extra private coin/parity records and costs differ; the equality concerns
the common current interface. Exact inversion, all 36 source probes and ten
continuing histories are reproduced by the packet. Neither successful
tomography nor synthetic data selects a physical law or its finer record
access contract.

## Coherent continuation and record access

A committed instrument does not determine what unread coherent records do
afterwards. The [coherent native continuation](observer/COHERENT_CONTINUATION.md)
embeds the binary shell's events in a 3,960-arc workspace. Its reflection
choice gives a wave recurrence for compressed amplitudes, and a second
continuation with the same first coherent event gives different Ramsey
predictions. Both survive all sixty proper carrier rotations, real
orthogonality, native transition support, a fixed uniform coherent state and
a common algebraic time reversal. A phase on every moving arc preserves every
committed C/E branch and changes the compressed trace cube by `−1/90`.

The [source-executed record test](observer/source_selection.py) realizes the
access distinction with existing source words: a complementary read before
commitment has plus probabilities `3/4` and `1/4` for two phase
implementations of the same committed instrument, while committing first
gives `1/2` for both. A selection argument for this layer must constrain the
cross-record phases and the continuation action; recovering a wave after
choosing the reflection does not.

## Response extensions: N2, N3 and N4

These candidates keep the native load repair and add a private six-level
factor at each port, with the existing skew-Hermitian response matrices `D_p`
from the [source construction](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex).
The response space is `(C⁶)^⊗12`. The matrices, their incidence-preserving
port identification, the twelve factors and the coupling `θ` are declared
inputs.

| ID | Law | Strongest established property | Principal missing connection |
| --- | --- | --- | --- |
| **N2** | The native current drives the same local response at both endpoints | Preserves the native load instrument and locality; inherits completion and finite binary-network convergence; protects local block populations | Full source grammar and selection against N3; all-input marginal preservation excludes feedback |
| **N3** | The current drives opposite, centred endpoint increments | Same preservation and locality; one global change of coordinates reduces it to the native law | A reason to select or reject its trivial closed-loop response |
| **N4** | One coherent load–response exchange conjugates the binary C/E laws | Private preparations change bare load probabilities; twelve central addresses preserve its quantum channel in a complete abstract scalar-seam grammar | Physical realization, preparation/readout, retained-record service and clock comparison; its label `j` is not a bare-load current |

### N2: the current drives the endpoint response

Replace each native branch by

```text
U_(e,j) = exp[θ j(D_b−D_a)],
K_r = L_r ⊗ U_(e,j)^(a) ⊗ U_(e,j)^(b),
```

with identity on the other factors. Both endpoints receive the same unitary.

| Property | Established scope |
| --- | --- |
| Load and records | For every correlated input, tracing out the response preserves C and each classical record block exactly. Coherences between distinct retained records can change. |
| Locality and symmetry | Disjoint seams commute; simultaneous proper port transport and response transport give finite carrier covariance. |
| Completion | A seam-dependent controlled unitary conjugates the update to C_e tensor identity. The same conjugation transfers E and the exact single-seam completion formula. |
| Network convergence | On the finite binary shell, fixed positive seam mixtures converge to the conditional expectation onto their common fixed algebra. If γ_C and γ_E are their gaps off that algebra, γ_E/2 ≤ γ_C ≤ γ_E. Its full algebra and absolute gap are uncomputed. |
| Persistent information | Each local projector diag(1,1,1,0,0,0) is protected. Arbitrary private coherences are not claimed to survive. |
| Response histories | Closed load loops can leave internal memory. For an input that factors into load and product private states, the record-conditioned response factors across ports and creates no private-to-private entanglement. |
| Backreaction | None: changing the internal response cannot change unconditional native load probabilities, by exact marginal preservation for all inputs. |

The convergence result concerns forgotten-record channels; it does not give
convergence of a growing coherent history or unlimited record storage.
Conserved current differences omit the constant central response direction,
so attaching `D_b−D_a` does not derive the full gauge structure or a
physical gauge field. Loop memory alone does not establish non-Abelian
curvature: a commuting-matrix control produces loop memory too. Exact native
source binding, covariance, locality and completion identities have finite
matrix checks and a complete three-event comparison with N3; see the
[detailed join](CONSTRUCTION.md#geometry-response-and-local-event-clocks).
A reciprocal extension must relax exact preservation of C on specified
internal inputs, because requiring it for every input rules out feedback.

### N3: the local flat alternative

With `J=Σ_p D_p` and `Dbar_p=D_p−J/12`:

```text
K_flat,r = L_r ⊗ exp(−θ j Dbar_a)^(a) ⊗ exp(+θ j Dbar_b)^(b).
```

For `W(x)=⊗_p exp(θ x_p Dbar_p)`, the same load-controlled unitary conjugates
every seam to C_e tensor identity. Every closed load loop therefore has
identity response and the full encoded internal algebra survives. Load
preservation, classical records, locality and proper covariance hold as for
N2, with no response-to-load feedback.

One local six-level state, one seam attempt and one local two-outcome effect
distinguish N2 and N3 with probability difference **0.049095 at θ=1/8**,
including both fair coin outcomes. A reference-assisted input that retains
all eight coin outcomes of a three-event history gives trace distance about
**0.137390 at θ=1/8**, about **0.069358 at θ=1/16**, and zero at θ=0, with no
postselection. Preparation, effect and coupling are supplied. Preservation
and convergence alone cannot select N2 over this simpler alternative.

### N4: a reciprocal coherent exchange

On the binary shell let `S_e` swap the endpoint occupations. With
`X=S_e(N_a−N_b)²`, `Y=iS_e(N_a−N_b)` and the same private factors:

```text
g=D_b−D_a, h=D_a+D_b−J/6, J=Σ_p D_p,
G=g^(a)+g^(b), H=h^(a)+h^(b),
W_e=exp[θ(Y⊗G+X⊗H)],     K_r=W_e(L_r⊗I)W_e†.
```

The exponent equals `θ(s_+⊗B−s_−⊗B†)` with `s_+=|10><01|` and `B=H−iG`: one
coherent exchange with no state variable beyond N2. The coupling, the equal
relative weighting, the exponential law and the quantum interpretation are
supplied. Both response directions omit the constant central generator.

**Why test it.** Any single-product conjugation `exp(θ A⊗G)` with local
number-preserving `A` leaves C's load marginal unchanged, which blocks
feedback for that whole family. N4 evades this scoped obstruction, and
neither summand alone has the tested feedback. The independent
[rational feedback control](CONSTRUCTION.md#retained-audit-controls) gives a
number-preserving CPTP escape with probability difference `144/625`; its
scalar-only private commutant excludes every single-product exponential
representation. N4 is a reciprocity experiment, and its fixed coefficient
ratio is an assumption even without a fitted parameter.

**Properties and result.** N4 is local, total-occupation conserving and
proper-carrier covariant, with protected 3+3 block populations. Same-seam
completion and finite binary-network convergence transfer by conjugacy; its
full common algebra and absolute gap are uncomputed, and it is restricted to
the binary shell. With fixed bare load `I_2/2`, supplied private preparations
change the endpoint occupation probability from **0.482020673 to 0.517979327
at θ=1/8**; the one-term controls give `1/2`.

**Main cost.** The record label `j` refers to the conjugated occupation, so
the native current-derived geometry and clock identities cannot simply be
imported. The [shared bare-occupation protocol](CONSTRUCTION.md#geometry-response-and-local-event-clocks)
reconstructs `K=αL_ico/396` at the full tracial reference with `α≥1/9` for
every coupling, retains the rank-three port geometry and certifies nonzero
two-attempt memory for every nonconstant carrier mode at θ=1/8 and 1/16. A
change of clock units cannot remove it. The
[central-address construction](../code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair)
preserves the full native quantum channel in twelve identical quantum blocks
and supplies a complete abstract scalar-seam grammar on the binary shell;
sampling an occupied port has a proved quantum information cost. The
[local-tick extension](local_tick/README.md) derives the null when the
observer counts incident attempts with all remote events retained, and keeps
a strictly positive certified N4 distinction at both couplings.

## S1: the canonical scalar benchmark

Weighted seam agreement supplies a scalar stiffness. Given positive repair
capacities `M`, grounded stiffness `K_D`, mass `m` and conjugate momentum
`Π`, the tested conservative completion is

```text
H = Πᵀ M^−1 Π/2 + fᵀ(K_D+m²M)f/2,
f_dot = M^−1 Π,       Π_dot = −(K_D+m²M)f.
```

Using the same capacity as inertia is a hypothesis. With the stated
linearity, velocity, energy and time-reversal premises it fixes this
completion, which oscillates rather than settling. The
[construction](CONSTRUCTION.md#retained-scalar-and-entropy-constructions)
reuses the free scalar continuum and quantum-instrument results under a
supplied canonical representation, state, controls and clock. An eight-site
construction composes preparation and retained observations on one finite
code with mathematical trace-distance error **<0.001237** before physical
service errors. This is the original q=3, d=8192 construction recipe: its
full source word has not been emitted or executed. The later two-site and
native histories below are executed on their own stated domains.

The [configuration-equilibrium control](CONSTRUCTION.md#retained-audit-controls)
preserves only the normalized classical configuration Gibbs marginal; its
full phase-space reference and quantum vacuum change, and nonzero kinetic
deformation violates the same-capacity inertia premise. The
[executed scalar adapter](observer/scalar_response/README.md) lets the first
native current prepare a signed field displacement, performs a binary Weyl
read and feeds only its result back: the mean second current is `5r_q/132`,
and the scalar probability error `ε_q` bounds the complete history distance
by `5ε_q/33`. The [two-site source implementation](observer/source_scalar/README.md)
compiles and executes a separate finite scalar action on the twelve-port
processor with both pointer outcomes retained. S1 supplies field benchmarks;
its dynamics are not identified with the native repair. The open connection
is to identify inertia and the Hamiltonian model parameter with observations
and compatible observer-relative clock readings from the same microscopic
process. The retained [kinetic alternatives](CONSTRUCTION.md#retained-scalar-and-entropy-constructions)
provide dynamical discriminators.

## Runnable observer protocols

These protocols run the native law C with an explicit finite record tape,
controller and stopping rule on the uniform 66-state binary source, and the
native implementation also checks the full 78-state extension. Two events
means two successive repair invocations.

The [native observer](observer/README.md) records `(c1,j1)`, issues a
supplied forecast from that record, chooses seam `(0,5)` when `j1>0` and
`(0,1)` otherwise, records the second outcome and halts. All 264 binary-source
paths are enumerated. The [capacity calculation](observer/CAPACITY.md) gives
11 terminal record words and 16 reachable prefixes, a four-bit classical
prefix encoding. The observer's action reads only the current, so it runs
unchanged for instruments satisfying the declared current-record interface;
its forecast also requires the coin. Unread-channel covariance alone does
not supply either actually committed record.

The [electromagnetic detector protocol](observer/p_response/README.md)
accepts the closure coordinate P (N enters no formula of the detector) and
computes the declared map's source-only Thomson-endpoint value `α(P)`. Its
charge assignment is `Q_v=g(P)N_v` with `g(P)=sqrt(4π α(P))`, and its
detector law is a supplied Aharonov–Bohm phase. A prepared path qubit and
calibrated loop exposure χ
produce a committed bit with `Pr(y|j1,P,χ)=[1+y sin(g(P)j1χ)]/2`; only y
reaches the next-seam controller. All 528 weighted histories are retained.
The response `E[g(P)j2]=(5g(P)/264)sin(g(P)χ)` vanishes under disconnected
and independently shuffled detector controls. Charge realization, flux
source, probe preparation and readout are specified hardware choices.

Both constructions preserve committed classical record values through
feedback and permit disturbance of coherences between record values. Neither
is an online learner or a complete physical universe. The detector reads the
native current; it is not the bare-occupation instrument that N4 requires.

## Excluded alternatives

| Alternative | Requirement it fails | Evidence |
| --- | --- | --- |
| Tie kept by the endpoint that has more | Load-independent tie-break; freezes the binary shell | Proposition 1 |
| Record of an endpoint's pre-repair load | Local number-phase covariance | [Record-law classification](record_laws/README.md) |
| N2 with one shared six-level processor | Disjoint-seam commutation | Declared locality requirement |
| Any single-product exponential response conjugation | Cannot feed internal state back to the load | One-term obstruction and the `144/625` control |
| Kinetic deformation `I+εA` of S1 | Same-capacity inertia | [Retained audit controls](CONSTRUCTION.md#retained-audit-controls) |

An odd-total integer seam reaches a nearest-balanced shell rather than exact
pathwise equality, and the full-carrier channel is not the binary conditional
expectation. The [identifiability checks](identifiability_checks/README.md)
and [construction counterexamples](CONSTRUCTION.md#retained-audit-controls)
state further domains excluded by exact controls.

## Maintaining this register

Use discriminating requirements, test deletion controls, count independent
assumptions and try the smallest diagnostic first, following the
[construction obligations](CONSTRUCTION.md#remaining-construction-work). A
shorter rewritten formula does not by itself reduce a law's physical inputs.
Keep each entry consistent with its law, state space, parameters, schedule,
record policy and physical interpretation, and give a scientifically distinct
alternative its own identifier. For each claimed recovery or agreement
result, record the legitimate set, preserved data, perturbation domain,
schedule class, nearest primary theorem, checked premise mapping and remaining
physical discriminator. The [bounded work products](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md#7-work-products-and-scientific-contribution)
fit the owners below.

[CP-00 #1135](https://github.com/FloatingPragma/observer-patch-holography/issues/1135)
owns the declared law/program boundary, and
[CP-01 #1136](https://github.com/FloatingPragma/observer-patch-holography/issues/1136)
owns complete source execution. They replace the coordination in #1076 and
#1077. The original q=3 scalar preparation and five retained-pointer
histories remain an unexecuted benchmark. CP-00/CP-13 must explicitly decide
whether to complete or retire it from the active path without relabeling
its status. Its large common-code word is unexecuted. Later complete
native and two-site programs establish their own finite milestones and do
not silently replace that deliverable. No entry closes a physical premise by
renaming it. Completion requires the simplest shared rule to satisfy the
required constraints on one common model, with discriminating evidence for
its physical identification.
