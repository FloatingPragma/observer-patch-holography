# A local action for native quantum observations

The native C/E instruments fix an interacting local generator on the binary
operator algebra. No field coupling or Gaussian action needs to be fitted to
obtain that generator. This packet also bounds how much an observer's **whole
local quantum and record output** changes when distant ports are omitted.
The bound includes arbitrary untouched references and interventions controlled
by the observer's own local records.

On one carrier, it also derives a current fluctuation action in one port's
recorded movement clock: the coefficient is `5/2`, with no fitted time rate.
The exact finite current blocks have correlated squares and non-Gaussian
fourth moments; the quadratic action describes their proved scaling limit.

It also identifies finite compatible neighborhoods of the actual spherical
screen source, without enlarging the elementary carrier. This adds a quantum
and complete-history extension to the existing
[density-only propagation bound](../federated_observer/SCALING.md). It reuses
the [signed-difference theorem](../observer/QUANTUM_MIXING.md) and
[actual federation instruments](../federated_quantum_sector/README.md);
the general block formula is reused. The attachment supplies
a local action representation, an observational truncation theorem and a
conditional local-limit construction for the full quantum process.

## What the declared inputs supply

There are binary factors at the ports of a finite graph, and the same
native endpoint repair acts on each seam. C and E remain two specified
candidate instruments; this calculation does not select between them.
The existing architecture supplies the twelve-port carriers, their seam
wiring and the allowed quantum interface. Uniform, independent seam
selection is a declared scheduling protocol. Coherent preparation and local
phase interventions in the finite controls are declared services.

For a graph with `M` seams, let `Phi_e` be one native repair with its records
discarded, and let `Phi=(1/M) sum_e Phi_e`. These maps are constructed from
the source endpoint function and the actual C/E record equality. Their
extension by identity on other ports is exactly the existing finite tensor
split. A federation needs no larger elementary carrier.

## The source fixes a local Liouville action

On the Hilbert space of operators with the Hilbert--Schmidt inner product,
define

```text
h_e = I-Phi_e,                 H_native = sum_e h_e,
Phi = I-H_native/M.
```

The endpoint calculation constructs orthogonal projectors such that
`h_e=P_1+(1/2)P_half` for C and `h_e=P_1` for E. Each term is positive and
acts only on the two endpoint operator factors. These facts hold on the
full binary tensor algebra, including coherences between number sectors.
The identity is a common zero vector, so the local operator is frustration
free. Its kernel need not be one dimensional: vacuum/full coherence remains
fixed on the unrestricted algebra. The fixed-number uniqueness
theorem has a narrower domain.

The signed-difference decomposition makes the action explicit. For a matrix
unit `|x><y|`, write `z=x-y`, and `Z={v:z_v=0}`. The common occupied sites
form a subset `S` of `Z`. On each fixed common-number block,

```text
H_native[z] = (1/2) L_exclusion(G[Z]) + V_L(z) I,
V_C(z) = (1/2) sum_edges |z_a-z_b|,
V_E(z) = sum_edges 1[z_a != z_b].
```

Thus its quadratic form on block amplitudes `f(S)` is

```text
<f,H_native[z]f>
 = (1/4) sum_(S,ab in G[Z]) |f(S)-f(S^ab)|^2
   + V_L(z) sum_S |f(S)|^2.
```

An exchange of equal occupations contributes zero. These coefficients come
from the native fair placement coin and record rule. The hard occupation
constraint and seam exchange couple distinct ports: for example the local
generator takes the population `|10><10|` to half the difference between
that population and `|01><01|`. It cannot be a sum of independent on-site
updates. Coherence differences additionally carry the explicitly derived
C/E killing terms. A Gaussian field with an independently chosen potential
is not inserted.

`H_native` is a positive **Liouville** operator acting on operators. It is
not a Hamiltonian acting on the source's state vectors. In particular the
native CP dynamics mixes some pure states and is not unitary. For the
original discrete scheduler, the exact transfer remains `I-H_native/M`;
replacing it by `exp(-H_native/M)` would change predictions. The logarithm
of a transfer, when defined, also differs from this local sum.

## An exact action representation of each finite observation history

The source instrument also has a **unitary history representation**, even
on one elementary carrier. Dissipation of the retained source state does not
exclude such a representation when its discarded record environment is kept.
The [exact construction](history_action.py) supplies a dimensionless Hermitian
action for it. This reuses finite instrument dilation and the composition
principle of [quantum networks](https://arxiv.org/abs/0904.4483); it is a concrete
attachment to the C/E endpoint maps, rather than a new general dilation theorem.
The [existing source compiler](../observer/source_discrimination/README.md)
separately realizes those committed instruments with charged coherent services.

Write each native seam Kraus operator as `K_alpha=W_alpha/sqrt(2)`, where
`W_alpha` is an integer matrix computed from the endpoint update and the C/E
record dictionary. Its label retains the fair coin and the hidden candidate
record as well as the visible current. Stacking these matrices gives
`W*W=2I`; therefore `V=W/sqrt(2)` is an isometry. On
`H_output direct_sum H_input`, set

```text
D = I-V V*,             U = [[D,V],[V*,0]],
P_minus = (I-U)/2,       A_0 = pi P_minus.
```

The exact identities `D*=D`, `D²=D`, `DV=0` and `V*V=I` give
`U*=U`, `U²=I` and `exp(-i A_0)=U`. An input in the second summand is mapped
to `V psi` in the first. Reading its current label and tracing the remaining
label gives exactly the native selective CP map, including its disturbance.
The completion has dimension 20 for C and 28 for E on the two active binary
ports. Spectator ports are extended by identity; no federation is required.
This direct sum is an encoding space, not a derived additional physical sector.

**Complete observation contract.** For an arbitrary density matrix on the
source and an untouched reference, projection onto each current label and
partial trace of its hidden label gives
`sum_(alpha:current=j) (K_alpha tensor I) rho (K_alpha* tensor I)`.
This is an identity of linear maps, not just of outcome probabilities.
Suppose complete unnormalized branches agree up to a given event. Apply
the same next intervention on each accessible record branch, retaining its
outcomes and any controller memory, and then the same local instrument.
Linearity gives equality of every resulting branch. Induction proves equality
of all final records and quantum outputs for every finite adaptive protocol
using those declared interventions. Fresh hidden labels are retained at each
event; resetting an old environment without accounting for its information
would change the contract. No protocol may gain access to those hidden labels
merely because a representation stores them.

The execution checks a two-repair protocol on three ports: the first current
chooses the next seam and a local phase operation. Local tensor addresses
`(0,1,2)` map to native carrier ports `(0,1,5)`. All three edges belong to
the source seam list, and all 48 local state/seam/coin branches agree with
their twelve-port restriction. All 64 input matrix units
and every complete current word agree exactly with a separate floor/ceiling
recursion. The shared runtime also reproduces the full terminal states with
an entangled untouched reference and a serialized checkpoint. Thus the
representation preserves sequential observations and feedback, not only one
unread channel. Arbitrary coherent interference between hidden record labels
is an additional experiment, outside this equality's accessible-record scope.

**What this action determines.** This construction needs no additional law to
reproduce the existing native event boundaries. It chooses a completion of
the source isometry. A finite implementation may choose such a completion;
uniqueness is not a requirement for implementing the instrument. However,
intermediate evolution in that completion is not determined by those event
observations. With `P_plus=(I+U)/2`, every
`A_m=pi P_minus+2 pi m P_plus`, for integer `m`, has the same endpoint `U`.
For an initial input-sector vector, the probability of returning to that
sector at interpolation parameter `theta` is

```text
[1+cos((pi-2 pi m) theta)]/2.
```

At `theta=1/3`, the positive actions `A_0` and `A_2` give respectively `3/4`
and `0`. This discrimination requires an intermediate read of the auxiliary
completion sector. Every allowed native-boundary observation, including
adaptive ones, remains identical. The example limits an inference to new
intermediate readouts; it does not reject a valid chosen implementation.

The dimensionless gate action, local dissipative Liouville form, path
surprisal and transfer logarithm are different objects. The construction
supplies no autonomous source Hamiltonian, lab seconds, `hbar` normalization,
Lorentzian spacetime action or continuum field limit. It also does not make
the necessary record storage, preparation and control services free. Such a
physical identification must preserve the declared observation and resource
contract. The existing nonpositive fractional transfer-step result constrains
that reduced transfer; it does not prohibit this enlarged history realization.

## A current action in an observer's own recorded clock

The [local-clock construction](current_clock_action.py) connects the existing
[current martingale](../event_circulation/README.md) directly to a clock made
from native records. The domain is the same 66-state binary total-two carrier
with uniform independent seam/coin attempts. The observer retains **all
carrier current records**, but its clock ticks only when occupation at port
`v=0` changes. This access is broader than a port seeing only incident records;
it needs no federation, second clock, independent replica or Poisson schedule.
It is a declared observation protocol, not a derivation of the source wiring
or choice of scheduler. The result concerns native record laws without
intervening controls. C and E share these laws; the result does not select
between their different coherence rules.

Let `Q0,Q1` be the global-movement transition kernels split according to whether
that movement leaves port `v` unchanged or ticks it. These are computed from
the endpoint branches, without fitted transition rates. Write `d(x)` for all
moving branches and `d_v(x)` for those ticking `v`. The global-movement
Palm law is `mu(x)=d(x)/600`. The local-tick Palm law and stopped kernel are

```text
nu_v(x) = d_v(x)/100,       F = (I-Q0)^(-1),       R = F Q1,
nu_v R = nu_v,             nu_v F = 6 mu.
```

Here Palm means the state law just after a stationary clock tick. The second
identity is an exact **stopped occupation measure**, stronger than replacing
random waiting times by their mean. In particular a local tick takes six
global movements on average under `nu_v`. This preparation is operational:
from `I66/66`, perform one native attempt and herald a tick at `v`. Its success
probability is `5/198` per attempt; failure `193/198` is part of the experiment.
The successful preparation has support 51, so it is not a faithful reference
on the full 66-state algebra. Starting with an arbitrary configuration or
waiting for the first tick from a uniform state does not supply this Palm law.

**Exact current moments.** For the source incidence matrix `B`, set
`P=I-B^T(BB^T)^+B`. Each global-movement cycle current `PJ` has zero conditional
mean. Let `T` be the number of global movements until the next local tick,
including its current. Direct enumeration
gives the uniform bound

```text
max_x Pr_x(T>3) = 2963/3200 = beta < 1,
Pr_x(T>3k) <= beta^k.
```

Thus `T` has moments of every order uniformly in the starting configuration.
Each stopped current block is bounded in norm by `T`; stopped martingale and
second-moment identities follow by truncation with this integrable bound.
The occupation identity `nu_v F=6 mu` then gives the exact block covariance
`P/5`. The local-tick chain preserves `nu_v`, so for its consecutive current
blocks `D_l`,

```text
E[D_l | complete past at preceding local tick] = 0,
E[D_l D_l^T] = P/5,       E[D_l D_k^T] = 0 for k != l.
```

For the native triangle `(0,1,5)`, the signed circulation weight has squared
norm three. Consequently `Var(sum_(l=1)^m triangle_current_l)=3m/5` exactly
under the declared Palm preparation. These finite covariance statements do
not depend on exchanging a limit with expectation. Gradient currents
telescope to bounded occupation differences at the stopped endpoints.

**Scaling action.** `R` is irreducible on its 51-state support. It has period
two: each tick flips the occupation at `v`. Ergodic averages of its conditional
block covariances nevertheless converge; aperiodic pointwise convergence is
unnecessary. The geometric waiting bound supplies uniformly bounded higher
moments and hence the martingale Lindeberg condition. The same imported
martingale functional limit theorem used by the current packet therefore gives

```text
n^(-1/2) sum_(l<=floor(n*u)) D_l -> W(u),
Cov(W(u),W(w)) = min(u,w) P/5.
```

On a fixed partition of rescaled **local tick count** `u`, the limiting
Gaussian increment density on the cycle subspace has action

```text
A[y] = (5/2) sum_l ||y(u_l)-y(u_(l-1))||^2/(u_l-u_(l-1)).
```

The scalar triangle action has coefficient `5/6`. The rank-19 cycle space
counts current coordinates, not spatial dimensions. This is a finite-partition
fluctuation action; no differentiable Brownian-path density or finite-count
Gaussian error estimate is asserted. Local occupation has exact multiplier
`-1` in this clock because the clock selects flips. This alternation is not
evidence of a propagating spatial wave and is not a modular phase.

**Finite records retain more information.** The exact resolvent also computes
the full first four stopped current moments. If `Qb^(r)` marks every branch
with its triangle score to power `r`, then for a terminal weight `b`,

```text
M_0 = F Q1 b,
M_r = F [Q1^(r) b + sum_(k<r) binom(r,k) Q0^(r-k) M_k].
```

The [retained receipt](current_clock_action_results.json) gives exact rational
values. A single block's triangle variance is `3/5`, but its fourth cumulant
is approximately `0.256840`, and the covariance of two consecutive squared
block currents is approximately `0.020723`. Thus finite blocks are neither
Gaussian nor independent; their zero signed cross-covariance does not erase
that dependence.

The same computed `M_2` also predicts nonstationary preparations: for any
declared initial configuration law `p`, the variance through `m` local ticks
is exactly `p (I+R+...+R^(m-1)) M_2`. Conditional zero mean removes all
cross terms. This retains the initial-state correction instead of applying
`3m/5` to a run that did not start in the local Palm law.

The shared quantum runtime executes the native herald and a **fixed two-attempt
window** for each C/E candidate, retaining all 98 words, the herald failure,
full terminal densities and serialized continuation. It continues through the
second attempt even when the first ticks. Conditional on herald success, the
probability of no further local tick in that window is exactly `2033/2250`.
Independent floor/ceiling branch counts check every terminal density. This
finite execution is not an execution of the infinite stopped quantum channel.
The classical stopped-current calculation does not replace that channel's
coherence or waiting-time effects.

A separate 49-word native experiment uses the first signed current to choose
the second seam. Its complete record law differs from the unchanged uniform
schedule by total variation `28/5625`; the exact joint generating-function
value `E[2^(g1+g2) 3^(tick1+tick2)]` is `153281/120000`. This is a controlled
negative example: stationary action coefficients cannot be carried through
arbitrary feedback without rechecking its law. Neither this source-derived
action nor a conversion of tick counts into units establishes laboratory
seconds, relativistic propagation, a Lorentzian field action or correspondence
with measured physical fields.

## Exact path weights and the optional continuous parameter

An optional bookkeeping construction gives each seam an independent
rate-one Poisson attempt process. Its parameter `s` counts expected
attempts **per seam**. The total count is `N~Poisson(Ms)`, and conditional
on `N=n` the ordered seam word is exactly the original uniform scheduler.
Consequently the unread channel is

```text
T(s) = sum_(n>=0) exp(-Ms) (Ms)^n/n! Phi^n
     = exp[-s H_native].
```

This is also a local Lindblad generator: writing each native Kraus map as
`K_(e,a)`, completeness gives
`-H_native(X)=sum_(e,a)[K_(e,a)XK_(e,a)^* - {K_(e,a)^*K_(e,a),X}/2]`.
No Hamiltonian term or fitted rate is required for this particular
scheduler extension. The Poisson count distribution itself is declared;
it is not derived as an observer's physical clock.

More explicitly, one ordered seam-and-coin path of length `n`, after
integrating its ordered attempt locations within `[0,s]`, has coefficient
`exp(-Ms)(s/2)^n/n!`. A matrix unit contributes only if its left and right
native paths give the same hidden Kraus label at every step. Visible
currents are then retained, while hidden alternatives are summed as CP
maps. This is a local path-weight action with hard record-compatibility
constraints. Complex input coherences and final quantum readouts also
contract with this kernel; the construction is not a classical
hidden-variable replacement for arbitrary quantum experiments.

For protocols depending on event order and records, every joint record law
at a deterministic native count remains available by conditioning this
auxiliary count construction on `N=n`. A protocol using externally assigned
Poisson locations needs those locations retained as additional inputs. Clock readings
formed from the native records can therefore be retained directly.
A physical duration, modular-flow identification, Gaussian measure and
Lorentzian action do not follow from this reparameterization.

## Complete local observation theorem

Let `A` be the fixed set of observed ports, `Delta` the maximum graph degree,
and `B_R(A)` its closed graph-radius-`R` ball. Compare the full process with
a simulation retaining this ball and its initial reduced state, including
any untouched reference. Delete seams crossing the ball boundary. In a
fixed-count discrete simulation, keep deleted attempts as identity holds
under the original scheduler, preserving its seam probabilities.

The observer may keep every native current/coin record of seams wholly
inside `A`, execute local instruments on `A`, and choose later such
instruments from those local records. It may also retain the final quantum
state on `A` jointly with all those records and an arbitrary reference.
Local memories storing these records are included. The additional intervention
protocol is finite, or locally finite and well posed on every bounded horizon;
an accumulation of infinitely many control operations is excluded. The controller may not
consult external records, act outside `A`, or observe the global attempt
tape as an extra signal. The native scheduler is state independent.

For the **unconditional complete output**, in trace distance, the error is
at most

```text
discrete n attempts:
  min{1, |A| binom(n,R+1) (Delta/M)^(R+1)},

optional per-edge Poisson parameter s:
  min{1, |A| (Delta*s)^(R+1)/(R+1)!}.
```

The binomial coefficient is zero when `n<R+1`. The same bounds control all
measurements of the retained output, and hence joint field and local-clock
readouts within the declared horizon. Selected rare outcomes can amplify
conditional errors; the bound does not silently become a postselected
bound. An unbounded stopping time requires an additional stopping-tail
estimate. Finite caps must retain their timeout outcomes.

**Proof.** Couple the full and truncated source schedules inside the ball.
For a fixed seam word, trace the dependency cone of `A` backwards. An edge
can enlarge that cone only when it meets the existing cone. Record memories
and their allowed feedback remain in `A`. If the cone never crosses the
retained boundary, the whole local instrument is supported within the ball.
The original initial state and its restriction then give identical outputs,
including arbitrary references. This is the actual endpoint instance of
[DependencyCone](../../Lean/ObserverPatchHolography/Locality/DependencyCone.lean)
and the nonselective finite tensor identity in
[SlotLocalitySurface](../../Lean/EventAlgebra/SlotLocalitySurface.lean).
The causal-access requirement follows the existing
[AdaptiveScheduler](../../Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean)
contract; no new Lean proof is claimed.

Crossing the boundary requires a temporally ordered walk of at least
`R+1` edges starting in `A`. There are at most `|A| Delta^(R+1)` possible
walk prefixes. In the uniform discrete schedule, a particular ordered
edge list occurs at specified `R+1` attempt positions with probability
`M^(-(R+1))`; union over their `binom(n,R+1)` choices. For independent
rate-one seam Poisson processes, the expected number of ordered occurrences
is `s^(R+1)/(R+1)!`, including repeated edges by the Poisson factorial
moment identity. A union bound gives the two displayed probabilities.
On every good schedule the complete outputs coincide; on a bad schedule
the trace distance of normalized outputs is at most one. Averaging proves
the result without any dimension factor or reference-size restriction.

The idea is related to general dissipative
[Lieb--Robinson locality](https://arxiv.org/abs/1003.3675).
Here the explicit native Poisson unraveling permits the elementary proof
above, including local records. Its bound does not assume that a
single-particle swap coupling also controls coherent histories.

## Finite executions and conditional federation limit

The runtime uses the actual source path `(0,1,82,72)` across two elementary
carriers, with its three seams uniformly selected for three attempts. This
is a declared three-seam experiment, not the full-federation scheduler.
All binary states on the four active ports are included; all other ports
are empty and are not touched by this declared schedule. The initial probe
is entangled with a reference qubit, and a remote port starts occupied.
Current records on the first seam control a local phase intervention.

For each candidate, all 27 schedules are executed both in full and with
the outermost seam replaced by identity. The runtime retains 125 complete
local record branches across these schedules, full terminal matrices and
a serialized checkpoint continuation. Exactly `7/27` of schedules allow
a boundary-crossing cone. All others agree to roundoff. An independent
rational test checks complete local output equality for **all 256 input
matrix units on every good schedule**, including the feedback operation.
The bad-schedule controls have a nonzero discrepancy, so truncation is not
assumed exact at arbitrary horizons.

A second executed control reads an external port and flips a local port
from that record. Two initial states agreeing on the radius-two ball then
produce perfectly distinguishable local outputs with zero native attempts.
That controller violates the local-access hypothesis. Local native gates
alone do not constrain information supplied through a global controller.

The receipt also evaluates the uniform bound for the existing 240, 960 and
3,840-port federations, whose maximum degree is six. For two observed ports,
`R=8` and `s=1/4`, the complete quantum-history error allowance is
`243/1146880 < 0.000212` at every size. The discrete attempt allowances and
actual retained ball sizes are reported separately. This is a simulation
error guarantee under the stated scheduler and local-access assumptions,
not an empirical error bar or a measured physical propagation speed.

For a **fixed bounded-degree infinite graph** with consistent local initial
states, use independent per-edge Poisson processes with the auxiliary parameter
`s`. A uniform choice among all edges of an infinite graph is unavailable.
The Poisson bounds make finite-ball quantum observation laws Cauchy for every
bounded `s` horizon. The two-ball comparison uses a common
inner ball and the sum of their truncation allowances. Positivity,
normalization and compatibility pass to the limit. This constructs a
continuing local quantum process and all specified local instrument
histories from federations of the same elementary factors. The limiting maps compose to give CP dynamics on the quasilocal algebra.
An individual observer subalgebra need not evolve autonomously: its reduced
description can retain the memory exhibited by the joint clock process.

This conditional thermodynamic construction does not require a new larger
elementary carrier. It does require the graph neighborhoods, source
instruments, initial restrictions and observer access maps to agree. The
finite screen identifications below supply these hypotheses on their
declared neighborhoods, but not on an entire infinite sequence. The direct
parent-map coarse-graining obstruction applies. The result supplies no
Lorentz symmetry, continuum dimension, or physical identification of `s`.
Those require a further observational correspondence and scaling argument.

## Compatible neighborhoods of the actual screen

The [screen computation](screen_compatibility.py) uses the unchanged pinned
`port_pair` federation builder through level five. A root carrier is the
central child of face zero at each subdivision, with cell address `4^L-1`.
It identifies every vertex and induced seam of its port-graph neighborhood,
preserving local port labels, elementary-carrier membership and native
instrument semantics. The retained [result](screen_results.json) contains
the complete source-address maps; agreement of hashes alone is not the test.

| Source levels | Global ports | Matched radius around all twelve root ports | Matched ports |
| --- | --- | --- | --- |
| 1 and 3 | 960 and 15,360 | 3 | 45 |
| 2 and 4 | 3,840 and 61,440 | 6 | 114 |
| 3 and 5 | 15,360 and 245,760 | 9 | 222 |

These are local observation identifications within the existing spherical
federations. They are not the coarse cell-parent map and do not make coarse
counts autonomous. A central subdivision reverses its local orientation;
adjacent levels do not preserve the same root port labels. The listed
comparisons use parity subsequences. At the next radius in each comparison,
the particular fixed-label identification fails, and its explicit differing
adjacency row is retained. That failure does not exclude a different port-frame
isomorphism or a controlled approximate comparison.

With matching initial density restrictions on the common ball, including any
untouched reference, the two complete root-observer output laws differ by at
most the sum of the two truncation bounds. At levels three and five, for
per-seam attempt intensity `s=1/4`, this gives

```text
trace distance <= 2 * 12 * (6/4)^10 / 10!
               = 2187/5734400 < 0.000382.
```

This controls the whole local quantum output, all allowed local records and
their feedback, under the local-access hypotheses above. A record from a seam
crossing the root's boundary enlarges the observation support and requires
recomputing the bound. Comparing equal global attempt counts would give
different per-seam intensities: the two screens have respectively 40,320 and
645,120 seams. The comparison therefore uses the stated independent per-seam
Poisson extension; it does not identify this parameter with observer time.

The executable also transports an actual four-seam C/E experiment between
source paths `(757,756,735,734)` and `(12277,12276,12243,12242)`. It includes an
entangled reference, a current-controlled local phase operation, all 81
complete current words and serialized continuation. Source seams are sorted
global addresses, so signed currents must be reversed where their orientation
differs from the common path. Independent rational tests compare every visible
outcome on all 256 four-port matrix units; hence the correspondence extends
to every input and untouched reference on this support. The runtime checks
the full branch density matrices. This finite tape holds every unselected
seam: it tests instrument transport, not a simulation of the full Poisson
screen process or a replacement for its analytical error bound.

### A nontrivial compatible state domain

Use the checked full binary tensor algebra, with no separate
total-two constraint on each elementary carrier. A declared density `0<rho<1`
gives the faithful product reference
`sigma_rho = tensor_ports diag(1-rho,rho)`. Its local restrictions agree exactly
under every retained address map. Every native unread seam map preserves it:
the reference is scalar within each total-number block, while each map is
unital and conserves number. The executable checks stationarity at `rho=1/3`
and branchwise number-charge preservation for both candidates. These
references have nonzero local fluctuations at every screen size. Their density
is supplied by preparation, not selected by repair.

On nested finite port sets, the embedding `A -> A tensor I` and restriction
`identity tensor reference expectation` are state preserving, completely
positive and compose exactly. This is the tensor-reference instance of the
paper's [reference multiresolution construction](../../paper/screen_microphysics_and_observer_synchronization.tex)
and the existing [CP restriction contract](../../Lean/QFT/CPRestrictionNet.lean).
The modular flow is generated by `log(rho/(1-rho))*N`. Every selective native
branch commutes with that flow because it preserves number charge; stationarity
is asserted only for the unread sum. The modular action is trivial on the
current POVMs and fixed-total observables. It does not furnish observable
clock ticks or equate modular parameter with Poisson intensity.

An arbitrary matching local preparation, including a root entangled with a
reference, can replace the product on the matched region. If the initial
restrictions differ by trace distance `delta`, add `delta` to the output bound.
The uniform **global total-two** preparations are an important different case:
their 222-port restrictions at levels three and five differ by exactly
`831558460593/30921625444352`, about `0.02689`. Their occupation law is
hypergeometric; for any fixed local region it approaches the vacuum as the
screen grows. Geometry compatibility alone cannot turn those dilute
preparations into a nontrivial finite-density limit.

### Exact local density-response predictions

The native first-moment law gives the common density uniformization
`Q=I-L/12`, where every port degree is at most six. Its root `12 x 12` blocks
agree exactly between levels three and five through `Q^22`; a root-port
response first differs at power 23. All integer numerators and the first
differing entry are retained. The radius-nine identification guarantees
agreement through power 18 by returning-path support; the larger exact
agreement range is an executed finite result.

For the Poisson screen scheduler, the actual density response is
`exp(-s L/2)=exp[6s(Q-I)]`. Thus the common coefficients give a tighter
bound on this particular readout: after preparing one particle at any root
port and leaving the other ports empty, its probability of being found at any
specified root port differs by less than `10^-17` at `s=1/4`. The receipt uses
the exact rational majorant
`(3/2)^23 / (23! * (1-(3/2)/24))` for the remaining Poisson tail. This is a
mean-density observation prediction under supplied preparation and final
occupation readout. It does not control the complete quantum history, assign
one `Q` step to one native global attempt, or compare with laboratory data.

### Comparing fields in the observer's own incident clock

The [incident-clock construction](INCIDENT_REFINEMENT.md) removes the auxiliary
Poisson parameter from the finite source-screen comparison. It counts ordered
causal seam words before the observer's kth incident attempt, including
zero-current attempts. Actual boundary-walk counts yield a complete local
quantum and retained-record output distance below `6.656e-9` after four ticks
and `4.296e-7` after eight ticks between levels three and five. The observation
support includes all fifteen endpoints of the root's 33 incident seams.

At Bernoulli density `1/3`, the same stopped source has a nonzero, certified
two-time root-field covariance after four ticks in
`[0.15332581621, 0.15332582288]` on either full screen. Directed integer
propagation charges both unresolved waits and arithmetic losses. Executed
full-screen histories retain native feedback, all clock marks, external-cap
timeouts and serialized continuation. Separate exact C/E instrument controls
check every four-port matrix unit and an entangled reference. These domains
are distinguished in the shared contract; a first-moment density reduction
does not replace a quantum history.

### What a limit requires

The nearest [tower theorem](../../Lean/Tower/ConsensusTower.lean) requires
compatible algebras, states and dynamics over its whole domain. Here the
finite tensor identifications, explicit compatible states, native instrument
transport and bounded local output comparison are supplied for the listed
source pairs. An all-level stabilization or approximation theorem is
required. A proof would need to control spherical subdivision geometry,
local port-routing decisions and state restrictions as the radius increases;
the pinned executable routing uses floating direction scores. Six finite
levels do not prove such an induction. These certificates therefore establish
finite source-screen compatibility with quantified observational error, not a
continuum limit, emergent physical dimension or relativistic spacetime.

```sh
python -m model.native_action --check model/native_action/results.json
python -m model.native_action.history_action --check model/native_action/history_action_results.json
python -m model.native_action.screen_compatibility --check model/native_action/screen_results.json
python -m model.native_action.incident_refinement --check model/native_action/incident_results.json
python -m pytest -q -W error model/native_action
```
