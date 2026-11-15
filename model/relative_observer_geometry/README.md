# Relative observer fields, records and clocks

This packet closes a specific observer-access gap: on a **known diagonal native
preparation**, absolute occupation fields can be read from retained signed
currents. They need not be simulator-only variables. Two regional observers
can reconstruct compatible parts of the same field history, continue from
their own checkpoints, and compare local event clocks on shared records.

The finite readouts and the simultaneous refinement of their incident clocks
are conditional on declared access and scheduling. They do not identify the
fields as physical matter or the local counters as physical time.

## What executes

The carrier is the existing binary total-two sector: 66 configurations of two
occupied ports among twelve. The native endpoint rule acts on an oriented
admitted seam `(a,b)`, commits its fair coin and signed current `j`, and changes
the loads by

\[
n_i(k)-n_i(k-1)=j_k(\mathbf1_{i=b_k}-\mathbf1_{i=a_k}).
\]

The rule, carrier, preparation service, record access and scheduler are premises. The packet first checks this equality on **all 3,960** native
configuration/seam/coin branches. Summing these one-step identities gives an
exact reconstruction for every finite prefix, without importing a field action:

\[
n_i(N)=n_i(0)+\sum_{k\le N}j_k(\mathbf1_{i=b_k}-\mathbf1_{i=a_k}).
\]

Each `LocalCounter` receives only its initial local loads and every incident
seam record. Boundary seams are essential, including seams whose other endpoint
is outside the region. Nonincident records are absent; incident zero-current
records remain. Common event identifiers authenticate correspondence and give a
source-order cut. They carry no duration, physical simultaneity or global clock.
Authenticity and completeness of the incident record service are premises:
omitting a boundary current can be locally undetectable.

The retained finite example supplies the seam tape
`(0,1), (0,1), (0,5), (0,5), (0,1)`. It enumerates all **2,112** histories from
66 equally weighted, **disclosed** initial configurations and 32 fair-coin words.
There is no successful-outcome selection. These are exact laws on a declared
tape, not the stationary uniform-seam experiment. Regional observers access
ports `0..6` and `5..11`; their overlap agrees on every history and their union
covers all twelve ports. At the same retained prefix their values therefore
glue uniquely. Local checkpoint restoration and continuation agree for every
history. No hidden load vector occurs in a local checkpoint.

Separately, the primitive CQ runtime executes the first three operations from
the disclosed preparation with ports 0 and 11 occupied. Each operation has six
coin/current outcomes. All **216** outcome words and full terminal density
blocks remain, including **208 zero-probability words**. Exact enumeration
checks the eight positive words, and runtime continuation checks the complete
serialized checkpoint state. The record-only `n_0` expectation equals the separately
audited terminal occupation expectation. Zero-probability impossible words are
kept with an arbitrary zero public-field value and are not declared physically
consistent trajectories. Floating runtime comparisons are diagnostics at
`2e-11`, not directed roundoff certificates; the diagonal enumeration and
continuity identities are exact rational/integer computations.

## Local fields and common state

For a known basis preparation, the current prefix determines each local field
value and its products. Averaging those record functions gives the same
occupation correlators as the diagonal source law. A randomly chosen initial
configuration is usable in this way only when its local preparation label is
disclosed. An unknown initial law instead permits conditional estimates; it
does not generally reveal the realized occupation. This packet does not
implement a Bayesian filter or extend the pathwise assertion to arbitrary
coherent preparations.

The existing effective-field construction can predict **current increments**
without this initial disclosure. Absolute occupation readback and increment
readback thus have different access requirements. No extra direct occupation
measurement is silently granted here.

At the classical record level, a local view is a deterministic map `v(h)` from
the retained common history. Its state is the pushforward of the common law.
Pullback gives the explicit finite GNS isometry

\[
(Jf)(h)=f(v(h)),\qquad
\sum_h p(h)|f(v(h))|^2=\sum_v p(v)|f(v)|^2.
\]

The complete finite enumeration records the local fibers and checks their
exact masses. Distinct indicator fibers are disjoint, so their off-diagonal
Gram entries vanish. This is a commutative record-algebra inclusion and a
compatible expectation, not a constructed quantum local net. In particular,
state marginalization does not supply the homomorphic algebra retractions in
[`FiniteCausalObserverNet`](../../Lean/QFT/FiniteCausalObserverNet.lean).
[`primitive_observer.algebra`](../primitive_observer/algebra.py) supplies the
actual CQ public-record expectation used for the runtime readback check.

Pathwise value gluing requires actual pointwise coverage and a shared retained
cut. It does **not** reconstruct a joint probability law from regional
marginals. The negative control uses the laws with equally likely occupied
pairs `(0,10)/(1,11)` and `(0,11)/(1,10)`: both regional marginal laws on `0..5`
and `6..11` agree, while `E[n_0 n_10]` is respectively `1/2` and `0`. This is
the kind of failure ruled out by an explicit descent receipt, not by isotony
alone; see [`ObserverNetDescent`](../../Lean/QFT/ObserverNetDescent.lean).

## A successful and a failed relative clock comparison

A local incident-event counter increments on each nonzero movement incident
to that observer. On the actual coin word `0,1,0,1,0`, shared seam `(0,1)`
movements occur at record identifiers 1, 2 and 5. The local `(port 0, port 1)`
tick readings are

| Shared event identifier | Port 0 ticks | Port 1 ticks |
|---|---:|---:|
| 1 | 1 | 1 |
| 2 | 2 | 2 |
| 5 | 5 | 3 |

The first two shared events fix the affine comparison `v=u`. The predesignated
third event has residual **−2** because port 0 also recorded the excursion
through port 5. Its probability is exactly **1/32 given the disclosed initial
preparation**. Local event counters need not be affine on every shared-event
domain; this does not preclude other valid comparisons or clock constructions.

A positive control reads the shared-seam movement count `u` and a second
declared display `v=2u+5`. The same three shared events give `(1,7)`, `(2,9)` and
`(3,11)` with zero third-event residual. The offset and scale are a supplied
display convention, not independently recovered physical clock calibration.

These computations use the exact two-anchor/third-event algebra isolated in
[`Time/ClockComparison`](../../Lean/Time/ClockComparison.lean). That Lean
interface additionally assumes timelike histories; no such geometric
assumption is instantiated here. The relation of these native local tick
records to observer-relative modular flow is open. A
preferred universal physical clock is not required.

## Inferring nonuniform scheduling from clock comparisons

[clock_reconstruction.py](clock_reconstruction.py) removes uniform seam
selection from the pair-clock geometry inference. Assume independent draws
from an unknown positive seam distribution on the connected twelve-port
graph. Every incident event ticks, including zero current. The comparison
service supplies the ordered pair tapes and terminal coincidences. Write
`a_i` for port i's incident probability and `p_ij` for the shared seam
probability, zero on a nonedge. Observable pair laws give

```text
mu_ij = E[number of j ticks through the next i tick] = a_j/a_i,
nu_ij = Pr[j also ticks at the terminal i event] = p_ij/a_i.
```

Each seam ticks exactly two ports, hence `sum a_i=2`. Set `r_0=1` and
`r_i=mu_0i`; then `a_i=2r_i/sum r` and `p_ij=a_i nu_ij` uniquely recover
the scheduler. Probability normalization supplies no elapsed-time unit.
The estimator checks all 132 oriented pair laws, opposite coincidences,
row completeness and connectivity. It receives no load preparation, seam
table, final state or coordinates. Positive support then recovers the
unweighted native rank-three Gram; unequal weights need not preserve that
band as the slow band of the actual dynamics.

For one pair with exposures `a,b` and shared weight `c`, its count generating
function is `(a-c+cz)/(a+b-c-(b-c)z)`. The inferred weights also determine the
marked density transfer `B(z)(I-A(z))^-1`, with the remote waiting tail
summed exactly. Under the declared stationary binary total-two preparation,
this predicts the complete joint law of the initial A occupation and the B
occupation after specified A ticks. Independent enumeration of all 66 states
and native endpoint branches reproduces four count/field comparisons for
uniform and biased schedules. Both recover the same structural Gram but
predict different joint fields; using uniform kinetics for the biased
clock records fails the control.

This identifies an IID scheduler from exact ensemble laws. A stationary
cyclic control has the same pair means and terminal coincidences, but its
zero-count probability is `4/5` rather than `4/9` for ports 0 and 1.
Temporal independence therefore is testable with further interval
records. Finite data also cannot distinguish an absent seam from an
arbitrarily rare one without an exposure bound. Local initial-bit disclosure and
incident-current readback are necessary for the field observation. A
physical clock or Lorentzian geometry is not inferred from this probability
normalization.

## Geometry connection and remaining boundary

The existing record-only native response calibration is executed again. Its
recovered unweighted graph gives exactly the **same thirty seam interfaces**
used by the field counters; its exact slow response Gram has rank three.
This reuses [`observer_reconstruction/native_geometry.py`](../observer_reconstruction/native_geometry.py),
including its supplied uniform-seam calibration and twelve preparations with
both loads at one port. The separate [single-carrier readout](intrinsic_geometry.py)
removes uniform seam exposure and fair calibration coins as requirements for
recovering this structural metric, as explained below. Those calibrations live in the larger total-two sector,
not the binary history sector. They connect two uses of the same endpoint rule
and labeled interfaces; they are not measurements on these same five-step
histories. No target coordinate array enters the reconstruction. The selected
slow-band unit Gram matrix `X=4P` supplies a finite Euclidean spatial interpretation up
to orthogonal isometry. Its 66 unordered port pairs fall into exact shells:

| Graph shell | Pair count | Squared chord length in the recovered unit Gram |
|---|---:|---:|
| Native seam | 30 | `2 - 2 sqrt(5)/5` |
| Two seams | 30 | `2 + 2 sqrt(5)/5` |
| Antipodal | 6 | `4` |

The executable checks every Gram pair against the graph-distance class and
confirms that all thirty native exchange seams are exactly the shortest
metric shell. The rank-three projector is symmetric and idempotent, so `X` is
a positive semidefinite rank-three Gram with unit diagonal and zero row sums.
Port relabeling preserves the exact shell census; a symmetric change to one
edge entry is rejected. Together with
[`native_effective_fields`](../native_effective_fields/README.md), this means
the source-derived hard-core exchange terms are local on this reconstructed
finite metric. The logarithm of the transfer can contain longer-range
configuration couplings.

This is an exact spatial interpretation of one finite twelve-port carrier,
not a measurement of the dimension of physical space. It supplies no
physical length unit, global spacetime, Lorentz cone, signal-speed law or
continuum limit. The time readouts below are record parameters, and the
metric reconstruction uses a separately disclosed calibration preparation.

### What one carrier establishes without a federation

One carrier suffices for the finite record algebra, local modular flow,
current/occupation identities, the exact response Gram and its spatial shells.
Two regional observers can also be restrictions of that **same carrier**;
their record agreement does not require a network of carrier copies. None of
these local results needs the spherical federation construction below.

The independently replayable [single-carrier result](intrinsic_results.json)
also reduces the calibration premises. Let `w(p,e)>0` be the recorded exposure
of seam `e=(a,b)` in preparation `2e_p`, and let `j` be its current. Every
endpoint difference is even, so the native rule gives
`j=1[p=a]-1[p=b]` for either coin. Dividing the joint record law by `w(p,e)`
gives the conditional current, and hence

\[
L_{\cdot p}=\sum_{e=(a,b)}\mathbb E[j\mid p,e](e_a-e_b).
\]

This reconstructs the same unweighted Laplacian and exact rank-three Gram
even when each preparation has a different seam schedule and its coins are
biased or deterministic. The executable runs 540 native branches with such
protocols and obtains exactly the original Gram. It uses no target coordinates
or source incidence table in the estimator.

The deterministic calibration current permits a further simplification:
**one native shot per preparation and seam**, 360 in total, recovers the same
Gram without probabilities in its input records. Conditional current variance
is exactly zero for these known even preparations, so empirical exposure
weights introduce no sampling error in the declared native model. Tests also
run every setting with coin zero alone and coin one alone, and with unequal
repetition counts. Known preparations, authenticated records and complete seam
coverage are required. A missing seam does not inherit the carrier
certificate; contradictory currents are rejected. This is a noiseless native
model result. Laboratory preparation and detector errors need a separate
error model; no experimental accuracy claim follows from the 360-shot result.
The seam labels disclose which endpoints interact. This is geometry recovered
from that accessible interface and its currents; it does not derive the
carrier architecture or its observation services from self-reference.

This premise reduction concerns the **structural readout**, not the dynamics
under every scheduler. Doubling the weight of one seam gives the actual mean
response `T'=I-(L+bbᵀ)/62`, with `b=e_a-e_b`. The correction `bbᵀ` is positive
and annihilates constants, so it cannot lower the first nonconstant eigenvalue
`5-sqrt(5)`. Exact elimination over `Q(sqrt(5))` shows that its multiplicity
falls from three to **two**. The two remaining modes are therefore slowest;
the normalized long-attempt response has rank two. This is checked for each
of the thirty possible doubled seams. Thus arbitrary positive exposure does
not establish a three-dimensional slow response. A native slow-band claim
needs the spectral certificate for its actual scheduler. Exposure removal
recovers the original structural rank-three Gram from those records.

The meaning of “spatial readout” also matters. All modes of the same recovered
response give the exact inverse-response Gram

\[
G=(L+\mathbf1\mathbf1^T/12)^{-1}-\mathbf1\mathbf1^T/12,
\qquad R_{ab}=G_{aa}+G_{bb}-2G_{ab}.
\]

The executable verifies `LG=I-11ᵀ/12` exactly. Since `L` is positive on the
eleven-dimensional constant-free subspace, `G` is a positive Gram of rank
**eleven**. The squared Euclidean distances `R` have shells `11/30`, `7/15`
and `1/2`, with the same pair counts 30, 30 and 6. Equivalently, solving
`Lv=e_a-e_b` with zero mean gives `v_a-v_b=R_ab`, checked on every pair.
This is an alternative response-derived metric, not an executed ruler
apparatus. It does not contradict the specified slow-band theorem: that
theorem selects its own spectral readout. It shows why the word “dimension”
must name the readout, rather than infer a physical bulk dimension from one
Gram rank.

A continuum claim needs a scale family or another limit construction with
compatible observable laws, effective locality and controlled errors. It
does **not** require this particular federation architecture, or even several
observers at each regulator. A single observer can access an enlarging family
of regions. The current federation is a conditional way to study extension,
joint records and boundary effects; it adds supplied inter-carrier routing
and does not by itself close the physical-space obligation. The smallest
decisive work can therefore remain on one carrier until an extension is
needed for a stated observable or limiting claim.

### Geometry and propagation on the binary quantum sector

The [binary spatial dynamics](binary_dynamics.py) construction removes the
need for doubled-load calibration when connecting geometry to the existing
quantum clocks. Its twelve known preparations are occupied pairs `(0,i)`,
for `i=1..11`, and `(1,2)`. Their occupation columns form an invertible
12-by-12 matrix `S`, although every preparation lies in the same 66-state
binary sector. The estimator receives only those preparation labels, seam
endpoints and **current records**; the placement coin stays hidden.

For each preparation/seam setting, equal endpoint occupations give current
zero. Unequal endpoints give current zero or the signed occupation difference,
each with probability one half under the declared native rule. Freshly
preparing the same known state and repeating until both current values occur
therefore identifies that setting's mean exactly. This uses the supplied fair
coin law; observed frequencies do not establish fairness or estimate the mean.
Addressed seams, preparation/reset and faithful current access are explicit
services.

There are 250 deterministic settings and 110 two-outcome settings. The latter
have a cap of 24 native attempts each. If every fresh coin is fair conditional
on the past, the probability of any timeout is at most
`110 * 2^(1-24) = 55/4194304`; at most 2,890 attempts execute. A failed setting
is an explicit timeout and supplies no calibrated geometry. The retained
seeded execution uses **578 attempts**, succeeds, and retains its entire
current transcript. The seed reproduces software branches; it does not certify
independent laboratory entropy. The estimator checks current support and
native consistency, not the authenticity of preparation services or reported
cap metadata. The execution driver enforces its cap.

The current means give the mean output columns `Y`; exact inversion recovers
`T=Y S^-1=I-L/60` and the same rank-three structural Gram as the integer-load
calibration. No exact probability table, target coordinates, doubled loads or
new quantum extension enter this binary protocol. The family of preparations
is a collection of reset trials, not one uninterrupted record history.

The operator connection holds on that same quantum carrier. Write
`N(f)=sum_v f_v N_v` for occupation observables. Each native C or E Kraus
branch is a partial injection of binary configurations. Its pullback of a
diagonal occupation has no off-diagonal entries, and complete branch summation
gives the exact Heisenberg identity

\[
\Phi^*N(f)=N(Tf).
\]

Thus the inferred first-moment propagation holds for **every initial binary
density**, including coherent preparations. This operator identity does not
grant pathwise hidden occupations for a coherent unknown initial state.
Known basis preparations and complete incident currents retain the existing
pathwise occupation readback. The C current instrument supplies the calibration
records; C and E share this unconditioned first-moment law while their other
quantum readouts can differ.

On the stationary uniform binary state, the centered occupation GNS inner
product and all attempt-lag covariances follow from the same source:

\[
C_0=\frac5{33}Q,\qquad
C_k=\frac5{33}T^kQ,\qquad Q=I-\mathbf1\mathbf1^T/12.
\]

The executable verifies the formula by propagating the complete 66-state
native kernel at lags 0, 1, 2, 4 and 8. The Heisenberg identity proves it for
every nonnegative integer lag. The constant-free occupation subspace of the
observer's GNS space has dimension eleven; its slow spectral subspace has
dimension three and carries the same spatial Gram. The stationary state is
used to compute this covariance reference, not required by the calibration.
These are exact source covariances, not a finite-data covariance estimator.

The recovered response also has an exact graph propagation front:
`(T^k)_ab=0` for `k<distance(a,b)` and is positive at the first allowed hop.
Every port pair is checked. Individual columns are response-kernel probes,
not fixed-total preparation interventions. An admissible zero-total contrast uses
the binary preparations `{0,3}` and `{1,3}`. Their mean occupations at port 2
differ by `0, 0, 1/1800, 41/27000` after zero through three native attempts;
independent full-state propagation verifies those numbers. The parameter is
uniformly scheduled native attempts, and distance is native seam hops.
Neither quantity is identified with modular phase, local movement count,
laboratory seconds or physical length by this result.

The same calculation states the remaining dynamical boundary precisely.
All eleven nonconstant mean modes have positive multipliers below one.
For example, the stationary GNS squared norm of `N_0-N_1` decreases from
`10/33` to `1459/5940` in one attempt. A state-preserving automorphism would
preserve that norm, so it cannot equal this record-forgotten native step on
these observables. A positive scalar clock rescaling likewise cannot turn
their decaying means into undamped wave phases. This concerns the direct
occupation/current readout. It does not rule out a Euclidean transfer
reconstruction, a larger reversible description with retained records, or
wave-like observations in a different justified sector.

### Temporal Hilbert space reconstructed from field records

The [temporal-record construction](temporal_records.py) recovers a finite
effective transfer from observer-accessible occupation correlations. Its
recovery function accepts only rational correlation values; it does not consult
the source graph, a hidden state vector, or a supplied Hamiltonian. The source
calculation separately establishes those exact inputs and checks predictions
beyond the reconstruction window.

The elementary protocol uses the uniform binary total-two preparation on the
same 66 configurations. Ports zero and one form the observer region. Each of
their nine incident seams counts as a tick, including a zero-current attempt;
remote repairs continue between ticks. Both initial bits are disclosed, and
their later values are reconstructed from signed incident currents. This is
the existing [incident-clock field protocol](../observer_time_fields/README.md),
without feedback. Preparation, independent uniform scheduling and fair native
coins are supplied services. No physical duration per repair is assigned.

Let `f=n_0-1/6`, and write `c_k=E[f(0)f(k)]`. The established positive-transfer
dictionary gives `c_k=<f,T^k f>`, with `0<T<=I`. Therefore the matrices

```text
G_ij = c_(i+j),       J_ij = c_(i+j+1)
```

are the inner-product and transfer matrices on the span of
`f,Tf,...,T^(r-1)f`. This is a reflection-time Gram construction from recorded
correlations. It is not the ordinary covariance matrix of a chronological
sample sequence, which would use time differences. When `G` is positive
definite, `R=G^(-1)J` is self-adjoint in its measured inner product. The exact
checks `J>0` and `G-J>0` make `R` a strict positive contraction.

Correlations at lags zero through fourteen determine a **seven-dimensional
cyclic field space**. Its eighth-vector Schur remainder is exactly zero. If
`a=G^(-1)(c_7,...,c_13)^T`, that equality is

```text
||T^7 f - sum_(j=0)^6 a_j T^j f||^2 = c_14 - a^T G a = 0.
```

Under the established positive-transfer representation, this proves
closure and the recurrence for every later lag. It is stronger than fitting
a finite window. The reconstructed transfer has characteristic polynomial

```text
(753 x^3 - 1681 x^2 + 1241 x - 303)
  * (18147 x^4 - 54265 x^3 + 60655 x^2 - 30039 x + 5562)
  / 13664691.
```

Exact real-root isolation and a rational logarithm remainder bound certify all
seven eigenvalues and `-log(lambda)` values in the
[receipt](temporal_results.json). The transfer eigenvalues range from about
`0.633967902` to `0.862715126`; the corresponding effective logarithm values
range from `0.147670740` to `0.455756954` per incident tick. They are
dimensionless model predictions, not measured energies in joules or a native
coherent evolution. Independent propagation of the complete 66-state source
agrees at every reconstructed lag and at held-out lags fifteen through twenty.

The two readable field combinations `n_0+n_1-1/3` and `n_0-n_1` split this
space into four symmetric and three antisymmetric modes. An independent
observability calculation on the twelve-port mean response finds the same
seven modes. Four of the eleven centered occupation directions are invisible
to these linear two-field correlations at every lag. This does not make them
invisible to every possible record protocol. The three antisymmetric temporal
modes have distinct eigenvalues; they are not the carrier's three-dimensional
spatial slow band. This result reconstructs the cyclic part seen by a declared
detector, rather than asserting that it is the full observer algebra or a
spacetime coordinate system.

Exact laws and finite sampled records have different guarantees. If every
used correlation has error at most `epsilon`, each `r`-dimensional Hankel
matrix has operator error at most `r epsilon`. With
`g=1/tr(G^(-1))<=lambda_min(G)`, the generalized eigenvalues have error at
most `2r epsilon/(g-r epsilon)` when `r epsilon<g`. This controls the finite
pencil; a noisy near-zero Schur remainder does not prove exact rank or an
infinite recurrence. The raw seven-mode Gram has `g` about `5.46e-23`, so
recovering all seven lines from ordinary sample precision is poorly
conditioned. The symmetric and antisymmetric Grams have the better certified
floors `1.09e-12` and `1.03e-6`, respectively, without eliminating their own
measurement requirements.

A separate execution uses 2,048 independent reset histories for each of lags
one through four. Later fields use only the actual current counters. All four
joint-one frequencies fall within `1/16` of their exact predictions; the
conditional simultaneous Hoeffding failure bound is below `9.01e-7` under
the declared independent source services. The fixed seed makes the execution
reproducible and does not certify IID randomness. These are source-model
finite-lag checks, not experimental validation or seven-line spectroscopy.

### Temporal representation on the actual source screens

The same Gram construction also applies to the distinct finite-density
[levels-three/five protocol](../observer_time_fields/REFINEMENT.md): Bernoulli
density `1/3`, the whole twelve-port root region, and 33 incident clock seams.
Its existing directed source evaluator and complete-output locality bounds
enclose the root covariance at lags zero through three on either full screen.
All outward export rounding is included before constructing the pencil.

These intervals certify `G>0`, `J>0` and `G-J>0` on the two-vector space
`span{f,Tf}`, with `f=n_0-1/3`. The lower Gram bound exceeds `0.0001655` and
the ordered eigenvalues of the compressed transfer are enclosed near
`0.88702` and `0.97079`, with error below `0.00002`. Thus a nontrivial local
temporal Hilbert subspace and its positive effective transfer survive the
actual finite-screen comparison, including the arithmetic and spatial
budgets. These two values are **compression eigenvalues**, not claims of two
isolated spectral lines or closure under later evolution. They do not equal
the elementary source's seven-mode spectrum. These finite temporal pencils
do not identify a physical spacetime metric or a laboratory duration.

### Incident clocks on the simultaneous density continuum scale

The [actual density snapshot limit](../observer_time_fields/REFINEMENT.md#actual-native-density-observations-at-several-continuum-times)
uses global native attempt count at scale `a_l=1260*16^l`. Its time parameter
also has a local incident-event realization. Fix a nonempty event set before
drawing the independent uniform seam schedule. A port or region ticks once
when a seam meets it, **including zero current**. If its event set contains
`d_l` of the `M_l=630*4^l` source seams, then

```text
p_l = d_l/M_l,                a_l p_l = 2 d_l 4^l,
N_l(n) = number of its ticks through attempt n,
C_l(t) = N_l(floor(a_l t))/(a_l p_l),
T_l(k) = first attempt with k ticks,   T_l(0)=0,
U_l(t) = T_l(floor(a_l p_l t))/a_l.
```

These are actual incident records and their inverse count, with a declared
dimensionless normalization. A primitive port has `d_l=5` or `6`; a whole
carrier has `d_l=33`. A complete base-face patch has
`d_l=(63*4^l+3*2^l)/2`. The executable checks these counts on the actual source
and ancestry at levels zero through two. The formula for every level follows
from thirty internal seams per carrier, three intercarrier seams per cell,
and `3*2^l` seams crossing the ancestral patch boundary.

For any `T,delta>0`, put `n_*=ceil(a_l(T+delta))`. The centered count
`Z_n=N_l(n)-np_l` is a martingale with variance `np_l(1-p_l)`.
Applying the maximal inequality to its squared absolute value gives

```text
Pr[max_(n<=n_*) |Z_n| > a_l p_l delta]
    <= n_* (1-p_l)/(a_l^2 p_l delta^2).
```

On the complementary event, uniformly for `0<=t<=T`,
`|C_l(t)-t|<=delta+1/a_l`. Also `N_l(n_*)>=a_l p_l T`, so every requested
inverse count is reached within the retained horizon. At its first occurrence,
`N_l(T_l(k))=k`; therefore `|U_l(t)-t|<=delta+1/(a_l p_l)`.
The receipt exports these rational finite bounds, capped at probability one.
For example, choosing `delta_l=2^(-floor(l/2))` proves uniform convergence
in probability of both clocks to the identity for **every nonempty** such
event set: `a_l p_l>=2*4^l` and, at `T=1`, the failure bound is at most
`2/2^l`. The displayed bound converges to zero on every fixed finite horizon.

Several observers use the same seam schedule. Their clocks are correlated;
their pair covariance per attempt is `p_(A intersect B)-p_A p_B`.
A union bound gives simultaneous clock convergence without an independence
assumption. On the simultaneous good event, B's normalized display at A's
inverse count differs from `t` by at most `delta_A+delta_B+1/(a_l p_A)`.
The shared count horizon must cover both requested inverse clocks. Complete
enumeration of the elementary native source also finds covariance
`-25/2376` between port zero's one-attempt tick and `n_0(0)n_0(1)`.
Thus clock–field dependence cannot be discarded at finite level.

The time-change step is pathwise. For a density path `X_l` in a normed state
space and its modulus `omega_(X_l,H)`, whenever `U_l([0,T])` lies in `[0,H]`,

```text
sup_(t<=T) ||X_l(U_l(t))-X_l(t)||
    <= omega_(X_l,H)(sup_(t<=T)|U_l(t)-t|).
```

The linked density construction proves a continuous path limit in
`C([0,T],H^-3(S2))` for linear interpolation of actual density snapshots and
in `C([0,T],H^-2(S2))` after their finite `G_l^(1/2)` readout filter. Its step
versions are uniformly close to these interpolations. The displayed modulus
bound therefore transfers the same limiting laws to observations made at any
fixed finite collection of incident clocks; the stopped step versions
converge in `D` to the same continuous laws. The joint convergence uses a
deterministic inverse-clock limit and adds no independence assumption between
clock marks and density records. Finite-dimensional density convergence
alone would not justify this random-time substitution.

There is also a joint fluctuation prediction. For finitely many clocks with
fixed bounded incidence counts `d_i` and pairwise shared-seam counts `c_ij`,
define

```text
W_(l,i)(t) = [N_(l,i)(floor(a_l t))-floor(a_l t)*d_i/M_l]
               / sqrt(a_l/M_l).
```

Its exact covariance is
`floor(a_l min(t,u))/a_l * (c_ij-d_i*d_j/M_l)`. Independent schedule
increments give a joint Gaussian finite-dimensional limit with covariance
`min(t,u)c_ij`. To check the triangular-array CLT directly, a fixed linear
combination of centered one-step counts has magnitude at most
`sum|theta_i| sqrt(M_l/a_l)`, which tends to zero, and its total variance
converges to the displayed quadratic form. Taylor's characteristic-function
remainder is bounded by this maximum times the total variance divided by six;
the quadratic logarithm remainder also vanishes. Disjoint time blocks are
independent. This proves the finite-dimensional Brownian law, including
degenerate combinations, without imposing independent port clocks.

For adjacent degree-five ports the limiting normalized correlation is `1/5`;
for the two degree-six endpoints of an intercarrier seam it is `1/6`;
nonadjacent port clocks have limiting correlation zero. The executable
checks actual witnesses and the finite common-scheduler correction at three
source levels. Each pattern exists at every level: nine of each carrier's
twelve ports have degree five, each port has five internal neighbors, and
every external seam joins two degree-six ports.

Under the same stationary preparation, the density proof conditions on the
entire seam and auxiliary-swap history. Its conditional characteristic
function converges in probability, and hence in `L1`, to the deterministic
Gaussian value. Multiplying by any bounded characteristic function of the
seam-only clock vector therefore factors the limiting joint characteristic
function. Thus these centered clock fluctuations and the density snapshots
are **asymptotically independent**, despite the finite two-time-field
dependence above. This statement concerns finite vectors; it does not extend
to clocks selected by accepted repairs, currents or other state-dependent
events. The single linear endpoint field has zero covariance with
the clock mark at finite level; its zero covariance alone would not prove
the stronger joint result.

The port can count its own ticks; this does not grant it the complete spatial
record access needed for an arbitrary field smearing. Initial local field
disclosure, complete incident currents, any admitted regional aggregation,
the uniform independent scheduler and the source refinement are explicit.
The normalization connects native incident clocks to the density continuum
parameter. It does not identify that parameter with modular flow, Lorentzian
proper time, or atomic-clock seconds.

### Predictable port access on the density continuum scale

The fixed-port requirement can be removed for the path correspondence.
Before attempt `j`, choose a port `v_j` using any earlier available records,
including density records. Its degree `d_j` is five or six and is disclosed.
The fresh seam is independent and uniform on the same `M_l` seams.
This is passive access: choosing the port changes neither the preparation
nor the repair law. The path and any access needed to change ports are
supplied; autonomous motion or a physical worldline is not constructed.

Let `I_j` mark every seam incident on `v_j`, including zero current, and put

```text
b_l = a_l/M_l = 2*4^l,
R_n = sum_(j<=n) (30/d_j) I_j,
C_l(t) = R_floor(a_l*t)/(30 b_l),
T_l(k) = inf{n:R_n>=k},   T_l(0)=0,
U_l(t) = T_l(floor(30 b_l t))/a_l.
```

The clock uses integer credits, six or five per incident event. Conditional
on the entire past, its next credit has mean `30/M_l` and variance
`900[1/(M_l d_j)-1/M_l^2]`. Thus `R_n-30n/M_l` is a martingale even when
the access decision depends on earlier field observations. For
`n_*=ceil(a_l(T+delta))`, the squared-martingale maximal inequality gives

```text
Pr[max_(n<=n_*) |R_n-30n/M_l| > 30 b_l delta]
    <= n_* (M_l-5)/(5 a_l^2 delta^2).                      (PC1)
```

On its complementary event, `|C_l(t)-t|<=delta+1/a_l` uniformly on `[0,T]`.
Also `R_(n_*)>=30 b_l T`, so every inverse threshold is reached. At the
first crossing of a positive integer `k`, `k<=R_(T_l(k))<=k+5`.
Combining this overshoot with the floor error gives
`|U_l(t)-t|<=delta+1/(5 b_l)`, including the initial zero threshold.
Taking `delta_l=2^(-floor(l/2))` proves uniform convergence in probability
on every fixed finite horizon. A union bound treats finitely many such
access protocols. The preceding pathwise modulus inequality then transfers
the same `H^-3` density and `H^-2` Green-field limits to their inverse clocks;
the stopped step paths converge in `D` to the same continuous laws.
No clock–field independence is required. This argument does **not** transfer
the fixed-incidence Gaussian clock-jitter or independence prediction to
field-adaptive access.

The finite control executes all 1,260 first native branches on the actual
240-port source. A nonzero current at port zero moves the next readout to
its adjacent degree-six port one. Five histories switch access; all 1,587,600
two-step histories are retained in the conditional clock enumeration.
The conditional credit mean stays `1/21`, while the bare tick rate changes
from `1/126` to `1/105`. Choosing the endpoint **after** seeing the current
seam instead gives credit mean `1787/315`; discarding zero currents gives
`1/42` for this preparation. Both violate the required compensator.
Exact finite prefix and inverse-crossing tests check (PC1) and its
deterministic inversion, including threshold overshoot.

This establishes a dimensionless clock for adaptive observations of the
same density continuum. Same-degree bare clocks share their mean rate;
the result contains no derived velocity dependence, local duration for
changing access, or atomic-frequency relation. Passive modular-flow
identification is separately restricted by the
[all-input native modular obstruction](../relative_modular_dynamics/README.md#passive-repair-cannot-act-as-a-modular-automorphism).

### Active coherent waves in a native incident clock

The [admitted coherent wave](../observer_time_fields/REFINEMENT.md#two-continuations-of-the-same-native-density-spectrum)
can also use a primitive incident counter. The timing carrier and protected
wave register are distinct apparatus regions. Each event at a fixed port
of degree `d=5` or `6`, **including zero current**, triggers one fixed walk
step `W_l`. The timing source pauses during atomic coherent services, as in
the [existing event-driven modular probe](../transfer_modular_probe/README.md#a-native-event-drives-the-modular-step).
The wave register receives no passive timing-source repair between steps.
Consequently its continuing operation is exactly `W_l^N` after `N` ticks;
interleaved passive repairs on that register would define a different law.
The coherent preparation, reflection, protected storage, common amplitude
and Ramsey readout are the stated apparatus contract. Service events
are charged separately from timing-carrier incidents.

Set `Q_l=M_l 4^l=630*16^l`, `D_l=4^l L_l` and
`beta=2d/sqrt(630)`. The native clock gives
`N_l(floor(a_l tau))/sqrt(Q_l)=beta C_l(tau)`.
On its preceding finite good event, let `eta=delta+1/a_l`. For the proved
spectral band `0<=D_l<=Lambda<4Q_l`, the coherent dispersion estimate gives

```text
e_l = Lambda/[6(4Q_l-Lambda)],
sup_(tau<=T) ||cos(N_l(floor(a_l tau)) arccos K_l)
                - cos(beta tau sqrt(D_l))||_band
    <= beta sqrt(Lambda) [eta+(T+eta)e_l].                 (WC1)
```

Here `K_l=I-L_l/(2M_l)` is the density transfer and the cosine is the
existing walk's compressed amplitude. Indeed its random wave parameter
differs from `beta tau` by at most `beta eta`; its frequency error costs at
most `beta(T+eta)sqrt(Lambda)e_l`. The cosine Lipschitz bound proves (WC1)
for each retained clock history. Independence of clocks and fields is
unnecessary. Along the same spatial subsequence, the proved bounded-band wave convergence and the clock LLN therefore give
`cos(beta tau sqrt(A))` in probability, uniformly on compact parameter
intervals. This closes the independent wave-counter requirement within
this active protocol; it does not select the continuation from passive repair.

On the good event at most `N_max=floor(a_l p_l(T+eta))` steps are triggered.
If each implemented step is unitary and differs from `W_l` in operator norm
by at most `epsilon_l`, telescoping adds at most `N_max epsilon_l` to the
amplitude bound. For raw arrays with norm at most `1+epsilon_l`, the bound
is `(1+epsilon_l)^N_max-1`, bounded by
`N_max epsilon_l/(1-N_max epsilon_l)` when the denominator is positive.
The error estimate vanishes when `sqrt(Q_l) epsilon_l -> 0`.
A fixed gate tolerance does not prove this growing limit. No source-pulse
compiler for the growing walk is executed here; preparation/readout errors
need their own allowances.

An exact elementary control retains all 237,600 two-attempt timing-carrier
histories and both Ramsey outputs for the existing 3,960-arc walk. The
protected input is `(6n_0-1)/sqrt(330)`. Its plus probabilities after zero,
one and two incident-triggered steps are `1`, `21/22` and `91/110`; the
complete two-attempt probability is `1297/1320`. Triggering only on nonzero
current gives `7901/7920`; repreparing before each triggered step gives `71/72`;
the opposite continuation reflection gives `391/396`. Thus the same native
clock history distinguishes these operational choices. These are finite
apparatus predictions, with no atomic-clock or observed Lorentzian-field
identification.

### Lorentzian density of the product wave

The [native-record product construction](../physical_space_time/OCCUPATION_FIELD.md#independent-native-records-give-a-radial-product-field)
has a sharper geometric consequence once its **comparison wave dynamics**
are specified. Write its spatial generator as
`B=-d_r^2+(A⊕0)` on the chosen interval times the source screen. The restored
angular constant mode is necessary here. The existing proof supplies a regular,
strongly local form with domain `H1` on the screen, equivalent to its ordinary
gradient energy. These hypotheses determine a measurable symmetric uniformly
elliptic coefficient tensor in each chart. Smoothness is not a hypothesis
needed for this measurable representation.

For completeness, let `Gamma` denote the form's energy measure, normalized
so `E(u)=integral dGamma(u)`. On a chart choose a smooth cutoff `v` and a
coordinate-linear function `w=xi.x`. Strong locality and the chain rule give

```text
E(v sin(nw)/n)+E(v cos(nw)/n)
    = integral v^2 dGamma(w)+E(v)/n^2.
```

Apply the proved upper and lower gradient-energy bounds to this identity
and let `n` tend to infinity. It follows that `Gamma(w)` is absolutely
continuous, with a density bounded above and below by positive multiples of
`|xi|^2`. Polarization gives the chart coefficients
`a_screen^ij=dGamma(x_i,x_j)/dmu`. The chain rule reconstructs the whole
form from these coefficients. This is the local application of
[Cipriani, Theorem 7.2 and Corollary 7.3](https://msp.org/pjm/2006/223-2/pjm-v223-n2-p03-p.pdf),
applied to the part form with `H1_0` domain in a bounded chart.
The model supplies the required form/domain hypotheses. Constants have zero
energy, and the radial product adds the coefficient `a_rr=1`. Chart changes
preserve the resulting measurable tensor.

Thus, in local product coordinates with measure `mu(x) d^3x`,
`E_B(phi)=integral mu a^ij partial_i phi partial_j phi d^3x`.
The chosen wave action has the coefficient density

```text
S(phi) = -1/2 integral C^{alpha beta} partial_alpha phi partial_beta phi d^4x,
C = mu diag(-1,a),
g = sqrt(-det C) C^-1 = mu sqrt(det a) diag(-1,a^-1),
sqrt(-det g) g^-1 = C.
```

The positive spatial tensor forces signature `(-,+,+,+)` almost everywhere.
The weak equation `partial_alpha(C^{alpha beta} partial_beta phi)=0`
is exactly the product wave equation. This reconstructs its full metric
density, including the conformal normalization, from its action. The
density-to-metric identity is established mathematics; see
[Visser, equations (23)–(26)](https://arxiv.org/pdf/gr-qc/9712010).
Its application here does not import fluid dynamics.

`wave_density_metric` checks this algebra exactly for rational coefficient
inputs. The anisotropic fixture tests the algebra; the native increment
readout below obtains finite averaged coefficients from source records.
Controls verify
spatial chart covariance and the change `t=kappa*tau` of wave-clock
coordinate. Changing the measure normalization preserves the characteristic
cone but changes the metric normalization; replacing `a` by a round tensor
merely because both have rank three changes the cone. The form readout,
measure and clock therefore carry information absent from the slow-band rank.

This establishes a **measurable Lorentz-signature description of the declared
product wave**. The radial domain, boundary condition and active wave
protocol are construction inputs. No smooth native coefficient field,
unique whole-sequence limit, curvature, Einstein dynamics, interacting QFT
or physical radial transport follows. If the **whole coefficient density**
`C=mu diag(-1,a)` has a smooth nondegenerate representative in a fixed
smooth product atlas, the displayed `g` is a smooth static Lorentz metric
and the same product wave is exactly its scalar wave equation, with the
declared radial boundary condition. This requires smooth positive `mu`
as well as smooth uniformly elliptic `a`. The static line element then
gives local proper-duration factor
`N(x)=sqrt(mu(x)) (det a(x))^(1/4)` relative to its wave coordinate; identifying
an actual clock apparatus with that line element is a separate measurement
contract. A faithful global thermal state's modular flow is likewise not
automatically the modular flow of every spatial patch.

The fixed native [radial-projection measure](../physical_space_time/OCCUPATION_FIELD.md#a-continuum-random-distribution-on-the-carrier-sphere)
does **not** satisfy this smooth-density premise in the standard smooth
sphere atlas. Let `n_F` be an outward unit face normal, `r` the polyhedron
inradius, `a0` one flat face's area, and `u` a spherical direction. On face
`F`, the source point is `x=r u/(n_F.u)`, so the retained Jacobian gives
`rho_F=dmu/domega=c/(n_F.u)^3`, `c=pi r^2/(5 a0)>0`.
At an open edge shared by distinct faces `F,G`, their support values
agree: `n_F.u=n_G.u=d>0`. Their density values agree, but their one-sided
spherical gradients differ by `-3c(n_F-n_G)/d^4`, a nonzero tangent vector.
The density is continuous and piecewise smooth but not `C1`. Changing its
values on the edges cannot make a smooth representative: continuous
representatives equal almost everywhere must agree everywhere. Since
`C^{tt}=-mu`, smoothness of `a` alone cannot repair this obstruction.
An alternate smooth observation dictionary needs its own measure and
operator proof; this is not a topological obstruction to every atlas.

**Smooth effective descriptions of the actual product readouts do exist.**
This statement does not require the preceding exact smooth-density premise.
Use the smooth reference volume `nu=dr/L dOmega/(4pi)` on the same compact
radial cylinder, and write `dmu=rho dnu`. The native `rho` is positive and
globally Lipschitz: it is a positive constant times
`(max_F n_F.u)^(-3)`, whose denominator is bounded away from zero.
Mollification and normalization in a fixed finite angular atlas give smooth positive `rho_epsilon`
converging uniformly to `rho`, with uniformly bounded gradients and
gradient convergence almost everywhere along a subsequence. Smooth
`a_epsilon=diag(1,a_screen_epsilon)` converge almost everywhere to `a`
with common ellipticity bounds. Both choices preserve the radial domain
and Dirichlet boundary condition, and their screen coefficients are independent of `r`. Set
`H_epsilon=L2(rho_epsilon dnu)`, `H=L2(rho dnu)` and identify them by
the explicit unitary `U_epsilon f=sqrt(rho_epsilon/rho) f`.
The corresponding smooth-theory smearing is `U_epsilon^-1 f`; this local
dictionary preserves support and converges to the original readout.

The transported energy forms Mosco-converge to the native product form.
Here are the measure-dependent details of the
[elliptic coefficient-convergence argument](https://arxiv.org/pdf/1412.0725).
Put `w_epsilon=sqrt(rho_epsilon/rho)`. Then `w_epsilon^±1->1` uniformly,
their gradients are bounded and tend to zero almost everywhere, and
multiplication by `w_epsilon^±1` converges strongly to the identity on
`H1`. The same multiplication error tends to zero in `H1` for a bounded
`H1` sequence after taking a Rellich-convergent `L2` subsequence: split
`f_epsilon grad w_epsilon` into
`(f_epsilon-f) grad w_epsilon+f grad w_epsilon`.
Thus transporting the forms introduces no limiting derivative term.
Dominated convergence gives recovery on each fixed form vector. For the
lower bound, use weak `H1` compactness and the dual energy formula with
`(rho_epsilon a_epsilon)^(-1)`; bounded ellipticity and almost-everywhere
convergence pass every fixed dual test to the limit. Zero radial-end
traces are preserved throughout.

Consequently `Btilde_epsilon=U_epsilon B_epsilon U_epsilon^-1` converges
in strong resolvent sense to `B`. This convergence is also in **norm
resolvent**: all shifted resolvents map the unit ball into one bounded
`H1` set, which is compact in `H`. If unit inputs converge weakly,
self-adjointness and strong resolvent convergence identify every compact
output limit with the limiting resolvent applied to that weak limit;
compactness of the limiting resolvent then gives the norm assertion.
The angular constant is retained, and the common Dirichlet radial gap
`B_epsilon,B >= (pi/L)^2` prevents an infrared singularity.

Functional calculus gives operator-norm convergence of
`q_b(Btilde_epsilon)` and of the ordered field covariance and sine
response, uniformly on compact wave-time intervals. These spectral
functions vanish at infinite energy. Cosine evolution converges strongly,
uniformly on such intervals for each fixed input. Momentum covariance is
unbounded; for fixed smearings in the `H1` form domain with zero radial-end
traces, its spectral tail above energy `K` is
bounded by a constant times their uniform form energy divided by `sqrt(K)`.
This yields convergence of the full smeared thermal CCR covariances and
symplectic forms. Every fixed finite word of Weyl readouts therefore
converges, as do finite adaptive Ramsey instruments formed from finite Weyl
sums and Weyl feedback. All their classical records and any fixed finite
detector/reference outputs are included; arbitrary sharp projectors or
trace-norm comparison of the entire surviving field are not asserted.

Each smooth coefficient density `C_epsilon=mu_epsilon diag(-1,a_epsilon)`
defines an exactly smooth static `1+3` Lorentz metric by the displayed
density identity. Its free scalar thermal theory consequently approximates
**any fixed finite menu** of these native-product field and clock readouts
at bounded wave times
to any requested tolerance, with one common sufficiently small `epsilon`
for the entire menu. The source preparation, radial design and active-clock
contracts are unchanged. Keep the same native clock tapes and controller;
their source, timing and timeout allowances add by the triangle inequality
to the uniform smooth-description error, without choosing a geometry per
readout. This is an existence and convergence result,
without a computable smoothing rate for the unknown limiting tensor.
It supplies neither one smooth metric exact for every observation nor
smoothness of the native limit, curvature or optimal-front convergence,
an interacting spatial continuum, or laboratory identification.

The [same two-mode quartic projection](../interacting_band_extension/README.md#interacting-readouts-on-the-same-radial-product)
also has a compatible smooth-background interpretation. On `(0,pi)` at
fixed covariant coupling `lambda`, its effective coupling is
`g_epsilon=lambda/pi integral mu_screen,epsilon^2 sqrt(det a_screen,epsilon) dx`.
Bounded ellipticity and dominated convergence give `g_epsilon->g`.
At each fixed occupation cutoff, the Hamiltonian, its Gibbs preparation
and bounded-time complete detector histories are continuous in this
coupling. First choose an occupation cutoff using the proved fixed-two-mode
limit, then one sufficiently small `epsilon` for that cutoff and the whole
observation menu. This diagonal comparison gives smooth-background
regulated interacting descriptions of those observations. Their state is
the Gibbs state of the **interacting** Hamiltonian; the free thermal state
belongs to the separate free-field protocol. The spatial two-mode cutoff
is retained, with no error bound against a full interacting local QFT.

### Native records determine the limiting tensor weakly

The tensor can be read from the same full-port density kernel, without
fitting an elliptic operator. Set `s=4^l`, `M_l=630s`, `N_l=2M_l s`, and
`K_l=I-L_l/(2M_l)`. A disclosed **total-one calibration** of the unchanged
native endpoint rule has exactly this transition kernel: each occupied
seam contributes one move and one hold among the `2M_l` equally likely
seam/coin records. Empty seams contribute two holds. Currents reconstruct
the terminal location from the disclosed start. This calibration is
distinct from tagging a particle in the fixed-number density experiment;
that experiment independently gives the same `K_l` in its two-point
covariance after the susceptibility and conserved constant are restored.

For cell probes `f_l,g_l`, repeat their values at every port, evolve for
`n_l=floor(t N_l)` native attempts, and average the twelve starting ports
in each cell. Denote the resulting positive kernel by `Pbar_l(t)`; it
preserves constants but is not a cell semigroup. For `t>0`, define

```text
Gamma_l,t(f,g) = [Pbar_l(t)(fg)-f Pbar_l(t)g-g Pbar_l(t)f+fg]/(2t).
```

This is the conditional mean of `Delta f Delta g/(2t)`, so every finite
matrix of probe increments is positive semidefinite. In a finite run use
the exact separation `t=n/N_l`; no fractional event is introduced. The
source full-port resolvent/heat theorem gives `Pbar_l(t)->exp(-t Abar)`
strongly along the existing spatial subsequence. For bounded `H1` probes
and bounded `H1` weight `phi`, cell averages converge in `L2`, as do their
products. Hence, first refining at fixed positive `t`, then taking `t->0`,

```text
lim_(t->0) lim_l integral phi Gamma_l,t(f,g) dmu
  = [E(f,phi g)+E(g,phi f)-E(fg,phi)]/2
  = integral phi a_screen^ij partial_i f partial_j g dmu.   (NR)
```

The first equality follows from `(I-exp(-t Abar))/t` in the quadratic-form
pairings; all three products belong to `H1`. No differentiation of a
rough coefficient is used. Choose cutoff coordinate probes `f_i=chi x_i`
with `chi=1` near `supp(phi)`. Dividing (NR) by `integral phi dmu` recovers
the weighted chart average of `a_screen^ij`. Shrinking nonnegative chart
weights then recovers the coefficient almost everywhere by Lebesgue
differentiation. This identifies the measurable tensor of each selected
limit, without a whole-sequence rate, pointwise finite-run accuracy,
smoothness or isotropy assumption.

The order of limits is essential. On an original face use barycentric
coordinates `(x,y)` with corners `(0,0),(1,0),(0,1)`. The unrelaxed
one-attempt increment tensor on its interior refined cells is exactly

```text
a_bare = [[1/36,-1/72],[-1/72,1/36]],
(3/7) a_bare <= a_screen(x) <= (15/19) a_bare   almost everywhere. (NB)
```

The lower inequality is the existing `H1` coercivity bound expressed in
this chart. For the upper inequality, sample a smooth compactly supported
probe at cell centroids. Its unweighted native trial energy converges by
the three centroid displacement directions to `integral grad f.a_bare.grad f`.
The source terminal-star conductances satisfy `w<=15/19`; their harmonic
lifts have no greater native energy and their cell means converge to the
same probe. Mosco liminf therefore bounds `E(f)` by `15/19` of this trial
energy. The cutoff sine/cosine localization argument above gives (NB)
pointwise. Taking the one-attempt derivative before refinement thus reads
a strictly different, unrelaxed tensor; it cannot determine `a_screen`.

`native_chart_increment` executes continuous cutoff coordinate probes on
source levels one and two, retaining all holds. Exact two-attempt endpoint
composition independently checks its polarized kernel calculation. A
cell-reset control assigns probability `1/3657830400` to a two-cross-seam
path whose true probability is zero: the source matching needs an internal
move between distinct cross seams. The mean-eliminated Schur generator
also has a positive off-diagonal entry and cannot replace this positive
kernel. One- and sixteen-attempt runs at levels one and two share
`t=1/20160`; these are finite comparisons, not evidence of an asymptotic rate.

The finite averaged tensor feeds `wave_density_metric` with declared
radial coefficient one and product measure `dr dmu=dr dx dy/10`. It yields
an exact Lorentz-signature **comparison metric for those averaged
coefficients**. A metric of an average is not an average of the metric.
The limiting metric uses (NR) and the subsequent shrinking-chart limit
before the nonlinear density-to-metric map. The radial domain and active
wave law are supplied; neither native density time nor this calibration
is identified with wave time, a modular parameter or laboratory seconds.

### Intrinsic distance is the product wave's first response

The recovered energy determines more than a comparison propagation bound.
For the declared product `X=(0,L) times screen`, retain the full angular
constant, radial Dirichlet endpoints and the same operator
`B=-d_r^2+(A⊕0)`. Write `Gamma_B(u)=grad u.a.grad u`, where
`a=diag(1,a_screen)`. For measurable sets `E,F` of positive finite measure,
define the **essential intrinsic set distance**

```text
d(E,F) = sup_u [essinf_F u - esssup_E u],
u bounded in the local form domain, Gamma_B(u)<=1 almost everywhere.
```

Use the local domain on compact interior exhaustions, not globally
boundary-zero test functions. This matters at Dirichlet endpoints.
The strongly local regular form and its uniform ellipticity are proved; (NR) determines its angular energy density. The radial factor is
the declared one. No smooth coefficients or pointwise heat kernel are needed.

The form chain rule gives
`||exp(lambda u) exp(-tau B) exp(-lambda u)||<=exp(lambda^2 tau)`
for every admissible `u`. Optimizing the resulting separated-support bound
over `lambda` and then `u` gives the Davies–Gaffney estimate
`||1_F exp(-tau B)1_E||<=exp(-d(E,F)^2/(4tau))`.
Its pairwise heat-to-wave implication
([Coulhon–Sikora, Theorem 3.4](https://arxiv.org/pdf/math/0609429))
gives `1_F cos(t sqrt(B))1_E=0` for `|t|<d(E,F)`; the same holds for
`1_F B^-1/2 sin(t sqrt(B))1_E`. The integrated short-time theorem of
[Ariyoshi–Hino, Theorem 2.7 and Remark 2.8](https://www.maths.tcd.ie/EMIS/journals/EJP-ECP/article/download/286/286-793-1-PB.pdf)
also applies to this symmetric local Dirichlet form, including its
nonconservative Dirichlet boundary. Apply their normalization to `E_B/2`
and set their time to `2tau`. Their essential set distance is exactly the
one above, so

```text
lim_(tau->0) -4tau log <1_F,exp(-tau B)1_E> = d(E,F)^2.     (OF1)
```

This proves that the propagation bound is optimal for complete local
readouts. Indeed spectral calculus gives
`exp(-tau B)=(pi tau)^-1/2 integral_0^infinity exp(-t^2/(4tau)) cos(t sqrt(B)) dt`.
If the cosine block vanished for all `0<t<R` with `R>d(E,F)`, this identity
and `||cos(t sqrt(B))||<=1` would bound the positive heat pairing by
`sqrt(m(E)m(F)) erfc(R/(2sqrt(tau)))`. Its logarithmic exponent contradicts
(OF1). Therefore, with an empty infimum interpreted as infinity,

```text
inf{t>0: ||1_F cos(t sqrt(B))1_E||>0}
 = inf{t>0: ||1_F B^-1/2 sin(t sqrt(B))1_E||>0}
 = d(E,F).                                                (OF2)
```

The sine equality follows by integrating the cosine and, conversely,
differentiating a zero sine block in the strong topology. An infimum need
not be attained. The result quantifies over the full local `L2` spaces;
one fixed signed probe and detector can cancel. In the free-field encoding,
the same sine block is its retarded commutator, so (OF2) is also a
state-independent onset of possible local linear response. Finite spectral
cutoffs need not have an exact support cone. A growing source apparatus
realizing all those local probes and modes is a separate implementation
obligation.

This distance is consistent with the reconstructed static metric:
`g_spatial/(-g_00)=a^-1` almost everywhere. For constant coefficients and
time weight `w`, the squared optical separation is `w dx.a^-1.dx`; its
dual covector saturates `xi.a.xi=w`. Changing the wave coordinate by
`t=kappa*tau` changes the onset to `d/kappa`. A modular clock satisfying
`sigma_s=alpha_(-b s)` for **this same field and algebra-state pair** gives
onset `|s|=d/b`; a different prepared marginal's modular flow supplies no
such identification. Neither statement calibrates seconds.

`native_optical_front_control` checks optical duality, chart covariance
and clock scaling on the finite native increment tensor. Its full-port
240-state eigenmode of eigenvalue six is zero outside one carrier, so
every remote detector misses that signed mode at all wave parameters.
A positive preparation crosses the same seam in one native attempt.
The continuous finite-matrix cosine also has a cross-seam term `t^2/2`;
its nonzero small-time tail cannot be substituted for the limiting front.
On the declared interval `(0,4)`, restored angular-constant probes on
`E=(1,3/2)` and `F=(2,5/2)` have distance `1/2`: for `0<=t<=1`, the
cosine pairing is `(t-1/2)_+/2` and the sine pairing `(t-1/2)_+^2/4`.
Independent Dirichlet spectral sums check these exact readouts with
explicit tail bounds. Removing the angular constant erases this control.

For rough tensors use the essential form distance, not arbitrary
representative-dependent curve lengths. Multiplying a constant tensor
by four only on a line changes no energy, heat or wave operator: that
line has volume measure zero. Yet a naive length integral along that
line is halved. This control excludes an unsupported identification with
smooth null geodesics. The achieved join is native weak energy readout,
intrinsic distance and optimal response of the **declared product wave**.
Radial geometry, active dynamics and clock/readout access are its contracts; fastest primitive transport and observed spacetime are not
identified by (OF2).

The existing angular record menu does not select that radial extension.
Let `L` be the full-port generator reconstructed from native endpoint
records, `R` any connected radial graph Laplacian with `R 1=0`, and
`J f=1_radial tensor f` with the radial constant normalized. Then

```text
B_R=I tensor L+R tensor I,     B_R J=J L,     J*J=I,
J* F(B_R) J=F(L)                                             (RS)
```

for bounded spectral functions `F`. The first equality on the invariant
subspace proves (RS) by the resolvent identity and spectral calculus.
Thus heat, cosine and positive shifted-resolvent readouts of this angular
sector agree for inequivalent connected local extensions. The continuum
version uses any connected radial Neumann interval or circle: their
normalized constant has eigenvalue zero, but their lengths and other
eigenvalues differ. These are radial-**average** observations, not ordinary
boundary Dirichlet-to-Neumann data. No equality is claimed for all radial
states, extra interventions or a full thermal state's normalization.

<a id="radial-nonidentifiability-under-complete-histories"></a>

The same intertwiner extends this statement to complete adaptive histories.
Put `U_R(t)=exp(-it sqrt(B_R))` and `U(t)=exp(-it sqrt(L))`. Spectral calculus
gives `U_R(t)J=JU(t)`. For any angular instrument or feedback operator K,
its lift `I_radial tensor K` satisfies `K_R J=J K`. Induction through each
branch word therefore gives

```text
W_R,h J = J W_h,
(W_R,h tensor I_ref) (J tensor I_ref) = (J tensor I_ref) (W_h tensor I_ref).
```

The branch h may contain native clock marks from the same clock instrument
in each extension, earlier measurement outcomes, record-selected durations
and feedback. For every angular/reference input
embedded by J, the complete classical-history and decoded angular/reference
output is identical. Keeping all stopped branches includes timeouts; extra
clock records do not select R. This statement requires angular-only controls
and the constant-radial input code. Nonconstant radial inputs, controls that
leave the code and a full-product Gibbs preparation are outside its scope.

The executable control uses the actual 240-port L, all 3600 two-attempt
native clock tapes, zero-current incidents and all 4900 tape/field outcomes.
The timer pauses during comparison services. Calculations group tapes only
when their conditional work state has the same incident count and field
bits; their weights sum the complete classical direct-sum output.
It propagates a two-column entangled reference, measures angular port parity,
uses a native-seam swap as first-bit feedback and chooses the second wave
duration `1/4` or `1/2` from that bit. Every stopped output agrees across the
three extensions below; deleting feedback changes the record law. These
spectral waves and angular instruments are comparison operations, not
additional executed source pulses or a physical clock identification.

`radial_nonidentifiability_control` uses the actual 240-port `L` and three
three-site radial graphs: a unit-weight path, a weight-four path and a
unit-weight cycle. Their radial spectra are `(0,1,3)`, `(0,4,12)` and
`(0,3,3)`, while all three exact 720-port intertwiners satisfy (RS).
A specified additional measurement distinguishes them: prepare the
normalized constant angular mode times radial contrast `(1,-1,0)/sqrt(2)`
and measure its cosine response in the **same angular-clock parameter**.
The second derivatives at zero are respectively `-5/2`, `-10` and `-3`.
At time `1/100`, the remainder `||R||^2 t^4/24` gives three disjoint exact
cosine-response intervals, so the separating probe works at finite time.
This identifies the missing information as access to a nonconstant radial
preparation and response, beyond the retained angular menu.

The original native histories and incident marks do not contain `R`;
these comparison extensions neither alter their law nor execute extra
native repair seams. Likewise the supplied length in the existing radial
record filter changes its decoded outputs, not its input current tape.
The result is underdetermination from this record menu, not a refutation
of the Dirichlet product or OPH. Rank-three carrier coordinates and `S2`
directional topology do not remove this ambiguity. Conditional on a
chosen product operator and compatible active clock, (OF2) does discharge
a separate propagation-speed postulate: the same form fixes its optimal
operator front and its conversion to that clock's parameter. Selecting
physical radial access and laboratory units requires further readouts.

### Native observer-centered radial records

The existing carrier also admits a radial readout without adding radial
sites. `recover_radial_records` receives labelled total-one preparations,
seams, coins and signed currents from all 720 native branches. The current
reconstructs the occupied endpoint; the declared uniform-seam/fair-coin
law fixes its row kernel `K=I-L/60`. Breadth-first distance from a chosen
observer port partitions the twelve ports into shells of sizes `(1,5,5,1)`.
No hidden final load or supplied radial graph enters this decoder. Complete
branch enumeration is exact source evidence, not a finite-sample noise bound;
preparation and complete current access are the same calibration contracts
as the existing structural readout.

Let `U` lift four shell values to the ports and `M=diag(1,5,5,1)`. Direct
native branch counts give

```text
       [ 5 -5  0  0 ]
R =    [-1  3 -2  0 ],       L U=U R,     U*U=M,
       [ 0 -2  3 -1 ]
       [ 0  0 -5  5 ]
J=U M^(-1/2),  J*J=I,  L J=J (M^(1/2) R M^(-1/2)).        (NRQ)
```

Every port in a shell has the same total transition weight to each other
shell. This is strong lumpability: multiplication with intermediate shell
projectors proves the complete shell-history law for every initial port,
not only its mean. The normalized radial generator is self-adjoint and has
eigenvalues `0,5-sqrt(5),6,5+sqrt(5)`. (NRQ) intertwines its heat, resolvent
and coherent cosine readouts with the **same native operator**. It is the
record-level radial specialization of the existing finite spectral theorem.
Exchanging a port from shell one with a port from shell two violates the
block-sum equalities: averaging such a partition would erase memory.

Mark every seam incident on the observer, including zero current. The
remote and incident substochastic row kernels also lump separately into
`Q0,Q1`; their row sums are `5/6,1/6`. Consequently every joint marked-shell
history closes, and the first-incident transition is
`T=(I-Q0)^(-1) Q1`, with timeout mass `(5/6)^N` after `N` attempts.
Its commutator satisfies `[T,R]_(0,1)=13/37`, excluding any scalar spectral
time change of `R`. Independent native two-attempt enumeration tests every
initial port, intermediate shell, final shell and pair of incident marks.

For the separately admitted protected coherent walk, one incident trigger
applies one step while the timing source pauses during service. The exact
720-arc native reversal/reflection on an input at the observer gives both
Ramsey outcomes; the plus probabilities after zero through three triggers
are `1,23/24,101/120,7271/10800`. These are `T_n(K)` compressions under the
[existing source walk](../observer/COHERENT_CONTINUATION.md), with its
preparation, reflection and readout services. The control executes rational
arc operations; it does not rerun source pulses for this particular input
or identify the passive stopped kernel with that active wave.

This radius describes shells **within the same twelve-port carrier**. Its
four-dimensional function space is a subspace of the twelve-dimensional
source space, not an independent factor that can be tensored with the full
angular field. It therefore supplies an operational finite radial readout
but does not select the extra interval in the product construction, a
continuum radial metric, or physical length units. Additional product-radial
observations in (RS) retain that distinct obligation.

### Transport through source refinement

The response-recovered carrier Gram can be transported into the actual level-3
and level-5 source screens without changing its twelve root-port labels. The
screens have exactly the same induced native-seam graph through nine hops from
the corresponding root carrier: 222 ports and 561 seams. The exact unitless
chord metric on the root twelve ports is therefore identical in both screens.
The root graph-hop ball counts at radii zero through ten are
`12, 15, 30, 45, 54, 84, 114, 132, 186, 222, 243`; the first induced adjacency
difference appears at radius ten. The separate
[`native_action` refinement packet](../native_action/INCIDENT_REFINEMENT.md)
bounds how much this finite boundary difference can affect the complete local
repair output at a declared incident-event horizon.

An exact breadth-first readout covers every retained federation level, 0
through 5. The largest root-port distances are 13, 31, 60, 121, 242 and 484
native seam hops. At level 5, the port-ball counts at radii 16, 32, 64, 128 and
256 are 576, 2,250, 8,805, 37,224 and 123,926; the corresponding exact
doubling ratios are `125/32`, `587/150`, `12408/2935` and `5633/1692`. This is
roughly quadratic growth across part of this finite range, followed by
saturation on the finite sphere. It is a property of native seam-hop balls,
not evidence for a calibrated three-dimensional physical volume law. The
rank-three metric applies to each twelve-port carrier; it has not been glued
into a bulk metric.

This transports an operational neighborhood and one local metric, not a
Euclidean coordinate chart over the 222-port ball. The existing federation
ancestry map, supplied by the source rather than recovered from observer
records, copies each local port label to every child carrier. It preserves
the carrier-cell adjacency, but it does not map every fine port seam to a
coarse port seam: between levels 0 and 1, 74 fine seams map to non-seams and
two collapse; between levels 1 and 2, 280 map to non-seams and four collapse.
The exact per-edge Laplacian intertwining defect is one at both steps. Thus
this map supplies neither a path-independent carrier-frame transport nor a
scale law for a spatial metric.

The federation builder also has an inter-carrier routing convention: it uses
the source sphere's cell-centre directions and a declared local tangent-frame
rule to choose which ports are joined. The builder marks that routing
`source_derived: false`; it is not inferred from native observer records. At
levels 3 and 5, respectively, only `113/160` and `5951/7680` of inter-carrier
endpoint pairs are antipodal-consistent with the common local port labels.
This convention constructs the seam list, but it does not supply a
record-derived frame-transition law or a scale factor. No path-independent
physical coordinates follow from it.

The obstruction is also visible without extrapolating the ancestry map. As a
non-derivation control, take an actual inter-carrier seam at the level-3 root
and additionally require its endpoint distance to equal the local
shortest-shell length. Two independent orientations of its neighboring
carrier can satisfy that same link-length constraint, preserve both carriers'
exact internal Gram matrices, while changing other cross-carrier Gram entries. This equality is used only in the control; it is not a model input.
The current repair records identify the seam endpoints but do not record the
missing cross-carrier inner products. So the substrate supplies exact
native-hop neighborhoods and a local unitless carrier metric; it does not supply a global physical space, inter-carrier ruler, or continuum geometry.

The source-order/frame packet in
[`Geometry/SourceOrderFrameCompatibilityPacket`](../../Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean)
needs its spatial placement, cone attachment and frame premises. This
packet adds no Lorentz cone, spacetime metric, signal-speed law or continuum
limit. Similarly, the finite representation in
[`QFT/SourceHistoryGNSDynamics`](../../Lean/QFT/SourceHistoryGNSDynamics.lean)
has its own supplied operator extension and flow parameter; record ordering
does not identify that parameter with these local ticks.

Absolute native fields have an observer-visible readback construction under
an explicit preparation/access contract, compatible local restrictions,
complete continuing records, and a testable clock comparison domain. The
single-carrier structural Gram is invariant under the admitted changes of
calibration exposure; the native dynamical slow band need not be. Connecting
these field and clock readouts to one effective spacetime description is an open obligation, with the joint observable laws and their premises held
fixed. These are finite model results without an empirical comparison.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry --check model/relative_observer_geometry/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry.intrinsic_geometry --check model/relative_observer_geometry/intrinsic_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry.binary_dynamics --check model/relative_observer_geometry/binary_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry.temporal_records --check model/relative_observer_geometry/temporal_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry.clock_reconstruction --check model/relative_observer_geometry/clock_reconstruction_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/relative_observer_geometry
```

Each result binds its producer and relevant native sources by SHA-256;
`results.json` also binds its runtime and cited Lean interfaces.
`manifest.json` binds every packet file.
