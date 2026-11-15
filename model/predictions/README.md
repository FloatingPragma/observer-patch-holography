# Fixed observer predictions and complete-count comparison

This package distinguishes two specified finite observer instruments using
an interference measurement. Both apply the same native repair to classical
loads. The current-record instrument C and contrast-record instrument E give
different probabilities when the source preparation is coherent. The settings,
sample sizes, error allowances and decision rule are fixed before comparison.

The package also audits the inherited P-dependent electromagnetic endpoint
against published data and replays a [public Bell-record comparison](BELL_RECORDS.md).
The latter constrains a local-factorized extension of the diagonal recorded
branch, under the experiment's stated setting and locality contract. These
retrospective calculations and the prospective finite-instrument protocol
have different evidential roles.

The [complete repair identification](../repair_identification/README.md)
extends the two-coordinate phase probe to all nine real instrument
coordinates. The [compiled-history forecast](../compiled_clock_history/README.md#calibrated-repair-law-and-complete-readout-predictions)
then carries calibration uncertainty through a complete two-tick experiment
with feedback and source-emitted apparatus. Under the stated synthetic
probability tubes, its two intervention settings retain a resolved joint
clock/detector contrast. The [resource-clock forecast](../laboratory_clock/README.md)
additionally predicts deadline-conditioned detector events under a joint
repair/pulse-duration calibration region. These are conditional source
predictions, not new discrepancies with SM + GR, externally registered
predictions or physical confirmations. The
[public-data audit](../../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
tests a direct event-count/optical-phase identification against public
Ar13+/Yb+ records. Its complete-block range test charges missing data and
calendar-window selection. It rejects that simple attachment under its
explicit bounded-error and readout assumptions; collective clocks and
interacting physical identification require their own dictionary. The
[regional extension](#collective-regional-optical-attachment) derives all
shared-event collective clocks on one carrier and tests every distinct pair
under the same direct readout hypothesis. No clock definition needs
external approval, and no atomic transition is derived by choosing a unit.

The [simultaneous source-clock limit](../relative_observer_geometry/README.md#incident-clocks-on-the-simultaneous-density-continuum-scale)
also predicts joint clock jitter and density observations. For fixed bounded
incidence patterns, centered incident counts have Gaussian finite-dimensional
covariance `min(t,u)*c_ij`, where `c_ij` counts shared seams. Adjacent degree-five
ports have limiting correlation `1/5`; the degree-six endpoints of an external
seam have `1/6`. Conditional native stirring proves asymptotic independence
from density snapshot vectors, while exact finite nonlinear clock–field
dependence is nonzero. Uniform inverse-clock bounds let the same local clocks
read the proved continuous density paths. These predictions use the declared
stationary preparation, IID seams and observation access. Public optical ratios
do not supply the joint density readouts or identify this source with atomic
physics; no additional empirical verdict follows.

## Run

From the repository root:

```sh
# Replay the retained exact predictions, statistical checks and endpoint audit.
python3 -m model.predictions --check model/predictions/results.json
python3 -m model.predictions.bell_records --check model/predictions/bell_results.json
python3 -m model.predictions.optical_clock --check model/predictions/optical_clock_results.json
python3 -m model.predictions.regional_optical --check model/predictions/regional_optical_results.json
python3 -m model.predictions.modular_optical --check model/predictions/modular_optical_results.json
python3 -m model.predictions.repair_response --check model/predictions/repair_response_results.json

# Produce an unfilled count template; the destination must not exist.
python3 -m model.predictions --template /tmp/observer-counts.json

# After independent acquisition supplies every count, evaluate the file.
python3 -m model.predictions --evaluate /tmp/observer-counts.json

# Exercise the decision rule with explicitly labeled synthetic data.
python3 -m model.predictions --demo C
python3 -m model.predictions --demo E

python3 -m pytest -q -W error model/predictions
```

The template contains `null` counts and cannot pass validation. Every setting
requires exactly 8,192 attempted trials, with all twelve common
`(coin,current,outcome)` bins and an abort count supplied explicitly. Counts
must be nonnegative integers, including explicit zeros. The evaluator rejects
missing or duplicated bins, missing trials, a different protocol hash and
undeclared fields. It never divides by the number of successful trials.

## Fixed apparatus and predictions

Use the [native observer's](../observer/README.md) 78-dimensional load carrier:
twelve nonnegative integer port loads with total load two. Ports use indices
0 through 11. Let `x={2,3}` mean one load at port 2 and one at port 3, and
let `z={1,2}` have the analogous meaning. Prepare

```text
|phi_theta> = (|x> + exp(i theta)|z>)/sqrt(2).
```

Apply one repair on seam `(0,1)`, retaining its record. Measure the binary
source projector `P_plus=|phi_0><phi_0|` and its complement. Readout acts
inside each committed-record block; it does not erase or interfere different
classical record values. Fresh preparations supply independent trials.

| Fixed preparation | C: probability of plus | E: probability of plus |
| --- | --- | --- |
| Phase 0 | `5/8` | `3/8` |
| Phase pi | `1/8` | `3/8` |
| Phase pi/2 | `3/8` | `3/8` |
| Equal incoherent mixture of x and z | `3/8` | `3/8` |

The phase-reversal contrast is `1/2` for C and zero for E. All four
preparations have identical classical load probabilities. Consequently a
classical endpoint comparison alone cannot distinguish these instruments.

The two columns are special cases of one formula. For every quantum record
law compatible with the native seam rule, the plus probability is
`3/8 + Re(g e^{-iθ})/2`, where `g` is the law's coherence between the empty
seam of `x` and the waiting single load of `z`; C has `g=1/2` and E has
`g=0`. The [record-law classification](../record_laws/README.md) proves this
and gives the other distinguished laws: recording only the current gives
`g=1/√2`, and recording the full seam occupation gives `g=0`. The aggregate
plus frequencies of the four settings therefore estimate `Re g` and `Im g`
for any admissible law, while the locked decision rule below compares the two
executed laws.
[DERIVATION.md](DERIVATION.md) gives the complete joint law and the exact
channel calculation.

The [instrument code](observer_instrument.py) derives these probabilities
from the native endpoint map. It exposes every record, outcome and
unnormalized terminal source density, including zero-probability events.
E retains its additional contrast label internally. The shared comparison
interface sums probabilities over that label; this does not recover C's
quantum coherence.

## Error allowances and verdicts

[protocol.json](protocol.json) binds the predictions to the implementation
hashes. For every setting the plan fixes:

- 8,192 independent attempted trials, with no optional stopping;
- a sampling half-width of `1/32`;
- an independently established apparatus probability-error budget of `1/64`;
- a familywise sampling-error requirement below `0.001`.

Four settings each contribute twelve joint-bin frequencies, the aggregate
plus frequency and the abort frequency: 56 statistics. An exact rational
certificate bounds the simultaneous sampling failure probability by less
than `0.000012606`. Each observed interval is compared with the candidate's
ideal probability widened by the apparatus budget. One disjoint interval
rejects that candidate under the declared contract. Testing C and E against
the same simultaneous intervals needs no additional multiplicity factor.

The apparatus budget covers preparation, repair, storage and readout errors
in each scored probability. It is a physical premise requiring evidence
independent of the comparison counts. If its calibration guarantee has
failure probability beta, add beta to the stated sampling bound for an
unconditional error guarantee. A count file cannot establish that calibration,
trial independence or its own provenance.

An abort is an attempted trial without a complete valid record and readout.
It counts in the denominator and has ideal probability zero. A statistically
resolved excess over its allowance returns `apparatus_budget_exceeded`, so
the evaluator cannot attribute the discrepancy to the ideal repair law.

Without supplied counts the verdict is `no_data`. Software demonstrations
are marked `synthetic`; external counts are marked `external_unverified`.
Both report `empirical_confirmation=false` and `provenance_verified=false`.
If one candidate survives and the other is rejected, the counts discriminate
the conditional instruments. If both survive, this test does not separate
them. If both are rejected, neither specified apparatus model explains the
counts within the declared premises. Survival means `not_rejected`, not proof.

## Physical scope

The native load interpretation, coherent preparation and source-projector
measurement are explicit apparatus choices. Both C and E are ordinary
quantum instruments. Their separation is a test of the chosen record law,
not a novel prediction distinguishing OPH from quantum mechanics.

Neither finite instrument result has a valid comparison with publicly
available physical data. The C/E protocol requires a specific coherent
12-port preparation, native fair-coin repair, committed coin/current record,
source-projector readout and complete attempted-trial counts; no matching
dataset or calibrated implementation is supplied. The
[incident-clock field construction](../observer_time_fields/README.md) predicts
local laws in dimensionless incident-repair counts, while the
[screen comparison](../native_action/INCIDENT_REFINEMENT.md) is conditional on
source port maps, matched initial local/reference states and a particular
incident-event record interface. There is no physical calibration connecting
those quantities to published laboratory time, position or field readouts.
The ideal modular-phase target is a supplied reference in a separate finite
model, and the CODATA endpoint comparison below is retrospective; neither is
an empirical test of these incident-clock results or the C/E protocol.

## Collective regional optical attachment

[regional_optical.py](regional_optical.py) extends the direct optical test to
every pair of distinct event sets generated by nonempty port regions on one
native carrier. The [collective clock theorem](../laboratory_clock/README.md#collective-clocks-on-the-same-carrier)
derives the full renewal law, its mean ratio `b/a`, and its residual fluctuations.
No independent copies, noise parameter or atomic model is introduced. The
tested attachment is

```
reported 300-second block ratio = gain*a*B(N)/(b*N) + error,
gain >= 1,    |error| <= 10^-12.
```

Here `N` is one fixed integer number of A ticks per block in the declared
`300*[642121496772644,642121496772646]` interval. The frequency bracket is a
hypothesis, not a hard bound inferred from a quoted standard uncertainty.
Each consecutive candidate calendar block, including every unobserved block,
must contain exactly N consecutive A intervals. This makes the block count
laws independent under the declared IID scheduler.
The normalization `a/b` follows from the count law and permits unequal native
clock rates. The clocks require ordered regional tape unions and coincidences.
The hypothesis includes no additional temporal filtering.

Rejection uses a point-probability bound, not variance alone. If the geometric
part has failure probability `q`, it is the sum of an independent
`Bernoulli(q/(1+q))` and twice a geometric variable with failure `q^2`.
The larger variance `v` of this Bernoulli and the terminal `Bernoulli(c/a)`
is strictly positive for every distinct event pair. Convolution and Fourier
inversion give `sup_k Pr(B(N)=k) <= sqrt(pi/(8*N*v))`: use
`|1-r+r*exp(it)|^N <= exp(-2*N*r*(1-r)*sin(t/2)^2)` and
`sin(|t|/2)>=|t|/pi` on `[-pi,pi]`, then enlarge the integral to the real line.
The code rounds outward using `pi<22/7` and rational arithmetic.

For an observed range `w`, the admissible B-count interval contains at most
`ceil((b/a)*N_upper*(w+2*10^-12)/gain_lower)+1` integers. All 2,848 incidence
classes are checked; the weakest bound occurs at `(a,b,c)=(5,30,5)`.
If any K of T candidate blocks lie in an interval of mass at most `p`, anchor
its smallest value and union over the other K−1 values. This gives
`T*choose(T-1,K-1)*p^(K-1)`. Arbitrary block removal and a fitted constant gain
or interval centre are allowed. The calculation additionally charges
`105120^3` for selecting calendar endpoints and K inside 2021, and all
11,782,056 ordered distinct event-set pairs. It covers neither searches in
other years nor other block durations.

Using complete UTC-grid summaries from the [public Ar13+/Yb+ data](https://zenodo.org/records/6901524),
the direct regional bounds are `10^-121` for 36Ar and `10^-22` for 40Ar.
A separate conditional extension averages eight **independent native replicas**
through a preselected causal, nonnegative, unit-integral linear filter with
**compact support `[0,30 s]`**. A ticks are equally spaced in the declared
calendar coordinate. Interior `[0,270 s]` ticks contribute their full filter
mass to one 300-second block. Conditioning on every noninterior increment
leaves independent block innovations with arbitrary fixed outside shifts;
the same anchor-and-subset bound applies. Binary geometric digits and a
split Fourier integral give primitive-port bounds `10^-123`, `10^-38` and
regional bounds `10^-123`, `10^-25` for 36Ar, 40Ar respectively. The regional
worst class is `(29,30,29)`. Filter selection from the measured values is not covered.

The [experiment](https://arxiv.org/html/2205.13053v2) reports eight servo
outputs and a response time constant around 30 seconds. Those facts do **not**
establish independent native replicas, equally spaced native ticks, or a
compact impulse response: a decaying filter can retain an infinite tail.
The extension rejects its stated conjunction, not the actual atomic servo.

A separate **moment-matching diagnostic** equates the unfiltered, zero-error
population variance to the sample Allan value. It requires at least
`2544097785171`, `2407037601663` independent primitive replicas, or
`53002037192`, `50146616702` at the finite carrier's regional floor, for
36Ar, 40Ar respectively. These are conditional algebraic targets, not
confidence bounds on physical resources; the sample Allan value is not a
certified population upper bound. No continuum refinement level follows.
The refining source has `M_l=630*4^l` seams and `a_l=1260*16^l` attempts per
density-coordinate unit, hence a fixed degree-d port has **expected**
`a_l*d/M_l=2*d*4^l` incidents. Its seconds calibration is free, and the
thirty-seam regional floor does not transfer to arbitrary refined regions.
Identical regional event sets satisfy `B(N)=N` and are outside the test.

The clock and field comparison also reuses the existing
[pair-incidence reconstruction](../relative_observer_geometry/README.md#inferring-nonuniform-scheduling-from-clock-comparisons)
in a directly testable moment form. Let `I` mark all ports incident to one
source event, `p=E[I]`, and `Sigma=Cov(I)`. Fair repair gives
`P=Id-diag(p)+(Sigma+p*p^T)/2` for the mean density transfer. Every event is
present, including zero current, and ordered coincident tapes supply the
index. The optical module checks this identity on all 3,960 native branches;
no stationary load preparation or second field rate is needed. For a fixed
window of N IID events, estimate `Q=E[I I^T]`, reusing its diagonal for p.
Chebyshev on its 78 distinct entries gives operator error at most `6*delta`
with failure at most `78/(4*N*delta^2)`. This is a finite sampling certificate
for the existing clock-to-field correspondence, not an independent physical
derivation. Public ratio files lack these joint clock/occupation records.

The protocol uses event order: preparation, one repair, then readout. It
requires no conversion of a solver counter to physical seconds. A laboratory
implementation must nevertheless implement and calibrate those operations.
No spatial continuum, gravity or universal physical clock is inferred from
this finite test.

## Recovering the modular clock from visible current records

[inverse_clock.py](inverse_clock.py) ties a current record to the coherent
readouts of the [same three-mode transfer](../transfer_modular_probe/README.md).
Its input is the complete movement-bit tape `B_t=|j_t|`, **including every
zero-current attempt**, under the uniform IID seam schedule and fair coin on
the binary total-two carrier. This is an admitted complete-attempt record
service, stronger than an isolated port's incident tape. Record position
defines its clock; no physical seconds conversion enters the construction.

Write `p=E[B_t]` and `c_n=E[B_t B_(t+n+1)]-p^2`. Thus `c_0` is the
covariance of consecutive records, not the same-time Bernoulli variance.
The actual marked strong lumping gives
`c_n=w_+ lambda_+^n+w_- lambda_-^n` with positive weights and distinct
positive roots. Four consecutive covariances invert this relation:

```
D = c_0*c_2-c_1^2 = w_+*w_-*(lambda_+-lambda_-)^2
S = (c_0*c_3-c_1*c_2)/D = lambda_++lambda_-
P = (c_1*c_3-c_2^2)/D  = lambda_+*lambda_-
lambda_± = (S ± sqrt(S^2-4P))/2
w_+ = (c_1-lambda_-*c_0)/(lambda_+-lambda_-)
w_- = (lambda_+*c_0-c_1)/(lambda_+-lambda_-)
```

Full 66-state marked path sums, computed independently of the forward
spectral formula, give `p=5/33`, `D=1/9702990000`, `S=5/3` and `P=311/450`.
Consequently the observed two-mode domain determines
`lambda_±=(25±sqrt(3))/30` and `E_±=-log(lambda_±)`. Their unit-free ratio
is `E_-/E_+≈2.2033420301466635`. The one-exponential approximation matching
`c_0,c_1` misses `c_2` by exactly `1/2673000`; matching an average rate or
one decay does not establish the joint prediction. This approximation is
a moment control, not a claimed alternative native instrument.

`identify` accepts rational probability intervals and encloses every
compatible model in this positive two-mode domain. It verifies positive
Hankel determinant, separated real roots in `[1/2,1)`, and positive residues;
failure produces an inconclusive result. Repeated uncertain parameters are
shared in the underlying model; outward interval arithmetic enlarges their
image. A successful enclosure does not certify that an arbitrary input box
contains a realizable model, or rule out additional hidden modes outside the
declared domain. It also does not distinguish C from E, whose classical
current tapes agree.

`predict_readout` uses those same root intervals to predict the existing
coherent code's X and Y quadratures: `q_X=(1+cos(s E_j))/2` and
`q_Y=(1+sin(s E_j))/2`. Preparation in a specified spectral mode, the
coherent phase operation and the detector basis are declared services;
passive repair is not identified with this unitary. The retained example
uses `s=4`, both modes, and all four outcomes including padding. Synthetic
source-centred input tubes of radius `1e-12` predict the four zero-address
probabilities within approximately `[.946797,.948659]`, `[.721645,.723507]`,
`[.762113,.764252]` and `[.924060,.926199]`; the compiled source outputs lie
inside the corresponding outward bounds. The same inference predicts the
unused lag `c_4`. No separate field frequency or clock scale is fitted.
A common unit convention rescales both frequencies together. The raw arrays
retain padding without normalization and receive twice their complete
half-diamond error as an event-weight allowance.

The probability tubes have a dependent-tape coverage contract. Use the same
`M` starting positions from `M+5` retained bits to estimate the mean and five
lag products, reserving the fifth product for validation. A bounded function
of `L` consecutive marked events has variance at most `1/4`. For separated
windows whose starting positions differ by `k>=L`, conditioning on their
endpoint states and the reversible marked lumping gives absolute covariance
at most `(9/10)^(k-L)/4`, since its nonconstant eigenvalues are below `9/10`.
Overlapping windows retain the `1/4` bound. Summing covariances gives
`Var(mean_L) <= (2L-1+2/(1-9/10))/(4M)`. For the six statistics (`L=1..6`),
Chebyshev and a union bound give simultaneous tube failure at most
`39/(M*delta^2)`, without independent-window or independent-validation data
assumptions. For an arbitrary initial binary density, discarding a prefix
of 221 native attempts adds at most `1e-3` by the existing
[unread warm-up theorem](../record_laws/README.md); no retained controller
record is used to select the subsequent tape.

At `delta=1e-12`, this conservative proof gives total failure at most `.01`
for `M=4333333333333333333333333334` windows. This is a sufficient bound,
not a minimal data requirement, a collected dataset, or a feasible
experimental protocol. Broad tubes can be inconclusive. The retained
receipt exercises algebraic inference, interval coverage and complete
source readouts; the public optical data below do not supply this marked
current tape or the coherent-mode measurement dictionary. An atomic or
continuum-field identification must preserve these joint constraints before
those observations can test this correspondence.

Replay with `python -m model.predictions.inverse_clock --check
model/predictions/inverse_clock_results.json`.

### Preparing the finite edge reference from source pulses

The same replay implements an active preparation for the
[native edge generator](../observer/COHERENT_CONTINUATION.md#a-finite-edge-generator-representation-from-native-repair).
Its density transfer is `K12=I-L/60`; `L` has eigenvalues
`0, 5-sqrt(5), 6, 5+sqrt(5)` and sector ranks `1,3,5,3`. A closed-history
trace represents a Gibbs functional but does not prepare a coherent density.
Here a supplied protected register in `I12/12`, a declared event preparation,
and an explicit herald produce the density itself.

Prepare each of the thirty even event addresses with probability `(2-p)/60`
and each odd address with probability `p/60`. The existing compiled reversible
shift fixes an even address and swaps the corresponding seam endpoints at an
odd address. Unprepare the event register and herald a valid site address,
event zero and clean scratch.
The success operator is exactly
`F_p=(1-p)I+pK12=I-pL/60`. The retained instance uses `p=1/2` and two
successful heralds, producing `rho=F_p^4/Tr(F_p^4)`. The source routines are
executed pulse by pulse on all columns, then reused on the fifteen-qubit
register for all twelve logical input columns. Every first-failure output,
including address padding and dirty scratch, contributes to the total
probability. The complete raw filter is reused for the second herald.
Fixed numerical error guards test those arrays; they are not interval or
physical precision certificates. An executed source rotation
`exp(i*pi*Y_01/4)` followed by all sixteen address outcomes reads the
normalized zero-minus-one contrast `260343/43825510`; dephasing the prepared
density reduces this contrast to zero. Sector projectors below are
mathematical readout predictions, not a compiled sector detector.

For a fixed preparation parameter `t>=0`, set `p=30t/n<=1`. After `n`
successful heralds,

```
rho_n = (I-tL/(2n))^(2n) / Z_n  ->  sigma_t = exp(-tL)/Z
p_R(n,t) = d_R (1-t lambda_R/(2n))^(2n) / Z_n
Pr(all n heralds succeed) = Z_n/12 >= 1/12.
```

These relations tie preparation and sector readouts to the same generator,
with no independently fitted sector temperatures. Since `lambda_max<8`,
`0 <= -2n log(1-t lambda/(2n))-t lambda <= e_n`,
where `e_n=16t^2/[n(1-4t/n)]`. The finite weights are the Gibbs weights
multiplied by factors in `[exp(-e_n),1]`. Rejection coupling gives trace
distance at most `min(1,e_n)`; the normalized log ratios give
`D(rho_n||sigma_t)<=e_n`. The full-register modular gap inferred from
`log[(p_S/d_S)/(p_R/d_R)]` differs from `t(lambda_R-lambda_S)` by at most
`e_n`; it is fixed by the same preparation and sector law. Thus the Gibbs variational identity
`S(rho)-t Tr(rho L)=log Z-D(rho||sigma_t)` supplies an explicit
free-entropy deficit and stationarity at the limit on the trace-one fibre.
This is not a derivation of bulk-plus-area entropy stationarity.

The density has nonzero off-diagonal entries: dephasing preserves site
probabilities but changes the prepared reference. Its sector weights use
`d_R`, not the regular-group multiplicities `d_R^2`; the represented algebra
has dimension 44, with no four-dimensional irrep. Full group sewing and a
compact-group limit do not follow. The edge modular generator is `tL` up to
a scalar on the full logical algebra `B(C^12)`, which includes operators
between sectors. On the represented algebra `M1+M3+M5+M3`, the state is
central and its modular flow is the identity. The full-register generator
is distinct from the wave generator `sqrt(L)`; one clock conversion
cannot identify their gaps across all bands. The ideal growing filter has
the proved limit, while source precision and service budgets for that entire
family require a separate bound. Biased coherent preparation, protected
storage, initial density and herald readout are admitted observer services,
not a state selected by passive repair. Public clock-ratio datasets do not
contain this joint preparation and sector measurement, so these are
conditional model predictions, not an empirical discriminator against SM
or GR.

### Positive subordination on the prepared register

The [positive native-lag theorem](../physical_space_time/OCCUPATION_FIELD.md#the-poisson-extension-as-a-native-lag-record-readout)
and the [semigroup Poisson formula](https://arxiv.org/html/0910.2569v2#S1.E9)
give a constructive route from heat to square-root Gibbs weights:
`exp(-b sqrt(L)) = integral mu_b(t) exp(-tL) dt`, with
`mu_b(t)=b exp(-b²/(4t))/(2 sqrt(pi) t^(3/2))`.
Mixing normalized states `exp(-tL)/Z_t` with weights `mu_b(t)` is incorrect.
The conditioned mixture requires weights proportional to `mu_b(t) Z_t`.
Sampling a requested filter history first and retaining its actual success
probability supplies this factor without a fitted partition function.

The executable reuses the same source filter `F=I-hL/2`, where `h=p/30`.
For `a_k=binom(2k,k)/(4^k(2k-1))`, independently sampling
`J~Poisson(R)` and summing J jumps of law `a_k` gives

```
sum_n w_n F^(2n) = exp(-R sqrt(I-F²))
                 = exp[-b sqrt(L-hL²/4)],   b=R sqrt(h).
```

This is positive discrete subordination, with a controlled filter bias
relative to `exp(-b sqrt(L))`. The finite implementation caps jumps at B
and counts at J, replaces `exp(-R)` by a rational lower Taylor bound, and
rounds each coefficient down to a dyadic grid. A fixed word of fresh native
placement-coin records selects the requested n. Every unused interval is a
recorded sampling rejection; the coefficients are never renormalized.
The decoder checks actual native histories at every interval boundary.
The fair IID coin law, protected working register and classical controller
are explicit services.

Let `delta=1-sum w_n` and let s be a rational upper bound for `sqrt(8h)`.
The omitted positive operator is at most `delta I`. With
`q=max(w_0,1-Rs)` and `eta=Rsh/(1-h)`, the conditioned state has trace
distance at most `min(1,delta/max(q,1/12)+eta)` from the square-root Gibbs
target. If `delta<q`, constant-to-band modular gaps differ from
`b sqrt(lambda_j)` by at most `delta/(q-delta)+eta`.
For odd Taylor degree k, split R into `ceil(R)` arguments at most one.
The missing-mass budget is

```
delta <= R binom(2B,B)/4^B + R^(J+1)/(J+1)!
         + 3 ceil(R)/(k+1)! + (BJ+1) 2^(-sampling_bits).
```

Thus caps, Taylor precision and coin-word length can be chosen before the
readout. For fixed b, taking `h=1/m²`, `R=bm`, and vanishing truncation
budgets gives the square-root Gibbs limit. The untruncated jump law has
infinite mean. Even stopping at the first failed herald leaves infinite
expected work for the exact sampler, because the constant sector survives
every filter with initial weight `1/12`. Every capped trial has finite cost.
This is a convergent family, not a bounded-cost exact preparation.

The complete source instance uses `p=R=1/2`, B=J=2, 64 native coin events,
and at most four reused raw filters. It retains each sampling label,
first-failure event and all sixteen successful interference outcomes.
Its conservative Gibbs trace-distance bound is below `.217`; the finite
spectral diagnostic is about `.01344`. Incorrectly averaging normalized
component states changes the result by about `.000626` in trace distance.
A separate unexecuted B=64, J=8 contract bounds Gibbs trace distance by
`.046`, with at most 512 filters and 58,404,864 native filter pulses.
Sampler/control and readout costs are distinguished from these pulse counts.
Numerical source guards do not certify physical precision for the growing
family.

One prepared full-register state also fixes joint readout predictions.
Choose one normalized vector `v_j` in each positive eigenspace, set
`X_j=|0><v_j|+adjoint`, and let Y_j be its Pauli partner. If p_0,p_j are
the probabilities per vector, then
`S_j=<X_j²>=p_0+p_j`, `<[X_j,Y_j]>/(2i)=A_j=p_0-p_j`, and
`S_j/A_j=coth(g_j/2)` with `g_j=log(p_0/p_j)`. For the declared
square-root clock, the two-point law is
`S_j cos(sqrt(lambda_j)t)-i A_j sin(sqrt(lambda_j)t)`.
Independent copies normalized by `sqrt(2N A_j)` give the existing joint
CCR limit, with variance `coth(g_j/2)/2`; cross-mode covariances vanish.
This replaces separately assigned band polarizations by one state. The
small capped instance does not have an exact common inverse temperature.

Nor is `sqrt(L)` the exact finite native-walk frequency:
`Omega_30=sqrt(30) acos(I-L/60)`. The receipt separately bounds this
dispersion error; its band gaps converge only under the stated scaling.
The continuation below compiles the spectral-pair detector and square-root
comparison unitary on this logical register. Positive-phase encoding of the
native arc walk is separate. This is neither a prepared
Fock state nor a radial product state. The ideal densities are central on
the represented group algebra, where modular flow is trivial; source
approximation errors have separate guards. Nontrivial gaps concern the
full logical algebra. Passive selection, physical clocks and a public-data
identification do not follow.

<a id="prepared-reference-clock"></a>
### One prepared register, incident clock and feedback history

[prepared_clock.py](prepared_clock.py) continues that actual source preparation
through a source-compiled clock, a current-selected intervention and a
complete detector. The work register is never replaced by its ideal target.
The [retained replay](prepared_clock_results.json) compares every quantum
output, including an untouched twelve-dimensional reference, with an
independent endpoint/spectral calculation.

The canonical finite instance is:

| Component | Executed contract |
| --- | --- |
| Carrier | One binary total-two native carrier, initially occupied at ports 10 and 11; one protected four-qubit work register |
| Preparation | Supplied `I12/12`, represented with its retained purification; the existing source filter `F=I-L/120`; the existing 64-coin positive request sampler, capped at four filters |
| Records | Every sampler request/rejection, first herald failure, two complete native seam/coin/current records, and all sixteen terminal addresses |
| Clock | Every seam incident to port 10 triggers one fixed source word, including zero-current events |
| Feedback | The first committed current triggers a pair kick exactly when port 10 loses its load; the second repair starts from that repair's actual native state |
| Continuation | A source-bound checkpoint carries the unnormalized live arrays and binds the completed failure log; replay continues the serialized arrays |
| Refinement | This finite cap and two-attempt horizon are executed; the inherited sampler/filter limit has its separate approximation and cost contract |

All 64 sampling repairs use the empty seam `(0,1)`, so every coin word leaves
the carrier fixed. The two subsequent IID uniform-seam repairs therefore
continue that same carrier without a hidden request/state correlation.
All 3,600 ordered native tapes are enumerated. Six groups share work
operations; their full tapes and exact multiplicities remain available.
The horizon is two attempts, not a wait for two ticks: zero- and one-tick
outputs remain in the law.

The native endpoint graph fixes `L` and its four spectral bands. Exact
eigenvector equations and orthogonality establish the twelve-dimensional
basis; its compiled extension retains all sixteen addresses. The actual
finite preparation, rather than an asymptotic Gibbs replacement, fixes

```text
M = sum_(n=0)^4 w_n F^(2n),    rho = M/Tr(M),
m0 = sum_n w_n,               H = -log(M/m0),
U = exp(-iH).
```

The positive real logarithm fixes the gaps without a phase fit or integer
phase aliases. [prepared_clock_source.py](prepared_clock_source.py) emits
the basis changes, diagonal phases, intersector kicks and X/Y detectors
using the existing source compiler. Every routine runs on all of its
columns; execution reuses its raw returned matrix, retains padding and
performs no polar correction. The `sqrt(L)` comparison clock is also
compiled, but it is not the exact modular clock of this capped state.

For the ideal successful work marginal on `B(C12)`,
`rho^(is) O rho^(-is) = U(-s)^* O U(-s)`. Thus forward incident count `N`
corresponds to modular parameter `s=-N`. This follows by substituting
`log(rho)=-H+constant`, on all work operators, and composes with the
declared instruments by ordinary CP composition. It does **not** identify
the modular flow of the retained request/work algebra: the request-zero
branch is tracial, whereas the averaged work state has a nonzero gap.
The pure purification and four zero-weight padding modes likewise are
outside that faithful-state identification. The archive retains all these
records; taking its work marginal is an explicit access map.

The first noncommuting kick makes the clock observable. Let `m6` be the
eigenvalue of `M` on one chosen vector in the five-dimensional eigenvalue-six
band, and `g=log(m0/m6)`. The native event “first attempt loses the port-10
load and second attempt is incident” has probability `1/90`. On that event,
the unconditional X-address contrast between feedback on and off is

```text
(m0-m6) sin(g)/(12*90) = 9.09849736... * 10^-7.
```

The source execution reproduces it. Removing the clock or the
noncommuting kick removes this phase signal. The preparation success
probability is about `0.80215045`; it has not been divided out of the
contrast. Independent tests check the complete per-request X/Y laws,
native continuation, and corrupted reference phases that preserve every
outcome probability.

The complete error comparison uses **numerically checked operator guards**,
not interval arithmetic or measured device fidelities. For a one-filter
guard `delta`, put `g_filter=(1+delta)^2`. Retaining every first failure gives
preparation half-diamond allowance
`Eprep <= sum_n w_n (g_filter^n-1)/2`. If successful continuation has channel
gain `Gc` and half-diamond allowance `Ec`, telescoping gives
`Ecomplete <= Gc*Eprep + Psuccess*Ec`. Raw amplification is included; no
success conditioning enters this bound. A separate `1e-10` numerical
arithmetic allowance is included. The combined allowance is below `6e-8`,
and direct complete quantum-output discrepancies are below `1e-10`.
Each single-setting raw event weight can differ by `2*Ecomplete`; the
two-setting signed contrast remains positive after `4*Ecomplete`.
Conditioned claims require their separately reported success-probability
penalty. A future interval certificate would strengthen these numerical
checks without changing the finite experiment.

At most four source filters and two clock steps are invoked. The maximum
is **538,182 source pulses**, with one load and export per pulse; numerical
cache reuse removes no service invocation. Peak encoded storage is nineteen
qubits, including the retained reference and filter environment/scratch,
besides the native carrier and six-dimensional source processor. The receipt
also counts the serialized source-word packet, its byte-address cursor,
the separate expanded-pulse ordinal, records, counters, the five 64-bit
sampler thresholds and early-failure resource expectations. The source-word
packet is a storage subtotal; it does not include the controller/interpreter
program or working memory. The two observation-setting bits are counted separately.
At most 66 native-carrier placement coins are used: 64 for the sampler and
two for subsequent repairs. At most two IID uniform selections among thirty seams
are separate supplied service calls, with ideal mean `2*Psuccess`. Exact
uniform selection among thirty possibilities has no fixed finite fair-bit
implementation, so no such worst-case bit budget is claimed.
Preparation blanks, protected storage, coherent source services,
herald/readout and the controller's paused-carrier schedule are supplied.
Physical duration, drift and service fidelity are not bounded.

This completes the finite preparation-to-clock-and-detector connection.
It supplies neither a growing nonconstant product-field encoder nor the
original eight-mode scalar benchmark, and it does not identify empirical
seconds, interacting fields or a relativistic spacetime.

```sh
python -m model.predictions.prepared_clock --check model/predictions/prepared_clock_results.json
python -m pytest -q -W error model/predictions/test_prepared_clock.py model/predictions/test_prepared_clock_source.py
```

<a id="shared-state-fluctuation-field"></a>
### A field of fluctuations from the same prepared register

Write the twelve eigenvalues of the ideal successful work state as `p_j`, with
constant-vector weight `p_0`, and put `delta_j=p_0-p_j>0`. The **same state**
defines eleven pairs
`Q_j=X_(0j)/sqrt(2 delta_j)`, `P_j=Y_(0j)/sqrt(2 delta_j)`.
Their expected commutators are `i delta_jk`; their symmetric covariance is
diagonal, with each Q/P variance `coth(g_j/2)/2`. The existing clock rotates
each pair at `g_j=log(p_0/p_j)`. Thus one preparation fixes noise, ordered
response and clock rates together:
`<Q_j(t)Q_k>=delta_jk[coth(g_j/2)cos(g_j t)-i sin(g_j t)]/2`, where t counts
the ideal incident-triggered clock steps and its declared continuous extension.
The finite matrices do not themselves
obey exact canonical commutation relations. These observables belong to
the full logical algebra; the reference state is central on the smaller
represented group algebra.

For independent identically prepared successful copies, let
`F_N(A)=sum_(r=1)^N A^(r)/sqrt(N)`. If chronological factors have generators
`A_1,...,A_m`, their ordered characteristic function is exactly `z_N^N`,
where `z_N=Tr(rho exp(i A_m/sqrt(N))...exp(i A_1/sqrt(N)))`. Put

```text
C = sum_k Tr(rho A_k^2)/2 + sum_(l>k) Tr(rho A_l A_k),
D = sum_k ||A_k||,
|z_N^N-exp(-C)| <= (D^4/24 + |C|^2/2)/N.
```

The parity that fixes the constant vector and negates the excited vectors
commutes with `rho` and anticommutes with every chosen generator. Odd
Taylor terms vanish. The fourth derivative of the ordered unitary product
has norm at most `D^4`; its remainder is therefore at most `D^4/(24N^2)`.
Since `Re C=Tr(rho(sum A_k)^2)/2>=0`, the exponential remainder is at most
`|C|^2/(2N^2)`. Telescoping between two contractions gives the bound above.
This directly proves the ordered Gaussian/CCR limit without assigning
separate polarized states to the modes. Independent copies and their
preparation are premises.

The executed witness uses representatives of all three nonconstant bands,
one prepared work register, and a source-compiled pointer. All 32 joint
pointer/work addresses, twelve reference columns, requests, rejections and
first herald failures are retained. Its chronological generators are
`Q_3/10,P_3/10,Q_1/50,P_2/40`, with pulse-angle parameter `N=4096`.
Here subscripts label the three band representatives. The derived
`z_N^N` is approximately `0.76998918+0.00384689i`; the Gaussian limit is
`0.76997108+0.00384989i`. The finite-copy allowance is below `0.000944`;
including source numerical guards gives an allowance below `0.001462`.
The positive ordered phase is resolved. Raw estimator amplification
is charged when taking its Nth power. The request-zero branch supplies a
tracial control with zero Gaussian symplectic form; its finite ordered
phase need not vanish. A 4096-register circuit and access to all eleven
directions have not been executed. Pulse/load/export costs are counted for
the one-register witness; collecting successful copies by retrying has
mean preparation-attempt count `N/Psuccess` and no finite worst-case count.

<a id="growing-shared-state-field"></a>
The same construction has a quantitative growing-field contract. For
positive frequencies `omega_j`, use **field** variables
`q_j=Q_j/sqrt(omega_j)` and `pi_j=sqrt(omega_j)P_j`. With `g_j=b omega_j`
and model-time evolution `H/b`, their ordered covariance and response are

```text
W_b(t)=coth(b sqrt(D)/2) cos(t sqrt(D))/(2 sqrt(D))
       -i sin(t sqrt(D))/(2 sqrt(D)),
D_ret(t)=sin(t sqrt(D))/sqrt(D),       modular s=-t/b.
```

The finite prepared register verifies this normalization with
`omega_j=g_j`, `b=1`, its raw source clock and its actual successful
state. Using unscaled Q instead gives response derivative `omega_j`,
instead of one. These finite frequencies are not silently replaced by
`sqrt(L)`; the positive-log bound below controls that separate limit.

For a fixed menu of m real smeared factors, put
`R_k^2=sum_j(f_kj^2/omega_j+omega_j h_kj^2)` and `R=sum R_k`.
Each factor is a star matrix connecting one anchor to one excited vector,
with exact squared norm `sum_j[(f_kj^2/omega_j+omega_j h_kj^2)/(2 delta_j)]`.
All products and derivatives preserve the span of that anchor and the m
excited vectors, including their declared observation times. Its state
weight is at most `(m+1)p_0`. Thus the fourth-order Taylor remainder above
improves to `tau (sum_k ||A_k||)^4/(24N^2)`, where `tau` is this span's weight. If
`g_j,b omega_j>=gamma>0`, let `theta=gamma/(1+gamma)` and `Z=1/p_0<=d`:

```text
finite-copy error <= R^4 [(m+1)Z/96+1/8]/(theta^2 N).
```

Consequently `N/Z -> infinity` suffices on a positive fixed band; `N/d ->
infinity` is a conservative choice. The executable counterfamily
`rho=diag(2,1,...,1)/(d+1)` proves that choosing `N=d` does not suffice:
its unit-smear characteristic tends to `exp(3[cos(1/sqrt(2))-1])`, not
the Gaussian `exp(-3/4)`. This family tests the bound, not a native encoder.

The [executable contract](prepared_clock.py) also charges preparation,
clock, spectral and cutoff errors. Complete preparation half-trace error
E and success floor p give conditional trace error at most `4E/p` when
`2E<p`. A one-copy word error epsilon gives `zeta=4E/p+epsilon`; raw
Nth-power amplification is at most `N zeta/(1-N zeta)` when `N zeta<1`.
No bound is reported outside those domains. For positive filter mass m0,
`p=m0 Z/d>=m0/d`, so vanishing E alone is insufficient: a conservative
requirement is `N d E -> 0`. A continuum partition-function bound cannot
replace a uniform **finite-source** bound without its own proof. The native
family has such a proof from its existing ancestry and internal-gap bounds.

<a id="finite-native-trace-and-gap"></a>
Put `s=4^l`, `D_l=s L_l`, and write any full-port vector as `v=Uu+q`, where
q has zero mean in each carrier. For `E=<v,D_l v>`, Schur minimization,
[K1 and the internal gap](../physical_space_time/OCCUPATION_FIELD.md#a-common-positive-spatial-generator-from-the-native-green-field)
give `||(I-P_k)u||^2<=624*4^(-k)E` and `||q||^2<=E/(2s)`.
The centered coarse range has rank `20*4^k-1`, so min-max gives
`lambda_(20*4^k)(D_l)>=4^k/(624+1/2)` for `k<l`. At `k=l` only q
remains and the bound is `lambda_(20s)>=2s`. The full-source gap
`lambda_1>=s/(4096s+1)>=1/4097` handles the first nineteen eigenvalues.
For intermediate indices choose `20*4^k<=n<80*4^k`; high indices satisfy
`n<240s`. Thus, **at every finite level**,

```text
lambda_n(D_l)>=n/81940,
Tr exp(-b sqrt(D_l))<=1+163880/b^2.
```

The trace bound follows by integrating the decreasing function
`exp(-b sqrt(x/81940))`; the exact constant mode contributes one.
For the native density step `K=I-L_l/(2M_l)`, `M_l=630s`, a lazy event
filter equals `F=(1-p)I+pK=I-hD_l/2` precisely when `p=630hs^2<=1`.
This implies `h||D_l||<=12hs<2`. Retained positive subordination therefore
gives `M<=exp[-b sqrt(D_l-hD_l^2/4)]<=exp[-b sqrt(D_l/2)]`, and

```text
Z=Tr(M)/m0 <= (1+327760/b^2)/m0,
m0/d <= Psuccess <= (1+327760/b^2)/d,     d=240s.
```

This derives a uniform partition bound without removing the herald cost.
It also gives prepared gaps `g>=b/sqrt(8194)+log(m0)`. With omitted mass
`delta=1-m0` satisfying `delta/(1-delta)<=b/182`, both prepared and target
positive gaps are at least `b/182`, since `sqrt(8194)<91`. Raw preparation
and clock errors retain their separate allowances. The independent controls
check the actual 240- and 960-port spectra and the min-max index boundaries.
Contraction alone is insufficient: `D=diag(0,4/h)` makes `F=diag(1,-1)`,
so every even filter power preserves a high-energy mode. This endpoint
violates the native lazy-event constraint.

On a common embedded band, a gap error e contributes at most
`R^2[e/gamma^2+2T e/(b theta)]/2` for model times `|t|<=T`.
A normalized two-point spectral error kappa adds `R^2 kappa/2`.
Target-field phase-space tails of total radius r add `(R r+r^2/2)/theta`,
using the same positive target gap floor. Finite-source observables are
band-filtered; no bound on uncontrolled high source modes is presumed.
Event-count interpolation retains its separate native-clock bound; a
model-time phase bound alone does not establish a clock-reading error.
These constants permit a diagonal limit using the existing common
[A/B operator limits](../physical_space_time/OCCUPATION_FIELD.md#thermal-product-covariance-and-the-shared-wave-operator).
One first proves each bounded-band convergence, then removes the cutoff
under uniform tails `||D^(-1/4)f_tail||+||D^(1/4)h_tail|| -> 0`.
Smooth local smears require that tail justification; finite-band agreement
alone is insufficient. The sample full-port schedule `d=240*4^l`, `N=d^2`
has explicit decreasing error requirements, not measured compiler rates.
Growing state preparation, mode access and independent copies retain their
implementation costs. The source precision schedule must be chosen
jointly with the copy schedule: a dimension-independent error of order `4^(-r)` needs
`r-3l -> infinity` under the conservative full-port herald floor and
`N=d^2`. Keeping the same refinement index for both does not meet that
requirement; the stronger precision and sampler caps carry their full cost.
The native contract below uses the proved partition bound and charges the
complete preparation error to choose a different compatible copy schedule.

For mean-zero angular A the native constant vector supplies the anchor.
The Dirichlet product B is strictly positive: using it requires a separate
zero-energy reference sector `0 direct-sum B`, distinct from the physical
constant angular modes with positive radial energy. This additional sector
has not been constructed by the angular source fixture. No stationary
massless zero-frequency oscillator, physical radial transport or observed
Lorentzian geometry follows from this field-limit contract.

<a id="prepared-positive-log-limit"></a>
The clock also has a conditional limit tied to its preparation. For
`F=I-hD/2`, `M=sum w_n F^(2n)`, `m0=sum w_n`, put `H=log(m0)I-log(M)`.
Assume F is a contraction on the **whole encoded domain** and `h Lambda<2`
on a fixed band `0<=D<=Lambda`. With `b=R sqrt(h)`, missing sampler mass
`delta`, and `0<a<=exp(-b sqrt(Lambda))`, `delta<a` gives

```text
||H-b sqrt(D)|| <= delta/(a-delta)
  + b h Lambda^(3/2)/(4(1+sqrt(1-h Lambda/4))).
```

Positivity bounds the omitted sampler operator before taking the logarithm;
both scalar log corrections have the same sign, so normalization adds no
factor two. If each raw filter has error at most `epsilon`, a cap of K
filters adds mixture error `xi=2K epsilon/(1-2K epsilon)` and log error
`xi/(a-delta-xi)`, provided the denominators are positive. Insufficient
precision or retained mass produces no bound. This raw-source log bound
also requires the entire encoded domain to lie in the certified band, or
the source implementation to preserve that band: small leakage toward
arbitrarily small outside eigenmasses defeats a low-band log bound. The
finite twelve-dimensional heralded preparation satisfies this requirement;
the reported growing-schedule examples declare it as a compiler obligation.
Keeping caps fixed while `Kh` tends to zero instead collapses H to zero.
Combining this bound with a compatible spatial operator limit requires
bounded-band convergence followed by spectral cutoff removal; it supplies
no unbounded-generator norm limit or physical spacetime identification.

<a id="native-growing-preparation-and-clock-contract"></a>
The full native construction uses the sampler's zero-request branch to
control the whole spectrum, including leakage outside the observation band.
Set `k=2l+c`, with `4^c>=b^2 max(630,Lambda)`, `R=2^k`, `h=b^2 4^(-k)`.
This makes the native event weight admissible. With jump cap `B=2^(6k)`,
count cap `J=6R+2k`, `K=BJ`, Taylor order `4k+3`, and sampler bits at least
`max(48,ceil(log2(K+1))+2k,2R+2)`, omitted mass is at most `4a`,
where `a=4^(-k)`. The unchanged zero-request weight satisfies
`w0>=3^(-R)-2^(-sampling_bits)>=2^(-2R)=f`; no identity mass is added.

The [source resource bounds](prepared_clock_source.py) count the uniform
Bell input, native event filter, full spectral encoder and phase words.
For `P=max(Pinput,Pfilter)`, require primitive plus load/export error at
most `f a/(64dKP)` in both operator norm and complete-instrument half-diamond
norm. Composition gives complete preparation half-trace error
`E<=f a/(32d)`, including every stopped failure. For the unnormalized
successful marginal `Sraw`, `Mraw=d Sraw` obeys
`||Mraw-M||<=2dE<=f a/16`; hence the global logarithm error is at most
`a/(16-a)`. Normalization adds only a scalar to its logarithm. This
discharges the reducing-band condition for this quantified construction.

The ordered detector has a separate cost: complex star rotations are
synthesized by phase-adjusted real Householders and a controlled pair
rotation. Its pulse total `Pdet` includes the common encoder, pointer
preparation/readout and both clock words around every timed factor.
Requiring `2k+ceil(log2(8Pdet))` bits of primitive service accuracy bounds
the combined X/Y characteristic error by `4a/(8-a)<=a`, including raw
composition gain. `native_growing_field_contract` joins these costs with
the derived partition, positive-gap and success bounds on the **same**
finite `D_l`. For fixed b and band, `N=2^l` gives finite-copy error
`O(2^(-l))`, preparation amplification `O(Na)` and detector error `O(Na)`;
the gap/clock error also tends to zero. The conservative constants do not
claim useful accuracy at the first few levels. The finite resource formulas
are evaluated; increasing words, independent copies and physical service
accuracy are obligations. Continuum embeddings/tails, radial architecture,
interactions and laboratory units are not supplied by this contract.

The fixed finite gaps also admit a narrow public-data test. All nine ordered
assignments of the three executed constant-to-band gaps to the
[published Ar13+/Yb+ ratios](https://arxiv.org/pdf/2205.13053v2#page=20) are
incompatible, even allowing absolute attachment error `10^-6` in each gap
and `10^-12` in the ratio. Exact rational log enclosures leave a ratio
separation above `0.057` for each isotope. One common frequency scale
cancels; no preparation parameter, offset, harmonic, phase alias or
species-dependent scale is fitted. The error budgets are explicit tested
hypotheses, not hard bounds inferred from quoted standard uncertainties or
numerical compiler checks. This retrospective exclusion concerns this
direct assignment only; it supplies no atomic model or exclusion of OPH.

## Native repair response from passive clock and field records

[repair_response.py](repair_response.py) predicts a finite response using the
same repair law, fields and incident clocks that define the passive process.
There is no fitted force, response coefficient or conversion to seconds.
For the equilibrium numerical forecast, use the uniform fixed-number binary
reference law; finite native preparation has the error budget below.
Disclose the occupations at
ports 0, 1 and 2. If the port-0 bit is one, apply one fair native repair on
seam `(0,1)`; otherwise apply no repair. The trigger, coin and signed current
are retained. Postpulse fields are reconstructed from the disclosed bits and
that current, without another state query. The pulse contributes one tick to
each incident observer clock when it occurs, including a zero-current pulse.

Write `a=n_0`, `b=n_1` for the initial occupations in a *passive continuation*
experiment. The corresponding post-intervention law has density

```text
w = 1 + (b-a)/2,
Pr_feedback(F) - Pr_passive(F) = E_passive[(b-a) 1_F]/2.
```

Here `F` is any fixed subsequent record event, with the same continuation
rule in both experiments, independent of the earlier pulse mark. It may
include fields, incident-clock comparisons, feedback on subsequent records,
stopping and timeouts. A continuation that consults the pulse mark must use
the separate marked weights below, not the summed density. The identity follows from detailed
balance of the local fair repair: for its self-adjoint operator `R`,
`R a=(a+b)/2`, and the pulse policy has operator `I + a(R-I)`. This is the
finite Markov adjoint identity applied to the actual native intervention,
not a claim of a new general fluctuation-response theorem. It is exact at
the specified pulse amplitude; it is not an infinitesimal-response fit.

The complete marked prediction uses these passive-record weights:

| Retained trigger, coin, current | Weight |
| --- | --- |
| No pulse (`0,-1,0`, with `-1` the unused-coin label) | `1-a` |
| Pulse (`1,0,1`) | `(1-a)*b/2` |
| Pulse (`1,1,0`) | `a/2` |
| Pulse (`1,0,0`) | `a*b/2` |
| Each other declared pulse mark | `0` |

These are importance densities for unnormalized joint events, not
conditional probabilities of drawing a mark from each passive sample.
Their sum need not equal one at each sample. Their equilibrium averages
produce a normalized complete law. No division by that sum is permitted.
Pre-pulse disclosed fields can be recovered from the postpulse fields and
retained current, so the prediction preserves their record as well.

The executed forecast stops at the first subsequent port-0 incident event,
with a supplied watchdog cap of one or two uniform native attempts. The
second observer counts incidents at port 1 or port 2. Every signed-current
record meeting the three read ports, zero-current events, pulse marks and
timeouts is included. The exact timeout probability is `(5/6)^cap`.
For the event "the stop succeeds, the final port-1 field is one, and the
second clock has at least one *subsequent* tick", the cap-two predictions are:

| Second clock | Passive probability | Feedback probability | Response from passive covariance |
| --- | --- | --- | --- |
| Port 1 | `71/5400` | `34/2475` | `7/11880` |
| Port 2 | `1/216` | `5/792` | `1/594` |

The cap-one response is exactly zero for both clocks. Total clock readouts
also retain the pulse tick: subtracting that known contribution in the
specified event merely selects equal continuation menus for comparison.
The watchdog attempt count supplies no physical duration. Independent
primitive-path expansion checks all one-attempt record labels; exact native
propagation checks the capped laws. An unconditional extra repair has zero
marginal response at equilibrium. A nonstationary preparation violates the
stationary weight formula in general; the retained control has predicted total mass
`1/2` against actual mass one and is rejected without renormalization.

**Local preparation removes the equilibrium premise from the response identity.**
Let `mu` be any binary fixed-total-two population law, and let `S` exchange
the pulse endpoints. Only the moving mark can differ between the actual
pulse law `A_mu` and the passive weighted measure `W_mu`:

```text
(A_mu-W_mu)((1,0,1),s) = [mu(Ss)-mu(s)]/2,  (s_0,s_1)=(0,1),
L1(A_mu-W_mu) = TV(mu,S mu)/2 <= TV(mu,uniform).
```

All other marked differences vanish. Thus endpoint-swap invariance is
necessary and sufficient for the complete marked identity; global
stationarity is unnecessary. One **unconditional fair native repair of the
pulse seam** prepares `(mu+S mu)/2` exactly. It preserves the rest of the
configuration and costs one incident tick at ports 0 and 1 in both arms.
Its coin/current record must be unavailable to the later protocol. The
receipt executes this preparation and retains a nonstationary example at
TV distance `32/33` from equilibrium whose pulse changes the port-1 mean
from `1/2` to `3/4`, exactly as its passive records predict. This derives
the local preparation; the pulse controller and field access are the declared observation protocol.

Without that local preparation, any common complete classical continuation
contracts the displayed L1 error, independently of its finite stopping cap.
This includes every pulse mark, field, zero-current incident tick and timeout.
`W_mu` can have mass different from one: its L1 distance must **not** be
halved or silently renormalized. A retained perturbation of size `1/132`
attains both the L1 and event bound; its two-attempt timeout discrepancy is
`25/4752`. Equal endpoint means and total forecast mass one alone do not
ensure the identity, as a second exact control shows.

For numerical forecasts relative to equilibrium, the
[native warm-up theorem](../record_laws/README.md#preparation-by-finite-native-warm-up)
provides `delta_n=min(1,33*r_alpha^n)`, with
`r_alpha=1-alpha*(1-19079/20000)`. The actual complete readout probability
law is within TV `delta_n` of its equilibrium law; the passive weighted
forecast using the **same warmed preparation** is within L1 `delta_n` of
the actual pulse readout. Uniform scheduling gives both errors below
`1/1000` after 221 repair attempts. The binary-sector support, exogenous
seam coverage, admitted record law and discarded warm-up tape are inherited
hypotheses. These are unconditional bounds; conditioning on a selected
preparation outcome or rare later success requires a separate bound.
The response identity uses classical occupation/current records; coherent
readout of the repair environment requires a separate calculation.
The executable wrapper does not extend this binary readout to collisions,
identify a physical duration, or infer a preferred quantum record law.

For an arbitrary finite graph with `V` ports, `M` seams, fixed total `k`,
`chi=k*(V-k)/(V*(V-1))`, and native density transfer `K=I-L/(2M)`, the same
identity gives the complete mean occupation response to a pulse on `(a,b)`:

```text
Delta <n_c(n)> = chi/2 * [(K^n)_(b,c) - (K^n)_(a,c)].
```

For any chosen spatial normalization `s>0`, let `r=2*M*s/(2*M*s+z)`, `z>0`.
Then `(1/(2*M*s))*sum_(n>=0) r^(n+1)*Delta<n(n)>` equals
`chi/2*(e_b-e_a)^T*(s*L+z*I)^-1`. Full 66-configuration inverses independently
verify this response on the native carrier. This makes the
[spatial Green operator](../observer_time_fields/REFINEMENT.md#native-temporal-correlations-and-the-spatial-green-field)
a response to an observer's explicit repair operation, as well as a passive
correlation. One microscopic pulse does not by itself have a nonzero
continuum normalization. The following controlled source supplies one.

### A controlled continuum mean with an executed-event clock

Fix a real port profile `h`, disclose the initial binary load and reconstruct
subsequent occupations from retained currents. In each selection slot choose
a uniform seam `(a,b)` and execute its unchanged fair native repair with
probability

```text
p_execute = [1 + epsilon*(h_b-h_a)*(n_a-n_b)]/2.
```

Every selection, idle/execution flag, repair coin and current is retained.
The source certifies `|epsilon*(h_b-h_a)|<=1/2`. The designated dynamical
clock counts executed repairs, including zero-current repairs. It is an
aggregate event clock; it is not a primitive port-incident clock. Selection
slots and controller service have separate counts. Native fair bits can be
read by repairing a disclosed equal-occupation seam: its load stays fixed,
but every service write is a charged native event. An equal pair exists
in every binary carrier triangle. These resource counts do not identify
seconds, an atomic transition or an autonomous force law.

For the uniform fixed-number preparation the execution subkernel has
stationary incoming mass `pi/2`, even at nonzero epsilon; its row activity
depends on the load. Dividing by that activity gives the actual kernel at
successive executed repairs. Differentiating this normalization at zero
produces the extra `K` in the exact response:

```text
executed repair count r: d< n >/d epsilon|_0 = 2*chi*K*(I-K^r)*h,
selection slot count q: d< n >/d epsilon|_0 = 2*chi*(I-Q^q)*h,
Q = (I+K)/2.
```

Full 66-configuration propagation checks both identities. Replacing `q`
by `2r` is not exact: the retained finite controls have different responses.
This correction retains the correlation between clock stopping and load.

On the actual level-`l` carrier set `s=4^l`, `N=240s`, `M=630s` and density
one third, with its canonical susceptibility `chi`. For a fixed mean-zero
flat coordinate `f`, put `h_l=cell_average(f)/sqrt(N*chi)`, copied to the
ports of each cell, and read
`eta_l(g)=sum_p cell_average(g)_p*(n_p-1/3)/sqrt(N*chi)`.
The [common native generator](../observer_time_fields/REFINEMENT.md#the-full-native-dynamics-have-the-same-limit)
then gives, along the same spatial subsequence, at
`r=floor(2*M*s*tau)=floor(1260*16^l*tau)` executed repairs,

```text
d E[eta_l(g)]/d epsilon|_0 --> 2 <g,(I-exp(-tau*A))f>.
```

The small-`tau` source is `2Af` in the dual of the closed form domain when
`f` lies in that domain. This derives a controlled mean response through
the same `A` as the passive field and heat limit. It does not identify a
physical force or quantum field interaction.

The statement also has a finite, level-independent amplitude. At `tau=1`
stop at `r` executed repairs or `q=4r` selection slots, retaining the actual
terminal field on timeout. Against the zero-bias law the stopped likelihood
is `L=product(1+epsilon*s_j)`, where `s_j` is the signed activation score
and is zero after stopping. Its conditional mean is zero and variance is
at most `v=h_l^T L_graph h_l/M`. The
[coordinate energy bound](../physical_space_time/OCCUPATION_FIELD.md#a-common-positive-spatial-generator-from-the-native-green-field)
gives `q*v<=8/chi<36`. Writing `S=sum s_j`, martingale conditioning yields

```text
E L^2 <= (1+epsilon^2*v)^q <= exp(epsilon^2*q*v),
E (L-1-epsilon*S)^2 <= (1+epsilon^2*v)^q-1-q*epsilon^2*v.
```

Indeed the remainder increments are `epsilon*s_j*(L_previous-1)`;
conditional centering removes their cross terms. When `epsilon^2*q*v<=1/2`,
Cauchy--Schwarz bounds the nonlinear mean error by
`36*epsilon^2*||g||_L2(mu)`, uniformly in refinement.

The cap comparison retains the stopping correlations. At zero bias the
uncapped slot count `T_r` has mean `2r` and
`P(T_r>4r)<=2^(-r/2)`. With `B` the maximum absolute terminal field, the
squared difference of capped and uncapped derivatives is at most
`8*B^2*r*v*P(T_r>4r)`: use the complete stopped score, whose second moment
is at most `2r*v`, and the fact that the fields differ only on timeout.
The biased timeout obeys
`P_epsilon(timeout)^2<=E_0[L^2]*P_0(timeout)`; it is not assigned the
unbiased watchdog probability. The executable propagates all likelihood,
score and field moments on the complete finite state space and verifies
these correlations, including a nonzero nonlinear mean correction.

For a flat coordinate of the unit carrier, `||f||<=1` and
`||f||^2>=3/16`. Since `A>=3/4096`,
`2<f,(I-exp(-A))f> >= 9/32792`. Consequently the fixed choice
`epsilon=10^-6` has the strictly positive bound

```text
liminf_l E_epsilon[eta_l(f)] >= epsilon*9/32792 - 36*epsilon^2
                            = 244359/1024750000000000 > 0.
```

This is a nonvanishing controlled continuum signal. It does not prove a
unique nonlinear law at fixed epsilon. A finite fair-bit implementation
rounds each activation probability down using
`b=ceil(log2(16*N*q^2))` bits. Coupling bounds the path discrepancy by
`q*2^-b`; the squared terminal-mean error is at most
`1/(64*q^2*N*chi)`, which vanishes. At most `b*q` native coin-service
writes are charged, plus the separately supplied uniform-seam selection
service. This implementation is a declared observer protocol, not a
claim that passive repair spontaneously chooses this drive or clock.

The [published Ar13+/Yb+ files](https://zenodo.org/records/6901524) contain
corrected optical frequency ratios. They do not contain these disclosed
occupations, trigger/current records and corresponding intervention arms.
Their [experimental description](https://arxiv.org/html/2205.13053v2) supplies
an atomic apparatus; this packet supplies no identification of its readouts
with the native pulse. These exact native forecasts therefore add no public
measurement agreement or discriminator against SM + GR. Choosing event-count
units is sufficient to define the model clocks; identifying an observed
physical response must preserve the joint relations above.

Replay with `python -m model.predictions.repair_response --check
model/predictions/repair_response_results.json`.

## Minimal modular-transition optical attachment

The source-derived [three-mode transfer generator](../transfer_modular_probe/README.md)
has energies `0`, `E_+=-log((25+sqrt(3))/30)` and
`E_-=-log((25-sqrt(3))/30)`. Its positive gaps are `E_+`, `E_-` and
`E_--E_+`. [modular_optical.py](modular_optical.py) tests the explicit
hypothesis that the measured Ar13+ and Yb+ optical transitions correspond
directly to two such gaps with **one common positive frequency scale**.
That scale cancels in their ratio. A common Gibbs inverse temperature also
cancels, since `-log(T3^beta/Tr(T3^beta))=beta H+constant`. The test uses the
real logarithm of the positive transfer; it does not add integer harmonics,
alternative logarithm branches, species-dependent gains or offsets.

All nine ordered assignments, including the three same-gap assignments,
are evaluated by exact outward intervals. Their distinct ratios are
approximately `0.453856`, `0.546144`, `0.831019`, `1`, `1.203342`, `1.831019`
and `2.203342`. King et al. report the dimensionless ratios
`R40=1.05776938758748094(11)` and `R36=1.05776646273518748(13)` in
[Extended Data Table 1 of the accepted manuscript](https://arxiv.org/pdf/2205.13053v2#page=20),
published in [Nature 611, 43–47 (2022)](https://doi.org/10.1038/s41586-022-05245-4).
These table values and their quoted standard uncertainties are retained in
[the small data record](data/modular_optical_ratios.json); this comparison
does not reanalyse the measured time traces.

Under the declared **total absolute ratio error at most `1e-12`**, every
assignment is disjoint from either observation. The closest ratio is one;
the certified absolute separation exceeds `0.057` for each isotope. The
error budget combines measurement and attachment errors and is part of
the tested hypothesis. Quoted standard uncertainties are not converted
into hard bounds or a significance claim, and the two observations are
compared separately without an independence assumption. A `0.06` total
ratio-error control does permit an assignment. Allowing separate arbitrary
frequency scales would also remove the test's constraint.

This is a retrospective failure of the minimal direct-transition assignment.
The source does not derive atomic structure, identify these ions with the
three-state code, or supply an optical detector mapping. It therefore neither
excludes all modular clocks nor supplies a new discriminator against SM + GR.
It records a concrete restriction on a proposed physical dictionary while
keeping the source spectrum and its observer readout independent of the
measured ratios.

## Inherited coupling comparison

[shared_parameter.py](shared_parameter.py) propagates that same inherited
structured-running branch's forward P enclosure through the conditional
`n_s=1-P/48` map before reading observations. Its
[receipt](shared_parameter_results.json) gives `n_s≈0.9660214146696063` and
replays all 379,175 retained ACT chain rows: both interval endpoints have
tail weight `3114/508985`, unchanged from the alpha-calibrated branch.
The coupling endpoint misses by `300.388218` ppm. This tests one
shared parameter across those two comparisons without selecting new
parameters or completing their physical maps. It supplies no joint
likelihood or physical theory-error bound. Replay with
`python -m model.predictions.shared_parameter --check`.

[endpoint_audit.py](endpoint_audit.py) compares the pinned
`thomson_structured_running` output used by the
[P-dependent observer](../observer/p_response/README.md) with the
[CODATA 2022 value published by NIST](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
`alpha^-1=137.035999177(21)`. The inherited output is approximately
`136.994835177413`, a deficit of **300.388218 parts per million** relative
to that central value. Exact rational arithmetic makes its narrow numerical
interval disjoint from the central value plus or minus five reported
standard uncertainties.

The [P-derivation guide](../../code/P_derivation/README.md) states the
incomplete physical transport obligations. The numerical interval does not
bound those missing physical contributions. The audit rejects the additional
hypothesis that this output is the complete physical endpoint with zero
missing correction. It assigns no theory sigma significance, performs no fit,
and does not recalibrate P from the target. This is a retrospective endpoint
gap, not a refutation of OPH as a whole or a prospective prediction success.
