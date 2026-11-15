# A native current counter and its measurement laws

The same executable binary observer predicts its current-count and
waiting-record distributions. No independent random clock or Hamiltonian is
added. Its 66 native configurations reduce exactly to a three-state process
for **every finite absolute-current record word**. This reduction makes longer
record predictions inexpensive and identifies which scheduler assumptions can
be removed.

The counter is `N(n) = sum(abs(j_k), k=1..n)`, where `j_k` is the
committed native current. Its readings are observable integer event counts.
Here `n` always means **native attempt count**, including the recorded zero
currents. A physical duration, proper time or relation to an experimental
clock has not been identified.

## Local field readout without an initial-bit disclosure

[`bit_readout.py`](bit_readout.py) reconstructs a port's binary occupation
from its **own signed incident currents**, including zero-current ticks.
Before the first nonzero current its bit is unknown. A local `+1` sets the
outgoing bit to one; `-1` sets it to zero. Subsequent currents update that
bit exactly. The initial bit can then be recovered retrospectively. This
supplies no earlier causal access to it; unsigned movement ticks alone do
not identify the bit.

The source audit checks the number-projector input/output supports of every
C and E Kraus branch. Zero-current and remote branches preserve local number;
a signed branch fixes its incoming and outgoing number. The identities compose
on arbitrary quantum inputs and untouched references. For coherent inputs,
the recovered initial label describes the resolved instrument domain, without
assigning a hidden classical bit before measurement.

An exact twelve-class marked reduction of all 66 source configurations
preserves every local signed-current word. With `A` the remote kernel and
`B_j` the incident kernels, `Q_j=(I-A)^-1 B_j` sums all remote waits. The
unresolved kernel satisfies `Q_0^16 1 <= (3/8)1`, giving

```
P(unresolved after m local incident ticks) <= (3/8)^floor(m/16).
```

Thus 128 local ticks give failure at most `6561/16777216 < 1/1000`, uniformly
over initial states. The exact worst-case expected count is `129960/8261`
local ticks. All twelve ports obey the same certificate. A union bound gives
failure below `1/1000` after 160 ticks at each port, without assuming independent
observers. Joint record access and a common comparison checkpoint are separate
services; the individual stopping times are not simultaneous measurements.

Every all-zero finite prefix retains an unresolved outcome. Conditioning on
resolution changes field statistics: among immediate resolutions the initial
bit probability is `1/2`, although its stationary prior is `1/6`. Keep this
outcome selection in any subsequent field/tensor comparison. If a complete
readout window has resolved initial-one mass `r_1` and unknown mass `u`, its
initial field mean lies in `[r_1,r_1+u]`. This uses all outcomes; finite-sample
confidence limits require a separate sampling contract. The schedule is
independent uniform seam selection with fair native placements; feedback or
coherent interventions require another waiting law. No stationary preparation,
global attempt counter or physical duration enters this readout certificate.

```sh
python -m model.native_clock.bit_readout --check
python -m pytest -q -W error model/native_clock/test_bit_readout.py
```

### Finite confidence from one continuing local history

[`trajectory_readout.py`](trajectory_readout.py) supplies the sampling contract
for stationary local moments using one continuing signed-current tape. No
stationary preparation, resets, independent replicas or initial-bit disclosure
are needed. The horizon is fixed before reading the tape, in **own incident
ticks**. Under the same IID uniform seam schedule and fair native placements,
the stationary occupation, consecutive-occupation product and movement means
are respectively `1/6`, `1/11` and `5/33`. The product differs from `1/36`:
consecutive readings are correlated.

For the exact stopped kernel `Q`, a transition reward `h` has conditional
mean `r`, stationary mean `mu`, and a checked rational solution
`(I-Q)g=r-mu`. On each source path,
`sum(h-mu)=g(X_0)-g(X_n)+sum D_k`, where `D_k` is conditionally centered.
Let `R=osc(g)` and let `W` bound the conditional range of
`h(X_k,X_(k+1))+g(X_(k+1))`. Iterating the bounded-variable exponential
inequality gives
`P(|mean(h)-mu|>(R+W*t)/n) <= 2 exp(-2*t^2/n)`.
The program chooses the integer `t=ceil(sqrt(n*(s+1)/2))`, so each failure
probability is at most `2^-s`, with no floating-point confidence calculation.
For occupation, `R=W=3`. All twelve ports have the same exact constants;
their observations need not be independent.

The first nonzero current retrospectively resolves every local bit in the
window. An all-zero window returns `[0,1]` for the two bit moments and
reports its exact zero movement count. Coverage includes these unknown
outcomes and applies from arbitrary binary total-two inputs, including
coherent states with untouched references. It is not conditional on resolution,
an arbitrary stopping-time guarantee, or a confidence interval for an initial
nonstationary occupation. Three simultaneous intervals have total failure
at most `3/4096` at the default setting. A retained 32,768-tick nonstationary
native execution satisfies all three predictions; repeated repair of one
empty incident seam fails the movement test. This control changes the
scheduler, showing why its hypothesis remains necessary. No physical dataset
is identified with this local tape.

```sh
python -m model.native_clock.trajectory_readout --check
```

### Preparing a regional clock with local incident counts

[`regional_preparation.py`](regional_preparation.py) prepares a faithful
reference from arbitrary binary total-two inputs using the four-port cap's
incident ticks, including zeros. If `p` is its incident seam fraction and
`r=19079/20000` is the native mixing bound, the stopped quantum channel has
reference-uniform half trace error at most `33*p^(-1/2)*r^n` after `n` local
ticks. A similarity transform proves this for both C and E; simply assuming
the stopped channel is self-adjoint would be incorrect.

For the declared cap, 496 local ticks suffice for the precision24 warm-up.
The **next specified incident tick** retains both movement and failure.
The uniform-input movement probability is `5/33`; finite preparation gives
an interval around it. Complement repair then selects a faithful regional
reference with minimum eigenvalue `1/120` and cap modular gaps `log(23/22)`
and `log(5/4)`. These differ from the globally heralded reference's gaps.
The receipt composes conditioning, finite complement repair and modular
generator errors, including arbitrary initially untouched references.

Warm-up and heralding use no global attempt counter. Regional record
collection, an IID seam schedule and a fixed local horizon are declared.
The complement controller counts its own repairs and reports completion;
those remote repairs produce no incident tick at the cap. Retained preparation records must be ignored by
the continuation; the bound does not assert independence from that tape.
The reference defines a modular clock, whose active implementation and
physical calibration are separate obligations.

```sh
python -m model.native_clock.regional_preparation --check
```

The [single-carrier clock comparison](COMPARISON.md) additionally derives
relative rates and explicit finite uncertainty for all twelve local current
clocks, without stationary preparation or a federation. Their asymptotic
count covariances recover all three native port-distance shells. These
unit-independent relations are testable after an operational apparatus
identification; choosing seconds is a convention, not the missing result.

The [measured host-clock attachment](HOST_CLOCK.md) records actual native
execution against a monotonic host counter, calibrates one local clock on
one block and checks a different clock on three untouched blocks. Their
duration errors are 6.29%, 2.29% and 7.76% against a declared 25% engineering
tolerance. A host-workload control preserves every native record and rejects
that calibration in all four blocks. This is an executable apparatus check;
it supplies no traceable lab calibration or identification with proper time.

The [joint clock bridge](CLOCK_BRIDGE.md) derives how incident interactions,
including zero currents, convert to local movement counts. Their stationary
ratio is `33/5`, but their exact correlations and finite uncertainty rule out
using that mean as a pathwise conversion. The bridge retains field reads and
extra feedback hits, and replays the existing host capture with its unchanged
calibration: all three baseline blocks remain inside the propagated interval,
while the corresponding workload controls reject it. This uses the same
total-two carrier as the capture, not the full binary preparation of the
field-refinement protocol.

## Active configuration readout and basis preparation

[`configuration.py`](configuration.py) supplies a different, active protocol:
a controller chooses seams from its recorded currents and resolves the whole
binary total-two configuration, without an IID sampler or an initial-bit
disclosure. It first scans for a nonzero current, which fixes two opposite
endpoint bits. It then grows a connected known region. To read an unknown
neighbour, it probes from a known endpoint; if no movement is seen, it routes
an opposite known bit to that endpoint and probes again. One of these two
probes has unequal endpoints. A moving current resolves the neighbour; an
exhausted retry budget produces failure, never a guessed bit. Routing itself
uses only witnessed native movements on known unequal endpoints.

Every primitive C/E map is a weighted partial injection. Once the current
tape resolves the configuration, the cumulative signed currents also resolve
its initial label `x`. Each successful fine path therefore has Kraus form
`a_h |y_h><x|`. Native feedback can restore `x`, or route to any chosen binary
total-two target `t`. Discarding the private fine tape gives, respectively,

```
restore:      rho -> sum_x |x><x|_record tensor s_x P_x rho P_x,
reset:        rho -> |t><t| tensor sum_x s_x <x|rho|x>.
```

Without restoration, the decoded initial label is the same but the terminal
work configuration can depend on the path. These formulas include arbitrary
untouched references. Here `s_x` is the
actual finite success probability for input `x`, not one. The complete
failure branch has the missing mass. Each adaptive node applies a complete
C/E instrument, and every finite terminal branch is retained, so induction
over the controller tree gives Kraus completeness. Comparing against ideal rank-one
measurement or replacement, with an empty ideal failure branch, gives
half-diamond error at most `max_x(1-s_x)`: the missing successful block and
the actual failure block each contribute half that mass. Conditional
success can steer a correlated reference; no independence is inferred from
postselection. Discarding a record is an explicit access choice, not erasure
of its physical cost.

The exact source replay checks the controller on every basis input and
retains its failure probability and maximum native-attempt budget. The
source configuration graph has diameter four. With `a` seed scans and `r`
attempts per mixed-edge retry, reset costs at most `30a+79r` attempts and
has failure at most `2^(-8a)+69*2^-r`. The ten growth stages need at most
`1+...+10=55` swaps and twenty probes, followed by four terminal swaps.
At most ten growth probes are mixed, giving 69 fallible blocks; every seed
scan has at least eight mixed seams. Four scans and 32 retries give
**2648 attempts and error at most `70/2^32`**. One initial 30-seam covering
sweep absorbs arbitrary full total-two collision inputs; keeping its signed
currents preserves retrodiction. Full-carrier reset therefore costs at most
2678 attempts. Binary restoration is allowed; recreating collision states
by these binary repairs is not.

`python -m model.native_clock.configuration` prints its receipt. `--check`
checks the retained bytes; `--write --output PATH` creates a new file and
refuses to overwrite existing evidence.

The current-only controller stores at most 2678 `(seam,current)` records,
78 candidate-original/current mappings and twelve known-bit entries. A
packed encoding of its tape, length counter, support map and known bits
uses at most 19,328 bits; this excludes routing scratch, immutable program
tables and host-language overhead. These are logical storage components,
not a measured hardware RAM or thermodynamic erasure cost. Discarding access
to a preparation tape does not make its production free.

The carrier, fair C/E instrument, seam control and record memory are supplied.
No coherent preparation follows from this classical protocol. Coarsening
its fine measurement to one port is **not** a Lüders port measurement: it
also destroys coherences between different configurations with the same
port value. The [auxiliary-pointer construction](../native_interface/README.md)
avoids that disturbance when a continuing work-bit instrument is required.

## Why the model does not predict seconds

The source contains integer attempt counts and dimensionless probabilities,
but no quantity with units of seconds. The executable
`duration_scale_identifiability` audit assigns illustrative values of `1/100`,
`1` and `100` seconds per attempt. They leave the exact one-attempt record law
`(28/33, 5/33)` unchanged, while implying different event rates per second
(`500/33`, `5/33` and `1/660`). Thus the native law fixes records per attempt;
it does not fix an attempt duration or a laboratory event rate.

No laboratory apparatus or public dataset has been operationally matched to
this current record. Under a separately validated constant attempt-rate map,
the minimal scale calibration is one nonzero native-attempt difference paired
with its independently measured lab interval. One pair alone does not justify
that the same rate holds over other intervals. Without such an apparatus map,
fitting generic lab event data would not test this model.

## What runs and what is measured

The program prepares the declared stationary mixture `I66/66` of the native
binary total-two carrier. At each attempt it selects one of the 30 source
seams uniformly, uses the native fair placement coin, and executes the existing
rule C. The readout retains `bit = abs(j)`: one for a moving current, zero for
a no-current outcome. Seam, coin and current direction are probability-summed
inside the complete instrument. The current counter is a function of those
retained bits; the observer never needs to inspect the hidden configuration or
its token-separation class.

The [producer](clock.py) executes the actual 120 native Kraus branches through
the shared [primitive observer runtime](../primitive_observer/runtime.py).
One five-attempt history retains all 32 bit words and their complete
unnormalized 66-dimensional terminal density matrices. Every word probability
and terminal density is compared with the exact native transition calculation.
A saved two-attempt checkpoint is serialized, independently replayed and
continued to the same complete output.

First-current laws are extracted from the executed prefixes: `1`, `01`,
`001`, `0001`, and the remaining `0000` tail. Their probabilities sum to one.
The analogous next-current law after a retained first `1` is explicitly
conditional on that first event, whose probability is `5/33`. The full parent
history, including first `0`, remains present. A physical apparatus that halts
its native execution at the first current is not implemented in this packet.
The derived stopping law describes those retained prefixes; it does not treat
later fixed-horizon terminal states as the stopped states.

## Exact source reduction

The occupied ports have three separation classes on the supplied icosahedral
carrier: adjacent, distance two, and antipodal, with 30, 30 and 6 configurations.
For **every** configuration in a class, summing each marked native transition
over each destination class gives the same row. The producer checks this
strong lumping property for all 3,960 primitive branch addresses. With rows
acting on class probability vectors, the exact marked matrices are

```text
        [13/15   0     0  ]          [1/15  1/15   0  ]
T0  =   [  0    5/6    0  ]    T1 =  [1/15  1/15  1/30]
        [  0     0    5/6 ]          [  0   1/6    0  ]
```

`T0 + T1` is stochastic. Its stationary class law is
`pi = (5/11, 5/11, 1/11)`. The full word law is
`P(b1...bn) = pi T_b1 ... T_bn 1`. This is an exact relation to the original
native kernel, not a fit to measured words. It also holds for any initial
class distribution; tests execute nonstationary basis preparations separately.
The reduction concerns this classical record family. It does not reduce the
full quantum observer algebra or choose between coherent laws C and E.

Counts follow from the coefficients of `pi (T0 + z T1)^n 1`. The retained
[receipt](results.json) includes **all** count probabilities at 1, 2, 8 and
32 attempts, including zero, rather than a conditional histogram of successes.
It also retains first-current probabilities through 16 attempts and the exact
remaining no-current tail.

## Exact model laws for new measurements

| Observable and preparation | Exact prediction |
| --- | --- |
| Moving bit at a stationary attempt | `5/33` |
| First next current, sampled at an arbitrary stationary attempt | Mean waiting count `147/22` attempts |
| Interval after an actual current, in the stationary event-conditioned law | Mean `33/5` attempts; variance `951/25` attempts squared |
| Two adjacent event-conditioned waiting intervals | Covariance `9/100` attempts squared |
| Next bit is one, given last bit one | `23/150` |
| Next bit is one, given last two bits `11` | `53/345` |

The distinction between the two waiting preparations is measurable. A
stationary sample encounters the more slowly moving states more often than a
sample taken immediately after a current. The latter, conventionally called
the *Palm* preparation, has class law
`nu = pi T1 / (5/33) = (2/5, 1/2, 1/10)`. It is obtained here by conditioning
an actual retained source event, not by assuming a fresh independent clock.

No-current outcomes leave the classical configuration unchanged. The next
current therefore has a geometric waiting law conditional on that configuration.
For integer `m >= 0`, the survival laws are exactly

```text
P(W > m | stationary sample) = (5/11)(13/15)^m + (6/11)(5/6)^m
P(W > m | after a current)   = (2/5) (13/15)^m + (3/5) (5/6)^m.
```

These formulas include the no-current tail at every finite horizon. Their
geometric sums yield the means and second moments above. They are not
independent waiting intervals: after moving, the class evolves according to
the derived jump matrix

```text
J = [1/2  1/2   0 ]
    [2/5  2/5  1/5]
    [ 0    1    0 ].
```

Writing `m = (15/2, 6, 6)` for the conditional mean waits, stationarity
`nu J = nu` gives the adjacent product moment `sum(nu_i m_i J_ij m_j)`.
Subtracting `(33/5)^2` gives the nonzero covariance `9/100`. Replacing this
process with independent equal-rate ticks loses that observation. Likewise,
the two different conditional next-bit probabilities in the table rule out
a first-order Markov description on the current bit alone. The three hidden
classes retain predictive memory without adding a clock module.

## Which scheduler premises disappear

The bare source graph supplies a narrower scheduler selection. Two checked
port automorphisms move any of its 30 seams to any other. Consequently a
normalized **one-step seam law invariant under those source symmetries** has
weight `1/30` on every seam. This is an exact symmetry implication, not a
claim that the adaptive two-event observer uses that law. It also says
nothing about independence across attempts. The
[symmetry/scheduling gate](symmetry_schedule.py) executes a stationary
countercontrol: repeat the last seam with probability `q`, otherwise redraw
uniformly. For every rational `0 <= q < 1`, each seam has uniform
one-step exposure, every fixed-seam native kernel preserves the same uniform
total-two configuration law, and all seams remain reachable. Yet the
probability of two successive moving currents is
`(23+52q)/990`: `23/990` for IID seams and `49/990` at `q=1/2`, while each
single-attempt moving rate remains `5/33`. Within this declared scheduler
family, the retained current bits determine `q=(990 P(11)-23)/52` from their
exact two-event law; finite-sample uncertainty needs a bound. Thus
source symmetry and stationarity reduce the uniform-marginal premise
but do not select the temporal scheduling law required by the complete field
and clock records.

The **stationary mean rate `5/33` does not require uniform seam choice**.
For each fixed seam exactly 20 of the 66 configurations have unequal endpoint
occupancy, and the fair coin produces a move with probability one half.
Every seam thus has the same stationary moving probability `20/(66*2)`.
Its native kernel is symmetric, so every fixed normalized mixture of these
seam kernels preserves the uniform stationary law. This proves the mean
identity for arbitrary fixed, state-independent seam probabilities. It does
not cover adaptive or state-dependent scheduling.

If all seam probabilities are positive, every configuration has a positive
escape probability `r_i`. Symmetry gives the event-conditioned configuration
law `nu_i = r_i/(66*p)`, where `p=5/33`. Consequently
`sum(nu_i/r_i) = 1/p = 33/5`. The observer's mean event interval is then fixed
without choosing uniform seam weights. Positivity matters: selecting only
one seam leaves 46 configurations forever inactive. That stationary mixture
has mean move rate `5/33`, but its event-conditioned mean interval is
`2`, not `33/5`. This exact inactive-component negative control is retained.

The full measurement law remains scheduler-sensitive. Doubling one seam's
relative weight, with weights `(2/31,1/31,...,1/31)`, leaves the stationary
mean and mean event interval unchanged but changes the observable two-current
probability from `23/990` to `1477/63426`. The stopped configuration transition
law changes too. Mean clock normalization therefore does not eliminate the
whole scheduling premise.

An additional declared idle choice with probability `1-q`, independent at
each attempt, gives `T0' = (1-q)I + q T0`, `T1' = q T1`. The stopped moving
kernel is unchanged because `T1'(i,j)/r_i' = T1(i,j)/r_i`; its stationary
event-conditioned law also stays fixed. This removes only global idling from
the stopped transition sequence. Zero records, waiting counts and attempt
costs change. The producer checks `q=1/2` on the full 66-state kernel;
the cancellation establishes the identity for any `0<q<=1`.

## Passive age from two current-derived endpoints

The [passive-age controller](passive_age.py) estimates a fixed unknown number
of unread native repair attempts. It requires no coherent phase operation,
reset, stationary preparation or independent copies. The same work and
auxiliary carriers continue across trials. The observer receives endpoint
current records and failures; the interval scheduler's records and stopping
gate are private. Giving those scheduler counts to the decoder would make
the estimation problem trivial.

First the [configuration controller](configuration.py) discovers a binary
configuration through actual signed currents. Choose the antipodal pair
containing its first witnessed occupied port. Its disclosed initial
occupation `Y0` is either one or two. After the unread interval, a second
configuration read retrodicts the endpoint occupation `Yn`. The second read
disturbs work; its final configuration is not substituted for its retrodicted
input.

For each of the six antipodal pairs the exact native Heisenberg identity is

```
P*(Y - 1/3) = (9/10)(Y - 1/3).
```

The integer source audit covers all 66 binary configurations and both C/E
instruments; partial-injection Kraus supports also establish the off-diagonal
identity. With ideal complete endpoint reads, the score
`Z=(Yn-1/3)/(Y0-1/3)` has conditional mean `(9/10)^n` for every initial
configuration and every past transcript. The
conclusive initial read supplies that configuration from arbitrary quantum
input and untouched references. No assumption that successive trials are
independent is used.

Every failed start read, scheduler selection or endpoint read contributes
score zero and remains in the denominator. The score lies in `[-1/2,5/2]`.
If each discovery fails with probability at most `delta`, each bounded
seam selector fails with probability at most `tau`, and `n<=W`, its complete
conditional mean differs from `(9/10)^n` by at most
`(7/2)delta+(5/2)W tau`. A martingale Hoeffding interval for the observed mean,
enlarged by this bias, gives a finite set of integer ages. The decoder takes
only complete record counts and the declared window. Empty or multiple
candidate sets are reported explicitly.

Uniform seams are chosen by actual coarse-coin operations on a reusable
auxiliary. The selector's exact `I/2` coin effect removes an independent
random-number service even for initially correlated carriers. Its attempts
and timeouts are counted separately from work attempts, including the costly
endpoint reads. The retained run executes 65,536 continuing trials per law
in the window zero through four; a separate analytic budget covers the
larger window zero through twelve. The latter is not an executed history.
For the executed design, the coverage failure is at most `1/4096` per law
and the mean radius is smaller than half the minimum separation of the
five candidate means. Its worst-case budget is 341,049,374 native attempts,
including endpoint discovery and auxiliary selection. These service costs
are not the unknown interval length.

Both the public receipt checker and the CLI own a fresh source replay. A
matching pair of supplied receipts cannot certify its own stream digest,
terminal carriers or operation counts. Structural checks additionally enforce
native support, current-prefix continuity and complete resource identities.

This is a destructive ensemble readout of a delimited interval, with a fixed
unknown length across trials. It does not generate that interval, supply an
individual monotone clock, or identify physical seconds. The fixed-length
condition cannot be dropped: a mixture of ages zero and two with probabilities
`9/19` and `10/19` has the same score mean as age one but mean age `20/19`.
This is a failure of this mean inverse, not an asserted equality of complete
record distributions. Joint current access, finite memory and the declared
native instrument are part of the observer architecture.

## Existing theorem attachment and limits

The existing
[`SourceRecurrenceClock.returnCount_block_additive`, `returnCount_mono` and `returnCount_le`](../../Lean/QFT/SourceRecurrenceClock.lean)
apply directly to the current record: take `i=1` and the two-valued trajectory
`omega(k)=abs(j_k)`. These count theorems accept any such record sequence; they
require no Markov assumption. They give adjacent-block additivity, monotonicity
and the bound `N(n)<=n` for this instantiated observer counter.

The nearest waiting-law result in that file, `meanReturn_hasSum`, belongs to
a different committed two-state chain. Its numerical return law is not
imported here. The advance is the checked native-kernel mapping, exact marked
reduction and directly calculated waiting/count distributions compared with full
executed CQ output. These are exact finite arithmetic checks and direct
mathematical derivations; no new Lean proof is claimed.

The [count-clock enclosure](../../Lean/Time/SourceCountClockEnclosure.lean)
assumes independently justified count-volume errors and a fourth-power
proper-time/volume relation. The present current counter supplies neither.
It therefore establishes no spatial dimension, Lorentzian geometry, proper
time calibration or physical frequency. The supplied carrier, native rule,
stationary preparation and chosen current observable are explicit.

The related [`relative_modular_dynamics`](../relative_modular_dynamics/README.md)
packet derives a periodic Y-readout response to the Palm state's modular
parameter from an existing coherent preparation. That gives a local modular
phase readout under the declared ideal phase target. Its exact boundary test
also shows why the current native rule cannot supply that phase evolution:
real coherent preparation and real repair branches leave the Y readout
unbiased at every integer native-attempt count. Native attempt counts and the
modular parameter are distinct model readouts, with no derived conversion
between them, and neither specifies physical duration or a lab clock.

The new laws are predictions of the executable source model. No laboratory
measurement has yet been operationally identified with this native current
record, so this packet makes no new empirical agreement claim and changes no
frozen anharmonic result. Its useful next target is a justified interpretation
of a measured event counter, with the joint count and waiting laws checked
together. Matching the mean alone would miss the distinguishing correlations.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.native_clock --check model/native_clock/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/native_clock
```

`--write PATH` creates a new receipt exclusively and refuses to overwrite an
existing one. Exact probabilities are rational strings. Floating CQ residuals
are replay diagnostics with tolerance `2e-11`, not outward-rounded numerical
certificates. Public source and theorem hashes bind the retained computation.
