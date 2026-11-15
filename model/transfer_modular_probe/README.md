# An executable modular clock for the native record-transfer description

The native current records determine a finite transfer description whose
modular flow can be executed through the existing coherent source compiler.
No oscillator action, mass or external clock calibration enters the phase
generator. A native incident record can trigger one fixed modular step through
the same compiler. The coherent code, preparation and event-to-probe coupling
are apparatus choices. This is a finite operational connection, without a
claim that native dissipative waiting is unitary evolution or that a continuum
spacetime has been recovered.

## From actual records to the reference and generator

The [native clock](../native_clock/README.md) proves marked strong lumping of
all 66 native configurations into three token-separation classes. Its exact
transition matrix and stationary class weights are

```text
P = [[14/15, 1/15, 0], [1/15, 9/10, 1/30], [0, 1/6, 5/6]],
pi = (5/11, 5/11, 1/11).
```

Reversibility makes `S=D^(1/2) P D^(-1/2)` symmetric, where `D=diag(pi)`.
The existing full transfer has exact bounds `(2/3)I <= K <= I`; its invariant
class-constant restriction is positive as well. The source characteristic
polynomial factors as

```text
(x-1)(x^2-5x/3+311/450).
```

Thus the three eigenvalues are `1`, `lambda_+=(25+sqrt(3))/30` and
`lambda_-=(25-sqrt(3))/30`. On this reconstructed three-dimensional Hilbert
space, `sigma=S/Tr(S)=3S/8` is faithful and

```text
H = -log S,                 log sigma = -H - log(8/3) I,
sigma^(is) A sigma^(-is) = exp(-isH) A exp(isH).
```

The dual state evolution is `exp(isH) rho exp(-isH)`. This fixes the orientation
convention; the source phase is `+s E_j`, with `E_j=-log(lambda_j)`.
The reference is derived from the positive native transfer. This does not
prove that OPH uniquely selects this reference among all admissible ones.

The reconstruction uses the finite reflection-positive transfer argument in
[the native transfer derivation](../primitive_observer/DERIVATION_TRANSFER.md),
the exact marked reduction, and finite inner modular flow as in
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean). The additional
work is an explicit compatible source program and its numerical certificate;
no new abstract modular theorem or Lean proof is claimed.

This `M3` is a representation of the reduced history, not the native `M66`
working-state algebra. Its reference differs from the diagonal Palm reference
used by [the native modular probe](../native_modular_reference/README.md).
An equality of their physical interpretations is not assumed.

## Complete source execution

The spectral basis is encoded in addresses 0, 1 and 2 of two work qubits.
Address 3 is explicit padding. The complete source target is
`diag(1, exp(is E_+), exp(is E_-), 1)`. It implements the modular action on the
whole effective `M3` algebra, with a declared padding extension. On the block
algebra `M3 direct_sum C`, the positive extension `S direct_sum 1`, normalized
by `11/3`, produces the same action. This is not a new spatial dimension.

Each probe starts from the declared blank code and executes:

1. An emitted balanced source word preparing `(ket(0)+ket(j))/sqrt(2)`.
2. The complete Walsh phase word, whose phases come from the source spectrum.
3. An emitted X or Y analysis word and a complete four-outcome computational
   read, followed by the emitted inverse analysis word.
4. For the continuing experiment, another modular segment and another X read
   on that same disturbed state, retaining both records and all padding.

The ideal single-read probabilities at address zero are

```text
p_X(0;s,j) = [1+cos(s E_j)]/2,
p_Y(0;s,j) = [1+sin(s E_j)]/2.
```

At modular parameter `s=8`, the X probabilities are approximately `0.801845656`
and `0.277055154` for the two modes. There is no physical-seconds conversion.
Both quadratures give a phase reading modulo a period; a single random
outcome does not determine elapsed time. These selected-interval experiments
use programmed phases. The event-driven protocol below has no variable phase
input during execution.

The continuing `(4,1)` experiment retains 16 two-read records and uses 430
primitive source pulses, or 860 proper-code service events. For
`q_k=[1+cos(s_k E_j)]/2`, its nonzero ideal record probabilities are
`q1*q2`, `q1*(1-q2)`, `(1-q1)*(1-q2)` and `(1-q1)*q2`. Repreparing the original
state between reads gives a different law and is a regression control.

Every operation is carried by the shared validated `runtime.Program`.
A scoped pure-state Kraus evaluator applies its original binary arrays to
vectors and retains every output, including tiny padding leakage. This avoids
subtractive cancellation at nominally zero branch weights. It neither repairs
nor relaxes the generic runtime's relative positivity guards. An independent
generic-density execution checks the complete unconditioned terminal state.
The scoped checkpoint serializes the exact binary vectors and resources;
resumption verifies them against a fresh complete prefix. It is not asserted
to be the generic runtime's checkpoint format.

## A native event drives the modular step

[`event_clock.py`](event_clock.py) couples a continuing native carrier to the
same four-address probe. At each attempt, a uniformly selected native seam
executes its fair-coin repair. An attempt incident to port 0 produces the local
seam/current record and triggers **the same fixed `U(1)=exp(iH)` word**, including
when the current is zero. Remote attempts update the carrier but do not advance
the probe. There is no elapsed-time or adjustable-phase argument to this
controller. The unit coupling per event fixes the modular parameter convention.
The reference `S/Tr(S)` remains the three-class *attempt-transfer* reference;
it is not silently replaced by a port-incident transfer reference.

The clock carrier and coherent probe/service code are distinct regions of the
declared apparatus. The carrier scheduler pauses during source preparation,
modular gates, pointer reads, analysis undo and interventions. Source service
pulses are charged separately and do not recursively become clock events.
This atomic scheduling/coupling contract implements an active finite observer;
it is not an autonomous free-running hardware realization or a theorem about
the passive repair law. No unrecorded backlog is discarded.

Each of two incident intervals has an attempt cap `c=1` or `2`. The supplied
watchdog controls termination only; its hidden attempt count does not select
a phase or become a local clock reading. A timeout is a retained output with
the terminal carrier and probe state, not a discarded trial. After the first
incident event, the probe is read and its analysis basis undone. A positive
native current also triggers the fixed balanced gate as a noncommuting
intervention. The second interval continues both disturbed systems without
reset. Coins and remote records are unread and summed by their complete
instruments; the retained interface is the local seam/current tape, pointer
outcomes and timeout. Terminal native populations are an execution-certificate
output, not a claim of primitive local access to the whole carrier.

Let `A` be the remote native branch-count matrix and `B_m` the count matrix
for visible incident mark `m`. Rows of `A` sum to 50 and rows of `sum_m B_m`
sum to 10, out of 60 seam/coin branches. The complete stopped maps are

```text
M_m(c) = sum_{r=0}^{c-1} A^r B_m / 60^(r+1),
M_timeout(c) = A^c / 60^c.
```

Their sum is stochastic. Independent floor/ceiling path recursion checks
every entry for every starting configuration. For each cap and probe mode,
the execution retains all 3,661 joint output labels and exact terminal native
populations, including zero-probability labels and raw padding leakage. Its
native preparation is the uniform diagonal state on the 66 binary
configurations. Independent density evolution using the closed spectral roots
checks the complete probe output. The first timeout has probability `(5/6)^c`;
both ticks have probability `[1-(5/6)^c]^2`. These are unconditional complete
output laws, not laws conditioned on success.

There are at most eight source gates, four native attempts and two pointer
reads. The gate budget is 512 primitive source pulses and 1,024 proper-code
service events; those costs are distinct from clock-carrier attempts. For each
seam, the complete native `C` Kraus family consists of recorded
partial permutations with squared column sum one; hence it is CPTP, including
coherent native inputs. Exact native instruments, classical controls and
complete reads contract the channel comparison. Applying the existing
four-column gate bound at each possible source use gives

```text
half-diamond error <= epsilon(2+epsilon)/2
                     * sum_{k=0}^7 (1+epsilon)^(2k) < 1e-10.
```

This covers the complete output, arbitrary joint carrier/probe inputs and
untouched references, including record-controlled interventions and timeouts.
The executed witness is diagonal in the native carrier; a full coherent
66-state trajectory is not numerically claimed. Gate propagation residuals
are diagnostics, whereas the bound encloses the original binary gate arrays.

The no-reset control removes the optional current intervention. For mode 2
and cap 2, the unconditional event “both ticks succeed and the second pointer
read is zero” differs by more than `1/1000` when the probe is freshly prepared
between reads. The ideal difference is approximately `0.00145105936`; the two
compiled complete protocols add at most `4e-10` event-weight error: raw maps
need not preserve trace exactly, so each event error is bounded by the full
trace norm, twice its half-trace allowance. The proof
uses the exact phase enclosure, not the displayed decimal.

## Exact phase uniqueness and finite reading precision

The two nonzero frequencies are rationally independent. Suppose
`lambda_+^m lambda_-^n=1` for integers `m,n`. Conjugation in `Q(sqrt(3))`
exchanges the roots. Dividing the original and conjugate equalities gives
`(lambda_+/lambda_-)^(m-n)=1`; since their ratio exceeds one, `m=n`.
Their product is `311/450`, strictly between zero and one, so `m=n=0`.
Thus `E_+/E_-` is irrational. Equality of both exact phases at two real modular
parameters would give integer winding numbers with that rational ratio;
the parameters must coincide. A probe with nonzero amplitudes in all three
spectral modes has the same injective ideal phase orbit.

Exact uniqueness is not perfect clock resolution. Outward rational interval
powering of the single `U(1)` target certifies that any two integer readings in
`0,...,64` have maximum complex-phase distance greater than `1/4`. For the
ideal X/Y ensemble readouts of both modes, four probability-mean errors at
most `1/32` give per-mode phase error at most `sqrt(2)/16 < 1/8`, sufficient
for unique nearest-pair decoding in that window. Independent repeated known
preparations with `R` shots per quadrature give the elementary Chebyshev and
union bound `min(1,1024/R)`. This is an ensemble reading bound, not a one-shot
measurement of a continuously disturbed clock. Repeating the fixed compiled
step, with complete preparation/analysis/undo, uses at most 69 gates per
quadrature in this window and has source half-diamond error below `1e-9`.
This window concerns unread probe powers before a final quadrature read;
the retained joint native-history execution above covers two ticks.
Physical service and sampling errors require their own bounds.

At 272 ticks, both phases return within `0.04504` of their initial values,
certified by exact intervals. More generally, partitioning the phase torus
into `Q^2` boxes shows that for every positive integer `Q`, some integer
`1 <= n <= Q^2` has both phase chords at most `2*pi/Q`. Arbitrarily close
recurrences therefore prohibit an unlimited-window finite-precision reading,
despite the absence of exact nonzero aliases.

## Connection to measured native correlations

The same two modes give the connected covariance of two actual native
movement bits separated by `n` unread attempts:

```text
C_E(n) = w_+ lambda_+^n + w_- lambda_-^n,
w_+ = 1/7260 + sqrt(3)/16335,
w_- = 1/7260 - sqrt(3)/16335.
```

The producer checks this against exact full 66-state marked path sums at
`n=0,1,4,8`. The new X/Y source readouts reconstruct
`C_M(s)=sum_j w_j exp(is E_j)`, the corresponding analytic continuation.
The marked transfer matrices are **not POVM effects**; their measurement
semantics remain the native current instrument. The new coherent probes
measure spectral phases in the reconstructed code. Analytic continuation
plus its executable simulation is a supported correspondence, while agreement
of that interpretation with independent physical observations is open.

## Clock and energy scale

The roots `lambda_+` and `lambda_-` describe decay per **native attempt**.
The declared active controller uses their gaps `E_j=-log(lambda_j)` as phases
per **incident event**. This choice does not identify an attempt, an incident
event, a modular parameter increment and a laboratory second. The
[clock-scale sanity check](clock_scale_sanity.py) exhibits the ambiguity
without a fitted target: for any positive assignment `tau` seconds per
incident event, `omega_j=E_j/tau` and `t_n=n*tau` give the same native phase
`omega_j*t_n=n*E_j`. Changing `tau` rescales both assigned frequencies and any
conditional energy `hbar*omega_j`, while the exact native word law and the
dimensionless gap ratio remain fixed. The two gaps are transfer/probe modes;
they are not particle rest masses. An absolute particle-mass result would
also need a source-selected particle operator, a justified pole/rest-energy
readout and a common laboratory clock identification. The attempt-to-incident
waiting law alone supplies none of those identifications.

## Error and scope

Exact rational square-root and atanh-series intervals enclose the source
energies. A real-argument complex Taylor remainder encloses every target phase.
Each original binary source matrix is compared on all four input columns with
the complete ideal unitary. Its Frobenius bound supplies an operator bound
`epsilon <= 1e-11`, half-diamond bound `epsilon(2+epsilon)/2` and channel-norm
bound `(1+epsilon)^2`. The certificate retains amplification through every
source use and arbitrary untouched reference systems. The full largest
seven-gate experiment has half-diamond allowance at most `1e-10` against the
declared ideal code experiment; the intervening exact computational read
contracts the comparison.

That certificate covers raw compiled gates. Subsequent binary state propagation
and scalar probability residuals are diagnostic, and physical apparatus
service errors are unbounded. Tiny raw normalization defects are retained as
record weights. No empirical counts, interacting continuum limit or Lorentz
metric are supplied here. A further physical identification must preserve
observer-accessible correlations and intervention responses, without requiring
one universal clock or fundamental fields.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.transfer_modular_probe --check model/transfer_modular_probe/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.transfer_modular_probe.event_clock --check model/transfer_modular_probe/event_clock_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/transfer_modular_probe
```
