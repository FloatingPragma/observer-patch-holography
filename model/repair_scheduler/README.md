# Preparation and scheduling supplied by native repair

The observer can prepare its native working state and choose random seams
without an external random-number service, an exactly stationary initial
state, an initially independent auxiliary, or repeated auxiliary resets.
This packet constructs those operations from the existing C/E instruments,
executes complete continuing observation histories, and accounts for their
errors and finite resources. The fair native instrument, two carriers,
classical control and the existing coherent mixer are specified inputs.

Two results do different jobs. A deterministic sweep prepares the working
carrier. A separate auxiliary carrier's coarse coin records implement the
uniform random seam choices needed during the existing phase experiment.
The sweep is used only for preparation: its dynamics are not substituted for
the random-seam clock, transfer or field laws.

## Native preparation without a random seam schedule

Use the actual twelve-port graph, its thirty seams, and total load two.
Visit every seam once in the fixed source order, discard all preparation
records, and repeat. A covering sweep sends the entire 78-dimensional
carrier into its 66-dimensional binary sector exactly. A collision at a
port stays there until its first incident repair splits it into two unit
loads. Every port is visited, and subsequent repairs cannot recreate a
collision. The conclusion holds for arbitrary input densities and untouched
references: zero population outside the binary sector forces the positive
output's corresponding rows and columns to vanish. No support is deleted
numerically and no output is renormalized.

On the binary sector let P be the product, in the specified order, of the
thirty classical heat-bath kernels obtained by calling the endpoint rule
with both placement coins. Each kernel is doubly stochastic. For m sweeps,
`classical_power` calculates P^m with integer numerators and denominator
2^(30m), and checks every row and column sum. Set J_ij=1/66 and

```
Q_m = P^m-J,
a_m = max_i sum_j |(Q_m)_ij|,
b_m = max_j sum_i |(Q_m)_ij|.
```

The singular norm obeys `||Q_m||_2 <= sqrt(a_m b_m)`. The implementation
encloses this square root upwards with exact rational arithmetic. P need
not be symmetric or reversible; a random-scan spectral gap is not imported.

For a matrix unit `|x><y|`, the signed load difference z=x-y is unchanged
by every surviving native branch. On a seam (a,b), when z_a=z_b=0 the
shared particles undergo a doubly stochastic lazy exchange; when
z_a=z_b is nonzero nothing changes. For C, a difference of one attenuates
by 1/2 and a difference of two kills the unit. For E, every nonzero
difference kills it. Every nonzero z meets a boundary. The surviving C
sectors have at least ten difference-one seams, so their full-sweep
Hilbert--Schmidt contraction is at most 2^-10. E removes every
off-diagonal sector in one sweep.

There is also a certificate uniform over the entire admitted
[current-record family](../record_laws/README.md). Positivity of each current
Gram matrix gives `|G_j(x,y)| <= sqrt(p_j(x)p_j(y))`, the existing C0
entrywise envelope. Difference-one and difference-two seams contribute
at most 1/sqrt(2) and 1/2 respectively. The weighted cut count
`c_1+2c_2` is at least ten on every nonzero difference sector, so a full
sweep contracts the envelope by at most 1/32. Its shared-particle factors
are doubly stochastic, and hence contract Hilbert--Schmidt norm even when
multiplied in different orders. Complex phases cannot exceed this positive
pathwise envelope. As in the existing law-uniform burn-in proof, classical
law-only feedback from earlier coarse currents is covered, provided the
seam order stays fixed and its controller and preparation tape are
discarded. That law-selecting controller starts fixed and independent of
the source/reference. Arbitrary coherent environment feedback is outside
this grammar; it is distinct from the later auxiliary coin construction.

Write c=2^-10 for C, c=0 for E, or c=1/32 for the whole family. For every
binary input rho_AR and arbitrary untouched reference R,

```
D((Sweep^m tensor id_R)(rho_AR), (I66/66) tensor rho_R)
 <= min(1, 33 max(sqrt(a_m b_m), c^m), b_m/2 + 33 c^m).       (1)
```

Here D is half the trace norm. For the second bound, the diagonal reference
blocks rho_jj are positive and their traces sum to one, giving a contribution
at most b_m/2. The off-diagonal operator is Hilbert--Schmidt orthogonal to
the diagonal. For a pure input its reference Schmidt rank is at most 66;
trace norm is at most 66 times Hilbert--Schmidt norm on the resulting
support. Convexity proves the same inequality for arbitrary mixed inputs
and reference dimensions. The first bound applies this argument to the
whole traceless part. For a subnormalized input of trace p, multiply the
bound by p. For a full 78-state input, first execute one absorbing sweep
and apply (1) to the remaining sweeps. With classical law feedback, apply
the bound separately conditional on the controller state after absorption,
then sum the subnormalized branches. Their reference marginals sum to the
original marginal, so the correlation created during absorption adds no
independence assumption.

**Seventeen sweeps, 510 native attempts, certify distance at most
0.000000259 from the binary uniform state, including arbitrary untouched
references and arbitrary initial collision population.** The exact retained
bound is about 2.58673903e-7. This is a sufficient preparation budget for the
actual fixed sweep, not a claim of the fastest mixing protocol. General
reference states, thermal Gibbs states and physical units are not selected
by this tracial preparation.

## A reusable native source of scheduler bits

Put an auxiliary carrier A alongside work and reference WR. On one fixed
auxiliary seam let M_c be the native C/E CP map that retains coin c and
forgets the current (or contrast). Source enumeration proves

```
M_c^*(I_A) = I_A/2,
Tr_A[(M_c tensor id_WR)(rho_AWR)] = rho_WR/2.                  (2)
```

Every fixed-record source map is a partial injection with coefficient
1/sqrt(2). For either fixed coin, each input column occurs once across the
hidden record groups; injectivity makes every off-diagonal effect zero.
This proves (2) on both the 78-state carrier and its binary restriction,
for every source seam. No initial factorization between A and WR is used.

The normalized map 2M_c is trace preserving on A. Apply (2) conditionally
at each subsequent invocation, including after work operations chosen from
earlier coins. Every length-k coin word has reduced WR map `2^-k id`,
interleaved with the declared work operations. The auxiliary may keep
arbitrary quantum memory and be reused indefinitely; it is never reset.
This is a property of the specified repeated instrument, not an inference
from one-call statistical fairness. A hidden fair bit reused forever is a
counterexample to the latter inference.

Read five coarse coins as a binary integer. Values 0 through 29 select the
corresponding seam; reject 30 and 31 and repeat, for at most M rounds.
For every input and untouched reference, the exact reduced CP weights are

```
selected seam e at round r: (1/16)^(r-1) / 32,
timeout:                   (1/16)^M,
selected seam e, all rounds: (1-16^-M)/30.                    (3)
```

Every accepted seam is exactly uniform. Timeout is a retained outcome with
unchanged work state, not a discarded failed draw. Each request costs at
most 5M auxiliary native attempts. Ordinary modulo-30 decoding would be
biased and is rejected by an exact control.

The auxiliary quantum output and its fine current/contrast records are
outside this coarse interface. Keeping or conditioning on them can steer
an initially correlated work/reference state and invalidates (2) as a
product-output claim. The auxiliary must also be disjoint from the work:
repairing the work itself generally changes it. An auxiliary whose fixed
seam is empty can supply fair coins while making no movements at all;
coin counts therefore do not supply a movement-clock rate or physical time.
The native fair instrument already contains stochasticity. This result does
not derive it from the deterministic endpoint function alone, or remove
the blank/service requirements of its separate source compilation.

## One continuing complete experiment

The implemented controller uses seventeen preparation sweeps before each
trial. It then executes the existing native signed-current herald, balanced
mixer, three unread repairs, second mixer and native coin/current pointer.
Each of the three random seams is selected by (3) with M=7. A selector
timeout halts that trial at its actual disturbed state; the next trial
starts from that state and applies the same preparation. Failed heralds
also are trials. No initial state is replaced between trials.

The public record is `(herald, stage, bit)`: failed heralds have stage zero;
stages one through three mark the corresponding selector timeout; stage
four contains the completed pointer bit. Both useful herald signs, both
pointer values, six timeout records and the impossible failed-herald/read-one
record are retained: twelve outcomes per trial and 144 for two trials.
Preparation tapes and successful scheduler scratch records are deliberately
forgotten. "Complete" refers to this specified public interface, not an
independence claim about every discarded microscopic record.

The producer executes two continuing trials for each C/E law from a complex
full-sector input with collision probability 9/25. It retains every joint
record probability and checks the conditional second-trial laws. A separate
full-sector input entangled with a reference qubit is propagated through
preparation and every complete trial output. The accepted work channels are
enumerated from actual native calls, and an independent existing source
enumerator checks their entire reference output. All auxiliary 32-word
trees at every rejection round are separately executed on one entangled
auxiliary without resets. Equation (2) justifies tracing it out in the
main histories; a dense two-carrier joint state is not claimed to have been
simulated throughout that calculation.

Trace-distance contraction makes (1) valid after **any common subsequent
complete CPTP instrument**, including adaptive feedback and quantum outputs.
For the present trial, coupling its bounded selector with the unlimited
uniform selector adds at most `3*16^-7`. Conditioning on each previous
public history gives the same bounds, because the next preparation accepts
every incoming state and the auxiliary identity accepts every correlation.
Telescoping therefore gives, for N trials, complete transcript error at most

```
N * (259/10^9 + 3/16^7).                                    (4)
```

No IID trial or reset-state premise is added. The known six-record ideal
C/E likelihood experiment is embedded by assigning zero ideal mass to the
new timeouts. Retain its fixed N=3927 and three-attempt menu, whose ideal
wrong-candidate bound is 1/200. Any timeout makes the actual classifier
return `scheduler_timeout_inconclusive`; an impossible record returns
`outside_both_declared_models`. All twelve nonnegative integer counts and
the fixed total are mandatory. The wrong-candidate probability is bounded by

```
1/200 + 3927*(259/10^9 + 3/16^7) < 0.006061.                 (5)
```

The scheduler-timeout probability has the separate union bound
`3927*3/16^7 < 0.000044`. The complete campaign is a proved resource/error
schedule, not 3927 executed or measured trials, and not a new frozen test.
The former design's minimum-cost claim is not transferred to this apparatus.

The worst-case campaign uses 2,434,740 native attempts, counting both
carriers, and the same 267,036 mixer pulses and 534,072 proper-code
load/export events as the existing menu. At most 510 working preparation
attempts, five working herald/probe/pointer attempts and 105 auxiliary
attempts occur per trial. Source mixer arrays are freshly emitted and
compared with their full-code target using the existing `1e-10` two-call
envelope. Equations (1)--(5) concern the exact native/ideal-mixer model;
raw-array comparisons and floating history residuals are separate from
physical service fidelity. No hardware accuracy or duration is certified.

## Which finite premises this supplies

| Existing interface | Concrete realization here | What remains |
| --- | --- | --- |
| Native state and quantum channel | Actual full total-two input, native Kraus maps, exact absorbing support map and reference-uniform preparation bound | Total-two architecture and selected native instrument |
| Stationary-preparation input to complete observation bounds | Finite deterministic sweep, (1), then arbitrary complete readout by trace-distance contraction | Preparation records must be discarded; the finite state is approximately stationary |
| Uniform independent seam choices in the phase experiment | Auxiliary coin effects (2), bounded selector (3), timeout correction (4) | Auxiliary access, fair native instrument and classical controller |
| Born readout, retained records and feedback | Complete CP branches and unchanged copies of earlier classical labels; subsequent Kraus maps act on work and read those labels | Record storage, coherent mixer and its proper-code services |
| Repeated-trial prediction without IID preparation | Conditional preparation and selector bounds, complete timeout-aware classifier, (4)--(5) | Apparatus accuracy and a physically identified campaign |

Appending a committed classical label writes a new register; all later
work maps are controlled by, and copy, that label. The corresponding
projector intertwines every later Kraus operator, so old evidence is
durable after writing. The tensor identity on an untouched reference or
unused buffer preserves these identities and error bounds. This supplies
these particular finite channel/record interfaces, not all premises in the
paper library or a spatial/continuum refinement theorem.

The empirical obligation is separate. The existing
[public-data audit](../../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
retains actual measured discrepancies and supplied physical dictionaries.
Its oscillator preparation uses a capped Fock mixture, not the native
uniform binary carrier prepared here. Replacing its calibration uncertainty
or declaring agreement with every known measurement would require a new,
justified apparatus connection. This packet supplies no such claim, no new
laboratory data, and no derivation of a thermal Hamiltonian, radial geometry
or observed physical units. It removes identified preparation and scheduler
inputs in one complete observer without relabeling the remaining inputs.

## Reproduce and challenge

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.repair_scheduler --check model/repair_scheduler/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/repair_scheduler model/native_equilibration/test_sweeps.py
```

The exact local matrix-unit, integer-kernel and coin-effect checks support
the analytic bounds. Independent native Kraus oracles, entangled reference
histories, impossible/timeout outcomes, fine-record steering, biased modulo
sampling, shared-bit memory and malformed-input controls challenge the
implementation. Exact counts, structures, source bindings and rational
bounds replay exactly; only floating diagnostics use the established
comparison tolerance. Separate domain checks reject negative probabilities,
nonzero mass for impossible records, and incomplete probability mass even
when individual entries would fall within that tolerance.
