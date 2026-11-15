# Source-compiled apparatus in a complete observer-clock history

This construction implements the coherent preparation, union intervention
and terminal detector of the
[quantum observer correspondence](../quantum_observer_correspondence/README.md)
by words emitted by the retained six-level source compiler. Their original
binary output arrays are propagated through the same stopped native history.
Every classical word and its full quantum state remain in the output,
including unresolved waits and the declared complement-abort outcomes.

The parent result proves the quantum correspondence, its
arbitrary-reference extension and a joint clock/detector response. This
connection is an executable apparatus attachment with a complete-history
error certificate. It does not select a repair law, derive the apparatus from
native repairs alone, or identify a physical duration. Native C, the uniform
seam scheduler, fair coin and the compiler's source-control services are
specified inputs.

## State, operations and records

This specializes the [canonical primitive observer contract](../primitive_observer/README.md):
unnormalized quantum states are indexed by committed classical records,
complete instruments append those records, and feedback reads the tape.
`history.py` evaluates that composition in sparse exact arithmetic for this
fixed protocol. Its serialized prefix binds the source files, parameters,
resource counts and every work/reference component. Replay validates that
prefix before continuation uses its deserialized states.

The ideal native domain consists of the 66 binary total-two states of the
twelve-port source. Observer A accesses `(0,1)` and B accesses `(0,5)`.
Their union is `(0,1,5)`. A scheduled attempt selects one of the thirty
source seams uniformly and a fair placement coin. An incident clock ticks
when the selected seam meets its region, including zero-current repairs.
A retained native mark contains both clock increments and the three signed
occupation increments on the union. Hidden seam and coin alternatives are
summed as CP alternatives, never as coherent amplitudes.

The nine-qubit chart has address

```text
local_bits + 8 * spectator_label.
```

Its low three bits are the physical occupations at `(0,1,5)`, in that order.
The other ports are `(2,3,4,6,7,8,9,10,11)`. Within each remote occupation
shell of total zero, one or two, their sorted binary tuples receive a
six-bit label. Local occupation total determines the applicable shell.
This gives an injective isometry `V:C^66 -> C^512`, checked on all 66
source states. It is a separate chart; it is not the collision-inclusive encoding
used by the separate source-discrimination construction.

Write `P=VV*` and `Q=I-P`. Each native instrument is extended to the whole
512-dimensional register by its embedded code Kraus operators and a
separately retained Kraus operator `Q` reporting complement abort. These
operators satisfy the complete instrument identity on all 512 inputs.
An abort ends native evolution. It does not renormalize the code component
or discard the complement state. The terminal detector is applied.
This extension specifies what happens to compiler errors outside the code;
it is not an assertion of additional physical native states.

The ideal local gates preserve local occupation total, so the shell chart
implements their physical union locality exactly. An erroneous raw matrix
entry need not preserve that total. It can enter padding or another valid
code with a different shell interpretation. Both effects are included in
the global map error; the complement test alone is not a locality or error
certificate.

The experiment is the following complete finite history.

1. From a declared blank, execute the source preparation of
   `(|ports 0,1>+|ports 0,5>)/sqrt(2)`. The reference setting instead prepares
   `(|ports 0,1>|0>+|ports 0,5>|1>)/sqrt(2)`. The reference participates in
   this preparation and is untouched thereafter. There is no initial
   occupation disclosure.
2. Request one union incident, with a cap of two scheduled attempts. Retain
   a first-interval timeout or complement abort as a terminal status.
3. On first completion, execute native C feedback on `(0,1)` if the first
   recorded change at port zero is positive, and on `(0,5)` otherwise.
   Retain its current and count its incident ticks. Then execute
   `U_gamma=R_(1,5)(-1)^(gamma*n0*n1)`, with fixed `gamma=0` or `1`.
   In the ordered pair basis `(|10>,|01>)`, the rotation is
   `[[3/5,-4/5],[4/5,3/5]]` and fixes `00,11`.
4. Request a second union incident with the same cap, retaining every
   completion, timeout and complement abort.
5. On every terminal status execute the compiled binary read of occupation
   at port one. Retain both outcomes and the full 512-dimensional work
   matrix, together with the untouched reference when present.

On the ideal code, this is the parent 1,442-word experiment. The extended
cap-two output has 1,656 words, including the explicit abort words, whose
ideal probabilities are zero. No zero entry, timeout, detector outcome or
padding component is removed. Source pulse invocations are apparatus costs;
they do not add native incident ticks.

## What the compiled execution means

The compiler emits the existing source primitives
`exp(sum_p c_p D_p)` and contracts them through its declared proper-code
load/export services. Their exact invariant-code identities are inherited
from that compiler; its numerical leakage guard is not a hardware error
certificate. Keeping all 512 chart addresses does not assert that a separate
six-level export-error environment has been simulated. The preparation uses
24 source pulses, or 27 when it also prepares the reference. The intervention uses 24 pulses at both
`gamma` settings, including the requested zero-angle control operations.
The detector uses a three-pulse source CNOT to a blank pointer and retains
both pointer slices.

Let `M_gamma` be the returned complex `8 by 8` intervention array and
`L_0,L_1` the returned detector arrays. These are the actual arrays from
the emitted words, not rounded versions of a separately supplied ideal
matrix. Every real and imaginary binary component is converted to its exact
rational value. The history then propagates those rational complex values
and the exact native probabilities. Consequently there is no subsequent
floating matrix multiplication remainder in that retained execution.

The binary arrays need not be exactly unitary or trace preserving. They are
used as returned, without polar projection, symmetrization, normalization or
probability clipping. Their raw complete output is a CP matrix calculation
whose possible mass defect is part of its distance from the normalized
ideal instrument. The certificate is about these executed numbers. It is
not an interval proof for an unperformed analogue device, nor a conversion
of the source compiler's diagnostic tolerances into hardware guarantees.

The source certificate compares the full local arrays with the exact
rational intervention and the ideal stacked detector
`W=(P_0,P_1)^T`. Tensoring the resulting operator bounds with the untouched
spectator or reference identity does not enlarge them. This avoids a
factor of 64 from the chart and a factor from reference dimension.

## Complete adaptive error bound

For any matrix `M` and ideal isometry `U`, suppose
`||M-U||_op <= eta`. Define

```text
e(eta) = eta*(2+eta)/2,
g(eta) = (1+eta)^2.
```

For every operator `X`, including an arbitrary untouched reference,

```text
M X M* - U X U* = (M-U) X M* + U X (M-U)*.
```

Schatten submultiplicativity gives half-diamond discrepancy at most
`e(eta)`, and the CP map with Kraus operator `M` has trace-norm gain at
most `g(eta)`. The same argument applies to a stacked instrument:
discarding its coherent outcome cross blocks is a CPTP map, so stacking
`L_0,L_1` bounds their complete classical/quantum output.

Use the certified intervention and detector bounds to define
`eta_u, e_u, g_u` and `eta_r, e_r, g_r`. Let `sigma` be the raw prepared
state, `rho` its normalized ideal code target, and

```text
epsilon_p >= ||sigma-rho||_1/2.
```

This preparation error includes all local, reference and padding entries;
`sigma` is not assumed normalized or confined to the code. On code, the
first interval succeeds with the state-independent probability

```text
s = 1-(3/5)^2 = 16/25.
```

On the full chart, the effect of all branches that reach the apparatus is
exactly `s P`: complement inputs abort before that point. The pre-apparatus
feedback instrument is complete and preserves the code. Thus the complete
history before its terminal read, using the raw intervention, has gain

```text
G = 1+s*(g_u-1).
```

Indeed, its effect is bounded by
`I-sP+s*g_u*P <= G I`. For a CP map this effect bound gives the same
trace-norm bound on Hermitian operators, also after reference extension.
Replacing the intervention on the ideal prepared history costs at most
`s e_u`, because the unnormalized pre-intervention branches have total
trace `s`. All later native continuation, timeouts and aborts form a
complete instrument and contract that discrepancy.

It follows that the full returned classical/quantum matrices satisfy

```text
||H_raw(sigma)-H_ideal(rho)||_1/2
    <= g_r * (G*epsilon_p + s*e_u) + e_r.
```

This is a bound on half the sum of the nuclear norms of all unnormalized output
block differences. It holds for mixed inputs and arbitrary finite reference
systems whenever the corresponding preparation discrepancy is supplied.
There is no multiplication by the number of classical words, no successful
branch normalization, and no assumption that the raw preparation remained
in the code. All deterministic feedback uses the committed record only.
The proof would fail if a branch were selected from an unread occupation
or if complement aborts were dropped.

In particular the bound controls every subsequent common CPTP readout or
continuation. It also bounds half the entrywise absolute difference of the
complete record masses. Since the raw matrices are not assumed normalized,
the trace defect is bounded by twice the displayed error; calling them
exact probabilities is not needed for the certificate.

## Connection to the effective clock and field laws

The inherited native remote map has trace effect `3/5`, and its incident
maps have total effect `2/5`. The exact correspondence uses
`R=(2/5)(I-A)^(-1)` and its decoder on the transported matrix cones.
The ideal local intervention, detector and feedback commute with the
remote source maps. Its established all-matrix-unit identity therefore
identifies `H_ideal` with the decoded effective history, including arbitrary
references and finite-cap outputs. The displayed budget closes the
connection from emitted apparatus arrays to that exact effective target.

The raw compiled intervention can fail exact number preservation and
remote commutation. The result does not silently apply the unchanged
effective dictionary to that raw matrix. It compares its actual native
execution with the ideal history for which the correspondence is proved.
An exact effective execution of the raw intervention would require its own
transported intervention dictionary.

The inherited joint clock/detector distinction is therefore meaningful
for this apparatus execution up to the certified complete error. A common
record coarse graining contracts the bound. For comparisons of two
settings, add their complete error bounds to the inherited ideal
distinguishability statement. The coherent interaction is a chosen
source-control operation; its successful compilation does not make it the
Hamiltonian or action selected by native repair.

### Cap refinement and the completed-clock target

The retained receipt uses cap two; the executor accepts caps one through
three, matching the inherited comparison oracle, with independent controls
at caps one and three. The same geometric argument gives
`s_h=1-(3/5)^h` for every positive finite cap and the composition formula
above with `s_h` replacing `16/25`. The ideal capped output has total timeout
mass `delta_h=1-s_h^2`. Comparison with the uncapped completed-clock target
on an alphabet retaining timeout flags therefore adds `delta_h` to the
apparatus error. At cap two that allowance is `369/625`. The `1e-11` result
certifies implementation of the declared capped experiment. The analytic
timeout allowance tends to zero with the cap; execution at arbitrarily
large caps, spatial refinement and a physical continuum limit are separate
claims. Pulse costs stay bounded because there is only one intervention;
the maximum scheduled-repair cost is `2h`.

## Calibrated repair law and complete readout predictions

The [law-calibration bridge](law_calibration.py) connects the nine-quadrature
[repair identification](../repair_identification/README.md) to this same
continuing experiment. Its full-carrier current instrument restricts to the
invariant 66-state binary sector; forgetting private records is a CP
coarse-graining. The bound retains the joint clocks, detector, conditional
quantum state, every timeout and an arbitrary untouched reference. It applies
to the common-current interface, not coherent access to private Kraus labels.
It assumes the same relabelled instrument law on every seam. A calibration
of one seam alone does not establish that homogeneity.

Let `delta` bound the half-diamond difference of one repair instrument from
nominal C. For cap `h`, let `p=2/5` and `q=1-(3/5)^h`. Each waiting interval
uses `q/p` attempts on average. The first success additionally triggers one
feedback repair and a second wait. The complete ideal output therefore has
half trace-distance error at most

```text
min(1, mu_h * delta),   mu_h = (1+q)*q/p + q.
```

To prove this, telescope the chronological CP instruments, weighting each
local discrepancy by the probability that its call is reached. The scheduler
fixes these probabilities independently of the state and admissible law;
later identical controls contract the complete norm. The feedback repair
contributes `q`. No successful-branch conditioning is used. At cap two,
`mu_h=408/125`, improving the worst-branch multiplier five.

Nine probability tubes with C as a feasible witness constrain the entire
admissible family. If the tubes are centered on C with radius `r`, the
instrument bound is `delta <= min(1,7r)`. This uses a fixed nominal witness,
not a projection of noisy estimates onto the CP family. The strict complete
count interface retains aborts, separately declared sampling design,
apparatus errors and simultaneous coverage assumptions.

For the declared radius `r=1/10000`, the complete ideal-law error is
`357/156250`. The unconditional probability of both waits completing and the
terminal detector reporting one is `124967/750000` at `gamma=0` and
`35131/150000` at `gamma=1`. The nominal contrast is `1056/15625`; its
calibrated interval is strictly positive at the declared tube radius. These are conditional source
forecasts with no measured calibration counts. They are not a discrepancy
with Standard Model or general-relativistic observations. The
[current receipt](law_calibration_results.json) retains both settings and
their probability ranges without refitting.

The raw comparison bounds the executed nominal-C compiler against an
alternative calibrated ideal history, by the triangle through ideal C.
It does not execute the raw compiler with the alternative repair law.
Its independently certified compiled-array discrepancy is added separately.
Each replay recertifies its emitted arrays against the retained `1e-11`
envelope; platform-dependent exact residuals are not portable constants.
Raw maps need not preserve trace exactly: a half trace-norm
error `epsilon` allows event-weight error `2*epsilon`. Neither a probability
normalization nor a cap of one is applied to this raw-map error. Conditioning
on successful records requires a positive success-probability lower bound
and its own enlarged uncertainty. Physical calibration additionally needs a
validated detector/clock dictionary and independent observation data.

```bash
python -m model.compiled_clock_history.law_calibration --check model/compiled_clock_history/law_calibration_results.json
```

## Resources and remaining obligations

The ideal two-request protocol has at most four scheduled native attempts
and one feedback repair. Their expectations are respectively `328/125`
and `16/25`. Apparatus pulse expectations are

```text
preparation_pulses + 24*(16/25) + 3,
```

namely `1059/25` without a reference and `1134/25` for the compiled Bell
preparation. The maximum pulse counts are 51 and 54. These are logical
source pulse invocations, not elapsed physical time. Compiling or caching a
matrix once does not remove its repeated invocation from a branch budget.
The exact raw output also permits mass-weighted resource diagnostics;
those are distinguished from normalized ideal probability expectations.

The work register has dimension 512; the reference execution retains work
and a two-state reference. Exact rational computation carries complete
conditional quantum matrices rather than replacing them by their traces.
The terminal pointer uses one additional fresh qubit. A fixed 14-bit record
buffer holds two four-bit scheduled marks, one two-bit feedback mark, a
three-bit terminal status and one pointer bit. The feedback code uses the
three-outcome menu determined by the saved first mark and its feedback seam;
it is not a code for the union of both seam menus. A reusable two-bit counter
enforces each cap. Code-complement checks are declared services before each
native attempt, with at most five invocations. Their abort does not tick a
native clock. None of these storage or service counts specifies physical
duration, drift or device fidelity.
The retained producer and tests supply source bindings, emitted-word
counts, exact certificates, complete-output comparisons and hostile
controls. Numerical values are reported by the receipt generated from
those executions.

Blank preparation, coherent load/export, native C, scheduling, record
delivery, pointer readout and physical control accuracy are declared
services. The native repairs themselves are evaluated by their established
exact instruments; this attachment does not compile the entire scheduler
and every repair into source pulses. It proves neither a physical clock
scale nor an all-level spatial limit, Lorentzian dynamics, continuum field
theory or selection between admissible repair laws.
