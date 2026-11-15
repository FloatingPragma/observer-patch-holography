# A native phase measurement on the stopped interacting work history

The [clocked common observer](../clocked_common_observer/README.md) retained
native quantum states but read only the work register. Its C and E candidates
therefore had identical work-record probabilities. This packet measures the
retained native phase on **that same history**, and predicts its joint law
with all four work records. It also replaces independently supplied coherent
preparation and phase effects by the existing native herald and coin/current
detector. The balanced operation they share has a
[34-pulse implementation](../compiled_native_interface/README.md) in the
existing coherent source compiler.

The useful preparation probability is `5/33`. All `28/33` failed preparations,
both successful herald signs, the stopping slot and current sign, all sixteen
work records, the final native phase bit, and both unresolved timeouts remain
in the output. There is no laboratory dataset or empirical-agreement claim.

## What actually executes

Prepare the native diagonal state `I/66` and the finite sixteen-dimensional
work vacuum independently, once each. Execute one fixed-seam native repair.
A zero current ends this trial as a failed herald; the prepared work vacuum
is retained. Either nonzero current heralds one of the two orientation
mixtures. Apply the balanced mixer to that **unnormalized** native branch.
The herald sign determines whether it prepared plus or minus phase.

For each succeeding controller slot, execute the native signed-current
instrument and one work waiting evolution of supplied model coordinate
`delta=1/4`. On the first nonzero current at slot `n`, latch the native
register, and let the evolved work execute its current-controlled
kick, four-read continuation and record-conditioned feedback. It consumes
the state after `n` waiting steps, without a new work preparation.

Apply the same mixer to the retained native register and execute the
fixed-seam native pointer. Its raw plus bit is determined by
`2*current+2*coin-1>0`. Reverse the interpreted bit for the minus herald, so
zero means “matches the prepared phase” for either sign. This native phase
operation commutes with the work suffix at fixed stopping slot and sign;
evaluating it first in the factorized simulation represents the same joint
instrument. Neither subsystem's later records control the other in this
packet.

If no current occurs before the cap, retain both continuing states without
performing a phase measurement. Measuring that timeout would disturb the
state needed for cap extension. The phase sentinel `-1` means **no phase
operation occurred**, not a third pointer outcome. The complete record is
`(herald, slot, current, work_bits, phase_bit)`, with one failed-herald record
and two unmeasured timeouts. A cap `N` has `128*N+3` output records, including
zero-probability outcomes.

The runtime executes actual native Kraus maps and actual source-compiled
work maps. It does not receive a normalized successful phase state or an
independently supplied phase POVM. Native and work terminal matrices remain
as tensor factors for every record; forming their large tensor products is
unnecessary. This is an exact representation of the declared independent
preparation and classical control, not physical copying or independent
experiments at each possible stopping slot.

## The precise source-implementation boundary

The finite native execution uses the numerically evaluated **ideal balanced
service** in its 66-dimensional invariant sector. Its original binary output
arrays are certified directly against exact rational targets. Separately,
the producer executes the complete source mixer on all **512 code addresses**
and attaches its directed operator/channel certificate. It does not compress
the source matrix, remove its small padding leakage, or replace its original
entries by a unitary projection.

The exact full ideal mixer has the same two-coordinate action at each of
64 spectator addresses. It preserves the native code, agrees there with the
66-dimensional mixer, and also acts on padding pairs. An explicit CPTP
extension of every intervening native operation first dephases the native
code versus its complement, applies the original native instrument on the
code, and places identity evolution of the complement in its zero-current
or bit-zero outcome. This extension is an additional operational convention;
it changes no native-sector target. The work operations have their own
complete finite CPTP target.

If the actual source matrix has operator error `epsilon` against the exact
ideal mixer, its half-diamond error is at most
`d=epsilon*(2+epsilon)/2`, and its CP trace-norm gain is at most
`g=(1+epsilon)^2`. At most two controlled uses, separated by exact CPTP
processing, cost **`d*(1+g)`**, including possible numerical norm amplification.
That full-carrier comparison also covers an untouched work/reference system.
It is charged without dividing by the successful-herald probability.

Thus the receipt contains two distinct bounds: the actual factorized binary
66-dimensional execution against the exact ideal-service experiment, and its
comparison to exact applications of the full binary source mixer with the
other ideal operations. A floating 512-dimensional history is not executed
or certified here. Binary source maps are not assumed exactly trace preserving;
their complete record weights are not renormalized. The compiled service
removes an independently supplied mixer **implementation**, while coherent
source services and the balanced target design are premises.

## Complete numerical certificate and finite theorem attachments

For every native outcome, exact rational arithmetic evaluates the declared
herald, first-current subchannels, mixer and native pointer. All reached
states are block diagonal between the active phase corner and its complement,
so both mixer factors cancel the square roots. A cross-corner coherence is
refused by this exact domain rather than discarded. Each original reported
66-by-66 matrix is compared entry by entry, including off-support entries and
non-Hermitian arithmetic defects. A directed Frobenius bound gives its half
nuclear-norm error `e_r` and reported norm at most `p_r+2*e_r`.

The existing all-input work certificates supply `w_r` for the complete work
menu at each stopping slot and current sign. Failure uses the once-prepared
vacuum; timeout uses the retained waiting state. The complete joint bound is

```text
sum_r e_r + sum_r (p_r + 2*e_r)*w_r.
```

No additional factor of sixteen is needed: each `w_r` covers the
complete work instrument. The preparation error, repeated numerical norm
gain, rounded waiting products and every returned work density were charged
by that certificate. The computed work budget is bound to the executed
histories; its public digest is not authentication of a caller-fabricated
proof object.

The exact unresolved population is

```text
(5/33)*[(2/5)*(13/15)^N + (3/5)*(5/6)^N].
```

It is added once for comparison with unlimited stopping in the same finite
work model. Cap extension means continuing both retained timeout states;
forgetting later records cannot undo later quantum operations. Effect-weight
intervals use the full nuclear error `2*e`, since the binary outputs need not
have exactly unit trace.

The [finite Kraus composition theorem](../../Lean/Dynamics/ChoiCPTP.lean)
applies to the ideal complete native and work maps, including untouched
references; the record direct sum can be embedded in a sufficiently large
finite tape. The native
[record-count identities](../../Lean/QFT/SourceRecurrenceClock.lean) apply to
the first-current counter, without importing that file's separate two-state
return law or a physical clock calibration. The
[native-interface reachability argument](../native_interface/README.md)
prohibits constructing coherence from diagonal inputs using the
committed partial injections alone. The coherent compiler services are the
additional resource. No further Lean theorem is claimed.

## Predictions and their meaning

At cap 64, the complete output contains **8,195 records**. The exact
unconditional probability of a successful stop followed by a matching phase
bit is approximately **0.0980095163 for C** and **0.0900157544 for E**. Their
difference is approximately **0.0079937619**. These include the failed heralds
and unresolved timeouts; they are not conditional probabilities promoted to
full-trial success rates.

The phase-match event attains the total variation distance of the complete
ideal classical record laws for this construction. At fixed herald, slot and
current sign, both candidates use exactly the same normalized conditional
work-record distribution. Appending those work records therefore leaves
that full classical distance unchanged. Work-only marginals are equal.
The gain is the connection to one continuing interacting history and its
joint native/work observables, not greater discrimination than the
retained native record alone allows.

The [receipt](results.json) additionally gives bounded probabilities for
phase match together with final gauge bit one, and for phase match together
with unequal smooth-read bits. It reports C-minus-E intervals and the
unnormalized correlation witness
`P(match,gauge)*P(success)-P(match)*P(gauge)`. This equals the squared success
probability times the conditional covariance, without discarding any trials.
The source-comparison half-L1 bound is labelled as a distance between
unnormalized binary record weights, rather than claiming exact normalization.

At cap 64, the complete ideal-service output error is at most `0.000003`.
The failed-herald probability is `28/33`, and the unresolved mass is
approximately `0.00000716142`. The full source-mixer comparison has the same
displayed finite error bound; its two-call contribution is below `1e-12`.
The C-minus-E record-weight difference for **phase match and final gauge one**
is in **`[0.005151, 0.005175]`**; for phase match and unequal smooth bits it
is in **`[0.002564, 0.002587]`**. Both compare the complete finite record laws;
these are weights of unnormalized binary record maps, not normalized
probabilities.
The covariance numerator is in `[-0.000035,-0.000030]` for C and
`[-0.000031,-0.000026]` for E, establishing phase/work correlation after
conditioning on a successful stop, with that conditioning's mass explicit.

Only the ideal invariant-sector process has the stated vanishing-tail
unlimited-stopping comparison. Under the specified padding extension, binary
source leakage can remain in a zero-current sector indefinitely. Its rotating
work state need not have a terminal-state limit. The source-mixer comparison
and its complete record-weight intervals therefore concern finite caps; the
receipt does not promote the ideal tail theorem to an unlimited binary-source
instrument.

Successful paths use `n+2` native attempts, two mixer calls, `n` work waits
and one complete work suffix. Timeouts use `N+1` native attempts and one mixer;
failed heralds use one native attempt and no mixer or wait. All branches pay
for one work preparation. The receipt keeps expected and maximum work-plus-
mixer source pulse budgets separate from the native channel attempts, which
are not compiled into source words here. Expected costs use the exact ideal
invariant-sector clock law; they are not asserted as expectations of the
slightly non-normalized binary-source comparator. The finite maximum path
cost is valid with the declared padding extension.

The remaining premises are the binary sector and diagonal initial state,
the chosen seam and uniform scheduler, the coherent compiler services,
the effective work Hamiltonian and finite occupation cutoff, classical
controllers, and work model coordinate per slot. Slots are native attempt
counts after the herald. The existing
[finite thermal-clock result](../common_source/THERMAL_CLOCK.md) aligns
unitary work evolution with observer-relative modular flow through
`s=-t/beta` for a supplied faithful Gibbs reference. This packet does not
establish the relation between those modular parameters, the native clock
records and work evolution, or select the Gibbs reference. No fundamental
field or preferred global clock is required by that open task. Physical
source errors, physical memory/halting, a selected effective observation law,
occupation/spatial refinement and an empirical comparison are open.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.clocked_quantum_readout --check model/clocked_quantum_readout/results.json
python -m pytest -q -W error model/clocked_quantum_readout
```

Generation uses `--write` and refuses to overwrite a file. Replay freshly
executes the source maps, computes the exact comparisons and checks every
source binding and portable outward receipt value.
