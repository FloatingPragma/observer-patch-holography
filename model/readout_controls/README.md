# Readout controls and a prospective native-law test

This packet separates two questions: which deposited residuals cannot be
removed by changing the oscillator action, and which fresh measurement
could distinguish the native C and E channels without fitting either
candidate. It executes both readouts, retains the separate empirical
scores, and provides no OPH-specific laboratory confirmation.

## An action-independent empirical comparison

The [McGarry experiment](https://arxiv.org/html/2603.04744v1), Section II.2.3,
identifies its binary spin means with real and imaginary characteristic
components. The [public deposit](https://zenodo.org/records/18752287) supplies
the aggregate measurements. This comparison uses authenticated local
copies of its two one-dimensional scans and their
[provenance](../anharmonic_observations/provenance.json). Appendix E describes
the characteristic scans and their nominal sampling. Accepted trial counts
and individual outcomes are not present in these files.

Every normalized state satisfies

```text
chi(0) = Tr(rho I) = 1,
chi(−beta) = conjugate(chi(beta)).
```

Consequently the ideal real-origin detector gives probabilities `(1,0)`;
the imaginary-origin detector gives `(1/2,1/2)`. At opposite imaginary
coordinates its means sum to zero. Changing the initial oscillator state,
anharmonic action, dephasing or elapsed duration cannot change these
identities. Postselection changes the state but preserves them whenever
the selected state is normalized and both probe settings address that state.
They are standard measurement identities, not distinctive predictions of OPH.

The producer freshly executes four native-compiled oscillator detectors:
real and imaginary origin, and imaginary reads at `beta=±0.25i`. Four
input controls include a nonzero signed response and complex coherence at
the occupation edge. Both pointer outcomes are retained. Exact rational
arithmetic on all entries of the returned binary Kraus matrices bounds the
effect-sum and identity defects by `1e−10`. The arbitrary-state origin
identity itself is independent of the eight-level demonstration cutoff.

The inherited compiler uses its guarded polar representation. This packet
certifies the matrices it receives; it does not certify the original pulse
product's representation error or floating density propagation. To state a
normalized mathematical comparison, it bounds the canonical completion
`K_b (sum K_b* K_b)^(-1/2)`. That completion is not executed, and no runtime
state or raw source matrix is normalized or repaired here. Physical coherent
services and readout errors are separate obligations.

The empirical rule retains all 46 origin measurements and all 240 nominal
opposite-coordinate pairs, without fitting an offset, gain or action. All
values and residuals are in [the receipt](results.json).

| Fixed comparison | Result |
| --- | --- |
| Origin measurements, spin-mean RMSE | `0.05008124` |
| Real-origin opposite-outcome frequencies at 0, 2 and 4 ms | `0.002`, `0`, `0.006` |
| Minimum probability RMSE across the 707 positive-duration means, from origin identities alone | `0.00629129` |
| Additional nominal opposite-probe constraints included | `0.03419381` |

These are **descriptive sample residuals**, not estimates of a hardware-error
rate or statistical rejection. A fair imaginary-origin detector has sampling
noise even when it is ideal. The nonzero real-origin frequencies are
incompatible with a noiseless deterministic detector *if* these entries
are literal raw frequencies; they cannot identify which apparatus or
processing step differs. The comparison does not reconstruct raw counts from
printed values or assign unauthenticated standard errors.

The lower bound follows without solving a dynamics model. For two measured
means `a,b` and an ideal odd prediction `z,−z`,

```text
(a−z)^2 + (b+z)^2 = (a+b)^2/2 + 2[z−(a−b)/2]^2.
```

The producer sums this exact floor and the origin residuals, then divides
by the original 707-point count. Probability RMSE is half spin-mean RMSE.
The paired-coordinate floor assumes that the nominal opposite coordinates
really mean opposite apparatus settings for the same preparation. Their
printed decimal sums differ from zero by at most `1.3e−16`; no bound on
arbitrarily high-energy states is inferred from that small coordinate
roundoff. The origin-only floor does not need this pairing premise.

This explains a part of the residual structure that action changes cannot
remove under the declared readout model. It does not explain the remaining
trajectory mismatch, prove that the dataset is unreliable, or justify
subtracting its residuals from the dynamical-model score. A contemporaneous null
and phase-reversal calibration, retaining accepted and rejected counts,
is a concrete next laboratory requirement.

## A complete native discrimination experiment

The second calculation gives a prospective test of the two defined
native candidate channels. C retains the coin and repair current; E
additionally resolves the pre-repair endpoint contrast. Their classical
configuration laws agree, while their preserved coherences differ. The
[native interface](../native_interface/README.md) and
[compiled mixer](../compiled_native_interface/README.md) supply the
preparation/readout construction.

Start with the declared stationary `I66/66`. A native fixed-seam current
heralds either of two orientation preparations. Use the same source-compiled
balanced mixer to prepare and analyze the phase, with four unread native
attempts between them. The final native coin/current record determines
whether the phase matched its preparation. The parameter is the count of
native attempts, not a duration in seconds. There is no fitted time scale.

The producer checks the actual source coherence blocks and freshly executes
the complete six-outcome runtime, including serialized checkpoint
continuation. The exact predictions are

| Quantity | C | E |
| --- | --- | --- |
| Useful-herald probability | `5/33` | `5/33` |
| Conditional matching probability | `1921/2592` | `12401/20000` |
| Failed-herald probability | `28/33` | `28/33` |

Each herald sign contributes half the useful mass. The impossible
`(failed,1)` record remains present with zero probability. The total
variation distance of the two complete distributions is
`12259/668250`, about `0.018345`. These probabilities were derived in the
interface construction; the result is a fixed-trial discrimination
rule and certified sample budget that counts every failed trial.

Let `P,Q` denote these six probabilities and
`A=sum_y sqrt(P(y)Q(y))`. For `N` independent fresh trials, select C when
its likelihood is at least E's. Under either candidate, Markov's inequality
on the square root likelihood ratio gives individual decision error at
most `A^N`. The producer encloses every square root and logarithm by exact
rational arithmetic. **3,579 total trials suffice for an error bound of
1% under either locked ideal candidate.** There is no conditioning on a
favorable success count or optional stopping. The future count classifier
accepts all six counts and rejects an outcome impossible under both models.

The short source mixer is freshly executed on its complete 512-dimensional
lift, with all padding retained. Its two-call comparison bound is at most
`1e−10`, including norm amplification of the original binary maps. The
six-record runtime is the bounded 66-dimensional ideal-target surrogate;
it is not a projection of the raw source matrix. The comparison extends
native operations as CPTP maps on the padding and retains a failed record
there. Raw rounded source matrices are not assumed exactly trace preserving.
The statistical design applies to the exact ideal candidate laws, while
source representation diagnostics and physical implementation errors are separate. For independent physical trials, a trace-preserving trial law within
total variation `epsilon` of its candidate gives the conservative bound
`0.01 + 3579*epsilon`. The same bound holds without independence if every
actual conditional trial law, given its past, is uniformly within `epsilon`
of the ideal candidate law. Marginal per-trial closeness alone is insufficient.
No such laboratory `epsilon` is supplied here. The
[continuing-carrier protocol](../closed_prediction_protocol/README.md) gives
finite preparation and sequential conditional-law accounting.

The coherent proper-code service, stationary preparation and uniform seam
schedule are declared apparatus premises. A working simulation of this
protocol would distinguish C from E for that apparatus; it would not by
itself show which law describes an unconstructed physical observer. No
public dataset with this native protocol is identified here. Matching a
different oscillator experiment by inventing a detector or clock mapping
would not close that gap.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.readout_controls --check model/readout_controls/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/readout_controls
```

The packet runs offline. Source files, protocol, original data and existing
attribution are hash-bound. Exact data fractions, candidate laws, trial
counts and bounds must replay exactly; executable floating diagnostics use
the shared numerical tolerance. The original data are attributed to
McGarry and coauthors under CC BY 4.0, as recorded in their existing
provenance. Existing empirical producers and receipts are unchanged.
