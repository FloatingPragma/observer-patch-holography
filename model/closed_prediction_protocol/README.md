# A finite sequential native-law prediction protocol

This protocol uses **finite unread native preparation on the same continuing
carrier** and bounds the resulting dependent record sequence. It combines the
mixing theorem, native phase interface and source mixer with complete
preparation, memory and decision-error accounting. Its
[ideal product-law comparison](../readout_controls/README.md) defines the
reference likelihoods; the actual trials need neither exact reset states
nor independence.

The locked design uses **379 discarded native attempts before each of
3,927 trials**, with three unread attempts inside each useful phase probe.
The exact-source model's probability of choosing the wrong C/E candidate is
at most **0.75%**. A separate, explicitly unauthenticated physical service
allowance could raise this to **1%**. No elapsed physical time, empirical
native-law identification or relativistic spacetime recovery is claimed.

## What runs from end to end

The supplied architecture is one twelve-port carrier in its binary
total-two sector, with 66 configurations. The supplied native law is C or
E. Preparation uses independent uniform choices among its thirty seams and
a fair placement coin. This schedule is an experimental input; it does not
follow from the production observer's adaptive scheduler.

1. Apply 379 actual unread native channels to whatever sector state the
   carrier has. Forget their preparation records. Retaining and
   conditioning on those records would define a different preparation.
2. Perform the existing fixed-seam native current herald. Retain both useful
   signs and failure. The useful reference yield is `5/33`.
3. On either useful branch, execute the declared source balanced mixer,
   three unread native attempts, a second mixer, and the native coin/current
   pointer. Interpret the pointer using the herald sign. On failure, retain
   the failed record and continue to the next trial.
4. Keep the resulting carrier state. Start the next trial with its own 379
   discarded attempts. Do not replace that state with `I66/66`.
5. After exactly 3,927 trials, compare the two locked likelihoods using all
   six counts. Failed heralds count as trials; there is no success-count
   stopping rule. An impossible record returns `outside_both_declared_models`.

The producer executes two consecutive trials for each law, retaining all
36 joint records and their full unnormalized 66-dimensional terminal states
in memory. The receipt retains every joint probability and complete-state
audit residual. It starts from a nonstationary complex superposition; the
second trial uses each actual first-trial terminal branch. It also executes
preparation of a carrier entangled with a two-level untouched reference.
These finite executions check the construction; they are not 3,927 sampled
experimental trials, and their computed probabilities are not observations.

Preparation enumerates all native partial-injection matrix-unit paths and
adds their exact integer multiplicities divided by 60. Every discarded
attempt is applied separately. This is a sparse execution of the source
Kraus sum, with no spectral exponential or finite-time reset. A dense
complex-reference check independently compares it with the shared runtime's
complete native Kraus application. Preparation checkpoints are serialized
and replay-validated; the herald checkpoint also uses the shared runtime's
independent full-prefix validation.

The coherent pulse executor separately runs the source mixer word and lifts
it to all 512 encoded addresses, including padding. The six-record runtime is the existing
66-dimensional ideal-target surrogate. The freshly recomputed full-source
two-mixer comparison is bounded by `1e-10`; raw source matrices are not
projected into a supposedly actual 66-dimensional gate. That comparison is
an exact bound on the returned floating arrays against the target, including
norm amplification. It is neither a physical gate-fidelity measurement nor
a bound on all floating state propagation. The numerical runtime residuals
are diagnostics, separate from the analytic decision guarantee.

## Why repeated trials need no exact reset or independence premise

The existing [native equilibration bound](../native_equilibration/README.md)
is uniform over arbitrary input states and untouched references. Within
this fixed sector,

```text
D((Phi^n tensor id_R)(rho), (I66/66) tensor rho_R)
    <= 33 (229/240)^n.
```

The gap comparison is exact: `(9/4)^2 > 5` gives
`(5-sqrt(5))/60 >= 11/240`. Exact rational powers verify that 379 attempts
give distance at most `epsilon = 1/1570800`. This is the smallest attempt
count certified by that particular conservative gap bound at this
allowance, not the measured minimum mixing time.

Condition on any previous visible history of nonzero probability. Its
incoming carrier state may depend arbitrarily on that history, and may be
entangled with a passive reference. Apply the same unread preparation and
then the same complete instrument. Trace-distance contraction gives

```text
TV(actual next-trial law given the past, ideal candidate law) <= epsilon.
```

The previous record can affect the incoming state; it does not select an
undeclared replacement instrument. The reference remains untouched during
the trial. An environment that feeds uncontrolled information back into the
carrier during its preparation or measurement is an additional service
error, not a case silently covered by this statement.

For a sequence of length `N`, telescope the conditional probability kernels
one at a time. Each replacement contributes at most `epsilon` in total
variation, so the actual transcript is within `N epsilon` of the ideal
product law. This argument does **not** assume the actual trials are IID.
The untouched-reference version also permits intertrial quantum correlations
before the next preparation. No finite equilibration is asserted to erase
the old classical record or to make the actual states exactly independent.

The ideal candidate laws have Hellinger affinity
`A = sum_y sqrt(P_C(y) P_E(y))`. The fixed-sample likelihood test has either
individual error at most `A^N`; exact outward square-root and logarithm
bounds give `A^3927 <= 1/200`. The inherited classifier counts all failures.
Combining the two existing ingredients therefore gives

```text
wrong-candidate probability <= 1/200 + 3927/1570800 = 3/400.
```

If each complete physical trial additionally differs from its declared
model by conditional total variation at most `eta`, uniformly over its
past and relevant references, add `3927 eta`. Meeting
`eta <= 1/1570800` suffices for the stated 1% target. This is a required
physical tolerance, **not an achieved apparatus certification**. It must
cover relevant preparation, control, native-law, storage and readout errors;
the compiler's floating-array certificate does not supply it.

## The locked predictions and their resource cost

The design rule was fixed before any future observations: search unread
probe counts one through eight and minimize the worst-case total native
attempt count, including every preparation and failed trial. It selects
three unread attempts. This is optimal only within that declared finite
menu, not among all possible observer experiments. The separate four-attempt ideal product-law design requires 3,579 trials;
it has a different error budget and does not include this preparation cost.

The useful-branch matching probabilities at three attempts are
`341/432` for C and `1343/2000` for E. The failed-herald probability is
`28/33` under either ideal law. Each sign has half of the useful mass; the
complete six-record laws differ in total variation by `1591/89100`.
These are finite conditional model predictions, with all failures retained.

A failed trial uses 380 native attempts. A useful trial uses 384 native
attempts and two 34-pulse mixers. The worst-case full design therefore needs
**1,507,968 native attempts**, **267,036 source pulses**, and **534,072
proper-code load/export events**. These count different operations and do
not add up to a physical duration. The ideal product-reference expected
native cost is 1,494,640; the actual dependent experiment need not have
exactly that expectation. The certificate executes only the two-trial
checks described above, not this entire sampled campaign.

## Derived consequences and interface inputs

Exact stationary preparation and exact trial independence are not
premises of this finite test. The native repair law supplies a quantitative
preparation operation, and the continuing-state protocol has a complete
sequential error bound. The source compiler supplies the chosen mixer once
its coherent proper-code services are admitted; committed native repairs
alone cannot create the needed coherence from diagonal inputs.

The fixed particle sector, native law candidate, random schedule, coherent
transport/control services and classical record implementation are
declared architecture/interface inputs. They are not claimed to follow from
the closure hypothesis, the no-outside rule or the agreement condition of the
[OPH basis](../../docs/BASIS.md). No public laboratory dataset implementing this native
protocol is identified here. A meaningful empirical comparison requires an
apparatus mapping and an authenticated complete-trial error bound; fitting
an unrelated oscillator trace does not provide those inputs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.closed_prediction_protocol --check model/closed_prediction_protocol/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/closed_prediction_protocol
```

The packet is offline, source-bound and reproducible. Exact laws, attempt
counts and statistical inequalities replay exactly; floating execution
diagnostics use the shared numerical tolerance.
