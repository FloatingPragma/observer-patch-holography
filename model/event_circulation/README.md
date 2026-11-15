# Circulation in the observer's recorded event count

The [native event clock](../native_clock/README.md) and
[one-history circulation](../native_circulation/README.md) describe the same
native process. Joining them gives an exact prediction without assigning a
duration to a native attempt:

```text
E[S_N] = 0,                  Var(S_N) = N/10.
```

Here `N` counts nonzero native movements and `S_N` sums their signed currents
around the chosen triangular face. A movement away from that face has score
zero but advances `N`. These are event counts, not seconds. The result
assumes the event-stationary preparation described below, the binary total-two
carrier, uniform seam selection, fair coins and access to seam/current records.
It supplies no laboratory observation.

## The exact join

The original 66 configurations have 8 allowed moving branch addresses when the
occupied ports are adjacent, and 10 otherwise. Write this degree as `d(x)`.
There are 60 equally weighted seam/coin addresses per attempt. Every no-current
branch leaves the diagonal configuration unchanged, so summing all waits gives

```text
h(x) = d(x)/60,
Q_event(x,y) = Q_move(x,y)/h(x),
mu(x) = d(x)/600.
```

The executed native rule, not an assumed random walk, supplies every branch.
Exact enumeration verifies all marked transition rows, normalization,
irreducibility and detailed balance for `mu`. The adjacent configurations carry
total event-stationary mass `2/5`, and the others `3/5`. A geometric waiting sum
is legitimate because `h(x)>=2/15`; every next event occurs almost surely.
At a finite cap `M` the exact Palm-start timeout probability is

```text
(2/5)(13/15)^M + (3/5)(5/6)^M.
```

These tails are retained, not normalized away.

One actual native attempt from the stationary diagonal state `I66/66`,
conditioned on a moving outcome at that specified attempt, prepares `mu`.
The preparation succeeds with probability `5/33`; its failure probability
`28/33` remains part of the experiment. Starting from `I66/66` and waiting
until its first movement generally gives a different law. The receipt and tests
include this negative control. No hidden configuration measurement is required
to herald the event-stationary preparation.

Let `B` be the oriented source incidence matrix and
`P=I-B^T(BB^T)^+B` its cycle projection, reused from the existing circulation
packet. Conditional on a moving event at configuration `x`, direct enumeration
gives

```text
E[J | x, movement] = -B^T x/d(x),
E_mu[J J^T] = I30/30,         E_mu[(PJ)(PJ)^T] = P/30.
```

Thus `PJ` remains a bounded martingale difference in the complete event history.
Unequal-event covariances vanish without independence. The triangular weights
are in `ker B` and have squared norm 3, giving the exact variance `N/10`.
Gradient currents telescope to a bounded endpoint load change.

The finite, irreducible event chain satisfies the ergodic hypothesis used in
the existing bounded-martingale limit argument. In rescaled event count,
`N^(-1/2) sum_(k<floor(Nu)) PJ_k` converges to Brownian motion on the cycle
subspace with covariance `min(u,v) P/30`. On a fixed partition, its limiting
increment density has quadratic action

```text
15 sum_k ||delta y_k||^2/(u_k-u_(k-1)).
```

The coefficient 15 is the event-parameter version of the existing attempt
coefficient 99: the stationary conversion is `5/33` events per attempt.
This reparameterizes the same source-derived fluctuation law. Its relation to
the observer's modular flow and effective interacting field description is
not established; no universal physical clock or fundamental field is required.
Rank 19 counts cycle
coordinates, not spatial dimensions. The imported martingale limit theorem is
the one discussed in the parent packet; no Lean limit proof or finite-count
Gaussian error estimate is claimed.

## What this simplifies, and what it does not

Adding identity holding with state-dependent activity `0<a(x)<=1`
cancels from every marked event transition:

```text
a(x) Q_move,g(x,y) / (a(x) h(x)) = Q_event,g(x,y).
```

Consequently, all recorded event-counter laws from the same starting state
are unchanged. The receipt checks a nonuniform rational activity on all 66
states. This removes idle-attempt frequency from these observables. It does
not remove relative seam weights, preparation or record access. Attempt costs
and no-event records change. It also does not apply when another work system
evolves during the waits, as in the
[native-clocked common observer](../clocked_common_observer/README.md).

The nearest existing finite results are
[`GreenKubo.lean`](../../Lean/Thermodynamics/GreenKubo.lean),
[`NativeRepairFourCycle.lean`](../../Lean/Thermodynamics/NativeRepairFourCycle.lean)
and [`SourceRecurrenceClock.lean`](../../Lean/QFT/SourceRecurrenceClock.lean).
The parent packets map their transport and count hypotheses. This
packet adds the actual marked event-chain attachment, its Palm preparation and
event covariance; it does not substitute a different finite witness or a
two-state waiting theorem for the actual 66-state source.

## Execution and retained predictions

The shared primitive runtime starts from `I66/66` and executes the native
movement herald. Success continues with three native attempts; failure retains
its full terminal state and fills inactive slots without further native attempts.
The four continuation outcomes distinguish no movement
from movement with triangular score `-1`, `0` or `+1`. All 128 complete words
and full 66-dimensional terminal density matrices are compared with independent
exact native counts. A serialized checkpoint is continued and compared with
uninterrupted execution. First-event prefixes give nine stopped records plus
the timeout state and preparation failure. The receipt gives both conditional
and unconditional masses. These are extracted from actual executed prefixes; no
physical halt apparatus is implemented.

Conditional on preparation success, longer predictions use an exact integer recursion for the joint counter and
terminal configuration, retaining every score at 1, 4, 16, 64 and 256 events.
The common row denominator is 40 and the initial denominator is 600. The
receipt stores complete counter laws, exact moments and hashes of the full
joint integer arrays. No trajectories are sampled. Runtime floating comparisons
are diagnostics, while the finite event laws and moment identities use exact
rational arithmetic. The diagonal derivation is not asserted to identify the
full stopped quantum instrument on coherent inputs or arbitrary references.

At 256 events the variance is `128/5` and the exact finite law has excess
kurtosis about `0.0957799`. Successive squared increments have mean product
`11/500`, whereas independence would give `1/100`. The Gaussian limit therefore
does not replace the retained finite-count distribution with independent
Gaussian steps.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.event_circulation --check model/event_circulation/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/event_circulation
```
