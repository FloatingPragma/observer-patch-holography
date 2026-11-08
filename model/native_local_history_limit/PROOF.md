# Complete local repair records at growing capacity

The result concerns the fair binary endpoint law, its independent uniform
seam selector and an explicitly declared local observation. It controls
complete visible C/E records and their quantum disturbance. The population
comparison in [SCALING.md](../federated_observer/SCALING.md) sums out the
native coin, and the [exchange-clock comparison](../native_exchange_clock/PROOF.md)
controls endpoint correlation operators. Neither is a complete-instrument
statement. The present proof uses the actual endpoint attachment and
backward-cone identity of [native_local_response](../native_local_response/PROOF.md).

## Source and observation

Let a finite simple graph have `V` binary ports, `M` seams and maximum
degree `Delta`. Each of `T` attempts independently selects one seam
uniformly, then applies its complete fair native instrument. A fixed C/E
mixture, including its tag, fair coin, impossible raw slots and current,
is allowed. Unobserved outcomes are summed, never conditioned away.
The graph family, conserved occupation and state are declared inputs.

Choose observed seams and final ports whose endpoints lie in a nonempty
fixed set `S`. Retain every raw record on those seams, its order and a
declared time-bin label, together with the final state on the chosen ports
and any external reference `R`. Trace all other record values and exact
global tick addresses. The native archive remains readable; this theorem
does not erase it or forbid an observer from using a richer observation.

The clock partition is `0=t_0<...<t_b=t`, with `M*t_j` integers, and
`T=M*t`. An event at tick `k` has bin label `j` when
`M*t_(j-1)<k<=M*t_j`. This clock is the per-seam **attempt** scale.
The population exchange time of the companion package is `t/2`.
Neither parameter is identified with physical proper time. A common
partition is required across capacities. The evaluator rejects unaligned
boundaries rather than silently rounding them.

The native memory convention is export-oldest, so a finite FIFO does not
halt the experiment. A host resource stop remains a distinct outcome.
The theorem concerns the ideal finite kernel; admission of a program and
expansion of its complete branch tree are separate checks. Record-adaptive
global selectors are outside the independent-selector hypothesis.

## 1. Localize the full quantum instrument

Let `B` contain the graph ball of radius `r` about `S`, and let `K` be the
number of seams internal to `B`. Replace every noninternal attempt by the
identity and discard its record. Keep the global attempt count and the
original internal seam maps. For a fixed selector word, the full and
window instruments agree whenever the backward dependency cone stays in
`B`. Starting the cone with all of `S` at the terminal time contains every
actual record root at its earlier event time. This also covers the random
set of observed indices: it is determined by the input-independent selector
word, before the quantum outcomes are sampled.

Escape requires a time-ordered path of `r+1` seam selections. For each
starting vertex there are at most `Delta^(r+1)` such paths, and each
specified sequence of distinct time indices has probability `M^(-r-1)`.
Thus the complete observation channels satisfy

```
(1/2)||Phi_native - Phi_window||_diamond
 <= a_N := min(1, |S| binomial(T,r+1) (Delta/M)^(r+1))
 <= min(1, |S| (Delta*t)^(r+1)/(r+1)!).
```

Conditional maps are trace preserving when all their raw outcomes are
retained. On selector words where escape occurs, their half-diamond
distance is at most one; elsewhere it is zero. Averaging proves the bound
for arbitrary initial coherences and an arbitrary entangled reference.
This is a full-instrument bound, not a hidden-swap-bit identification.

## 2. A marked Poisson comparison with the same local records

On the finite window put independent rate-one Poisson **attempt** clocks
on its `K` internal seams. At each marked event apply the same complete
native endpoint instrument. This defines a comparison process; it is not
an added native scheduler or a unitary evolution law.

In a native tick interval of length `1/M`, the window selector is empty
with probability `1-lambda`, or selects each internal seam with probability
`1/M`, where `lambda=K/M`. The Poisson selector has empty probability
`exp(-lambda)`, exactly one event on each seam with probability
`exp(-lambda)/M`, and all multiple-event words. The common probability mass
is exactly

```
1-lambda + lambda*exp(-lambda).
```

Consequently the one-tick total variation is
`lambda*(1-exp(-lambda)) <= lambda^2`. Independently couple each of the
`T` tick intervals. When their marked words match, use the identical
ordered native quantum instrument, including all visible coins and raw
values. Trace distance contracts under this common instrument, even with
an arbitrary reference. This gives

```
(1/2)||Phi_window - Phi_Poisson,B||_diamond
 <= d_N := min(1, T*K^2/M^2) = min(1,t*K^2/M).
```

On a matching singleton, its continuous event time and the native tick
have the same time-bin label because every boundary lies on the grid.
Within-bin order is retained. Poisson intervals with two or more events
are included in the failure allowance; they are not deleted or
renormalized. The proof also permits fixed common local interventions at
bin boundaries when their supports are included in `S`.

Exact timestamps cannot be substituted: a native event time lies on a
countable grid, whereas a Poisson event time has a continuous law. Their
conditional timestamp laws have total variation one. Complete global
record lengths also distinguish the deterministic and Poisson schedulers.
The observation map is therefore an essential premise, not a cosmetic
change in units.

## 3. State restriction and the capacity bound

If the initial joint `B,R` states have half-trace distance at most `epsilon`,
contractivity and the two preceding comparisons give

```
D(output_native, output_common_Poisson_window) <= epsilon + a_N + d_N,
```

capped at one. No factor depending on the number of output records is
needed: the bound controls the entire classical-quantum output.

For the fair fixed-occupation state `tau_(V,q)`, the exact canonical
restriction theorem gives distance at most `v(v-1)/(2V)` from the
`v=|B|`-port Bernoulli product with density `q/V`. Couple each product bit
to one with a declared limiting density `p`. A union bound adds
`v*|q/V-p|`. Therefore

```
epsilon <= min(1, v(v-1)/(2V) + v*|q/V-p|).
```

This specialization permits an independent common reference. A correlated
reference requires the actual joint `B,R` bound; agreement of the system
marginal alone is insufficient. The native law conserves `q` and does not
select its dependence on capacity or the density `p`.

For two source graphs with isomorphic rooted windows, transport all seam
addresses, endpoint orientations and observed records by that isomorphism.
Their distance is at most the sum of their two `epsilon+a_N+d_N` bounds.
An endpoint reversal also requires the native coin/current relabeling;
the executable controls use orientation-preserving maps instead.

## 4. Local process limit and its hypotheses

Suppose `M_N,V_N` grow, rooted windows eventually agree at every fixed
radius, the endpoint law and observation agree under their identifications,
and initial joint restrictions converge consistently. For each fixed
window, `d_N` tends to zero. The continuous backward-cone estimate is

```
a_infinity(r,t) <= min(1, |S| (Delta*t)^(r+1)/(r+1)!).
```

It follows either by the factorial moment measure of the independent
marked clocks, or by the same ordered-path argument with ordered time
integrals of volume `t^(r+1)/(r+1)!`. The bound tends to zero as `r` grows.
It constructs the local limiting instrument by a Cauchy comparison of
finite windows and controls all finite-time local raw records and their
backaction. For the canonical specialization, the distance to a matching
infinite local process is bounded by `epsilon+a_N+d_N+a_infinity`.

Bounded degree and a finite alphabet give finitely many rooted radius-`r`
graphs. A diagonal subsequence therefore has locally stable rooted balls
after declared address relabelings. This is a subsequence compactness
statement, not a proof that the selected source refinement has a unique
rooted limit. The fixed Bernoulli-density initial family is consistent on
such balls. Arbitrary initial/reference families require their own
compatibility or convergent subsequence.

The actual evaluated fair federations have `V=240*4^level`,
`M=630*4^level` and degree at most six. Their cell-zero radius-one windows
at levels zero, one and two are explicitly isomorphic: twelve original
ports, three attached halo ports, and thirty-three internal seams. The
attachment ports are `(1,3,10)`. The executable checks their edge sets,
endpoint orientation and transported addresses, then composes the state,
selector and exact four-state halo-escape bounds on one common clock
partition. This finite check supplies no all-radius convergence theorem
for the floating routed source constructor.

## 5. Finite evaluation without silent postselection

The finite-window Poisson event count has mean `K*t`. Truncating its
enumeration after `J` events requires a separate retained overflow outcome.
Its mass is bounded by `min(1,(K*t)^(J+1)/(J+1)!)`, since the falling
factorial of order `J+1` is at least `(J+1)!` above that cutoff. Any chosen
trace-preserving replacement on that overflow branch changes the output
by at most this allowance. This is an error bound, not successful native
execution of infinitely many branches.

The tests enumerate complete raw maps on small endpoint graphs, retain
an entangled reference, check that within-bin chronology matters, and
verify exact selector and canonical-state controls. The owning source
attachment checks every H66 native seam map and the four federation
endpoint patterns. Arbitrary-size rational bounds use integer census data;
they do not execute a cosmological carrier.

This theorem makes a fixed local experiment controllable as its declared
surrounding source grows, under explicit state and geometric hypotheses.
It selects no physical capacity, unique density, three-dimensional space,
Lorentzian causal cone, interacting chiral action, particle or measured
mass. All eight physical source-to-spectrum requirements retain their
separate identification obligations.
