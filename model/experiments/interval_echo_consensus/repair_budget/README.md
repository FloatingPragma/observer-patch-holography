# How many pair requests until the material settles?

With fresh uniform requests, a finite observer can predict the entire distribution of the number of successful repairs by a deadline. The prediction uses the initial vacancy and double counts and the number of ports. The angular phases of the material can change in flight without changing this forecast. The results apply to the [routed strict-write apparatus](../routed_counts/), with a protected blank routing tree, enough writable chords and a request drawn after the preceding free flight. One attempt includes a refusal. A finite risk budget refers to a number of requests. Physical clock time depends on how they are served.

## Exact observer probability

Let N be the number of source ports, v the number holding zero, d the number holding two, and C=N(N−1)/2. A routed request at a blank chord registers a write precisely for the two address orders 02 and 20. Summing the click projectors over all C unordered address pairs gives vd times the identity on the entire (v,d) sector. Each sharp string contributes exactly vd, and every projector is diagonal in that basis. The identity therefore holds for coherent superpositions, mixed states and entanglement with a retained reference. The full 99-state free seam conserves both load histograms. The commutator tests cover blank and written sectors. An arbitrary preceding flight, even one chosen in response to published results, leaves the success probability unchanged if the pair selection is fresh conditional on that history. These claims concern the material census. Address-resolved signs or other records can remain history dependent.

A write consumes one vacancy and one double. A quiet offer changes neither count. Put m=min(v,d) and a=|v−d|. At the stage with i pairs of defects left, a request succeeds with probability

\[
p_i=i(i+a)/C .
\]

Conditional on the enrolled initial census, the time to the m-th write is the sum of independent geometric waiting times with these parameters, each counting its successful request. Independence follows by multiplying the constant conditional success and refusal probabilities for every sequence of public results. For a finite deadline the probability of each possible number of writes follows by multiplying the finite transition matrix, retaining the completed state. For an initial mixture over inventories the observer averages the resulting laws with the preparation weights. Conditioning on a success alone would discard the requested timing information.

## Mean and finite budget

The exact mean is \(C\sum_{i=1}^m 1/[i(i+a)]\). For i≥2 the summand is bounded by 1/[(i−1)i], whose sum telescopes. The mean is strictly less than **N(N−1) requests**. On a 20-owner board with 240 ports and just one vacancy and one double, the mean is exactly **28,680 requests**. That last pair already requires order N² attempts in this uniform selection protocol. At most two arbitrary port transpositions carry the chosen pair to the chord. Each uses fewer than 2N adjacent tree swaps. Undoing the route leaves an O(N) upper bound on elementary routed swaps per attempt. These are finite construction costs for this service. No search over complete material histories or complexity bound for other repair programs enters the argument.

Write q=(1+a)/C for the last success probability. Separate the last geometric wait from the sum V of the preceding ones. Conditional on V, the survival probability at t is at most (1−q)^(t−V), also when t<V. For each other stage the geometric generating function gives E[(1−q)^(−G_i)]=p_i/(p_i−q). Thus, for m≥1 and 0<q<1,

\[
\Pr(T>t)\le \frac{m(a+2)}{m+a+1}(1-q)^t .
\]

The finite product telescopes because i(i+a)−(1+a)=(i−1)(i+a+1). At m=1 the tail is exact. Arbitrarily long quiet prefixes have positive probability at the last stage whenever q<1. Thus almost-sure completion has a quantitative risk budget, while a deterministic finite deadline requires a different schedule. At m=0 completion is immediate. The two-port q=1 case takes one request. [Lean arithmetic proofs](../../../../Lean/Screen/RoutedAttemptBudget.lean)

The same apparatus also admits a deterministic schedule under a different flight condition. In one protected scan offer each of the C pairs once, undoing each route before the next offer. Sharp surviving vacancies and doubles keep their addresses, and no write creates a new vacancy or double. If an opposite pair survived the scan, its own offer would have accepted it. This proves completion in at most C requests. The basis-column argument extends to the normalized output on coherent inputs. Free intervals can run before that scan, or between complete scans. The random-offer bound above permits a flight between every request, including one chosen from the earlier record. [Census theorem and native route](../routed_counts/)

## Native counts and controls

Use four physical ports forming a spanning tree and two spare native chords inside one cell. Prepare two doubles and two vacancies either sharply or coherently in two orders. Between interrogations evolve the load density matrix on the three tree seams, changing the duration with the stage and attempt. After a successful interrogation park its written chord. Independently construct every signed jump directly from the 99-state native echo columns. Sum the channel over the six actual pair addresses and over each retained sign outcome. This calculation supplies integer synthetic cohorts through the stated Born readout. The forecast uses only the initial census and number of ports. [Protocol](PROTOCOL.md) · [Program](measure.py)

At five requests the coherent run recorded **1,032 / 139,252 / 121,860** starts with zero, one and two writes, respectively, out of 262,144. All 16 native count banks agree with the independently calculated distribution within a simultaneous 99% half-width of 0.00449. An orientation-symmetric weakening of the strict write gives 47,699 completions at the same appointment. Parking the material and requesting the same pair repeatedly gives zero completions. All three controls, including the disabled write, disagree with the native forecast at all seven nonzero appointments for both preparations. The large-board entries in the receipt are analytic predictions. [All counts and forecasts](counts.json) · [Independent tests](test_budget.py)

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.repair_budget.measure
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/repair_budget
(cd Lean && lake build RoutedAttemptBudget)
```
