# Can a random elapsed time explain the source response of a written clock?

An observer writes the order of two loads at a fixed seam. With an ordinary excess released one link away, the observer sometimes reads the opposite order on the same two ports. The two clock modes that store this order have identical energy on the filled reference screen. Waiting any amount of time on the reference preserves the order exactly. A positive distribution over reference-clock elapsed times therefore cannot account for all the source-conditioned clock records, even with an unbounded distribution or a correlation between elapsed time and the final source location. The result identifies an internal material response in this finite clock. [Protocol](PROTOCOL.md) · [Held replication](replication.json)

## What observers measure

The native write maps a sharp (2,0) or (0,2) at the addressed seam to two unit-filled ports with a definite contrast sign in a symmetric coin. The observer lets the excess move freely on the 240-port federation. At the appointment, the observer applies the same native echo and measures the two local loads. A (2,0) or (0,2) records the recovered order. When the excess occupies an endpoint, the read gives (2,1) or (1,2). Both unresolved events and all other outcomes enter the denominator. The preparation, release appointments, terminal pulse, sharp detector and Born sampling rule are supplied. A separate material dimer measures the exposure on its calibrated first branch. The assay reads only the local material ports. [Native column checks and pilot acquisition](measure.py) · [Independent sharp-basis acquisition](replicate.py)

Each entry below counts reversals in 1,073,741,824 independent simulated preparations per orientation and appointment. The finite generator gives the bin probabilities. Each multinomial draw samples the complete count bank without stepping through a billion trajectories. The independent bank uses a fresh seed and a 960-state generator assembled in sharp coin, sign and excess-position coordinates.

| Dimer wait | Reversals from (2,0) | Reversals from (0,2) |
| ---: | ---: | ---: |
| 1 | 157 | 147 |
| 2 | 455 | 393 |
| 4 | 666 | 642 |
| 8 | 527 | 531 |
| 12 | 220 | 243 |

At wait 4 the observer starting with (2,0) counted 1,073,416,782 returns, 666 reversals and 324,376 unresolved endpoint occupations. The simultaneous 99% exact binomial interval for reversal was 0.525 to 0.727 per million. The blank upper bound was 0.0105 per million. Both starting orders cleared the exact-binomial comparison at every held wait. The distribution-free Hoeffding bands were too wide to resolve any pair. The blank and parked one-link controls recorded zero reversals, as did the changed law with contrast-two sign exchange removed. The latter preserves ordinary excess hopping and clock sensitivity to endpoint occupation. The dimer intervals cover each programmed duration on the calibrated branch. [Complete counts and intervals](replication.json)

The pilot used 1,048,576 launches per setting. It recorded zero reversals at all five held one-link waits in both orientations. The endpoint launch resolved a reversal against blank at six of the eight appointments in each orientation under simultaneous exact binomial intervals. At wait 4 the counts were 958 and 964. Two-link launches gave up to six reversals per bank and never separated from the blank bound. Neither distance three nor fourteen resolved a reversal. Parking the source on the clock gave an unresolved endpoint read on every launch. The independent replication increased the one-link count after the pilot's sensitivity was known and held the wait menu fixed. [Pilot counts](counts.npz) · [Decisions by setting](results.json)

## Why an arbitrary time distribution has an exact null

In the symmetric coin, let the two clock states be the even and odd superpositions of the signed contrast-two writes. On the filled reference both energies are exactly zero. Either prepared sign lies entirely in this degenerate subspace. For each elapsed time τ, applying the blank clock evolution to either sign returns it. Averaging over any positive distribution of τ preserves both signs. Even if a trial's elapsed time specifies where its excess lands, the final echo either recovers the original order or reports one of the two unresolved endpoint occupancies. A reversed order has exactly zero probability. This adds a test on a coherence between the two contrast parities. The [two separate fast clock readings](../free_source_clock/MIXTURE.md) continue to admit positive elapsed-time mixtures when examined by themselves.

Let D project the excess onto the two written endpoints. On the even and odd branches the exact spatial generators are H + 32D and H + 33D respectively. H contains every other ordinary blank seam and the complete degree-reference correction. For a sharp launch at j write the two propagated waves as u and v. The local opposite-order probability is the sum of |u(x) − v(x)|²/4 over x outside the endpoints. At an endpoint x the unresolved load outcome instead has probability (|u(x)|² + |v(x)|²)/2. These formulas include the effect of the clock on the source. Removing sign exchange at magnitude two makes the two spatial generators identical, so the changed law predicts zero reversals at all times. [Reconstruction and checks](measure.py)

A path calculation gives a second check, including launches whose effect the counted detector cannot resolve. For a launch at distance d > 0, consider paths returning the source to its launch port with reversed order. In the difference of the two propagator series the first contributing walk uses d hops to reach an endpoint, one on-site insertion there and d hops back. If m(x) counts shortest paths from the launch to endpoint x of length d, every such product has the same sign and the difference of the two insertion strengths has magnitude one. The joint return-and-reversal probability starts at

$$
\left[\frac{\sum_{x}m(x)^2}{2\,4^d(2d+1)!}\right]^2 t^{4d+2}.
$$

At one link the coefficient is 1/576 at power six. For an endpoint launch, demanding a final excess at a particular adjacent non-endpoint gives coefficient 1/64 at power four. These witnesses lower-bound the inclusive opposite-order probability at sufficiently short positive times. They supply no minimum detectable rate at finite appointments. An independent integer recurrence checks the leading coefficient and all lower vanishing powers on the full board for every enrolled launch, including power 58 at fourteen links. [Exact coefficients](results.json) · [Independent checks](test_degenerate.py)

For [OW-27](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L131), the observation excludes pure retiming as an explanation of this apparatus's complete response on the stated source. A coarse-grained geometric delay could coexist with the interaction. Separating such a delay requires a calibrated causal signal or another transported clock on the same source. The [two-belt radar experiment](../EXTERIOR_PROTOCOL.md) uses a different preparation. No gravitational redshift, light bending or Shapiro delay is identified by this test.

Run from the repository root.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_27_metric.degenerate_clock.measure
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_27_metric.degenerate_clock.replicate
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_27_metric/degenerate_clock
```
