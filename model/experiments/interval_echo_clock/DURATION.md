# Can a clock discrepancy disappear in a longer observation?

The answer depends on the motion, the measured quantity and the detector. The two trial hands use the same moving patch’s local history. One responds to every positive primary. The other responds to the elapsed event the overlap repair posts. Percentages below are raw-minus-accepted counts divided by expected home counts. The five-bit clock, reciprocal guide, primary opportunities and π/16 oscillator pulse are fixed throughout. [Frozen comparison](DURATION_PROTOCOL.md) · [All samples and controls](duration.json)

| Motion | Observation window | Expected extra raw ticks | Gap as a fraction of elapsed home time | Terminal oscillator at 8,192 preparations per analyzer |
| --- | ---: | ---: | ---: | --- |
| One 2-out/6-back journey | 8 opportunities | 0.17593 | 4.3983% | Unresolved |
| Repeat that journey 60 times | 480 | 10.55601 | 4.3983% | Resolved |
| One journey followed by 472 opportunities at rest | 480 | 0.17593 | 0.0733% | Unresolved |
| Repeat a slower 31-out/31-back journey seven times | 434 | 0.11293 | 0.0520% | Unresolved |
| Stationary throughout | 480 | 0 | 0 | Unresolved |

## Absolute lag and fractional lag

For a completed interval of q opportunities, let s(q) be the probability of suppressing a positive closing event. The mean number of omitted marks is s(q)/2. The two rules respond equally on ordinary rest events. In a window of T opportunities with n(q) completed intervals of each size, the mean absolute gap is Σ n(q)s(q)/2, while its fraction of the expected home count is Σ n(q)s(q)/T. Repeating one motion pattern therefore preserves its fractional gap. Following one trip with rest leaves the accumulated gap fixed, and the fraction falls as 1/T. Nothing in that dilution makes the traveller regain missed ticks. Expected quantities average over preparations. A single trip can accumulate an integer gap of zero or more ticks. A second material check follows the moving patch, hand, oscillator and all physical event bits for all 480 opportunities on the 20-cell screen. The repeated and parked branches share the exact first eight local records at each matched seed. Each of the 472 parked services posts its own ordinary primary, and both observers recover the same final elapsed word. Both fixed-seed first trips omitted zero ticks. The parked histories ended at zero and the 60-trip histories omitted seven and eight ticks respectively. The ensemble can have a positive mean even when both individual short trips read zero. The statistical panels sample the local conditional instrument. [Material histories](rest.json) · [Physical hold protocol](REST_PROTOCOL.md) No window reaches the nine-bit elapsed counter’s capacity. [Implementation and finite trajectory check](duration.py)

## What the detector can distinguish

A final rail measurement reads a periodic phase. Under an independent event with posting probability p, its mean complex fringe is multiplied by 1−p+p exp(iπ/16). Independent events multiply these factors. The angle wraps, and averaging variable tick counts reduces fringe visibility. In repeated 2/6 travel the raw visibility falls from 0.962 at eight opportunities to 0.099 at 480. Even a persistent event-rate difference can become arbitrarily hard to see with a fixed terminal-only detector and a fixed number of preparations. The positive accumulated record or a more closely spaced series of readings answers a different measurement question.

At 31/31, all three tested windows gave an unresolved terminal-only comparison at both 128 and 8,192 preparations per analyzer. At 434 opportunities the larger of the two predicted differences in terminal click probability was 0.135 percentage points. The observed four-standard-error resolution margins were 3.12 and 3.10 percentage points at the two analyzer settings. The paired repair-record measurement, which retains the individual omitted marks, estimated the fractional elapsed gap as 0.0507% ± 0.0017 percentage points (one standard error). The four-standard-error criterion in the frozen protocol applies separately at each exposure. It specifies this simulated apparatus’s sensitivity. A statement about physical clocks would require identifying which actual clock transitions couple to which repair records, and the resolution of an experiment that compares them. [Sample counts and visibility](duration.json) · [Clock instrument](README.md)

## Comparison with physical clocks

Chou and colleagues resolved ordinary kinematic time dilation in a comparison of two optical atomic clocks at relative speeds below 10 m/s. To test for an additional Interval Echo effect in such data, the model would have to identify the repair observable driving each atomic transition and predict the residual shift along the experimental trajectory after ordinary time dilation is accounted for. The present fractional gaps compare two deliberately chosen couplings inside the simulation. [Optical clock experiment](https://www.nist.gov/publications/relativity-and-optical-clocks)

```sh
python -m model.experiments.interval_echo_clock.duration
python -m model.experiments.interval_echo_clock.rest
python -m model run --law interval_echo --clock-trip --rest-after 472 2 6
```
