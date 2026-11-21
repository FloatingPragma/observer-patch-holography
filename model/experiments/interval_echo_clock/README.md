# A second clock at Interval Echo reunion

Two finite patches can bring a coherent seam oscillator and the registered Interval Echo hand to the same return. When the oscillator advanced on every positive local repair, its moving/home elapsed ratio was consistent with one on all three journeys. The registered hand read about 0.956 or 0.968. Which local repair events advance the clock makes a measurable difference. [Frozen protocol](PROTOCOL.md) · [Retained measurements](results.json)

| Moving oscillator advances on | Out/back service opportunities | Registered traveller/home | Oscillator traveller/home | Oscillator minus registered, ±1 SE | Frozen equivalence decision |
| --- | --- | ---: | ---: | ---: | --- |
| Positive primary repair | 2/6 | 0.95620 | 0.99937 | +0.04317 ± 0.00159 | Disagrees |
| Positive primary repair | 6/2 | 0.95630 | 1.00047 | +0.04417 ± 0.00156 | Disagrees |
| Positive primary repair | 4/4 | 0.96826 | 1.00120 | +0.03294 ± 0.00164 | Disagrees |
| Accepted Interval Echo tick | 2/6 | 0.95596 | 0.95580 | −0.00016 ± 0.00137 | Agrees |
| Accepted Interval Echo tick | 6/2 | 0.95585 | 0.95482 | −0.00104 ± 0.00123 | Agrees |
| Accepted Interval Echo tick | 4/4 | 0.96805 | 0.96589 | −0.00216 ± 0.00126 | Unresolved |

The registered word and the oscillator readout share a repair trajectory in each preparation. Ratios in different rows come from separate ensembles. Agreement requires the entire observed difference ±4 block standard errors to fit in ±0.006. The 4/4 accepted-tick interval extends to −0.00718. The exact finite-branch expectation for that difference is approximately −0.00010. The finite sample therefore leaves its frozen equivalence question unresolved. Across all rows the measured visibility of both hands was above 0.960. The separate reciprocal-optics study gave elapsed forecasts of 0.95595, 0.95602 and 0.96817 for these schedules. Its supplied optics and embeddings have their own scope. [One-link optical observations](../../evaluation/ONE_LINK_CLOCK.md)

## Independent replication of the uniform journey

A separately frozen experiment repeated the 4/4 journey with 8,388,608 fresh preparations per analyzer and clock rule. It reused the same gates and tolerances, enrolled fifteen fresh material witnesses and had one fixed statistical look. The accepted-tick dimer differed from the registered hand by **+0.000046 ± 0.000497** (one block standard error). Its four-error interval [−0.001944, +0.002035] fits within the ±0.006 equivalence band. On the matched positive-primary control, the difference was **+0.032042 ± 0.000444**. Its interval [0.030267, 0.033817] clears the upper limit. The larger sample resolves the uniform accepted-event comparison while independently reproducing the raw-primary discrepancy. [Replication protocol](REPLICATION.md) · [Independent receipt](replication.json)

The ordinary material clock now places the accepted event on paired load rails. Both enrolled hands read those rails through native local contacts before advancing. The separate publication study tested every overwrite and both terminal readouts on the three spherical boards. [Physical event publication and controls](PUBLICATION.md)

## What the patches do

An oscillator begins with a double on one end of an isolated native seam and an empty opposite rail. Its next local positive primary is produced by ordinary signed repair. In the primary instrument each such record opens a π/16 pulse of the same native Hamiltonian. In the accepted-tick instrument a local closing repair may suppress this pulse under Interval Echo’s existing rule. The home oscillator advances at each positive home primary in either case. The oscillator’s state stores its accumulated phase physically. Initial alignment, the pulse duration, the protected seam, transport windows, the guide and cadence are supplied. This probes two possible local event interfaces for a clock. The accepted-tick oscillator and the registered hand share the same elapsed-event filter. [Native seam law](../../candidate_laws/contrast_memory_spherical_v1/law.py) · [Material clock](../../candidate_laws/interval_echo/clock.py)

The moving oscillator occupies spare rails 6 and 7 of a bank in the traveller’s finite patch. The home oscillator occupies a private bank. On every addressed transfer, the ordinary π contact swaps all four possible sharp binary inputs with the same phase. Linearity propagates both coherent rail amplitudes through the permutation even when its path crosses another occupied rail. During stationary services the probe remains isolated except for its local exposure. Forty-five material records exercised five seeds, three schedules and the three spherical boards of 20, 80 and 320 cells. The experiment tracked the two coherent rail labels and their read vacancy through every elementary swap and verified their return. Each conditional detector branch can be published as a signed bit at the reunion seam and read by both owners. [Route, pointer and publication checks](experiment.py)

At reunion a local read of rail 6 against vacant port 10 measures a cosine. Repeating the preparation with an additional π/2 analysis pulse gives the other quadrature. The argument of the resulting mean complex signal gives the ensemble clock reading. All eight opportunities fit within an unwrapped phase interval for these windows. Under raw-primary gating the home and traveller undergo the same positive-primary probability at every opportunity, so their ensemble phase ratio is exactly one. Under accepted-tick gating the characteristic function of the local accepted count predicts oscillator ratios 0.955883, 0.955883 and 0.968145. The small difference between the mean phase and the mean elapsed count follows from averaging a periodic pointer over random marks.

For each of the six rows and two analyzer choices the retained data contain 1,048,576 simulated preparations in 32 blocks. The cohort uses the exact joint law of local primary and accepted repair words conditional on the material handoff records and samples both terminal probes and both elapsed words. Identical handoff sequences on different boards share a cohort. The raw receipt also records the local words along every material witness, the conditional endpoint probabilities and the signed readback for both possible bits. It includes all quiet outcomes. A dimer with its pulse switched off has zero phase. [Block sums and native witness records](results.json)

## Reproduce

Run from the RER repository root.

```sh
python -m model.experiments.interval_echo_clock
python -m model.experiments.interval_echo_clock.replication
python -m pytest -q model/experiments/interval_echo_clock
```

The first two commands recompute the original experiment and the independent replication and compare their complete results with the stored receipts. The independent tests check the signed pointer directly on the complete native seam, compare the retained block sums with the measured ratios, and check the zero-response and stationary limits.


## Duration, sensitivity and interval capacity

Repeated 2/6 travel preserves a 4.3983% expected difference between raw and accepted timing. After a single trip followed by rest the absolute tick gap persists while its fraction of total elapsed time tends to zero. At slower sustained 31/31 motion the expected difference is 0.0520%. The terminal oscillator did not resolve it at 8,192 preparations per analyzer. The paired repair histories did. [Long-window controls](DURATION.md)

The precision trip uses five-bit interval fields. Exact finite-walk enumeration bounds the mean-rate discrepancy from the reciprocal unit-link target by 0.000488401 over all positive integer intervals. At the three-bit capacity boundary, enlarging the bank also repairs the physical 8/8 optical comparison on the same guide. [Finite capacity certificate](CAPACITY.md) · [Held optical comparisons](../interval_echo_optics/README.md)
