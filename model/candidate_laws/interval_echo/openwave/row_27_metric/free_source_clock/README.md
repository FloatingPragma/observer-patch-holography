# Can one elapsed time explain two clocks beside a freely moving source?

**One mobile excess changes both written clocks, and their records exclude a single deterministic change of elapsed time as the sole explanation at three launch positions.** The two clocks have independently calibrated blank frequencies 32 and 33. Their faces are ordinary binary seam coins. Each observer records whether the final coin differs from the starting coin. The same observer can also detect whether the excess reached either clock endpoint. The complete 240-port host evolves freely with every native interaction open. [Preparations and readouts](PROTOCOL.md) · [Count banks](counts.npz) · [Native experiment](measure.py)

The clock sits on a filled seam at ports 0 and 5. Its prepared contrast parity selects the clock standard. A sharp excess, supplied independently after preparation, may move past, through or away from the clock's ports. Both starts of the clock coin count equally. At each of 193 matched appointments each clock supplies 1,048,576 terminal records. A separate material dimer times the exposure at one eighth duration. One port per initial native separation is enrolled. The count of links labels a preparation and supplies no proper-distance scale. The source's separately measured formation energy on the filled host is 51.00 exchange units at the degree-six launch and 66.50 at each degree-five launch. These energies and launch geometry vary together in the one-link comparison. [Source energy and clock fits](results.json) · [Independent calibration of the two lines](../../row_03_clock/written_coin/README.md)

| Initial separation (port) | Largest change in even flip probability from blank | Largest change in odd flip probability from blank | Matched appointments excluding a common elapsed time |
| --- | ---: | ---: | ---: |
| At the clock (0) | 98.79 percentage points | 99.09 percentage points | 98 / 193 |
| 1 link (1) | 5.59 | 5.44 | 125 / 193 |
| 2 links (2) | 17.05 | 16.30 | 132 / 193 |
| 3 links (13) | 0.16 | 0.15 | 0 / 193 |
| 14 links (156) | 0.17 | 0.18 | 0 / 193 |
| Source absent | Reference | Reference | 0 / 193 |

The conservative simultaneous 99% counting radius for each source-minus-blank flip difference is 0.574 percentage points. Overlap with the far control gives a limit at the measured exposure and precision. The study makes no radial-field fit. At zero contrast the even clock gives zero flips for every launch and wait. The odd clock has blank frequency 0.99989, and an endpoint launch changes its flip probability by up to 56.31 percentage points. This control preserves the alphabet, geometry, preparation and exchanges.

## A time comparison that eliminates the clock coordinate

A single deterministic proper time τ accumulated by the two calibrated clocks would give even and odd flip probabilities sin²(32τ/2) and sin²(33τ/2). Hence the even records predict the same value of cos(1056τ) as the odd records. In terms of the Chebyshev polynomials this requires

$$T_{33}(1-2p_{\rm even}) = T_{32}(1-2p_{\rm odd}).$$

The condition applies to each individual exposure, including a time-dependent lapse between the start and that exposure. There is no fitted frequency, trajectory or winding count. At the one-link launch, appointment 108 gives an even-clock interval **[0.9508, 0.99995]** and an odd-clock interval **[−1, −0.9671]** for the common quantity. These disjoint intervals account for sampling variation in all 4,246 registered probability estimates together. The rule was designed after the first frequency fits. It uses inversion of the exact binomial distribution and includes all interior extrema of the two polynomials. Removing the source, or holding it at any of the five enrolled positions, gives overlapping intervals at every appointment. A source parked at an endpoint stops both clocks, satisfying the equality at zero elapsed phase. The freely moving source produces a different result. [Derivation and scope](COMMON_TIME.md) · [Per-setting intervals](common_time.json) · [Calculation](common_time.py)

There is also an exact result at the two launches where the count intervals overlap. Expand the difference of the two Chebyshev expressions around zero exposure in exchange units. Its first nonzero coefficient occurs at orders 6, 6, 8, 10 and 32 for the five launches in the table. For one link it is −513040 times t⁶. For 14 links it is −2058299 times t³² divided by 185940969061374140350464000000. The coefficients follow by multiplying integer powers of the two native position generators and composing the resulting even polynomials using rational arithmetic. Thus the ideal records depart from deterministic common retiming for each enrolled launch. The counting experiment resolves that departure at the three closest launches. [Exact expansion](series.py) · [How distributions of common times compare with these counts](MIXTURE.md)

A fit of a single sinusoid to the first half of the even-coin records misses the second half by RMS 0.0343 at one link and 0.0363 at two links (design tolerance 0.025). Launching on the clock itself suppresses the fitted excursion below 0.006 for both standards. The blank frequencies are 32.0041 and 32.9976. At three links, one cross-standard frequency forecast gives RMS 0.02525, narrowly above the tolerance, while both single-standard fits pass and the source-minus-blank differences lie inside the counting envelope. The paired identity supplies the separate coordinate-free test in the table.

## What the excess changes

In fixed contrast parity let E be the blank gap, D the projector onto source positions 0 and 5, and H the excess generator from the ordinary blank contacts and the host degree references. The exact native columns on the symmetric-coin and antisymmetric-coin branches give H + E D and H + E I. With both endpoints filled, only the antisymmetric coin costs E. With a double load at either endpoint, both coins cost E. The source can therefore distinguish the two clock branches simply by scattering. The sharp coin flip probability is one quarter the squared norm of the difference of these two freely propagated source waves. Coarse reading of the final coin sums over the complete source position record.

Separate cohorts select one coin-energy branch on a fresh written seam before installing the source. At one link, appointment 7, the excess arrives at an endpoint with conditional probabilities **0.0219% and 4.1313%** for the two arms of the even clock. Every rejected preparation and selected miss is retained. Changing how the observer starts the clock can therefore change where the source arrives. An independent check evolves all four sharp coin/contrast labels with their original native contact matrices. [Operator and count checks](test_measure.py)

This experiment excludes one deterministic common local time acting alone on these two *ideal* clock faces. A shared gravitational contribution on top of the observed material scattering requires an independent prediction that separates the two. The result assumes matched exposure shutters and supplies neither an autonomous astronomical source nor an exchanged light signal.

Run from the repository root:

```sh
python -m model.candidate_laws.interval_echo.openwave.row_27_metric.free_source_clock.measure
python -m model.candidate_laws.interval_echo.openwave.row_27_metric.free_source_clock.common_time
python -m model.candidate_laws.interval_echo.openwave.row_27_metric.free_source_clock.series
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_27_metric/free_source_clock
```
