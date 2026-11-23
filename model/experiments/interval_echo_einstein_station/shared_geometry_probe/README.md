# Do two native probes see one source-dependent geometry?

A vacancy and an excess load have equal free propagation on an intercell dimer. An independent ordinary-source energy read on a second dimer predicts the terminal responses of **both** loads at nine bridge placements. The source can nevertheless polarize a probe launched as an equal mixture of its two ports. This observation distinguishes exchange of quantum information with the source from evolution in a source-dependent classical geometry on the complete two-port propagation space. A second experiment finds this distinction near a freely generated seam correlation when a fresh probe can traverse the **entire federation**. [Two-dimer protocol](PROTOCOL.md) · [Free-flight protocol](FLIGHT_PROTOCOL.md)

## What did the ordinary source predict?

Each port initially carries one unit. The source carries one additional or missing unit on its own blank dimer. Independent exchange-parity and inverse-echo seam measurements resolve two stationary source orbitals a unit apart in energy, with identical port populations and zero current. Calibrate the map from input orientation to the three terminal pointer probabilities using one ordinary source vacancy and one probe vacancy at the first bridge of the 20-cell board. Hold the calibration fixed for the opposite load signs, a quarter-mixture of source orbitals, further input orientations and the remaining bridges on 20, 80 and 320 cells. The experiment supplies coherent launches and calibrated coherent filters as instruments. Sharp endpoint counts give a separate check below. An independent source-free dimer clock reads each exposure with maximum phase error 0.000095 in native time units. [Recorded counts and addresses](results.json) · [Source meter and registration](measure.py)

| Source preparation | Held binary frequencies | Maximum absolute prediction error | Outside simultaneous envelope |
| --- | ---: | ---: | ---: |
| Lower ordinary orbital | 4,620 | 0.0001331 | 0 |
| Higher ordinary orbital | 4,620 | 0.0001548 | 0 |
| Quarter mixture | 4,752 | 0.0001418 | 0 |
| Symmetric orbital with odd-reversal marks | 4,752 | 0.9766295 | 4,680 |

The simultaneous 99% prediction radius is 0.0008374 for these simulated counts. On the calibration bridge a sharp probe vacancy reaches the receiving gateway after two exposures with probability 0.32374 (lower ordinary source), 0.34252 (higher ordinary source) and 0.32843 (quarter mixture). The collar-mark preparation raises the independently metered regional energy by one unit at fixed orbital. Its recorded arrival probability is 0.32438, and the largest error in the table occurs on a phase-sensitive pointer. Stress components or pressure and strain have yet to be measured on this family. The experiment compares sources of equal regional energy and equal mean load current.

## What distinguishes a quantum contact from propagation on a common geometry?

A completely mixed input on two ports has probability one half for each normalized binary analyzer. Every unitary on those ports preserves that state, including a classical mixture of source-controlled unitaries with an arbitrary wait or Hamiltonian. After four exposures, for *every* placement and sign, the largest measured component exceeds 0.42674 in absolute displacement from one half for either pure ordinary source. Its conservative simultaneous lower bound is 0.42656. Switching the sharp seam coins disconnects the contact on the registered law. Independently switching off the contrast term gives a changed-law control. Every enrolled mixed-probe component in both controls agrees with one half within 0.000177.

The equal responses of opposite loads follow directly from the native 99-state seam restriction. With probe endpoints a,b and source endpoints c,d its traceless operator is the traceless part of −(Xₚ + Xₛ)/2 − g |b,c⟩⟨b,c|, with g = 32 for both ordinary orbitals and 33 for the odd-reversal marks. The compensating sharp seams cancel the one-body position offset. Tracing out the dynamical source leaves an interacting probe with a nonunital channel, meaning a completely mixed launch acquires polarization. A geometry with explicit source recoil, additional traced modes, or a separate geometric contribution alongside contact remains compatible with that result.

## Does it survive on an unconstrained federation?

A freely moving load prepares two initially independent sharp seam heads at four exposures. Remove that load without postselection. On each new copy, inject a fresh vacancy or excess uniformly over all 240 or 960 ports. Compare its six endpoint and phase measurements on the armed bridge and on a bridge connecting owners two and three hops away. Any unitary on the whole position space, and every classical mixture thereof, predicts probability 1/N for each normalized detector. The product control obtains the two actual marginal head states from independent source copies. A fully mixed four-head preparation provides a second, exactly uniform control. [Terminal histograms and simultaneous intervals](flight.json) · [Full propagation and detector projections](flight.py)

| Cells | Generated collar, largest resolved near departure from 1/N | Correlation alone, largest resolved near departure from product |
| ---: | ---: | ---: |
| 20 | 0.0039119 | 0.0000013380 |
| 80 | 0.00015344 | 0.0000042674 |

The values are simultaneous 99% lower bounds over all 864 independent simulated cohorts. Every one of the 144 near readings on the generated collar resolves a displacement from uniform. All 48 sharp endpoint settings alone resolve the near departure from uniform, and 36 of those settings also distinguish the two actual head states from the independent product. Of the 144 comparisons to the product of the same head marginals, 115 resolve an additional correlation contribution. Neither comparison resolves a displacement at the distant station. All 288 settings with fully mixed heads include 1/N. The result measures a local action of the dynamically generated collar on matter. Determining proper curvature and universal gravitational coupling requires independently calibrated null stress and proper-distance observations at the same physical refinement.

## Reproduce the readouts

Run from the repository root with NumPy and SciPy. Both executable entry points require a fresh output filename. The forward-density regression checks all 960 input basis states of the 20-cell propagation, including complex phase analyzers. [Independent checks](test_probe.py)

```sh
python -m model.experiments.interval_echo_einstein_station.shared_geometry_probe.measure --output /tmp/ordinary-source.json
python -m model.experiments.interval_echo_einstein_station.shared_geometry_probe.flight --output /tmp/formed-collar.json
python -m pytest -q model/experiments/interval_echo_einstein_station/shared_geometry_probe
```
