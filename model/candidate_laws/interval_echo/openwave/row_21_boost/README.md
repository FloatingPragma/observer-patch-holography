# A material pair shortens as it moves on the scheduled guide

Two equal loads on a native route acquire a momentum-dependent separation in a spectrally selected state. On the 80-cell board the pair spans **9.81 port spacings at rest and 4.38 at 0.879 of the measured signal-front speed**. A separate center-fringe measurement tracks its motion. The observed contraction differs measurably from the Lorentz length formula. Even fitting the limiting speed cannot account for the three moving profiles on a given board. The local finding is a contracting material configuration under a supplied laboratory schedule. OpenWave lists the M10 Lorentz criterion as untested. [Criterion](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L125) · [Protocol and prior inspection](PROTOCOL.md) · [Raw outcomes](results.json)

## How the pair gets its shape

Use the two staggered matchings of the existing one-load Klein–Gordon guide. The three closed routes visit 110, 214 and 246 actual ports and cross 15, 36 and 41 cell boundaries. Two equal doubled loads start on adjacent ports. A translation-phase source sets their total winding. On a matching containing only one double the exact 99-state seam gate acts on the two possible endpoint positions. If both doubles occupy the same seam it acts on the distinct native input (2,2). Together these rules determine the interacting pair. The unchanged sharp source offers only unit separation at every winding. The terminal distribution of internal distances comes from the two-load evolution under a 16,384-period coherent frequency filter.

In each momentum sector choose the isolated line with the largest summed overlap on the two adjacent input orientations. This selects a branch of the scheduled material using model-informed spectroscopy. The nearest gap at each setting, admission, fidelity and all rejects are recorded. An independent sample counts both occupied ports at the end of a period. The observable is their shorter arc separation in port spacings. The rest and fastest profiles are also measured 16 periods later. Translation of the center or spreading of its wave packet leaves this relative coordinate unchanged. Disagreement with the contraction forecast can therefore be detected even when the center's location is completely delocalized. The laboratory supplies both matching durations, the spatial phase source, frequency filter, shutters and simultaneous reads. The uniform timing in both directions is part of this apparatus. [Original guide and one-load dispersion](../row_30/README.md)

## Speed, size and a held common-speed test

At the beginning, middle and end of one period, independent copies report the number of doubles on the odd ports, between zero and two. The right endpoint of every active seam has odd parity in the first layer and even parity in the second. Thus the center velocity in two-port cells per period is ((2\langle n_o\rangle_{1/2}-\langle n_o\rangle_0-\langle n_o\rangle_1)/4). The count averages include blocked contacts occupied by both endpoints. Sending one sharp double in either direction reaches exactly six then twelve cells at those numbers of periods on the *same* detuned guide. The observed front probabilities agree with (cos^{2t}(0.2)). A signal moving beyond one cell per period is excluded by the two local layers. This is a causal front speed along the driven guide.

| Board | Moving winding magnitude | Measured speed / front speed | Measured mean separation / rest | Lorentz prediction using front speed |
| --- | ---: | ---: | ---: | ---: |
| 20 cells | 2, 3, 7 | 0.517, 0.667, 0.886 | 0.840, 0.717, 0.437 | 0.856, 0.745, 0.464 |
| 80 cells | 3, 7, 13 | 0.423, 0.729, 0.879 | 0.893, 0.650, 0.446 | 0.906, 0.685, 0.477 |
| 320 cells | 4, 8, 15 | 0.475, 0.727, 0.879 | 0.863, 0.652, 0.445 | 0.880, 0.687, 0.476 |

All eighteen signed comparisons separate the observed ratio from the front-speed forecast under the simultaneous 99% counting bounds. The intervals condition on the actual number of accepted trials. Each setting starts (2^{32}) simulated copies; 2.5–5.1% pass the mixed-source tag on the registered material. The onsets, waveform selection and three speed ranges were explored before the count checks. The probabilities predicted by the chosen microscopic law are the sampling model. A changed law with the contrast coupling removed leaves the same speed front and contracts by a different amount. On 80 cells its fastest measured ratio is 0.375.

At each speed independently invert the same formula for a possible limiting speed. On 80 cells the lowest moving winding requires **0.9253–0.9568**, the middle **0.9547–0.9631**, and the highest **0.9797–0.9845** guide cells per period. The low and high intervals are also disjoint on the other two boards. This excludes an exact common-speed Lorentz-contracted family in the enrolled finite regime. A long-wavelength limit or a different material pair requires its own preparation and frame test.

## Does the center actually travel?

Coherently admit adjacent windings on either side of a central one and measure the midpoint of the shorter occupied arc. Its spatial fringe shifts under the same scheduled evolution. At the three 80-cell settings the second source orientation gives speed magnitudes **0.4196, 0.7271 and 0.8784**. Each agrees with its independently computed two-line phase difference divided by wave-number difference. Across the three boards that source has 4–10% admission, visibility 0.46–0.50 on registered material, and conditional selected-line fidelity above 0.999999. Finite sidelobes and sampling of both quadratures enter the secant speed bounds. The single-winding current measures the central velocity; the fringe uses two distinct windings and agrees more closely with that current as their spacing decreases.

The first orientation has very little overlap with the selected branch. On the 80- and 320-cell registered guides, its accepted fraction is only about (5\times10^{-8}) to (8\times10^{-8}). Sidelobes and small counts defeat its phase-speed read on all three boards under the stated bounds. The receipt retains those attempts and the targeted second-orientation measurements. The finite clocked guide exposes an underlying preferred period. Equal-time endpoint reads occur in one frame. Observer-independent spacetime covariance would need transformations between physical frames and a common invariant dynamics over the domain being claimed.

## Replay

From the repository root using the model's scientific Python environment:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_21_boost.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_21_boost
```

The [program](measure.py) recreates all six material and control runs. The [checks](test_boost.py) multiply the full two-body native gates on a small ring, compare all translation blocks, and independently reconstruct velocities and separations from the sealed counts.
