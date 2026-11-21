# What a mode-selective preparation measures on the free screen

**Six native pulse programs isolate six distinct standing energy modes. A coherent pair transfers occupation across a cell boundary after release.** The independently prepared zero-wave-number state splits across several energies, and its traveling pulse gives different peak speeds at successive receivers. These results resolve the temporal spectrum of this preparation. They leave the massive Klein–Gordon dispersion in OpenWave row 30 open. [Protocol](PROTOCOL.md) · [Raw primary counts](results.json)

## Free spectroscopy and transfer

The 20-cell sphere carries 630 written collars with sharp ordered arrows in even sign superpositions. One unit starts at port 0. All physical hops allowed from this start span 32 joint configurations on 12 ports in three cells. Each move reverses a physical collar. For preparation we alternate 96 short free intervals with eight independently addressed native degree-reference phases. Every seam acts freely from release to reception. Terminal observers recover the load from incident arrowheads and the locally known initial records. Phase control, the ordered collars and the native exposure intervals are supplied apparatus. This experiment times its readings by native exposure. [Generator and source](source.py) · [Frozen preparation settings](controls.json) · [Detector identity](../row_30_arrows/spectrum.py)

The six selected energies had the largest launch overlap in the pilot. With all phases disabled the largest energy weight is 33.56%. With the compiled controls the six single-mode states and their five coherent pairs each have preparation infidelity below 1.7 × 10⁻¹¹ by numerical replay. No trials are discarded. The isolated states' root-occupation probabilities change by less than 7.4 × 10⁻⁷ over 16 native exposure units. The oscillations below occur when a second prepared level provides a phase reference.

| Level paired with level 1 | Frequency from first 257 exposures | Held R² on next 256 | Observed range of occupation outside the home cell |
| ---: | ---: | ---: | ---: |
| 4 | 5.69179 | 0.999979 | 27.31 percentage points |
| 12 | 14.30217 | 0.999976 | 4.29 percentage points |
| 0 | 10.09269 | 0.999946 | 2.45 percentage points |
| 31 | 25.75996 | 0.999920 | 0.23 percentage points |
| 22 | 23.44706 | 0.999898 | 2.00 percentage points |

The decoder scanned angular frequencies from 0.5 to 50 using home clicks alone. All five coherent pairs passed. All five corresponding incoherent mixtures missed the held criterion. Each exposure has 262,144 simulated terminal trials, including every quiet outcome. For the strongest pair, 14,148 loads were already outside the birth cell at release. At exposure 0.5625, the sampled appointment nearest the first half-cycle, that count was **84,762**. This is coherent intercell redistribution in a finite standing pair. [Readout and controls](experiment.py) · [Time series](results.json)

## Why there is no rest-frequency or speed calibration here

We gave one distance unit to each allowed path step across a physical seam. The resulting graph Laplacian defines an independent squared wave number. For each of its 32 spatial basis vectors we compiled a native phase sequence. Reversing the order of its elementary pulses realizes the transposed sequence. Applying it to a fresh input and counting returns to the unique root measures the squared projection on the desired real spatial vector. This readout uses positive service durations. Each of the 32 analyzers was tested on a separate ensemble, including its quiet outcomes. [Spatial analyzer](analyzer.py) · [Frozen settings](analyzers.json)

- In the six isolated energy modes the upper bounds on the largest single wavelength weight run from 16.48% to 32.87%. These simultaneous bounds include count uncertainty and analyzer preparation error. The registered gate was 90%.
- Preparing spatial zero independently gives **29,978 of 262,144** clicks in level 1 and **43,453 of 262,144** in level 0 on separate ensembles. Their weights each exceed 10% even after the simultaneous error allowance. The six chosen level analyzers together account for about one third of this source. Its calculated free return probability falls to 0.000148.
- Across every energy eigenspace, including optimal combinations in degeneracies, the sole 90% match has energy 4 and squared wave number 1. Its support is entirely on ports 9 and 11 of one cell. Projecting instead onto harmonics of the full 240-port sphere with an equal coherent weight across histories at each port gives at most 27.16% weight in any one harmonic eigenspace. This full-screen alternative is a numerical projection.
- The independently released sharp load gives count-peak times 0.14453, 0.24219 and 0.25391 on consecutive physical ports 1, 82 and 83. The successive speeds are **10.24** and **85.33** native seams per exposure. Exact-probability peaks also miss the 10% agreement gate. These short-range maxima cannot set a common speed for this packet.

The uniform rest has a simple obstruction. At a joint path configuration, the generator acting on a constant amplitude multiplies it by κ times the difference between the blank-reference count and the allowed path degree. That coefficient ranges from −48 to 0. The energy variance is 137.75 in squared native units. After the primary experiment we tested a staggered convention, reversing the spatial phase at every step from the root. This transforms the generator into a multiple of the path Laplacian plus the native diagonal degree-reference potential. The staggered zero mode has smaller variance, 12.0. Five energy eigenspaces now meet the 90% gate at nonzero wave numbers, including one two-dimensional eigenspace. Staggered rest still gives **148,892** and **48,390** clicks out of 262,144 in two separate energy analyzers. Neither convention supplies an isolated rest mode. [Follow-up implementation](staggered.py) · [Frozen staggered controls](staggered_controls.json) · [Counts and all 31 energy eigenspaces](staggered.json)

## What this supports for OpenWave

The public row asks for emergent massive dispersion of a matter wave. The platform's M10 review distinguishes a kinetic source with its dispersion supplied from evolution that produces the measured curve. The M5 review identifies the frequencies of perturbations measured under time evolution as the relevant evidence there. Our five held predictions are transition frequencies for five prepared pairs at later times. A direct test of the row would instead use a transported free branch's measured rest frequency and speed to predict its angular frequency at **other wavelengths**. The rest and velocity calibration gates both failed on this background. Estimating a mass or comparing held wavelength fits against a parabola would exceed these records. [Public OpenWave matrix](https://github.com/openwave-labs/openwave/blob/main/MODELS.md) · [Public row review](https://github.com/openwave-labs/openwave/pull/741)

The separate driven-guide experiment provides a positive control for resolving wavelength dependence. Its detector counts predicted 48 held signed momenta within 1.05% across three screens. The guide's alternating contact durations supply its gap. The present free experiment establishes what the same candidate's ordered written-arrow preparation retains when every native seam runs after release. [Driven guide and independent calibrations](../row_30/README.md)

A separately prepared 194-cell static channel now gives the affirmative free comparison. Its frequency error is at most 1.19% across four held wavelengths. An independently calibrated parabola misses by 6.68% on the fourth. The independent review reanalyzes its raw terminal counts, the physical energy reference, the independent speed calibration and nearby competing dispersions. [Static-channel review](CHANNEL_REVIEW.md) · [Count uncertainty and alternatives](channel_review.json)

## Replay

From the RER repository root with NumPy and SciPy installed, use the frozen pulse files for the complete replay. The first two commands compare the new counts and spectral audits with their archived JSON receipts. The tests independently check the transposed analyzer on a complex input, exact free propagation, disabled phase controls and the two rest variances.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_selective.experiment
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_selective.staggered
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_selective
```

To compile fresh settings, run the `source` and `analyzer` modules with `--output NEW.json` separately or the `staggered` module with `--compile NEW.json`. Numerical optimizers can generate different pulse sequences with comparable preparation fidelity. The two replay modules accept `--output NEW.json` to preserve an independent receipt.
