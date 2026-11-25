# Can observers identify a gravitational source in a collar correlation?

The source observer records two head readings and their joint parity. The second observer records local clock ticks and where a pulse is received, including missed receptions. A third reads energy transfer using a phase kick and a binary energy measurement. These are three independent instruments on copies of the same preparation. The hardware addresses route the pulses. Directions, image areas and focusing below refer to quantities inferred from the observers’ instruments. The flagship shows how the predictions of a finite admissible observation menu on bounded wave times can have a smooth relativistic 3+1 description to specified accuracy. The gravitational question is whether one description predicts how the observed focusing changes when only the joint correlation changes. [Executable model’s measurement principle](../../../../docs/OPH_EXECUTABLE_HYPOTHESIS.md#how-the-model-is-tested) · [Flagship proposition on smooth effective field descriptions](../../../../flagship/from_observer_consensus_to_standard_physics.tex)

## A source selected through joint measurements

Hold the two individual head distributions and the entire parked-energy distribution fixed. Measure the nine pairs of Pauli components separately in each parity sector. The eighteen signed means give an orthonormal coordinate system for the allowed joint variations. Hold the four-station survey volume fixed to first order. Within these settings, the component of the modular readout along the allowed directions selects the largest positive modular change per unit Euclidean change in the joint readings. That normalization is a search choice. No camera response enters it. Both endpoint preparations have positive weights on the executable preparation menu. [Selection protocol](PROTOCOL.md) · [Explicit recipes, source counts and directional predictions](results.json)

The measured slopes refer to the calibrated modular reader in nats. Its weights and the public preparation settings are inputs to the counting intervals. The expected difference between the two source readings is +0.00022946 nats. The first and precision count banks gave the following results. The intervals have simultaneous 99% coverage over the Pauli settings in each separate bank, conditional on the simulator’s probability law.

| Simulated trials per setting | Counted difference, nats | Simultaneous half-width, nats |
| ---: | ---: | ---: |
| 67,108,864 | +0.00008805 | 0.00106743 |
| 68,719,476,736 | +0.00022940 | 0.00003336 |

The two endpoints give the same conditional one-head readings. The odd-energy outcome probability is 0.20243684 in both. At a decorrelated reference reconstructed from the two marginals the modular derivative vanishes. This distinguishes a joint source record from changes in the individually visible material or its parked energy.

## What the other observers record

An ordinary intervention calibrates the relation between energy current and the difference of forward and reverse optical focusing on the vacancy probe. At the four-tick appointment this calibration predicts an interval from −442.40 to −292.88 for the collar. The camera observer obtains +7.08 to +16.10. These are focusing changes per unit preparation change, in the same calibrated clock and detector units. A decoder using only the reception, clock and phase-kick counts and instrument settings gives both intervals. The total allocated counting risk is 0.009. The finite appointment spacing is 1/64. Using the analytic zero-spacing and infinitesimal-source limits gives −350.80 predicted and +11.06 observed. The sign discrepancy holds in the finite records themselves. [Raw clock, meter and camera counts](observer_records.json) · [Count decoder](observer_records.py)

The comparison that uses the sum of forward and reverse focusing also has the opposite sign to the calibrated work prediction on both probes. Those are calculated response derivatives. The claim about the vacancy difference above has an independent, counted acquisition. The choice of exchange clock, material pulse and identification of measured bridge flux with directional stress are part of the hypothesis under test.

## A seam that carries actual interaction energy

A distinct four-word seam provides a concrete nonparticle energy candidate. Two preparations with identical separately measured marginals have energies 12 and 4 on the isolated native contact. A thermometer and joint counting meter measure their difference as 7.95234 ± 0.02032 energy units (one counting standard error). With the coin phase reversed, the energy difference vanishes. With uniform contact energy coupled to a supplied homogeneous gravitational dynamics, observers’ radar records give a volume-weighted curvature change of 12.79997–12.80000 at three held times. A comparison gravitational law reproduces the ordinary calibration and gives zero collar curvature on the same matter records. These calculations identify exactly which inference depends on assigning the seam energy to gravity. The two-head experiment above uses a different preparation and tests its own connection to clock and reception records. [Native energy, measured curvature and alternative coupling](../../../../code/matched_einstein_limit/seam_overlap/README.md) · [Equal-volume optical rims](../../../../code/matched_einstein_limit/seam_overlap/rim.json)

## Comparison with astronomical records

An astronomical identification would require using the same collar prescription, its distribution and the ordinary coupling to predict independently measured orbital velocities or lensing. The measured pair statistic in nats has no established conversion to a galactic stress profile. One existing extension of the compact-source law to a gas disk predicts at least 62.97 km/s at 9.95 kpc in AGC 114905. The measured, pressure-corrected speed is 26.94 ± 1.83 km/s. This is a failed prediction of that disk extension. The observer-record contrast above tests the proposed finite source identification at a different scale. [Disk calculation and data provenance](../../../../code/cosmology/rar_deep_regime/cut_spacing/README.md) · [Mancera Piña et al., observations of AGC 114905](https://arxiv.org/abs/2404.06537)

Run from the RER root with the project’s Python environment.

```sh
PYTHONPATH=code:. python -m openwave.interval_echo.row_19.selected_collar.measure
PYTHONPATH=code:. python -m openwave.interval_echo.row_19.selected_collar.observer_records --counts-only
PYTHONPATH=code:. python -m openwave.interval_echo.row_19.selected_collar.search --counts-only
```

The optical count decoders use published instrument settings and complete count banks. Omitting `--counts-only` regenerates the simulated events. Each high-count experiment samples batches from the computed probability laws. The uncertainty intervals describe simulated counting uncertainty.
