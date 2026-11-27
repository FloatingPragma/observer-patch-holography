# A freely arriving wall and its rest-frequency test

The [OpenWave de Broglie clock criterion](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#simplest-test-per-criterion) asks for a self-starting, frequency-rigid internal oscillation at $\omega=mc^2/\hbar$. M10 carries the upstream development mark 🚧. A free endpoint read on two native massive strands detects arrival of a wall. Neither its first nor its seventh physical dimer subsequently oscillates at the independently reconstructed bulk rest rate of 15. Exact algebra excludes that frequency for every initial state and every fixed terminal measurement on the 13-, 29- and 49-dimer open strands. The full local criterion is Partial. [Prospective protocol](PROTOCOL.md) · [Proof, including the eight-bond counterexample](PROOF.md)

## Arrivals and terminal records

Prepare one wall in a sine envelope on bonds 1 through 15, with zero amplitude at the endpoint. After two free exchange-clock radians read the first dimer. Its single orientation announces whether the endpoint contains an unsatisfied bond. Preparing the reflected envelope controls for walls created near the detector by free pair production. The near input increases the arrival probability by 0.002661 and 0.002676 on the two lengths, with a sampler standard error of 0.000016 on each difference. A particular click need not come from the prepared wall.

After either answer, allow further free flight and read the first and seventh physical dimer orientations together. At each of 2,048 appointments the eight joint words are seeded multinomial draws from their calculated finite-law Born probabilities. Each appointment enrolls $2^{24}$ copies. The hand, phase and density totals are independent binomial or multinomial draws from their respective probabilities. A two-quadrature dimer hand measures the mean time step as 0.009999992 exchange radians. [Physical probabilities](apparatus.py) · [Counts](counts.npz) · [Counts-only summary](results.json)

| Source | Arrival near | Arrival mirrored | Amplitude at 15, first port, miss / hit | Amplitude at 15, seventh port, miss / hit |
| --- | ---: | ---: | ---: | ---: |
| 29 dimers, native | 0.003412 | 0.000751 | 0.000156 / 0.000187 | 0.000056 / 0.000103 |
| 49 dimers, native | 0.003416 | 0.000740 | 0.000155 / 0.000157 | 0.000048 / 0.000122 |
| 49 dimers, half contrast | 0.005528 | 0.002684 | 0.000134 / 0.000111 | 0.000127 / 0.002569 |
| 49 dimers, planted oscillator after trigger | 0.003420 | 0.000743 | 1.000017 / 1.000081 | 0.368343 / 0.702624 |

The native ceiling was 0.01 and the planted first-port minimum was 0.97. The planted arm sets the post-trigger contrast to zero and the exchange strength to 7.5. Finite-window leakage, hand noise and shot noise enter the tabulated Fourier amplitudes. The width of the interval excluded in the exact native word laws follows from the generator. Any fixed function of at most seven simultaneous bond-wall bits on any finite open matched strand has no line in $(14,18)$, including after either endpoint outcome. On the 49-dimer strand an eighth-port many-body transition falls between 15.23386 and 15.23402. This exact counterexample concerns a different, coherent many-body preparation. An irreducible reflection polynomial on each of the three named open strands also bounds every *rational* frequency of any fixed reader by one in magnitude. [Algebra and scope](PROOF.md)

## Independent weak-kick response of the same source

Fresh copies of the near envelope receive opposite spatial phase gradients. Adjacent-wall phase quadratures measure the kicks. Separate bond clicks at four hand times determine the difference in the first moments of wall occupation. A single copy contributes its measured bond outcome at the assigned address. We retain every outcome, including additional walls during flight. For a measured half-gradient $k$ and a half-difference of velocities $v$, $k/v$ estimates a finite kinetic mass in dimer and clock units with $\hbar=1$. A separately calibrated material coefficient is $c^2=16.00156$ (95% sampler bootstrap interval 15.99361–16.00977). [Coupling meter](../../row_30_collective/coupling_meter.json) · [Independent ring rest spectroscopy](../../submission/row_30_cautions/README.md#the-ring-source)

| Strand and kick | Measured $k$ | Measured $v$ | $k/v$ | $(k/v)c^2$ ± sampler s.e. |
| ---: | ---: | ---: | ---: | ---: |
| 29 dimers, smaller kick | 0.09976 | 0.10374 | 0.96166 | 15.388 ± 0.084 |
| 29 dimers, larger kick | 0.20014 | 0.20648 | 0.96927 | 15.510 ± 0.042 |
| 49 dimers, smaller kick | 0.10009 | 0.10415 | 0.96101 | 15.378 ± 0.085 |
| 49 dimers, larger kick | 0.19919 | 0.20645 | 0.96483 | 15.439 ± 0.043 |

These four recoil products are inside the exact silent interval. Each exceeds the independently reconstructed ring rest rates, 15.000082 and 14.999985. The tabulated standard errors propagate the kick, density and hand sampling at fixed coupling calibration. The finite envelope, boundary, momentum, flight time and changes in total wall number affect the measured response. The infinite-strand zero-momentum inverse curvature is $15/16$. The flight uses unheralded copies; measuring the endpoint to trigger a clock changes the state.

The source of motion, the read appointment and the material are commissioned. The free arrival supplies a real candidate initiation event. Identifying an electron requires separate particle and charge evidence. Every orientation carries one unit per complete dimer, so this stock measurement assigns zero excess to the wall relative to the background. A spacetime light speed common to other signals also requires an independent measurement. Exact arithmetic for these strands and the finite-resolution exclusion for the named readers do not decide the criterion on other Interval Echo sources. [Pinned decisions and official cell](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#results-by-model)

From the repository root, decode the archived counts and run the direct physical and algebraic checks with

~~~sh
python -m model.candidate_laws.interval_echo.openwave.row_03_clock.wall_arrival.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock/wall_arrival
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_03_clock/wall_arrival/WallBand.lean)
~~~

The acquisition option `--acquire` regenerates the simulated counts on the installed NumPy sampler. Portable verification reads the archived integers and checks all eight word probabilities against the direct orientation calculation on a small strand.
