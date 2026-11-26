# A slow signed clock after a single occupied write

In the Interval Echo simulation, one sharp load enters a blank seam. Its native echo writes a positive direction. After release, with every screen contact open, the writing owner's local memory almost completely reverses and then recovers its sign. On the 320-cell board its calculated positive probability falls from 1 to 0.005588 at tick 200 and returns to 0.993868 at tick 400. A separate exchange dimer supplies one tick per half-turn. The signed clock's fitted angular frequency is 0.005108 radians per exchange unit. The smaller board gives 0.005109. [Frozen protocol](PROTOCOL.md) · [Complete calculated probabilities and seeded samples](records.json.gz)

The source poses a direct test of one proposed de Broglie interpretation. On fresh copies a component instrument gives mean formation energy 42.25 relative to the empty blank host. For every single unit load with one naturally writable seam, including coherent combinations of both direction parities, the full formation spectrum lies in [40,63] on each registered board. Every single-copy internal angular frequency is consequently at most 23. The observed freely running sign therefore cannot tick at the formation energy for any preparation or fixed terminal reader in this class, in units with ℏ = 1. **Full OW-03 remains Partial locally.** The source has a prepared initial write and a stationary memory location. Identifying a particle rest energy requires independently measured inertia on the same excitation. The upstream M10 entry remains untested. [Native sector proof](PROOF.md) · [Official criterion and M10 assessment](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L437)

## What the observer can read

The enrolled contacts are (0,1) on 20 cells and (2,3) on 320. A sharp amount-one load at the first endpoint and one native echo prepare a definite positive-direction memory. Releasing it leaves the mobile load entangled with the recording coin. The observer supplies a separate shuttered exchange dimer, counts its complete endpoint alternations and reads the original core's signed direction at the registered appointment. Each cohort starts fresh and has no intermediate direction read. The timer's 401 sharp endpoint reports and its calibration at zero, quarter cycle and half cycle are in the archive. The programmed half-turn exposure, choice of seam, one preparatory echo and terminal record access are supplied. The slow frequency arises from propagation under the all-open native contact law. [Timer, estimator and record sampler](measure.py)

| Cells | Quantity | Tick 0 | Tick 200 | Tick 400 |
| ---: | --- | ---: | ---: | ---: |
| 20 | Native positive probability | 1 | 0.005599 | 0.993845 |
| 20 | Native sampled positive counts | 1,048,576 | 5,896 | 1,042,227 |
| 20 | Zero-contrast positive probability | 1 | 0.476495 | 0.565515 |
| 320 | Native positive probability | 1 | 0.005588 | 0.993868 |
| 320 | Native sampled positive counts | 1,048,576 | 5,683 | 1,042,168 |
| 320 | Zero-contrast positive probability | 1 | 0.447630 | 0.518096 |

Counts are PCG64-seeded multinomial samples from the calculated probabilities, with 1,048,576 simulated starts per setting. Both directions and the impossible complement are in every bank. Reversing the incoming endpoint swaps the entire positive and negative laws. An equal mixture of the two sharp directions gives exactly one half throughout under both flight laws. A separately planted sign-flip Hamiltonian produces the unit-frequency reference response. The simultaneous 99% Hoeffding allowance of 0.002274 per bounded frequency describes Monte Carlo sampling resolution. Replay checks the dyadic sampling weights against the fresh native-law calculation and reproduces every sampled count. [Simulated cohorts](records.json.gz)

The earlier parity-resolved orbital spectra informed the preregistered observation window and search band [0.003,0.008]. The direction time series was reserved. The fit used ticks 0, 80, 160, 240 and 400, allowing a constant offset and both sinusoidal quadratures. The four held ticks were 40, 120, 200 and 320. On 20 cells the fitted frequency was 0.00510886 and the largest held probability error was 0.001083. On 320 cells the corresponding values were 0.00510821 and 0.001196. Both pass the frozen error threshold of 0.025. Sampling at these separated appointments alone admits frequency aliases. An independent full-operator audit gives ground-to-ground differences of 0.005105404 and 0.005105447. The two ground states contribute a signed cosine coefficient above 0.9948, and the remaining terms together have absolute amplitude below 0.148 at every time. This ensures a persistent line on each finite board while allowing weaker beats. [Operator and detector derivation](PROOF.md) · [Prior parity-specific spectra](../../row_31_orbitals/occupied_write/README.md)

## Energy on the same source

Separate fresh cohorts read the energy of each contributing actual seam. An inverse echo, a magnitude read and a reversal-parity comparison resolve the local generator's energy. Each private contact reads the squared endpoint load. A complete empty-host comparison gives zero. Summed component outcomes give sampled means 42.2614 and 42.2566 on 20 and 320 cells respectively, with simultaneous sampler allowance ±0.2024. The exact summed mean of the calculated instrument probabilities is 42.25 at both sites. This measurement uses no direction-clock fit. The upper bound of 23 on internal frequencies holds for sharp or mixed states, both recording coins, any fixed accessible effect and any core seam. A separately stocked magnitude-two collar lies outside the single natural-write domain. An independent motional mass and a physical unit calibration remain open. [Source and class-wide proof](PROOF.md) · [Finite support inequality](InternalBand.lean)

From the repository root, use one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_03_clock.released_direction
nice -n 15 .venv/bin/python -B -O -W error -m "$p.measure"
nice -n 15 .venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock/released_direction
(cd Lean && nice -n 15 lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_03_clock/released_direction/InternalBand.lean)
```
