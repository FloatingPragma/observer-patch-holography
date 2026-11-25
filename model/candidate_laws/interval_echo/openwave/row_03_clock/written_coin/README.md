# A neutral written seam ticks at its independently measured excitation gap

An observer can start a clock by reading one coin on a written seam and then opening every contact. The same coin changes at **31.995** or **33.001** radians per exchange-clock unit, according to the prepared contrast parity. A second detector on fresh copies of the identical preparation resolves the upper excitation line at **32.007** or **32.990** above the filled host. The preparation occupies that line on half its launches. Its measured mean formation energy is **15.982** or **16.467**, respectively. The result establishes a free internal gap clock at a fixed address. The mass of a particle moving through the federation requires its own identification. [Design](PROTOCOL.md) · [Records](records.json)

## What the observers do

On the 20-, 80- and 320-cell boards the selected seam joins ports 0 and 5 of cell zero. Supply a filled host with one excess and one vacancy installed at that seam. Half exchange followed by a local reference phase on port 5 prepares the even or odd superposition of the two orientations. The native echo converts both orientations to filled endpoints and writes their signed contrasts on the same collar. Read just the coin, preserving coherence between contrast signs. Either outcome starts the clock. A separate comparison uses a sharp, unsplit excess and vacancy, whose written contrast sign is known. Installing the neutral pair, addressing the preparation gates and supplying the coarse coin read are part of the laboratory. The clock runs without a driving pulse after release.

At each of 193 appointments the observer compares the final coin with its initially recorded coin. It retains equal, flipped and blank results. A separate native dimer supplies swap counts over one eighth of the unknown exposure. Unfolding those counts gives the observer's time in dimer-gap units over the enrolled range. The prepared subdivisions and shutters are apparatus. A frequency fitted to the first 97 appointments forecasts the last 96. The three board embeddings share this exact invariant local sector and hence one count acquisition per setting. They are separate geometry checks of a single clock law. [Preparation and acquisition](measure.py)

| Free contrast | Source parity | Frequency from coin | Held RMS error | Upper energy bin | Mean formation energy |
| ---: | --- | ---: | ---: | ---: | ---: |
| 8 | even | 31.9952 | 0.0140 | 32.0067 | 15.9818 |
| 8 | odd | 33.0013 | 0.0128 | 32.9900 | 16.4671 |
| 4 | even | 16.0004 | 0.0064 | 15.9052 | 8.0026 |
| 4 | odd | 16.9995 | 0.0054 | 17.0847 | 8.4850 |

Each coin appointment has 1,048,576 launches, including both possible starting coins. An independent timing bank has 67,108,864 starts per appointment. The 512-bin energy pointer has spacing about 0.196 in these units. At native contrast the component meter on the even-parity source returned **524,884 zero-energy and 523,692 occupied even outcomes**. On the odd-parity source it returned **525,332 zero-energy and 523,244 occupied odd outcomes**. Each difference between the sampled mean and half the upper energy is within the sampling uncertainty. A centered equal-weight coin has exact mean energy 16 or 16.5.

The sharp, unsplit source produces both rapid lines. A single fitted frequency misses its held native records with probability RMS error 0.268. Erasing the coherence between energy branches gives exactly one-half flip probability at every appointment. Removing contrast during flight stops the even clock and leaves the odd clock ringing at 1.0000 in exchange units. A blank, uniformly filled source always gives the blank outcome. Empty or leakage-only spectral windows on these controls have no interpretation as detected upper lines.

## Can the two clocks distinguish a material change from a change in elapsed time?

A common change in elapsed time relative to the independent reference multiplies the two clock rates by the same factor. A change in contrast strength shifts these two material lines by the same **amount**. The first three readings, the two native frequencies and the changed even frequency, suffice to predict the changed odd line under each hypothesis. The common factor predicts **16.480–16.527**, while an equal shift predicts **16.968–17.046**. The other 96 appointments on the odd clock give **16.996–17.004**. The clock construction distinguishes these effects under this intervention. The intervention changes a coefficient in the law. A gravity experiment instead compares clocks at different locations in the same law.

This is a retrospective comparison of archived count banks. The three reference intervals use their first 97 appointments only. Each coin reading uses all starts, and each duration comes from the companion dimer count on its calibrated first timing branch. The calculation inverts both binomial readings, allows every oscillation branch below the original search limit of 46 and takes an outer bound on each frequency. The joint intervals have at least 99% coverage for the specified four panels under the free, sharp-coin response and the reference-clock model. A synthetic control with a common factor of one half admits the common-factor forecast. [Count analysis](rate_comparison.py) · [Control tests](test_rate_comparison.py)

## Why the fast line survives on an open board

Let the two signed contrast labels at magnitude two be the second two-state system, and let the first system be the coin. Write X for sign reversal and P for the projector onto the antisymmetric coin. With both port loads equal to one, the native Hamiltonian on this four-dimensional written subspace is exactly

$$H=P\otimes[4\kappa I+(I-X)/2].$$

The blank filled state and the symmetric-coin state both have zero energy relative to the filled host. To see the formula, the two oriented inputs of the echo write into the symmetric coin with opposite contrast labels. A second echo sends that part back to blank memory, where the memory generator vanishes. The antisymmetric coin is orthogonal to every such write and the echo fixes it. Its native generator measures contrast square four and reverses the sign. Every neighboring blank seam has equal, unit-filled endpoints and annihilates these states. All degree-reference contacts add the same scalar on the packet and the filled reference. The proof therefore covers free all-open flight with every collar initially blank except the one prepared mark. [Exact native columns and independent count checks](check.py)

In a fixed parity, a sharp coin carries equal amplitudes on the symmetric zero-energy state and the antisymmetric excited state. Write E = 4κ for even contrast parity or E = 4κ+1 for odd. The probability of reading the opposite coin after free exposure t is exactly sin²(Et/2). A native echo and a sharp coin read prepare both amplitudes; the archived reference-phase exclusion applies to a different, zero-contrast starting subspace. Here contrast magnitude two is a physical stored record.

Two energy instruments use separate cohorts and consume neither the coin frequencies nor their fits. The Fourier instrument coherently controls how long the active local term runs against a parked branch and reads the register in 512 Fourier bins. A separate dimer measures the clock increment. Idle and active branches have the same background and static degree-reference contribution. The component instrument reverses the echo, then reads whether the memory is blank or has magnitude two and which reversal parity it carries. A controlled full exchange against idle measures the latter. The supplied pulse switch, phase standard and known local energy coefficients enter these measurements. The upper energy is the cost of the antisymmetric-coin branch relative to the zero-energy host. Both branches and every unsuccessful comparison remain in the denominator. The study uses exchange-clock units with reduced Planck constant one.

## When can the mean equal the frequency?

Within either parity, write p for the excited-state population of the released collar. Its mean formation energy is pE. Any off-diagonal coherence in the energy basis has magnitude at most the square root of p(1−p). Thus setting the mean equal to E forces p=1, which gives constant coin probabilities. Source signs, starting coins and any initial density operator on that parity obey this constraint.

Allowing both parities can manufacture a numerical coincidence. At κ=8, mix one part freely ringing even sharp coin with sixteen parts stationary odd upper state. The resulting mean is (16+16×33)/17 = 32. Its coin-flip oscillation at 32 has only one seventeenth of the original contrast, with an independent initial coin label for the stationary cohort. With fixed mixing weights and κ=4 the mean is (8+16×17)/17 = 280/17 while the fast even line is 16. This is an arithmetic mixture chosen after seeing the energies. Preparing the odd upper cohort without destroying its stationarity requires its own energy selection; a fresh sharp coin measurement on it would change the state. For any fixed state with populations p₀, p₊, p₋, the mean is 4κ(p₊+p₋)+p₋. Equality with either fast line on an interval of κ forces full weight on just that line. It then supplies no fast internal beat.

The location of this clock is fixed by the same invariant columns. For any chosen seam, its written four-dimensional subspace stays at that seam and every other blank, filled seam acts as zero relative to the host. The free matrix elements between clocks written on distinct seams therefore vanish, including seams sharing an endpoint. Superposing their locations cannot give this clock a propagating dispersion. Relocating it would require an additional interaction or a different preparation. Its gap alone supplies no independently measured inertial rest mass. OW-03 has a local internal-frequency measurement on this specified source. The upstream row on the pinned model is untested. [Official criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L103)

From the RER root, reproduce the count banks and the separate exact polynomial and distribution checks with

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.written_coin.measure
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.written_coin.check
```
