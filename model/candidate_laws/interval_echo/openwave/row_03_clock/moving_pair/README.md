# Does a moving pair keep time across a cell boundary?

Two different loads start at opposite ends of a seam. The clock hand is which end holds the larger load. On an isolated seam this hand flips back and forth at one radian per exchange-time unit. After every contact opens, observers in the next cell can still read the phase set at launch when both loads reach their detector together. On all three tested bridges the arriving hand drifts away from the isolated clock. The visits are rare, with peak rates of **0.48, 0.49 and 1.25 per million launches**. [Acquisition and follow-up design](PROTOCOL.md) · [Complete count banks](records.json)

## What arrives?

The sender starts with one amount-one and one amount-two port in a blank, empty 20-cell sphere. Running the source seam alone for one of four phases sets its hand. The supplied contacts and their shutters prepare the phase. Nothing drives the 57,360 ordered positions during flight. The receiver reads both local loads across the named seam in its own cell, counting the prescribed order, the reverse order and all misses. This measurement can establish simultaneous arrival with a phase signal. It does not record the paths between launch and arrival. On this preparation the energy bound below keeps the two amounts adjacent with probability at least **27/32** at every free exposure.

The table gives the fractional arrival rate per original launch. The phase drift is measured between two separately reacquired appointments relative to a fresh isolated-bond clock. A constant offset on the route cancels.

| Bridge | Sender and receiver cells | Peak intact arrivals per million | Held appointments | Phase drift on repeat launches |
| ---: | --- | ---: | --- | ---: |
| 0 | 0 to 6 | 0.482 | 34 and 44 | −1.516 radians |
| 10 | 5 to 15 | 0.486 | 34 and 44 | −1.537 radians |
| 20 | 10 to 15 | 1.250 | 42 and 44 | +1.603 radians |

Let f and r be the amplitudes to collect the original and reversed orders after a sharp launch. Swapping the two species everywhere commutes with the relative generator. The even and odd transfer amplitudes are therefore a = f + r and b = f − r. For a launch exchange phase φ, the two collection probabilities are |a ± b exp(−iφ)|²/4. Their sum A is independent of the phase. Their signed difference D depends on it. The four phase settings give the linear statistic z = (D₀ − Dπ + i(Dπ/₂ − D₃π/₂))/2, whose expectation is a* b. It includes every miss. Thus nonzero z requires both order parities to reach the very same receiving seam coherently. A random mixture of opposite source phases erases it. An isolated-bond cohort measures the reference phase with the same four scans; a separate native exchange timer measures exposure in 32 equal subdivisions.

The quarter-turn and three-quarter-turn inputs have identical initial probabilities for both sharp load orders. Their relative coherences have opposite signs. The difference of their two remote order imbalances equals twice the imaginary part of z. The fresh counts resolve this difference at every one of the six enrolled appointments. Its magnitude ranges from 0.132 to 0.380 per million original launches. Thus a classical mixture retaining only the starting load order loses a measured part of the signal.

The first scan of 49 appointments used 2³² simulated launches at each phase setting. Each bridge resolved a nonzero hand. Every phase-erased cohort was compatible with zero. The zero-contrast control retained all exchange paths and reached the same receivers more often, peaking at 0.134% to 0.159% per original launch. The phase signal alone consequently identifies no contrast-bound particle. The first scan's simultaneous uncertainty was too large to resolve the native phase-only discrepancy. Its selected appointments were then fixed for an independent repeat with 2³⁸ simulated launches per phase. The independent drift at each bridge excludes even a fixed phase offset with arbitrary visibility at every appointment. A freely varying phase delay or a different physical clock remains possible.

For the test, multiply the measured z at each appointment by the conjugate of the independent bench hand. A common phase offset forces these complex numbers onto a fixed ray. For two measured values u,v, the imaginary part of u times the conjugate of v must vanish even if both signs on the ray are allowed. On the fresh counts the absolute imaginary parts exceed their simultaneous probability-error allowances by factors 6.29, 6.27 and 5.19. Each raw bin count has an exact two-sided binomial interval with tail probability 10⁻¹⁵ per side. Union and triangle bounds cover the real and imaginary combinations, with total counting-error probability below 2 × 10⁻¹¹ across the recorded marginal intervals. The finite bank and the procedure for selecting the follow-up times are public. All counts are simulated.

A local impulse on the expanded pair gives a flow predicted by an independent kinetic pointer. The pair retains at least 84.31% adjacency probability at every freely chosen time after the impulse. [Recoil, clock and binding measurements](recoil/README.md)

Each reported confidence bound uses direct finite-binomial tails with outward rounding. At bridge 10, appointment 45, the reverse-quarter preparation gave 999 reversed orders in 2³² starts, with a two-sided marginal interval [1.7887904 × 10⁻⁷, 2.9620525 × 10⁻⁷] at tail probability 10⁻¹⁵ per side. A separate recurrence of binomial term ratios verifies both rare-event limits. The independently repeated phase-drift tests use 2³⁸ starts per phase. [Tail inversion](counting.py) · [All count bounds and rare-event recurrence](test_counting.py)

## Independent energy and the scope of the clock

The sharp readings on the two amounts, the previously measured uniform single-load energies 48 and 192, and a separate full-exchange parity measurement fix the initial mean energy. The two adjacent amounts subtract 32 units from their separate uniform energies, giving 208. A uniformly chosen one of their ten unequal bonds has exchange-odd probability one half on each of the four preparations. Its controlled full swap reads this parity against idle. Summing the bond energies adds five. The twelve independent count banks give **212.988 to 213.010**, against the exact 213. These are empty-host formation energies in native units. [Two-load sector, native energy reference and binding experiment](../../row_31_orbitals/README.md)

In this sector the all-open generator referred to the separate uniform loads is H = K − 32P, with K half the Laplacian of legal swaps and P the adjacent-pair projector. Since K is positive, the conserved mean −27 implies that the adjacent probability is at least 27/32. Each configuration has at most twelve unequal seams. The Laplacian row-sum bound gives 0 ≤ K ≤ 12, hence every energy relative to the empty host lies between 208 and 252. This uses only the degree bound shared by the 20-, 80- and 320-cell boards. Every fixed terminal effect on this entire two-amount sector has Bohr frequencies of at most 44. An internal beat at the measured host-referenced mean near 213 is therefore excluded for these preparations. The same input could have a separately defined inertial rest energy. Determining it requires a dispersion or acceleration measurement on that state with a calibrated ruler.

The propagator uses the full ordered configuration graph. Its Chebyshev series omits at most 2 × 10⁻¹⁵ in operator norm in exact arithmetic, by summing the factorial bounds on the Bessel coefficients. The stated counting intervals apply to simulated probabilities. A separate test assembles all microscopic seam columns on an isolated cell and checks a full-board arrival against an independent matrix exponential.

From the repository root, replay the records and their independent checks with

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.measure
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock/moving_pair
```
