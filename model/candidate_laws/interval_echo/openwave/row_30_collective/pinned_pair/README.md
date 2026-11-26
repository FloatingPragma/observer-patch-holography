# Which waves can the eleven-contact pair launch?

A patch puts a relative phase between the two occupied seam words of the admitted M10.4 preparation. Its own thirty seam detectors then read the freely released pair. Reversing the relative sign produces up to 729 times the probability of leaving the two source words. Every phase gives the same time profile after this change in weight. These phase impulses cannot select a family of dispersing matter waves on the admitted pair. The separately commissioned dimer strand has its own source and wavelength receivers. [Frozen experiment](PROTOCOL.md) · [Complete count banks](counts.json) · [Native acquisition](measure.py)

Write X for the pair on (156,157) and Z for the pair on (156,166). The three occupied preparation reads give the accepted amplitude proportional to −13X + 14Z. All eight unsuccessful attempts have total probability (364/729)⁸. A full exchange swaps the second tip onto the degree reference at 156. A delay on that reference followed by the same exchange gives the multiplier exp(−iφ) on Z. The two addressed operations use a separate preparation interval. During the wait all contacts of the prepared board evolve freely. The changed-law control sets the contrast to zero only during that wait. [Original native operations and receiver](../../row_04/README.md)

The antisymmetric word D = (X−Z)/√2 has exact energy 1 + 3κ for any real contrast κ. Its orthogonal partner S = (X+Z)/√2 can escape. The conditional bright weight after the pulse is

$$w(\phi)=(365-364\cos\phi)/730.$$

For any receiver outcome annihilating D, the conditional detection probability at time t is exactly w(φ) times the probability on an initial S. This follows by evolving the two components and applying the same detector to each. It also holds for sums of such outcomes. The full spectral measure of the kicked word is the mixture of a stationary line at 1 + 3κ and the spectral measure of S with that same weight w(φ). Both supports are independent of φ. Even an arbitrary diagonal phase pattern on the original ports contributes only the relative phase between the two occupied words. A phase gradient here therefore changes the amount of a common spectrum that a receiver sees. An actual matter-wave dispersion measurement requires independently distinguishable wavelengths with their own frequencies. [Integer eigenvector verification](../../row_04/topology_census.py) · [Lean receiver identity](PhaseRead.lean)

The acquisition starts on the original 240-port carrier and evolves all 28,680 distinct pair configurations. Each bank contains 2²⁷ simulated starts, with thirty source-cell seam addresses, accepted quiet and source exhaustion. The table gives the all-start probability of an accepted result outside X and Z. This aggregate includes other seam words in the source cell. It is a change of word and need not be an intercell arrival.

| Exchange time | No phase impulse | π impulse | π impulse, exchange-only flight |
| ---: | ---: | ---: | ---: |
| 2 | 0.0000054 | 0.003848 | 0.776083 |
| 3 | 0.0000251 | 0.018329 | 0.901672 |
| 5 | 0.0000132 | 0.009758 | 0.793135 |
| 8 | 0.0000338 | 0.024662 | 0.747156 |
| 13 | 0.0000180 | 0.013052 | 0.994470 |
| 21 | 0.0000206 | 0.015043 | 0.949434 |

At each wait the π bank predicted the other four banks. The largest residual across the 29 applicable individual bins and their union was 0.361 of the allowed simultaneous sampling error, rounded up. Dephasing the preparation removes the response to the kick at every wait. Individual detections on X and Z still allow interference with D. The phase pulse, station clock and eleven written contacts are specified parts of this test.

From the repository root, run the following commands with one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_30_collective.pinned_pair
.venv/bin/python -m "$pkg.measure"
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_collective/pinned_pair
```
