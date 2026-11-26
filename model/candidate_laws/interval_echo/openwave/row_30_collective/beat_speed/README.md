# How fast do neighboring matter waves move?

Two adjacent wavelengths of a freely released neutral island produce a fringe that moves around its material ring. On the ten-dimer ring, increasing the outward momentum first raises the measured fringe speed from about 0.01504 to 0.03944 dimer steps per exchange time. Near the zone edge it falls from about 0.03960 to 0.01515. The same reversal occurs on the twelve-dimer ring, in both directions, and when the contrast is halved during flight. **Every choice of continuum Klein–Gordon mass, speed and energy origin fails on the full measured momentum ray in the dimer-step ruler.** The first increments need a separate experimental range to establish how accurately they approximate a continuum dispersion near rest. [Protocol frozen at 7ff0474f0](PROTOCOL.md) · [Complete counts](counts.json.gz) · [Velocities and bounds](results.json)

| Ring dimers | Flight contrast | Minimum increase between first two speeds | Minimum decrease between last two speeds |
| ---: | ---: | ---: | ---: |
| 10 | 8, native | 0.02437 | 0.02444 |
| 12 | 8, native | 0.02195 | 0.02204 |
| 10 | 4, control | 0.04761 | 0.04810 |
| 12 | 4, control | 0.04278 | 0.04346 |

Each entry rounds down the smaller lower bound across the two launch directions, in dimer steps per exchange time. The material uses the same native bonds on the five- and six-owner rings of the 80-cell board. A conserved sharp stock commissions the dimer array and its connectors. A local mixing and phase service prepares each pair of carriers before the contacts evolve freely. Every terminal census keeps all 4,194,304 simulated starts, counting individual single-upper words and one inclusive complement. Native contrast is eight. The four-contrast comparison changes the law during the flight. [Material ancestry](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/README.md) · [Addressed source circuit](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/source.py)

## What the patches measure

Separate fresh preparations interfere each source character with the uniform character at relative phases zero and π/2. The two all-start spatial histograms determine its signed Fourier index. The uniform-only control is flat to the registered tolerance. The Nyquist character is the alternating sequence and has just one sampled sign. Every release with two adjacent characters has a resolved spatial beat at their measured difference. Observers read each fringe again on fresh copies at fourteen nonzero exchange times, ending at 8192. All 44 velocity settings returned a unique clock winding.

Let the single-upper return of the measured character $a$ be $A_a(t)$, and let $k_a$ denote its spatial wave number. Translation covariance gives the single-upper probability at dimer $x$ as

$$p_x(t)=\frac{|e^{ik_ax} A_a(t)+e^{ik_bx} A_b(t)|^2}{2n}.$$

The nonconstant Fourier component at $k_b-k_a$ therefore has phase $\arg(A_a^* A_b)$. Its change in time divided by the spatial step gives the speed of the fringe, equal in magnitude to the energy secant on these adjacent momenta. A component of principal weight $w>1/2$ stays within $\arcsin((1-w)/w)$ of its line's phase. The read combines the calculated bound for each of the two sharp-word preparations with simultaneous 99% Hoeffding bounds on the measured complex moments. It uses the existing polar material's broad initial recoil range of ±1 exchange unit to begin the phase unwrapping. An independent test evolves every orientation word for both zone-edge directions at each size and contrast and compares the terminal probabilities with the archived counts. [Acquisition and decoder](measure.py) · [Direct word-space checks](test_speed.py)

## Why the curvature decides the comparison

A positive continuum Klein–Gordon energy is the Euclidean length of $(m,ck)$. By the triangle inequality, twice the length at a midpoint cannot exceed the sum of the lengths at its two neighbors. On equally spaced measured wave numbers, each successive energy increment must therefore be at least the preceding one. A fixed energy offset cancels in each increment. A linear drift adds equally to successive increments and cancels in their difference. The negative branch reverses the inequality. The positive low-momentum and negative high-momentum curvatures on each ray eliminate both branches even if each ring and direction has its own calibration. [Lean proof of the adjacent-secant inequality](Curvature.lean)

The previous island experiment used a measured rest and first recoil to predict fifteen further increments. The present sign test uses neither calibration. It determines the speed directly from neighboring nonzero carriers at the zone edge. Folding a different momentum coordinate into the spatial ruler changes the hypothesis being tested. The material's near-rest inertial response and the separate lattice Klein–Gordon dispersion of a commissioned wall are compatible with this measurement. [Earlier calibrated polar test](../polar_recoil/README.md)

## Independent wall-matter comparison

On the open 29- and 49-dimer wall strands, existing adjacent-wall phase, clock and local displacement counts give 28 finite-packet velocities over both native and halved contrast. The first two measured velocities have positive order and the last two have negative order on all four size/material panels. Independent local coupling records supply the mass and speed for the ordinary continuum group derivative $c^2k/\sqrt{m^2+c^2k^2}$. Averaging it over each prepared finite packet's Fourier power on the principal momentum interval gives agreement with **9/12** measured velocity intervals for $k\le0.50$, and **0/16** at greater momenta. The corresponding lattice forecast agrees on **12/12** and **15/16**. These retrospective coverage counts use point calibrations and the finite-width packet model. The confidence bands belong to the packet instrument. [Original measurements and uncertainties](../packet_speed/README.md) · [Reconstruction from counts](wall.py)

OW-30 retains its finite demonstration and partial full-criterion grade. The two polar rings exclude a full-band continuum assignment. The wall result identifies a low-momentum window for a different native preparation. The official strand submission has a separate admission and future-draw requirement. [Row decision](../../submission/ROW_30_DECISION.md)

Run from the RER root with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_30_collective.beat_speed
.venv/bin/python -m "$p.measure"
.venv/bin/python -m "$p.wall"
.venv/bin/python -m pytest -q --pyargs "$p.test_speed"
(cd Lean && lake env lean ../model/candidate_laws/interval_echo/openwave/row_30_collective/beat_speed/Curvature.lean)
~~~
