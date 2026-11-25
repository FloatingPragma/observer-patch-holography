# A lossless orbital sign on one charged face

A simulated local detector can send a signed defect to the dark interferometer exit by taking two coherent paths around a degeneracy of three material ports. Moving the control loop off the degeneracy sends the defect back to the bright exit. The same source gives positive phase under a complete spatial turn and under opposite orders of perpendicular half-turns. The dark exit identifies an orbital geometric phase under local scalar exposures. This distinction matters when using an interference sign to assign electron spin. [Derivation](PROOF.md) · [Source and acquisition](measure.py) · [Complete records](records.json)

| Intervention | Vacancy at the dark exit | Excess at the dark exit |
| --- | ---: | ---: |
| Enclose the degeneracy once | 99.816% | 99.912% |
| Displace the same loop outside the degeneracy | 0.194% | 0.043% |
| Enclose it twice | 1.699% | 0.014% |
| Park the direction of the local exposure | 0% | 0% |
| Rotate the complete cell template through 360° | 0% | 0% |
| Compare opposite orders of perpendicular template half-turns | 0% | 0% |
| Suppress contrast during the enclosing drive | 0% | 0% |

Each entry is a simulated frequency from 4,194,304 launches. Every launch yields an exit and an address for the vacancy or excess. The driven values use 1,024 steps per arm. All banks at 64, 256 and 1,024 steps, including leakage between orbital lines, enter the archive. The simultaneous 99% Hoeffding radius for each bin and dark marginal is 0.118 percentage points, rounded up. Finite-speed orbital transitions explain the difference from the adiabatic limit. The 99% dark criterion and the displaced-loop and spatial controls distinguish a parameter-space sign from a physical half-turn-order sign on this source. The host, laboratory controls and ideal readout define the observation domain. This experiment assigns neither electromagnetic charge units nor an absolute rest energy. [Tests and full-word replay](test_cycle.py)

## Preparation and readout

Start from the unit-filled, blank-collar 20-cell sphere. Select the first face for which all three ports have an actual degree-five private contrast reference. Its addresses are 0, 5 and 11. Change the middle port to zero or two, giving exactly one missing or additional unit. Isolate the face. A native half-exchange on the last two ports and a private exposure of π/16 on the last port produce their reflection-odd orbital with unit probability for either sign. The exchange, private reference and three open face contacts come from registered Interval Echo. The shutters, nonnegative private service durations, coherent arm, initial stock and balanced recombiner are supplied. Neither intermediate measurements nor conditional survival enter the cycle.

Label the ports by angles φ equal to 0, 2π/3 and 4π/3. In one step, expose each private reference for half its scheduled dwell, open all three face contacts together for h = 0.2, then complete the second half of each dwell. The total private dwell for port j is

\[
t_j=\frac{0.35h}{8}\bigl(3+\cos(\theta-\phi_j)+b\cos\phi_j\bigr).
\]

The total duration of private exposures is independent of θ and b. In the two arms use the opposite midpoint sequences θ = ±πw(k+1/2)/N, k = 0,…,N−1. The enclosing loop uses (b,w) = (0,1), the displaced loop (2,1), the double winding (0,2), and the parked drive (0,0). Recombine once. Each count records the exit and the full ternary occupation of the face, which fixes both the signed load and its address. In the separate rotation banks, both arms admit all twelve local addresses. A changed-drive comparison sets the common seam and private contrast strength to zero throughout the drive while retaining the prepared orbital, routing, durations and readers. Every such bank gives bright interference. A high enclosing dark fraction accompanied by an equally high displaced fraction would defeat the proposed winding discriminator.

The finite material words and computation of ideal probabilities are auditable separately from the simulated detector samples. A test evolves all 27 ternary words independently of the three-column acquisition. Another identifies the tangent splitting by exact symbolic arithmetic. Lean proves the mirror relation for histories of arbitrary length under the verified single-step intertwining hypothesis. A negative spatial half-turn-order sign and a matching exchange measurement on a freely mobile source would distinguish an intrinsic spin candidate. The native freely propagating orbital and the bound-trimer exchange have independent source records. [Mirror proof](../../../../../Lean/SourceLaw/OrbitalMirrorCycle.lean) · [Freely transported orbital](../README.md) · [Native exchange](../mobile_core_exchange/README.md)

## Reproduction

Run from the repository root with its NumPy, SciPy and SymPy environment. A single worker suffices.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.conical_cycle.measure --output /tmp/oph-orbital-cycle.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/conical_cycle
(cd Lean && LEAN_NUM_THREADS=1 lake env lean --threads=1 SourceLaw/OrbitalMirrorCycle.lean)
```
