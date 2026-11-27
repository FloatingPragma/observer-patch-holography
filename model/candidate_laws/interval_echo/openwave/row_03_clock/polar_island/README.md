# The polar island and its neighboring clock hand

On the 80-cell Interval Echo screen an observer can prepare a sharp upper dimer, release it around a closed ring, and read a neighboring exchange hand whenever a randomly addressed dimer reports upper. The independent rest-frequency hand reads **39.91317** at native contrast. An all-lower launch gives a stronger local rest line than either loaded launch. The ground level and the lowest mobile island are the unique pair of energies with exactly this separation anywhere in the complete ring orientation space. A local tick at that exact frequency therefore measures coherence with the host. OW-03 remains Partial. [Prospective design](PROTOCOL.md) · [Original-denominator records](records.json) · [Full-space enclosures](enclosures.json)

## What the observer measures

The ring uses the two existing 80-cell mounts. Static coherent connector collars and terminal apparatus are supplied. Each upper orientation reverses the two endpoint loads of one dimer. One sharp upper at site zero, a coherent distribution of one upper on two neighboring sites, and an all-lower sham evolve freely. The reference comes from two separate return instruments on sharp source stocks and the prior remote motion and recoil reads use this same moving band. The two-site control tests a change of launch profile. [Mount and reduced generator](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/PROOF.md) · [Independent rest and recoil counts](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/fringes.json)

At the terminal appointment the observer chooses a uniform dimer index. Its signed endpoint contrast resolves orientation. On the right neighbor an exchange-parity read resolves the symmetric and antisymmetric orientations. Equivalently, native addressed exchange and private-contrast quarter rotations followed by an orientation measurement realize that basis. The two dimers are disjoint. A classical random sign with expectation equal to the cosine of the reference phase labels each terminal record. This sign uses the known waiting time and supplies no excitation during flight. All four joint outcomes and both reference signs are retained. The reported signed contrast assigns zero to either lower-tag outcome and divides the signed upper counts by **all** launches. Thus a rare tag cannot inflate the measured amplitude. The delay average uses two independent integer draws and includes every sampled wait through the exact triangular probability law. The records are simulated multinomial draws from these probabilities. [Native local operation and read](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/PROOF.md) · [Acquisition](measure.py)

The last two controls remove only ground coherences or select the one-island spectral band. Both use analytical spectral operations in the comparison calculation. The sharp physical launches and all counts on their free runs precede that interpretation. Each of the 216 banks represents 2^38 independent simulated launches. The conservative standard error for the signed contrast is at most 1.91 parts per million.

| Dimers | Contrast | Independent rest hand | Island bandwidth | Upper launch, center (ppm) | Upper launch, ground dephased (ppm) | All-lower sham, center (ppm) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 8 | 39.913169 | 0.099323 | −0.03 | −0.28 | −22.37 |
| 12 | 8 | 39.913166 | 0.099323 | 0.75 | −1.25 | −22.39 |
| 10 | 4 | 19.830148 | 0.194811 | 9.62 | −5.93 | −170.67 |
| 12 | 4 | 19.830138 | 0.194811 | 6.66 | −7.68 | −170.48 |

The table gives the longest delay window. Its first zero occurs about 0.00153 angular units from the center. The ground-dephased sham is below 0.2 ppm in magnitude at all three probe frequencies. The loaded runs fail the frozen eight-standard-error center-versus-sides rule after dephasing on both lengths and contrasts. The exact spectral decomposition retains a small excited-sector response in the finite window, including −6.61 and −7.09 ppm for the upper launch at half contrast. Nearby many-body lines and finite-window leakage can contribute there. Within the isolated mobile band every address-averaged joint effect is time independent. These are distinct experimental and spectral conclusions. [Finite proofs](PROOF.md) · [Lean lemmas](PassiveRead.lean)

## Replay

From the RER root, choose one numerical thread. The enclosure command checks every energy, including pairs with unequal momentum. The independent small-ring test enumerates every physically allowed delay and both neighbor projectors directly in word space.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
lane=model.candidate_laws.interval_echo.openwave.row_03_clock.polar_island
folder=model/candidate_laws/interval_echo/openwave/row_03_clock/polar_island
.venv/bin/python -B -W error -m "$lane.measure" --output /tmp/polar-clock-records.json
.venv/bin/python -B -W error -m "$lane.exclusion" --output /tmp/polar-clock-enclosures.json
.venv/bin/python -m pytest -q "$folder"
(cd Lean && lake env lean -j1 "../$folder/PassiveRead.lean")
```

The [acquisition laws](acquisition_laws.json) hold the binary64 probabilities of all eight outcomes in every bank. The record test recomputes each law, requires agreement with the stored law within 10⁻¹², and redraws the seeded counts from the stored law. The spectral summaries are compared within binary64 roundoff. A fresh acquisition through the output option draws from the recomputed laws, and the laws option writes them.
