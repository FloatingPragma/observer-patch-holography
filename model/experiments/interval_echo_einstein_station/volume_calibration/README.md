# Measuring a region by inward and outward crossings

A reciprocal crossing experiment measures how much of the accessible launch space lies inside the home cell. That fraction stays fixed for all four independent perturbations of the two matter registers that preserve both marginals. The aperture and launch volume fail a separate three-dimensional calibration on the blank board. Each refinement quadruples the launch volume of the same addressed parent region and doubles its short-hand aperture. In a similarity-preserving three-dimensional calibration, quadrupling volume would multiply area by about 2.52. [Protocol](DESIGN.md) · [Recorded probabilities and simulated counts](results.json)

| Blank screen | Launch volume of parent | Short-hand aperture | Relative A³/V² in short-hand limit | Relative A³/V² at hand 2 |
| ---: | ---: | ---: | ---: | ---: |
| 20 cells | 12 | 3 | 1 | 1 |
| 80 cells | 48 | 6 | 0.500 | 0.502 |
| 320 cells | 192 | 12 | 0.250 | 0.252 |

The shortest appointment, at hand 0.5, gives relative ratios 1, 0.499992 and 0.249999. One isotropic conversion factor per board cancels out of this comparison. The parent addresses and uniform launch protocol are apparatus inputs. The result excludes interpreting this pair of observables as the volume and surface area of geometrically similar three-dimensional regions with common microscopic calibration. A radial or other degree of freedom absent from this measurement could define a different bulk ruler.

## What the inward preparation measures

With $P$ the projector on the region's sharp addresses, $Q=I-P$, and $U$ the native one-load flight, the outward probability is $x=\|QUP\|_F^2/N_R$. The inward preparation gives $y=\|PUQ\|_F^2/N_C$. Unitarity makes the two numerators equal because each is $N_R-\|PUP\|_F^2$. Consequently $y/(x+y)=N_R/(N_R+N_C)$ whenever there is a crossing. The finite volume result follows exactly from unitarity within the uniform source ensemble. Matter can change the crossing probability without changing this volume. The occupied home has 9 reachable sharp addresses and the three complete mobile components contain 237, 957 and 3837, respectively.

An observer asks whether a sharp load started and ended on opposite sides and retains quiet outcomes. Separate simulated binomial draws for inward and outward preparations give the stored finite-sample volume fractions. The smallest board includes the full complementary population from an independent dense flight, for the blank board and all nine occupied words. On the two larger boards the full inward population uses unitarity. Independently propagated single exterior preparations check the reciprocal predictions for a direct neighbour and an exterior port one additional hop away. The aperture and its derivatives use complete endpoint probabilities. The material hand uses the station's two complementary endpoint analyses. [Acquisition and checks](measure.py)

For finite hand reading $t$, the calibrated aperture is $A_t=4N_R x/t^2$. Expanding $QUP=-itQHP+O(t^2)$ gives $A_0=4\|QHP\|_F^2$. A blank exchange of magnitude one half contributes one unit. With three-dimensional length unit $h$, the putative area and volume would be proportional to $h^2 A$ and $h^3 V$, respectively. Their cubed-area to squared-volume ratio is independent of $h$. Counts on the blank descendants instead satisfy $A_0^2/V=3/4$ on all three boards, consistent with counting an interface on a sheet. The test leaves open a different choice of spatial regions or geometric observables.

## Regional entropy at fixed launch volume

The resident, hand, protected tag and two persistent collision registers use the established occupied-station preparation. The reference is the nine-word distribution with first row $(4,6,2)/36$ and cyclic subsequent rows. Each tangent moves probability among four joint words while preserving both marginal distributions and hence the same mean stored energy. The four independent tangent entropy slopes are −0.007991, −0.049771, +0.041780 and −0.007991 nats per unit displacement. Independent joint-word counts measure these responses with approximately 0.000275 nats standard error per tangent. [Occupied source preparation](../refinement.py) · [Earlier clock-distance comparison](../regional_geometry/README.md)

All nine word Hamiltonians differ on diagonal contact energies only. Their short-hand apertures are exactly one, so the short-hand area derivative vanishes for every preparation. For finite appointments the four-vector of area derivatives has norm of order $10^{-10}$ to $10^{-9}$. At each appointment, fit a single unrestricted coefficient in $\delta S+\eta\,\delta A$ to all four directions. The relative residuals in the declared basis are 99.43%, 99.25% and 99.40% at hand 0.5, 1 and 2, respectively, on each board. At hand 2 the best unrestricted coefficient is negative. The independent dense complementary flight on 20 cells reproduces all four finite aperture gradients within $1.4\times10^{-17}$. The tiny area derivatives are calculations from the complete simulated record laws. Finite-shot area measurements at the deposited sample sizes do not resolve them. Restoring the occupied energies makes all collision words yield the same aperture. [Population and finite-sample results](results.json)

An exploratory follow-up calibrated a single marginal-preserving source tangent to leave all three finite apertures unchanged to first order on 20 cells. Its predicted entropy slope is **+0.00765437 nats per unit displacement**. Fresh simulated memory acquisitions gave +0.00752 ± 0.000275, +0.00793 ± 0.000275 and +0.00741 ± 0.000275 on 20, 80 and 320 cells. With the original tangent the largest transferred aperture slope had magnitude $8.9\times10^{-18}$, below the numerical discrepancy of the smallest board's independent flight implementations. The smaller board's vanishing finite slopes are fitted. This tangent distinguishes matter entropy from the three declared transmission apertures on the tested family. G4 still requires a calibrated proper volume, a physical area, and equilibrium of the appropriate total entropy at an independently identified vacuum. [Follow-up design](FOLLOWUP.md)

## Reproduce

Run both commands from the repository root in its Python environment.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_einstein_station.volume_calibration.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_einstein_station/volume_calibration
```
