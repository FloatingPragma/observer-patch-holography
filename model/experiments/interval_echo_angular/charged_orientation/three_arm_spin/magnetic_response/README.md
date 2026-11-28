# Stationary bond currents on the arriving winding lines

On each of the three physical guides, observers can identify the lowest extended negative vacancy line and compare the two outputs of a local current meter. For the isolated winding doublet, the signed difference is **exactly zero on every active contact**. The measured winding phase can distinguish the two circular inputs even though neither sources stationary circulation through these bonds. A head-dependent phase control on the same chart changes the measured port probabilities. [Study terms](PROTOCOL.md) · [Complete detector probabilities](results.json)

## Finite implication from the pulse records

Let $A$ be the vacancy transfer operator from an oriented seam's source port to its destination, including every admissible head flip. A local observer exposes only that contact. The native two-rail blocks have equal diagonal energies and off-diagonal magnitudes $a$. Apply exposures $\pi/(4a)$ and $3\pi/(4a)$ on fresh copies. Their difference in destination probability is the expectation of $J/a$, where $J=i(A^\dagger-A)$. Blocked inputs are counted at their actual positions and cancel between the two banks. The source and destination counts, elsewhere, and both preparation refusal categories enter the primary probability law. The calculation verifies the two effects directly on the blank and written contact columns. [Native pulse check and likelihoods](study.py)

The real free generator preserves the three-arm permutation group. For a seam in one arm, exchange the other two arms. This preserves both the oriented seam transfer and its clock, including all sixty-four head words. Within one isolated real standard doublet this reflection has one even and one odd direction. The current matrix is purely imaginary and skew after division by $i$. Its off-diagonal entries vanish under commutation with the reflection, and both diagonal entries vanish by skew symmetry. The conclusion holds for either circular label, mixtures, and a retained reference. The standard-row spectral enclosure excludes another copy at the same line in the prepared symmetry sector. For observers who verify admission to that line, this proves the null current on every seam without a further current measurement. [Matrix lemma](CurrentReflection.lean) · [Original source and clock selection](../recoil_search/THRESHOLD.md)

| Physical links per arm | Separate bond meters | Largest conditional score at the four tested waits | Certified native ceiling at every free wait | Largest change under the matched phase control |
| ---: | ---: | ---: | ---: | ---: |
| 19 | 63 | $1.203\times10^{-9}$ | $2.152\times10^{-6}$ | 0.001423 |
| 20 | 66 | $8.204\times10^{-10}$ | $2.339\times10^{-6}$ | 0.001242 |
| 21 | 69 | $9.191\times10^{-10}$ | $2.552\times10^{-6}$ | 0.001089 |

The evaluated native waits are zero, one, four and sixteen exchange exposures. The control comparisons use waits one, four and sixteen. The control phase depends on the local triangle bit and follows the same rule on all three physical arms, preserving their permutation symmetry. Both experiments use the original prepared input and native meter. For each finite-clock circular input the cyclic symmetry and the exact balance across each longitudinal cut imply the ceiling in the table at all free waits. The bound includes leakage outside the isolated line. The corresponding ceilings on the written contacts are $2.017\times10^{-7}$, $2.193\times10^{-7}$ and $2.392\times10^{-7}$. [Finite-clock derivation and replay](FINITE_CLOCK.md)

## Which applied phases can these observers distinguish?

A head-independent hop phase on a marked edge can be absorbed into the phase convention of its flipped bit. Blank edges act at fixed head word. Their graph is three paths between common hubs. Removing vertex phases leaves two independent through-arm phase differences. The receipt verifies the identity in integer coefficient arithmetic on every allowed transition, in both directions by conjugation. This is a phase-reference equivalence with the source and coherent instruments transformed together.

The two independent blank-route interferometers give unit bright probability and one-half quadrature probability at zero inserted phase. A supplied $\pi/8$ phase gives 0.96194 and 0.69134 on each route comparison. Both leave the full head record unchanged. In the changed flight, a blank loop comparing two left head bits of opposite parity has a genuine relative phase of $\pi/4$ in magnitude. The control thus exercises a physically detectable phase and the local current decoder. It supplies no autonomous magnetic source.

**OpenWave disposition.** Row 5 stays on HOLD. Any static dipole construction sourcing its field only from the stationary bond currents in this guide would have zero source on the tested doublet. A separate field degree of freedom would require its own coupling and readout on the same branch. The finite measurements supply no electromagnetic moment calibration, angular-momentum scale, or Larmor response for an electron.

From the OPH root, replay the source and the independent checks.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.magnetic_response
.venv/bin/python -B -m "$pkg.study"
.venv/bin/python -B -m pytest -q --pyargs "$pkg.test_study"
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/magnetic_response/CurrentReflection.lean)
```
