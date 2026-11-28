# A charge and its compensator on the written two-fan fixture

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

At the 14–24 gateway, a vacuum-referenced interferometer reads a low vacancy line at **3.548015** in measured ring-clock units. The first observer predicts the line for a separated moving vacancy and excess from this read. An independent detector and clock give **7.096021**, against **7.096031** predicted. A second held read tests the archived orbital inputs at the single-vacancy frequency. The lowest line supports only one source direction and causes the large, time-selected doublet arrival to disappear from almost all launches. These are simulated records on a finite prepared fixture. [Protocol](PROTOCOL.md) · [Counts and bounds](records.json)

## Neutral comparison and finite implication

The source uses the full archived written four-face fixture, with 160 configurations in each single-charge sector. In the pair arm an addressed feeder between two copies reverses one head, changing the unit loads $(1,1)$ to $(0,2)$. Exact checks on all 256 neutral input words preserve both endpoint head balances and total signed charge. The ordinary free written contact has zero matrix element for this conversion. The declared preparation pulse implements the transfer and inverse. The feeder is off during the free wait and is excluded from the archived degree potential. The vacuum arm keeps both copies unit-filled and the same local clock running.

For one vacancy read, park the other copy's positive charge at its degree-six source. Its neutral rings stay active in the vacuum and its occupied endpoint has zero degree potential. Let $s_-$ and $s_+$ be the two dressed gateway launch preparations and $H_-$ and $H_+$ their free Hamiltonians. The return contrast is measured on the reference qubit after both arms are recombined. For the vacuum launch, its signed complex contrast is

$$
C_{\mathrm{park}}(t)=\langle s_-|e^{-itH_-}|s_-\rangle.
$$

With the positive charge released on the second copy, the dwell operator factors because the feeder is closed. Complementing all head orientations is an exact permutation of the two free Hamiltonians and maps the vacuum launch to its opposite-sign counterpart. Thus

$$
C_{\mathrm{move}}(t)=\langle s_-|e^{-itH_-}|s_-\rangle\,
\langle s_+|e^{-itH_+}|s_+\rangle=C_{\mathrm{park}}(t)^2.
$$

The joint spectrum contains the sum of independently resolved contributing frequencies. Exact symmetric integer inertia gives counts $(0,1,1)$ below the cuts $(3.547,3.549,4.65)$, respectively. A nonzero ground line read on the parked panel therefore entails a twice-frequency line of nonzero weight on the moving panel, under the tested preparation, clock and closed-feeder hypotheses. The held detector checks this consequence with separate trials.

The neutral clock reads two signed interference quadratures of a local ring prepared in an equal coherent superposition of its vacuum and one odd mode. The clock detector uses the actual ring coherence. A second ring clock uses the other outer-ring parity on an independent shot stream. The simultaneous 99% interval for the observed pair-to-clock ratio minus the predicted ratio is **[−0.001325, +0.001309]**, within the fixed ±0.002 resolution. Removing the positive copy's hopping while retaining the complete contact diagonals yields a held line at **19.547992** in the same clock units. [Acquisition and replay](measure.py)

## Held arrival and the energy zero

A phase-calibrated coherent dwell filter is applied to the two exact archived odd launch preparations. It accepts about **0.01211%** of first-ring-odd launches and **1.09864%** of last-ring-odd launches. The receiver's gateway charge-and-local-doublet counts are below at the archived dwell 27.0625. Every row has $2^{22}$ launches, including every declined filter.

| Launch | Full release | Phase filter |
| --- | ---: | ---: |
| Vacuum | 171,752 | 6,846 |
| First outer ring odd | 445,474 | 5 |
| Last outer ring odd | 1,244,071 | 451 |

The numerical filter has a finite response away from the ground line. The exact isolated eigenspace has rank one, so its transfer matrix has a kernel on the two-dimensional archived source. No choice of waiting time makes an exact filter at this line transmit an arbitrary orbital input. The finite accepted vacuum launch gives nearly half its charge clicks in either cell, including in the read directly after filtering. These readings describe an energy-selected excitation delocalized over the two fans and a timed unfiltered orbital transfer on the same apparatus. Establishing a confined two-input branch on an extended freely active medium requires a test on that medium. [Exact parked-threshold comparison](../binding_audit/README.md)

A common energy term proportional to total signed charge cancels in the neutral comparison. A per-defect energy term is different. Adding $\mu=1/2$ to each occupied defect changes only a scalar within the flight blocks, so every fixed-defect detector probability is exactly invariant at all times. The same vacuum-to-pair Ramsey read moves by **0.999994** model frequency units because the pair arm carries two defects. Thus the Ramsey read measures an energy difference that ordinary charge-arrival records leave undetermined. The assignment of half the moving-pair gap to a common single-branch gap uses the declared charge-complement symmetry and the selected neutral reference. The ring transition supplies internal units. A model-to-laboratory clock identification and a localized freely propagating branch are missing from a 511 keV rest-energy inference.

From the repository root, reproduce all count banks and certificates with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_recoil/measure.py
```
