# A rotating material probe beside a persistent current

A vacancy circulates on a five-port ring in either of two degenerate stationary states. They have opposite currents and the same density at every port. On a second ring, a standing-wave probe can tell the source directions apart once the rings interact. A calibrated pulse schedule lets the probe precess while the circulating source retains its original mode with high probability. Reversing the material charge of the probe preserves the precession sign. The measured effect comes from the two-ring interaction. [Protocol, physical ports and controls](PROTOCOL.md) · [Compact source on the intact federation](../row_05_dipole/README.md)

## Which interaction transmits the direction?

The rings are the even and odd sites of the native ten-port icosahedral belt, with ten interring contacts. One vacancy occupies the first ring. The second has a vacancy or an excess. The guide holds the population of each ring fixed. Intraring exchange and the diagonal contact energies remain active. Projecting the native 99-state seam gives a connected energy of −17 for two vacancies on contacting ports and +15.5 for a vacancy beside an excess. The two signed source states give the same average energy against each sharp probe position.

The observer measures both quadratures of the probe's angular coherence, every source momentum and every probe miss. With all contacts active, reversing the source direction changes the sine quadrature by −1.56167 for the vacancy probe at time 9 and −1.72570 for the excess at time 12.875. The corresponding source-mode return probabilities are 78.94% and 89.22%. Both signals contain 15 distinct nonzero absolute frequencies above the stated numerical resolution. With a single contact, reflecting the source ring about the contacting port interchanges the sources and preserves the entire probe instrument: every probe reading is exactly identical for opposite source directions. With the source continuously held in one momentum, its density is uniform, so the two projected probe generators agree even with all ten contacts. [Complete spectra and geometry](results.json) · [Integer joint records](counts.npz)

## Precession with a small source disturbance

Use a standing probe in modes ±2 and issue separate source-rim, probe-rim and interring pulses of durations 1/4, 1/8, and 1/512 or 1/1024. All pulses use the same native contact. The timing, Fourier-mode preparation, population guards and coherent analyzers are supplied. Two contacts through intermediate pair configurations give a direction-sensitive phase. At both contact durations second-order perturbation from the isolated rings misses the 0.10-radian tolerance over a predicted turn. [Predictions and held tests](dispersive.json)

A second interferometer compares separately prepared two-ring product modes while the cycle acts. Reading their bare return phases calibrates the beat without using the standing-wave flight data. A fresh standing probe is then read after 20,000 to 96,000 cycles for both source signs. All 40 held rotations agree with the calibrated phase to 0.017 radians. The lowest observed probe visibility is 96.74%. For the shorter contact, the original source mode returns in at least 99.28% of the enrolled probes. The per-cycle spectral weights imply a conservative return bound of 98.15% throughout the sequence at cycle boundaries, including unsampled cycles. The two branches of the standing wave have distinct total fifth-turn momenta, so their source-return probabilities add.

| Probe sign | Contact duration | Measured beat per cycle, source +1 | Independent finite spectrum |
| --- | ---: | ---: | ---: |
| Vacancy | 1/512 | +0.00044282 | +0.00044281 |
| Excess | 1/512 | +0.00009339 | +0.00009335 |
| Vacancy | 1/1024 | +0.00008862 | +0.00008868 |
| Excess | 1/1024 | +0.00003814 | +0.00003812 |

Reversing source circulation reverses the precession. Holding its direction fixed and reversing the probe charge preserves the sign at both durations. This sign test distinguishes the material interaction from a coupling proportional to the product of the two signed orbital charges. Assigning an electron magnetic moment requires an independently sourced exterior field and mass and angular-momentum standards for that same charge. [Calibration and 40 held readings](beat.json)

From the repository root, regenerate the native joint reads, both weak-contact comparisons, and the independently reconstructed long-time propagators.

```sh
python -m model.candidate_laws.interval_echo.openwave.row_05_counter.experiment
python -m model.candidate_laws.interval_echo.openwave.row_05_counter.dispersive
python -m model.candidate_laws.interval_echo.openwave.row_05_counter.beat
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_05_counter
```
