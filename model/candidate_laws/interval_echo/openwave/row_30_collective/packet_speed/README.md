# How fast does the wall travel?

The same material has a long-wavelength Klein–Gordon coefficient corresponding to **4 dimer spacings per exchange-clock unit** and a maximum group speed of **1**. On a 49-dimer strand, the fastest of seven tested phase gradients gave a signed wall-density speed of **0.9903** [0.9784, 1.0023]. Opposite preparations have the same initial wall probabilities. A receiver sees their density move in opposite directions. The speed is half the difference in their first-moment slopes, retaining every click and quiet reading. [Acquisition design](PROTOCOL.md) · [All 28 comparisons](results.json) · [Detector counts](counts.npz)

## Two ways to obtain a speed

Let J be the coupling between neighbors and h the dimer exchange. The meter in the same commissioned material gives J = 8.00029 and h = 0.500031 in exchange units. The positive band has squared frequency ω² = 4(J−h)² + 8Jh(1−cos k). Its low-momentum coefficient is c² = 4Jh, and its velocity is v(k) = 4Jh sin(k)/ω(k). For J ≥ h the exact identity

$$4h^2\omega^2-(4Jh\sin k)^2=16h^2(J\cos k-h)^2$$

bounds |v| by 2h with equality at cos k = h/J. Exchanging J and h covers the other ordering. Thus the maximum is 2 min(J,h), whereas the coefficient c is 2√(Jh). These quantities coincide at the gapless point. The maximum is a group speed for this bulk band. A causal cone for all signals in the federation requires a separate measurement. [Native reduction and coupling measurements](../README.md)

The experiment uses two open strands prepared before this acquisition on the 80- and 320-cell boards. Local services set the phase gradients of one-wall wave packets. Independent neighboring-wall quadratures determine their actual wave numbers. At four appointments the receiver reads individual wall counts on fresh copies. Two calibrated exchange-clock quadratures measure the flight durations in separate cohorts. The initial sources have finite width, the detector reads the bare wall density, and the strands have ends. The full finite-chain propagation includes all three effects. The confidence intervals invert the binomial counts, include the timer uncertainty and cover the enrolled readings simultaneously at 99% or greater under the specified simulation.

| Strand | Contrast during flight | Coefficient c | Bulk maximum | Largest packet speed in scan |
| --- | ---: | ---: | ---: | ---: |
| 29 dimers | 8 | 4.0002 | 1.0001 | 0.9862 |
| 49 dimers | 8 | 4.0002 | 1.0001 | 0.9903 |
| 29 dimers | 4 | 2.8286 | 1.0001 | 0.9791 |
| 49 dimers | 4 | 2.8286 | 1.0001 | 0.9830 |

The changed-law forecast uses the known halving of the contrast term with the native coupling-meter calibration. Across all 14 native velocities the bulk group-velocity forecast averaged over the declared packet envelope lies inside the counting intervals. In 13 of 14 changed-law comparisons it does likewise. At k = 1.5083 on the shorter changed strand the bulk forecast of 0.9904 exceeds the observed interval [0.9694, 0.9889]. This comparison uses a point forecast with the coupling estimates fixed. It tests the conversion from an extended band to a finite-packet detector.

## Response near rest

At zero momentum the exact band curvature is ∂²ω/∂k² = c²/ω(0). With reduced Planck constant one, its reciprocal is a kinetic mass in the strand's own distance and time units. The native values give 15/16 = 0.9375. The 49-dimer experiment reads k and velocity independently. At k = 0.1000 and 0.2000 the ratios k/v are **0.9478** and **0.9545**, respectively. These finite-width, finite-k responses are compatible with the band curvature at the precision of these counts. They concern a phase-prepared collective wall in commissioned material. A mass in laboratory units, a freely formed electron and a signal speed common to different materials each require their own identification.

## Does the packet keep its launch velocity?

An independent instrument reads the change in preparation energy under a phase gradient. Within the initially occupied one-wall subspace all diagonal energies agree. A flip of dimer j connects the wall at j to the wall at j+1. For a normalized, real envelope a, the prepared mean energy thus has the form E(k) = constant − 2h cos(k) Σ aⱼaⱼ₊₁. Reading local dimer flips on copies with gradients k+π/2 and k−π/2 gives its exact central derivative, [E(k+π/2)−E(k−π/2)]/2. The position moment is diagonal in wall words and preserves this subspace. Its initial Heisenberg velocity is consequently the same derivative. The full Hamiltonian couples the bare wall to other configurations after release.

For k = 0.800 on the 49-dimer strand, this initial velocity is in **[0.71161, 0.71192]**, whereas the four-appointment flight gives **[0.73042, 0.75278]** at native contrast. The initial and flight intervals are disjoint in 10 of 14 native settings and 12 of 14 half-contrast settings. The energy experiment was designed with the flight counts available. Independent local flip counts cover its 14 derivatives simultaneously at 99% by a bounded-sum inequality. Joint coverage with the separate 99% flight intervals is at least 98%. The finite-time wall response therefore contains dynamics beyond the prepared energy slope. [Energy read](work_meter.py) · [All comparisons](work.json) · [Energy counts](work_counts.npz)

From the repository root replay the acquisition and the small-chain physical check with

~~~sh
.venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_30_collective.packet_speed.measure
.venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_30_collective.packet_speed.work_meter
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_collective/packet_speed
~~~
