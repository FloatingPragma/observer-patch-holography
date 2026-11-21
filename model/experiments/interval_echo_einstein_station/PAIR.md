# The regional modular charge and the pair’s measured energy

A single gain calibrated on current port 0 gives the wrong sign for the central energy change on current port 9. This comparison uses the same four releases, two gateways, elapsed exposure and central energy meter as the finite energy-balance experiment. We used their equally weighted source mixture as the common modular reference. This choice keeps the modular reconstruction independent of the energy meter. [Protocol](PAIR_PROTOCOL.md) · [Four releases and numerical support checks](pair.json) · [Energy experiment](../interval_echo_null_balance/PROTOCOL.md)

| Signed current change at native time 4 | Regional modular charge difference, nats | Central energy difference, native units | Energy predicted using the port-0 gain |
| --- | ---: | ---: | ---: |
| Port 0, calibration | +0.00379039 | +0.00965558 | +0.00965558 |
| Port 9, held energy | −0.00459299 | +0.09979154 | −0.01170011 |

Both pulse signs at both addresses enter the unlabelled reference preparation. Port 9’s energy is held out of the gain calibration. In the four-way mixture σ, shifting a small weight from the minus to the plus release at a given port produces the tangent Δρ = ρ₊ − ρ₋. For |ε| ≤ 1/4 the state σ + εΔρ remains a mixture of the same physical preparations. The first-law slope on this support is tr(Δρ · (−log σ)). The independently measured energy of each preparation is the sum of the cell’s thirty actual internal seam energies during propagation. The ratios ΔE/ΔK have opposite signs across the two addresses. The calibration gives +2.54738 energy units per nat, and the held prediction misses by −0.111492 energy units. Keeping only the diagonal local word probabilities also fails and changes the slope of the calibration itself. [Signed differences and diagonal-only control](pair.json)

**What the patch could read.** On the 20-cell sphere, the bound pair propagates through cells 0, 6 and 4. At the terminal appointment a reader on cells 0, 1, 2 and 5 measures the state of cell 0 using its previously closed third gateway. It distinguishes every zero-, one- and two-load word on the twelve ports. Fresh copies receive 4,422 alternative real or imaginary two-word interference settings and one direct word setting. We rebuilt the native configuration gates and the terminal readout paths. With exact outcome probabilities the reader recovered all four reduced density matrices to operator norm error below 2 × 10⁻¹⁶. The reference has rank 17 on its 79-word space with sector ranks 1, 12 and 4. The smallest positive eigenvalue is 6.88 × 10⁻⁹. All four states lie inside its support up to numerical error below 3.4 × 10⁻¹⁴ in Frobenius norm. The three registered spectral thresholds give the same ranks and slopes.

**What the finite survey resolved.** With 262,144 simulated preparations per setting for the common reference, unconstrained linear tomography returned 36 negative eigenvalues, the smallest −0.00461. Its logarithm gives no physical modular generator at that exposure. The signed charge comparisons above use exact local outcome probabilities. The finite count archive retains every outcome, including quiet. A positivity-constrained estimator or a different measurement budget would need its own uncertainty and near-null support analysis. [All finite reference counts](pair_reads.npz)

G2 requires the *geometrically normalized* cap modular charge and directional stress of propagating matter on one common reference. The result rules out equating this particular regional source-mixture charge with the ordinary central energy by a single gain on these two tangents. The independent energy meter sums the unweighted internal terms. A geometrically weighted stress integral, a physical vacuum, cap shape and modular boost normalization remain open. The two current ports select different perturbations of one pair, and have no assigned timelike observer directions. [Einstein premises](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

From the RER root the first command verifies the source, apparatus probabilities and finite draws. The second also rebuilds the native tomography compiler in the enclosing LAW research workspace.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.pair
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.pair --verify-compiler
```
