# Rotation-selected port current

Interval Echo can optionally equip one complete cell with an addressed port pulse. Its bound pair responds to twelve independent currents on the specified two-band instrument. At native time 16, switching the sign of the pulse at port 0 changes the distribution of the patch’s ordinary seam readings by **0.111596 in total variation**. The receiver reads all thirty possible occupied seams and quiet. The calculation gives predicted probabilities for fresh preparations. [Pulse and detector](currents.py) · [Every outcome and control](currents.json)

## How a current is selected

The service uses the oriented incidence of the shared echosahedron. A unit at port $p$ supplies the one-load source. The corresponding two-load source places a pair uniformly on the five edges touching $p$. In both cases the prepared excitation is the normalized difference between this source and its antipodal image. For each energy present in that source, the detector asks whether both units occupy one native seam. It retains levels with conditional agreement exceeding one half. The one-load control admits every occupied port. At $\kappa=8$ the two-load source reaches four levels. The two retained triplets have agreement probabilities 0.991773 and 0.998393. The energies use the contrast law's sector-dependent choice of zero. [Native pair Hamiltonian](../contrast_memory_spherical_v1/sectors.py) · [Source spectrum](currents.json)

Sum the permutation matrices for the proper rotations which turn the neighbours of a fixed port by one edge around its pentagon. Twelve rotations meet that incidence test. Let $P$ be the positive spectral projector of their sum. The mean class-sum eigenvalues on the selected triplets are 6.472136 and −2.472136. For each port the controller uses the band-diagonal part of $P n_p P+J_p$, where $n_p$ measures occupation and $J_p$ circulates a unit along the oriented ring of five neighbours. Gain is one. The two-load version lifts both one-load operators to the hard-core pair sector. The ideal service exponentiates that operator in the selected modes and acts as the identity on the complement. The local preparations, ring orientation, protected spectral service, gain, signed pulse strength and exposure are supplied settings. The implementation uses the native rotations and derives both bands from the prepared source. [Shared rotation action](../../executable_model/symmetry.py) · [Operator construction](currents.py)

## What the control distinguishes

| Tested menu | Independent input directions | Lie closure | Visible after removing a shared phase |
| --- | ---: | ---: | ---: |
| Rotation-selected port services | 12 | 12 | 12 |
| Same services plus the band-diagonal occupation at port 0 | 13 | 18 | 17 |

The results agree for the one-load control and the two-load bound excitation. On each sector all 66 brackets lie in the narrow span and 720 port-rotation conjugation checks pass, with residual norms below $10^{-11}$. These matrix checks establish covariance on one cell. The wider menu quantifies how a separately permitted physical operation changes the closure. Completeness over all accessible controls, accepted intercell rotations under refinement and the gauge theorem's matter premises remain open. [Restricted and extended menu calculations](currents.json)

The detector experiment starts the odd two-load source at port 1, applies $+0.41$ or $-0.41$ at each of the twelve ports, and then freely evolves the intact cell for 16 native time units. The whole input, including its weak unselected components, enters each terminal reading. Across the twelve addresses the total-variation distance ranges from 0.005372 to 0.111596. All 24 signed distributions are archived. This isolated-cell exposure has an externally selected duration. [Predicted patch records](currents.json)

The earlier workspace gate compiler built the unit-gain two-load pulses on all 66 basis inputs from native exchanges, echo, addressable written-contrast phases and protected control windows. Its separate synthesis receipt is `LAW/current-interval-clock/theorem-calibration/search/current/sensitivity.json` in the research workspace, SHA-256 `68d4d36f7227624e11b81c0df548510599fe032dee248edfaa66ad83762f446e`. The replay here executes the ideal guarded operator and its terminal detector. Call `pulse(port_number, strength, loads=2)` to obtain that unitary on the complete 66-state binary pair sector. The one-load version uses `loads=1`.

The one-load colour triplet also constrains propagation. A sharp gateway has exactly one-quarter weight in this band. On all three tested boards every everywhere-singlet single-load trajectory has zero probability current between cells. A colour filter at just the source and receiver permits coherent orbital transfer of signed defects under the specified extension of the one-load meter. [Colour and orbital transport](../../experiments/interval_echo_gauge/singlet_flight/README.md)

From the RER repository root, print the full receipt with the first command. The second reruns the independent detector comparison, including the return pulse and a separate construction of the two-load Hamiltonian.

~~~sh
python -m model.candidate_laws.interval_echo.currents
python -m pytest -q model/candidate_laws/interval_echo/test_currents.py
~~~
