# Does the winding qubit travel at one rest energy?

On the source with one odd seam, the receiving owner can read a qubit at one fixed material port. An exact energy test on this same sixteen-state reachable component gives a different result. Its integer characteristic polynomial has no repeated roots. For each energy $E$, let $P_E=|u_E\rangle\langle u_E|$ and let $V$ embed the source code. Selecting this line maps any input density $\rho$ to
$$
P_EV\rho V^*P_E=\langle u_E|V\rho V^*|u_E\rangle P_E.
$$
The input changes the probability. Given this outcome the material occupies a fixed state, so even a receiving owner with arbitrary instruments cannot obtain an unknown qubit. Forgetting the line outcome mixes these measure-and-prepare channels and cannot restore entanglement. The receipt gives a nonzero arrival contribution from every line. For comparison, keeping all eight positive formation energies **coherently** gives a negative partially transposed channel eigenvalue of −0.0076433. This positive-subspace projection is a calculation; the finite-wait test below uses native free evolution and local detectors.

On each fresh launch, wait a random time uniformly between zero and $W$, discard that setting and let the charge fly for the original sixteen units. The preparation, arrival herald and three local analyzing pulses are identical to the fixed-port experiment. In the energy basis the forgotten wait multiplies an off-diagonal with separation $d$ by $e^{-idW/2}\operatorname{sinc}(dW/2)$, with $\operatorname{sinc}(x)=\sin(x)/x$. The simulated banks retain every success, both analyzer outputs, and every failed arrival.

| Width $W$ | Mean fixed-port arrival | Smallest eigenvalue of the partially transposed Choi matrix |
|---:|---:|---:|
| 0 | 0.338820 | −0.008108 |
| 1 | 0.380030 | −0.004009 |
| 2 | 0.403195 | +0.005407 |
| 8 | 0.336834 | +0.032739 |

The simultaneous 99% error in the last column is at most 0.002230 from 16,777,216 original starts per setting and 72 settings. At widths zero and one the accepted channel can carry entanglement. At widths two and eight the accepted two-qubit Choi matrix has positive partial transpose, which is sufficient for separability in dimension two by two. Providing the failure flag adds only a classical outcome. After an ideal sharp energy measurement with the outcome forgotten, the exact acceptance is 0.326623 and the minimum is +0.055480. The finite windows tested are individual appointments, without an assumption of monotonicity at intermediate widths. The observed qubit uses coherence between distinct formation energies in this source.

[Executable acquisition](rest_read.py) · [All energy lines and counts](rest_read.json) · [Native source and local analyzers](measure.py)
