# How precisely must the observer time the charged qubit's flight?

The two tip orbitals still transmit entanglement when the observer adds an unrecorded waiting time drawn uniformly from 0 to 2 native units before the established 131,072-unit flight. With a width of 64, about **30.38%** of mixed starts still arrive, but the accepted reference–orbital state's negativity is at most **0.00111**, compared with at least **0.02832** with synchronized timing. At width 64 the dominant, much slower charge-transfer beat retains **99.999993%** of its visibility. Thus timing can resolve the orbital encoding and the passage of charge quite differently. [Recorded detector banks](waiting.json) · [Physical source and receiver](../README.md)

The new experiment uses the matched even head parity. On every fresh preparation, draw $u$ uniformly in $[0,W]$, evolve the native material freely for $u$, then run the usual flight for $T=131072$. The controller supplies timing; the receiver is given neither $u$ nor an energy outcome. All eight final head-coin words contribute to the receiving channel. On six prepared source states and three Pauli readout axes, retain both coincidence outcomes and every failed arrival. The simulator samples $2^{30}$ starts per setting from the exact marginalized outcome law. This is physical free waiting with a forgotten schedule. It requires no band projector or global energy-meter implementation. [Executable](waiting.py)

If $K_{j,c}$ is the arriving two-tip amplitude at energy $E_j$ for final coin word $c$, the averaged Choi matrix includes each pair of energy amplitudes with the factor
$$
 e^{-i(E_j-E_k)(T+W/2)}
 \frac{\sin[(E_j-E_k)W/2]}{(E_j-E_k)W/2}.
$$
The factor is one at zero energy difference before the flight phase is included. Taking the expectation before sampling detector outcomes performs the same average as physically drawing and subsequently erasing the schedule for each trial. An independent test integrates explicit evolution of the full native sector at positive quadrature times for $W=1/4$. [Spectral amplitudes and gap certificate](README.md) · [Independent test](test_waiting.py)

| Timing width $W$ | Measured mixed arrival | Measured smallest partial-transpose eigenvalue | Dominant charge beat visibility |
| --- | ---: | ---: | ---: |
| 0 | 0.421166 | −0.0122623 | 1 |
| 0.25 | 0.364286 | −0.0049715 | >0.999999999 |
| 1 | 0.293203 | −0.0014392 | >0.999999999 |
| 2 | 0.293540 | −0.0006309 | >0.999999999 |
| 4 | 0.302834 | −0.0003145 | >0.999999999 |
| 64 | 0.303753 | −0.0000015 | 0.999999934 |
| 4096 | 0.307515 | −0.0000034 | 0.999729 |
| 319147.19 | 0.165028 | −0.0000001 | 0.0000362 |

The last width is a full period inferred from the independent detector frequency in the previous experiment. The joint 99% error allowance is **0.000333** on the eigenvalue and **0.000064** on the mixed arrival probability, across all reported tomography settings and controls. The witness remains below zero through width 2. The width-4 interval contains zero; it does not measure an exact transition. For this two-qubit reference and receiver, a partial transpose has at most one negative eigenvalue. Indeed any two-dimensional subspace of two-qubit vectors contains a product vector: the homogeneous determinant of its two generating coefficient matrices has a projective root over the complex numbers. A partial transpose has nonnegative expectation on every product vector, which excludes a two-dimensional strictly negative subspace. Its absolute magnitude divided by the successful arrival probability gives the negativity conditioned on arrival. At width 2 the resulting simultaneous interval is 0.00101–0.00329, and at width 64 it is 0–0.00111. These intervals characterize the accepted Choi state formed from an initially maximally entangled pair.

A synchronized control uses the same mean wait as the width-64 experiment, exactly 32 extra units. Its measured arrival is 0.409288 and its partial-transpose eigenvalue is −0.011554. The accepted negativity is at least 0.02741 there. The difference isolates the effect of timing dispersion from a shift of the mean appointment. At zero total flight time the receiver has zero arrivals. Independent unconditional charge measurements on the same physical preparation give the signed current through the bridge. [Local current and beat measurements](README.md)

In the infinite-width limit every distinct-energy cross term vanishes. The eigenvalues on the reachable component are all simple; each retained energy has just one input bra, even when all coin words are included. The limiting channel therefore separates the input from its reference exactly. Its mixed arrival is 0.395053. Large finite waits need not approach this value monotonically, and the smallest splittings are much slower than the main charge beat. A finite-resolution observer can use the synchronized trimer as a coherent orbital carrier under the settings measured here. Neither that success nor the timing scale supplies a common degenerate energy or a physical spin-½ rotation.

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.live_trimer.energy.waiting
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/live_trimer/energy/test_waiting.py
```
