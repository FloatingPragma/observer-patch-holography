# Local return after a disturbance of the cold ring

On stations 2532 and 2533 the symmetric and antisymmetric single-load counts measure the exchange bond energy. A phase reversal at 2532 raises this energy from \(-0.318319\) to \(+0.318319\). After 96 identical propagation steps, the same detector reads \(-0.311108\) in expectation. A synthetic bank of 262,144 launches per setting gives \(-0.318121,+0.318907,-0.311129\) for the reference, disturbed and late readings. Every empty and doubly occupied outcome is included. [Finite record and counts](results.json) · [Replay](measure.py)

The background is the prepared half-filled sea on the 246-station ring. After the initial local intervention the fixed native program repeats the same one-step unitary \(U=\exp(-iH)\), where \(H=-\tfrac12\sum X_{j,j+1}\). Its exact group law permits a terminal read at fractional action time as well. The propagation uses the neighboring stations, phases and routing of the addressed preparation. It has no feedback or cooling during flight. [Action, source construction and all-open controls](../README.md)

**Uniform local-return theorem.** Let an unconditional quantum channel act on a fixed finite set \(S\) of ring stations, with Kraus operators preserving its total load. Start from any stationary background at the chosen filling. For the occupation of any normalized mode \(f\), the change after action time \(t\) obeys
\[
\left|\delta\langle n(f)\rangle_t\right|\le 2\eta,\qquad
\eta=\|P_S e^{iHt}f\|.
\]
Here \(H\) in the one-particle expression denotes the hopping matrix. For \(\eta<1\) project the evolved mode onto the complement of \(S\) and normalize to obtain \(g\). Its occupation commutes with all local Kraus operators, so its expectation is preserved by the intervention. The difference between the two one-particle rank-one projectors has eigenvalues \(\pm\eta\). On fermion Fock space their occupation difference also has norm \(\eta\), by diagonalizing the two modes. Applying this bound before and after the channel proves the inequality. At \(\eta=1\) the probability range suffices. This proof permits local measurements with their outcomes discarded, arbitrary local phases and entangling operations on the enrolled patch.

On an infinite chain the propagator at fixed separation \(d\) is \(i^d J_d(t)=O_d(t^{-1/2})\), by the Fourier integral of \(\exp(it\cos k)\). Hence every fixed finite local intervention returns the expectations of all fixed local number-preserving quadratic observables to their background values. Polarization gives the off-diagonal quadratic observables from occupations. Finite-ring measurements use the following bound before a disturbance can circle the apparatus. [Bessel integral representations](https://dlmf.nist.gov/10.9) · [Large-argument asymptotics](https://dlmf.nist.gov/10.17)

Take \(S\) to be both receiving stations and \(f_\pm\) the symmetric and antisymmetric modes. On the infinite chain both norms \(\eta_\pm\) equal \((J_0(t)^2+J_1(t)^2)^{1/2}\). The derivative of the sum of squares is \(-2J_1(t)^2/t\). On the finite ring the diagonal and nearest-bond propagators include images at separations of magnitude at least 245. For \(t\le96\) each entry's image sum has absolute bound
\[
\epsilon=\frac{2(48)^{245}}{245!}\frac1{1-48/246}<5.79\times10^{-69}.
\]
This follows by summing \(|J_m(t)|\le(t/2)^{|m|}/|m|!\), obtained from the Poisson integral for integer order. Thus \(\eta_\pm\le(J_0(s)^2+J_1(s)^2)^{1/2}+2\epsilon\) throughout \([s,96]\). The bond energy is half the difference of the two occupations. Rational alternating-series enclosures at \(s=16,32,64\) give the following absolute limits for **every** charge-preserving intervention on the two stations and every time in each window.

| Action-time window | Largest allowed bond-energy change |
| --- | ---: |
| 16 through 96 | 0.393757974 |
| 32 through 96 | 0.281231502 |
| 64 through 96 | 0.200011116 |

A sharper exact calculation covers any channel on the one binary station \(s\) preserving its load. Such a channel multiplies its off-diagonal coherence by a scalar \(z\) in the closed unit disk. Write \(C\) for the Fermi projection, \(u=e^{-itH}|s\rangle\), \(v=e^{-itH}(C-I/2)|s\rangle\), and let \(h\) be the measured bond. For one choice of phase convention the energy change is \(2\operatorname{Re}[(z-1)\langle v|h|u\rangle]\). At ages \(0,16,64,96\) the maximum absolute changes over this entire disk are \(0.636637,\ 0.043824,\ 0.010356,\ 0.007210\). The phase reversal attains the maxima at those appointments. A direct calculation in the six-station fixed-load word basis independently checks the formula. [Local controls and fixed-load test](test_relaxation.py)

The cold input is supplied away from the disturbance. A uniform mixture over the same fixed-filling sector stays stationary at bond energy zero under the same evolution. The pulse train protects the ring during propagation. Releasing every federation contact gives a different physical generator and causes immediate leakage into the exterior. The independent all-open two-load experiment also has stationary interior entropy biases at identical total-energy distributions. Local return establishes stability of a specified field background against confined disturbances. Selection of that background from an extensive hot preparation, the proper-area response and the Einstein premises require additional inputs. [All-open stationary records](../../autonomous_selection/README.md)
