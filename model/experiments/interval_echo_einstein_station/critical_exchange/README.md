# Interval entropy and local energy on one Interval Echo ring

Owners of consecutive stations can infer an individual bond's signed energy change by varying the boundaries of the region on which they measure entropy. A second reader measures the same bond by interference. On a 246-station ring the two measurements give both signs of a coherent energy redistribution. Across three source variations, two target bonds and four choices of the opposite endpoint, their calculated Born responses differ by at most 0.744 percent. The ring realizes a prepared one-dimensional field experiment with the active contact law. A direct free-release test distinguishes this source from an equilibrium of the unaddressed federation. [Protocol](PROTOCOL.md) · [Executable and controls](measure.py) · [Born rates and synthetic counts](results.json)

## Independent action on the native source

On binary loads and a blank seam, the full 99-dimensional Interval Echo term restricts to
\[
h=(\kappa+\tfrac12)P_{\ne}-\tfrac12 X,\qquad
X=|10\rangle\langle01|+|01\rangle\langle10|,\qquad \kappa=8.
\]
The available one-port parity \(Z=(-1)^n\) preserves the first term and reverses \(X\). The two terms commute. The entire blank contact has integer spectrum. With \(g(t)=\exp(-ith)\), the exact physical word
\[
g(t)Zg(2\pi-t)Z=\exp(itX)
\]
holds for \(0\le t\le2\pi\), including both equal-load words and any spectator reference. Native duration-\(\pi\) swaps transport data to the degree-five parity service and restore the transit ports. The registered phase writer supplies arbitrary one-port phases. Together with the exchange word these generate every adjacent two-mode rotation. Removing exchange at the same timing gives zero transfer. Halving the contrast changes the parity service and halves the test's transfer probability. [Port-phase apparatus](../../interval_echo_gauge/port_phases/README.md) · [Closed-route selection](../../../candidate_laws/interval_echo/openwave/row_30/experiment.py)

In the ordered binary occupation basis a nearest-neighbor \(X\) is fermion hopping. The closing bond acquires the usual sign \((-1)^{M-1}\) at total filling \(M\). The independently specified local energy is \(h_j=-X_{j-1,j}/2\). The sum \(H=\sum_j h_j\) has one-particle dispersion \(-\cos k\) and Fermi velocity one in station and action-time units. A Fourier rotation, one phase per mode and the inverse rotation implement \(\exp(-itH)\) exactly on this fixed-filling sector. The readout of \(h_j\) is half the probability of an antisymmetric single-load result minus half the probability of a symmetric result. Empty and double outcomes contribute zero. An addressed balanced splitter and a port phase resolve the two alternatives.

The source has \(M=123\) and occupies the negative-energy orbitals. Their preparation from a sharp filling uses adjacent rotations and phases. The executable replays 30,135 rotations on actual neighboring stations along the route and verifies every output orbital to \(7.5\times10^{-15}\). An occupation-mode rotation followed by native load detection measures each modal probability. The quadratic description computes the output rates of these instructions on arbitrary superpositions at this filling. The logical speed uses periods of the programmed action. Its relation to moving material clocks in three spatial directions is a separate physical identification.

## Exact common modes across overlapping intervals

Let \(N\) be the even ring length, \(a\) a left boundary, \(b=a+L\), and \(\alpha=\pi/N\). Station \(a\) is the first enclosed site and \(b-1\) the last. The Fermi projection restricted to these sites has entries
\[
c_d=\frac{\sin(\pi d/2)}{N\sin(\alpha d)},\qquad c_0=\tfrac12,\qquad C_{ij}=c_{i-j}.
\]
This follows by summing the occupied plane waves. For \(L\le N/2\), the occupied and empty Fourier submatrices each have full row rank by the Vandermonde determinant. All eigenvalues of \(C\) therefore lie strictly between zero and one. The vacuum on the interval factors into independent occupied/vacant alternatives in its correlation eigenbasis, by the Slater determinant rule. Its one-particle modular generator is \(K=\log[(I-C)C^{-1}]\). For a normalized differentiable source variation, \(S'=\operatorname{tr}(K C')\). This identity also applies to a mixture variation whose nearby states are not Gaussian, since \(-\log\rho_I\) at the reference is quadratic.

Weight the independently specified bonds with the circle's conformal interval profile,
\[
f_L(s)=\frac{\sin(\alpha s)\sin(\alpha(L-s))}{\alpha\sin(\alpha L)},\qquad
D_I=2\pi\sum_{s=1}^{L-1} f_L(s)\,h_{a+s}.
\]
The factor uses the velocity fixed by the dispersion. On one-particle amplitudes the off-diagonal bond entries are \(t_s=-\pi f_L(s)\), with \(t_0=t_L=0\). With site indices relative to \(a\) and \(d=i-j\),
\[
[D_I,C]_{ij}=(t_i-t_{j+1})c_{d-1}+(t_{i+1}-t_j)c_{d+1}.
\]
The two differences have the same prefactor times \(\sin(\alpha(d-1))\) and \(\sin(\alpha(d+1))\), respectively. The remaining sum is proportional to \(\sin(\pi(d-1)/2)+\sin(\pi(d+1)/2)=0\). For a zero index the original difference is zero. Thus \([D_I,K]=0\) for every contiguous interval in the faithful range, including translated and nested intervals. A cut outside the interval gives the same argument on a ring with a parity twist. The tridiagonal matrix has simple spectrum because its interior bonds are nonzero. Every modular mode can consequently be read in the basis specified by the action. This applies the free-ring commuting operator to a source and meters built from Interval Echo. [Eisler and Peschel, finite chains](https://arxiv.org/abs/1805.00078)

Commutation determines eigenvectors. At four stations the ratios of the two positive eigenvalues of \(K\) to the corresponding eigenvalues of \(D_I\) are \(0.987411\) and \(1.014914\). For a direct certificate, in the reflection-even basis the two matrices \(C,D_I\) have rows \((1/2+c_3,c_1),(c_1,1/2+c_1)\) and \((0,t_1),(t_1,t_2)\), respectively. Directed interval arithmetic bounds the two ratios in the disjoint intervals \((0.987,0.988)\) and \((1.014,1.016)\). A common multiplier cannot give the finite modular operator on every perturbation. The interval of 24 includes a reference mode probability of \(1.32\times10^{-17}\). A separate 45-decimal diagonalization of the Fermi covariance checks its rare probability and its entropy response. [Tests and enclosures](test_exchange.py)

## Signed regional and local tests

For each reference mode \(u\) count its occupancy \(p\) and vacancy \(1-p\). Its entropy contribution is \(\log[(1-p)/p]\,\delta p\). On different copies measure the symmetric and antisymmetric bond outcomes. The sources are a plane-wave replacement, a coherent two-wave replacement, and a reversible rotation from an occupied into an empty orbital. For the rotation, the interference phase places an energy node between the two target bonds.

The common bond energy can also be reconstructed entirely from neighboring apertures. Write \(D(a,b)\) for the weighted expectation change and \(e_b\) for the energy change of bond \((b-1,b)\). Put \(F(a,b)=\sin(\alpha(b-a))D(a,b)/(2N)\). The sine recurrence cancels every interior bond in
\[
F(a,b+1)-2\cos\alpha\,F(a,b)+F(a,b-1)=
\sin\alpha\,\sin(\alpha(b-a))\, e_b.
\]
The exact energy inverse is independent of the chosen left endpoint. Replacing the three energy responses by the three entropy derivatives provides a local prediction with no bond reading as input. Comparing it with the independently read \(e_b\) gives the errors below. Each column reports the maximum over all declared centres or target bonds and all three sources.

| Aperture length | Regional entropy vs. weighted energy | Local energy inferred from entropy vs. direct bond |
| ---: | ---: | ---: |
| 4 | 0.1940% | 0.7440% |
| 8 | 0.1432% | 0.1572% |
| 16 | 0.0403% | 0.0569% |
| 24 | 0.0192% | 0.0370% |

The reversible rotation gives \(+1.05888\times10^{-5}\) on bond 51 and \(-1.05888\times10^{-5}\) on bond 55. All four left endpoints recover both signs. For the eight integrated apertures the computed action/covariance commutator norm is below \(1.4\times10^{-13}\). Uniform weights with the same sum give a norm between 0.277 and 1.938. At four sites an independent calculation using the full 16-state reduced density matrices checks the entropy derivative of the replacement mixture.

The separate synthetic count bank uses \(2^{30}\) copies per setting on four stations. It gives an entropy response of 0.00439 with simultaneous 99 percent interval \([0.00292,0.00587]\), and a bond-energy prediction of 0.00442 with interval \([0.00240,0.00644]\). Both measurements resolve a positive response. These count intervals include the finite difference between the exact Born responses, 0.00441186 and 0.00440332. The tighter entries in the table use calculated Born probabilities.

## Free release and the Einstein premises

The full, unaddressed, blank-memory federation supplies a direct control. Every one of the 975 seams touching the prepared ring has disagreement probability at least \(1/2\). For a seam ending outside, this is the mean occupation of the ring endpoint. For two ring sites it is \(1/2+2|C_{ij}|^2\) by the Slater rule. The positive native exchange and degree-reference terms give the full free energy a lower bound \(8(975)/2=3900\). Filling ports 0 through 122 at the same total load instead gives 32 unequal seams and degree deficit 92. Its exact free expectation is \(8.5(32)+8(92)=1008\). The prepared sea therefore lies above the free fixed-filling ground energy. Exterior port 710 has exactly one ring neighbor, 697. The only first-order route into its initially empty load slot has amplitude operator \(-\sigma^-_{697}/2\). Its receiver probability begins at \(t^2/8+O(t^3)\) in native all-open time. This control tests the physical distinction between the addressed propagation and a passive release.

For the continuum massless chiral field there is a corresponding operator identity. On a circle set \(z(x)=\log[\sin(\alpha(x-a))/\sin(\alpha(b-x))]\), so \(z'=1/f\). Under the unitary half-density change from \(x\) to \(z\), the nonlocal part of the chiral vacuum covariance is proportional to \(1/\sinh[(z-z')/2]\). Its Fourier multiplier, including the diagonal half-identity, is \(1/(1+\exp(\pm2\pi p_z))\), with the sign set by chirality. The modular generator is consequently \(\pm2\pi p_z\). In station coordinates this gives the circle-weighted chiral energy \( \pm2\pi[-i(f\partial_x+f'/2)]\). The native dispersion has the two corresponding linear branches at \(k=\pm\pi/2\). A uniform estimate for the finite modular logarithm requires control near occupancies zero and one. The four aperture measurements are a finite comparison, with a fixed mesh and ring circumference. A three-dimensional Einstein conclusion additionally requires physical proper-area and volume records and a common stress on the theorem's frame and radius family. [Casini and Huerta, massless Dirac modular flow](https://arxiv.org/abs/0903.5284) · [Eisler, Tonni and Peschel, continuum and finite-chain comparison](https://arxiv.org/abs/1902.04474)

From the repository root, run the following with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_einstein_station.critical_exchange.measure
python -m pytest -q model/experiments/interval_echo_einstein_station/critical_exchange
~~~

## How observers select and disturb the reference

Work comparisons select the half-filled reference from an arbitrary input with the same total load. The controller observes the hopping energies before it measures any interval entropy. Once the reference is cold, fixed guided propagation restores its local bond readings after a bounded disturbance, on a pre-return time window of the finite guide. The original all-open federation gives a different evolution. [Independent work selection](selection/README.md) · [Local return and control](relaxation/README.md)
