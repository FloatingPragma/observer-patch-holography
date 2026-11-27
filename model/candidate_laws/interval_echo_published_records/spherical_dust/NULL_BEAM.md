# Light-bundle focusing at the ordinary coupling

The homogeneous interior of the falling dust sphere obeys physical null balance on the same archive preparation. An observer measures transverse angular separation of outgoing light, normalizes its affine parameter using phase and dust-clock observations, and compares the resulting focusing with a separate rest-action and frequency meter. The held source uses the coupling previously measured by releasing an ordinary exterior probe. At affine distance \(1/8\), the held optical and matter readings are 1.247056797 and 1.246635999. At fixed affine length \(1/4\), the integrated optical deficit and matter prediction differ by \(1.17\times10^{-9}\). [Finite sensor tapes](null_beam.json) · [Acquisition and independent decoders](null_beam.py)

## The pointwise and integrated statements

Use the interior scale \(a(\eta)=(1-h\eta/2)^2\), \(h^2=2GM_g\), normalized to one at the outgoing vertex. \(M_g\) is the stock to which the shell motion responds. For the common-source law it is the total pressureless rest-action stock \(M_m=M_o+M_c\) initially within unit radius. A radial conformal ray has coordinate distance \(\eta\). Emit it with unit frequency in the vertex dust clock. At another dust clock the observed frequency is \(\omega=1/a\), and \(d\lambda=a^2d\eta\) defines the affinely normalized parameter. The proper transverse separation per unit angular chord is \(D=a\eta\). Local spatial isotropy gives screen matrix \(D I_2\) and zero shear. Angular Jacobi propagation then measures \(R_{kk}=-2D''/D\) at positive affine distance. At the vertex \(D(0)=0,D'(0)=1\).

Put \(q=1-h\eta/2=(1-5h\lambda/2)^{1/5}\). Directly differentiating the physical separations with \(d/d\lambda=-(h/2)q^{-4}d/dq\) gives
\[
 D=\frac2h(q^2-q^3),\qquad D'=3q^{-2}-2q^{-3},\qquad
 -2D''/D=3h^2q^{-10}.
\]
For zero stock the continuous limits give \(D=\lambda\) and zero focusing. Independently, the comoving worldlines have \(\rho=3M_m/(4\pi a^3)\), zero pressure and \(T_{kk}=\rho\omega^2\). Thus the forecast \(8\pi G T_{kk}=6GM_m q^{-10}\) equals the screen reading for \(M_g=M_m\). The same conclusion holds locally in every direction and at every regular interior event by the spatial symmetries and affine rescaling. The vertex choice fixes units for the finite comparison. [Shell dynamics and worldline stress](PROOF.md)

Twice integrating the angular Jacobi equation gives an additional finite-length comparison:
\[
 2[L-D(L)]=\int_0^L (L-\lambda)D(\lambda)\,8\pi G T_{kk}(\lambda)\,d\lambda .
\]
The left uses the endpoint screen alone; the right uses intermediate weights and an independently normalized matter flux. The endpoint carries zero weight on the right. The cross-sectional area is \(4\pi D(L)^2\) for unit-sphere solid angle and its signed deficit relative to \(4\pi L^2\) has leading coefficient \((2\pi/3)R_{kk}(0)L^4\) in \(4\pi[L^2-D(L)^2]\). This is fixed-affine cone transport. The null premise of the tensor theorem is the pointwise \(R_{kk}=8\pi G T_{kk}\); a ball-entropy modular generator and its KMS normalization are different observables. The interior proof excludes the material surface discontinuity, any exterior shear and the claim that native Interval Echo produces this spherical constitutive law. [Einstein theorem from common null balance](../../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

## What each detector contributes

At subdivisions \(d=32,64,128\), acquire targets \(\lambda_i=i/d\) through \(L=1/4\) and use tick \(\epsilon=d^{-10}\). The emitter sets the phase to the nearest tick. The affine reader receives four raw dust clock marks at conformal phase \(-1/4,0,1/4,1/2\). They determine the cubic proper-clock polynomial \(\tau(\eta)\) and the rate \(a=d\tau/d\eta\). It reconstructs \(\lambda=\int_0^\eta a^2\,dz/a(0)\) without the mass or coupling. Four further independent dust timestamps at phase offsets \((-2,-1,1,2)/d^2\) read the local rate for the photon-frequency meter. Their five-point divided difference is exact on a cubic before quantization. This phase/clock convention can be checked at every acquired event. The held finest tape's maximum affine reconstruction error is \(3.96\times10^{-22}\).

Three transverse chord meters span the \(x,y,45^\circ\) directions, with chord factors \(14/25,10/13,16/17\). Each records an integer proper separation. Dividing by its own chord gives \(D\); the diagonal channel tests the off-diagonal screen distortion. At the finest held midpoint the integer readings are 82530180366150873941, 113365632371086365303 and 138706185489329199900, in units \(128^{-10}\). The largest discrepancy in normalized angular distance on that tape is \(1.14\times10^{-21}\). Isotropy proves the all-direction claim; the three readings test that stated screen model.

The matter receiver weighs the conserved shell tickets between labels \(1/4\) and \(1/2\), and independently measures their two proper radii \(r_i,r_o\) at the appointment. It reads the density from \(\Delta M/[4\pi(r_o^3-r_i^3)/3]\). For the finest held midpoint the two radius counts are 274257154389995663771 and 548514308779991327543, the shell mass is 0.062973485, and \(\omega=1.076172127\). The ordinary release at radius two gives \(G_{\rm obs}=0.25000002649\). The work forecast is \(6G_{\rm obs}\Delta M\omega^2/(r_o^3-r_i^3)\). Neither the source masses nor any fitted screen coefficient enter the optical decoder. [Ordinary release](body.py) · [Mass-ticket protocol](survey.py)

| Stock, finest grid | Optical \(R_{kk}\) | \(8\pi G_{\rm obs}\rho\omega^2\) | \(2[L-D(L)]\) | Weighted work integral |
| --- | ---: | ---: | ---: | ---: |
| Ordinary | 0.65707713 | 0.65692596 | 0.00173123267 | 0.00173123248 |
| Ordinary plus held coincidence | 1.24705680 | 1.24663600 | 0.00332121029 | 0.00332120912 |
| Half geometric response to coincidence | 0.94169716 | 1.19247678 | 0.00249398275 | 0.00315904412 |
| Held geometry with coincidence omitted from work | 1.24705680 | 0.72173663 | 0.00332121029 | 0.00192280528 |

The held-minus-ordinary changes are 0.58997967 and 0.58971004 in the point comparison, and 0.00158997762 and 0.00158997664 in the integrated comparison. These are independent physical readouts inside the declared candidate. They do not identify the coincidence stock with astronomical dark matter.

## Finite error estimates

The point reader differences three measured screen distances about \(\lambda=1/8\). All stocks satisfy \(h\le11/20\), \(M_m\le3/5\), \(q\ge0.9\) along the ray, and \(q\ge0.95\) throughout the point stencil. Each screen distance differs from its nominal-affine value by at most \(3\epsilon\): phase rounding shifts \(\lambda\) by at most \(\epsilon/2\), \(|D'|\le2\), and every chord and length tick is known. For spacing \(s=1/d\) the second-difference truncation is at most \(s^2\sup|D^{(4)}|/12\), with \(D^{(4)}=-h^3(624-504q)/(8q^{18})\). Propagating the three screen tick errors and the denominator through \(-2D''/D\) gives the optical bound implemented in the decoder.

Independent errors enter the material forecast monotonically. Each radius has error at most \(\epsilon/2\); the four local clock marks give rate error \(3 d^2\epsilon/4\). The four calibration marks imply initial rate error below \(10\epsilon\), and uniform rate error below \(100\epsilon\) for \(0\le\eta\le1/3\) by the displayed cubic interpolation in the code. The affine schedule error test \(120\epsilon\) includes rate squaring, initial normalization and the rounded target phase. A material forecast moved from its rounded phase to the nominal affine event costs less than \(30\epsilon\) using the derivative of \(6GM_m q^{-10}\). The ordinary exterior release starts at radius two. Twice integrated acceleration is a positive weighted average of \(GM_o/r(t)^2\). Its two measured endpoint radii therefore enclose the entire falling arc, yielding a two-sided bound on the ordinary gain from those observations, including the radius tick. The held point residual decreases from 0.006745 to 0.001684 to 0.000421. Its finest deterministic bound is 0.000621.

Composite Simpson integration uses the positive weights at the uniformly spaced affine appointments and the interval bounds on each measured \(D\) and \(T_{kk}\). With \(J=8\pi G D T_{kk}=(12GM_m/h)(q^{-8}-q^{-7})\), the fourth derivative of the integrand is \((L-\lambda)J^{(4)}-4J^{(3)}\). Each derivative follows from \((d/d\lambda)^n q^p=(-h/2)^n q^{p-5n}\prod_{j=0}^{n-1}(p-5j)\). The Simpson remainder is bounded by \(Ls^4\sup|[(L-\lambda)J]^{(4)}|/180\). Adding endpoint, source, clock and ordinary-calibration uncertainties gives a finest held bound \(8.87\times10^{-8}\). Both altered mechanisms violate their own point and integral envelopes. The optical and matter records can also be corrupted independently in the regression tests. [Bounds and corruption tests](test_null_beam.py)
