# A common optical reading of the hidden two-head phase

Two adjacent launch ports can reveal a phase that every sharp-port leading-arrival clock misses. On the same certified stationary source, the vacancy and the excess recover identical leading interference intensities. Both see the signed, volume-preserving phase change. Using an ordinary source to calibrate each aperture's threshold length nevertheless predicts the wrong sign for that change. [Stationary source](../../full_algebra/README.md) · [Sharp-port comparison and phase direction](../README.md) · [Aperture readings and counts](results.json)

## A material aperture on one blank seam

Use the triples $(x,y;z)=(0,5;39)$ and $(39,40;0)$. In both cases $x,y$ share an unmarked native seam and both are three graph hops from $z$. On a fresh source of either charge, a blank half exchange prepares $(|x\rangle+i|y\rangle)/\sqrt2$ up to a global phase. Leave the two parked heads untouched. Release all native services for an exposure $t$, then read the occupation of $z$. On another fresh source, launch from $z$ and apply the same blank half exchange before a sharp $x$ reading. Its bra has the same complex coefficients as the forward ket, corresponding to the conjugate aperture state. Choose the two orientations with equal probabilities for each single charge. Every failure to reach the terminal or the selected aperture output counts as a miss.

For one charge let $U,V$ be the real integer third-hop blocks from the respective launch ports to $z$, formed with minus twice the native off-diagonal Hamiltonian. The quarter-phase amplitude is proportional to $A=U+iV$. The reversed amplitude uses $A^T$. The symmetrized effect for $f=\lim_{t\to0}P(t)/t^6$ is
$$ E=\frac{A^\dagger A+\overline A A^T}{4\cdot 2^6(3!)^2}. $$
The denominator includes both the normalized aperture preparation and the orientation lottery. Switching charge transposes $U$ and $V$ on both parity blocks, so it interchanges the two positive effects. This identity follows from integer multiplication on the actual 960-port board. It holds on arbitrary two-head preparations in the enrolled parity blocks. Replacing the coherent launch by a classical mixture, or averaging positive and negative quarter phases, removes the imaginary cross term. Switching off the squared contrast on the same board makes both leading phase signals vanish. [Native blocks and apparatus replay](measure.py)

In the odd block let $D_1=Y\otimes Z/2$, $D_2=Z\otimes Y/2$, with zero perturbation in the even block. The two optical responses are
$$
\frac{\partial(f_L,f_R)}{\partial(D_1,D_2)}
=\begin{pmatrix}-1/8&0\\0&-9/64\end{pmatrix}.
$$
The Euclidean response norm is at least $1/8$ on every unit tangent in that plane. The independently volume-matched tangent at the numerical center has coefficients $(-0.9996945,-0.0247174)$, so the leading signed intensity changes are $+0.1249618$ and $+0.0034759$ per unit tangent. The exact source enclosure preserves both signs.

## Ordinary calibration and held source changes

Define the dimensionless threshold length $\ell=(2304 f)^{-1/6}$. This is the ratio of times needed to attain the same vanishingly small click probability. The denominator reference is the blank unique three-hop link $(13,792)$, whose leading coefficient is $1/2304$. The source-volume adjustment always comes from the independent six-unit four-station survey. The reference $\chi$ fixes the full head charge $K_{\rm full}=-\log\chi$ and the collar charge $K_{\rm collar}=\log(\chi_L\otimes\chi_R)-\log\chi$. The following values are population slopes at the numerical center.

| Volume-matched source | Full charge, nats | Collar charge, nats | Threshold slope $(0,5;39)$ | Threshold slope $(39,40;0)$ |
| --- | ---: | ---: | ---: | ---: |
| Ordinary exchange calibration | +0.409453 | +0.000345446 | +0.0723494 | +3.26651 |
| Held unpairing | +0.231106 | +0.231106 | +0.0395963 | +0.988636 |
| Held two-head phase | +0.000324171 | +0.000324171 | −0.0195700 | −0.000543722 |

At the equilibrium the two threshold lengths are 0.327932 and 0.327877 blank-link units. Fit one coefficient per aperture using only the ordinary charge and threshold slopes. For full head charge the ratios of observed to predicted unpairing are 0.96964 and 0.53622. Independently propagating source, volume and effect uncertainties encloses these ratios within $[0.9542,0.9847]$ and $[0.5266,0.5454]$. The phase ratios are negative and lie within $[-346.70,-336.61]$ and $[-0.2135,-0.2070]$. For collar charge the ordinary response lies between $0.0001884$ and $0.0005025$ nats. Its held unpairing ratios lie within $[0.000439,0.001209]$ and $[0.000242,0.000670]$, and both phase ratios again have the wrong sign. The angle adjustment preserves both modular readings by unitary invariance. The verifier replays native stationary-radar effects, encloses the volume tangents and the local marginal rotation, and bounds both logarithms. [Interval calculation](certificate.py) · [Ratios at the exact enclosed equilibrium](certificate.json)

## A finite exposure and its counting cost

At $t=1/64$ both species give positive normalized phase responses on both apertures. The pooled-orientation sampled values are $0.124531$ (vacancy) and $0.124312$ (excess) on the left, and $0.003942$ and $0.003955$ on the right. The simultaneous 99% sampling half-width is below $0.000211$ for these four measurements. At $1/128$ all four lower bounds also exceed zero. These finite-exposure readings can differ between charges. The coefficient identity concerns the limit.

The simulation uses $2^{68}$ binary trials per setting, sampled in binomial batches, and records all click and miss counts. The 16 source settings use native product Pauli preparations and their unread-sign counterparts. For a given product projector $R$, twice the response on $(R/2,R/2)$ minus the response on $(R,0)$ gives the odd-block response on $R$. Taking signed sums over the four Pauli eigenvalue pairs yields $D_1$ and $D_2$. Linearity recovers the local signed optical slope on the selected equilibrium without amplifying a finite difference between nearly identical sources. The independently determined volume matching weights are applied afterward. For each binomial count $X$ the simultaneous bound on its expected count is $\sqrt{2X s}+2s$ with $s=\log(2\cdot128/0.01)$, propagated through the signed sums. The count bank is synthetic. It illustrates the statistical price of a sixth-order optical signal. [Count bank and exact exposure settings](results.json)

The threshold compares an extended material aperture with a terminal. It measures a directed material response to phase and pairing. The experiment supplies neither spatial geodesics between sharp points nor the common spacetime null curvature required for Einstein's equation. Its sign test excludes a scalar charge-to-threshold calibration on these sources. A tensor stress, a different geometric instrument or a physical refinement requires its own transfer calculation.

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.shared_radar.aperture.measure
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.shared_radar.aperture.certificate
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/shared_radar/aperture
```
