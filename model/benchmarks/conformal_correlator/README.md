# Local correlations in an interacting four-dimensional conformal theory

This evaluator computes the protected two-point functions and the first interaction correction to a local four-point function in Euclidean SU(2) $\mathcal N=4$ supersymmetric Yang–Mills theory. The inserted fields are dimension-two, gauge-invariant scalar bilinears in the stress-tensor multiplet. The gauge theory, its vacuum, spacetime coordinates and coupling are inputs. The calculation gives a numerical target for an observer instrument. [Four-point function, equations (2.11), (2.18), (2.22) and (3.13)](https://arxiv.org/abs/1108.3557) · [Independent box-integral formula](https://www.sciencedirect.com/science/article/pii/0370269393911187)

Choose two orthogonal pairs of complex null polarizations with $Y_1\cdot Y_2=Y_3\cdot Y_4=1$ and all cross-pair contractions zero. Write $d_{ij}=|x_i-x_j|^2$ and $a=g_{\rm YM}^2 N_c/(4\pi^2)$. The two surviving two-point functions are $C_{ij}=(N_c^2-1)/(2(4\pi^2)^2d_{ij}^2)$. Free connected Wick graphs vanish in this polarization. The published loop correction reduces to

$$
\frac{\langle O_1 O_2 O_3 O_4\rangle}{C_{12}C_{34}}
= 1 - \frac{2 a\,d_{12}d_{34}}{N_c^2-1}
  \int_{\mathbb R^4}\frac{d^4y}{\pi^2\prod_i |y-x_i|^2}
  + O(a^2).
$$

At $a=0.01$, the normalized values **through one loop** are **0.984374** for the regular tetrahedron and **0.977909** for the skew configuration. At zero coupling both are exactly one. The box integral was computed from all twelve ordered Feynman-parameter sectors, then checked against the independent cross-ratio dilogarithm. The relative difference was below $10^{-12}$ on both fixtures. Simultaneous permutations of positions and polarizations, scaling and inversion preserve the expected weights. These numerical errors concern the coefficient of $a$. Higher perturbative orders have their own contributions. [Coordinates, two-point functions and controls](one_loop.json) · [Evaluator](one_loop.py)

The flat-space conformal vacuum has excitations at arbitrarily low energy. An exact dilation scales each nonzero energy by an arbitrary positive factor. A positive physical mass gap therefore requires a different regime or sector. The pure SU(2) lattice calculation uses its own action and displays gap bounds at finite volume. [Pure-gauge Hamiltonian and bounds](../../../code/yang_mills/README.md#finite-volume-gap-for-pure-su2-in-three-spatial-dimensions)

## Electron mass in the conformal vacuum and on the Coulomb branch

Suppose a nonempty charge sector has an internal charge that commutes with dilations. Let $m_{\rm th}$ be the infimum of its positive-energy invariant masses. Every dilation by $\lambda>0$ maps this sector to itself and multiplies each invariant mass by $1/\lambda$. Hence $m_{\rm th}=m_{\rm th}/\lambda$ for every $\lambda$, forcing $m_{\rm th}=0$. This argument permits a continuum of charged states and requires no isolated particle pole. It excludes the observed electron as a lowest massive charged excitation of this flat, dilation-invariant vacuum. Its measured rest energy is $0.51099895069(16)$ MeV. The two scalar four-point ratios are unchanged when every separation is multiplied by 1/2, 2 or 4. This numerical check probes the scalar insertions. The charged-sector conclusion follows from the symmetry. [CODATA mass calibration, Table XXXV](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf) · [Scalar evaluator](one_loop.py)

On the Coulomb branch of the same gauge theory, choose a scalar expectation $vT^3$, where $T^3=\operatorname{diag}(1,-1)/2$. The raising and lowering matrices obey $[T^3,T^\pm]=\pm T^\pm$. The quadratic adjoint mass operator has eigenvalues $g^2v^2,0,g^2v^2$ on $T^-,T^3,T^+$. For each charged root, unbroken supersymmetry gives a massive multiplet with spin-1, spin-1/2 and spin-0 states of the same mass. Setting $|gv|$ from the electron's laboratory mass calibrates every member to that mass. This gives a charged fermion in a specified vacuum, with charged scalar and vector partners and an adjustable scale. It supplies no electron identification in Interval Echo. [Coulomb-branch multiplets](https://arxiv.org/abs/1902.07205) · [Spectrum, Eq. (6.46)](https://doi.org/10.1007/JHEP09(2018)125)

Run from the repository root.

```sh
python -m model.benchmarks.conformal_correlator.one_loop --output /tmp/echo-cft.json
python -m pytest -q model/benchmarks/conformal_correlator
```
