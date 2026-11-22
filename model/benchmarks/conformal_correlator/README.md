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

Run from the repository root.

```sh
python -m model.benchmarks.conformal_correlator.one_loop --output /tmp/echo-cft.json
python -m pytest -q model/benchmarks/conformal_correlator
```
