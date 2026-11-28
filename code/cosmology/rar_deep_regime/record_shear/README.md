# A lensing prediction from a specific record-formation process

If the stationary mark-walk response supplies the extra mass density, the measured excess surface density must be **concave as a function of inverse squared radius** outside the projected mark injectors. This prediction uses no fitted acceleration, mark rate or galaxy mass scale. We compared all 20 three-bin inequalities in each of two published isolated-galaxy shear profiles. The largest departures in the forbidden direction were **0.40 standard errors for KiDS-bright and 1.46 for GAMA**. These measurements permit the predicted shape. A common inverse-square density and other ordinary mass distributions can also pass it. [Local test design](PROTOCOL.md) · [Published shear profiles](https://kids.strw.leidenuniv.nl/sciencedata.php) · [Brouwer et al.](https://doi.org/10.1051/0004-6361/202040108)

## From record formation to a sign on the sky

On the published comparison chart a mark stays at its address with probability 1/4 and takes each of six adjacent steps with probability 1/8. With survival $s$ and compact injection $J$, the stationary visits obey $u=J+sPu$. On the native triangular source cells, the sparse, paired, sign-resolved record reader has leading response proportional to the sum of squared differences of the visit amplitudes attached to the two owners of each shared sign. Other registered examples include a squared boundary average when the majority sign is retained. These follow from the chosen correlation-weakening intervention with its positive reference and the sparse-admission limit. [Stationary transport](../../../matched_einstein_limit/collar_identifiability/refinement/transport.py) · [Boundary-sign reader](../../../matched_einstein_limit/collar_identifiability/refinement/boundary_sources.py) · [Scope of the native records and chart](../../../matched_einstein_limit/collar_identifiability/refinement/README.md)

Write the density candidate as $q=\sum_e v_e^2$ with each $v_e$ a fixed linear combination of same-transport fields. This also covers any fixed positive semidefinite quadratic form. On the cubic chart, with unit grid spacing and $j_e$ the same combination of injectors,
$$
\Delta_3 v_e=8(s^{-1}-1)v_e-8j_e/s,\qquad
\Delta_3 v_e^2=2v_e\Delta_3v_e+\sum_{y\sim x}(v_e(y)-v_e(x))^2.
$$
It follows that $\Delta_3 q\ge0$ off the sources. For the projected column $\Sigma=\sum_z q$, the sum of the interior vertical second differences leaves the negative of the two boundary-face values. Moving them to the other side adds nonnegative terms to $\Delta_2\Sigma$. This works for arbitrary source shape, even when some of the signed link fields change sign. At an injection the additional term $-16v_ej_e/s$ can reverse the sign.

A 31³ absorbing chart used the incidence matrix of the four cells descended from one face of the 80-cell federation. All **957** sky columns without an injector passed for survival 1 and 7/8. The numerical balance residual was below $3\times10^{-14}$. Controls included a column with injection (negative curvature), reversal of the modular sign (negative curvature away from injection), an independent continuum calculation for a displaced point source, and a deliberately convex shear law. [Calculation](measure.py) · [Independent checks](test_measure.py) · [Finite results](results.json)

In a continuum transport chart where $\nabla^2 v_e=\mu^2v_e$ away from injections, $\nabla^2 q=2\sum_e(|\nabla v_e|^2+\mu^2 v_e^2)\ge0$. Write $\Sigma(R,\phi)=\int q\,dz$. For a complete line of sight with vanishing end gradients, $\nabla^2_\perp\Sigma\ge0$ wherever that sightline misses every injector. This step uses a transport chart; the projected density interpretation also requires a common physical ruler and ordinary thin-lens coupling with negligible stress corrections. A time-dependent, advected or spatially heterogeneous transport law can change the sign.

Let $x=\ln R$, $S=\Delta\Sigma=\bar\Sigma(<R)-\langle\Sigma(R)\rangle_\phi$ and $t=R^{-2}$. Angular averaging the two-dimensional Laplacian removes the angular derivative by periodicity, giving $\langle\Sigma\rangle_{xx}\ge0$ on a complete source-free annulus. Differentiating the enclosed mean gives $S_x=-2S-\langle\Sigma\rangle_x$. Hence
$$
4t^2 S_{tt}= S_{xx}+2S_x=-\langle\Sigma\rangle_{xx}\le0.
$$
A compact ordinary baryonic mass adds a term linear in $t$, and a uniform mass sheet contributes zero to $S$. Both cancel in every three-point curvature test. A finite return budget, a sharp imposed truncation, injection throughout the annulus, or a spatially varying conversion of visits into stress needs its own source calculation.

## Comparison with measured tangential shear

The measured horizontal coordinate in the release is the estimated point-source baryonic acceleration $G M_b/R^2$ for each lens. Multiplying inverse squared radius by a positive constant preserves concavity for each object. A fixed positive mixture of such lens profiles preserves it as well. This application assumes that the selected bins sample a common mixture of isolated lenses, that the ordinary baryons are compact inside each tested annulus, and that the annuli are free of additional mark injectors in projection. The public binned profiles do not locate injectors or fix the lens-by-lens weights independently in each bin. [Lens selection and binning](https://doi.org/10.1051/0004-6361/202040108)

For every ordered triple in the six enrolled bins, let $t_1<t_2<t_3$ and $\lambda=(t_3-t_2)/(t_3-t_1)$. Concavity requires $C=\lambda S_1+(1-\lambda)S_3-S_2\le0$. The error on $C$ uses the full released covariance after multiplicative shear calibration. A common positive weighting of fractional positions within all logarithmic bins gives another concave function of their centres. For arbitrary positive weights of a continuous common profile, each bin has some effective position between its edges. The extremal chord weights then bracket all possibilities. A conservative violation requires both edge-weight contrasts to exceed the threshold. No weights are adjusted to make a triple pass or fail.

| Stack | Largest centre chord excess | Standard error | Centre excess / error | Largest violation robust to unknown bin weights |
| --- | ---: | ---: | ---: | ---: |
| KiDS-bright | 0.875 | 2.179 | 0.402 | −1.003 standard errors |
| GAMA | 13.836 | 9.500 | 1.456 | 1.214 standard errors |

Excess and error are in the release’s $h_{70}\,M_\odot\,{\rm pc}^{-2}$. The strongest centre contrasts both concern the bins at $4.04,7.13,12.58\times10^{-13}\ {\rm m\,s}^{-2}$. The one-sided 95% family threshold is 2.807 standard errors across 20 triples under the reported Gaussian covariance. The cross-shear controls gave maxima of 0.318 and 0.715 standard errors using the same covariance. The two samples share imaging and are evaluated separately. Overall shear measurements on these galaxies were known before the shape test was specified. [Published data release](https://kids.strw.leidenuniv.nl/sciencedata.php) · [All measurements and covariance matrices](results.json)

A resolved violation on common lenses with independently measured source-free annuli would reject the stated record-to-lensing assignment at any positive normalization. Passing the shape test leaves the charge per cut, its gravitational coupling and the total dark abundance undetermined.

Run from the repository root with one BLAS thread and `PYTHONPATH=code:.`.

~~~sh
python -Bm cosmology.rar_deep_regime.record_shear.measure
python -Bm cosmology.rar_deep_regime.record_shear.measure --archive /path/to/brouwer2021_rar.tar
python -Bm pytest -q code/cosmology/rar_deep_regime/record_shear
~~~

The archive mode verifies the authors’ file hash and reconstructs both profiles and covariance matrices. The default reuses the archived six-bin inputs and recomputes the graph and observational tests.
