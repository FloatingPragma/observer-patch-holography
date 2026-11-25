# Does the same orbital correction predict shear around different galaxy types?

Brouwer and colleagues measured different weak-lensing signals around early- and late-type galaxies after matching their stellar-mass distributions. Their Figure 8 release includes the covariance between the two populations. Here we ask how large that difference is in units of the deep orbital correction calibrated on SPARC, whether unknown sampling within the bins can account for it, and what a radial-pressure explanation would require at a finite outer boundary. [KiDS measurement, especially §5.4](https://doi.org/10.1051/0004-6361/202040108) · [Released shear and covariance](https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar) · [Frozen orbital calibration](../receipts/lensing_profile.json)

## What the shear measures

At baryonic acceleration $b=GM_b/R^2$, a compact baryonic lens plus a spherical scale-free envelope has excess surface density $Y=A(b)+qZ(b)$, where $A=4b/(\pi\,4G)$ and $Z=\sqrt{a_0 b}/(4G)$. The point-mass identity is $\Delta\Sigma=M_b/(\pi R^2)$. The unit envelope is the weak-field, zero-slip optical prediction for the orbital excess $\sqrt{a_0 b}$, with the full-sample SPARC value $a_0=8.7504\times10^{-11}\,{\rm m\,s^{-2}}$. The ratio $q$ is dimensionless. Measure shear directly in the released excess-surface-density units. Inferring an orbital acceleration from shear would already assume the optical transfer under test. [Unsplit radial-profile experiment](../lensing_profile.py)

We use the six bins lying wholly above the analysis isolation floor of $10^{-13}\,{\rm m\,s^{-2}}$ and wholly below the deep cutoff of $1.2\times10^{-11}\,{\rm m\,s^{-2}}$. The two releases split on rest-frame $u-r=2.5$ and Sérsic index $n=2$. In each split the high group denotes redder or more centrally concentrated lenses. We restore the full 30 by 30 matrix, apply the author's multiplicative corrections to shear and covariance, verify both marginal standard errors and then select the twelve required measurements. This preserves cross-population covariance. The two splits reuse the survey and their significances cannot be combined as independent evidence. [Measured bins and source hashes](profiles.json) · [KiDS science products](https://kids.strw.leidenuniv.nl/sciencedata.php)

Let $d=DY$ with $D=[-I,I]$ acting on the two stacked profiles, $V=DCD^T$, and $w=V^{-1}Z/(Z^TV^{-1}Z)$. At the published bin centres any shared radial signal cancels in $d$. The measured contrast $w^Td$ is **0.586 ± 0.121 by colour** and **0.467 ± 0.122 by Sérsic index**. Applying the same weights to cross shear gives −0.016 and −0.068, using the tangential covariance for scale. The squared covariance distances from equal profiles are 27.15 and 17.73 in six dimensions. Fitting one difference in envelope amplitude leaves distances of 3.54 and 3.00 in five residual dimensions. These are retrospective conditional statistics. The independently calibrated orbital scale supplies their units. Only a microscopic source that predicts the amplitudes without using these shear measurements would supply an OPH prediction. [Covariance calculation and results](results.json)

### What if the red and blue stacks sample different parts of a bin?

For any common $q\ge0$ and any positive normalized pair weights, the predicted mean is between $A(l)+qZ(l)$ and $A(u)+qZ(u)$ in a bin $[l,u]$. Let $k=\sqrt{u/l}$, which is constant across these logarithmic bins. Even if the high group samples the upper edge and the low group samples the lower edge,

$$
Y_{\rm high}-kY_{\rm low}\;\le\;A(u)-kA(l).
$$

This cancels the unknown common $q$. With the minimum-variance unit-envelope weights for this contrast, all six weights are positive, so the inequality survives averaging. The measurements exceed that upper bound by **1.51 statistical errors** in colour and **0.50** in Sérsic index. Matching the stellar-mass distributions does not establish equality of the source–lens pair distributions inside acceleration bins. Publishing the weighted bin means of $b$ and $\sqrt b$ by type would directly fix the prediction for this law.

The whole twelve-bin calculation allows every stack and bin its own positive pair distribution. Every value between the two edges is attainable as an edge mixture. At frozen $q=1$ the minimum squared covariance distances are 16.81 and 16.24. Allowing one common $q\ge0$ gives 11.30 and 9.03, at fitted values 1.166 and 1.236. The code brackets each global minimum between a feasible squared distance and a lower bound from a supporting plane. In detail, if every allowed mean $m$ satisfies $v^Tm\le H$, completing the square gives $(Y-m)^TC^{-1}(Y-m)\ge2(v^TY-H)-v^TCv$. Each interval is narrower than $10^{-6}$. These distances to a constrained set have a boundary-dependent sampling distribution. Baryonic calibration, extended gas, projection, companions and population-dependent environmental structure remain physical sources of uncertainty. The KiDS authors discuss gas and environmental explanations for the original type dependence. [Fitted and frozen distances](results.json) · [Brouwer et al.](https://arxiv.org/abs/2106.11677)

## What radial pressure would require

Assume ordinary weak-field geodesic optics in a static spherical metric, in units $G=c=1$, with gravitational potential $\Phi$ and mass function $m$. The leading orbital and Weyl mass profiles obey $m_{\rm RC}=r^2\Phi'$ and $m_{\rm lens}=(m_{\rm RC}+m)/2$. On a scale-free exterior $m=hr$, $\rho=h/(4\pi r^2)$, and $p_r=w_r\rho$, the radial Einstein equation gives $s=m_{\rm RC}/r=h(1+w_r)$ to leading order. Hold the independently measured orbital amplitude $s>0$ fixed. Where the scale-free optical projection applies, the lens amplitude relative to the zero-pressure orbit prediction then obeys

$$
q=\frac{2+w_r}{2(1+w_r)},\qquad
h=s(2q-1),\qquad
w_r=\frac{2(1-q)}{2q-1}.
$$

Thus $q>1$ requires radial tension for this explanation. With the nominal compact baryon term, the red stack at bin centres gives $q=1.413\pm0.080$ and the blue stack $q=0.827\pm0.090$ in separate marginal fits. Those conditional amplitudes compare two surveys' populations. KiDS has no matched independent rotation measurement in this calculation. Finite line-of-sight support, gas and departures from spherical symmetry can also change the shear. A stress inferred by solving this equation for $w_r$ has already used the light deflection it is meant to predict. [Rotation and lensing mass derivation by Faber and Visser](https://arxiv.org/abs/astro-ph/0512213) · [Finite-edge projection](../finite_halo/README.md)

### A lower mass bound for returning radial tension to vacuum

Let a static, spherically symmetric effective stress have $p_r(r_0)<0$ and $p_r(R_t)=0$ in a smooth transition to vacuum. Assume $R_t>r_0>0$, nonnegative density, $\rho+p_r\ge0$, $p_t\le\rho$ and attractive potential $\Phi'\ge0$ throughout. The mass function has $m'=4\pi r^2\rho$. The radial conservation equation is exact, with no weak-field expansion,

$$
(r^2p_r)'=2rp_t-r^2(\rho+p_r)\Phi'.
$$

Integration and the stated sign conditions give

$$
r_0^2|p_r(r_0)|\le 2\int_{r_0}^{R_t}rp_t\,dr
\le 2\int_{r_0}^{R_t}r\rho\,dr
\le \frac{m(R_t)-m(r_0)}{2\pi r_0}.
$$

Consequently the transition must contain at least $2\pi r_0^3|p_r(r_0)|$ of mass. If the inner excess has $m_A(r_0)=hr_0$ and $p_r(r_0)=w_r h/(4\pi r_0^2)$, this costs at least $|w_r|/2$ of that inner excess mass. If the density throughout also obeys $\rho(r)\le h/(4\pi r^2)$, the second integral is bounded by $(h/2\pi)\log(R_t/r_0)$, so $R_t/r_0\ge\exp(|w_r|/2)$. For illustration, assigning the entire nominal red-stack enhancement to the scale-free pressure gives $w_r\simeq-0.452$, a minimum transition mass of 22.6% of the inner excess, and a minimum radius ratio of 1.25 under the additional density cap. These are necessary costs conditional on that assignment, including scale-free optical projection. A layer with surface stress, a violation of the energy conditions, exchange with another stress component or a time-dependent source requires its own matching calculation. [Finite spherical halo and conservation convention](../finite_halo/README.md)

The OPH test would need an independently read collar mass, radial and tangential stress, outer boundary, and corresponding orbital and lensing profiles from one physical source. Differences in distance to Earth are not measurements of local OPH record density. The type contrasts and transition bound specify quantities to check if a native source supplies those records.

## Reproduce

From this directory, run the offline comparison and regression tests. Supply the original archive to verify source hashes and all extracted bins against the release.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python measure.py
python -m unittest discover -s . -p 'test_*.py'
curl -fsSL https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar -o /tmp/brouwer2021_rar.tar
python measure.py --archive /tmp/brouwer2021_rar.tar
```
