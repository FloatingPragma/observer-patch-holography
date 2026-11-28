# A stellar check on M31’s outer orbits

Keck spectroscopy of red-giant stars in M31’s dynamically hot halo gives a line-of-sight dispersion of **98.0 km/s**, with a published 16th–84th percentile interval of **91.4–105.2 km/s** in the overlapping 24–40 kpc bin. At the measured median radius our frozen finite-survival model predicts **105.5 km/s** for isotropic orbits. The lossless model gives **121.6 km/s** and baryons alone give **58.9 km/s**. Radially biased stellar orbits give **98.4, 117.5 and 53.1 km/s**, respectively, at the preselected anisotropy $\beta=0.5$. This comparison supports the finite model over the lossless continuation on the stated dynamical assumptions. [Stellar spectroscopy, Table 3 and Figure 7](https://doi.org/10.3847/1538-4357/aa9f26) · [Frozen source and previous gas comparison](README.md) · [All sensitivity cases](tracer.json)

## How the stellar measurement enters

We selected SPLASH on its tracer identity and radial coverage before opening its dispersion table. Only the 24–40 kpc bin overlaps the disputed H I rings. It contains 382 spectra left after satellite rejection, including 240 stars whose membership probabilities favor M31. The authors fit a velocity mixture containing Milky Way contaminants and M31 tidal debris. The reported dispersion is the width of the hot component about its fitted mean. Its posterior interval describes that mixture fit. The authors estimate that measurement and frame-transformation errors typically inflate the measured dispersion by only a few tenths of a km/s. [Selection record](TRACER_PROTOCOL.md) · [Stellar spectroscopy, §§II–V](https://arxiv.org/html/1711.02700)

The authors plot the median of all stars in the bin, with the actual sample range as a horizontal bar. The archived vector figure places that bar at **27.06–38.18 kpc** with median **33.20 kpc** for their 783 kpc M31 distance. We rescale the angular positions to the 785 kpc normalization of the frozen galaxy. This gives 27.13–38.27 kpc and a median of 33.29 kpc. The measured speed dispersion is unchanged by this coordinate conversion. We also bound the forecast across the entire populated range, since the aggregate hot-star estimate does not specify the effective radial weights of its likelihood. [Figure 7 and distance convention](https://arxiv.org/html/1711.02700) · [Primary table and vector-coordinate checks](tracer.py)

## Forward projection

For the hot stellar tracer we adopt $\nu(r)\propto r^{-\gamma}$ and $\gamma=3.2$. Under spherical deprojection the density exponent is one steeper than the projected brightness exponent. Gilbert et al. measure a projected slope of $-2.2\pm0.3$ after statistically removing cold substructure, over 9–90 kpc. This source of tracer weights shares some spectra and debris classification with the velocity analysis. Its slope uses star counts and was fixed separately from the present comparison. We test $\gamma=2.9,3.2,3.5$. Extending the measured outer slope is an assumption. Using the paper’s cored profile with the same asymptotic slope and core radius 5.2 kpc changes the isotropic finite forecast to 105.95 km/s. Their photometry allows spherical symmetry and also admits a mildly prolate halo. [Stellar surface brightness, §IV and Table 3](https://arxiv.org/html/1210.3362#S4.SS1)

The stars are pressure-supported tracers of the frozen total field. We use the prior spherical baryonic force, the same registered dark source and the same finite or infinite survival. For a stationary spherical tracer with zero mean streaming, constant $\beta=1-\sigma_\theta^2/\sigma_r^2$ and radial pressure vanishing at infinity, the Jeans solution and its line-of-sight projection are
$$
\sigma_r^2(r)=\int_0^1 t^{\gamma-2\beta-1}\,v_c^2(r/t)\,dt,
\qquad
\sigma_{\rm los}^2(R)=
\frac{\int_0^{\pi/2}\cos^{\gamma-2}\theta\,(1-\beta\cos^2\theta)\,
\sigma_r^2(R/\cos\theta)\,d\theta}
{\int_0^{\pi/2}\cos^{\gamma-2}\theta\,d\theta}.
$$
The gravitational field from all source shells is integrated to infinity, including the analytic outer tail. For the central finite model **84.4%** of the predicted variance at the median comes from the force at true radii beyond the outermost measured H I ring at 38.09 kpc. The contributions above 90 and 175 kpc are 10.6% and 1.56%. These are variance weights in the Jeans equation. They differ from the fraction of stars located outside those radii. [Numerical projection and equations](tracer.py) · [Independent single-integral and analytic checks](test_tracer.py)

## What the stars can distinguish

| Halo or gas input | Isotropic prediction at 33.29 kpc | Constraint across the populated bin |
| --- | ---: | ---: |
| Baryons alone | 58.9 | 55.0–65.2 |
| Frozen finite mark lifetime | 105.5 | 101.9–110.8 |
| Infinite mark lifetime | 121.6 | 119.7–124.4 |
| Chemin H I, one favorable quoted error per speed | at least 123.0 | at least 116.7 |
| Corbelli H I, one favorable quoted error per speed | at least 98.9 | at least 92.2 |

Speeds here are predicted **stellar dispersions** in km/s. The H I comparisons use only a lower bound on the mass required under circular disk motion and a nonnegative spherical dark density. For each observed H I ring define $G M_i=r_i[(v_i-k\,\delta v_i)^2-v_{b,\rm disk}^2(r_i)]$. At true radius $s$ impose $GM_A(s)\ge\max(0,\max_{r_i\le s}GM_i)$. This gives a positive stepwise enclosed dark mass. Use the spherical ordinary mass for the halo stars and project this minimum dark field. Any further positive dark mass raises the stellar dispersion under the same Jeans assumptions. This construction allows the most rapid permitted decline beyond the gas and keeps the galaxy’s stellar input fixed. It was added after opening the stellar data. Simultaneous one-error shifts are sensitivity cases and have no stated joint probability. The corresponding primary comparison under a *constant* continuation of the outer gas circular speed throughout all sampled forces would give much larger dispersions of 153.8 and 134.2 km/s. Measuring a last gas speed does not specify that continuation. [H I data and ordinary-force subtraction](results.json) · [Positive mass bounds](tracer.json)

The table covers unknown radial weights. For a power-law tracer at constant anisotropy, inserting $v_c^2(s)=GM(s)/s$ into both positive projection integrals shows $R\sigma_{\rm los}^2(R)$ is nondecreasing when $M$ is nondecreasing. Enclosing that product by its values at adjacent grid edges bounds $\sigma_{\rm los}(R)$ continuously within each interval. Any positively weighted average of the variances stays inside those bounds when the component mean is common across the bin. The published mixture likelihood may have additional selection dependence.

On isotropy and the central tracer slope, even the minimum mass required by Chemin’s outer disk gives a stellar dispersion above the published 84th percentile across the entire bin. Corbelli’s broader errors permit a lower envelope inside that interval. The result changes with orbital anisotropy and tracer density. With $\gamma=3.5$, $\beta=+0.5$ and the same favorable gas shifts, the lower bound across the bin falls to **101.3 km/s** for the Chemin profile. Thus the non-H I data favor a lower sustained halo force under isotropy and constrain the projected pressure. One dispersed stellar bin cannot resolve the local circular-speed rise without stronger knowledge of stellar anisotropy and angular sampling. [Stellar error model](https://arxiv.org/html/1711.02700#S5) · [All assumptions and sensitivity combinations](tracer.json)

Run from the repository root with one BLAS thread and `PYTHONPATH=code:.`. The fetch option also verifies the primary spectroscopy table, the source PDF plot coordinates and the hot-tracer photometric slope. The default uses the archived measurements.

~~~sh
python -Bm cosmology.rar_deep_regime.m31_transport.tracer --fetch
python -Bm pytest -q -W error code/cosmology/rar_deep_regime/m31_transport
~~~
