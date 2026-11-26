# Two telescopes test stellar orbital anisotropy in J0946+1006

Keck measured the central stellar dispersion of the foreground galaxy as $287\pm11$ km/s. The same galaxy's VLT/MUSE observation places the second moment in a separate bin near 1 arcsec below 247 km/s after outward rounding. For the [frozen Chabrier orbital field](../README.md), the central Keck aperture can reach at most **265.4 km/s** under isotropy, even allowing arbitrary attractive orbital gravity outside $S=2.07$ arcsec. The two readings remain incompatible for every radial tracer profile $0\le\beta(r)\le r^2/(r^2+a^2)$ with **$a\ge3.11$ arcsec** (11.1 kpc). At the outer ring the limiting envelope already allows $\beta=0.307$, and it allows asymptotically radial orbits. These are sufficient exclusion bounds. Smaller $a$ permits a wider set of models.

The inherited angular distance is 3.57692 kpc per arcsec (flat matter-plus-lambda distance with Hubble constant 70 km/s/Mpc and matter parameter 0.3 at lens redshift 0.222). This is also the distance used in the Newtonian part of the force. Turner et al. use 3.69 kpc per arcsec in their own cosmology. The tables and bin supports are compared in angular coordinates.

The comparison assumes the fixed sphericalized stellar light and its role as a common spectroscopic tracer, the individual stellar mass and SPARC acceleration ranges, equilibrium spherical raw second moments, and the aperture and seeing responses below. Flattening or an additional interior acceleration changes the prediction. Stellar orbital anisotropy and the radial/tangential components of the **gravitating source** are separate quantities. The comparison constrains the field used by the [two-ring pressure test](../../README.md); the force carrier is undetermined.

## Measurements and numerical results

Sonnenfeld et al. publish 11 separately tabulated Keck/LRIS rectangular slit measurements, each with its own mean velocity and dispersion. Their two observing nights in December 2006 had 0.8 arcsec seeing. Their aperture is 1.00 by 0.21 arcsec and the reported trace-centroid accuracy is about 0.02 arcsec. We allow a Euclidean displacement of 0.03 arcsec, assume a uniformly transmitting rectangle and circular Gaussian seeing, and integrate numerator and brightness over the rectangle separately. For each measured signed velocity and dispersion we use the conservative raw-second-moment endpoint $\sqrt{[|v|-\Delta v]_+^2+(\sigma-\Delta\sigma)^2}$. The central result uses only $\sigma=287\pm11$; it remains valid with arbitrary systemic velocity because $v_{\rm rms}\ge\sigma$. All 11 rows are archived with the independent spectra. [Sonnenfeld et al., Tables 6 and 7](https://arxiv.org/html/1111.4215) · [Measured values](keck.json)

The independent outer second moment comes from the black marker nearest one arcsecond in Turner et al. Figure 12, measured with VLT/MUSE and 0.5 arcsec FWHM. Its whole adopted footprint is $[0.64,1.37]$ arcsec. The MUSE central bin is reserved to check the strong-radial-orbit example below. Its measured upper bar is 280.93 km/s, increased by two original raster pixels and rounded up to 282. Different point-spread functions, angular radii and observation dates enter the calculation. Neither dataset supplies a joint confidence level for the adopted endpoints. Sonnenfeld et al.'s Hubble photometry, stellar-population masses and light scales enter through the existing calibration; the observed Keck slit velocities never set these SPS masses. [Turner et al., Figure 12](https://arxiv.org/html/2401.08771v1) · [MUSE extraction and source provenance](../README.md)

| Orbit assumption | Keck central upper or exclusion | With matched MUSE luminosity centroid |
|---|---:|---:|
| Isotropy | 265.4 km/s | 256.7 km/s |
| Arbitrary radial changes of either sign, $|\beta|\le b$ | all $b\le0.018$ excluded | all $b\le0.036$ excluded |
| Any $0\le\beta(r)\le r^2/(r^2+a^2)$ | all $a\ge3.11$ arcsec excluded | all $a\ge2.51$ arcsec excluded |

The first result column accepts any nonnegative MUSE bin weights on the previously adopted support. The optional last column adds a measurement condition. The quoted luminosity-weighted MUSE bin position must use the same effective light weights as its extracted second moment. Its centre has radius 1.010 arcsec. Allowing 0.010 arcsec in raster position and 0.142 arcsec for arbitrary within-spaxel light displacement implies $\mathbb E R\ge|\mathbb E\mathbf R|\ge0.85$ arcsec. If those weights differ, use the support-only column. The general sign-changing bound concerns **spherical** $\beta=1-\langle v_\theta^2\rangle/\langle v_r^2\rangle$ for uncentered velocities relative to systemic; spherical symmetry gives equal tangential components. Turner's fitted **cylindrical** Jeans anisotropy and its posterior do not supply a bound on this quantity. [All numerical bounds and slit offsets](results.json)

## Why one stellar bin controls the exterior gravity

Keep the stellar tracer density $\nu$ fixed. For a radial source shell at $s$, define

$$
 W(r,s)=\exp\!\left(2\int_r^s\beta(t)\,dt/t\right),\qquad
 K(R,s)=\frac1s\int_0^{\sqrt{s^2-R^2}}
 W(r,s)\left(1-\beta(r)R^2/r^2\right) dz,\quad r^2=R^2+z^2 .
$$

Set $K=0$ for $R\ge s$. Integrating the spherical Jeans equation with its zero outer integration constant yields $j_r(r)=\int_r^\infty W(r,s)\nu(s)v_c^2(s)\,ds/s$. Projection and exchange of the positive integrals give half the second-moment numerator $J(R)=\int\nu(s)v_c^2(s)K(R,s)\,ds$. The isotropic kernel is $K_0(R,s)=\sqrt{(1-R^2/s^2)_+}$.

**Anisotropy of either sign.** If $|\beta|\le b<1/2$, set $z=y\sqrt{s^2-R^2}$. Since $r/s\ge y$, the integrating factor lies between $y^{2b}$ and $y^{-2b}$. Since the angular factor lies between $1-b$ and $1+b$, integrating in $y$ gives $L K_0\le K\le U K_0$ with $L=(1-b)/(1+2b)$ and $U=(1+b)/(1-2b)$. For $s\ge S$ take exterior budget $B=U\int_S^\infty\nu v_c^2\,ds$. At any sky position the exterior contribution lies between $B(L/U)K_0(R,S)$ and $B$.

**Anisotropy increasing towards radial orbits.** Suppose $0\le\beta(r)\le r^2/(a^2+r^2)$, without requiring that $\beta$ itself be monotone. For the interior shells,
$K_0/(1+s^2/a^2)\le K\le(1+s^2/a^2)K_0$. For the exterior, $W(z,s)$ decreases with $z$, and for $r^2=z^2+R^2$,

$$
 \frac{W(r,s)}{W(z,s)}\ge\frac{a^2+z^2}{a^2+r^2},
 \qquad 1-\beta(r)R^2/r^2\ge\frac{a^2+z^2}{a^2+r^2}.
$$

If $Z=\sqrt{s^2-R^2}$, monotonicity gives $\int_0^Z W(z,s) dz\ge (Z/s)\int_0^s W(z,s) dz$. Also $W(r,s)\le W(z,s)$ and the angular factor is at most one. Combining these inequalities, for every $s\ge S$,

$$
 Q(R)=\frac{\sqrt{(1-R^2/S^2)_+}}{(1+R^2/a^2)^2}
 \ \le\ \frac{K(R,s)}{K(0,s)}\ \le1 .
$$

Here $B=\int_S^\infty\nu v_c^2 K(0,s) ds$. This handles populations tending to completely radial orbits, including unbounded *ratios* of far-shell kernels to the isotropic kernel. Finite outer-bin speed forces a finite $B$ because the blurred $Q$ has strictly positive weight on that bin. The bounds require no exterior force amplitude, slope, cut radius or mass model.

**Combining the two telescope responses.** After MUSE seeing, every point in the outer footprint satisfies $J_{\rm out}\ge j_-(x)+B q(x)$, with $j_-$ the blurred lower interior bound, $q$ the blurred exterior factor. The two mass endpoints for the lower interior estimate are $(10.79,11.19)$ in $\log_{10}(M_\star/M_\odot)$; the upper uses $(10.94,11.32)$. Both acceleration endpoints come from the same independent SPARC receipt as in the parent test. Writing $I$ for the blurred brightness, the upper measured 247 km/s implies
$B\le\max_x[(247^2 I(x)-j_-(x))/q(x)]$. The central Keck squared speed is bounded above by $(\langle j_+\rangle+B)/\langle I\rangle$ with that instrument's seeing and rectangular aperture.

The optional centroid constraint is a finite extremum over luminosity-weighted radial bins. If $A_i$ and $q_i$ are the lower bounds on $j_-/I$ and $q/I$ in each cell and $u_i$ its upper radius, maximize $(247^2-\sum p_i A_i)/\sum p_i q_i$ over $p_i\ge0$, $\sum p_i=1$ and $\sum p_i u_i\ge0.85$. A linear-fractional functional with positive denominator attains a maximum at a vertex. These vertices are single bins or pairs bracketing 0.85 with the mean at 0.85. An independent linear program verifies the two-bin calculation.

The numerical implementation takes upper/lower endpoints on positive source-shell sums, sky annuli after Gaussian blurring, and a uniform rectangular grid covering each slit bin. Angular monotonicity holds for the **bounding** interior fields because they are positive mixtures of decreasing isotropic kernels. Gaussian convolution preserves radial monotonicity. This permits independent upper and lower bounds for each rectangular cell under every permitted center displacement. The code uses ordinary double precision with outward rounding for the reported thresholds; it does not certify special-function roundoff. [Calculation](measure.py) · [Independent tests](test_slit.py)

## A strong radial example and its held measurement

For an Osipkov–Merritt tracer with $a=0.4$ arcsec, take the upper Chabrier masses and acceleration, extend the frozen force along the full ray, and add the Newtonian force of an $8\times10^{10}\,M_\odot$ spherical shell at the outer ring. Inside that ring the added force vanishes. These illustrative parameters were selected while exploring the Keck central and MUSE outer moments. For this profile the angular integrand is $(s^2+a^2)(z^2+a^2)/[s(z^2+R^2+a^2)^2]$. Its closed-form integral decreases with $R$. The projection gives **287.25 km/s** for the Keck centre and **241.19 km/s** at the outer MUSE bin's luminosity-weighted centre. Both measurements bracket these values using the published centre-evaluation convention for the outer bin.

The MUSE centre provides a held check. This same tracer predicts **291.06 km/s** at the origin and a numerical lower sum of **290.37 km/s** across its entire 0–0.16 arcsec footprint. The upper observed bar is **282 km/s** after outward rounding. A finer source integration with a larger outer cutoff changes each reported speed by less than 0.05 km/s. This is a degeneracy of two projected second moments and a failure on the third. A physical stellar distribution function and a formation mechanism for the added shell are separate requirements. [Explicit projected control](strong_radial.py)

## Reproduction

From the repository root, with the normal dependencies installed:

```sh
export PYTHONPATH=code OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m cosmology.rar_deep_regime.double_ring_pressure.resolved_kinematics.slit_anisotropy.measure
python -m pytest -q code/cosmology/rar_deep_regime/double_ring_pressure/resolved_kinematics/slit_anisotropy
```
