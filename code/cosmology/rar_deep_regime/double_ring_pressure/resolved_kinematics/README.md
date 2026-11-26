# Resolved stellar spectra test the orbital field of the double-ring galaxy

In SDSS J0946+1006 the inner lens galaxy has both Hubble arcs and resolved MUSE stellar spectra. The fixed spherical Chabrier continuation used in the [two-ring pressure calculation](../README.md) can deliver too little stellar motion at the centre even with free exterior gravity. Assuming an equilibrium isotropic stellar tracer, fixed published light scales and the footprints below, the conservative upper sum is **262 km/s** at the centre, against a conservative **270 km/s** lower plotted reading. The lighter mass branch therefore fails to cover both stellar error bars under these assumptions. A more massive stellar calibration gives 278 km/s at the centre in its nominal full continuation. These spectra independently test the force used to infer a signed pressure contrast. They do not identify a collar stress.

## Observations and frozen inputs

Turner, Smith and Collett analysed 5.2 hours of VLT/MUSE integral-field spectroscopy (ESO 0102.A-0950), fitting the foreground stellar spectrum together with the luminous background arc. The black Figure 12 symbols are the stellar root-mean-square line-of-sight speeds, including resolved rotation through the raw second moment relative to systemic velocity. The colored curves are model fits, for which the authors used lensing normalization. We use only the black observations. The stellar spectra constitute an independent readout of the same galaxy. Their imposed lensing mass never enters this test. [Paper and Figure 12](https://arxiv.org/html/2401.08771v1) · [Journal publication](https://doi.org/10.1093/mnras/stae263) · [Original figure archive](https://arxiv.org/src/2401.08771v1)

[digitize.py](digitize.py) checks the SHA-256 of the 654 × 498 original figure, fits both coordinate scales from axis ticks, reads 41 identifiable black marker cores and verifies the two registered bar endpoints. All positions are in [muse.json](muse.json). The authors measured 53 bins and excluded nine in their outer-region dynamical analysis. Positions in one radial plot collapse different sky azimuths. Unresolved overlapping marker cores and bins outside the extracted window are excluded. The median comparisons below use this visible subset only. The outer comparison uses the marker nearest 1 arcsec by radius, independently of its speed.

| Stellar bin | Marker centre (arcsec, km/s) | Plotted radial bar (arcsec) | Plotted speed bar (km/s) | Allowed footprint (arcsec) |
|---|---:|---:|---:|---:|
| Central | 0.000, 275.98 | 0.000, 0.000 | 271.03–280.93 | 0–0.16 |
| Nearest 1 arcsec | 1.010, 239.48 | 0.797–1.214 | 233.71–245.46 | 0.64–1.37 |

Each speed endpoint is widened by two original raster pixels (0.413 km/s) and rounded outwards, with a further upward allowance on the outer bin to 247 km/s. Spatially we add two original raster pixels (0.0095 arcsec) and the half-diagonal of a 0.2 arcsec square MUSE spaxel (0.1415 arcsec) to the plotted width. As Figure 12 describes bin widths, using them to bound *all contributing spaxel centres* is an explicit footprint hypothesis. The two numerical bounds also apply if the bin masks, independently of how the bars were drawn, satisfy the displayed windows. The seeing model is the authors' circular Gaussian with FWHM 0.5 arcsec. The published error bars are compared individually. There is no new joint confidence level, remeasurement of the spectra, PSF uncertainty estimate or covariance extraction. [Observing and analysis methods](https://arxiv.org/html/2401.08771v1) · [ESO MUSE manual](https://www.eso.org/sci/facilities/paranal/instruments/muse/doc/ESO-261650_MUSE_User_Manual_8p2.pdf)

The two sphericalized F814W components have core radii (0.066, 0.082) arcsec, truncation radii (0.50, 6.05) arcsec and magnitudes (18.75, 17.15). The stellar mass, angular-distance convention and independent SPARC acceleration are inherited from the [parent calculation](../measure.py). Its frozen distance scale is 3.57692 kpc per arcsec; Turner et al. report 3.69 kpc per arcsec for their own cosmology. Our measured bins and profile inputs use angles throughout, with the frozen lens distance entering the Newtonian force. We use the frozen inner force $v_c^2=GM_*/r+\sqrt{G M_*a_0}$, with $M_*$ the spherical stellar mass enclosed at $r$. The light scales are fixed at the central published values for this test. The upper Chabrier case sets both logarithmic stellar masses to their upper endpoints (10.94, 11.32) and the independent acceleration to its upper bootstrap endpoint. The nominal masses are (10.85, 11.27) and (11.13, 11.52) for the respective Chabrier and Salpeter calibrations. Both terms increase with mass and acceleration, so the upper interior force covers all smaller component masses and accelerations at these fixed shapes. [Photometry and stellar populations](https://arxiv.org/abs/1111.4215) · [SPARC acceleration calibration](../../receipts/sparc_deep_regime_diagnostic.json)

## Bounding every nonnegative exterior force

Let $\nu(s)$ be the spherical tracer luminosity density. Write $I(R)$ for half its projected surface brightness, and $J(R)$ for half the line-of-sight projected second-moment numerator. In a steady spherical isotropic tracer with vanishing radial pressure at infinity, $j(r)=\int_r^\infty \nu(s)v_c^2(s)\,ds/s$. Projecting $j$ along the line of sight and exchanging the two nonnegative integrals gives

$$
 J(R)=\int_0^\infty \nu(s)v_c^2(s)K_s(R)\,ds,\qquad
 K_s(R)=\sqrt{(1-R^2/s^2)_+}.
$$

The measured quantity is $v_{\rm rms}^2=J/I$ after convolving both fields separately with the circular Gaussian. Write $I_\sigma$, $J_\sigma$ for these fields and $S=2.07$ arcsec for the cutoff. For $s\ge S$ we have $1-R^2/S^2\le K_s(R)\le1$. Gaussian seeing at sky position $x$ has $\mathbb E R^2=x^2+2\sigma^2$. Set $q(x)=1-(x^2+2\sigma^2)/S^2$. If $B=\int_S^\infty\nu v_c^2\,ds$, the blurred exterior contribution at $x$ lies between $q(x)B$ and $B$. Here $q$ is strictly positive over the entire outer stellar bin. Consequently an upper measured outer speed $U$ bounds $B$ independently of the exterior force profile,

$$
 B\le U^2\max_{x\in\mathrm{outer}}\frac{I_\sigma(x)}{q(x)},\qquad
 v^2_{\rm central}\le
 \max_{x\in\mathrm{centre}}\frac{J_{<,\sigma}^{\rm upper}(x)+B}{I_\sigma(x)} .
$$

These inequalities hold with arbitrary nonnegative aperture weights on the stated footprints. The observed moment is a ratio of two equally weighted integrals of $J_\sigma$ and $I_\sigma$. For instance, $J_\sigma\ge qB$ on the outer footprint implies $U^2\int I_\sigma\ge B\int q$, giving the first bound. The inner field in the numerator takes the two maximal component masses and acceleration. The interior term in the outer observable is discarded using its nonnegative sign. Exterior mass, slope and parametrization can all vary while the inequality persists. If exterior force were to yield a divergent $B$, the finite measured outer moment and $q>0$ would already exclude it.

Both $I$ and $J_<$ decrease with projected radius: the prescribed light decreases, and each positive Jeans shell has decreasing $K_s(R)$. Convolution with a decreasing circular Gaussian preserves this order (decompose the two profiles into concentric disks, whose intersection decreases with separation). On a radial subinterval the code bounds the numerator by its inner endpoint and each decreasing denominator by its outer endpoint.

[enclosure.py](enclosure.py) independently bounds the intrinsic moment by positive sums over 3,001 shells. On each shell it takes the largest enclosed mass and kernel, the largest tracer density, and the smallest radius in the Newtonian denominator. The central density supplies the finite first-shell bound. For a PSF at $x$, the squared displaced radius in units of $\sigma^2$ has a noncentral chi-square distribution with two degrees of freedom. Integrating annulus weights against left or right field endpoints, including the light tail beyond 6 arcsec, gives the numerical upper sum **261.74281 km/s**. [measure.py](measure.py) independently integrates the Gaussian Bessel kernel and the projected Jeans moment. Its finer estimate using the same sky-footprint partition is **260.303 km/s**. The mathematical enclosure uses monotonicity; the implementation evaluates standard special functions in double precision and has no machine-verified rounding-error certificate.

## Minimum missing interior orbital force

The same two bins set a positive target for any proposed correction. Let $a$ be kpc per arcsec and $D=\sup_{0<s<S}as\,[\delta v_c^2(s)]_+/G$ be the largest additional spherical enclosed orbital mass equivalent above the maximum frozen interior force, at any radius. This quantity bounds the positive increase in $v_c^2$ at angle $s$ by $GD/(as)$. The seeing-averaged $K_s$ is at most the probability that a photon falls inside the disk of radius $s$. The peak of the Gaussian bounds this probability by $\min(1,s^2/(2\sigma^2))$, even for a point mass at the origin. Thus, relative to the maximum interior numerator plus the same model-independent exterior bound, the positive addition to the central squared moment cannot exceed

$$
 \frac{GD}{a I_\sigma(0.16'')}
 \int_0^S \frac{\nu(s)}s\min\left(1,\frac{s^2}{2\sigma^2}\right) ds .
$$

To reach the 270 km/s lower reading the addition must cover at least $270^2-262^2$ in squared speed. The minimum blurred brightness over the central footprint and an upper sum on the integral give **$D\ge 8.43\times10^8 M_\odot$** (rounded further down to $8.4\times10^8$ for the paper). The same bound applies to the maximum enclosed mass of any supplementary spherical Newtonian source. For a nonlinear gravity law the quantity is the mass equivalent of the *total added orbital acceleration*, including any response to baryons. Stellar mass, ordinary compact objects and collar contributions are not decomposed by these two measurements. There is no determination of a unique mass profile or central black-hole mass.

## Measured comparisons

| Stellar mass and acceleration | Central full continuation | Inner Einstein ring, 1.43 arcsec | Median absolute residual (41 markers) | Central upper with free exterior, direct quadrature |
|---|---:|---:|---:|---:|
| Nominal Chabrier | 207.96 | 179.85 | 66.77 | 242.94 |
| Upper Chabrier | 228.49 | 195.75 | 48.62 | 260.30 |
| Nominal Salpeter | 278.09 | 230.67 | 8.33 | 303.81 |

All entries are km/s. The first three columns evaluate the continuation along the full line of sight at the plotted bin centres. The last column allows all exterior forces and the stated finite bin footprints. The comparison uses no fitted scale in this lens. The heavier mass calibration covers the central plotted bar; the full map contains azimuthal differences and overlap between bins, so the median residual provides no acceptance statistic. Tracer orbital anisotropy, spatial flattening, different light scales, an added central mass, a different seeing response or a remeasured bin mask require recalculation. Stellar tracer isotropy in Jeans' equation is separate from the anisotropy of the **total gravitating stress** in the [conservation constraint](../stress_conservation/README.md). A measured geometry and stellar motion constrain effective gravitational fields. Their attribution to an OPH collar additionally requires an independent material stress prediction at the ordinary coupling.

## Replay

From the repository root, with its Python dependencies installed:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code
p=cosmology.rar_deep_regime.double_ring_pressure.resolved_kinematics
curl -L https://arxiv.org/src/2401.08771v1 -o /tmp/turner-2401.08771v1.tar
python -m "$p.digitize" /tmp/turner-2401.08771v1.tar
python -m "$p.measure"
python -m "$p.enclosure"
python -m pytest -q code/cosmology/rar_deep_regime/double_ring_pressure/resolved_kinematics
```

[Measured markers and endpoints](muse.json) · [Numerical results](results.json) · [Independent Jeans and PSF checks](test_kinematics.py)
