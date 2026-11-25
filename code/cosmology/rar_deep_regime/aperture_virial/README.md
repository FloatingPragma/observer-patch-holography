# Radial stress from an Einstein ring and stellar spectra

A static spherical metric makes the difference between two mass readings a direct measure of the radial pressure integrated along a line of sight. The light reading is the mass inside the cylinder drawn by an Einstein ring. The orbital reading is the projection of the mass inferred from the radial force. We compare the two on 53 individual SLACS lenses while allowing the mix of radial and tangential stellar motion to vary freely with position. The calculation tests the continuation of the SPARC-calibrated square-root acceleration into bright stellar bodies. A native collar-stress prediction would have to supply the signed pressure column independently. [Matched input galaxies](sample.json) · [Bolton et al., Tables 4 and 5](https://doi.org/10.1086/589327)

## The pressure column

Use areal radius and the static spherical metric with $g_{tt}=-e^{2\Phi/c^2}$ and $g_{rr}=(1-2Gm/(rc^2))^{-1}$. To leading weak-field order the orbital mass is $M_o=r^2\Phi'/G=m+4\pi r^3p_r/c^2$. Light responds to the Weyl mass $M_w=(M_o+m)/2$. These formulas permit pressure comparable to the local density times $c^2$ while the metric perturbations stay small. The pressure is a component of the total physical stress. Assigning it to a collar requires subtracting other sources in the same metric. [Weak-field derivation by Faber and Visser](https://arxiv.org/abs/astro-ph/0512213)

For a spherical mass profile $M$, the cylinder of projected radius $R$ contains all of a shell with $r<R$ and the fraction $1-\sqrt{1-R^2/r^2}$ of an exterior shell. Integration by parts, including a central point mass if present, gives
$$
  M_{\rm cyl}(R)=\int_0^{\pi/2}M(R\sec\theta)\cos\theta\,d\theta.
$$
This requires a vanishing outer boundary term, for example $M(r)=O(r)$, and locally convergent integrals. Apply the linear projection to $M_o-M_w=2\pi r^3p_r/c^2$. With $z=R\tan\theta$ and $\overline\Sigma=M_{\rm cyl}/(\pi R^2)$,
$$
 P(R)\equiv\int_{-\infty}^\infty p_r(\sqrt{R^2+z^2})\,dz
 = c^2\,[\overline\Sigma_o(<R)-\overline\Sigma_w(<R)].
$$
An excess of optical over orbital mass therefore measures a negative pressure column under these conditions. A single circular speed at $R$ cannot fix $\overline\Sigma_o(<R)$ because stars and light sample different paths. Our spectral calculation supplies an explicit radial-force family.

Since $d\overline\Sigma/dR=-2\Delta\Sigma/R$, a difference of pressure columns uses the shear difference on a finite observed annulus:
$$
 P(R_1)-P(R_2)=2c^2\int_{R_1}^{R_2}
       [\Delta\Sigma_o(R)-\Delta\Sigma_w(R)]\,\frac{dR}{R}.
$$
An additive mass sheet in an inferred surface density drops out of this contrast. The measured ellipticity must be converted to shear with the appropriate source distances and reduced-shear correction. If the pressure decays and the Abel integrals converge, a resolved profile gives
$$
 p_r(r)=\frac{2c^2}{\pi}\int_r^\infty
 \frac{\Delta\Sigma_o(R)-\Delta\Sigma_w(R)}{R\sqrt{R^2-r^2}}\,dR .
$$
The strong-lens comparison below gives one pressure column per object conditional on the mass calibration, foreground model and external convergence. It does not invert a radial pressure profile from one ring.

## Stellar speeds with unrestricted orbital anisotropy

Assume steady, nonrotating spherical starlight with deprojected luminosity density $\nu$. Let $j_r=\nu\langle v_r^2\rangle$ and $j_t=\nu\langle v_\theta^2\rangle=\nu\langle v_\phi^2\rangle$ be the stellar *tracer* stresses. They weight the stellar spectrum. The unobserved stress $p_r$ in the metric equation above has different units and meaning. The Jeans equation in the given force field $v_c^2=r\Phi'$ is
$$
  (r^2j_r)'=2rj_t-r\nu v_c^2,\qquad
  0\le j_r(r)\le J(r)=r^{-2}\int_r^\infty s\nu(s)v_c^2(s)\,ds.
$$
The inequality uses $j_t\ge0$ and $r^2j_r\to0$ at infinity. [Spherical aperture dynamics](https://arxiv.org/abs/1401.4462)

Consider a circular aperture of radius $a$. At radius $r$ it accepts the two polar caps with $|\cos\theta|\ge u=\sqrt{\max(0,1-a^2/r^2)}$. Put $A=(1-u^3)/3$ and $B=1-u-A$. The line-of-sight kinetic numerator divided by $4\pi$ is $N=\int r^2(Aj_r+Bj_t)dr$. Eliminate $j_t$ with Jeans and integrate the $j'_r$ term by parts. For $r<a$ its coefficient after integration is zero. For $r>a$ the coefficient is $a^2/(2u)$. The two pieces agree in the integration-by-parts boundary term at $r=a$ and the singularity in $1/u$ is integrable. Thus
$$
 N=N_{\rm circ}+\frac{a^2}{2}\int_0^\infty j_r(\sqrt{a^2+z^2})\,dz,\qquad
 N_{\rm circ}=\frac12\int r^2B\nu v_c^2\,dr .
$$
The enclosed luminosity is positive and fixed. The lower moment takes $(j_r,j_t)=(0,\nu v_c^2/2)$ and the upper takes $(J,0)$. Both solve the moment equations. Convex mixtures give all intermediate moments. A nonnegative stellar distribution function can impose further restrictions. These are necessary orbital bounds encompassing radially varying velocity anisotropy.

A centred circular fibre convolved with an isotropic decreasing Gaussian has a decreasing radial acceptance function, as follows by superposing concentric disks and using monotonicity of their overlap with offset. Its layer-cake representation is a positive mixture of top hats. The same lower and upper bounds therefore apply with the measured-aperture weighting, including atmospheric blurring. We use the 1.5 arcsec SDSS fibre radius, a common seeing FWHM of 1.4 arcsec, the reported intermediate-axis effective radii and a spherical de Vaucouleurs light profile. The numerator is linear in $v_c^2$ for each bounding stress. A separate isotropic Jeans projection reproduces the 53 archived aperture checks. [Numerical implementation](measure.py) · [Independent tests](test_virial.py) · [Earlier isotropic prediction](../stellar_lensing.py)

## Comparing the same galaxies

The first test assigns a flat total circular-speed curve. Given the lens and source distances, a constant Weyl-to-orbit weight $q_{\rm tot}$ gives $\theta_E=2\pi q_{\rm tot}v_c^2 D_{ls}/(c^2 D_s)$ in radians. Each speed lies between the two tracer bounds. At nominal spectral speeds the lower bound from one lens is 1.242 while another lens requires at most 0.921. Allowing each speed two quoted standard errors and each isothermal-ellipsoid ring radius 2% gives the common necessary interval **[0.926, 1.160]**. All 53 admit a unit total weight in this enlarged box.

The anomalous-component test uses the independently calibrated $a_0=8.7504\times10^{-11}\,\mathrm{m\,s^{-2}}$, $v_*^2=G M_*(<r)/r$ and $v_A^2=\sqrt{G M_*(<r)a_0}$. The projected stellar mass follows the light. Writing stellar mass as $10^{11}m\,M_\odot$, both bounds on the stellar squared speed take the form $k_*m+k_A\sqrt m$. Invert the upper spectral moment at the measured lower speed and the lower moment at the measured upper speed to bracket $m$. Intersect this with an independently supplied stellar mass range wherever available.

Let $B$ and $A$ be the projected orbital masses per unit $m$ and $\sqrt m$ at the observed ring. The ring gives $M_E=c^2R_E\theta_E D_s/(4GD_{ls})$. Independently ordinary lensing for the stars and a constant optical weight for the anomaly give
$$
 q_A(m)=\frac{M_E-Bm}{A\sqrt m}.
$$
This strictly decreases for $m>0$. Each admissible mass interval yields its complete interval in the relaxed moment model. For variable pressure $q_A$ at one ring is an aperture-effective ratio, with the positive orbital projection providing the weights along the ray. The pressure-column identity continues to apply. We use the frozen spherical square-root prescription throughout the stellar body and along the ray. This extension, the stellar mass-to-light shape, isolation, projection and fixed distance calibration are physical assumptions. The comparison supplies no inner-to-outer transfer for a law specified only on a deep, baryon-enclosing annulus.

* With free stellar normalization and each quoted speed allowed two standard errors, the 53 spectra and isothermal-ellipsoid rings give a shared necessary $q_A$ interval of **[0.955, 1.441]**. The light-traces-mass reconstruction of the same arcs gives **[0.925, 1.342]**. Unit optical response meets every interval in both reductions.
* For the 18 exact matches to the gravitating-mass experiment with archived stellar-population estimates, allow each mass between the Chabrier estimate minus one quoted log error and the Salpeter estimate plus one. Spectra narrow ten individual mass boxes. The shared interval is **[1.038, 1.233]** from isothermal-ellipsoid imaging and **[0.943, 1.174]** from light-traces-mass imaging. The two active endpoints here come from the population brackets. Spectra and arcs alone on those same 18 give [0.955, 1.755] and [0.925, 1.715]. The two image reductions share the arcs. [Stellar population measurements](https://arxiv.org/abs/0911.2471) · [Source values](../data/slacs_weak_lenses.dat)
* J2300+0022 has a conditional projected radial stress $P/c^2$ between **−1872 and −52 $M_\odot\,\mathrm{pc}^{-2}$** at 4.534 kpc in the isothermal-ellipsoid reduction. At 4.424 kpc the alternative reduction permits **−1796 to +79** in the same units. Imaging and stellar calibration determine whether these data require a nonzero pressure column. Nonstellar baryonic pressures, external convergence, gas, non-spherical backgrounds and intervening matter would enter a source attribution.
* In the 18-lens mass-box comparison, scaling the nominal isothermal-ellipsoid ring radii down by 2% gives [0.960, 1.167]. Using the upper SPARC scale endpoint gives [0.953, 1.132]. Freezing the anomalous orbital mass beyond ten effective radii gives [1.072, 1.271] while preserving the enclosed anomalous mass at the edge. These are separate sensitivity cases. The chosen boxes have no joint confidence level. [Results for each lens and all seven controls](results.json)

The catalogue selection and two full-table digests are shared with the isotropic comparison. Every included galaxy is a single early-type grade-A system with a successful image fit and measured effective radius, stellar dispersion and error. Detection and follow-up of strong lenses affect this sample. Double the radial grid and increase the angular quadrature order from 96 to 160 to reproduce all common weight endpoints to $6\times10^{-5}$. The spectra measure luminosity-weighted stellar dispersions. Rotation, galaxy-specific seeing, profile variation and observational covariances can alter this comparison.

From the repository root, run

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m cosmology.rar_deep_regime.aperture_virial.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/cosmology/rar_deep_regime/aperture_virial
~~~

Supplying `--catalogue /path/to/bolton.tsv` also verifies the archived 53-lens input against the checked source tables. [VizieR source tables](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J%2FApJ%2F682%2F964&-out.all&-out.max=unlimited)
