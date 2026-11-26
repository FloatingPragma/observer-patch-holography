# Pressure between the two Einstein rings of J0946+1006

The two sources behind this galaxy sample the primary lens at 5.115 and 7.404 kpc. On a centered spherical weak-field continuation of the SPARC acceleration law, the lighter of two published stellar-population models requires the inner line of sight to carry a **smaller integrated radial pressure** than the outer one. Even removing the entire halo outside the second ring cannot eliminate the deficit. The heavier stellar model admits zero pressure. Both statements depend on interpreting the observed partial arcs as aligned circular critical curves and on which radial orbital force the stars generate. [Input measurements and primary source fits](https://arxiv.org/abs/0801.1555) · [Stellar light and masses](https://arxiv.org/abs/1111.4215) · [Outer-source spectroscopy](https://doi.org/10.1093/mnras/stab1399)


Static radial stress conservation turns the two aperture contrasts into a conditional lower bound on the weighted excess of tangential over radial pressure. Its aligned Chabrier value divided by $c^2$ is 0.007570 solar masses per cubic parsec. [Anisotropy bounds, stellar controls and assumptions](stress_conservation/README.md)

## Two lens planes

Let $C_i=c^2D_i/(4\pi G D_lD_{li})$ be the primary critical surface density for each source. Write $R_i=D_l\theta_i$ and $b_{12}$ for the **reduced radial deflection of the second-source ray by the intervening first source**, in the same angular units as $\theta_2$. It includes the $D_{12}/D_2$ factor. Allow an unknown constant surface-density sheet $S$ on the primary lens plane external to the galaxy dynamics. The aligned thin-lens equations are

$$
\bar\Sigma_{w1}+S=C_1,\qquad
\bar\Sigma_{w2}+S=C_2(1-b_{12}/\theta_2).
$$

A centered positive secondary gives $b_{12}\ge0$ if the outer ray reaches it on the same side of the axis. A monotonically decreasing primary mean Weyl density between the rings suffices to keep the ray on this side. Off-axis sources, external shear, variable foreground convergence and further lens planes can change the relations. The spectroscopy fixes three redshifts. No assignment of a primary optical power law is required.

Let $\Delta\bar\Sigma=\bar\Sigma(<R_1)-\bar\Sigma(<R_2)$ and $P_i=\int p_r(\sqrt{R_i^2+z^2})\,dz$. The static spherical weak-field relation $P_i/c^2=\bar\Sigma_{oi}-\bar\Sigma_{wi}$ uses the projection of the radial **orbital** mass $M_o=rv_c^2/G$. Subtracting the two measurements eliminates $S$ exactly:

$$
 (P_1-P_2)/c^2=\Delta\bar\Sigma_o-(C_1-C_2)-C_2b_{12}/\theta_2
 \ \le\ \Delta\bar\Sigma_o-(C_1-C_2).
$$

A negative upper bound excludes both zero radial pressure and any radial pressure decreasing outwards on the full lines of sight. It need not make either individual pressure column negative. The primary sheet cancellation concerns one physical sheet at one redshift. Other changes to the inferred deflections need independent calibration. [Pressure-column derivation](../aperture_virial/README.md)

## Why the far halo cannot rescue the force law

A spherical unit-mass shell at $r$ puts the fraction $f_i(r)=1-\sqrt{\max(0,1-R_i^2/r^2)}$ inside cylinder $i$. Its weight in $\Delta\bar\Sigma_o$ is $W=f_1/(\pi R_1^2)-f_2/(\pi R_2^2)$. For $r\ge R_2$, $f_i/R_i^2=1/[r^2(1+\sqrt{1-R_i^2/r^2})]$, so $W<0$. On the interior $W\le W_0=(R_1^{-2}-R_2^{-2})/\pi$, with equality at $r\le R_1$. Thus the upper limit at fixed interior profile over **every nonnegative exterior orbital mass measure** is achieved by holding the enclosed mass fixed at $M_o(R_2)$. This positivity concerns the effective orbital mass density $M'_o/(4\pi r^2)$ and must be checked if stresses are allowed to produce a negative effective density. If this truncated upper bound falls short by $L$, any added, nonnegative effective orbital mass profile that restores a nonnegative pressure contrast must put at least $L/W_0$ inside $R_2$. Additional baryons would modify both the direct force and the square-root term.

Set $s=\sqrt{1-(R_1/R_2)^2}$ and $r(x)=[(1-x^2)/R_1^2+x^2/R_2^2]^{-1/2}$. Projection by parts of the truncated profile gives

$$
 D_{\max}=\frac{s}{\pi}\left[
 \frac{1}{R_1^2}\int_0^1 M_o(r(x))\,dx
 -\frac{M_o(R_2)}{R_2^2(1+s)}\right].
$$

Only the interior force enters. We use $M_o(r)=B(r)+r\sqrt{a_0 B(r)/G}$ out to $R_2$, where $B$ is the cumulative stellar mass. The compact exterior theorem alone does not specify this interior continuation.

## Bounds from the observed profiles

The two components of the lens have sphericalized truncated pseudo-isothermal profiles from Table 3 of Sonnenfeld et al. For angular core and cutoff radii $c,t$ their spherical mass fractions at angular radius $v$ are $F(v)=2[t\arctan(v/t)-c\arctan(v/c)]/[\pi(t-c)]$. The corresponding projected fraction is $1-[\sqrt{t^2+v^2}-\sqrt{c^2+v^2}]/(t-c)$. Each fraction increases with $v$ and decreases with either scale. For $F$, the numerator divided by $t-c$ averages the derivative of $h(w)=w\arctan(v/w)$, with $h''(w)=-2v^3/(w^2+v^2)^2<0$.

The continuous interval bound uses three additional facts.

1. $D_{\max}$ decreases with $R_1$ and increases with $R_2$ on this source class. For nonincreasing positive stellar density $\rho$ put $h=d\log B/d\log r$. The derivative of the anomalous density, divided by $\sqrt{a_0B/G}/(4\pi r^3)$, is $-1-(h-2)^2/4+(h/2)d\log\rho/d\log r<0$. The truncated combined density is radially nonincreasing, so its mean projected density decreases with the inner radius. At fixed $R_1$ the derivative with respect to $R_2$ is $[2 M_o-R_2 M_o's/(1+s)]/(\pi R_2^3)>0$ because $R_2 M_o'\le 3M_o$.
2. The mass coefficient of each population is positive. At the extremal radii the verified bounds give $F_i(R_1)/F_i(R_2)>R_1/[R_2(1+s)]$. For the square-root contribution use $\partial_{m_i}\sqrt B=F_i/(2\sqrt B)$. On the integrated segment bound $r\ge R_1$, $F_i(r)\ge F_i(R_1)$ and $B(r)\le B(R_2)$. The stated inequality then proves the derivative of the truncated contrast is positive. It also implies positivity for the direct stellar term. The two individual maximum masses thus bound every mixture in the box.
3. The mass integrand increases in $x$. The implementation uses 192 right endpoints and lower bounds on the subtracted boundary term. It partitions each cutoff radius into two intervals, evaluates the monotone core and cutoff envelopes, bounds the acceleration term by the upper SPARC bootstrap scale, and uses directed interval arithmetic. For the flat $H_0=70$, $\Omega_m=0.3$ distance calibration the trapezoid error uses $|d^2(1/E)/dz^2|<60$ on $0\le z\le2.035$. Computed endpoints are rounded outwards. [Implementation](measure.py) · [Independent shell and two-plane checks](test_pressure.py)

The redshifts are $(0.222,0.609,2.035)$. The published double-ring radii are $1.43\pm0.01$ and $2.07\pm0.02$ arcsec. In arcseconds the two stellar core radii are $0.066\pm0.010,\ 0.082\pm0.010$, and the cutoff radii are $0.50\pm0.05,\ 6.05\pm0.10$. The Chabrier logarithmic component masses range independently within $[10.79,10.94]$ and $[11.19,11.32]$ in solar units, from Table 8's individual uncertainties. Shape and mass ranges, treated as rectangular, are not a joint confidence region. We use the full independently calibrated SPARC scale interval $[7.3902,10.3793]\times10^{-11}\ {\rm m\,s^{-2}}$. These inputs give $C_1-C_2=1055.58\ M_\odot\,{\rm pc}^{-2}$ at the stated distances.

| Ring selection | Upper $(P_1-P_2)/c^2$, $M_\odot\,{\rm pc}^{-2}$ | Additional interior effective mass if $P_1\ge P_2$ |
| --- | ---: | ---: |
| Each quoted ring error | $-114$ | $\ge1.72\times10^{10}\ M_\odot$ |
| Twice each ring error | $-76$ | $\ge1.10\times10^{10}\ M_\odot$ |
| Include the two Bolton inner-ring fits down to $1.38$ arcsec, double outer error | $-10$ | $\ge1.4\times10^9\ M_\odot$ |

The alternate inner angles are $1.38$ and $1.39$ arcsec from two reductions of the same Hubble images. If departures from aligned circular lensing shift the two dimensionless lens equations by at most $\epsilon$, the pressure contrast uncertainty is at most $(C_1+C_2)\epsilon$ before other systematics. A differential geometry error at the percent level can therefore matter. The stellar components themselves have fitted projected axis ratios 0.66 and 0.71 with different orientations. The arc data in this calculation do not independently bound the error of their spherical approximation. [Matched catalogue](../aperture_virial/sample.json)

Using the nominal Salpeter component masses instead, $\log_{10}(m_i/M_\odot)=(11.13,11.52)$, and continuing the model along the full ray gives an explicit zero-pressure solution. The inner ring requires a primary sheet of $0.1452\,C_1$. The second ray reaches the intermediate plane at $+0.656$ arcsec and needs a positive reduced bend of $0.10775$ arcsec. A centered singular isothermal secondary would give this bend at $80.76$ km/s. Gavazzi et al. infer $59\pm13$ km/s from the secondary luminosity and an externally calibrated Tully–Fisher relation. This is a model-dependent comparison, separate from a direct dispersion spectrum. The primary does have an independent SDSS fibre measurement, $263\pm21$ km/s. Weighting the two measured F814W luminosities as tracers, the untruncated Salpeter force gives an isotropic spherical Jeans dispersion of 259.36 km/s. The full aperture bounds with arbitrary radial variation of stellar anisotropy are [227.06,307.78] km/s. The nominal Chabrier field gives 196.92 km/s under isotropy and [170.98,236.90] km/s with unrestricted anisotropy. These calculations use the same 1.5 arcsec fibre radius and common 1.4 arcsec seeing FWHM as the 53-galaxy comparison. The pressure test constrains a stellar calibration and an effective force law together. It supplies no microscopic collar-pressure measurement. [All values and controls](results.json)

## Optical mass in two circular apertures of an elliptical lens

Collett and Auger fit the two sets of Hubble arcs with pixelated sources, an elliptical power-law primary, external shear and an isothermal secondary. Their Table 1 gives the primary intermediate-axis Einstein scale $b=1.397\pm0.001$ arcsec for the **first** source plane, projected density slope $\eta=1.027^{+0.023}_{-0.025}$, and projected mass axis ratio $q=0.946^{+0.009}_{-0.005}$. These describe the total optical mass fit. The flatter stellar components belong to the orbital model. The fit allowed the relative lensing-distance factor to vary. Multiplication by $C_1$ uses the spectroscopic first-source and primary redshifts and the distance calibration above. [Hubble image model, Table 1 and note 4](https://doi.org/10.1093/mnras/stu1190)

For elliptical radius $s^2=qx^2+y^2/q$, the fitted convergence is $\kappa=(2-\eta)(b/s)^\eta/2$. Its mean inside the ellipse $s<b$ is unity. Integration over the **circular** sky aperture $x^2+y^2<\theta^2$ instead gives the primary optical contrast

$$
W=\Delta\bar\Sigma_w=C_1 Q(q,\eta)\left[(b/\theta_1)^\eta-(b/\theta_2)^\eta\right],\qquad
Q=\frac1{2\pi}\int_0^{2\pi}(q\cos^2\phi+q^{-1}\sin^2\phi)^{-\eta/2}\,d\phi .
$$

A uniform primary sheet cancels in the difference of means and the external shear has zero convergence. The reconstruction includes the secondary, so the inferred primary $W$ already accounts for its fitted deflection. The nominal value at the discovery radii is $1165.79\ M_\odot\,{\rm pc}^{-2}$. Form a rectangle of the three quoted marginal parameter ranges. Jensen's inequality gives $Q\ge[(q+q^{-1})/2]^{-\eta/2}$, bounded below using the smallest $q$ and largest $\eta$. The square bracket increases with $b$ and $\theta_2$ and decreases with $\theta_1$. At these extremal radii it increases with $\eta$ at fixed $b=b_{\min}$, since $t e^{-\eta t}$ increases for $0\le t=\log(\theta/b)<1/\eta$. These bounds share the same apertures with the spherical orbital upper limit $D=\Delta\bar\Sigma_o$.

| Angular range | Lower $W$ | Upper $D$ | Lower $W-D$ |
| --- | ---: | ---: | ---: |
| Quoted ring errors | 1090 | 941 | 149 |
| Doubled ring errors | 1040 | 979 | 61 |
| Extended inner range, $1.38$ to $1.41$ arcsec, doubled outer error | 1143 | 1045 | 98 |
| Extended inner range, $1.41$ to $1.45$ arcsec, doubled outer error | 1040 | 979 | 61 |

All entries are outward-rounded surface densities in $M_\odot\,{\rm pc}^{-2}$. The result compares an observed-image optical fit with the spherical Chabrier force continuation. Interpreting $D-W$ as a radial-pressure column difference additionally requires the spherical stress identity. The parameter rectangle has no assigned joint confidence level. Both optical reductions use the same arcs. A physical lensing bound requires a lower limit on the primary mass contrast $W$ at matched circles after allowing for the primary profile, secondary mass, lens centering and line-of-sight degeneracies. A general mass-sheet transformation also rescales the primary. A fit with separate stellar and halo components favors heavier Salpeter-like stellar masses. [Composite lens model and parameter uncertainties](https://arxiv.org/abs/2602.20889) · [Elliptical aperture calculation](elliptical.py) · [Replayed bounds](elliptical.json)

For the aligned calculation, write dimensionless errors $e_i$ on the right sides of its two lens equations, giving $C_1(1+e_1)$ and $C_2(1-b_{12}/\theta_2+e_2)$. The change that can erase a deficit is $H=C_2 e_2-C_1 e_1$. If the tabulated guaranteed aligned deficit at zero secondary bend is $L$, a nonnegative pressure contrast requires $H\ge L+C_2 b_{12}/\theta_2$. On the box $|e_i|\le\epsilon$, the exact maximum of $H$ is $(C_1+C_2)\epsilon$, attained at $e_1=-\epsilon$, $e_2=+\epsilon$. Thus the smallest per-plane error envelopes capable of erasing the **guaranteed margins** of 114, 76 and 10 are respectively 1.7532%, 1.1688% and 0.1538% at $b_{12}=0$. A positive reduced radial bend raises the threshold. The actual force residual can exceed the conservative guaranteed margin. A physical sheet produces $C_2e_2=C_1e_1$ and cancels.

From the repository root, run

~~~sh
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m cosmology.rar_deep_regime.double_ring_pressure.measure
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m cosmology.rar_deep_regime.double_ring_pressure.elliptical
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/cosmology/rar_deep_regime/double_ring_pressure
~~~
