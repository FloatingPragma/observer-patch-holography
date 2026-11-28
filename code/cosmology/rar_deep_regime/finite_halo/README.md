# Orbits, clocks and bending on a finite spherical halo

Consider the OPH conditional exterior density proportional to 1/r² and compact baryonic mass M_b. In geometric units write the mass length as m=b+hr inside an outer edge R_t, where b=GM_b/c² and h=√(GM_b a₀)/c². Matching to vacuum outside the edge gives a finite total mass. A specific static choice is zero radial pressure. In the metric with g_tt=−exp(2Φ) and g_rr=(1−2m/r)⁻¹ the radial Einstein equation gives (v_circ/c)²=rΦ′=m/(r−2m). Conservation fixes each tangential pressure to ρ(v_circ/c)²/2, where ρ is the energy density. Exterior clock ratios and geodesic deflections use that same metric. The [executable metric and shear reader](measure.py) treats the inner baryonic region as outside its optical domain.

This static exterior is an effective constitutive example. The registered collar records supply no radial pressure law, outer edge, inner matching mass or cosmological evolution. On a scale-free annulus with p_r=wρ, the radial Einstein equation instead gives (v_circ/c)²=h(1+w)/(1−2h). The enclosed mass alone cannot resolve w. Zero radial pressure is closely related to an Einstein cluster at the level of the stress tensor. No constituent particle population is assigned in this comparison. [Lake](https://arxiv.org/abs/gr-qc/0607057) · [Faber and Visser](https://arxiv.org/abs/astro-ph/0512213) · [Checked source identities](../../../../Lean/ObserverPatchHolography/EinsteinBranch/HaloStress.lean)

## Which radial distributions do the measured orbits allow?

The source-counted pressureless sphere has eight joint records per eleven ordinary ones. Applying its total factor of 19/11 at every measured tracer radius fails in 40 of 136 SPARC galaxies at five quoted speed errors. Even if each galaxy chooses independent nonnegative weights for gas, disk and bulge, their three radial shapes fail in 23 galaxies. Widening each ring's speed interval by another 5 km/s leaves seven exclusions. At 10 km/s NGC 1003 still requires at least 2.36% of its outer radial acceleration beyond these templates, compared with 34.51% at the original five-error widths. These are catalogue interval tests conditional on the supplied tracer templates and geometry. They locate a deficit a record-derived source would have to supply. [Catalogue reader](tracer_shapes.py) · [Exact witnesses for all 2,176 comparisons](orbit_certificates.json)

## How much can a finite outer edge change the shear?

For the weak-field projection of a sharply truncated excess A/r², with y=R/R_t, the ratio of differential surface density to the infinite-profile value πA/R is

$$
F(y)=
\begin{cases}
[4y/(1+\sqrt{1-y^2})+2\arccos(y)]/\pi,&0\le y\le1,\\
4/(\pi y),&y\ge1 .
\end{cases}
$$

The maximum is 4/π at the edge. Positive mixtures of cutoffs obey the same upper bound. In six isolated bins per stack, the frozen full-sample SPARC scale of 8.7504 × 10⁻¹¹ m/s² and the published bin centres give tangential measurements 4.00 and 2.78 conditional statistical errors above the maximum for a single common distribution of edge accelerations on KiDS-bright and GAMA respectively. Maximizing over unknown positive within-bin weights as well reduces the gaps to 0.98 and 2.09 errors. The cross-shear controls are −0.29 and −0.04 errors from zero. Mass calibration and within-bin sampling affect this comparison. [Frozen protocol](PROTOCOL.md) · [Retrospective common-edge calculation](COMMON_EDGE_PROTOCOL.md) · [Six-bin readings](results.json) · [KiDS public data release](https://kids.strw.leidenuniv.nl/sciencedata.php)

## Can one profile fit the dwarf and the lensing stacks?

The outer observed ring of AGC 114905 is at 9.95 kpc, with circular speed 26.94 ± 1.83 km/s and total baryonic mass (1.47 ± 0.14) × 10⁹ M☉. Assign all its observed squared speed to the positive spherical anomalous component. This gives a conservative upper anomalous mass fraction of 1.142 at the measured radius for the catalogue central values. A universal sharp edge then requires g_t≥6.71 × 10⁻¹¹ m/s² at the frozen SPARC scale. In both lensing samples all enrolled projected radii are beyond that edge, even at the lower scale and the generous settings below. Outside a spherical edge the column contribution of the halo itself vanishes and its enclosed mass obeys ΔΣ=M/(πR²). The bound on this mass ratio transfers independently of a₀. [Transfer protocol](TRANSFER_PROTOCOL.md) · [Published dwarf measurement](https://arxiv.org/abs/2404.06537)

A broader test allows any shared normalized enclosed profile f(x)=M_A/M_b, x=R/√(GM_b), for which f is nondecreasing and f/x is nonincreasing across the compared radii. A positive thin spherical shell outside projected R contributes in units of its mass divided by πR² the amount −(1−c)²/(2c), where c=√(1−R²/r²). Thus ΔΣ≤M(<R)/(πR²), even when some material lies outside the lensing aperture. The dwarf bounds f(x) by f_max max(1,x/x_d). Combining the two bounds needs no circular-speed asymptote, a₀ or chosen edge. It assumes a universal mass-normalized profile over these radii, ordinary weak-field stress optics, the stated compact baryonic lens estimate and nonnegative baryonic support in the measured disk. [Profile protocol](PROFILE_PROTOCOL.md) · [Cross-data calculation](transfer.py) · [Bounds by bin](transfer.json)

| Survey | Measured weighted ΔΣ | Central shared-profile upper bound | Generous sharp-edge upper bound | Generous shared-profile upper bound |
|---|---:|---:|---:|---:|
| KiDS-bright | 17.30 ± 0.69 | 5.04 | 9.94 | 11.36 |
| GAMA | 22.73 ± 2.99 | 4.97 | 9.74 | 11.18 |

Units are M☉/pc². The upper bounds allow arbitrary nonnegative sampling weights within each published acceleration bin. The generous calculation lowers the dwarf mass by two quoted errors, raises distance and speed by two quoted errors, also permits inclination 27° instead of 31°, and increases each lens mass by 0.2 dex. Those settings do not represent a joint probability interval. The quoted errors propagate the published shear covariance at fixed sample membership. The two lens samples use KiDS imaging. The dwarf ring lies inside the measured baryonic disk. These results disfavor extending the same spherical sublinear profile to that ring and the outer lens apertures under the tested matter and optical calibrations. The model may have a different matching radius, mass distribution or environmental response. The inputs do not identify a physical collar mechanism for such a change.

The [regression checks](test_halo.py) differentiate the metric to recover its stress, test Schwarzschild and scale-free light deflection, independently integrate spherical columns, and replay the data bounds. With the public tar archive downloaded from [KiDS](https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar), run from the repository root as follows.

```sh
PYTHONPATH=code python -m cosmology.rar_deep_regime.finite_halo.measure /path/to/brouwer2021_rar.tar --output /tmp/halo-shear.json
PYTHONPATH=code python -m cosmology.rar_deep_regime.finite_halo.transfer --output /tmp/halo-transfer.json
python -W error -m pytest -q code/cosmology/rar_deep_regime/finite_halo
```

## Reverse calibration from shear

With the SPARC normalization withheld, the six-bin shear mean at catalogue masses and bin centres gives untrimmed spherical scales of 1.530 × 10⁻¹⁰ m/s² for KiDS-bright and 2.813 × 10⁻¹⁰ m/s² for GAMA. At the measured outer AGC 114905 ring, assigning all the rotation to an untrimmed compact-source anomaly gives a central upper scale of 2.699 × 10⁻¹² m/s². Incorporating the measured disk force in the local acceleration prescription instead gives 2.990 × 10⁻¹³ m/s². With the generous dwarf speed and mass allowances the spherical upper scale becomes 9.195 × 10⁻¹² m/s². Even allowing the shear the upper bin and mass settings, its maximum edge enhancement and a two-statistical-error downward shift, the respective lower scales are 3.139 × 10⁻¹¹ and 3.649 × 10⁻¹¹ m/s². Independent outer-edge changes can suppress the dwarf speed for a given scale; the profile transfer above tests that alternative for a common family. [Calibration protocol](CALIBRATION_PROTOCOL.md) · [Shear inverse and dwarf readings](calibration.json)

The microscopic dictionary a₀=G n² c_var identifies one product of proper cut density n and per-source variance coefficient c_var. Scaling n by a positive factor and dividing c_var by its square preserves all the spherical profile predictions. A calibration of that product from orbits or shear does not count or locate observer patches. Independent regional readings must distinguish those microscopic quantities and their coupling to gravity. [Finite nonidentifiability theorem](../../../../Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean)
