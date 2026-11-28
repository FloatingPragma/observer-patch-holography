# Deep-regime comparison for the anomalous-density profile

`rar_deep_regime_diagnostic.py` compares the dark-matter paper's conditional
deep-regime laws, `a_A = sqrt(a_b a_0)` and `v^4 = G M_b a_0`, with the SPARC
rotation-curve sample (snapshot under `data/`, see `data/PROVENANCE.md`).

A shared-threshold extension of the native binary source gives a profile-based test on ten galaxies and an amplitude comparison with four lensing stacks. Direct stellar distances, outer speeds and gas maps constrain the same two galaxies. The mass required by lensing exceeds one layer of sharp written records even at full delivery. [Resolved boundary test](boundary_source/README.md) · [Native single-layer mass calculation](../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/README.md)

The law has one dimensional constant `a_0` that the OPH source does not fix.
The diagnostic therefore reports

* the anomalous-excess exponent and the fitted `a_0` in the total model
  `g_model = g_bar + sqrt(g_bar a_0)` on fixed absolute subsets
  `g_bar < f (1.2e-10 m/s^2)`, for `f` in `{0.3, 0.1, 0.03}`;
* the baryonic Tully-Fisher exponent and the `a_0` implied by its
  normalisation;
* galaxy-bootstrap intervals for both constants, which the conditional
  density profile requires to agree;
* the scatter of the additive all-gradient extension
  `g = g_bar + sqrt(g_bar a_0)` and its Solar-System anomaly, which the
  paper's high-gradient branch is required to remove.

Status: `diagnostic_postdiction_fixed_absolute_cuts_one_fitted_constant`,
`physical_claim: false`,
`source_derived_output: false`. The exponents are the parameter-free content;
`a_0` is a fitted comparison value. The result is a tension under fixed
mass-to-light ratios and equal retained-point weighting, not a full likelihood
or a verdict on the generic OPH dark-sector theorem. The BTFR side uses the
catalogued finite-radius `Vflat` as an asymptotic-speed proxy without a
baryonic subtraction; resolving that known bias is part of the required joint
likelihood.

```bash
python3 code/cosmology/rar_deep_regime/rar_deep_regime_diagnostic.py \
  --output code/cosmology/rar_deep_regime/receipts/sparc_deep_regime_diagnostic.json
python3 -m pytest code/cosmology/rar_deep_regime
```

## Observer protocol for the dark-gravity prediction

The target is a joint distribution of recorded orbital speeds, image shear and cosmological observables. The physical hypothesis is that gravitational measurements respond to stress carried by imperfectly recovered overlap records. A catalogue entry for distance, the number of telescope exposures, and a hidden simulation clock supply no measurement of the recovery defect. OpenWave row 19 has a separate acceptance test for a stable neutral massive state. [Dark-sector paper](../../../cosmology/oph_dark_matter_paper.pdf) · [Official state criterion](https://github.com/openwave-labs/openwave/blob/f8fd02b2e0d96cfa0d868259daec2c9822375d34/MODELS.md#L121)

A mechanism-neutral comparison asks whether one process predicts both the magnitude and the location of the excess on held systems. Two separate written collars transmit a measurable correlation to a flying receiver at fixed ordinary occupancy. Their connected logarithmic charge has opposite signs under two faithful references. The local projected acceleration has an exact initial null, and six stations' conditional return times fail a three-dimensional distance test even in the vacant control. The source test therefore requires an independently prepared physical reference for the cut. The geometry experiment requires directions, clocks and rulers validated in the reference region before comparing weak source perturbations. Match ordinary and collar stress profiles, retain unsuccessful probes, and measure area at fixed proper volume under refinement. [Two-write measurements and proofs](../../../model/experiments/interval_echo_einstein_station/README.md) · [Observer gravity benchmark](../../../model/evaluation/EVIDENCE.md#patch-10)

The small-ball normalization gives a direct consistency test. Let $E(ℓ)$ be a signed collar-energy variation in one frame and $ℓ$ the proper radius of the independently measured local ball. The proposed stress is proportional to $E(ℓ)/ℓ^4$. If the whole collar variation is an error of order $o(ℓ^4)$, it gives zero limiting stress. A nonzero source requires a leading coefficient that agrees across translated observer frames, while the residual after subtracting that coefficient is $o(ℓ^4)$. Calibrate proper radius, boost reference, ordinary stress and the gravitational coupling on the same preparations. A mean return time that fails its vacuum geometry control cannot supply $ℓ$. The galactic profile radius and the local ball radius enter different tests. [Einstein derivation](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) · [Dark-source criterion](../../../cosmology/oph_dark_matter_paper.pdf)

1. **Specify the patch statistic.** Prepare or generate the reference and the record-recovery map in finite patches. Have overlapping observers read the same retained cut in a common frame and count it once. Measure the signed mean, its full variance and the variance between physically recorded histories. Repeat with an independent tag and with full recovery. Record which statistic the stress law couples to. Two ensembles with identical untagged states can have different between-history variances. A source that uses that difference requires the history record in the physical source. [Tagged source and complete readouts](../../../model/experiments/interval_echo_einstein_station/backreaction.json)
2. **Calibrate the response with ordinary matter.** Use material clocks, held null directions and an area and proper-volume instrument on one refinement family. Vary ordinary matter to fix the coupling and vary the collar statistic independently at fixed ordinary energy. Predict what a separate geometry instrument reads, including its quiet and failure outcomes. Repeat at different radii with two independent and two correlated sources. Measure the number of distinct physical cuts per proper length and the per-cut stress variance in mass units. The deep branch then predicts the normalization from $a_0=G n^2 c$. In a simulation the source can be fully specified from independent inputs even if a telescope cannot directly resolve one microscopic cut. [Einstein requirements and local tests](../../../model/evaluation/EVIDENCE.md#g6)
3. **Test the field seen by the same galaxies.** Before looking at comparison residuals, specify the outer radial domain, the baryonic mass distribution, allowable geometry and kinematic uncertainties, environmental selection, and how the two metric potentials govern orbits and light. Forward-model line-of-sight velocities and tangential shear at matched physical radii with the measured lens–source pair weights and their joint covariance. Use the same tracer and radial selection rule for every eligible galaxy. An inner velocity field may constrain a different disk orientation than the outer gas. Cross shear and independently observed velocity fields control measurement systematics. For a spherical $1/r^2$ excess with negligible slip extending along the line of sight, the ratio below equals one without fitting $a_0$ or the inner mass. A finite outer boundary has a calculable correction. Equal ratios are also possible for another mechanism producing the same field. [Prospective annular prediction and projection controls](receipts/lensing_profile.json) · [Disjoint rotation test](receipts/held_galaxies.json)
4. **Test source attribution across scales.** Freeze the repair law, preparation and readout, conversions, nuisance distributions, intrinsic scatter and rejection rule before opening a fresh observational sample. Require the same source to predict held rotation curves and shear, merger lensing relative to galaxies and hot plasma, and the gravitational effects on CMB spectra and structure growth. Compare predictions with a collisionless-matter baseline on identical data, nuisance information and uncertainty budgets. A successful test of the field does not identify its carrier if both models give the same observer records. Source attribution needs an observable on which their predictions differ. [Cluster shear and gas measurements](https://arxiv.org/abs/astro-ph/0608407) · [Planck temperature, polarization and lensing](https://arxiv.org/abs/1807.06209)

For a single spherical excess, write $v_A^2=v_{\rm obs}^2-v_{\rm bar}^2$ and $D(R)=\Delta\Sigma_{\rm obs}(R)-\Delta\Sigma_{\rm bar}(R)$. At $0<R_0<R$ the matched observable is

$$
\frac{4G\,[R^2D(R)-R_0^2D(R_0)]}
     {R v_A^2(R)-R_0 v_A^2(R_0)}=1.
$$

When the denominator is poorly resolved, use the difference of the numerator and denominator with their joint covariance. For compact spherical baryons the annular differences also cancel their inner mass. Extended gas or disks require their own orbital and projected-light models. With a sharp outer cutoff $L$, radii $R_0=L/4$, $R=L/2$ give **1.014799** under the same zero-slip hypothesis. The released stacked acceleration bins and the disjoint SPARC catalogue measure different lenses. The matched ratio requires its own observation. [Numerical projection and cutoff checks](lensing_profile.py) · [Annular lensing statistic](https://arxiv.org/abs/0911.4973)

## Predictions for held galaxies

Five folds defined by galaxy name hold out every radius of a galaxy together.
On the historical deep cut, each training fold fits one acceleration scale.
Across 115 held galaxies (960 radii), the conditional deep law predicts total
radial acceleration with an RMS log error of **0.15055 dex**. Each galaxy has
equal weight. The empirical radial acceleration relation gives 0.15100 dex,
a fitted constant multiplier of baryonic gravity gives 0.17525 dex, and
baryons alone give 0.69277 dex. The conditional law improves on the fitted
multiplier by 0.02471 dex (paired galaxy bootstrap range 0.01083–0.03867 dex).
Its difference from the empirical relation is −0.00046 dex (−0.00203 to
+0.00117 dex). Both intervals condition on the five training fits. Removing
the fractional velocity-error cut admits eight nonpositive excesses and gives
RMS errors of 0.19301 and 0.19434 dex for the conditional and empirical laws.

At the 97 matched outer rings, the ratio of the flat-speed and **single-ring**
radial estimates of `a_0` equals `r v_bar²/(G M_b)`. Raising each measured
outer speed by 20% after selection moves the flat-speed estimate by a median
of 0.382 dex and preserves the ratio to within 2 × 10⁻¹⁵ dex. The published
radial fit also includes interior radii. These retrospective tests establish
transfer of a fitted acceleration law within SPARC. The catalogue provides
rotation and baryonic data. Identifying the modular source requires an
independent observable of collar charge and a calibrated map to acceleration.
[Held-galaxy receipt](receipts/held_galaxies.json) · [SPARC observations](https://astroweb.case.edu/SPARC/)

```bash
python3 code/cosmology/rar_deep_regime/held_galaxies.py
```

### Can a single exterior profile fit individual rotation curves?

The proposed density gives a constant anomalous contribution to squared
orbital speed outside a compact baryonic source. We test this consequence
using total photometric mass and the catalogue's radial Newtonian force
from gas and stars. We select 151 deep-regime rings in 41 galaxies with no
catalogued bulge contribution, beyond both four stellar disk scale lengths
and the H I isodensity radius. These cuts are independent of measured speed.
Gas can extend beyond that radius. Applying the compact-source asymptote
here is a physical hypothesis. [SPARC tables and radius conventions](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/AJ/152/157?format=html&tex=true)

| Distance and inclination allowance | Disk mass-to-light ratio | Lower bound on largest speed residual, in quoted errors |
| --- | --- | ---: |
| Catalogue values | 0.5 | 33.01 |
| Within one quoted error | 0.4–0.6 | 6.00 |
| Within two quoted errors | 0.4–0.6 | 3.32 |

At the simultaneous 95% velocity threshold of 3.590 quoted errors, the
one-error geometry box requires $a_0 \ge 1.373\times10^{-10}$ from NGC 0024
and $a_0 \le 0.844\times10^{-10}\,\mathrm{m\,s^{-2}}$ from NGC 0801. In the
two-error box, $a_0=1.146\times10^{-10}\,\mathrm{m\,s^{-2}}$ and one shared
distance, inclination and disk mass-to-light ratio per galaxy reproduce all
151 speeds within 3.590 quoted errors. The archive records those 41
assignments and reprojects every speed. This is a feasibility result with
a scale chosen from the data. The velocity threshold assumes marginal
Gaussian errors and conditions on the chosen nuisance ranges. It gives no
joint probability for the assigned distances, inclinations or stellar masses.
[Galaxy assignments and bounds](receipts/held_galaxies.json)

A second test fits the deep radial-acceleration scale on the **74 other**
qualifying SPARC galaxies (584 radii) and transfers the result to all 41
exterior galaxies. The fitted scale is $7.60\times10^{-11}\,\mathrm{m\,s^{-2}}$.
Resampling the 74 training galaxies gives a 95% range of $(6.19,9.22)\times10^{-11}\,\mathrm{m\,s^{-2}}$ at catalogue inputs.
With the fitted scale fixed, the largest exterior speed residual is at least
6.25 quoted errors even in the two-error geometry box. This discrepancy
tests transfer between the deep radial formula and its compact photometric
exterior approximation on disjoint galaxies. [Fixed-scale transfer test](receipts/held_galaxies.json)

### Is the outer discrepancy caused by replacing the disk force with a compact mass?

Two predictions use the same 151 rings and the same acceleration scale fitted on
74 other galaxies. One adds $\sqrt{G M_{\rm phot}a_0}$ to the baryonic squared speed.
The other adds $\sqrt{R V_{\rm bar}^2 a_0}$ using the tabulated gas and
disk force at each ring. Both terms scale linearly with distance at fixed
angular size and flux. The archive lists both nominal speeds for every ring.
At NGC 24's outermost ring (11.27 kpc), the compact and local-force values are
81.11 and 83.12 km/s. SPARC records 110.00 ± 2.83 km/s. Taking the most
favorable distance, inclination and mass-to-light ratio inside the two-error
catalogue boxes raises the two predictions to 92.30 and 94.47 km/s.
The respective residuals are at least 6.254 and 5.488 quoted speed errors.
For the local-force prescription only that ring is incompatible with the
3.590-error simultaneous comparison. A single parameter set per galaxy
covers all other 150 rings. This diagnostic keeps the incompatible ring in
the primary score. Even the upper end of the training-scale bootstrap
interval gives at most 98.29 km/s in the original geometry boxes. Repeating
the training calibration using only galaxies with zero catalogued bulge contribution or
removing the training velocity-error cut lowers its scale.
[Predictions, training variants and shared galaxy assignments](receipts/held_galaxies.json)

The underlying H I publication gives 109.7 ± 2.0 km/s at 320 arcseconds and
a kinematic inclination of 64° ± 3°. Its neutral-gas and optical isophotes
give 76° and 78°. The same galaxy's inner Hα velocity field gives 75° ± 3°.
An alternative stellar tip distance compilation gives 7.67 ± 0.32 Mpc.
At the central Hα inclination, central alternative distance and disk
mass-to-light ratio 0.5, the local-force prediction is 91.57 km/s in the
SPARC speed frame. At 8.31 Mpc, 78° and 0.6, the three enrolled ring speeds
are 105.02, 102.37 and 100.10 km/s, with a maximum residual of 3.499 quoted
errors. This point uses the observed discrepancy to choose galaxy inputs.
The H I publication and SPARC use the same outer observations. The stellar
distance is independent of those speeds; the distance compilations can
share underlying imaging. [Chemin et al., §5.1 and Table 7](https://arxiv.org/pdf/astro-ph/0609148)
· [Dicaire et al., Table 4 and Appendix A](https://arxiv.org/pdf/0801.2971)
· [Schombert et al., §2.3 and Table 2](https://doi.org/10.3847/1538-3881/ab9d88)

A separate sensitivity check applies the same catalogue rule to all 41
galaxies. Six have alternative TRGB distances in Schombert et al. Table 2
and two have Hα inclinations in Dicaire et al. Table 4. Angles are restricted
to the physically admissible range. All other geometry inputs keep their
SPARC values. Changing only the distance boxes leaves one ring at least
4.890 quoted errors from the held prediction. Changing only the inclination
boxes leaves one at least 3.783 errors away. Combining both sources admits
one assignment per galaxy for **all 151 rings**, with a maximum residual of
3.5834 quoted errors. Those assignments consume the held speeds. The result
quantifies sensitivity to projection and baryon geometry at a fixed scale.
An independent prediction requires distances, inclination as a function of
radius, and deprojected gas and stellar maps fixed before reading the outer
speeds. These force substitutions supply no map from the observer's repair
history to the gravitational source. [Complete geometry controls](receipts/held_galaxies.json)

For comparison, at the catalogue geometry a freely varying, nonnegative and
nondecreasing spherical enclosed mass per galaxy accommodates all 151 speeds
within 0.366 quoted errors. Including bulges extends that result to 237
rings in 58 galaxies. The rotation data constrain the shape and normalization
of the excess. Attribution to a modular source requires an independent
charge observable linked to the same radii.

## Source-bound full-RAR calibration

`sparc_full_rar_calibration.py` is deliberately separate from the conditional
OPH deep law above. It authenticates the committed CDS/VizieR snapshot and
fits the published empirical interpolation
`g_obs = g_bar / (1 - exp(-sqrt(g_bar/a_0)))`. The receipt reproduces the
153-galaxy parent selection, retains 2696 strict-cut points in 147 galaxies,
and records the remaining three-point difference from the published census.
It also records fixed mass-to-light and simple weighting sensitivity.

Status: `retrospective_same_data_calibration_diagnostic`,
`oph_specific_empirical_evidence: false`,
`independent_held_out_data: false`. The interpolation, acceleration scale, and
SPARC data predate this replay, `a_0` is fitted, and the nuisance likelihood is
incomplete. This producer strengthens source custody and reproducibility; it
does not test an OPH-derived value.

```bash
python3 code/cosmology/rar_deep_regime/sparc_full_rar_calibration.py --check
python3 -m pytest -q \
  code/cosmology/rar_deep_regime/test_sparc_full_rar_calibration.py
```

## Quadrupole cross-check (uncertified diagnostic)

`qumond_quadrupole_crosscheck.py` computes the Solar-System external-field
quadrupole `Q2` of a local nonlinear field law from first principles (inner
multipole expansion of the phantom potential) and compares it with Blanchet
and Novak (2011) for the simple function (2 percent) and with the Park et al.
(2026) benchmark for the radial-acceleration function (`8e-5`). The central
values are not a local numerical certificate: the integrations emit thousands
of warnings and lack singularity subtraction, radial-tail bounds, and an
independent high-precision replay. It also
records the footprint of the density formulation used in the paper: zero
quadrupole for an isotropic ambient density, tidal scale `1.9e-31 s^-2`, and
`5e-15` solar masses inside ten astronomical units. The historical field-law
assessment and the current certification boundary are documented in
`DARK_SECTOR_FIELD_LAW_POSTMORTEM.md`.

```bash
python3 code/cosmology/rar_deep_regime/qumond_quadrupole_crosscheck.py \
  --output code/cosmology/rar_deep_regime/receipts/qumond_quadrupole_crosscheck.json
```

## Where can the lensing peaks occur in a cluster merger?

The western Bullet subcluster tests every strict maximum allowed by a total convergence made of two circular, radially nonincreasing profiles centred on the X-ray plasma. Moving toward the line segment between the centres shortens the distance to both. The western peak lies 36.61 arcsec outside that segment, about 162 kpc at the authors' plate scale. Their bootstrap offset significance refers to the gas centre. The later JWST analysis reports about 150 kpc for the western peak. A dynamical source should predict its sky position before comparison with the lensing map. Replay with `python bullet_cluster.py --check` from this directory. [Calculation](bullet_cluster.py) · [Receipt](receipts/bullet_cluster.json) · [Clowe et al., §3 and Table 2](https://doi.org/10.1086/508162) · [Cha et al., §IV.2](https://arxiv.org/abs/2503.21870)

## A rotation-calibrated prediction for an independent shear measurement

Brouwer et al. release tangential and cross shear as excess surface densities in the same baryonic-acceleration bins, with multiplicative corrections and a full covariance for each stack. We use seven bin centres in the reliable-isolation, deep regime. For compact spherical baryons and negligible gravitational slip, projecting the paper’s proposed density gives $4G\Delta\Sigma=(4/\pi)g_b+\sqrt{a_0g_b}$. The scale comes from the existing deep rotation fit to 115 SPARC galaxies. [Public data and acknowledgments](https://kids.strw.leidenuniv.nl/sciencedata.php) · [Brouwer et al., Fig. 4](https://doi.org/10.1051/0004-6361/202040108)

| Lens sample | Fixed projected profile χ² | Compact baryons χ² | Cross shear against zero χ² |
| --- | ---: | ---: | ---: |
| KiDS-bright | 71.10 | 739.66 | 4.32 |
| GAMA | 16.66 | 58.02 | 6.43 |
| KiDS-bright, hot-gas mass alternative | 43.00 | 764.71 | 5.77 |

The primary projection underpredicts six of the seven KiDS-bright tangential bins. The lowest bin straddles the isolation floor. Restriction to the six bins wholly above it gives χ² = 40.67 for the fixed KiDS-bright profile. The statistic conditions on the published mass conventions, the bin-centre approximation and the supplied statistical covariance. The hot-gas column changes the baryonic binning while retaining the compact-source projection. An extended gas distribution needs its own lensing calculation. The stacks share KiDS source data. Each result stands separately. Curvature sourced by a measured modular charge, and the relativistic relation between the potential measured by orbits and the one measured by light, remain open. [Calibrated measurements, covariance and predictions](receipts/lensing_profile.json)

**Mass normalization and finite bins.** Brouwer et al. report a ±0.2 dex stellar-mass systematic. With published memberships fixed, we scan that range as a common shift of the entire compact baryonic mass. For the six KiDS-bright bins the lowest χ² is 27.01 at the bin centres and 27.10 with uniform lens–source pair density per unit sky area. The endpoints of the independent SPARC scale interval give 56.23 and 29.11 at nominal masses and bin centres. If each bin can have its own unrestricted nonnegative pair weights, a common +0.135 dex shift permits χ² = 4.93. That construction puts all weight at the upper edges of the first two bins and the lower edge of the last. The published summaries give aggregate weights. Establishing the actual profile requires pair-level sampling and stellar/gas fractions. [Survey mass systematics, §3.3 and Fig. 4](https://doi.org/10.1051/0004-6361/202040108) · [Numerical controls](receipts/lensing_profile.json)

**A sign test independent of amplitude.** Divide each bin's excess surface density by the square root of its central baryonic acceleration in units of $10^{-12}\,\mathrm{m\,s^{-2}}$. For the same fractional sampling distribution in each log bin, the isothermal contribution becomes constant and the positive compact-baryon contribution rises with acceleration. The lower-three mean minus the upper-three mean must therefore be nonpositive for any positive constant amplitudes. KiDS-bright gives $+4.25\pm1.75\,h_{70} M_\odot\,\mathrm{pc}^{-2}$ from the supplied covariance. GAMA gives $+15.67\pm7.63$ in the same units. The retrospective contrast tests the common shape. Distinguishing source pressure, gravitational slip, mass-model errors and sampling requires independent measurements. A direct test would measure orbital tracers and shear around the same lenses at matched radii, applying the model to each weighted pair. Strong-lens and stellar-dynamical measurements already constrain slip in the inner region of ESO 325-G004. Applying that constraint to these outer haloes would require a radial model. [Collett et al.](https://arxiv.org/abs/1806.08300)

To replay from the verified 2.3 MB public archive, run from the RER root.

```sh
curl -fsSLo /tmp/brouwer2021_rar.tar https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar
python code/cosmology/rar_deep_regime/lensing_profile.py /tmp/brouwer2021_rar.tar
```

### Shear from a specific record-formation law

For stationary, spatially uniform mark transport, the sparse sign-aligned record density implies that excess surface density is concave in inverse squared radius on projected annuli free of mark injectors. This assumes ordinary thin-lens coupling and gives a shape test at any positive normalization. The largest departures in the forbidden direction across 20 measured chords were 0.40 standard errors for KiDS-bright and 1.46 for GAMA. The one-sided threshold was 2.81 for each sample. The published stacks do not independently establish the locations of mark injectors or the common lens mixture across bins. [Derivation, finite graph and survey replay](record_shear/README.md)

### Stellar motion and gravitational lensing in the same galaxies

For 53 single early-type SLACS galaxies with well-measured spectra and modeled arcs, the stars supply a prediction of the Einstein angle. Assume a spherical total field with constant circular speed, equal weak-field metric potentials and an isotropic de Vaucouleurs stellar distribution. The Jeans equation then relates the circular speed to the dispersion measured through the seeing-convolved SDSS fibre. At observed Einstein radii of 1.14–6.92 kpc the median predicted/measured angle is 1.101. The median absolute fractional error is 16.1%. Comparing in the measured dispersion coordinate gives 39 of 53 within two assigned errors and χ² = 163.2. Trying the survey's alternate mass model for the same arc images gives a 16.8% median absolute angle error. Galaxy structure, orbital anisotropy and selection contribute additional uncertainty. [Bolton et al., Tables 3–5](https://doi.org/10.1086/589327) · [Per-galaxy predictions](receipts/stellar_lensing.json)

For all 22 galaxies in the SLACS weak-lensing stack, published stellar-population masses supply the compact baryonic sources. Stellar photometry and background shear can share imaging systematics. The SPARC acceleration scale sets the isothermal excess. Requiring each lens to be deep and at least ten stellar effective radii out leaves two bin centres. The control extends each galaxy's spectroscopically inferred circular speed flat to the outer measurement. Predictions use complete-annulus pair weights with uniform source density and shape variance. Surface densities in the receipt are in $h M_\odot\,\mathrm{pc}^{-2}$, with $h=0.7$. [Gavazzi et al., Eq. 10 and Tables 1–2](https://doi.org/10.1086/519237) · [Auger et al., Table 4](https://doi.org/10.1088/0004-637X/705/2/1099) · [Matched stellar inputs](data/slacs_weak_lenses.dat)

The complete bin comparisons, including both stellar-mass prescriptions, the spectroscopic control and measured cross shear, are archived with their units and conditional uncertainties. [Per-bin observations and forecasts](receipts/stellar_lensing.json)

The nominal stellar-source predictions fall below both tangential measurements. The same two bins give χ² = 13.36 for the cross-shear zero control. The actual pair weights and extended gas are unmeasured here. At the outer bin, take each stellar mass up to its Salpeter estimate plus one quoted error and the SPARC scale up to its upper bootstrap endpoint. Allow arbitrary positive weights over galaxies and arbitrary radii within the ordered radial bin. If the weighted mean source efficiency for each lens in that bin equals its catalogue value, the conservative upper bound is 58.02, 0.12 quoted errors below the measured central value. Allowing arbitrary source efficiencies between zero and one raises the upper bound to 96.42. A tighter comparison needs the weighted pair radii, an independent outer shear control, and gas and environmental mass constraints. Both data sets compare different tracers of the same total metric. Resolving the proposed annular ratio additionally requires orbital measurements at the shear radii. [Calculation](stellar_lensing.py) · [Measurements and numerical bounds](receipts/stellar_lensing.json)

The 35 inner lenses outside the weak-lensing stack give a median observed-to-predicted Einstein angle of 0.855. Applying this factor to the 18 held inner lenses raises their median absolute angle error from 11.6% to 15.8%. Applied to the stellar speeds extended flat in radius, it gives outer surface-density contrasts of 44.0 and 25.2 against measurements of $52\pm17$ and $60\pm17$ in $h\,M_\odot\,\mathrm{pc}^{-2}$. The conditional two-bin χ² rises from 3.23 to 4.41. The split uses galaxy identities. Sample selection, radial profile shape and outer shear uncertainties affect its interpretation. [Disjoint normalization test](receipts/stellar_lensing.json)

A published satellite-velocity and weak-shear comparison uses common central galaxies, but selects satellites with radius and velocity cuts derived from lensing fits. Its radial velocity distributions therefore cannot serve directly as an independent orbital forecast for this annular test. A new reduction would select hosts and satellites without using shear, model orbital anisotropy and interlopers, and carry calibrated source weights into the same host and radius bins. [Ma et al., §3.2.1](https://arxiv.org/abs/2409.13329)

Replay from the RER root using the public SLACS catalogue. The program verifies the scientific columns against pinned digests.

```sh
curl -fsSLo /tmp/bolton_slacs.tsv 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J%2FApJ%2F682%2F964&-out.all&-out.max=unlimited'
OPENBLAS_NUM_THREADS=1 python code/cosmology/rar_deep_regime/stellar_lensing.py /tmp/bolton_slacs.tsv
```

A static spherical branch with vanishing radial pressure gives exact circular, clock and bending measurements after an outer edge is specified. The public rotation and shear records also test the larger class of shared normalized profiles with sublinear enclosed-mass growth. At generous mass, inclination and bin-weight settings, the outer speed of AGC 114905 bounds the weighted excess surface densities at 11.36 and 11.18 M☉/pc² against measured KiDS-bright and GAMA values of 17.30 ± 0.69 and 22.73 ± 2.99. The physical identification of the measured dwarf radius with this branch is a hypothesis. [Relativistic source and cross-data bounds](finite_halo/README.md) · [External gas-rich galaxy sample](udg_transfer/README.md)

Nine early-type galaxies have measured hydrogen orbits enclosing at least 90% of their modelled baryons beyond five effective radii. At nominal stellar masses and distances NGC 4262 and NGC 3626 give incompatible excess masses. At lower acceleration, a common linear profile fits three outer orbits from the 80% selection and both survey mean shears within the stated intervals. The six-bin KiDS-bright shear profile misses at nominal bin centres even after fitting its scale. Unknown pair weights within those bins materially weaken this test. [Outer-orbit comparison and uncertainties](outer_orbit_transfer/README.md)

M31 supplies a distributed source from photometry and 100 resolved hydrogen rings. Its coherent stationary mark field can increase the enclosed dark mass by at most a factor of 1.443 from 28.17 to 38.09 kpc for any constant positive gain and any lifetime. The circular speeds inferred by Chemin et al. require a factor of 2.727, or 2.126 when both speeds shift by two quoted errors in the favorable directions. A separate survey by Corbelli et al. permits the model within its wider outer velocity errors. [M31 orbits, source assumptions and all-lifetime bound](m31_transport/README.md)

A closed redistribution with unchanged gravitational mass per unit energy must debit the surroundings for any net mass inside a galaxy aperture. At 291 kpc, four isolated-galaxy lensing stacks require weighted mean exterior deficits of 0.40–3.42 trillion solar masses even after two favorable statistical errors, the deprojection allowance and twice the maximum ordinary mass in each bin. A donor at the cosmic baryon density at redshift 0.5 would need to extend at least 1.66–3.39 Mpc. [Compensation, direct shear and finite-shell control](compensation/README.md)

For ten nearby galaxies, the perimeter of the measured baryonic surface-density field gives a median held outer-speed error of 35 km/s. Independent stellar distances for a preselected pair give an excess-speed-squared ratio of 1.29. Its measured optical boundaries predict 5.75. Two-dimensional maps of stars, H I and molecular gas predict 4.09, while subtracting the ordinary force of those same mapped components gives an observed ratio of 1.54. The two gas phases alone predict 1.56 on this pair. Calibrated on NGC 3198, their full maps miss three of the other eight outer speeds by 76–99 km/s. [Resolved boundary source and limits](boundary_source/README.md)
