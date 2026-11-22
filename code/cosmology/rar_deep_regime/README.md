# Deep-regime comparison for the anomalous-density profile

`rar_deep_regime_diagnostic.py` compares the dark-matter paper's conditional
deep-regime laws, `a_A = sqrt(a_b a_0)` and `v^4 = G M_b a_0`, with the SPARC
rotation-curve sample (snapshot under `data/`, see `data/PROVENANCE.md`).

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
