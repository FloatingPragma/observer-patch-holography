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

To replay from the verified 2.3 MB public archive, run from the RER root.

```sh
curl -fsSLo /tmp/brouwer2021_rar.tar https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar
python code/cosmology/rar_deep_regime/lensing_profile.py /tmp/brouwer2021_rar.tar
```
