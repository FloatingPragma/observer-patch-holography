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
