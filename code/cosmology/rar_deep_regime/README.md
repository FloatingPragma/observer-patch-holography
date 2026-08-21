# Deep-regime comparison for the repair-charge condensate

`rar_deep_regime_diagnostic.py` compares the dark-matter paper's deep-gradient
galaxy law, `a_R = sqrt(a_b a_0)` and `v^4 = G M_b a_0`, with the SPARC
rotation-curve sample (snapshot under `data/`, see `data/PROVENANCE.md`).

The law has one dimensional constant `a_0` that the OPH source does not fix.
The diagnostic therefore reports

* the deep radial-acceleration exponent and the fitted `a_0` on the
  self-consistent deep subset `g_bar < f a_0`, for `f` in `{0.3, 0.1, 0.03}`;
* the baryonic Tully-Fisher exponent and the `a_0` implied by its
  normalisation;
* whether the two constants agree within the combined scatter, which the
  condensate law forces;
* the scatter of the additive all-gradient extension
  `g = g_bar + sqrt(g_bar a_0)` and its Solar-System anomaly, which the
  paper's high-gradient branch is required to remove.

Status: `diagnostic_postdiction_one_fitted_constant`, `physical_claim: false`,
`source_derived_output: false`. The exponents are the parameter-free content;
`a_0` is a fitted comparison value.

```bash
python3 code/cosmology/rar_deep_regime/rar_deep_regime_diagnostic.py \
  --output code/cosmology/rar_deep_regime/receipts/sparc_deep_regime_diagnostic.json
python3 -m pytest code/cosmology/rar_deep_regime
```

## Quadrupole cross-check (documentary)

`qumond_quadrupole_crosscheck.py` computes the Solar-System external-field
quadrupole `Q2` of a local nonlinear field law from first principles (inner
multipole expansion of the phantom potential) and compares it with Blanchet
and Novak (2011) for the simple function (2 percent) and with the Park et al.
(2026) benchmark for the radial-acceleration function (`8e-5`). It also
records the footprint of the density formulation used in the paper: zero
quadrupole for an isotropic ambient density, tidal scale `1.9e-31 s^-2`, and
`5e-15` solar masses inside ten astronomical units. The reason the field-law
formulation was retired is written up in
`docs/DARK_SECTOR_FIELD_LAW_POSTMORTEM.md`.

```bash
python3 code/cosmology/rar_deep_regime/qumond_quadrupole_crosscheck.py \
  --output code/cosmology/rar_deep_regime/receipts/qumond_quadrupole_crosscheck.json
```
