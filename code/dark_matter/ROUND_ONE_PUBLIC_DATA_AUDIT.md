# Round-One OPH Public-Data Audit

This note records the first audited no-new-observer-simulation comparisons for
the OPH cosmology continuation. It supplements the older generated empirical
scorecard. It does not alter that scorecard's invalidated compressed-cosmology
status and does not define a combined OPH score.

## Results

| Surface | Audited result | Status |
| --- | ---: | --- |
| Fixed Z6 RAR | `0.1328335133 dex` RMS on 2,693 rows; same-data diagonal-chi-square optimum has `0.1329467473 dex` RMS | Positive retrospective fixed-formula check; conditional/calibrated coefficient. |
| Unit-collar RAR | `0.1342081979 dex` RMS | Conditional Jensen-band endpoint. |
| BTFR slope | error-aware `3.8456543 +/- 0.0858172`; OPH slope 4 is `1.79854 sigma` high | Compatible and inside the published systematic range 3.5–4.0. |
| Z6 BTFR normalization | data-minus-model `-0.1347626 dex` at 100 km/s; `-6.4707 sigma` statistical-only | Pressure point; global galaxy nuisances are not marginalized. |
| Cassini Z6 quadrupole | `3.6201778e-26 s^-2`; `19.2232 sigma` fixed-input residual | Universal/full-source static extension excluded; settled-galaxy applicability theorem missing. |
| Cassini unit endpoint | `3.4021876e-26 s^-2`; `18.0122 sigma` fixed-input residual | The current Jensen collar band does not rescue universal applicability. |
| Analytic scalar tilt | `n_s=1-P/48=0.9660214956`; `+0.267 sigma` from Planck summary | Positive conditional arithmetic check. |
| Analytic tilt CAMB/Planck TT | `chi2/bin=0.9544982`, baseline `0.9444957`, total diagonal `Delta chi2=+0.83020` over 83 bins | Positive conditional transfer diagnostic; not official likelihood. |

The previous unweighted BTFR slope `3.4863` and `0.9469 dex` intercept
mismatch are superseded. They discarded uncertainties in both coordinates and
compared unequal-slope intercepts at 1 km/s, far outside the data.

## Reproducibility Surface

The maintained implementations live in the sibling
[`oph-physics-sim`](https://github.com/muellerberndt/oph-physics-sim)
repository:

- `oph_fpe/cosmology/btfr_likelihood.py`
- `oph_fpe/cosmology/rar_fixed_comparison.py`
- `oph_fpe/cosmology/cassini_external_field.py`
- `oph_fpe/cosmology/public_data_comparisons.py`
- `tools/best_of_public_data.py`
- `tools/analytic_p48_planck_check.py`
- `tests/test_public_data_comparisons.py`
- `schemas/cosmology/best_of_public_data_comparisons.schema.json`

From that repository:

```bash
python3 -m pytest -q tests/test_public_data_comparisons.py

python3 tools/best_of_public_data.py \
  --primary-run runs/oph_universe_64k_final_audited_20260711 \
  --baseline-run runs/oph_universe_64k_final_audited_20260711 \
  --out runs/best_of_public_data/current \
  --strict

python3 tools/analytic_p48_planck_check.py \
  --out runs/best_of_public_data/analytic_p48_planck
```

The comparison bundle binds the public table hashes, exact source modules,
selected run directory, and result sidecars. There are currently zero frozen
physical-prediction receipts.

## Required Next Calculations

1. Solve the conservative disk PDE rather than only the algebraic RAR proxy,
   then test curl-field residuals against disk geometry and environment.
2. Run a hierarchical joint RAR/BTFR likelihood with distance, inclination,
   stellar mass-to-light, gas, intrinsic-scatter, and covariance nuisances.
3. Forward-project the conditional no-slip branch into KiDS ESD vectors and
   their full covariance.
4. Execute the analytic P/48 branch with official Planck and ACT likelihoods.
5. Derive, rather than fit, a source applicability/coarse-graining rule that
   preserves the galaxy branch while satisfying Cassini.

## Claim Boundary

- RAR and BTFR are conditional continuation checks, not recovered-core
  predictions.
- The BTFR table is a separate public observable drawn from an overlapping
  SPARC sample, not a wholly independent galaxy population.
- Cassini excludes universal application of the displayed static equation. It
  does not by itself exclude OPH's explicitly old-settled-galaxy scope.
- The favorable analytic tilt and unfavorable finite-clock tilt are distinct
  branches and may not be substituted after inspecting residuals.
- The archived compressed-cosmology row remains invalidated because it is
  hard-coded, lacks an attached covariance, and inherits a rejected neutrino
  input.

