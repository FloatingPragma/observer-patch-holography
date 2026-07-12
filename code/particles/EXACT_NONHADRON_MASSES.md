# Public Non-Hadron Mass Outputs

Generated: `2026-07-12T03:22:10Z`

This bundle gives numeric non-hadron mass outputs that are not target-anchored witness rows.
Target-anchored exact fits and compare-only absolute attachments are withheld from this prediction table and kept in `EXACT_FITS_ONLY.md` for audit/debug use.
Charged-lepton and quark exact same-family/current-family surfaces remain available as exact-fit audit witnesses, but their numeric values are not public mass outputs under the strict non-circularity policy.
Absolute neutrino masses are likewise withheld here while the absolute attachment remains compare-only; dimensionless PMNS and mass-splitting-ratio comparisons stay in `RESULTS_STATUS.md`.
Photon, gluon, and graviton labels are not emitted as `0 GeV` particle predictions. Their zero hard quadratic parameters are tracked separately as conditional classical/perturbative carrier modes in `CARRIER_MODE_ACCEPTANCE.md`; the quantum-particle gate remains open.

| Particle | Mass Output | Kind | Scope | Source |
| --- | ---: | --- | --- | --- |
| Higgs Boson | `125.1995304097179 GeV` | `conditional_declared_surface_higgs_top_candidate` | `declared_d10_d11_running_matching_threshold_surface_only` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |

## Separated Classical Carrier Modes

| Carrier | Hard quadratic mass parameter squared | Classical mode gate | Quantum particle gate |
| --- | ---: | --- | --- |
| Electromagnetic carrier | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` |
| Color gauge carrier | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` |
| Einstein tensor carrier | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` |
