# Conditional Candidate Outputs

Generated: `2026-07-12T11:21:03Z`

This surface renders the numeric conditional candidates that the strict
source-only promotion policy withholds from the public OPH-value columns.
Each row carries its recorded claim label and its open gates. No row is
promotable, and this surface must not feed back into promotion decisions.
The public claim boundary stays in `RESULTS_STATUS.md`; target-anchored
exact fits stay in `EXACT_FITS_ONLY.md`.

The electroweak row is the target-free D10 value law. It is the OPH
derivation for the pair: no measured mass enters, and the value follows as
an exact implication of the five quotient-transport premises. The
diagnostic base-running chart and the target-anchored freeze-once adapter
stay in the technical audit README and off this surface.

## Conditional GeV Chart Candidates

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `m_H` | `125.1995304097179 GeV` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, qt1_qt5_quotient_path_source_entailment, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |
| `m_t` | `172.3523553288312 GeV` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, qt1_qt5_quotient_path_source_entailment, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion, top_threshold_control` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |
| `M_W` | `80.37700001539531 GeV` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `M_Z` | `91.18797807794321 GeV` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |

## Conditional Scale-Free Chart Ratios

These ratios cancel the chart scale. The value-law triple shares one
emitted quintet, so the three coordinates stand or fall together. The
Higgs/top ratio shares the declared-surface split pair.

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `M_W / M_Z` | `0.8814429457652062` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `sin2_theta_W_eff` | `0.22305833336075578` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `alpha_em_eff_inv` | `132.7524088191668` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `m_H / m_t` | `0.726416126839982` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |

## Conditional Dimensionless Candidates

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `alpha_inv_no_hadron_near_endpoint` | `137.035959500817279952949994585787193262029` | `deterministic_pixel_root_plus_unified_width_term` | `source_spectral_endpoint, same_scheme_hadronic_remainder, full_interval_source_certificate` | `code/particles/runs/status/final_end_to_end_predictions.json` |
| `v_over_E_star` | `2.0199803239725553e-17` | `closed_form_hierarchy_law_on_public_endpoint_packet` | `e_star_physical_normalization, clock_scale_certificate, frozen_running_threshold_convention` | `code/particles/hierarchy/computations/hierarchy_numeric_witness.json` |

Supersession note (2026-07-14): the recorded
`alpha_inv_no_hadron_near_endpoint` value is a mixed-provenance display packet
(the inner value evaluated at `P_fwd` plus `alpha_U` evaluated at the SL-3
pixel) from a run that predates the converged precision-100 reruns; it is not
a fixed point of any single declared map and stays on record as a display
packet only. The certified self-consistent gauge-width fixed point is
`alpha^-1 = 137.035660136946577...`, and the certified source root is
`alpha^-1 = 136.994835177413...` (CL-2/CL-6 in `../../CLOSURE_LEDGER.md`).
