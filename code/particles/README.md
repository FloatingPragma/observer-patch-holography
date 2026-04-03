# OPH Particle Derivation Code

This directory is the canonical particle-spectrum code path for OPH inside
`reverse-engineering-reality`.

## Scope

The goal is to keep one auditable derivation surface from the OPH inputs
already derived in the papers to the emitted particle-spectrum artifacts used by
the particle paper:

- electroweak calibration
- Higgs/top
- charged leptons
- quarks
- neutrinos
- hadrons
- public status rendering

Historical Oracle batches, worker logs, and transient handoff material are not
part of this canonical tree.

## Active Layout

- [calibration](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/calibration)
- [flavor](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/flavor)
- [leptons](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/leptons)
- [neutrino](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/neutrino)
- [hadron](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron)
- [qcd](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/qcd)
- [runs](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs)
- [scripts](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/scripts)
- [RESULTS_STATUS.md](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/RESULTS_STATUS.md)
- [particle_mass_derivation_graph.svg](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/particle_mass_derivation_graph.svg)
- [task_tracker.yaml](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/task_tracker.yaml)

## Current High-Level Chain

- electroweak calibration:
  `calibration/derive_d10_ew_observable_family.py ->
  calibration/derive_d10_ew_source_transport_pair.py ->
  calibration/derive_d10_ew_population_evaluator.py ->
  calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py ->
  calibration/derive_d10_ew_source_transport_readout.py`
- Higgs/top:
  `calibration/derive_d11_forward_seed.py ->
  calibration/derive_d11_forward_seed_promotion_certificate.py`
- charged leptons:
  support and scale artifacts feeding
  `leptons/derive_lepton_excitation_gap_map.py ->
  leptons/derive_lepton_log_spectrum_readout.py ->
  leptons/build_forward_charged_leptons.py`
- quarks:
  `flavor/derive_quark_sector_mean_split.py ->
  flavor/derive_quark_sector_descent.py ->
  flavor/build_forward_yukawas.py`
- neutrinos:
  `neutrino/derive_neutrino_scale_anchor.py ->
  neutrino/derive_family_response_tensor.py ->
  neutrino/derive_majorana_holonomy_lift.py ->
  neutrino/derive_majorana_phase_pullback_metric.py ->
  neutrino/build_forward_majorana_matrix.py ->
  neutrino/build_forward_splittings.py ->
  neutrino/derive_neutrino_weighted_cycle_repair.py ->
  neutrino/derive_neutrino_bridge_rigidity_theorem.py ->
  neutrino/derive_neutrino_absolute_attachment_theorem.py ->
  neutrino/export_forward_neutrino_closure_bundle.py`
- hadrons:
  `qcd/derive_lambda_msbar_descendant.py ->
  hadron/derive_full_unquenched_correlator.py ->
  hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py ->
  hadron/derive_stable_channel_cfg_source_measure_payload.py ->
  hadron/derive_stable_channel_sequence_evaluation.py ->
  hadron/derive_stable_channel_groundstate_readout.py`
- rendered public surface:
  `scripts/build_results_status_table.py`

## Main Outputs

- status table:
  [RESULTS_STATUS.md](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/RESULTS_STATUS.md)
- exact-fits-only diagnostic surface:
  [EXACT_FITS_ONLY.md](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/EXACT_FITS_ONLY.md)
- exact non-hadron mass bundle:
  [EXACT_NONHADRON_MASSES.md](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/EXACT_NONHADRON_MASSES.md)
- machine-readable status:
  [results_status.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/results_status.json)
- machine-readable exact-fits-only surface:
  [exact_fits_only.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/exact_fits_only.json)
- machine-readable exact non-hadron mass bundle:
  [exact_nonhadron_masses.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/exact_nonhadron_masses.json)
- frozen status artifact:
  [status_table_forward_current.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/status/status_table_forward_current.json)
- frozen exact-fits-only artifact:
  [exact_fits_only_current.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/status/exact_fits_only_current.json)
- frozen exact non-hadron mass bundle:
  [exact_nonhadron_masses_current.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/status/exact_nonhadron_masses_current.json)
- derivation graph:
  [particle_mass_derivation_graph.svg](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/particle_mass_derivation_graph.svg)

## Current Quark Boundary

- maximal theorem-emitted package artifact:
  [quark_maximal_theorem_emitted_package.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/flavor/quark_maximal_theorem_emitted_package.json)
- present-premise no-go artifact:
  [quark_current_premise_no_go_theorem.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/flavor/quark_current_premise_no_go_theorem.json)
- minimal-extension closure artifact:
  [quark_minimal_extension_closure_theorem.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/flavor/quark_minimal_extension_closure_theorem.json)

These artifacts fix the current quark theorem boundary. On the present premise
set, the maximal theorem-emitted package consists of the D12 mass ray, the
negative selector `sigma_ref`, and the restricted-scope affine mean package
with `g_ch = 0.9231656602589082` on `shared_budget_only` and
`(g_u, g_d) = (0.7797392875757557, 0.12172551081512113)` on
`current_family_only`. The stronger physical closure objects
`Theta_ud^mass`, `Theta_ud^phys`, and `Theta_ud^abs` are absent on that
ledger. The exact minimal extension triple above that package is
`H_mass : ell_ud = log(c_d / c_u)`,
`H_phys : s_ud^phys : M_ud^{CR,phys} -> Sigma_ud^phys`, and
`H_abs : A_q^phys : Sigma_ud^phys -> R`.

## Typical Rebuild

From `reverse-engineering-reality/code/particles`:

```bash
python3 calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py
python3 calibration/derive_d10_ew_source_transport_readout.py
python3 calibration/derive_d10_ew_exactness_audit.py
python3 calibration/derive_d11_reference_exact_adapter.py
python3 neutrino/derive_neutrino_weighted_cycle_repair.py
python3 neutrino/derive_neutrino_bridge_rigidity_theorem.py
python3 neutrino/derive_neutrino_absolute_attachment_theorem.py
python3 neutrino/export_forward_neutrino_closure_bundle.py
python3 neutrino/derive_neutrino_two_parameter_exact_adapter.py
python3 hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py
python3 hadron/derive_stable_channel_sequence_evaluation.py
python3 hadron/derive_current_hadron_lane_audit.py
python3 scripts/build_results_status_table.py
python3 scripts/build_exact_fit_surface.py
python3 scripts/build_exact_nonhadron_mass_bundle.py
python3 scripts/generate_mass_derivation_svg.py
```

## One-Shot CLI Table

For a disposable runtime rebuild that re-runs the active D10/D11/UV builders,
stages the current canonical flavor/lepton/neutrino public-surface artifacts,
and prints the resulting particle status table directly in the terminal:

```bash
python3 compute_current_output_table.py
```

Useful flags:

```bash
python3 compute_current_output_table.py --show-paths
python3 compute_current_output_table.py --with-hadrons --show-paths
python3 compute_current_output_table.py --no-print-table --show-paths
python3 compute_current_output_table.py --verbose
python3 compute_current_output_table.py --format markdown
python3 compute_current_output_table.py --format json
python3 compute_current_output_table.py --color always
```

## Focused Verification

```bash
python3 -m pytest \
  calibration/test_d10_ew_w_anchor_neutral_shear_factorization.py \
  calibration/test_d10_ew_source_transport_readout_artifact.py \
  calibration/test_d10_ew_exactness_audit.py \
  calibration/test_d10_current_carrier_frontier_split.py \
  hadron/test_runtime_schedule_receipt_n_therm_and_n_sep.py \
  hadron/test_stable_channel_sequence_evaluation.py \
  hadron/test_current_hadron_lane_audit.py \
  test_results_status_candidate_policy.py \
  test_results_status_quark_promotion_policy.py \
  test_results_status_structural_rows.py \
  test_predictive_builders_reference_free.py
```

## Paper Surface

The code here feeds the particle paper:

- [deriving_the_particle_zoo_from_observer_consistency.tex](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/paper/deriving_the_particle_zoo_from_observer_consistency.tex)
- [deriving_the_particle_zoo_from_observer_consistency.pdf](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/paper/deriving_the_particle_zoo_from_observer_consistency.pdf)
