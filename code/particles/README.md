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
  neutrino/build_forward_splittings.py`
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
- machine-readable status:
  [results_status.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/results_status.json)
- frozen status artifact:
  [status_table_forward_current.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/status/status_table_forward_current.json)
- derivation graph:
  [particle_mass_derivation_graph.svg](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/particle_mass_derivation_graph.svg)

## Typical Rebuild

From `reverse-engineering-reality/code/particles`:

```bash
python3 calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py
python3 calibration/derive_d10_ew_source_transport_readout.py
python3 calibration/derive_d10_ew_exactness_audit.py
python3 hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py
python3 hadron/derive_stable_channel_sequence_evaluation.py
python3 hadron/derive_current_hadron_lane_audit.py
python3 scripts/build_results_status_table.py
python3 scripts/generate_mass_derivation_svg.py
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
