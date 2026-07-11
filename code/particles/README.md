# OPH Particle Derivation Code

This directory is the canonical particle-spectrum code path for OPH inside
`reverse-engineering-reality`.

## Scope

The goal is to keep one auditable derivation surface from the OPH inputs
derived in the papers to the emitted particle-spectrum artifacts used by
the particle paper:

- electroweak calibration
- Higgs/top
- charged leptons
- quarks
- neutrinos
- hadrons
- public claim rendering

Historical Oracle batches, worker logs, and transient handoff material are not
part of this canonical tree.

## Active Layout

- [calibration](calibration)
- [flavor](flavor)
- [leptons](leptons)
- [neutrino](neutrino)
- [hadron](hadron)
- [qcd](qcd)
- [jwst](jwst)
- [fractional](fractional)
- [uhe](uhe)
- [compact_transients](compact_transients)
- [hierarchy](hierarchy)
- [runs](runs)
- [scripts](scripts)
- [RESULTS_STATUS.md](RESULTS_STATUS.md)
- [DERIVATION_GAP_LEDGER.md](DERIVATION_GAP_LEDGER.md)
- [GAP_BUNDLE_CAMPAIGN.md](GAP_BUNDLE_CAMPAIGN.md)
- [campaigns/gap_bundle](campaigns/gap_bundle)
- [particle_mass_derivation_graph.svg](particle_mass_derivation_graph.svg)

## High-Level Chain

- electroweak calibration:
  `calibration/derive_d10_ew_observable_family.py ->
  calibration/derive_d10_ew_source_transport_pair.py ->
  calibration/derive_d10_ew_population_evaluator.py ->
  calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py ->
  calibration/derive_d10_ew_source_transport_readout.py`
- electroweak hierarchy certificate:
  `hierarchy/certificates/DAG_U.json ->
  hierarchy/certificates/R_U_interval_certificate.json ->
  hierarchy/certificates/R_U_krawczyk_certificate.json ->
  hierarchy/computations/hierarchy_numeric_witness.json ->
  hierarchy/certificates/R_EW_global_capacity_certificate.json ->
  hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json ->
  hierarchy/certificates/R_pixel_screen_resonance_summary.json ->
  hierarchy/issue_332_rg_naturality_certificate.json`
- Higgs/top:
  `calibration/derive_bw_higgs_carrier_bridge.py` records the Borel-Weil
  one-Higgs carrier bridge `H_OPH = H^0(CP1,O(1)) ~= C^2` as representation
  and symmetry-breaking geometry only, before the quantitative D11 lane
  `calibration/derive_d11_declared_calibration_surface.py ->
  calibration/derive_d11_forward_seed.py ->
  calibration/derive_d11_forward_seed_promotion_certificate.py ->
  calibration/derive_d11_fixed_ray_no_go_theorem.py ->
  calibration/derive_d11_live_exact_split_pair_theorem.py`
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
  The weighted-cycle output is a target-informed template candidate, not a
  prediction. NuFIT 6.1 rejects its correlated theta23-delta point at the
  declared 3σ gate. `neutrino/score_neutrino_nufit61.py` records the official
  profile result, and `neutrino/audit_neutrino_pmns_conventions.py` finds no
  admissible convention rescue. The bridge and absolute attachment remain
  compare-only.
- hadrons:
  `qcd/derive_lambda_msbar_descendant.py ->
  hadron/derive_full_unquenched_correlator.py ->
  hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py ->
  hadron/derive_stable_channel_cfg_source_measure_payload.py ->
  hadron/derive_stable_channel_sequence_evaluation.py ->
  hadron/derive_stable_channel_groundstate_readout.py`
- JWST compact-object source-release bridge:
  `jwst/build_compact_object_source_release_receipts.py`
  mirrors the JWST compact-object simulator claim ladder in this paper-stack
  code tree. It is diagnostic only, not a particle prediction and not an OPH
  confirmation claim.
- Fractional quotient-sector sandbox:
  `fractional/build_fractional_quotient_receipts.py`
  mirrors the fractional exciton/FQAH simulator claim ladder. It is
  diagnostic only and remains blocked at the material-specific Hamiltonian proof.
- High-energy messenger coefficient emission:
  `uhe/build_uhe_coefficient_emission_receipts.py`
  mirrors the source-only UHE coefficient-emitter receipt ladder. It freezes the
  quotient, source law, compact-engine source loads, baseline, feature map,
  moment targets, solver, no-UHE-data-use DAG, and common-source lock without
  analyzing neutrino, cosmic-ray, or gamma event data.
- Compact-transient receipt scaffold:
  `compact_transients/build_compact_transient_receipts.py`
  mirrors the compact record-surface ladder for FRBs, old-host compact
  sources, and black-hole recycling. It is diagnostic by default at
  `CR2_CONDITIONAL_PHENOMENOLOGY`, blocks physical-prediction promotion until
  control, refinement, freeze, and held-out likelihood receipts exist, and fails
  closed if a black-hole generation prior reads ringdown residuals.
- rendered public surface:
  `scripts/build_results_status_table.py`

The live Higgs lane is a conditional split candidate on the declared
electroweak running, matching, and threshold surface. It emits
`m_H = 125.1995304097179 GeV` and a companion top coordinate
`m_t = 172.3523553288312 GeV`.
The row is not promoted as a strict source-only particle prediction until the
D10 target-free repair closes.
The exact public running-top row uses the PDG 2025 cross-section entry
`Q007TP4`.
The auxiliary direct-top average `Q007TP` is compare-only; #207 is closed as a
corpus-limited codomain no-go.

The hierarchy proof bundle is a separate audit lane from the rounded `1.63094`
calibration carrier. It records the public endpoint branch
`P_C = 1.630968209403959324879279847782648941`,
`alpha_U(P_C) = 0.041124336195630495`, and
`v/E_star = 2.0199803239725553e-17`, together with a source-audit branch that
keeps the public Thomson endpoint out of the upstream solve. Its Krawczyk
certificate proves a unique source zero for the declared `R_U` formula stack
inside the supplied interval. The local/global bridge closes on the exact
capacity value
`N_CRC^EW = 3.5323546226929906511187512962330547600462e122`,
with zero bridge residual, the 12-port screen sieve, the 24-slot oriented
repair register, and `epsilon_H = 0` on the selected source-to-Higgs branch.
The pixel-screen receipt records the same selected `(P_*,N_CRC^EW)` pair as an
equal-area screen chart:
`K_cell = 4*N_CRC^EW/P_*`,
`K_cell*(P_*/4)=N_CRC^EW`,
`Lambda_CRC*l_star^2=3*pi/N_CRC^EW`, and
`Lambda_CRC*a_cell=3*pi*P_*/N_CRC^EW=12*pi/K_cell`.
It is a summary receipt of the existing certificates, not an SI Lambda or
primitive-carrier promotion.
Full SI gravity remains gated by the no-G clock stack.

## Main Outputs

- claim table:
  [RESULTS_STATUS.md](RESULTS_STATUS.md)
- systematic open-gap ledger after the compressed `P`-trunk simplification:
  [DERIVATION_GAP_LEDGER.md](DERIVATION_GAP_LEDGER.md)
- bundled closure campaign for addressing the open blockers as coupled
  packets:
  [GAP_BUNDLE_CAMPAIGN.md](GAP_BUNDLE_CAMPAIGN.md)
  and [campaigns/gap_bundle](campaigns/gap_bundle)
- exact-fits-only diagnostic surface:
  [EXACT_FITS_ONLY.md](EXACT_FITS_ONLY.md)
- exact non-hadron mass bundle:
  [EXACT_NONHADRON_MASSES.md](EXACT_NONHADRON_MASSES.md)
- machine-readable claim table:
  [results_status.json](results_status.json)
- machine-readable exact-fits-only surface:
  [exact_fits_only.json](exact_fits_only.json)
- machine-readable exact non-hadron mass bundle:
  [exact_nonhadron_masses.json](exact_nonhadron_masses.json)
- frozen claim-table artifact:
  [status_table_forward_current.json](runs/status/status_table_forward_current.json)
- frozen exact-fits-only artifact:
  [exact_fits_only_current.json](runs/status/exact_fits_only_current.json)
- frozen exact non-hadron mass bundle:
  [exact_nonhadron_masses_current.json](runs/status/exact_nonhadron_masses_current.json)
- derivation graph:
  [particle_mass_derivation_graph.svg](particle_mass_derivation_graph.svg)
- electroweak hierarchy proof bundle:
  [hierarchy](hierarchy)

## Quark Closure Surface

- maximal theorem-emitted package artifact:
  [quark_maximal_theorem_emitted_package.json](runs/flavor/quark_maximal_theorem_emitted_package.json)
- target-free mass bridge artifacts:
  [light_quark_overlap_defect_value_law.json](runs/flavor/light_quark_overlap_defect_value_law.json)
  and
  [quark_d12_t1_value_law.json](runs/flavor/quark_d12_t1_value_law.json)
- physical-sheet contract artifact:
  [quark_lane_closure_contract.json](runs/flavor/quark_lane_closure_contract.json)
- algebraic-collapse artifact:
  [quark_absolute_readout_algebraic_collapse.json](runs/flavor/quark_absolute_readout_algebraic_collapse.json)
- exact `current_family_only` PDG witness:
  [quark_current_family_exact_pdg_theorem.json](runs/flavor/quark_current_family_exact_pdg_theorem.json)
- restricted current-family common-refinement transport-frame sector-attached lift:
  [quark_current_family_transport_frame_sector_attached_lift.json](runs/flavor/quark_current_family_transport_frame_sector_attached_lift.json)
- restricted current-family common-refinement transport-frame physical sigma lift theorem:
  [quark_current_family_transport_frame_physical_sigma_lift_theorem.json](runs/flavor/quark_current_family_transport_frame_physical_sigma_lift_theorem.json)
- restricted current-family common-refinement transport-frame strengthened physical sigma lift theorem:
  [quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json](runs/flavor/quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json)
- restricted current-family common-refinement transport-frame absolute sector readout theorem:
  [quark_current_family_transport_frame_absolute_sector_readout_theorem.json](runs/flavor/quark_current_family_transport_frame_absolute_sector_readout_theorem.json)
- restricted current-family common-refinement transport-frame exact PDG completion:
  [quark_current_family_transport_frame_exact_pdg_completion.json](runs/flavor/quark_current_family_transport_frame_exact_pdg_completion.json)
- restricted current-family common-refinement transport-frame exact forward Yukawas:
  [quark_current_family_transport_frame_exact_forward_yukawas.json](runs/flavor/quark_current_family_transport_frame_exact_forward_yukawas.json)
- restricted current-family common-refinement end-to-end exact PDG derivation chain:
  [quark_current_family_end_to_end_exact_pdg_derivation_chain.json](runs/flavor/quark_current_family_end_to_end_exact_pdg_derivation_chain.json)
- missing source-only sigma selector gate:
  [quark_sigma_source_datum_no_target_leak_required.json](runs/flavor/quark_sigma_source_datum_no_target_leak_required.json)
- selected-class physical sigma-datum descent witness:
  [quark_public_physical_sigma_datum_descent.json](runs/flavor/quark_public_physical_sigma_datum_descent.json)
- selected-class conditional end-to-end Yukawa witness:
  [quark_public_exact_yukawa_end_to_end_theorem.json](runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json)
- public exact Yukawa promotion frontier:
  [quark_public_exact_yukawa_promotion_frontier.json](runs/flavor/quark_public_exact_yukawa_promotion_frontier.json)
- selected-class public closure summary:
  [quark_public_strengthened_physical_sigma_lift_frontier.json](runs/flavor/quark_public_strengthened_physical_sigma_lift_frontier.json)

These artifacts fix the quark claim boundary on the local code surface. The
maximal theorem-emitted package consists of the D12 mass ray, the negative
selector `sigma_ref`, and the restricted-scope affine mean package with
`g_ch = 0.9231656602589082` on `shared_budget_only` and
`(g_u, g_d) = (0.7797392875757557, 0.12172551081512113)` on
`current_family_only`. A separate target-free mass bridge internalizes
`Delta_ud^overlap = (1/6) * log(c_d / c_u)`, equivalently
`quark_d12_t1_value_law`, on the emitted D12 ray. The exact six-mass witness
closes on `current_family_only` and matches the official PDG 2025 API running
quark reference sextet exactly on that declared target surface; its top coordinate is
the PDG 2025 cross-section mass entry. The auxiliary direct-top entry is
compare-only; [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207)
is closed as a corpus-limited codomain no-go. A separate restricted theorem chain on the explicit
  `current_family_common_refinement_transport_frame_only` carrier emits a
sector-attached `Sigma_ud^phys` element, the exact physical sigma datum, the
restricted absolute sector readout, the same sextet, and explicit exact
forward Yukawas `Y_u` and `Y_d`. The declared end-to-end exact chain closes on
that carrier. A separate direct public descent witness proves representative
independence on the selected public physical quark frame class chosen by `P`,
but the physical sigma datum is target-derived on the current public surface.
Promotion requires `QUARK_SIGMA_SOURCE_SELECTOR` and
`NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE`. Once a source-only sigma datum is
supplied, the affine mean law emits `(g_u, g_d)` algebraically, and the exact
forward construction emits the same sextet together with explicit exact forward
Yukawas `Y_u` and `Y_d`. Until then this is selected-class conditional support
only, not a public source-only quark mass prediction. It does not claim a global
classification of all quark frame classes. The upstream generation-bundle
transfer route is an alternative source-sigma route, not a closed theorem.

## Typical Rebuild

From `reverse-engineering-reality/code/particles`:

```bash
python3 calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py
python3 calibration/derive_d10_ew_source_transport_readout.py
python3 calibration/derive_d10_ew_exactness_audit.py
python3 calibration/derive_d11_declared_calibration_surface.py
python3 calibration/derive_d11_forward_seed.py
python3 calibration/derive_d11_forward_seed_promotion_certificate.py
python3 calibration/derive_d11_live_exact_higgs_promotion.py
python3 calibration/derive_d11_live_exact_split_pair_theorem.py
python3 calibration/derive_d11_reference_exact_adapter.py
python3 neutrino/derive_neutrino_weighted_cycle_repair.py
python3 neutrino/derive_neutrino_bridge_rigidity_theorem.py
python3 neutrino/derive_neutrino_absolute_attachment_theorem.py
python3 neutrino/export_forward_neutrino_closure_bundle.py
python3 neutrino/derive_neutrino_two_parameter_exact_adapter.py
python3 hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py
python3 hadron/derive_stable_channel_sequence_evaluation.py
python3 hadron/derive_current_hadron_lane_audit.py
python3 jwst/build_compact_object_source_release_receipts.py
python3 fractional/build_fractional_quotient_receipts.py
python3 scripts/build_results_status_table.py
python3 scripts/build_derivation_gap_ledger.py
python3 scripts/build_exact_fit_surface.py
python3 scripts/build_exact_nonhadron_mass_bundle.py
python3 scripts/generate_mass_derivation_svg.py
```

## One-Shot CLI Table

For a disposable runtime rebuild that re-runs the active D10/D11/UV builders,
stages the canonical flavor/lepton/neutrino public-surface artifacts,
and prints the resulting particle claim table directly in the terminal:

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
  jwst/test_compact_object_source_release_receipts.py \
  fractional/test_fractional_quotient_receipts.py \
  hierarchy/test_hierarchy_bundle.py \
  uv/test_oph_bd_threshold_spectrum_receipts.py \
  test_results_status_candidate_policy.py \
  test_results_status_quark_promotion_policy.py \
  test_results_status_structural_rows.py \
  test_predictive_builders_reference_free.py
```

## Paper Surface

The code here feeds the particle paper:

- [deriving_the_particle_zoo_from_observer_consistency.tex](../../paper/deriving_the_particle_zoo_from_observer_consistency.tex)
- [deriving_the_particle_zoo_from_observer_consistency.pdf](../../paper/deriving_the_particle_zoo_from_observer_consistency.pdf)

## License And Patent Policy

This particle code surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../PATENTS.md).
