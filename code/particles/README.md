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
  and group-action geometry only. The receipt distinguishes the projective
  ray's two-torus stabilizer from the nonzero vacuum vector's `U(1)_Q`
  stabilizer before the quantitative D11 lane
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
  compare-only. The shared-basis identity cancels the charged-lepton matrix by
  construction and does not supply the missing physical charged-lepton basis.
  Earlier intrinsic builders also exported left SVD vectors instead of the
  Majorana Takagi matrix; that implementation error is corrected. This lane
  emits neither a physical PMNS matrix nor absolute neutrino masses.
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

The W/Z lane has three distinct numeric surfaces:

- the selected-carrier chart emits
  `80.38629169244275 GeV` and `91.18290444674243 GeV`;
- the archived D10 value-law candidate emits
  `80.37700001539531 GeV` and `91.18797807794321 GeV`;
- the frozen `80.377 GeV`, `91.18797809193725 GeV` pair is a
  comparison-only adapter.

The archived value law is an exact implication of the QT1--QT5 quotient-path
certificate. The repository does not emit that certificate. A proof-producing
enumeration must construct the quotient canonicalizer and finite path lists,
verify the incidence counts, exact color weighting convention, a
representation-derived `Z_6` invariant-rank fraction (not merely group order),
the fibre Gram form and residual pairing, and prove a
positive MAR gap over an explicit deformation class. The executable QT receipt
checks only the downstream algebra and chart Jacobian; it fails closed on every
source-entailment field.

The live Higgs lane is a conditional split candidate on the declared
electroweak running, matching, and threshold surface. It emits
`m_H = 125.1995304097179 GeV` and a companion top coordinate
`m_t = 172.3523553288312 GeV`.
The row is not a strict source-only particle prediction. Its inherited gates
include the final source root, an independently physical `E_star`, QT1--QT5,
the concrete RG/matching/scheme receipt, the DS1--DS5 D11 split-character and
rigidity certificate, and top-threshold
control, a complex-pole/uncertainty certificate, and a prospective source DAG.
The full proof-producing implementation contract is documented in
[`calibration/WZH_SOURCE_CLOSURE_CAMPAIGN.md`](calibration/WZH_SOURCE_CLOSURE_CAMPAIGN.md).
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
These are dimensionless hierarchy/naturality statements. They do not turn
`v/E_star` into a GeV mass until the physical meaning and normalization of
`E_star` are independently source-closed.
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
- target-anchored `current_family_only` mixed-convention audit witness:
  [quark_current_family_exact_pdg_theorem.json](runs/flavor/quark_current_family_exact_pdg_theorem.json)
- restricted current-family common-refinement transport-frame sector-attached lift:
  [quark_current_family_transport_frame_sector_attached_lift.json](runs/flavor/quark_current_family_transport_frame_sector_attached_lift.json)
- restricted current-family common-refinement transport-frame physical sigma lift theorem:
  [quark_current_family_transport_frame_physical_sigma_lift_theorem.json](runs/flavor/quark_current_family_transport_frame_physical_sigma_lift_theorem.json)
- restricted current-family common-refinement transport-frame strengthened physical sigma lift theorem:
  [quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json](runs/flavor/quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json)
- restricted current-family common-refinement transport-frame absolute sector readout theorem:
  [quark_current_family_transport_frame_absolute_sector_readout_theorem.json](runs/flavor/quark_current_family_transport_frame_absolute_sector_readout_theorem.json)
- restricted current-family common-refinement transport-frame target-audit completion:
  [quark_current_family_transport_frame_exact_pdg_completion.json](runs/flavor/quark_current_family_transport_frame_exact_pdg_completion.json)
- restricted current-family common-refinement transport-frame dimensionful mass textures:
  [quark_current_family_transport_frame_exact_forward_yukawas.json](runs/flavor/quark_current_family_transport_frame_exact_forward_yukawas.json)
- restricted current-family common-refinement end-to-end target-audit chain:
  [quark_current_family_end_to_end_exact_pdg_derivation_chain.json](runs/flavor/quark_current_family_end_to_end_exact_pdg_derivation_chain.json)
- theorem-grade source-spread non-identifiability obstruction:
  [quark_sigma_source_nonidentifiability_obstruction.json](runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json)
- source-only spread gate projected from that obstruction:
  [quark_sigma_source_datum_no_target_leak_required.json](runs/flavor/quark_sigma_source_datum_no_target_leak_required.json)
- running-mass scheme and physical-Yukawa normalization obstruction:
  [quark_running_mass_scheme_convention_obstruction.json](runs/flavor/quark_running_mass_scheme_convention_obstruction.json)
- selected-class physical sigma-datum descent witness:
  [quark_public_physical_sigma_datum_descent.json](runs/flavor/quark_public_physical_sigma_datum_descent.json)
- selected-class target-audit mass-texture wrapper:
  [quark_public_exact_yukawa_end_to_end_theorem.json](runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json)
- public exact Yukawa promotion frontier:
  [quark_public_exact_yukawa_promotion_frontier.json](runs/flavor/quark_public_exact_yukawa_promotion_frontier.json)
- selected-class public closure summary:
  [quark_public_strengthened_physical_sigma_lift_frontier.json](runs/flavor/quark_public_strengthened_physical_sigma_lift_frontier.json)

These artifacts fix the quark claim boundary on the local code surface. After
all target rows, exact target witnesses, fitted spreads, and residuals against
them are removed, the source equations determine two ordered profile rays but
leave their endpoint spans independent. The compatible spread fiber is exactly
`(R_{>0})^2`. Its free rescaling action preserves the source identities and
changes the affine mass readout, so it is a physical non-identifiability rather
than a gauge redundancy. The selected-class descent proves representative
independence only after a spread datum is attached; it does not select the
datum. The edge-statistics candidate also leaves two unfixed correction
coefficients and begins from a hand-written family-kernel template.

The stronger no-extra-axiom result is emitted in
[`quark_axiom_level_yukawa_moduli_nonidentifiability.json`](runs/flavor/quark_axiom_level_yukawa_moduli_nonidentifiability.json).
MAR orders the structural tuple `(chi_cpl, N_nonab, N_c, N_g)` and contains no
Yukawa eigenvalue coordinate. Independent positive rescalings of the two
centered quark profiles therefore remain physically distinct equal-score MAR
minima. The positive family has no smallest member, so MAR cannot be read as a
hidden numerical Yukawa selector under its stated definition.

The same-family and common-refinement artifacts reproduce their chosen target
coordinates after target inversion. They are audit surfaces, not source-only
predictions. Their packet combines light-quark `MSbar` coordinates at 2 GeV,
charm and bottom `MSbar` coordinates at self-scale, and a separate top pole
extraction. The stored GeV-valued matrices therefore certify mass textures,
not physical dimensionless Yukawa matrices. A physical Yukawa construction
would require source-emitted RG trajectories, common-scale transport with
threshold matching, a top conversion, the running Higgs expectation value in
the same scheme, and `y_q(mu) = sqrt(2) m_q(mu) / v(mu)`. Numeric public quark
rows remain withheld. Reopening the source theorem requires a new OPH source
observable that breaks the independent positive-rescaling action without a
dependency path to target data.

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
