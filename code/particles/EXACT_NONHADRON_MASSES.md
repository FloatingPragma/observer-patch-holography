# Exact Non-Hadron Masses

Generated: `2026-05-06T07:44:09Z`

This bundle gives one exact mass output for every non-hadron particle on the declared OPH surfaces.
It records exact-output surfaces rather than one uniform theorem tier.
For quarks, the exact theorem surface matches the official PDG 2025 API running-quark target surface on the selected public physical quark frame class chosen by `P`.
The same selected-class theorem emits explicit exact forward Yukawas `Y_u` and `Y_d`, and the same sextet is also realized on `current_family_only` and on the restricted current-family common-refinement transport-frame carrier.
For charged leptons, this bundle records the exact same-family witness surface. The theorem surface also contains the live same-label `q_e` readback, a source-side determinant character for a fixed formal source multiplicity vector, a conditional determinant-line lift, and an algebraic charged-mass readout once theorem-grade `A_ch(P)` is given. Issue #201 is closed as a corpus-limited no-go: the available corpus does not attach that source-side character to the physical charged determinant line.
The top coordinate uses the PDG 2025 cross-section mass entry `Q007TP4`. The auxiliary direct-top entry `Q007TP` is compare-only; [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go by `code/particles/runs/calibration/direct_top_bridge_contract.json`.

| Particle | Exact Mass | Kind | Scope | Source |
| --- | ---: | --- | --- | --- |
| Photon | `0.0 GeV` | `structural_zero` | `structural` | `structural_gauge_redundancy_surface` |
| Gluon | `0.0 GeV` | `structural_zero` | `structural` | `structural_color_gauge_surface` |
| Graviton | `0.0 GeV` | `structural_zero` | `structural` | `structural_diffeomorphism_redundancy_surface` |
| W Boson | `80.377 GeV` | `exact_frozen_target_compare_only_adapter` | `frozen_authoritative_target_surface` | `code/particles/runs/calibration/d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json` |
| Z Boson | `91.18797809193725 GeV` | `exact_frozen_target_compare_only_adapter` | `frozen_authoritative_target_surface` | `code/particles/runs/calibration/d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json` |
| Higgs Boson | `125.1995304097179 GeV` | `exact_source_only_higgs_top_split_calibration_theorem` | `declared_d10_d11_running_matching_threshold_surface_only` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |
| Electron | `0.0005109989499999994 GeV` | `exact_target_anchored_current_family_witness` | `current_family_only` | `code/particles/runs/leptons/lepton_current_family_exact_readout.json` |
| Muon | `0.10565837550000004 GeV` | `exact_target_anchored_current_family_witness` | `current_family_only` | `code/particles/runs/leptons/lepton_current_family_exact_readout.json` |
| Tau | `1.7769324651340912 GeV` | `exact_target_anchored_current_family_witness` | `current_family_only` | `code/particles/runs/leptons/lepton_current_family_exact_readout.json` |
| Up Quark | `0.0021600000000000005 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Charm Quark | `1.2729999999999992 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Top Quark | `172.35235532883115 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Down Quark | `0.004699999999999999 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Strange Quark | `0.09349999999999999 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Bottom Quark | `4.182999999999994 GeV` | `selected_class_theorem_grade_exact_forward_quark_closure` | `selected_public_physical_quark_frame_class_only` | `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json` |
| Electron Neutrino | `0.017454720257976796 eV` | `theorem_grade_weighted_cycle_absolute_attachment` | `weighted_cycle_bridge_rigid_absolute_family` | `code/particles/runs/neutrino/neutrino_absolute_attachment_theorem.json` |
| Muon Neutrino | `0.019481987935919015 eV` | `theorem_grade_weighted_cycle_absolute_attachment` | `weighted_cycle_bridge_rigid_absolute_family` | `code/particles/runs/neutrino/neutrino_absolute_attachment_theorem.json` |
| Tau Neutrino | `0.05307522145074924 eV` | `theorem_grade_weighted_cycle_absolute_attachment` | `weighted_cycle_bridge_rigid_absolute_family` | `code/particles/runs/neutrino/neutrino_absolute_attachment_theorem.json` |
