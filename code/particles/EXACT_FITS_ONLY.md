# Exact Fits Only

Generated: `2026-05-06T07:44:09Z`

This surface lists exact target matches on declared OPH carriers. It separates theorem-grade selected-class outputs from compare-only and carrier-restricted exact surfaces.
For quarks, the selected-class theorem and its supporting exact carriers coincide with the official PDG 2025 API running-quark target surface.

## Electroweak Frozen-Target Exact Pair

- Fit kind: `exact_frozen_target_compare_only_adapter`
- Scope: `frozen_authoritative_target_surface`
- Promotable: `false`
- Source artifact: `code/particles/runs/calibration/d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json`
- Max absolute residual: `0.0`
- Note: Exact on the frozen-authoritative D10 repair surface. The public theorem surface is the target-free source-only emission, which is separate and differs only at the `1e-8 GeV` scale.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m_W` | `80.377` | `80.377` |
| `m_Z` | `91.18797809193725` | `91.18797809193725` |

## Higgs/Top Reference Exact Adapter

- Fit kind: `exact_target_anchored_compare_only_inverse_slice`
- Scope: `compare_only_inverse_slice`
- Promotable: `false`
- Source artifact: `code/particles/runs/calibration/d11_reference_exact_adapter.json`
- Max absolute residual: `0.0`
- Note: Exact only as a compare-only inverse slice on the D11 Jacobian. The live D11 theorem lane uses the source-only split theorem on the declared D10/D11 surface. That theorem emits `m_H = 125.1995304097179 GeV` and a companion top coordinate `m_t = 172.3523553288312 GeV`. The Higgs row lands on the PDG 2025 Higgs average at quoted precision. The exact public running-top row uses the PDG 2025 cross-section entry `Q007TP4`. The auxiliary direct-top average `Q007TP` is compare-only; [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go by `code/particles/runs/calibration/direct_top_bridge_contract.json`.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m_H` | `125.1995304097179` | `125.1995304097179` |
| `m_t` | `172.3523553288312` | `172.3523553288312` |

## Charged Current-Family Exact Witness

- Fit kind: `exact_target_anchored_current_family_witness`
- Scope: `current_family_only`
- Promotable: `false`
- Source artifact: `code/particles/runs/leptons/lepton_current_family_exact_readout.json`
- Max absolute residual: `1.1102230246251565e-15`
- Note: Exact on the ordered charged eigenvalue triple, with a closed ordered-three-point readout theorem inside `current_family_only`, and with the scoped affine coordinate `A_ch_current_family` closed on that same exact family. The charged theorem lane does not emit a theorem-grade absolute anchor; [#201](https://github.com/FloatingPragma/observer-patch-holography/issues/201) is closed as a corpus-limited no-go by `code/particles/runs/leptons/charged_end_to_end_impossibility_theorem.json`.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m_e` | `0.0005109989499999994` | `0.0005109989499999999` |
| `m_mu` | `0.10565837550000004` | `0.10565837550000001` |
| `m_tau` | `1.7769324651340912` | `1.77693246513409` |

## Quark Selected-Class Exact Theorem

- Fit kind: `selected_class_theorem_grade_exact_forward_quark_closure`
- Scope: `selected_public_physical_quark_frame_class_only`
- Promotable: `true`
- Source artifact: `code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json`
- Max absolute residual: `5.684341886080802e-14`
- Note: Exact theorem on the selected public physical quark frame class chosen by `P`. `oph_quark_public_physical_sigma_datum_descent` makes the exact physical sigma datum target-free public on that selected class, and `oph_quark_public_exact_yukawa_end_to_end_theorem` emits the exact PDG 2025 running-quark sextet together with explicit exact forward Yukawas `Y_u` and `Y_d`. The top coordinate uses PDG summary `Q007TP4`. The auxiliary direct-top entry `Q007TP` is compare-only; [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go by `code/particles/runs/calibration/direct_top_bridge_contract.json`. This is selected-class closure only. It does not claim a global classification of all quark frame classes.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m_u` | `0.0021600000000000005` | `0.00216` |
| `m_c` | `1.2729999999999992` | `1.273` |
| `m_t` | `172.35235532883115` | `172.3523553288312` |
| `m_d` | `0.004699999999999999` | `0.0047` |
| `m_s` | `0.09349999999999999` | `0.0935` |
| `m_b` | `4.182999999999994` | `4.183` |

## Quark Current-Family Exact Witness

- Fit kind: `exact_target_anchored_current_family_witness`
- Scope: `current_family_only`
- Promotable: `false`
- Source artifact: `code/particles/runs/flavor/quark_current_family_exact_readout.json`
- Max absolute residual: `1.1368683772161603e-13`
- Note: Exact on the official PDG 2025 API running-quark target surface on the ordered three-point quark family witness, with the internal same-family quadratic readout closed on the fixed carrier and the selected-sheet exact closure packaged on `sigma_ref`. The top coordinate uses PDG summary `Q007TP4`. The auxiliary direct-top entry `Q007TP` is compare-only; [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go. The declared scope is `current_family_only`. A separate restricted theorem chain emits a sector-attached `Sigma_ud^phys` element on the explicit `current_family_common_refinement_transport_frame_only` carrier, and the merged transport-frame theorem reconstructs the same running sextet exactly on `current_family_common_refinement_transport_frame_only`. The declared transport-frame chain also closes explicit exact forward Yukawas `Y_u` and `Y_d` with certification status `forward_matrix_certified`, and the full declared-carrier chain is recorded in `oph_quark_current_family_end_to_end_exact_pdg_derivation_chain`. A separate target-free mass bridge closes `Delta_ud_overlap = (1/6) * log(c_d / c_u)`, equivalently `quark_d12_t1_value_law`, on the emitted D12 ray. A separate public theorem closes on the selected public physical quark frame class chosen by `P`: `oph_quark_public_physical_sigma_datum_descent` makes the exact physical sigma datum target-free public on that selected class, and `oph_quark_public_exact_yukawa_end_to_end_theorem` emits the same exact sextet together with explicit exact forward Yukawas `Y_u` and `Y_d`. This entry remains an exact-fit surface rather than that public selected-class theorem.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m_u` | `0.0021599999999999996` | `0.00216` |
| `m_c` | `1.2729999999999995` | `1.273` |
| `m_t` | `172.3523553288311` | `172.3523553288312` |
| `m_d` | `0.004700000000000002` | `0.0047` |
| `m_s` | `0.09349999999999999` | `0.0935` |
| `m_b` | `4.182999999999994` | `4.183` |

## Neutrino Two-Parameter Exact Adapter

- Fit kind: `exact_two_observable_compare_only_segment_adapter`
- Scope: `compare_only_two_parameter_segment_adapter`
- Promotable: `false`
- Source artifact: `code/particles/runs/neutrino/neutrino_two_parameter_exact_adapter.json`
- Max absolute residual: `4.0657581468206416e-20`
- Note: Exact compare-only fit to both representative PDG central splittings by moving along the explicit positive selector segment and then rescaling with one positive `lambda_nu`. It is diagnostic-only after the emitted weighted-cycle bridge-rigidity and absolute-attachment theorems. On that same exact compare-only branch, the explicit bridge coordinates are `B_nu = 6.69675975` and `C_nu = 0.99952948`, but they remain sidecars and must not feed back into theorem state.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m1_eV` | `0.01745663294772044` | `n/a` |
| `m2_eV` | `0.019484199595350048` | `n/a` |
| `m3_eV` | `0.053081390655025595` | `n/a` |
| `delta_m21_sq_eV2` | `7.490000000000005e-05` | `7.49e-05` |
| `delta_m31_sq_eV2` | `0.0025129` | `0.0025129` |
| `delta_m32_sq_eV2` | `0.002438` | `0.002438` |

## Neutrino Atmospheric Only Exact Adapter

- Fit kind: `exact_single_observable_compare_only_adapter`
- Scope: `compare_only`
- Promotable: `false`
- Source artifact: `code/particles/runs/neutrino/neutrino_compare_only_scale_fit.json`
- Exact matched observable: `Delta m32^2`
- Note: Exact only for one splitting observable on the repaired weighted-cycle family. The same artifact states that no single `lambda_nu` hits both central splittings exactly.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m1_eV` | `0.017456756479999103` | `n/a` |
| `m2_eV` | `0.019484260653687455` | `n/a` |
| `m3_eV` | `0.053081413067295344` | `n/a` |
| `delta_m21_sq_eV2` | `7.489806641884242e-05` | `7.49e-05` |
| `delta_m31_sq_eV2` | `0.0025128980664188426` | `n/a` |
| `delta_m32_sq_eV2` | `0.002438` | `0.002438` |

## Neutrino Solar Only Exact Adapter

- Fit kind: `exact_single_observable_compare_only_adapter`
- Scope: `compare_only`
- Promotable: `false`
- Source artifact: `code/particles/runs/neutrino/neutrino_compare_only_scale_fit.json`
- Exact matched observable: `Delta m21^2`
- Note: Exact only for one splitting observable on the repaired weighted-cycle family. The same artifact states that no single `lambda_nu` hits both central splittings exactly.

| Observable | Value | Reference |
| --- | ---: | ---: |
| `m1_eV` | `0.017456981811834547` | `n/a` |
| `m2_eV` | `0.01948451215654942` | `n/a` |
| `m3_eV` | `0.05308209824224456` | `n/a` |
| `delta_m21_sq_eV2` | `7.49e-05` | `7.49e-05` |
| `delta_m31_sq_eV2` | `0.0025129629398205804` | `n/a` |
| `delta_m32_sq_eV2` | `0.0024380629398205803` | `0.002438` |
