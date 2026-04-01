#!/usr/bin/env python3
"""Build the `/particles`-native public prediction surface.

Chain role: assemble the live per-sector outputs from the active local
derivation chain into one public candidate-or-gap table.

Mathematics: this file does not derive masses itself; it applies promotion
policy, ledger mapping, residual reporting, and surface provenance rules.

OPH-derived inputs: the local `/particles` calibration, flavor, neutrino, and
hadron artifacts only. No legacy ancillary predictor surface is imported here.

Output: `RESULTS_STATUS.md`, `results_status.json`, and the machine-readable
public surface snapshot used for audits and progress tracking.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
LEDGER_YAML = ROOT / "particles" / "ledger.yaml"
DEFAULT_MD_OUT = ROOT / "particles" / "RESULTS_STATUS.md"
DEFAULT_JSON_OUT = ROOT / "particles" / "results_status.json"
DEFAULT_FORWARD_OUT = ROOT / "particles" / "runs" / "status" / "status_table_forward_current.json"
UV_BW_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_internalization_scaffold.json"
UV_BW_PRELIMIT_SYSTEM = ROOT / "particles" / "runs" / "uv" / "bw_realized_transported_cap_local_system.json"
UV_BW_FIXED_LOCAL_COLLAR_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
UV_BW_CARRIED_SCHEDULE = ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"
UV_BW_CAP_PAIR_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_scaling_limit_cap_pair_extraction_scaffold.json"
UV_BW_RIGIDITY_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_ordered_cut_pair_rigidity_scaffold.json"
FORWARD_YUKAWAS = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
QUARK_SECTOR_MEAN_SPLIT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
D10_SOURCE_TRANSPORT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
D11_FORWARD_SEED = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
FORWARD_CHARGED_LEPTONS = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
FORWARD_NEUTRINO_BUNDLE = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
NEUTRINO_EXACT_BLOCKERS = ROOT / "particles" / "runs" / "neutrino" / "exact_blocking_items.json"
NEUTRINO_WEIGHTED_CYCLE_REPAIR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
NEUTRINO_TWO_PARAMETER_EXACT_ADAPTER = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
NEUTRINO_LAMBDA_BRIDGE_CANDIDATE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"
PUBLIC_SURFACE_KIND = "particles_native_candidate_or_gap_surface"
P_DEFAULT = 1.63094
LOG_DIM_H_DEFAULT = 1.0e122


GROUP_ORDER = ["Bosons", "Leptons", "Quarks", "Hadrons"]
NEUTRINO_OSCILLATION_SOURCE_URL = "https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf"
NEUTRINO_OSCILLATION_REFERENCE_LABEL = "PDG 2025 NO reference"
NEUTRINO_PDG_2025_NO_CENTRAL = {
    "theta12_deg": 33.68,
    "theta23_deg": 43.3,
    "theta13_deg": 8.56,
    "delta_deg": 212.0,
    "delta_m21_sq_eV2": 7.49e-5,
    "delta_m32_sq_eV2": 2.438e-3,
}
NEUTRINO_PDG_2025_NO_1SIGMA = {
    "theta12_deg": {"plus": 0.73, "minus": 0.70},
    "theta23_deg": {"plus": 1.0, "minus": 0.9},
    "theta13_deg": {"plus": 0.11, "minus": 0.11},
    "delta_deg": {"plus": 26.0, "minus": 41.0},
    "delta_m21_sq_eV2": {"plus": 0.19e-5, "minus": 0.20e-5},
    "delta_m32_sq_eV2": {"plus": 0.021e-3, "minus": 0.019e-3},
}
D10_MASS_PAIR_NOTE = (
    "Derived from the D10 calibration chain "
    "`derive_d10_ew_observable_family.py -> derive_d10_ew_source_transport_pair.py -> "
    "derive_d10_ew_population_evaluator.py -> derive_d10_ew_exact_closure_beyond_current_carrier.py -> "
    "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py -> "
    "derive_d10_ew_tau2_current_carrier_obstruction.py -> derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py -> "
    "derive_d10_ew_exact_mass_pair_chart_current_carrier.py -> derive_d10_ew_repair_branch_beyond_current_carrier.py -> "
    "derive_d10_ew_repair_target_point_diagnostic.py -> derive_d10_ew_w_anchor_neutral_shear_factorization.py -> "
    "derive_d10_ew_target_free_repair_value_law.py -> derive_d10_ew_source_transport_readout.py`. "
    "Calibration here means that the shared pixel scale `P` is first fixed on the declared D10 running/matching surface, which in turn fixes the D10 source basis "
    "`(alpha2_mz, alphaY_mz, eta_source, v_report)`. "
    "The live forward transmutation certificate makes that order explicit on disk: the same source-only basis reconstructs `alpha_U`, the unified diffusion parameter `t_U = 4*pi^2*alpha_U`, and the transmutation exponent `t_tr = 2*pi / ((N_c + 1) * alpha_U)` without reading them back from measured couplings. "
    "The selected current-carrier chart is closed and remains explicit on disk, but the active public electroweak surface is now the target-free source-only theorem `EWTargetFreeRepairValueLaw_D10`. "
    "That theorem emits the repaired chart `(tau2_tree_exact, delta_n_tree_exact)` from the D10 source basis alone using `lambda_EW = eta_source^2 / (4 * beta_EW)`, then emits one coherent electroweak quintet from one repaired coupling pair. "
    "On the paper-facing theorem lane the transmutation factor is `beta_transmutation_EW = N_c + 1`; older overloaded beta ratios survive only on compare-only validation readouts. "
    "So the public D10 W/Z values are no longer freeze-once rows. The older freeze-once coherent repair law is retained only as compare-only validation and agrees with the target-free theorem to machine scale: about `+1.54e-08` GeV on `W` and `-1.40e-08` GeV on `Z`. "
    "That frozen-target repair pair is also surfaced separately on the exact-hit diagnostic side as `oph_d10_ew_w_anchor_neutral_shear_factorization`, where it hits the canonical `W/Z` references exactly on one frozen authoritative repair surface. "
    "This closes the electroweak mass-side lane on the Phase II calibration tier; the earlier source-only underdetermination theorem, minimal conditional route through `ColorBalancedQuadraticRepairDescent_D10`, and former candidate `EWTargetEmitter_D10` remain on disk only as historical scaffolding beneath the promoted theorem."
)
D11_NOTE = (
    "Derived from `derive_d11_forward_seed.py -> derive_d11_forward_seed_promotion_certificate.py`, which propagates the D10 gauge core into the compact D11 forward seed, certifies the fixed-ray forward path, and reads out the D11 mass row from the Jacobian surface. "
    "A separate exact-hit sidecar is now also on disk as `oph_d11_reference_exact_adapter`: it solves the linear D11 Jacobian against the canonical Higgs/top reference pair and therefore hits those references exactly, but only as a compare-only inverse slice. "
    "The live public D11 rows remain the reference-free forward-seed outputs, not the inverse adapter."
)
CHARGED_CONTINUATION_NOTE = (
    "No public value is emitted yet on the theorem lane. A separate exact same-family witness is already on disk: `oph_lepton_current_family_exact_readout` reproduces the charged reference triple exactly on the same ordered eigenvalue family, its target-anchored ordered-three-point readout chain is closed within `current_family_only` by `oph_lepton_current_family_quadratic_readout_theorem`, and the scoped same-family affine coordinate is closed on that same witness by `oph_lepton_current_family_affine_anchor_theorem`; those closures do not promote the live charged theorem lane. The active charged path is "
    "`derive_charged_sector_local_current_support_obstruction_certificate.py -> "
    "derive_charged_sector_local_minimal_source_support_extension_emitter.py -> "
    "derive_charged_sector_local_support_extension_completion_law.py -> "
    "derive_charged_sector_local_support_extension_source_scalar_pair_readback.py -> "
    "derive_charged_d12_continuation_followup.py -> "
    "derive_charged_sector_local_support_extension_eta_source_readback.py -> "
    "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py -> "
    "derive_charged_lepton_absolute_scale_coordinate_shell.py -> "
    "derive_lepton_excitation_gap_map.py -> derive_lepton_log_spectrum_readout.py -> "
    "build_forward_charged_leptons.py`; the live same-carrier scalar order is "
    "`eta_source_support_extension_log_per_side` and then "
    "`sigma_source_support_extension_total_log_per_side`, with the smaller ordered source-scalar pair readback now explicit on disk. "
    "A representation-consistent absolute-scale shell is also explicit: future charged scale code must emit either "
    "`mu_e_absolute_log_candidate` or `g_e_linear_candidate` and convert exactly once via `g_e = exp(mu_e_absolute_log_candidate)`. "
    "But the present charged theorem still fixes only the centered charged log class modulo a common shift, so the absolute scale `g_e` remains unresolved on the live theorem lane. "
    "At theorem level, the exact waiting set is sharper than a standalone eta/sigma fit: the charged sector-response object is still only the latent candidate `C_hat_e^{cand}`, not a declared theorem-grade operator. Promoting that candidate is blocked by the upstream theorem `oph_generation_bundle_branch_generator_splitting`, reduced further to the clause `compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split`. The local corpus proves neither exact vanishing nor uniform quadratic smallness of that descended commutator yet. The exact minimal operator-side extension is already packaged on disk as `central_split_quadratic_commutator_transfer`, but `current_corpus_contains_theorem = false`: it has not been internalized into the live theorem corpus. On the absolute side, the charged equalizer route is now an explicit no-go under common-shift symmetry: the current theorem emits only the quotient class of charged logs modulo `(1,1,1)`, so no theorem-grade `g_e` or `Delta_e_abs` exists yet. The layered frontier is now explicit on disk as `oph_charged_absolute_frontier_factorization`: on the current surface the missing affine object is `A_ch`, while conditional on future theorem-grade `C_hat_e` promotion the post-promotion burden sharpens further to the descended physical affine scalar `mu_phys(Y_e)`, carried by `oph_charged_mu_physical_descent_reduction`, with `C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I`, `s_det(Y_e) = 3 * mu_phys(Y_e)`, and `A_ch(Y_e) = mu_phys(Y_e)`. "
    "The sharpest constructive route is therefore two-layered and still only an extension candidate, not a current-corpus closure. First, if the actual centered compressed generator factors through centered Schur-type `P->Q->P` feedback with a refinement-uniform middle-factor bound, then internalizing that packaged extension would close the transfer gap exactly when the descended commutator vanishes and otherwise only quadratically. On the current local certificate, the proxy margin would survive such an internalization whenever the uniform bound satisfies about `M < 119.5600535277701`. That route can promote only the centered proxy `C_hat_e^{cand}`. Beyond that promotion step, the post-promotion single slot is the descended physical scalar carried by `oph_charged_mu_physical_descent_reduction`: once the uncentered lift is refinement-stable on theorem-grade physical `Y_e`, the family cocycle is forced to vanish on same-`Y_e` refinement pairs and the remaining missing object is one `charged_physical_affine_scalar_mu`, from which the uncentered lift, determinant-line section, and affine charged anchor follow canonically, with `A_ch(Y_e) = (1/3) log det(Y_e)` or equivalently `A_ch(Y_e) = (1/3) tr(log Y_e)`. "
    "A D12 continuation bridge exists under the extra assumptions A1-A3 and gives eta = -6.729586682888832 and sigma = 8.154061112725994 with near-exact centered-log shape closure, "
    "but the theorem-grade lane still lacks emitted eta, sigma, and absolute scale. On that continuation bridge the compare-only absolute target would be `g_e* = 0.04577885783568762`, equivalently `Delta_e_abs* = 3.003986333402356`, and that target is kept strictly non-promotable until a theorem-grade absolute anchor `A_ch` exists on the live branch."
)
QUARK_CONTINUATION_NOTE = (
    "A separate exact same-family witness is already on disk too: `oph_quark_current_family_exact_readout` reproduces the six quark reference masses exactly on the same ordered three-point family, its internal current-family readout chain is closed through `oph_quark_current_family_quadratic_readout_theorem`, and the selected-sheet exact completion on `sigma_ref` is packaged as `oph_quark_current_family_selected_sheet_exact_closure`; its recorded scope field is still `current_family_only`, so it does not resolve the wrong-branch D12 CKM no-go or emit `quark_d12_t1_value_law`. Derived from the local quark chain "
    "`derive_quark_sector_mean_split.py -> derive_quark_sector_descent.py -> "
    "build_forward_yukawas.py -> derive_quark_d12_overlap_transport_law.py -> "
    "derive_quark_quadratic_even_transport_scalar.py -> derive_generation_bundle_same_label_physical_invariant_bundle.py -> "
    "derive_quark_scalarized_continuation_bundle.py -> derive_quark_d12_mass_branch_and_ckm_residual.py`, using the active reference-free forward Yukawa "
    "candidate on the `/particles` public surface. These rows are D12 "
    "continuation-level only: the March 28, 2026 consolidation against the OPH "
    "tier ledger shows that a nonzero light-quark pure-`B` source selector is "
    "not available at recovered-core tier. The active builder still waits on "
    "the pure-`B` payload pair "
    "`source_readback_u_log_per_side_and_source_readback_d_log_per_side`, but the broader D12 continuation branch is now scalarized: "
    "the mass side reduces to `Delta_ud_overlap` together with the quadratic-even scalar `eta_Q_centered`, and the mixing side closes on the same D12 branch because the forward Yukawa step already emits the same-label transport unitary `V_CKM^fwd = U_u^dagger U_d`, whose principal logarithm yields the gauge-fixed generator `K_CKM`. "
    "On the current same-family continuation branch, the exact mass-side reduction sharpens to the emitted ray `D12_ud_mass_ray`, with `Delta_ud_overlap = ray_modulus / 5` and `eta_Q_centered = -((1 - x2^2) / 27) * ray_modulus`; that is a real structural reduction on the current sheet. "
    "The strongest current exact-mean specialization gives "
    "`alpha_u = 1.0007763698011345`, `alpha_d = 1.008463281557513`, "
    "`Delta_ud_overlap = 0.14049991320632976`, `eta_Q_centered = -0.018104730181494357`, "
    "with RMS log-mass error about `5.21e-05`, but those values remain compare-derived rather than OPH-emitted. "
    "The retained same-family numerical point on the emitted ray is sample-only: `t1_sample = ray_modulus = 0.6695617711471163`, giving `Delta_ud_overlap = 0.13391235422942327` and `eta_Q_centered = -0.018155152181872827`. "
    "But the current D12 sheet is now an explicit no-go for the physical CKM shell on the emitted same-label left-handed solver surface. The selector is already fixed to `sigma_ref`, but the current theorem-grade rigidity result is only same-sheet rephasing rigidity; the corpus does not yet emit a theorem-grade CKM rigidity law under `ray_modulus` motion on `D12_ud_mass_ray`. The emitted angles still undershoot the comparison shell by factors about `29.68`, `34.02`, and `117.06`. "
    "The current local surface is no longer formally insufficient at the selector level: `oph_quark_sigma_ud_singleton_uniqueness_theorem` closes the emitted same-label left-handed orbit to the singleton `{sigma_ref}`, and `quark_relative_sheet_selector` now emits the theorem-grade selector value `sigma_ref` on the local solver surface. That is a real but negative closure: the singleton is still the current D12 reference sheet, so the physical CKM-shell mismatch survives unchanged on the selected branch. The only finite local scan now on disk is a 4001-point same-sheet `Delta_ud_overlap` scan against `reference_targets`; it remains comparison-only and cannot change the selected branch on the theorem lane. "
    "A smaller already-local finite basis orbit has also been exposed and ruled out: swapping in right or conjugate-right singular bases can move CKM moduli, but that orbit leaves the ordered same-label left-eigenframe domain and is therefore diagnostic-only. A second representative-level transport-frame diagnostic orbit can be extracted from the common-refinement line-lift data, and its gauge-invariant self-overlap `F0^dagger F1` improves the compare-only CKM shell loss from about `46.62` to `2.21`; but the recorded residual objectwise-U(1) quotient acts only by diagonal rephasings `M -> D_u^dagger M D_d` and therefore preserves all entrywise moduli `|M_ij|`. That transport-frame signal is still not an emitted physical `Sigma_ud` element, but for the smaller exact reason that no sector-attached lift/equivalence identifies the common-refinement matrix with an emitted same-label left-handed `Sigma_ud` element carrying `{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM}`. "
    "So the exact next quark object is no longer another branch selector: on the emitted local solver surface the selector is already fixed to `sigma_ref`, and the next theorem-side burden is the one-scalar value law `quark_d12_t1_value_law` on the already-emitted `D12_ud_mass_ray`, with `intrinsic_scale_law_D12` retained only as the derived wrapper, independent of target masses and independent of CKM/CP. The public quark rows therefore remain continuation-level and still live on the selected wrong-branch D12 sheet."
)
NEUTRINO_CONTINUATION_NOTE = (
    "No public value is emitted yet on the theorem lane, but a stronger exact sidecar is already on disk: the two-parameter positive-segment adapter hits the representative central splittings exactly while staying compare-only because the theorem lane still waits on `C_nu`. The active neutrino path is `derive_neutrino_scale_anchor.py -> "
    "derive_family_response_tensor.py -> derive_majorana_holonomy_lift.py -> derive_majorana_phase_pullback_metric.py -> "
    "build_forward_majorana_matrix.py -> build_forward_splittings.py -> export_forward_neutrino_closure_bundle.py`, and the "
    "new exact intrinsic branch is `derive_same_label_scalar_certificate.py -> build_intrinsic_neutrino_mass_eigenstate_bundle.py -> "
    "derive_intrinsic_neutrino_exact_eta_map.py -> derive_intrinsic_neutrino_exact_mixing_law_validation.py -> "
    "derive_neutrino_weighted_cycle_repair.py`. "
    "The old isotropic intrinsic ansatz is now explicitly ruled out on the live branch by the exact cap `max |Delta m^2| <= 8 a rho`, which is far below the atmospheric scale. "
    "The strongest branch on disk is therefore the shared-basis weighted-cycle lift: using the live same-label scalar certificate `(gap_e, overlap_sq_e)`, the live flavor cocycle invariants `gamma` and `eps`, the live selector phases, and the standard-math-fixed balanced transport-load selector on the positive segment between `chi = 1 + eps` and `1 + gamma_half`, it emits `D_nu = (chi + 1 + gamma_half) / 2` and the edge law `w_e = q_e^(1 + gamma + eps / D_nu)`. "
    "That current weighted-cycle branch lands in the PDG 2025 oscillation window at the dimensionless pattern level, with `theta12 = 34.2259 deg`, `theta23 = 49.7228 deg`, `theta13 = 8.68636 deg`, `delta = 305.581 deg`, and `J = -0.02753`; its predicted splitting ratio is `Delta m21^2 / Delta m32^2 = 0.03072111`. "
    "No hidden discrete branch remains on that weighted-cycle lane: the current attached stack is still constant on a one-parameter positive amplitude orbit. "
    "On the sharpened bridge scaffold, the paper-facing amplitude parameterization is still `B_nu = lambda_nu * q_mean^p_nu / m_star_eV`, but after factoring out the already-emitted positive proxy `P_nu = I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1` the exact remaining object is smaller: one reduced bridge-correction invariant `C_nu = B_nu / P_nu`. "
    "The sharper formulation of the remaining gap is an absolute-attachment theorem from the D10 amplitude sector to that weighted-cycle scale-free normal form; the naive direct attachment `lambda_nu = m_star_eV` is diagnostic-only and undershoots the atmospheric scale badly. "
    "The sharpest current local interface for that attachment is the Majorana overlap-defect scalar route. The realized same-label gap/defect readback beneath that route is already complete from live flavor overlap and spectral-gap certificates, so the remaining neutrino gap is no longer a missing pullback payload. The normalized overlap-defect weight section `oph_same_label_overlap_defect_weight_normalizer` is already closed from the live same-label scalar certificate, and the selector-centered common-refinement bundle descent now feeds the weighted-cycle bridge scaffolds rather than an isotropic live branch. Because the proxy `P_nu` is already internal to that emitted stack, the current corpus emits `B_nu` if and only if it emits `C_nu`, so the irreducible exact attachment object is the reduced correction scalar `C_nu` rather than a raw bridge amplitude above the normalizer. "
    "So that residual bridge law remains open on the live theorem lane, and absolute neutrino masses and absolute `Delta m^2` values are still compare-only unless an external atmospheric anchor is supplied. "
    "The sharpened no-go statement is therefore that the current weighted-cycle branch emits only the one-parameter absolute family `m_i = lambda_nu * mhat_i`, `Delta m^2_ij = lambda_nu^2 * Delta_hat_ij`; all PMNS data and dimensionless hierarchy ratios are unchanged across that family, so no unique theorem-grade absolute spectrum exists yet. The attached irreducibility theorem goes one step further: the present attached stack factors exactly through `q_e = q_mean * qbar_e` and still leaves a one-parameter positive amplitude orbit, so a qbar-only collapse law does not follow from the current corpus. The worker-validated correction is sharper still: the current selected-point scalar `I_nu^(wc) = 0.5 * sum_e qbar_e * |z_e(psi_wc) - 1|^2` is already a function of the emitted `qbar_e`, `psi_wc`, and `psi*` data, so it cannot be the missing bridge-external scalar. The exact remaining object is therefore the reduced correction invariant `C_nu`, while the raw paper-facing amplitude remains `B_nu = (I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1) * C_nu`, equivalently `A_nu / m_star_eV`. "
    "A cross-route compare-only corridor now sharpens that same bridge further: the direct `C_nu` audit yields a target-containing `C_nu` window `[0.9994296, 1.0015784]` with live compare-only value `0.9995473560`, and that window induces `B_nu in [6.69600, 6.71040]` with midpoint `6.70320` and relative half-width about `0.107%`, improving on the older three-route bridge corridor `[6.69315, 6.71259]`. A narrower shortlist-consensus window `[6.70484, 6.71259]` is also on disk, but it remains diagnostic-only because it does not contain the live compare-only target. "
    "A compare-only single-scale fit shows how sharp that remaining gap is: fitting `lambda_nu` against PDG 2025 central splittings gives `lambda_nu = 1.72390`, `m = (0.01746, 0.01948, 0.05308) eV`, `Delta m21^2 = 7.48982e-5 eV^2`, and `Delta m32^2 = 2.43801e-3 eV^2`, with residuals about `-9.02e-4 sigma` on `Delta m21^2` and `+2.84e-4 sigma` on `Delta m32^2`. "
    "A stronger exact compare-only adapter is now also on disk: moving along the already-emitted positive load segment `D_tau = (1-tau_nu) * chi + tau_nu * (1 + gamma_half)` and then rescaling by one positive `lambda_nu` hits both representative PDG 2025 central splittings exactly, with `tau_nu = 0.49717`, `lambda_nu = 1.72396`, and `m = (0.01746, 0.01948, 0.05308) eV`. That two-parameter segment adapter remains strictly non-promotable because the theorem lane still waits on the reduced bridge-correction invariant `C_nu`. "
    "The older atmospheric-only and solar-only exact one-observable adapters remain on disk as narrower special slices, and the same repair artifact still states that no single `lambda_nu` fits both central splittings exactly on the live promoted midpoint branch. "
    "The sharpest current constructive interface is still the Majorana overlap-defect scalar route, but after exact q_mean-factorization the remaining scalar is parameterized as `lambda_nu = (m_star_eV / q_mean^p_nu) * (I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1) * C_nu`, not as a closed `F_nu(qbar, I_nu^(wc))` law. The best closed constructive object beneath that bridge is the defect-weighted same-label family `q_e = sqrt(g_e * d_e)` together with its induced `mu_e = mu_nu * exp(eta_e) / mean_f(exp(eta_f))`; it is the first spectrum-moving object below the reduced correction scalar, but it still does not emit `C_nu` or `B_nu` itself. A compare-only residual-amplitude audit shows that `A_nu / m_star` is numerically much cleaner than the raw `lambda_nu / m_star` bridge factor. On the core residual scalar pool the strongest low-complexity clue is `sqrt(I_nu) * sqrt(Delta_hat_21 / Delta_hat_32) / sum_defect`. On the extended pool that also includes defect-weighted `mu_e` family scales normalized by `m_star`, no stronger low-complexity clue survives once the trivial closed normalizer `sum_qbar = 3` is removed from the search, so the extended pool does not beat the same core clue. Both remain compare-only and not theorem-grade. A tempting closed-form numerology `lambda_nu = gamma / sqrt(Delta_hat_21 / Delta_hat_32)` is explicitly ruled out as theorem-grade by the same positive-rescaling no-go, because both `gamma` and `Delta_hat_21 / Delta_hat_32` are orbit invariants on the scale-free family. "
    "A separate compare-only load-segment audit records nearby alternative selectors on the same positive segment; the geometric mean gives the smallest raw representative central-ratio error, but the promoted theorem branch is the arithmetic midpoint because it is the unique balanced/least-distortion selector on that one-dimensional affine segment. "
    "With the hard-separated compare-only atmospheric anchor `Delta m32^2 = 2.438e-3 eV^2`, the same branch gives `m = (0.01746, 0.01948, 0.05308) eV`, `Delta m21^2 = 7.48981e-5 eV^2`, and `Delta m31^2 = 2.51290e-3 eV^2`; relative to the representative PDG 2025 central values, the solar and atmospheric residuals are only about `-9.02e-4 sigma` and `+2.84e-4 sigma`. These rows therefore remain unpromoted as exact OPH mass predictions because the reduced bridge-correction invariant `C_nu` is still missing."
)
HADRON_CONTINUATION_NOTE = (
    "Rows are suppressed by default because hadrons are execution-contract-frozen on the current branch rather than paper-derived outputs. The active hadron path is `derive_lambda_msbar_descendant.py -> "
    "derive_full_unquenched_correlator.py -> derive_stable_channel_cfg_source_measure_payload.py -> "
    "derive_runtime_schedule_receipt_n_therm_and_n_sep.py -> derive_stable_channel_sequence_population.py -> "
    "derive_hadron_production_geometry_summary.py -> derive_stable_channel_sequence_evaluation.py -> "
    "derive_stable_channel_groundstate_readout.py`, and a separate diagnostic-only surrogate bridge "
    "`derive_hadron_surrogate_execution_bridge_status.py` records that the full receipt/writeback/evaluation/convergence/systematics path "
    "has been closed on a surrogate HMC/RHMC kernel. The operational barrier is also lower-friction now: `run_production_backend_writeback.py` executes the backend-export -> receipt -> dump -> payload -> evaluation -> closure-report path in one command once a real production export exists. The production geometry is explicit: 3 seeded 2+1 ensembles, 6 cfg total, naive raw gauge storage about "
    "`2.80071464105088e14` bytes for all cfg, and a backend correlator dump of `195264` float64 bytes. "
    "Public hadron rows still require one production backend export bundle on the seeded family with publication-complete manifest provenance and real `pi_iso`, `N_iso_direct`, and `N_iso_exchange` correlator arrays, followed by production continuum/volume/chiral/statistical systematics; the first local derivative after that bundle lands is the normalized production dump "
    "`backend_correlator_dump.production.json`."
)
INVENTORY: List[Dict[str, Any]] = [
    {
        "particle_id": "photon",
        "label": "gamma",
        "group": "Bosons",
        "prediction_key": "m_gamma",
        "ledger_id": "structural.massless.photon",
        "note": "Structural massless gauge sector.",
    },
    {
        "particle_id": "gluon",
        "label": "g (8 color states)",
        "group": "Bosons",
        "prediction_key": "m_gluon",
        "ledger_id": "structural.massless.gluons",
        "note": "Structural massless color gauge sector.",
    },
    {
        "particle_id": "graviton",
        "label": "graviton",
        "group": "Bosons",
        "prediction_key": "m_graviton",
        "ledger_id": "structural.massless.graviton",
        "note": "Structural massless spin-2 sector from the OPH dynamical-metric and diffeomorphism branch.",
    },
    {
        "particle_id": "w_boson",
        "label": "W",
        "group": "Bosons",
        "prediction_key": "mW_run",
        "ledger_id": "calibration.d10.electroweak",
        "note": D10_MASS_PAIR_NOTE,
    },
    {
        "particle_id": "z_boson",
        "label": "Z",
        "group": "Bosons",
        "prediction_key": "mZ_run",
        "ledger_id": "calibration.d10.electroweak",
        "note": D10_MASS_PAIR_NOTE,
    },
    {
        "particle_id": "higgs",
        "label": "H",
        "group": "Bosons",
        "prediction_key": "crit_mH_tree",
        "ledger_id": "secondary.d11.higgs_top",
        "note": D11_NOTE,
    },
    {
        "particle_id": "electron",
        "label": "e",
        "group": "Leptons",
        "prediction_key": "m_e",
        "ledger_id": "continuation.flavor.charged_leptons",
        "note": CHARGED_CONTINUATION_NOTE,
    },
    {
        "particle_id": "muon",
        "label": "mu",
        "group": "Leptons",
        "prediction_key": "m_mu",
        "ledger_id": "continuation.flavor.charged_leptons",
        "note": CHARGED_CONTINUATION_NOTE,
    },
    {
        "particle_id": "tau",
        "label": "tau",
        "group": "Leptons",
        "prediction_key": "m_tau",
        "ledger_id": "continuation.flavor.charged_leptons",
        "note": CHARGED_CONTINUATION_NOTE,
    },
    {
        "particle_id": "electron_neutrino",
        "label": "nu_e",
        "group": "Leptons",
        "prediction_key": "m_nu_e",
        "ledger_id": "continuation.neutrinos.d6_estimate",
        "note": NEUTRINO_CONTINUATION_NOTE,
    },
    {
        "particle_id": "muon_neutrino",
        "label": "nu_mu",
        "group": "Leptons",
        "prediction_key": "m_nu_mu",
        "ledger_id": "continuation.neutrinos.d6_estimate",
        "note": NEUTRINO_CONTINUATION_NOTE,
    },
    {
        "particle_id": "tau_neutrino",
        "label": "nu_tau",
        "group": "Leptons",
        "prediction_key": "m_nu_tau",
        "ledger_id": "continuation.neutrinos.d6_estimate",
        "note": NEUTRINO_CONTINUATION_NOTE,
    },
    {
        "particle_id": "up_quark",
        "label": "u",
        "group": "Quarks",
        "prediction_key": "m_u",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
    },
    {
        "particle_id": "down_quark",
        "label": "d",
        "group": "Quarks",
        "prediction_key": "m_d",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
    },
    {
        "particle_id": "strange_quark",
        "label": "s",
        "group": "Quarks",
        "prediction_key": "m_s",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
    },
    {
        "particle_id": "charm_quark",
        "label": "c",
        "group": "Quarks",
        "prediction_key": "m_c",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
    },
    {
        "particle_id": "bottom_quark",
        "label": "b",
        "group": "Quarks",
        "prediction_key": "m_b",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
    },
    {
        "particle_id": "top_quark",
        "label": "t",
        "group": "Quarks",
        "prediction_key": "crit_mt_pole",
        "ledger_id": "secondary.d11.higgs_top",
        "note": D11_NOTE,
        "extra_prediction_keys": ["m_t"],
    },
    {
        "particle_id": "proton",
        "label": "p",
        "group": "Hadrons",
        "prediction_key": "m_p",
        "ledger_id": "simulation.hadrons.current_lane",
        "note": HADRON_CONTINUATION_NOTE,
    },
    {
        "particle_id": "neutron",
        "label": "n",
        "group": "Hadrons",
        "prediction_key": "m_n",
        "ledger_id": "simulation.hadrons.current_lane",
        "note": HADRON_CONTINUATION_NOTE,
    },
    {
        "particle_id": "neutral_pion",
        "label": "pi0 proxy",
        "group": "Hadrons",
        "prediction_key": "m_pi",
        "ledger_id": "simulation.hadrons.current_lane",
        "note": HADRON_CONTINUATION_NOTE,
    },
    {
        "particle_id": "rho_770_0",
        "label": "rho(770)0 proxy",
        "group": "Hadrons",
        "prediction_key": "m_rho",
        "ledger_id": "simulation.hadrons.current_lane",
        "note": HADRON_CONTINUATION_NOTE,
    },
]


def _canonical_artifact_ref(path: pathlib.Path | str) -> str:
    path = pathlib.Path(path)
    try:
        rel = path.relative_to(ROOT)
    except ValueError:
        return str(path)
    return f"code/{rel.as_posix()}"


def _effective_hadron_profile(*, with_hadrons: bool, hadron_profile: str) -> str:
    return hadron_profile if with_hadrons else "suppressed"


def format_gev(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    if value == 0.0:
        return "0"
    abs_value = abs(value)
    if abs_value < 1.0e-9 or abs_value >= 1.0e4:
        return f"{value:.6e}"
    if abs_value < 1.0e-4:
        return f"{value:.12g}"
    if abs_value < 1.0:
        return f"{value:.10g}"
    return f"{value:.9g}"


def format_scalar(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    if value == 0.0:
        return "0"
    abs_value = abs(value)
    if abs_value < 1.0e-9 or abs_value >= 1.0e4:
        return f"{value:.6e}"
    if abs_value < 1.0e-4:
        return f"{value:.12g}"
    if abs_value < 1.0:
        return f"{value:.8g}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def format_observable_value(value: Optional[float], unit: str) -> str:
    if value is None:
        return "n/a"
    return f"{format_scalar(value)} {unit}".rstrip()


def format_observable_reference(value: float, unit: str, *, err_plus: float, err_minus: float) -> str:
    if abs(err_plus - err_minus) <= max(abs(err_plus), abs(err_minus), 1.0) * 1.0e-12:
        return f"{format_scalar(value)} +- {format_scalar(err_plus)} {unit}".rstrip()
    return f"{format_scalar(value)} +{format_scalar(err_plus)} -{format_scalar(err_minus)} {unit}".rstrip()


def format_observable_delta(pred_value: float, reference_value: float, unit: str) -> str:
    delta = pred_value - reference_value
    rel = None if reference_value == 0 else delta / reference_value
    display = format_observable_value(delta, unit)
    if rel is None:
        return display
    return f"{display} ({rel:+.3e})"


def format_reference(entry: Dict[str, Any]) -> str:
    if entry.get("value_gev") is not None:
        value_gev = float(entry["value_gev"])
        if entry.get("reference_kind") == "upper_limit":
            return f"<{format_gev(value_gev)} GeV"
        err_plus = entry.get("error_plus_gev")
        err_minus = entry.get("error_minus_gev")
        if err_plus is not None and err_minus is not None:
            err_plus = float(err_plus)
            err_minus = float(err_minus)
            if abs(err_plus - err_minus) <= max(err_plus, err_minus, 1.0) * 1.0e-12:
                return f"{format_gev(value_gev)} +- {format_gev(err_plus)} GeV"
            return f"{format_gev(value_gev)} +{format_gev(err_plus)} -{format_gev(err_minus)} GeV"
        return f"{format_gev(value_gev)} GeV"
    display = entry.get("display")
    if display:
        return str(display)
    return "n/a"


def format_delta(pred_value: Optional[float], reference: Dict[str, Any]) -> str:
    ref_kind = reference.get("reference_kind")
    ref_value = reference.get("value_gev")
    if pred_value is None:
        return "n/a"
    if ref_kind == "upper_limit" and ref_value is not None:
        return "within limit" if pred_value <= ref_value else f"+{format_gev(pred_value - ref_value)} above limit"
    if ref_kind != "value" or ref_value is None:
        return "n/a"
    delta = pred_value - float(ref_value)
    rel = None if ref_value == 0 else delta / float(ref_value)
    if rel is None:
        return format_gev(delta)
    return f"{format_gev(delta)} ({rel:+.3e})"


def build_note(
    row_spec: Dict[str, Any],
    reference: Dict[str, Any],
    prediction: Dict[str, Any],
    ledger_entry: Dict[str, Any],
) -> str:
    pieces: List[str] = [row_spec["note"]]
    if row_spec["particle_id"] == "top_quark" and prediction.get("d12_m_t_sidecar_gev") is not None:
        pieces.append(f"D12 sidecar value: {format_gev(float(prediction['d12_m_t_sidecar_gev']))} GeV.")
    ref_note = reference.get("comment")
    if ref_note:
        ref_note_text = str(ref_note).strip()
        if ref_note_text.lower().startswith("of "):
            ref_note_text = "Reference is " + ref_note_text[3:]
        pieces.append(ref_note_text)
    if row_spec["particle_id"] in {"up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark"}:
        pieces.append("PDG quark references are running masses, not direct free-particle pole masses.")
    if row_spec["particle_id"] in {"electron_neutrino", "muon_neutrino", "tau_neutrino"}:
        pieces.append("Comparison remains partly qualitative because the current weighted-cycle oscillation branch still waits on the reduced bridge-correction invariant `C_nu` above the emitted proxy.")
    if row_spec["group"] == "Hadrons":
        pieces.append("Use this only when explicitly debugging the hadron pipeline.")
    return " ".join(piece for piece in pieces if piece)


def load_reference_entries(path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload["entries"]


def load_ledger_entries(path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {entry["id"]: entry for entry in payload["entries"]}


def _d10_public_mass_pair_allowed(readout: Dict[str, Any]) -> bool:
    mass_pair = dict(readout.get("mass_pair_predictive_candidate", {}))
    return bool(readout.get("public_surface_candidate_allowed", False)) and all(
        key in mass_pair for key in ("MW_pole", "MZ_pole")
    )


def _d11_public_seed_allowed(seed: Dict[str, Any]) -> bool:
    mass_readout = dict(seed.get("mass_readout", {}))
    return bool(seed.get("public_surface_candidate_allowed", False)) and all(
        key in mass_readout for key in ("mH_gev", "mt_pole_gev")
    )


def _quark_public_forward_allowed(forward: Dict[str, Any], mean_split: Dict[str, Any]) -> bool:
    return (
        bool(forward.get("public_surface_candidate_allowed", False))
        and forward.get("source_mode") == "factorized_descent"
        and mean_split.get("active_candidate") != "current_family_exact_witness"
    )


def _charged_public_candidate_allowed(forward: Dict[str, Any]) -> bool:
    return bool(forward.get("public_surface_candidate_allowed", False))


def _neutrino_public_candidate_allowed(bundle: Dict[str, Any]) -> bool:
    return bool(bundle.get("public_surface_candidate_allowed", False))


def _neutrino_repaired_branch_waiting_absolute_scale(blockers: Dict[str, Any]) -> bool:
    status = dict(blockers.get("live_continuation_branch_status", {}))
    exact_blockers = list(blockers.get("exact_blockers") or [])
    blocker_names = {item.get("name") for item in exact_blockers}
    blocker_kinds = {item.get("kind") for item in exact_blockers}
    repaired_status = status.get("status") in {
        "physically_repaired_up_to_one_positive_scale",
        "physically_repaired_up_to_one_reduced_bridge_correction_invariant",
    }
    repaired_blocker_surface = not exact_blockers or (
        "one_positive_neutrino_bridge_correction_invariant" in blocker_names
        and "reduced_bridge_correction_invariant" in blocker_kinds
    )
    return (
        repaired_status
        and repaired_blocker_surface
        and bool(status.get("same_label_scalar_certificate_present"))
        and bool(status.get("shared_charged_left_basis_present"))
        and bool(status.get("repair_artifact_present"))
    )


def build_surface_state(*, with_hadrons: bool) -> Dict[str, Any]:
    d10_active = False
    d11_active = False
    charged_active = False
    neutrino_active = False
    neutrino_repaired_branch = False
    quark_active = False

    if D10_SOURCE_TRANSPORT_READOUT.exists():
        readout = json.loads(D10_SOURCE_TRANSPORT_READOUT.read_text(encoding="utf-8"))
        d10_active = _d10_public_mass_pair_allowed(readout)

    if D11_FORWARD_SEED.exists():
        seed = json.loads(D11_FORWARD_SEED.read_text(encoding="utf-8"))
        d11_active = _d11_public_seed_allowed(seed)

    if FORWARD_CHARGED_LEPTONS.exists():
        charged = json.loads(FORWARD_CHARGED_LEPTONS.read_text(encoding="utf-8"))
        charged_active = _charged_public_candidate_allowed(charged)

    if FORWARD_NEUTRINO_BUNDLE.exists():
        bundle = json.loads(FORWARD_NEUTRINO_BUNDLE.read_text(encoding="utf-8"))
        neutrino_active = _neutrino_public_candidate_allowed(bundle)
    if NEUTRINO_EXACT_BLOCKERS.exists():
        blockers = json.loads(NEUTRINO_EXACT_BLOCKERS.read_text(encoding="utf-8"))
        neutrino_repaired_branch = _neutrino_repaired_branch_waiting_absolute_scale(blockers)

    if FORWARD_YUKAWAS.exists() and QUARK_SECTOR_MEAN_SPLIT.exists():
        forward = json.loads(FORWARD_YUKAWAS.read_text(encoding="utf-8"))
        mean_split = json.loads(QUARK_SECTOR_MEAN_SPLIT.read_text(encoding="utf-8"))
        quark_active = _quark_public_forward_allowed(forward, mean_split)

    return {
        "public_surface_kind": PUBLIC_SURFACE_KIND,
        "surface_policy": "local_candidate_or_gap_only",
        "active_local_public_candidates": {
            "d10_mass_pair": d10_active,
            "d11_forward_seed": d11_active,
            "charged_local_candidate": charged_active,
            "neutrino_local_candidate": neutrino_active,
            "neutrino_repaired_branch": neutrino_repaired_branch,
            "quark_forward_candidate": quark_active,
            "hadrons_enabled": with_hadrons,
        },
    }


def apply_local_candidate_overrides(prediction: Dict[str, Any]) -> Dict[str, Any]:
    updated = dict(prediction)
    if prediction.get("m_t") is not None and updated.get("d12_m_t_sidecar_gev") is None:
        updated["d12_m_t_sidecar_gev"] = float(prediction["m_t"])

    updated.setdefault("m_gamma", 0.0)
    updated.setdefault("m_gluon", 0.0)
    updated.setdefault("m_graviton", 0.0)

    if D10_SOURCE_TRANSPORT_READOUT.exists():
        readout = json.loads(D10_SOURCE_TRANSPORT_READOUT.read_text(encoding="utf-8"))
        mass_pair = dict(readout.get("mass_pair_predictive_candidate", {}))
        if _d10_public_mass_pair_allowed(readout):
            updated.update(
                {
                    "mW_run": float(mass_pair["MW_pole"]),
                    "mZ_run": float(mass_pair["MZ_pole"]),
                }
            )

    if D11_FORWARD_SEED.exists():
        seed = json.loads(D11_FORWARD_SEED.read_text(encoding="utf-8"))
        if _d11_public_seed_allowed(seed):
            mass_readout = dict(seed.get("mass_readout", {}))
            updated.update(
                {
                    "crit_mH_tree": float(mass_readout["mH_gev"]),
                    "crit_mt_pole": float(mass_readout["mt_pole_gev"]),
                }
            )

    if FORWARD_YUKAWAS.exists() and QUARK_SECTOR_MEAN_SPLIT.exists():
        forward = json.loads(FORWARD_YUKAWAS.read_text(encoding="utf-8"))
        mean_split = json.loads(QUARK_SECTOR_MEAN_SPLIT.read_text(encoding="utf-8"))
        if _quark_public_forward_allowed(forward, mean_split):
            u = [float(x) for x in forward["singular_values_u"]]
            d = [float(x) for x in forward["singular_values_d"]]
            updated.update(
                {
                    "m_u": u[0],
                    "m_c": u[1],
                    "m_d": d[0],
                    "m_s": d[1],
                    "m_b": d[2],
                }
            )

    return updated


def build_neutrino_oscillation_comparison_rows(surface_state: Dict[str, Any]) -> List[Dict[str, Any]]:
    active = dict(surface_state["active_local_public_candidates"])
    if not active.get("neutrino_repaired_branch"):
        return []
    if not NEUTRINO_WEIGHTED_CYCLE_REPAIR.exists():
        return []

    repair = json.loads(NEUTRINO_WEIGHTED_CYCLE_REPAIR.read_text(encoding="utf-8"))
    pmns = dict(repair.get("pmns_observables", {}))
    anchored = dict(repair.get("compare_only_atmospheric_anchor", {}))
    two_parameter_adapter = (
        json.loads(NEUTRINO_TWO_PARAMETER_EXACT_ADAPTER.read_text(encoding="utf-8"))
        if NEUTRINO_TWO_PARAMETER_EXACT_ADAPTER.exists()
        else None
    )
    ratio_value = repair.get("dimensionless_ratio_dm21_over_dm32")
    ratio_reference = (
        NEUTRINO_PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"] / NEUTRINO_PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"]
    )

    def _row(
        *,
        observable_id: str,
        observable: str,
        status: str,
        prediction_value: float,
        reference_value: float,
        err_plus: float,
        err_minus: float,
        unit: str,
        note: str,
    ) -> Dict[str, Any]:
        return {
            "observable_id": observable_id,
            "observable": observable,
            "status": status,
            "prediction_value": float(prediction_value),
            "prediction_display": format_observable_value(float(prediction_value), unit),
            "reference_value": float(reference_value),
            "reference_display": format_observable_reference(
                float(reference_value),
                unit,
                err_plus=float(err_plus),
                err_minus=float(err_minus),
            ),
            "delta_display": format_observable_delta(float(prediction_value), float(reference_value), unit),
            "note": note,
            "reference_source_url": NEUTRINO_OSCILLATION_SOURCE_URL,
            "reference_label": NEUTRINO_OSCILLATION_REFERENCE_LABEL,
            "unit": unit,
        }

    rows = [
        _row(
            observable_id="theta12_deg",
            observable="theta12",
            status="weighted_cycle_dimensionless",
            prediction_value=float(pmns["theta12_deg"]),
            reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["theta12_deg"],
            err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["theta12_deg"]["plus"],
            err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["theta12_deg"]["minus"],
            unit="deg",
            note="Direct PMNS angle from the current weighted-cycle branch; no absolute mass anchor is used.",
        ),
        _row(
            observable_id="theta23_deg",
            observable="theta23",
            status="weighted_cycle_dimensionless",
            prediction_value=float(pmns["theta23_deg"]),
            reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["theta23_deg"],
            err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["theta23_deg"]["plus"],
            err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["theta23_deg"]["minus"],
            unit="deg",
            note="Direct PMNS angle from the current weighted-cycle branch; this lands inside the PDG 3sigma NO window near its upper edge.",
        ),
        _row(
            observable_id="theta13_deg",
            observable="theta13",
            status="weighted_cycle_dimensionless",
            prediction_value=float(pmns["theta13_deg"]),
            reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["theta13_deg"],
            err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["theta13_deg"]["plus"],
            err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["theta13_deg"]["minus"],
            unit="deg",
            note="Direct PMNS angle from the current weighted-cycle branch; no absolute mass anchor is used.",
        ),
        _row(
            observable_id="delta_cp_deg",
            observable="delta_CP",
            status="weighted_cycle_dimensionless",
            prediction_value=float(pmns["delta_deg"]),
            reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_deg"],
            err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_deg"]["plus"],
            err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_deg"]["minus"],
            unit="deg",
            note="Direct PMNS phase from the current weighted-cycle branch; inside the PDG 3sigma NO window but displaced from the NO best fit.",
        ),
    ]

    if ratio_value is not None:
        rows.append(
            {
                "observable_id": "delta_m21_sq_over_delta_m32_sq",
                "observable": "Delta m21^2 / Delta m32^2",
                "status": "weighted_cycle_dimensionless",
                "prediction_value": float(ratio_value),
                "prediction_display": format_scalar(float(ratio_value)),
                "reference_value": float(ratio_reference),
                "reference_display": format_scalar(float(ratio_reference)),
                "delta_display": format_observable_delta(float(ratio_value), float(ratio_reference), ""),
                "note": "Dimensionless hierarchy ratio from the current weighted-cycle branch; this is the cleanest split comparison because it does not depend on the missing absolute normalization scalar.",
                "reference_source_url": NEUTRINO_OSCILLATION_SOURCE_URL,
                "reference_label": NEUTRINO_OSCILLATION_REFERENCE_LABEL,
                "unit": "",
            }
        )

    if two_parameter_adapter:
        exact_dm2 = dict(two_parameter_adapter["exact_outputs"]["delta_m_sq_eV2"])
        rows.extend(
            [
                _row(
                    observable_id="delta_m21_sq_eV2",
                    observable="Delta m21^2",
                    status="compare_only",
                    prediction_value=float(exact_dm2["21"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="Exact compare-only central solar splitting from the two-parameter positive-segment neutrino adapter; this remains non-promotable because the live theorem lane still waits on the reduced bridge-correction invariant C_nu.",
                ),
                _row(
                    observable_id="delta_m32_sq_eV2",
                    observable="Delta m32^2",
                    status="compare_only",
                    prediction_value=float(exact_dm2["32"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="Exact compare-only central atmospheric splitting from the same two-parameter positive-segment neutrino adapter; the older one-parameter atmospheric anchor remains on disk only as a narrower diagnostic slice.",
                ),
            ]
        )
    elif "delta_m21_sq_eV2" in anchored and "delta_m32_sq_eV2" in anchored:
        rows.extend(
            [
                _row(
                    observable_id="delta_m21_sq_eV2",
                    observable="Delta m21^2",
                    status="compare_only",
                    prediction_value=float(anchored["delta_m21_sq_eV2"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="Absolute solar splitting after compare-only anchoring with the atmospheric Delta m32^2 input; this is not yet a promoted theorem-grade OPH output.",
                ),
                _row(
                    observable_id="delta_m32_sq_eV2",
                    observable="Delta m32^2",
                    status="compare_only_anchor",
                    prediction_value=float(anchored["delta_m32_sq_eV2"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="This is the external atmospheric anchor used to put the repaired dimensionless branch on an eV scale, so it is shown only as compare-only context, not as an independent prediction.",
                ),
            ]
        )

    return rows


def prediction_surface_for_row(row_spec: Dict[str, Any], surface_state: Dict[str, Any], *, with_hadrons: bool) -> str:
    active = dict(surface_state["active_local_public_candidates"])
    particle_id = row_spec["particle_id"]
    if particle_id in {"photon", "gluon", "graviton"}:
        return "particles_structural_massless"
    if particle_id in {"w_boson", "z_boson"} and active.get("d10_mass_pair"):
        return "local_d10_public_mass_pair_candidate"
    if particle_id in {"higgs", "top_quark"} and active.get("d11_forward_seed"):
        return "local_d11_forward_seed_candidate"
    if particle_id in {"electron", "muon", "tau"} and active.get("charged_local_candidate"):
        return "local_charged_public_candidate"
    if particle_id in {"electron_neutrino", "muon_neutrino", "tau_neutrino"} and active.get("neutrino_local_candidate"):
        return "local_neutrino_public_candidate"
    if particle_id in {"electron_neutrino", "muon_neutrino", "tau_neutrino"} and active.get("neutrino_repaired_branch"):
        return "local_neutrino_repaired_branch_waiting_absolute_scale"
    if particle_id in {
        "up_quark",
        "down_quark",
        "strange_quark",
        "charm_quark",
        "bottom_quark",
    } and active.get("quark_forward_candidate"):
        return "local_quark_public_forward_candidate"
    if row_spec["group"] == "Hadrons" and not with_hadrons:
        return "suppressed"
    return "particles_gap"


def build_rows(
    prediction: Dict[str, Any],
    reference_entries: Dict[str, Dict[str, Any]],
    ledger_entries: Dict[str, Dict[str, Any]],
    *,
    with_hadrons: bool,
    surface_state: Dict[str, Any],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for row_spec in INVENTORY:
        if row_spec["group"] == "Hadrons" and not with_hadrons:
            continue
        reference = reference_entries[row_spec["particle_id"]]
        ledger_entry = ledger_entries[row_spec["ledger_id"]]
        pred_value = prediction.get(row_spec["prediction_key"])
        pred_value = None if pred_value is None else float(pred_value)
        rows.append(
            {
                "particle_id": row_spec["particle_id"],
                "particle": row_spec["label"],
                "group": row_spec["group"],
                "status": ledger_entry["tier"],
                "status_label": ledger_entry["label"],
                "prediction_key": row_spec["prediction_key"],
                "prediction_value_gev": pred_value,
                "prediction_display_gev": format_gev(pred_value),
                "reference_kind": reference["reference_kind"],
                "reference_display": format_reference(reference),
                "reference_value_gev": reference.get("value_gev"),
                "delta_display": format_delta(pred_value, reference),
                "prediction_surface": prediction_surface_for_row(row_spec, surface_state, with_hadrons=with_hadrons),
                "note": build_note(row_spec, reference, prediction, ledger_entry),
                "reference_source_url": reference["source"]["url"],
            }
        )
    return rows


def build_premise_boundaries() -> Dict[str, Any]:
    uv_boundary = json.loads(UV_BW_SCAFFOLD.read_text(encoding="utf-8"))["public_status_boundary"]
    witness_chain = []
    for item in uv_boundary.get("local_intermediate_witness_chain", []):
        normalized_item = dict(item)
        artifact = normalized_item.get("artifact")
        if artifact:
            normalized_item["artifact"] = _canonical_artifact_ref(artifact)
        witness_chain.append(normalized_item)
    if witness_chain:
        uv_boundary["local_intermediate_witness_chain"] = witness_chain
    uv_boundary["canonical_scaffold_artifacts"] = [
        _canonical_artifact_ref(UV_BW_PRELIMIT_SYSTEM),
        _canonical_artifact_ref(UV_BW_FIXED_LOCAL_COLLAR_DATUM),
        _canonical_artifact_ref(UV_BW_CARRIED_SCHEDULE),
        _canonical_artifact_ref(UV_BW_CAP_PAIR_SCAFFOLD),
        _canonical_artifact_ref(UV_BW_RIGIDITY_SCAFFOLD),
    ]
    uv_boundary["prelimit_system_artifact"] = _canonical_artifact_ref(UV_BW_PRELIMIT_SYSTEM)
    uv_boundary["remaining_missing_emitted_witness_artifact"] = _canonical_artifact_ref(UV_BW_CARRIED_SCHEDULE)
    uv_boundary["smaller_remaining_raw_datum_artifact"] = _canonical_artifact_ref(UV_BW_FIXED_LOCAL_COLLAR_DATUM)
    if NEUTRINO_LAMBDA_BRIDGE_CANDIDATE.exists():
        uv_boundary["neutrino_local_bridge_candidate_context"] = _canonical_artifact_ref(
            NEUTRINO_LAMBDA_BRIDGE_CANDIDATE
        )
    return {
        "uv_bw_internalization": uv_boundary,
    }


def render_markdown(
    *,
    rows: List[Dict[str, Any]],
    comparison_rows: List[Dict[str, Any]],
    generated_utc: str,
    P: float,
    log_dim_H: float,
    loops: int,
    with_hadrons: bool,
    hadron_profile: str,
    reference_payload: Dict[str, Any],
    surface_state: Dict[str, Any],
    premise_boundaries: Dict[str, Any],
) -> str:
    groups_present = [group for group in GROUP_ORDER if any(item["group"] == group for item in rows)]
    hadron_profile_display = _effective_hadron_profile(with_hadrons=with_hadrons, hadron_profile=hadron_profile)
    lines: List[str] = [
        "# Particle Results Status",
        "",
        f"Generated: `{generated_utc}`",
        "",
        f"Inputs: `P={P}` | `log_dim_H={log_dim_H}` | `loops={loops}` | `with_hadrons={with_hadrons}` | `hadron_profile={hadron_profile_display}`",
        "",
        f"Public Surface: `{surface_state['public_surface_kind']}`",
        "",
        f"Surface Policy: `{surface_state['surface_policy']}`",
        "",
        "Active Local Public Candidates: "
        f"`D10={surface_state['active_local_public_candidates']['d10_mass_pair']}` | "
        f"`D11={surface_state['active_local_public_candidates']['d11_forward_seed']}` | "
        f"`charged={surface_state['active_local_public_candidates']['charged_local_candidate']}` | "
        f"`neutrinos={surface_state['active_local_public_candidates']['neutrino_local_candidate']}` | "
        f"`neutrino_repaired={surface_state['active_local_public_candidates']['neutrino_repaired_branch']}` | "
        f"`quarks={surface_state['active_local_public_candidates']['quark_forward_candidate']}` | "
        f"`hadrons_enabled={surface_state['active_local_public_candidates']['hadrons_enabled']}`",
        "",
        "This table is a `/particles`-native audit surface. If a sector has no live local public candidate yet, the value is reported as `n/a`; legacy fallback predictors are not used.",
        "",
        "Hadron rows are intentionally suppressed by default because the hadron lane is execution-contract-frozen: promotable rows require a real production backend export bundle plus production systematics, not just further symbolic derivation. Re-enable them only for explicit hadron debugging with `--with-hadrons`.",
        "",
        f"Measured/reference values are pinned from the official {reference_payload['source']['label']} {reference_payload['source']['edition']} machine-readable surface where available, with explicit manual structural-context entries for non-PDG rows such as gluons, graviton, and flavor neutrinos: {reference_payload['source']['api_info_url']}.",
        "",
    ]

    uv_boundary = premise_boundaries.get("uv_bw_internalization")
    if uv_boundary:
        lines.extend(
            [
                "## Premise Boundaries",
                "",
                f"- `uv_bw_internalization`: `{uv_boundary['status']}`",
                f"- First exact object: `{uv_boundary['remaining_object']}`",
                f"- Second exact object: `{uv_boundary['follow_on_object']}`",
                f"- Remaining split: `{uv_boundary['remaining_objects'][0]}` + `{uv_boundary['remaining_objects'][1]}`",
                f"- Internalized scope: {uv_boundary['current_internalized_scope']}",
                f"- Why still open: {uv_boundary['reason_current_corpus_fails']}",
                f"- First theorem object: {uv_boundary['statement']}",
                f"- Second theorem object: {uv_boundary['follow_on_statement']}",
                f"- Candidate extension status: `{uv_boundary['candidate_extension_status']}`",
                f"- Filled witnesses already on disk: `{uv_boundary['filled_contract_witnesses'][0]}`, `{uv_boundary['filled_contract_witnesses'][1]}`, `{uv_boundary['filled_contract_witnesses'][2]}`",
                f"- Remaining emitted witness: `{uv_boundary['remaining_missing_emitted_witness']}`",
                f"- Remaining witness formula: `{uv_boundary['remaining_missing_emitted_witness_formula']}`",
                f"- Smaller raw datum beneath that witness: `{uv_boundary['smaller_remaining_raw_datum']}`",
                f"- Smaller raw datum artifact: `{uv_boundary.get('smaller_remaining_raw_datum_artifact', str(UV_BW_FIXED_LOCAL_COLLAR_DATUM))}`",
                f"- Smallest exact blocker: `{uv_boundary.get('smallest_exact_blocker', 'n/a')}`",
                f"- Smallest exact blocker formula: `{uv_boundary.get('smallest_exact_blocker_formula', 'n/a')}`",
                f"- Candidate extension route: {uv_boundary['candidate_extension_route']}",
                f"- Candidate extension target: `{uv_boundary['candidate_extension_target']}`",
                "",
            ]
        )

    for group in groups_present:
        lines.extend(
            [
                f"## {group}",
                "",
                "| Particle | Status | OPH value (GeV) | Measured / reference | Delta | Note |",
                "| --- | --- | ---: | --- | --- | --- |",
            ]
        )
        for row in [item for item in rows if item["group"] == group]:
            lines.append(
                f"| {row['particle']} | {row['status']} | {row['prediction_display_gev']} | {row['reference_display']} | {row['delta_display']} | {row['note']} |"
            )
        lines.append("")
        if group == "Leptons" and comparison_rows:
            lines.extend(
                [
                    "## Neutrino Oscillation Comparison",
                    "",
                    "| Observable | Status | OPH value | PDG 2025 NO reference | Delta | Note |",
                    "| --- | --- | --- | --- | --- | --- |",
                ]
            )
            for row in comparison_rows:
                lines.append(
                    f"| {row['observable']} | {row['status']} | {row['prediction_display']} | {row['reference_display']} | {row['delta_display']} | {row['note']} |"
                )
            lines.append("")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current `/particles` results status table.")
    parser.add_argument("--P", type=float, default=P_DEFAULT, help="Metadata-only pixel constant.")
    parser.add_argument("--log-dim-H", type=float, default=LOG_DIM_H_DEFAULT, help="Metadata-only screen-capacity constant.")
    parser.add_argument("--loops", type=int, default=4, choices=[1, 2, 3, 4], help="Metadata-only loop-order tag.")
    parser.add_argument("--with-hadrons", dest="with_hadrons", action="store_true", help="Include the current debug hadron lane.")
    parser.add_argument("--no-hadrons", dest="with_hadrons", action="store_false", help="Skip hadron computation and suppress hadron rows.")
    parser.add_argument("--hadron-profile", default="serious", choices=["demo", "quick", "serious"], help="Hadron profile for optional comparison rows.")
    parser.add_argument("--reference-json", default=str(REFERENCE_JSON), help="Pinned reference JSON path.")
    parser.add_argument("--ledger-yaml", default=str(LEDGER_YAML), help="Claim ledger path.")
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT), help="Markdown output path.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT), help="JSON output path.")
    parser.add_argument("--forward-out", default=str(DEFAULT_FORWARD_OUT), help="Forward artifact output path.")
    parser.set_defaults(with_hadrons=False)
    args = parser.parse_args()

    with_hadrons = bool(args.with_hadrons)
    surface_state = build_surface_state(with_hadrons=with_hadrons)
    premise_boundaries = build_premise_boundaries()
    reference_payload = json.loads(pathlib.Path(args.reference_json).read_text(encoding="utf-8"))
    reference_entries = reference_payload["entries"]
    ledger_entries = load_ledger_entries(pathlib.Path(args.ledger_yaml))
    prediction = apply_local_candidate_overrides({})
    rows = build_rows(
        prediction,
        reference_entries,
        ledger_entries,
        with_hadrons=with_hadrons,
        surface_state=surface_state,
    )
    comparison_rows = build_neutrino_oscillation_comparison_rows(surface_state)
    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    effective_hadron_profile = _effective_hadron_profile(
        with_hadrons=with_hadrons,
        hadron_profile=str(args.hadron_profile),
    )

    markdown = render_markdown(
        rows=rows,
        comparison_rows=comparison_rows,
        generated_utc=generated_utc,
        P=float(args.P),
        log_dim_H=float(args.log_dim_H),
        loops=int(args.loops),
        with_hadrons=with_hadrons,
        hadron_profile=str(args.hadron_profile),
        reference_payload=reference_payload,
        surface_state=surface_state,
        premise_boundaries=premise_boundaries,
    )

    markdown_out = pathlib.Path(args.markdown_out)
    markdown_out.write_text(markdown + "\n", encoding="utf-8")

    forward_out = pathlib.Path(args.forward_out)
    forward_out.parent.mkdir(parents=True, exist_ok=True)
    forward_payload = {
        "artifact": "oph_status_table_forward_current",
        "generated_utc": generated_utc,
        "status": "particles_native_candidate_or_gap_surface",
        "inputs": {
            "P": float(args.P),
            "log_dim_H": float(args.log_dim_H),
            "loops": int(args.loops),
            "with_hadrons": with_hadrons,
            "hadron_profile": effective_hadron_profile,
        },
        "surface_state": surface_state,
        "premise_boundaries": premise_boundaries,
        "rows": rows,
        "comparison_rows": comparison_rows,
    }
    forward_out.write_text(json.dumps(forward_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    json_out = pathlib.Path(args.json_out)
    json_out.write_text(
        json.dumps(
            {
                "generated_utc": generated_utc,
                "inputs": {
                    "P": float(args.P),
                    "loops": int(args.loops),
                    "log_dim_H": float(args.log_dim_H),
                    "with_hadrons": with_hadrons,
                    "hadron_profile": effective_hadron_profile,
                },
                "surface_state": surface_state,
                "premise_boundaries": premise_boundaries,
                "reference_source": reference_payload["source"],
                "rows": rows,
                "comparison_rows": comparison_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"saved: {markdown_out}")
    print(f"saved: {json_out}")
    print(f"saved: {forward_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
