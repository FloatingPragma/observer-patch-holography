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
QUARK_SHARED_ABSOLUTE_NORM_BINDING = ROOT / "particles" / "runs" / "flavor" / "quark_shared_absolute_norm_binding.json"
QUARK_RELATIVE_SHEET_SELECTOR = ROOT / "particles" / "runs" / "flavor" / "quark_relative_sheet_selector.json"
QUARK_PUBLIC_SIGMA_DESCENT = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
QUARK_PUBLIC_EXACT_YUKAWA_THEOREM = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
D10_SOURCE_TRANSPORT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
D11_FORWARD_SEED = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
FORWARD_CHARGED_LEPTONS = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
FORWARD_NEUTRINO_BUNDLE = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
NEUTRINO_BRIDGE_RIGIDITY_THEOREM = ROOT / "particles" / "runs" / "neutrino" / "neutrino_bridge_rigidity_theorem.json"
NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_theorem.json"
NEUTRINO_EXACT_BLOCKERS = ROOT / "particles" / "runs" / "neutrino" / "exact_blocking_items.json"
NEUTRINO_WEIGHTED_CYCLE_REPAIR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
NEUTRINO_TWO_PARAMETER_EXACT_ADAPTER = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
NEUTRINO_EXACT_ADAPTER_BRIDGE_COORDINATE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_exact_adapter_bridge_coordinate.json"
NEUTRINO_LAMBDA_BRIDGE_CANDIDATE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"
NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM = (
    ROOT / "particles" / "runs" / "neutrino" / "neutrino_physical_majorana_phase_theorem.json"
)
QUARK_D12_INTERNAL_BACKREAD_CLOSURE = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_continuation_closure.json"
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
    "The D10 mass-side surface fixes the `W/Z` pair from that source trunk, and the freeze-once coherent pair serves as compare-only validation on the same family. "
    "The electromagnetic row is not taken from the compact hypercharge slice. Its source anchor is `a0 = alpha_em^-1(m_Z^2) = 128.30576920234813` on the same running family, and the physical fine-structure constant is read on the Ward-projected `U(1)_Q` transport family with Thomson endpoint `alpha^-1(0) = 137.035999177`. "
    "`EWTransportKernel_D10`, `EWTransportReadoutCoherence_D10`, and `EWScalarProvenanceEquality_D10` supply the shared-kernel and shared-provenance gate for that readout. "
    "The compact anti-diagonal carrier slice and the source-only target-free emitter artifacts sit on disk as carrier or compare surfaces; they do not define the public electromagnetic theorem."
)
D11_NOTE = (
    "Derived from `derive_d11_forward_seed.py -> derive_d11_forward_seed_promotion_certificate.py`, which propagates the D10 gauge core into the compact D11 forward seed, closes the emitted one-scalar forward branch on the fixed ray, and reads out the D11 mass row from the Jacobian surface. "
    "A separate exact-hit sidecar is on disk as `oph_d11_reference_exact_adapter`: it solves the linear D11 Jacobian against the canonical Higgs/top reference pair and therefore hits those references exactly, but only as a compare-only inverse slice. "
    "The live public D11 rows remain the reference-free forward-seed outputs, not the inverse adapter."
)
_NEUTRINO_EXACT_BRIDGE_COORDINATE = (
    json.loads(NEUTRINO_EXACT_ADAPTER_BRIDGE_COORDINATE.read_text(encoding="utf-8"))
    if NEUTRINO_EXACT_ADAPTER_BRIDGE_COORDINATE.exists()
    else None
)
_NEUTRINO_BRIDGE_RIGIDITY = (
    json.loads(NEUTRINO_BRIDGE_RIGIDITY_THEOREM.read_text(encoding="utf-8"))
    if NEUTRINO_BRIDGE_RIGIDITY_THEOREM.exists()
    else None
)
_NEUTRINO_ABSOLUTE_ATTACHMENT = (
    json.loads(NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.read_text(encoding="utf-8"))
    if NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.exists()
    else None
)
_QUARK_D12_INTERNAL_BACKREAD = (
    json.loads(QUARK_D12_INTERNAL_BACKREAD_CLOSURE.read_text(encoding="utf-8"))
    if QUARK_D12_INTERNAL_BACKREAD_CLOSURE.exists()
    else None
)
_QUARK_SECTOR_MEAN_SPLIT = (
    json.loads(QUARK_SECTOR_MEAN_SPLIT.read_text(encoding="utf-8"))
    if QUARK_SECTOR_MEAN_SPLIT.exists()
    else None
)
_QUARK_SHARED_NORM_BINDING = (
    json.loads(QUARK_SHARED_ABSOLUTE_NORM_BINDING.read_text(encoding="utf-8"))
    if QUARK_SHARED_ABSOLUTE_NORM_BINDING.exists()
    else None
)
_QUARK_RELATIVE_SELECTOR = (
    json.loads(QUARK_RELATIVE_SHEET_SELECTOR.read_text(encoding="utf-8"))
    if QUARK_RELATIVE_SHEET_SELECTOR.exists()
    else None
)
_QUARK_PUBLIC_SIGMA_DESCENT = (
    json.loads(QUARK_PUBLIC_SIGMA_DESCENT.read_text(encoding="utf-8"))
    if QUARK_PUBLIC_SIGMA_DESCENT.exists()
    else None
)
_QUARK_PUBLIC_EXACT_YUKAWA = (
    json.loads(QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.read_text(encoding="utf-8"))
    if QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.exists()
    else None
)
_QUARK_D12_INTERNAL_BACKREAD_NOTE = ""
if _QUARK_D12_INTERNAL_BACKREAD is not None:
    _quark_closed = _QUARK_D12_INTERNAL_BACKREAD["closed_mass_side_package"]
    _QUARK_D12_INTERNAL_BACKREAD_NOTE = (
        " A separate continuation-only internal backread sidecar is on disk as "
        "`oph_quark_d12_internal_backread_continuation_closure`: using the emitted reference-free "
        "forward light-quark pair together with the explicit D12 backread assumptions, it fixes "
        f"`Delta_ud_overlap = {_quark_closed['Delta_ud_overlap']:.14f}`, "
        f"`t1 = {_quark_closed['t1']:.14f}`, "
        f"`eta_Q_centered = {_quark_closed['eta_Q_centered']:.14f}`, and "
        f"`kappa_Q = {_quark_closed['kappa_Q']:.14f}` on that continuation surface. "
        "That sidecar does not replace the public theorem frontier and does not repair the wrong-sheet CKM boundary."
    )
CHARGED_CONTINUATION_NOTE = (
    "No public charged value is emitted on the theorem lane. The repo contains an exact same-family witness on "
    "`current_family_only`, a conditional determinant-line lift on theorem-grade physical charged data, and the "
    "downstream algebraic readout from theorem-grade `A_ch(P)`. The single open charged theorem is the source landing "
    "`B_D10(P) -> Y_e(P)` or equivalently `B_D10(P) -> s_det(P)`."
)
_QUARK_GCH = _QUARK_SHARED_NORM_BINDING["g_ch"] if _QUARK_SHARED_NORM_BINDING is not None else 0.9231656602589082
_QUARK_SHARED_SCOPE = _QUARK_SHARED_NORM_BINDING.get("theorem_scope", "shared_budget_only") if _QUARK_SHARED_NORM_BINDING is not None else "shared_budget_only"
_QUARK_SIGMA_U = _QUARK_SECTOR_MEAN_SPLIT["sigma_u_total_log_per_side"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 5.5905
_QUARK_SIGMA_D = _QUARK_SECTOR_MEAN_SPLIT["sigma_d_total_log_per_side"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 3.3049
_QUARK_SIGMA_SEED = _QUARK_SECTOR_MEAN_SPLIT["sigma_seed_ud_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 4.4477
_QUARK_ETA_UD = _QUARK_SECTOR_MEAN_SPLIT["eta_ud_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 1.1428
_QUARK_A_UD = _QUARK_SECTOR_MEAN_SPLIT["A_ud_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 0.2467442927388686
_QUARK_B_UD = _QUARK_SECTOR_MEAN_SPLIT["B_ud_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 0.8125616987157023
_QUARK_GU = _QUARK_SECTOR_MEAN_SPLIT["g_u_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 0.7797392875757557
_QUARK_GD = _QUARK_SECTOR_MEAN_SPLIT["g_d_candidate"] if _QUARK_SECTOR_MEAN_SPLIT is not None else 0.12172551081512113
_QUARK_MEAN_SCOPE = _QUARK_SECTOR_MEAN_SPLIT.get("theorem_scope", "current_family_only") if _QUARK_SECTOR_MEAN_SPLIT is not None else "current_family_only"
_QUARK_SELECTOR_TOKEN = (
    _QUARK_RELATIVE_SELECTOR["quark_relative_sheet_selector"]["value"]["canonical_token"]
    if _QUARK_RELATIVE_SELECTOR is not None
    else "D12::same_label_left::reference_sheet"
)

QUARK_CONTINUATION_NOTE = (
    "Selected-class exact quark theorem surface. "
    f"`{_QUARK_PUBLIC_SIGMA_DESCENT['artifact'] if _QUARK_PUBLIC_SIGMA_DESCENT else 'oph_quark_public_physical_sigma_datum_descent'}` "
    "fixes the exact physical sigma datum on the public quark frame class chosen by `P`, and "
    f"`{_QUARK_PUBLIC_EXACT_YUKAWA['artifact'] if _QUARK_PUBLIC_EXACT_YUKAWA else 'oph_quark_public_exact_yukawa_end_to_end_theorem'}` "
    "emits the exact PDG 2025 running-quark sextet together with explicit exact forward Yukawas `Y_u` and `Y_d`. "
    "Supporting exact surfaces: `oph_quark_current_family_exact_readout` on `current_family_only` and "
    "`oph_quark_current_family_transport_frame_exact_pdg_completion` plus "
    "`oph_quark_current_family_transport_frame_exact_forward_yukawas` on the declared common-refinement transport-frame carrier. "
    "The D12 mass bridge is target-free on the emitted ray, and the sextet uses the PDG 2025 cross-section top entry rather than the auxiliary direct-top entry. "
    "Scope: selected-class closure only; no global classification of quark frame classes. "
    "Synchronization anchor: [#198](https://github.com/FloatingPragma/observer-patch-holography/issues/198)."
)
NEUTRINO_CONTINUATION_NOTE = (
    "Derived from `derive_neutrino_weighted_cycle_repair.py -> "
    "derive_neutrino_bridge_rigidity_theorem.py -> derive_neutrino_absolute_attachment_theorem.py -> "
    "export_forward_neutrino_closure_bundle.py`. The isotropic intrinsic branch is excluded by the exact atmospheric cap, "
    "and the weighted-cycle / Majorana-holonomy branch fixes the full scale-free PMNS/hierarchy shape with "
    "`theta12 = 34.2259 deg`, `theta23 = 49.7228 deg`, `theta13 = 8.68636 deg`, `delta = 305.581 deg`, "
    "`J = -0.02753`, and `Delta m21^2 / Delta m32^2 = 0.03072111`. "
    "The weighted-cycle bridge rigidity theorem emits the physical reduced bridge invariant "
    "`C_nu = sum_gap^2 * prod_qbar * solar_response_over_mstar^-0.5 = 0.9994295999075177`, "
    "and the emitted proxy is `P_nu = 6.699825740519345`. "
    "The absolute attachment theorem then emits `B_nu = P_nu * C_nu = 6.696004159297337`, "
    "`lambda_nu = (m_star_eV / q_mean^p_nu) * P_nu * C_nu = 1.7237014208357415`, and therefore one absolute family "
    "`m_i = lambda_nu * mhat_i`, `Delta m^2_ij = lambda_nu^2 * Delta_hat_ij`. "
    "The proof-facing neutrino lane therefore runs through the emitted theorem pair rather than through the compare-only adapter or the bridge corridor. "
    "The two-parameter exact adapter, the bridge corridor, and the reduced-correction candidate audit remain on disk only as diagnostic surfaces beneath the theorem lane."
)
HADRON_CONTINUATION_NOTE = (
    "Rows are suppressed by default because hadrons are execution-contract-frozen on the current branch rather than paper-derived outputs. The active hadron path is `derive_lambda_msbar_descendant.py -> "
    "derive_full_unquenched_correlator.py -> derive_stable_channel_cfg_source_measure_payload.py -> "
    "derive_runtime_schedule_receipt_n_therm_and_n_sep.py -> derive_stable_channel_sequence_population.py -> "
    "derive_hadron_production_geometry_summary.py -> derive_stable_channel_sequence_evaluation.py -> "
    "derive_stable_channel_groundstate_readout.py`, and a separate diagnostic-only surrogate bridge "
    "`derive_hadron_surrogate_execution_bridge_status.py` records that the full receipt/writeback/evaluation/convergence/systematics path "
    "has been closed on a surrogate HMC/RHMC kernel. The operational barrier has lower friction: `run_production_backend_writeback.py` executes the backend-export -> receipt -> dump -> payload -> evaluation -> closure-report path in one command once a real production export exists. The production geometry is explicit: 3 seeded 2+1 ensembles, 6 cfg total, naive raw gauge storage about "
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
        "prediction_key": "m_t",
        "ledger_id": "continuation.flavor.quarks",
        "note": QUARK_CONTINUATION_NOTE,
        "extra_prediction_keys": ["crit_mt_pole"],
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
    if not path.is_absolute():
        return path.as_posix()
    try:
        rel = path.relative_to(ROOT)
    except ValueError:
        return path.as_posix()
    return f"code/{rel.as_posix()}"


def _effective_hadron_profile(*, with_hadrons: bool, hadron_profile: str) -> str:
    return hadron_profile if with_hadrons else "suppressed"


def _deprogress_text(text: str) -> str:
    replacements = {
        "already packages": "packages",
        "already latent": "latent",
        "already-emitted": "emitted",
        "already schedule-independent": "schedule-independent",
        "already downstream of": "downstream of",
        "now exactly checked": "checked explicitly",
        "now integrated": "integrated",
        "now gives": "gives",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    return text


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


def _quark_public_exact_theorem_allowed(payload: Dict[str, Any]) -> bool:
    return bool(payload.get("public_promotion_allowed", False)) and payload.get(
        "proof_status"
    ) == "closed_target_free_public_exact_yukawa_end_to_end_theorem"


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
        neutrino_repaired_branch = neutrino_repaired_branch or neutrino_active
    if NEUTRINO_EXACT_BLOCKERS.exists():
        blockers = json.loads(NEUTRINO_EXACT_BLOCKERS.read_text(encoding="utf-8"))
        neutrino_repaired_branch = neutrino_repaired_branch or _neutrino_repaired_branch_waiting_absolute_scale(blockers)

    if QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.exists():
        theorem = json.loads(QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.read_text(encoding="utf-8"))
        quark_active = _quark_public_exact_theorem_allowed(theorem)
    elif FORWARD_YUKAWAS.exists() and QUARK_SECTOR_MEAN_SPLIT.exists():
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

    if NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.exists():
        theorem = json.loads(NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.read_text(encoding="utf-8"))
        if theorem.get("public_surface_candidate_allowed", False):
            masses_eV = [float(x) for x in theorem["outputs"]["masses_eV"]]
            updated.update(
                {
                    "m_nu_e": masses_eV[0] * 1.0e-9,
                    "m_nu_mu": masses_eV[1] * 1.0e-9,
                    "m_nu_tau": masses_eV[2] * 1.0e-9,
                }
            )

    if QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.exists():
        theorem = json.loads(QUARK_PUBLIC_EXACT_YUKAWA_THEOREM.read_text(encoding="utf-8"))
        if _quark_public_exact_theorem_allowed(theorem):
            exact = theorem["public_exact_outputs"]["exact_running_values_gev"]
            updated.update(
                {
                    "m_u": float(exact["u"]),
                    "m_c": float(exact["c"]),
                    "m_d": float(exact["d"]),
                    "m_s": float(exact["s"]),
                    "m_b": float(exact["b"]),
                    "m_t": float(exact["t"]),
                }
            )

    elif FORWARD_YUKAWAS.exists() and QUARK_SECTOR_MEAN_SPLIT.exists():
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
    absolute_theorem = (
        json.loads(NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.read_text(encoding="utf-8"))
        if NEUTRINO_ABSOLUTE_ATTACHMENT_THEOREM.exists()
        else None
    )
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

    if absolute_theorem and absolute_theorem.get("public_surface_candidate_allowed", False):
        theorem_dm2 = dict(absolute_theorem["outputs"]["delta_m_sq_eV2"])
        rows.extend(
            [
                _row(
                    observable_id="delta_m21_sq_eV2",
                    observable="Delta m21^2",
                    status="theorem_grade",
                    prediction_value=float(theorem_dm2["21"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m21_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="Theorem-grade absolute solar splitting from the emitted weighted-cycle bridge rigidity and absolute attachment theorems.",
                ),
                _row(
                    observable_id="delta_m32_sq_eV2",
                    observable="Delta m32^2",
                    status="theorem_grade",
                    prediction_value=float(theorem_dm2["32"]),
                    reference_value=NEUTRINO_PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"],
                    err_plus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["plus"],
                    err_minus=NEUTRINO_PDG_2025_NO_1SIGMA["delta_m32_sq_eV2"]["minus"],
                    unit="eV^2",
                    note="Theorem-grade absolute atmospheric splitting from the emitted weighted-cycle bridge rigidity and absolute attachment theorems.",
                ),
            ]
        )
    elif two_parameter_adapter:
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
                    note="Absolute solar splitting after compare-only anchoring with the atmospheric Delta m32^2 input; this is not a promoted theorem-grade OPH output.",
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


def build_majorana_phase_surface_rows(surface_state: Dict[str, Any]) -> List[Dict[str, Any]]:
    active = dict(surface_state["active_local_public_candidates"])
    if not active.get("neutrino_repaired_branch"):
        return []

    theorem = (
        json.loads(NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM.read_text(encoding="utf-8"))
        if NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM.exists()
        else None
    )
    if (
        theorem
        and theorem.get("status") == "theorem_grade_emitted"
        and theorem.get("public_surface_candidate_allowed", False)
    ):
        emitted = dict(theorem.get("emitted_parameters") or theorem.get("candidate_parameters") or {})
        alpha21_note = (
            "Theorem-grade physical Majorana phase on the repaired shared-basis weighted-cycle surface. "
            "The readout uses the canonical Takagi congruence of the emitted symmetric cycle matrix together with the electron-row gauge `U_e1 in R_{>0}`, so it is insensitive to computational column rephasings of the intermediate unitary."
        )
        alpha31_note = (
            "Same theorem surface as `alpha21^(Maj)`: emitted by the canonical Takagi congruence readout on the repaired shared-basis weighted-cycle matrix."
        )
        return [
            {
                "observable_id": "alpha21_majorana",
                "observable": "alpha21^(Maj)",
                "status": "theorem_grade",
                "prediction_value": float(emitted["alpha21_deg_0_to_360"]),
                "prediction_display": format_observable_value(float(emitted["alpha21_deg_0_to_360"]), "deg"),
                "unit": "deg",
                "note": alpha21_note,
            },
            {
                "observable_id": "alpha31_majorana",
                "observable": "alpha31^(Maj)",
                "status": "theorem_grade",
                "prediction_value": float(emitted["alpha31_deg_0_to_360"]),
                "prediction_display": format_observable_value(float(emitted["alpha31_deg_0_to_360"]), "deg"),
                "unit": "deg",
                "note": alpha31_note,
            },
        ]

    if theorem is None:
        alpha21_note = (
            "Still absent on the current public surface, not declared unphysical. "
            "The weighted-cycle branch already emits PMNS-type observables, masses, and splittings, "
            "but no declared public Majorana readout artifact is available yet."
        )
        alpha31_note = "Same current boundary as `alpha21^(Maj)`."
    else:
        blocker = str(
            theorem.get("public_promotion_blocker")
            or "The repaired weighted-cycle branch has not yet been represented explicitly on the closed shared basis for public Majorana promotion."
        )
        alpha21_note = (
            "Still absent on the current public surface, not declared unphysical. "
            "A canonical Takagi readout candidate exists on the repaired weighted-cycle matrix, "
            f"but public promotion remains open. {blocker}"
        )
        alpha31_note = (
            "Same current boundary as `alpha21^(Maj)`: still absent on the current public surface while the weighted-cycle "
            f"candidate remains below the public-promotion blocker. {blocker}"
        )

    return [
        {
            "observable_id": "alpha21_majorana",
            "observable": "alpha21^(Maj)",
            "status": "still_absent",
            "prediction_value": None,
            "prediction_display": "n/a",
            "unit": "",
            "note": alpha21_note,
        },
        {
            "observable_id": "alpha31_majorana",
            "observable": "alpha31^(Maj)",
            "status": "still_absent",
            "prediction_value": None,
            "prediction_display": "n/a",
            "unit": "",
            "note": alpha31_note,
        },
    ]


def prediction_surface_for_row(row_spec: Dict[str, Any], surface_state: Dict[str, Any], *, with_hadrons: bool) -> str:
    active = dict(surface_state["active_local_public_candidates"])
    particle_id = row_spec["particle_id"]
    if particle_id in {"photon", "gluon", "graviton"}:
        return "particles_structural_massless"
    if particle_id in {"w_boson", "z_boson"} and active.get("d10_mass_pair"):
        return "local_d10_public_mass_pair_candidate"
    if particle_id in {"electron", "muon", "tau"} and active.get("charged_local_candidate"):
        return "local_charged_public_candidate"
    if particle_id in {"electron_neutrino", "muon_neutrino", "tau_neutrino"} and (
        active.get("neutrino_local_candidate") or active.get("neutrino_repaired_branch")
    ):
        return "local_neutrino_weighted_cycle_absolute_attachment"
    if particle_id in {
        "up_quark",
        "down_quark",
        "strange_quark",
        "charm_quark",
        "bottom_quark",
        "top_quark",
    } and active.get("quark_forward_candidate"):
        return "selected_public_quark_exact_yukawa_theorem_surface"
    if particle_id in {"higgs", "top_quark"} and active.get("d11_forward_seed"):
        return "local_d11_forward_seed_candidate"
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
    majorana_rows: Optional[List[Dict[str, Any]]] = None,
) -> str:
    majorana_rows = majorana_rows or []
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
        "This table is a `/particles`-native audit surface. If a sector has no live local public candidate, the value is reported as `n/a`; legacy fallback predictors are not used.",
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
                f"- Why open: {_deprogress_text(uv_boundary['reason_current_corpus_fails'])}",
                f"- First theorem object: {_deprogress_text(uv_boundary['statement'])}",
                f"- Second theorem object: {_deprogress_text(uv_boundary['follow_on_statement'])}",
                f"- Candidate extension status: `{uv_boundary['candidate_extension_status']}`",
                f"- Filled witnesses on disk: `{uv_boundary['filled_contract_witnesses'][0]}`, `{uv_boundary['filled_contract_witnesses'][1]}`, `{uv_boundary['filled_contract_witnesses'][2]}`",
                f"- Remaining emitted witness: `{uv_boundary['remaining_missing_emitted_witness']}`",
                f"- Remaining witness formula: `{uv_boundary['remaining_missing_emitted_witness_formula']}`",
                f"- Smaller raw datum beneath that witness: `{uv_boundary['smaller_remaining_raw_datum']}`",
                f"- Smaller raw datum artifact: `{uv_boundary.get('smaller_remaining_raw_datum_artifact', str(UV_BW_FIXED_LOCAL_COLLAR_DATUM))}`",
                f"- Smallest exact blocker: `{uv_boundary.get('smallest_exact_blocker', 'n/a')}`",
                f"- Smallest exact blocker formula: `{uv_boundary.get('smallest_exact_blocker_formula', 'n/a')}`",
                f"- Candidate extension route: {_deprogress_text(uv_boundary['candidate_extension_route'])}",
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
        if group == "Leptons" and majorana_rows:
            lines.extend(
                [
                    "## Majorana Phase Surface",
                    "",
                    "| Observable | Status | OPH value | Note |",
                    "| --- | --- | --- | --- |",
                ]
            )
            for row in majorana_rows:
                lines.append(
                    f"| {row['observable']} | {row['status']} | {row['prediction_display']} | {row['note']} |"
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
    majorana_rows = build_majorana_phase_surface_rows(surface_state)
    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    effective_hadron_profile = _effective_hadron_profile(
        with_hadrons=with_hadrons,
        hadron_profile=str(args.hadron_profile),
    )

    markdown = render_markdown(
        rows=rows,
        comparison_rows=comparison_rows,
        majorana_rows=majorana_rows,
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
        "majorana_rows": majorana_rows,
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
                "majorana_rows": majorana_rows,
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
