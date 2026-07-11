#!/usr/bin/env python3
"""Guard explicit public-surface policy helpers in the claim table builder."""

from __future__ import annotations

import importlib.util
import pathlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_results_status_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_results_status_table", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_d10_and_d11_public_candidate_policy_is_explicit() -> None:
    module = _load_module()
    d10 = {
        "public_surface_candidate_allowed": True,
        "prediction_promotion_allowed": True,
        "mass_pair_predictive_candidate": {"MW_pole": 1.0, "MZ_pole": 2.0},
    }
    d11 = {
        "public_surface_candidate_allowed": True,
        "mass_readout": {"mH_gev": 3.0, "mt_pole_gev": 4.0},
    }
    assert module._d10_public_mass_pair_allowed(d10) is True
    assert module._d11_public_seed_allowed(d11) is True
    d10["prediction_promotion_allowed"] = False
    assert module._d10_public_mass_pair_allowed(d10) is False


def test_neutrino_repaired_branch_policy_is_explicit() -> None:
    module = _load_module()
    blockers = {
        "exact_blockers": [
            {
                "name": "one_positive_neutrino_bridge_correction_invariant",
                "kind": "reduced_bridge_correction_invariant",
            }
        ],
        "live_continuation_branch_status": {
            "status": "physically_repaired_up_to_one_reduced_bridge_correction_invariant",
            "same_label_scalar_certificate_present": True,
            "shared_charged_left_basis_present": True,
            "repair_artifact_present": True,
        }
    }
    assert module._neutrino_repaired_branch_waiting_absolute_scale(blockers) is True


def test_quark_public_forward_policy_uses_explicit_surface_gate() -> None:
    module = _load_module()
    forward = {
        "public_surface_candidate_allowed": True,
        "source_mode": "factorized_descent",
        "promotion_blockers": [],
    }
    mean_split = {"active_candidate": "predictive_candidate"}
    assert module._quark_public_forward_allowed(forward, mean_split) is True


def test_surface_state_exposes_particles_native_policy_and_local_public_candidates() -> None:
    module = _load_module()
    state = module.build_surface_state(with_hadrons=False)
    assert state["public_surface_kind"] == "particles_native_candidate_or_gap_surface"
    assert state["surface_policy"] == "local_candidate_or_gap_only"
    assert "active_local_public_candidates" in state
    active = state["active_local_public_candidates"]
    assert set(active) == {
        "d10_mass_pair",
        "d11_forward_seed",
        "charged_local_candidate",
        "neutrino_local_candidate",
        "neutrino_repaired_branch",
        "quark_forward_candidate",
        "hadrons_enabled",
    }
    assert active["charged_local_candidate"] is False
    assert active["neutrino_local_candidate"] is False
    assert active["neutrino_repaired_branch"] is False
    assert active["hadrons_enabled"] is False


def test_neutrino_rows_get_repaired_branch_surface_before_absolute_scale() -> None:
    module = _load_module()
    surface = module.prediction_surface_for_row(
        {
            "particle_id": "electron_neutrino",
            "group": "Leptons",
        },
        {
            "active_local_public_candidates": {
                "d10_mass_pair": False,
                "d11_forward_seed": False,
                "charged_local_candidate": False,
                "neutrino_local_candidate": False,
                "neutrino_repaired_branch": True,
                "quark_forward_candidate": False,
                "hadrons_enabled": False,
            }
        },
        with_hadrons=False,
    )
    assert surface == "rejected_target_informed_neutrino_candidate"


def test_premise_boundaries_use_repo_stable_artifact_refs() -> None:
    module = _load_module()
    uv = module.build_premise_boundaries()["uv_bw_internalization"]
    assert uv["prelimit_system_artifact"] == "code/particles/runs/uv/bw_realized_transported_cap_local_system.json"
    assert uv["remaining_missing_emitted_witness_artifact"] == (
        "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json"
    )
    assert uv["smaller_remaining_raw_datum_artifact"] == (
        "code/particles/runs/uv/bw_fixed_local_collar_markov_faithfulness_datum.json"
    )
    assert uv["neutrino_local_bridge_candidate_context"] == (
        "code/particles/runs/neutrino/neutrino_lambda_nu_bridge_candidate.json"
    )


def test_top_note_uses_preserved_sidecar_value() -> None:
    module = _load_module()
    updated = module.apply_local_candidate_overrides(
        {
            "m_t": 174.49070203607485,
            "crit_mt_pole": 171.1,
            "crit_mH_tree": 126.5,
        }
    )
    note = module.build_note(
        {
            "particle_id": "top_quark",
            "label": "t",
            "group": "Quarks",
            "prediction_key": "crit_mt_pole",
            "ledger_id": "secondary.d11.higgs_top",
            "note": "test",
        },
        {"reference_kind": "value", "value_gev": 0.0, "comment": None},
        updated,
        {"blocked_by": [], "tier": "secondary_quantitative", "label": "secondary_quantitative"},
    )
    assert "174.490702" in note
