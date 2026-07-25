#!/usr/bin/env python3
"""Audit the conditional absolute-amplitude freedom of the weighted-cycle candidate.

Chain role: preserve a downstream scale diagnostic without mistaking it for the
first physical gap in the neutrino lane.

Mathematics: conditional on the rejected weighted-cycle ansatz, one obtains a scale-free normal
form `m_i = lambda_nu * mhat_i`, `Delta m^2_ij = lambda_nu^2 * Delta_hat_ij`.
The current D10 scale anchor `m_star` is dimensionful, but no source-closed theorem
attaches it to that candidate normal form. This script records the conditional
scale freedom and direct-attachment diagnostic; it does not close the upstream
operator, physical-basis, or mass-label problems.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_REPAIR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_SCALE_ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_absolute_amplitude_bridge.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the repaired neutrino absolute-amplitude bridge audit.")
    parser.add_argument("--repair", default=str(DEFAULT_REPAIR))
    parser.add_argument("--scale-anchor", default=str(DEFAULT_SCALE_ANCHOR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    repair = _load_json(Path(args.repair))
    scale_anchor = _load_json(Path(args.scale_anchor))

    if repair.get("artifact") != "oph_neutrino_weighted_cycle_repair":
        raise SystemExit("repair artifact mismatch")
    if scale_anchor.get("artifact") != "oph_neutrino_scale_anchor":
        raise SystemExit("scale-anchor artifact mismatch")

    mhat = [float(x) for x in repair["scale_free_mass_normal_form"]["masses"]]
    dh = {key: float(val) for key, val in repair["scale_free_dm2_normal_form"]["dm2"].items()}
    m_star_eV = float(scale_anchor["anchors"]["m_star_gev"]) * 1.0e9

    direct_masses = [m_star_eV * x for x in mhat]
    direct_dm2 = {key: m_star_eV * m_star_eV * val for key, val in dh.items()}
    compare_only_anchor = dict(repair.get("compare_only_atmospheric_anchor") or {})
    compare_anchor_dm32 = float(compare_only_anchor.get("delta_m32_sq_input_eV2") or 0.0)
    direct_shortfall = None
    if compare_anchor_dm32 > 0.0 and direct_dm2["32"] > 0.0:
        direct_shortfall = compare_anchor_dm32 / direct_dm2["32"]

    payload = {
        "artifact": "oph_neutrino_weighted_cycle_absolute_amplitude_bridge",
        "generated_utc": _timestamp(),
        "status": "conditional_scale_audit_on_rejected_source_open_candidate",
        "theorem_status": "not_established_physical_branch_upstream_gates_open",
        "public_promotion_allowed": False,
        "weighted_cycle_source_closure_status": repair.get("source_closure_status"),
        "source_artifacts": {
            "neutrino_weighted_cycle_repair": str(Path(args.repair)),
            "neutrino_scale_anchor": str(Path(args.scale_anchor)),
        },
        "fixed_live_objects": {
            "scale_free_mass_normal_form": repair.get("scale_free_mass_normal_form"),
            "scale_free_dm2_normal_form": repair.get("scale_free_dm2_normal_form"),
            "pmns_observables": repair.get("pmns_observables"),
            "m_star_gev": scale_anchor["anchors"]["m_star_gev"],
            "m_star_formula": scale_anchor["anchors"]["m_star_formula"],
        },
        "no_go_statement": (
            "Conditional on the rejected weighted-cycle ansatz, the stored normal form retains one positive scale. "
            "The live D10 scale anchor does not supply a source-closed attachment to that candidate. More importantly, "
            "the physical lane is blocked earlier by the source-open neutrino operator, charged basis, and mass-label rule."
        ),
        "direct_scale_anchor_attachment_diagnostic": {
            "candidate_rule": "lambda_nu = m_star_eV",
            "m_star_eV": m_star_eV,
            "absolute_masses_eV": direct_masses,
            "absolute_dm2_eV2": direct_dm2,
            "compare_only_atmospheric_shortfall_factor": direct_shortfall,
            "status": "diagnostic_only_nonpromotable",
        },
        "remaining_object": "neutrino_weighted_cycle_absolute_amplitude_bridge",
        "remaining_object_kind": "dimensionful_attachment_theorem",
        "first_physical_missing_object": "source_closed_neutrino_operator_basis_and_mass_label_contract",
        "remaining_object_contract": {
            "must_emit": "A_nu > 0",
            "role": "dimensionful amplitude that upgrades the repaired weighted-cycle normal form to absolute masses and splittings",
            "absolute_mass_rule": "m_i = A_nu * mhat_i",
            "absolute_splitting_rule": "Delta m^2_ij = A_nu^2 * Delta_hat_ij",
            "no_external_inputs": [
                "measured_oscillation_targets",
                "external_eV_anchor",
            ],
            "specialization_requirement": "On the isotropic specialization, the bridge must reduce to the live D10 scale anchor or an emitted equivalent attachment law.",
        },
        "equivalent_scalar": {
            "name": "lambda_nu",
            "relation_to_bridge": "lambda_nu = A_nu in the eV-valued normalization used by the repaired weighted-cycle artifact",
        },
        "notes": [
            "The dimensionless selector, topology, exponent, basis placement, and ordering are not source-derived; this artifact must not imply otherwise.",
            "Only after a source-closed physical branch exists would its residual absolute-scale attachment become a physical theorem target.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
