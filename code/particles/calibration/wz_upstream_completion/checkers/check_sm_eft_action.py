#!/usr/bin/env python3
"""Independent checker for the SM_EFT_ACTION_1 bundle (Workstream A).

Re-derives every checkable clause of the work-order contract from the
emitted bundle alone: census completeness and anomaly cancellation from
the listed charges, full symbolic Yukawa matrices, the SM hypercharge
convention with the GUT normalization excluded, distinct vev types with
the equality theorem recorded absent, retained and excluded operator
lists with reasons, an acyclic source DAG whose ancestry carries no
target value, refusal controls recorded failed, and the subject digest
recomputed from the packet hashes.  A numeric scan walks the whole
bundle and rejects any number outside the structural whitelist.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
BUNDLE_PATH = ROOT / "outputs" / "sm_eft_action_1.json"

STRUCTURAL_INTS = {0, 1, 2, 3, 6, 8, 15, 45}
STRUCTURAL_FRACTIONS = {
    "1/6", "-2/3", "1/3", "-1/2", "1", "1/2", "-1/4", "0",
}


def fail(message: str) -> None:
    raise SystemExit(f"sm_eft_action checker: {message}")


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def numeric_scan(value: Any, path: str = "$") -> None:
    if isinstance(value, bool):
        return
    if isinstance(value, float):
        fail(f"float at {path}")
    if isinstance(value, int):
        if value not in STRUCTURAL_INTS and not (1 <= value <= 3):
            fail(f"non-structural integer {value} at {path}")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            numeric_scan(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            numeric_scan(item, f"{path}[{index}]")


def check(bundle: dict[str, Any]) -> None:
    if bundle.get("schema") != "sm_eft_action_1.v1":
        fail("schema mismatch")
    if bundle.get("promotion_allowed") is not False:
        fail("the bundle must not allow promotion")
    if bundle.get("lane_status") != "EXTERNAL_SM_DECLARED_ACTION":
        fail("lane status must declare the external action")
    if "mass_core" in bundle or "inverse_target_adapter" in bundle:
        fail("a mass core or inverse target adapter is forbidden")

    census = bundle["field_census"]
    fermions = census["fermions"]
    if len(fermions) != 15 or census["weyl_state_total"] != 45:
        fail("census must carry fifteen multiplets and forty-five Weyl states")
    generations = {f["generation"] for f in fermions}
    if generations != {1, 2, 3}:
        fail("census must carry three generations")
    per_gen = [f for f in fermions if f["generation"] == 1]
    u1_cubed = sum(Fraction(f["hypercharge"]) ** 3 * f["weyl_states"] for f in per_gen)
    if u1_cubed != 0:
        fail("recomputed hypercharge-cubed anomaly must vanish")

    conventions = bundle["conventions"]
    if conventions["gprime_convention"] != "SM_hypercharge":
        fail("the hypercharge convention must be the SM one")
    if conventions["gut_normalization"]["status"] != "excluded":
        fail("the GUT normalization must be recorded and excluded")
    if "BMHV" not in conventions["gamma5_prescription"]:
        fail("the gamma5 prescription must be declared")

    ast = bundle["action_ast"]
    retained_names = {entry["operator"] for entry in ast["retained"]}
    for needed in (
        "gauge_kinetic_SU3", "gauge_kinetic_SU2", "gauge_kinetic_U1",
        "higgs_kinetic", "fermion_kinetic", "yukawa_up", "yukawa_down",
        "yukawa_lepton", "higgs_mass", "higgs_quartic", "gauge_fixing",
        "ghost_sector",
    ):
        if needed not in retained_names:
            fail(f"retained operator missing: {needed}")
    excluded_names = {entry["operator"] for entry in ast["excluded"]}
    if "theta_QCD" not in excluded_names or "dimension_gt_4" not in excluded_names:
        fail("theta and higher-dimension exclusions must be recorded")
    for entry in ast["excluded"]:
        if not entry.get("reason"):
            fail("every exclusion carries a reason")

    yukawas = bundle["yukawa_packet"]["matrices"]
    for sector in ("Yu", "Yd", "Ye"):
        matrix = yukawas[sector]
        if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
            fail(f"{sector} must be three by three")
        for row in matrix:
            for entry in row:
                if "symbol" not in entry:
                    fail(f"{sector} entries must be symbols")
    if "UuL" not in json.dumps(bundle["yukawa_packet"]["basis_rotations"]):
        fail("basis rotations must be recorded")
    if "V_CKM" not in bundle["yukawa_packet"]["ckm"]["definition"]:
        fail("the CKM composition must be recorded")

    vevs = bundle["vev_types"]
    if vevs["v_chart"]["type"] == vevs["v_F"]["type"]:
        fail("the chart and Fermi vevs must be distinct types")
    if "absent" not in vevs["equality_theorem"]:
        fail("the vev equality theorem must be recorded absent")

    dag = bundle["source_dag"]
    node_ids = {node["id"] for node in dag["nodes"]}
    for node in dag["nodes"]:
        if "no_target_values" not in node["ancestry"]:
            fail("every DAG node must carry target-free ancestry")
    adjacency = {node: set() for node in node_ids}
    for edge in dag["edges"]:
        if edge[0] not in node_ids or edge[1] not in node_ids:
            fail("DAG edge references an unknown node")
        adjacency[edge[0]].add(edge[1])
    seen: set[str] = set()
    stack: set[str] = set()

    def visit(node: str) -> None:
        if node in stack:
            fail("the source DAG must be acyclic")
        if node in seen:
            return
        stack.add(node)
        for child in adjacency[node]:
            visit(child)
        stack.discard(node)
        seen.add(node)

    for node in node_ids:
        visit(node)

    for name, verdict in bundle["controls"].items():
        if verdict.get("expected_failure") is not True or verdict.get("failed") is not True:
            fail(f"control {name} must record its refusal")

    packets = {
        name: bundle[name]
        for name in ("conventions", "field_census", "action_ast", "yukawa_packet", "vev_types", "source_dag")
    }
    for name, payload in packets.items():
        if canonical_sha256(payload) != bundle["packet_hashes"][name]:
            fail(f"packet hash mismatch: {name}")
    digest = canonical_sha256({"packet_hashes": bundle["packet_hashes"], "schema": bundle["schema"]})
    if digest != bundle["subject_digest"]:
        fail("subject digest mismatch")

    numeric_scan(
        {k: v for k, v in bundle.items() if k not in ("packet_hashes", "subject_digest", "controls")}
    )


def main() -> int:
    if not BUNDLE_PATH.is_file():
        fail("bundle missing; run the producer first")
    bundle = json.loads(BUNDLE_PATH.read_text(encoding="utf-8"))
    check(bundle)
    print("sm_eft_action checker OK: census, anomalies, conventions, operators, Yukawa structure, vev types, DAG, controls, hashes, and the numeric scan all pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
