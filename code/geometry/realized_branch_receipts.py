#!/usr/bin/env python3
"""Realized-branch receipt evaluation for the Einstein branch-entry gate (#503).

Evaluates the decidable incidence receipts of the compact paper's
geometry-producer packet (Definition `def:spherical-incidence-receipt`) on the
currently exported realized consensus artifacts, and emits the machine status
consumed by Remark `rem:realized-branch-status`:

* the verified rooted-tree packet-net domain (issue #238 export) is the only
  repair domain with a full correctness export to date; its overlap structure
  is a tree, so the support-visible incidence complex has no 2-simplices,
  Euler characteristic 1, and fails the spherical incidence receipt --- the
  realized geometric branch is NOT certified nonempty by it;
* a designed cellulated-sphere mesh (icosahedron) passes the receipt but is
  declared geometry, not a consensus product, and is labeled as such.

Run:
    python3 code/geometry/realized_branch_receipts.py
writes code/geometry/runs/realized_branch_receipt_report.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from quotient_cap_readout import (  # noqa: E402
    IncidenceComplex,
    classify_surface,
    euler_characteristic,
    icosahedron,
    incidence_complex,
    is_closed_surface,
    is_connected,
    orient,
    spherical_incidence_receipt,
)

ROOT = HERE.parents[1]
TREE_ARTIFACT = ROOT / "code" / "consensus" / "runs" / "verified_tree_packet_net_domain.json"
CYCLIC_ARTIFACT = HERE / "runs" / "cyclic_cap_net_run_domain.json"
MODULAR_ARTIFACT = HERE / "runs" / "modular_clock_instrumentation_report.json"
NULLNET_ARTIFACT = HERE / "runs" / "null_net_receipt_report.json"
EVENT_ARTIFACT = HERE / "runs" / "realized_event_receipt_report.json"
BULK_ARTIFACT = HERE / "runs" / "bulk_depth_receipt_report.json"
REPORT_PATH = HERE / "runs" / "realized_branch_receipt_report.json"


def tree_packet_net_records() -> list[tuple[int, ...]]:
    """Overlap records of the verified tree packet-net domain: one record per
    oriented interface edge (the artifact's own graph payload)."""
    with open(TREE_ARTIFACT) as f:
        artifact = json.load(f)
    graph = artifact["domain"]["graph"]
    vertices = {name: idx for idx, name in enumerate(graph["vertices"])}
    records = []
    for edge in graph["oriented_edges"]:
        a, b = edge.split("-")
        records.append(tuple(sorted((vertices[a], vertices[b]))))
    return records


def receipt_report_for(records: list[tuple[int, ...]]) -> dict:
    k = incidence_complex(records)
    return {
        "n_vertices": len(k.vertices),
        "n_edges": len(k.edges),
        "n_triangles": len(k.triangles),
        "euler_characteristic": euler_characteristic(k),
        "connected": is_connected(k),
        "closed_surface": is_closed_surface(k),
        "orientable_coherent": orient(k) is not None if k.triangles else False,
        "spherical_incidence_receipt": spherical_incidence_receipt(k),
        "surface_classification": classify_surface(k),
    }


def cyclic_run_summary() -> dict | None:
    """Summary of the cyclic cap-net repair run artifact, if present."""
    if not CYCLIC_ARTIFACT.exists():
        return None
    with open(CYCLIC_ARTIFACT) as f:
        run = json.load(f)
    return {
        "source": "code/geometry/runs/cyclic_cap_net_run_domain.json",
        "provenance": run["provenance"],
        "receipts_witnessed": run["receipts_witnessed"],
        "receipts_pending": run["receipts_pending"],
    }


def modular_instrumentation_summary() -> dict | None:
    """Summary of the modular-clock instrumentation artifact, if present."""
    if not MODULAR_ARTIFACT.exists():
        return None
    with open(MODULAR_ARTIFACT) as f:
        run = json.load(f)
    return {
        "source": "code/geometry/runs/modular_clock_instrumentation_report.json",
        "scope": run["scope"],
        "receipts_witnessed": run["receipts_witnessed"],
        "receipts_pending": run["receipts_pending"],
        "kms_final_median_residual": run["verdicts"]["kms_final_median_residual"],
        "crossratio_final_relative_error":
            run["verdicts"]["crossratio_final_relative_error"],
    }


def build_report() -> dict:
    tree_records = tree_packet_net_records()
    tree = receipt_report_for(tree_records)

    ico_tris, _ = icosahedron()
    designed = receipt_report_for(list(ico_tris))

    cyclic = cyclic_run_summary()
    modular = modular_instrumentation_summary()
    nullnet = None
    if NULLNET_ARTIFACT.exists():
        with open(NULLNET_ARTIFACT) as f:
            nn = json.load(f)
        nullnet = {
            "source": "code/geometry/runs/null_net_receipt_report.json",
            "scope": nn["scope"],
            "receipts_witnessed": nn["receipts_witnessed"],
            "receipts_pending": nn["receipts_pending"],
        }
    events = None
    if EVENT_ARTIFACT.exists():
        with open(EVENT_ARTIFACT) as f:
            ev = json.load(f)
        events = {
            "source": "code/geometry/runs/realized_event_receipt_report.json",
            "scope": ev["scope"],
            "receipts_witnessed": ev["receipts_witnessed"],
            "receipts_pending": ev["receipts_pending"],
        }
    bulk = None
    if BULK_ARTIFACT.exists():
        with open(BULK_ARTIFACT) as f:
            bd = json.load(f)
        bulk = {
            "source": "code/geometry/runs/bulk_depth_receipt_report.json",
            "scope": bd["scope"],
            "receipts_witnessed": bd["receipts_witnessed"],
            "caveats": bd["caveats"],
            "timelike_eigenvalues_across_seeds":
                bd["timelike_eigenvalues_across_seeds"],
            "countermodel_signature":
                bd["countermodel_global_coupling"]["signature"],
        }

    evaluations = {
        "verified_tree_packet_net_domain": {
            "source": str(TREE_ARTIFACT.relative_to(ROOT)),
            "provenance": "consensus_product (issue #238 verified export)",
            **tree,
        },
        "designed_cellulated_sphere_icosahedron": {
            "source": "code/geometry/quotient_cap_readout.py:icosahedron()",
            "provenance": "declared_geometry (NOT a consensus product)",
            **designed,
        },
    }
    if cyclic is not None:
        evaluations["cyclic_cap_net_repair_run"] = cyclic
    if modular is not None:
        evaluations["modular_clock_instrumentation"] = modular
    if nullnet is not None:
        evaluations["null_net_receipts"] = nullnet
    if events is not None:
        evaluations["realized_event_receipts"] = events
    if bulk is not None:
        evaluations["bulk_depth_receipts"] = bulk

    # full nonemptiness needs ALL finite receipt families on one realized
    # tower; the cyclic run witnesses the D1 + incidence + mesh + naturality
    # families with explicit branch selection, and the boundary-collar
    # instrumentation witnesses the cross-ratio and 2pi-KMS families on the
    # declared Gaussian MaxEnt states; the cap-interior, null-net, event,
    # and physical-identification families remain pending, so the gate
    # stays open.
    topology = bool(cyclic is not None
                    and all(cyclic["receipts_witnessed"].values()))
    modular_ok = bool(modular is not None
                      and all(modular["receipts_witnessed"].values()))
    nullnet_ok = bool(nullnet is not None
                      and all(nullnet["receipts_witnessed"].values()))
    events_screen_ok = bool(
        events is not None and all(
            v for k, v in events["receipts_witnessed"].items()
            if k != "e3_rank_four_bulk_depth"
        )
    )
    bulk_depth_ok = bool(
        bulk is not None
        and all(bulk["receipts_witnessed"].values())
    )
    if topology and modular_ok and nullnet_ok:
        event_clause = (
            "the realized-event instrumentation witnesses screen "
            "population, certified separation, the Moebius chart cocycle, "
            "and an intrinsic Lorentzian (1,2) ancestry cone on the "
            "screen sheet; " if events_screen_ok else
            "realized-event evaluation incomplete; "
        )
        bulk_clause = (
            "the bulk-depth channel witnesses E3 rank four, bulk PCA "
            "dimension four, and a seed-stable intrinsic ancestry cone of "
            "Lorentzian signature (1,3) with the modular clock timelike, "
            "on the corpus-consistent inverse-system dynamics, with the "
            "narrow-margin/window caveats explicit and the "
            "strong-coupling (2,2) countermodel preserved; "
            if bulk_depth_ok else
            "the bulk-depth channel is not witnessed; "
        )
        status = (
            "OPEN: the cyclic cap-net repair run realizes the D1 repair "
            "clauses with the spherical-incidence, mesh, and naturality "
            "receipts (explicit branch selection); the free-fermion "
            "boundary-collar instrumentation witnesses the modular "
            "cross-ratio and geometric 2pi-KMS receipts; the null-net "
            "instrumentation witnesses NTI, weak additivity, separating "
            "faithfulness, mixed-GNS Cauchy, one-particle HSM compression, "
            "and percent-level modular Lie closure (rate clause open); "
            + event_clause + bulk_clause +
            "remaining pending: cone-margin convergence with commit "
            "density, Cyc limit clause, Lie-closure rate, second-quantized "
            "MI, cap-interior data, and the UC/VR/scale "
            "physical-identification receipts; #503 stays open on the "
            "full nonemptiness clause."
        )
    elif topology and modular_ok:
        status = (
            "OPEN: topology/mesh and boundary-collar modular families "
            "realized; null-net receipt evaluation incomplete."
        )
    elif topology:
        status = (
            "OPEN: topology/mesh families realized by the cyclic run; "
            "modular receipt families pending."
        )
    else:
        status = (
            "OPEN: no realized consensus artifact carries the spherical "
            "incidence receipt."
        )
    return {
        "artifact": "einstein_branch_realized_receipt_evaluation",
        "issue": 503,
        "paper_anchor": "rem:realized-branch-status",
        "evaluations": evaluations,
        "realized_geometric_branch_certified_nonempty": False,
        "topology_mesh_families_realized_with_branch_selection": topology,
        "boundary_collar_modular_families_witnessed": modular_ok,
        "null_net_families_witnessed_one_particle": nullnet_ok,
        "screen_event_families_witnessed": events_screen_ok,
        "bulk_depth_channel_witnessed": bulk_depth_ok,
        "status": status,
    }


def main() -> None:
    report = build_report()
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT_PATH, "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
