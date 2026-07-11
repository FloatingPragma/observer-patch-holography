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


def build_report() -> dict:
    tree_records = tree_packet_net_records()
    tree = receipt_report_for(tree_records)

    ico_tris, _ = icosahedron()
    designed = receipt_report_for(list(ico_tris))

    realized_nonempty = bool(tree["spherical_incidence_receipt"])
    return {
        "artifact": "einstein_branch_realized_receipt_evaluation",
        "issue": 503,
        "paper_anchor": "rem:realized-branch-status",
        "evaluations": {
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
        },
        "realized_geometric_branch_certified_nonempty": realized_nonempty,
        "status": (
            "OPEN: no realized consensus artifact carries the spherical "
            "incidence receipt; the composed branch-entry theorem is "
            "conditional and #503 remains open on its nonemptiness clause."
            if not realized_nonempty else
            "REALIZED: a consensus artifact carries the receipt."
        ),
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
