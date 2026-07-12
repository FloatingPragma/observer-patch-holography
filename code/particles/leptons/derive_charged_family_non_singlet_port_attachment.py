#!/usr/bin/env python3
"""Conditional icosahedral moment theorem for charged-family attachment.

The live hierarchy certificate supplies an equal invariant load and therefore
does not populate this theorem's non-singlet input.  Given a future
refinement-stable, source-emitted twelve-port record, this module constructs
ID-independent first and quadrupole moments and tests the simple-spectrum gate
that canonically yields three family lines and two centered shape coordinates.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HIERARCHY = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_screen_sieve_icosahedral_certificate.json"
)
DEFAULT_PACKET = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_12_port_non_singlet_record.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_family_non_singlet_port_attachment.json"
)


def icosahedron_vertices() -> np.ndarray:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices: list[tuple[float, float, float]] = []
    for a in (-1.0, 1.0):
        for b in (-phi, phi):
            vertices.append((0.0, a, b))
            vertices.append((a, b, 0.0))
            vertices.append((b, 0.0, a))
    points = np.asarray(vertices, dtype=float)
    return points / np.linalg.norm(points[0])


def port_moments(weights: Sequence[float], vertices: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    points = icosahedron_vertices() if vertices is None else np.asarray(vertices, dtype=float)
    w = np.asarray(weights, dtype=float)
    if points.shape != (12, 3) or w.shape != (12,):
        raise ValueError("expected twelve port weights and twelve three-dimensional vertices")
    centered_w = w - np.mean(w)
    first = np.einsum("i,ij->j", centered_w, points)
    identity = np.eye(3)
    quadrupole = sum(
        centered_w[i] * (np.outer(point, point) - identity / 3.0)
        for i, point in enumerate(points)
    )
    return first, quadrupole


def attachment_readout(weights: Sequence[float], x2: float, tolerance: float = 1.0e-10) -> dict[str, Any]:
    first, quadrupole = port_moments(weights)
    eigenvalues, eigenvectors = np.linalg.eigh(quadrupole)
    gaps = np.diff(eigenvalues)
    simple = bool(np.min(gaps) > tolerance)
    sigma = float(eigenvalues[2] - eigenvalues[0])
    eta = float((1.0 + x2) * sigma - 2.0 * (eigenvalues[1] - eigenvalues[0]))
    return {
        "first_moment": first.tolist(),
        "quadrupole": quadrupole.tolist(),
        "trace_residual": float(np.trace(quadrupole)),
        "ordered_eigenvalues": eigenvalues.tolist(),
        "eigenvalue_gaps": gaps.tolist(),
        "simple_spectrum": simple,
        "family_lines_columns_sign_ambiguous": eigenvectors.tolist(),
        "sigma_source_support_extension_total_log_per_side": sigma if simple else None,
        "eta_source_support_extension_log_per_side": eta if simple else None,
    }


def build_artifact(hierarchy: dict[str, Any], packet: dict[str, Any] | None) -> dict[str, Any]:
    if int(hierarchy["orbit_stabilizer"]["orbit_size"]) != 12:
        raise ValueError("attachment theorem requires the declared twelve-port orbit")
    packet_ok = bool(
        packet
        and packet.get("artifact") == "oph_charged_12_port_non_singlet_record"
        and packet.get("source_only") is True
        and packet.get("frozen_before_charged_mass_comparison") is True
        and len(packet.get("port_weights", [])) == 12
    )
    readout = attachment_readout(packet["port_weights"], float(packet["x2"])) if packet_ok else None
    simple = bool(readout and readout["simple_spectrum"])
    promotion = packet_ok and simple
    blockers: list[str] = []
    if not packet_ok:
        blockers.append("source_emitted_frozen_12_port_non_singlet_record")
    if packet_ok and not simple:
        blockers.append("simple_quadrupole_spectrum")
    return {
        "artifact": "oph_charged_family_non_singlet_port_attachment",
        "status": "CLOSED_CONDITIONAL_THEOREM_LIVE_INPUT_OPEN",
        "public_promotion_allowed": promotion,
        "live_input_status": "available" if packet_ok else "absent",
        "blockers": blockers,
        "conditional_theorem": {
            "id": "charged_family_non_singlet_port_attachment",
            "statement": (
                "Given a refinement-stable source-only scalar record w_i on the twelve "
                "icosahedral ports, form the centered first moment m and traceless second "
                "moment Q. These tensors are invariant under port relabeling and covariant "
                "under the icosahedral action. If Q has simple spectrum, its three spectral "
                "projectors define three ID-independent unoriented family lines. Ordering "
                "the eigenvalues defines two independent adjacent gaps spanning the centered "
                "three-family plane. Their invertible change of coordinates emits the live "
                "charged support-extension scalars sigma and eta."
            ),
            "moments": {
                "first": "m = sum_i (w_i - mean(w)) p_i",
                "quadrupole": "Q = sum_i (w_i - mean(w)) (p_i p_i^T - I/3)",
            },
            "shape_map": {
                "sigma": "lambda_3 - lambda_1",
                "eta": "(1+x2)*(lambda_3-lambda_1) - 2*(lambda_2-lambda_1)",
                "jacobian_determinant_from_adjacent_gaps": -2.0,
            },
            "id_independence": (
                "Simultaneous permutation of port weights and coordinates leaves m and Q "
                "unchanged; eigenline signs are physically irrelevant projective choices."
            ),
        },
        "live_hierarchy_boundary": (
            "The current screen certificate contains twelve equal unit curvature defects and "
            "only the invariant X/12 load. It emits no nonuniform refinement-stable w_i record."
        ),
        "no_target_leak_contract": {
            "required_packet_artifact": "oph_charged_12_port_non_singlet_record",
            "required_flags": ["source_only=true", "frozen_before_charged_mass_comparison=true"],
            "forbidden_ancestry": ["m_e", "m_mu", "m_tau", "PDG", "CODATA", "charged target ratios"],
        },
        "readout": readout,
        "next_exact_missing_object": "refinement_stable_source_only_nonuniform_12_port_record",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hierarchy", type=Path, default=HIERARCHY)
    parser.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    hierarchy = json.loads(args.hierarchy.read_text(encoding="utf-8"))
    packet = json.loads(args.packet.read_text(encoding="utf-8")) if args.packet.exists() else None
    artifact = build_artifact(hierarchy, packet)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"], "blockers": artifact["blockers"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
