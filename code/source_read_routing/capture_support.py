"""Capture full shared-vertex W12 wiring from the public simulator primitives.

The unpublished support_wiring.py is not used.  Level three must reproduce
the existing fixture exactly before either larger support is exported.
This optional capture needs a simulator checkout; replay needs only fixtures.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def capture(simulator, revision):
    paths = ["oph_fpe/core/icosahedral.py", "oph_fpe/core/screen_ports.py",
             "oph_fpe/gauge/covariant_overlap.py", "oph_fpe/finite_groups.py"]
    # A source revision is an actual content pin, not merely worktree context.
    for path in paths:
        committed = subprocess.check_output(["git", "-C", str(simulator), "show", f"{revision}:{path}"])
        if committed != (simulator/path).read_bytes():
            raise ValueError("uncommitted simulator source differs from its revision: "+path)
    sys.path.insert(0, str(simulator))
    from oph_fpe.core.icosahedral import build_geodesic_icosahedral_tower
    from oph_fpe.core.screen_ports import assign_echosahedral_ports

    old = json.loads((ROOT / "code/source_routing/support_w12_l3.json").read_text(encoding="utf-8"))
    source = {"repository": "https://github.com/muellerberndt/oph-physics-sim",
              "revision": revision,
              "files": {p: hashlib.sha256((simulator / p).read_bytes()).hexdigest() for p in paths}}
    tower = build_geodesic_icosahedral_tower(5)
    for level in (3, 4, 5):
        mesh = tower.levels[level]
        incidence = defaultdict(list)
        for cell, face in enumerate(mesh.faces):
            for vertex in face:
                incidence[int(vertex)].append(cell)
        pairs = sorted({(a, b) for cells in incidence.values()
                        for a in cells for b in cells if a < b})
        endpoints = np.asarray(pairs, dtype=np.int64)
        points = mesh.vertices[mesh.faces].sum(axis=1)
        points /= np.linalg.norm(points, axis=1, keepdims=True)
        ports = assign_echosahedral_ports(endpoints[:, 0], endpoints[:, 1], len(mesh.faces), points=points)
        if ports.overflow_count:
            raise ValueError("twelve-port assignment overflow")
        glued = np.column_stack((endpoints[:, 0], ports.left_port,
                                 endpoints[:, 1], ports.right_port)).tolist()
        if level == 3:
            if mesh.faces.tolist() != old["faces"] or glued != old["glued_pairs"]:
                raise ValueError("public primitive reconstruction differs from captured W12 L3")
            print("Public primitives reproduce all captured L3 faces and glued port pairs.", flush=True)
            continue
        packet = {"schema": "oph.source_read_routing.support.v1", "level": level,
                  "carriers": len(mesh.faces), "faces": mesh.faces.tolist(),
                  "intra_carrier_seams": old["intra_carrier_seams"], "glued_pairs": glued,
                  "source": source, "l3_fixture_sha256": hashlib.sha256(
                      (ROOT / "code/source_routing/support_w12_l3.json").read_bytes()).hexdigest(),
                  "scope": "Declared shared-vertex W12 geometric port assignment; no source selection or physical geometry claim."}
        path = HERE / f"support_w12_l{level}.json"
        path.write_bytes(canonical(packet))
        print(f"{path.name}: {packet['carriers']} carriers, {len(glued)} glued seams", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulator", type=Path)
    parser.add_argument("--revision")
    args = parser.parse_args()
    revision = args.revision or subprocess.check_output(
        ["git", "-C", str(args.simulator), "rev-parse", "HEAD"], text=True).strip()
    capture(args.simulator.resolve(), revision)
