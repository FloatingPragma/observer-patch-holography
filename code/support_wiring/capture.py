"""Capture public geometry and committed cell transport; optional producer input."""
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
REVISION = "7faa47b5cf00b42f6bf6b3e95ed4f7eb64f4239f"
FILES = ("oph_fpe/core/icosahedral.py", "oph_fpe/core/screen_ports.py",
         "oph_fpe/gauge/covariant_overlap.py", "oph_fpe/finite_groups.py")


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def capture(simulator: Path, destination: Path):
    pins = {}
    for name in FILES:
        raw = subprocess.check_output(["git", "-C", str(simulator), "show", f"{REVISION}:{name}"])
        if raw != (simulator / name).read_bytes():
            raise ValueError("working source differs from pinned revision: " + name)
        pins[name] = hashlib.sha256(raw).hexdigest()
    sys.path.insert(0, str(simulator))
    from oph_fpe.core.icosahedral import build_geodesic_icosahedral_tower
    from oph_fpe.core.screen_ports import assign_echosahedral_ports
    old_path = ROOT / "code/source_routing/support_w12_l3.json"
    old = json.loads(old_path.read_text())
    tower = build_geodesic_icosahedral_tower(5)
    destination.mkdir(parents=True, exist_ok=True)
    manifest = {"schema": "oph.support_wiring.geometry.v1", "repository":
                "https://github.com/muellerberndt/oph-physics-sim", "revision": REVISION,
                "source_sha256": pins, "l3_original_sha256": hashlib.sha256(old_path.read_bytes()).hexdigest(),
                "levels": {}}
    for level in (3, 4, 5):
        mesh = tower.levels[level]
        centres = mesh.vertices[mesh.faces].sum(axis=1)
        centres /= np.linalg.norm(centres, axis=1, keepdims=True)
        vertices, edges = defaultdict(list), defaultdict(list)
        for c, face in enumerate(mesh.faces):
            for i in range(3):
                vertices[int(face[i])].append(c)
                edges[tuple(sorted((int(face[i]), int(face[(i + 1) % 3]))))].append(c)
        data = {"vertices": mesh.vertices, "faces": mesh.faces, "centres": centres,
                "areas": mesh.spherical_face_areas,
                "seams": np.asarray(old["intra_carrier_seams"], dtype=np.int64)}
        detail = {"geometry": mesh.receipt()}
        for wiring, incidence in (("w12", vertices), ("w3", edges)):
            pairs = np.asarray(sorted({(a, b) for cells in incidence.values()
                                       for a in cells for b in cells if a < b}), dtype=np.int64)
            port = assign_echosahedral_ports(pairs[:, 0], pairs[:, 1], len(mesh.faces), points=centres)
            assert port.overflow_count == 0
            glued = np.column_stack((pairs[:, 0], port.left_port, pairs[:, 1], port.right_port))
            data[wiring] = glued
            detail[wiring] = {"routing_mode": port.routing_mode, "local_frame_hash": port.local_frame_hash}
            if level == 3 and wiring == "w12":
                assert mesh.faces.tolist() == old["faces"]
                assert glued.tolist() == old["glued_pairs"]
        if level > 3:
            mapping = tower.cell_refinements[level - 1]
            data["parent"] = mapping.child_to_parent
            data["children"] = np.asarray(mapping.children_by_parent)
            data["expectation_weights"] = mapping.conditional_expectation_weights
            detail["refinement"] = mapping.receipt()
        path = destination / f"geometry_l{level}.npz"
        np.savez_compressed(path, **data)
        detail["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest["levels"][str(level)] = detail
        print(f"captured L{level}: {len(mesh.faces)} carriers, {len(data['w12'])} W12 glued seams", flush=True)
    (destination / "geometry.json").write_bytes(json_bytes(manifest))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("simulator", type=Path)
    p.add_argument("--output", type=Path, default=HERE / "geometry")
    a = p.parse_args()
    capture(a.simulator.resolve(), a.output)
