"""Infer binary spatial directions from anonymous material station measurements."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from .clocked_apparatus import clock_laws, patches
from .clocked_field import stopping_law
from .compiler import Compiler
from .gate import require

HERE = Path(__file__).resolve().parent
BOOK = HERE / "blind_chart_protocol.json"


def infer(calibration):
    n = len(calibration)
    symmetric = calibration + calibration.T
    ordered = np.sort(symmetric[np.triu_indices(n, k=1)])
    gaps = np.diff(ordered)
    split = np.argmax(gaps)
    require(
        gaps[split] > 0 and np.count_nonzero(gaps == gaps[split]) == 1,
        "unresolved edge class",
    )
    edges = symmetric > ordered[split]
    np.fill_diagonal(edges, False)
    degree = edges.sum(axis=1)
    require(min(degree) > 0 and np.all(degree == degree[0]), "irregular station menu")
    hop = shortest_path(edges, directed=False, unweighted=True)
    require(np.isfinite(hop).all(), "disconnected stations")
    cuts = set()
    for a, b in zip(*np.where(np.triu(edges, 1)), strict=True):
        require(np.all(hop[a] != hop[b]), "edge has an equidistant station")
        bit = hop[a] < hop[b]
        cuts.add(tuple((bit ^ bit[0]).astype(int).tolist()))
    xyz = np.array(sorted(cuts)).T
    hamming = np.count_nonzero(xyz[:, None, :] != xyz[None, :, :], axis=2)
    require(
        np.array_equal(hop, hamming) and n == 2 ** len(cuts),
        "inferred cuts do not give a complete binary cube",
    )
    self_weight = np.trace(calibration) / calibration.sum()
    compressed = np.eye(n) * self_weight + (1 - self_weight) * edges / degree[0]
    return xyz, compressed


def terminal(panels, stopping):
    """Complete source rows: matched receiver, unmatched, clock timeout."""
    success = sum(
        weight * panel for weight, panel in zip(stopping, panels, strict=True)
    )
    miss = stopping.sum() - success.sum(axis=1)
    require(miss.min() > -1e-10, "scalar detector over-normalized")
    return np.column_stack(
        (success, np.maximum(miss, 0), np.full(len(success), 1 - stopping.sum()))
    )


def run():
    design = json.loads(BOOK.read_text())
    spec = json.loads((HERE / "clocked_protocol.json").read_text())["fixture"]
    host, data, _, front, bank = patches(spec["level"])
    weights, _ = clock_laws(host, front, bank, spec)
    comp = Compiler(data, spec["marker"], spec["shape"])
    n = int(np.prod(comp.shape[:3]))
    size = len(comp.slots) // n
    rng = np.random.default_rng(design["seed"])
    permutation = rng.permutation(n)
    # Only the apparatus sees slot group membership. The inference takes a count matrix.
    emit = np.repeat(np.eye(n), size, axis=0)[:, permutation] / np.sqrt(size)
    cap = spec["alarm_after_periods"]
    calibration = None
    physical = {}
    for arm in ("nominal", "shift_only"):
        columns = comp.columns(comp.build(arm))
        step = columns[comp.slots]
        require(
            np.linalg.norm(step.conj().T @ step - np.eye(len(step)), ord=2) < 1e-10,
            "marked-sector leakage",
        )
        state = emit
        panels = []
        for k in range(cap + 1):
            panels.append((abs(emit.T @ state) ** 2).T)
            if k == 1:
                position = (
                    (abs(state.reshape(n, size, n)) ** 2).sum(axis=1)[permutation].T
                )
                if arm == "nominal":
                    reference = position
                    calibration = np.array(
                        [
                            rng.multinomial(design["shots"], row / row.sum())
                            for row in position
                        ]
                    )
                else:
                    require(
                        np.max(abs(reference - position)) < 1e-12,
                        "first-step control differs",
                    )
            state = step @ state
        physical[arm] = panels
    coords, matrix = infer(calibration)
    matrices = [np.eye(n), matrix]
    for _ in range(2, cap + 1):
        matrices.append(2 * matrix @ matrices[-1] - matrices[-2])
    predicted = [abs(m) ** 2 for m in matrices]
    rows = []
    alias_floors = np.zeros(coords.shape[1])
    for r in (1, 2, 3):
        pmf, _ = stopping_law(r, cap)
        require(
            np.array_equal(pmf, weights[r]), "clock differs from physical first passage"
        )
        forecast = terminal(predicted, pmf)
        for arm, panels in physical.items():
            law = terminal(panels, pmf)
            counts = np.array(
                [rng.multinomial(design["shots"], f / f.sum()) for f in law]
            )
            if arm == "nominal":
                # Dropping a cut identifies two terminals. Any common prediction errs by at least half their difference.
                for axis in range(coords.shape[1]):
                    kept = np.delete(coords, axis, axis=1)
                    for u, v in zip(
                        *np.where(
                            np.triu(np.all(kept[:, None] == kept[None, :], axis=-1), 1)
                        ),
                        strict=True,
                    ):
                        alias_floors[axis] = max(
                            alias_floors[axis], max(abs(law[:, u] - law[:, v])) / 2
                        )
            rows.append(
                dict(
                    clock=r,
                    arm=arm,
                    detector_counts=counts.tolist(),
                    maximum_probability_error=float(np.max(abs(forecast - law))),
                    maximum_count_error=float(
                        np.max(abs(forecast - counts / design["shots"]))
                    ),
                )
            )
    return dict(
        protocol_sha256=sha256(BOOK.read_bytes()).hexdigest(),
        source_sha256={
            name: sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                "blind_chart.py",
                "compiler.py",
                "clocked_protocol.json",
                "clocked_apparatus.py",
                "clocked_field.py",
                "clocked_sources.json",
            )
        },
        calibration_counts=calibration.tolist(),
        inferred_coordinates=coords.tolist(),
        fitted_self_weight=float(matrix[0, 0]),
        removed_direction_error_lower_bounds=alias_floors.tolist(),
        held=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    path.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
