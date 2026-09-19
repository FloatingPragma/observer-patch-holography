"""Canonical dyadic means and authenticated, two-writer provenance on W12."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARCHIVE = ROOT / "evidence/support_wiring_776"
STEPS = (1, 5, 30, 100, 300)


def canonical(x):
    return (json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, data):
    Path(path).write_bytes(canonical(data))


def preparation(n):
    return np.fromiter((hashlib.sha256(f"oph776:20260919:{i}".encode()).digest()[0] % 8
                        for i in range(n)), dtype=np.int64, count=n)


def carrier(seams):
    adjacency = np.zeros((12, 12), dtype=int)
    for a, b in seams:
        adjacency[a, b] = adjacency[b, a] = 1
    lap = 5 * np.eye(12) - adjacency
    w, v = np.linalg.eigh(lap)
    p = v[:, abs(w - (5 - np.sqrt(5))) < 1e-8]
    p = p @ p.T
    gram = 4 * p
    basis = next(idx for idx in itertools.combinations(range(12), 3)
                 if np.linalg.det(gram[np.ix_(idx, idx)]) > 1e-6)
    factor = np.linalg.cholesky(gram[np.ix_(basis, basis)])
    frame = np.linalg.solve(factor, gram[list(basis)]).T
    assert np.max(abs(frame @ frame.T - gram)) < 1e-12
    return lap, p, frame


def phases(geometry, wiring="w12"):
    used = [set() for _ in range(12)]
    groups = []
    for a, b in sorted(map(tuple, geometry["seams"])):
        color = next(k for k in itertools.count() if k not in used[a] | used[b])
        while len(groups) <= color:
            groups.append([])
        groups[color].append((a, b))
        used[a].add(color)
        used[b].add(color)
    offsets = 12 * np.arange(len(geometry["faces"]))
    result = [("intra", (np.asarray(pairs)[None, :, :] + offsets[:, None, None]).reshape(-1, 2))
              for pairs in groups]
    if wiring != "isolated":
        g = geometry[wiring]
        result.append(("glued", np.column_stack((12 * g[:, 0] + g[:, 1], 12 * g[:, 2] + g[:, 3]))))
    return result


def limbs(values):
    """Lossless little-endian three-limb nonnegative integers; fail on overflow."""
    assert all(0 <= v < 2 ** 192 for v in values)
    mask = (1 << 64) - 1
    return np.column_stack([np.asarray((values >> shift) & mask, dtype=np.uint64)
                            for shift in (0, 64, 128)])


def unlimbs(values):
    return (values[:, 0].astype(object) + (values[:, 1].astype(object) << 64)
            + (values[:, 2].astype(object) << 128))


def root_binding():
    return {"schema": "oph.support_wiring.trace.v1", "specification_sha256": sha(HERE / "SPECIFICATION.md"),
            "geometry_sha256": sha(HERE / "geometry/geometry.json"), "sweeps_per_level": 4,
            "levels": [3, 4, 5], "numerator_bits": 192, "law": "exact_pair_mean",
            "events": "two simultaneous destination writes per mean action"}


def build_trace(output):
    output.mkdir(parents=True, exist_ok=True)
    binding = root_binding()
    chain = hashlib.sha256(canonical(binding)).hexdigest()
    chunks = []
    event_id, exponent, phase_id = 0, 0, 0
    state = writers = None
    for level in (3, 4, 5):
        g = dict(np.load(HERE / f"geometry/geometry_l{level}.npz"))
        _, _, frame = carrier(g["seams"])
        n = 12 * len(g["faces"])
        if level == 3:
            state = preparation(n).astype(object)
            incoming = np.full(n, -1, dtype=np.int32)
            kind = "prepare"
        else:
            source = (g["parent"][:, None] * 12 + np.arange(12)).reshape(-1)
            state = state[source].copy()
            incoming = writers[source]
            kind = "copy"
        writers = np.arange(event_id, event_id + n, dtype=np.int32)

        def emit(kind, sweep, data, count):
            nonlocal event_id, phase_id, chain
            data["positions"] = (np.asarray(state, dtype=float).reshape(-1, 12) / 2.0 ** exponent) @ frame
            data["port0_writers"] = writers[::12]
            filename = f"phase_{phase_id:03d}.npz"
            np.savez_compressed(output / filename, **data)
            item = {"file": filename, "kind": kind, "level": level, "sweep": sweep,
                    "phase": phase_id, "exponent": exponent, "first_event": event_id,
                    "event_count": count, "sha256": sha(output / filename), "previous": chain}
            chain = hashlib.sha256(canonical(item)).hexdigest()
            item["chain"] = chain
            chunks.append(item)
            event_id += count
            phase_id += 1

        emit(kind, 0, {"parents": incoming, "values": limbs(state)}, n)
        for sweep in range(1, 5):
            for kind, endpoints in phases(g):
                a, b = endpoints.T
                read_from = np.column_stack((writers[a], writers[b]))
                result = state[a] + state[b]
                state *= 2
                state[a] = result
                state[b] = result
                exponent += 1
                writers[a] = event_id + 2 * np.arange(len(a), dtype=np.int32)
                writers[b] = event_id + 2 * np.arange(len(a), dtype=np.int32) + 1
                emit(kind, sweep, {"endpoints": endpoints.astype(np.int32), "parents": read_from,
                                   "values": limbs(result)}, 2 * len(a))
            print(f"trace L{level} sweep {sweep}: {event_id:,} writes, denominator 2^{exponent}", flush=True)
    manifest = {"binding": binding, "chunks": chunks, "final_chain": chain, "events": event_id}
    save_json(output / "trace.json", manifest)
    return manifest


def laplacian(g, wiring):
    ep = np.concatenate([e for _, e in phases(g, wiring)])
    n = 12 * len(g["faces"])
    a, b = ep.T
    adjacency = sparse.coo_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()
    return sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel()) - adjacency, 2 * len(ep) / len(g["faces"])


def controls(output):
    rows, matrices = [], {}
    for level in (3, 4, 5):
        g = dict(np.load(HERE / f"geometry/geometry_l{level}.npz"))
        n = len(g["faces"])
        initial = preparation(12 * n).astype(float)
        _, projector, _ = carrier(g["seams"])
        degree = np.bincount(g["w12"][:, [0, 2]].ravel(), minlength=n)
        defect, regular = np.flatnonzero(degree == 11), np.flatnonzero(degree == 12)
        samples = sorted(set([int(defect[0]), int(defect[-1]), int(regular[0]), int(regular[-1])]))
        for wiring in ("isolated", "w3", "w12"):
            ordered = phases(g, wiring)
            mean = (np.repeat(initial.reshape(-1, 12).mean(axis=1), 12) if wiring == "isolated"
                    else np.full(len(initial), initial.mean()))
            endpoints, diagnostics = [], []
            for reverse in (False, True):
                x = initial.copy()
                v = [float(np.sum((x - mean) ** 2))]
                for _ in range(4):
                    for _, ep in (list(reversed(ordered)) if reverse else ordered):
                        a, b = ep.T
                        y = (x[a] + x[b]) / 2
                        x[a] = x[b] = y
                    v.append(float(np.sum((x - mean) ** 2)))
                endpoints.append(x)
                diagnostics.append({"reverse": reverse, "squared_distance_by_sweep": v,
                                    "max_distance_to_component_mean": float(np.max(abs(x - mean))),
                                    "sum_drift": float(x.sum() - initial.sum()),
                                    "max_component_sum_drift": float(np.max(abs((x-initial).reshape(-1, 12).sum(axis=1))))
                                    if wiring == "isolated" else float(abs(x.sum()-initial.sum()))})
            row = {"level": level, "wiring": wiring, "confluence": diagnostics,
                   "finite_schedule_max_difference": float(np.max(abs(endpoints[0] - endpoints[1]))), "kernels": []}
            lap, denominator = laplacian(g, wiring)
            operator = sparse.eye(12 * n, format="csr") - lap / denominator
            if wiring == "isolated":
                operator = operator[:12, :12]
            q = np.eye(12) - np.ones((12, 12)) / 12
            for cell in ([0] if wiring == "isolated" else samples):
                block = slice(12 * cell, 12 * cell + 12)
                y = np.zeros((operator.shape[0], 12))
                y[block] = q
                for step in range(1, 2 * max(STEPS) + 1):
                    y = operator @ y
                    scale = np.max(abs(y))
                    if scale:
                        y /= scale
                    if step % 2 == 0 and step // 2 in STEPS:
                        kernel = q @ y[block]
                        kernel = (kernel + kernel.T) / 2
                        kernel *= 12 / np.trace(kernel)
                        key = f"l{level}_{wiring}_c{cell}_n{step//2}"
                        matrices[key] = kernel
                        row["kernels"].append({"cell": cell, "n": step // 2, "matrix": key,
                                               "slow_share": float(np.trace(projector @ kernel) / 12),
                                               "eigenvalues": np.linalg.eigvalsh(kernel)[::-1].tolist()})
                print(f"kernels L{level} {wiring} cell {cell} done", flush=True)
            row["denominator"] = denominator
            rows.append(row)
    np.savez_compressed(output / "kernels.npz", **matrices)
    save_json(output / "controls.json", rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["trace", "controls"])
    parser.add_argument("--output", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mode == "trace":
        build_trace(args.output / "trace")
    else:
        controls(args.output)
