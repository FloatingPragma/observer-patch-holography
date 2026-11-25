"""Impulse delivered by a bridge to each internal center-of-energy velocity."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from .. import bridge

HERE = Path(__file__).resolve().parent
UNIT = 4096
RULERS = ("hops", "front")


def marked_block(route, slot, q, odd):
    blocks = 2 * bridge.native.iso.ranging.parity_blocks(
        q, 8, -1 if odd and slot == 0 else 1
    )
    np.testing.assert_array_equal(blocks, blocks.astype(np.int64))
    return blocks.astype(np.int64)


def coordinates(board, route):
    """Two frozen charts, in 1/4096 blank-link units, from each gateway inward."""
    st = bridge.action.preparations()["equilibrium"]
    charts = {}
    marks = (route[:2], route[2:])
    for slot, seam in enumerate(marks):
        tag = ("left", "right")[slot]
        gateway = route[1 + slot]
        offset = (gateway // 12) * 12
        start = gateway - offset
        obs = []
        for odd in (False, True):
            t = marked_block(route, slot, 1, odd)[1, 2]
            t = np.kron(t, np.eye(2)) if slot == 0 else np.kron(np.eye(2), t)
            obs.append((t @ t.T + t.T @ t) / 2)
        coefficient = np.einsum("bij,bji->", obs, st).real
        # Trace-norm source radius 6e-7 and the operator-norm bound by row sum.
        err = max(np.linalg.norm(x, ord=np.inf) for x in obs) * 6.01e-7
        lo, hi = UNIT / np.sqrt(coefficient + err), UNIT / np.sqrt(coefficient - err)
        label = int(np.rint(lo))
        assert label == int(np.rint(hi)) and lo > label - 0.499 and hi < label + 0.499
        edge = tuple(sorted((p - offset for p in seam)))
        g = (
            board.adjacency()[offset : offset + 12, offset : offset + 12]
            .toarray()
            .astype(np.int64)
            * UNIT
        )
        costs = {}
        for ruler in RULERS:
            if ruler == "front":
                g[edge] = g[edge[::-1]] = label
            x = shortest_path(g, directed=False, indices=start)
            np.testing.assert_array_equal(x, x.astype(np.int64))
            costs[ruler] = x.astype(np.int64)
        charts[tag] = dict(
            gateway=gateway,
            offset=offset,
            marked_length_numerator=label,
            marked_length_enclosure=[float(lo / UNIT), float(hi / UNIT)],
            x=costs,
        )
    return charts


def moment(board, route, cell, slot, q, odd, x):
    """Integer matrix of 4*UNIT*B. Share every internal seam equally at its two ports."""
    offset = 12 * (route[1 + slot] // 12)
    endpoints = sorted((route[:2], route[2:])[slot])
    marker = marked_block(route, slot, q, odd)
    units = np.repeat(x, 4)
    result = units[:, None] * cell + cell * units[None, :]
    identity = np.eye(4, dtype=np.int64)
    for a, b in board.seams:
        if a // 12 != offset // 12 or b // 12 != offset // 12:
            continue
        for p, t in ((a, b), (b, a)):
            diag = 17 * identity
            if [a, b] == endpoints:
                u = endpoints.index(p) + 1
                part = marker[u, u] - marker[0, 0]
                eye = np.eye(2, dtype=np.int64)
                diag = np.kron(part, eye) if slot == 0 else np.kron(eye, part)
            i, j = p - offset, t - offset
            result[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] += (x[j] - x[i]) * diag
    np.testing.assert_array_equal(result, result.T)
    return result


def setup(board, route, q, odd, charts):
    ports = [p for c in charts.values() for p in range(c["offset"], c["offset"] + 12)]
    idx = bridge.beam.indices(ports)
    b = bridge.bridges(board, [route[1:3]])[idx, :][:, idx].toarray()
    h = bridge.operators(board, route, q, odd)
    zero = np.zeros((48, 48), dtype=np.int64)
    families = {}
    for slot, (tag, chart) in enumerate(charts.items()):
        cell = h[chart["offset"] // 12]
        act = block_diag(cell, zero) if slot == 0 else block_diag(zero, cell)
        for name, x in chart["x"].items():
            inner = moment(board, route, cell, slot, q, odd, x)
            dipole = block_diag(inner, zero) if slot == 0 else block_diag(zero, inner)
            k = act @ dipole - dipole @ act
            force = k @ b - b @ k
            np.testing.assert_array_equal(force, force.T)
            # Shifting the energy origin by constant*cell energy fixes this momentum and force.
            translated = moment(board, route, cell, slot, q, odd, x + 3 * UNIT)
            np.testing.assert_array_equal(translated - inner, 6 * UNIT * cell)
            families[tag, name] = (force, b, act, dipole)
    return ports, families


def impulse_effects(w, mats, step):
    """Richardson mixed difference of interrupted readings with forward pulses."""
    _, b, h, dipole = mats
    ev, basis = np.linalg.eigh(h / 2)
    eb, ub = np.linalg.eigh(b / 2)
    d = dipole / (4 * UNIT)
    out = np.zeros((4, 4), dtype=complex)
    for size, weight in ((step, -1.0), (step / 2, 2.0)):
        energy_step = (basis * np.exp(-1j * size * ev)) @ basis.T
        bridge_step = (ub * np.exp(-1j * size * eb)) @ ub.T
        # Subtract the independently interrupted readings before and after
        # the internal pulse to avoid subtracting four large state means.
        delta = energy_step.conj().T @ d @ energy_step - d
        response = bridge_step.conj().T @ delta @ bridge_step - delta
        out += (
            weight
            * sum(arm.conj().T @ response @ arm for arm in np.hsplit(w, 10))
            / (10 * size**2)
        )
    return out


def run():
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    charts = coordinates(board, route)
    rings = bridge.beam.screens(board, route)
    states = bridge.tangents()
    target = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    results = {}
    for q in (-1, 1):
        observations = {}
        for odd in (False, True):
            ports, families = setup(board, route, q, odd, charts)
            h = bridge.native.generator(board, route, q, odd)
            ht = -2j * h
            start = np.zeros((h.shape[0], len(target)))
            start[target] = np.eye(len(target))
            w = expm_multiply(ht, start, traceA=ht.diagonal().sum())[
                bridge.beam.indices(ports)
            ]
            for key, mats in families.items():
                force = mats[0] / (16 * UNIT)
                arms = [arm.conj().T @ force @ arm for arm in np.hsplit(w, 10)]
                bucket = observations.setdefault(
                    "_".join(key),
                    {"commutator": [], "pulse_1_512": [], "pulse_1_1024": []},
                )
                bucket["commutator"].append(np.mean(arms, axis=0))
                for n in (512, 1024):
                    bucket[f"pulse_1_{n}"].append(impulse_effects(w, mats, 1 / n))
        result = {
            key: {
                method: bridge.evaluate(blocks, states)
                for method, blocks in sample.items()
            }
            for key, sample in observations.items()
        }
        results[str(q)] = result
    for chart in charts.values():
        chart["x"] = {k: v.tolist() for k, v in chart["x"].items()}
    return dict(unit=UNIT, time=4, rulers=charts, force=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
