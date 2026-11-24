"""Independently prepared Ramsey and mobile qutrit probes at the poles of a wall belt."""

import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt
from model.executable_model import carrier
from model.experiments.interval_echo_encounter.conjugation import free_seam

HERE = Path(__file__).resolve().parent
SHOTS = 262144
SEED = 1211248
SHORT = (1 / 64, 1 / 32)
LONG = (1.0, 2.0)
WALL = {3: (0, 4, 7), 4: (0, 1, 4, 7)}


def geometry():
    board = carrier.federation(0)
    for poles, sites in sorted(belt.belts().items()):
        bridges = [(p, v) for p in poles for v in board.neighbours[p] if v >= 12]
        if not bridges:
            continue
        pole, far = bridges[0]
        other = next(p for p in poles if p != pole)
        near = [i for i, v in enumerate(sites) if v in board.neighbours[pole]]
        opposite = [i for i, v in enumerate(sites) if v in board.neighbours[other]]
        assert len(near) == 5 and set(near) == set(range(10)) - set(opposite)
        assert set((i + 1) % 10 for i in near) == set(opposite)
        assert set(board.neighbours[far]).isdisjoint(sites)
        assert board.degree[pole] == board.degree[far] == 6
        return dict(
            belt=list(sites), near=pole, opposite=other, far=far, near_indices=near
        )
    raise ArithmeticError("No polar gateway in first cell")


def contact(strength, sign):
    """One fresh probe, relative to a host contact, from the registered seam matrix."""
    seam = free_seam(strength)

    def energy(host, probe):
        index = law.SRC.index((host, probe))
        return float(seam[index, index])

    empty_shift, filled_shift = [energy(h, 1 + sign) - energy(h, 1) for h in (1, 2)]
    swap = float(seam[law.SRC.index((1, 1 + sign)), law.SRC.index((1 + sign, 1))])
    return empty_shift, filled_shift - empty_shift, swap


def translated(word, offset):
    return tuple(sorted((v + offset) % 10 for v in word))


@lru_cache(maxsize=None)
def sector(n):
    return BinarySector(belt.ring(10), n)


def evolve(h, state, time):
    return expm_multiply((-1j * time) * h, state) if time else state.copy()


def preparation(kind, n, h, shift=0, sense=1):
    box = sector(n)
    word = WALL[n] if kind in {"sharp", "flight"} else tuple(range(n))
    state = np.zeros(len(box.states), complex)
    state[box.index[translated(word, shift)]] = 1
    if kind == "flight":
        state[box.index[translated(word, shift + 3)]] = sense * 1j
        state = evolve(h, state / np.sqrt(2), 1)
    return state


def detector(n):
    box = sector(n)
    labels = np.full(len(box.states), 10)
    if n in WALL:
        for j in range(10):
            labels[box.index[translated(WALL[n], j)]] = j
    return labels


def group(labels, probability):
    return np.bincount(labels, weights=probability, minlength=11)


def rounded(values):
    return np.round(values, 12).tolist()


def record(weights, rng):
    probability = np.asarray(weights, dtype=float)
    if np.min(probability) < -1e-13 or abs(probability.sum() - 1) > 1e-10:
        raise ArithmeticError("Instrument lost probability")
    prob = np.maximum(probability, 0)
    draws = rng.multinomial(SHOTS, (prob / prob.sum()).ravel())
    return dict(probabilities=rounded(prob), counts=draws.reshape(prob.shape).tolist())


def common(n, gain, exposure, rng):
    # Ordered arms (reference, sample); their visibility is <reference|sample>.
    phase = np.exp(-1j * gain * n * exposure)
    return [record([(1 + v) / 2, (1 - v) / 2], rng) for v in (phase.real, phase.imag)]


def differential(state, h, occ, n, gain, duration, labels, rng):
    arms = [
        evolve(h + sparse.diags(gain * count), state, duration)
        for count in (occ, n - occ)
    ]
    visibility = np.vdot(*arms)
    frames = []
    for combiner in (1, -1j):
        probabilities = [
            group(labels, abs((arms[0] + s * combiner * arms[1]) / 2) ** 2)
            for s in (1, -1)
        ]
        frames.append(record(probabilities, rng))
    return dict(
        visibility=rounded([visibility.real, visibility.imag]), quadratures=frames
    )


def rail(state, h, occ, gain, swap, exposure, bridge, labels, rng):
    # The near and far endpoint have the same five-contact host self energy.
    # Keeping the bridge diagonal during a Zeno hold adds only a scalar on this sector.
    size = len(state)
    off = sparse.eye(size, format="csr") * (swap if bridge else 0.0)
    joint_h = sparse.bmat([[h + sparse.diags(gain * occ), off], [off, h]], format="csr")
    outgoing = evolve(joint_h, np.concatenate((np.zeros(size), state)), exposure)
    probabilities = [group(labels, abs(a) ** 2) for a in outgoing.reshape(2, size)]
    return record(probabilities, rng)


def inputs():
    for n in (0, 1):
        yield n, "calibration", 0, 1
    for n in (3, 4):
        yield n, "compact", 0, 1
        yield n, "sharp", 0, 1
        for shift in (0, 1):
            for sense in (-1, 1):
                yield n, "flight", shift, sense


def run():
    rng = np.random.default_rng(SEED)
    geom = geometry()
    rows = []
    for strength in (8, 0):
        for n, kind, shift, sense in inputs():
            box = sector(n)
            h = sparse.csr_matrix(box.hamiltonian(strength))
            state = preparation(kind, n, h, shift, sense)
            occ = np.array(
                [len(set(w).intersection(geom["near_indices"])) for w in box.states]
            )
            labels = detector(n)
            entry = dict(
                kappa=strength,
                loads=n,
                kind=kind,
                shift=shift,
                sense=sense,
                arrival_p=rounded(group(labels, abs(state) ** 2)),
                probes=[],
            )
            for sign in (-1, 1):
                self_shift, gain, swap = contact(strength, sign)
                probe = dict(
                    sign=sign,
                    host_self=self_shift,
                    incremental_contact=gain,
                    bridge_swap=swap,
                    common=[common(n, gain, t, rng) for t in SHORT],
                    differential=[
                        differential(state, h, occ, n, gain, t, labels, rng)
                        for t in SHORT
                    ],
                    mobile=[
                        rail(state, h, occ, gain, swap, t, opened, labels, rng)
                        for t in LONG
                        for opened in (False, True)
                    ],
                )
                entry["probes"].append(probe)
            rows.append(entry)
    return dict(
        seed=SEED,
        shots=SHOTS,
        short_times=SHORT,
        long_times=LONG,
        geometry=geom,
        rows=rows,
        bulk_reference=dict(input_cell=[2, 1, 1], host=1, excess_per_port=1 / 3),
    )


def infer(common_read):
    cos, sin = [np.diff(frame["counts"][::-1])[0] / SHOTS for frame in common_read]
    return np.arctan2(sin, cos)


def metrics(result):
    def find(k, n, kind, shift=0, sense=1):
        return next(
            r
            for r in result["rows"]
            if (r["kappa"], r["loads"], r["kind"], r["shift"], r["sense"])
            == (k, n, kind, shift, sense)
        )

    summaries = []
    for qindex, sign in enumerate((-1, 1)):
        for time_index, t in enumerate(SHORT):

            def phase(row):
                return infer(row["probes"][qindex]["common"][time_index])

            base = phase(find(8, 0, "calibration"))
            unit = phase(find(8, 1, "calibration")) - base
            for n in (3, 4):
                flight = find(8, n, "flight")
                summaries.append(
                    dict(
                        sign=sign,
                        exposure=t,
                        loads=n,
                        flight_bulk_excess=(phase(flight) - base) / unit - 10 / 3,
                        compact_bulk_excess=(phase(find(8, n, "compact")) - base) / unit
                        - 10 / 3,
                    )
                )
    measures = []
    for n in (0, 1, 3, 4):
        kinds = ("calibration",) if n < 3 else ("sharp", "flight", "compact")
        for kind in kinds:
            row = find(8, n, kind)
            for probe in row["probes"]:
                mobile = probe["mobile"]
                measures.append(
                    dict(
                        loads=n,
                        kind=kind,
                        sign=probe["sign"],
                        near_when_open=rounded(
                            [np.sum(mobile[i]["probabilities"][0]) for i in (1, 3)]
                        ),
                        near_when_held=rounded(
                            [np.sum(mobile[i]["probabilities"][0]) for i in (0, 2)]
                        ),
                    )
                )
    return dict(calibrated_common=summaries, near_endpoint=measures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    result["summary"] = metrics(result)
    path = HERE / "results.json"
    if args.write:
        path.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    elif json.loads(path.read_text()) != json.loads(json.dumps(result)):
        raise SystemExit("Native polar-probe replay differs from archived record")
    print(
        f"{len(result['rows'])} preparations verified, {SHOTS} trials per detector setting"
    )


if __name__ == "__main__":
    main()
