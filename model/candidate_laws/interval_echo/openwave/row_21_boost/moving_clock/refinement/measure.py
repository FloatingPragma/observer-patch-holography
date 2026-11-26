"""Material readout on three geometrically selected guides at the same velocity."""

import argparse
import json
from pathlib import Path

import numpy as np

from .. import apparatus as lab, preparation, readback

HERE = Path(__file__).resolve().parent
SHOTS = 2**34
SIZES = (27, 54, 108)
PHASES = (1, 1j, -1, -1j)


def guides(board):
    path, pos, seen = [0], {0: 0}, {0}
    stack = [iter(board.neighbours[0])]
    found = {}
    while stack:
        v = next(stack[-1], None)
        if v is None:
            stack.pop()
            pos.pop(path.pop())
        elif v not in seen:
            seen.add(v)
            pos[v] = len(path)
            path.append(v)
            stack.append(iter(board.neighbours[v]))
        elif v in pos:
            route = path[pos[v] :]
            n = len(route) // 2
            if n not in SIZES or 2 * n != len(route) or n in found:
                continue
            if len({x // 12 for x in route}) < 3:
                continue
            for j in range(len(route)):
                cycle = route[j:] + route[:j]
                outs = [
                    sorted(set(board.neighbours[x]) - set(cycle)) for x in cycle[:2]
                ]
                if all(outs):
                    found[n] = dict(ports=cycle, vacant=[z[0] for z in outs])
                    break
            else:
                raise ValueError("First geometry-selected route lacks origin receivers")
    if set(found) != set(SIZES):
        raise ValueError("Incomplete refinement triple")
    return found


def clock(arms, layers, stations, rail, rng, pilot=False):
    counts, mixed, pilots = [], [], []
    for t, where in enumerate(stations):
        a, b = arms[:, where, rail]
        counts.append(
            [int(rng.binomial(SHOTS, abs(a + q * b) ** 2 / 2)) for q in PHASES]
        )
        if pilot:
            p = (abs(a) ** 2 + abs(b) ** 2) / 2
            mixed.append(rng.binomial(SHOTS, p, size=4).tolist())
            if t in (0, len(stations) - 1):
                pilots.append(
                    [
                        lab.draw(abs(arm) ** 2, SHOTS, rng).reshape(arm.shape).tolist()
                        for arm in arms
                    ]
                )
        arms = lab.step(arms, layers)
    return dict(occupied=counts, mixture=mixed, pilot_full=pilots)


def panel(n, h, s, rng):
    dt = 10 * n // 27
    d = int(np.sign(h)) * 4 * dt // 5
    delta = (3 / 4) * (6 * np.pi / n)
    holds = (np.pi, np.pi - 2 * delta)
    layers = lab.gates(s, 8, holds=holds)
    rail = int(h < 0)
    stations = np.rint(np.linspace(0, d, dt + 1)).astype(int) % n
    result = {}
    source_error = 0.0
    for plane in (True, False):
        arms = lab.arms(n, h, layers, s, 8, plane=plane, holds=holds)
        for state in [*arms, *((arms[0] + q * arms[1]) / np.sqrt(2) for q in PHASES)]:
            source_error = max(source_error, preparation.build(state, s)[1])
        if plane:
            rails = [
                [
                    lab.draw((abs(ket) ** 2).sum(axis=0), SHOTS, rng).tolist()
                    for ket in (arm, arm @ layers[0].T, lab.step(arm, layers))
                ]
                for arm in arms
            ]
        result["plane" if plane else "packet"] = clock(
            arms, layers, stations, rail, rng, pilot=not plane
        )
    echo = lab.radar(n, 0, d, layers, SHOTS, rng, turn=dt)
    hand = lab.reference_hand(n, [echo[x] for x in ("emission", "turn", "receipt")])
    return dict(
        n=n,
        h=h,
        exchange=s,
        delta=delta,
        dt=dt,
        d=d,
        stations=stations.tolist(),
        rail=rail,
        source_error=source_error,
        arm_current=rails,
        echo=echo,
        hand=hand,
        **result,
    )


def acquire():
    rng = np.random.default_rng(211944)
    board = lab.federation(2)
    rings = guides(board)
    for route in rings.values():
        ports = route["ports"]
        assert all(
            b in board.neighbours[a]
            for a, b in zip(ports, ports[1:] + ports[:1], strict=True)
        )
        route["private_phase_port"] = next(p for p in ports if board.degree[p] == 5)
    return dict(
        seed=211944,
        shots=SHOTS,
        guides=rings,
        rows=[
            panel(n, h, s, rng) for n in SIZES for s in (1.0, 0.8) for h in (0, 3, -3)
        ],
    )


def ratio(a, ref):
    if not a["resolved"] or not ref["resolved"]:
        return None
    x, dx = a["elapsed"], a["phase_error"]
    y, dy = ref["elapsed"], ref["phase_error"]
    if min(x - dx, y - dy) <= 0:
        return None
    return [float((x - dx) / (y + dy)), float(x / y), float((x + dx) / (y - dy))]


def read(data):
    shots = data["shots"]
    rows, rest, shortfalls, contractions = [], {}, {}, []
    for row in data["rows"]:
        n, h, s = row["n"], row["h"], row["exchange"]
        e, t, r = row["hand"]
        if [e, t, r] != [t - abs(row["d"]), row["dt"], t + abs(row["d"])]:
            raise ValueError("Unresolved radar dates")
        echo = row["echo"]
        if any(
            sum(echo[key]) != shots for key in ("forward_then_return", "parked_turn")
        ):
            raise ValueError("Incomplete echo trials")
        observed = np.sqrt(e * r) / t
        raw = {
            kind: readback.phase_bins(row[kind]["occupied"], shots)
            for kind in ("plane", "packet")
        }
        if h == 0:
            rest[n, s] = raw
        compared = {kind: ratio(raw[kind], rest[n, s][kind]) for kind in raw}
        mixed = readback.phase_bins(row["packet"]["mixture"], shots)
        pilot = np.asarray(row["packet"]["pilot_full"])
        if np.any(pilot.sum(axis=(2, 3)) != shots):
            raise ValueError("Incomplete pilot")
        centers = [readback.center(p.sum(axis=(0, 2))) for p in pilot]
        c = np.asarray(row["arm_current"], dtype=float)
        if np.any(c.sum(axis=2) != shots):
            raise ValueError("Incomplete current arms")
        v = (c[:, 2, 0] - c[:, 1, 0]).mean() / shots
        dv = readback.uncertainty(c[:, [2, 1], 0], shots).sum() / 2
        rows.append(
            dict(
                n=n,
                h=h,
                exchange=s,
                radar=observed,
                ratios=compared,
                centers=centers,
                current=[v - dv, v, v + dv],
                mixture_resolved=mixed["resolved"],
                echo=echo["forward_then_return"][-1],
                parked=echo["parked_turn"][-1],
                early=echo["early_clicks"],
                source_error=row["source_error"],
                clock_resolved={k: raw[k]["resolved"] for k in raw},
            )
        )
        if h and compared["plane"]:
            lo, mid, hi = compared["plane"]
            gap = [observed - hi, observed - mid, observed - lo]
            last = shortfalls.get((h, s))
            if last is not None and min(*gap, *last) > 0:
                contractions.append(
                    dict(
                        n=n,
                        h=h,
                        exchange=s,
                        bounds=[gap[0] / last[2], gap[1] / last[1], gap[2] / last[0]],
                    )
                )
            shortfalls[h, s] = gap
    return dict(rows=rows, contractions=contractions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", type=Path)
    parser.add_argument("--read", type=Path)
    args = parser.parse_args()
    if args.record:
        args.record.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    elif args.read:
        print(
            json.dumps(
                read(json.loads(args.read.read_text())), indent=1, allow_nan=False
            )
        )
