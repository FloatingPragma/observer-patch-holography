"""Measure the neighbouring receiver with its material hand clear of the gateway."""

import argparse
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.experiments.material_inertia.neutral_defect import ordered_index
from .measure import (
    HERE,
    NEAR,
    FAR,
    OTHER,
    acquire,
    branched,
    flight,
    native_gates,
    sample,
)


def endpoint_swap(car, x, y, a, b):
    route = list(range(car.n))
    route[a], route[b] = route[b], route[a]
    ports = np.asarray(route)
    return ordered_index(ports[x], ports[y], car.n)


def range_bank(h, car, y, tick, initial, turn, rng):
    raw = branched(h, car, y, tick, [initial], ranging=True, turn=turn)
    out = {}
    for deadline, bank in raw.items():
        keys = sorted(bank)
        probs = [float(bank[key][0]) for key in keys]
        counts = sample(probs, rng)
        certified = [
            i
            for i, key in enumerate(keys)
            if key.startswith("visited:") and not key.endswith((":12", ":timeout"))
        ]
        arrived = [i for i, key in enumerate(keys) if key.startswith("visited:")]
        out[str(deadline)] = dict(
            keys=keys,
            counts=counts.tolist(),
            returns=int(counts[certified].sum()),
            remote_marks=int(counts[arrived].sum()),
            return_probability=round(sum(probs[i] for i in certified), 13),
        )
    return out


def run():
    native_gates()
    car, x, y, gateway_tick, k, v = flight()
    h = k + v
    dist = shortest_path(car.adjacency(), directed=False, unweighted=True)[NEAR, :12]
    hand = int(np.argmax(dist))
    mate = min(i for i in car.neighbours[hand] if i < 12)
    tick = ((x == hand) | (x == mate)) & ~((y == hand) | (y == mate))
    rng = np.random.default_rng(305220)
    sources = [(hand, FAR), (hand, OTHER), (hand, NEAR), (hand, OTHER)]
    f = branched(h, car, y, tick, sources)
    result = dict(hand=[hand, mate], station_distance=int(dist[hand]))
    for name, cols in [("passive", [0, 1]), ("relay", [2, 3]), ("sham", [0, 0])]:
        result[name] = acquire(
            {d: {key: p[cols] for key, p in bank.items()} for d, bank in f.items()}, rng
        )
    result["half_exchange"] = acquire(
        branched(k * 0.5 + v, car, y, tick, sources[:2]), rng
    )
    for label, initial, ticks in [
        ("separated", (hand, NEAR), tick),
        ("gateway", (NEAR, 0), gateway_tick),
    ]:
        swap = endpoint_swap(car, x, y, initial[1], FAR)
        # The gateway calibration uses S(0,1) S(1,82) S(0,1). The direct
        # handshake at the second station uses the native pulse S(1,82).
        if FAR not in car.neighbours[NEAR] or 0 not in car.neighbours[NEAR]:
            raise ArithmeticError("missing turnaround contact")
        for service, turn in [("passive", None), ("turn", swap)]:
            result[f"range_{label}_{service}"] = range_bank(
                h, car, y, ticks, initial, turn, rng
            )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=str)
    args = parser.parse_args()
    data = run()
    if args.record:
        with open(args.record, "x") as f:
            json.dump(data, f, separators=(",", ":"))
            f.write("\n")
    elif data != json.loads((HERE / "separated.json").read_text()):
        raise SystemExit("separate-hand receipts differ")
    for name, rows in data.items():
        if isinstance(rows, dict):
            print(
                name,
                [
                    (
                        d,
                        b.get("correct"),
                        b.get("detected"),
                        b.get("returns"),
                        b.get("return_probability"),
                    )
                    for d, b in rows.items()
                ],
            )
