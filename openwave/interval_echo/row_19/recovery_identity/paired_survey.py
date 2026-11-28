"""Clocked radar triangles on the two written heads after unselected delivery."""

import argparse
from hashlib import sha256
from itertools import permutations
import json

from mpmath import iv
import numpy as np

from . import paired_relay as pair, relay, relay_precision as precision

RECEIPT = relay.HERE / "paired_survey.json"
SHOTS = 2**74
SEED = 191204
ROUTES = tuple(permutations(range(3), 2))
MENU = (*relay.NAMES, "blank_a", "blank_b")
# 10 return thresholds per route, four combinations of courier and probe.
L = np.log(4 * (len(MENU) * 6 * 10 * 4 + 2 * 12) / 0.01)


def effects(board, sites, paths, q):
    h = relay.generator(board, [(p[:2], 1) for p in paths], q)
    u = relay.RANGING.kernels(h, sites, 0.5, 11)
    first = {r: relay.RANGING.first_hits(u, *r) for r in ROUTES}
    panels = []
    for i, j in ROUTES:
        arms = []
        for n in range(2, 12):
            a = np.array([first[j, i][n - t] @ first[i, j][t] for t in range(1, n)])
            arms.append(np.einsum("taj,tak->jk", a.conj(), a))
        arms.append(np.eye(4) - sum(arms))
        panels.append(arms)
    panels = np.array(panels)
    if np.linalg.eigvalsh(panels).min() < -2e-12:
        raise ValueError("Nonpositive return or timeout effect")
    return panels


def acquire():
    delivery = json.loads(pair.RECEIPT.read_text())
    board = relay.FORMATION.federation(1)
    paths = delivery["paths"]
    sites = [*paths[0], paths[1][0]]
    rng = np.random.default_rng(SEED)
    clock, _ = relay.RANGING.clocks(0.5, 12)
    hand = [[pair.draws((x, 1 - x), SHOTS, rng) for x in row] for row in clock[1:]]
    banks = {}
    for q in (-1, 1):
        panel = effects(board, sites, paths, q)
        for i, sgn in enumerate((-1, 1)):
            states = [
                *pair.unpack(delivery["installed"][i]),
                np.diag([1, 0, 0, 0]),
                np.diag([1, 0, 0, 0]),
            ]
            laws = np.einsum("rbij,sji->srb", panel, states).real
            banks[f"{sgn},{q}"] = [
                [pair.draws(dist, SHOTS, rng) for dist in stock] for stock in laws
            ]
    data = dict(
        protocol_sha256=sha256(
            (relay.HERE / "paired_relay_protocol.md").read_bytes()
        ).hexdigest(),
        delivery_sha256=sha256(pair.RECEIPT.read_bytes()).hexdigest(),
        shots=SHOTS,
        seed=SEED,
        stocks=MENU,
        clock=hand,
        banks=banks,
    )
    data["measurements"] = decode(data)
    return data


def interval(lo, hi):
    return iv.mpf([lo, hi])


def limits(x):
    return [float(x.a), float(x.b)]


def length(bank):
    if (
        len(bank) != 11
        or any(type(m) is not int or m < 0 for m in bank)
        or sum(bank) != SHOTS
    ):
        raise ValueError("Incomplete radar histogram")
    # Twice the range in appointment units is the capped mean return bin.
    # A return on bin 2 is the earliest. Misses score a timeout at bin 12.
    total = 2 * SHOTS
    remaining = SHOTS
    error = 0.0
    for m in bank[:-1]:
        remaining -= m
        total += remaining
        error += precision.radius(remaining, L)
    return (iv.mpf(total) / SHOTS + interval(-error, error)) / 4


def clock_check(bank):
    angles = []
    counts = pair.checked(bank, (12, 2, 2), SHOTS)
    for tick, row in enumerate(counts, 1):
        moments = []
        for hits, _ in row:
            e = precision.radius(hits, L)
            moments.append(2 * (iv.mpf(hits) / SHOTS + interval(-e, e)) - 1)
        cos, sin = -moments[0], moments[1]
        angle = iv.atan2(sin, cos)
        if tick >= 7:
            angle += 2 * iv.pi
        if tick / 2 not in angle:
            raise ValueError("Registered clock step failed")
        angles.append(limits(angle))
    return angles


def area(bank):
    lengths = {r: length(b) for r, b in zip(ROUTES, bank, strict=True)}
    sides = [(lengths[i, j] + lengths[j, i]) / 2 for i, j in ((0, 1), (1, 2), (0, 2))]
    a, b, c = sides
    sq = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16
    if sq.a <= 0:
        raise ValueError("Triangle unresolved")
    return iv.sqrt(sq)


def decode(data):
    delivery = json.loads(pair.RECEIPT.read_text())
    if (
        data["delivery_sha256"] != sha256(pair.RECEIPT.read_bytes()).hexdigest()
        or data["protocol_sha256"] != delivery["protocol_sha256"]
        or data["shots"] != SHOTS
        or list(data["stocks"]) != list(MENU)
        or set(data["banks"]) != {"-1,-1", "-1,1", "1,-1", "1,1"}
    ):
        raise ValueError("Unexpected instrument or source")
    times = clock_check(data["clock"])
    result = {}
    for key, states in data["banks"].items():
        if len(states) != len(MENU):
            raise ValueError("Missing stock")
        a = [area(s) for s in states]
        ordinary = (a[1] - a[0]) / (2 * delivery["half_steps"]["ordinary"])
        held = (a[3] - a[2]) / (2 * delivery["half_steps"]["held"])
        ratio = None if 0 in ordinary else held / ordinary
        result[key] = dict(
            area_at_ordinary_minus=limits(a[0]),
            ordinary=limits(ordinary),
            held=limits(held),
            ratio=None if ratio is None else limits(ratio),
            source_compatible=None
            if ratio is None
            else 0 in ratio - iv.mpf(delivery["source_interval"]),
            formed_differences=[limits(a[i + 1] - a[i]) for i in (4, 6)],
            blank_control=limits(a[9] - a[8]),
        )
    return dict(angles=times, surveys=result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write", action="store_true")
    m.add_argument("--counts-only", action="store_true")
    args = p.parse_args()
    if args.counts_only:
        d = json.loads(RECEIPT.read_text())
        relay.tidal.verify(decode(d), d["measurements"])
    else:
        d = acquire()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(d, f, separators=(",", ":"), allow_nan=False)
                f.write("\n")
        else:
            relay.tidal.verify(d, json.loads(RECEIPT.read_text()))
    for name, read in d["measurements"]["surveys"].items():
        print(name, read["ratio"], read["source_compatible"], read["blank_control"])
