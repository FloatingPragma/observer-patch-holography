"""Read native contact work and ordinary-calibrated radar on the same source stocks."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import log

from mpmath import iv
import numpy as np

from . import paired_relay as pair, paired_survey as radar, relay

HERE = relay.HERE
PROTOCOL = HERE / "bridge_radar_protocol.md"
RECEIPT = HERE / "bridge_radar.json"
SHOTS = radar.SHOTS
SEED = 191207
STOCKS = list(relay.NAMES) + ["blank_a", "blank_b"]
# Ten survival counts for each directed radar, two work bins, two clock axes.
EVENTS = 2 * 2 * 2 * 10 * 6 * 10 + 2 * 10 * 2 + 12 * 2
L = log(4 * EVENTS / 0.01)
LABELS = ("ordinary", "paired", "formed_vacancy", "formed_excess")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def bridge_effects(w, ports):
    a, b = (w[i].reshape(-1, 4) for i in ports)
    out = []
    for amplitude in (a + b, a - b):
        out.append(amplitude.conj().T @ amplitude / 2)
    return np.array([*out, np.eye(4) - sum(out)])


def source_laws(board, route, paths, charge, stocks):
    laws = np.zeros((len(stocks), 3))
    for parity, wave in enumerate(pair.flights(board, route, paths, charge)):
        effects = bridge_effects(wave, route[1:3])
        laws += np.einsum("mij,sji->sm", effects, stocks[:, parity]).real
    return laws


def acquire():
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    paths = pair.paths(board, route, relay.tidal.chart(board, route)[0])
    delivery = json.loads(pair.RECEIPT.read_text())
    if [list(path) for path in paths] != delivery["paths"]:
        raise ValueError("Different receiving seams")
    stocks, steps, _, _ = relay.sources(board, route)
    rng = np.random.default_rng(SEED)
    clock, _ = relay.RANGING.clocks(0.5, 12)
    counts = [[pair.draws((x, 1 - x), SHOTS, rng) for x in row] for row in clock[1:]]
    contact, banks = {}, {}
    for i, courier in enumerate((-1, 1)):
        source = source_laws(board, route, paths, courier, stocks)
        contact[str(courier)] = [
            pair.draws(p, SHOTS, rng) for p in [*source, source[0], source[0]]
        ]
    for ix, path in enumerate(paths):
        sites = [*path, paths[1 - ix][0]]
        for probe in (-1, 1):
            panel = radar.effects(board, sites, paths, probe)
            for i, courier in enumerate((-1, 1)):
                states = [
                    *pair.unpack(delivery["installed"][i]),
                    np.diag([1, 0, 0, 0]),
                    np.diag([1, 0, 0, 0]),
                ]
                laws = np.einsum("rbij,sji->srb", panel, states).real
                banks[f"{ix},{courier},{probe}"] = [
                    [pair.draws(p, SHOTS, rng) for p in s] for s in laws
                ]
    data = dict(
        protocol=digest(PROTOCOL),
        delivery=digest(pair.RECEIPT),
        seed=SEED,
        shots=SHOTS,
        paths=paths,
        steps=steps,
        stocks=STOCKS,
        clock=counts,
        contact=contact,
        radar=banks,
    )
    data["results"] = decode(data)
    return data


def interval(count, factor=L):
    center = iv.mpf(int(count)) / SHOTS
    radius = radar.precision.radius(count, factor)
    return center + iv.mpf([-radius, radius])


def ranges(bank, factor=L):
    lengths = {}
    for ports, row in zip(radar.ROUTES, bank, strict=True):
        remaining = SHOTS
        time = iv.mpf(2)
        for m in row[:-1]:
            remaining -= m
            time += interval(remaining, factor)
        lengths[ports] = time / 4
    return lengths


def area(bank, factor=L):
    lengths = ranges(bank, factor)
    a, b, c = ((lengths[i, j] + lengths[j, i]) / 2 for i, j in ((0, 1), (1, 2), (0, 2)))
    h = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16
    if h.a <= 0:
        raise ValueError("Nonpositive radar triangle")
    return iv.sqrt(h)


def slopes(values, steps):
    scales = (2 * steps["ordinary"], 2 * steps["held"], 1, 1)
    return [
        (values[i + 1] - values[i]) / s
        for i, s in zip((0, 2, 4, 6), scales, strict=True)
    ]


def det3(rows):
    (a, b, c), (d, e, f), (g, h, i) = rows
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def ratio(a, b):
    return None if 0 in b else radar.limits(a / b)


def calibrated(source, geometry):
    products = [geometry[k] * source[0] - geometry[0] * source[k] for k in (1, 2, 3)]
    resolved = 0 not in source[0] and 0 not in geometry[0]
    return dict(
        source_slopes=list(map(radar.limits, source)),
        source_paired_ratio=ratio(source[1], source[0]),
        determinant=list(map(radar.limits, products)),
        compatible=[(0 in x) if resolved else None for x in products],
    )


def decode(data):
    delivery = json.loads(pair.RECEIPT.read_text())
    keys = {f"{i},{s},{p}" for i in (0, 1) for s in (-1, 1) for p in (-1, 1)}
    if (
        data["protocol"] != digest(PROTOCOL)
        or data["delivery"] != digest(pair.RECEIPT)
        or data["stocks"] != STOCKS
        or data["seed"] != SEED
        or data["shots"] != SHOTS
        or data["steps"] != delivery["half_steps"]
        or list(map(list, data["paths"])) != delivery["paths"]
        or set(data["radar"]) != keys
        or set(data["contact"]) != {"1", "-1"}
    ):
        raise ValueError("Wrong calibration, population or geometry")
    # The archived clock check uses a smaller union bound. Passing that
    # narrower interval also passes this acquisition's larger error bound.
    radar.clock_check(data["clock"])
    predictions = {}
    for key in sorted(keys):
        _, courier, _ = key.split(",")
        raw = pair.checked(data["contact"][courier], (10, 3), SHOTS)
        sym, antisym = ([interval(n) for n in raw[:, i]] for i in (0, 1))
        work = [8 * s + 9 * a for s, a in zip(sym, antisym, strict=True)]
        s, a, w = [slopes(values, data["steps"]) for values in (sym, antisym, work)]
        packets = pair.checked(data["radar"][key], (10, 6, 11), SHOTS)
        areas = [area(packet) for packet in packets]
        g = slopes(areas, data["steps"])
        minors = {}
        for inds in combinations(range(4), 3):
            minor = det3([(s[i], a[i], g[i]) for i in inds])
            minors[",".join(LABELS[i] for i in inds)] = radar.limits(minor)
        pred = iv.mpf(delivery["source_interval"])
        comparison = calibrated(w, g)
        predictions[key] = dict(
            geometry_slopes=list(map(radar.limits, g)),
            geometry_paired_ratio=ratio(g[1], g[0]),
            fixed_material=comparison,
            exchange_only=calibrated(a, g),
            scalar_compatible=(0 in g[1] - pred * g[0]) if 0 not in g[0] else None,
            symmetric_meter_minors=minors,
            symmetric_meter_excluded=any(
                lo > 0 or hi < 0 for lo, hi in minors.values()
            ),
            blank_area=radar.limits(areas[-1] - areas[-2]),
            blank_work=radar.limits(work[-1] - work[-2]),
        )
    return predictions


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--write", action="store_true")
    g.add_argument("--counts-only", action="store_true")
    args = p.parse_args()
    if args.counts_only:
        d = json.loads(RECEIPT.read_text())
        relay.tidal.verify(decode(d), d["results"])
    else:
        d = acquire()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(d, f, separators=(",", ":"), allow_nan=False)
                f.write("\n")
        else:
            relay.tidal.verify(d, json.loads(RECEIPT.read_text()))
    for k, v in d["results"].items():
        print(
            k,
            v["fixed_material"]["source_paired_ratio"],
            v["geometry_paired_ratio"],
            v["fixed_material"]["compatible"],
            v["symmetric_meter_excluded"],
        )
