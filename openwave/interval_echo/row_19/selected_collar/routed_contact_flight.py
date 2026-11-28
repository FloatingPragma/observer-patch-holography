"""Release two material probes across an exterior edge of a recorded-walk source."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import routed_contacts as routes

local = routes.local
HERE = local.HERE
N = 2**56
SEED = 28122892


def instrument():
    profile = next(
        v
        for v in local.inputs("routed_contacts")["profiles"]
        if v["stations"] == 20 and v["stop_bits"] == 5
    )
    board = local.survey.bulk.federation(0)
    steps, _, _ = routes.layout(board)
    match = routes.receiver_codes(board, steps)
    return native_instrument(board, profile, match)


def native_instrument(board, profile, match):
    u, v = profile["held_edge"]
    base = set(local.survey.routes.codebook()[0][0])
    ports = [12 * c + i for c in (u, v) for i in range(12)]
    loc = {z: i for i, z in enumerate(ports)}
    a, b = [loc[z] for z in profile["probe_ports"]]
    s1, s2 = [loc[z] for z in profile["ordinary_partners"]]
    degrees = np.array(board.degree)[ports]
    sites = (s1, s2)
    edges = []
    for i, j in board.seams:
        if i not in loc or j not in loc:
            continue
        if i // 12 != j // 12:
            if {loc[i], loc[j]} != {a, b}:
                continue
        elif (i % 12, j % 12) in base:
            continue
        edges.append(tuple(sorted((loc[i], loc[j]))))
    second = [
        {tuple(sorted((loc[12 * c + i], loc[12 * c + j]))) for i, j in match[c]}
        for c in (u, v)
    ]
    return profile, (u, v), ports, (a, b), sites, degrees, edges, second


def flight(serial):
    profile, (u, v), ports, (a, b), sites, degrees, edges, second = instrument()
    relay = local.inputs("shared_contact")["relay"]["serial" if serial else "0"]
    table = {}
    stock = {}
    raw = profile["joint_parity_and_stop"]
    rng = np.random.default_rng(SEED)
    for tag, writer, mass in (
        ("no_reports", "fresh", 0),
        ("fresh", "fresh", 0),
        ("shared", "shared", 0),
        ("ordinary_one", "fresh", 1),
        ("ordinary_two", "fresh", 2),
    ):
        p = F(relay[writer]["probability"])
        if tag == "no_reports":
            p = F(0)
        weights = defaultdict(F)
        for i, j, k in np.ndindex(2, 2, 3):
            for flip, r in enumerate((1 - p, p)):
                near, far = (
                    (i * flip, j * flip) if serial else (1 - i * flip, 1 - j * flip)
                )
                m = mass if k < 2 else 0
                key = (near, far, m, k if m else 0)
                weights[key] += F(raw[i][j][k]) * r
        assert sum(weights.values()) == 1
        bank = []
        deriv = [F(0)] * 5
        for key, w in weights.items():
            if not w:
                continue
            near, far, m, k = key
            if key not in stock:
                allowed = [
                    e
                    for e in edges
                    if not (near and e in second[0] or far and e in second[1])
                ]
                stock[key] = [
                    local.spectrum(ports, allowed, degrees, (sites[k], a, b), m, q)
                    for q in (1, 2)
                ]
            pair = stock[key]
            bank.append((w, pair))
            A, B = pair
            d = B["clock_diagonal"] - A["clock_diagonal"]
            for h, z in enumerate(
                (
                    A["curvature"],
                    B["curvature"],
                    F(int(d[0] - d[1]), 2),
                    A["plus_curvature"],
                    B["plus_curvature"],
                )
            ):
                deriv[h] += w * z
        table[tag] = {
            "derivatives": list(map(str, deriv[:3])),
            "in_phase_derivatives": list(map(str, deriv[3:])),
            "reads": local.timed_reads(bank, rng, N),
        }
    inc, doses = local.calibrations(table)
    gradient, outer = map(F, profile["exact_odd_on_held_edge"])
    du = gradient - outer
    response = F(1, 144) if serial else F(1, 72)
    assert inc == [
        F(17, 4) * du * response,
        F(65, 4) * du * response,
        -24 * du * response,
    ]
    for q in (0, 1):
        assert (
            F(table["shared"]["in_phase_derivatives"][q])
            - F(table["fresh"]["in_phase_derivatives"][q])
            == inc[q]
        )
    for ref, scale in (
        ("ordinary_one", (F(17, 2), F(65, 4), F(-31, 2))),
        ("ordinary_two", (F(65, 4), F(65, 2), F(-65, 2))),
    ):
        expected = (
            du
            * F(relay["fresh"]["probability"])
            * F(profile["exact_return_green"])
            / 32
        )
        if serial:
            d = F(profile["exact_return_green"])
            expected = du * ((2 * d - 1) - d * F(relay["fresh"]["probability"])) / 32
        assert [F(x) for x in doses[ref]["reference_change"]] == [
            expected * z for z in scale
        ]
    return {
        "seed": SEED,
        "starts_per_setting": N,
        "sites": [u, v],
        "source": 0,
        "ports": profile["probe_ports"],
        "partners": profile["ordinary_partners"],
        "arms": table,
        "overlap_change": list(map(str, inc)),
        "ordinary_doses": doses,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "ROUTED_CONTACTS.md",
                "routed_contacts.json",
                "shared_contact.json",
            )
        },
    }


def run():
    return {"separate_reports": flight(False), "serial_agreement": flight(True)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "routed_contact_flight.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    for name, row in z.items():
        print(
            name,
            row["sites"],
            list(map(float, map(F, row["overlap_change"]))),
            row["ordinary_doses"]["ordinary_one"],
        )
