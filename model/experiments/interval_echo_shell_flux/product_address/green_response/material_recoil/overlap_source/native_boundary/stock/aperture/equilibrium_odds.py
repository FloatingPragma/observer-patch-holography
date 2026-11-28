"""Infer the interaction of a stored head pair from independent ordinary energy flags."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json

from mpmath import iv
import numpy as np

from . import occupied_equilibrium as heat

flight = heat.flight
lab = flight.lab
HERE = heat.HERE
bounds = lab.interval_output
SHOTS = 2**50


def head_lines(blocks):
    """Return (filled head energy, cost of a defect at the seam), with multiplicity."""
    lines = {}
    for mag, (empty, shift) in blocks.items():
        vectors = (
            [np.array([1, sign]) for sign in (1, -1)]
            if mag == 0
            else [
                np.array([1, sign, coin, sign * coin])
                for coin, sign in product((1, -1), repeat=2)
            ]
        )
        np.testing.assert_array_equal(
            np.array(vectors) @ np.array(vectors).T, np.eye(len(empty)) * len(empty)
        )
        counts = Counter()
        for w in vectors:
            energies = []
            for m in (empty + shift, -shift):
                val = int((m @ w)[0])
                np.testing.assert_array_equal(m @ w, val * w)
                assert val >= 0 and val % 4 == 0
                energies.append(val // 4)
            counts[tuple(energies)] += 1
        lines[mag] = counts
    assert all(a == 0 for _, a in lines[0])
    return lines


def counted(prob, rng):
    """Binary ordinary energy-filter trials, including unattached or rejected modes."""
    assert 0 <= prob <= 1 + 3e-14
    k = int(rng.binomial(SHOTS, min(1, prob)))
    # Simultaneous 99% Hoeffding coverage over at most 64 new binary banks.
    half = float(np.sqrt(np.log(128 / 0.01) / (2 * SHOTS))) + 1e-12
    return k, iv.mpf([max(0, k / SHOTS - half), min(1, k / SHOTS + half)])


def odds(z):
    return z[0, 0] * z[2, 2] / (z[0, 2] * z[2, 0])


def measure_board(twice, marks, radii, lines, beta, rng, native):
    costs = sorted({a for bank in lines.values() for _, a in bank})
    # Same Gershgorin lower energy for every positive contact insertion.
    lower = int(np.min(2 * twice.diagonal() - abs(twice).sum(1))) / 2
    prices, energies, flag_means, flag_bounds = {}, {}, {}, {}
    raw = {}
    for a, b in product(costs, repeat=2):
        h = twice + np.diag(2 * (a * marks[0] + b * marks[1]))
        _, energy, trace = heat.heat_shells(h, radii, include_trace=True)
        prices[a, b], energies[a, b] = trace, energy
        levels = np.linalg.eigvalsh(h / 2)
        p = np.exp(-float(beta.mid) * (levels - lower)).mean()
        k, band = counted(p, rng)
        raw[f"{a},{b}"] = k
        expected = trace * iv.exp(beta * lower) / len(twice)
        assert abs(float(expected.mid) - p) < 5e-14
        flag_bounds[a, b], flag_means[a, b] = band, p

    predictions, internal, observed, joint_counts = {}, {}, {}, {}
    for i, j in product((0, 2), repeat=2):
        z = zu = iv.mpf(0)
        mean = 0.0
        for (e, a), n in lines[i].items():
            for (f, b), m in lines[j].items():
                weight = n * m * iv.exp(-beta * (e + f))
                cell = weight * prices[a, b]
                z += cell
                zu += cell * (energies[a, b] + e + f)
                mean += n * m * np.exp(-float(beta.mid) * (e + f)) * flag_means[a, b]
        predictions[i, j], internal[i, j] = z, zu / z
        dimension = lines[i].total() * lines[j].total()
        k, band = counted(mean / dimension, rng)
        joint_counts[f"{i},{j}"] = k
        observed[i, j] = band * dimension
    ratio = odds(predictions)
    log_odds = iv.log(ratio)
    cross_energy = internal[0, 0] - internal[0, 2] - internal[2, 0] + internal[2, 2]
    out = dict(
        ordinary_reference_lower=lower,
        sharp_filter_accepts=raw,
        joint_filter_accepts=joint_counts,
        sharp_filter_bands={f"{a},{b}": bounds(v) for (a, b), v in flag_bounds.items()},
        predicted_odds=bounds(ratio),
        counted_odds=bounds(odds(observed)),
        joint_free_energy=bounds(-log_odds / beta),
        joint_internal_energy=bounds(cross_energy),
        internal_energy_per_2200_repair=bounds(cross_energy / 14),
        joint_price_required_beyond_spectrum=bounds((log_odds - native) / beta),
    )
    return out, prices, flag_bounds


def replica_parameters(bank):
    """One-tracer bounds before summing the shared hidden head lines."""
    active = sorted({a for a, _ in bank if a > 0})
    u = [x / bank[0, 0] for a in active for x in (bank[a, 0], bank[0, a])]
    ratio = [
        bank[a, b] * bank[0, 0] / (bank[a, 0] * bank[0, b])
        for a, b in product(active, repeat=2)
    ]
    return max(x.b for x in u), min([iv.mpf(1), *(x.a for x in ratio)])


def every_inventory(banks, lines, beta):
    """Bound all nonnegative integer counts of distinguishable independent replicas."""
    weights = [n * iv.exp(-beta * e) for (e, _), n in lines[2].items()]
    inactive = sum(
        n * iv.exp(-beta * e) for (e, a), n in lines[2].items() if a == 0
    ) / sum(weights)
    params = [replica_parameters(v) for v in banks]
    assert all(0 < r <= 1 and 0 < t < 1 for t, r in params)
    # r^k >= t for each material; inactive arm pairs have ratio exactly one.
    k = 1
    while k < 2**20 and all(r ** (2 * k) >= t for t, r in params):
        k *= 2
    assert all(r**k >= t for t, r in params)
    # The global loss T is the product of the single-arm attenuations.
    # Each split gives a rigorous bound on both sides T >= cut and T <= cut.
    best = (iv.mpf(0), 0)
    for j in range(1, 33):
        cut = iv.mpf(2) ** (-j)
        frac = (1 - inactive) * cut / (inactive + (1 - inactive) * cut)
        bound = min((cut ** (iv.mpf(1) / k)).a, (1 - frac**2).a)
        if bound > best[0]:
            best = (bound, j)
    return dict(
        material_upper_attenuation_lower_odds=[
            (bounds(t)[1], bounds(r)[0]) for t, r in params
        ],
        inactive_head_fraction=bounds(inactive),
        power=k,
        split_binary_exponent=best[1],
        odds_floor=bounds(best[0])[0],
    )


def run():
    iv.dps = 65
    beta = iv.log(iv.mpf(8) / 7)
    board, old, edges, radii, _ = flight.chart()
    log_native = iv.mpf(lab.contact.overlap_bands(old["source"])["reference_log_odds"])
    rng = np.random.default_rng(3769201)
    mats = lab.contact.blocks()
    out = {}
    for tag, blocks in mats.items():
        lines = head_lines(blocks)
        panels, precise, measured = {}, [], []
        for load, species in ((0, "vacancy"), (2, "excess")):
            twice, marks = flight.apparatus(board, edges, load, -1, 1)
            if tag == "exchange":
                # On blank seams exchange has unit off-diagonal (twice native
                # energy); on a sharp occupied reference only unequal loads pay.
                np.fill_diagonal(twice, 0)
                np.fill_diagonal(twice, -twice.sum(axis=1))
                _, _, _, raw = lab.contact.stock.kernels()
                idx = [3 * load + 1, 3 + load]
                np.testing.assert_array_equal(raw[np.ix_(idx, idx)], [[2, -2], [-2, 2]])
                assert raw[4, 4] == 0
            row, sharp, flags = measure_board(
                twice, marks, radii, lines, beta, rng, log_native
            )
            panels[species] = row
            precise.append(sharp)
            measured.append(flags)
        arbitrary = every_inventory(measured, lines, beta)
        missing_price = (iv.log(iv.mpf(arbitrary["odds_floor"])) - log_native) / beta
        out[tag] = dict(
            required_joint_price_floor_all_inventories=bounds(missing_price)[0],
            head_lines={
                str(k): [[e, a, n] for (e, a), n in v.items()] for k, v in lines.items()
            },
            materials=panels,
            any_replica_inventory=every_inventory(precise, lines, beta),
            any_inventory_from_filters=arbitrary,
        )
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("EQUILIBRIUM_ODDS.md", "pair_contact.json")
        },
        beta=bounds(beta),
        native_log_odds=bounds(log_native),
        native_odds=bounds(iv.exp(log_native)),
        trials_per_filter=SHOTS,
        results=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "equilibrium_odds.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    for law, case in result["results"].items():
        print(
            law,
            {s: r["predicted_odds"] for s, r in case["materials"].items()},
            case["any_inventory_from_filters"],
        )
