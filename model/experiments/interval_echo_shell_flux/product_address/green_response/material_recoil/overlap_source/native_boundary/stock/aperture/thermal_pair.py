"""Measured connected charge, thermal position work and a separately released double."""

import argparse
from hashlib import sha256
from itertools import product
import json
import math

from mpmath import iv
import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import pair_probe as probe
from . import thermal_kernel as kernel

contact = probe.contact
motion = probe.motion
HERE = contact.HERE
BATCH = 2**46
BATCHES = 1024
SHOTS = BATCH * BATCHES
LABELS = ("00", "01", "10", "11")
STRENGTHS = (0, 32, 33)


def shells(values, i, j, rad):
    return (
        np.bincount(rad[i], values, minlength=max(rad) + 1)
        + np.bincount(rad[j], values, minlength=max(rad) + 1)
    ) / 2


def parity_weights(a, bit):
    if not bit:
        return int(a == 0)
    return 0.5 if a == 0 else 0.25


def arm_weights(cell, parities):
    return [
        parity_weights(a, cell // 2) * parity_weights(b, cell % 2) for a, b in parities
    ]


def interval_output(bound):
    return [
        float(np.nextafter(float(bound.a), -np.inf)),
        float(np.nextafter(float(bound.b), np.inf)),
    ]


def work_bands(base, counts, changes, base_error, cert, heat, source):
    """Cancel the shell-independent Gibbs normalization before bounding radial work."""
    iv.dps = 40
    interval = iv.mpf

    def positive(v, e):
        return interval([max(0.0, v - e), v + e])

    def increment(v, n):
        # Roundoff bound on nonnegative sums, plus the positive Taylor tail per address.
        return positive(v, 2 * v * cert["rounding_relative"] + n * cert["tail"])

    x = [
        {
            u: [increment(v, n) for v, n in zip(values, counts, strict=True)]
            for u, values in row.items()
        }
        for row in changes
    ]
    base = [
        positive(v, n * base_error + 1e-10 * v)
        for v, n in zip(base, counts, strict=True)
    ]
    z = sum(base)
    sites = sum(counts)
    both = heat["upper_each_diagonal"]
    keys = list(product(STRENGTHS, repeat=2))
    table = [
        [interval(int(k)) / 448 for k in (48, 176, 112, 112)],
        list(map(interval, [1, -1, -1, 1])),
    ]
    raw = np.array([arm_weights(cell, keys) for cell in range(4)]).T
    shifts = [{u: sum(bins) for u, bins in row.items()} for row in x]

    def weights(weights):
        values = []
        for (a, b), coeffs in zip(keys, raw, strict=True):
            norm = (
                z
                + shifts[0].get(a, 0)
                + shifts[1].get(b, 0)
                + (interval([0, sites * both]) if a * b else 0)
            )
            values.append(
                sum(float(v) * w for v, w in zip(coeffs, weights, strict=True)) / norm
            )
        return values

    p, m = map(weights, table)
    ap, am = sum(p), sum(m)
    delta = interval(source.tolist())
    aq = ap + delta * am
    geom = [interval(0) for _ in base]
    t = [interval(0) for _ in base]
    for side, row in enumerate(x):
        for u, bins in row.items():
            cp = sum(v for v, k in zip(p, keys, strict=True) if k[side] == u)
            cm = sum(v for v, k in zip(m, keys, strict=True) if k[side] == u)
            change = delta * (cm * ap - cp * am) / (ap * aq)
            for h, (v, b) in enumerate(zip(bins, base, strict=True)):
                t[h] += v / b * cp / ap
                geom[h] += v / b * change
    g = []
    for h, b in enumerate(base):
        maxres = counts[h] * both / b
        # The conditional residual is a convex mean for each of the two source arms.
        grow = interval([-float(maxres.b) * 1.00001, float(maxres.b) * 1.00001])
        g.append(
            iv.log1p(
                (geom[h] + grow) / (1 + t[h] + interval([0, float(maxres.b) * 1.00001]))
            )
        )

    def work(a, b):
        return interval_output(-(g[b] - g[a]) / kernel.BETA)

    return dict(
        relative_to_first=[work(0, h) for h in range(len(base))],
        successive=[work(h, h + 1) for h in range(len(base) - 1)],
        enrolled={f"{a}_{b}": work(a, b) for a, b in ((2, 4), (4, 6), (2, 6))},
    )


def sample(rng, probabilities):
    p = np.maximum(np.asarray(probabilities), 0)
    return rng.multinomial(BATCH, p / p.sum(), size=BATCHES).sum(0)


def meter_bank(values, distributions, rng):
    # One terminal sharp pair address determines this calibrated observable.
    levels, index = np.unique(values, return_inverse=True)
    banks = {
        tag: sample(rng, np.bincount(index, weights=p, minlength=len(levels))).tolist()
        for tag, p in distributions.items()
    }
    return dict(values=levels.tolist(), counts=banks)


def meter_intervals(meter, numerical):
    vals = np.array(meter["values"])
    margin = np.ptp(vals) * (np.sqrt(contact.RISK / (2 * SHOTS)) + numerical)
    out = {}
    for tag, bank in meter["counts"].items():
        assert sum(bank) == SHOTS
        out[tag] = vals @ (np.asarray(bank) / SHOTS) + np.array([-margin, margin])
    return out


def decode(raw, delta, gain):
    releases = {}
    for cut, meters in raw["release"].items():
        terms = {}
        for tag, meter in meters.items():
            b = meter_intervals(meter, raw["probability_error"])
            mixed = b["00"] + b["11"] - b["01"][::-1] - b["10"][::-1]
            terms[tag] = [
                probe.mix_bands(delta, mixed),
                b["product"] - b["native"][::-1],
            ]
        releases[cut] = {
            key: (
                terms["kinetic"][i] + probe.mix_bands(gain, terms["contrast"][i])
            ).tolist()
            for i, key in enumerate(("from_source", "held"))
        }
    shell = raw["shells"]
    e = np.sqrt(contact.RISK / (2 * SHOTS)) + raw["probability_error"]
    freq = {a: np.asarray(row) / SHOTS for a, row in shell.items()}
    assert all(sum(row) == SHOTS for row in shell.values())
    forward = probe.mix_bands
    mixed = freq["00"] - freq["10"] - freq["01"] + freq["11"]
    signed = [forward(delta, v + np.array([-4 * e, 4 * e])).tolist() for v in mixed]
    held = np.column_stack(
        (
            freq["product"] - freq["native"] - 2 * e,
            freq["product"] - freq["native"] + 2 * e,
        )
    )

    def ratio(a, b):
        return np.log([(a - e) / (b + e), (a + e) / (b - e)])

    odds = [ratio(a, b) for a, b in zip(freq["product"], freq["native"], strict=True)]
    direct = {
        f"{a}_{b}": ((odds[a] - odds[b][::-1]) / kernel.BETA).tolist()
        for a, b in ((0, 1), (2, 4), (4, 6), (2, 6))
    }
    filter_shift = (raw["one_filter"]["product"] - raw["one_filter"]["native"]) / SHOTS
    single = (
        np.asarray(raw["one_unit"]["product"]) / SHOTS
        - np.asarray(raw["one_unit"]["native"]) / SHOTS
    )
    return dict(
        release=releases,
        shells_from_source=signed,
        shells_held=held.tolist(),
        work_held=direct,
        one_filter_change=[filter_shift - 2 * e, filter_shift + 2 * e],
        one_unit_changes=np.column_stack((single - 2 * e, single + 2 * e)).tolist(),
    )


def single_control(board, edges, rad, rng):
    adjacency = board.adjacency().toarray().copy()
    for a, b in edges:
        adjacency[a, b] = adjacency[b, a] = 0
    marked = np.bincount(np.asarray(edges).ravel(), minlength=board.n)
    h = (np.diag(adjacency.sum(1)) - adjacency) / 2 - np.diag(motion.K * marked)
    energies, basis = np.linalg.eigh(h)
    p = (basis * basis) @ np.exp(-kernel.BETA * energies)
    bins = np.bincount(rad, p)
    return {tag: sample(rng, bins).tolist() for tag in ("native", "product")}


def run():
    contact.blocks()  # Native calibration and zero response for every single-unit input.
    _, table, _, _ = contact.source_counts()
    edges = contact.stock.carrier()[1][:2]
    board = motion.cold.federation(0)
    h, sites, i, j, _ = motion.hamiltonian(board, edges)
    rad = motion.cold.distances(board, sorted({x for e in edges for x in e}))
    distance = shortest_path(h != 0, directed=False, indices=sites[0], unweighted=True)[
        sites[1]
    ]
    heat = kernel.mixed_bounds(h, sites, distance)
    w, error = kernel.baseline(h)
    updates, cert = kernel.single_updates(h, sites)
    assert min(w) > error and min(x.min() for row in updates for x in row.values()) >= 0
    base = shells(w, i, j, rad)
    sizes = shells(np.ones_like(w), i, j, rad)
    shells_by_mark = [
        {u: shells(a, i, j, rad) for u, a in row.items()} for row in updates
    ]
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    delta = probe.log_factor(old["source"])
    gain = motion.timing.lab.ordinary_gain(tape) / 4
    work = work_bands(base, sizes, shells_by_mark, error, cert, heat, delta)
    sectors = {}
    acceptance = {}
    floor = kernel.bounds(h)[0] - 33
    for a, b in product(STRENGTHS, repeat=2):
        pos = w + updates[0].get(a, 0) + updates[1].get(b, 0)
        sectors[f"{-a},{-b}"] = pos / pos.sum()
        acceptance[f"{-a},{-b}"] = float(pos.mean() * np.exp(kernel.BETA * floor))
    cond = motion.binary_responses(sectors)
    native = table.ravel() / table.sum()
    comparison = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    distributions = dict(
        zip(LABELS, cond, strict=True), native=native @ cond, product=comparison @ cond
    )
    rng = np.random.default_rng(917611)
    prob = {tag: shells(v, i, j, rad) for tag, v in distributions.items()}
    raw = dict(
        shells={tag: sample(rng, v).tolist() for tag, v in prob.items()}, release={}
    )
    abs_error = len(w) * (
        error + 2 * cert["tail"] + heat["upper_each_diagonal"]
    ) + 2 * cert["rounding_relative"] * sum(
        max(v.sum() for v in row.values()) for row in updates
    )
    raw["probability_error"] = float(abs_error / (sum(w) - abs_error))
    f = motion.binary_responses(acceptance)
    raw["one_filter"] = {
        tag: int(rng.binomial(BATCH, float(dist @ f), size=BATCHES).sum())
        for tag, dist in (("native", native), ("product", comparison))
    }
    raw["one_unit"] = single_control(board, edges, rad, rng)
    exact = {}
    for r in motion.RADII:
        k, v = probe.operators(board, rad <= r)
        raw["release"][str(r)] = {
            tag: meter_bank(op.diagonal(), distributions, rng)
            for tag, op in (("kinetic", k), ("contrast", v))
        }
        exact[str(r)] = float(
            (distributions["product"] - distributions["native"])
            @ (k.diagonal() + motion.K * v.diagonal())
        )
    dec = decode(raw, delta, gain)
    assert dec["one_filter_change"][0] <= 0 <= dec["one_filter_change"][1]
    assert all(a <= 0 <= b for a, b in dec["one_unit_changes"])
    for r, target in exact.items():
        assert all(lo <= target <= hi for lo, hi in dec["release"][r].values())
    for k in ("shells_from_source", "shells_held"):
        for (lo, hi), target in zip(
            dec[k], prob["product"] - prob["native"], strict=True
        ):
            assert lo <= target <= hi
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "THERMAL_PAIR.md",
                HERE / "pair_contact.json",
                HERE / "mediated_release.json",
            )
        },
        edges=[list(e) for e in edges],
        pair_addresses=len(w),
        heat=heat,
        polynomial_degree=kernel.DEGREE,
        baseline_per_address_error=error,
        increment_error=cert,
        sites_per_shell=np.bincount(rad).tolist(),
        ordinary_base_shells=base.tolist(),
        single_increments=[
            {str(u): v.tolist() for u, v in row.items()} for row in shells_by_mark
        ],
        thermometer=kernel.BETA,
        source_delta=delta.tolist(),
        certified_work=work,
        shell_probabilities={tag: v.tolist() for tag, v in prob.items()},
        model_release=exact,
        one_unit_source_change=0,
        energy_filter_floor=floor,
        single_filter_gain_upper=float(
            delta[1] * math.exp(kernel.BETA * floor) * heat["upper_each_diagonal"] / 4
        ),
        shots=SHOTS,
        banks=raw,
        decoded=dec,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    group = cli.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "thermal_pair.json"
    if args.counts_only:
        result = json.loads(path.read_bytes())
        old = json.loads((HERE / "pair_contact.json").read_bytes())
        tape = json.loads((HERE / "mediated_release.json").read_bytes())
        contact.timing.lab.lab.cold.apparatus.check(
            result["decoded"],
            decode(
                result["banks"],
                probe.log_factor(old["source"]),
                motion.timing.lab.ordinary_gain(tape) / 4,
            ),
        )
    else:
        result = run()
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            contact.timing.lab.lab.cold.apparatus.check(
                json.loads(path.read_bytes()), result
            )
    print(result["certified_work"]["enrolled"], result["decoded"]["release"])
