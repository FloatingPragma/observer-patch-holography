"""Compare overlap and clamped ordinary sources at the same thermal instruments."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json

from mpmath import iv
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from . import thermal_pair as exp

HERE = exp.HERE
kernel = exp.kernel
probe = exp.probe
motion = exp.motion


def check_clamp():
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6)]
    size = 7
    fixed = 2
    _, _, raw, _ = exp.contact.stock.kernels()
    live = np.zeros((size, size), int)
    for i, j in edges:
        live[i, j] = live[j, i] = 1

    def compressed(number, occupied):
        allowed = [
            tuple(x)
            for x in combinations([k for k in range(size) if k != fixed], number)
        ]
        ix = {x: i for i, x in enumerate(allowed)}
        out = np.zeros((len(ix), len(ix)), int)
        for k, free in enumerate(allowed):
            word = np.zeros(size, int)
            word[list(free)] = 1
            word[fixed] = occupied
            for a, b in edges:
                col = 3 * word[a] + word[b]
                for r in np.flatnonzero(raw[:, col]):
                    assert r < 9
                    after = word.copy()
                    after[a], after[b] = divmod(int(r), 3)
                    if after[fixed] == occupied:
                        after[fixed] = 0
                        out[ix[tuple(np.flatnonzero(after))], k] += raw[r, col]
        return out, allowed

    bare = int((compressed(0, 1)[0] - compressed(0, 0)[0]).item())
    occ, words = compressed(2, 1)
    diff = occ - compressed(2, 0)[0] - bare * np.eye(len(occ), dtype=int)
    expected = (
        -4 * (2 * motion.K + 1) * np.array([sum(live[fixed, list(s)]) for s in words])
    )
    np.testing.assert_array_equal(diff, np.diag(expected))
    return dict(
        ordinary_pairs=len(words),
        checked_entries=diff.size,
        bare_four_energy=bare,
        permitted_shifts=sorted(set((expected // 4).tolist())),
    )


def quotient(a, b):
    if b[0] <= 0 <= b[1]:
        return None
    return exp.interval_output(iv.mpf(a) / iv.mpf(b))


def empirical_work(rng, prob, numerical):
    e = np.sqrt(exp.contact.RISK / (2 * exp.SHOTS))
    raw = {k: exp.sample(rng, p).tolist() for k, p in prob.items()}
    p = {k: np.asarray(v) / exp.SHOTS for k, v in raw.items()}
    work = {}
    for i, j in ((0, 1), (2, 6)):

        def band(tag, idx):
            return iv.mpf(
                [p[tag][idx] - e - numerical[tag], p[tag][idx] + e + numerical[tag]]
            )

        work[f"{i}_{j}"] = exp.interval_output(
            -iv.log(
                band("occupied", j)
                * band("blank", i)
                / (band("blank", j) * band("occupied", i))
            )
            / kernel.BETA
        )
    return raw, work


def run():
    check = check_clamp()
    board = motion.cold.federation(0)
    edges = exp.contact.stock.carrier()[1][:2]
    h0, sites, i, j, _ = motion.hamiltonian(board, edges)
    rad = motion.cold.distances(board, sorted({x for e in edges for x in e}))
    fixed = int(np.flatnonzero(rad == 1)[0])
    adj = board.adjacency().toarray()
    assert fixed not in np.ravel(edges)
    keep = (i != fixed) & (j != fixed)
    pair_i, pair_j = i[keep], j[keep]
    h = h0[keep][:, keep]
    rev = np.full(len(i), -1, int)
    rev[keep] = np.arange(keep.sum())
    sites = rev[sites]
    assert min(sites) >= 0
    near = adj[fixed, pair_i].astype(int) + adj[fixed, pair_j].astype(int)
    placed = h - sparse.diags((2 * motion.K + 1) * near, format="csr", dtype=float)
    w, err = kernel.baseline(h)
    occupied, err1 = kernel.baseline(placed)
    changes, cert = kernel.single_updates(h, sites)
    distance = shortest_path(h != 0, directed=False, indices=sites[0], unweighted=True)[
        sites[1]
    ]
    heat = kernel.mixed_bounds(h, sites, distance)

    def bins(values):
        return exp.shells(values, pair_i, pair_j, rad)

    base, one, counts = bins(w), bins(occupied), bins(np.ones_like(w))
    sums = [{u: bins(v) for u, v in row.items()} for row in changes]
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    gain = motion.timing.lab.ordinary_gain(tape) / 4
    delta = probe.log_factor(old["source"])
    work = exp.work_bands(base, counts, sums, err, cert, heat, delta)
    sectors = {}
    for a, b in product(exp.STRENGTHS, repeat=2):
        v = w + changes[0].get(a, 0) + changes[1].get(b, 0)
        sectors[f"{-a},{-b}"] = v / v.sum()
    cond = motion.binary_responses(sectors)
    table = np.array(old["accepted_weight"])
    native = table.ravel() / table.sum()
    indep = np.outer(table.sum(1), table.sum(0)).ravel() / float(table.sum()) ** 2
    distr = dict(
        zip(exp.LABELS, cond, strict=True),
        native=native @ cond,
        product=indep @ cond,
        blank=w / w.sum(),
        occupied=occupied / occupied.sum(),
    )
    exp_error = len(w) * (
        err + heat["upper_each_diagonal"] + 2 * cert["tail"]
    ) + 2 * cert["rounding_relative"] * sum(
        max(v.sum() for v in row.values()) for row in changes
    )
    local_error = exp_error / (w.sum() - exp_error)
    ordinary_error = len(w) * err1 / (occupied.sum() - len(w) * err1)
    numerical = {"blank": local_error, "occupied": ordinary_error}
    rng = np.random.default_rng(735793)
    raw, held = empirical_work(rng, {"blank": base, "occupied": one}, numerical)

    def weight(values, err, h):
        return iv.mpf([values[h] - err * counts[h], values[h] + err * counts[h]])

    tests = {}
    for a, b in ((0, 1), (2, 6)):
        f = exp.interval_output(
            -iv.log(
                weight(one, err1, b)
                * weight(base, err, a)
                / (weight(base, err, b) * weight(one, err1, a))
            )
            / kernel.BETA
        )
        overlap = work["relative_to_first"][1] if a == 0 else work["enrolled"]["2_6"]
        tests[f"work_{a}_{b}"] = dict(
            ordinary=f,
            ordinary_held=held[f"{a}_{b}"],
            overlap=overlap,
            ratio=quotient(overlap, f),
        )
    banks = {}
    for r in motion.RADII:
        k, v = probe.operators(board, rad <= r)
        terms = {}
        banks[str(r)] = {}
        for tag, op in (("kinetic", k), ("contrast", v)):
            meter = exp.meter_bank(op.diagonal()[keep], distr, rng)
            banks[str(r)][tag] = meter
            b = exp.meter_intervals(meter, local_error)
            b["occupied"] = exp.meter_intervals(meter, ordinary_error)["occupied"]
            mixed = b["00"] + b["11"] - b["10"][::-1] - b["01"][::-1]
            terms[tag] = [
                b["occupied"] - b["blank"][::-1],
                probe.mix_bands(delta, mixed),
                b["product"] - b["native"][::-1],
            ]
        normal, overlap, held_release = [
            (terms["kinetic"][t] + probe.mix_bands(gain, terms["contrast"][t])).tolist()
            for t in range(3)
        ]
        tests[f"release_{r}"] = dict(
            ordinary=normal,
            overlap=overlap,
            overlap_held=held_release,
            ratio=quotient(overlap, normal),
        )
    ratios = [row["ratio"] for row in tests.values() if row["ratio"] is not None]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "THERMAL_ORDINARY.md",
                HERE / "thermal_pair.json",
                HERE / "mediated_release.json",
            )
        },
        compression=check,
        clamp=fixed,
        pair_addresses=len(w),
        distance=int(distance),
        ordinary_baseline_radius=err1,
        numerical_probability_error={k: float(v) for k, v in numerical.items()},
        source_work=work,
        work_counts=raw,
        release_counts=banks,
        tests=tests,
        common_ratio=bool(max(r[0] for r in ratios) <= min(r[1] for r in ratios)),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "thermal_ordinary.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    print(result["clamp"], result["tests"], result["common_ratio"])
