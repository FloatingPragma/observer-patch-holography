"""Count first exits across a fixed surface using its measured mean link periods."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import shortest_path

from . import sector_register_balls as balls
from . import sector_register_pairs as pairs

HERE = balls.HERE
DEADLINES = (1, 2, 4, 8, 16, 32, 64, 128, 256)
GRID = 36**8 * 8192


def first_exits(W, inside):
    restricted = W[inside][:, inside].tocsr()
    # Each row of the killed integer kernel sums to at most 36. Subtract
    # survival from one. Eight steps have exact denominator dividing grid.
    grid = GRID
    assert max(len(inside), 36) * grid < 2**63
    v = np.full(restricted.shape[0], grid, dtype=np.int64)
    out = []
    for n in range(1, max(DEADLINES) + 1):
        numer = restricted @ v
        if n <= 8:
            assert np.all(numer % 36 == 0)
        v = (numer + 18) // 36
        if n in DEADLINES:
            out.append((len(v) * grid - int(v.sum())) / (len(v) * grid))
    return out


def interval_product(a, b):
    p = np.multiply.outer(a, b)
    return np.array([p.min(), p.max()])


def contrast(arr, delta, err, n):
    weights = delta.ravel()
    mid = float((arr - arr.mean()) @ weights)
    tol = (
        np.ptp(arr) * err
        + (abs(weights).sum() + err) * max(0, n - 8) / (2 * GRID)
        + 1e-15
    )
    if np.ptp(arr) == 0 and n <= 8:
        return [0.0, 0.0]
    return [mid - tol, mid + tol]


def run():
    (ports, _, adj, uv, cfg, _, periods), _, _, delta, err = pairs.apparatus()
    previous = json.loads((HERE / "sector_register_balls.json").read_text())
    cert = previous["all_register_size_certificate"]
    center, volume = cert["center"], cert["volume"]
    witness = next(
        r
        for r in previous["examples"]
        if r["center"] == center and r["volume"] == volume
    )
    distance = shortest_path(
        csr_matrix(adj), unweighted=True, indices=np.searchsorted(ports, center)
    )
    inside = np.zeros(len(ports), bool)
    inside[np.argsort(distance, kind="stable")[:volume]] = True
    u, v = uv.T
    assert (
        ports[uv[inside[u] != inside[v]]].tolist()
        == cert["constant_one_period_boundary_links"]
    )
    clocks = np.array([periods[tuple(c)] for c in cfg])
    times, group = np.unique(
        (clocks[:, None] + clocks).reshape(-1, len(uv)), axis=0, return_inverse=True
    )
    tables = {"interior_exit": [], "exterior_entry": []}
    exterior_kernel = None
    for t in times:
        assert np.all(12 % t == 0)
        W = csr_matrix(
            (12 // np.tile(t, 2), (np.r_[u, v], np.r_[v, u])), shape=adj.shape
        )
        wait = 36 - np.asarray(W.sum(1)).ravel()
        assert min(wait) >= 0
        W = (W + diags(wait, dtype=np.int64)).tocsr()
        ext = W[~inside][:, ~inside]
        if exterior_kernel is None:
            exterior_kernel = ext
        else:
            assert (ext != exterior_kernel).nnz == 0
        for tag, domain in (("interior_exit", inside), ("exterior_entry", ~inside)):
            tables[tag].append(first_exits(W, domain))
    # The word preparation controls mixture frequencies. Identify its joint
    # configurations from the observed home-occupancy classes.
    words = pairs.chart.source.thermal.apparatus()[0]
    lookup = {tuple(c): i for i, c in enumerate(cfg)}
    native_ids = [lookup[tuple(w[:4] > 0)] for w in words]
    matched = np.zeros((len(cfg), len(cfg)))
    for sign, arm in zip(
        [1, -1],
        previous["matched_local_clock_control"]["word_pairs_by_arm"],
        strict=True,
    ):
        for i, j in arm:
            matched[native_ids[i], native_ids[j]] += sign / 2
    ordinary = witness["ordinary"]["pooled_ordinary"][1]
    source = witness["cross"]["pooled_ordinary"][1]
    rows = {}
    for k, table in tables.items():
        rows[k] = []
        all_times = np.array(table)[group]
        for n, arr in zip(DEADLINES, all_times.T, strict=True):
            o, c = (
                np.array(contrast(arr, delta[tag], err[tag], n))
                for tag in ("ordinary", "cross")
            )
            matched_band = contrast(arr, matched, 0.0, n)
            if k == "exterior_entry":
                # Every killed transition and holding probability is identical.
                o, c, matched_band = np.zeros(2), np.zeros(2), [0.0, 0.0]
            row = {
                "matched_together_minus_crossed_change": matched_band,
                "deadline": n,
                "ordinary_change": o.tolist(),
                "source_change": c.tolist(),
            }
            if not o[0] <= 0 <= o[1]:
                gain = interval_product(-np.array(ordinary), 1 / o[::-1])
                residual = source + interval_product(gain, c)
                row.update(
                    ordinary_entropy_coefficient=gain.tolist(),
                    source_entropy_residual_nats=residual.tolist(),
                )
            rows[k].append(row)
    # Both comparisons can start from half together and half independent
    # native members. Switch a fraction in only one half of that baseline.
    # Its complete one-member laws then coincide in all three preparations.
    qc = -np.array(
        previous["matched_local_clock_control"]["relative_modular_contrast_nats"]
    )[::-1]
    ret = np.array(witness["cross"]["pooled_retained"][1])
    ex = np.array(witness["cross"]["pooled_exchanged"][1])
    qn = ret - ex[::-1]
    target = 0.00005
    last = rows["interior_exit"][-1]
    contrasts = {}
    for label, q, flux in (
        (
            "controlled",
            qc,
            -np.array(last["matched_together_minus_crossed_change"])[::-1],
        ),
        ("distinct_native", qn, np.array(last["source_change"])),
    ):
        fraction = 2 * target / q[::-1]
        assert 0 < min(fraction) <= max(fraction) < 1
        contrasts[label] = {
            "source_charge_nats": q.tolist(),
            "fraction_of_its_own_half": fraction.tolist(),
            "exit_change_at_equal_charge": interval_product(
                fraction / 2, flux
            ).tolist(),
        }
    assert (
        contrasts["controlled"]["exit_change_at_equal_charge"][0]
        > 0
        > contrasts["distinct_native"]["exit_change_at_equal_charge"][1]
    )
    return {
        "common_baseline_equal_mean_charge": {
            "target_increase_nats": target,
            "deadline": DEADLINES[-1],
            "interventions": contrasts,
        },
        "input_sha256": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("SECTOR_REGISTER_ESCAPE.md", "sector_register_balls.json")
        },
        "center": center,
        "volume": volume,
        "exterior_first_entry_identical_killed_kernel": exterior_kernel is not None,
        "count_intervals": rows,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dst = HERE / "sector_register_escape.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), result
        )
    print(result["count_intervals"]["interior_exit"])
