"""Flagged recovery of a cut record by occupied neighbouring cells."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import sqrt

import numpy as np

from cosmology.record_abundance import stationary
from . import paired_relay as pair, relay, retention_budget

HERE = relay.HERE
RECEIPT = HERE / "access.json"
N, REPEATS, SEED = 8192, 64, 191206


def neighbourhoods(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    rings = []
    for a, b in edges:
        ends = {a, b}
        one = (adj[a] | adj[b]) - ends
        two = (one | set.union(*(adj[j] for j in one))) - ends
        rings.append([sorted(one), sorted(two)])
    return [[row[k] for row in rings] for k in (0, 1)]


def meter_banks(failed, kept, stocks, basis, rng):
    bank = []
    for state in stocks:
        moment = pair.meter.head.station.pair(basis, state)
        rows = []
        for i, x in enumerate(moment):
            w = state[i // 9].trace().real
            p = np.array([w / 2 + x, w / 2 - x, 1 - w])
            # In a retained branch, the recovery reference equals the target:
            # the flagged score reads zero. Fresh directions read the same heads.
            rows.append(
                [rng.multinomial(REPEATS * n, p).tolist() for n in (failed, kept)]
            )
        bank.append(rows)
    return bank


def source_read(banks, failed, kept, c):
    a = np.array(banks, dtype=np.int64)
    totals = REPEATS * np.array([failed, kept])
    if a.shape != (2, 18, 2, 3) or a.min() < 0 or not np.all(a.sum(axis=3) == totals):
        raise ValueError("Incomplete conditional source bank")
    contrast = int(totals.sum())
    signed = a[:, :, 0, 0] - a[:, :, 0, 1]
    v = float(c @ (signed[0] - signed[1]) / (2 * contrast))
    # Each failed-source outcome lies in [-1/2,1/2]. All six panels share
    # the event for the 216 means. The access counts are conditioned upon.
    err = 2 * np.sum(abs(c)) * sqrt(totals[0] * np.log(2 * 216 / 0.01) / 2) / contrast
    return [v - err, v + err]


def summary(counts, panel):
    ncut = counts[:, 0]
    failed = counts[:, 1 + panel]
    mean = failed.sum() / ncut.sum()
    # Standard error of the ratio from independent whole-board replicates.
    se = np.std(failed - mean * ncut, ddof=1) / sqrt(len(ncut)) / np.mean(ncut)
    return float(mean), float(se)


def routed_score(p, r, success):
    """Trace the discarded wire after routing a single right head into storage."""
    restored = []
    marginal = pair.meter.head.marginals(r)[1]
    for rho, right in zip(r, marginal, strict=True):
        storage = right / np.trace(right)
        wires = np.kron(rho, storage).reshape((2,) * 6)
        if success:
            wires = wires.transpose(0, 2, 1, 3, 5, 4)
        restored.append(np.trace(wires, axis1=1, axis2=4).reshape(4, 4))
    delta = pair.meter.head.station.log_blocks(
        np.asarray(restored)
    ) - pair.meter.head.station.log_blocks(r)
    return float(pair.meter.head.station.pair(delta, p - r))


def homogeneous_ceiling(histogram, cosmic_occupation, required_gain):
    """Largest local occupancy compatible with the required dark/baryon ratio."""
    from scipy.optimize import brentq

    def failure(p):
        return sum(n * (1 - p) ** int(m) for m, n in histogram.items()) / sum(
            histogram.values()
        )

    required_local = required_gain * failure(cosmic_occupation)
    if required_local > 1:
        return None
    return brentq(lambda p: (1 - p) * failure(p) - required_local, 0, 1)


def run():
    board = relay.FORMATION.federation(1)
    edges = np.asarray(board.external) // 12
    rings = neighbourhoods(board.cells, edges)
    groups = stationary.screen.observations.schedules(board)[0]
    z = json.loads(pair.RECEIPT.read_text())["source"]
    basis = pair.meter.head.joint_basis()
    c = np.array(z["basis_coefficients"])
    p0 = np.array([3 * np.eye(4) / 16, np.eye(4) / 16])
    r = p0 + pair.H * np.einsum("i,ibjk->bjk", z["direction"], basis)
    closed = routed_score(p0, r, False)
    opened = routed_score(p0, r, True)
    rng = np.random.default_rng(SEED)
    panels = {}
    for p in (0.08, 0.16, 0.32):
        cells = rng.random((N, board.cells)) < p
        initial = cells.sum(axis=1)
        for t in range(11):
            a, b = groups[t % len(groups)][:2]
            swap = rng.integers(0, 2, (N, len(a)), dtype=np.uint8).astype(bool)
            stationary.screen.paired_moves(cells, a, b, swap)
        np.testing.assert_array_equal(initial, cells.sum(axis=1))
        cuts = cells[:, edges[:, 0]] != cells[:, edges[:, 1]]
        counts = [cuts.sum(axis=1)]
        for neighbours in rings:
            fails = np.stack(
                [~cells[:, region].any(axis=1) for region in neighbours], axis=1
            )
            counts.append((cuts & fails).sum(axis=1))
        cohort = np.stack(counts, axis=1)
        settings = []
        for k, neighbours in enumerate(rings):
            n = int(cohort[:, 0].sum())
            f = int(cohort[:, k + 1].sum())
            measured, se = summary(cohort, k)
            expected = float(np.mean([(1 - p) ** len(region) for region in neighbours]))
            bank = meter_banks(f, n - f, (p0, r), basis, rng)
            score = source_read(bank, f, n - f, c)
            forecast = measured * z["expected_charge"]
            settings.append(
                dict(
                    cuts=n,
                    failed=f,
                    population_failure=expected,
                    observed_failure=measured,
                    board_standard_error=se,
                    banks=bank,
                    predicted_conditional_score=forecast,
                    observed_score=score,
                    closed_route_score=closed,
                    guaranteed_route_score=opened,
                )
            )
        unique, freq = np.unique(cohort, axis=0, return_counts=True)
        panels[str(p)] = dict(
            boards=[[*map(int, x), int(y)] for x, y in zip(unique, freq)],
            settings=settings,
        )
    gain = retention_budget.run()["galaxy_min_failed_fraction_gain"]
    cosm = panels["0.16"]["settings"]
    return dict(
        protocol_sha256=sha256((HERE / "access_protocol.md").read_bytes()).hexdigest(),
        seed=SEED,
        witness_counts=[
            dict(Counter(str(len(region)) for region in ring)) for ring in rings
        ],
        boards_per_setting=N,
        repeats=REPEATS,
        panels=panels,
        galaxy_required_gain=gain,
        maximal_failure_gains=[1 / s["population_failure"] for s in cosm],
        required_local_failure=[gain * s["population_failure"] for s in cosm],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    d = run()
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(d, f, separators=(",", ":"), allow_nan=False)
            f.write("\n")
    else:
        relay.tidal.verify(d, json.loads(RECEIPT.read_text()))
    print(d["witness_counts"], d["maximal_failure_gains"], d["required_local_failure"])
    for p, s in d["panels"].items():
        print(
            p,
            [
                (v["observed_failure"], v["board_standard_error"], v["observed_score"])
                for v in s["settings"]
            ],
        )
    print(
        "Maximum homogeneous local occupancies",
        [
            homogeneous_ceiling(h, 0.16, d["galaxy_required_gain"])
            for h in d["witness_counts"]
        ],
    )
