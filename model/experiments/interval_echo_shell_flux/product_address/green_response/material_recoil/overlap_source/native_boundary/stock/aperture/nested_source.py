"""Test regional modular readings and measure a boost-weighted ordinary spectral line."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import closed_balls as balls

HERE = balls.HERE
PRIME = 1009
SHOTS = 2**24


def minor(rows):
    """Select an integer square minor certified invertible over a prime field."""
    a = np.asarray(rows, dtype=np.int64) % PRIME
    ids = np.arange(len(a))
    det = 1
    for j in range(a.shape[1]):
        candidates = np.flatnonzero(a[j:, j])
        if not len(candidates):
            raise ValueError("These comparisons do not certify full column rank")
        k = int(candidates[0]) + j
        a[[j, k]] = a[[k, j]]
        ids[[j, k]] = ids[[k, j]]
        p = int(a[j, j])
        det = det * p % PRIME
        factors = a[j + 1 :, j] * pow(p, -1, PRIME) % PRIME
        a[j + 1 :, j:] = (a[j + 1 :, j:] - factors[:, None] * a[j, j:]) % PRIME
    # Determinant of the selected original rows in their returned pivot order.
    return ids[: a.shape[1]], det


def weights(hops, r):
    return np.maximum(0, r * r - hops * hops)


def comparisons(hops, distances, start):
    rows, ids = [], []
    for r in (start, start + 1):
        same = np.all((distances < r) == (distances < r + 1), axis=1)
        proper = hops.max(axis=1) >= r + 1
        for center in np.flatnonzero(same & proper):
            rows.append(
                r * weights(hops[center], r + 1) - (r + 1) * weights(hops[center], r)
            )
            ids.append([int(center), r])
    rows = np.array(rows)
    chosen, det = minor(rows)
    confirm, second = minor(rows[chosen])
    assert np.array_equal(confirm, np.arange(len(chosen))) and det == second
    return dict(
        rows=len(rows),
        rank=len(chosen),
        modulus=PRIME,
        determinant_mod_prime=det,
        minor_rows=[ids[i] for i in chosen],
    )


def spectral_bank(rng, counts, distances, radius):
    # Each new preparation gets one fair endpoint and two uniform deadlines.
    # Resolve the native single-bond lines after those assignments. Record
    # every energy even on a timeout so the apparatus can be audited.
    lines = [0, 8, 9, 32, 33]
    counts_by_gate = np.zeros((2, len(lines)), dtype=np.int64)
    cut = np.maximum.outer(np.arange(radius), np.arange(radius))
    for n, count in enumerate(counts):
        settings = rng.multinomial(
            int(count), np.full(2 * radius**2, 1 / (2 * radius**2))
        )
        settings = settings.reshape(2, radius, radius)
        lo = [0, 1, 3][n]
        first = settings if n == 0 else rng.binomial(settings, 0.5)
        for endpoint, d in enumerate(distances):
            for gate in (0, 1):
                mask = (cut >= d) == bool(gate)
                counts_by_gate[gate, lo] += first[endpoint][mask].sum()
                if n:
                    counts_by_gate[gate, lo + 1] += (settings - first)[endpoint][
                        mask
                    ].sum()
    assert counts_by_gate.sum() == SHOTS
    return counts_by_gate, (radius / 2) * float(counts_by_gate[1] @ lines) / SHOTS


def spectra(w, cells, hops, h):
    # Independent check of the five lines and parity weights on the source bond.
    hb = h[:9, :9] / 4
    for n, edges in enumerate(([0], [8, 9], [32, 33])):
        idx = sorted({n, 3 * n})
        v, u = np.linalg.eigh(hb[np.ix_(idx, idx)])
        np.testing.assert_allclose(v, edges)
        np.testing.assert_allclose(u[0] ** 2, 1 / len(idx))
        assert not np.any(hb[np.ix_(idx, sorted(set(range(9)) - set(idx)))])
    old = json.loads((HERE / "closed_energy.json").read_bytes())["forecast_from_source"]
    rng = np.random.default_rng(1951229)
    trials = (np.ones(len(w)), np.asarray(balls.energy.least_change(cells), float))
    results = []
    for center in (20, 12):
        for r in (2, 4, 8):
            distances = hops[center, [13, 20]]
            factor = Q(int(weights(distances, r).sum()), 4 * r)
            reads, means = [], []
            for freq in trials:
                cnt = rng.multinomial(SHOTS, freq / freq.sum())
                by_n = np.bincount(w[:, 1], weights=cnt, minlength=3).astype(int)
                bank, mean = spectral_bank(rng, by_n, distances, r)
                reads.append(bank.tolist())
                means.append(mean)
            measured = means[1] - means[0]
            # For a difference of two independent sample means in [0,33R/2],
            # Hoeffding with a union bound over all six reported comparisons.
            error = 33 * r / 2 * np.sqrt(np.log(2 * 6 / 0.01) / SHOTS)
            obs = [measured - error, measured + error]
            model = float(factor * Q(31, 56))
            predicted = [float(factor) * p for p in old["bond_energy"]]
            assert max(obs[0], predicted[0]) <= min(obs[1], predicted[1])
            # Given the middle reading, a registered nonzero load resolves two
            # equiprobable lines. Its conditional surprise is exactly log(2).
            # A closed gate stores a common timeout. The independent deadline
            # labels themselves have the same distribution in both arms.
            parity = [sum(arm[1][1:]) / SHOTS for arm in reads]
            parity_change = (parity[1] - parity[0]) * np.log(2)
            parity_error = np.log(2) * np.sqrt(np.log(2 * 6 / 0.01) / SHOTS)
            alpha = 2 * factor / r
            parity_target = -float(alpha) * np.log(2) / 28
            parity_band = [parity_change - parity_error, parity_change + parity_error]
            assert parity_band[0] <= parity_target <= parity_band[1] < 0
            results.append(
                dict(
                    center=center,
                    radius=r,
                    endpoint_hops=distances.tolist(),
                    lines_by_gate_and_arm=reads,
                    weight=str(factor),
                    expected=str(factor * Q(31, 56)),
                    forecast_from_source=predicted,
                    conditional_line_tangent_nats=parity_target,
                    conditional_line_tangent_measured_99=parity_band,
                    measured=measured,
                    simultaneous_99=obs,
                )
            )
            assert obs[0] > 0 and obs[0] < model < obs[1]
    return results


def ordinary_control():
    # Start on blank records, with an ordinary double at the middle home or
    # an empty home. Measure reference and endpoints on separate trials.
    words = np.zeros((2, 7), dtype=int)
    words[1, 1] = 2
    _, _, _, _, _, energy = balls.energy.returned(words)
    rng = np.random.default_rng(513938)
    counts = [int(rng.binomial(SHOTS, p)) for p in (0.25, 0.25, 0.75)]
    ref, before, after = np.array(counts) / SHOTS
    eps = np.sqrt(np.log(6 / 0.01) / (2 * SHOTS))

    def odds(p):
        return np.log((1 - p) / p)

    interval = [
        (after - before - 2 * eps) * odds(ref + eps),
        (after - before + 2 * eps) * odds(ref - eps),
    ]
    tangent = 0.5 * np.log(3)
    assert 0 < interval[0] < tangent < interval[1]
    return dict(
        shots_per_bank=SHOTS,
        double_counts=counts,
        tangent_nats=tangent,
        measured_tangent_99=interval,
        original_to_final_whole_energy=str(Q(int(energy[1] - energy[0]), 8)),
    )


def inventory_spans():
    _, classes = balls.energy.pair.source.classes()
    result = {}
    for k, words in enumerate(classes):
        assert np.all(words[:, :4].sum(1) == 2 * k)
        b, home, seams, _, _, diagonal = balls.energy.returned(words)
        m = balls.energy.compressed_hamiltonian(words, b, home, seams, diagonal)
        radius = np.abs(m).sum(axis=1) - np.abs(diagonal)
        low, high = (int((diagonal - radius).min()), int((diagonal + radius).max()))
        if k == 0:
            assert not np.any(m)
        result[str(k)] = dict(
            words=len(words),
            energy_band=[str(Q(low, 4)), str(Q(high, 4))],
            span_upper=str(Q(high - low, 4)),
            span_per_double_upper=str(Q(high - low, 4 * k)) if k else None,
        )
    return result


def run():
    _, w, cells, _, _, h, *_ = balls.energy.pair.setup()
    board, homes, seams, terms, private, _ = balls.energy.returned(w)
    _, hops, _ = balls.radar.chart(board, seams)
    distances = np.column_stack(
        (hops[:, homes], *[hops[:, [u, v]].max(1) for u, v in seams])
    )
    rank = {str(r): comparisons(hops, distances, r) for r in (2, 4, 6)}
    mid = homes[1]
    remaining = np.delete(distances[mid], 1)
    last = int(min(remaining))
    assert last < hops[mid].max()
    # Any repair with the specified marginals and mean has the same tangent.
    counts = np.bincount(w[:, 1], minlength=3)
    assert counts.tolist() == [24, 400, 24]
    a = np.log(counts[1] / counts[0]) / 14
    z = balls.energy.least_change(cells)
    port = balls.energy_at_ports(board, terms, private, z)
    # Assemble a second time in rational arithmetic for exact weighted work.
    delta = (z - 1) / len(w)
    exact = delta @ private / 4
    for (u, v), price in zip(board.seams, terms @ delta / 8, strict=True):
        exact[u] += price
        exact[v] += price
    np.testing.assert_allclose(np.array(exact, float), port)
    return dict(
        sha256={
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (
                HERE / "NESTED_SOURCE.md",
                HERE / "closed_energy.json",
                HERE / "closed_balls.json",
            )
        },
        equal_record_comparisons=rank,
        middle=dict(
            center=mid,
            first_radius=2,
            last_radius=last,
            counts=counts.tolist(),
            tangent_nats=a,
            minimum_added_tangent_nats=(last / 2 - 1) * a,
        ),
        shots_per_arm_per_setting=SHOTS,
        spectral_readings=spectra(w, cells, hops, h),
        whole_work_at_middle={
            str(r): str(weights(hops[mid], r) @ exact / (2 * r))
            for r in range(2, last + 1)
        },
        total_work=str(sum(exact)),
        inventory_spans=inventory_spans(),
        ordinary_control=ordinary_control(),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "nested_source.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    print(out["middle"], out["inventory_spans"])
