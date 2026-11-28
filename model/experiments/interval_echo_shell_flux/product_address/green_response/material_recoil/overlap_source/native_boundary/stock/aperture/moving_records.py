"""Calibrate time-ordered material records and independent messages across regional cuts."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import moving_service as lab
from . import moving_current
from . import branch_balance
from . import thermal_radar

HERE = lab.energy.HERE
SOURCE_N = 2**29
METER_N = 2**36
RADII = (2, 4, 6, 8, 10, 12)


def ensemble(words, weights, extra=0, single=0):
    dist = defaultdict(Q)
    for w, p in zip(words, weights, strict=True):
        dist[lab.signature(w, extra, single)] += p
    return dict(dist)


def summary(pre):
    bulk, material, tv = lab.readers()
    bins = (pre @ bulk.T) / 8
    return dict(
        bins=bins,
        entropy=lab.bits(bins),
        material=lab.bits(pre @ material.T / 8),
        signals=pre @ tv.T / 16,
    )


def report(x):
    return dict(
        bits=x["entropy"].tolist(),
        material_given_start_bits=x["material"].tolist(),
        surprise_bins=[[str(p) for p in row] for row in x["bins"]],
        signal_advantage=[[str(p) for p in row] for row in x["signals"]],
    )


def meter(pre, rng):
    full, half = lab.gates()
    out = []
    for row in pre:
        prob = np.asarray(row, float)
        starts = rng.multinomial(METER_N, prob / prob.sum())
        surprise = np.zeros(8, dtype=np.int64)
        for col in np.flatnonzero(starts):
            counts = rng.multinomial(int(starts[col]), half[:, col] / 8)
            good = half[:, col] > 0
            np.add.at(surprise, half[good, col] - 1, counts[good])
        advantages = []
        for t in (1, 2):
            # Choose the sender bit and orientation independently in each
            # fresh trial, then perform the full native gate and decode.
            correct = 0
            for col, total in enumerate(rng.multinomial(METER_N, prob / prob.sum())):
                if not total:
                    continue
                m, ab = divmod(col, 9)
                a, b = divmod(ab, 3)
                for direction, parts in enumerate(
                    rng.multinomial(int(total), [0.25] * 4).reshape(2, 2)
                ):
                    inputs = (
                        [9 * m + b, 9 * m + 3 * t + b]
                        if direction == 0
                        else [9 * m + 3 * a, 9 * m + 3 * a + t]
                    )
                    columns = full[:, inputs].reshape(11, 3, 3, 2).sum(0).sum(direction)
                    assert np.all(columns.sum(0) == 4)
                    decode = (columns[:, 1] > columns[:, 0]).astype(int)
                    for bit, qty in enumerate(parts):
                        freq = rng.multinomial(int(qty), columns[:, bit] / 4)
                        correct += int(freq[decode == bit].sum())
            advantages.append(correct)
        # Fresh full-pulse draws compare the returned two loads and head.
        returns = 0
        for col, n in enumerate(rng.multinomial(METER_N, prob / prob.sum())):
            if n:
                counts = rng.multinomial(int(n), full[:, col] / 4)
                returns += int(counts[col])
        out.append(
            dict(
                half_surprise_bins=surprise.tolist(),
                full_correct=advantages,
                returned_word=returns,
            )
        )
    return out


def sweep_interval(bins, counts):
    from mpmath import iv

    iv.dps = 40
    val = iv.mpf(0)
    for i, q in enumerate(counts @ bins, 1):
        val += (
            iv.mpf(q.numerator) / int(q.denominator) * iv.log(iv.mpf(8) / i) / iv.log(2)
        )
    return [
        float(np.nextafter(float(val.a), -np.inf)),
        float(np.nextafter(float(val.b), np.inf)),
    ]


def measured(banks):
    ent = np.array(
        [lab.bits(np.asarray(row["half_surprise_bins"]) / METER_N) for row in banks]
    )
    sig = 2 * np.array([row["full_correct"] for row in banks]) / METER_N - 1
    return ent, np.column_stack(
        (sig, [row["returned_word"] / METER_N for row in banks])
    )


def history_run(board, homes, marks, word):
    # A finite trajectory begins at the actual returned preparation. Sampling
    # an equilibrium response is a separate calculation below.
    rng = np.random.default_rng(7212301)
    loads = np.zeros(board.n, dtype=int)
    loads[homes] = word[:4]
    head = np.zeros(len(board.seams), int)
    for mark, j in zip(marks, word[4:], strict=True):
        head[board.seams.index(mark)] = j
    initial = np.bincount(loads, minlength=3)
    kernel = lab.gates()[1]
    gates = np.cumsum(kernel, axis=0)
    digest = sha256()
    moved = written = 0
    for order in range(16):
        for idx in rng.permutation(len(board.seams)):
            a, b = board.seams[idx]
            col = 9 * head[idx] + 3 * loads[a] + loads[b]
            row = int(np.searchsorted(gates[:, col], rng.integers(8), side="right"))
            digest.update(bytes((a, b, col, row)))
            moved += int(row % 9 != col % 9)
            written += int(row // 9 != col // 9)
            head[idx], ab = divmod(row, 9)
            loads[a], loads[b] = divmod(ab, 3)
        np.testing.assert_array_equal(np.bincount(loads, minlength=3), initial)
    return dict(
        appointments=(order + 1) * len(board.seams),
        changed_load_pairs=moved,
        changed_heads=written,
        ordered_records_sha256=digest.hexdigest(),
    )


def surveys(hops, edges, types, orig, source, ordinary, single):
    eo, es = (x["entropy"] - orig["entropy"] for x in (ordinary, source))
    ao, ass = (x["signals"] - orig["signals"] for x in (ordinary, source))
    single_e = single["entropy"] - orig["entropy"]
    single_a = single["signals"] - orig["signals"]
    result = {}
    examples = []
    gains = []
    coefficients = [
        lab.log_coefficients(x["bins"] - orig["bins"])
        for x in (ordinary, source, single)
    ]
    for r in RADII:
        stations = {}
        for i, h in enumerate(hops):
            inside = h < r
            ends = inside[edges]
            interior = np.bincount(types[ends.all(1)], minlength=4)
            cut = np.bincount(types[ends[:, 0] ^ ends[:, 1]], minlength=4)
            entropy = np.array([interior @ eo, interior @ es, interior @ single_e])
            area = np.vstack((cut @ ao, cut @ ass, cut @ single_a))
            stations[i] = (interior, cut, entropy, area)
        panels = {}
        for t in (0, 1, 2):
            silent = sum(
                (a[0, t] == 0 and abs(e[0]) > 1e-12) for _, _, e, a in stations.values()
            )
            active = [(i, e, a) for i, (_, _, e, a) in stations.items() if a[0, t] != 0]
            # Stationarity with a positive area entropy would imply opposite signs.
            pos = [(i, e, a) for i, e, a in active if float(a[0, t]) * e[0] < 0]
            errors = [abs(e[1] - e[0] * float(a[1, t] / a[0, t])) for _, e, a in pos]

            def same(i, held):
                cnt, _, _, areas = stations[i]
                cross = (cnt @ coefficients[held]) * areas[0, t] - (
                    cnt @ coefficients[0]
                ) * areas[held, t]
                return all(x == 0 for x in cross)

            ordinary_errors = [
                abs(e[2] - e[0] * float(a[2, t] / a[0, t])) for _, e, a in pos
            ]
            if t == 2:
                gains.extend(-e[0] / float(a[0, t]) for _, e, a in pos)
            panels[("single_message", "double_message", "word_return")[t]] = dict(
                exact_held_source_matches=sum(same(i, 1) for i, _, _ in pos),
                exact_held_single_matches=sum(same(i, 2) for i, _, _ in pos),
                held_single_residual_above_1e_8=sum(
                    bool(v > 1e-8) for v in ordinary_errors
                ),
                zero_ordinary_area_nonzero_history=int(silent),
                ordinary_calibration_resolved=len(active),
                positive_ordinary_gain=len(pos),
                same_gain_residual_above_1e_8=sum(bool(v > 1e-8) for v in errors),
            )
        result[str(r)] = panels
        inner, cut, e, a = stations[20]
        examples.append(
            dict(
                center=20,
                radius=r,
                interior=inner.tolist(),
                boundary=cut.tolist(),
                ordinary_source_history_bits=e.tolist(),
                ordinary_source_signals=[[str(x) for x in row] for row in a],
            )
        )
    return dict(
        radii=result,
        middle=examples,
        positive_return_gain_range=[min(gains), max(gains)],
    )


def composition_certificate(w, cells, z, observed):
    from scipy.optimize import linprog
    from sympy import Matrix, linsolve, symbols

    n = lab.V
    m = (w[:, :4] == 2).sum(1)
    score = 8 * n - 32 + (24 - 2 * n) * m - 6 * m * m
    denom = n * (n - 1)
    base = np.full(len(w), Q(1, len(w)), object)
    delta = (z - base) @ score / denom
    assert delta == Q(44, 7 * denom)

    def empirical(arm):
        bank = observed[arm]
        return sum(
            int(k) * (n * n - sum(v * v for v in pop))
            for k, (pop, _) in zip(bank["counts"], bank["classes"], strict=True)
        ) / (SOURCE_N * denom)

    obs = empirical("source") - empirical("original")
    score_range = (int(score.max()) - int(score.min())) / denom
    error = score_range * np.sqrt(np.log(200) / SOURCE_N)
    assert obs - error > 0 and abs(obs - float(delta)) < error
    _, _, _, _, codes, h, ids, groups = lab.energy.pair.setup()
    total = lab.energy.returned(w)[-1]
    rows = np.vstack(
        ([cells == i for i in range(4)], w[:, :4].T, h.diagonal()[codes], total)
    )
    target = rows @ (z * len(w))
    A = lab.energy.pair.group_matrix(rows, ids, groups)
    cost = lab.energy.pair.group_matrix(score[None, :], ids, groups)[0]
    counter = {}
    for sign, name in ((1, "low"), (-1, "high")):
        opt = linprog(
            sign * cost,
            A_eq=A,
            b_eq=np.asarray(target, float),
            bounds=(0, None),
            method="highs",
        )
        assert opt.success
        keep = np.flatnonzero(opt.x > 1e-9)
        unknown = symbols(f"x:{len(keep)}")
        (sol,) = linsolve((Matrix(A[:, keep]), Matrix(target)), unknown)
        weights = np.full(groups, Q(0), dtype=object)
        weights[keep] = [Q(int(y.p), int(y.q)) for y in sol]
        assert all(weights >= 0) and all(A @ weights == target)
        shift = (cost @ weights - int(score.sum())) / (len(w) * denom)
        alternate = lab.before(ensemble(w, weights[ids] / len(w)))
        native = lab.before(ensemble(w, base))
        four_way = (alternate - native) @ branch_balance.classify()
        counter[name] = dict(
            delta=str(shift),
            group_weights=list(map(str, weights)),
            four_way_changes=list(map(str, four_way)),
        )
    assert Q(counter["low"]["delta"]) < 0 < Q(counter["high"]["delta"])
    assert any(sum(map(Q, row["four_way_changes"])) < 0 for row in counter.values())
    return dict(
        delta_m=str((z - base) @ m),
        delta_m_squared=str((z - base) @ (m * m)),
        blank_activity_delta=str(delta),
        estimate_from_four_home_counts=obs,
        measured_99=[obs - error, obs + error],
        sharp_score_range=score_range,
        ideal_accepted_per_arm=int(
            np.ceil(4 * score_range**2 * np.log(100) / float(delta) ** 2)
        ),
        same_charge_and_energy_witnesses=counter,
    )


def run():
    _, w, cells, *_ = lab.energy.pair.setup()
    board, home, marks, *_ = lab.energy.returned(w)
    assert board.n == lab.V
    edges = np.array(board.seams)
    _, hops, _ = thermal_radar.chart(board, marks)
    assert hops.max() < board.n
    types = np.array(
        [1 + marks.index(edge) if edge in marks else 0 for edge in board.seams]
    )
    numbers = np.bincount(types, minlength=4)
    uniform = np.full(len(w), Q(1, len(w)), object)
    z = lab.energy.least_change(cells) / len(w)
    source = ensemble(w, z)
    original = ensemble(w, uniform)
    ordinary = ensemble(w, uniform, extra=1)
    preparations = dict(
        original=original,
        source=source,
        ordinary=ordinary,
        single=ensemble(w, uniform, single=1),
    )
    rng = np.random.default_rng(41000321)
    certs = {}
    meters = {}
    predictions = {}
    # Source census: independent complete returned words. The invariant
    # composition and head components predict every later stationary meter.
    for arm, dist in preparations.items():
        p = np.array(list(dist.values()), float)
        counts = rng.multinomial(SOURCE_N, p / p.sum())
        empirical = {
            k: Q(int(cnt), SOURCE_N) for k, cnt in zip(dist, counts, strict=True)
        }
        census = summary(lab.before(empirical))
        before = lab.before(dist)
        certs[arm] = summary(before)
        meters[arm] = meter(before, rng)
        predictions[arm] = dict(
            counts=counts.tolist(),
            classes=[[list(pop), [list(c) for c in m]] for pop, m in dist],
            entropy=census["entropy"].tolist(),
            signals=np.asarray(census["signals"], float).tolist(),
        )
    records = {}
    # Three measured-arm differences for four edge types and four
    # readouts, within each of the source-census and terminal-meter families.
    # Allocate half of the total error to each bank family.
    for tag in ("source", "ordinary", "single"):
        change = certs[tag]["entropy"] - certs["original"]["entropy"]
        sig = certs[tag]["signals"] - certs["original"]["signals"]
        ent, s = measured(meters[tag])
        e0, s0 = measured(meters["original"])
        held = np.column_stack((ent - e0, s - s0))
        predicted = np.column_stack(
            (
                np.array(predictions[tag]["entropy"])
                - predictions["original"]["entropy"],
                np.array(predictions[tag]["signals"])
                - predictions["original"]["signals"],
            )
        )
        truth = np.column_stack((change, np.asarray(sig, float)))
        # Conditional-entropy per initial code is at most three bits. Terminal
        # advantage is a +/-1 score. Source-side expected advantage is in [0,1].
        halfwidth = np.array([3, 2, 2, 1]) * np.sqrt(np.log(2 * 48 / 0.005) / METER_N)
        forecastwidth = np.array([3, 1, 1, 1]) * np.sqrt(
            np.log(2 * 48 / 0.005) / SOURCE_N
        )
        assert np.all(abs(held - truth) <= halfwidth) and np.all(
            abs(predicted - truth) <= forecastwidth
        )
        enclosed = sweep_interval(
            certs[tag]["bins"] - certs["original"]["bins"], numbers
        )
        records[tag] = dict(
            held=held.tolist(),
            forecast=predicted.tolist(),
            expected=truth.tolist(),
            measured_halfwidth=halfwidth.tolist(),
            forecast_halfwidth=forecastwidth.tolist(),
            sweep_bits=sum(enclosed) / 2,
            sweep_bits_interval=enclosed,
        )
    spectra = lab.support_bounds(lab.energy.pair.source.classes()[1], numbers[0].item())
    base = Q(2 * len(edges), lab.V)
    # With one ordinary double in a blank board, the entire path entropy
    # is one bit whenever that double lies on the selected contact.
    blank = np.zeros((1, 7), int)
    b0, b1 = [summary(lab.before(ensemble(blank, [Q(1)], extra=e))) for e in (0, 1)]
    np.testing.assert_allclose(numbers @ (b1["entropy"] - b0["entropy"]), float(base))
    assert np.all(b1["signals"][:, :2] == b0["signals"][:, :2])
    return dict(
        sha256={
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (
                HERE / "MOVING_RECORDS.md",
                HERE / "nested_source.json",
                HERE / "closed_pairing.json",
            )
        },
        edges_by_type=numbers.tolist(),
        marks=marks,
        source_count=SOURCE_N,
        meter_count=METER_N,
        blank_double_bits=str(base),
        spectrum=spectra,
        sharp_trajectory=history_run(board, home, marks, w[0]),
        source_census=predictions,
        branch_balance=branch_balance.run(
            preparations, numbers, edges, types, hops, certs, records, RADII
        ),
        composition=composition_certificate(w, cells, z, predictions),
        currents=moving_current.run(
            numbers, lab.before(original), lab.before(source), board
        ),
        meters=meters,
        stationary_laws={k: report(x) for k, x in certs.items()},
        changes=records,
        regional=surveys(
            hops,
            edges,
            types,
            certs["original"],
            certs["source"],
            certs["ordinary"],
            certs["single"],
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opts = cli.parse_args()
    out = run()
    path = HERE / "moving_records.json"
    if opts.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        thermal_radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    print({k: v["sweep_bits"] for k, v in out["changes"].items()}, out["spectrum"])
