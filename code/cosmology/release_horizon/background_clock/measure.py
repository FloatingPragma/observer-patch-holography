"""Read pulse periods, clock coincidences, census changes and thermal events on one screen."""

import argparse
from fractions import Fraction as Q
import json
from math import log, sqrt, prod
from pathlib import Path

import numpy as np

if __package__:
    from .material import HERE, clock_gates, instructions, screen, stream
else:
    from material import HERE, clock_gates, instructions, screen, stream
from horizon import upper, trigger_histories
from relay import gates

CONFIG = json.loads((HERE / "protocol.json").read_text())


def scale(level, t):
    return min(2**level, 1 if t <= 8 else 2 if t <= 16 else 4)


def census(level, qs, flip, fixed=False, repeat=False, mute=None):
    group = 4**level
    cells = np.arange(20 * group)
    marker = np.zeros(len(cells), dtype=np.int8)
    elapsed = [np.zeros(len(cells), dtype=int)]
    emitted = []
    for layer, q in enumerate(qs, 1):
        stride = group if fixed else group // (q * q)
        enrolled = cells[cells % stride == stride - 1]
        counts = np.zeros(20, dtype=int)
        local_elapsed = np.zeros(len(cells), dtype=int)
        for stroke in range(q):
            previous = marker.copy()
            marker = np.take(flip, previous)
            if mute and (layer, stroke) == mute[:2]:
                marker[mute[2]] = previous[mute[2]]
            changed = marker != previous
            local_elapsed += changed
            np.add.at(counts, enrolled // group, changed[enrolled].astype(int))
            if repeat:
                reread = marker.copy()
                np.add.at(counts, enrolled // group, (reread != marker)[enrolled])
        elapsed.append(elapsed[-1] + local_elapsed)
        emitted.append(counts.tolist())
    return np.stack(elapsed), emitted


def crest_times(trace, owner):
    return np.flatnonzero(trace[:, owner]).tolist()


def signal_readings(trace, receivers, elapsed):
    times = crest_times(trace, 0)
    records = {}
    for v in receivers:
        arrivals = crest_times(trace, v)
        local = elapsed[arrivals, v].tolist()
        source = elapsed[times, 0].tolist()
        periods = [
            str(Q(b - a, d - c))
            for a, b, c, d in zip(local, local[1:], source, source[1:])
        ]
        records[str(v)] = dict(
            arrival_layers=arrivals,
            local_clocks=local,
            period_ratios=periods,
            received=len(arrivals),
            missing=4 - len(arrivals),
        )
    return dict(
        emission_layers=times,
        emitter_clock=elapsed[times, 0].tolist(),
        receivers=records,
    )


def compare(trace, target, seam, elapsed, half, flip, mode):
    # Separate comparison bank. Its neighbor uses a free dimer and leaves both
    # intervening half steps coherent. Counts come from successive sharp pointers.
    left = right = 0
    lpos = rpos = 0
    for _ in range(2):
        nxt = int(flip[lpos])
        left += int(nxt != lpos)
        lpos = nxt
    for _ in range(3):
        nxt = int(flip[rpos])
        right += int(nxt != rpos)
        rpos = nxt
    state = np.eye(2, dtype=complex)[:, rpos]
    seen = []
    attempts = []
    quiet = 0
    invalid = 0
    skip = 0
    drift_pending = False
    total = 0
    for layer, bit in enumerate(trace[:, target]):
        steps = elapsed[layer] - elapsed[layer - 1] if layer else 0
        for _ in range(steps):
            if drift_pending and total % 2 == 0:
                skip = 2
                drift_pending = False
            total += 1
            nxt = int(flip[lpos])
            left += int(nxt != lpos)
            lpos = nxt
            if skip:
                skip -= 1
            else:
                state = half @ state
            if total % 2 == 0:
                obs = int(np.argmax(abs(state)))
                np.testing.assert_allclose(
                    abs(state) ** 2, np.eye(2)[:, obs], atol=4e-12
                )
                right += int(obs != rpos)
                rpos = obs
                state = np.eye(2, dtype=complex)[:, obs]
        signs = instructions(*seam)[1][int(mode != "cut"), int(bit)]
        attempts.append([layer, left, right, *signs])
        if not signs:
            invalid += 1
        elif signs == [0, 0]:
            quiet += 1
        elif signs == [1, 1]:
            # A crest between interrogations reads the retained completed-period count.
            seen.append([layer, left, right])
            if mode == "drift" and len(seen) == 2:
                drift_pending = True
        else:
            raise ValueError("Unpaired common record")
    fit = None
    if len(seen) >= 2:
        (_, x, y), (_, xx, yy) = seen[:2]
        rate = Q(yy - y, xx - x)
        bias = y - rate * x
        fit = dict(
            rate=str(rate),
            offset=str(bias),
            residual=[str(b - rate * a - bias) for _, a, b in seen[2:]],
        )
    return dict(
        events=seen, attempts=attempts, quiet=quiet, invalid=invalid, calibration=fit
    )


def thermal(bins, cells, elapsed, originals, seed, cadence):
    rng = np.random.default_rng(seed)
    words = np.arange(256)
    permutations = [
        words ^ ((((words >> k) ^ (words >> (k + 1))) & 1) * (3 << k))
        for k in (0, 2, 4, 6)
    ]
    init = [
        prod(p if w & (1 << j) else 1 - p for j, p in enumerate((1 / 5, 1 / 3) * 4))
        for w in words
    ]
    terminal = sum(((words >> (2 * j)) & 1) << j for j in range(4))
    output = np.zeros(16, np.int64)
    means = [Q(0)] * 3
    mass = Q(1)
    remaining = originals
    by_age = []
    for t, p, trace in bins:
        active = np.ones(4, bool)
        offers = np.zeros(4, int)
        for s in range(1, len(trace)):
            offers += active * (
                1 if cadence == "layer" else elapsed[s, cells] - elapsed[s - 1, cells]
            )
            active &= ~trace[s, cells].astype(bool)
        conditional = [upper(int(n)) for n in offers]
        a, b = conditional[0] - conditional[1], conditional[2] - conditional[3]
        means = [x + p * y for x, y in zip(means, (a, b, a * b), strict=True)]
        draws = remaining if mass == p else rng.binomial(remaining, float(p / mass))
        state = rng.multinomial(draws, init)
        for n, perm in zip(offers, permutations, strict=True):
            for _ in range(n):
                accepted = rng.binomial(state, 0.25)
                state = state - accepted + accepted[perm]
        final = np.zeros(16, np.int64)
        np.add.at(final, terminal, state)
        output += final
        by_age.append(
            dict(
                trigger=t,
                probability=str(p),
                offers=offers.tolist(),
                misses=active.tolist(),
                energy_words=final.tolist(),
            )
        )
        mass -= p
        remaining -= draws
    assert mass == remaining == 0 and output.sum() == originals
    codes = np.arange(16)
    a = (codes & 1) - ((codes >> 1) & 1)
    b = ((codes >> 2) & 1) - ((codes >> 3) & 1)
    observed = [float(output @ z) / originals for z in (a, b, a * b)]
    width = sqrt(2 * log(2 * 6 / 0.01) / originals)
    assert all(abs(float(x) - y) <= width for x, y in zip(means, observed, strict=True))
    return dict(
        ages=by_age,
        energy_words=output.tolist(),
        covariance=str(means[2] - means[0] * means[1]),
        observed=observed[2] - observed[0] * observed[1],
        simultaneous_99pct_halfwidth=width * (1 + abs(observed[0]) + abs(observed[1]))
        + width**2,
    )


def run():
    gates()
    half, flip, park = clock_gates()
    records = []
    for level in CONFIG["levels"]:
        board, levels, radii, seams, neighbors, parent = screen(level)
        horizon = CONFIG["deadline"]
        qs = [scale(level, k) for k in range(1, horizon + 1)]
        elapsed, face_counts = census(level, qs, flip)
        flat, flat_counts = census(level, [1] * horizon, flip)
        fixed_clock, fixed = census(level, qs, flip, fixed=True)
        np.testing.assert_array_equal(fixed_clock, elapsed)
        parked, silent = census(level, qs, park)
        assert not np.any(parked) and not np.any(silent)
        rep, rep_counts = census(level, qs, flip, repeat=True)
        np.testing.assert_array_equal(rep, elapsed)
        assert rep_counts == face_counts
        d = int(radii.max())
        receivers = [
            int(np.flatnonzero(radii == r)[0])
            for r in sorted({min(i, d) for i in (4, 8, 16)})
        ]
        target = receivers[len(receivers) // 2]
        witness = min(neighbors[target])
        seam = seams[target, witness]
        traces = {}
        for t in range(horizon + 2):
            states = stream(levels, board.cells, horizon, t)
            expected = np.isin(
                np.arange(horizon + 1)[:, None] - t - radii[None, :], (0, 2, 4, 6)
            )
            assert np.array_equal(states, expected)
            traces[t] = states
        pilot = traces[0]
        cut = stream(levels, board.cells, horizon, 0, cut=True)
        assert not np.any(cut[:, 1:])
        muted, _ = census(
            level, qs, flip, mute=(int(radii[receivers[-1]]) + 1, 0, receivers[-1])
        )
        cases = {
            name: compare(
                pilot, target, seam, elapsed[:, target].tolist(), half, flip, name
            )
            for name in ("native", "drift", "cut")
        }
        assert cases["native"]["calibration"]["residual"] == ["0", "0"]
        assert cases["drift"]["calibration"]["residual"][-1] == "-1"
        assert cases["cut"]["calibration"] is None
        receipts = []
        for t, p in trigger_histories(horizon):
            receipts.append(
                dict(
                    trigger=t,
                    probability=str(p),
                    **signal_readings(traces[t], receivers, elapsed),
                )
            )
        row = dict(
            cells=board.cells,
            receivers=receivers,
            radii=[int(radii[i]) for i in receivers],
            comparison_seam=[12 * target + seam[0], 12 * witness + seam[1]],
            comparison=cases,
            clock=elapsed[:, 0].tolist(),
            counts_per_face=face_counts,
            static_counts_per_face=flat_counts[0],
            fixed_enrollment_per_layer=[v[0] for v in fixed],
            parked_census_count=0,
            muted_receiver=signal_readings(pilot, receivers, muted),
            commanded=signal_readings(pilot, receivers, elapsed),
            static=signal_readings(pilot, receivers, flat),
            cut=signal_readings(cut, receivers, elapsed),
            trigger_mixture=receipts,
        )
        if level == 2:
            archive = json.loads((HERE.parent / "results.json").read_text())
            reference = archive["boards"][2]["deadlines"][0]["local_differences"]
            cells = reference["cells"]
            deadline = 25
            bins = [
                (t, p, traces[t][: deadline + 1])
                for t, p in trigger_histories(deadline)
            ]
            # At the tail t=26 no pulse can have entered a t<=25 trace.
            row["thermal"] = dict(
                cells=cells, deadline=deadline, marker=int(elapsed[deadline, cells[0]])
            )
            for idx, cadence in enumerate(("layer", "microstep")):
                row["thermal"][cadence] = thermal(
                    bins,
                    cells,
                    elapsed,
                    CONFIG["contrast_originals"],
                    CONFIG["seed"] + idx,
                    cadence,
                )
            assert (
                row["thermal"]["layer"]["covariance"]
                == reference["banks"]["earlier_release"]["covariance"]
            )
            # Propagate the actual supported record dependencies at a comparison
            # origin two layers before the reading. Each target also retains its history.
            past = [{i} for i in range(board.cells)]
            for _ in range(2):
                past = [old | past[parent[i]] for i, old in enumerate(past)]
            # Independent interventions. Each column seeds a single material
            # crest at a different owner; the parked root adds no later crest.
            # Record every hit on either endpoint, including at the reset.
            pulses = np.eye(board.cells, dtype=np.int8)
            delivered = pulses.astype(bool)
            for step in range(2):
                for parents, children, tables in levels:
                    old_parents = pulses[parents, :]
                    old_children = pulses[children, :]
                    pulses[children, :] = tables[
                        np.arange(len(children))[:, None], 1, old_parents, old_children
                    ]
                pulses[0, :] = 0
                delivered |= pulses.astype(bool)
            measured = [np.flatnonzero(delivered[v]).tolist() for v in cells]
            assert measured == [sorted(past[v]) for v in cells]
            row["thermal"]["reset_layer"] = 23
            row["thermal"]["reset_marker"] = int(elapsed[23, cells[0]])
            row["thermal"]["pulse_inputs_per_probe"] = measured
            a, b, c, e = cells
            row["thermal"]["shared_postreset_owners"] = sorted(
                (past[a] | past[b]) & (past[c] | past[e])
            )
        records.append(row)
    return dict(boards=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), separators=(",", ":")) + "\n"
    if args.out:
        args.out.write_text(data)
    elif json.loads(data) != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("Archive mismatch")
