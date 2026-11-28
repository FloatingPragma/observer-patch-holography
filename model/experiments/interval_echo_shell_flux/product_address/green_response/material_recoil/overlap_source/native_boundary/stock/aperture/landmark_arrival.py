"""Two independent beacons read ordinary material at first arrival on the native walk."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.sparse.csgraph import laplacian, shortest_path
from sympy import Matrix, ZZ
from sympy.polys.matrices import DomainMatrix

from . import closed_passage as passage
from . import collective_source as source
from . import landmark_release as ranging
from . import landmark_preparation
from . import landmark_calibration

HERE = source.tidal.HERE


def blank_exits(adj, pins):
    n = len(adj)
    interior = np.array([i for i in range(n) if i not in pins])
    L = laplacian(adj)

    def matrix(x):
        return DomainMatrix.from_list(x.tolist(), ZZ).to_sparse()

    A = matrix(L[np.ix_(interior, interior)])
    B = matrix(-L[np.ix_(interior, pins)])
    X, den = A.solve_den(B)
    assert A * X == B * den
    F = np.zeros((n, 4), object)
    F[interior] = X.to_list()
    F[pins, np.arange(4)] = int(den)
    assert np.all(F >= 0) and np.all(F.sum(1) == den)
    return F, int(den)


def eventual(table, pins, F, den):
    # The last two pins enclose the only head-changing material crossing.
    n = len(F)
    A = np.eye(22, dtype=int).astype(object) * 24 * den
    rhs = np.zeros(22, object)
    for m in range(11):
        for j, x in enumerate(pins[2:]):
            i = 2 * m + j
            for dest in table[m * n + x]:
                h, y = divmod(int(dest), n)
                rhs[i] += F[y, 1]
                A[i, 2 * h : 2 * h + 2] -= F[y, 2:]
    solution = Matrix(A.tolist()).inv() * Matrix(rhs.tolist())
    heads = np.array([Q(int(p.p), int(p.q)) for p in solution], object).reshape(11, 2)
    values = (F[:, 1] + heads @ F[:, 2:].T) / Q(den)
    assert np.all(values >= 0) and np.all(values <= 1)
    assert np.all(values[:, pins[0]] == 0) and np.all(values[:, pins[1]] == 1)
    live = np.tile(~np.isin(np.arange(n), pins[:2]), 11)
    flat = values.ravel()
    assert np.all(flat[table[live]].sum(1) == 24 * flat[live])
    return values


def held(table, slots, pins, seed, starts, shots=2**15):
    n = len(slots)
    a, b, u, v = pins
    k = 10
    block = np.eye(n, dtype=np.int64)
    stopped = slots.copy()
    stopped[pins] = np.array(pins)[:, None]
    for _ in range(k):
        block = block[stopped].sum(1)
    tickets = 6**k
    cdf = (block.cumsum(1) + np.arange(n)[:, None] * tickets).ravel()
    assert np.all(block.sum(1) == tickets)
    rng = np.random.default_rng(seed)
    result = {}
    for x in starts:
        label = np.repeat(np.arange(11), shots)
        state = label * n + x
        wins = np.zeros(11, int)
        losses = np.zeros(11, int)
        for _ in range(12000):
            if not len(state):
                break
            head, site = np.divmod(state, n)
            at = (site == u) | (site == v)
            if np.any(at):
                state[at] = table[state[at], rng.integers(24, size=int(at.sum()))]
            free = ~at
            if np.any(free):
                ticket = site[free] * tickets + rng.integers(
                    tickets, size=int(free.sum())
                )
                state[free] = (
                    head[free] * n + np.searchsorted(cdf, ticket, side="right") % n
                )
            site = state % n
            win = site == b
            lose = site == a
            wins += np.bincount(label[win], minlength=11)
            losses += np.bincount(label[lose], minlength=11)
            live = ~(win | lose)
            state, label = state[live], label[live]
        pending = shots - wins - losses
        err = np.sqrt(np.log(2 * 22 / 0.003) / (2 * shots))
        result[str(x)] = dict(
            b=wins.tolist(),
            a=losses.tolist(),
            pending=pending.tolist(),
            intervals=np.column_stack(
                (
                    np.maximum(0, wins / shots - err),
                    np.minimum(1, (wins + pending) / shots + err),
                )
            ).tolist(),
        )
    return dict(trials_per_head=shots, seed=seed, max_blocks=12000, launches=result)


def pulse_inputs():
    words, weights, _ = source.records()
    board, _, marks, *_ = source.tidal.lab.energy.returned(words)
    transfers = passage.full_pulses(words[:, 4:])
    for m in np.unique(words[:, 4:6]):
        assert all(transfers[str(t)][m] == [0, 0] for t in (1, 2))
    _, half, _, _ = passage.pair.stock.kernels()
    gate = half // 8 - 4 * np.eye(99, dtype=int)
    tables = {t: passage.pulse_table(board, marks, t, gate, False) for t in (1, 2)}
    blind = passage.pulse_table(board, marks, 2, gate, True)
    # On the changed law every head and both load values communicate by
    # the same ordinary permutation, including on every written seam.
    _, free = passage.contacts(board, [])
    for m, rows in enumerate(blind.reshape(11, board.n, 24)):
        assert np.array_equal(rows, m * board.n + np.repeat(free, 4, axis=1))
    mass = {
        key: np.array([sum(w[words[:, -1] == m]) for m in range(11)], object)
        for key, w in weights.items()
    }
    return board, marks, tables, mass


def check_two_layers(table, slots, edge):
    # Collapse the four partially transmitting labels by which end can
    # send across the active contact. Check all 24 tickets at every port.
    n = len(slots)
    u, v = edge
    labels = {7: 0, 9: 0, 2: 1, 4: 1}
    for m, sheet in labels.items():
        for x in range(n):
            target = table[m * n + x]
            layer = np.array([labels[int(dest // n)] for dest in target])
            actual = np.bincount(n * layer + target % n, minlength=2 * n)
            expected = np.bincount(n * sheet + slots[x], minlength=2 * n) * 4
            if (sheet, x) in ((0, u), (1, v)):
                expected[n * sheet + x] -= 2
                expected[n * (1 - sheet) + (v if x == u else u)] += 2
            assert np.array_equal(actual, expected)


def landmark_ratios(adj, contact, da, db):
    # Ground the second beacon, then impose the first as a Dirichlet zero.
    # For beacon pair A,B the resulting Green entries have denominator
    # 4*d*Rraw[A,B], where d is the raw bypass resistance denominator.
    H, d = ranging.grounded_inverse(adj)
    raw, d = ranging.ranges((H, d), [], [])
    u, v = contact
    ratios = {}
    silent = 0
    for b in range(len(adj)):
        if b in contact:
            continue
        for a in range(b):
            if a in contact:
                continue
            length = raw[a, b]
            au = raw[u, b] + length - raw[u, a]
            av = raw[v, b] + length - raw[v, a]
            uu = 4 * length * raw[u, b] - au * au
            vv = 4 * length * raw[v, b] - av * av
            uv = 2 * length * (raw[u, b] + raw[v, b] - raw[u, v]) - au * av
            den = 4 * d * length
            assert min(uu, vv, uv) >= 0
            weight = Q(den + uu + vv - 2 * uv, 2 * (2 * den + uu + vv))
            assert 0 < weight < Q(1, 2)
            if au == av:
                silent += 1
            else:
                ratios[a, b] = 1 + db * weight / da
    lo = min(ratios, key=ratios.get)
    hi = max(ratios, key=ratios.get)
    return dict(
        silent_beacon_pairs=silent,
        sensitive_pairs=len(ratios),
        minimum_ratio=str(ratios[lo]),
        maximum_ratio=str(ratios[hi]),
        minimum_pair=list(lo),
        maximum_pair=list(hi),
    ), ratios


def finite_arrivals(tables, pins, x, measures):
    # Hold the physical word through all appointments. Both beacons absorb
    # on first arrival. The complement contains all trips in transit.
    out = {}
    for t, table in tables.items():
        size = len(table) // 11
        F = np.zeros((len(table), 2), object)
        F[np.arange(11) * size + pins[0], 0] = 1
        F[np.arange(11) * size + pins[1], 1] = 1
        terminal = F.sum(1) == 1
        T = table.copy()
        T[terminal] = np.flatnonzero(terminal)[:, None]
        rows = {}
        for tick in range(1, 65):
            F = F[T].sum(1)
            if tick in (4, 8, 16, 32, 64):
                at = F[np.arange(11) * size + x]
                assert np.all(at.sum(1) <= 24**tick)
                contrast = measures @ at / Q(24**tick)
                rows[str(tick)] = list(map(str, [*contrast, -sum(contrast)]))
        out[str(t)] = rows
    assert out["1"] == out["2"]
    return out


def parked_profiles(adj, slots, marks, eligible, a, da):
    H, d = ranging.grounded_inverse(adj)
    closed, _ = ranging.ranges((H, d), [], [])
    u, v = marks[-1]
    Z = H[:, u] - H[:, v]
    resistance = Q(Z[u] - Z[v], d)
    endpoints = np.ravel(marks)
    center = int(shortest_path(adj, unweighted=True, indices=endpoints).max(0).argmin())
    # Mean change of range in one monitored six-slot material tick. Eligible
    # blank stars have the same local ticks on all the fixed-head preparations.
    strength = da / (6 * (1 + resistance))
    change = [
        -strength * Q(sum((Z[y] - Z[x]) ** 2 for y in row), d * d)
        for x, row in enumerate(slots)
    ]
    profiles = {}
    for anchor in (center, u, v):
        radii = [
            Q(s, d) - a * Q((Z[anchor] - z) ** 2, d * d) / (1 + resistance)
            for s, z in zip(closed[anchor], Z, strict=True)
        ]
        profiles[str(anchor)] = ranging.shell_profile(radii, eligible, change, None)
    return profiles


def summarize(v, eligible):
    vals = v[eligible]
    return dict(
        launches=len(vals),
        toward_b=int((vals > 0).sum()),
        toward_a=int((vals < 0).sum()),
        silent=int((vals == 0).sum()),
        minimum=str(min(vals)),
        maximum=str(max(vals)),
    )


def run():
    board, marks, tables, mass = pulse_inputs()
    adj, slots = passage.contacts(board, marks)
    distance = shortest_path(
        adj, directed=False, unweighted=True, indices=np.ravel(marks)
    )
    eligible = distance.min(0) >= 2
    delta = mass["source"] - mass["original"]
    assert sum(delta) == 0
    rng = np.random.default_rng(995025)
    shots = 2**24
    acquired = {k: rng.multinomial(shots, np.array(p, float)) for k, p in mass.items()}
    # Total variation of the difference bounds its action on [0,1].
    # With zero-sum mass errors the half L1 norm suffices.
    mass_err = 11 * np.sqrt(np.log(2 * 22 / 0.003) / (2 * shots))
    measured = (acquired["source"] - acquired["original"]) / shots
    active = list(source.PARTIAL)
    da, db = delta[0], sum(delta[active])
    assert db / da == -Q(9, 7)
    assert all(len(set(weights[active])) == 1 for weights in mass.values())
    check_two_layers(tables[1], slots, marks[-1])
    survey, ratios = landmark_ratios(adj, marks[-1], da, db)
    matched, accepted_weights = landmark_preparation.prepare()
    prepared = np.array(
        list(map(Q, matched["preparations"]["entropy_stationary"]["head_fractions"])),
        object,
    )
    matched_delta = prepared - mass["original"]
    matched_a = matched_delta[0]
    # A directly recorded bounded score for each end of the interval
    # [0,1/2] for the partially transmitting contribution.
    margin = np.sqrt(np.log(4 / 0.003) / shots)
    scores = [measured[0] + w * sum(measured[active]) for w in (0, 0.5)]
    effective_bounds = [min(scores) - margin, max(scores) + margin]
    assert effective_bounds[0] > 0
    panels = {}
    for a, b in ((25, 36), (104, 116)):
        pins = [a, b, *marks[-1]]
        F, den = blank_exits(adj, pins)
        field = {t: eventual(tab, pins, F, den) for t, tab in tables.items()}
        live = eligible & ~np.isin(np.arange(board.n), [a, b])
        exits = {t: delta @ v for t, v in field.items()}
        residual = exits[1] - exits[2]
        ratio = ratios[a, b]
        assert np.all(exits[1] == ratio * exits[2])
        common = matched_delta @ field[2]
        assert np.all(common == matched_delta @ field[1])
        assert np.all(common == matched_a * (field[2][0] - field[2][1]))
        # Independent native first-step solution versus the two-sheet formula.
        weight = (ratio - 1) * da / db
        blank, closed = field[2][0], field[2][1]
        assert np.all(field[2][1:] == closed)
        assert np.all(field[1][0] == blank)
        assert np.all(field[1][[1, 3, 5, 6, 8, 10]] == closed)
        assert np.all(field[1][active].sum(0) / 4 == closed + weight * (blank - closed))
        panel = dict(
            ordinary_exit_changes={
                str(t): summarize(x, live) for t, x in exits.items()
            },
            opposite_sign=int(np.sum((exits[1] * exits[2] < 0) & live)),
            equal_changes=int(np.sum((exits[1] == exits[2]) & live)),
            single_change_relative_to_double_bearing=summarize(residual, live),
            double_relative_to_own_bearing=0,
            head_blind_change=0,
            matched_partial_common_read=summarize(common, live),
            matched_partial_fixed_stations={
                str(x): str(common[x]) for x in (25, 36, 104, 116) if x not in pins[:2]
            },
            partial_fraction_of_open_response=str(weight),
            single_to_double_exit_ratio=str(ratio),
        )
        sample = {}
        for x in [25, 36, 104, 116]:
            if x not in pins[:2]:
                sample[str(x)] = {str(t): str(exits[t][x]) for t in tables}
        panel["fixed_stations"] = sample
        controls = {}
        for key, row in matched["preparations"].items():
            d = np.array(list(map(Q, row["head_fractions"])), object) - mass["original"]
            assert np.all(d @ field[1] == d @ field[2])
            assert np.all(d @ field[2] == Q(row["blank_change"]) * (blank - closed))
            controls[key] = {
                str(x): str((d @ field[1])[x])
                for x in (25, 36, 104, 116)
                if x not in pins[:2]
            }
        panel["modular_and_energy_matched_controls"] = controls
        if (a, b) == (25, 36):
            records = {}
            for t, table in tables.items():
                bank = held(table, slots, pins, 466118 + t, [104])
                intervals = np.array(bank["launches"]["104"]["intervals"])
                truth = np.array(field[t][:, 104], float)
                assert np.all(intervals[:, 0] <= truth) and np.all(
                    truth <= intervals[:, 1]
                )
                mix = measured[:, None] * intervals
                band = [
                    float(mix.min(1).sum() - mass_err),
                    float(mix.max(1).sum() + mass_err),
                ]
                assert band[0] <= float(exits[t][104]) <= band[1]
                bank["measured_repair_change"] = band
                # At fixed enumerated enrollment only the held endpoint counts vary.
                mix = np.asarray(delta, float)[:, None] * intervals
                bank["change_at_exact_enrollment"] = [
                    float(mix.min(1).sum()),
                    float(mix.max(1).sum()),
                ]
                records[str(t)] = bank
            panel["held"] = records
            panel["matched_finite_arrivals"] = finite_arrivals(
                tables, pins, 104, matched_delta
            )
            # The held double provides the open-minus-closed transfer on this
            # same native chart. Unresolved trips enter both endpoint bounds.
            counts = records["2"]["launches"]["104"]
            total = records["2"]["trials_per_head"]
            edge = []
            for heads in ([0], list(range(1, 11))):
                N = len(heads) * total
                e = np.sqrt(np.log(4 / 0.003) / (2 * N))
                hits = sum(counts["b"][j] for j in heads)
                pending = sum(counts["pending"][j] for j in heads)
                edge.append([max(0, hits / N - e), min(1, (hits + pending) / N + e)])
            transfer = np.array([edge[0][0] - edge[1][1], edge[0][1] - edge[1][0]])
            outer = np.outer(transfer, effective_bounds)
            inferred = [float(outer.min()), float(outer.max())]
            assert all(inferred[0] <= exits[t][104] <= inferred[1] < 0 for t in (1, 2))
            panel["held_transfer"] = dict(
                open_minus_closed=transfer.tolist(),
                precursor_change_both_probes=inferred,
                combined_coverage_at_least=0.994,
            )
            outer = np.outer(transfer, matched["measured_changes"]["blank"])
            band = [float(outer.min()), float(outer.max())]
            assert band[0] <= float(common[104]) <= band[1] < 0
            panel["held_transfer"]["matched_common_change"] = band
        panels[f"{a}/{b}"] = panel
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "LANDMARK_ARRIVAL.md", HERE / "closed_passage.json")
        },
        scores_blank_and_half_partial=list(map(float, scores)),
        effective_precursor_bounds=list(map(float, effective_bounds)),
        all_beacon_pairs=survey,
        **landmark_calibration.measured_equivalence(accepted_weights),
        ordinary_energy_matched=matched,
        parked_profiles=parked_profiles(
            adj, slots, marks, eligible, mass["original"][0], matched_a
        ),
        shots_per_source=shots,
        head_counts={k: v.tolist() for k, v in acquired.items()},
        third_head_populations={k: list(map(str, v)) for k, v in mass.items()},
        panels=panels,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = HERE / "landmark_arrival.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    for key, panel in result["panels"].items():
        print(key, panel["equal_changes"], panel["opposite_sign"])
