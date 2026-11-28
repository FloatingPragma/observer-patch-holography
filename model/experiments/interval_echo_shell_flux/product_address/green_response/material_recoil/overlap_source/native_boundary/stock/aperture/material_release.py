"""Predict a remote material census from two shared contact-current readings."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
import json

import numpy as np

from . import moving_records as records

lab = records.lab
HERE = records.HERE
LAUNCHES = ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2))
N = 2**46


def matching_layers(n, blanks):
    layers = []
    for a, b in sorted(blanks):
        for partner in layers:
            if partner[a] == a and partner[b] == b:
                break
        else:
            partner = np.arange(n)
            layers.append(partner)
        partner[a], partner[b] = b, a
    return layers


def forward(values, layers):
    for partner in layers:
        values = (values + values[partner]) * Q(1, 2)
    return values


def ancestry(ends, layers):
    """Exact two-address backwards paths; a shared bond tosses just one coin."""
    counts = {tuple(ends): 1}
    for partner in reversed(layers):
        new = defaultdict(int)
        for (i, j), w in counts.items():
            p, q = int(partner[i]), int(partner[j])
            if p == j:
                new[i, j] += 2 * w
                new[p, q] += 2 * w
            else:
                for a in (i, p):
                    for b in (j, q):
                        new[a, b] += w
        counts = new
    result = defaultdict(Q)
    A, B = ends

    def kind(i):
        return 0 if i == A else 1 if i == B else 2

    for (i, j), w in counts.items():
        assert i != j
        result[kind(i), kind(j)] += Q(w, 4 ** len(layers))
    assert sum(result.values()) == 1
    return dict(result)


@lru_cache(maxsize=None)
def native(col, hold_material=False):
    half = lab.gates()[1]
    law = defaultdict(Q)
    for i in np.flatnonzero(half[:, col]):
        row = int(i)
        if hold_material:
            row = 9 * (row // 9) + col % 9
        law[row] += Q(int(half[i, col]), 8)
    return list(law.items())


def kernel_audit(dist, marks):
    for a in range(3):
        for b in range(3):
            law = defaultdict(Q)
            law[a, b] += Q(1, 2)
            law[b, a] += Q(1, 2)
            actual = {(k % 9 // 3, k % 3): p for k, p in native(3 * a + b)}
            assert actual == law and all(k < 9 for k, _ in native(3 * a + b))
    heads = [set() for _ in marks]
    for _, labels in dist:
        for group, part in zip(heads, labels, strict=True):
            group.update(part)
    for labels in heads[:2]:
        for m in labels:
            for ab in range(9):
                assert all(
                    k % 9 == ab and k // 9 in labels for k, _ in native(9 * m + ab)
                )


def background(pop, i, j, a, b):
    """Original A/B are parked without selecting on their values."""
    if i < 2 and j < 2:
        return [((a, b)[i], (a, b)[j], Q(1))]
    size = sum(pop)
    if i == j == 2:
        return [
            (u, v, Q(c * (d - int(u == v)), size * (size - 1)))
            for u, c in enumerate(pop)
            for v, d in enumerate(pop)
            if c and d > int(u == v)
        ]
    return [
        ((a, b)[i] if i < 2 else t, (a, b)[j] if j < 2 else t, Q(c, size))
        for t, c in enumerate(pop)
        if c
    ]


def acquire(pre, hold_material):
    post, jumps = defaultdict(Q), defaultdict(Q)
    for col, prob in pre.items():
        for row, w in native(col, hold_material):
            post[row] += w * prob
            jumps[row % 3 - col % 3] += w * prob
    assert sum(post.values()) == sum(jumps.values()) == 1
    return post, jumps


def marginals(post):
    values = np.full((2, 3), Q(0), object)
    for col, p in post.items():
        values[0, col % 9 // 3] += p
        values[1, col % 3] += p
    return values


def trial(dist, ends, layers, ancestors, launch, hold_material=False):
    A, B = ends
    n = len(layers[0])
    start = defaultdict(Q)
    resident = defaultdict(Q)
    stationary = np.full(3, Q(0), object)
    for (pop, labels), weight in dist.items():
        stationary += np.asarray(pop, object) * (weight / n)
        for m in labels[-1]:
            col = 9 * m + 3 * launch[0] + launch[1]
            prob = weight / len(labels[-1])
            start[col] += prob
            for row, w in native(col, hold_material):
                resident[pop, row] += w * prob
    pre2 = defaultdict(Q)
    for (pop, row), prob in resident.items():
        head, ab = divmod(row, 9)
        a, b = divmod(ab, 3)
        for (i, j), c in ancestors.items():
            for u, v, w in background(pop, i, j, a, b):
                pre2[9 * head + 3 * u + v] += prob * c * w
    post1, jump1 = acquire(start, hold_material)
    post2, jump2 = acquire(pre2, hold_material)
    initial = np.tile(stationary, (n, 1))
    initial[[A, B]] = np.eye(3, dtype=int)[list(launch)]
    # Independent computation using all three categorical terminal records.
    census = initial.copy()
    census[[A, B]] = marginals(post1)
    census = forward(census, layers)
    assert np.array_equal(census[[A, B]], marginals(pre2))
    census[[A, B]] = marginals(post2)
    census = forward(census, layers[::-1])
    J = np.array(
        [sum(k * p for k, p in jump.items()) for jump in (jump1, jump2)], object
    )
    dipole = np.full(n, Q(0), object)
    dipole[A], dipole[B] = -1, 1
    green = np.array(
        [forward(forward(dipole, layers), layers[::-1]), forward(dipole, layers[::-1])],
        object,
    )
    baseline = forward(forward(initial @ np.arange(3), layers), layers[::-1])
    assert np.array_equal(census @ np.arange(3), baseline + J @ green)
    assert all(census.sum(axis=1) == 1) and sum(J @ green) == 0
    return census, (jump1, jump2), J, green


def mix(arms, weights):
    out = defaultdict(Q)
    for name, weight in weights.items():
        for k, p in arms[name].items():
            out[k] += Q(weight) * p
    return dict(out)


def samples(laws, rng, scores, n=N):
    counts = []
    for law in laws:
        prob = np.array([float(law.get(int(s), 0)) for s in scores])
        counts.append(rng.multinomial(n, prob / prob.sum()))
    return np.array(counts)


def whole_board_check(board, marks, layers, dist, launch, prediction, rng):
    """Direct sharp runs, with the full board present through every appointment."""
    n = 8192
    batch = np.arange(n)
    head = np.zeros((n, len(marks)), dtype=np.int8)
    loads = np.zeros((n, board.n), dtype=np.int8)
    selected = rng.choice(len(dist), n, p=np.array(list(dist.values()), float))
    for i, (pop, labels) in enumerate(dist):
        ids = batch[selected == i]
        if not len(ids):
            continue
        # A uniform random permutation of the original material. Park both
        # endpoints *after* this draw, independently of their values.
        order = np.argsort(rng.random((len(ids), board.n)), axis=1)
        for t, sites in [
            (1, order[:, : pop[1]]),
            (2, order[:, pop[1] : pop[1] + pop[2]]),
        ]:
            loads[ids[:, None], sites] = t
        for j, label in enumerate(labels):
            head[ids, j] = rng.choice(label, len(ids))
    A, B = marks[-1]
    parked = loads[:, [A, B]].copy()
    initial_inventory = loads.sum(1) + sum(launch)
    loads[:, [A, B]] = launch
    half = lab.gates()[1]
    tickets = np.stack([np.repeat(np.arange(99), half[:, col]) for col in range(99)])
    local_counts = []
    for sequence in (layers, layers[::-1]):
        for j in (2, 0, 1):
            u, v = marks[j]
            col = 9 * head[:, j] + 3 * loads[:, u] + loads[:, v]
            row = tickets[col, rng.integers(8, size=n)]
            if j == 2:
                local_counts.append(np.bincount(row % 3 - loads[:, v] + 2, minlength=5))
            head[:, j] = row // 9
            loads[:, u], loads[:, v] = row % 9 // 3, row % 3
        for partner in sequence:
            u = np.flatnonzero(partner > np.arange(board.n))
            v = partner[u]
            # Read the full native ticket, without using the averaging formula.
            col = 3 * loads[:, u] + loads[:, v]
            row = tickets[col, rng.integers(8, size=col.shape)]
            assert np.all(row < 9)
            loads[:, u], loads[:, v] = row // 3, row % 3
    assert np.array_equal(loads.sum(1) + parked.sum(1), initial_inventory)
    counts = np.stack([(loads == t).sum(0) for t in range(3)], axis=1)
    census, jumps = prediction[:2]
    delta = np.max(abs(counts / n - np.asarray(census, float)))
    for obs, law in zip(local_counts, jumps, strict=True):
        expected = [float(law.get(t, 0)) for t in range(-2, 3)]
        delta = max(delta, np.max(abs(obs / n - expected)))
    # Two schedules, 240 three-bin reads and two five-bin local reads each.
    band = np.sqrt(np.log(2 * 2 * (board.n * 3 + 10) / 0.01) / (2 * n))
    assert delta < band
    return dict(
        sharp_runs=n,
        appointments_per_run=2 * len(board.seams),
        launch=launch,
        maximum_frequency_error=float(delta),
        simultaneous_99_halfwidth=float(band),
        local_jump_counts=[x.tolist() for x in local_counts],
        terminal_load_counts=counts.tolist(),
    )


def experiment(board, marks, hops, arms, weights, reverse):
    A, B = marks[-1]
    layers = matching_layers(board.n, set(board.seams) - set(marks))
    if reverse:
        layers.reverse()
    ancestors = ancestry((A, B), layers)
    radii = hops[219]
    contrasts = dict(
        matched=weights,
        repair=[{"source": "1"}, {"original": "1"}],
        ordinary_single=[{"single": "1"}, {"original": "1"}],
        ordinary_double=[{"ordinary": "1"}, {"original": "1"}],
    )
    results = {name: {} for name in contrasts}
    controls = {}
    for launch in LAUNCHES:
        _, _, J, _ = trial(arms["source"], (A, B), layers, ancestors, launch, True)
        controls[str(launch)] = list(map(str, J))
        assert all(J == 0)
    rng = np.random.default_rng(671239 + int(reverse))
    cache = {}
    for name, mixtures in contrasts.items():
        exact, local, held = [], [], []
        for arm in mixtures:
            key = tuple(sorted(arm.items()))
            if key not in cache:
                dist = mix(arms, arm)
                cache[key] = [
                    trial(dist, (A, B), layers, ancestors, launch)
                    for launch in LAUNCHES
                ]
            jrows, observed, terminal = [], [], []
            for census, jumps, J, green in cache[key]:
                jrows.append(J)
                if name == "matched":
                    freqs = samples(jumps, rng, range(-2, 3))
                    observed.append((freqs @ np.arange(-2, 3) / N).tolist())
                    held_counts = samples(
                        [dict(enumerate(x)) for x in census], rng, range(3)
                    )
                    terminal.append(held_counts @ np.arange(3) / N)
            exact.append(np.array(jrows, object))
            local.append(np.array(observed))
            held.append(np.array(terminal))
        diff = exact[0] - exact[1]
        for i, launch in enumerate(LAUNCHES[1:], 1):
            amount = sum(launch)
            dj = (diff[i] - diff[0]) / amount
            # Subtract the no-probe run within each preparation. Across the
            # two preparations the two baseline differences cancel exactly.
            x = dj @ green
            radial = x @ radii
            tag = ",".join(map(str, launch))
            report = dict(
                current_per_unit=list(map(str, dj)),
                mean_radius_change=str(radial),
                aperture_load=[str(x[radii < r].sum()) for r in records.RADII],
                stations_changed=sum(p != 0 for p in x),
            )
            if name == "matched":
                loc = local[0] - local[1]
                forecast = (loc[i] - loc[0]) @ np.asarray(green @ radii, float) / amount
                term = (held[0][i] - held[0][0] - held[1][i] + held[1][0]) / amount
                # The total contrast is exactly zero; stratify by terminal
                # station and subtract a common constant from the radius score.
                centered = radii - radii.mean()
                arrived = float(term @ centered)
                # Two schedules, four releases and two banks. Each contrast
                # uses four independent preparations. Local jump range is 4,
                # terminal load range is 2. Union bound over the 16 intervals.
                factor = np.sqrt(2 * np.log(2 * 16 / 0.01) / N) / amount
                wf = float(
                    4 * np.linalg.norm(np.asarray(green @ radii, float)) * factor
                )
                wt = float(2 * np.linalg.norm(centered) * factor)
                assert (
                    abs(float(radial) - forecast) <= wf
                    and abs(float(radial) - arrived) <= wt
                )
                report.update(
                    local_radius_estimate=float(forecast),
                    local_99_halfwidth=wf,
                    terminal_radius_estimate=arrived,
                    terminal_99_halfwidth=wt,
                )
            results[name][tag] = report
        for amount in (1, 2):
            left = diff[LAUNCHES.index((amount, 0))] - diff[0]
            right = diff[LAUNCHES.index((0, amount))] - diff[0]
            assert all(left + right == 0)
    first = tuple(sorted(weights[0].items()))
    index = LAUNCHES.index((2, 0))
    direct = whole_board_check(
        board,
        marks,
        layers,
        mix(arms, weights[0]),
        LAUNCHES[index],
        cache[first][index],
        np.random.default_rng(665713 + int(reverse)),
    )
    return dict(
        sharp_board_replay=direct,
        held_material_control_currents=controls,
        layer_sizes=[int((p != np.arange(board.n)).sum() // 2) for p in layers],
        ancestry={str(k): str(v) for k, v in ancestors.items()},
        contrasts=results,
    )


def equilibrium_release(dist, n):
    """Late material readings after the measured, parked intervention."""
    panels = []
    for launch in LAUNCHES:
        marginal = np.full(3, Q(0), object)
        for (pop, _), w in dist.items():
            for old_a, old_b, prob in background(pop, 2, 2, 0, 0):
                counts = list(pop)
                for t in (old_a, old_b):
                    counts[t] -= 1
                for t in launch:
                    counts[t] += 1
                assert min(counts) >= 0 and sum(counts) == n
                marginal += np.asarray(counts, object) * (w * prob / n)
        panels.append(marginal)
    shifts = {}
    for launch, law in zip(LAUNCHES[1:], panels[1:], strict=True):
        t = sum(launch)
        exact = np.eye(3, dtype=int)[t] - np.eye(3, dtype=int)[0]
        assert np.array_equal(n * (law - panels[0]), exact)
        shifts[str(launch)] = list(map(str, law - panels[0]))
    return dict(
        per_station_empty_launch=list(map(str, panels[0])),
        per_station_release_changes=shifts,
    )


def run():
    _, words, cells, *_ = lab.energy.pair.setup()
    board, _, marks, *_ = lab.energy.returned(words)
    _, hops, _ = records.thermal_radar.chart(board, marks)
    center = 219
    uniform = np.full(len(words), Q(1, len(words)), object)
    arms = dict(
        original=records.ensemble(words, uniform),
        source=records.ensemble(words, lab.energy.least_change(cells) / len(words)),
        single=records.ensemble(words, uniform, single=1),
        ordinary=records.ensemble(words, uniform, extra=1),
    )
    A, B = marks[-1]
    assert hops[center, A] == hops[center, B] == 7
    # Blanks deliver a perfect message in one full pulse, the first two
    # written contacts cannot signal at all with either load alphabet.
    # Lower even the third contact's cost to zero: this bounds every
    # possible nonnegative time assignment on that contact.
    zero_cost = np.minimum.reduce(
        [hops[center], hops[center, A] + hops[B], hops[center, B] + hops[A]]
    )
    assert np.array_equal(zero_cost, hops[center])
    for dist in arms.values():
        kernel_audit(dist, marks)
        signals = records.summary(lab.before(dist))["signals"][:, :2]
        assert np.all(signals[0] == 1) and np.all(signals[1:3] == 0)
        assert np.all(signals[3] > 0)
    frozen = json.loads((HERE / "moving_records.json").read_bytes())
    weights = frozen["branch_balance"]["first_ball_enclosing_third"][
        "history_and_volume_matched"
    ]["arm_weights"]
    # Each scheduled local kernel is doubly stochastic and can hold.
    # A sweep can enact any one allowed transition while all other visits
    # hold. Connected blank swaps permute every fixed material inventory.
    # Consequently every recurrent class is uniform over its possible
    # material permutations, also when the release changes the head class.
    assert hops.max() < board.n and np.all(lab.gates()[1].diagonal() > 0)
    return dict(
        equilibrium={
            name: equilibrium_release(dist, board.n) for name, dist in arms.items()
        },
        preparation_sha256=sha256(
            (HERE / "moving_records.json").read_bytes()
        ).hexdigest(),
        protocol_sha256=sha256((HERE / "MATERIAL_RELEASE.md").read_bytes()).hexdigest(),
        marks=marks,
        station=center,
        matched_mixtures=weights,
        radar_endpoint_distances=hops[center, [A, B]].tolist(),
        stations_invariant_under_zero_cost_shortcut=int(
            np.sum(zero_cost == hops[center])
        ),
        aperture_radii=list(records.RADII),
        ideal_draws_per_histogram=N,
        schedules={
            order: experiment(board, marks, hops, arms, weights, rev)
            for order, rev in [("forward_reverse", False), ("reverse_forward", True)]
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = HERE / "material_release.json"
    if args.write:
        path.write_text(receipt)
    else:
        assert path.read_text() == receipt
    print(
        "Both schedules and all 240 terminal readings reproduce the signed-current prediction."
    )
