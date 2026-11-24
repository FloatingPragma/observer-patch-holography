"""Remote polarization of five native electrodes around the joined three-wall source."""

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.stats import binom

from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from ...row_13_walls.energy import is_wall
from ...row_20_electric import experiment as electric
from . import belt, upper

SHOTS = 2**22
SEED = 121314


def geometry():
    board = carrier.federation(1)
    gateways = {
        (u // 12, v // 12): u for edge in board.external for u, v in (edge, edge[::-1])
    }
    network = nx.Graph(gateways.keys())
    axes = belt.belts()

    def charts(owner, adjoining, indices, outside_poles):
        for poles, order in axes.items():
            if sum(board.degree[owner * 12 + p] == 6 for p in poles) != outside_poles:
                continue
            for direction, origin in product((1, -1), range(10)):
                ports = [
                    12 * owner + order[(origin + direction * j) % 10] for j in range(10)
                ]
                if all(
                    ports[t] == gateways[owner, o]
                    for t, o in zip(indices, adjoining, strict=True)
                ):
                    yield (
                        ports,
                        [
                            owner * 12 + p
                            for p in poles
                            if board.degree[12 * owner + p] == 6
                        ],
                    )

    for middle in sorted(network):
        for left, right in combinations(sorted(network[middle]), 2):
            if network.has_edge(left, right):
                continue
            menu = [
                list(charts(*args))
                for args in (
                    (left, [middle], [4], 2),
                    (middle, [left, right], [1, 5], 1),
                    (right, [middle], [0], 2),
                )
            ]
            if not all(menu):
                continue
            layout = [m[0] for m in menu]
            ports = [p for order, _ in layout for p in order]
            active = {
                tuple(sorted((ports.index(a), ports.index(b))))
                for a, b in board.seams
                if a in ports and b in ports
            }
            inside = {
                (10 * i + a, 10 * i + b)
                for i in range(3)
                for a, b in belt.ring(10).seams
            }
            bridges = sorted(active - inside)
            assert inside <= active and bridges == [(4, 11), (15, 20)]
            assert all(
                board.degree[p] - sum(i in e for e in active) == 1
                for i, p in enumerate(ports)
            )
            return board, dict(
                owners=[left, middle, right],
                ports=ports,
                poles=[p for _, p in layout],
                bridges=bridges,
            )
    raise ValueError("No equivalent three-belt geometry with five exterior poles")


def routes(board, layout, dimers=6):
    graph = nx.Graph(board.seams)
    owners = set(layout["owners"])
    distances = nx.multi_source_dijkstra_path_length(
        graph, layout["ports"], weight=None
    )
    result = []
    for i, group in enumerate(layout["poles"]):
        for pole in group:
            start = next(p for p in board.neighbours[pole] if p // 12 != pole // 12)
            best = None

            def extend(path, blocked):
                nonlocal best
                if len(path) == 2 * dimers:
                    score = (
                        distances[path[-1]],
                        len({p // 12 for p in path}),
                        tuple(-p for p in path),
                    )
                    if best is None or score > best[0]:
                        best = (score, path)
                    return
                for a in board.neighbours[path[-1]]:
                    if a // 12 in owners or a in blocked:
                        continue
                    for b in board.neighbours[a]:
                        if (
                            b // 12 in owners
                            or b in blocked
                            or board.degree[a] != board.degree[b]
                        ):
                            continue
                        guards = set()
                        if len(path) == 2:
                            guards = set(board.neighbours[start]) - {pole, a}
                        if b not in guards:
                            extend(path + [a, b], blocked | {a, b} | guards)

            extend([pole, start], {pole, start})
            path = best[1]
            references = sorted(set(board.neighbours[start]) - {pole, path[2]})
            assert len(references) == 4 and not set(path) & set(references)
            sides = [
                j
                for j, p in enumerate(layout["ports"][10 * i : 10 * i + 10])
                if p in board.neighbours[pole]
            ]
            assert len(sides) == 5
            result.append(
                dict(
                    belt=i,
                    pole=pole,
                    ports=path,
                    references=references,
                    observed_sites=sides,
                    last_distances=[distances[path[2 * n - 2]] for n in (4, 6)],
                    receiver_distances=[distances[path[2 * n - 4]] for n in (4, 6)],
                )
            )
    return result


def read(n, source, kappa, times, broken=False, ordered=True):
    """Read one distant orientation with either detector preparation."""
    j, h = electric.couplings(kappa)
    word = np.arange(1 << n)
    z = 1 - 2 * ((word[:, None] >> np.arange(n)) & 1)
    bonds = z[:, :-1] * z[:, 1:]
    if broken:
        bonds[:, (n - 1) // 2] = 0
    generator = np.diag(j * (bonds.sum(axis=1) - source * z[:, 0]))
    for i in range(n):
        generator[word, word ^ (1 << i)] -= h
    e, v = np.linalg.eigh(generator)
    phase = np.exp(-1j * (e - e[0])[:, None] * np.asarray(times))
    if ordered:
        alternating = sum(1 << i for i in range(1, n, 2))
        prob = sum(
            abs(v @ (v[i, :, None] * phase)) ** 2 / 2
            for i in (alternating, len(word) - 1 - alternating)
        )
        return z[:, -2] @ prob
    amplitudes = v @ ((v.sum(axis=0) / np.sqrt(len(word)))[:, None] * phase)
    return z[:, -1] @ (abs(amplitudes) ** 2)


def instrument(n, kappa, rng):
    times = np.arange(129) / 8
    unit = read(n, 1, kappa, times)
    idx = int(np.argmax(abs(unit)))
    t = float(times[idx])
    menu = {s: float(read(n, s, kappa, [t])[0]) for s in range(-5, 6)}
    assert abs(menu[0]) < 1e-11
    assert all(abs(menu[s] + menu[-s]) < 1e-11 for s in menu)
    calibrated = {}
    for s in (-1, 0, 1, 2):
        plus = int(rng.binomial(SHOTS, (1 + menu[s]) / 2))
        calibrated[str(s)] = dict(
            plus=plus, minus=SHOTS - plus, mean=2 * plus / SHOTS - 1
        )
    broken = float(read(n, 1, kappa, [t], True)[0])
    assert abs(broken) < 1e-11
    blind_repeats = 129
    port_error = float(
        binom.cdf(blind_repeats // 2, blind_repeats, 0.5 + abs(menu[1]) / 4)
    )
    last = read(n, 1, kappa, times, ordered=False)
    return dict(
        shot_repeats_per_bit=blind_repeats,
        bit_misread=port_error,
        all_visible_bit_error_upper=min(1.0, 25 * port_error),
        independent_x_input_last_dimer_peak=float(np.max(abs(last))),
        dimers=n,
        kappa=kappa,
        time=t,
        signal=menu,
        calibration=calibrated,
        broken=broken,
        broken_plus=int(rng.binomial(SHOTS, (1 + broken) / 2)),
    )


def chart(word):
    if not is_wall(word):
        return None
    occupied = tuple(i for i in range(10) if (word >> i) & 1)
    return belt.wall(occupied, 10)


def shoulders(word, radius):
    _, _, centre = chart(word)
    windows = np.flatnonzero(
        abs(belt.signed_distance(np.arange(10), centre, 10)) <= radius
    )
    weights = np.zeros(10, dtype=int)
    for u in windows:
        for delta in (-1, 0, 1):
            weights[(u + delta) % 10] += 1
    return weights, len(windows)


def exterior_classes():
    wall = [w for w in range(1024) if is_wall(w)]
    pentagon = sum(1 << i for i in range(0, 10, 2))
    result = {}
    for loads in (3, 4):
        pairs = []
        for word in wall:
            if word.bit_count() != loads:
                continue
            copies = [
                v
                for v in range(1024)
                if v.bit_count() == loads and v & pentagon == word & pentagon
            ]
            others = [v for v in copies if not is_wall(v)]
            assert others
            pairs.append(
                dict(
                    wall=word,
                    other=min(others),
                    same_view=len(copies),
                    same_view_walls=len(copies) - len(others),
                )
            )
        compact = sum(1 << i for i in range(loads))

        def rotate(w, i):
            return ((w << i) | (w >> (10 - i))) & 1023

        def views(source, segmented):
            patterns = [rotate(source, i) & pentagon for i in range(10)]
            return Counter(patterns if segmented else [v.bit_count() for v in patterns])

        full = views(pairs[0]["wall"], False)
        assert full == views(compact, False)
        w, c = (views(s, True) for s in (pairs[0]["wall"], compact))
        tv = sum(abs(w[k] - c[k]) for k in set(w) | set(c)) / 20
        result[str(loads)] = dict(
            indistinguishable=pairs,
            whole_pole_counts=dict(full),
            segmented_single_pole_variation=tv,
        )
    return result


def joint(words, weights, panel, rng):
    """A uniformly chosen window-covered port and a one-in-three independent reference."""
    unit = panel["signal"][1]
    seen = np.array([bool(chart(w)) for w in range(1024)])
    pieces = np.stack([(words >> (i * 10)) & 1023 for i in range(3)], axis=1)
    accepted = seen[pieces].all(axis=1)
    records = []
    # Each source shot addresses a belt port uniformly and compares it with a
    # random bulk phase. Use the terminal wall centre only in the analysis.
    for side in (0, 2):
        for radius in (2, 3):
            # Terminal read groups copies by wall type and number of windows.
            # Counts keep trials outside the aperture and all failures.
            bins = {}
            for w in np.unique(pieces[accepted, side]):
                at = accepted & (pieces[:, side] == w)
                p = float(weights[at].sum())
                loads = int(w).bit_count()
                coverage, nwin = shoulders(int(w), radius)
                for site, multiplicity in enumerate(coverage):
                    sign = (int(w) >> site & 1) - 1 / 3
                    label = (loads, nwin, int(multiplicity))
                    mass, contrast = bins.get(label, (0.0, 0.0))
                    bins[label] = (mass + p / 10, contrast + p * unit * sign / 10)
            labels = sorted(bins)
            probs = [
                v
                for key in labels
                for v in (
                    (bins[key][0] + bins[key][1]) / 2,
                    (bins[key][0] - bins[key][1]) / 2,
                )
            ]
            outside = float(weights[~accepted].sum())
            # Resolve the remote bright/dark readings for non-three-wall outcomes too.
            mean_bits = (
                np.array([int(w).bit_count() for w in pieces[~accepted, side]]) / 10
            )
            outside_contrast = unit * float(weights[~accepted] @ (mean_bits - 1 / 3))
            draws = rng.multinomial(
                SHOTS,
                np.clip(
                    probs
                    + [
                        (outside + outside_contrast) / 2,
                        (outside - outside_contrast) / 2,
                    ],
                    0,
                    1,
                ),
            )
            doses = {int(load_class): [0, 0] for load_class, _, _ in labels}
            for (load_class, _, m), plus, minus in zip(
                labels, draws[:-2:2], draws[1:-2:2], strict=True
            ):
                doses[load_class][0] += m * int(plus - minus)
                doses[load_class][1] += int(plus + minus)
            gain = (
                panel["calibration"]["1"]["mean"] - panel["calibration"]["-1"]["mean"]
            ) / 2
            estimates = {
                load_class: 10 * diff / (3 * count * gain)
                for load_class, (diff, count) in doses.items()
            }
            records.append(
                dict(
                    side=side,
                    radius=radius,
                    labels=labels,
                    plus_minus=list(map(int, draws)),
                    by_load_charge=estimates,
                )
            )
    return dict(three_wall_probability=float(weights[accepted].sum()), samples=records)


def pooled(sectors, panels):
    """Conditional on the observed group sizes, every covered wall has the same target mean."""
    result = {}
    # Two gain quadratures plus two wall classes for each of two instruments.
    alpha = 0.01 / 8
    for key, panel in panels.items():
        gain = (
            panel["calibration"]["1"]["mean"] - panel["calibration"]["-1"]["mean"]
        ) / 2
        gain_error = np.sqrt(2 * np.log(2 / alpha) / SHOTS)
        data = {}
        for loads in (3, 4):
            total = 0
            numerator = 0
            for sector in sectors:
                for sample in sector["receivers"][key]["samples"]:
                    counts = sample["plus_minus"]
                    for i, (label, _, multiplicity) in enumerate(sample["labels"]):
                        if label == loads:
                            plus, minus = counts[2 * i : 2 * i + 2]
                            total += plus + minus
                            numerator += multiplicity * (plus - minus)
            mean = numerator / total
            # Single-shot weighted orientation lies in [-3,3].
            margin = 3 * np.sqrt(2 * np.log(2 / alpha) / total)
            interval = [
                10 * x / (3 * y)
                for x, y in product(
                    (mean - margin, mean + margin),
                    (gain - gain_error, gain + gain_error),
                )
            ]
            data[str(loads)] = dict(
                accepted_trials=total,
                estimated_charge=10 * mean / (3 * gain),
                simultaneous_99_interval=[min(interval), max(interval)],
            )
        result[key] = data
    return result


def run():
    board, layout = geometry()
    sensors = routes(board, layout)
    rng = np.random.default_rng(SEED)
    panels = {str(n): instrument(n, candidate.KAPPA, rng) for n in (4, 6)}
    controls = {str(n): instrument(n, 0.0, rng) for n in (4, 6)}
    sectors = []
    for n in (10, 11):
        (words, _, _, _), (energy, amplitudes, _) = upper.trial_state(
            upper.bins(), layout["bridges"], n
        )
        p = abs(amplitudes) ** 2
        sectors.append(
            dict(
                loads=n,
                source_energy=energy,
                receivers={
                    key: joint(words, p, panel, rng) for key, panel in panels.items()
                },
            )
        )
    return currents._portable(
        dict(
            shots=SHOTS,
            seed=SEED,
            geometry=layout,
            routes=sensors,
            registered=panels,
            exchange_only=controls,
            middle_ambiguity=exterior_classes(),
            sectors=sectors,
            pooled=pooled(sectors, panels),
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as f:
        f.write(json.dumps(run(), indent=2, allow_nan=False) + "\n")
