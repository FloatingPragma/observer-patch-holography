"""Native coincidence support for two initially filled, charge-neutral receiving owners."""

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import ceil, factorial
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.openwave.row_01_charge.measure import (
    sharp_contrast,
)
from model.executable_model import federation

RECEIPT = Path(__file__).with_name("results.json")


def local_certificate():
    """Integer support on all 99 columns, separately for exchange and contrast."""
    a = np.zeros((90, 9), dtype=np.int64)
    for i, (left, right) in enumerate(law.SRC):
        for coin in (0, 1):
            record = law.LABELS.index((coin, left - right))
            repaired = law.native_event(left, right, coin)
            a[9 * record + law.SRC.index(repaired), i] = 1
    b = 2 * np.eye(90, dtype=np.int64) - a @ a.T
    magnitude = [-1] + [abs(d) for _, d in law.LABELS]
    flip = [law.LABELS.index((c, -d)) for c, d in law.LABELS]
    twice_g = (
        np.eye(10, dtype=np.int64) - np.eye(10, dtype=np.int64)[flip],
        2 * np.diag(np.square(magnitude[1:])),
    )
    histogram = np.array(
        [tuple(pair.count(i) for i in range(3)) for pair in law.SRC] * 11
    )
    changes_inventory = np.any(histogram[:, None] != histogram[None, :], axis=2)
    labels = np.repeat(magnitude, 9)
    changes_collar = labels[:, None] != labels[None, :]
    blocks = []
    for g, actual in zip(twice_g, law.local_term_parts(), strict=True):
        g = np.kron(g, np.eye(9, dtype=np.int64))
        assert not np.any(a.T @ g @ b)
        octuple = np.zeros((99, 99), dtype=np.int64)
        octuple[:9, :9] = 2 * a.T @ g @ a
        octuple[9:, 9:] = b @ g @ b
        np.testing.assert_allclose(octuple / 8, actual, atol=1e-13, rtol=0)
        blocks.append(
            dict(
                nonzero_entries=int(np.count_nonzero(octuple)),
                inventory_changing=int(np.count_nonzero(octuple[changes_inventory])),
                collar_class_changing=int(np.count_nonzero(octuple[changes_collar])),
            )
        )
        assert not np.any(octuple[changes_inventory | changes_collar])
    return dict(
        dimension=99,
        exchange=blocks[0],
        contrast=blocks[1],
        terminal_port_minus_zero=[sharp_contrast(i, 0) for i in range(3)],
    )


def fixture(board):
    """First two external gateways of owner zero; one distinct inner address each."""
    gates = [(a, b) for a, b in board.external if a // 12 == 0]
    used = set(a for a, _ in gates)
    arms = []
    for g, r in gates[:2]:
        inner = next(i for i in board.neighbours[g] if i < 12 and i not in used)
        outer = next(i for i in board.neighbours[r] if i // 12 == r // 12)
        used.add(inner)
        arms.append((inner, g, r, outer))
    assert len({r // 12 for _, _, r, _ in arms}) == 2
    spectator = next(i for i in range(12) if i not in used)
    return arms, spectator


def swaps(state, board):
    """Generate each different-load edge swap once, with native amplitude -1/2."""
    holes, doubles = state
    for i in holes:
        for j in board.neighbours[i]:
            if j in holes:
                continue
            dst_holes = tuple(sorted(j if p == i else p for p in holes))
            if j in doubles:
                yield dst_holes, tuple(sorted(i if p == j else p for p in doubles))
            else:
                yield dst_holes, doubles
    for i in doubles:
        for j in board.neighbours[i]:
            if j not in holes and j not in doubles:
                yield holes, tuple(sorted(j if p == i else p for p in doubles))


def shortest_layers(start, board):
    layers = [{start: 1}]
    seen = {start}
    for _ in range(4):
        layer = defaultdict(int)
        for state, ways in layers[-1].items():
            for reached in swaps(state, board):
                if reached not in seen:
                    layer[reached] += ways
        seen.update(layer)
        layers.append(layer)
    return layers


def full_board(level):
    board = federation(level)
    arms, spectator = fixture(board)
    initial = (
        tuple(sorted([spectator, *(i for i, _, _, _ in arms)])),
        tuple(sorted(g for _, g, _, _ in arms)),
    )
    terminal = (
        tuple(sorted([spectator, *(r for _, _, r, _ in arms)])),
        tuple(sorted(s for _, _, _, s in arms)),
    )
    source = shortest_layers(initial, board)
    target = shortest_layers(terminal, board)
    for length in range(9):
        left = source[length // 2]
        right = target[length - length // 2]
        paths = sum(left[s] * right[s] for s in left.keys() & right.keys())
        if paths:
            break
    else:
        raise AssertionError("The exhibited two separate four-swap paths were lost")
    coefficient = Fraction(paths, 2**length * factorial(length))
    return dict(
        cells=board.cells,
        arms=arms,
        remaining_vacancy=spectator,
        source_holes_doubles=initial,
        observed_holes_doubles=terminal,
        first_nonzero_amplitude_order=length,
        shortest_paths=paths,
        leading_full_word_probability_coefficient=str(coefficient**2),
    )


def arm_response(board, ports, conversion_strength, cut):
    """Full four-port alphabet, including every possible contrast word at arrival."""
    words = list(product(range(3), repeat=4))
    word_index = {v: i for i, v in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    local = law.seam_term(8)
    pair = law.SRC.index((1, 1))
    for state in ((0, 2), (2, 0)):
        other = law.SRC.index(state)
        local[other, pair] += conversion_strength / np.sqrt(2)
        local[pair, other] += conversion_strength / np.sqrt(2)
    for col, word in enumerate(words):
        h[col, col] += 8 * sum(
            (6 - board.degree[p]) * n**2 for p, n in zip(ports, word)
        )
        for edge in range(3):
            if cut and edge == 1:
                continue
            src = law.SRC.index(word[edge : edge + 2])
            for dst in np.flatnonzero(local[:, src]):
                target = word[:edge] + law.SRC[dst] + word[edge + 2 :]
                h[word_index[target], col] += local[dst, src]
    _, components = connected_components(h != 0, directed=False)
    times = [0.25, 1, 4, 16, 64]
    responses = {}
    for seeded, word in ((False, (1, 1, 1, 1)), (True, (0, 2, 1, 1))):
        ix = np.flatnonzero(components == components[word_index[word]])
        tag = [j for j, i in enumerate(ix) if words[i][2:] in ((0, 2), (2, 0))]
        e, v = eigh(h[np.ix_(ix, ix)])
        amps = (
            np.exp(-1j * np.outer(times, e)) * v[ix.tolist().index(word_index[word]), :]
        ) @ v.T
        np.testing.assert_allclose(np.sum(abs(amps) ** 2, axis=1), 1, atol=1e-12)
        responses[seeded] = np.sum(abs(amps[:, tag]) ** 2, axis=1)
    return times, responses


def bench():
    """Two four-port source/receiver routes, with the surplus vacancy held at rest."""
    board = federation(0)
    arms, spectator = fixture(board)
    rows = []
    shots = 2**20
    rng = np.random.default_rng(113024)
    for strength, cut in ((0, False), (0, True), (8, False)):
        times, left = arm_response(board, arms[0], strength, cut)
        _, right = arm_response(board, arms[1], strength, cut)
        for seeds in ((False, False), (True, False), (True, True)):
            panels = []
            for t, left_prob, right_prob in zip(
                times, left[seeds[0]], right[seeds[1]], strict=True
            ):
                prob = [
                    (1 - left_prob) * (1 - right_prob),
                    left_prob * (1 - right_prob),
                    (1 - left_prob) * right_prob,
                    left_prob * right_prob,
                ]
                # All accepted and rejected records; disjoint routes multiply on this bench only.
                count = rng.multinomial(shots, prob)
                panels.append(
                    dict(
                        time=t,
                        none_left_right_both_counts=count.tolist(),
                        both_probability=float(f"{left_prob * right_prob:.8g}"),
                    )
                )
            rows.append(
                dict(
                    conversion_strength=strength,
                    bridge_cut=cut,
                    seeded_neutral_pairs=sum(seeds),
                    panels=panels,
                )
            )
    return dict(
        arms=arms, fixed_vacancy=spectator, shots_per_setting=shots, settings=rows
    )


def formation_receiver_bounds():
    """Apply the blank-host spectral envelope to two neutral, exhaustive terminal patches."""
    path = RECEIPT.parent.parent / "excited_bounds.json"
    spectrum = json.loads(path.read_text())
    certificate = spectrum["any_three_gateways"]
    lower = spectrum["verified_negative_sector_floor"]
    slope, tilt = map(Fraction, certificate["per_defect_and_charge_coefficients"])
    rows = {}
    for ratio, reported in certificate["necessary_dressing_sizes"].items():
        # E_parent >= ratio * E_light >= ratio * lower; Q=-1, D=2*p+1.
        needed = lower * Fraction(ratio)
        pairs = ceil((needed + tilt - slope) / (2 * slope))
        assert slope * (2 * pairs - 1) - tilt < needed <= slope * (2 * pairs + 1) - tilt
        assert 2 * pairs + 1 == reported["minimum_defects_in_fixed_composition"]
        # One-vacancy daughter has no doubles. All other doubles enter the
        # two receivers exactly when the remaining material returns to its host.
        larger = (pairs + 1) // 2
        rows[ratio] = dict(
            minimum_defect_ports=2 * pairs + 1,
            minimum_pairs_in_two_exhaustive_receivers=pairs,
            minimum_pairs_in_one_of_the_two_receivers=larger,
            minimum_ports_in_that_neutral_receiver=2 * larger,
            minimum_doubles_outside_two_single_cell_receivers=max(0, pairs - 12),
        )
    return rows


def run():
    arrays = []
    for state in product(range(3), repeat=4):
        # Independent receiver check on two two-port output words.
        if all(sum(pair) == 2 and pair != (1, 1) for pair in (state[:2], state[2:])):
            arrays.append(list(state))
    assert len(arrays) == 4 and all(w.count(2) == 2 and w.count(0) == 2 for w in arrays)
    return dict(
        local=local_certificate(),
        conditional_one_vacancy_daughter=formation_receiver_bounds(),
        two_dimer_neutral_records=arrays,
        all_open_sharp_preparations=[full_board(i) for i in range(3)],
        two_arm_bench=bench(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2) + "\n"
    if args.write:
        RECEIPT.write_text(encoded)
    else:
        stored = json.loads(RECEIPT.read_text())
        current = json.loads(encoded)
        for actual, expected in zip(
            current["two_arm_bench"]["settings"],
            stored["two_arm_bench"]["settings"],
            strict=True,
        ):
            for a, e in zip(actual["panels"], expected["panels"], strict=True):
                np.testing.assert_allclose(
                    a.pop("both_probability"),
                    e.pop("both_probability"),
                    atol=0,
                    rtol=2e-5,
                )
        if current != stored:
            raise SystemExit(
                "Neutral-ejecta certificate or observer counts differ from stored receipt"
            )
        print(
            "Reproduced the native selection rule, three full-board witnesses and changed-law control."
        )
