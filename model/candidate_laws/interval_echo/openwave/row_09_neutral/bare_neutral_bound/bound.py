"""An exact blank-exterior neutral formation-energy test on physical cells."""

import argparse
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from ......executable_model import cell, federation
from .....contrast_memory_spherical_v1 import law
from ..passive_decay.measure import compare_records


HERE = Path(__file__).resolve().parent
KAPPA = 8
LABELS = 32768
PERIOD = 512


def microscopic_check():
    """Independently rebuild the blank restriction from the full native write."""
    source = [(a, b) for a in range(3) for b in range(3)]
    assert source == list(law.SRC)
    write = np.zeros((90, 9), dtype=int)
    flip = np.zeros((90, 90), dtype=int)
    cost = np.zeros((90, 90), dtype=int)
    for j, (a, b) in enumerate(source):
        s = a + b
        for coin in (0, 1):
            y = s // 2 + coin * (s % 2), s // 2 + (1 - coin) * (s % 2)
            write[9 * law.LABELS.index((coin, a - b)) + source.index(y), j] = 1
    for j, (coin, contrast) in enumerate(law.LABELS):
        for i in range(9):
            flip[9 * j + i, 9 * law.LABELS.index((coin, -contrast)) + i] = 1
            cost[9 * j + i, 9 * j + i] = contrast**2
    twice = np.eye(90, dtype=int) - flip + 2 * KAPPA * cost
    np.testing.assert_array_equal(write.T @ write, 2 * np.eye(9, dtype=int))
    np.testing.assert_array_equal(
        (2 * np.eye(90, dtype=int) - write @ write.T) @ twice @ write, 0
    )
    exact = write.T @ twice @ write
    expected = 2 * np.eye(9, dtype=int)
    for j, (a, b) in enumerate(source):
        expected[j, j] += 4 * KAPPA * (a - b) ** 2
        expected[source.index((b, a)), j] -= 2
    np.testing.assert_array_equal(exact, expected)
    exchange, contrast = law.local_term_parts()
    np.testing.assert_allclose(
        4 * (exchange + KAPPA * contrast)[:9, :9], expected, atol=3e-14, rtol=0
    )
    np.testing.assert_allclose(
        (exchange + KAPPA * contrast)[9:, :9], 0, atol=1e-14, rtol=0
    )
    return exact.tolist()


def census(shape):
    words = (np.arange(3**12, dtype=int)[:, None] // (3 ** np.arange(12))) % 3 - 1
    p = np.sum(words == 1, axis=1)
    m = np.sum(words == -1, axis=1)
    q = p - m
    internal = sum((words[:, a] - words[:, b]) ** 2 for a, b in shape.seams)
    corrections = words**2 + 2 * words
    private = np.sort(corrections, axis=1)[:, :9].sum(axis=1)
    score = 2 * (internal + private) - 3 * q
    sectors = {}
    for plus in range(13):
        for minus in range(13 - plus):
            mask = (p == plus) & (m == minus)
            sectors[plus, minus] = int(min(score[mask]))
    assert len(sectors) == 91 and sectors[0, 0] == 0
    assert min(v for k, v in sectors.items() if k != (0, 0)) == 7
    # For both signs, verify the equality conditions that force at least
    # one added unit of contact contrast on an intact whole federation.
    allowed = sorted(combinations(range(12), 3))
    equality = {}
    for signed, total in ((-1, 11), (0, 18), (1, 7)):
        selected = np.flatnonzero((q == signed) & (p + m > 0) & (score == total))
        assert len(selected)
        cases = 0
        for j in selected:
            assert (int(p[j]), int(m[j])) == {-1: (0, 1), 0: (1, 1), 1: (1, 0)}[signed]
            for outside in allowed:
                if (
                    int(corrections[j].sum() - sum(corrections[j, a] for a in outside))
                    != private[j]
                ):
                    continue
                cases += 1
                assert all(
                    (value == 1) == (a in outside)
                    for a, value in enumerate(words[j])
                    if value
                )
        equality[signed] = cases
    return sectors, equality


def subset_check(shape, expected):
    neighbors = [sum(1 << j for j in shape.neighbours[i]) for i in range(12)]
    universe = (1 << 12) - 1
    sets = [[j for j in range(12) if word >> j & 1] for word in range(universe + 1)]
    cuts = [
        sum((neighbors[j] & (universe ^ s)).bit_count() for j in bitset)
        for s, bitset in enumerate(sets)
    ]
    observed = dict.fromkeys(expected, 10000)
    visited = 0
    for positive, ports in enumerate(sets):
        vacant = complement = universe ^ positive
        plus = len(ports)
        while True:
            minus = vacant.bit_count()
            charge = plus - minus
            crossing = sum((neighbors[j] & vacant).bit_count() for j in ports)
            private = 3 * max(0, plus - 3) - min(9, minus)
            score = (
                2 * (cuts[positive] + cuts[vacant] + 2 * crossing + private)
                - 3 * charge
            )
            observed[plus, minus] = min(observed[plus, minus], score)
            visited += 1
            if not vacant:
                break
            vacant = (vacant - 1) & complement
    assert visited == 3**12 and observed == expected


def region_bounds(sectors):
    nonempty = {k: v for k, v in sectors.items() if k != (0, 0)}
    singles = {k[0]: v for k, v in nonempty.items() if k[0] == k[1]}
    pairs = {
        (a, b): ea + eb
        for a, ea in nonempty.items()
        for b, eb in nonempty.items()
        if a[0] + b[0] == a[1] + b[1]
    }
    one_cell = min(singles.values())
    two_cells = min(pairs.values())
    assert one_cell == two_cells == 18
    assert [p for p, v in singles.items() if v == 18] == [1]
    assert all(a[0] + b[0] == 1 for (a, b), v in pairs.items() if v == 18)
    # An intact board always has an external edge from the unique excess
    # to a unit-filled port for every word at the lower equality. The one
    # vacancy at equality is private, hence cannot be its exterior neighbor.
    high_inventory = min(
        min(v for p, v in singles.items() if p >= 2),
        min(v for (a, b), v in pairs.items() if a[0] + b[0] >= 2),
        3 * min(nonempty.values()) + 1,  # even once total charge is zero
    )
    assert high_inventory >= 20
    return {
        "boundary_shut_floor_contrast_units": 9,
        "intact_board_floor_contrast_units": 10,
        "two_or_more_pairs_floor_contrast_units": high_inventory // 2,
    }


def mechanics(shape, locations, references, coupling):
    """Construct H minus the cold formation reference on a terminal cell."""
    index = {where: j for j, where in enumerate(locations)}
    h = np.zeros((len(index), len(index)))
    for j, where in enumerate(locations):
        valid = [v for v in where if v is not None]
        assert len(valid) == len(set(valid))
        word = np.ones(12, dtype=int)
        vacancy, excess = where
        if vacancy is not None:
            word[vacancy] = 0
        if excess is not None:
            word[excess] = 2
        value = coupling * (
            sum((word[a] - word[b]) ** 2 for a, b in shape.seams)
            + references @ (word**2 - 1)
        )
        h[j, j] += value
        for a, b in shape.seams:
            after = tuple(b if v == a else a if v == b else v for v in where)
            if after != where:
                h[j, j] += 0.5
                h[index[after], j] -= 0.5
    np.testing.assert_array_equal(h, h.T)
    return h


def clock_banks(shape, references, cutoff):
    """Apply the same coarse Fourier energy aperture to four signed load banks."""
    bank = {
        "cold_filled": [(None, None)],
        "vacancies": [(i, None) for i in range(12)],
        "excesses": [(None, i) for i in range(12)],
        "neutral_nontrivial": list(permutations(range(12), 2)),
    }
    energies = np.arange(1, int(cutoff * LABELS / PERIOD) + 1) * PERIOD / LABELS
    energies = energies[energies >= 0.25]
    record = {}
    for k in (8, 0):
        arms = {}
        for name, where in bank.items():
            h = mechanics(shape, where, references, k)
            spectrum = eigh(h, eigvals_only=True)
            assert max(abs(spectrum)) < PERIOD / 2
            phase = (energies[:, None] - spectrum) / PERIOD
            kernel = (np.sinc(phase * LABELS) / np.sinc(phase)) ** 2
            response = np.sum(kernel, axis=0)
            assert max(response) <= 1 + 1e-12
            arms[name] = {
                "material_words": len(where),
                "reported_vacancy_count": int(where[0][0] is not None),
                "reported_excess_count": int(where[0][1] is not None),
                "energy_lowest_and_highest": [float(spectrum[0]), float(spectrum[-1])],
                "computed_bright_per_uniform_sharp_enrollment": float(
                    np.mean(response)
                ),
                "computed_maximum_bright_over_sharp_energy_lines": float(max(response)),
            }
        record[str(k)] = arms
    assert record["8"]["neutral_nontrivial"]["energy_lowest_and_highest"][0] >= 72
    assert (
        record["8"]["neutral_nontrivial"][
            "computed_maximum_bright_over_sharp_energy_lines"
        ]
        < 1e-5
    )
    assert (
        record["0"]["neutral_nontrivial"][
            "computed_bright_per_uniform_sharp_enrollment"
        ]
        > 0.01
    )
    return record


def checked_federations():
    result = []
    for depth in range(3):
        board = federation(depth)
        shape = cell()
        ports = board.n
        external = [edge for edge in board.seams if edge[0] // 12 != edge[1] // 12]
        ends = [v for edge in external for v in edge]
        assert len(ends) == len(set(ends))
        assert all(sum(v // 12 == c for v in ends) == 3 for c in range(board.cells))
        assert set(board.seams) - set(external) == {
            (12 * c + a, 12 * c + b) for c in range(board.cells) for a, b in shape.seams
        }
        assert connected_components(board.adjacency(), directed=False)[0] == 1
        vacancy = next(p for p in range(ports) if p not in ends)
        excess = next(p for p in sorted(ends) if p not in board.neighbours[vacancy])
        assert board.degree[vacancy] == 5 and board.degree[excess] == 6
        result.append(
            {
                "cells": board.cells,
                "ports": ports,
                "sharp_pair_potential_floor_witness": [vacancy, excess],
                "exchange_only_first_positive_neutral_line_upper": (1 / 8)
                / (1 - 1 / board.cells),
            }
        )
    return result


def run():
    shape = cell()
    native = microscopic_check()
    sectors, equality = census(shape)
    subset_check(shape, sectors)
    bounds = region_bounds(sectors)
    b = federation(0)
    ref = 6 - np.array(b.degree[:12])
    assert ref.sum() == 9
    parent_receipt = HERE.parent / "passive_decay/records.json"
    passive = json.loads(parent_receipt.read_text())
    parent = passive["source_levels"][3]
    cutoff = parent / 10
    assert cutoff < 4 < KAPPA * min(bounds.values())
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "independent_parent_receipt_sha256": sha256(
            parent_receipt.read_bytes()
        ).hexdigest(),
        "native_blank_4H_micro_columns": native,
        "words_checked": 3**12,
        "local_counts_plus_minus_and_twice_shifted_floor": [
            [*k, v] for k, v in sectors.items()
        ],
        "intact_board_low_equality_placements_per_signed_sector": {
            str(k): v for k, v in equality.items()
        },
        "all_nonvacuum_neutral_formation_bounds": bounds,
        "boards": checked_federations(),
        "parent_formation_and_one_tenth_cutoff": [parent, cutoff],
        "finite_512_unit_period_32768_label_joint_charge_and_energy_read": clock_banks(
            shape, ref, cutoff
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    data = run()
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "results.json").read_text()), data)
        print(
            "Verified the exact all-stock native neutral bounds and the changed-law clock."
        )


if __name__ == "__main__":
    main()
