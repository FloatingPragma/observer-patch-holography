"""Independent word, local instrument and finite radius checks."""

from fractions import Fraction as F
import json
from itertools import combinations

import networkx as nx
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

from . import measure, reader, words


def receipt():
    return json.loads((measure.HERE / "results.json").read_text())


def test_all_possible_first_crossings_and_minimal_paths():
    for level in range(3):
        board = federation(level)
        exits = measure.audit_board(level)["exits"]
        external = dict(exits)
        outputs = set()
        maximum_columns = 0
        for plus in combinations(range(12), 2):
            for minus in combinations(sorted(set(range(12)) - set(plus)), 2):
                active = (set(plus) | set(minus)) & external.keys()
                for site in active:

                    def cross(pair):
                        return tuple(
                            sorted(external[site] if i == site else i for i in pair)
                        )

                    word = cross(plus), cross(minus)
                    assert word not in outputs
                    outputs.add(word)
                maximum_columns = max(maximum_columns, len(active))
        assert maximum_columns == 3
        graph = nx.Graph(board.seams)
        owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
        hops = nx.single_source_shortest_path_length(owners, 0)
        distances = nx.multi_source_dijkstra_path_length(graph, range(12))
        for i, distance in distances.items():
            assert distance >= 2 * hops[i // 12] - 1


def test_native_and_moving_changed_contact_against_full_seam():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    for ratio in (F(0), F(1, 2), F(1), F(2)):
        for strength in measure.STRENGTHS:
            unitary = expm(
                -1j * np.pi * float(ratio) * (float(strength) * exchange + 8 * contrast)
            )
            p = measure.transfer_probability(ratio, strength)
            for a, b in ((0, 1), (2, 1)):
                column = unitary[:, law.SRC.index((a, b))]
                expected = np.zeros(99)
                expected[law.SRC.index((a, b))], expected[law.SRC.index((b, a))] = (
                    1 - p,
                    p,
                )
                np.testing.assert_allclose(abs(column) ** 2, expected, atol=1e-13)


def test_opposite_full_flight_responses_on_legal_words():
    planted = receipt()["planted_full_flight"]
    for row, sign in zip(planted, (1, -1), strict=True):
        trial = words.release(row["plus"], row["minus"])
        measure.verify(trial, row)
        t = F(1, trial["denominator"])
        sources = [[F(sign)] * 2, [F(sign)] * 2, [F(0)] * 2]
        lo, hi = measure.exact_bounds(t, F(1), 5, sources)["signed_squared_radius"]
        value = trial["signed_radius_per_squared_time"]
        err = trial["truncation_radius_per_squared_time"]
        assert float(lo) < value - err < value + err < float(hi)
        assert 0.499 < value * sign < 0.501


def test_recorded_slope_and_count_only_clock():
    out = receipt()
    assert (
        reader.census()["decoded_local_words"]
        == out["terminal_reader"]["decoded_local_words"]
    )
    with pytest.raises(ValueError):
        reader.decode([0])
    for size, laws in out["full_flight_bounds_per_squared_exchange_time"].items():
        for t, native in laws["1"].items():
            half = laws["1/2"][t]["signed_squared_radius"]
            assert all(x == 0 for x in laws["0"][t]["signed_squared_radius"])
            assert 0 < half[0] <= half[1] < native["signed_squared_radius"][0]
            assert native["slope_in_squared_wave_number"][1] < 0
            if size == "20":
                clock = out["local_clock_checks"][t]
                lo, hi = clock["duration_interval"]
                assert lo < float(F(t)) < hi
                assert clock["full_native_radius_per_squared_time"][0] > 0


def test_complete_banks_and_sharp_source_independent_of_scatter():
    report = receipt()["bridge_only"]
    direct = measure.sharp_source(json.loads(measure.PORTS.read_text()))
    for (counts, _), distribution in zip(
        direct, report["incoming_vacancy_unit_excess"], strict=True
    ):
        np.testing.assert_allclose(
            np.array(counts) / sum(counts), distribution, atol=0.0002
        )
    for panel in report["banks"].values():
        for row in panel.values():
            assert sum(row["counts"]) == report["starts"]
            assert sum(row["grid"]) == report["grid"]
            assert sum(row["probabilities"]) == pytest.approx(1)
    closed = report["banks"]["exchange=0, duration_over_pi=1"][
        "uniformly_chosen_first_receiver"
    ]
    assert closed["counts"][0] == closed["counts"][2] == 0
    with pytest.raises(ValueError):
        measure.verify([1, 2], [2, 1])
