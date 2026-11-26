"""Independent readout of the material selector, support graphs and counted launches."""

from functools import cache
from math import log
import numpy as np

from .apparatus import federation
from .selector import published, SHOTS, SEED
from .test_count import word


@cache
def bank():
    return published("selector")


def test_routed_selector_word_is_measured_and_reversible():
    result = bank()
    cfg = published("checkpoint")["config"]
    board = federation(0)
    reserved = set(cfg["owners"] + cfg["auditors"] + [cfg["boundary"][1] // 12])
    edges = result["dimers"]
    selected = {p for pair in edges for p in pair}
    assert len(selected) == 16 and not reserved & {p // 12 for p in selected}
    assert all(
        b in board.neighbours[a]
        and a // 12 == b // 12
        and board.degree[a] == board.degree[b] == 5
        for a, b in edges
    )
    for station in result["station_routes"]:
        assert (
            station["sharp"] == [1, 0, 2]
            and {r["port"] for r in station["paths"]} == selected
        )
        for r in station["paths"]:
            path = r["path"]
            assert (
                path[0] == r["port"]
                and path[-1] == station["head"]
                and set(path).isdisjoint(station["references"])
            )
            assert all(y in board.neighbours[x] for x, y in zip(path, path[1:]))
    for j, r in enumerate(result["selector_truth_table"]):
        assert r["raw"][0] == r["raw"][1]
        pairs = list(zip(r["raw"][0][::2], r["raw"][0][1::2]))
        assert all(pair in ((0, 1), (1, 0)) for pair in pairs)
        measured = sum((pair == (1, 0)) << i for i, pair in enumerate(pairs))
        assert measured == r["decoded"] == j
    assert result["plate"] == {
        "parked": [[1, 0], [0, 1]],
        "half": [[0.5, 0.5], [0.5, 0.5]],
        "full": [[0, 1], [1, 0]],
    }


def test_independent_conditional_counts_and_hidden_vs_indexed_grades():
    result = bank()
    checkpoint = published("checkpoint")
    words = [word(r, 8)[0] for r in checkpoint["raw"]]
    assert words == result["public_inputs"]
    rng = np.random.Generator(np.random.PCG64(SEED))
    report_words = np.array(words)
    total_bins = 0
    ell = log(2 * 262144 / 0.01)
    for name, panel in result["arms"].items():
        width = panel["width"]
        reported = np.asarray(panel["terminal_by_measured_pattern"])
        counts = np.asarray(panel["counts_by_measured_pattern"])
        terminal = np.asarray(panel["terminal_counts"])
        # Reconstruct every receiving outcome independently using tree bridge
        # inventory coordinates and the separate prior sharp read on its final word.
        patterns = list(range(1 << width))
        shifts = (
            [((1 << (j - 1)) if j else 0) for j in patterns]
            if name == "one"
            else patterns
        )
        expect = report_words[np.arange(256)[:, None] ^ shifts]
        np.testing.assert_array_equal(reported, expect)
        # Each start is one PCG64 integer below 2**40; its top `width` bits
        # pick one of the equally weighted selector patterns.
        reproduced = np.array(
            [
                np.bincount(
                    rng.integers(0, 1 << 40, size=SHOTS, dtype=np.int64)
                    >> (40 - width),
                    minlength=len(shifts),
                )
                for _ in range(256)
            ]
        )
        np.testing.assert_array_equal(counts, reproduced)
        hidden = np.zeros((256, 256), int)
        np.add.at(hidden, (np.arange(256)[:, None], expect), counts)
        np.testing.assert_array_equal(hidden, terminal)
        assert np.all(counts > 0) and np.all(terminal.sum(axis=1) == SHOTS)
        allowed = np.zeros((256, 256), bool)
        np.put_along_axis(allowed, expect, True, axis=1)
        np.testing.assert_array_equal(terminal > 0, allowed)
        prob = allowed / len(shifts)
        radius = np.sqrt(2 * SHOTS * prob * (1 - prob) * ell) + 2 * ell / 3
        loss = np.abs(terminal - SHOTS * prob) / radius
        assert (
            loss.max() <= 1 and abs(loss.max() - panel["max_allowance_fraction"]) < 1e-8
        )
        total_bins += terminal.size
        ports = panel["certificate"]
        admissible = [set(row) for row in expect.tolist()]
        assert ports["edges"] == sum(
            not admissible[i].isdisjoint(other)
            for i in range(256)
            for other in admissible[i + 1 :]
        )
        for initial, code_words in (
            (range(128), ports["fixed_class_code_words"]),
            (range(256), ports["code_words"]),
        ):
            code = [words.index(w) for w in code_words]
            union = set().union(*(admissible[i] for i in initial))
            assert all(i in initial for i in code)
            assert set().union(*(admissible[i] for i in code)) == union
            assert len(union) == len(code) * len(shifts)
        # A visible diagnostic log distinguishes the permutation, making
        # (reported selector pattern, checkpoint word) reversible for every input.
        assert all(len(set(output)) == 256 for output in expect.T)
        if name == "one":
            assert ports["known_indexed_one_swap_sizes"] == [128, 256]
    assert total_bins == result["checked_terminal_bins"] == 262144


def test_subset_draw_covers_the_complete_finite_word_group():
    result = bank()
    groups = published("continuations")["families"]
    row_order = np.argsort(result["public_inputs"])
    for arm, title in (("internal", "repeat_inside"), ("boundary", "repeat_boundary")):
        traces = np.asarray(result["arms"][arm]["terminal_by_measured_pattern"])[
            row_order
        ]
        assert {tuple(column) for column in traces.T} == {
            tuple(r["terminal"]) for r in groups[title]["routes"]
        }
    assert [
        (
            len(p["certificate"]["fixed_class_code_words"]),
            len(p["certificate"]["code_words"]),
        )
        for p in result["arms"].values()
    ] == [(16, 32), (1, 2), (1, 1), (128, 256)]
