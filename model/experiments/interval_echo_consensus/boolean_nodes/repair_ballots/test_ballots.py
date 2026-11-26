"""Compare every material trace with the canonical single-site replacement definition."""

import json
from functools import lru_cache
from itertools import combinations, product

import numpy as np

from .measure import HERE, Apparatus

GRAPHS = json.loads((HERE / "results.json").read_text())["graphs"]


def neighborhood(edges, n):
    return [[b if a == i else a for a, b in edges if i in (a, b)] for i in range(n)]


def repairs(bits, adjacent):
    changed = {}
    for i, neighbors in enumerate(adjacent):
        for new in (0, 1):
            if (
                new != bits[i]
                and any(bits[i] != bits[j] for j in neighbors)
                and all(new == bits[j] for j in neighbors)
            ):
                changed[i] = bits[:i] + (new,) + bits[i + 1 :]
    return changed


def encoded(flags):
    return "".join("1" if x else "0" for x in flags)


def test_all_material_attempts_against_canonical_replacements():
    for name, bank in GRAPHS.items():
        adj = neighborhood(bank["edges"], len(bank["runs"][0]["input"]))
        for row in bank["runs"]:
            bits = tuple(map(int, row["input"]))
            mode = row["mode"]

            def survey(word):
                return [
                    encoded(word[i] != word[j] for j in ns) for i, ns in enumerate(adj)
                ]

            truth = survey(bits)
            sensor = ["0" * len(ns) for ns in adj] if mode == "parked" else truth
            assert row["local_survey"] == sensor
            flags = {i for i, s in enumerate(sensor) if s and "0" not in s}
            initial = sum(bits[i] != bits[j] for i, j in bank["edges"])
            assert (
                len(row["attempts"])
                == row["budget"]
                == (0 if mode == "short" else initial)
            )
            for t, (site, observed, guard, wrote, withdrawn) in enumerate(
                row["attempts"]
            ):
                legal = repairs(bits, adj)
                scan = []
                order = [
                    (t + row["direction"] * k) % len(bits) for k in range(len(bits))
                ]
                expected = order[0]
                if mode != "blind":
                    for i in order:
                        scan.append(i in flags)
                        if i in flags:
                            expected = i
                            break
                assert (expected, encoded(scan)) == (site, observed)
                expected_guard = (
                    "0" * len(adj[site]) if mode == "parked" else survey(bits)[site]
                )
                assert expected_guard == guard
                assert bool(wrote) == (mode != "parked" and site in legal)
                if mode == "native":
                    assert flags == set(legal)
                    assert bool(legal) == bool(wrote)
                if wrote:
                    bits = legal[site]
                    if name == "sphere":
                        assert len(adj[site]) == 3
                changed = (
                    [site] if mode == "stale" else ([site, *adj[site]] if wrote else [])
                )
                assert withdrawn == [i for i in changed if i in flags]
                flags.difference_update(changed)
            assert row["output"] == encoded(bits)
            assert row["latched"] == encoded(i in flags for i in range(len(bits)))
            assert row["final_survey"] == survey(bits)
            assert row["service_probes_exchanges"][0] >= 2 * len(bank["edges"])
            assert row["audit_probes_exchanges"][0] == 3 * len(bits) + 2 * len(
                bank["edges"]
            )
            if mode == "native":
                assert not repairs(bits, adj) and not flags


def test_every_possible_terminal_word_on_small_carriers():
    for name in ("path", "pentagon"):
        graph = GRAPHS[name]
        n = len(graph["runs"][0]["input"])
        adj = neighborhood(graph["edges"], n)

        @lru_cache(None)
        def terminals(word):
            children = repairs(word, adj).values()
            if not children:
                return {word}
            return set.union(*(terminals(child) for child in children))

        competing = 0
        reaches_unanimity = 0
        forces_unanimity = 0
        for word in product((0, 1), repeat=n):
            lit = set(repairs(word, adj))
            predicted = set()
            for count in range(len(lit) + 1):
                for subset in combinations(lit, count):
                    selected = set(subset)
                    if any(selected.intersection(adj[i]) for i in selected):
                        continue
                    if any(selected.isdisjoint(adj[i]) for i in lit - selected):
                        continue
                    predicted.add(
                        tuple(b ^ int(i in selected) for i, b in enumerate(word))
                    )
            assert terminals(word) == predicted
            assert (len(predicted) > 1) == any(
                i in lit and j in lit for i, j in graph["edges"]
            )
            competing += len(predicted) > 1
            uniform = any(len(set(end)) == 1 for end in terminals(word))
            reaches_unanimity += uniform
            forces_unanimity += uniform and len(predicted) == 1
        assert competing == graph["inputs_with_multiple_terminal_words"]
        assert reaches_unanimity == graph["inputs_with_a_unanimous_terminal_word"]
        assert forces_unanimity == graph["inputs_forcing_a_fixed_unanimous_word"]


def test_full_board_census_by_complement_pairs():
    graph = GRAPHS["sphere"]
    adj = neighborhood(graph["edges"], 20)
    counts = {}
    competing_total = 0
    # Each pair of complementary words has exactly one member with owner 19 dark.
    for start in range(0, 2**19, 4096):
        words = np.arange(start, start + 4096, dtype=np.uint32)
        bits = ((words[:, None] >> np.arange(20, dtype=np.uint32)) & 1).astype(bool)
        enabled = np.column_stack(
            [
                np.all(bits[:, ns] == ~bits[:, i : i + 1], axis=1)
                for i, ns in enumerate(adj)
            ]
        )
        left, right = np.array(graph["edges"]).T
        competing_total += 2 * int(
            np.any(enabled[:, left] & enabled[:, right], axis=1).sum()
        )
        quiet = bits[~np.any(enabled, axis=1)]
        rank = (quiet[:, left] != quiet[:, right]).sum(axis=1)
        for k, v in zip(*np.unique(rank, return_counts=True), strict=True):
            counts[str(k)] = counts.get(str(k), 0) + 2 * int(v)
    assert counts == graph["quiescent_by_residual_mismatches"]
    assert sum(counts.values()) == 133630 and counts["0"] == 2
    assert competing_total == graph["inputs_with_multiple_terminal_words"] == 474412


def test_empty_incidence_is_an_eligibility_stutter():
    row = Apparatus(2, []).run((0, 1), 1, "native")
    assert row["output"] == "01" and row["latched"] == "00" and row["budget"] == 0


def test_sphere_agreement_basin_by_independent_sets():
    graph = GRAPHS["sphere"]
    adj = neighborhood(graph["edges"], 20)

    @lru_cache(None)
    def count_independent(vertices):
        if not vertices:
            return 1
        i = min(vertices)
        without = vertices - {i}
        return count_independent(without) + count_independent(
            without.difference(adj[i])
        )

    # The carrier contains a pentagon. Its two color classes cannot both be independent.
    assert all((i + 1) % 5 in adj[i] for i in range(5))
    count = count_independent(frozenset(range(20)))
    assert graph["inputs_with_a_unanimous_terminal_word"] == 2 * count == 11656
    assert graph["inputs_forcing_a_fixed_unanimous_word"] == 7736
