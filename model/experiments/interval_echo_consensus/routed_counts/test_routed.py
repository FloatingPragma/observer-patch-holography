from itertools import combinations, permutations, product
import json

import numpy as np

from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import federation
from .flight import open_generator
from .material import Routing, attempt, certify
from .measure import HERE


def test_route_is_exact_pullback_on_the_actual_tree_and_can_hit_dirty_endpoints():
    board = federation(0)
    net = Routing(board)
    certify()
    chord = net.chords[0]
    domain = [*chord, 0, 83, 157]
    for i, j in permutations(dict.fromkeys(domain), 2):
        labels = list(range(board.n))
        legs = net.legs((i, j), chord)
        for a, b in legs:
            assert b in board.neighbours[a]
            assert {a, b} != set(chord)
            labels[a], labels[b] = labels[b], labels[a]
        assert [labels[v] for v in chord] == [i, j]


def test_all_words_and_local_forks_with_nonuniform_background():
    board = federation(0)
    net = Routing(board)
    sites = (0, 14, 82, 83)
    ordered = list(combinations(sites, 2))
    for word in product(range(3), repeat=len(sites)):
        # Unrelated unit loads can be present on every routing and chord endpoint.
        expected = min(word.count(0), word.count(2))
        for rev in (False, True):
            full = [1] * board.n
            for s, x in zip(sites, word, strict=True):
                full[s] = x
            used = 0
            for pair in ordered[::-1] if rev else ordered:
                chord = net.chords[-1 - used] if rev else net.chords[used]
                outcome, _ = attempt(full, net, pair, chord)
                used += int(outcome)
            assert used == expected
            assert full.count(2) == word.count(2) - expected
            assert full.count(0) == word.count(0) - expected


def test_complete_open_columns_by_direct_physical_seam_action():
    b = federation(0)
    h, c = open_generator(b)
    for v, j in [(0, 1), (12, 13), (0, 12), (0, 82), (137, 96)]:
        loads = [1] * b.n
        loads[v], loads[j] = 0, 2
        # Add independent native nine-state matrix elements in a word,
        # including each actual private port service and vacuum offset.
        kin = {}
        energy = 0
        for a, e in b.seams:
            x, y = loads[a], loads[e]
            energy += KAPPA * (x - y) ** 2
            if x != y:
                nv = e if v == a else a if v == e else v
                nj = e if j == a else a if j == e else j
                other = int(nv * b.n + nj)
                kin[other] = kin.get(other, 0) - 0.5
        energy += sum(
            KAPPA * (6 - int(deg)) * (load**2 - 1)
            for deg, load in zip(b.degree, loads, strict=True)
        )
        k = v * b.n + j
        assert c[k] == energy
        kin[k] = energy - sum(kin.values())
        col = h[:, k].tocoo()
        assert dict(zip(col.row, col.data, strict=True)) == kin


def test_counts_cover_every_launch_and_coarse_and_address_claims_separate():
    r = json.loads((HERE / "results.json").read_bytes())
    for g in r["sharp_boards"]:
        assert 2 * g["chords"] >= g["ports"]
        for t in g["banks"]:
            word = t["input_word"]
            assert len(t["positive_at"]) == len(t["marks"])
            if not t["exchange_disabled"]:
                assert len(t["marks"]) == min(word.count("0"), word.count("2"))
                assert t["quiet_remaining_pairs_when_normal"] is not None
    flight = r["coherent_flight"]
    starts = flight["starts_per_bank"]
    # Simultaneous 99% Hoeffding bound for at most 100 scored signed/bin indicators.
    margin = np.sqrt(np.log(2 * 100 / 0.01) / (2 * starts))
    with np.load(HERE / "counts.npz") as archive:
        for t in flight["cohorts"]:
            name = f"{t['mode']}_{int(t['reverse'])}_{t['reference']}"
            a = archive[name]
            assert a.shape == (flight["ports"] * (flight["ports"] - 1) // 2 + 1, 2)
            assert a.sum() == starts
            np.testing.assert_array_equal(
                a[flight["target_row"]], t["selected_pair_reference"]
            )
            np.testing.assert_array_equal(
                a[:-1].sum(axis=0), t["complete_marker_reference"]
            )
            np.testing.assert_array_equal(a[-1], t["quiet"])
            assert abs(t["marker"] / starts - 2 / 3) < margin
            assert max(abs(a[:-1].sum(axis=0) / starts - 1 / 3)) < margin
            if t["reference"] == "X":
                assert t["quiet"][1] == 0
                sep = abs(np.diff(t["selected_pair_reference"])[0]) / starts
                if t["mode"] == "native":
                    assert sep > 10 * margin
                else:
                    assert sep < 2 * margin


def test_dense_scan_uses_exactly_the_available_material_and_all_offers():
    saved = json.loads((HERE / "dense.json").read_bytes())
    word = saved["initial_word"]
    native = [b for b in saved["banks"] if not b["exchange_disabled"]]
    assert len(saved["occupied_memory_labels"]) == 11
    dirty = {tuple(m["chord"]) for m in saved["occupied_memory_labels"]}
    v, d = word.count(0), word.count(2)
    total = sum(word)
    for b in saved["banks"]:
        success = len(b["accepted"])
        assert b["quiet"] + success == saved["total_offers"] == 28680
        marks = {tuple(chord) for _, chord in b["accepted"]}
        assert not dirty.intersection(marks)
        assert len(marks) == success
        assert b["cost"] == saved["initial_cost"] - 2 * success
        if not b["exchange_disabled"]:
            assert b["cost"] == total + 2 * max(0, total - len(word))
            assert success == min(v, d)
    assert native[0]["accepted"] != native[1]["accepted"]
