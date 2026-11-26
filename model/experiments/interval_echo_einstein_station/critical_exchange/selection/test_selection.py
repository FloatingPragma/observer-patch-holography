from itertools import combinations
import json
from math import comb

import numpy as np
import pytest

from . import measure
from ..measure import ring_modes


def test_absorption_rate_on_every_word_including_degenerate_energies():
    n = 10
    k, _, filled = ring_modes(n)
    energies = -np.cos(k).round(12)
    choices = measure.work_choices([(e, e) for e in energies])
    upper = set(np.flatnonzero(~filled))
    for word in combinations(range(n), n // 2):
        word = set(word)
        m = len(upper & word)
        change = []
        for i, j in combinations(range(n), 2):
            if j in word:
                i, j = j, i
            if i in word and j not in word and choices[i][j]:
                change.append(int(i in upper) - int(j in upper))
        assert change.count(1) == m**2
        assert all(x in (0, 1) for x in change)
        if m == 0:
            assert not change


def test_alternating_sharp_site_word_matches_the_measured_orbital_distribution():
    n = 6
    _, phase, _ = ring_modes(n)
    for word in combinations(range(n), n // 2):
        probability = abs(np.linalg.det(phase[::2, list(word)]))**2
        separate_pairs = len({i % (n // 2) for i in word}) == n // 2
        assert probability == pytest.approx(2**(-n // 2) if separate_pairs else 0)


def test_same_calibration_selects_a_cold_or_hot_source_without_lost_starts():
    saved = json.loads((measure.HERE / "results.json").read_text())
    audit = measure.build(saved["counts_symmetric_antisymmetric_quiet"])
    measure.audit(audit, saved)
    assert audit["independent_pair_count"] == comb(246, 2)
    assert 0.022 < audit["resolved_gap_lower"] < audit["exact_gap"]
    for name, record in audit["trajectories"].items():
        assert not record["unfinished"]
        assert record["final_upper_count"] == (123 if name == "reversed" else 0)
        lo, hi = record["work_from_count_intervals"]
        assert lo < record["released_hopping_work"] < hi
        assert hi < 0 if name == "reversed" else lo > 0


def test_a_calibration_with_no_resolved_gap_is_rejected():
    with pytest.raises(ValueError, match="misses|separate"):
        measure.build([[measure.SHOTS // 2, measure.SHOTS // 2, 0]] * 246)
