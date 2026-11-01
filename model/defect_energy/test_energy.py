"""Independent native ordered-particle and arbitrary-reference checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

import pytest

from tools import strict_json
from ..observer import record_recovery as records
from ..defect_protection import protection
from ..exchange_sector import statistics as exchange
from . import energy
from .__main__ import OUTPUT, build, check


def ordered_step(law, ordered, seam, coin):
    """Track two distinguishable labels directly through the native update."""
    state = tuple(int(p in ordered) for p in range(12))
    out, record = records.recorded_step(law, state, seam, coin)
    before, after = set(ordered), {p for p, n in enumerate(out) if n}
    removed, added = before-after, after-before
    if removed:
        assert len(removed) == len(added) == 1
        old, new = removed.pop(), added.pop()
        ordered = tuple(new if p == old else p for p in ordered)
    return ordered, record


@pytest.mark.parametrize("law", ("C", "E"))
def test_all_complete_lifted_columns_against_native_label_tracking(law):
    expected = [Counter() for _ in range(132)]
    for r, pair in product((0, 1), exchange.PAIRS):
        ordered = pair if r == 0 else pair[::-1]
        source = r*66+exchange.INDEX[pair]
        for seam, coin in product(records.seams(), (0, 1)):
            final, _ = ordered_step(law, ordered, seam, coin)
            target_pair = tuple(sorted(final))
            target = int(final != target_pair)*66+exchange.INDEX[target_pair]
            expected[source][target] += 1
    actual, _ = energy.transfer_counts(law)
    assert actual == tuple(tuple(sorted(row.items())) for row in expected)
    action = energy.action_certificate(law)
    assert action["source_columns_checked"] == 3960
    assert action["nonzero_parity_columns"] == 208


@pytest.mark.parametrize("law", ("C", "E"))
@pytest.mark.parametrize("word", (energy.EVEN_WORD, energy.ODD_WORD,
                                  energy.MATCHED_EVEN_WORD, energy.MATCHED_ODD_WORD))
def test_every_word_outcome_against_native_ordered_history(law, word):
    actual = energy.word_instrument(law, energy.INITIAL_PAIR, word)
    assert len(actual) == 2**len(word)
    assert sum(row["probability"] for row in actual) == 1
    for row in actual:
        ordered = energy.INITIAL_PAIR
        movements = 0
        for seam, coin in zip(word, row["coins"], strict=True):
            following, _ = ordered_step(law, ordered, seam, coin)
            movements += ordered != following
            ordered = following
        assert row["final_pair"] == tuple(sorted(ordered))
        assert row["memory_flip"] == int(ordered != tuple(sorted(ordered)))
        assert row["movements"] == movements
        assert row["memory_flip"] == row["orientation_toggles"] % 2


@pytest.mark.parametrize("law", ("C", "E"))
def test_closed_action_rechart_and_exact_matched_resources(law):
    value = energy.matched_cost_certificate(law)
    left, right = value["selected_nonrenormalized_witnesses"]
    canonical = lambda row: Counter((a, tuple(b), c, d) for a, b, c, d in row["history"])
    assert canonical(left) == canonical(right)
    assert value["state_dependent_orientation_toggle_counts"] == [2, 1]
    assert value["common_compiled_resources"]["orientation_feedback_pulses"] == 27908
    assert value["common_compiled_resources"]["native_pulses"] == (
        3561005 if law == "C" else 3561026)
    lookup = {(r.seam, r.coin, r.raw): r for r in protection.branches(law)}
    for witness in (left, right):
        # All single-vertex chart changes generate every Z2 vertex rechart.
        for changed in range(66):
            chart = tuple(int(i == changed) for i in range(66))
            initial = exchange.INDEX[energy.INITIAL_PAIR]
            target = initial
            for _, seam, coin, raw in witness["history"]:
                target = dict(protection.lifted_pairs(
                    lookup[tuple(seam), coin, raw], chart))[target]
            assert target % 66 == initial
            assert target // 66 == witness["memory_flip"]


@pytest.mark.parametrize("law", ("C", "E"))
def test_exact_reference_safe_one_step_memory_channel(law):
    # Independently push one half of a Bell state through the native ordered
    # process.  The exact Choi state verifies the entire memory channel,
    # including off-diagonal operators and any untouched reference.
    choi = [[F() for _ in range(4)] for _ in range(4)]
    for pair, seam, coin in product(exchange.PAIRS, records.seams(), (0, 1)):
        labels = []
        for r in (0, 1):
            final, _ = ordered_step(law, pair if r == 0 else pair[::-1], seam, coin)
            labels.append(2*int(final != tuple(sorted(final)))+r)
        for a, b in product(labels, repeat=2):
            choi[a][b] += F(1, 66*60*2)
    p = F(26, 495)
    expected = [[F() for _ in range(4)] for _ in range(4)]
    for indices, weight in (((0, 3), 1-p), ((1, 2), p)):
        for a, b in product(indices, repeat=2):
            expected[a][b] += weight/2
    assert choi == expected
    one, two = energy.memory_channel(law, 1), energy.memory_channel(law, 2)
    assert one["memory_flip_probability"] == p
    assert two["memory_flip_probability"] == F(733, 7425)
    assert two["Y_Z_eigenvalue"]-one["Y_Z_eigenvalue"]**2 == F(398, 245025)
    assert 0 < p < F(1, 2)  # distance to every pure state is at least p
    for n in (0, 1, 2, 3, 6, 12, 64):
        assert energy.memory_channel(law, n)["expected_native_movements"] == F(5*n, 33)


@pytest.mark.parametrize("attempts", (-1, 65, True, 1.0, "1"))
def test_bounded_exact_domain_rejects_invalid_times(attempts):
    with pytest.raises(ValueError):
        energy.memory_channel("C", attempts)


@pytest.mark.parametrize("law,pair,word", (("other", (0, 11), ()),
    ("C", (True, 11), ()), ("C", (0, 0), ()), ("C", (0, 11), ((0, 1),)*9),
    ("C", (0, 11), ((0, 2),)), ("C", (0, 11), ([0, 1],))))
def test_word_domain_rejects_non_native_inputs(law, pair, word):
    with pytest.raises(ValueError):
        energy.word_instrument(law, pair, word)


def test_receipt_exact_replay_and_tamper_rejection():
    expected = build()
    actual = strict_json.load(OUTPUT)
    check(actual, expected)
    actual["laws"]["C"]["one_step_distance_from_every_unitary_lower_bound"] = "0"
    with pytest.raises(ValueError, match="exact replay"):
        check(actual, expected)
