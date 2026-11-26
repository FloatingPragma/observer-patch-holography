import json
import itertools

import numpy as np
import pytest

from . import measure as m
from . import analyzer, apparatus


def test_native_exterior_action_in_every_local_sector():
    for size in (1, 2, 3):
        f = analyzer.fourier(size)
        true = np.array(
            [
                [analyzer.minor(f, y, x) for x in range(1 << size)]
                for y in range(1 << size)
            ]
        )
        np.testing.assert_allclose(apparatus.programmed(f, 8), true, atol=3e-13, rtol=0)


def test_all_start_replay_and_exact_one_and_two_reference_caps():
    saved = json.loads((m.HERE / "records.json").read_text())
    assert json.loads(json.dumps(m.acquire())) == saved
    scores, width = m.score(saved)
    for key, ideal in zip(
        m.MENUS, (0.5, 0.75, (1 + 1 / np.sqrt(2)) / 2, 0.5, 0.5), strict=True
    ):
        assert abs(scores[key] - ideal) < width
    lo, hi = [scores["native", r] for r in (1, 2)]
    assert hi - lo > 2 * width and lo - 0.5 > width
    for L in range(7):
        A = np.diag([0.5] * L, 1) + np.diag([0.5] * L, -1) if L else np.zeros((1, 1))
        assert np.isclose(np.linalg.eigvalsh(A)[-1], np.cos(np.pi / (L + 2)))
    saved["banks"].pop()
    with pytest.raises(ValueError, match="Missing preparation"):
        m.score(saved)


def test_receiver_uses_material_report_and_preserves_ambiguous_ends():
    for bit in (0, 1):
        assert m.decode_word(apparatus.announce(bit)) == 1 + bit
    assert m.decode_word("11") is None
    for size, station in itertools.product((1, 2, 3), (0, 1)):
        for mask in range(1 << size):
            assert m.decode_word(apparatus.word_read(station, size, mask)) == mask
    for size in (1, 2, 3):
        for extreme in (0, (1 << size) - 1):
            _, guess = analyzer.receiver(size, extreme)
            assert set(guess.values()) == {int(extreme == 0)}


def test_recorder_keeps_all_physical_refusals(monkeypatch):
    word_read, announce = apparatus.word_read, apparatus.announce
    rng = np.random.default_rng(m.SEED)
    # Each of the four detector stages can record an erasure or invalid code.
    for stage, failure in itertools.product(
        ("tag", "near", "bus", "far"), (None, "11")
    ):

        def broken_word(station, size, mask):
            if stage == "tag" and station == 0 and size == 1:
                return failure
            if stage == "near" and station == 0 and size == 3:
                return None if failure is None else failure * size
            if stage == "far" and station == 1 and size == 3:
                return None if failure is None else failure * size
            return word_read(station, size, mask)

        monkeypatch.setattr(apparatus, "word_read", broken_word)
        if stage == "bus":
            monkeypatch.setattr(
                apparatus,
                "announce",
                lambda bit: None if failure is None else failure * 2,
            )
        else:
            monkeypatch.setattr(apparatus, "announce", announce)
        bank = m.bank("path_read", 2, 0, rng)
        hits, refusals, hist = m.count_bank(bank)
        assert (hits, refusals) == (0, m.SHOTS)
        assert hist.sum() == (0 if stage in ("tag", "near") else m.SHOTS)
        assert all(guess is None for _, _, _, _, guess, _ in bank["records"])


def test_bob_uses_received_word_even_if_a_valid_bit_is_changed(monkeypatch):
    announce = apparatus.announce
    monkeypatch.setattr(apparatus, "announce", lambda bit: announce(1 - bit))
    bank = m.bank("native", 0, 0, np.random.default_rng(m.SEED))
    _, refusals, _ = m.count_bank(bank)
    assert refusals == 0
    for _, near, payload, _, guess, _ in bank["records"]:
        delivered = m.received_word(payload, 1)
        assert delivered != m.decode_word(near, 1)
        assert guess == int(delivered == 0)


def test_missing_counts_forged_guess_and_malformed_words_are_distinct_from_refusal():
    saved = json.loads((m.HERE / "records.json").read_text())
    bank = saved["banks"][0]
    row = next(row for row in bank["records"] if row[-1] > 0)
    for index, bad, reason in (
        (-1, -1, "Invalid enrolled count"),
        (-1, 0.5, "Invalid enrolled count"),
        (-1, True, "Invalid enrolled count"),
        (-1, row[-1] - 1, "Incomplete enrolled bank"),
        (-2, 1 - row[-2], "Published guess"),
        (-2, bool(row[-2]), "Published guess"),
        (1, "0", "Malformed detector word"),
        (2, [], "Malformed message packet"),
        (3, "0x", "Malformed detector word"),
    ):
        before = row[index]
        row[index] = bad
        with pytest.raises(ValueError, match=reason):
            m.count_bank(bank)
        row[index] = before


def test_signwise_forecast_cannot_be_hidden_by_correct_pooled_control():
    saved = json.loads((m.HERE / "records.json").read_text())
    plus, minus = saved["banks"][-2:]
    assert plus["source_sign"] == 1 and minus["source_sign"] == 0
    # Shift the same number of starts to the same measured plus outcomes in both
    # banks. Opposite signs yield opposite accuracy errors and cancel when pooled.
    for bank in (plus, minus):
        source = next(row for row in bank["records"] if row[-2] == 0 and row[-1] > 2000)
        sink = next(row for row in bank["records"] if row[-2] == 1 and row[-1] > 2000)
        source[-1] -= 2000
        sink[-1] += 2000
    with pytest.raises(m.ForecastMismatch, match="Per-sign accuracy"):
        m.score(saved)
