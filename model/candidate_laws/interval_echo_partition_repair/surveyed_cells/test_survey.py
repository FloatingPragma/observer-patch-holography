"""Independent checks of the release law, overlapping cuts and held predictions."""

from fractions import Fraction as F
from itertools import combinations, product
import json
from math import isclose, log

import numpy as np
import pytest

from ..regional_check import check_saved
from . import cropped
from . import measure as m


def test_repair_and_all_pair_held_marginals():
    for biased in (False, True):
        weights = dict(
            zip(m.TAGS, m.tag_weights("ordinary", biased=biased), strict=True)
        )
        for axis in range(3):
            for tag in m.TAGS:
                for dest, p in m.controller_step(tag, axis, biased).items():
                    assert (
                        weights[tag] * p
                        == weights[dest] * m.controller_step(dest, axis, biased)[tag]
                    )
        if not biased:
            for start in m.TAGS:
                row = {start: F(1)}
                for axis in range(3):
                    nxt = {}
                    for x, mass in row.items():
                        for y, p in m.controller_step(x, axis).items():
                            nxt[y] = nxt.get(y, F(0)) + p * mass
                    row = nxt
                assert row == weights
    for axes in combinations(range(3), 2):
        for value in product((0, 1), repeat=2):
            for t in (-m.STEP, m.STEP):
                assert sum(
                    p
                    for k, p in zip(m.TAGS, m.tag_weights("held", t), strict=True)
                    if tuple(k[i] for i in axes) == value
                ) == F(1, 4)


def test_overlapping_restrictions_and_fiber_sizes():
    small = m.CUTS["overlap_pair"]
    for t in (-m.STEP, m.STEP):
        weights = m.tag_weights("held", t)
        expected = m.histogram(small, weights)
        for big in [m.CUTS["ordinary_line"], m.CUTS["held_sheet"]]:
            restricted = {}
            for (cells, _, word), p in m.histogram(big, weights).items():
                sub = tuple(cells[big.index(x)] for x in small)
                large_bonds = m.interior_registers(cells)
                sub_word = sum(
                    ((word >> large_bonds.index(e)) & 1) << i
                    for i, e in enumerate(m.interior_registers(sub))
                )
                key = sub, (), sub_word
                restricted[key] = restricted.get(key, F(0)) + p
            assert restricted == expected
    for cut in m.CUTS.values():
        ref = m.histogram(cut, m.tag_weights("held"))
        c = len({(cells, ctl) for cells, ctl, _ in ref})
        for (cells, _, _), p in ref.items():
            bonds = m.interior_registers(cells)
            a, v = m.rod_survey(cells)
            assert 2 * len(bonds) == 6 * v - a
            assert p == F(1, c * 2 ** len(bonds))


def test_proof_on_arbitrary_nonuniform_law():
    # At fixed cell volume any observed distribution, including a sharp word,
    # satisfies the entropy gap; the condition applies on its full support.
    for cut in m.CUTS.values():
        ref = m.histogram(cut, m.tag_weights("held"))
        keys = list(ref)
        p = np.arange(1, len(keys) + 1, dtype=float)
        p /= p.sum()
        pi = np.array(list(map(float, ref.values())))
        sizes = {c: m.rod_survey(c)[0] for c in {x for x, _, _ in keys}}
        areas = np.array([sizes[cells] for cells, _, _ in keys])

        def entropy(w):
            return -sum(w * np.log(w))

        change = entropy(p) - entropy(pi) + log(2) / 2 * (p - pi) @ areas
        assert isclose(change, -sum(p * np.log(p / pi)), rel_tol=1e-12)


def test_sector_control_and_frozen_banks():
    result = json.loads((m.HERE / "results.json").read_text())
    replay = m.run()
    for name, row in replay["regions"].items():
        assert row["records"] == result["regions"][name]["records"]
    check_saved(replay, result)
    assert result["ordinary_eta_99"][0] < log(2) / 2 < result["ordinary_eta_99"][1]
    assert result["held_residual_99"][0] < 0 < result["held_residual_99"][1]
    for name, expected in [
        ("ordinary_line", -1),
        ("held_sheet", 1),
        ("overlap_pair", 0.5),
    ]:
        cut = m.CUTS[name]
        arm = "ordinary" if name == "ordinary_line" else "held"
        si, ai, v, dg = m.projection_slopes(cut, arm)
        assert ai == expected and v == 0
        assert isclose(si, -log(2) * ai / 2, abs_tol=3e-15)
        sc, ac, _, _ = m.projection_slopes(cut, arm, True)
        assert isclose(sc + log(2) * ac / 2, -float(dg) * log(2), abs_tol=3e-15)
        observed = result["regions"][name]["measured_secant_99"]
        assert observed["entropy"][0] < si < observed["entropy"][1]
        assert observed["area"][0] < ai < observed["area"][1]


def test_replay_rebuilds_all_derived_fields_without_using_a_sampler(monkeypatch):
    saved = json.loads((m.HERE / "results.json").read_text())

    def unexpected_generator(*args):
        raise AssertionError("frozen bank replay used a sampler")

    monkeypatch.setattr(m.np.random, "Generator", unexpected_generator)
    only_counts = {key: saved[key] for key in ("seed", "trials_per_bank")}
    only_counts["regions"] = {name: {"records": row["records"]} for name, row in saved["regions"].items()}
    rebuilt = m.replay(only_counts)
    check_saved(rebuilt, saved)
    saved["held_residual_99"][0] += 1e-8
    with pytest.raises(ValueError, match="held_residual_99"):
        check_saved(rebuilt, saved)


def test_balanced_count_tampering_fails_even_if_statistics_barely_change():
    saved = json.loads((m.HERE / "results.json").read_text())
    counts = saved["regions"]["ordinary_line"]["records"]["1"][m.BULK]
    counts[0] -= 1
    counts[1] += 1
    assert sum(counts) == m.TRIALS
    with pytest.raises(ValueError, match="SHA-256"):
        m.replay(saved)


def test_model_checks_word_order_detector_counts_and_exact_integers(monkeypatch):
    cut, arm, sign = m.CUTS["ordinary_line"], "ordinary", 1
    saved = json.loads((m.HERE / "results.json").read_text())
    bank = saved["regions"]["ordinary_line"]["records"]["1"]
    words, sector, energy = m.terminal_probabilities(cut, arm, sign)
    with monkeypatch.context() as patch:
        patch.setattr(m, "terminal_probabilities", lambda *args: (tuple(reversed(words)), sector, energy))
        with pytest.raises(ValueError, match="bulk word"):
            m.verify_model_bank(cut, arm, sign, bank)
    for channel in ("rod_sector_hits", "ordinary_energy_hits"):
        with pytest.raises(ValueError, match="probability bound"):
            m.verify_model_bank(cut, arm, sign, dict(bank, **{channel: 0}))
    with pytest.raises(ValueError, match="integer counts"):
        m.verify_model_bank(cut, arm, sign, dict(bank, rod_sector_hits=float(bank["rod_sector_hits"])))


def test_fresh_seed_generates_independent_records_and_obeys_the_model_bounds():
    saved = json.loads((m.HERE / "results.json").read_text())
    fresh = m.generate(m.SEED + 1)
    assert fresh["seed"] == m.SEED + 1
    for name, cut in m.CUTS.items():
        banks = fresh["regions"][name]["records"]
        assert banks != saved["regions"][name]["records"]
        for sign in (-1, 1):
            m.verify_model_bank(cut, "ordinary" if name == "ordinary_line" else "held", sign, banks[str(sign)])


def test_actual_smooth_ball_cropping_depends_on_admission():
    for r, h in [(F(1), F(1, 8)), (F(1, 4), F(1, 32)), (F(1, 8), F(1, 512))]:
        a = 4 * r * r
        assert cropped.cut_area_per_pi(r, h) == 3 * a / 2 - h * h / 2
        assert cropped.midpoint_count_per_pi(r, h) == 4 * r**3 / h**3
