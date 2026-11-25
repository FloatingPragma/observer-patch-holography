from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import log, prod, sqrt

import numpy as np
import sympy as sp

from . import measure as m


def test_native_energy_and_three_way_source():
    data = json.loads((m.HERE / "records.json").read_text())
    assert json.loads(json.dumps(m.run())) == data
    census = data["released_energy"]
    assert census["host_seams"] == 630
    for a in census["arms"]:
        s, b = Fraction(a["s"]), Fraction(a["b"])
        assert Fraction(a["host_sum_mean"]) == (2795 - 65 * b) / 4
        assert Fraction(a["ordinary_edge_mean"]) == 65 * (1 - b) / 4
        meter = a["half_pulse_meter"]
        assert meter["arrivals"] + meter["refusals"] == meter["launches"]
        assert Fraction(meter["predicted_arrival"]) * 65 == Fraction(
            a["ordinary_edge_mean"]
        )
        assert abs(meter["arrivals"] / meter["launches"] - float((1 - b) / 4)) < sqrt(
            log(1000) / (2 * meter["launches"])
        )
        assert meter["parked_arrivals"] == 0
        expected = []
        for p in m.CORNERS:
            x, y, z = [2 * v - 11 for v in p]
            expected.append(int(4 * (1 + s * x * y * z + b * x)))
        assert expected == a["weights_per_32"]
    for a in census["arms"][:3]:
        assert (
            a["prob_each_energy_32_and_33_per_64"]
            == census["arms"][0]["prob_each_energy_32_and_33_per_64"]
        )
        # Every ordered port pair has the same reduced classical state,
        # including off-diagonal choices of stored source correlations.
        for i, j in combinations(range(3), 2):
            marg = Counter()
            for p, w in zip(m.CORNERS, a["weights_per_32"], strict=True):
                marg[p[i], p[j]] += w
            assert marg == dict.fromkeys(product((5, 6), repeat=2), 8)
    independent_entropies = [
        -sum(w / 32 * log(w / 32) for w in arm["weights_per_32"])
        for arm in census["arms"][:3]
    ]
    assert (
        abs(
            independent_entropies[1]
            - independent_entropies[0]
            - data["released_entropy_gain_from_minus_half_to_zero_nats"]
        )
        < 2e-15
    )


def test_stopped_exit_probabilities_and_all_launches():
    # Build the eight-corner killed chain independently by graph adjacency.
    corners = list(product((5, 6), repeat=3))
    q = sp.Matrix(
        8,
        8,
        lambda i, j: sp.Rational(
            2
            if i == j
            else sum(abs(a - b) for a, b in zip(corners[i], corners[j], strict=True))
            == 1,
            8,
        ),
    )
    # Each corner has three different exterior moves; each keeps all three signs.
    exit_scores = sp.Matrix(
        [sp.Rational(9, 8) * prod(2 * x - 11 for x in p) for p in corners]
    )
    assert list((sp.eye(8) - q).inv() * exit_scores) == [
        prod(2 * x - 11 for x in p) for p in corners
    ]
    data = json.loads((m.HERE / "records.json").read_text())
    for arm in data["angular_banks"]:
        s = Fraction(arm["s"])
        scores = np.array([prod(2 * x - 11 for x in p) for p in arm["collectors"]])
        exits = np.array(arm["exit_counts_by_tick_and_collector"])
        assert exits.sum() == sum(arm["observed_launch_counts"]) == arm["shots"]
        assert set(scores) == {-3, 3}
        prop = exits[:, scores > 0].sum() / arm["shots"]
        assert Fraction(arm["predicted_positive_score"]) == (1 + s / 3) / 2

        # Nine simultaneous comparisons (eventual and first two exit times).
        def within(sample, mean, n):
            return abs(sample - mean) <= sqrt(log(18 / 0.01) / (2 * n))

        assert within(prop, float((1 + s / 3) / 2), arm["shots"])
        for i in range(2):
            n = exits[i].sum()
            assert within(
                exits[i, scores > 0].sum() / n,
                float((1 + s * Fraction(-1, 5) ** i) / 2),
                n,
            )


def test_modular_reference_and_held_terminal_read():
    data = json.loads((m.HERE / "records.json").read_text())
    ref = data["reference_modular"]
    per_sign = Counter()
    for corner, n in zip(m.CORNERS, ref["reference_counts_by_corner"], strict=True):
        per_sign[m.angular_score(corner)] += n
    n = ref["reference_shots"]
    assert sum(per_sign.values()) == n
    # Three bounded Bernoulli estimates at joint confidence at least 99%.
    qtol = sqrt(log(6 / 0.01) / (2 * n))
    assert abs(per_sign[1] / n - 1 / 4) < qtol
    coefficient = log(per_sign[-1] / per_sign[1]) / 2
    first, last = data["angular_banks"][0], data["angular_banks"][-1]
    score_change = (
        6 * (last["positive_score_count"] - first["positive_score_count"]) / n
    )
    for entry in (first, last):
        assert (
            abs(
                entry["positive_score_count"] / n
                - (1 + float(Fraction(entry["s"])) / 3) / 2
            )
            < qtol
        )
    assert (
        abs(
            coefficient * score_change
            - ref["observed_fixed_reference_modular_change_nats"]
        )
        < 1e-14
    )
    reference, final = [
        data["released_energy"]["arms"][i]["weights_per_32"] for i in (0, 2)
    ]
    modular_change = -sum(
        (v - w) / 32 * log(w / 32) for w, v in zip(reference, final, strict=True)
    )
    entropy_change = sum(
        w / 32 * log(w / 32) - v / 32 * log(v / 32)
        for w, v in zip(reference, final, strict=True)
    )
    assert (
        abs(modular_change - ref["predicted_fixed_reference_modular_change_nats"])
        < 1e-14
    )
    assert abs(entropy_change - ref["predicted_held_entropy_change_nats"]) < 1e-14


def test_directed_ranges_and_fired_clock():
    from .. import measure as apparatus

    for panel in m.ranges():
        pos, k = panel["points"], panel["step"]
        for i, j in combinations(range(4), 2):
            measured = panel["outward_times"][f"{i}{j}"]
            # Return to the emitting station using the reversed native π contacts.
            assert m.active_range(pos[j], pos[i]) == measured
            assert measured == max(
                abs(a - b) for a, b in zip(pos[i], pos[j], strict=True)
            )
        assert (panel["outward_times"]["02"], panel["outward_times"]["13"]) == (
            2 * k,
            2 * k,
        )
    clock = m.population_clock()
    assert clock["shell_totals_targets_and_eight_times_s_and_b_slopes"] == [
        [(2 * k) ** 3 - (2 * k - 2) ** 3, 0, 0] for k in range(1, 7)
    ]
    signed = clock["targets_and_eight_times_nonzero_correlation_slope"]
    assert len(signed) == 48 and sum(row[-1] for row in signed) == 0
    assert m.harmonic_check()[-1]["interior_addresses"] == 1000
    assert all(apparatus.release_corner(p) == p for p in m.CORNERS)


def test_quadratic_null_fit_on_each_body_diagonal_shell():
    clock = m.population_clock()
    signs = list(product((-1, 1), repeat=3))
    modes = [i for n in (0, 1, 2) for i in combinations(range(3), n)]
    tensor = sp.Matrix([[prod(t[i] for i in mode) for mode in modes] for t in signs])
    assert tensor.T * tensor == 8 * sp.eye(7)
    for row in clock["eight_diagonal_tensor_tests"]:
        n = row["shell"]
        slopes = {
            tuple(-1 if x < 6 else 1 for x in (a, b, c)): sp.Rational(h, 8)
            for a, b, c, h in clock["targets_and_eight_times_nonzero_correlation_slope"]
            if max(abs(2 * x - 11) for x in (a, b, c)) == 2 * n - 1
        }
        measured = sp.Matrix([slopes[t] for t in signs])
        assert tensor.row_join(measured).rank() == 8
        assert tensor.T * measured == sp.zeros(7, 1)
        assert max(abs(v) for v in measured) == sp.Rational(
            row["best_uniform_tensor_error"]
        )
