"""Audit the transported preparation on actual physical port words and local effects."""

from itertools import product
from math import inf, pi

import numpy as np
import pytest

from . import apparatus as app
from . import chart


def test_common_charts_independent_coordinates_and_pair_mean():
    report = chart.scalar_chart_report()
    a, b = report["charts"]
    assert a["joint_sharp_work_words"] == 1024
    assert b["joint_sharp_work_words"] == 4096
    assert [r["requested"] for r in a["coordinates"]] == [
        r["requested"] for r in b["coordinates"]
    ]
    assert a["coordinates"][0]["sharp_data_support"] == [0]
    assert a["coordinates"][-1]["sharp_data_support"] == [2]
    assert all(r["sharp_data_support"] == [0, 2] for r in b["coordinates"])
    labels = np.array(list(product((0, 2), repeat=12)))
    for config in (a, b):
        lo, hi = config["endpoints"]
        rows = config["coordinates"]
        local = np.array([r["ternary_read"] for r in rows])
        probs = local[np.arange(12), labels].prod(axis=1)
        np.testing.assert_allclose(probs.sum(), 1, atol=8e-13)
        expectations = probs @ labels
        requested = [r["requested"] for r in rows]
        np.testing.assert_allclose(
            lo + (hi - lo) * expectations / 2, requested, atol=4e-12
        )
        np.testing.assert_allclose(
            [r["decoded_scalar"] for r in rows], requested, atol=4e-12
        )
        # Compare the predicted scalar mean against all native output columns
        # of a prospective overlap, using the actual chart preparation.
        gate = app.pulse(pi / 2)
        for i, j in app.SEAMS:
            pair = np.diag(np.kron(local[i], local[j]))
            migrated = np.diag(gate @ pair @ gate.conj().T).real
            digits = np.array(app.LAW.SRC)
            recovered = lo + (hi - lo) * (migrated @ digits) / 2
            np.testing.assert_allclose(
                recovered, (requested[i] + requested[j]) / 2, atol=4e-12
            )


def test_stage_act_unstage_factors_entangled_spectators_on_all_routes():
    lay = chart.stock_layout()
    board = app.federation(0)
    h0, h1 = lay["head"]
    head_gate = app.pulse(0.937)
    blank = app.LAW.SRC.index((0, 2))
    for i, stock in enumerate(lay["stock_ports"]):
        paths = lay["paths_to_head"][i]
        # Use a transit port and a head port as the two entangled spectator
        # coordinates. They are restored even if they move during staging.
        intermediate = next(
            port
            for path in paths
            for port in path[1:-1]
            if port not in (i, stock, h0, h1)
        )
        operations = [edge for path in paths for edge in chart.exchange_word(path)]
        assert all(b in board.neighbours[a] for a, b in operations)
        actual = {}
        expected = {}
        for level, phase in ((0, 1 / np.sqrt(2)), (2, 1j / np.sqrt(2))):
            word = [(u + 1) % 3 for u in range(board.n)]
            word[i], word[stock] = 0, 2
            word[h0], word[intermediate] = level, 2 - level
            original = word.copy()
            for u, v in operations:
                word[u], word[v] = word[v], word[u]
            input_col = app.LAW.SRC.index((word[h0], word[h1]))
            for row, (x, y) in enumerate(app.LAW.SRC):
                if abs(head_gate[row, input_col]) < 1e-13:
                    continue
                observed = word.copy()
                observed[h0], observed[h1] = x, y
                for u, v in reversed(operations):
                    observed[u], observed[v] = observed[v], observed[u]
                actual[tuple(observed)] = phase * head_gate[row, input_col]
            for row in (app.LAW.SRC.index((2, 0)), blank):
                x, y = app.LAW.SRC[row]
                end = original.copy()
                end[i], end[stock] = x, y
                expected[tuple(end)] = phase * head_gate[row, blank]
        assert actual.keys() == expected.keys()
        np.testing.assert_allclose(
            list(actual.values()), [expected[k] for k in actual], atol=4e-12
        )


@pytest.mark.parametrize(
    "values,lo,hi",
    [
        ([0] * 11, -1, 1),
        ([0] * 11 + [1.00001], -1, 1),
        ([-1.00001] + [0] * 11, -1, 1),
        ([0] * 12, 1, 1),
        ([0] * 12, 1, -1),
        ([float("nan")] * 12, -1, 1),
        ([0] * 12, -inf, 1),
    ],
)
def test_outside_or_degenerate_chart_declined_before_operation(
    values, lo, hi, monkeypatch
):
    def forbidden(*unused):
        raise AssertionError("An invalid chart reached the physical apparatus")

    monkeypatch.setattr(app, "pulse", forbidden)
    with pytest.raises(ValueError):
        chart.prepare_chart(values, lo, hi)
