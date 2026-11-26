import json
from copy import deepcopy

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_encounter.conjugation import defect_generator
from . import measure as m
from . import weak


def test_principal_single_defect_and_quiet_orbital():
    board = m.geometry()[0]
    order = list(m.PORTS)
    for sign in (-1, 1):
        _, h8 = m.sector((sign,), 8, True)
        full = defect_generator(board, sign, 8, 0)
        np.testing.assert_array_equal(
            h8.toarray(), 8 * full[order, :][:, order].toarray()
        )
        for ps in (-1, 1):
            m.closed_certificate(sign, ps, 8)


def test_actual_and_deformed_sectors_obey_the_registered_gate_tests():
    raw = json.loads((m.HERE / "counts.json").read_text())
    assert all(
        all(v["resolved"] == (row["bridge"] == "open") for v in row["exterior"])
        for row in m.decode(raw)["studies"]
    )
    for strength, opened, ss, ps, t, bins in raw["banks"]:
        if ps == 0:
            assert sum(np.array(bins).reshape(3, 9)[:, 4]) == m.STARTS
    invalid = deepcopy(raw)
    invalid["banks"] += invalid["banks"][:1]
    with pytest.raises(ValueError, match="duplicate"):
        m.decode(invalid)


def test_weak_opening_against_independent_two_sided_differentiation():
    ss, ps = 1, -1
    h0, v, initial, bins = weak.source(ss, ps)
    step = 0.0005
    t = 2
    out = []
    for g in (-step, 0, step):
        psi = expm_multiply(t * (h0 + g * v), initial)
        out.append(np.bincount(bins % 9, weights=abs(psi) ** 2, minlength=9))
    np.testing.assert_allclose(
        (out[2] + out[0] - 2 * out[1]) / (2 * step**2),
        weak.quadratic(ss, ps)[t][1],
        atol=3e-6,
        rtol=0,
    )


def test_weak_registration_and_remainder():
    raw = json.loads((m.HERE / "weak_counts.json").read_text())
    summary = weak.report(raw)
    assert summary["quadratic_with_remainder"] == summary["held_bins"] == 216
    invalid = deepcopy(raw)
    invalid["calibration"][1][-1] = weak.STARTS
    with pytest.raises(ValueError, match="intensity"):
        weak.report(invalid)
