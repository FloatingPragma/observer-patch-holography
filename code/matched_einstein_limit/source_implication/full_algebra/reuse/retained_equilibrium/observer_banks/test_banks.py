"""Verify every terminal marginal, all spectral brackets and independently decoded secants."""

from copy import deepcopy
import json

import numpy as np
import pytest

from . import audit, measure


def archive(name):
    return json.loads((audit.HERE / name).read_text())


def test_all_integer_marginals_and_independent_finite_sample_bounds():
    raw, cert = archive("records.json"), archive("certificate.json")
    assert cert["denominator"] == audit.GRID
    assert audit.marginal_check(raw) > 0
    out = audit.audit(raw, cert["dyadic_root_boxes"])
    assert out == cert["secant_intervals_millionths"]
    assert out["entropy"][0] > 0 > out["area"][1]
    for key in ("volume", "ordinary_residual"):
        assert out[key][0] < 0 < out[key][1]
    assert out["double_gain_control"][1] < -380000
    other = measure.decode(raw["counts"])["secant"]
    for key, (lo, hi) in out.items():
        assert lo / 1e6 < other[key]["value"] < hi / 1e6


def test_changed_terminal_law_is_rejected():
    raw = deepcopy(archive("records.json"))
    row = raw["counts"][0]["clock"][0]
    # Leave the shot count invariant but move a macroscopic fraction of the true timeout bank.
    loss = int(row[-1]) // 4
    row[-1] = str(int(row[-1]) - loss)
    row[0] = str(int(row[0]) + loss)
    with pytest.raises(AssertionError):
        audit.marginal_check(raw)


def test_source_law_and_existing_finite_difference_target():
    raw = archive("records.json")
    for current, reference in zip(measure.prepare(), raw["tickets"], strict=True):
        for key in ("pauli", "clock"):
            assert np.max(np.abs(np.array(current[key]) - reference[key])) <= 1
    paired = archive("../results.json")["held"]["finite"][0]
    cert = archive("certificate.json")["secant_intervals_millionths"]
    for key in ("entropy", "area", "volume"):
        low, high = cert[key]
        assert low * 1e-6 < paired[key] < high * 1e-6
