import json

import numpy as np
import pytest

from . import bracket, compiler, measure


def test_full_frozen_class_and_exact_bracket():
    candidates = bracket.ballot()
    assert sum(x["candidates"] for x in candidates.values()) == 23396
    assert {name: len(row["matched"]) for name, row in candidates.items()} == {
        "P0": 0,
        "P1": 1,
        "P2": 1,
        "F1": 0,
        "F2": 1,
        "G1": 0,
        "G2": 1,
    }
    assert candidates["P1"]["matched"][0] == [[["A1(2)", 1]], [["A1(1)", 1]]]
    _, frame, matrices, pivots = bracket.sample()
    assert len(pivots) == 3
    np.testing.assert_allclose(np.eye(3), frame.T @ frame, atol=1e-14)
    assert np.linalg.matrix_rank(matrices.reshape(12, 9), tol=1e-10) == 3


def test_twelve_address_compilation_and_excluded_dimensions():
    n = 12
    fourier = np.exp(2j * np.pi * np.arange(n)[:, None] * np.arange(n) / n) / np.sqrt(n)
    np.testing.assert_allclose(
        compiler.execute(compiler.compile_unitary(fourier)), fourier, atol=4e-12, rtol=0
    )
    assert compiler.to_gate(11, 9)
    with pytest.raises(ValueError):
        compiler.compile_unitary(np.eye(13))
    with pytest.raises(ValueError):
        compiler.compile_unitary(np.eye(0))


def test_replay_with_frame_center_spin_order_and_refusals():
    r = measure.run()
    assert r == json.loads((measure.HERE / "records.json").read_text())
    assert all(
        sum(counts) == r["starts_per_bank"]
        for b in r["banks"]
        for counts in b["counts"]
    )
    for k, row in r["covered"].items():
        if k.endswith(("/native", "/remounted")):
            assert row[0] == 63
    assert r["covered"]["ordered/cut_bridge"][0] == 0
    assert r["covered"]["ordered/half_contrast"][0] < 20
    assert r["covered"]["ordered/remounted"][2] < 20
    assert r["covered"]["ordered/native"][3] < 50
    rows = {b["word"] + "/" + b["setting"]: np.array(b["counts"]) for b in r["banks"]}
    bright = r["sources"].index(["pair", 0, 2, 0.0]) * 7 + r["analyzers"].index(
        ["pair", 0, 2, 0.0]
    )
    assert list(rows["center/native"][bright]) == [8192, 0, 0, 0]
    assert list(rows["turn/native"][bright]) == [0, 0, 8192, 0]
    assert np.all(rows["ordered/cut_bridge"][:, 3] == 8192)
    largest_difference = (
        np.max(abs(rows["ordered/native"] - rows["reversed/native"])) / 8192
    )
    assert largest_difference > 2 * r["simultaneous_radius"]
