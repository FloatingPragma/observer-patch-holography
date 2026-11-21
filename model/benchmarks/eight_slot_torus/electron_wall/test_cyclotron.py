import json

import numpy as np
import pytest

from model.benchmarks.eight_slot_torus.electron_wall import cyclotron as c
from model.benchmarks.eight_slot_torus.electron_wall import walk as w


def test_magnetic_grades_against_literal_update():
    n, length = 7, 3
    shape = (n, 1, 1)
    frame = np.array([1, 1j])
    z = np.array([1, -1])
    rng = np.random.default_rng(15)
    for q in (-1, 1):
        conn = w.links(shape)
        conn[1, :, 0, 0] = np.exp(2j * np.pi * np.arange(n) / n)
        for grade in (1, -1):
            seed = rng.normal(size=(n, length, 2)) + 1j * rng.normal(size=(n, length, 2))
            def encode(v):
                a = v * frame / np.sqrt(2)
                return np.stack([a, grade * a[:, ::-1] * z], axis=1)[:, None, None]
            state = w.project(w.step(w.prepare(encode(seed)), conn, q))
            matrix = c.magnetic_block(n, length, q, grade)
            right = matrix @ seed.ravel()
            np.testing.assert_allclose(12 * state, encode(right.reshape(seed.shape)), atol=1e-12)
            eigs = np.linalg.eigvalsh(matrix.toarray())
            # The finite table contains every solver-returned pole, even when both probes miss it.
            measured = c.modes(128, length, q, grade, 6)
            for weights in measured["weights"]:
                assert len(weights) == 7 and sum(weights) == pytest.approx(1.)
            assert (abs(eigs) < 12).all()


def test_registered_failure_and_weak_series():
    base = json.loads((c.HERE / "cyclotron.json").read_text())
    previous = json.loads((c.HERE / "cyclotron_lines.json").read_text())
    weak = json.loads((c.HERE / "weak_lines.json").read_text())
    target = w.rest_properties(base["private_length"])["g"]
    for archive, periods in ((previous, base["periods"]), (weak, [4096, 8192, 16384])):
        assert [(r["period"], r["charge"]) for r in archive] == [(n, q) for n in periods for q in base["charges"]]
        c.check_reversal(archive)
    assert not c.initial_gate(previous, target)
    assert c.weak_gate(weak, previous, target)
    for missing in ([], weak[:-1], [r for r in weak if r["charge"] == 1], weak[:2] + weak):
        assert not c.weak_gate(missing, previous, target)
    nonfinite = [dict(r, g=float("nan")) for r in weak]
    assert not c.weak_gate(nonfinite, previous, target)
    # CI regenerates both charges on the first torus. The CLI replays all larger carriers.
    fresh = c.experiment(dict(base, periods=[base["periods"][0]]))
    c.compare(previous[:2], fresh)
    bad = json.loads(json.dumps(previous[:2]))
    bad[0]["grades"][0]["weights"] = 0
    with pytest.raises(ValueError):
        c.compare(bad, fresh)
    bad = json.loads(json.dumps(previous[:2]))
    bad[0]["grades"][0]["assigned"][0] += 1e-12
    with pytest.raises(ValueError, match="integer label"):
        c.compare(bad, fresh)
