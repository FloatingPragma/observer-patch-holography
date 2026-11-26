import json
from itertools import product

import numpy as np
import pytest

from . import collisions, measure


def test_native_blank_block_and_reference():
    for k in (0, 0.25, 0.5, 4, 8):
        np.testing.assert_allclose(
            measure.local_certificate(k), [2 * k - 1, 2 * k], atol=1e-13, rtol=0
        )
    np.testing.assert_allclose(measure.parity_effect(), np.full((2, 2), 0.5))
    n = np.arange(3)
    np.testing.assert_equal(
        n**2 - 1, np.maximum(1, n) ** 2 - 1 + np.minimum(1, n) ** 2 - 1
    )


def test_three_site_compression_against_full_native_matrices():
    # Path with independent referenced endpoints. Column compression includes
    # a swap into a colliding configuration only on the separated side.
    for k in (0.25, 8):
        contact = measure.contact(k)
        h = np.kron(contact, np.eye(3)) + np.kron(np.eye(3), contact)
        words = list(product(range(3), repeat=3))
        assert list(product(range(3), repeat=2)) == list(measure.SRC)
        h += np.diag([k * (a * a + c * c - 2) for a, b, c in words])
        positive = [words.index(tuple(max(1, t) for t in w)) for w in words]
        negative = [words.index(tuple(min(1, t) for t in w)) for w in words]
        product_energy = h[np.ix_(positive, positive)] * np.equal.outer(
            negative, negative
        ) + h[np.ix_(negative, negative)] * np.equal.outer(positive, positive)
        expected = np.zeros_like(h)
        for col, w in enumerate(words):
            for i, j in ((0, 1), (1, 2)):
                if abs(w[i] - w[j]) == 2:
                    swapped = list(w)
                    swapped[i], swapped[j] = swapped[j], swapped[i]
                    expected[col, col] += 2 * k - 0.5
                    expected[words.index(tuple(swapped)), col] -= 0.5
        np.testing.assert_allclose(h - product_energy, expected, atol=1e-13)
        idx = [i for i, w in enumerate(words) if w.count(0) == w.count(2) == 1]
        if k == 8:
            a = [i for i, w in enumerate(words) if w.count(2) == 1 and 0 not in w]
            b = [i for i, w in enumerate(words) if w.count(0) == 1 and 2 not in w]
            floor = sum(np.linalg.eigvalsh(h[np.ix_(z, z)])[0] for z in (a, b))
            assert np.linalg.eigvalsh(h[np.ix_(idx, idx)])[0] > floor


def test_held_contact_records():
    data = json.loads((measure.HERE / "counts.json").read_text())
    rows = measure.decode(data)
    assert all(
        abs(r["measured"] - r["forecast"]) < r["bound"] + r["forecast_bound"]
        for r in rows
    )
    native = [r for r in rows if r["k"] in (4, 8)]
    assert all(r["measured"] > r["bound"] for r in native)
    even = [r for r in rows if r["k"] == 0.25 and r["phase"] == "even"]
    assert len(even) == 4 and all(r["measured"] < -r["bound"] for r in even)
    data["banks"][1] = data["banks"][0]
    with pytest.raises(ValueError):
        measure.decode(data)


def test_independent_fragment_collisions():
    data = json.loads((measure.HERE / "collisions.json").read_text())
    assert [r["counts"] for r in collisions.acquire()["banks"]] == [
        r["counts"] for r in data["banks"]
    ]
    assert all(
        r["probability_lower"] > 0 and r["gap"] > 0 for r in collisions.decode(data)
    )
