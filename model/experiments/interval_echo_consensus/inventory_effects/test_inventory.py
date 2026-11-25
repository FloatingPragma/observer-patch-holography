import json
from itertools import product

import numpy as np
import pytest
from copy import deepcopy
from scipy.linalg import expm

from . import certificate as cert
from . import circuit as c
from . import measure as m
from . import symmetry


def test_blank_formula_and_all_spectator_separation():
    record = cert.certify()
    assert record["checked_distinct_positive_gaps_per_species"][-1] == [8, [729] * 3]
    swap = np.eye(9)[:, [3 * b + a for a, b in product(range(3), repeat=2)]]
    np.testing.assert_allclose(expm(-1j * np.pi * c.local_h(c.K)), swap, atol=3e-14)
    # Enumerate a full three-port Hamiltonian in the tensor basis and isolate
    # matrix units by looking at the ordered-double-commutator eigenvalues.
    words = list(product(range(3), repeat=3))
    for a, b in cert.PAIRS:
        q = np.array([(131 + 1) * w[1] ** 2 + cert.diagonal(w[1], w[2]) for w in words])
        edges = [
            (i, words.index((b, a, w[2])))
            for i, w in enumerate(words)
            if w[:2] == (a, b)
        ]
        assert len({abs(q[i] - q[j]) for i, j in edges}) == 3
        for i, j in edges:
            assert abs(q[j] - q[i]) == cert.detuning(a, b, [words[i][2]])


def test_native_reader_on_both_columns_and_all_three_effects():
    y = np.array([[0, -1j], [1j, 0]])
    h = np.eye(2) / 4 + y / 8
    native_input = np.eye(729, dtype=complex)[:, c.INDEX[:2]]
    for menu, desired in (
        ("fine", [h, h, np.eye(2) - 2 * h]),
        ("union", [2 * h, np.eye(2) - 2 * h, 0 * h]),
        ("single", [h, np.eye(2) - h, 0 * h]),
    ):
        branches = c.play(c.compile_reader(c.target(menu)), native_input, c.K)
        outside = np.delete(branches, c.INDEX, axis=0)
        np.testing.assert_allclose(outside, 0, atol=1e-12)
        for j, target in enumerate(desired):
            v = branches[c.INDEX[2 * j : 2 * j + 2], :]
            np.testing.assert_allclose(v.conj().T @ v, target, atol=1e-12)


def test_reader_counts_and_cubic_control():
    data = json.loads((m.HERE / "records.json").read_bytes())
    for key, record in data["law"].items():
        for bank in record["banks"].values():
            assert all(sum(bank[src]) == n for src, n in m.N.items())
        verdict = m.score(record)
        assert verdict == record["readback"]
        flags = [
            d["distinguished"]
            for row in verdict["comparisons"].values()
            for d in row.values()
        ]
        assert any(flags) == (float(key) == 0)
        if float(key) == c.K:
            assert verdict["historical_single_against_cubic"]["distinguished"]
            assert verdict["historical_union_against_cubic"]["distinguished"]
            first = (
                np.array(record["banks"]["fine"]["plus_y"])[:2].sum() / m.N["plus_y"]
            )
            second = (
                np.array(record["banks"]["fine"]["minus_y"])[:2].sum() / m.N["minus_y"]
            )
            assert abs(first - 3 / 8) < 0.003 and abs(second - 1 / 8) < 0.003


def test_sharp_symmetry_and_count_forecast():
    out = symmetry.run()
    assert out == json.loads((m.HERE / "symmetry.json").read_bytes())
    for mode in ("balanced", "swapped", "unbalanced"):
        for src, n in (("zero", m.N["zero"]), ("one", m.N["one"])):
            raw = np.array(out["banks"][mode][src]).reshape(7, 8)
            assert raw.sum() == n and not raw[2:].any()
            if mode == "unbalanced":
                assert raw[int(src == "one")].sum() == n
            else:
                assert max(abs(raw[:2] / n - 1 / 16).ravel()) < 0.002


def test_missing_bank_and_minimum_disagreement():
    old = json.loads((m.HERE / "records.json").read_bytes())["law"]["8.0"]
    assert m.score(old)["maximum_cross_context_disagreement_upper_99"] < 1 / 256
    for damage in ("source", "lost", "bins"):
        row = deepcopy(old)
        if damage == "source":
            del row["banks"]["single"]["one"]
        if damage == "lost":
            row["banks"]["fine"]["one"][0] -= 1
        if damage == "bins":
            row["source_reads"]["zero"].append(0)
        with pytest.raises(ValueError):
            m.score(row)


def test_counts_determine_positive_density_and_predict_noncommuting_readers():
    from . import tomography as t

    rec = json.loads((m.HERE / "tomography.json").read_bytes())
    assert np.linalg.norm(np.cross(*t.DIRECTIONS)) > 0.1
    for k, rows in rec["law"].items():
        for row in rows.values():
            r = t.readback(row["banks"])
            assert r == row["state_readback"]
            assert r["positive_at_99"]
            assert all(x["passes"] == (float(k) == c.K) for x in r["forecasts"])
            broken = deepcopy(row["banks"])
            broken["X"][0][0] += 1
            with pytest.raises(ValueError):
                t.readback(broken)
