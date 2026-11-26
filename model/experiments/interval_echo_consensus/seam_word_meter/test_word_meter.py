import gzip
import json
from fractions import Fraction
from itertools import product
from math import pi

import numpy as np
import sympy as sp

from . import apparatus as app
from .measure import HERE, series

RECEIPT = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def test_exact_ternary_effects_and_dirty_collar():
    # Exact Gaussian-rational propagator after the integer contrast phase.
    eye = sp.eye(9)
    swap = eye[:, [app.LAW.SRC.index((b, a)) for a, b in app.LAW.SRC]]
    gate = (1 - sp.I) / 2 * eye + (1 + sp.I) / 2 * swap
    np.testing.assert_allclose(
        np.array(gate).astype(complex), app.pulse(pi / 2), atol=app.TOL
    )
    projectors = [sp.diag(*[int(a == k) for a, b in app.LAW.SRC]) for k in range(3)]
    for p in projectors:
        mate = swap * p * swap
        for obj in (p, mate):
            reduced = sum(
                (j * gate.H * obj * gate * j for j in projectors), sp.zeros(9)
            )
            assert reduced.applyfunc(sp.simplify) == (p + mate) / 2
    assert RECEIPT["local"]["initial_first_port_then_two_echoes"] == [0, 1, 0]
    assert RECEIPT["local"]["phase_inputs_initial_unread_read"] == {
        "1": [0.5, 0, 0.5],
        "-1": [0.5, 1, 0.5],
    }


def test_every_word_against_independently_labeled_classical_paths():
    # Given sharp input, a resolved position record followed by a half-exchange
    # has the same subsequent position statistics as a labeled path with a fair swap.
    # The coherent branch conditioned on all reads is checked separately.
    for idx, word in enumerate(product(range(30), repeat=2)):
        actual = np.asarray(RECEIPT["two_offers"][idx])
        expected = np.zeros_like(actual)
        for start in range(12):
            for actions in product((0, 1), repeat=2):
                site = start
                reads = 0
                for edge, do_swap in zip(word, actions, strict=True):
                    a, b = app.SEAMS[edge]
                    reads = 2 * reads + int(site == a)
                    if do_swap and site in (a, b):
                        site = a + b - site
                expected[reads, site, start] += 1
        np.testing.assert_array_equal(actual, expected)


def test_dependent_calendar_and_coherent_return():
    held = RECEIPT["held"]
    returned = dict((tuple(w), p) for w, p in held["returned_word_law"])
    assert set(returned) == set(series())
    for w, p in series().items():
        assert abs(returned[w] - float(p)) < 1e-10
    assert app.schedule({(): Fraction(1)}) == ({(): 1}, {(): 1.0}, [])
    for p in (0, 1, np.sqrt(2) / 3):
        np.testing.assert_allclose(app.calendar_bit(p), [1 - p, p], atol=app.TOL)
    for row in held["banks"].values():
        for case in row["modes"].values():
            for branch in case["conditional_records"]:
                records = np.array(
                    [v for _, v in branch["paired_read_transcript_and_port"]]
                )
                assert abs(records.sum() - 1) < 2e-9
    assert (
        abs(held["banks"]["localized"]["joint_minus_independent_TV"] - 1159 / 27000)
        < 1e-11
    )
    for key in ("quarter_phases", "localized"):
        modes = held["banks"][key]["modes"]
        assert not np.allclose(
            modes["native"]["all_launch_terminal"],
            modes["unread"]["all_launch_terminal"],
        )


def test_arbitrary_coherent_reference_and_fixed_real_calibration():
    rng = np.random.default_rng(14147)
    v = rng.normal(size=(9, 5)) + 1j * rng.normal(size=(9, 5))
    v /= np.linalg.norm(v)
    joint = np.einsum("ij,ij->i", v.conj(), v).real.reshape(3, 3)
    score = np.array([-7.5, 1.5, np.sqrt(19)])
    expected = (score @ joint.sum(axis=0) + score @ joint.sum(axis=1)) / 2
    before = v.T.conj() @ v
    after = np.zeros_like(before)
    values = np.zeros(2)
    for p in app.projectors():
        u = app.pulse(pi / 2) @ np.kron(p, np.eye(3)) @ v
        after += u.conj().T @ u
        weights = (abs(u) ** 2).sum(axis=1).reshape(3, 3)
        values += score @ np.stack([weights.sum(axis=0), weights.sum(axis=1)], axis=1)
    np.testing.assert_allclose(after, before, atol=app.TOL)
    np.testing.assert_allclose(values, [expected, expected], atol=app.TOL)


def test_uniform_mean_band_and_orientation_multiplicity():
    row = RECEIPT["repeated_seam_response"]
    assert row["distinct_maps"] == 30
    levels = (
        1 - np.array([5 + np.sqrt(5)] * 3 + [6] * 5 + [5 - np.sqrt(5)] * 3 + [0]) / 60
    )
    np.testing.assert_allclose(row["eigenvalues"], levels, atol=1e-10)


def test_deep_words_and_rare_positive_branches():
    long_word = (0,) * 1400
    assert app.schedule({long_word: 1}) == ({long_word: 1}, {long_word: 1.0}, [])
    rare = Fraction(1, 10**400)
    w0, w1 = long_word + (0,), long_word + (1,)
    exact, numeric, nodes = app.schedule({w0: 1 - rare, w1: rare, (5,) * 1401: 0})
    assert exact == {w0: 1 - rare, w1: rare} and set(numeric) == {w0, w1}
    assert nodes[0]["native_binary64_probs"] is None and numeric[w1] is None
    x = sp.sqrt(2) / 3
    ideal, observed, _ = app.schedule({(0,): 1 - x, (29,): x})
    assert ideal == {(0,): 1 - x, (29,): x}
    assert abs(observed[(29,)] - float(x)) < app.TOL
