import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from . import apparatus, bound, measure

HERE = Path(__file__).resolve().parent


def independent_chain(n, holes):
    """Blank physical-edge Hamiltonian on sharp words, built from actual route adjacency."""
    arms = apparatus.fixture()[:2]
    r = arms[0] if n == 8 else arms[1]
    board = federation(0)
    sealed = set(apparatus.fixture()[2])
    edges = [
        (r.index(a), r.index(b))
        for a, b in board.seams
        if a in r and b in r and (a, b) not in sealed
    ]
    assert len(edges) == n - 1
    words = []
    for x in range(n):
        for zs in combinations([i for i in range(n) if i != x], holes):
            words.append(tuple(2 if i == x else 0 if i in zs else 1 for i in range(n)))
    ix = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    for j, w in enumerate(words):
        for left, r in edges:
            z = list(w)
            z[left], z[r] = z[r], z[left]
            h[j, j] += 0.5 + 8 * (w[left] - w[r]) ** 2
            h[ix[tuple(z)], j] -= 0.5
    return words, h


def test_native_local_meter_and_closed_route():
    record = apparatus.geometry_record()
    assert record["receiver_owners"] == [6, 4]
    assert len(record["silent_seams"]) == 44
    blocks, frames = apparatus.local_frames(8)
    np.testing.assert_array_equal(blocks[0], [[8.5, -0.5], [-0.5, 17]])
    np.testing.assert_array_equal(blocks[1], [[49.5, -0.5], [-0.5, 49.5]])
    np.testing.assert_allclose(frames[1][:, 0], [2**-0.5] * 2)
    # The same physical private-quarter and half-exchange decode either
    # antisymmetric load order for the independent parity calibration.
    for order in [(0, 1), (0, 2)]:
        f = apparatus.rotate_frame(np.array([[0.5, -0.5], [-0.5, 0.5]]), 8, order)
        np.testing.assert_allclose(
            abs(f.T @ [2**-0.5, -(2**-0.5)]) ** 2, [0, 1], atol=1e-12
        )


def test_two_free_wings_independent_matrix_exponential():
    saved = json.loads((HERE / "results.json").read_text())
    frames, _ = measure.spectra()
    factors = []
    for n, holes in [(8, 2), (6, 1)]:
        words, h = independent_chain(n, holes)
        native_words, hn = apparatus.sector(n, holes, 1, 8)
        ix = [native_words.index(w) for w in words]
        np.testing.assert_array_equal(h, hn[np.ix_(ix, ix)])
        a = [(0, 1), (3, 4)] if n == 8 else [(1, 2)]
        d = [(0, 1), (5, 6)] if n == 8 else [(3, 4)]
        f = frames if n == 8 else [frames[1]]
        s = apparatus.frame_columns(words, a, f)[0][:, 1]
        target = apparatus.frame_columns(words, d, f)[0][:, 1]
        step = expm(-1j * h / 2)
        y = np.linalg.matrix_power(step, 256) @ s
        p = abs(target @ y) ** 2
        factors.append(p)
    expected = saved["settings"][0]["records"][4]
    assert abs(np.prod(factors) - expected["matching_probability"]) < 1e-10
    assert abs(expected["matching_count"] / measure.SHOTS - np.prod(factors)) < 0.001
    assert 0.905 < factors[0] < 0.907
    assert 0.921 < factors[1] < 0.923


def test_inventory_calibrations_and_all_time_gap():
    result = measure.run()
    archived = json.loads((HERE / "results.json").read_text())

    def check(a, b):
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            for k in a:
                check(a[k], b[k])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for v, u in zip(a, b, strict=True):
                check(v, u)
        elif isinstance(a, float):
            assert abs(a - b) < 2e-9
        else:
            assert a == b

    check(result, archived)
    for case in result["settings"]:
        for row in case["records"]:
            assert np.array(row["three_read_counts"]).sum() == measure.SHOTS
    cut = result["settings"][4]["records"][4]
    assert cut["matching_count"] == 0
    assert np.array(cut["three_read_counts"])[:, :, 1].sum() > 3_800_000
    raised = np.array(result["settings"][2]["records"][4]["three_read_counts"])
    assert raised[0, :2, :].sum() == 0
    assert np.array(result["no_seeded_pair"])[:, :2, :2].sum() == 0
    b = bound.certificate()
    assert 0.03378 < b["all_time_lower_charge_and_either_neutral_upper"] < 0.0338
    surplus = result["first_wing_filled_surplus_control"]["probabilities"]
    assert 0.408 < surplus[1] < 0.409


def test_inclusive_receiver_envelope():
    panel = bound.receiver_envelope()
    assert panel["neutral_upper_fraction"] == "4001/80"
    assert 0.02555 < panel["all_three_port_words_probability_upper"] < 0.025563
    assert panel["any_length_single_boundary_upper"] < 0.032679
    for n in (6, 8, 10):
        words, h = apparatus.sector(n, 1, 1, 8)
        row_bounds = np.diag(h) + np.sum(abs(h), axis=1) - abs(np.diag(h))
        for word, bound_at_word in zip(words, row_bounds, strict=True):
            touching = abs(word.index(0) - word.index(2)) == 1
            assert bound_at_word <= (51 if touching else 36)
    at_128 = json.loads((HERE / "results.json").read_text())[
        "raised_charge_all_receiver_words"
    ][4]
    assert sum(at_128["six_words_and_miss"][:6]) == 0
    assert 2.55e-8 < at_128["emission_probability"] < 2.58e-8
