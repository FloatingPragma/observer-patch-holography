import json
from hashlib import sha256
from math import comb

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.executable_model import Carrier, cell, federation
from model.experiments.interval_echo_consensus.sign_graph import PrivateRecordWalk
from . import candidate
from .probe import geometry, pair_source
from .release import HERE, certified_window, check_receipt, fourth_derivative_bound, joint_read, word_transform


def test_probe_and_source_on_every_board():
    for board in (cell(), *(federation(level) for level in range(3))):
        gram, contacts, rest, bits = geometry(board)
        assert set(gram.data) == {8}
        assert set(contacts) == set(board.seams)
        assert bits >= len(board.seams)*8
        actual = candidate.generator(board) - board.laplacian()/2
        np.testing.assert_array_equal(actual.diagonal(), rest)
    # Both sharp pointer outcomes of one interrogation at the exchange pi exposure.
    x = np.array([[0, 1], [1, 0]])
    for active in (0, 1):
        u = expm(-1j*np.pi*active*(np.eye(2)-x)/2)
        np.testing.assert_allclose(u, x if active else np.eye(2), atol=1e-14)


def test_every_sector_and_joint_read_against_the_full_graph():
    board = Carrier("cycle", 4, [(0, 1), (0, 3), (1, 2), (2, 3)])
    edge = (0, 1)
    walk = PrivateRecordWalk(board, "primitive")
    pairs, odd, _, active, width, _, bits = pair_source(board, edge)
    size = 2**bits
    tags = [j for j, (a, b, _) in enumerate(walk.offers) if (a, b) == edge]
    assert len(tags) == width == 2 and size == 256
    mask = sum(1 << j for j in tags)
    positions = list(zip(map(int, pairs.aa), map(int, pairs.bb), strict=True))
    rows, cols, vals = [], [], []
    for j, position in enumerate(positions):
        for bank in range(size):
            origin = j*size+bank
            neighbors = list(walk.neighbors(position, bank))
            rows.append(origin)
            cols.append(origin)
            vals.append(len(neighbors)/2)
            for dest, final in neighbors:
                rows.append(pairs.index(*dest)*size+final)
                cols.append(origin)
                vals.append(-.5)
    literal = sparse.csr_matrix((vals, (rows, cols)), shape=(size*pairs.P,)*2)
    start = np.zeros((pairs.P, size))
    for bank in range(size):
        if not bank & mask:
            start[pairs.index(*edge), bank] = (-1)**bank.bit_count()/np.sqrt(size/2**width)
    home = [pairs.index(*seam) for seam in board.seams]
    psi = np.eye(pairs.P)[pairs.index(*edge)]
    for t in (.7, 7):
        full = abs(expm_multiply(-1j*t*literal, start.ravel()).reshape(start.shape))**2
        groups = np.zeros((width+1, pairs.P))
        for bank in range(size):
            groups[(bank & mask).bit_count()] += full[:, bank]
        expected = np.column_stack((groups[:, home], np.delete(groups, home, axis=1).sum(axis=1)))
        states = np.array([expm(-1j*t*(odd.toarray()-np.diag((width-m)*active))) @ psi for m in range(width+1)])
        np.testing.assert_allclose(joint_read(states, home), expected, atol=1e-12, rtol=0)
        states[-1] *= -1
        assert np.max(abs(joint_read(states, home)-expected)) > .01
    changed = np.array([bool(bank & mask) for bank in range(size)])
    first = literal @ start.ravel()
    second = (literal @ first).reshape(start.shape)
    assert np.linalg.norm(first.reshape(start.shape)[:, changed]) == 0
    assert np.linalg.norm(second[:, changed])**2/4 == pytest.approx(2*width/64)
    # In the union rule, both occupied ends enable precisely the extra tagged
    # flips on the launch. Apply these columns directly to the sharp input.
    rival_first = first.reshape(start.shape).copy()
    launch = pairs.index(*edge)
    for j in tags:
        rival_first[launch] += (start[launch] - start[launch, np.arange(size) ^ (1 << j)])/2
    assert np.linalg.norm(rival_first[:, changed])**2 == pytest.approx(width/4)


def test_archived_joint_marginals_and_historical_source():
    from model.experiments.interval_echo_consensus.record_release import source
    from model.experiments.interval_echo_consensus.bounded_clock_flight import DetectorEvolution

    result = json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    spec = json.loads((HERE/"protocol.json").read_text(encoding="utf-8"))
    board = federation(spec["board_level"])
    pairs, odd, common, active, width, _, _ = pair_source(board, tuple(spec["launch"]))
    ref, launch, home, _ = source("primitive")
    assert not (odd != ref).nnz and launch == pairs.index(*spec["launch"])
    assert not (odd != pairs.hamiltonian(8)).nnz
    archived = (HERE.parents[2]/"openwave/interval_echo/row_04/row_04_receipt.json").read_bytes()
    assert sha256(archived).hexdigest() == "7f6e458fd0b7cf9b9e50fb15a5c4b77f944a2457a9c53cdcca88efdf4d4c3bdd"
    historical = json.loads(archived)
    assert historical["source"] == result["launch"] and historical["addresses"] == [*result["seams"], "quiet"]
    for arm in historical["arms"]:
        assert arm["strength"] in (0, 8)
        for j, (t, record) in enumerate(zip(spec["exposures"], arm["readings"], strict=True)):
            actual = (result["appointments"][str(t)]["odd_home_and_quiet"] if arm["strength"] == 8
                      else result["controls"]["all_even"]["home_and_quiet"][j])
            assert record["exposure"] == t
            np.testing.assert_allclose(actual, record["probabilities"], rtol=0, atol=5.1e-10)
    assert result["first_changed_word_quartic"] == 9*width/64
    times = spec["exposures"]
    probs = []
    for m in range(width+1):
        h = odd-sparse.diags((width-m)*active, dtype=float)
        probs.append(DetectorEvolution(h.tocsr(), times).read(launch, home))
    weights = np.array([comb(width, m)/2**width for m in range(width+1)])
    for j, t in enumerate(times):
        row = result["appointments"][str(t)]
        joint = np.array(row["weight_and_home_or_quiet"])
        assert joint.shape == (width+1, len(home)+1) and np.min(joint) >= 0
        np.testing.assert_allclose(joint.sum(axis=0), weights @ np.array(probs)[:, j],
                                   atol=1e-8, rtol=0)
        assert probs[-1][j, :-1].sum() == pytest.approx(row["odd_home"], abs=1e-9)
        assert joint[1:].sum() == pytest.approx(row["changed_word"], abs=1e-9)
    assert common[launch] == width and active[launch] == 0
    assert result["window"]["odd"]["gate"] == "certified"
    windows = result["window"]["by_odd_count"]
    assert len(windows) == width+1 and all(w["gate"] == "certified" for w in windows)
    assert result["window"]["all_tagged_preparations_lower"] == min(w["continuous_lower"] for w in windows)
    assert result["window"]["all_tagged_preparations_lower"] > .94099
    transform, counts = word_transform(width)
    np.testing.assert_array_equal(transform[0], counts/2**width)


def test_direct_sharp_word_readout_first_appointment():
    from .check_words import check
    for times in ([], [7], [1, 1], [True]):
        with pytest.raises(ValueError):
            check(times)
    assert check([1])[1] < 1e-9


def test_certified_interpolation_interior_minimum_and_receipt_rejections():
    h = sparse.csr_matrix([[.5, -.5], [-.5, .5]])
    bound = fourth_derivative_bound(h, 0)
    # p(t) = (1 + cos(t))/2 has a minimum strictly between two supplied nodes.
    times = np.arange(14)/4
    p = (1+np.cos(times))/2
    window = certified_window(p, -.5*np.sin(times), bound, 4, 1e-12)
    assert p.min() > 0 and window["continuous_lower"] == 0 and window["gate"] == "refuted"
    for value in (np.uint64(2**63+1), 2+2j, float("inf")):
        with pytest.raises(ValueError):
            fourth_derivative_bound(sparse.csr_matrix([[value]]), 0)
    with pytest.raises(ValueError):
        certified_window(p+2, p, bound, 4, 1e-12)
    for home in ([], [-1], [0, 0], [2]):
        with pytest.raises(ValueError):
            joint_read(np.ones((1, 2)), home)
    for wrong in ({"p": [float("nan")]}, {"p": [.5, 0.]}, {"p": [False]}, {"p": [.5], "extra": 1}):
        with pytest.raises(ValueError):
            check_receipt({"p": [.5]}, wrong)
