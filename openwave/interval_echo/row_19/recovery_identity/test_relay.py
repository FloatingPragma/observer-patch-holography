from copy import deepcopy
from fractions import Fraction
import json

import pytest

from openwave.interval_echo.row_19.recovery_identity import relay_precision as repeat


def test_banks_and_prediction():
    archived = json.loads(repeat.RECEIPT.read_text())
    panels = repeat.decode(archived)
    assert panels == archived["results"]
    for method in ("release", "finite"):
        lo, hi = panels["even"][method]["held_over_ordinary"]
        assert lo < hi < 0 < archived["source_interval"][0]
        assert panels["even"][method]["compatible"] is False
        assert all(lo < hi < 0 for lo, hi in panels["even"][method]["formed"])
    weights = repeat.allowable_weights(panels, archived["source_interval"])
    assert 0.49 < weights["release"][0] < weights["release"][1] < 0.51
    assert 0.83 < weights["finite"][0] <= weights["finite"][1] == 1
    assert panels["vacancy"]["finite"]["compatible"]
    assert panels["excess"]["finite"]["compatible"] is False
    for q, sign in (("vacancy", 1), ("excess", -1)):
        for lo, hi in panels[q]["finite"]["formed"]:
            assert sign * lo > 0 and sign * hi > 0


def test_incomplete_or_swapped_counts_rejected():
    archived = json.loads(repeat.RECEIPT.read_text())
    broken = deepcopy(archived)
    broken["shell"][0][4][1][1][0] += 1
    with pytest.raises(ValueError, match="Incomplete"):
        repeat.decode(broken)
    broken = deepcopy(archived)
    broken["operators_times_8"][0].reverse()
    with pytest.raises(AssertionError):
        repeat.decode(broken)


def test_tail_sum_includes_every_count_and_handles_zero_calibration():
    n = repeat.SHOTS
    mean, half = repeat.shell_read([n // 2, 0, n // 4, n // 4], 12)
    assert mean == Fraction(5, 4)
    assert half > 0
    panel = repeat.intervals([1] * 8, [0] * 8, {"ordinary": 1, "held": 1}, (1, 2))
    assert panel["compatible"] is None


def test_three_contact_delivery_from_prior_two_head_solver():
    """Attach the third restriction to the previous full two-head solver."""
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import expm_multiply

    from openwave.interval_echo.row_19.recovery_identity import relay

    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    archive = json.loads(relay.RECEIPT.read_text())
    path = archive["stations"][0]["path"]
    edge = sorted(path[:2])
    stocks = relay.sources(board, route)[0]
    direct = np.zeros((2, 2), complex)
    # One held stock populates both source parity blocks independently.
    for sector in (0, 1):
        h2 = relay.RANGING.generator(board, route, 1, 8, bool(sector))
        h3 = sparse.kron(h2, sparse.eye(2), format="csr")
        local = relay.FORMATION.local_parts(1, 8)
        for i in range(2):
            for j in range(2):
                head = local[i + 1, j + 1].copy()
                if i == j:
                    head -= local[0, 0] + 8.5 * np.eye(2)
                else:
                    head += 0.5 * np.eye(2)
                address = sparse.csr_matrix(
                    ([1], ([edge[i]], [edge[j]])), shape=(board.n, board.n)
                )
                h3 += sparse.kron(address, np.kron(np.eye(4), head), format="csr")
        h3 -= h3.diagonal().mean() * sparse.eye(h3.shape[0], format="csr")
        general = relay.generator(
            board, [(route[:2], 1 - 2 * sector), (route[2:], 1), (edge, 1)], 1
        )
        np.testing.assert_allclose((general - h3).data, 0, atol=4e-12, rtol=0)
        weights, vectors = np.linalg.eigh(stocks[3, sector])
        start = np.zeros((8 * board.n, 4), complex)
        start[8 * route[0] + 2 * np.arange(4)] = vectors
        wave = expm_multiply(-8j * h3, start, traceA=0).reshape(-1, 2, 4)
        for i, weight in enumerate(weights):
            component = wave[:, :, i]
            direct += weight * component.T @ component.conj()
    target = archive["stations"][0]["appointments"]["8"]["bloch"][3]
    read = [np.trace(direct @ sigma).real for sigma in relay.FORMATION.PAULI]
    np.testing.assert_allclose(read, target, atol=4e-12, rtol=0)
