"""Finite coherent motion in the two-beacon coordinate supplied by exterior ranging."""

import numpy as np
from scipy.linalg import eigh


def doubled_generator(board, active, raw, load):
    # Pulse only active native bonds. The physical private contacts stay at
    # their original sites when the observer gates off the other bonds.
    h = raw[:9, :9]
    pair = np.array([3 * load, load])
    np.testing.assert_array_equal(
        h[np.ix_(pair, pair)], [[32 * load**2 + 2, -2], [-2, 32 * load**2 + 2]]
    )
    adj = board.adjacency()[active][:, active].toarray().astype(np.int64)
    degree = adj.sum(1)
    return (
        np.diag(
            (16 * load**2 + 1) * degree
            + 16 * load**2 * (6 - board.degree[active])
            - 96 * load**2
        )
        - adj
    )


def derivatives(J):
    J2 = J @ J
    J3 = J2 @ J
    J4 = J3 @ J
    # Derivatives of |exp(-itJ/2)|^2. On the diagonal, products involving
    # the identity select the initial address. All numerators are integers.
    second = 2 * (J * J - np.diag(np.diag(J2)))
    fourth = 2 * np.diag(np.diag(J4)) - 8 * J * J3 + 6 * J2 * J2
    assert np.all(second.sum(1) == 0) and np.all(fourth.sum(1) == 0)
    return second, fourth


def compare(src, ordinary, err, oe, ids, ports):
    # Signed changes outside their error bands test every nonnegative
    # common conversion on precisely the same release and two-beacon chart.
    opposite = ((src > err) & (ordinary < -oe)) | ((src < -err) & (ordinary > oe))
    pick = (
        int(np.argmax(np.minimum(abs(src) - err, abs(ordinary) - oe) * opposite))
        if opposite.any()
        else int(np.argmax(abs(src)))
    )
    return {
        "range": [float(src.min()), float(src.max())],
        "ordinary_range": [float(ordinary.min()), float(ordinary.max())],
        "source_error_upper": err,
        "ordinary_error_upper": oe,
        "opposite_sign_count": int(opposite.sum()),
        "witness_port": int(ports[ids[pick]]),
        "witness_change": float(src[pick]),
        "witness_ordinary": float(ordinary[pick]),
    }


def run(board, ports, clearance, ab, field, ordinary, source_err, ordinary_err, raw):
    active = np.flatnonzero(clearance > 0)
    physical = ports[active]
    from scipy.sparse.csgraph import shortest_path

    d = shortest_path(board.adjacency(), unweighted=True)
    far = (clearance[active] >= 5) & (d[np.ix_(physical, ports[ab])].min(1) >= 5)
    ids = active[far]
    src, cal = field[active], ordinary[active]
    fullJ = doubled_generator(board, np.arange(board.n), raw, 1)
    full2, full4 = derivatives(fullJ)
    summary = {}
    for load in (1, 2):
        J = doubled_generator(board, physical, raw, load)
        d2, d4 = derivatives(J)
        for local, full in ((d2, full2), (d4, full4)):
            outside = np.setdiff1d(np.arange(board.n), physical)
            assert np.all(full[np.ix_(physical[far], outside)] == 0)
            np.testing.assert_array_equal(
                local[far], full[np.ix_(physical[far], physical)]
            )
        assert np.max(abs(d2[far] @ src) / 4) < 1e-9
        k4 = d4[far] / 16
        width = abs(k4).sum(1).max()
        out = {
            "fourth_derivative": compare(
                k4 @ src, k4 @ cal, width * source_err, width * ordinary_err, ids, ports
            )
        }
        energy, V = eigh(J / 2)
        for t in (0.25, 1.0, 2.0, 4.0):
            unitary = (V[far] * np.exp(-1j * t * energy)) @ V.T
            K = abs(unitary) ** 2
            assert np.max(abs(K.sum(1) - 1)) < 1e-12
            K -= np.eye(len(active))[far]
            x, y = K @ src, K @ cal
            width = abs(K).sum(1).max()
            out[str(t)] = compare(
                x,
                y,
                width * source_err + 1e-12,
                width * ordinary_err + 1e-12,
                ids,
                ports,
            )
        summary[str(load)] = out
    return {
        "clean_ports": physical.tolist(),
        "launches": ports[ids].tolist(),
        "loads": summary,
    }
