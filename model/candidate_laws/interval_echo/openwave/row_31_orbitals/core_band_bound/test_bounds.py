from itertools import pairwise

import numpy as np

from . import measure


def test_native_empty_record_and_uniform_gap():
    assert len(measure.native_columns()) == 5
    for a, threshold in ((1, -3), (2, -27)):
        bands = measure.bands(range(7), a, 6)
        assert bands[-2][-1] == threshold
        assert all(u[3] < v[2] for u, v in pairwise(bands))
        assert bands[0][2] == -48 * a * a


def test_full_native_term_compression_for_overlapping_core():
    law = measure.candidate.native_law
    a0, a1 = law.local_term_parts()
    h = a0 + 8 * a1
    record = (law.echo()[:, 0]).reshape(11, 9)[:, 0]
    board = measure.cell()
    written = {tuple(board.seams[0]), tuple(board.seams[1]), tuple(board.seams[-1])}
    for amount in (1, 2):
        matrix = np.diag((6 - board.degree) * 8.0 * amount**2)
        for x, y in board.seams:
            anc = record if (x, y) in written else np.eye(11)[0]
            src = np.eye(9)[:, [law.SRC.index((amount, 0)), law.SRC.index((0, amount))]]
            frame = np.kron(anc[:, None], src)
            native = frame.T @ h @ frame
            np.testing.assert_allclose(h @ frame, frame @ native, atol=2e-13)
            matrix[np.ix_([x, y], [x, y])] += native
        e = np.linalg.eigvalsh(matrix - 48 * amount**2 * np.eye(board.n))
        offset = 0
        for _, n, lo, hi in measure.audit_cell(board, written, amount):
            assert e[offset] >= lo - 1e-12 and e[offset + n - 1] <= hi + 1e-12
            offset += n
