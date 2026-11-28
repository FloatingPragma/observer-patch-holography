"""Free propagation on a closed native board with its original written edge."""

import argparse
import json
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply
from model.executable_model import Carrier, federation
from . import source
from ..closed_region.measure import check, layout

TIMES = (0, 1, 4, 16, 64)


def operators(board):
    seam = min(board.external)
    a, b = seam
    remote = np.ones(board.n)
    remote[[a, b]] = 0
    law = source.LAW
    h = source.formation.emission.four_h(8)
    memories = [i + 1 for i, (_, d) in enumerate(law.LABELS) if abs(d) == 1]
    words = [(0, 0), (1, 0), (0, 1)]
    ix = np.array([m * 9 + law.SRC.index(word) for word in words for m in memories])
    assert not h[np.ix_(np.setdiff1d(np.arange(99), ix), ix)].any()
    local = h[np.ix_(ix, ix)] / 4
    assert not local[:4, 4:].any()
    written = sparse.kron(sparse.diags(remote), local[:4, :4], format="csr")
    loc = np.r_[4 * a : 4 * a + 4, 4 * b : 4 * b + 4]
    rows, cols = np.nonzero(local[4:, 4:])
    written += sparse.csr_matrix(
        (local[4 + rows, 4 + cols], (loc[rows], loc[cols])), shape=written.shape
    )
    edges = np.array([edge for edge in board.seams if edge != seam])
    left, right = edges.T
    adjacency = sparse.csr_matrix(
        (np.ones(2 * len(edges)), (np.r_[left, right], np.r_[right, left])),
        shape=(board.n, board.n),
    )
    degree = 8 * (6 - board.degree)
    blank = (
        sparse.diags(8.5 * np.bincount(edges.ravel(), minlength=board.n) + degree)
        - adjacency / 2
    )
    matrix = sparse.kron(blank, sparse.eye(4), format="csr") + written
    branches = np.flatnonzero(law.echo()[:, law.SRC.index((0, 1))])
    indices = [
        (memories.index(int(i // 9)), (a if law.SRC[i % 9][0] else b)) for i in branches
    ]
    starts = np.zeros((4 * board.n, 4))
    weights = []
    for j, (memory, _) in enumerate(indices):
        for k, (_, port) in enumerate(indices):
            starts[4 * port + memory, 2 * j + k] = 1
            weights.append(-0.25 if j == k else 0.25)
    dist = (
        shortest_path(
            board.adjacency(), directed=False, unweighted=True, indices=list(seam)
        )
        .min(axis=0)
        .astype(int)
    )
    return matrix, written, edges, degree, seam, dist, starts, np.array(weights)


def read(states, weights, parts):
    mat, written, edges, degree, seam, dist, _, _ = parts
    # Columns carry signed ensemble weights, including degenerate energy blocks.
    psi = states.reshape(len(dist), 4, -1)
    occupancy = (abs(psi) ** 2).sum(axis=1) @ weights
    left, right = edges.T
    interference = np.einsum(
        "imk,imk,k->i", psi[left].conj(), psi[right], weights, optimize=True
    ).real
    en = 8.5 * (occupancy[left] + occupancy[right]) - interference
    per = degree * occupancy
    for endpoint in (left, right):
        np.add.at(per, endpoint, en / 2)
    written_energy = float(
        (np.sum(states.conj() * (written @ states), axis=0).real) @ weights
    )
    per[list(seam)] += written_energy / 2
    total = float((np.sum(states.conj() * (mat @ states), axis=0).real) @ weights)
    np.testing.assert_allclose(per.sum(), total, atol=3e-10)
    np.testing.assert_allclose([total, occupancy.sum()], [17 / 8, 0], atol=3e-9)
    return dict(
        written_energy=written_energy,
        shell_energy=np.bincount(dist, weights=per).tolist(),
        enclosed_energy=np.cumsum(np.bincount(dist, weights=per)).tolist(),
        enclosed_load=np.cumsum(np.bincount(dist, weights=occupancy)).tolist(),
    )


def aperture8(parts, radius):
    """Exact integer matrix for eight times the endpoint-shared energy."""
    h, written, edges, degree, seam, dist, _, _ = parts
    mask = (dist <= radius).astype(np.int64)
    left, right = edges.T
    shares = mask[left] + mask[right]
    diagonal = 8 * mask * degree
    for endpoint in (left, right):
        np.add.at(diagonal, endpoint, 34 * shares)
    exchange = sparse.csr_matrix(
        (-2 * np.tile(shares, 2), (np.r_[left, right], np.r_[right, left])),
        shape=(len(mask), len(mask)),
    )
    op = sparse.kron(
        sparse.diags(diagonal, dtype=np.int64) + exchange, sparse.eye(4), format="csr"
    )
    op += 4 * sum(mask[list(seam)]) * written
    integer = op.astype(np.int64)
    assert (integer != op).nnz == 0
    if mask.all():
        assert (integer != 8 * h).nnz == 0
    return integer


def curvature(parts):
    h, _, _, _, _, dist, starts, weights = parts
    scaled = 4 * h
    four = scaled.astype(np.int64)
    assert (four != scaled).nnz == 0
    base = starts.astype(np.int64)
    one = four @ base
    two = four @ one
    signs = (4 * weights).astype(np.int64)
    values = []
    for radius in range(max(dist) + 1):
        a = aperture8(parts, radius)
        numerator = ((one * (a @ one) - two * (a @ base)).sum(axis=0)) @ signs
        values.append(Fraction(2 * int(numerator), 4 * 4**2 * 8))
    return [str(x) for x in values]


def run_case(board):
    parts = operators(board)
    h, _, _, _, seam, dist, starts, weights = parts
    states = starts.copy()
    dt = 0
    center = float(h.diagonal().mean())
    generator = -1j * (h - center * sparse.eye(h.shape[0], format="csr"))
    rows = {}
    for t in TIMES:
        if t > dt:
            states = expm_multiply((t - dt) * generator, states, traceA=0.0)
        rows[str(t)] = read(states, weights, parts)
        dt = t
    if board.cells == 2:
        es, vs = np.linalg.eigh(h.toarray())
        for tol in (1e-8, 1e-11):
            groups = np.split(np.arange(len(es)), np.where(np.diff(es) > tol)[0] + 1)
            projected = np.hstack([vs[:, ix] @ (vs[:, ix].T @ starts) for ix in groups])
            mean = read(projected, np.tile(weights, len(groups)), parts)
            if "infinite_mean" in rows:
                check(rows["infinite_mean"], mean)
            rows["infinite_mean"] = mean
    return dict(
        owners=board.cells,
        dimension=h.shape[0],
        seam=seam,
        ports_at_each_radius=np.bincount(dist).tolist(),
        enclosed_energy_second_derivative_at_zero=curvature(parts),
        times=rows,
    )


def run():
    _, edges, _, _ = layout(2)
    boards = (
        Carrier("closed two-owner assembly", 24, edges),
        federation(0),
        federation(1),
    )
    return [run_case(board) for board in boards]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = json.loads(json.dumps(run()))
    path = source.HERE / "free.json"
    if p.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
