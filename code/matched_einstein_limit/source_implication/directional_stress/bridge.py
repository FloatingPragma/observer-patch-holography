"""Energy flux into the receiving cell at the four-tick screen appointment."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..action_response import measure as action
from ..beam_response import measure as beam

HERE = Path(__file__).resolve().parent
native = action.native


def operators(board, route, charge, odd):
    """Twice the internal cell actions, each on its own 48-dimensional block."""
    eye = np.eye(4, dtype=np.int64)
    coin = np.eye(2, dtype=np.int64)
    cells = np.zeros((board.cells, 12, 4, 12, 4), dtype=np.int64)
    marked = [tuple(sorted(pair)) for pair in (route[:2], route[2:])]
    assert 16 * abs(2 * charge + charge**2) * max(6 - board.degree) <= 2 * 64
    for p, degree in enumerate(board.degree):
        i = p % 12
        cells[p // 12, i, :, i] += 16 * (2 * charge + charge**2) * (6 - degree) * eye
    for a, b in board.seams:
        if a // 12 != b // 12:
            continue
        h = cells[a // 12]
        i, j = a % 12, b % 12
        if (a, b) in marked:
            slot = marked.index((a, b))
            blocks = 2 * native.iso.ranging.parity_blocks(
                charge, 8, -1 if odd and slot == 0 else 1
            )
            np.testing.assert_array_equal(blocks, blocks.astype(np.int64))
            blocks = blocks.astype(np.int64)
            np.testing.assert_array_equal(blocks[0, 0], blocks[0, 0, 0, 0] * coin)
            contact = np.block(
                [
                    [
                        blocks[i + 1, j + 1] - int(i == j) * blocks[0, 0]
                        for j in range(2)
                    ]
                    for i in range(2)
                ]
            )
            assert np.max(np.abs(contact).sum(axis=1)) <= 2 * 64
            for u, dst in enumerate((i, j)):
                for v, src in enumerate((i, j)):
                    part = blocks[1 + u, 1 + v] - int(u == v) * blocks[0, 0]
                    h[dst, :, src] += (
                        np.kron(part, coin) if slot == 0 else np.kron(coin, part)
                    )
        else:
            h[i, :, i] += 17 * eye
            h[j, :, j] += 17 * eye
            h[i, :, j] -= eye
            h[j, :, i] -= eye
    return cells.reshape(board.cells, 48, 48)


def bridges(board, edges):
    aa, bb = np.asarray(edges).T
    n = len(aa)
    values = np.repeat([17, 17, -1, -1], n)
    x, y = np.r_[aa, bb, aa, bb], np.r_[aa, bb, bb, aa]
    matrix = sparse.csr_matrix((values, (x, y)), shape=(board.n, board.n))
    return sparse.kron(matrix, sparse.eye(4, dtype=np.int64), format="csr")


def masked(cells, owners):
    allowed = set(owners)
    return sparse.block_diag(
        [
            sparse.csr_matrix(block)
            if i in allowed
            else sparse.csr_matrix((48, 48), dtype=np.int64)
            for i, block in enumerate(cells)
        ],
        format="csr",
    )


def current(board, cells, edge):
    a, b = edge
    diff = masked(cells, [b // 12]) - masked(cells, [a // 12])
    bridge = bridges(board, [edge])
    return bridge @ diff - diff @ bridge  # physical current = i times this / 8


def audit(board, cells, route, generator):
    complete = masked(cells, range(board.cells)) + bridges(board, board.external)
    delta = complete - generator
    np.testing.assert_array_equal(delta.diagonal(), delta.diagonal()[0])
    assert (delta - sparse.diags(delta.diagonal(), dtype=np.int64)).count_nonzero() == 0
    receiving = route[2] // 12
    incident = [
        (a, b) if b // 12 == receiving else (b, a)
        for a, b in board.external
        if receiving in (a // 12, b // 12)
    ]
    # [2H,4E_receiver] is 8 times the physical time derivative divided by i.
    region4 = 2 * masked(cells, [receiving]) + bridges(board, incident)
    total = sum(current(board, cells, e) for e in incident)
    assert (complete @ region4 - region4 @ complete - total).count_nonzero() == 0
    return [list(e) for e in incident]


def tangents():
    states = action.preparations()
    rows = json.loads((HERE.parent / "signed_transfers.json").read_text())["tests"]
    rotation = (states["plus"] - states["minus"]) / 2
    out = {"reference": states["equilibrium"]}
    for row in rows:
        name = "pairing" if row["source"] == "collar" else "ordinary"
        out[name] = (
            states[name] - states["equilibrium"] + row["angle_per_unit"] * rotation
        )
    return out


def evaluate(effects, states):
    return {
        k: float(np.einsum("bij,bji->", effects, v).real) for k, v in states.items()
    }


def report():
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    rings = beam.screens(board, route)
    states = tangents()
    nodes = list(route[1:3]) + list(rings[0]) + list(rings[1])
    rows = {}
    for q in (-1, 1):
        ensembles = [[], [], []]
        for odd in (False, True):
            h = native.generator(board, route, q, odd)
            cells = operators(board, route, q, odd)
            inc = audit(board, cells, route, h)
            full = current(board, cells, route[1:3])
            support = beam.indices(nodes)
            k = full[support, :][:, support].toarray()
            assert np.count_nonzero(k) == full.count_nonzero()
            pos = beam.indices(list(rings[0]) + list(rings[1]))
            w = np.zeros((h.shape[0], len(pos)))
            w[pos] = np.eye(len(pos))
            ht = -2j * h
            w = expm_multiply(ht, w, traceA=ht.diagonal().sum())[support]
            joint = 1j * k / 8
            b = bridges(board, [route[1:3]])[support, :][:, support].toarray() / 2
            owners = [p // 12 for p in route[1:3]]
            diff = masked(cells, [owners[1]]) - masked(cells, [owners[0]])
            d = diff[support, :][:, support].toarray() / 2
            # Exact Fourier degree <=9 for the isolated bridge. Check its
            # independent phase-kick read recovers the full current operator.
            energy, basis = np.linalg.eigh(b)
            np.testing.assert_array_equal(np.rint(energy), np.round(energy, 7))
            assert set(np.rint(energy)) == {0, 8, 9}
            read = np.zeros_like(joint)
            for m in range(19):
                theta = 2 * np.pi * m / 19
                u = (basis * np.exp(-1j * theta * energy)) @ basis.T
                freq = np.arange(-9, 10)
                weight = np.real(np.sum(1j * freq * np.exp(-1j * theta * freq))) / 19
                read += (weight / 2) * (u.conj().T @ d @ u)
            np.testing.assert_allclose(read, joint, atol=1e-12, rtol=0)
            head = np.array([v.conj().T @ joint @ v for v in np.hsplit(w, 10)])
            for j, subset in enumerate((head[:5], head[5:], head)):
                ensembles[j].append(subset.mean(axis=0))
        groups = {
            k: evaluate(np.array(e), states)
            for k, e in zip(
                ("forward_launch", "reverse_launch", "equal_tagged_launch"), ensembles
            )
        }
        rows[str(q)] = groups
    return dict(
        time=4,
        forward_ports=rings[0],
        reverse_ports=rings[1],
        oriented_edge=list(route[1:3]),
        receiving_edges=inc,
        charges=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "bridge.json").write_text(data)
    else:
        print(data, end="")
