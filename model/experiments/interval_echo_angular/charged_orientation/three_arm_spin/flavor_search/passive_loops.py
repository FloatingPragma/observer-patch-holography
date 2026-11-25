"""Directed single-lap transfer of a unit defect with every native head allowed."""

from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, expm
from scipy.sparse import csr_matrix, kron

from ...flux_doublet.measure import MOD, face_edges, faces, links
from ...live_trimer.measure import LAW, blocks
from ..apparatus import fixture

LABELS = ("blank", (0, 1), (1, -1), (1, 1), (2, -1), (2, 1))
TIMES = (0, 0.5, np.pi / 9, np.pi / 8, np.pi)


def local_transfers(sign, forward):
    """Check the invariant components of 8H exactly and the projected exp(-itH)."""
    endpoints = (1 + sign, 1) if forward else (1, 1 + sign)
    index = LAW.SRC.index(endpoints)
    other = LAW.SRC.index(endpoints[::-1])
    enabled = int((sign == 1) == forward)
    hops = []
    for label in LABELS:
        h8 = blocks()[label]
        source = np.arange(index, len(h8), 9)
        dest = np.arange(other, len(h8), 9)
        expected = np.zeros((len(source), len(source)), dtype=int)
        if label == "blank":
            expected[0, 0] = -4
        elif label[0] == 1:
            expected[1 - enabled, enabled] = -32 if label[1] == 1 else -36
        np.testing.assert_array_equal(h8[np.ix_(dest, source)], expected)
        # Each input column either stays at its original word or participates
        # in one isolated two-state exchange, so support holds for all t.
        for j in range(len(source)):
            initial = source[j]
            allowed = [initial] + dest[np.flatnonzero(expected[:, j])].tolist()
            outside = np.setdiff1d(np.arange(len(h8)), allowed)
            assert not h8[np.ix_(outside, allowed)].any()
        for time in TIMES:
            unitary = expm(-1j * time * h8 / 8)
            transfer = unitary[np.ix_(dest, source)]
            assert np.max(abs(transfer[expected == 0]), initial=0.0) < 1e-13
            for k, j in zip(*np.nonzero(expected), strict=True):
                initial, final = source[j], dest[k]
                assert h8[initial, initial] == h8[final, final]
                # The isolated symmetric dimer gives the transfer analytically.
                phase = np.exp(-1j * time * h8[initial, initial] / 8)
                dim_amp = -1j * phase * np.sin(time * expected[k, j] / 8)
                np.testing.assert_allclose(transfer[k, j], dim_amp, atol=9e-14)
        hops.append(expected)
    head = block_diag(*hops)
    projected = np.zeros((11, 11), dtype=int)
    projected[0, 0] = 1
    np.testing.assert_array_equal(head @ head, 16 * projected)
    return head


def check_triangle(hops):
    """Multiply sparse integer matrices on all 11^3 head words of a receiver face."""
    board = fixture()
    g, y, z = board["legs"][0][-3:]
    path = ((g, y), (y, z), (z, g))
    assert path[0] in board["blank"] and path[1] in board["blank"]
    assert path[2] in board["marked"]
    rows = []
    for sign in (-1, 1):
        # Match each local head to its listed contact direction.
        w = csr_matrix([[1]], dtype=np.int64)
        for a, b in path:
            forward = (a, b) in board["marked"] or (a, b) in board["blank"]
            w = kron(w, csr_matrix(hops[sign, forward]), format="csr")
        twice = w @ w
        assert twice.nnz == 1 and twice[0, 0] == (-4) ** 6
        assert w[0].nnz == w[:, [0]].nnz == 1 and w[0, 0] == (-4) ** 3
        # Receiver preparation restricts the first two heads to blank and
        # the third to its live magnitude-one even-parity sector (indices 5,6).
        even = w[[5, 6], :][:, [5, 6]]
        assert even.nnz == 1 and (even @ even).nnz == 0
        rows.append(
            dict(
                sign=sign,
                one_lap_nonzero=w.nnz,
                two_laps_nonzero=twice.nnz,
                prepared_one_lap_rank=1,
                prepared_two_laps_rank=0,
            )
        )
    return dict(
        ports=[g, y, z], full_head_dimension=11**3, rows=rows, integer_scale_per_step=8
    )


def supplied_phase_control():
    """A separately imposed U(1) connection gives a detectable blank-loop phase."""
    native_blank = blocks()["blank"] / 8
    at_pi = expm(-1j * np.pi * native_blank)
    b = at_pi[LAW.SRC.index((1, 0)), LAW.SRC.index((0, 1))]
    np.testing.assert_allclose(b, 1, atol=1e-13)
    rows = []
    for q in (0, 1, 10):
        phase_table = links(q)
        for face in faces():
            steps = [phase_table[s, d] for s, d in face_edges(face)]
            rotated = np.prod([b * np.exp(2j * np.pi * step / MOD) for step in steps])
            relative = rotated / b**3
            np.testing.assert_allclose(
                relative, np.exp(2j * np.pi * q / MOD), atol=1e-13
            )
        angle = 2 * np.pi * q / MOD
        rows.append(
            dict(
                supplied_flux_units=q,
                faces_verified=len(faces()),
                bright_x=round((1 + np.cos(angle)) / 2, 12),
                bright_y=round((1 + np.sin(angle)) / 2, 12),
            )
        )
    return rows


def run():
    directions = list(product((-1, 1), (False, True)))
    hops = {(s, f): local_transfers(s, f) for s, f in directions}
    for s in (-1, 1):
        np.testing.assert_array_equal(hops[s, False], hops[s, True].T)
    return dict(
        local_dimension=sum(len(blocks()[sector]) for sector in LABELS),
        appointments=[float(t) for t in TIMES],
        oriented_transitions=[
            dict(
                sign=s,
                forward=f,
                minus_eight_hop=(-hops[s, f][hops[s, f] != 0]).tolist(),
            )
            for s, f in directions
        ],
        receiver_triangle=check_triangle(hops),
        supplied_flux_control=supplied_phase_control(),
    )


if __name__ == "__main__":
    Path(__file__).with_name("passive_loops.json").write_text(
        json.dumps(run(), indent=2) + "\n"
    )
