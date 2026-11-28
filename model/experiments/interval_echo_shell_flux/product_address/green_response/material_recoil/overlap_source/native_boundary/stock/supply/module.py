"""Direct contraction of the native 99-state operators on six physical ports."""

from collections import defaultdict

import numpy as np

from .. import measure

EDGES = [(2, 3), (2, 1), (3, 4), (1, 0), (4, 5)]


def apply(state, gate, edge, edges=EDGES, ports=6):
    a, b = edges[edge]
    result = defaultdict(complex)
    for word, amp in state.items():
        xa, xb, old = word[a], word[b], word[ports + edge]
        col = 9 * old + 3 * xa + xb
        for row in np.flatnonzero(abs(gate[:, col]) > 1e-12):
            new = list(word)
            new[a], new[b] = measure.LAW.SRC[row % 9]
            new[ports + edge] = int(row // 9)
            result[tuple(new)] += gate[row, col] * amp
    return {word: v for word, v in result.items() if abs(v) > 1e-12}


def inner(left, right):
    return sum(v.conjugate() * right.get(k, 0) for k, v in left.items())


def native_microscope():
    _, _, h4, _ = measure.kernels()
    h = h4 / 4
    echo = measure.LAW.echo()
    energy, vectors = np.linalg.eigh(h)
    swap = (vectors * np.exp(-1j * np.pi * energy)) @ vectors.T
    np.testing.assert_allclose(swap[:9, :9], measure.LAW._load_ops()[0], atol=1e-13)
    # Destructive public readout: undo capture, phase one donor endpoint
    # at its empty degree-five wing's ordinary reference, apply a half
    # exchange, and count public ordinary occupations on donor endpoints.
    one = [3, 1]
    half = (vectors * np.exp(-1j * np.pi * energy / 2)) @ vectors.T
    analyzer = half[np.ix_(one, one)] @ np.diag([-1j, 1])
    for ix, expected in [
        (0, np.array([[1, -1], [-1, 1]]) / 2),
        (1, np.ones((2, 2)) / 2),
    ]:
        row = analyzer[ix, :]
        np.testing.assert_allclose(np.outer(row.conj(), row), expected, atol=1e-13)
    _, hcontrast = measure.LAW.local_term_parts()
    np.testing.assert_allclose(
        hcontrast[np.ix_([0, 3], [0, 3])], np.diag([0, 1]), atol=1e-13
    )
    magnitudes = measure.MAG[np.arange(99) // 9]
    assert not np.any(h4[magnitudes[:, None] != magnitudes[None, :]])
    assert not np.any(h4[9:, :9])
    zero = [
        9 * (1 + j) + s
        for j, (_, d) in enumerate(measure.LAW.LABELS)
        if d == 0
        for s in range(9)
    ]
    assert not np.any(h4[:, zero])
    # Integer certificate for the separated-head energy bound. On empty
    # endpoints, (4h - 32I)/4 is a projection in the magnitude-one code.
    ix = [9 * (i + 1) for i, (_, d) in enumerate(measure.LAW.LABELS) if abs(d) == 1]
    gram = h4[np.ix_(ix, ix)] - 32 * np.eye(4, dtype=int)
    assert not h4[np.ix_(ix, np.setdiff1d(np.arange(99), ix))].any()
    np.testing.assert_array_equal(gram @ gram, 4 * gram)
    np.testing.assert_array_equal(gram, gram.T)

    def cold(port):
        word = [0] * 11
        if port is not None:
            word[port] = 1
        return {tuple(word): 1}

    def run(state):
        rows = []

        def take():
            np.testing.assert_allclose(inner(state, state), 1, atol=1e-13)
            rows.append(
                [inner(state, apply(state, h, i)).real for i in range(len(EDGES))]
            )

        take()
        for gates in [
            [(echo, 0)],
            [(swap, 1), (swap, 2)],
            [(swap, 3), (swap, 4)],
            [(echo, 1), (echo, 2)],
        ]:
            for gate, edge in gates:
                state = apply(state, gate, edge)
            take()
        return rows, state

    # For a closed single module the pre-write path with all five heads
    # blank is 2H=17D-A in the one-ordinary-unit sector. A tridiagonal
    # determinant recurrence uses integers, without choosing eigenvalue bins.
    dets = []
    for e in [16, 17, 18]:
        prior, det = 0, 1
        for degree in [1, 2, 2, 2, 2, 1]:
            prior, det = det, (17 * degree - 2 * e) * det - prior
        dets.append(det)
    assert dets == [1248, 288, 1656]
    nominal = {}
    states = []
    for name, state in [("vacuum", cold(None)), ("left", cold(2)), ("right", cold(3))]:
        es, s = run(state)
        states.append(s)
        nominal[name] = np.round(np.sum(es, axis=1), 12).tolist()
        for word in s:
            assert not any(word[1:5])
            assert measure.MAG[word[6]] == int(name != "vacuum")
            assert sum(word[:6]) == int(name != "vacuum")
    # Includes all diagonal and off-diagonal matrix elements: isometry and
    # energy spectral projectors on vacuum, one-left and one-right inputs.
    gram = np.array([[inner(a, b) for b in states] for a in states])
    reduced = np.array([[inner(a, apply(b, h, 0)) for b in states] for a in states])
    np.testing.assert_allclose(gram, np.eye(3), atol=1e-13)
    np.testing.assert_allclose(reduced, h[np.ix_([0, 3, 1], [0, 3, 1])], atol=1e-13)
    for sign, name in [(1, "even"), (-1, "odd")]:
        state = {**cold(2), **{k: sign * v for k, v in cold(3).items()}}
        es, s = run({k: v / np.sqrt(2) for k, v in state.items()})
        # An actual eigenvector on the record after transfer. This proves
        # that the full local spectral measure, not just the mean, transfers.
        acted = apply(s, h, 0)
        expected = 8 + (sign < 0)
        assert (
            max(
                abs(acted.get(k, 0) - expected * s.get(k, 0))
                for k in acted.keys() | s.keys()
            )
            < 1e-12
        )
        nominal[name] = np.round(np.sum(es, axis=1), 12).tolist()
    return dict(
        six_port_stage_energies=nominal,
        all_99_free_blank_to_written_columns_zero=True,
        all_99_magnitude_blocks_invariant=True,
        both_zero_labels_dark_on_all_loads=True,
        transferred_even_odd_eigenvalues=[8, 9],
        ordinary_parity_readout_checked=True,
        initial_2H_minus_2E_determinants_at_terminal_energies=dets,
    )
