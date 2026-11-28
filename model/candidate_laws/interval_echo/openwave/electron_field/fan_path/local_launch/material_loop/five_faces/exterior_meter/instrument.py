"""Native single-contact pulses followed by a read of the inherited 1045 head."""

import numpy as np

from ..apron import OUT, instrument
from ..wheel import (
    EDGES,
    RIM,
    WHEEL,
    edge,
    face_gate,
    flip,
    head,
    sides,
    wheel_matrices,
)


SHARED = (1053, 1052)
PRECEDING = (1048, 1053)


def swap(basis, side):
    """An isolated full material exchange includes the untouched words."""
    k = WHEEL.index(edge(*side))
    ends = WHEEL[k]
    ids = {word: i for i, word in enumerate(basis)}
    dest = []
    for port, word in basis:
        if port in ends and head(ends, word, k) == port:
            row = (sum(ends) - port, word ^ (1 << k))
        else:
            row = (port, word)
        dest.append(ids[row])
    u = np.eye(len(basis))[:, dest]
    np.testing.assert_array_equal(u @ u, np.eye(len(basis)))
    return u


def exterior(basis, route):
    closed, src, flags, edges = instrument(basis, route, 0)
    opened, src1, flags1, edges1 = instrument(basis, route, 0.5)
    assert src == src1 and edges == edges1
    np.testing.assert_array_equal(flags, flags1)
    end = route["ordered_two_word_trace"][-1]
    for _, word in end:
        for link, expected in (((1045, 1052), 1045), ((1045, 1053), 1053)):
            assert head(link, word, EDGES.index(link)) == expected
    # One unit of the native face Laplacian is I - R, even on disabled words.
    full_face = np.eye(len(flags)) - (opened - closed) / 0.5
    np.testing.assert_array_equal(full_face.T @ full_face, np.eye(len(flags)))
    arrival = full_face[:, src]
    click = arrival[flags].T @ arrival[flags]
    k = WHEEL.index(edge(*SHARED))
    eligible = [head(WHEEL[k], w, k) == SHARED[1] for _, w in basis]
    np.testing.assert_array_equal(click, np.diag(eligible))
    # Charge conjugation complements both prepared outside heads and every
    # wheel head; vacancy transfer is enabled at the end opposite the head.
    all_edges = [*WHEEL, *edges]
    outside_bits = sum(
        (head(e, end[0][1], EDGES.index(e)) == e[0]) << j
        for j, e in enumerate(edges, len(WHEEL))
    )
    gate = face_gate(all_edges, OUT)
    mask = (1 << len(all_edges)) - 1
    read_bit = 1 << all_edges.index((1045, 1052))
    for (port, w), open_word in zip(basis, eligible, strict=True):
        normal = w | outside_bits
        opposite = normal ^ mask
        assert flip(opposite, gate) == flip(normal, gate) ^ mask
        assert bool((flip(opposite, gate) ^ opposite) & read_bit) == open_word
        if port in SHARED:
            assert (head(WHEEL[k], w, k) == port) == (
                head(WHEEL[k], opposite, k) != port
            )
    identity = np.eye(len(basis))
    shared = swap(basis, SHARED)
    preceding = swap(basis, PRECEDING)
    # exp(-i pi (I-S)/4) and its complex conjugate use positive native waits.
    short = ((identity + shared) - 1j * (identity - shared)) / 2
    long = short.conj()
    np.testing.assert_array_equal(short @ short @ short, long)
    np.testing.assert_array_equal(short.conj().T @ short, identity)

    def after(pulse):
        return pulse.conj().T @ click @ pulse

    plus, minus = after(short), after(long)
    score = plus - minus
    _, currents = wheel_matrices(basis, 0.5)
    ix = sides(RIM).index(SHARED)
    # Every entry is dyadic: these comparisons are exact on the full space.
    np.testing.assert_array_equal(score, 1j * (click @ shared - shared @ click))
    np.testing.assert_array_equal(4 * score, currents[ix])
    np.testing.assert_array_equal(
        score @ score, np.diag((shared != identity).any(axis=0))
    )
    effects = {
        "bare": click,
        "shared": after(shared),
        "preceding_then_shared": after(shared @ preceding),
        "agreement": (identity + score) / 2,
    }
    setup = {
        "rim_side": list(SHARED),
        "outside_face": list(OUT),
        "read_head": [1045, 1052],
        "space_before_pulse": len(basis),
        "space_with_outside_face": len(flags),
        "eligible_outside_inputs": int(sum(eligible)),
        "active_material_pairs": int(np.linalg.matrix_rank(score)) // 2,
        "charge_conjugate_words_checked": len(basis),
        "material_rate": 4,
        "short_long_and_idle_in_units_of_pi": [1 / 16, 3 / 16, 1 / 8],
        "outside_full_wait_in_units_of_pi": 1,
        "closed_face_click": 0,
        "no_quadrature_arm_contrast": 0,
    }
    return effects, plus, minus, currents, setup
