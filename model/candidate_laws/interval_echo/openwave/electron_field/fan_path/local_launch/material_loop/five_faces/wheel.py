"""Route the two actual arriving words into a five-face native component."""

from collections import deque
from itertools import pairwise

import numpy as np

from model.executable_model import federation
from ....fixture import terminal_fan


A, LAND, CENTER = 1045, 1052, 1047
RIM = (1046, 1048, 1053, 1052, 1050)
FAN = (1045, 1049, 1053, 1052, 1051, 1044)
FEED = (1045, 1052)
ROUTE = ((1045, 1049, 1053), (1045, 1052, 1053), (1047, 1052, 1053), FEED)


def edge(a, b):
    return (min(a, b), max(a, b))


def sides(face):
    return list(pairwise((*face, face[0])))


FACES = tuple((CENTER, a, b) for a, b in sides(RIM))
WHEEL = sorted({edge(a, b) for f in FACES for a, b in sides(f)})
LOCAL = sorted({edge(a, b) for j in range(1, 5) for a, b in sides((0, j, j + 1))})
WRITTEN = [edge(FAN[i], FAN[j]) for i, j in LOCAL]
EDGES = WRITTEN + sorted(set(WHEEL) - set(WRITTEN))


def head(e, word, i):
    return e[0] if word & (1 << i) else e[1]


def bit_for(e, h, i):
    return int(e[0] == h) << i


def face_gate(edges, face):
    order = sides(face)
    mask = sum(1 << edges.index(edge(a, b)) for a, b in order)
    direction = sum(bit_for(edge(a, b), b, edges.index(edge(a, b))) for a, b in order)
    return mask, direction


def flip(word, gate):
    mask, cw = gate
    return word ^ mask if word & mask in (cw, mask ^ cw) else word


def raw_word(w):
    return sum(
        bit_for(e, FAN[i] if w >> j & 1 else FAN[k], j)
        for j, ((i, k), e) in enumerate(zip(LOCAL, WRITTEN, strict=True))
    )


def project_word(word):
    return sum(
        bit_for(e, head(e, word, EDGES.index(e)), i) for i, e in enumerate(WHEEL)
    )


def in_sector(x, w):
    hs = [head(e, w, i) for i, e in enumerate(WHEEL)]
    return (
        x in RIM
        and hs.count(CENTER) == 4
        and all(hs.count(p) == 1 + int(x == p) for p in RIM)
    )


def geometry():
    board = federation(2)
    faces = {frozenset(f) for f in board.faces()}
    assert all(frozenset(f) in faces for f in (*ROUTE[:3], *FACES))
    wheels = {h: terminal_fan(board, h)[1:] for h in range(1044, 1056)}
    equal = [
        h for h, r in wheels.items() if len({int(board.degree[p]) for p in r}) == 1
    ]
    assert equal == [CENTER] and wheels[CENTER] == RIM
    assert all(board.degree[p] == 5 for p in RIM)
    return wheels


def preparation():
    geometry()
    gates = [face_gate(EDGES, f) for f in ROUTE[:3]]
    feeder = EDGES.index(FEED)

    def change(row, pulse):
        x, w = row
        if pulse < 3:
            return x, flip(w, gates[pulse])
        if x in FEED and head(FEED, w, feeder) == x:
            return (FEED[1] if x == FEED[0] else FEED[0]), w ^ (1 << feeder)
        return row

    # Search each supplied word and stop at its first shortest preparation.
    candidates = []
    for background in range(1 << (len(EDGES) - 9)):
        start = tuple((A, raw_word(w) | (background << 9)) for w in (334, 365))
        todo = deque([(start, ())])
        visited = {start}
        while todo:
            pair, pulses = todo.popleft()
            (x, w), (y, v) = pair
            if (
                x == y == LAND
                and w ^ v == gates[2][0]
                and in_sector(x, project_word(w))
                and in_sector(y, project_word(v))
            ):
                candidates.append((len(pulses), background, pulses, start, pair))
                break
            for i in range(4):
                new = tuple(change(z, i) for z in pair)
                if new not in visited:
                    todo.append((new, (*pulses, i)))
                    visited.add(new)
    if not candidates:
        raise ValueError(
            "No coherent face transfer in the enrolled three-face apparatus"
        )
    _, background, pulses, start, pair = min(candidates)
    trace = [start]
    for p in pulses:
        trace.append(tuple(change(z, p) for z in trace[-1]))
    # Complement each bit and reverse the charge sign; every ring and
    # material condition is unchanged. Check at every pulse for both words.
    for old, new, idx in zip(trace, trace[1:], pulses):
        for (x, w), (y, v) in zip(old, new, strict=True):
            mask = (1 << len(EDGES)) - 1
            if idx < 3:
                assert flip(mask ^ w, gates[idx]) == mask ^ v and x == y
            else:
                enabled = x in FEED and head(FEED, mask ^ w, feeder) != x
                assert enabled == (x != y) and (mask ^ v) == (
                    (mask ^ w) ^ (1 << feeder if enabled else 0)
                )
    return {
        "center": CENTER,
        "positive_rim_order": list(RIM),
        "input_plus_new_edge_order": [list(e) for e in EDGES],
        "compatible_background_count": len(candidates),
        "chosen_background": background,
        "native_address_schedule": [list(ROUTE[i]) for i in pulses],
        "ordered_two_word_trace": [[list(v) for v in pair] for pair in trace],
        "wheel_input_words": [project_word(w) for _, w in pair],
    }


def wheel_basis():
    basis = [(p, w) for p in RIM for w in range(1 << len(WHEEL)) if in_sector(p, w)]
    assert len(basis) == 30
    return basis


def rotation(basis):
    cycle = {x: y for x, y in sides(RIM)}
    cycle[CENTER] = CENTER
    perm = []
    for x, w in basis:
        v = 0
        for i, e in enumerate(WHEEL):
            dest = edge(cycle[e[0]], cycle[e[1]])
            j = WHEEL.index(dest)
            v |= bit_for(dest, cycle[head(e, w, i)], j)
        perm.append(basis.index((cycle[x], v)))
    return np.array(perm)


def wheel_matrices(basis, g, flux=0, t=4):
    n = len(basis)
    ids = {z: i for i, z in enumerate(basis)}
    h = np.diag(np.full(n, 8, complex))
    currents = np.zeros((5, n, n), complex)
    gates = [face_gate(WHEEL, f) for f in FACES]
    for i, (p, w) in enumerate(basis):
        for gate in gates:
            v = flip(w, gate)
            if v != w:
                h[ids[p, v], i] -= g
                h[i, i] += g
        for a, (u, v) in enumerate(sides(RIM)):
            e = edge(u, v)
            b = WHEEL.index(e)
            if p not in e or head(e, w, b) != p:
                continue
            dest = v if p == u else u
            j = ids[dest, w ^ (1 << b)]
            sign = 1 if p == u else -1
            h[j, i] -= t * np.exp(1j * sign * flux / 5)
            h[i, i] += t
            currents[a, j, i] = 1j * t * sign
    np.testing.assert_allclose(h, h.conj().T, atol=2e-15, rtol=0)
    return h, currents
