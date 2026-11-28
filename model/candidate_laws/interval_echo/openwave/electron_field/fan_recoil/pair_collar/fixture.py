"""Enroll the native collar, two induced paths and their blank signed sectors."""

from functools import lru_cache
from itertools import pairwise

import numpy as np
from scipy.sparse import diags

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation

SOURCE = (1010, 1045)
CELLS = 48
RATE = 0.5


def geometry():
    board = federation(2)
    assert SOURCE[1] in board.neighbours[SOURCE[0]]
    partner = {x: y for a, b in board.external for x, y in ((a, b), (b, a))}

    @lru_cache(None)
    def choices(a, stage):
        assert board.degree[a] == 6
        owner = a // 12
        return tuple(
            (a, m, q, partner[q])
            for m in sorted(board.neighbours[a])
            if m // 12 == owner and board.degree[m] == 5
            for q in sorted(board.neighbours[m])
            if q != a
            and q // 12 == owner
            and q in partner
            and (stage or q not in board.neighbours[a])
        )

    def find(start, excluded, first=None):
        def walk(a, seen, path):
            if len(path) == CELLS + 1:
                return path
            for x, m, q, after in choices(a, not path):
                triple = (x, m, q)
                if not path and first is not None and triple != first:
                    continue
                if len(path) < CELLS and after // 12 in seen:
                    continue
                result = walk(after, seen | {after // 12}, (*path, triple))
                if result is not None:
                    return result
            return None

        return walk(start, excluded | {start // 12}, ())

    forced = (1045, 1052, 1047)

    def owners(path):
        return {row[0] // 12 for row in path}

    baseline = find(SOURCE[1], {SOURCE[0] // 12}, forced)
    for protected in range(2, CELLS + 2):
        left = find(SOURCE[0], owners(baseline[:protected]))
        if left is None:
            continue
        right = find(SOURCE[1], owners(left), forced)
        if right is not None:
            paths = [
                [int(p) for triple in arm for p in triple] for arm in (left, right)
            ]
            bulk = {p for arm in paths for p in arm[3:]}
            edges = {tuple(sorted(e)) for arm in paths for e in pairwise(arm[3:])}
            measured = {tuple(sorted((a, b))) for a, b in board.seams if {a, b} <= bulk}
            if measured == edges:
                break
    else:
        raise ValueError("No two induced native bulk paths from the pair collar")
    for arm in paths:
        np.testing.assert_array_equal(board.degree[arm], np.tile([6, 5, 6], CELLS + 1))
        assert all(b in board.neighbours[a] for a, b in pairwise(arm))
        assert len({p // 12 for p in arm[::3]}) == CELLS + 1
    assert not (owners(left) & owners(right))
    return board, paths, protected


def native_numbers():
    blank = law.seam_term(8)
    n = law.SRC.index
    assert blank[n((1, 1)), n((1, 1))] == 0
    for sign in (-1, 1):
        a, b = (1 + sign, 1), (1, 1 + sign)
        assert blank[n(a), n(a)] == blank[n(b), n(b)] == 8.5
        assert blank[n(a), n(b)] == blank[n(b), n(a)] == -0.5
    k0, k1 = law.local_term_parts()
    j = np.eye(law.DIM)[:, :9]
    np.testing.assert_allclose((k0 + 8 * k1) @ j, j @ blank, atol=5e-14, rtol=0)
    return k0 + 8 * k1


def source_vectors():
    generator = native_numbers()
    echo = law.echo()
    n = law.SRC.index
    raw = np.eye(law.DIM)[:, [n((0, 2)), n((2, 0)), n((1, 1))]]
    chambers = echo @ raw
    # Both words recorded by the source observer are 11 before the pulse.
    for j in range(3):
        for idx in np.flatnonzero(abs(chambers[:, j]) > 1e-14):
            assert idx >= 9 and law.SRC[idx % 9] == (1, 1)
    np.testing.assert_allclose(chambers.T @ chambers, np.eye(3), atol=2e-15, rtol=0)
    np.testing.assert_allclose(generator @ chambers, 0, atol=5e-14, rtol=0)
    created = echo @ chambers
    np.testing.assert_allclose(created, raw, atol=2e-15, rtol=0)
    np.testing.assert_allclose(echo @ created, chambers, atol=2e-15, rtol=0)
    np.testing.assert_allclose(
        generator @ raw,
        raw @ np.array([[32.5, -0.5, 0], [-0.5, 32.5, 0], [0, 0, 0]]),
        atol=5e-14,
        rtol=0,
    )
    return generator, echo, raw, chambers


def signed_generator(board, ports, sign, reference, cut=False):
    assert sign in (-1, 1) and reference in (0, 1)
    degrees = np.asarray(board.degree)[ports]
    # All physical incident contrasts remain. External fixed unit loads
    # suppress only their swap. The private comparison uses the exact D^2 term.
    virtual = 6 - degrees
    onsite = 8 * (
        degrees + virtual * ((1 + sign - reference) ** 2 - (1 - reference) ** 2)
    )
    np.testing.assert_array_equal(onsite, 48 + 16 * sign * (1 - reference) * virtual)
    links = np.full(len(ports) - 1, RATE)
    if cut:
        links[2::3] = 0
    diag = onsite + np.r_[0, links] + np.r_[links, 0]
    return diags([-links, diag, -links], [-1, 0, 1], format="csr")


def assert_joint_columns(board, paths, reference):
    # Independent physical diagonal check over all two-charge positions.
    # Cross terms cannot connect the disjoint induced active paths.
    left, right = [a[3:] for a in paths]
    for sign_left in (-1, 1):
        a = signed_generator(board, left, sign_left, reference)
        b = signed_generator(board, right, -sign_left, reference)
        for ip, p in enumerate(left):
            for iq, q in enumerate(right):
                assert q not in board.neighbours[p]
                first = len(board.neighbours[p])
                second = len(board.neighbours[q])
                v = 8 * (first + second)
                for port, sgn in ((p, sign_left), (q, -sign_left)):
                    v += (
                        8
                        * (6 - board.degree[port])
                        * ((1 + sgn - reference) ** 2 - (1 - reference) ** 2)
                    )
                # The only surviving exchange diagonals are the within-arm
                # incident contacts. On all other links the guard parks the neighbour.
                v += RATE * (
                    int(ip > 0)
                    + int(ip < len(left) - 1)
                    + int(iq > 0)
                    + int(iq < len(right) - 1)
                )
                assert a[ip, ip] + b[iq, iq] == v
