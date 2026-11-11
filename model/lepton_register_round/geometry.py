"""Exact carrier geometry over Q(sqrt5): faces, the face third-turn, screen orientations and the hexagon of a face axis.

Coordinates are the public source chart of `model.common_source.slow_geometry`; a number a + b sqrt5 is the pair (a, b).
"""
from fractions import Fraction as F
from functools import lru_cache

from model.common_source import slow_geometry

FACE = (4, 5, 11)
C_SEAM, E_SEAM = (2, 3), (2, 10)


def _add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def _sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def _mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def _sign(x):
    a, b = x
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    if a * a > 5 * b * b:
        return 1 if a > 0 else -1
    return 1 if b > 0 else -1


def _cross(u, v):
    return (_sub(_mul(u[1], v[2]), _mul(u[2], v[1])), _sub(_mul(u[2], v[0]), _mul(u[0], v[2])),
            _sub(_mul(u[0], v[1]), _mul(u[1], v[0])))


def _dot(u, v):
    out = (F(0), F(0))
    for a, b in zip(u, v):
        out = _add(out, _mul(a, b))
    return out


@lru_cache(maxsize=1)
def carrier():
    source = slow_geometry.source()
    seams = tuple(tuple(int(v) for v in e) for e in source.seams())
    rotations = tuple(tuple(int(v) for v in p) for p in source.rotations())
    points = tuple(tuple((F(a), F(b)) for a, b in row) for row in slow_geometry.raw_frame())
    neighbours = {p: frozenset(q for e in seams for q in e if p in e and q != p) for p in range(12)}
    return seams, rotations, points, neighbours


def compose(outer, inner):
    return tuple(outer[inner[x]] for x in range(len(inner)))


def order(p):
    identity, current, n = tuple(range(len(p))), p, 1
    while current != identity:
        current, n = compose(p, current), n + 1
    return n


@lru_cache(maxsize=1)
def faces():
    seams, _, _, nb = carrier()
    return tuple(sorted({tuple(sorted((a, b, c))) for a in range(12) for b in nb[a] for c in nb[b]
                         if c in nb[a] and len({a, b, c}) == 3}))


def antipodal_face(face):
    _, _, _, nb = carrier()
    return next(g for g in faces() if not set(g) & set(face) and all(not (nb[x] & set(g)) for x in face))


def outward(face):
    """The cyclic order of the face that runs counterclockwise seen from outside."""
    _, _, x, _ = carrier()
    p, q, r = face
    normal = _cross(tuple(_sub(a, b) for a, b in zip(x[q], x[p])), tuple(_sub(a, b) for a, b in zip(x[r], x[p])))
    centre = tuple(_add(_add(a, b), c) for a, b, c in zip(x[p], x[q], x[r]))
    return (p, q, r) if _sign(_dot(normal, centre)) > 0 else (p, r, q)


def third_turn(face):
    """The rotation of order three fixing the face and turning it counterclockwise seen from outside."""
    _, rotations, _, _ = carrier()
    p, q, r = outward(face)
    found = [h for h in rotations if h[p] == q and h[q] == r and h[r] == p]
    if len(found) != 1:
        raise ArithmeticError('face third-turn is not unique')
    return found[0]


def half_turn(seam):
    _, rotations, _, _ = carrier()
    a, b = seam
    found = [h for h in rotations if order(h) == 2 and h[a] == b and h[b] == a]
    if len(found) != 1:
        raise ArithmeticError('seam half-turn is not unique')
    return found[0]


def screen_sign(face, a, b):
    """Sign of n_f . (x_a x x_b): +1 when a -> b runs in the sense of the face third-turn about its axis."""
    _, _, x, _ = carrier()
    normal = tuple(_add(_add(u, v), w) for u, v, w in zip(*(x[p] for p in face)))
    return _sign(_dot(normal, _cross(x[a], x[b])))


@lru_cache(maxsize=None)
def hexagon(face=FACE, start=2):
    """Ports of the hexagon around the axis of the face, in the screen sense from `start`."""
    seams, _, _, _ = carrier()
    ports = set(range(12)) - set(face) - set(antipodal_face(face))
    following = {}
    for a, b in seams:
        if a in ports and b in ports:
            if screen_sign(face, a, b) > 0:
                following[a] = b
            else:
                following[b] = a
    loop = [start]
    while len(loop) < len(ports):
        loop.append(following[loop[-1]])
    if following[loop[-1]] != start:
        raise ArithmeticError('hexagon is not a closed loop')
    return tuple(loop)


def hexagon_seams(face=FACE, start=2):
    loop = hexagon(face, start)
    return tuple((loop[i], loop[(i + 1) % len(loop)]) for i in range(len(loop)))


def summary():
    """Geometry of the register: the face fixed by the two seam half-turns, its hexagon and its orientations."""
    seams, rotations, _, _ = carrier()
    g_c, g_e = half_turn(C_SEAM), half_turn(E_SEAM)
    g = compose(g_e, g_c)
    fixed = [f for f in faces() if {g[i] for i in f} == set(f)]
    rho = third_turn(FACE)
    loop = hexagon()
    position = {p: i for i, p in enumerate(loop)}
    advance = {(position[rho[p]] - position[p]) % 6 for p in loop}
    ring = hexagon_seams()
    e_type = {tuple(sorted(E_SEAM))}
    current = E_SEAM
    for _ in range(2):
        current = (rho[current[0]], rho[current[1]])
        e_type.add(tuple(sorted(current)))
    c_type = {tuple(sorted(s)) for s in ring} - e_type
    orbit_a = sorted({2, rho[2], rho[rho[2]]})
    return dict(
        seams=len(seams), rotations=len(rotations), faces=len(faces()),
        half_turn_product_order=order(g), faces_fixed_by_half_turn_product=[list(f) for f in fixed],
        face=list(FACE), antipodal_face=list(antipodal_face(FACE)), outward_order=list(outward(FACE)),
        third_turn_is_inverse_of_half_turn_product=rho == compose(g, g),
        hexagon_in_screen_sense=list(loop), third_turn_advances_hexagon_by=sorted(advance),
        c_seam_screen_orientation=[a for a in C_SEAM] if screen_sign(FACE, *C_SEAM) > 0 else list(C_SEAM[::-1]),
        e_seam_screen_orientation=[a for a in E_SEAM] if screen_sign(FACE, *E_SEAM) > 0 else list(E_SEAM[::-1]),
        e_type_seams=sorted(list(s) for s in e_type), c_type_seams=sorted(list(s) for s in c_type),
        third_turn_orbit_of_port_2=orbit_a)
