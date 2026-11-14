"""Portable twelve-port carrier and canonical ring torus.

The labelled carrier is the pinned RER icosahedron, reproduced here as finite
combinatorial data. Ring ports are the vertices adjacent to two vertices of
face (0, 11, 5), ordered around that face. No archive import is required.
The torus has two carriers per honeycomb cell; this is a 2D screen geometry.
"""
import operator

from .construction import Construction

CARRIER_SEAMS = (
    (0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7),
    (1, 8), (1, 9), (2, 3), (2, 4), (2, 6), (2, 10), (2, 11),
    (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9), (4, 11),
    (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10),
    (8, 9), (10, 11),
)
REFERENCE_FACE = (0, 11, 5)
RING_PORTS = (1, 10, 4)


def carrier(mass=0.0):
    return Construction(12, CARRIER_SEAMS, mass)


def flat_ring(lx, ly, mass=0.0):
    """Return the periodic ring construction, including the native seam order.

    Positive integer cell counts include thin tori; distinct ring ports keep
    their seam events separate even when carrier neighbours coincide.
    """
    lx, ly = operator.index(lx), operator.index(ly)
    if lx < 1 or ly < 1:
        raise ValueError("positive integer cell counts required")

    def cid(i, j, sub):
        return 2 * ((i % lx) * ly + (j % ly)) + sub

    count = 2 * lx * ly
    seams = [(12*c + a, 12*c + b) for c in range(count) for a, b in CARRIER_SEAMS]
    for i in range(lx):
        for j in range(ly):
            a = cid(i, j, 0)
            neighbours = (cid(i, j, 1), cid(i-1, j, 1), cid(i, j-1, 1))
            seams.extend((12*a + port, 12*b + port)
                         for b, port in zip(neighbours, RING_PORTS))
    return Construction(12*count, tuple(seams), mass)
