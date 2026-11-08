"""Exact corner-root topology, independent of floating port choices."""
from functools import lru_cache
from model.physical_space_time import uniform_green as affine


def integer(value, lower=0, upper=None):
    if type(value) is not int or value < lower or (upper is not None and value > upper):
        raise ValueError('an exact integer in the declared domain required')
    return value


class Topology:
    """Evaluate addressed source simplices without constructing the sphere.

    Internal vertex identifiers are exact shared-midpoint names. They are
    not the allocated source array indices; cell addresses are literal.
    """
    def __init__(self):
        self.vertices = {i: (0, None) for i in range(12)}
        self.midpoints = {}
        self.next_vertex = 12

    def midpoint(self, level, left, right):
        key = (level, *sorted((left, right)))
        if key not in self.midpoints:
            self.midpoints[key] = self.next_vertex
            self.vertices[self.next_vertex] = (level, (left, right))
            self.next_vertex += 1
        return self.midpoints[key]

    @lru_cache(None)
    def raw_child(self, level, cell):
        a, b, c = self.cell(level-1, cell//4)
        ab = self.midpoint(level, a, b)
        bc = self.midpoint(level, b, c)
        ca = self.midpoint(level, c, a)
        return ((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca))[cell % 4]

    @lru_cache(None)
    def cell(self, level, cell):
        integer(level); integer(cell, upper=20*4**level-1)
        return affine.FACES[cell] if level == 0 else self.raw_child(level, cell)

    @lru_cache(None)
    def neighbour(self, level, cell, side):
        integer(side, upper=2)
        face = self.cell(level, cell)
        edge = {face[side], face[(side+1) % 3]}
        if level == 0:
            candidates = range(20)
        else:
            raw = self.raw_child(level, cell)
            raw_side = next(i for i in range(3)
                if {raw[i], raw[(i+1) % 3]} == edge)
            child, parent = cell % 4, cell//4
            if child == 3:
                other = 4*parent+(raw_side+1) % 3
            elif raw_side == 1:
                other = 4*parent+3
            else:
                parent_side = child if raw_side == 0 else (child+2) % 3
                endpoint = self.cell(level-1, parent)[child]
                parent_other, _ = self.neighbour(level-1, parent, parent_side)
                other = 4*parent_other+self.cell(level-1, parent_other).index(endpoint)
            candidates = (other,)
        for other in candidates:
            if other == cell:
                continue
            other_face = self.cell(level, other)
            for other_side in range(3):
                if {other_face[other_side], other_face[(other_side+1) % 3]} == edge:
                    return other, other_side
        raise ArithmeticError('the addressed source edge has no matching neighbour')

    def ball(self, level, radius):
        integer(level); integer(radius)
        seen, frontier = {0}, {0}
        for _ in range(radius):
            frontier = {self.neighbour(level, c, side)[0]
                        for c in frontier for side in range(3)}-seen
            seen |= frontier
        return tuple(sorted(seen))


def address(level, cell):
    integer(level); integer(cell, upper=20*4**level-1)
    return divmod(cell, 4**level)


def minimum_level(radius):
    """A sufficient collar: 2**level > radius+2."""
    integer(radius)
    return (radius+2).bit_length()


def signature(level, radius):
    """Canonical corner fan, including every seam incident to its ball."""
    integer(level); integer(radius)
    if level < minimum_level(radius):
        raise ValueError('the corner ball needs its proved original-face collar')
    topology = Topology()
    cells = topology.ball(level, radius)
    return tuple((address(level, c), tuple(address(level,
        topology.neighbour(level, c, k)[0]) for k in range(3))) for c in cells)


def infinite_ball(radius):
    """An exact finite ball of the explicitly identified corner topology.

    This supplies no port routing or coordinates at infinite refinement.
    """
    return signature(minimum_level(radius), radius)


def canonical_window_size(radius):
    return 12*len(infinite_ball(radius))
