"""The rotations of the echosahedron and the incidence of its ports, seams and wedges.

A rotation is a permutation of the twelve ports of the fixed cell that maps seams to seams and is realized by a
proper rotation of the icosahedron. A wedge is an ordered triple ``(a, b, c)``: the two seams ``(a, b)`` and
``(b, c)`` that meet at port ``b``.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from .carrier import automorphisms, cell, icosahedron_from_coordinates, isomorphism


@lru_cache(maxsize=None)
def rotations() -> tuple[tuple[int, ...], ...]:
    """The sixty proper rotations of the cell, as port permutations in the labelling of ``CELL_SEAMS``."""
    c = cell()
    pts, seams = icosahedron_from_coordinates()
    iso = isomorphism(c.seams, seams)
    if iso is None:
        raise ValueError('the cell is not the icosahedron')
    X = pts[iso]
    out = []
    for g in automorphisms(c.seams):
        image = X[list(g)]
        R, *_ = np.linalg.lstsq(X, image, rcond=None)
        if np.abs(X @ R - image).max() < 1e-9 and np.linalg.det(R) > 0:
            out.append(tuple(int(v) for v in g))
    if len(out) != 60:
        raise ValueError(f'expected sixty rotations, found {len(out)}')
    return tuple(out)


def ports_of(label) -> frozenset[int]:
    """The ports a register touches: ('port', p), ('seam', (a, b)) or ('wedge', (a, b, c))."""
    kind, where = label
    if kind == 'port':
        return frozenset((int(where),))
    if kind in ('seam', 'wedge'):
        return frozenset(int(v) for v in where)
    raise ValueError(f'unknown register kind {kind!r}')


def incident(first, second, seams) -> bool:
    """Whether two registers of one cell may be coupled by a local term.

    A register is coupled to itself; two ports when a seam joins them; a port to a seam or wedge that contains
    it; a seam to a wedge that contains it; two wedges or two seams when they share a seam or a port of a seam.
    """
    if first == second:
        return True
    seamset = seams if isinstance(seams, (set, frozenset)) else {tuple(sorted(e)) for e in seams}
    a, b = ports_of(first), ports_of(second)
    kinds = {first[0], second[0]}
    if kinds == {'port'}:
        return tuple(sorted(a | b)) in seamset
    if 'port' in kinds:
        return (a <= b) or (b <= a)
    if kinds == {'seam', 'wedge'}:
        seam, wedge = (first, second) if first[0] == 'seam' else (second, first)
        x, y, z = wedge[1]
        return tuple(sorted(seam[1])) in {tuple(sorted((x, y))), tuple(sorted((y, z)))}
    return bool(a & b) and len(a | b) <= 4
