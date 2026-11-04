"""Independent orientation, discrete degree and energy-edge controls."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import shells


@pytest.mark.parametrize('radius', [1, 2, 3])
def test_closed_oriented_triangulation_and_actual_source_edges(radius):
    triangles, source = shells.shell(radius)
    oriented = Counter()
    all_edges = set()
    for t in triangles:
        a, b, c = t
        face_normal = shells.cross(tuple(y-x for x, y in zip(a, b, strict=True)),
                                  tuple(y-x for x, y in zip(a, c, strict=True)))
        assert shells.dot(a, face_normal) > 0
        for x, y in ((a, b), (b, c), (c, a)):
            oriented[(x, y)] += 1
            all_edges.add(shells.edge(x, y))
        assert shells.edge(a, b) in source and shells.edge(a, c) in source
    assert all(count == oriented[(y, x)] == 1 for (x, y), count in oriented.items())
    vertices = {x for t in triangles for x in t}
    assert len(vertices)-len(all_edges)+len(triangles) == 2
    assert len(source) == 48*radius**2
    assert len(all_edges)-len(source) == 24*radius**2  # interpolation diagonals


@pytest.mark.parametrize('center', [(0, 0, 0), (3, -2, 7)])
def test_exact_degree_has_no_rounded_solid_angles(center):
    triangles, _ = shells.shell(2, center)
    directions = shells.radial_directions(triangles, center)
    forward = shells.degree(triangles, directions)
    backward = shells.degree(tuple((a, c, b) for a, b, c in triangles), directions)
    reversed_target = {p: tuple(-v for v in n) for p, n in directions.items()}
    assert forward['degree'] == 1 and backward['degree'] == -1
    assert shells.degree(triangles, reversed_target)['degree'] == -1
    assert shells.degree(triangles, {p: (1, 0, 0) for p in directions})['degree'] == 0
    for target in ((1, 13, 170), (-1, 13, 170)):
        assert shells.degree(triangles, directions, target)['degree'] == 1


def test_rational_positive_rescaling_does_not_change_cone_degree():
    triangles, _ = shells.shell(1)
    directions = shells.radial_directions(triangles)
    rescaled = {p: tuple(F(i+1, i+2)*x for x in n) for i, (p, n) in enumerate(directions.items())}
    assert shells.degree(triangles, rescaled)['degree'] == 1


def test_spherical_area_bound_against_independent_octant_integral():
    # Three positive coordinate axes bound an octant of the unit sphere.
    theta, phi = s.symbols('theta phi', real=True)
    area = s.integrate(s.integrate(s.sin(theta), (theta, 0, s.pi/2)), (phi, 0, s.pi/2))
    assert area == s.pi/2
    vertices = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    chord_sum = sum(sum((a-b)**2 for a, b in zip(vertices[0], v, strict=True)) for v in vertices[1:])
    assert chord_sum == 4 and area < chord_sum
    assert 2*s.atan(abs(shells.det(*vertices))) == area


def test_edge_support_disjointness_is_checked_without_volume_assumptions():
    shells_by_radius = [shells.shell(r)[1] for r in range(1, 5)]
    assert sum(map(len, shells_by_radius)) == len(set().union(*shells_by_radius))
    assert not shells.shell(2, (-3, 0, 0))[1].intersection(shells.shell(2, (3, 0, 0))[1])
    assert shells.shell(2, (-2, 0, 0))[1].intersection(shells.shell(2, (2, 0, 0))[1])


def test_nonacute_shell_has_an_exact_nearest_neighbor_cost():
    triangles, edges = shells.shell(1)
    directions = {p: (1, 0, 0) for t in triangles for p in t}
    directions[(-1, -1, -1)] = (-1, 0, 0)
    with pytest.raises(ValueError, match='nonnegative'):
        shells.degree(triangles, directions)
    source_sum = sum(sum((a-b)**2 for a, b in zip(directions[x], directions[y], strict=True))
                     for x, y in edges)
    assert source_sum == 12 > 1


@pytest.mark.parametrize('damage', ['omit', 'reverse', 'duplicate'])
def test_exact_degree_requires_a_complete_oriented_closed_shell(damage):
    triangles, _ = shells.shell(1)
    directions = shells.radial_directions(triangles)
    broken = list(triangles)
    if damage == 'omit':
        broken.pop()
    elif damage == 'reverse':
        a, b, c = broken[0]
        broken[0] = (a, c, b)
    else:
        broken.append(broken[0])
    with pytest.raises(ValueError, match='closed triangulation'):
        shells.degree(broken, directions)


@pytest.mark.parametrize('damage', ['nonacute', 'zero', 'missing', 'float', 'target'])
def test_degree_rejects_nonadmissible_or_ambiguous_inputs(damage):
    triangles, _ = shells.shell(1)
    directions = shells.radial_directions(triangles)
    target = None
    if damage == 'nonacute':
        p = next(iter(directions))
        directions[p] = tuple(-x for x in directions[p])
    elif damage == 'zero':
        directions[next(iter(directions))] = (0, 0, 0)
    elif damage == 'missing':
        directions.pop(next(iter(directions)))
    elif damage == 'float':
        directions[next(iter(directions))] = (1., 0, 0)
    else:
        target = (1, 0, 0)
    with pytest.raises(ValueError):
        shells.degree(triangles, directions, target)
