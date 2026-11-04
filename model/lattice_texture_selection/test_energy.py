"""Independent bond operators, deformation paths and topology crossing."""
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import energy


def test_covariance_energy_equals_direct_rotated_two_spin_operator():
    sigma = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    a, b = -s.Rational(1, 4), -s.Rational(1, 8)
    covariance = (a, b, b)
    rho = s.eye(4)/4+sum((c*s.kronecker_product(p, p) for c, p in zip(covariance, sigma, strict=True)), s.zeros(4))
    assert all(value >= 0 for value in rho.eigenvals()) and s.trace(rho) == 1
    H = sum((s.kronecker_product(p, p)/4 for p in sigma), s.zeros(4))-s.eye(4)/4
    R, Q = energy.rotation((1, 2, 0, -1)), energy.rotation((2, -1, 1, 3))
    rotated = sum((s.kronecker_product(sum((R[i, j]*sigma[j] for j in range(3)), s.zeros(2)),
                                      sum((Q[i, j]*sigma[j] for j in range(3)), s.zeros(2)))/4
                   for i in range(3)), s.zeros(4))-s.eye(4)/4
    actual = s.simplify(s.trace(rho*(rotated-H)))
    assert actual == energy.bond_excess(R, Q, a, b)
    assert 0 <= actual <= 1
    assert s.trace(rho*H) == a+2*b-s.Rational(1, 4)
    # This positive two-spin control validates the operator identity; its
    # covariance is not assigned to the infinite cubic quantum ground.


def test_lexicographic_cut_against_direct_edge_set_count():
    for side in range(1, 5):
        points = sorted(product(range(side), repeat=3))
        edges = [(x, y) for i, x in enumerate(points) for y in points[i+1:]
                 if sum(abs(a-b) for a, b in zip(x, y, strict=True)) == 1]
        direct = []
        for prefix in range(len(points)+1):
            done = set(points[:prefix])
            direct.append(sum((a in done) != (b in done) for a, b in edges))
        assert tuple(direct) == energy.cut_sizes(side)


def test_explicit_rotation_reset_path_on_all_eight_cube_sites():
    points = list(product(range(2), repeat=3))
    quaternions = {p: (1, i % 3, (i+1) % 3, (i+2) % 3) for i, p in enumerate(points)}
    initial = {p: energy.rotation(q) for p, q in quaternions.items()}
    edges = set()
    for p in points:
        for axis in range(3):
            for sign in (-1, 1):
                other = list(p)
                other[axis] += sign
                edges.add(tuple(sorted((p, tuple(other)))))
    a, b = -s.Rational(1, 4), -s.Rational(1, 8)
    def total(rotations):
        return sum(energy.bond_excess(rotations.get(x, s.eye(3)), rotations.get(y, s.eye(3)), a, b)
                   for x, y in edges)
    initial_energy = total(initial)
    rotations = dict(initial)
    for point in points:
        q = quaternions[point]
        for t in (s.S.Zero, s.Rational(1, 3), s.Rational(2, 3), s.S.One):
            path_q = [(1-t)*q[0]+t]+[(1-t)*v for v in q[1:]]
            rotations[point] = energy.rotation(path_q)
            assert rotations[point].T*rotations[point] == s.eye(3)
            assert 0 <= total(rotations) <= initial_energy+energy.deformation_extra_bound(2)
        assert rotations[point] == s.eye(3)
    assert total(rotations) == 0


def test_exact_convex_hull_crossing_and_absolute_floor():
    u, v, w = s.Matrix([0, 1, 0]), s.Matrix([1, 0, 0]), s.Matrix([-1, 0, 0])
    assert (v+w)/2 == s.zeros(3, 1)
    assert (u-v).norm()**2+(u-w).norm()**2 == 4
    result = energy.certificate()
    assert result['crossing_two_source_chord_squares_lower'] == 4
    assert not result['crossing_floor_is_barrier_above_initial_energy']
    assert energy.shell_cost([1, -1, 0, 2], 3, F(1, 5)) == F(12, 5)


def test_nonacute_diagonal_cannot_hide_its_actual_edge_cost():
    u = s.Matrix([1, 0, 0])
    v, w = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]), s.Matrix([s.Rational(3, 5), -s.Rational(4, 5), 0])
    assert u.dot(v) > 0 and u.dot(w) > 0 and v.dot(w) < 0
    diagonal = (v-w).dot(v-w)
    source_chords = (u-v).dot(u-v)+(u-w).dot(u-w)
    assert 2 < diagonal <= 2*source_chords
    assert source_chords == s.Rational(8, 5) > 1
    assert energy.mixed_shell_cost(2, [1, -1, 0, 2], 3, F(1, 5)) == {
        'constant': F(3, 5), 'coefficient_of_pi': F(12, 5)}


def test_nonacute_shell_count_does_not_accept_boolean_or_negative_values():
    for count in (True, -1, .5):
        with pytest.raises(ValueError):
            energy.mixed_shell_cost(count, [])


@pytest.mark.parametrize('value', [0, -1, 1., True])
def test_finite_energy_domain_rejects_unsupported_parameters(value):
    with pytest.raises(ValueError):
        energy.deformation_extra_bound(value)
