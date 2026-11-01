"""Independent source-edge, full-band and texture-degree controls."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from ..field_symmetry.symmetry import hidden_spin_matrix
from ..full_fock_field.field import FieldGeometry
from . import algebra as a


def exact(matrix):
    if hasattr(matrix, 'toarray'):
        matrix = matrix.toarray()
    return s.ImmutableMatrix([[s.Rational(float(z.real))+s.I*s.Rational(float(z.imag))
                               for z in row] for row in np.asarray(matrix)])


def test_exhaustive_edges_against_original_open_cube_hopping():
    geometry = FieldGeometry((2, 2, 2), speed=1/8, interaction=1.)
    h = exact(geometry.single_h())
    for index in range(4):
        coefficients = tuple(int(j == index) for j in range(4))
        potential = s.diag(*(a.onsite(x, coefficients) for x in geometry.sites))
        assert h*potential+potential*h == s.zeros(16)
    # An untextured physical spin field fails on a source link in its own axis.
    wrong = s.diag(*(a.SIGMA[0] for _ in geometry.sites))
    assert h*wrong+wrong*h != s.zeros(16)
    report = a.exhaustion_certificate()
    assert (report['real_unknowns'], report['constraint_rank'], report['solution_dimension']) == (32, 28, 4)
    assert a.source_anchor_certificate()['actual_source_constraint_rank'] == 28


def test_source_anchor_rejects_wrong_hopping_axis(monkeypatch):
    from scipy import sparse
    original = FieldGeometry.single_h

    def changed(geometry):
        h = original(geometry).toarray()
        h[0:2, 8:10] = 1j*np.diag([1, -1])/16
        h[8:10, 0:2] = h[0:2, 8:10].conj().T
        return sparse.csr_matrix(h)

    monkeypatch.setattr(FieldGeometry, 'single_h', changed)
    with pytest.raises(ArithmeticError, match='actual source'):
        a.source_anchor_certificate()


@pytest.mark.parametrize('site', [(-3, 2, 5), (0, 0, 0), (4, 7, 2)])
def test_edge_recursion_on_unbounded_integer_sites(site):
    coefficients = (F(2, 3), F(-7, 5), 11, F(1, 9))
    source = a.onsite(site, coefficients)
    for j, sigma in enumerate(a.SIGMA):
        neighbour = tuple(x+int(i == j) for i, x in enumerate(site))
        assert a.onsite(neighbour, coefficients) == -sigma*source*sigma
        second = tuple(x+2*int(i == j) for i, x in enumerate(site))
        assert a.onsite(second, coefficients) == source


def test_spin_textures_are_source_hidden_frame_neel_order():
    for site in product(range(3), repeat=3):
        for j in range(3):
            expected = (-1)**sum(site)*exact(hidden_spin_matrix(site, j))
            assert a.onsite_basis(site)[j+1] == expected


def test_folded_mass_keeps_every_valley_and_characteristic_polynomial():
    report = a.clifford_certificate()
    assert report['folded_dimension'] == 16 and report['weyl_valleys'] == 8
    k, masses = (F(1, 2), F(2, 3), F(3, 5)), (F(7, 4), 1, 2, 2)
    h = a.folded_hamiltonian(k, masses)
    x = s.Symbol('x')
    r2 = sum(v*v for v in masses[1:])
    k2 = sum(v*v for v in k)
    # Direct determinant knows neither the claimed eigenvalue formula nor gap.
    expected = ((x*x-s.Rational(k2)-s.Rational(masses[0]**2+r2))**2
                -4*s.Rational(masses[0]**2*r2))**4
    assert h.charpoly(x).as_expr().expand() == expected.expand()
    assert sum(a.spectrum(k, masses).values()) == 16
    assert s.expand(s.prod((x-e)**n for e, n in a.spectrum(k, masses).items())) == expected.expand()


def test_four_dirac_blocks_are_complete_irreducible_and_rotatable():
    report = a.dirac_decomposition_certificate()
    assert report['dirac_block_count'] == 4
    assert sorted(sum(report['block_rows'], [])) == list(range(16))
    gamma, beta = a.folded_generators()
    for i, j in ((1, 2), (2, 3), (3, 1)):
        rotation_generator = beta[i]*beta[j]
        assert rotation_generator.H == -rotation_generator
        assert all(rotation_generator*g == g*rotation_generator for g in gamma)


@pytest.mark.parametrize('mass', [(0, 0, 0, 0), (0, 1, 2, 2), (5, 0, 0, 0),
                                  (3, 1, 2, 2), (-3, 1, 2, 2), (5, 1, 2, 2)])
def test_exact_zero_momentum_spectrum_including_gap_closing(mass):
    actual = a.folded_hamiltonian((0, 0, 0), mass).eigenvals()
    assert actual == a.spectrum((0, 0, 0), mass)
    assert min(e*e for e in actual) == a.bulk_gap_squared(mass)


def test_three_components_and_no_accidental_spin_sign_wall():
    assert a.phase((2, 0, 0, 1)) == 'positive_scalar'
    assert a.phase((-2, 0, 0, 1)) == 'negative_scalar'
    assert a.phase((1, 0, 0, 2)) == 'spin_texture'
    assert a.phase((1, 0, 0, -2)) == 'spin_texture'
    assert a.phase((1, 0, 0, 1)) == a.phase((0, 0, 0, 0)) == 'gapless'
    for t in (F(j, 12) for j in range(25)):
        mass = a.wall_rotation(t, F(7, 3))
        assert mass[0] == mass[3] == 0
        assert sum(x*x for x in mass[1:]) == F(49, 9)
        assert a.bulk_gap_squared(mass) == s.Rational(49, 9)
        assert a.phase(mass) == 'spin_texture'
    assert a.wall_rotation(0) == (0, 1, 0, 0)
    assert a.wall_rotation(1) == (0, 0, 1, 0)
    assert a.wall_rotation(2) == (0, -1, 0, 0)


def test_spatial_rotation_is_not_silently_an_exact_uniform_mass_identity():
    geometry = FieldGeometry((2, 2, 2), speed=1., interaction=0.)
    h = exact(geometry.single_h())
    potential = s.diag(*(a.onsite(x, a.wall_rotation(F(x[0], 2))) for x in geometry.sites))
    assert all(a.bulk_gap_squared(a.wall_rotation(F(x[0], 2))) == 1 for x in geometry.sites)
    assert h*potential+potential*h != s.zeros(16)


@pytest.mark.parametrize('mass', [(5, 1, 2, 2), (-5, 1, 2, 2), (1, 1, 2, 2),
                                  (F(-1, 5), F(1, 3), 0, 0), (100, -50, 0, 0)])
def test_constructive_component_retractions_preserve_the_exact_gap(mass):
    original_phase = a.phase(mass)
    for t in (F(i, 16) for i in range(17)):
        changed = a.component_retraction(mass, t)
        assert a.phase(changed) == original_phase
        assert a.bulk_gap_squared(changed) > 0
    endpoint = a.component_retraction(mass, 1)
    assert a.component_retraction(endpoint, F(3, 7)) == endpoint
    if original_phase == 'spin_texture':
        assert endpoint[0] == 0 and sum(abs(x) for x in endpoint[1:]) == 1
    else:
        assert endpoint == (1 if mass[0] > 0 else -1, 0, 0, 0)


@pytest.mark.parametrize('direction', [(1, 2, 3), (-2, 3, -5), (7, -11, 13)])
def test_oriented_degree_uses_exact_generic_ray_and_all_faces(direction):
    positive = a.texture_degree(a.VERTICES, direction)
    negative = a.texture_degree(tuple(tuple(-x for x in v) for v in a.VERTICES), direction)
    assert positive['degree'] == 1 and negative['degree'] == -1
    assert positive['faces_checked'] == negative['faces_checked'] == 8
    assert len(positive['intersections']) == len(negative['intersections']) == 1
    assert positive['fermionic_zero_mode_proved'] is False


def test_no_winding_when_all_target_faces_lie_in_one_hemisphere():
    translated = tuple(tuple(x+y for x, y in zip(v, (10, 20, 40))) for v in a.VERTICES)
    assert a.texture_degree(translated)['degree'] == 0


def test_rational_perturbation_homotopy_preserves_degree():
    changed = tuple(tuple(F(x)+F((i+2*j)%5-2, 50) for j, x in enumerate(v))
                    for i, v in enumerate(a.VERTICES))
    stable = a.hedgehog_stability(changed)
    assert stable['degree'] == a.texture_degree(changed)['degree'] == 1
    assert F(stable['face_separation_lower']) > 0
    anti = tuple(tuple(-x for x in v) for v in changed)
    assert a.hedgehog_stability(anti, -1)['degree'] == a.texture_degree(anti)['degree'] == -1
    with pytest.raises(ValueError):
        a.hedgehog_stability(anti)


@pytest.mark.parametrize('bad', [True, 1., float('nan'), float('inf'), 1j, '1/2', None, s.Rational(1, 2)])
def test_inexact_or_unowned_coefficient_types_are_rejected(bad):
    with pytest.raises(ValueError):
        a.onsite((0, 0, 0), (0, 1, bad, 0))


def test_mutated_fraction_and_returned_certificate_cannot_forge_next_result():
    bad = F(1, 2)
    bad._denominator = -2
    with pytest.raises(ValueError):
        a.rational(bad)
    result = a.exhaustion_certificate()
    result['basis_columns'][0][0] = 1000
    assert a.exhaustion_certificate()['basis_columns'][0][0] == 1


@pytest.mark.parametrize('bad', [(), [1, 2], [0]*5, np.zeros(4), iter((0, 1, 2, 3))])
def test_complete_coefficient_domain_required(bad):
    with pytest.raises(ValueError):
        a.phase(bad)


@pytest.mark.parametrize('site', [(True, 0, 0), (0., 0, 0), [0, 0, 0], (0, 0)])
def test_site_domain_is_exact(site):
    with pytest.raises(ValueError):
        a.onsite_basis(site)


def test_degree_cannot_admit_empty_degenerate_or_boundary_data():
    for values in ((), a.VERTICES[:-1], [(0, 0, 0)]*6):
        with pytest.raises(ValueError):
            a.texture_degree(values)
    for direction in ((0, 0, 0), (1, 0, 0)):
        with pytest.raises(ValueError):
            a.texture_degree(a.VERTICES, direction)
    changed = list(a.VERTICES)
    changed[0] = (F(4, 3), 0, 0)
    with pytest.raises(ValueError, match='strict margin'):
        a.hedgehog_stability(changed)
    with pytest.raises(ValueError):
        a.hedgehog_stability(a.VERTICES, True)
    for t, magnitude in ((-1, 1), (3, 1), (1, 0), (1, -1)):
        with pytest.raises(ValueError):
            a.wall_rotation(t, magnitude)
    for mass, t in (((0, 0, 0, 0), 0), ((1, 1, 0, 0), 1), ((1, 0, 0, 0), 2)):
        with pytest.raises(ValueError):
            a.component_retraction(mass, t)
