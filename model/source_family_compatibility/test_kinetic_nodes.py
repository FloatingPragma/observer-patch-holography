"""Independent source Fourier, projector flux and valley-mixing controls."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from model.source_family_compatibility import kinetic_nodes as k
from model.full_fock_field.field import FieldGeometry
from model.native_mass_channels import algebra as masses


def test_literal_source_link_symbol_and_all_eight_jacobians():
    geometry = FieldGeometry((3, 3, 3), speed=.125, interaction=1.)
    h = geometry.single_h().toarray()
    origin = geometry.sites.index((1, 1, 1))
    for b in product((0, 1), repeat=3):
        q = np.array([.17, -.23, .31])
        p = np.pi*np.array(b)+q
        action = np.zeros((2, 2), complex)
        for j in range(3):
            for direction in (-1, 1):
                target = [1, 1, 1]; target[j] += direction
                other = geometry.sites.index(tuple(target))
                action += h[2*origin:2*origin+2, 2*other:2*other+2]*np.exp(1j*direction*p[j])
        expected = sum((-geometry.speed*np.sin(p[j])*np.array(masses.SIGMA[j], complex)
                        for j in range(3)), np.zeros((2, 2), complex))
        np.testing.assert_allclose(action, expected, atol=2e-17)
        row = next(row for row in k.node_rows() if row['source_node_pi_units'] == list(b))
        assert row['kinetic_winding'] == -(-1)**sum(b)


def test_projector_flux_direct_eigenvector_sphere_triangulation():
    # Gauge-invariant oriented triangle phases, using eigenvectors of the
    # actual nonlinear sin symbol. The mesh is independent of the proof's
    # linear homotopy and catches the Chern convention and valley signs.
    north = np.array([0., 0., 1.]); south = -north
    rings = [[north]]
    for theta in np.linspace(0, np.pi, 13)[1:-1]:
        rings.append([np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                      for phi in np.arange(24)*2*np.pi/24])
    rings.append([south])
    triangles = []
    for j in range(24):
        triangles.append((north, rings[1][j], rings[1][(j+1)%24]))
        triangles.append((rings[-2][j], south, rings[-2][(j+1)%24]))
    for first, second in zip(rings[1:-2], rings[2:-1]):
        for j in range(24):
            z = (j+1)%24
            triangles.extend(((first[j], second[j], second[z]), (first[j], second[z], first[z])))
    sigma = tuple(np.array(x, complex) for x in masses.SIGMA)
    for row in k.node_rows():
        signs = np.array(row['source_Jacobian_over_speed'])
        for band in (-1, 1):
            phase = 0.
            def vector(x):
                h = sum((value*matrix for value, matrix in zip(signs*np.sin(.6*x), sigma)), np.zeros((2, 2), complex))
                return np.linalg.eigh(h)[1][:, int(band == 1)]
            for a, b, c in triangles:
                u, v, w = vector(a), vector(b), vector(c)
                phase += np.angle(np.vdot(u, v)*np.vdot(v, w)*np.vdot(w, u))
            assert abs(phase/(2*np.pi)-band*row['kinetic_winding']) < 2e-14


def test_gap_bound_and_strict_stability_condition():
    assert k.sphere_gap(F(1), F(3, 2)) == F(5, 4)
    assert k.winding_stability(1, F(5, 6))['degree_preserved_by_straight_homotopy'] is False
    assert k.winding_stability(1, F(4, 5))['degree_preserved_by_straight_homotopy'] is True
    rng = np.random.default_rng(187)
    for radius in (F(1, 20), F(1, 2), F(1)):
        points = rng.normal(size=(200, 3))
        points *= float(radius)/np.linalg.norm(points, axis=1)[:, None]
        assert np.linalg.norm(np.sin(points), axis=1).min() >= float(k.sphere_gap(radius))


@pytest.mark.parametrize('value', [True, 0, -1, 1., F(3, 2), float('nan')])
def test_invalid_spheres_rejected(value):
    with pytest.raises(ValueError):
        k.sphere_gap(value)


def test_actual_onsite_fourier_transfer_pairs_opposite_nodes():
    # The literal source onsite fields, transformed on the complete 2^3 cell,
    # reproduce every block of the folded mass matrices.
    _, beta = masses.folded_generators()
    for channel in range(4):
        for left in k.VALLEYS:
            for right in k.VALLEYS:
                block = sum((((-1)**(sum((x-y)*z for x, y, z in zip(right, left, site)) % 2))*masses.onsite_basis(site)[channel]
                             for site in k.VALLEYS), s.zeros(2))/8
                i, j = 2*k.VALLEYS.index(left), 2*k.VALLEYS.index(right)
                assert block == beta[channel][i:i+2, j:j+2]
    row = k.mass_mixing()
    assert len(row['channels']) == 4
    assert all(row['opposite_chirality_only'] for row in row['channels'])


def test_certificate_keeps_free_and_interacting_scope_separate():
    row = k.certificate()
    assert row['total_kinetic_winding'] == 0
    assert len(row['nodes']) == 8
    assert all(value['physical_particle_assignment'] is None for value in row['nodes'])
    assert row['free_field_scope']['default_U_over_c_eight_in_this_class'] is False
    assert row['scope']['default_interacting_Green_function_winding_derived'] is False
    row['nodes'][0]['kinetic_winding'] = 100
    assert k.node_rows()[0]['kinetic_winding'] == -1
