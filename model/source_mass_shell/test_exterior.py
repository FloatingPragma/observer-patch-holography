"""Literal radial source, full halo, compact packets and countercontrols."""
from fractions import Fraction as F
from itertools import product
from math import fsum

import numpy as np
import pytest
import sympy as s

from model.bounded_defect_index.finite import source_matrix, sites
from model.native_mass_channels.algebra import SIGMA
from . import bounds, exterior


def test_packet_action_against_independent_existing_radial_source_matrix():
    value = exterior.packet(1, 2, (F(1, 3), F(-1, 4), F(1, 2)),
                            spacing=F(1, 2), velocity=F(3, 2), mass=2)
    actual = exterior.act(value)
    points = sites(4); index = {point: j for j, point in enumerate(points)}
    lifted = np.zeros(2*len(points), complex)
    target = np.zeros_like(lifted)
    for local in product(range(4), repeat=3):
        point = tuple(local[j]-1+value['anchor'][j] for j in range(3))
        row = 2*index[point]
        lifted[row:row+2] = value['state'][local]
        target[row:row+2] = actual[local]
    old_source = source_matrix(4, spacing=F(1, 2), velocity=F(3, 2), mass=2)
    assert np.linalg.norm(old_source@lifted-target) < 1e-12
    # The source sees the boundary rows outside the packet support.
    support = np.linalg.norm(value['state'], axis=-1) > 1e-15
    assert np.linalg.norm(actual[~support]) > .1


@pytest.mark.parametrize('energy_sign', [-1, 1])
def test_all_valley_packet_and_residual_bounds(energy_sign):
    for width in (2, 4, 8):
        value = exterior.packet(width, 1000, (F(1, 4), F(-2, 3), F(3, 5)),
                                spacing=F(1, 2), velocity=1, mass=3, energy_sign=energy_sign)
        psi = value['state']
        assert abs(fsum(float(z.real)**2+float(z.imag)**2 for z in psi.flat)-1) < 1e-12
        assert all(np.linalg.norm(value['fiber_spinor'][2*j:2*j+2]) > 1e-4 for j in range(8))
        actual, frozen = exterior.act(value), exterior.act(value, freeze_direction=True)
        audit = value['bounds']
        assert np.linalg.norm((frozen-value['energy']*psi).ravel()) <= float(F(audit['kinetic_envelope_error_upper']))
        assert np.linalg.norm((actual-frozen).ravel()) <= float(F(audit['radial_texture_error_upper']))
        assert np.linalg.norm((actual-value['energy']*psi).ravel()) <= float(F(audit['lattice_eigenvector_residual_upper']))


def test_envelope_and_texture_errors_decrease_independently():
    envelope = [exterior.diagnostic(width, 10000, (F(1, 3), F(1, 4), F(2, 5)))
                for width in (4, 8, 16)]
    errors = [row['kinetic_envelope_residual'] for row in envelope]
    assert all(1.5 < x/y < 2.05 for x, y in zip(errors, errors[1:]))
    texture = [exterior.diagnostic(4, distance, (F(1, 3), 0, F(1, 2)))
               for distance in (200, 400, 800)]
    errors = [row['texture_residual'] for row in texture]
    assert all(1.9 < x/y < 2.05 for x, y in zip(errors, errors[1:]))


def test_distant_anchor_does_not_corrupt_the_common_Bloch_phase():
    args = (4,)
    momentum = (F(1, 3), F(1, 5), F(3, 7))
    near = exterior.packet(*args, 100, momentum, spacing=F(1, 2))
    far = exterior.packet(*args, 10**14, momentum, spacing=F(1, 2))
    assert np.array_equal(near['state'], far['state'])
    assert np.array_equal(exterior.act(near, freeze_direction=True),
                          exterior.act(far, freeze_direction=True))
    with pytest.raises(ValueError, match='numerical displacement'):
        exterior.packet(2, 2**54, momentum)
    assert bounds.exterior_packet_bound(1, 1, 1, 2, 2**54)['anchor_lattice_coordinates'][2] == 2**54


def test_exact_hopping_mass_gradient_cancellation_and_translation_obstruction():
    def onsite(point):
        norm = s.sqrt(sum(x*x for x in point))
        return (sum((s.Rational(x)*(-1)**(x % 2)*SIGMA[j]/norm
                     for j, x in enumerate(point)), s.zeros(2)) if norm else s.zeros(2))
    for n in ((0, 0, 0), (1, 2, 3), (-2, 1, 4)):
        for j in range(3):
            y = tuple(x+int(j == q) for q, x in enumerate(n))
            un = s.Matrix(n)/s.sqrt(sum(x*x for x in n)) if any(n) else s.zeros(3, 1)
            uy = s.Matrix(y)/s.sqrt(sum(x*x for x in y))
            direct = s.I*SIGMA[j]*onsite(y)+onsite(n)*s.I*SIGMA[j]
            gradient = s.I*sum(((-1)**(n[l] % 2)*(un[l]-uy[l])*SIGMA[l]*SIGMA[j]
                                for l in range(3)), s.zeros(2))
            assert (direct-gradient).applyfunc(s.simplify) == s.zeros(2)
    for shift in ((1, 0, 0), (2, 0, 0), (2, -4, 6)):
        difference = onsite(tuple(-x for x in shift))-onsite((0, 0, 0))
        assert (difference*difference-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)


def test_arbitrarily_small_certified_residual_does_not_execute_huge_packet():
    previous = None
    for n in (2, 4, 8, 16):
        result = bounds.exterior_packet_bound(F(1, n), 1, 3, n**3, 2*n**8, (1, 1, 1))
        residual = F(result['continuum_eigenvector_residual_upper'])
        if previous is not None:
            assert residual < previous/3
        previous = residual
    large = bounds.certificate()['same_radial_source_exterior_witness']
    assert F(large['continuum_eigenvector_residual_upper']) < F(61, 100000)
    assert large['full_packet_rows'] == 16*10**18
    assert large['eigensolver_or_large_packet_execution_used_for_bound'] is False


@pytest.mark.parametrize('width,distance', [(True, 100), (0, 100), (2, 3), (2, False), (2, 0)])
def test_packet_domain_rejects_wrong_cell_or_staggering_frame(width, distance):
    with pytest.raises(ValueError):
        bounds.exterior_packet_bound(1, 1, 1, width, distance)
