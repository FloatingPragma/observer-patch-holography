from fractions import Fraction as F
from itertools import product
from math import atan2, pi

import numpy as np
import pytest
from scipy import linalg as la

from model.lattice_texture_selection.shells import shell
from model.certified_winding_core.spectrum import certify
from . import certification, spectral, topology


@pytest.mark.parametrize('a', [1, 2, 3])
def test_exact_degree_matches_independent_spherical_area(a):
    value = topology.topology(a)
    for row in value['local_shells']:
        triangles = shell(1, tuple(row['center']))[0]
        total = 0.
        for triangle in triangles:
            u, v, w = [np.array(topology.direction(p, a), float) for p in triangle]
            u, v, w = (x/la.norm(x) for x in (u, v, w))
            total += 2*atan2(np.linalg.det(np.array([u, v, w])),
                             1+u@v+v@w+w@u)
        assert abs(total/(4*pi)-row['degree']) < 1e-12
    assert value['exterior_shell']['degree'] == 0
    assert value['local_shell_edges_disjoint'] == (a > 1)


@pytest.mark.parametrize('bad', [0, -1, 1.5, True])
def test_invalid_geometry_rejected(bad):
    with pytest.raises(ValueError):
        topology.topology(bad)


def test_literal_hopping_matches_independent_tensor_assembly():
    radius, mass, a = 3, 16, 1
    value = spectral.source(radius, mass, a)
    side = 2*radius+1
    derivative = np.diag(np.full(side-1, .5j), 1)+np.diag(np.full(side-1, -.5j), -1)
    independent = np.zeros(value['h'].shape, complex)
    for axis in range(3):
        factors = [derivative if j == axis else np.eye(side) for j in range(3)]
        independent += np.kron(np.kron(np.kron(*factors[:2]), factors[2]), spectral.SIGMA[axis])
    for i, point in enumerate(product(range(-radius, radius+1), repeat=3)):
        n = np.array(topology.direction(point, a), float)
        if la.norm(n):
            n /= la.norm(n)
        independent[2*i:2*i+2, 2*i:2*i+2] += mass*sum(n[j]*spectral.SIGMA[j] for j in range(3))
    assert np.array_equal(independent, value['h'].toarray())
    assert spectral.norm(value['h']@value['C']+value['C']@value['h'].conjugate()) < 1e-12


@pytest.mark.parametrize('a', [1, 2])
def test_complete_core_coupling_gram_and_energy_support(a):
    value = spectral.source(a+2, 16, a)
    core = value['core_indices']
    exterior = [i for i in range(value['h'].shape[0]) if i not in core]
    b = value['h'][exterior][:, core].toarray()
    eigenvalues = la.eigvalsh(b.conj().T@b)
    assert np.allclose(eigenvalues, [1.25, 1.25, 1.75, 1.75] if a == 1 else [1.5]*4)
    assert np.count_nonzero(value['h'][core][:, core].toarray()) == 0
    bound = topology.topology(a)
    assert bound['spin_energy_upper_coefficient_of_J'] == 3*(2*a+3)**2*(2*a+4)


def test_boundary_error_decreases_and_domains_fail_closed():
    bounds = [topology.spectral_bounds(64, 1, radius) for radius in (3, 4, 5)]
    errors = [F(x['ordered_band_energy_infinite_volume_error_upper']) for x in bounds]
    assert errors[2] < errors[1] < errors[0] < F(1, 10**7)
    for mass, a, radius in [(7, 1, 3), (64., 1, 3), (64, 0, 3), (64, 1, 2)]:
        with pytest.raises(ValueError):
            topology.spectral_bounds(mass, a, radius)


def test_rational_source_encloses_independent_high_precision_directions():
    import sympy as s
    hr, hi, error = certification.rounded_source()
    exact = spectral.source()['h'].toarray()
    rounded = (hr.toarray()+1j*hi.toarray())/2**26
    assert la.norm(exact-rounded, 'fro')**2 <= float(error)*(1+1e-5)
    for point in ((1, 2, 1), (2, -1, 0), (0, 0, 2)):
        f = topology.direction(point)
        norm = s.sqrt(sum(x*x for x in f))
        for x in f:
            bound = certification.iv.divide((F(x), F(x)), certification.iv.sqrt_interval(sum(y*y for y in f)))
            assert s.N(s.Rational(bound[0])-x/norm, 80) <= 0
            assert s.N(s.Rational(bound[1])-x/norm, 80) >= 0


def test_rational_certificate_cannot_accept_missing_or_corrupt_basis():
    hr, hi, error = certification.rounded_source()
    n = hr.shape[0]
    with pytest.raises(ValueError, match='complete'):
        certify(hr, hi, np.zeros(n, dtype=np.int64), np.zeros((n, 4), dtype=np.int64),
                np.zeros((n, 4), dtype=np.int64), error, 26)
    with pytest.raises(ArithmeticError, match='invertible'):
        certify(hr, hi, np.zeros(n, dtype=np.int64), np.zeros((n, n), dtype=np.int64),
                np.zeros((n, n), dtype=np.int64), error, 26)


@pytest.mark.parametrize('a', [1, 2, 3])
def test_exact_leading_schur_cancellation(a):
    value = topology.leading_schur(a)
    expected = {'-1/4': 2, '1/4': 2} if a == 1 else {'0': 4}
    assert value['eigenvalue_multiplicities'] == expected
    assert value['all_single_core_diagonal_terms_cancel']
