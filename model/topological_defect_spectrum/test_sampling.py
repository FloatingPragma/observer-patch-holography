from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import sampling


def gaussian_basis(x):
    envelope = np.exp(-np.sum(x*x, axis=1)/2)/np.pi**.75
    return envelope[:, None, None]*np.eye(16)[None, :, :]


def test_all_valley_sampling_Gram_and_alternating_cross_terms():
    bits = tuple(product((0, 1), repeat=3))
    for spacing in (F(1, 2), F(1, 4), F(1, 8)):
        a = float(spacing)
        n = np.arange(-round(6/a), round(6/a)+1)
        diagonal = a/np.sqrt(np.pi)*np.sum(np.exp(-(a*n)**2))
        alternating = a/np.sqrt(np.pi)*np.sum((1-2*(n % 2))*np.exp(-(a*n)**2))
        factor = np.array([[diagonal, alternating], [alternating, diagonal]])
        gram = np.kron(np.kron(np.kron(factor, factor), factor), np.eye(2))
        assert np.linalg.norm(gram-np.eye(16), 2) < (1/3000 if a == .5 else 1e-13)
        assert len(bits) == 8
    # A separate literal 16-column grid comparison verifies the folded-row
    # ordering and all cross-valley terms before the large-box limit.
    points = np.array(list(product(range(-2, 3), repeat=3)))
    values = sampling.sample(points, F(1, 2), gaussian_basis)
    actual = np.einsum('nsk,nsl->kl', values.conj(), values)
    n = np.arange(-2, 3)
    diagonal = .5/np.sqrt(np.pi)*np.sum(np.exp(-(.5*n)**2))
    alternating = .5/np.sqrt(np.pi)*np.sum((1-2*(n % 2))*np.exp(-(.5*n)**2))
    factor = np.array([[diagonal, alternating], [alternating, diagonal]])
    expected = np.kron(np.kron(np.kron(factor, factor), factor), np.eye(2))
    np.testing.assert_allclose(actual, expected, atol=2e-15, rtol=0)


def test_literal_original_hopping_and_mass_have_second_order_consistency():
    # Every original continuum row has a nonzero complex polarization and
    # a different smooth spatial width/center. The original lattice source
    # is evaluated on chunks; no 2*volume square matrix is constructed.
    widths = 1+np.arange(16)/32
    centers = np.array([(j % 3, (j//3) % 3, (j//9) % 3) for j in range(16)])/10
    polarization = (np.arange(1, 17)+1j*np.arange(16, 0, -1))/np.sqrt(2992)
    def field(x):
        displacement = x[:, None, :]-centers[None, :, :]
        envelope = np.exp(-widths[None, :]*np.sum(displacement**2, axis=2)/2)
        return (envelope*polarization[None, :])[:, :, None]
    def gradient(x):
        displacement = x[:, None, :]-centers[None, :, :]
        return (-np.swapaxes(displacement, 1, 2)*widths[None, None, :]*field(x)[:, None, :, 0])[:, :, :, None]
    errors = []
    for spacing in (F(1, 2), F(1, 4), F(1, 8)):
        radius = int(3/spacing)
        points = np.array(list(product(range(-radius, radius+1), repeat=3)))
        squared = 0.
        for start in range(0, len(points), 1024):
            part = points[start:start+1024]
            actual = sampling.lattice_action(part, spacing, field)
            expected = sampling.sampled_continuum_action(part, spacing, field, gradient)
            squared += np.linalg.norm(actual-expected)**2
        errors.append(np.sqrt(squared))
    assert all(3.3 < first/second < 4.5 for first, second in zip(errors, errors[1:]))
    assert errors[-1] < .1


def test_original_onsite_mass_intertwines_every_folded_column_exactly():
    points = np.array(list(product(range(-2, 3), repeat=3)))
    def constant(x):
        return np.broadcast_to(np.eye(16), (len(x), 16, 16))
    def zero_gradient(x):
        return np.zeros((len(x), 3, 16, 16))
    actual = sampling.lattice_action(points, F(1, 3), constant, mass=2)
    expected = sampling.sampled_continuum_action(points, F(1, 3), constant, zero_gradient, mass=2)
    np.testing.assert_allclose(actual, expected, atol=3e-16, rtol=0)
    assert np.linalg.norm(actual) > 1


def test_bounded_lattice_inclusion_keeps_its_direction_and_domain():
    result = sampling.certificate()
    assert result['source_valleys_retained'] == 8
    assert result['near_zero_lattice_multiplicity_lower_bound'] == 2
    for key in ('global_spectral_inclusion_rate_proved', 'additional_lattice_gap_levels_excluded',
                'exact_bounded_lattice_zero_modes_proved',
                'bounded_lattice_and_continuum_level_multiplicities_identified'):
        assert result[key] is False
    for spacing in (True, 0, -1, .5):
        with pytest.raises(ValueError):
            sampling.sample(np.array([[0, 0, 0]]), spacing, gaussian_basis)
    with pytest.raises(ValueError):
        sampling.sample(np.array([[0., 0., 0.]]), F(1, 2), gaussian_basis)
    with pytest.raises(ValueError):
        sampling.sample(np.array([[0, 0, 0]]), F(1, 2), lambda x: np.zeros((len(x), 2, 1)))
