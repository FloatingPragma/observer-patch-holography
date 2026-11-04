from fractions import Fraction as F

import numpy as np
import pytest
from scipy.sparse import diags, eye
from scipy.sparse.linalg import eigsh, norm

from . import bounds, finite


def test_exact_localization_schedule_and_small_spacing_domain():
    for a, v, m in ((F(1), F(1), F(1)), (F(1, 8), F(1), F(1)),
                    (F(1), F(1, 50), F(1)), (F(2, 3), F(5, 7), F(11, 13))):
        value = bounds.localization(a, v, m)
        radius = value['reference_radius_sites']
        assert radius >= 1 and radius >= 12*v/(a*m)
        assert value['exterior_squared_gap_lower'] >= m*m/2
        assert value['weighted_tail_prefactor_after_radius_2R'] <= 1
        assert 2*a*radius <= value['unperturbed_radius_2aR_upper']
        for distance in (0, 1, 3, 9):
            assert bounds.tail_bound(distance, a, v, m) == value['weight_per_site']**(-2*distance)
    assert bounds.localization(F(1, 8))['small_spacing_physical_decay_domain'] is True
    assert bounds.localization(1, F(1, 50))['small_spacing_physical_decay_domain'] is False
    assert bounds.localization(core_support_radius=20)['reference_radius_sites'] == 20


def test_exterior_compression_retains_boundary_terms_and_actual_gap():
    # A finite outer box has the same exterior-block identity; all matrix
    # operations below use literal original-source hopping and mass rows.
    a, v, m = F(1), F(1, 20), F(1)
    value = bounds.localization(a, v, m)
    radius = value['reference_radius_sites']
    points = np.array(finite.sites(3))
    kept = np.repeat(np.max(abs(points), axis=1) > radius, 2)
    h = finite.source_matrix(3, a, v, m)[kept][:, kept]
    mass = finite.source_matrix(3, a, v, m, hopping_fraction=0)[kept][:, kept]
    kinetic = h-mass
    anticommutator = kinetic@mass+mass@kinetic
    difference = h@h-kinetic@kinetic-float(m*m)*eye(h.shape[0])-anticommutator
    assert norm(difference) < 2e-14
    actual_squared_gap = eigsh(h@h, k=1, which='SA', return_eigenvectors=False,
                              v0=np.linspace(.1, 1., h.shape[0]), tol=1e-12)[0]
    assert actual_squared_gap >= float(value['exterior_squared_gap_lower'])-1e-12
    distance = np.repeat(np.maximum(0, np.max(abs(points), axis=1)-radius), 2)[kept]
    weights = float(value['weight_per_site'])**distance
    changed = diags(weights)@h@diags(1/weights)-h
    actual_weight_error = np.linalg.norm(changed.toarray(), 2)
    assert actual_weight_error <= float(value['weight_conjugation_error_upper'])+1e-12


def test_independent_complete_finite_zero_vectors_obey_the_tail_bound():
    a, v, m = F(1), F(1, 50), F(1)
    h = finite.source_matrix(4, a, v, m)
    # A deterministic block power filter avoids platform-dependent basis
    # choices for an exactly repeated complex eigenvalue. Both origin rows
    # and every intermediate original-source matrix row are retained.
    vectors = np.zeros((h.shape[0], 2), complex)
    origin = finite.sites(4).index((0, 0, 0))
    vectors[2*origin:2*origin+2] = np.eye(2)
    for _ in range(60):
        vectors -= h@(h@vectors)/2
    assert np.min(np.linalg.svd(vectors, compute_uv=False)) > .9
    vectors = np.linalg.qr(vectors)[0]
    assert np.linalg.norm(h@vectors) < 1e-12
    np.testing.assert_allclose(vectors.conj().T@vectors, np.eye(2), atol=2e-13, rtol=0)
    radius = bounds.localization(a, v, m)['reference_radius_sites']
    distances = np.repeat(np.max(abs(np.array(finite.sites(4))), axis=1), 2)
    for steps in (1, 2):
        tail = np.sum(abs(vectors[distances >= 2*radius+steps])**2, axis=0)
        assert np.max(tail) <= float(bounds.tail_bound(steps, a, v, m))
        assert np.max(tail) < 1e-7


def test_invalid_domains_and_nonlocal_perturbation_scope():
    for inputs in ((0, 1, 1), (True, 1, 1), (1., 1, 1), (1, -1, 1)):
        with pytest.raises(ValueError):
            bounds.localization(*inputs)
    for radius in (True, -1, 1.):
        with pytest.raises(ValueError):
            bounds.localization(core_support_radius=radius)
    for distance in (True, -1, 1.):
        with pytest.raises(ValueError):
            bounds.tail_bound(distance)
    assert bounds.certificate()['localization_extended_to_arbitrary_nonlocal_compact_perturbations'] is False
