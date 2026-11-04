"""Independent source, topology and spectral-exhaustion controls."""
import numpy as np
import pytest
import sympy as s

from ..bounded_defect_index import finite as original
from . import algebra, search


def test_complete_callias_intertwiner_and_degree_sign():
    value = algebra.index_certificate()
    assert value['all_six_generators_and_grading_intertwined_exactly']
    assert value['Callias_4row_blocks_in_positive_to_negative_grading'] == 2
    for d in (-4, -3, -1, 1, 2, 3, 4):
        assert algebra.degree_certificate(d)['source_continuum_chiral_index'] == 2*d


def test_rational_map_matches_original_hedgehog_and_antipodal_parity():
    for point in original.sites(2):
        value = search.direction(point, 1)
        expected = np.array(point)/np.linalg.norm(point) if point != (0, 0, 0) else np.zeros(3)
        np.testing.assert_allclose(value, expected, atol=2e-15)
        for d in (1, 2, 3, 4):
            n = search.direction(point, d)
            np.testing.assert_allclose(search.direction(tuple(-x for x in point), d),
                                       [(-1)**d*n[0], (-1)**d*n[1], -n[2]], atol=3e-15)
            if point != (0, 0, 0):
                assert abs(np.linalg.norm(n)-1) < 3e-15


def test_literal_links_and_all_sites_match_parent_for_degree_one():
    for radius in (1, 2):
        np.testing.assert_allclose(search.source_matrix(radius, 1, mass=4).toarray(),
                                   original.source_matrix(radius, mass=4).toarray(), atol=2e-15)


@pytest.mark.parametrize('degree', [1, 2, 3, 4])
def test_opposite_degree_keeps_all_rows_and_particle_hole_spectrum(degree):
    points = original.sites(1)
    ds = np.repeat([(-1)**(sum(p) % 2) for p in points], 2)
    h = search.source_matrix(1, degree, mass=4).toarray()
    opposite = search.source_matrix(1, -degree, mass=4).toarray()
    np.testing.assert_allclose(opposite, -ds[:, None]*h*ds[None, :], atol=1e-15)
    np.testing.assert_allclose(np.linalg.eigvalsh(h), np.linalg.eigvalsh(opposite), atol=2e-14)
    assert h.shape == (54, 54)


@pytest.mark.parametrize('degree', [-3, -1, 1, 3, 5])
def test_strong_mass_index_does_not_count_winding_as_families(degree):
    result = algebra.strong_mass_certificate(degree, mass=4)
    assert result['bounds']['kernel_dimension_exact'] == 2
    finite = search.finite_spectrum(1, degree, mass=4)
    assert len(finite['all_eigenvalues']) == finite['dimension'] == 54
    assert finite['numerical_zero_count'] == 2
    assert finite['window_rank'] == 2


def test_invalid_degree_and_even_degree_inversion_claim_rejected():
    for d in (0, True, 1.0, s.Integer(1)):
        with pytest.raises(ValueError):
            algebra.degree_certificate(d)
    with pytest.raises(ValueError, match='odd degree'):
        algebra.strong_mass_certificate(2, mass=4)
