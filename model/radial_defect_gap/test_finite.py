import numpy as np
import pytest
from scipy.linalg import eigh, solve

from . import bounds, finite


@pytest.mark.parametrize('mass', [3, 4, 64])
def test_complete_original_spectrum_and_origin_schur_identity(mass):
    value = finite.kernel_frame(1, mass=mass)
    h = value['H'].toarray()
    d, v = value['D'].toarray(), value['V']
    gap = float(value['estimate']['one_body_gap_lower'])
    energies = np.linalg.eigvalsh(h)
    assert np.count_nonzero(abs(energies) < 1e-12) == 2
    assert np.min(abs(energies[abs(energies) > 1e-12])) >= gap-1e-12
    np.testing.assert_allclose(v.conj().T@v, 1.5*np.eye(2), atol=1e-15, rtol=0)
    assert np.linalg.norm(v.conj().T@solve(d, v, assume_a='her')) < 1e-13
    assert np.linalg.norm(h@value['frame']) < 1e-12
    np.testing.assert_allclose(value['frame'].conj().T@value['frame'], np.eye(2), atol=1e-14, rtol=0)
    for energy in (-gap/2, gap/2):
        inverse_v = solve(d-energy*np.eye(len(d)), v, assume_a='her')
        zero_inverse_v = solve(d, v, assume_a='her')
        bracket = np.eye(2)+v.conj().T@solve(d-energy*np.eye(len(d)), zero_inverse_v, assume_a='her')
        schur = -energy*np.eye(2)-v.conj().T@inverse_v
        np.testing.assert_allclose(schur, -energy*bracket, atol=2e-13, rtol=0)
        assert np.min(np.linalg.eigvalsh(bracket)) >= 1-1e-13


@pytest.mark.parametrize('sign', [-1, 1])
def test_every_local_zero_frame_and_filled_sea_origin_weight(sign):
    value = finite.kernel_frame(2, sign=sign)
    frame = value['frame'].reshape(-1, 2, 2)
    weights = value['density']
    np.testing.assert_allclose(frame[::-1], frame, atol=2e-14, rtol=0)
    for row, weight in zip(frame, weights, strict=True):
        np.testing.assert_allclose(row.conj().T@row, weight*np.eye(2), atol=2e-14, rtol=0)
        np.testing.assert_allclose(row@row.conj().T, weight*np.eye(2), atol=2e-14, rtol=0)
        assert abs(abs(np.linalg.det(row))**2-weight**2) < 1e-14
    assert abs(weights.sum()-1) < 1e-14
    origin = value['points'].index((0, 0, 0))
    rho = weights[origin]
    assert rho >= float(value['estimate']['origin_weight_lower'])
    assert np.sum(weights**2) >= float(value['estimate']['ipr_lower'])
    energies, vectors = eigh(value['H'].toarray())
    negative = vectors[:, energies < -1]
    sea = negative@negative.conj().T
    rows = value['P_rows']
    np.testing.assert_allclose(sea[np.ix_(rows, rows)], (1-rho)*np.eye(2)/2, atol=2e-14, rtol=0)
    deviations = []
    for occupation in ((0, 0), (1, 0), (0, 1), (1, 1)):
        covariance = sea+value['frame']@np.diag(occupation)@value['frame'].conj().T
        local = covariance[np.ix_(rows, rows)]
        deviations.append((1-np.trace(local)+2*np.linalg.det(local)).real)
    np.testing.assert_allclose(deviations, [(1+rho*rho)/2, (1-rho*rho)/2,
                                          (1-rho*rho)/2, (1+rho*rho)/2], atol=2e-14, rtol=0)


def test_local_walk_truncation_and_exact_support_prefix():
    value = finite.kernel_frame(2)
    for terms in (0, 1, 2, 3, 6):
        approximation = finite.neumann_frame(2, terms)
        error_squared = np.linalg.norm(approximation-value['raw_frame'], 2)**2
        assert error_squared <= float(bounds.neumann_error_squared(terms))+1e-14
    for radius in (1, 2):
        small = finite.neumann_frame(radius, radius)
        big = finite.neumann_frame(radius+1, radius)
        np.testing.assert_allclose(finite.embed(radius, radius+1, small), big, atol=1e-15, rtol=0)
    beyond_small = finite.embed(1, 2, finite.neumann_frame(1, 2))
    beyond_big = finite.neumann_frame(2, 2)
    assert np.linalg.norm(beyond_small-beyond_big) > .01


@pytest.mark.parametrize('mass', [4, 64])
def test_volume_projector_convergence_and_localization_bound(mass):
    small, big = finite.kernel_frame(1, mass=mass), finite.kernel_frame(2, mass=mass)
    embedded = finite.embed(1, 2, small['frame'])
    actual = np.linalg.norm(embedded@embedded.conj().T-big['frame']@big['frame'].conj().T, 2)**2
    assert actual <= float(bounds.finite_volume_projector_error_squared(1, mass=mass))+1e-14
    distances = np.max(abs(np.array(big['points'])), axis=1)
    for distance in (1, 2):
        probability = big['density'][distances >= distance].sum()
        assert probability <= float(bounds.tail_probability_upper(distance, mass=mass))+1e-14
    assert big['density'][distances >= 1].sum() <= 1-float(big['estimate']['origin_weight_lower'])+1e-14


def test_countercontrols_reject_an_unproved_strong_mass_domain():
    with pytest.raises(ValueError):
        finite.kernel_frame(1, mass=1)
    with pytest.raises(ValueError):
        finite.neumann_frame(1, True)
    with pytest.raises(ValueError):
        finite.kernel_frame(0)
