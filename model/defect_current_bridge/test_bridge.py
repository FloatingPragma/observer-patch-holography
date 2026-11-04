"""Independent coefficient, carrier, symmetry and counterexample controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse as sp
import sympy as s

from model.common_source.chiral_candidate import action
from model.conditional_standard_model.spectrum import ActionInputs, mass_kernel
from . import algebra, construction as c, controls


def test_exact_local_identity_certificate():
    result = algebra.certificate()
    assert result['exact_dynamics_and_current_covariance_of_state_map'] is True
    assert result['normalized_trace_complement_is_vacuum'] is False


def test_independent_oriented_native_Higgs_convention():
    m = s.Matrix([[s.Rational(1, 101), s.Rational(2, 103), 0],
                  [s.Rational(2, 103), s.Rational(3, 107), 0], [0, 0, s.Rational(5, 109)]])
    inputs = ActionInputs(1, 1, 1, 1, 1,
                          s.diag(s.Rational(1, 6**6), s.Rational(1, 6**3), 1),
                          s.diag(s.Rational(1, 6**6), s.Rational(1, 6**4), s.Rational(1, 6**2)),
                          s.sqrt(2)*m, ('exact independent tensor-convention control',))
    actual = c.native_kernel(np.array(m, float)).toarray()
    expected = np.array(mass_kernel(inputs), complex)
    np.testing.assert_allclose(actual, expected, atol=1e-15, rtol=1e-15)


def test_native_representation_anti_Hermitian_convention():
    generators = action.generators()
    assert len(generators) == 12
    assert all(g.H == -g for g in generators)


@pytest.mark.parametrize('m', [0, F(1), F(61, 2)])
def test_rational_band_isolation(m):
    result = algebra.spectral_bound(m)
    assert F(result['band_separation_lower']) > 0
    assert result['full_many_body_band_isolation'] is False


@pytest.mark.parametrize('m,d', [(31, 62), (32, 62), (-1, 62), (1., 62), (1, 0)])
def test_band_bound_rejects_missing_domain(m, d):
    with pytest.raises(ValueError):
        algebra.spectral_bound(m, d)


@pytest.mark.parametrize('invalid', [np.zeros((2, 3)), [[0, 1], [0, 0]], [[float('nan')]]])
def test_pairing_requires_finite_symmetric_kernel(invalid):
    with pytest.raises(ValueError):
        c.assemble(invalid)


def test_complex_K_full_source_spectrum_and_broken_stagger_control():
    result = controls.independent_control()
    assert result['independent_full_Nambu_rows'] == 216
    assert result['finite_diagonalization_is_an_exact_interval_certificate'] is False


def test_zero_kernel_preserves_all_zero_directions():
    value = c.assemble(np.zeros((2, 2)))
    assert value['core'].shape == (8, 8)
    assert value['core'].nnz == 0
    assert c.norm(value['bdg']@value['inclusion']) < 1e-11


def test_localized_state_site_distribution_and_position():
    _, w, _, source = c.radial()
    spin = np.array([1, 2j])/np.sqrt(5)
    state = (w@spin).reshape(-1, 2)
    density = np.sum(abs(state)**2, axis=1)
    np.testing.assert_allclose(density, source['density'], atol=1e-15)
    np.testing.assert_allclose(density@np.array(source['points']), 0, atol=1e-15)
    assert sum(density) == pytest.approx(1)


def test_current_labels_without_native_action_do_not_intertwine():
    value = c.assemble(c.native_kernel(np.diag([.01, .02, .03])))
    g = np.array(action.generators()[0], complex)
    wrong = g.copy()
    wrong[0, 0] += 1j
    full = sp.kron(sp.eye(54), sp.csr_matrix(g))
    small = sp.kron(sp.eye(2), sp.csr_matrix(wrong))
    w = sp.kron(value['W'], sp.eye(45))
    assert c.norm(full@w-w@small) > 1


def test_no_unstated_quark_sector_when_disabled():
    kernel = c.native_kernel(np.diag([.01, .02, .03]), quark_design=False)
    assert np.linalg.matrix_rank(kernel.toarray(), tol=1e-12) == 6


def test_record_pairing_linearity_and_each_record_core_invariance():
    from model.common_source import native_family_record_memory as memory
    projectors = memory._projectors(False)
    weights = [(2, .02), (3, .03), (10, .05)]
    mean = sum((weight*np.array(projectors[index], float) for index, weight in weights), np.zeros((3, 3)))
    actual = c.native_kernel(mean, quark_design=False)
    expected = sp.csr_matrix((45, 45), dtype=complex)
    for index, weight in weights:
        k = c.native_kernel(np.array(projectors[index], float), quark_design=False)
        expected += weight*k
        value = c.assemble(k)
        assert c.residuals(value)['complete_Nambu_intertwiner'] < 1e-10
    assert c.norm(actual-expected) < 1e-14
    service = algebra.pulse_service_certificate()
    assert service['each_record_preserves_same_defect_CAR_subalgebra'] is True
    assert service['uniform_in_volume_continuum_bound'] is False
