from fractions import Fraction as F
import math

import numpy as np
import pytest

from model.observer.source_scalar import compiler
from model.particle_species_detector import native_fft as f


def test_complete_3D_Fourier_plan_matches_all_axes_and_untouched_references():
    rng = np.random.default_rng(179)
    state = rng.normal(size=(1024, 5))+1j*rng.normal(size=(1024, 5))
    expected = np.fft.fftn(state.reshape(8, 8, 8, 2, 5), axes=(0, 1, 2), norm='ortho').reshape(1024, 5)
    np.testing.assert_allclose(f.apply(state), expected, atol=1e-14, rtol=0)
    np.testing.assert_allclose(f.apply(f.apply(state), inverse=True), state, atol=1e-14, rtol=0)


def test_native_butterflies_retain_vacuum_and_two_particle_scalar_phases():
    for phase_pi in (F(1), F(3, 4), F(1, 2), F(1, 4)):
        phase = np.exp(1j*math.pi*float(phase_pi))
        one = np.array([[1, -phase], [1, phase]])/math.sqrt(2)
        expected = np.eye(4, dtype=complex)
        expected[np.ix_([1, 2], [1, 2])] = one
        expected[3, 3] = phase
        for inverse in (False, True):
            actual = compiler.execute(f.butterfly_word(phase_pi, inverse), 2)
            np.testing.assert_allclose(actual, expected.conj().T if inverse else expected, atol=2e-14, rtol=0)


def test_every_explicit_bit_reversal_route_is_the_claimed_mode_transposition():
    for operation in f.plan():
        if operation.kind != 'swap': continue
        indices = list(range(1024))
        for left in f.adjacent_route(operation): indices[left], indices[left+1] = indices[left+1], indices[left]
        expected = list(range(1024)); p, q = operation.first, operation.second
        expected[p], expected[q] = expected[q], expected[p]
        assert indices == expected


def test_exact_complete_source_bound_and_all_output_dictionary_are_explicit():
    result = f.certificate()
    assert result['N1_dimension'] == 1024
    assert result['two_mode_butterflies'] == 4608
    assert result['occupation_pointer_complete_final_labels_each'] == 6
    assert F(result['preparation_plus_detection_complete_source_error']) < F(1, 10000)
    assert result['full_Fock_native_state_array_executed'] is False
    assert result['native_map_assumed_exactly_Gaussian_or_number_preserving'] is False
    assert result['exact_Fourier_factorization']['exact_unnormalized_entries_checked'] == 64
    assert result['exact_Fourier_factorization']['exact_Gram_entries_checked'] == 64


def test_truncated_or_masked_one_particle_inputs_are_rejected():
    with pytest.raises(ValueError): f.apply(np.zeros(1023))
    with pytest.raises(ValueError): f.apply(np.ma.array(np.zeros(1024)))
