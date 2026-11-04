"""Independent source, spectral-instrument and localization controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import mpmath as mp
import numpy as np
import pytest
from scipy import sparse

from ..particle_species_detector import detector, native_fft
from ..periodic_momentum.momentum import PeriodicGeometry, one_particle
from . import native, spectrum as search


def test_exact_source_matches_existing_periodic_action():
    real, imag, _ = detector.integer_operators(4)
    original = np.asarray(one_particle(PeriodicGeometry((4, 4, 4))), complex)
    np.testing.assert_array_equal((real+1j*imag).toarray(), 16*original)
    # Fourier transform of the actual oriented source link i*sigma/16.
    import sympy as s
    root = (1+s.I)/s.sqrt(2)
    for k, (a, b) in enumerate(search.TWICE_SINE):
        assert s.simplify(s.I*(root**k-root**(-k))+(a+b*s.sqrt(2))) == 0


def test_all_exact_multiplicities_and_contact_annihilation():
    result = search.exact_certificate()
    assert result['complete_count'] == 1024
    assert [r['multiplicity_per_energy'] for r in result['shells']] == [16, 48, 120, 160, 120, 48, 8]
    assert sum(r['signed_branches']*r['multiplicity_per_energy'] for r in result['shells']) == 1024
    # Direct occupation addresses for every one-particle mode; no spin or
    # position output is omitted. A local doublon is a strict countercontrol.
    for occupied in range(1024):
        state = 1 << occupied
        assert sum(bool(state & (1 << (2*x))) and bool(state & (1 << (2*x+1)))
                   for x in range(512)) == 0
    assert (3 & 1) and (3 & 2)


def test_complete_detector_diagonalizes_every_original_source_column():
    columns = search.apply(np.eye(1024), inverse=True)
    np.testing.assert_allclose(detector.original_hamiltonian()@columns,
                               columns*search.energy_labels(), atol=2e-14, rtol=0)
    law = abs(search.apply(columns))**2
    np.testing.assert_allclose(law, np.eye(1024), atol=2e-13, rtol=0)
    np.testing.assert_allclose(law.sum(axis=0), np.ones(1024), atol=2e-13, rtol=0)


def test_omitting_spin_rotation_misidentifies_generic_energy_modes():
    label = detector.address((1, 0, 0), 0)
    basis = np.zeros(1024); basis[label] = 1
    energy_mode = search.apply(basis, inverse=True)
    incomplete = abs(native_fft.apply(energy_mode))**2
    assert incomplete[label] == pytest.approx(.5)
    assert incomplete[label+1] == pytest.approx(.5)
    complete = abs(search.apply(energy_mode))**2
    assert complete[label] == pytest.approx(1)


def test_local_projector_bound_and_saturating_state():
    # Build spectral projectors from an independent dense real-space eigensolver.
    h = detector.original_hamiltonian().toarray()
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    for n in range(7):
        energy = np.sqrt(2*n)/16
        selected = abs(eigenvalues-energy) < 1e-12
        v = eigenvectors[:, selected]
        count = [16, 48, 120, 160, 120, 48, 8][n]
        assert v.shape[1] == count
        for site in (0, 1, 127, 511):
            block = v[2*site:2*site+2]
            np.testing.assert_allclose(block@block.conjugate().T,
                                       count/1024*np.eye(2), atol=1e-13, rtol=0)
        state = v@v[0].conjugate()/np.sqrt(count/1024)
        assert np.vdot(state[:2], state[:2]).real == pytest.approx(count/1024, abs=1e-13)
        np.testing.assert_allclose(h@state, energy*state, atol=1e-13, rtol=0)
    theorem = search.localization_certificate()
    assert F(theorem['stationary_one_site_occupancy_upper']) == F(5, 32)
    assert F(theorem['six_site_exterior_occupancy_lower']) == F(1, 16)


def test_localized_input_is_not_stationary_and_zero_speed_is_counterexample():
    local = np.zeros(1024); local[0] = 1
    h = detector.original_hamiltonian()
    assert np.linalg.norm(h@local) > .1
    rho = sparse.csr_matrix(np.outer(local, local))
    assert sparse.linalg.norm(h@rho-rho@h) > .1
    np.testing.assert_array_equal((0*h)@local, np.zeros(1024))
    proof = search.laurent_certificate()
    assert proof['zero_energy_included'] is True
    assert proof['c_zero_excluded_from_no_go'] is True


def test_spin_rotation_enclosures_against_high_precision_radicals():
    with mp.workdps(100):
        for k in product(range(8), repeat=3):
            n = search.shell(k)
            x, y, z = (mp.mpf(a)+b*mp.sqrt(2) for a, b in search.direction(k))
            if not n or (x == y == 0 and z > 0):
                expected = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
            elif x == y == 0:
                expected = (mp.mpf(0), mp.mpf(1), mp.mpf(0))
            else:
                e = mp.sqrt(2*n); norm = mp.sqrt(2*e*(e+z))
                expected = ((e+z)/norm, x/norm, -y/norm)
            for interval, value in zip(native.target_entries(k), expected, strict=True):
                assert mp.mpf(interval.lo)/native.iv.SCALE <= value <= mp.mpf(interval.hi)/native.iv.SCALE


def test_complete_native_dilation_and_error_cost():
    result = native.certificate()
    assert result['spin_rotations'] == 492
    assert result['spin_rotation_native_pulses'] == 12792
    assert result['detector_native_pulses'] == 32925688
    assert F(result['complete_detector_operator_error']) < F(4155, 10**12)
    assert F(result['complete_output_unhalved_L1_error_upper']) == 2*F(result['complete_detector_operator_error'])
    assert native.local_execution()['distinct_complete_local_blocks_executed'] == 125
    assert result['higher_N_outcomes_are_interacting_energy_labels'] is False
    assert result['complete_source_error_projected_to_N1'] is False


def test_corrupted_native_angle_exceeds_certificate(monkeypatch):
    original = native.rotation_word
    def wrong(k):
        word, alpha, time = original(k)
        if word:
            word = tuple({**event, 'angle': event['angle']+.01} if index == 0 else event
                         for index, event in enumerate(word))
        return word, alpha, time
    native.rotation_certificate.cache_clear()
    monkeypatch.setattr(native, 'rotation_word', wrong)
    try:
        with pytest.raises(ArithmeticError, match='declared ceiling'):
            native.certificate()
    finally:
        native.rotation_certificate.cache_clear()


def test_record_rejects_scope_type_and_multiplicity_mutations():
    from . import packet
    value = packet.build()
    for key, replacement in (('observed_species_identifications', 1),
                             ('source_number_sector', True),
                             ('entire_completion_gate_accepted', True)):
        bad = deepcopy(value)
        bad['scope'][key] = replacement
        with pytest.raises(ValueError, match='differs from replay'):
            packet.validate(bad, value)
    bad = deepcopy(value)
    bad['exact_complete_spectrum']['shells'][3]['multiplicity_per_energy'] -= 1
    with pytest.raises(ValueError, match='differs from replay'):
        packet.validate(bad, value)


@pytest.mark.parametrize('momentum', [(0, 0), (0, 0, 8), (True, 0, 0), (0., 0, 0), '000'])
def test_invalid_momentum_rejected(momentum):
    with pytest.raises(ValueError):
        search.rotation(momentum)


@pytest.mark.parametrize('value', [np.zeros(1023), np.zeros((1024, 0)),
                                   np.full(1024, np.nan), np.ma.array(np.zeros(1024))])
def test_invalid_detector_input_rejected(value):
    with pytest.raises(ValueError):
        search.apply(value)
