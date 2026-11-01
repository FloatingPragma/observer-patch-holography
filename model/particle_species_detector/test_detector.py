from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest

from ..periodic_momentum.momentum import PeriodicGeometry, one_particle
from . import detector as d, native_fft


def test_integer_generator_is_bound_to_the_existing_original_source():
    real, imag, _ = d.integer_operators(4)
    expected = one_particle(PeriodicGeometry((4, 4, 4)))
    np.testing.assert_array_equal((real+1j*imag).toarray(), 16*np.asarray(expected, complex))


def test_exact_orthonormal_spectral_modes_and_unresolved_complement():
    result = d.exact_species_certificate()
    assert result['identified_orthonormal_modes'] == 16
    assert result['distinct_operational_Weyl_flavors'] == 8
    assert result['unresolved_single_particle_outcomes'] == 1008
    assert not result['z_current_commutes_with_H_on_all_N1_states']
    assert result['z_current_stationary_on_selected_modes']


def test_all_eight_flavors_have_two_opposite_energy_and_helicity_branches():
    catalog = d.mode_catalog()
    assert len({row['detector_address'] for row in catalog}) == 16
    for flavor in sorted({row['flavor'] for row in catalog}):
        rows = [row for row in catalog if row['flavor'] == flavor]
        assert {row['energy_sqrt2_coefficient'] for row in rows} == {'-1/16', '1/16'}
        assert {row['helicity'] for row in rows} == {'-1/2', '1/2'}
        assert all(row['model_number'] == 1 and row['CAR_parity'] == -1 for row in rows)
        assert all(row['continuum_mass'] == '0' for row in rows)
        assert not any(row['electric_charge_assigned'] or row['physical_species_name_assigned'] for row in rows)
    positive = [row for row in catalog if row['energy_branch'] == '+']
    assert sum(row['helicity'] == '1/2' for row in positive) == 4


@pytest.mark.parametrize('time', [F(), F(1, 2), F(8), F(-8)])
def test_actual_original_H_propagation_and_complete_detector_outputs(time):
    result = d.simulate_detector(time)
    probabilities = np.array(result['output_probabilities_all_1024_rows'])
    assert probabilities.shape == (1024, 16)
    assert np.max(abs(probabilities.sum(axis=0)-1)) < 1e-13
    np.testing.assert_allclose(result['identified_outcome_confusion_matrix'], np.eye(16), atol=1e-13, rtol=0)
    assert result['complete_ideal_output_vector_difference_diagnostic'] < 2e-13
    assert not result['probabilities_normalized_or_postselected']


def test_wrong_Fourier_direction_fails_all_chosen_momentum_labels():
    vectors = d.exact_mode_vectors()
    correct = abs(native_fft.apply(vectors))**2
    wrong = abs(native_fft.apply(vectors, inverse=True))**2
    addresses = [row['detector_address'] for row in d.mode_catalog()]
    assert np.min(correct[addresses, np.arange(16)]) > 1-1e-13
    assert np.max(wrong[addresses, np.arange(16)]) < 1e-25


def test_no_generic_energy_assignment_is_made_to_spin_z_outcomes():
    # A momentum with a nonzero x component mixes spin-z: this is one of
    # the retained1008 outcomes, not silently assigned the family formula.
    basis = np.zeros(1024, complex)
    label = d.address((1, 0, 0), 0)
    basis[label] = 1
    vector = native_fft.apply(basis, inverse=True)
    h_vector = d.original_hamiltonian()@vector
    assert abs(np.vdot(vector, h_vector)) < 1e-14
    assert np.linalg.norm(h_vector) > .08
    assert label not in {row['detector_address'] for row in d.mode_catalog()}


def test_complete_seeded_count_illustration_and_no_postselection():
    result = d.sample_counts()
    assert sum(result['counts_all_1024_single_particle_addresses']) == 65536
    assert all(row['count'] > 0 for row in result['counts_by_identified_mode'])
    assert F(result['simultaneous_failure_probability_upper']) < F(1, 10**11)
    assert result['counts_are_ideal_model_simulation']
    assert not result['conditioned_on_one_particle_or_success']


def test_finite_periodic_refinement_energies_are_enclosed_and_converge():
    table = d.refinement_table()
    with mp.workdps(100):
        previous = mp.mpf(0)
        for row in table:
            a = F(row['spacing'])
            convert = lambda value: mp.mpf(F(value).numerator)/F(value).denominator
            exact = mp.sin(convert(a)*mp.pi/4)/(8*convert(a))
            lo, hi = map(convert, row['positive_energy_interval'])
            assert lo <= exact <= hi
            assert previous < exact < mp.pi/32
            previous = exact
            assert row['side']*a == 8
            assert row['continuum_mass'] == '0'


@pytest.mark.parametrize('index', [0, 1, 15, 511, 1023])
def test_detector_address_round_trip(index):
    site, spin = d.decode(index)
    assert d.address(site, spin) == index


def test_invalid_input_domains_are_rejected():
    calls = [lambda: d.address((0, 0, 8), 0), lambda: d.address((0, 0), 0),
             lambda: d.decode(True), lambda: d.decode(-1), lambda: d.integer_operators(2),
             lambda: d.fft_reference(np.zeros(1023)), lambda: d.simulate_detector(.5),
             lambda: d.sample_counts(trials=100), lambda: d.sample_counts(seed=True)]
    for call in calls:
        with pytest.raises(ValueError): call()
