"""Independent source, support, accounting and mutation controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
import sys

import numpy as np
import pytest

from . import certificate as c


def test_source_counts_are_actual_complete_endpoint_enumeration():
    independent = np.zeros((66, 66), dtype=np.int64)
    for i, state in enumerate(c.STATES):
        for seam in c.quantum_mixing.seams():
            for coin in (0, 1):
                final, _ = c.native.native(state, seam, coin)
                independent[i, c.INDEX[final]] += 1
    assert np.array_equal(independent, c.counts())
    assert np.all(independent.sum(axis=0) == 60)


def test_spectrum_retains_every_mode_and_independent_trace_check():
    value = c.transfer_spectrum()
    assert value['all_66_Newton_trace_coefficients_checked']
    assert sum(row['multiplicity'] for row in value['eigenvalue_intervals']) == 66
    assert value['distinct_eigenvalues'] == 16
    assert value['physical_masses_identified'] is False


def test_commutant_is_operator_equations_not_rank_matching():
    value = c.transfer_sectors()
    assert value['ranks'] == [1, 11, 54]
    assert value['matrix_entry_orbit_sizes'] == [66, 1320, 2970]
    assert value['nonzero_subspace_invariant_under_every_seam_requires_all_12_ports']
    assert not value['full_occupation_effect_algebra_preserves_these_sectors']


def test_every_port_subset_is_classified():
    value = c.compact_search()
    assert sum(row['core_count'] for row in value['classification']) == 4096
    assert sum(row['orbit_size'] for row in value['representatives']) == 4096
    assert all(row['dimension'] == 0 for row in value['classification'] if row['ports'] <= 9)
    assert sum(row['orbit_size'] for row in value['representatives']
               if len(row['core']) == 10 and row['compact_invariant_dimension']) == 6


def test_identity_hold_is_a_positive_control_for_support_algorithm():
    # Changing the dynamics to an identity hold preserves every core state.
    dimension, _ = c.compact_dimension_upper((0, 1, 2, 3), np.eye(66, dtype=np.int64))
    assert dimension == 6
    assert c.compact_dimension_upper((0, 1, 2, 3))[0] == 0


@pytest.mark.parametrize('core', [(1, 0), (0, 0), (True,), (-1,), (12,), [0, 1]])
def test_invalid_cores_are_rejected(core):
    with pytest.raises(ValueError, match='native port tuple'):
        c.compact_dimension_upper(core)


def test_exact_compact_modes_and_native_correlation_without_quantum_readout():
    core, modes = c.compact_modes()
    for eigenvalue, vector in modes:
        assert all(set(c.PAIRS[i]) <= set(core) for i, x in enumerate(vector) if x)
        propagated = vector.astype(object)
        for n in range(5):
            assert int(vector @ propagated) == 20*eigenvalue**n
            assert Q(int(vector @ propagated), 66*60**n) == Q(10, 33)*Q(eigenvalue, 60)**n
            propagated = c.counts().astype(object) @ propagated


def test_no_native_raw_outcome_is_dropped():
    rows = c.recorded_maps()
    assert len(rows) == 720
    assert sum(not row[-1] for row in rows) == 420
    for law in ('C', 'E'):
        counts = np.zeros(66, dtype=np.int64)
        for name, _, _, _, mapping in rows:
            if name == law:
                for source, _ in mapping:
                    counts[source] += 1
        assert np.all(counts == 60)


def test_full_records_destroy_transfer_protection_and_selector_breaks_confinement():
    value = c.mode_and_instrument_certificate()
    for row in value['modes']:
        assert Q(row['all_positive_nonuniform_schedule_control']['one_attempt_exterior_amplitude_norm_squared']) == Q(1, 80)
        for law in ('C', 'E'):
            assert row['full_recorded_instrument'][law]['complete_probability'] == '1'
            assert 0 < Q(row['full_recorded_instrument'][law]['correlation_sector_weight_after_one_attempt']) < 1
    assert value['ground_reference_exterior_trace_distance'] == '7/22'
    assert not value['coherent_mode_preparation_or_native_measurement_implemented']


def test_recorded_channel_result_independent_of_sector_quadratic_shortcut():
    _, modes = c.compact_modes()
    vector = modes[0][1]
    density_numerator = np.outer(vector, vector)
    for law in ('C', 'E'):
        output = np.zeros((66, 66), dtype=np.int64)
        for name, _, _, _, mapping in c.recorded_maps():
            if name == law and mapping:
                inputs, outputs = np.array(mapping, dtype=int).T
                output[np.ix_(outputs, outputs)] += density_numerator[np.ix_(inputs, inputs)]
        weight = Q(int(np.trace(c.sector_projectors()[2] @ output)), 330*60*20)
        expected = c.mode_and_instrument_certificate()['modes'][0]['full_recorded_instrument'][law]
        assert weight == Q(expected['correlation_sector_weight_after_one_attempt'])


def test_local_phase_has_same_exterior_positive_energy_and_complete_history_decay():
    value = c.local_phase_excitation()
    assert value['strict_positive_energy_excess_interval'] == ['1/1008', '1/672']
    assert value['full_11_port_exterior_exactly_equal']
    assert value['full_native_record_plus_exterior_equal_at_every_finite_history']
    assert value['native_decay_comparison']['C']['unread_difference_multiplier'] == '5/6'
    assert value['native_decay_comparison']['E']['unread_difference_multiplier'] == '7/10'
    assert not value['state_prepared_by_native_phase_control']


def test_local_energy_bound_agrees_with_independent_spectral_evaluation():
    matrix = c.counts().astype(float)/60
    eigenvalues, vectors = np.linalg.eigh(matrix)
    hamiltonian = (vectors * -np.log(eigenvalues)) @ vectors.T
    reference = matrix/np.trace(matrix)
    signs = 1-2*np.array(c.STATES)[:, 0]
    excited = signs[:, None]*reference*signs[None, :]
    energy = np.trace(hamiltonian @ (excited-reference))
    assert float(Q(1, 1008)) < energy < float(Q(1, 672))
    assert np.linalg.eigvalsh(excited).min() > 0
    assert np.allclose(np.diag(reference), np.diag(excited))


def test_local_phase_full_exterior_partial_trace_is_exact():
    matrix = c.counts()
    signs = 1-2*np.array(c.STATES)[:, 0]
    difference = signs[:, None]*matrix*signs[None, :]-matrix
    partial = {}
    for i, left in enumerate(c.STATES):
        for j, right in enumerate(c.STATES):
            if left[0] == right[0]:
                key = left[1:], right[1:]
                partial[key] = partial.get(key, 0)+int(difference[i, j])
    assert set(partial.values()) == {0}


def test_receipt_replay_and_scope():
    value = c.report()
    assert c.strict_json.load(c.HERE/'results.json') == value
    assert not value['scope']['noncompact_local_density_perturbations_excluded']
    assert not value['scope']['all_refinements_or_history_observers_exhausted']


@pytest.mark.parametrize('field', ['physical_particle_or_mass_identified', 'native_unitary_time_identified'])
def test_resealed_scope_promotion_is_rejected(tmp_path, monkeypatch, field):
    value = c.report()
    changed = deepcopy(value)
    changed['scope'][field] = True
    (tmp_path/'results.json').write_text(json.dumps(changed))
    monkeypatch.setattr(c, 'HERE', tmp_path)
    monkeypatch.setattr(sys, 'argv', ['native-defects', '--check'])
    with pytest.raises(ValueError, match='differs from exact replay'):
        c.main()
