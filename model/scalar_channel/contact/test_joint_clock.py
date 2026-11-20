"""Native reductions, complete pointers and no-fit joint record forecasts."""

import copy
import math

import numpy as np
import pytest
from scipy import sparse as sp
from scipy.linalg import expm

from ...contrast_memory_spherical_v1 import law as native
from ..cage import construction as cage
from . import joint_clock as clock
from .__main__ import HERE
from .joint_source import ENDPOINTS, MARK, MODES, arm_sources, radar_source, source_certificate
from .law import collar_states


@pytest.fixture(scope='module')
def protocol():
    return dict(endpoints=list(ENDPOINTS), mark=MARK, kappa=8, dt=.5, ticks=12,
                shots=1048576, seed=91731, joint_times=[.5, 2, 6],
                joint_signs=[1, -1], joint_phases=list(clock.PHASES),
                joint_drive_strength=64, couplings=[0, 1, -1],
                joint_familywise_alpha=1e-6, exact_forecast_tolerance=5e-11)


@pytest.fixture(scope='module')
def experiment(protocol):
    return clock.run(protocol)


def test_complete_native_memory_reduces_to_shared_source(protocol):
    D = radar_source(protocol, 0)
    scalar_eye = sp.eye(MODES, format='csr')
    memory_eye = sp.eye(native.DIM, format='csr')
    exchange, square = native.local_term_parts()
    marked = sp.csr_matrix(([1.], ([MARK], [MARK])), shape=(MODES, MODES))
    for sign in (1, -1):
        arms = collar_states(sign)
        embedding = sp.kron(sp.csr_matrix(arms), scalar_eye, format='csr')
        drive = arms[:, [0]] @ arms[:, [1]].T + arms[:, [1]] @ arms[:, [0]].T
        for coupling in (0, 1, -1):
            complete = (sp.kron(memory_eye, sp.csr_matrix(D))
                        + sp.kron(sp.csr_matrix(exchange + 8 * square + 64 * drive), scalar_eye)
                        + coupling * sp.kron(sp.csr_matrix(square), marked))
            H0, H1 = arm_sources(protocol, coupling, sign)
            reduced = sp.bmat([[sp.csr_matrix(H0), 64 * scalar_eye],
                              [64 * scalar_eye, sp.csr_matrix(H1)]], format='csr')
            residual = complete @ embedding - embedding @ reduced
            assert abs(residual).max() < 5e-13
    witness = source_certificate(protocol)
    assert witness['invariant_subspace_max_error'] < 1e-13
    assert witness['scalar_mapping_max_error'] < 1e-13
    assert witness['pointer_outcomes'] == 300
    assert abs(witness['removed_common_endpoint_energy'] - 48) < 1e-12


def test_marked_controls_preserve_all_six_written_stations(protocol):
    base = cage.law(cage.cell())
    for coupling in (0, 1, -1):
        expected = base.copy()
        expected[list(ENDPOINTS), list(ENDPOINTS)] += 1
        expected[MARK, MARK] += coupling
        np.testing.assert_allclose(radar_source(protocol, coupling), expected, atol=1e-13, rtol=0)
    for field, replacement in (('mark', ENDPOINTS[0]), ('kappa', 0), ('endpoints', list(ENDPOINTS)[::-1])):
        changed = {**protocol, field: replacement}
        with pytest.raises(ValueError):
            radar_source(changed, 1)
    with pytest.raises(ValueError):
        radar_source(protocol, True)


def test_joint_instrument_and_intervened_histories_are_complete(protocol, experiment):
    summary, arrays = experiment
    assert arrays['joint_probabilities'].shape == (3, 2, 7, 3, 4, 2, 150)
    for name, axes in (('joint', (-2, -1)), ('arm_path', (-3, -2, -1)),
                       ('pure_arm', (-2, -1)), ('dephased', (-2, -1)),
                       ('independent_clock', (-1,)), ('independent_scalar', (-1,))):
        np.testing.assert_allclose(arrays[name + '_probabilities'].sum(axis=axes), 1, atol=2e-12)
        np.testing.assert_array_equal(arrays[name + '_counts'].sum(axis=axes), protocol['shots'])
    np.testing.assert_array_equal(arrays['coarse_counts'], clock.coarse_grain(arrays['joint_counts']))
    assert np.min(arrays['coarse_probabilities'][..., -1].sum(axis=-1)) > 0
    assert summary['independent_forecast_max_error'] < protocol['exact_forecast_tolerance']
    assert summary['retained_arm_path_max_error'] < protocol['exact_forecast_tolerance']
    assert summary['four_phase_path_max_error'] < protocol['exact_forecast_tolerance']
    assert summary['sampling_within_bound']


def test_decoupled_source_factorizes_and_marked_source_has_joint_information(experiment):
    _, arrays = experiment
    p = arrays['joint_probabilities']
    product = p.sum(axis=-1, keepdims=True) * p.sum(axis=-2, keepdims=True)
    np.testing.assert_allclose(p[0], product[0], atol=2e-12, rtol=0)
    assert np.max(clock._tv(p[1], product[1])) > 1e-3
    assert np.max(clock._tv(p[2], product[2])) > 1e-3
    np.testing.assert_allclose(arrays['dephased_probabilities'],
                               np.repeat(arrays['arm_path_probabilities'].sum(axis=-3)[..., None, :, :],
                                         4, axis=-3), atol=1e-13)
    assert np.max(clock._tv(p[1], arrays['dephased_probabilities'][1])) > .05


def test_pure_arm_controls_use_the_identical_propagator(protocol, experiment):
    _, arrays = experiment
    for ci, coupling in enumerate(clock.COUPLINGS):
        for arm in (0, 1):
            H = radar_source(protocol, coupling if arm else 0)
            U = expm(-2j * H)
            expected = abs(U[:, clock.SEEDS].T) ** 2
            np.testing.assert_allclose(arrays['pure_arm_probabilities'][ci, 0, :, 1, arm, arm],
                                       expected, rtol=0, atol=2e-12)
            np.testing.assert_allclose(arrays['pure_arm_probabilities'][ci, 0, :, :, arm, 1-arm],
                                       0, rtol=0, atol=1e-14)


def test_source_is_retained_during_both_finite_pulses(protocol, experiment):
    _, arrays = experiment
    H0, H1 = arm_sources(protocol, 1, 1)
    I = np.eye(MODES)
    source_off_pulse = expm(-1j * math.pi / 4 * np.block([[0 * I, I], [I, 0 * I]]))
    seed = np.eye(2 * MODES)[:, MARK]
    initial = source_off_pulse @ seed
    waiting = np.r_[expm(-2j * H0) @ initial[:MODES], expm(-2j * H1) @ initial[MODES:]]
    actual = abs(source_off_pulse @ waiting) ** 2
    source_on = arrays['joint_probabilities'][1, 0, 6, 1, 0].ravel()
    assert .5 * abs(actual - source_on).sum() > 1e-4


def test_terminal_pointer_preserves_the_external_reference_on_average():
    rng = np.random.default_rng(751)
    joint = rng.normal(size=(2 * MODES, 3)) + 1j * rng.normal(size=(2 * MODES, 3))
    joint /= np.linalg.norm(joint)
    reference_before = joint.T @ joint.conj()
    # Each instrument branch retains its unnormalized quantum reference.
    # The orthogonal record can be copied to both observers without copying
    # an unknown quantum state of the scalar or the reference.
    retained = np.einsum('ir,is->irs', joint, joint.conj())
    np.testing.assert_allclose(retained.sum(axis=0), reference_before, atol=1e-14)
    assert abs(np.trace(retained.sum(axis=0)).real - 1) < 1e-14


def test_selection_uses_only_the_declared_training_fold(protocol, experiment):
    summary, arrays = experiment
    altered = dict(arrays)
    changed = arrays['joint_counts'].copy()
    # Reassign every held-out sixth-station count to one complete outcome.
    changed[1, :, 5] = 0
    changed[1, :, 5, :, :, 0, 0] = protocol['shots']
    altered['joint_counts'] = changed
    scored = clock._scores(protocol, altered, summary['simultaneous_sampling_total_variation_bound'])
    original = summary['held_out_model_selection']
    assert scored['folds']['selection'] == original['folds']['selection']
    assert scored['selected_models'] == original['selected_models']
    assert scored['folds']['sixth_station'] != original['folds']['sixth_station']
    assert set(scored['folds']['selection']['models']) == set(clock.MODELS)


def test_archive_analysis_rejects_dropped_records_and_false_forecasts(protocol, experiment):
    _, arrays = experiment
    for name in ('joint_counts', 'forecast_probabilities', 'dephased_probabilities',
                 'independent_clock_probabilities', 'coarse_counts'):
        changed = dict(arrays)
        changed[name] = arrays[name].copy()
        changed[name].flat[0] += 1
        with pytest.raises(ValueError):
            clock.analyze(protocol, changed)
    wrong = copy.deepcopy(protocol)
    wrong['joint_times'][0] = 1
    with pytest.raises(ValueError):
        clock.analyze(wrong, arrays)


def test_sampling_lattice_retains_all_probability_mass():
    p = np.array([[.25, .25, .25, .25], [1e-16, .6, .4 - 1e-16, 0]])
    weights = clock.canonical_sampling_weights(p, 1)
    np.testing.assert_array_equal(weights.sum(axis=1), clock.SAMPLING_LATTICE)
    assert np.max(.5 * abs(weights / clock.SAMPLING_LATTICE - p).sum(axis=1)) <= 4 / (2 * clock.SAMPLING_LATTICE)
    assert weights[1, 0] == 0
    # A tied maximum uses the first index as the fixed remainder recipient.
    three = np.full((1, 3), 1 / 3)
    np.testing.assert_array_equal(clock.canonical_sampling_weights(three, 1),
                                  [[333333333334, 333333333333, 333333333333]])


def test_archived_clock_authenticates_all_six_families_and_scores_original_counts(protocol):
    with np.load(HERE / 'clock_weights.npz', allow_pickle=False) as weights, \
         np.load(HERE / 'radar_counts.npz', allow_pickle=False) as counts:
        report, arrays = clock.run(protocol, acquisition=(counts, weights))
        validation = report['original_acquisition']
        assert report['sampling_within_bound']
        assert validation['weight_digest'] == '546d1450118a6d66397e7db5ad25914cf5174da3615252d518ea1557a2a64c8d'
        assert set(validation['families']) == {name for name, _ in clock.FAMILIES}
        assert sum(v['rows'] * v['bins_per_row'] for v in validation['families'].values()) == 530208
        for name, _ in clock.FAMILIES:
            np.testing.assert_array_equal(arrays[name + '_counts'], counts[name + '_counts'])
        accepted = clock.authenticated_weights(weights)
    for name in accepted:
        changed = accepted[name].copy()
        changed.flat[0] += 1
        changed.flat[1] -= 1
        with pytest.raises(ValueError, match='original integer clock weights differ'):
            clock.authenticated_weights(dict(accepted, **{name: changed}))
    name = 'joint_sampling_weights'
    with pytest.raises(ValueError, match='family'):
        clock.authenticated_weights({k: v for k, v in accepted.items() if k != name})
    with pytest.raises(ValueError, match='signed int64'):
        clock.authenticated_weights(dict(accepted, **{name: accepted[name].astype(float)}))
    # Shape is checked separately because the original digest commits only to names and bytes.
    with pytest.raises(ValueError, match='shape'):
        clock.analyze(protocol, arrays, acquisition_weights=dict(accepted, **{name: accepted[name].ravel()}))
    bad_counts = arrays['independent_clock_counts'].copy()
    row = bad_counts.reshape(-1, 2)[0]
    rare = row.argmin()
    row[:] = 0
    row[rare] = protocol['shots']
    with pytest.raises(ValueError, match='archived counts fail'):
        clock.analyze(protocol, dict(arrays, independent_clock_counts=bad_counts), acquisition_weights=accepted)


def test_frozen_clock_weights_tolerate_rounding_boundary_crossings_and_reject_changed_source():
    p = np.array([[.40, .35, .25]])
    delta = np.array([[1, -1, 0]]) / clock.SAMPLING_LATTICE
    weights = clock.canonical_sampling_weights(p + .49 * delta, 1)
    check = clock.check_acquisition_law(p + .51 * delta, weights, 1)
    assert check['locally_requantized_bins_different'] == 2
    with pytest.raises(ValueError, match='tolerance'):
        clock.check_acquisition_law(p + 5 * delta, weights, 1)
    with pytest.raises(ValueError, match='total'):
        clock.check_acquisition_law(p, weights + 1, 1)
