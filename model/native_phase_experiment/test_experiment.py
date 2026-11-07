"""Independent source and receipt regression tests for the phase experiment."""
from copy import deepcopy
from fractions import Fraction
import subprocess
import sys

import numpy as np
import pytest

from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import experiment
from .independent import compute
from .verify import verify, _close
from .provenance import source_graph, source_hashes
from .__main__ import _persist


@pytest.fixture(scope='module')
def executed():
    return experiment.report(keep_arrays=True)


@pytest.fixture(scope='module')
def expected():
    return compute(experiment.ROOT, compact=False)


def test_complete_source_execution_matches_independent_counts(executed, expected):
    result = verify(executed[0], expected)
    assert result['all_pass']
    assert result['raw_records_compared'] == 1440
    assert result['maximum_raw_readout_error'] < 1e-12


def test_source_preparation_and_phase_retain_full_occupation_register(executed):
    result, arrays = executed
    assert arrays['prepared'].shape == (4096, 366)
    assert arrays['twisted'].shape == (4096, 366)
    for setting in ('reference', 'twisted'):
        assert result['preparation']['densities'][setting]['outside_total_two_probability'] < 1e-24
    runs = experiment.source_runs()
    pair = np.zeros(8)
    pair[[5, 6]] = -1 / np.sqrt(2)
    np.testing.assert_allclose(runs['pair'].success[:, 0], pair, atol=2e-12, rtol=0)
    np.testing.assert_allclose(runs['twisted_phase'].success, np.diag([1., -1.]), atol=2e-12, rtol=0)


def test_selector_retains_exhaustion_probability(executed, expected):
    result = executed[0]
    timeout = float(Fraction(expected['finite_selector']['timeout_probability']))
    assert abs(result['selector']['timeout_probability'] - timeout) < 2e-20
    assert result['selector']['timeout_probability'] > 0
    assert len(result['selector']['slot_probabilities']) == 1680


@pytest.mark.parametrize('mutation', ('omit_raw', 'duplicate_raw', 'omit_pointer', 'change_score',
                                     'omit_timeout', 'change_density', 'change_source', 'nonfinite',
                                     'move_mass_to_abort'))
def test_independent_verifier_rejects_changed_output(executed, expected, mutation):
    changed = deepcopy(executed[0])
    if mutation == 'omit_raw':
        changed['trials'][1]['records'].pop()
    elif mutation == 'duplicate_raw':
        changed['trials'][1]['records'][-1] = deepcopy(changed['trials'][1]['records'][0])
    elif mutation == 'omit_pointer':
        changed['trials'][1]['records'][0]['outcomes'].pop()
    elif mutation == 'change_score':
        changed['pointer_records'][1]['score'] = 1
    elif mutation == 'omit_timeout':
        changed['selector']['timeout_probability'] = 0.
    elif mutation == 'change_density':
        changed['preparation']['densities']['twisted']['real'][0][0] += 0.01
    elif mutation == 'change_source':
        changed['source_sha256']['model/observer/native.py'] = '0' * 64
    elif mutation == 'nonfinite':
        changed['trials'][0]['score'] = float('nan')
    else:
        changed['trials'][0]['records'][0]['outcomes'][0] -= 0.0001
        changed['trials'][0]['records'][0]['outcomes'][14] += 0.0001
    with pytest.raises((ValueError, ArithmeticError)):
        verify(changed, expected)


def test_exact_lineage_classifier_distinguishes_loss_and_hold(executed):
    lineage = executed[0]['protected_obstruction_reuse']
    assert lineage['C']['survival'] == '5/6'
    assert lineage['E']['survival'] == '7/10'
    assert lineage['C']['classification']['almostSure']
    assert lineage['E']['classification']['almostSure']
    assert lineage['hold_control']['closedTrap']


def test_archive_implicit_suffix_preserves_existing_file(tmp_path):
    archive = tmp_path / 'saved.npz'
    archive.write_bytes(b'existing archive')
    receipt = tmp_path / 'receipt.json'
    command = [sys.executable, '-B', '-m', 'model.native_phase_experiment',
               '--write', str(receipt), '--archive', str(tmp_path / 'saved')]
    result = subprocess.run(command, cwd=experiment.ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert archive.read_bytes() == b'existing archive'
    assert not receipt.exists()


def test_exclusive_archive_reservation_rolls_back_receipt(tmp_path):
    archive = tmp_path / 'saved.npz'
    archive.write_bytes(b'concurrent archive')
    receipt = tmp_path / 'receipt.json'
    with pytest.raises(FileExistsError):
        _persist({'result': True}, {'value': np.array([1.])}, receipt, archive)
    assert archive.read_bytes() == b'concurrent archive'
    assert not receipt.exists()


def test_exclusive_receipt_preserves_concurrent_file(tmp_path):
    receipt = tmp_path / 'receipt.json'
    receipt.write_bytes(b'concurrent receipt')
    with pytest.raises(FileExistsError):
        _persist({'result': True}, {}, receipt, None)
    assert receipt.read_bytes() == b'concurrent receipt'


def test_exact_rational_comparison_does_not_round_the_target_first():
    with pytest.raises(ValueError):
        _close(0.3, Fraction(3, 10), 'exact rational target', tolerance=1e-20)


def test_source_closure_includes_transitive_runtime_and_package_inputs():
    graph = source_graph()
    paths = set(graph['files'])
    assert 'model/native_phase_experiment/__init__.py' in paths
    assert 'model/native_phase_experiment/provenance.py' in paths
    assert 'model/reference_recovery/recovery.py' in paths
    assert 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/finite_groups.py' in paths
    assert set(source_hashes()) == paths


@pytest.mark.parametrize('mutation', (
    'component_weight', 'component_ports', 'component_slots', 'success_zero',
    'success_two', 'success_nan', 'success_missing', 'bit_distribution',
    'selector_diagnostic_nan', 'lineage_pair_count', 'lineage_kernel',
    'hold_extra_kernel', 'resource_missing', 'resource_zero', 'source_extra_run',
    'target_boolean_port', 'raw_boolean_coin'))
def test_complete_verifier_rejects_profile_metadata_changes(executed, expected, mutation):
    changed = deepcopy(executed[0])
    if mutation == 'component_weight':
        changed['preparation']['components'][0]['weight'] += 0.01
    elif mutation == 'component_ports':
        changed['preparation']['components'][0]['ports'][0] = 11
    elif mutation == 'component_slots':
        changed['preparation']['components'][0]['selector_slots'].clear()
    elif mutation.startswith('success_'):
        if mutation == 'success_missing':
            del changed['preparation']['success_probability']
        else:
            changed['preparation']['success_probability'] = {
                'success_zero': 0., 'success_two': 2., 'success_nan': float('nan')}[mutation]
    elif mutation == 'bit_distribution':
        changed['selector']['bit_record_probabilities'][0] += 0.1
    elif mutation == 'selector_diagnostic_nan':
        changed['selector']['normalization_error'] = float('nan')
    elif mutation == 'lineage_pair_count':
        changed['protected_obstruction_reuse']['C']['pair_count'] = 0
    elif mutation == 'lineage_kernel':
        changed['protected_obstruction_reuse']['C']['kernel'][0] = [60, 0]
    elif mutation == 'hold_extra_kernel':
        changed['protected_obstruction_reuse']['hold_control']['kernel'] = [[0, 0], [0, 0]]
    elif mutation == 'resource_missing':
        del changed['resource_envelope']
    elif mutation == 'resource_zero':
        changed['resource_envelope']['maximum_trial_source_pulses']['C'] = 0
    elif mutation == 'source_extra_run':
        changed['source_execution']['unexecuted'] = {}
    elif mutation == 'target_boolean_port':
        changed['target']['analyzer_ports'][1] = True
    else:
        changed['trials'][1]['records'][0]['coin'] = False
    with pytest.raises(ValueError):
        verify(changed, expected)


def test_density_admission_uses_exact_global_frobenius_bound(executed, expected):
    changed = deepcopy(executed[0])
    values = ('0x1.fb1fb2061bdf6p-7', '0x1.38138ee4deee9p-12')
    for column, value in zip((0, 6), values, strict=True):
        changed['preparation']['densities']['reference']['real'][0][column] = float.fromhex(value)
    q = Fraction(expected['finite_selector']['accepted_probability'])
    target = {(row, column): q * Fraction(value)
              for row, column, value in expected['preparation']['rho']['entries']}
    density = changed['preparation']['densities']['reference']
    exact_error = sum((Fraction(density['real'][row][column]) - target.get((row, column), 0)) ** 2
                      + Fraction(density['imag'][row][column]) ** 2
                      for row in range(66) for column in range(66))
    assert exact_error > Fraction(str(experiment.TOLERANCE)) ** 2
    with pytest.raises(ValueError, match='reference full density'):
        verify(changed, expected)


def test_density_inside_the_declared_bound_is_admitted(executed, expected):
    changed = deepcopy(executed[0])
    changed['preparation']['densities']['reference']['real'][0][0] += 1e-10
    assert verify(changed, expected)['all_pass']


def test_retained_receipt_is_bound_to_current_source_and_observations(executed, expected):
    retained = strict_json.load(experiment.ROOT / 'model/native_phase_experiment/results.json')
    assert verify(retained, expected)['all_pass']
    compare_receipt(retained, executed[0])


@pytest.mark.parametrize('control', ('loss', 'hold'))
def test_lineage_sanity_is_not_only_same_provider_comparison(executed, expected, monkeypatch, control):
    changed = deepcopy(executed[0])
    lineage = changed['protected_obstruction_reuse']
    if control == 'loss':
        lineage['C']['classification']['almostSure'] = False
    else:
        lineage['hold_control']['closedTrap'] = False
    monkeypatch.setattr(experiment, 'classify_phase_loss', lambda root: deepcopy(lineage))
    with pytest.raises(ValueError):
        verify(changed, expected)


@pytest.mark.parametrize('mutation', ('recipe_mass', 'balanced_recipe_density', 'bit_record_mass'))
def test_joint_profile_rejects_accumulated_scalar_errors(executed, expected, mutation):
    changed = deepcopy(executed[0])
    if mutation == 'recipe_mass':
        for component in changed['preparation']['components']:
            component['weight'] += 1.5e-10
    elif mutation == 'balanced_recipe_density':
        points = [item for item in changed['preparation']['components'] if item['kind'] == 'point']
        for index, component in enumerate(points):
            component['weight'] += 1.5e-10 if index < 33 else -1.5e-10
    else:
        changed['selector']['bit_record_probabilities'] = [
            value + 1.5e-10 for value in changed['selector']['bit_record_probabilities']]
    with pytest.raises(ValueError):
        verify(changed, expected)


def test_joint_profile_admits_errors_inside_unchanged_global_bounds(executed, expected):
    changed = deepcopy(executed[0])
    for component in changed['preparation']['components']:
        component['weight'] += 1e-13
    assert verify(changed, expected)['all_pass']
