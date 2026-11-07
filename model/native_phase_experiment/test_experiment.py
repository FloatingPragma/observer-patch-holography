"""Independent source and receipt regression tests for the phase experiment."""
from copy import deepcopy
from fractions import Fraction
import subprocess
import sys

import numpy as np
import pytest

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
