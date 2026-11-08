import copy
import json
import subprocess
import sys
from fractions import Fraction as F

import numpy as np
import pytest

from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_actual_blank_preparation_and_full_record_domains(result):
    pulse = result['native_blank_preparation_and_current_pulse']
    assert pulse['actual_input_is_blank'] is True
    assert pulse['ground_preparation_pulses'] == 1065
    assert pulse['pulses'] == 3441
    assert pulse['full_work_rows'] == 65536
    assert F(pulse['decoded_regional_entropy_change_upper']) < F(-4, 10**6)
    assert pulse['failure_and_offcode_decoder_state'] == 'I_4/4'
    assert pulse['decoded_entropy_decrease_diagnostic'] > 4e-6
    assert pulse['decoded_minimum_eigenvalue_diagnostic'] > 1/100
    assert pulse['decoded_trace_residual'] < 2e-11
    for row in packet.instrument_rows(result):
        assert len(row['native_first_failure_weights']) == row['pulses']
        assert row['postselected_or_renormalized'] is False
        assert row['complete_output_distance_diagnostic'] < 2e-11
    detector = result['native_complete_detector']
    assert len(detector['settings']) == 16
    assert all(row['retained_exterior_columns'] == 64 for row in detector['settings'])
    identity = detector['settings'][0]
    assert identity['pauli'] == 'II'
    assert identity['ideal_pointer_probabilities'][1] == 0
    assert result['native_modular_flow']['constant_modular_phase_retained'] is True
    contrast = result['native_complete_current_readout']
    assert F(contrast['plus_probability_increase_lower']) > F(17, 10**6)
    assert contrast['plus_probability_increase_diagnostic'] > float(F(contrast['plus_probability_increase_lower']))
    for row in contrast['experiments']:
        assert row['input_is_blank'] is True
        assert row['full_work_rows'] == 131072
        assert sum(row['pointer_success_probabilities'])+row['complete_first_failure_probability'] == pytest.approx(1, abs=2e-11)


def test_retained_receipt(result):
    packet.compare(json.loads((packet.HERE/'results.json').read_text()), result)


@pytest.mark.parametrize('key', [k for k, v in packet.claims().items() if not v])
def test_false_physical_premises_cannot_be_promoted(result, monkeypatch, key):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    value['claims'][key] = True
    with pytest.raises(ValueError): packet.validate(value)


@pytest.mark.parametrize('change', ['zero_outcome', 'failure_record', 'floor', 'clock', 'hash', 'typed_bool'])
def test_mutations_reject_incomplete_or_changed_source_claim(result, monkeypatch, change):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    if change == 'zero_outcome': value['native_complete_detector']['settings'].pop(0)
    if change == 'failure_record': value['native_blank_preparation_and_current_pulse']['native_first_failure_weights'].pop()
    if change == 'floor': value['source_reference']['actual_floor'] = '1/4'
    if change == 'clock': value['claims']['same_complete_eight_spin_source_reference_and_H_B_energy_convention'] = False
    if change == 'hash': value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    if change == 'typed_bool': value['all_pass'] = 1
    with pytest.raises(ValueError): packet.validate(value)


def test_collectively_tiny_forged_failures_fail_global_amplitude_check(result, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    row = value['native_blank_preparation_and_current_pulse']
    row['native_first_failure_weights'] = (np.asarray(row['native_first_failure_weights'])+1e-20).tolist()
    with pytest.raises(ValueError, match='amplitude'): packet.validate(value)


def test_actual_runtime_imports_bound_in_clean_process():
    script = '''
import json,pathlib,sys
from model.spin_modular_geometry import packet
packet.produce()
bound=packet.source_bindings()
missing=[]
for module in list(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None: continue
    path=pathlib.Path(path).resolve()
    if path.suffix=='.py' and path.is_relative_to(packet.ROOT):
        name=path.relative_to(packet.ROOT).as_posix()
        if name not in bound: missing.append(name)
assert not missing, sorted(set(missing))
'''
    run = subprocess.run([sys.executable, '-W', 'error', '-c', script], cwd=packet.ROOT,
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout+run.stderr


def test_explicit_receipt_path_is_consumed_before_replay(tmp_path):
    path = tmp_path/'missing.json'
    run = subprocess.run([sys.executable, '-m', 'model.spin_modular_geometry', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0 and str(path) in run.stderr
