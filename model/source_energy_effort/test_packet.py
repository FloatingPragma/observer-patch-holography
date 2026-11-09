import copy
import json
import math
import subprocess
import sys
from fractions import Fraction as F

import numpy as np
import pytest

from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_complete_native_energy_and_independent_clock(result):
    measurement = result['native_energy_measurement']
    assert len(measurement['settings']) == 64
    assert measurement['creation']['pulses'] == 1092
    assert measurement['full_work_rows_per_energy_instrument'] == 262144
    assert F(measurement['complete_protocol_TV_error_upper']) < F(665, 10**15)
    assert measurement['complete_probability_residual'] < 2e-11
    assert len([row for row in measurement['settings'] if row['coefficient'] == '0']) == 20
    for row in measurement['settings']:
        assert len(row['charge_energy_success_probabilities']) == 2
        assert all(len(q) == 2 for q in row['charge_energy_success_probabilities'])
        assert row['input_prefix_success_was_not_renormalized'] is True
    phase = result['native_clock_action']
    assert phase['pulses'] == 4907 and phase['controlled_scalar_phase_retained'] is True
    assert F(phase['energy_derivative_error_upper']) < F(359, 10**6)
    for index, row in enumerate(measurement['charge_energy_readouts']):
        assert abs(phase['unnormalized_mean_energy_from_sine_derivative'][index]-row['unconditioned_excitation_moment']) < float(F(phase['energy_derivative_error_upper']))
    for row in packet.instruments(result):
        assert len(row['native_first_failure_weights']) == row['pulses']
        assert row['postselected_or_renormalized'] is False


def test_retained_exact_receipt_and_readable_energy_table(result):
    packet.compare(json.loads((packet.HERE/'results.json').read_text()), result)
    assert (packet.HERE/'ENERGY_READOUT.md').read_text() == result['energy_table_markdown']


@pytest.mark.parametrize('key', [k for k, v in packet.claims().items() if v is False])
def test_physical_promotions_are_rejected(result, monkeypatch, key):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    changed = copy.deepcopy(result)
    changed['claims'][key] = True
    with pytest.raises(ValueError): packet.validate(changed)


@pytest.mark.parametrize('kind', ['setting', 'zero', 'failure', 'phase', 'energy', 'hash', 'typed_bool'])
def test_complete_record_and_source_mutations_fail(result, monkeypatch, kind):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    changed = copy.deepcopy(result)
    if kind == 'setting': changed['native_energy_measurement']['settings'].pop()
    if kind == 'zero': changed['native_energy_measurement']['settings'][-1]['charge_energy_success_probabilities'][0].pop()
    if kind == 'failure': changed['native_clock_action']['native_first_failure_weights'].pop()
    if kind == 'phase': changed['native_clock_action']['controlled_scalar_phase_retained'] = False
    if kind == 'energy': changed['exact_same_source_work']['channels'][0]['mean_excitation_lower'] = '2'
    if kind == 'hash': changed['source_bindings'][next(iter(changed['source_bindings']))] = '0'*64
    if kind == 'typed_bool': changed['all_pass'] = 1
    with pytest.raises(ValueError): packet.validate(changed)


def test_global_failure_amplitudes_cannot_be_forged(result, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    changed = copy.deepcopy(result)
    row = changed['native_clock_action']
    row['native_first_failure_weights'] = (np.asarray(row['native_first_failure_weights'])+1e-20).tolist()
    with pytest.raises(ValueError, match='amplitude'): packet.validate(changed)


def test_source_import_closure_in_clean_full_replay():
    code = '''
import pathlib,sys
from model.source_energy_effort import packet
packet.produce()
bindings=packet.source_bindings()
missing=[]
for module in list(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name is None: continue
    path=pathlib.Path(name).resolve()
    if path.suffix=='.py' and path.is_relative_to(packet.ROOT):
        key=path.relative_to(packet.ROOT).as_posix()
        if key not in bindings: missing.append(key)
assert not missing, sorted(set(missing))
'''
    run = subprocess.run([sys.executable, '-W', 'error', '-c', code], cwd=packet.ROOT,
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout+run.stderr


def test_check_reads_exact_requested_path(tmp_path):
    path = tmp_path/'missing.json'
    run = subprocess.run([sys.executable, '-m', 'model.source_energy_effort', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0 and str(path) in run.stderr


def test_conditioned_clock_comparison_propagates_original_probability_error():
    retained = json.loads((packet.HERE/'results.json').read_text())
    changed = copy.deepcopy(retained)
    clock = changed['native_clock_action']
    time = float(F(clock['time']))
    row = clock['joint_charge_clock_success_probabilities'][0]
    row[0] += 9e-12
    row[1] -= 9e-12
    clock['unnormalized_mean_energy_from_sine_derivative'][0] = (row[0]-row[1])/time
    difference = abs(clock['unnormalized_mean_energy_from_sine_derivative'][0]
                     - retained['native_clock_action']['unnormalized_mean_energy_from_sine_derivative'][0])
    assert difference > 2e-11
    assert difference < 2*2e-11/time
    packet.compare(changed, retained)
    packet.compare(retained, changed)


@pytest.mark.parametrize('mutation', ['score', 'probability', 'clock', 'failure', 'physical_bound', 'negative_score'])
def test_conditioned_clock_comparison_rejects_inconsistent_or_unadmitted_records(mutation):
    retained = json.loads((packet.HERE/'results.json').read_text())
    changed = copy.deepcopy(retained)
    clock = changed['native_clock_action']
    if mutation == 'score':
        # Smaller than the propagated allowance, but not computed from its records.
        clock['unnormalized_mean_energy_from_sine_derivative'][0] += 1e-9
    elif mutation == 'probability':
        row = clock['joint_charge_clock_success_probabilities'][0]
        row[0] += 3e-11
        clock['unnormalized_mean_energy_from_sine_derivative'][0] = (row[0]-row[1])/float(F(clock['time']))
    elif mutation == 'clock':
        clock['time'] = '1/2048'
        clock['unnormalized_mean_energy_from_sine_derivative'] = [
            (row[0]-row[1])/float(F(clock['time']))
            for row in clock['joint_charge_clock_success_probabilities']]
    elif mutation == 'failure':
        clock['native_first_failure_weights'].pop()
    elif mutation == 'negative_score':
        clock['joint_charge_clock_success_probabilities'][0].reverse()
        clock['unnormalized_mean_energy_from_sine_derivative'][0] *= -1
        # Even a self-consistent negative score retains the old type/domain guard.
        retained = copy.deepcopy(changed)
    else:
        clock['energy_derivative_error_upper'] = '1/1000'
    with pytest.raises(ValueError):
        packet.compare(changed, retained)


def test_conditioned_clock_comparison_rejects_one_ulp_score_forgery():
    retained = json.loads((packet.HERE/'results.json').read_text())
    changed = copy.deepcopy(retained)
    scores = changed['native_clock_action']['unnormalized_mean_energy_from_sine_derivative']
    scores[0] = math.nextafter(scores[0], math.inf)
    with pytest.raises(ValueError, match='own original probabilities'):
        packet.compare(changed, retained)
