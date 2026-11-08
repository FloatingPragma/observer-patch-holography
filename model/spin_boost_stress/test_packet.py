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


def test_actual_complete_weighted_flow(result):
    row = result['native_source_candidate']
    assert row['input_is_blank'] is True and row['complete_work_rows'] == 131072
    assert row['pulses'] == 11821
    assert len(row['native_first_failure_weights']) == 11821
    assert F(row['native_unconditioned_minus_pointer_probability_lower']) > F(5, 10**9)
    assert row['native_unconditioned_minus_pointer_probability'] > float(F(row['native_unconditioned_minus_pointer_probability_lower']))
    total = row['native_unconditioned_minus_pointer_probability']+row['native_plus_pointer_probability']+row['complete_first_failure_probability']
    assert total == pytest.approx(1, abs=2e-11)
    assert row['reserved_six_level_pointer_probabilities'] == [0.0]*4
    assert row['postselected_or_renormalized'] is False


def test_retained_receipt(result):
    packet.compare(json.loads((packet.HERE/'results.json').read_text()), result)


@pytest.mark.parametrize('key', [k for k, v in packet.claims().items() if not v])
def test_false_physical_attachments_cannot_be_promoted(result, monkeypatch, key):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    value['claims'][key] = True
    with pytest.raises(ValueError): packet.validate(value)


@pytest.mark.parametrize('kind', ['force', 'dipole', 'failure', 'speed', 'hash', 'typed_bool'])
def test_exact_source_or_record_mutation_fails(result, monkeypatch, kind):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    if kind == 'force': value['source_forced_stress'][0]['total_force_nonzero_as_source_operator'] = False
    if kind == 'dipole': value['coarse_energy_and_dipole_transport']['maps_with_nonzero_subcell_energy_dipole'] = 0
    if kind == 'failure': value['native_source_candidate']['native_first_failure_weights'].pop()
    if kind == 'speed': value['candidate_speed_controls'][0][0]['physical_speed_selected_by_this_normalization'] = True
    if kind == 'hash': value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    if kind == 'typed_bool': value['all_pass'] = 1
    with pytest.raises(ValueError): packet.validate(value)


def test_collective_failure_amplitude_mutation(result, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    value = copy.deepcopy(result)
    row = value['native_source_candidate']
    row['native_first_failure_weights'] = (np.asarray(row['native_first_failure_weights'])+1e-20).tolist()
    with pytest.raises(ValueError, match='amplitude'): packet.validate(value)


def test_source_import_coverage_on_clean_full_producer():
    script = '''
import pathlib,sys
from model.spin_boost_stress import packet
packet.produce()
bindings=packet.source_bindings(); missing=[]
for module in list(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name is None: continue
    path=pathlib.Path(name).resolve()
    if path.suffix=='.py' and path.is_relative_to(packet.ROOT):
        key=path.relative_to(packet.ROOT).as_posix()
        if key not in bindings: missing.append(key)
assert not missing, sorted(set(missing))
'''
    run = subprocess.run([sys.executable, '-W', 'error', '-c', script], cwd=packet.ROOT,
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout+run.stderr


def test_explicit_check_path(tmp_path):
    path = tmp_path/'missing.json'
    run = subprocess.run([sys.executable, '-m', 'model.spin_boost_stress', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0 and str(path) in run.stderr
