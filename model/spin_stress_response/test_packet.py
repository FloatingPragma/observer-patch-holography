import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_actual_native_complete_response(result):
    assert result['all_pass'] is True
    for row in result['complete_native_runs']:
        assert row['pulses'] == 144
        assert row['native_first_failure_slots'] == 510
        assert len(row['native_first_failure_weights']) == 144
        assert all(len(x) == 4 for x in row['native_first_failure_weights'])
        assert row['reserved_six_level_pointer_probabilities'] == [0.0]*4
        assert sum(row['source_pointer_probabilities'])+row['complete_first_failure_probability'] == pytest.approx(1, abs=2e-11)
        assert row['success_postselected_or_renormalized'] is False
        assert row['complete_output_distance_diagnostic'] < 1e-12
    p = [r['source_pointer_probabilities'][0] for r in result['complete_native_runs']]
    assert p[0] > p[1] > p[2]


def test_retained_receipt(result):
    retained = json.loads((packet.HERE/'results.json').read_text())
    packet.compare(retained, result)


@pytest.mark.parametrize('key', [k for k, v in packet.claims().items() if v is False])
def test_reject_claim_promotion(result, monkeypatch, key):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    altered = copy.deepcopy(result)
    altered['claims'][key] = True
    with pytest.raises(ValueError):
        packet.validate(altered)


@pytest.mark.parametrize('change', ['missing_failure', 'normalized_success', 'floor', 'stale_source', 'typed_boolean'])
def test_reject_corrupt_complete_or_exact_data(result, monkeypatch, change):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    altered = copy.deepcopy(result)
    if change == 'missing_failure':
        altered['complete_native_runs'][0]['native_first_failure_weights'].pop()
    elif change == 'normalized_success':
        altered['complete_native_runs'][0]['source_pointer_probabilities'][0] += 0.01
    elif change == 'floor':
        altered['same_source_modular_response']['uniform_regional_floor'] = '1/4'
    elif change == 'stale_source':
        altered['source_bindings'][next(iter(altered['source_bindings']))] = '0'*64
    else:
        altered['all_pass'] = 1
    with pytest.raises(ValueError):
        packet.validate(altered)


def test_joint_tiny_failure_forgery_is_rejected(result, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: result)
    altered = copy.deepcopy(result)
    values = np.asarray(altered['complete_native_runs'][0]['native_first_failure_weights'])
    altered['complete_native_runs'][0]['native_first_failure_weights'] = (values+1e-20).tolist()
    with pytest.raises(ValueError, match='amplitude'):
        packet.validate(altered)


def test_binding_covers_runtime_imports_in_fresh_process():
    script = '''
import json, pathlib, sys
from model.spin_stress_response import packet
packet.produce()
bound=set(packet.source_bindings())
missing=[]
for module in list(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None: continue
    path=pathlib.Path(path).resolve()
    if path.suffix == '.py' and path.is_relative_to(packet.ROOT):
        relative=path.relative_to(packet.ROOT).as_posix()
        if relative not in bound: missing.append(relative)
print(json.dumps(sorted(set(missing))))
assert not missing
'''
    completed = subprocess.run([sys.executable, '-c', script], cwd=packet.ROOT,
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout+completed.stderr


def test_check_consumes_explicit_path(tmp_path):
    missing = tmp_path/'missing.json'
    completed = subprocess.run([sys.executable, '-m', 'model.spin_stress_response',
                                '--check', str(missing)], cwd=packet.ROOT,
                               capture_output=True, text=True)
    assert completed.returncode != 0
    assert str(missing) in completed.stderr
