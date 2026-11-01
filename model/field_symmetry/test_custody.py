"""Fresh public evidence cannot authenticate a previously imported producer."""
import os
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from model.field_symmetry import packet


def test_preimport_python_change_uses_current_code_even_with_valid_stale_pyc(tmp_path):
    checkout = tmp_path/'checkout'
    for name in packet.source_bindings():
        target = checkout/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((packet.ROOT/name).read_bytes())
    environment = dict(os.environ)
    environment.pop('PYTHONPYCACHEPREFIX', None)
    environment.pop('PYTHONDONTWRITEBYTECODE', None)
    environment['PYTHONUTF8'] = '1'
    environment['OPENBLAS_NUM_THREADS'] = '1'
    environment['OMP_NUM_THREADS'] = '1'
    # Preserve size and timestamp deliberately: a fresh interpreter using the
    # default bytecode cache would also execute the old producer here.
    script = '''
import os
from pathlib import Path
from model.field_symmetry import readout
path = Path(readout.__file__)
assert Path(readout.__cached__).is_file()
stat = path.stat()
old = path.read_bytes()
new = old.replace(b'10**11', b'10**9 ')
assert old != new and len(old) == len(new)
path.write_bytes(new)
os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
from model.field_symmetry import packet
result = packet.report()
assert all(row['numerical_primitive_ceiling'] == '1/1000000000'
           for row in result['complete_source_experiments']), 'stale Python producer accepted'
assert result['source_bindings'][path.relative_to(Path.cwd()).as_posix()] == packet.source_bindings()['model/field_symmetry/readout.py']
'''
    result = subprocess.run([sys.executable, '-c', script], cwd=checkout,
                            env=environment, capture_output=True, text=True,
                            encoding='utf-8', timeout=180)
    assert result.returncode == 0, result.stdout+result.stderr


@pytest.fixture(scope='module')
def worker_payload():
    return packet.report()


@pytest.mark.parametrize('corruption', ['exit', 'malformed', 'extra_stdout', 'duplicate',
                                     'nonfinite', 'missing', 'extra_key', 'binding', 'after',
                                     'empty_experiments', 'empty_symmetry', 'empty_signal',
                                     'promote_physics', 'wrong_campaign', 'bool_count'])
def test_fresh_worker_refuses_broken_transport_and_custody(worker_payload, monkeypatch, corruption):
    result = deepcopy(worker_payload)
    if corruption == 'nonfinite':
        result['complete_source_experiments'][0]['complete_output_error'] = float('nan')
    elif corruption == 'missing':
        del result['symmetry']
    elif corruption == 'extra_key':
        result['unexpected'] = True
    elif corruption == 'binding':
        result['source_bindings'][next(iter(result['source_bindings']))] = '0'*64
    elif corruption == 'empty_experiments': result['complete_source_experiments'] = []
    elif corruption == 'empty_symmetry': result['symmetry'] = {}
    elif corruption == 'empty_signal': result['signal'] = {}
    elif corruption == 'promote_physics': result['physical_species_identified'] = True
    elif corruption == 'wrong_campaign': result['campaign_id'] = 'invented'
    elif corruption == 'bool_count': result['physical_comparison_count_added'] = False
    output = json.dumps(result)
    if corruption == 'malformed': output = '{'
    elif corruption == 'extra_stdout': output = 'unowned output\n'+output
    elif corruption == 'duplicate': output = output[:-1]+', "schema": "ambiguous"}'
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **k:
                        subprocess.CompletedProcess(a[0], 1 if corruption == 'exit' else 0,
                                                    stdout=output, stderr='worker failure'))
    monkeypatch.setattr(packet, '_produce', lambda: pytest.fail('no in-process fallback'))
    if corruption == 'after':
        before = dict(worker_payload['source_bindings'])
        after = dict(before); after[next(iter(after))] = '0'*64
        calls = iter((before, after))
        monkeypatch.setattr(packet, 'source_bindings', lambda: next(calls))
    with pytest.raises(ValueError):
        packet.report()


def test_public_replay_ignores_warm_parent_mutation_and_callers_cwd(monkeypatch, tmp_path):
    monkeypatch.setattr(packet.symmetry, 'certificate', lambda: pytest.fail('warm parent producer used'))
    monkeypatch.chdir(tmp_path)
    result = packet.report()
    assert result['symmetry']['full_carrier_checks'][-1]['complete_Fock_dimension'] == 65536


def test_retained_primitive_gate_uses_fresh_child_ceiling(worker_payload, monkeypatch):
    fresh = deepcopy(worker_payload)
    fresh['complete_source_experiments'][0]['numerical_primitive_ceiling'] = '1/1000000000000'
    retained = deepcopy(fresh)
    retained['complete_source_experiments'][0]['raw_primitive_residual'] = 5e-12
    # The older imported producer would accept this residual. It is not the
    # authority for the current child's stricter numerical contract.
    monkeypatch.setattr(packet.readout, 'PULSE_CEILING', packet.F(1, 10**11))
    with pytest.raises(ValueError, match='fresh numerical admission'):
        packet.compare(retained, fresh)
