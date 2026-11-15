"""Signed-source charge, rest-energy and frame results reproduce exactly."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location('signed_source_' + name, HERE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_verify_recomputes_three_results():
    result = load('verify').verify()
    charge = result['recomputed']['charge']
    assert charge['face_smith'] == [1] * 19 + [0] and charge['boundary_rank'] == 11
    assert charge['native_trials_absorbed'] == 264 and charge['changed_law_controls_absorbed'] == 0
    energy = result['recomputed']['rest_energy']
    assert abs(energy['E6'] - 6.014331228348673) < 1e-12
    assert abs(energy['pair_keV'] - 1024.4508487023415) < 1e-9
    frame = result['recomputed']['frame']
    assert frame['branches'] == 720 and frame['half_turn_returns_covariance'] and frame['half_turn_reverses_mean']


def test_changed_graph_is_rejected(tmp_path):
    recompute = load('recompute')
    raw = json.loads(recompute.GRAPH.read_text())
    raw['seams'] = raw['seams'][1:]
    bad = tmp_path / 'graph.json'
    bad.write_text(json.dumps(raw))
    recompute.GRAPH = bad
    try:
        recompute.load_graph()
    except ValueError as error:
        assert 'sealed graph' in str(error)
    else:
        raise AssertionError('changed graph accepted')


def test_claim_ledger_is_current():
    out = subprocess.run([sys.executable, str(ROOT / 'tools/build_openwave_signed_source_claim.py'), '--check'],
                         capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stdout + out.stderr
