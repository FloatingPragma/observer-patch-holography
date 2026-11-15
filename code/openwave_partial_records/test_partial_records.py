"""Partial OpenWave results on native records reproduce from their bound evidence."""
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location('partial_records_' + name, HERE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_verify_recomputes_every_row():
    result = load('verify').verify()
    assert result['archived_files_verified'] == 24 and result['providers'] == 11
    rec = result['recomputed']
    assert rec['clock']['lambda'] == '209/210' and rec['clock']['triangle'] == [0, 1, 5]
    assert rec['clock']['position_envelope'] == list(range(12)) and rec['clock']['within_one_twelve_port_cell']
    assert rec['spin']['crossing_sign'] == 1 and rec['spin']['spin_commutes_with_every_branch']
    assert rec['contact']['contact_branches_removed'] == 2520 and rec['contact']['pointer_outcomes'] == 17
    assert rec['residence']['exponent'] == 1.035144 and rec['residence']['profile95'] == [0.983376, 1.086373]
    assert rec['orbital']['shell_multiplicity'] == {'1': 1, '2': 4, '3': 9, '4': 16}
    assert result['sealed']['row_8']['mean_contact_energy'] == 11.051928582384567


def test_changed_archive_is_rejected(tmp_path):
    verify = load('verify')
    copy = tmp_path / 'evidence'
    shutil.copytree(verify.EVIDENCE, copy)
    receipt = copy / 'compact_internal_clock' / 'RECEIPT.json'
    receipt.write_text(receipt.read_text().replace('"209/210"', '"208/210"'))
    verify.EVIDENCE = copy
    with pytest.raises(ValueError, match='archived bytes differ'):
        verify.custody()


def test_unbound_file_is_rejected(tmp_path):
    verify = load('verify')
    copy = tmp_path / 'evidence'
    shutil.copytree(verify.EVIDENCE, copy)
    (copy / 'row20_residence' / 'extra.json').write_text('{}')
    verify.EVIDENCE = copy
    with pytest.raises(ValueError, match='exactly the bound files'):
        verify.custody()


def test_changed_graph_is_rejected(tmp_path):
    recompute = load('recompute')
    raw = json.loads(recompute.GRAPH.read_text())
    raw['seams'] = raw['seams'][1:]
    bad = tmp_path / 'graph.json'
    bad.write_text(json.dumps(raw))
    recompute.GRAPH = bad
    with pytest.raises(ValueError, match='sealed graph'):
        recompute.load_graph()


def test_claim_ledger_is_current():
    out = subprocess.run([sys.executable, str(ROOT / 'tools/build_openwave_partial_records_claim.py'), '--check'],
                         capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stdout + out.stderr
    ledger = json.loads((ROOT / 'claims/openwave_partial_records.json').read_text())
    assert [r['row'] for r in ledger['rows']] == [3, 4, 6, 8, 20, 31]
    assert all(r['status'] == 'partial' and r['official_status'] == 'untested' for r in ledger['rows'])
    assert ledger['official_OpenWave_awards_added'] == 0
