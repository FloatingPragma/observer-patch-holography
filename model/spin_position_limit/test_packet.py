"""Receipt path, complete output and claim-boundary checks."""
from copy import deepcopy
import json
import subprocess
import sys
import pytest
from . import packet


@pytest.fixture(scope='module')
def result():return packet.build()


def test_complete_retained_result_matches_producer(result):
    assert packet.canonical(result)==packet.OUTPUT.read_text()


def test_scope_promotion_rejected(tmp_path,result):
    changed=deepcopy(result)
    changed['scope']['mass_shell_or_empirical_species_identification']=True
    p=tmp_path/'promoted.json';p.write_text(packet.canonical(changed))
    with pytest.raises(ValueError):packet.verify(p)


def test_cli_reads_supplied_path(tmp_path):
    p=tmp_path/'invalid.json';p.write_text('{"schema": "other"}')
    run=subprocess.run([sys.executable,'-m','model.spin_position_limit','--check',str(p)],
                       cwd=packet.ROOT,text=True,capture_output=True)
    assert run.returncode!=0
    assert 'receipt' in run.stderr


def test_duplicate_keys_rejected_before_replay(tmp_path):
    p=tmp_path/'duplicates.json';p.write_text('{"schema":"x","schema":"x"}')
    with pytest.raises(ValueError):packet.verify(p)
