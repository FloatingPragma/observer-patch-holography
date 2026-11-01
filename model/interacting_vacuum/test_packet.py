import sys

import pytest

from model.interacting_vacuum import packet


def test_supplied_missing_receipt_path_fails_before_expensive_replay(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['interacting_vacuum', '--check', str(tmp_path/'missing.json')])
    monkeypatch.setattr(packet, 'report', lambda: pytest.fail('missing receipt must fail before replay'))
    with pytest.raises(FileNotFoundError): packet.main()


def test_receipt_binds_core_proofs_source_and_no_float_admission():
    bindings = packet.source_bindings()
    for path in ('model/interacting_vacuum/ground.py',
                 'model/interacting_vacuum/global_ground.py',
                 'model/interacting_vacuum/charged.py',
                 'model/cube_eta_protection/cube.py',
                 'model/cube_eta_protection/repair.py',
                 'paper/tex_fragments/INTERACTING_VACUUM.tex'):
        assert path in bindings
    assert packet.canonical({'a': True}) != packet.canonical({'a': 1})
