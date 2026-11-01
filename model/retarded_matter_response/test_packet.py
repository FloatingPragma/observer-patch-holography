import sys

import pytest

from model.retarded_matter_response import packet


def test_missing_supplied_receipt_fails_before_replay(tmp_path,monkeypatch):
    monkeypatch.setattr(sys,'argv',['retarded','--check',str(tmp_path/'missing.json')])
    monkeypatch.setattr(packet,'report',lambda:pytest.fail('must validate input path first'))
    with pytest.raises(FileNotFoundError):packet.main()


def test_prior_global_gap_and_complete_source_binding_are_authenticated():
    assert packet.prior_certificates()['full_global_ground_gap_authenticated']=='7/1024'
    bindings=packet.source_bindings()
    for path in ('model/interacting_vacuum/global_ground.py',
                 'model/interacting_vacuum/results.json',
                 'model/thermodynamic_ground/retarded.py',
                 'model/retarded_matter_response/static.py',
                 'paper/tex_fragments/RETARDED_MATTER_RESPONSE.tex',
                 'model/observer_particles/geometric_matter_join.py'):
        assert path in bindings
    assert packet.canonical({'x':True})!=packet.canonical({'x':1})
