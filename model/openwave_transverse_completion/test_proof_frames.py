"""Fixed complete proof coordinates preserve every certified rest interval."""
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import pytest
from . import source as s
from .proof_frames import fixed_eigen,WITNESSES


def vacuum_matrix(cap):
    n=np.arange(-cap,cap+1,dtype=np.int64)
    return np.diag(4*n*n+4)+np.diag(np.full(2*cap,-2,dtype=np.int64),1)+np.diag(np.full(2*cap,-2,dtype=np.int64),-1)


def test_all_retained_rest_certificates_replay_without_a_new_eigenframe(monkeypatch):
    def no_frame(*args,**kwargs):raise AssertionError('A retained proof must not request a fresh eigenframe')
    monkeypatch.setattr(s,'eigh',no_frame)
    here=Path(__file__).parent
    transverse=json.loads((here/'results.json').read_text())
    assert s.energy_certificate()==transverse['infinite_loop_charged_rest_certificate']
    dirac=json.loads((here.parent/'openwave_dirac_maxwell_refinement/results.json').read_text())
    for record in dirac['complete_executions']:
        assert s.field_floor(record['cap'])==record['field_vacuum_certificate']


def test_complete_proof_frame_rejects_malformed_or_wrong_operator_inputs():
    record=json.loads(WITNESSES.read_text());A=vacuum_matrix(24)
    changes=[lambda w:w['integer_frame'][0].__setitem__(0,True),
             lambda w:w['integer_frame'][0].__setitem__(0,2**80),
             lambda w:w['integer_frame'].pop(),
             lambda w:w['integer_eigenvalue_centres'].reverse()]
    for mutate in changes:
        bad=deepcopy(record);mutate(bad['frames']['vacuum_24'])
        with pytest.raises(AssertionError):fixed_eigen(A,'vacuum_24',bad)
    wrong=A.copy();wrong[0,0]+=1
    with pytest.raises(AssertionError):fixed_eigen(wrong,'vacuum_24',record)
    bad=deepcopy(record)
    for row in bad['frames']['vacuum_24']['integer_frame']:row[0]=0
    eigen=fixed_eigen(A,'vacuum_24',bad)
    with pytest.raises(AssertionError):s.exact_cap_enclosure(A,eigen)


def test_complete_proof_frame_inventory_is_pinned():
    from .__main__ import pins
    assert 'model/openwave_transverse_completion/complete_frames.json' in pins()
