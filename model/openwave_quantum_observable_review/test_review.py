import copy,json
import pytest
from .review import ROOT,covariance_review,trial_review,packet_guard,cube
import numpy as np


def test_covariance_offdiagonal_omission_rejected():
    r=json.loads((ROOT/'model/openwave_capacity_stress/results.json').read_text())
    assert len(covariance_review(r))==3
    bad=copy.deepcopy(r);bad['comparisons'][0]['full_covariance_rational'][0][1]='0'
    with pytest.raises(AssertionError):covariance_review(bad)


def test_omitted_failure_and_work_rejected():
    r=json.loads((ROOT/'model/openwave_transverse_gravity/results.json').read_text())['full_source_executions'][0]
    trial_review(r)
    bad=copy.deepcopy(r);bad['retained_no_arrival_and_photon_vacuum']['total_probability']=0
    with pytest.raises(AssertionError):trial_review(bad)
    bad=copy.deepcopy(r);bad['complete_history_work']['terminal_and_native_readout']=0
    with pytest.raises(AssertionError):trial_review(bad)
    bad=copy.deepcopy(r);bad['complete_preparation_and_work']['mean_preparation_work']=0
    with pytest.raises(AssertionError):trial_review(bad)


def test_longitudinal_packet_and_non_normalized_packet_rejected():
    G,C=cube(2);v=G[:,0].astype(complex);v/=np.linalg.norm(v)
    with pytest.raises(AssertionError):packet_guard(G,v)
    with pytest.raises(AssertionError):packet_guard(G,np.zeros(len(v)))
    assert not np.any(C@G)
