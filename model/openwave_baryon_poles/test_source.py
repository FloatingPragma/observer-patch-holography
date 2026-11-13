from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
import pytest
import sympy as sy
from model.openwave_baryon_poles import source as s


def test_actual_trial_poles_and_complete_spectral_weights():
    r=s.run()
    assert len(r['momentum_fibre_records'])==12
    for row in r['momentum_fibre_records']:
        p=row['pole'];z=F(p['prepared_trial_pole_weight_lower'])
        assert z>F(99,100)
        assert z+F(p['prepared_trial_other_spectral_weight_upper'])==1
        assert F(p['original_input_pole_event_probability_lower'])==z/5
        assert len(row['all_five_low_source_GeV_intervals'])==5
        assert F(p['complete_ground_gap_lower_source'])>0
    assert r['mass_refit_performed']is False
    assert r['physical_3plus1_Lorentz_particle_pole_identification']is False
    assert all(F(x['pole']['prepared_trial_pole_weight_lower'])>F(9999997,10000000)for x in r['rest_source_poles'])


def test_residue_enclosure_against_independent_full_coupled_spectrum():
    A=sy.diag(3,7);G=sy.eye(2);v=sy.Matrix([10,1])
    _,mean,var=s.trial_moments(A,G,v)
    bound=s.residue_bounds([(F(3),F(3)),(F(7),F(7))],F(30),F(1,25),F(1,575),mean,var)
    B=np.diag([.1,.2]);D=np.diag([30.,35.]);H=np.block([[np.array(A,float),B],[B,D]])
    e,U=np.linalg.eigh(H);trial=np.array([10.,1.,0.,0.])/np.sqrt(101)
    z=abs(trial@U[:,0])**2
    assert z>=float(F(bound['prepared_trial_pole_weight_lower']))
    assert np.linalg.norm(U[:2,0])**2>=float(F(bound['five_field_pole_residue_trace_interval'][0]))
    weights=abs(trial@U)**2
    assert abs(weights.sum()-1)<1e-14
    for point in(-1.,0.,1.):
        direct=trial@np.linalg.solve(point*np.eye(4)-H,trial)
        pole=weights[0]/(point-e[0]);complement=sum(weights[j]/(point-e[j])for j in range(1,4))
        assert abs(direct-pole-complement)<1e-14


def test_complete_filter_preserves_nonselective_energy_but_not_branch_energy():
    H=np.array([[2.,.3,.4],[.3,5.,.2],[.4,.2,12.]])
    rho=np.diag([.5,.5,0.]);v=np.array([3.,4.,0.])/5;K=np.outer(v,v);J=np.eye(3)-K
    post=K@rho@K+J@rho@J
    assert np.max(abs(post-rho))<1e-15
    assert abs(np.trace(H@(post-rho)))<1e-14
    assert abs(np.trace(K@rho@K)-.5)<1e-15
    assert abs(v@H@v-np.trace(H@rho))>.1
    assert np.max(abs(K.T@K+J.T@J-np.eye(3)))<1e-15


def test_witness_and_gap_mutations_fail_closed():
    witness=json.loads((Path(s.__file__).with_name('trial_witnesses.json')).read_text())[s.key(False,F(1),F(0))]
    witness['real_integer_coefficients'][0]='01'
    with pytest.raises(AssertionError):s.vector(witness)
    with pytest.raises(AssertionError):s.residue_bounds([(F(3),F(4)),(F(4),F(5))],F(30),F(1),F(1),F(3),F(0))
    with pytest.raises(AssertionError):s.residue_bounds([(F(3),F(3)),(F(7),F(7))],F(30),F(1),F(1),F(3),F(1000))
