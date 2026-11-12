from fractions import Fraction as F
import numpy as np
import pytest
from scipy.linalg import expm,sqrtm
from model.finite_area_response.response import Interval
from model.openwave_energy_lapse.source import native,exp_interval
from .observations import clocks,wave_matrix,wave,metric_and_curvature,preparation

def test_all_native_clocks_share_actual_complete_unitary_ratio():
    rows=clocks();m=np.asarray(native()[0],float);e,u=np.linalg.eigh(m);K=m-e[0]*np.eye(8)
    for row in rows['rows']:
        k=row['level'];v=(u[:,0]+u[:,k])/np.sqrt(2);plus=(u[:,0]+1j*u[:,k])/np.sqrt(2)
        for q,expected in zip((F(rows['prepared_q_near']),F(rows['prepared_q_far'])),row['plus_Y_probability_intervals']):
            final=expm(-1j*np.exp(float(q))*K)@v;p=abs(np.vdot(plus,final))**2
            lo,hi=map(F,expected);assert abs(p-float((lo+hi)/2))<3e-15
        joint=row['complete64_joint_records'];assert len(joint)==64
        lower=sum(F(r['probability_interval'][0])for r in joint);upper=sum(F(r['probability_interval'][1])for r in joint)
        assert lower<=1<=upper and upper-lower<F(1,10**18)
    ratios={tuple(r['near_to_far_frequency_ratio_interval'])for r in rows['rows']};assert len(ratios)==1

def test_uniform_wave_scaling_and_independent_edge_action():
    q=np.full(7,-.0002);matrix=wave_matrix(q)
    L=np.diag(np.full(7,2.))+np.diag(np.full(6,-1.),1)+np.diag(np.full(6,-1.),-1)
    np.testing.assert_allclose(matrix,np.exp(4*q[0])*L,atol=1e-15,rtol=0)
    q=np.array([-.0001,-.0002,0,0,-.0003,0,0]);S=np.diag(np.exp(2*q))
    np.testing.assert_allclose(wave_matrix(q),S@L@S,atol=1e-15,rtol=0)

def test_complete_wave_records_and_measurement_work():
    q=np.array([-.0001,-.0002,0,0,-.0003,0,0]);r=wave(q,3);h=sqrtm(wave_matrix(q));x=np.arange(7)
    psi=np.exp(-(x-7)**2/8+1j*.7*x);psi/=np.linalg.norm(psi);v=expm(-3j*h)@psi
    expected=[.5,*list(.5*abs(v)**2)]
    np.testing.assert_allclose(r['complete_vacuum_and_position_probabilities'],expected,atol=2e-14,rtol=0)
    actual=.5*np.dot(abs(v)**2,np.diag(h))-.5*np.vdot(v,h@v).real
    assert abs(actual-r['mean_detector_work'])<1e-14
    assert abs(r['source_energy_before']-r['source_energy_after'])<1e-14

def test_linear_geometry_balance_and_gamma_countercontrol():
    r=metric_and_curvature();assert r['gamma_zero_countercontrol_G11_at_3_0_0']=='2/27'
    assert r['linearized_spatial_Einstein_tensor']=='zero'
    assert not r['regional_entropy_stationarity_premises_discharged']
    for row in r['finite_residual_rows']:
        lo,hi=map(F,row['G00_minus_2rho_over_kappa_interval']);assert max(abs(lo),abs(hi))<F(1,10**10)

def test_nonpositive_domain_and_flat_control():
    with pytest.raises(ValueError):wave_matrix([0,.1,0])
    with pytest.raises(ValueError):wave_matrix([0,0])
    _,q=preparation();curved=wave(q,12);flat=wave([F(0)]*31,12)
    assert np.linalg.norm(np.array(curved['complete_vacuum_and_position_probabilities'])-flat['complete_vacuum_and_position_probabilities'])>1e-7
