from fractions import Fraction as F
import numpy as np
import pytest
from scipy.integrate import quad
from model.openwave_normalizable_particle_records import source as s


def test_normalizable_packet_and_complete_record_bounds():
    r=s.run();p=r['normalizable_packet']
    assert p['exact_L2_norm_squared']=='1'and p['position_variance_cells_squared_exact']=='655360'
    assert not p['finite_ring_or_delta_momentum_substitution']
    for row in r['members']:
        assert F(row['uniform_band_projector_probability_lower'])>F('0.9999614')
        assert F(row['position_record']['band_conditioned_inside_probability_lower'])>F('0.984025')
        assert F(row['complete_position_readout_mean_work_absolute_bound_source'])>0
    assert not r['complete_protocol']['finite_local_compilation_of_exact_band_filter_proved']
    for row in r['uniform_native_weak_packet_records'].values():
        assert F(row['conditional_band_input_current_and_proton_band_probability_interval'][0])>0


def test_independent_packet_integral_and_lattice_Fourier_normalization():
    kap=float(s.KAPPA);A=np.sqrt(15/(16*kap))
    density=lambda k:A*A*(1-(k/kap)**2)**2
    assert abs(quad(density,-kap,kap,epsabs=1e-12)[0]-1)<1e-12
    cosine=quad(lambda k:np.cos(k)*density(k),-kap,kap,epsabs=1e-12)[0]
    lo,hi=s.cosine_mean();assert abs(cosine-float(lo))<2e-15 and abs(cosine-float(hi))<2e-15
    # Exact analytical Fourier coefficients, independently summed on a large
    # nonperiodic position interval; the retained outside probability is bounded.
    xs=np.arange(1,20001,dtype=float);y=kap*xs
    prob=15/(2*np.pi*kap**5*xs**6)*(np.sin(y)-y*np.cos(y))**2
    total=5*kap/(6*np.pi)+2*prob.sum()
    assert 0<1-total<float(F(5)/(2*s.KAPPA**2*20000**2))


def test_complete_band_filter_and_gentle_position_bound_independent_matrix():
    psi=np.array([np.sqrt(.9),np.sqrt(.1),0.]);P=np.diag([1.,0.,1.]);p=psi@P@psi
    phi=P@psi/np.sqrt(p);rho=np.outer(psi,psi);out=np.outer(phi,phi)
    distance=np.linalg.norm(out-rho,ord='nuc')/2
    assert abs(distance-np.sqrt(1-p))<1e-14
    assert np.max(abs(P@P+(np.eye(3)-P)@(np.eye(3)-P)-np.eye(3)))<1e-15
    E=np.diag([0.,1.,0.]);assert abs(np.trace(E@(out-rho)))<=distance
    with pytest.raises(ValueError):s.moment(3)
    with pytest.raises(ValueError):s.moment(-2)
