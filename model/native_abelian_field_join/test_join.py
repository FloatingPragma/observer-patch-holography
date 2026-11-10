import pytest
import sympy as s
from fractions import Fraction as F
from . import model as m

@pytest.fixture(scope='module')
def result():return m.calculate()

def test_actual_native_reference_and_complete_dual_cost(result):
    local=result['local_reference'];globalref=result['global_reference']
    assert local['complete_reference_dimension']==2**31+1
    assert local['actual_ec_one_reference_gauge_factorization_checks']==12
    assert local['proper_cube_covariance_controls']==24
    assert len(local['complete_Fock_sector_exchange_controls'])==30
    assert globalref['full_compensator_modes']==52
    assert globalref['charge_exchange_preserves_energy_on_all_number_sectors']

def test_whole_bin_unnormalized_weight_and_all_six_records(result):
    row=result['whole_bin_coupling'];Q,p,w,a=s.symbols('Q p w a')
    probabilities=[1-Q,Q*(1-p),Q*p/2,Q*p*(s.Rational(1,2)-w),Q*p*w*a,Q*p*w*(1-a)]
    assert s.expand(sum(probabilities))==1
    assert row['conditional_native_bin_probability']=='2w_J'
    assert row['literal_native_CAR_pump_success_probability']=='1/2'
    assert not row['extra_factor_one_half_or_one_quarter_inserted']
    assert F(row['uniform_unconditional_native_arrival_lower'])>F('9.8e-6')
    for setting in row['sharper_actual_field_amplitude_transfers']:
        assert F(setting['native_unconditional_arrival_lower'])>0
        assert F(setting['native_unconditional_arrival_and_nonvacuum_lower'])>0
    assert not result['intervals']['unpenalized_full_charged_ground_identified']
    assert F(result['intervals']['alternative_three_active_charge_margin'])>0

def test_complete_kinematic_gauss_and_boundary_zeros(result):
    row=result['full_kinematic_charge_guards']
    assert row['full_kinematic_nonzero_hop_Gauss_checks']>1000
    assert row['full_kinematic_charge_pump_Gauss_checks']>1000
    assert row['retained_boundary_zero_hops']>0 and row['retained_boundary_zero_pumps']>0
    assert not row['physical_input_assumed_by_diagonal_charge_checks']

def test_literal_nonstationary_bin_dynamics_and_inactive_complement():
    import numpy as np
    from scipy.linalg import expm
    # A separate finite control retains two distinct native-bin energies and
    # a lower out-of-bin level. Only the bin vector transports; actual h does
    # not preserve that vector. No spectral result is transferred from this.
    center=.69;epsilon=.001;eta=1/24;t=.25;weight=.42
    psi=np.array([np.sqrt(1/3),np.sqrt(2/3),0.])
    project=np.outer(psi,psi);Pbin=np.diag([1.,1.,0.])
    native=np.diag([center-epsilon,center+epsilon,.1])
    averaged=np.diag([center,center,.1])
    assert np.linalg.norm(native@psi-center*psi)>0
    a=np.diag(np.sqrt(np.arange(1,4)),1)
    W=expm(.3*(a.T-a));Hph=np.diag(1.4*np.arange(4))
    def build(local):
        H=np.zeros((36,36))
        for y in range(3):
            sl=slice(12*y,12*(y+1))
            H[sl,sl]=np.kron(local+(4+.1*y)*np.eye(3),np.eye(4))+np.kron(np.eye(3),Hph)
        for u,v in ((0,1),(1,2)):
            su=slice(12*u,12*(u+1));sv=slice(12*v,12*(v+1))
            H[sv,su]=-eta*np.kron(project,W);H[su,sv]=-eta*np.kron(project,W.T)
        return H
    H=build(native);H0=build(averaged)
    assert np.linalg.norm(H-H0,2)<=epsilon+1e-14
    inside=np.zeros(36);inside[:12]=np.kron(psi,[1.,0.,0.,0.])
    outside=np.zeros(36);outside[8]=1
    assert np.allclose(H@outside,(4+.1)*outside)
    pump=np.sqrt(weight)*inside+np.sqrt(.5-weight)*outside
    assert abs(np.dot(pump,pump)-.5)<1e-15
    actual=expm(-1j*t*H);comparison=expm(-1j*t*H0)
    assert np.linalg.norm(actual-comparison,2)<=epsilon*t+1e-13
    target=slice(12,24)
    amplitude=np.linalg.norm((comparison@inside)[target])
    actual_probability=np.linalg.norm((actual@pump)[target])**2
    assert abs(actual_probability-weight*np.linalg.norm((actual@inside)[target])**2)<1e-16
    assert actual_probability>=weight*max(0,amplitude-epsilon*t)**2-1e-15
    whole_bin=np.kron(np.eye(3),np.kron(Pbin,np.eye(4)))
    assert np.linalg.norm(H@whole_bin-whole_bin@H)<1e-14
    assert abs(np.linalg.norm((np.eye(36)-whole_bin)@actual@pump)**2-(.5-weight))<1e-15
    assert np.linalg.eigvalsh(H)[0]<center+4-.1


def test_public_scientific_comparison_and_construction_scope():
    import copy
    from .certificate import raw_certificate,compare_research,certificate
    value=raw_certificate()
    assert compare_research(value)['all_scientific_fields_equal']
    wrong=copy.deepcopy(value)
    wrong['whole_bin_coupling']['literal_native_CAR_pump_success_probability']='1/4'
    with pytest.raises(ValueError):compare_research(wrong)
    public=certificate();scope=public['scope_boundaries']
    assert scope['each_original_intrablock_H_preserved']
    assert not scope['original_uniform_native_lattice_H_preserved']
    assert not scope['unpenalized_full_charged_ground_identified']
    assert scope['penalty_alternative_is_separate_Hamiltonian']
    assert public['whole_bin_coupling']['P0_projects_all26_copies_including_partners']
