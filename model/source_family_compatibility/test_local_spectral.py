"""Source, completeness, probability and locality controls for the local clock."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import local_spectral as m


@pytest.fixture(scope='module')
def certificate():
    return m.certificate()


@pytest.fixture(scope='module')
def source_clock():
    h,rho,a,projections=m.dimer_source()
    h=np.array(h,dtype=float);rho=np.array(rho,dtype=float);a=np.array(a,dtype=float)
    energies,vectors=np.linalg.eigh(h)
    ground=vectors[:,0];n=32
    unitaries=[(vectors*np.exp(-1j*t*energies))@vectors.T for t in range(n)]
    labels=[];effects=np.zeros((16,16),complex);kraus=[]
    for label,op in (('addition',a.T),('removal',a)):
        histories=[u@op@u.conj().T for u in unitaries]
        for j in range(n):
            k=sum(np.exp(2j*np.pi*j*t/n)*v for t,v in enumerate(histories))/n
            effects+=k.conj().T@k;kraus.append(k);labels.append((label,j))
    return h,rho,a,energies,vectors,ground,labels,kraus,effects


def test_actual_source_complete_every_input_and_both_charges(source_clock):
    *_,labels,kraus,effects=source_clock
    assert len(kraus)==len(labels)==64
    assert np.linalg.norm(effects-np.eye(16),ord=2)<2e-14
    for col in range(16):
        assert abs(sum(np.linalg.norm(k[:,col])**2 for k in kraus)-1)<2e-14
    # One charge or the wrong Fourier normalization destroys completeness.
    plus=sum(k.conj().T@k for k in kraus[:32])
    assert np.linalg.norm(plus-np.eye(16),ord=2)>.49
    assert np.linalg.norm(32*effects-np.eye(16),ord=2)>30


def test_actual_ground_records_match_certified_original_weights(source_clock,certificate):
    h,rho,a,energies,vectors,ground,labels,kraus,_=source_clock
    rows=certificate['actual_dimer_complete_clock']['all_phase_records']
    all_weights=[]
    for (label,j),k in zip(labels,kraus,strict=True):
        probability=np.linalg.norm(k@ground)**2
        lo,hi=map(F,rows[j]['each_charge_original_probability'])
        assert float(lo)-2e-14<=probability<=float(hi)+2e-14
        all_weights.append(probability)
        op=a.T if label=='addition'else a
        e=energies-energies[0]
        scalar=sum(np.exp(2j*np.pi*j*t/32)*np.exp(-1j*t*e)for t in range(32))/32
        spectral=(vectors*scalar)@vectors.T@op@ground
        assert np.linalg.norm(k@ground-spectral)<2e-14
    assert abs(sum(all_weights[:32])-.5)<1e-14
    assert abs(sum(all_weights[32:])-.5)<1e-14
    assert abs(sum(all_weights)-1)<2e-14


def test_excited_input_not_given_false_ground_relative_spectrum(source_clock):
    h,rho,a,energies,vectors,ground,labels,kraus,_=source_clock
    vacuum=np.eye(16)[:,0];differences=[]
    for (label,j),k in zip(labels,kraus,strict=True):
        op=a.T if label=='addition'else a
        e=energies-energies[0]
        scalar=sum(np.exp(2j*np.pi*j*t/32)*np.exp(-1j*t*e)for t in range(32))/32
        wrong=(vectors*scalar)@vectors.T@op@vacuum
        differences.append(abs(np.linalg.norm(k@vacuum)**2-np.linalg.norm(wrong)**2))
    assert sum(differences)>.5


def test_phase_count_does_not_multiply_the_complete_boundary_error():
    a=m.local_filter_bound(9,F(1,8),64,0,1,F(1,8))
    b=m.local_filter_bound(17,F(1,16),64,0,1,F(1,8))
    assert a['source_duration']==b['source_duration']==1
    assert a['complete_output_halftrace_upper']==b['complete_output_halftrace_upper']
    assert a['boundary_operator_error']==b['boundary_operator_error']
    assert not a['error_has_phase_count_prefactor']


def test_complete_dilation_parseval_controls_joint_reference():
    # Genuine source clocks at two durations test the Fourier equality itself.
    h,rho,a,energies,vectors,ground,_,_,_=source_clock.__wrapped__()
    n=8;dt=F(1,8)
    def unitary(t):return (vectors*np.exp(-1j*float(t)*energies))@vectors.T
    differences=[];time_differences=[]
    for op in (a.T,a):
        histories=[]
        for t in range(n):
            u=unitary(t*dt);v=unitary(t*dt/2)
            histories.append(u@op@u.conj().T-v@op@v.conj().T)
        time_differences.extend(histories)
        differences.extend(sum(np.exp(2j*np.pi*j*t/n)*x for t,x in enumerate(histories))/n for j in range(n))
    lhs=sum(x.conj().T@x for x in differences)
    rhs=sum(x.conj().T@x for x in time_differences)/n
    assert np.linalg.norm(lhs-rhs)<1e-14
    # Tensoring an untouched two-level reference preserves the operator identity.
    assert np.linalg.norm(np.kron(lhs-rhs,np.eye(2)))<2e-14


def test_default_positive_event_is_original_and_volume_independent(certificate):
    x=certificate['growing_default_source_local_positive_event']
    assert x['source_U_over_c']=='8' and x['phase_records']==[0,1025]
    assert x['Fejer_leakage_upper']=='11/202'
    assert F(x['published_each_cut_group_probability_lower'])==F(25,200704)
    assert F(x['each_cut_group_probability_lower'])>F(25,200704)
    p=x['local_cut_plan']
    assert p['phase_count']==4096 and p['charge_phase_records']==8192
    assert p['cut_radius']==3120 and p['cut_sites']==6241**3
    assert F(p['complete_output_halftrace_upper'])<=F(25,401408)
    assert not x['cut_ground_substituted'] and not x['shrinking_unit_hopping_energy_floor']
    assert not p['ground_marginal_prepared_by_planner']


def test_nonwrapped_window_preserves_mass_and_adds_complete_error(certificate):
    rows=certificate['actual_dimer_complete_clock']['all_phase_records']
    low=sum(F(rows[j]['each_charge_original_probability'][0])for j in (2,3))
    x=m.window_from_records(32,1,2,3,F(1,2),3,low,F(1,2),0,F(1,1000))
    expected=max(F(0),low-F(1,1000)-m.leakage(32,1,F(1,2))/2)
    assert x['original_true_window_probability_lower']==expected>0
    assert x['true_positive_energy_window'][0]==0
    assert x['true_positive_energy_window'][1]>1
    assert x['original_true_window_probability_lower']<F(1,2)
    # A missing UV tail would certify an aliased response as low energy.
    y=m.window_from_records(32,1,2,3,F(1,2),3,low,F(1,2),F(1,2))
    assert y['original_true_window_probability_lower']==0
    assert not x['sharp_projection_implemented']


def test_phase_alias_control_requires_actual_UV_input():
    n=32;e=.5
    def probability(energy,j):
        return abs(sum(np.exp(1j*t*(2*np.pi*j/n-energy))for t in range(n))/n)**2
    assert max(abs(probability(e,j)-probability(e+2*np.pi,j))for j in range(n))<1e-14
    with pytest.raises(ValueError):m.window_from_records(32,1,0,4,F(1,4),4,F(1,2),F(1,2),0)


def test_source_clock_homogeneity_preserves_error_and_weights():
    # Ordinary source duration is reduced by c under K(c,lambda*c)=c*K(1,lambda).
    a=m.local_filter_bound(32,1,100,0,1,F(1,8))
    b=m.local_filter_bound(32,F(1,8),100,0,1,1)
    assert a['complete_output_halftrace_upper']==b['complete_output_halftrace_upper']
    h,_,_=m.exact_source(2,((0,1,1),),F(1,8),F(1))
    k,_,_=m.exact_source(2,((0,1,1),),F(1),F(8))
    assert k==8*h


def test_cache_result_has_no_mutable_reference_alias(certificate):
    expected=repr(certificate)
    value=m.certificate();value['actual_dimer_complete_clock']['all_phase_records'][0]['each_charge_original_probability'][0]='0'
    value['growing_default_source_local_positive_event']['source_U']='999'
    assert repr(m.certificate())==expected


@pytest.mark.parametrize('args',[(True,1,0,1,F(1,8),F(1,10)),(2,0,0,1,F(1,8),F(1,10)),
    (2,1,0,2,F(1,8),F(1,10)),(2,1,0,1,0,F(1,10)),(2,1,0,1,F(1,8),0)])
def test_planner_rejects_invalid_source_inputs(args):
    with pytest.raises(ValueError):m.plan_filter(*args)


def test_physical_and_state_boundaries(certificate):
    assert certificate['physical_mass_prediction']is None
    assert certificate['specified_ground_marginal_error_required']
    assert not certificate['infinite_ground_preparation_or_subsequence_modulus_derived']
    assert not certificate['shrinking_unit_hopping_spectral_weight_proved']
    assert not certificate['spatial_Walsh_label_identified_as_momentum']


def test_original_lipschitz_score_bound_and_period_guard():
    x=m.spectral_score_bound(4096,F(1,16),F(1,4),24,F(1,2),F(1,1000),1,1,F(1,1000))
    assert x['original_spectral_score_error_upper']==F(1,8)+F(11,202)+F(1,250)
    assert not x['conditional_renormalization_applied']
    with pytest.raises(ValueError):m.spectral_score_bound(32,1,1,3,F(1,2),0)
    # A positive energy near zero can produce negative signed clock centres;
    # clamping them preserves rather than discards their original probability.
    energy=.02;n=64
    probs=[abs(sum(np.exp(1j*t*(2*np.pi*j/n-energy))for t in range(n))/n)**2 for j in range(n)]
    assert sum(probs[n//2+1:])>0
    score=sum(min(1,max(0,2*np.pi*(j if j<=n//2 else j-n)/n))*p for j,p in enumerate(probs))
    bound=m.spectral_score_bound(n,1,F(1,2),1,1,0)['original_spectral_score_error_upper']
    assert abs(score-energy)<=float(bound)


def test_fourth_taper_complete_time_normalization():
    for m0 in (1,2,4,16):
        coefficients,square=m.fourth_taper(m0)
        expected=np.array([1],dtype=object)
        for _ in range(4):expected=np.convolve(expected,np.ones(m0,dtype=object))
        assert coefficients[:len(expected)]==tuple(expected)
        assert all(x==0 for x in coefficients[len(expected):])
        assert sum(F(x*x,square)for x in coefficients)==1
        # Fourier Parseval for every pair of time positions makes the full
        # source CAR clock complete, without postselection on a time record.
        n=len(coefficients)
        for i in range(n):
            for j in range(n):
                scalar=sum(np.exp(2j*np.pi*k*(i-j)/n)for k in range(n))/n
                assert abs(scalar-(i==j))<2e-14


def test_fourth_taper_inverse_window_keeps_original_UV_and_locality(certificate):
    x=certificate['default_complete_inverse_energy_window'];p=F(25,100352)
    assert x['phase_count']==8192 and x['charge_phase_records']==16384
    assert x['phase_records']==[0,1665]
    assert F(x['original_total_UV_tail_above_20'])<p/16
    assert F(x['taper_leakage_upper'])==F(4096,7)*F(11,224)**7<p/16
    assert F(x['each_original_true_window_probability_lower'])>3*p/4
    assert F(x['published_each_original_true_window_probability_lower'])==3*p/4
    assert F(x['each_local_measured_group_probability_lower'])>3*p/4
    assert not x['native_tapered_clock_executed'] and x['physical_mass_assigned']is None


def test_taper_kernel_tail_directly_on_finite_circles():
    for m0 in (4,16,64):
        coeff,square=m.fourth_taper(m0);n=len(coeff)
        for phase in (0.,.2,.9):
            angle=2*np.pi*np.arange(n)/n-phase
            amp=np.polynomial.polynomial.polyval(np.exp(1j*angle),coeff)
            probs=abs(amp)**2/(n*square)
            assert abs(sum(probs)-1)<3e-13
            dist=abs((angle+np.pi)%(2*np.pi)-np.pi)
            tail=m.fourth_taper_leakage(n,1,1)
            assert sum(probs[dist>=1])<=float(tail)+3e-13


def test_tapered_instrument_on_every_actual_dimer_input(source_clock,certificate):
    h,rho,a,energies,vectors,ground,_,_,_=source_clock
    coefficients,square=m.fourth_taper(8);n=len(coefficients)
    unitaries=[(vectors*np.exp(-1j*t*energies))@vectors.T for t in range(n)]
    total=np.zeros((16,16),complex);weights=[]
    records=certificate['actual_dimer_tapered_complete_clock']['all_phase_records']
    for op in (a.T,a):
        histories=[u@op@u.conj().T for u in unitaries]
        for j in range(n):
            k=sum(coef*np.exp(2j*np.pi*j*t/n)*history
                  for t,(coef,history)in enumerate(zip(coefficients,histories,strict=True)))/np.sqrt(n*square)
            total+=k.conj().T@k
            p=float(np.linalg.norm(k@ground)**2);weights.append(p)
            lo,hi=map(F,records[j]['each_charge_original_probability'])
            assert float(lo)-3e-14<=p<=float(hi)+3e-14
            e=energies-energies[0]
            scalar=sum(coef*np.exp(2j*np.pi*j*t/n)*np.exp(-1j*t*e)
                       for t,coef in enumerate(coefficients))/np.sqrt(n*square)
            assert np.linalg.norm(k@ground-(vectors*scalar)@vectors.T@op@ground)<2e-14
    assert np.linalg.norm(total-np.eye(16),ord=2)<3e-14
    assert abs(sum(weights[:n])-.5)<1e-14 and abs(sum(weights[n:])-.5)<1e-14
