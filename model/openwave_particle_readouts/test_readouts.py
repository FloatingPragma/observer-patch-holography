"""Independent finite-matrix controls for the prepared-line bounds."""
import numpy as np
import json
import pytest
import sympy as s
from scipy.linalg import expm

from .readouts import prepared_line_bounds, native_rotation_controls, disturbance_bound


def test_complete_current_receipt_and_tamper_controls():
    from .readouts import HERE, certificate, validate_pins

    saved = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    validate_pins(saved)
    actual = certificate()
    assert actual == saved
    assert actual['prepared_line_motion']['source_clock_settings'] == ['0', '1', '3', '6']
    assert actual['same_source_hypercharge_background']['full_source_controls']['actual_hopping_coefficient'] == '1/17179869184'
    # Fresh numerical replay and source binding are separate checks.
    import copy
    corrupted = copy.deepcopy(saved)
    name = next(iter(corrupted['input_sha256']))
    corrupted['input_sha256'][name] = '0'*64
    with pytest.raises(ValueError, match='source pin mismatch'):
        validate_pins(corrupted)


@pytest.mark.parametrize('seed', [7, 13, 29, 41])
def test_independent_projector_and_time_evolution(seed):
    rng=np.random.default_rng(seed)
    a=np.array([0., .1, .3, .5, .8])
    h0=np.diag(a)
    raw=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
    perturb=(raw+raw.conj().T)/2
    perturb*=0.0001/np.linalg.norm(perturb,2)
    h=h0+perturb
    pulse=np.array([.1,.7,.3,.2,.1],dtype=complex)
    pulse/=np.linalg.norm(pulse)
    values,vectors=np.linalg.eigh(h)
    selected=(values>.08)&(values<.12)
    chi=vectors[:,selected]@(vectors[:,selected].conj().T@pulse)
    p=np.linalg.norm(chi)**2
    psi=chi/np.sqrt(p)
    weight=float(abs(pulse[1])**2)
    # Covers selected leakage and all unwanted atomic components separately.
    distances=[.08,.02,.18,.38,.68]
    eps=.0001*sum(1/x for x in distances)
    for t in (0,1,100,1000):
        result=prepared_line_bounds((s.Rational(str(weight)),)*2,
                                   s.Rational(str(eps)),s.Rational('0.0001'),t)
        evolved=expm(-1j*t*h)@psi
        nonorigin=1-abs(evolved[1])**2
        survival=abs(np.vdot(psi,evolved))**2
        assert nonorigin<=float(result['outside_origin_including_screened_cloud_upper'])+1e-12
        assert survival>=float(result['initial_state_return_probability_lower'])-1e-12
        assert float(result['probability'][0])<=p<=float(result['probability'][1])
    # Removing the actual mixing error would claim perfect localization.
    broken=prepared_line_bounds((s.Rational(str(weight)),)*2,0,0,0)
    assert float(broken['outside_origin_including_screened_cloud_upper'])==0
    assert 1-abs(psi[1])**2>1e-10


def test_arbitrary_vector_in_same_bin_is_not_admitted():
    # A degenerate energy bin includes a remote/cloud state even though the
    # specified local preparation has no amplitude there.
    h=np.diag([1.,1.,2.]);origin=np.diag([1.,0.,0.])
    cloud=np.array([0.,1.,0.]);prepared=np.array([1.,0.,0.])
    assert np.array_equal(h@cloud,cloud)
    assert np.vdot(cloud,origin@cloud)==0
    assert np.vdot(prepared,origin@prepared)==1
    with pytest.raises(ValueError):
        prepared_line_bounds((s.Rational(1,4),)*2,s.Rational(1,2),0,0)


def test_probability_records_retain_all_failures_and_cloud():
    q,p_failure,p_outside,p_line,origin=s.symbols('q f o p w')
    complete=(1-q)+q*p_failure+q*p_outside+q*p_line*origin+q*p_line*(1-origin)
    assert s.expand(complete.subs(p_outside,1-p_failure-p_line))==1
    broken=complete-q*p_line*(1-origin)
    assert s.expand(broken.subs(p_outside,1-p_failure-p_line))!=1


def test_full_source_spin_sign_and_reference_control():
    value=native_rotation_controls()
    assert value['all_occupation_count_controls']==121**2
    assert value['two_pi_pointer_probabilities_on_line']==['0','1']
    assert value['four_pi_pointer_probabilities_on_line']==['1','0']
    assert value['joint_source_reference_rotation_pointer_probabilities']==['1','0']
    assert value['broken_joint_rotation_rejected']


def test_bounds_reject_invalid_inputs():
    for weight,error,v,t in [((0,0),0,0,0),((1,2),0,0,0),((1,1),0,-1,0),((1,1),0,0,-1)]:
        with pytest.raises(ValueError):prepared_line_bounds(weight,error,v,t)


def test_disturbance_control_covers_leakage_and_loses_power_under_large_drive():
    probability=s.Rational(1,10000)
    protected=disturbance_bound(probability,s.Rational(1,1000),s.Rational(1,10**9),10**6)
    assert protected>s.Rational(999,1000)
    assert disturbance_bound(probability,s.Rational(1,1000),s.Rational(1,100),10**6)==0
    with pytest.raises(ValueError):disturbance_bound(2,0,0,0)
