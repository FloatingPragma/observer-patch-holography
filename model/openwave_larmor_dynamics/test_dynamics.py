import numpy as np
from fractions import Fraction
from .dynamics import stationary,quasistatic,midpoint_step,observe,diagonalize,evolve
from .experiment import PARAMS,contraction_certificate,current_controls,frozen_control
from .certificate import equivalent


def test_global_field_contraction_bound():
    c=contraction_certificate()
    assert Fraction(c['fixed_point_contraction_upper'])<Fraction(61,100)
    assert Fraction(c['absolute_kernel_Hessian_upper'])>0


def test_signed_current_and_complete_moment_instruments():
    c=current_controls();assert c['positive_branch_only_current_error']>1e-4


def test_independent_field_start_and_time_reversal():
    state,A,_=stationary(PARAMS)
    rng=np.random.default_rng(998);changed=A+.1*rng.normal(size=A.shape)
    af=quasistatic(state,A,PARAMS,.05)[0];bf=quasistatic(state,changed,PARAMS,.05)[0]
    assert np.linalg.norm(af-bf)<2e-11
    advanced,mid,_,_=midpoint_step(state,af,PARAMS,.05,5.)
    returned,_,_,_=midpoint_step(advanced,mid,PARAMS,.05,-5.)
    assert np.linalg.norm(returned-state)<2e-10
    obs,_=observe(advanced,mid,PARAMS,.05,True)
    assert obs['continuity_residual']<1e-11 and abs(obs['energy_record_total_probability']-1)<1e-11
    assert all(abs(sum(p)-1)<1e-14 for p in obs['spin_measurement_probabilities'])


def test_frozen_linear_propagator_and_source_pins():
    r=frozen_control(100.);assert r['direct_matrix_exponential_error']<2e-12
    assert max(x['self_field_defect'] for x in r['rows'])>1e-6
    assert not equivalent({'input_sha256':{'a':'1'}},{'input_sha256':{'a':'2'}})
