from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from . import algebra,bounds,packet


def test_exact_shell_interpolation_and_stress_frames():
    result=algebra.certificate()
    assert result['shell_moment_A']=='1/420'
    assert result['shell_moment_B']=='1/630'
    assert result['tracefree_stress_frame_determinant']=='13824/15625'
    assert result['boosted_ball_response_determinant']=='24407490807/2000000000'


def test_missing_or_repeated_direction_loses_tomographic_rank():
    assert algebra.response_matrix(algebra.DIRECTIONS[:-1]).rank()==8
    repeated=algebra.DIRECTIONS[:-1]+algebra.DIRECTIONS[:1]
    assert algebra.response_matrix(repeated).det()==0


def test_wave_and_hole_sign_are_both_needed():
    n=s.Matrix([s.Rational(3,5),s.Rational(4,5),0])
    p=s.Matrix([-1,*n]);j=-p;eta=s.diag(1,-1,-1,-1)
    classical=(p*j.T+j*p.T)/2
    assert classical[0,0]==-1
    steered=-eta*classical*eta
    assert steered[0,0]==1 and s.trace(eta*steered)==0


def test_full_car_wick_steering_and_local_versus_global_midpoint():
    # Two modes, one filled orbital; retained environment is not traced early.
    c=[]
    for mode in range(2):
        op=np.zeros((4,4))
        for col in range(4):
            if col&(1<<mode):op[col^(1<<mode),col]=(-1)**((col&((1<<mode)-1)).bit_count())
        c.append(op)
    psi=np.array([0.,1.,1.,0.])/np.sqrt(2)
    rho=np.outer(psi,psi);number=[x.T@x for x in c]
    B=number[1]-.5*np.eye(4)
    assert np.linalg.norm(B@c[0]-c[0]@B)==0
    derivative=np.trace(rho@B@number[0])
    covariance=np.array([[psi@a.T@b@psi for b in c] for a in c])
    assert derivative==pytest.approx(-abs(covariance[0,1])**2)
    states=[]
    for sign in (-1,1):
        root=np.diag(np.sqrt(np.diag(np.eye(4)+sign*B)))
        vector=root@psi;states.append(np.outer(vector,vector))
    midpoint=(states[0]+states[1])/2
    assert np.trace(midpoint@number[0])==pytest.approx(.5)
    assert np.trace(midpoint@c[0])==pytest.approx(0)
    assert np.linalg.norm(midpoint-rho)>.01
    for theta in (-.5,-.1,.1,.5):
        p=.5-theta/4
        relative=p*np.log(2*p)+(1-p)*np.log(2*(1-p))
        assert relative<=float(bounds.entropy_remainder(F(str(abs(theta)))))


def test_uniform_response_invertibility_including_entropy_difference():
    budget=bounds.neumann_budget(bounds.MAX_RADIUS_RATIO,finite_entropy_difference=True)
    assert 0<budget<F(1,4)
    assert bounds.neumann_budget(2*bounds.MAX_RADIUS_RATIO)>F(1,4)
    matrix=np.array(algebra.response_matrix()).astype(float)
    inverse=np.linalg.norm(np.linalg.inv(matrix),ord=np.inf)
    assert inverse==pytest.approx(68/5)
    # Adverse rank-one entry error at the exact uniform envelope.
    entry=float((bounds.BALL_ENTRY_ERROR+bounds.ENTROPY_ENTRY_ERROR)*bounds.MAX_RADIUS_RATIO**2)
    disturbed=matrix+entry*np.ones((9,9))
    assert np.linalg.norm(np.linalg.inv(disturbed),ord=np.inf)<272/15


@pytest.mark.parametrize('args',[(0,),(-1,),(F(1,4),),(True,),(0.001,)])
def test_invalid_radius_rejected(args):
    with pytest.raises(ValueError):bounds.neumann_budget(*args)


@pytest.mark.parametrize('args',[(True,),(-1,),(F(3,4),F(1,2)),(F(1,2),1)])
def test_invalid_entropy_domain_rejected(args):
    with pytest.raises(ValueError):bounds.entropy_remainder(*args)


@pytest.fixture(scope='module')
def receipt():return packet.build()


@pytest.mark.parametrize('key',['global_pure_vacuum_affine_two_sided_state_family_claimed',
    'native_shell_operations_implemented','joint_lattice_collar_limit_supplied',
    'area_Newton_or_generalized_stationarity_supplied','geometry_stress_null_balance_supplied',
    'complete_Einstein_register_inhabited','ten_component_trace_direction_supplied'])
def test_no_unproved_geometry_or_global_state_promotion(receipt,key):
    changed=deepcopy(receipt);changed['theorem'][key]=True
    with pytest.raises(ArithmeticError):packet.validate(changed)


def test_nested_types_and_current_custody_are_strict(receipt):
    for section,key in [('theorem','complete_Einstein_register_inhabited'),
                         ('bounds','global_pure_vacuum_affine_two_sided_perturbation_asserted'),
                         ('prior_custody','prior_full_source_replay')]:
        changed=deepcopy(receipt);changed[section][key]=0
        with pytest.raises(ArithmeticError):packet.validate(changed)
    changed=deepcopy(receipt);changed['algebra']['stress_covariant_sign']='negative'
    with pytest.raises(ArithmeticError):packet.validate(changed)
