"""Independent finite controls for the quantum Higgs and regulator join."""
from fractions import Fraction as F
import json

import numpy as np
from scipy.linalg import expm
import pytest
import sympy as s

from . import interacting_chiral as m


def jordan_creator(n, i):
    # Independent Jordan-Wigner tensor order: mode0 is the rightmost factor.
    z = s.diag(1, -1); raising = s.Matrix([[0, 0], [1, 0]])
    factors = [s.eye(2) if j > i else raising if j == i else z for j in reversed(range(n))]
    return s.ImmutableSparseMatrix(s.kronecker_product(*factors))


@pytest.fixture(scope='module')
def ops():
    return m.quantum_higgs_operators()


def test_literal_CAR_signs():
    for i in range(6):
        assert m.creator(6, i) == jordan_creator(6, i)
    a = [jordan_creator(6, i) for i in range(6)]
    for i in range(6):
        for j in range(6):
            assert a[i]*a[j].H+a[j].H*a[i] == (s.eye(64) if i == j else s.zeros(64))


def test_full_quantum_scalar_norm_from_independent_operators(ops):
    a = [jordan_creator(6, i) for i in range(6)]
    pair = [a[0]*a[5]-a[1]*a[4], a[2]*a[5]-a[3]*a[4]]
    gram = sum((b.H*b for b in pair), s.zeros(64))
    # Exact annihilating polynomial with all roots in[0,4], plus an attained4.
    polynomial = s.eye(64)
    for root in range(5):
        polynomial = polynomial*(gram-root*s.eye(64))
    assert polynomial == s.zeros(64)
    vacuum = s.eye(64)[:, 0]
    assert gram*vacuum == 4*vacuum
    assert sum((b.H*b for b in ops['pairs']), s.zeros(64)) == gram
    assert ops['W']*ops['W'] == s.zeros(192)
    assert ops['V'] == ops['V'].H


def test_full_current_invariance_and_missing_Higgs_negative_control(ops):
    v = ops['V']; nf = ops['N_F']; nh = ops['N_H']
    assert nf*v-v*nf == 2*(nh*v-v*nh)
    assert nh*v-v*nh != s.zeros(192)
    for g in ops['generators']:
        assert m._zero(g*v-v*g)
    # The last actual generator is q6Y. Removing its scalar action breaks it.
    fermion_charge = s.diag(*(sum((-3 if j < 4 else 6) for j in range(6) if mask >> j & 1)
                             for mask in range(64)))
    broken = s.ImmutableSparseMatrix(s.kronecker_product(s.eye(3), s.I*fermion_charge))
    assert broken*v-v*broken != s.zeros(192)


def test_noncommuting_common_H0_and_complete_full_input_pulses(ops):
    h0 = ops['N_H']+ops['N_F']/7
    assert h0*ops['V']-ops['V']*h0 != s.zeros(192)
    weights = [F(55,72), F(11,72), F(5,72), F(1,72)]
    matrices = []
    # A one-step finite gauge-invariant instrument; all input columns retained.
    for active in (0,0,1,1):
        h = np.array(h0+active*ops['V'], complex)
        matrices.append(expm(-1j*h/32))
    complete = sum((float(w)*u.conj().T@u for w,u in zip(weights,matrices)), np.zeros((192,192),complex))
    assert np.linalg.norm(complete-np.eye(192),ord=2) < 2e-14
    q = np.array(ops['N_F']-2*ops['N_H'],complex)
    assert max(np.linalg.norm(q@u-u@q,ord=2) for u in matrices) < 2e-14


def test_complex_family_rank_one_is_two_independent_rotations():
    u = s.Matrix([s.Rational(3,5), s.I*s.Rational(4,5)])
    p = u*u.H
    assert p*p == p and p != p.T
    # The L creator has coefficientsu, the e creator has conjugateu.
    a = [jordan_creator(4,i) for i in range(4)]
    left = u[0]*a[0]+u[1]*a[1]
    right = s.conjugate(u[0])*a[2]+s.conjugate(u[1])*a[3]
    direct = sum((p[i,j]*a[i]*a[2+j] for i in range(2) for j in range(2)),s.zeros(16))
    wrong = sum((p[i,j]*a[i]*a[2+j] for i in range(2) for j in range(2)),s.zeros(16))
    assert direct == left*right
    assert wrong != left*(u[0]*a[2]+u[1]*a[3])


def test_feedback_exact_cost_and_no_H0_factor():
    out=m.quantum_feedback_budget(1000000,1,F(1,12))
    assert out['vertex_norm_bound'] == 2
    assert out['uniform_fixed_memory_bound'] == F(216,11000000)
    assert out['archive_distance_squared_upper'] == F(432,11000000)
    assert out['common_H0_norm_enters_bound'] is False
    larger=m.quantum_feedback_budget(1000000,1,F(1,12),gain=3,cells=2)
    assert larger['joint_nonreset_probability'] == 1-F(11,12)**2
    assert larger['vertex_norm_bound'] == 12
    assert larger['uniform_fixed_memory_bound'] > out['uniform_fixed_memory_bound']


@pytest.mark.parametrize('kwargs',[{'cells':0},{'cells':True},{'gain':-1},{'p':1},{'p':True},{'steps':0}])
def test_feedback_domain_rejections(kwargs):
    args={'steps':100,'time':1,'p':F(1,12)};args.update(kwargs)
    with pytest.raises(ValueError):m.quantum_feedback_budget(**args)


def test_link_mean_and_bidirectional_residual_exact_on_mixed_states():
    u=s.Matrix([[0,0,0],[1,0,0],[0,0,0]])
    vectors=[s.Matrix([1,1,0]),s.Matrix([1,s.I,1]),s.Matrix([2,1,2])]
    rho=sum((v*v.H/(v.H*v)[0]/3 for v in vectors),s.zeros(3))
    mean=s.trace(rho*u);q=s.trace(rho*(u.H*u+u*u.H))
    residual=s.trace(rho*((u-mean*s.eye(3)).H*(u-mean*s.eye(3))
                         +(u.H-s.conjugate(mean)*s.eye(3)).H*(u.H-s.conjugate(mean)*s.eye(3))))
    assert s.simplify(residual-(q-2*mean*s.conjugate(mean))) == 0
    assert s.simplify(residual-2*mean*s.conjugate(mean)) > 0
    assert mean != 0
    plus=s.Matrix([1,1,0])/s.sqrt(2)
    pure=plus*plus.H
    assert s.trace(pure*u)==s.Rational(1,2)
    assert s.trace(pure*(u.H*u+u*u.H))-2*s.Rational(1,2)**2 == 2*s.Rational(1,2)**2


def test_identity_link_state_error_has_no_dark_escape():
    u=s.zeros(60);u[37,0]=1
    gram=(s.eye(60)-u).H*(s.eye(60)-u)
    values=gram.eigenvals()
    assert values[s.Integer(1)]==58
    assert values[(3-s.sqrt(5))/2]==1
    assert values[(3+s.sqrt(5))/2]==1
    assert m.identity_link_control()['mean_lambda_relative_two_direction_RMS_lower']=='1'


def test_dressed_path_preserves_each_local_generator():
    result=m.dressed_path_control()
    assert result['source_current_checks']==36
    assert result['full_CAR_bilinear_identity'] is True
    assert result['single_charged_endpoint_is_gauge_invariant'] is False


def test_gapped_quantum_projector_is_not_a_CAR_functor():
    result=m.quantum_projector_car_control()
    assert result['kernel_gap_exact']=='1'
    assert result['Ginsparg_Wilson_identity_exact'] is True
    assert result['projected_CAR_defect_spectrum']['-1/4']==2
    assert result['projected_CAR_defect_spectrum']['1/4']==1
    assert result['control_is_actual_cubic_Wilson_kernel'] is False


def test_naive_mirror_count_from_independent_symbol_derivatives():
    k=s.symbols('k0:3',real=True)
    jac=s.diag(*[-s.cos(x) for x in k])
    rows=m.sine_regulator_control()['nodes']
    assert len(rows)==8
    assert sum(row['orientation'] for row in rows)==0
    for row in rows:
        assert jac.det().subs(dict(zip(k,[s.pi*x for x in row['corner_pi_units']])))==row['orientation']


def test_actual_raising_link_Wilson_gap_has_no_free_zero():
    result=m.raising_link_wilson_control()
    assert [row['dimension'] for row in result['rows']]==[4,8,16,32,64]
    assert result['overlap_D_squared_singular_value_lower']=='7/3'
    # Independent numerical sign reconstruction of the complete active64 block.
    x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1,-1]);eye=np.eye(2)
    gamma=[np.kron(x,eye),np.kron(y,x),np.kron(y,y),np.kron(y,z)]
    g5=np.kron(np.kron(z,eye),np.eye(16));u=(x-1j*y)/2
    dw=4*np.eye(64,dtype=complex)
    for j in range(4):
        ops=[u if k==j else eye for k in range(4)]
        link=ops[0]
        for op in ops[1:]:link=np.kron(link,op)
        dw-=np.kron((np.eye(4)-gamma[j])/2,link)+np.kron((np.eye(4)+gamma[j])/2,link.conj().T)
    h=g5@(dw-np.eye(64));ev,vec=np.linalg.eigh(h)
    assert min(abs(ev))>=1-1e-13
    sign=(vec*np.sign(ev))@vec.conj().T;d=np.eye(64)+g5@sign
    assert min(np.linalg.svd(d,compute_uv=False))**2>=7/3-1e-12
    assert np.linalg.norm(g5@d+d@g5-d@g5@d)<1e-12


def test_public_operator_map_detached(ops):
    fresh=m.quantum_higgs_operators();fresh['V']=None
    assert m.quantum_higgs_operators()['V']==ops['V']


def test_json_certificate_complete_duties_and_mutation_detached():
    result=m.certificate();original=json.dumps(result,sort_keys=True)
    assert [row['original_mean_trace'] for row in result['memory_duties']]==['1','1/12','1/6','1/72']
    assert all(row['coarse_rows']==48 for row in result['memory_duties'])
    assert result['physical_particle_or_mass_assigned'] is False
    assert result['interacting_Lorentz_continuum_proved'] is False
    result['quantum_Higgs']['rank_one_family_vertex_norm_exact']='999'
    assert json.dumps(m.certificate(),sort_keys=True)==original
