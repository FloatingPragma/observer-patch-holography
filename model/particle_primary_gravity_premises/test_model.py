"""Independent source, finite entropy and boundary-domain controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import model as g


def _log(matrix):
    values,vectors=np.linalg.eigh(matrix)
    assert min(values)>0
    return (vectors*np.log(values))@vectors.conj().T


def _entropy(matrix):
    values=np.linalg.eigvalsh(matrix)
    return -sum(float(x*np.log(x)) for x in values if x>1e-14)


def _relative(left,right):
    return np.trace(left@(_log(left)-_log(right))).real


def _polynomial_density(H,beta,upper,order,scale):
    y=beta*(upper*np.eye(len(H))-H)/(2*scale)
    T=sum((np.linalg.matrix_power(y,k)/factorial(k) for k in range(order+1)),np.zeros_like(H))
    Q=np.linalg.matrix_power(T,scale)
    return Q@Q/np.trace(Q@Q)


def test_exact_polynomial_modular_derivative_and_monotone_defect():
    y=s.symbols('y',nonnegative=True)
    for m in (1,2,3,11):
        T=sum(y**k/s.factorial(k) for k in range(m+1))
        previous=sum(y**k/s.factorial(k) for k in range(m))
        assert s.diff(T,y)==previous
        assert s.expand(m*T-y*s.diff(T,y))==sum((m-k)*y**k/s.factorial(k) for k in range(m))
        for x in (F(),F(1,8),F(1,4),F(1,2)):
            ratio=F(s.diff(T,y).subs(y,s.Rational(x.numerator,x.denominator))/T.subs(y,s.Rational(x.numerator,x.denominator)))
            bound=g.polynomial_slope(1,0,1,m,1)
            assert bound['modular_slope_lower']<=ratio<=bound['modular_slope_upper']
    for bad in (0,-1,True):
        with pytest.raises(ValueError):g.polynomial_slope(1,0,1,bad,1)


def test_modular_work_sandwich_all_permutations_and_coherent_mixings():
    energies=np.array([-2.,-.5,1.,3.]);H=np.diag(energies)
    beta=.6;rho=_polynomial_density(H,beta,3,2,4);K=-_log(rho)
    bound=g.polynomial_slope('3/5',-2,3,2,4)
    lower=float(bound['modular_slope_lower']);upper=float(bound['modular_slope_upper'])
    unitaries=[np.eye(4)[list(p)] for p in permutations(range(4))]
    unitaries.append(np.array([[1j**(j*k)/2 for k in range(4)] for j in range(4)]))
    unitaries.append(expm(1j*np.array([[0,.3,0,.1],[.3,0,.4,0],[0,.4,0,.5],[.1,0,.5,0]])))
    for U in unitaries:
        changed=U@rho@U.conj().T
        work=np.trace((changed-rho)@H).real
        cost=np.trace((changed-rho)@K).real
        assert work>=-1e-14
        assert lower*work-1e-13<=cost<=upper*work+1e-13
        assert abs(cost-_relative(changed,rho))<1e-13
    slopes=np.diff(np.diag(K))/np.diff(energies)
    assert max(slopes)-min(slopes)>.001  # This control is not exact Gibbs.


def test_global_modular_clock_comparison_with_entangled_reference():
    H=np.diag([-2.,-.5,1.,3.]);rho=_polynomial_density(H,.6,3,2,4);K=-_log(rho)
    data=g.polynomial_slope('3/5',-2,3,2,4)
    defect=float(data['global_modular_affine_defect_range_upper'])
    residual=np.diag(K-.6*H);center=(min(residual)+max(residual))/2
    assert max(residual)-min(residual)<=defect+1e-14
    for t in (-2.,.3,1.):
        modular=expm(-1j*t*K);source=expm(-1j*t*.6*H)
        assert np.linalg.norm(modular-np.exp(-1j*t*center)*source,2)<=abs(t)*defect/2+1e-14
        # Complete output comparison on a maximally entangled source/reference.
        overlap=np.trace(modular.conj().T@source)/4
        half_trace=np.sqrt(max(0.,1-abs(overlap)**2))
        assert half_trace<=abs(t)*defect/2+1e-13


def test_regional_first_law_work_bound_and_complete_coin_entropy():
    H=np.diag([0.,1.,2.,5.]);beta=.5
    rho=_polynomial_density(H,beta,5,3,4)
    u=np.array([[1,1],[1,-1]])/np.sqrt(2);U=np.kron(u,np.eye(2))
    changed=U@rho@U.T
    region=lambda x:np.trace(x.reshape(2,2,2,2),axis1=1,axis2=3)
    sigma=region(rho);sigma_u=region(changed)
    work=np.trace((changed-rho)@H)
    full=_relative(changed,rho);local=_relative(sigma_u,sigma)
    assert 0<local<full<beta*work
    t=.25;mixture=(1-t)*sigma+t*sigma_u
    entropy_difference=_entropy(mixture)-_entropy(sigma)
    remainder=_relative(mixture,sigma)
    assert abs(entropy_difference+remainder-t*local)<1e-13
    assert 0<entropy_difference+remainder<=t*beta*work
    coin=np.zeros((4,4));coin[:2,:2]=(1-t)*sigma;coin[2:,2:]=t*sigma_u
    binary=-t*np.log(t)-(1-t)*np.log(1-t)
    assert abs(_entropy(coin)-binary-_entropy(sigma))<1e-13
    assert abs(_entropy(coin)-_entropy(mixture))>.1


def test_quantum_link_control_stationarity_and_actual_local_work_bound():
    control=g.stress.quantum_control();P=control['physical_projector']
    indices=[i for i in range(48) if P[i,i]==1]
    H=np.array(control['H'].extract(indices,indices),dtype=complex)
    A=np.array(control['h'][0].extract(indices,indices),dtype=complex)
    rho=_polynomial_density(H,1/16,257/30,3,2)
    assert np.linalg.norm(H@rho-rho@H)<1e-13
    k=np.linalg.norm(A,2);c=np.linalg.norm(H@A-A@H,2)
    slope=g.polynomial_slope('1/16','-257/30','257/30',3,2)
    for angle in (-.3,.1,.4):
        U=expm(1j*angle*A);changed=U@rho@U.conj().T
        work=np.trace((changed-rho)@H).real;cost=_relative(changed,rho)
        assert 0<work<=angle**2*k*c+1e-13
        assert float(slope['modular_slope_lower'])*work-1e-13<=cost<=work/16+1e-13
    # Actual Gauss-restricted control; no control spectrum is used as primary data.
    assert len(indices)==6
    assert not control['same_as_selected_spatial_primary']


def test_idle_purifier_cost_is_infinite_while_source_cost_is_finite():
    psi=s.Matrix([s.sqrt(s.Rational(2,3)),0,0,s.sqrt(s.Rational(1,3))])
    U=s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(2))
    out=U*psi;projector=psi*psi.H
    # Failure of support inclusion is the exact finite relative-entropy obstruction.
    assert s.simplify((out.H*(s.eye(4)-projector)*out)[0])>0
    rho=np.diag([2/3,1/3]);changed=np.diag([1/3,2/3])
    assert abs(_relative(changed,rho)-np.log(2)/3)<1e-14


def test_full_fock_spin_flip_caps_are_not_one_particle_caps():
    result=g.flip_control()
    assert len(result['directional_bond_controls'])==3
    levels=(-1,-1,0,0,0,0,1,1)
    spectrum=[sum(levels[i] for i in range(8) if mask&(1<<i)) for mask in range(256)]
    assert min(spectrum)==-2 and max(spectrum)==2
    assert max(levels)==1
    for current in g.channels.native.native_currents():
        assert s.trace(current)==0  # Filled M15 rows are actual gauge singlets.


@pytest.mark.parametrize('side',[2,3])
def test_open_boundary_coframe_identity_and_physical_extremizer(side):
    value=g.open_shift_control(side)
    assert value['missing_directional_face_slots']==3*side**2
    assert value['native_filled_negative_beta_expectation']==-45*side**2
    assert value['full_native_physical_norm_at_unit_spacing']==45*side**2


def test_open_boundary_identity_survives_nonconstant_unitary_link_values():
    side=2;cells=list(product(range(side),repeat=3));index={x:i for i,x in enumerate(cells)}
    T=s.zeros(8)
    for cell in cells:
        if cell[0]==0:
            other=(1,cell[1],cell[2]);phase=s.I**(1+cell[1]+2*cell[2])
            T[index[cell],index[other]]=phase
    I=s.eye(8);K=I-T;missing=I-T.H*T
    assert K.H*K==2*I-T-T.H-missing
    assert s.trace(missing)==4
    assert missing*missing==missing


@pytest.fixture(scope='module')
def member():return g.contract.member(0)


def test_authenticated_primary_modular_work_and_boundary_bounds(member):
    receipt=g.bind_member(member);g.verify_bound(member,receipt)
    result=receipt['result'];slope=result['modular_work_slope'];flip=result['nonzero_global_source_flip']
    assert F(slope['beta'])==22272
    assert F(slope['relative_modular_slope_defect_upper'])<F(1116,10**17)
    assert F(slope['global_modular_affine_defect_range_upper'])<F(325,10**8)
    assert F('18.05')<F(flip['cyclic_source_work_lower'])<F('18.06')
    assert F(flip['changed_hopping_norm_upper'])==48
    assert F(flip['cyclic_source_work_upper'])==144
    assert F(flip['complete_positive_work_record_probability_lower'])>F(9,16)
    assert F(flip['source_modular_relative_entropy_lower'])>400000
    assert flip['full_retained_purifier_relative_entropy']=='positive_infinity'
    assert F(result['local_regional_entropy_work']['adaptive_regional_modular_slope_upper'])<=F(1,64)
    assert result['open_coframe_boundary_obstruction']['full_Fock_and_physical_norm_exact']=='180'
    assert not result['proper_geometry_or_Einstein_equation_proved']
    assert not result['gravity_required_for_finite_mass_readout']


def test_local_entropy_bound_scales_quadratically_and_is_not_refinement_claim(member):
    first=g.local_entropy_work_bound(member,'1/1000');half=g.local_entropy_work_bound(member,'1/2000')
    assert first['cyclic_source_work_upper']==4*half['cyclic_source_work_upper']
    assert len(first['generator_owner_support'])==4
    assert first['flux_energy_index_region']==[[0,0,0]]
    assert not first['regional_modular_generator_equals_beta_local_energy']
    later=g.local_entropy_work_bound(g.contract.member(1))
    assert later['generator_norm_upper']>first['generator_norm_upper']
    assert later['adaptive_source_angle']<first['adaptive_source_angle']


def test_source_state_reference_clock_and_receipt_mutations_rejected(member):
    receipt=g.bind_member(member)
    altered=deepcopy(receipt);altered['result']['proper_geometry_or_Einstein_equation_proved']=True
    with pytest.raises(ValueError):g.verify_bound(member,altered)
    for section,key,value in (('clock','proper_time',0),('reference','idle_H',1),('state','beta',0)):
        changed=deepcopy(member);changed[section][key]=value
        with pytest.raises(ValueError):g.bind_member(changed)


def test_premise_map_symbols_exist_and_do_not_admit_primary_geometry():
    mapping=g.premise_map();rows=mapping['theorem_inputs']
    assert len(rows)==7
    assert len([r for r in rows if r['status']=='proved_on_primary'])==2
    assert mapping['canonical_target_node']=='geometric_gravity'
    assert not mapping['gravity_required_for_finite_mass_readout']
    assert all((g.contract.ROOT/name).is_file() for name in g.SOURCE_INPUTS)
