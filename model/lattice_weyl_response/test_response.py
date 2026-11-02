from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from . import algebra,bounds,controls,packet


def test_exact_shell_moments_and_H4_boundary():
    result=algebra.sobolev()
    assert result['H1_H2_H3_H4_integer_majorants']==[5,36,352,5036]
    assert result['exact_normalized_momentum_moments']==[
        '90678696/4288237','5320580800/4288237','529844446656/4288237','108752906109312/4288237']
    r=s.symbols('r');phi=r*r*(r-1)**4*(2-r)**4
    assert all(s.diff(phi,r,m).subs(r,x)==0 for m in range(4) for x in (1,2))
    assert s.diff(phi,r,4).subs(r,1)!=0 # fifth weak derivative has a boundary distribution


def test_moment_integration_against_independent_polynomial_antiderivative():
    # The ell=1 radial Laplacian squared, independently assembled.
    r=s.symbols('r',positive=True);phi=r*r*(r-1)**4*(2-r)**4
    lap=lambda f:s.diff(f,r,2)+2*s.diff(f,r)/r-2*f/r**2
    density=s.expand(r*r*lap(lap(phi))**2)
    integrated=s.integrate(density,(r,1,2))
    assert algebra.laurent_integral(density,r)==F(int(s.numer(integrated)),int(s.denom(integrated)))
    with pytest.raises(ArithmeticError):algebra.laurent_integral(1/r,r)


def test_no_gap_is_required_on_lattice_side_of_sylvester():
    rng=np.random.default_rng(941)
    for dimension in (2,3,5):
        for gap in (.001,.25,2.):
            # Diagonal bases are arbitrary; R need not commute with either generator.
            left=np.arange(dimension,dtype=float) # exact zero mode is allowed
            right=gap+np.arange(dimension,dtype=float)/2
            R=rng.normal(size=(dimension,dimension))+1j*rng.normal(size=(dimension,dimension))
            X=R/(left[:,None]+right[None,:])
            assert np.linalg.norm(left[:,None]*X+X*right[None,:]-R)<1e-12
            assert np.linalg.norm(X,2)<=np.linalg.norm(R,2)/gap*(1+1e-12)
    # Deleting epsilon from the estimate is false, even in dimension one.
    assert 1/F(1,100)>1


def test_central_difference_residual_at_all_eight_nodes():
    rng=np.random.default_rng(942)
    sigma=[np.array(x).astype(complex) for x in algebra.pauli()]
    for bits in np.ndindex(2,2,2):
        for a in (.01,.1,.7):
            p=rng.normal(size=3)*3
            discrete=sum(((-1)**bits[j]*np.sin(a*p[j])/a)*sigma[j] for j in range(3))
            continuum=sum(((-1)**bits[j]*p[j])*sigma[j] for j in range(3))
            assert np.linalg.norm(discrete-continuum,2)<=a*a*sum(abs(p)**3)/6+1e-12


def test_full_Fock_wick_and_noncommuting_detector_controls():
    result=controls.certificate()
    assert result['exact_nonzero_Wick_response']=='-3/8'
    assert result['support_intrusion_rejected'] and result['unsymmetrized_detector_rejected']


def test_all_valleys_and_exact_response_frame():
    result=algebra.certificate()
    assert result['frames']['Weyl_fields']==8
    assert result['response_inverse_norm']=='68/5'
    assert result['response_determinant']=='24407490807/2000000000'
    assert result['off_diagonal_momentum_symmetric_stress_contraction_verified']


def test_joint_refinements_have_uniform_finite_inverse():
    for m in (0,1,5,100):
        rho=F(1,2**28*(m+1));q=bounds.admissible_spacing(rho)
        error=bounds.errors(q,rho)
        assert 0<error['Neumann_budget']<F(1,4)
        assert error['normalized_lattice_entry_upper']<rho**3
        assert error['normalization_defect_upper']<F(1,4)
    # Taking a -> 0 without a rate relative to R does not imply this theorem.
    assert bounds.errors(F(1,100),F(1,2**28))['Neumann_budget']>1


@pytest.mark.parametrize('value',[0,-1,True,1.,F(1,8)])
def test_bad_joint_radius_rejected(value):
    with pytest.raises(ValueError):bounds.admissible_spacing(value)


@pytest.mark.parametrize('q,rho',[(0,F(1,16)),(True,F(1,16)),(F(1,10),F(1,16)),
                                 (F(1,100),0),(F(1,100),.01)])
def test_bad_error_domain_rejected(q,rho):
    with pytest.raises(ValueError):bounds.errors(q,rho)


@pytest.fixture(scope='module')
def receipt():return packet.build()


@pytest.mark.parametrize('key',['nonzero_U_stress_convergence_claimed',
    'lattice_regional_modular_log_convergence_claimed','boosted_modular_caps_realized',
    'native_gate_program_for_shell_detector_supplied','area_Newton_or_stationarity_supplied',
    'complete_Einstein_register_inhabited','arbitrary_all_valley_state_total_stress_claimed'])
def test_unproved_endpoint_promotion_rejected(receipt,key):
    changed=deepcopy(receipt);changed['theorem'][key]=True
    with pytest.raises(ArithmeticError):packet.validate(changed)


def test_proof_constants_source_bindings_and_nested_type_mutations_rejected(receipt):
    variants=[]
    changed=deepcopy(receipt);changed['theorem']['source_interaction_U']=False;variants.append(changed)
    changed=deepcopy(receipt);changed['controls']['support_intrusion_rejected']=1;variants.append(changed)
    changed=deepcopy(receipt);changed['bounds']['constants']['K4']='5035';variants.append(changed)
    changed=deepcopy(receipt);changed['prior_custody']['prior_full_source_replay']=True;variants.append(changed)
    changed=deepcopy(receipt);changed['source_bindings'].pop(next(iter(changed['source_bindings'])));variants.append(changed)
    for value in variants:
        with pytest.raises(ArithmeticError):packet.validate(value)
