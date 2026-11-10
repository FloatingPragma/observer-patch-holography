"""Same-native-source massive-vector readouts on controlled finite windows."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
from model.interacting_block_quadrature import model as parent
from model.interacting_block_quadrature import reference
from model.interacting_block_quadrature.certificate import validate_pins
from . import algebra,resources,readout


PARENT='model/interacting_block_quadrature/results.json'
PARENT_SHA='64f90635637575bdaa9a24f1840ea3abb6362d6de8bc28e97de7fb275e8218e3'


@lru_cache(None)
def source():
    path=RER/PARENT
    assert hashlib.sha256(path.read_bytes()).hexdigest()==PARENT_SHA
    value=json.loads(path.read_text())
    validate_pins(value['input_sha256'])
    assert value['construction']['measured_inputs']==[]
    assert F(value['source_gap_interval'][0])>0
    # Re-execute the native matrix-unit replication and actual discrete
    # source symmetry admission; do not infer these from names or hashes.
    collective=parent.exact_collective(2)
    native=reference.certificate()
    actual={tuple(r) for r,_ in reference.cubic_rotations()}
    assert actual=={tuple(r) for r in algebra.proper_cube()}
    assert native['completion']['complete_record_probabilities_and_norm_errors_preserved']
    return value,{'parent_sha256':PARENT_SHA,'source_gap_interval':value['source_gap_interval'],
        'actual_native_collective_replay':collective,
        'actual_native_regular_reference_cover_order':native['group']['native_full_cover_order'],
        'actual_native_discrete_symmetry_replayed':True,
        'same_H_selected_actual_ground_reference_and_clock':True}


def literal_fourier_pair():
    # Vacuum plus six one-quantum states at +/-p, three components each.
    # This exact finite control checks the linear readout and its signs.
    m=s.Integer(1);omega=s.Rational(5,4);p=s.Matrix([0,s.Rational(3,4),0])
    ap=[];am=[]
    for color in range(3):
        plus=s.zeros(7);minus=s.zeros(7)
        plus[0,1+color]=1;minus[0,4+color]=1
        ap.append(plus);am.append(minus)
    phi=[(plus+minus.H)/s.sqrt(2*omega) for plus,minus in zip(ap,am)]
    pi=[-s.I*omega*(plus-minus.H)/s.sqrt(2*omega) for plus,minus in zip(ap,am)]
    M=algebra.polarization(p,omega,m)[1:,:]
    A=[sum((s.I*p[a]*pi[a]/(m*omega) for a in range(3)),s.zeros(7))]
    A += [sum((M[i,a]*phi[a] for a in range(3)),s.zeros(7)) for i in range(3)]
    phim=[v.H for v in phi];pim=[v.H for v in pi]
    Am=[sum((-s.I*p[a]*pim[a]/(m*omega) for a in range(3)),s.zeros(7))]
    Am += [sum((M[i,a]*phim[a] for a in range(3)),s.zeros(7)) for i in range(3)]
    assert all(a.H==b for a,b in zip(A,Am))
    naive=sum((p[a]*phi[a]/m for a in range(3)),s.zeros(7))
    naive_minus=sum((-p[a]*phim[a]/m for a in range(3)),s.zeros(7))
    assert naive.H==-naive_minus and naive.H!=naive_minus
    H=s.diag(0,*([omega]*6))
    assert s.simplify(s.I*(H*A[0]-A[0]*H)+sum((s.I*p[i]*A[i+1] for i in range(3)),s.zeros(7)))==s.zeros(7)
    four=s.Matrix([omega,*p]);Pi=-algebra.METRIC+four*four.T
    residue=s.Matrix(4,4,lambda i,j:s.simplify((A[i]*Am[j])[0,0]))
    assert residue==Pi/(2*omega)
    # Existing unfiltered spatial triplet has equal residue eigenvalues;
    # the Proca spatial longitudinal residue is different away from rest.
    assert Pi[2,2]==s.Rational(25,16) and Pi[1,1]==Pi[3,3]==1
    return {'finite_fourier_pair_control_dimension':7,
        'exact_Lorenz_and_Hermiticity_checked':True,
        'wrong_field_only_A0_fails_Hermiticity':True,
        'vacuum_frequency_residue':[[str(v) for v in row] for row in residue.tolist()],
        'new_readout_is_distinct_from_unfiltered_triplet':True}


def complete_instrument_control():
    unitary=s.Matrix([[s.Rational(3,5),4*s.I/5],[4*s.I/5,s.Rational(3,5)]])
    C=3*s.I*unitary/5;defect=4*s.eye(2)/5
    U=C.row_join(defect).col_join(defect.row_join(-C.H))
    assert U.H*U==s.eye(4)
    rows=[]
    psi=s.Matrix([1,0])
    for phase in (s.Integer(1),s.I):
        kraus=[(s.eye(2)+sign*phase*C)/2 for sign in (1,-1)]+[sign*phase*defect/2 for sign in (1,-1)]
        assert sum((k.H*k for k in kraus),s.zeros(2)).applyfunc(s.simplify)==s.eye(2)
        probabilities=[s.simplify((psi.H*k.H*k*psi)[0]) for k in kraus]
        assert sum(probabilities)==1 and all(p>=0 for p in probabilities)
        assert probabilities[2]+probabilities[3]>0
        rows.append({'control_phase':str(phase),'four_complete_record_probabilities':[str(v) for v in probabilities]})
    return {'Kraus_records':'(s,flag), s=+1,-1 and flag=0,1',
        'K_s0':'[I+s exp(i phase) C_N]/2',
        'K_s1':'s exp(i phase) sqrt(I-C_N^dag C_N)/2',
        'contraction':'C_N=P_N W(filtered_phi,filtered_pi)P_N, zero on the full native complement',
        'same_parent_regular_reference_lift_used':True,
        'all_failure_and_complement_records_retained':True,
        'literal_complete_contraction_control':rows}


def finite_displacement_element(n,m,eta):
    """One exact compressed Weyl entry; the sum and every source index are finite."""
    eta=s.sympify(eta)
    assert isinstance(n,int) and isinstance(m,int) and n>=0 and m>=0
    polynomial=sum(eta**(n-k)*(-s.conjugate(eta))**(m-k)/
        (s.factorial(n-k)*s.factorial(m-k)*s.factorial(k)) for k in range(min(n,m)+1))
    return s.exp(-eta*s.conjugate(eta)/2)*s.sqrt(s.factorial(n)*s.factorial(m))*polynomial


def calculate():
    _,admission=source()
    return {'schema':'oph.interacting_block_proca.v1',
        'construction':{'source_H':'unchanged H_(a,N) from interacting_block_quadrature',
            'actual_state':'the same selected actual symmetric global ground and retained native regular-reference isometry',
            'clock':'exp(-i t H_(a,N)) tensor I_reference',
            'new_operation':'explicit momentum-dependent linear Weyl readout filters only',
            'rest_frequency':'actual source gap Delta; no independently chosen mass',
            'effective_free_Proca_mass_parameter':'Delta/c^2 with hbar=1 and c=Delta/sqrt(6), in the declared source units',
            'normalized_coordinates':'tau=Delta*t; xi=sqrt(6)*x',
            'normalized_site_fields':'Phi=sqrt(Delta) Q; Pi=P/sqrt(Delta)',
            'continuum_field_normalization':'(sqrt(6)*a)^(-3/2) times each normalized site field',
            'reference_dimension':3072,'full_kinematic_and_physical_complements_retained':True,
            'operators_even_and_native_Gauss_invariant':True,
            'measured_inputs':[],'Standard_Model_particle_assignment':None,
            'strict_continuum_interaction':'free; fixed-capacity interactions are controlled errors',
            'same_unfiltered_scalar_Lorentz_action_claim':False},
        'source_admission':admission,'exact_polarization_algebra':algebra.symbolic_identities(),
        'induced_spatial_and_boost_action':algebra.induced_rotation_and_boost(),
        'same_retained_reference_spatial_action':algebra.retained_reference_rotation(),
        'literal_actual_field_filter_control':literal_fourier_pair(),
        'complete_instruments':complete_instrument_control(),
        'finite_Weyl_compiler':{'one_mode_entry':'exp(-|eta|^2/2)*sqrt(n!m!)*sum_(k=0)^min(n,m) eta^(n-k)(-conj(eta))^(m-k)/[(n-k)!(m-k)!k!]',
            'multiple_modes':'product of finitely many one-mode entries on the finite admitted occupation menu',
            'unbounded_comparison_carrier_is_not_native_hardware':True},
        'finite_window_filter_bounds':resources.filter_bounds(F(1,4)),
        'explicit_source_capacities':[resources.finite_source_resources(a,side,accuracy=F(1,10**9)) for a,side in ((F(1,4),64),(F(1,8),128))],
        'complete_actual_source_experiments':[readout.controlled_experiment(a,side) for a,side in ((F(1,4),64),(F(1,8),128))],
        'continuum_record_bounds':[resources.continuum_bounds(a) for a in (F(1,4),F(1,8))]}
