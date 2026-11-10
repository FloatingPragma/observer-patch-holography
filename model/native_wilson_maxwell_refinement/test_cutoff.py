"""Independent finite matrix controls and exact source-bound checks."""
from fractions import Fraction as F
import importlib.util
import math
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_maxwell_cutoff_producer',HERE/'producer.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
b=p.bounds


def test_rational_enclosures_and_nonzero_parent_transfer():
    for value in (F(2),F(6),F(1,7),F(24)):
        assert b.root(value,False)**2<=value<=b.root(value)**2
    exp=b.exp_upper(F(3,4))
    partial=sum((F(3,4)**j/math.factorial(j) for j in range(120)),F(0))
    assert partial<exp and float(exp-partial)<1e-70
    saved=p.source();cut=b.cutoff_bound(3)
    result=p.record_transfer(saved,cut)['exact']
    assert F(cut['state_vector_error_upper'])<F(1863,10**12)
    assert F(result['comparison_actual_target_probability_lower'])>F(475,10**10)
    assert F(result['comparison_actual_binary_connected_lower'])>F(471,10**10)


def test_completed_H_and_literal_exponential_compression_countercontrol():
    # Exact single-mode block; boundary completion is independent of mode frame.
    for cap in range(1,7):
        a=s.zeros(cap+1)
        for n in range(1,cap+1):a[n-1,n]=s.sqrt(n)
        N=s.diag(*range(cap+1));top=s.zeros(cap+1);top[cap,cap]=1
        assert (a*a.H+a.H*a)/2+top*(N+s.eye(cap+1))/2==N+s.eye(cap+1)/2
    q=s.zeros(4)
    for n in range(1,4):q[n-1,n]=q[n,n-1]=s.sqrt(n)
    assert (q**6)[0,0]==15
    assert (q**8)[0,0]==81
    assert (q**8)[0,0]!=105  # Full vacuum eighth moment, hence exponentials differ.


def occupations(modes,cap):
    def recurse(left,slots,prefix=()):
        if slots==1:
            yield prefix+(left,);return
        for n in range(left+1):yield from recurse(left-n,slots-1,prefix+(n,))
    return tuple(state for number in range(cap+1) for state in recurse(number,modes))


def coordinates(cap,coefficients):
    states=occupations(len(coefficients[0]),cap);index={n:i for i,n in enumerate(states)}
    result=[]
    for c in coefficients:
        rows=[];cols=[];values=[]
        for col,state in enumerate(states):
            for j,n in enumerate(state):
                if n:
                    target=list(state);target[j]-=1;row=index[tuple(target)]
                    value=c[j]*math.sqrt(n)
                    rows.extend((row,col));cols.extend((col,row));values.extend((value,value.conjugate()))
        result.append(coo_matrix((values,(rows,cols)),shape=(len(states),len(states))).tocsr())
    return states,result


def test_actual_five_mode_noncommuting_word_control():
    V,nu=p.algebra.field.frame()
    c=np.array([[complex(s.N(V[e,j]/s.sqrt(2*nu[j]),18)) for j in range(5)] for e in (0,1,5,9)])
    g=F(1,64);cap=3;large=8
    smallstates,qs=coordinates(cap,c);largestates,ql=coordinates(large,c)
    vs=np.zeros(len(smallstates),complex);vl=np.zeros(len(largestates),complex);vs[0]=vl[0]=1
    for e,sign in ((0,1),(1,-1),(2,1),(3,1),(0,-1)):
        vs=expm_multiply(sign*1j*float(g)*qs[e],vs)
        vl=expm_multiply(sign*1j*float(g)*ql[e],vl)
    embedded=np.zeros_like(vl);lookup={state:i for i,state in enumerate(largestates)}
    for index,state in enumerate(smallstates):embedded[lookup[state]]=vs[index]
    observed=la.norm(vl-embedded)
    r=5*b.displacement_radius(g)
    theorem=float(b.coherent_word_error(cap,r)+b.coherent_word_error(large,r))
    assert observed>1e-11 and observed<theorem
    # A full five-mode total cutoff, with nonparallel consecutive link directions.
    assert len(smallstates)==56 and len(largestates)==1287


def test_retained_number_measurement_and_controlled_displacement_history():
    def run(cap):
        a=np.diag(np.sqrt(np.arange(1,cap+1)),1);q=a+a.T
        W=la.expm(.07j*q);dimension=cap+1
        # Ancilla 0/1 coherently controls a first displacement; a second
        # register stores even/odd photon number without postselection.
        state=np.zeros((2,2,dimension),complex);state[0,0,0]=state[1,0,0]=1/math.sqrt(2)
        state[1,0]=W@state[1,0]
        for n in range(dimension):
            if n%2:state[:,1,n],state[:,0,n]=state[:,0,n].copy(),state[:,1,n].copy()
        for control in range(2):
            for record in range(2):state[control,record]=la.expm((-.09j if record else .04j)*q)@state[control,record]
        return state
    small=run(3);large=run(16);embedded=np.zeros_like(large);embedded[:,:,:4]=small
    R=F(7,100)+F(9,100)
    bound=b.instrument_word_error(3,R)+b.instrument_word_error(16,R)
    assert la.norm(large-embedded)<float(bound)
    assert abs(np.linalg.norm(small)-1)<1e-14
    assert abs(np.linalg.norm(large)-1)<1e-14


def test_exact_cube_locality_and_current_sign():
    data=p.algebra.calculate()
    assert data['exact_transverse_plus_longitudinal_identity']
    assert s.Rational(1,2)-s.sympify(data['electric_vacuum_variance'])>0
    assert 2-s.sympify(data['magnetic_vacuum_variance'])>0
    # Literal two-mode CAR current on a single directed edge.
    a=s.Matrix([[0,1],[0,0]]);z=s.diag(1,-1)
    tail=s.kronecker_product(a,s.eye(2));head=s.kronecker_product(z,a)
    Fhop=head.H*tail;n0=tail.H*tail;n1=head.H*head
    assert n0*Fhop-Fhop*n0==-Fhop and n1*Fhop-Fhop*n1==Fhop
    H=Fhop+Fhop.H;flow=-s.I*(Fhop-Fhop.H)
    assert s.I*(H*n0-n0*H)==-flow
    assert s.I*(H*n1-n1*H)==flow


def test_cutoff_refinement_and_common_complete_alphabet():
    errors=[F(b.cutoff_bound(k,order=k+24)['state_vector_error_upper']) for k in range(1,9)]
    assert all(a>b_ for a,b_ in zip(errors,errors[1:]))
    assert errors[-1]<F(4,10**20)
    for cap in range(1,8):
        finite=290*(cap+1)+1;comparison=290*(cap+2)+1
        assert comparison-finite==290
    assert 290*4+1==1161


def test_fresh_certificate_without_numerical_eigensolvers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('numerical eigen service forbidden')
    for namespace in (np.linalg,la):
        for name in ('eig','eigh','eigvals','eigvalsh'):monkeypatch.setattr(namespace,name,forbidden)
    value=p.certificate()
    assert value['parent_finite_source_sha256']==p.PARENT_SHA
    assert value['exact_local_comparison_algebra']['exact_Gauss_and_electric_energy_decomposition']
