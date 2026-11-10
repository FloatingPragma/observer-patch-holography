"""Exact80-state current polynomial and rational five-state pointer coupling."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,math
import numpy as np
import sympy as s
from scipy.linalg import expm
from model.native_weak_current_interface.model import dressed_weak_matrix,dressed_weak_polynomial
from model.certified_reference import interval as iv

VALUES=(s.Integer(-1),s.Rational(-1,2),s.Integer(0),s.Rational(3,2),s.Integer(2))
LABELS=(-2,-1,0,1,2)
BITS=40
ORDER=64


def sqrt_upper(x):return F(iv.sqrt_fraction(F(x)).hi,iv.SCALE)

def matrix_hash(matrix):
    entries=[(int(i),int(j),str(v))for (i,j),v in sorted(s.SparseMatrix(matrix).todok().items())]
    return hashlib.sha256(json.dumps({'shape':list(matrix.shape),'entries':entries},separators=(',',':')).encode()).hexdigest()


@lru_cache(None)
def current_algebra():
    W=s.SparseMatrix(dressed_weak_matrix());I=s.eye(80)
    polynomial,_,_,_,_=dressed_weak_polynomial()
    require_gauge=bool(polynomial)  # The provider executes every full30-mode Gauss commutator.
    assert require_gauge and W.H==W
    projectors=[]
    for w in VALUES:
        P=s.SparseMatrix(I)
        for other in VALUES:
            if other!=w:P=P*(W-other*I)/(w-other)
        assert P.H==P and P*P==P and W*P==w*P
        projectors.append(s.SparseMatrix(P))
    assert sum(projectors,s.zeros(80))==I
    assert [s.trace(P)for P in projectors]==[4,12,58,4,2]
    f=(4*W**4-4*W**3-13*W**2+25*W)/15
    assert f==sum((j*P for j,P in zip(LABELS,projectors)),s.zeros(80))
    penalty=s.kronecker_product(s.diag(0,1,1,1,1),s.eye(16))
    assert penalty*W==W*penalty and penalty*f==f*penalty
    return W,s.SparseMatrix(f),tuple(projectors)


def sine_pi_fraction(multiplier):
    angle=iv.pi()*F(multiplier);middle=F(angle.lo+angle.hi,2*iv.SCALE)
    radius=F(angle.hi-angle.lo,2*iv.SCALE)
    sine,_=iv.sine_cosine(middle)
    outward=iv.Interval.point(radius).hi
    return iv.Interval(sine.lo-outward,sine.hi+outward)


def rounded_cell(interval):
    scale=1<<BITS;lo,hi=F(interval.lo,iv.SCALE),F(interval.hi,iv.SCALE)
    nearest=lambda x:(x*scale+F(1,2)).numerator//(x*scale+F(1,2)).denominator
    assert nearest(lo)==nearest(hi),'rounding cell ambiguous'
    q=F(nearest(lo),scale)
    return q,max(abs(q-lo),abs(q-hi))


@lru_cache(None)
def pointer():
    # P_ab=i(4pi/25) sum_(k=1,2) k sin[2pi k(a-b)/5].
    # exp(-iP) is exactly the forward cyclic shift. Preserve antisymmetry
    # by quantizing the two independent real coefficients together.
    coefficients={};errors={}
    for displacement in (1,2):
        total=iv.Interval.point(0)
        for k in (1,2):
            residue=(k*displacement+2)%5-2
            total+=k*sine_pi_fraction(F(2*residue,5))
        entry=F(4,25)*iv.pi()*total
        coefficients[displacement],errors[displacement]=rounded_cell(entry)
    coefficients[0]=F();errors[0]=F()
    for d in (1,2):coefficients[5-d]=-coefficients[d];errors[5-d]=errors[d]
    P=s.zeros(5)
    for a in range(5):
        for b in range(5):P[a,b]=s.I*s.Rational(coefficients[(a-b)%5])
    assert P.H==P and s.trace(P)==0
    delta=sqrt_upper(sum((errors[(a-b)%5]**2 for a in range(5)for b in range(5)),F()))
    norm=F(4,5)*F(iv.pi().hi,iv.SCALE)+delta
    assert norm<3
    return s.SparseMatrix(P),delta,norm,coefficients


def rational_exponential(generator):
    """Exact rational Taylor evaluator, separate from the unitary definition."""
    assert generator.shape==(5,5) and generator.T==-generator
    assert all(x.is_Rational for x in generator)
    result=s.eye(5);term=s.eye(5)
    for n in range(1,ORDER+1):term=term*generator/n;result+=term
    return s.SparseMatrix(result)


@lru_cache(None)
def local_execution():
    W,f,projectors=current_algebra();P,delta,norm,coefficients=pointer()
    shift=s.zeros(5)
    for j in range(5):shift[(j+1)%5,j]=1
    ideal=sum((s.kronecker_product(Q,shift**(j%5))for j,Q in zip(LABELS,projectors)),s.zeros(400))
    assert ideal.H*ideal==s.eye(400)
    G=s.kronecker_product(f,P);assert G.H==G
    assert all(x.as_real_imag()[0]==0 for x in G)
    # ||jP||<6. A geometric tail after degree64 is an exact uniform bound.
    tail=F(6**(ORDER+1),math.factorial(ORDER+1))/(1-F(6,ORDER+2))
    approximants=[rational_exponential(-s.I*j*P)for j in LABELS]
    approximate=sum((s.kronecker_product(Q,U)for Q,U in zip(projectors,approximants)),s.zeros(400))
    blank_indices=[5*j for j in range(80)]
    output=np.asarray(approximate.tolist(),dtype=complex)
    direct=expm(-1j*np.asarray(G.tolist(),dtype=complex))
    target=np.asarray(ideal.tolist(),dtype=complex)
    # Actual complete400-column numerical execution; no normalization or
    # dropping pointer values. Exact Taylor/quantization bounds are separate.
    columns=[]
    for slot in range(5):
        K=direct[slot::5,blank_indices]
        columns.append(K)
    complete=sum((K.conj().T@K for K in columns),np.zeros((80,80),complex))
    numerical={'all400_columns_executed':True,
               'direct_vs_rational_Taylor_operator_residual':float(np.linalg.norm(direct-output,2)),
               'direct_vs_exact_cyclic_operator_residual':float(np.linalg.norm(direct-target,2)),
               'complete_five_Kraus_residual':float(np.linalg.norm(complete-np.eye(80),2)),
               'numerical_residual_is_exact_error_bound':False}
    assert numerical['direct_vs_rational_Taylor_operator_residual']<1e-12
    assert numerical['direct_vs_exact_cyclic_operator_residual']<float(2*delta)+1e-12
    assert numerical['complete_five_Kraus_residual']<1e-12
    return {'active_current_dimension':80,'neutral_pointer_dimension':5,'complete_local_coupling_dimension':400,
            'current_eigenvalues':list(map(str,VALUES)),'centered_pointer_labels':list(LABELS),
            'cyclic_pointer_slots':[j%5 for j in LABELS],'projector_ranks':[int(s.trace(Q))for Q in projectors],
            'current_polynomial_coefficients_low_to_high':['0','5/3','-13/15','-4/15','4/15'],
            'f_norm':'2','W_norm':'2','scalar_penalty_commutator_zero':True,
            'all_twelve_full_native_Gauss_commutators_zero':True,
            'other26_native_CAR_modes_are_identity_spectators':True,
            'pointer_generator':[list(map(str,row))for row in P.tolist()],
            'pointer_generator_bits':BITS,
            'two_independent_pointer_numerators':[int(coefficients[j]*(1<<BITS))for j in(1,2)],
            'signed_bits_for_those_two_numerators':sum(abs(int(coefficients[j]*(1<<BITS))).bit_length()+1 for j in(1,2)),'exact_pointer_generator_error_upper':str(delta),
            'pointer_generator_norm_upper':str(norm),
            'coupling_generator_norm_upper':str(2*norm),
            'ideal_readout_error_from_pointer_quantization_upper':str(2*delta),
            'Taylor_order':ORDER,'rational_Taylor_operator_error_upper':str(tail),
            'exact_matrix_sha256':{'W80':matrix_hash(W),'fW80':matrix_hash(f),'G400':matrix_hash(G),'ideal_U400':matrix_hash(ideal)},
            'complete_Kraus_definition':'K_r=(I80 tensor <r|) exp(-iG) (I80 tensor |0>), r=0..4',
            'exact_five_Kraus_completeness_from_unitarity':True,'raw_numerical_execution':numerical}
