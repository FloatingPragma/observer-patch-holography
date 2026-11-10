"""Open Wilson matter, exact filled-band projector control and vacuum current."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from .current import matrices
from . import field
from model.native_block_dirac.algebra import pauli

HERE=Path(__file__).resolve().parent


def fraction(x):
    x=s.sympify(x)
    assert x.is_Rational
    return F(int(x.p),int(x.q))


def sqrt_bound(x,upper=False,bits=160):
    x=F(x);scale=2**bits;k=isqrt(x.numerator*scale**2//x.denominator)
    if upper and k*k*x.denominator<x.numerator*scale**2:k+=1
    return F(k,scale)


def outward(x,upper=False,bits=160):
    x=F(x);scale=2**bits
    return F(-(-x.numerator*scale//x.denominator) if upper else x.numerator*scale//x.denominator,scale)


@lru_cache(None)
def cube_matter():
    data=field.cube();vertices=data['vertices'];edges=data['edges'];N=len(vertices)
    alpha,beta=matrices();H=s.zeros(4*N);currents=[];degree=[0]*N
    for tail,head,axis in edges:degree[tail]+=1;degree[head]+=1
    for vertex in range(N):
        H[4*vertex:4*vertex+4,4*vertex:4*vertex+4]=(1+s.Rational(degree[vertex],2))*beta
    for tail,head,axis in edges:
        B=-(beta+s.I*alpha[axis])/2
        H[4*tail:4*tail+4,4*head:4*head+4]=B
        H[4*head:4*head+4,4*tail:4*tail+4]=B.H
        J=s.zeros(4*N)
        # +phase belongs to head^dagger tail for D=head-tail.
        J[4*tail:4*tail+4,4*head:4*head+4]=-s.I*B
        J[4*head:4*head+4,4*tail:4*tail+4]=s.I*B.H
        assert J.H==J and J**3==J and s.trace(J*J)==4
        currents.append(s.ImmutableSparseMatrix(J))
    charge_conjugation=s.kronecker_product(s.eye(N),s.I*beta*alpha[1])
    assert charge_conjugation*charge_conjugation.H==s.eye(4*N)
    assert charge_conjugation*H.conjugate()*charge_conjugation.H==-H
    time_reversal=s.kronecker_product(s.eye(N),s.diag(s.I*pauli()[1],s.I*pauli()[1]))
    assert time_reversal*H.conjugate()*time_reversal.H==H
    assert all(time_reversal*J.conjugate()*time_reversal.H==-J for J in currents)
    x=s.Symbol('x'); factors=[4*x*x+sign*4*x-c for c in (11,27,51) for sign in (-1,1)]
    minimal=s.Poly(s.prod(factors)/4**6,x)
    characteristic=s.prod(factor**(4 if c==27 else 2)
        for c in (11,27,51) for factor in (4*x*x-4*x-c,4*x*x+4*x-c))/4**16
    assert s.expand(H.charpoly(x).as_expr()-characteristic)==0
    domain_H=DomainMatrix.from_Matrix(H).convert_to(s.QQ_I)
    power=DomainMatrix.eye(H.rows,s.QQ_I);powers=[power.to_Matrix()]
    for _ in range(minimal.degree()):
        power=power*domain_H;powers.append(power.to_Matrix())
    assert sum((minimal.nth(j)*powers[j] for j in range(13)),s.zeros(H.rows)).applyfunc(s.expand)==s.zeros(H.rows)
    return H,tuple(currents),minimal,tuple(powers)


def sign_polynomial(minimal,terms=512):
    """Exact binomial polynomial, reduced modulo the executed annihilator."""
    degree=minimal.degree(); q=[fraction(minimal.nth(j)) for j in range(degree)]
    def reduce(values):
        values=list(values)
        for k in range(len(values)-1,degree-1,-1):
            value=values[k]
            if value:
                for j in range(degree):values[k-degree+j]-=value*q[j]
        return values[:degree]+[F(0)]*max(0,degree-len(values))
    def times_r(values):
        out=list(values)+[F(0),F(0)]
        for j,value in enumerate(values):out[j+2]-=value/25
        return reduce(out)
    term=[F(1)]+[F(0)]*(degree-1);total=term.copy()
    for n in range(terms):
        term=[x*F(2*n+1,2*n+2) for x in times_r(term)]
        total=[x+y for x,y in zip(total,term)]
    coefficients=reduce([F(0)]+[x/5 for x in total])
    assert all(coefficients[j]==0 for j in range(0,degree,2))
    # |h|>=1 and ||h||<5: binomial coefficients<=1, r<=24/25.
    sign_error=25*F(24,25)**(terms+1)
    return coefficients,sign_error


def cadd(z,w):return (z[0]+w[0],z[1]+w[1])
def csub(z,w):return (z[0]-w[0],z[1]-w[1])
def cmul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
ZERO=(F(0),F(0))


@lru_cache(None)
def filled_current(terms=512):
    H,currents,minimal,powers=cube_matter();coeff,error=sign_polynomial(minimal,terms)
    projector=s.eye(H.rows)/2
    for n,c in enumerate(coeff):
        if c:projector-=s.Rational(c.numerator,c.denominator)*powers[n]/2
    P=[[ZERO for _ in range(H.cols)] for _ in range(H.rows)]
    for row in range(H.rows):
        for col in range(row,H.cols):
            exact=tuple(map(fraction,s.expand(projector[row,col]).as_real_imag()))
            real,imag=(outward(value,bits=128) for value in exact)
            if row==col:assert imag==0
            P[row][col]=(real,imag);P[col][row]=(real,-imag)
    Js=[];As=[];Bs=[];supports=[]
    for current in currents:
        entries={(i,j):tuple(map(fraction,s.expand(value).as_real_imag()))
                 for (i,j),value in current.todok().items()}
        support=sorted({i for i,j in entries});supports.append(support);Js.append(entries)
        product={}
        for (i,j),z in entries.items():
            for row in range(H.rows):
                key=(row,j)
                product[key]=cadd(product.get(key,ZERO),cmul(P[row][i],z))
        As.append(product)
        Bs.append({key:csub(entries.get(key,ZERO),value) for key,value in product.items()})
    weights=[field.cube()['P4'],field.cube()['P6']];answers=[]
    for weight in weights:
        answer=ZERO
        for e in range(len(currents)):
            for f in range(len(currents)):
                coefficient=fraction(weight[e,f])
                if not coefficient:continue
                value=ZERO
                for i in supports[f]:
                    for j in supports[e]:
                        value=cadd(value,cmul(As[e].get((i,j),ZERO),Bs[f].get((j,i),ZERO)))
                answer=cadd(answer,(coefficient*value[0],coefficient*value[1]))
        assert answer[1]==0;answers.append(answer[0])
    assert all(x>0 for x in answers)
    a,b=answers
    approximate_low=a/4+b/(2*sqrt_bound(6,True))
    approximate_high=a/4+b/(2*sqrt_bound(6))
    # Each complex entry changes by at most2*2^-128, so the complete
    #32-dimensional operator rounding error is at most64*2^-128.
    eps=error/2+F(64,2**128)
    # ||J_e||_HS=2; sum_ef |Cov(Q)_ef| <=144/4=36.
    transfer=144*eps*(2+eps)
    lower=outward(approximate_low-transfer);upper=outward(approximate_high+transfer,True)
    assert lower>0
    return {'one_body_dimension':H.rows,'all_many_body_matter_modes':H.rows,
        'positive_free_levels_over_m':['sqrt(3)-1/2','sqrt(3)+1/2','sqrt(7)-1/2','sqrt(7)+1/2','sqrt(13)-1/2','sqrt(13)+1/2'],
        'positive_level_multiplicities':[2,2,4,4,2,2],
        'executed_complete_annihilator':str(minimal.as_expr()),
        'full_one_body_charge_conjugation_and_Hermiticity_checked':True,
        'sign_binomial_last_index':terms,'sign_uniform_error_upper':str(outward(error,True)),
        'projector_uniform_error_upper':str(outward(eps,True)),
        'Hermitian_projector_entry_rounding_bits':128,
        'complete_operator_rounding_error_upper':str(F(64,2**128)),
        'current_square_transfer_error_upper':str(outward(transfer,True)),
        'pair_one_photon_V1_vacuum_norm_squared_over_m_squared_interval':[str(lower),str(upper)],
        'pair_one_photon_V1_vacuum_norm_over_m_lower':str(sqrt_bound(lower)),
        'diagnostic_current_square_interval':[float(lower),float(upper)],
        'filled_vacuum_mean_current_all_twelve_edges':'0',
        'full_matter_Fock_dimension':2**H.rows,
        'full_Fock_Hamiltonian_matrix_numerically_enumerated':False,
        'numerical_one_body_eigenvectors_or_eigenvalues_used':False}


if __name__=='__main__':
    import json
    print(json.dumps(filled_current(),indent=2))
