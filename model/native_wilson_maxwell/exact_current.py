"""Independent exact radical projector for the cube Wilson current."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from .field import cube
from model.native_block_dirac.algebra import pauli
from model.affine_gauss_mobile.model import ivalue,rounded
import mpmath as mp

def complex_fraction(value):
    out=[]
    for c in s.expand(value).as_real_imag():
        assert c.is_Rational;out.append(F(int(c.p),int(c.q)))
    return tuple(out)

def multiply(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def add(z,w):return z[0]+w[0],z[1]+w[1]

@lru_cache(None)
def literal_matter():
    """Independent Kronecker assembly, using the actual oriented cube."""
    sigmas=pauli();beta=s.kronecker_product(sigmas[2],s.eye(2))
    alpha=[s.kronecker_product(sigmas[0],sigma) for sigma in sigmas]
    transfers=[s.zeros(8) for _ in range(3)];currents=[]
    for tail,head,axis in cube()['edges']:
        T=s.zeros(8);T[tail,head]=1;transfers[axis]+=T
        J=s.I*s.kronecker_product(T-T.T,beta)/2-s.kronecker_product(T+T.T,alpha[axis])/2
        assert J==J.H and s.trace(J*J)==4
        currents.append(s.ImmutableSparseMatrix(J))
    H=s.Rational(5,2)*s.kronecker_product(s.eye(8),beta)
    for j,T in enumerate(transfers):
        H-=s.kronecker_product(T+T.T,beta)/2+s.I*s.kronecker_product(T-T.T,alpha[j])/2
    assert H==H.H
    x=s.Symbol('x');minimal=s.Poly(s.prod(((x+s.Rational(sign,2))**2-r)
      for r in (3,7,13) for sign in (-1,1)),x)
    dm=DomainMatrix.from_Matrix(H).convert_to(s.QQ_I)
    power=DomainMatrix.eye(32,s.QQ_I);powers=[power.to_Matrix()]
    for _ in range(12):power=power*dm;powers.append(power.to_Matrix())
    assert sum((minimal.nth(j)*powers[j] for j in range(13)),s.zeros(32))==s.zeros(32)
    return s.ImmutableMatrix(H),tuple(currents),minimal,tuple(powers)

@lru_cache(None)
def calculate():
    H,Js,minimal,powers=literal_matter();z=s.Symbol('z');parameters=(3,7,13)
    factors={r:s.Poly(z*z-(2*r+s.Rational(1,2))*z+(r-s.Rational(1,4))**2,z) for r in parameters}
    Q=s.prod(factors.values());matrices=[];projectors=[]
    for r in parameters:
        other=Q.exquo(factors[r]);E=(other*s.invert(other,factors[r])).rem(Q)
        assert (E*E-E).rem(Q).is_zero
        numerator=s.Poly((3*r+s.Rational(1,4)-z)/(2*(r-s.Rational(1,4))),z)
        R=(E*numerator).rem(Q)
        assert (s.Poly(z,z)*R*R-r*E).rem(Q).is_zero
        M=sum((R.nth(k)*powers[2*k+1] for k in range(6)),s.zeros(32))
        P=sum((E.nth(k)*powers[2*k] for k in range(6)),s.zeros(32))
        assert M==M.H and P==P.H and s.trace(M)==0
        matrices.append(s.ImmutableMatrix(M));projectors.append(s.ImmutableMatrix(P))
    assert sum(projectors,s.zeros(32))==s.eye(32)
    assert all(s.trace(M*J)==0 for M in matrices for J in Js)
    domains=[DomainMatrix.from_Matrix(M).convert_to(s.QQ_I) for M in matrices]
    for i,M in enumerate(domains):
        for j,N in enumerate(domains):
            expected=parameters[i]*projectors[i] if i==j else s.zeros(32)
            assert (M*N).to_Matrix()==expected
    # sign(H)=sum_r M_r/sqrt(r). The six positive roots are sqrt(r)±1/2,
    # and interpolation gives H f(H²)=+1 there and -1 on opposite roots.
    entries=[{key:complex_fraction(value) for key,value in J.todok().items()} for J in Js]
    rational_matrices=[[[complex_fraction(M[i,j]) for j in range(32)] for i in range(32)] for M in matrices]
    answers=[];all_coefficients=[]
    for weight in (cube()['P4'],cube()['P6']):
        contractions=[]
        for r in range(3):
            line=[]
            for t in range(3):
                answer=(F(0),F(0))
                for e in range(12):
                    for f in range(12):
                        coefficient=F(int(weight[e,f].p),int(weight[e,f].q))
                        if not coefficient:continue
                        value=(F(0),F(0))
                        # Direct four-index contraction, independent of the
                        # interval producer's rounded P*J support algorithm.
                        for (j,k),Je in entries[e].items():
                            for (l,i),Jf in entries[f].items():
                                term=multiply(multiply(rational_matrices[r][i][j],Je),
                                              multiply(rational_matrices[t][k][l],Jf))
                                value=add(value,term)
                        answer=add(answer,(coefficient*value[0],coefficient*value[1]))
                assert answer[1]==0;line.append(s.Rational(answer[0].numerator,answer[0].denominator))
            contractions.append(line)
        assert s.Matrix(contractions)==s.Matrix(contractions).T
        # The terms linear in sign(H) cancel after the symmetric e/f sum.
        value=s.trace(weight)-sum(contractions[r][t]/(4*s.sqrt(parameters[r]*parameters[t]))
          for r in range(3) for t in range(3))
        answers.append(s.simplify(value));all_coefficients.append(contractions)
    sigma=s.simplify(answers[0]/4+answers[1]/(2*s.sqrt(6)))
    mp.iv.dps=80;enclosure=rounded(ivalue(sigma),places=60)
    assert F(enclosure[0])>F(2043,2000)**2
    return {'exact_sign_H':'sum_(r=3,7,13) M_r/sqrt(r); rational M_r M_s=delta_rs r E_r',
      'exact_filled_projector':'(I-sign_H)/2','rank':16,
      'complete_CRT_projector_ranks':[int(s.trace(P)) for P in projectors],
      'exact_transverse_current_contractions':all_coefficients,
      'P4_current_square_contraction':answers[0],'P6_current_square_contraction':answers[1],
      'exact_sigma_squared_over_m_squared':sigma,'outward_interval':enclosure,
      'strict_sigma_over_m_lower':'2043/2000','all_twelve_filled_vacuum_current_means':'0',
      'diagnostic_sigma_squared_over_m_squared':float(sigma),
      'binomial_sign_approximation_or_numerical_eigensolver_used':False}

if __name__=='__main__':
    import json
    print(json.dumps(calculate(),indent=2,default=str))
