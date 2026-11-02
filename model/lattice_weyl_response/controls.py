"""Independent finite full-CAR and noncommuting detector checks."""
from functools import lru_cache
import sympy as s
from model.weyl_geometric_modular.algebra import pauli


def annihilators(modes):
    result=[]
    for mode in range(modes):
        matrix=s.zeros(2**modes)
        for column in range(2**modes):
            if column&(1<<mode):
                matrix[column^(1<<mode),column]=(-1)**((column&((1<<mode)-1)).bit_count())
        result.append(matrix)
    return result


def second_quantize(matrix,c):
    return sum((matrix[i,j]*c[i].adjoint()*c[j]
                for i in range(len(c)) for j in range(len(c))),s.zeros(c[0].rows))


@lru_cache(maxsize=1)
def certificate():
    # Two occupied orbitals in the complete four-mode Fock space.
    U=s.Matrix([[1,1],[1,-1],[1,s.I],[1,-s.I]])/2
    if U.adjoint()*U != s.eye(2):
        raise ArithmeticError('finite Slater orbitals are not orthonormal')
    psi=s.zeros(16,1)
    for i in range(4):
        for j in range(i+1,4):
            psi[(1<<i)+(1<<j)]=U.extract([i,j],[0,1]).det()
    c=annihilators(4);P=U*U.adjoint();f=s.eye(4)[:,3];F=P*f
    Q=s.Matrix([[2,s.I,1,0],[-s.I,-1,2*s.I,0],[1,-2*s.I,3,0],[0,0,0,0]])
    B=c[3].adjoint()*c[3]-P[3,3]*s.eye(16)
    observable=second_quantize(Q,c)
    direct=s.simplify((psi.adjoint()*B*observable*psi)[0])
    predicted=s.simplify(-(F.adjoint()*Q*F)[0])
    if direct != predicted or direct == 0 or B*observable != observable*B:
        raise ArithmeticError('full-Fock disjoint-support Wick response failed')
    # Intruding into the shell invalidates the negative rank-one formula.
    intrusion=s.simplify((psi.adjoint()*B*(c[3].adjoint()*c[3])*psi)[0])
    wrong=s.simplify(-abs(P[3,3])**2)
    if intrusion == wrong:
        raise ArithmeticError('missing-support-separation mutation was not detected')
    # A genuine derivative/spin/weight noncommuting example.
    p=s.Matrix([[0,-s.I,0],[s.I,0,-s.I],[0,s.I,0]])/2
    sx,_,sz=pauli();h=s.kronecker_product(p,sx)
    A=s.kronecker_product(s.eye(3),s.Rational(5,4)*s.eye(2)-s.Rational(3,4)*sz)
    D=s.Rational(5,4)*h-s.Rational(3,4)*s.kronecker_product(p,s.eye(2))
    w=s.kronecker_product(s.diag(1,2,0),s.eye(2))
    raw=w*A*D;correct=(raw+D*w*A)/2
    if correct != correct.adjoint() or raw == raw.adjoint():
        raise ArithmeticError('symmetric stress detector control failed')
    return {'full_Fock_modes':4,'full_Fock_dimension':16,
            'exact_nonzero_Wick_response':str(direct),
            'support_intrusion_rejected':True,'unsymmetrized_detector_rejected':True,
            'symmetric_detector_verified':True,
            'finite_control_is_continuum_proof':False}
