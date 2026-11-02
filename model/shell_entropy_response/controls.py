"""Finite numerical cross-checks; certified constants live in bounds.py."""
import numpy as np
from scipy.integrate import quad
from model.shrinking_modular_entropy.algebra import finite_example


def finite_scalar(value):
    result=float(value)
    if not np.isfinite(result):raise ArithmeticError('nonfinite entropy diagnostic')
    return result


def log_matrix(matrix):
    a=np.asarray(matrix,dtype=complex)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not np.isfinite(a).all() or not np.allclose(a,a.conj().T,rtol=0,atol=1e-13):
        raise ValueError('finite Hermitian matrix required')
    values,vectors=np.linalg.eigh(a)
    if not np.isfinite(values).all() or not np.isfinite(vectors).all():
        raise ArithmeticError('nonfinite eigendecomposition')
    if np.min(values)<=0:raise ValueError('positive definite matrix required')
    result=(vectors*np.log(values))@vectors.conj().T
    if not np.isfinite(result).all():raise ArithmeticError('nonfinite matrix logarithm')
    return result


def relative_entropy(left,right):
    return finite_scalar(np.trace(left@(log_matrix(left)-log_matrix(right))).real)


def log_hessian(matrix,direction):
    if not np.isfinite(matrix).all() or not np.isfinite(direction).all():
        raise ValueError('finite Hessian inputs required')
    values,vectors=np.linalg.eigh(matrix)
    if not np.isfinite(values).all() or not np.isfinite(vectors).all():
        raise ArithmeticError('nonfinite eigendecomposition')
    if np.min(values)<=0:raise ValueError('positive definite matrix required')
    d=vectors.conj().T@direction@vectors
    total=0.0
    for i,x in enumerate(values):
        for j,y in enumerate(values):
            slope=1/x if abs(x-y)<1e-13 else (np.log(x)-np.log(y))/(x-y)
            total+=abs(d[i,j])**2*slope
    return finite_scalar(total)


def finite_checks():
    e=finite_example();A=np.array(e['A'],dtype=float);f=np.array(e['F'],dtype=float).ravel()
    X=np.outer(f,f);I=np.eye(2);theta=e['theta'];rho0=np.array(e['rho'].subs(theta,0),dtype=float)
    tangent=np.array(e['rho'].diff(theta),dtype=float);m=float(f@f)
    def scalars(t):
        return finite_scalar(f@np.linalg.solve(A+t*I,f)),finite_scalar(f@np.linalg.solve(I-A+t*I,f))
    def chi(t0):
        def integrand(t):
            q,r=scalars(t)
            return (q/(1-t0*q))**2+(r/(1+t0*r))**2
        return finite_scalar(quad(integrand,0,np.inf,epsabs=2e-12,epsrel=2e-12)[0])
    chi0=chi(0)
    if not 4*m*m-1e-12<=chi0<=m+1e-12:raise ArithmeticError('natural Hessian bound')
    for t0 in (-0.75,-0.5,-0.125,0.125,0.5,0.75):
        At=A-t0*X;rho=np.array(e['rho'].subs(theta,t0),dtype=float)
        from_fock=log_hessian(rho,tangent)
        from_covariance=log_hessian(At,X)+log_hessian(I-At,X)
        from_integral=chi(t0)
        if max(abs(from_fock-from_integral),abs(from_covariance-from_integral))>2e-11:
            raise ArithmeticError('full Fock and one-body Hessians disagree')
        full=relative_entropy(rho,rho0)
        one=finite_scalar(np.trace(At@(log_matrix(At)-log_matrix(A))+(I-At)@(log_matrix(I-At)-log_matrix(I-A))).real)
        integral=finite_scalar(t0*t0*quad(lambda u:(1-u)*chi(u*t0),0,1,epsabs=2e-11)[0])
        if max(abs(full-one),abs(full-integral))>2e-11:raise ArithmeticError('finite relative entropy disagreement')
        b=abs(t0)
        if not 2*m*m*t0*t0-1e-12<=full<=chi0*t0*t0/(2*(1-b)**2)+1e-12:
            raise ArithmeticError('finite amplitude bounds fail')
        if abs(full-chi0*t0*t0/2)>chi0*abs(t0)**3/(3*(1-b)**3)+1e-12:
            raise ArithmeticError('finite amplitude remainder fails')
    # Finite CAR exhaustion is monotone for the actual conditioned states.
    for t0 in (-0.5,0.5):
        rho=np.array(e['rho'].subs(theta,t0),dtype=float)
        rho_small=np.diag([rho[0,0]+rho[2,2],rho[1,1]+rho[3,3]])
        base_small=np.diag([rho0[0,0]+rho0[2,2],rho0[1,1]+rho0[3,3]])
        if not -1e-12<=relative_entropy(rho_small,base_small)<=relative_entropy(rho,rho0)+1e-12:
            raise ArithmeticError('finite CAR restriction monotonicity')
    # Wrong equal signs in particle/hole denominators give a different answer.
    wrong=finite_scalar(quad(lambda t:sum((v/(1-0.5*v))**2 for v in scalars(t)),0,np.inf)[0])
    if abs(wrong-chi(0.5))<1e-4:raise ArithmeticError('sign negative control ineffective')
    raw_norm=9/16;raw=np.sqrt(raw_norm)*f
    if np.linalg.norm(np.outer(raw,raw)/raw_norm-X)>1e-14:raise ArithmeticError('cell normalization')
    if np.linalg.norm(np.outer(raw,raw)/float(e['occupation'])-X)<1e-3:raise ArithmeticError('normalization negative control ineffective')
    # Spectator copies remain present; the response direction has support in
    # one chosen block, without deleting the other seven blocks.
    full_A=np.kron(np.eye(8),A);full_X=np.zeros((16,16));full_X[:2,:2]=X
    if abs(log_hessian(full_A,full_X)+log_hessian(np.eye(16)-full_A,full_X)-chi0)>2e-11:
        raise ArithmeticError('spectator direct-sum control')
    return {'full_Fock_Hessian_and_entropy_checked':True,'noncommuting_amplitudes_checked':6,
            'integral_CAR_and_covariance_forms_agree':True,
            'opposite_resolvent_sign_mutation_rejected':True,
            'occupation_in_place_of_cell_norm_rejected':True,
            'eight_block_spectator_control':True,
            'finite_CAR_restriction_monotonicity_checked':True,
            'absolute_diagnostic_tolerance':'2e-11','floating_diagnostics_used_as_proof':False}


def endpoint_checks():
    for exponent in (2,4,8,12):
        x=10.0**(-exponent);m=x*(1-x)
        hessian=m*m*(1/x+1/(1-x))
        if abs(hessian-m)>1e-14:raise ArithmeticError('natural endpoint control')
    # Fixed norm without the natural factorization does not have this bound.
    x=1e-12;m=0.01
    if not m*m*(1/x+1/(1-x))>1e6:raise ArithmeticError('unnatural endpoint control ineffective')
    return {'endpoint_exponents_checked':[2,4,8,12],
            'natural_Hessian_stays_bounded':True,'arbitrary_fixed_norm_probe_extension_rejected':True}


def certificate():
    return {'finite':finite_checks(),'endpoints':endpoint_checks()}
