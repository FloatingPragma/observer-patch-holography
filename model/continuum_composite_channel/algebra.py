"""Exact full-Fock identities for the nonlinear canonical contact channel."""
import sympy as s


def annihilator(mode, modes=2):
    result=s.zeros(1<<modes)
    for mask in range(1<<modes):
        if mask>>mode&1:
            result[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return result


def anti(a,b):
    return a*b+b*a


def certificate():
    c=[annihilator(j) for j in range(2)]
    n=[a.adjoint()*a for a in c]
    b=[n[1-j]*c[j] for j in range(2)]
    d=[2*b[j]-c[j] for j in range(2)]
    identity=s.eye(4);zero=s.zeros(4);contact=n[0]*n[1]
    phase=identity-2*contact
    for j in range(2):
        assert phase*c[j]*phase==-d[j]
        assert contact*c[j]-c[j]*contact==-b[j]
        for k in range(2):
            assert anti(d[j],d[k])==zero
            assert anti(d[j],d[k].adjoint())==(identity if j==k else zero)
            assert anti(b[j],b[k].adjoint())==anti(b[j],c[k].adjoint())
    # Independent symbolic complex spinors: no reality restriction on the smears.
    x=s.symbols('u0r u0i u1r u1i v0r v0i v1r v1i',real=True)
    u=s.Matrix([x[0]+s.I*x[1],x[2]+s.I*x[3]])
    v=s.Matrix([x[4]+s.I*x[5],x[6]+s.I*x[7]])
    def smear(fields,w):return sum((s.conjugate(w[j])*fields[j] for j in range(2)),zero)
    cu,cv,bu,du=smear(c,u),smear(c,v),smear(b,u),smear(d,u)
    pairing=(u.adjoint()*v)[0]
    q=pairing*s.eye(2)-v*u.adjoint()
    epsilon=s.Matrix([[0,1],[-1,0]])
    assert s.simplify(q-(epsilon*s.conjugate(u))*(epsilon*s.conjugate(v)).adjoint())==s.zeros(2)
    Q=sum((q[j,k]*c[j].adjoint()*c[k] for j in range(2) for k in range(2)),zero)
    P=s.conjugate(u[1]*v[0]-u[0]*v[1])*c[0]*c[1]
    assert s.simplify(anti(bu,cv.adjoint())-Q)==zero
    assert s.simplify(anti(bu,cv)-P)==zero
    assert s.simplify(anti(du,cv.adjoint())-(2*Q-pairing*identity))==zero
    assert s.simplify(anti(du,cv)-2*P)==zero
    assert s.simplify(s.trace(q.adjoint()*q)-(u.adjoint()*u)[0]*(v.adjoint()*v)[0])==0
    return {'onsite_Fock_dimension':4,'complex_spinor_real_parameters':8,
            'phase_identity':'(I-2 n_up n_down)c_sigma(I-2 n_up n_down)=-d_sigma',
            'composite':'d_sigma=(2 n_opposite-1)c_sigma',
            'composite_CAR_exact':True,'mixed_normal_and_anomalous_formulas_exact':True,
            'rank_one_quadratic_coefficient_exact':True,
            'contact_equation':'[H,c(u)]=-c((h+U/2)u)-(U/2)d(u)',
            'microscopic_mode_added':False}


def finite_covariance_control():
    """A non-diagonal pure Slater covariance and a full 16-state replay."""
    c=[annihilator(j,4) for j in range(4)]
    hadamard=s.Matrix([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2
    orbitals=hadamard[:,:2];pi=orbitals*orbitals.adjoint()
    vacuum=s.zeros(16,1);vacuum[0]=1
    creators=[sum((orbitals[j,k]*c[j].adjoint() for j in range(4)),s.zeros(16)) for k in range(2)]
    state=creators[0]*creators[1]*vacuum
    assert (state.adjoint()*state)[0]==1
    u=[s.Matrix([1,s.I])/2,s.Matrix([1,-1])/2]
    v=[s.Matrix([1,1])/2,s.Matrix([s.I,1])/2]
    q=s.diag(*[((a.adjoint()*b)[0]*s.eye(2)-b*a.adjoint()) for a,b in zip(u,v)])
    Q=sum((q[j,k]*c[j].adjoint()*c[k] for j in range(4) for k in range(4)),s.zeros(16))
    mean=(state.adjoint()*Q*state)[0]
    variance=s.simplify((state.adjoint()*Q.adjoint()*Q*state)[0]-s.conjugate(mean)*mean)
    contraction=s.simplify(s.trace(pi*q.adjoint()*(s.eye(4)-pi)*q))
    assert variance==contraction
    I4=s.simplify(s.trace(q.adjoint()*q));assert 0<=variance<=I4
    return {'Fock_dimension':16,'offdiagonal_covariance':True,
            'exact_variance':str(variance),'one_particle_trace_variance':str(contraction),
            'Hilbert_Schmidt_upper_bound':str(I4),'exact_full_Fock_replay':True}
