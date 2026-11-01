"""Independent finite-Fock diagnostics for the degree-controlled Dyson bound."""
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from . import algebra,bounds


def test_nonlinear_contact_degree_spectral_bound():
    modes=4;dim=1<<modes
    c=[np.array(algebra.annihilator(j,modes)).astype(complex) for j in range(modes)]
    D=sum((c[j].conj().T@c[j]@c[j+1].conj().T@c[j+1] for j in (0,2)),np.zeros((dim,dim),complex))
    labels=np.diag(D).real
    rng=np.random.default_rng(1827)
    for degree in range(1,7):
        A=np.eye(dim,dtype=complex)
        for _ in range(degree):
            coefficient=rng.normal(size=modes)+1j*rng.normal(size=modes)
            field=sum((coefficient[j]*c[j] for j in range(modes)),np.zeros((dim,dim),complex))
            A=A@(field if rng.integers(2) else field.conj().T)
        active=np.abs(A)>1e-10
        assert np.max(np.abs(labels[:,None]-labels[None,:])[active],initial=0)<=degree
        assert np.linalg.norm(D@A-A@D,2)<=degree*np.linalg.norm(A,2)+1e-10


def test_noncommuting_full_fock_dyson_remainder():
    modes=4;dim=1<<modes
    c=[np.array(algebra.annihilator(j,modes)).astype(complex) for j in range(modes)]
    hopping=(c[0].conj().T@c[2]+c[2].conj().T@c[0]
             +0.7j*(c[1].conj().T@c[3]-c[3].conj().T@c[1]))
    D=sum((c[j].conj().T@c[j]@c[j+1].conj().T@c[j+1] for j in (0,2)),np.zeros((dim,dim),complex))
    U=0.5;T=0.3;N=4;A=(c[0]+1j*c[3])/np.sqrt(2)
    assert np.linalg.norm(hopping@D-D@hopping)>0
    initial=np.zeros((N,dim,dim),complex);initial[0]=A
    def derivative(_,flat):
        terms=flat.reshape(N,dim,dim);out=np.empty_like(terms)
        for n in range(N):
            out[n]=1j*(hopping@terms[n]-terms[n]@hopping)
            if n:out[n]+=1j*U*(D@terms[n-1]-terms[n-1]@D)
        return out.ravel()
    solution=solve_ivp(derivative,(0,T),initial.ravel(),rtol=1e-11,atol=1e-13)
    assert solution.success
    approximation=solution.y[:,-1].reshape(N,dim,dim).sum(axis=0)
    unitary=expm(1j*T*(hopping+U*D));exact=unitary@A@unitary.conj().T
    error=np.linalg.norm(exact-approximation,2)
    assert 1e-9<error<float(bounds.dyson_remainder(1,N,F(1,2),F(3,10)))
    # Omitting the actual contact evolution fails the same high-order envelope.
    free=expm(1j*T*hopping);free_error=np.linalg.norm(exact-free@A@free.conj().T,2)
    assert free_error>10*float(bounds.dyson_remainder(1,N,F(1,2),F(3,10)))


def test_positive_quadratic_memory_can_escape_in_energy():
    c=[np.array(algebra.annihilator(j)).astype(complex) for j in range(2)]
    E0=1;k=F(1,2);T=F(3,2);Z=2
    previous=None
    for L in (8,32,128):
        data=bounds.ultraviolet_channel(E0,k,L,T,Z)
        h=np.array([[E0,float(k)],[float(k),L]],complex)
        H=sum((h[j,l]*c[j].conj().T@c[l] for j in range(2) for l in range(2)),np.zeros((4,4),complex))
        assert np.min(np.linalg.eigvalsh(H))>=-1e-12
        assert np.allclose(H@c[0]-c[0]@H,-E0*c[0]-float(k)*c[1])
        W=expm(1j*float(T)*H)
        actual=np.linalg.norm(W@c[0]@W.conj().T-np.exp(-1j*E0*float(T))*c[0],2)
        assert actual<=float(data['finite_time_field_error_upper'])
        z=1j
        G=np.linalg.inv(z*np.eye(2)-h)[0,0]
        assert abs(1/G-(z-E0-float(k*k)/(z-L)))<1e-12
        assert abs(float(k*k)/(z-L))<=float(data['selfenergy_upper_on_disk'])
        assert data['memory_Gram']==F(1,4)
        if previous is not None:assert actual<previous
        previous=actual
