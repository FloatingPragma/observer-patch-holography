"""Independent full-matrix derivative and continuity checks of the mean-field source."""
from pathlib import Path
from itertools import product
import hashlib,json
import numpy as np
import sympy as s
from scipy.linalg import eigh,expm
from model.openwave_larmor_dynamics.dynamics import diagonalize,current

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PAULI=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.])]
BETA=np.diag([1.,1.,-1.,-1.]);GAMMA=[np.kron(np.array([[0,1],[1,0]]),x) for x in PAULI]


def independent_dirac(field,p,external):
    n=p['n'];a=p['a'];N=p['field_n'];coords=list(product(range(n),repeat=3))
    lookup={x:i for i,x in enumerate(coords)};v=n**3;off=(N-n)//2
    out=np.zeros((4*v,4*v),complex)
    for x,i in lookup.items():
        pos=(np.array(x)-(n-1)/2)*a;ix=np.arange(4)*v+i
        out[np.ix_(ix,ix)]+=BETA*(p['mass']+p['trap']*(pos@pos)+3*p['wilson']/a)
        for mu in range(3):
            y=list(x);y[mu]+=1
            if y[mu]>=n:continue
            j=lookup[tuple(y)];jx=np.arange(4)*v+j
            midpoint=pos.copy();midpoint[mu]+=a/2
            avec=np.cross(np.array([0.,0.,external]),midpoint)/2
            phase=np.exp(1j*a*(field[tuple(np.array(x)+off)+(mu,)]+avec[mu]))
            block=(1j*GAMMA[mu]-p['wilson']*BETA)*phase/(2*a)
            out[np.ix_(ix,jx)]+=block;out[np.ix_(jx,ix)]+=block.conj().T
    return out


def abs_matrix(d):
    w,u=eigh(d);return (u*abs(w))@u.conj().T


def run():
    p=dict(n=3,a=1.,field_n=9,mass=8.,trap=.01,wilson=.5,kappa=.1)
    rng=np.random.default_rng(71234);field=rng.normal(scale=.02,size=(9,9,9,3));B=.03
    matrix=independent_dirac(field,p,B);w,u,ph=diagonalize(field,p,B)
    error=np.linalg.norm(matrix-(u*w)@u.conj().T,2);assert error<2e-12
    kernel=abs_matrix(matrix);max_derivative=0.;max_continuity=0.
    for _ in range(3):
        psi=rng.normal(size=108)+1j*rng.normal(size=108);psi/=np.linalg.norm(psi)
        J=current(psi,w,u,ph,p)
        for index in [(3,3,3,0),(4,3,4,1),(4,4,3,2)]:
            eps=2e-5;plus=field.copy();minus=field.copy();plus[index]+=eps;minus[index]-=eps
            derivative=np.vdot(psi,(abs_matrix(independent_dirac(plus,p,B))-abs_matrix(independent_dirac(minus,p,B)))@psi).real/(2*eps)
            max_derivative=max(max_derivative,abs(derivative+J[index]))
        dt=1e-4
        f=expm(-1j*dt*kernel)@psi;b=expm(1j*dt*kernel)@psi
        rate=((abs(f)**2-abs(b)**2)/(2*dt)).reshape(4,3,3,3).sum(axis=0)
        dadj=sum(np.roll(J[...,mu],1,axis=mu)-J[...,mu] for mu in range(3))
        max_continuity=max(max_continuity,float(np.linalg.norm(rate+dadj[3:6,3:6,3:6])))
    assert max_derivative<3e-9 and max_continuity<2e-8
    M=s.Rational(1403,100);m=s.Integer(8);r=s.Rational(1,2)
    C=((1+(M/m)**2)*2*(1+r*r)+M*(1+r)/2)/m
    lower=4*(s.Rational(1,3)-s.Rational(22,63)**3/6)**2
    q=s.Rational(1,10)*C/lower
    assert C==s.Rational(5302169,2048000) and q<s.Rational(608135,1000000)
    files=[HERE/'independent_review.py',*(HERE.parent/p for p in ['PROOF.md','dynamics.py','experiment.py','results.json'])]
    return {'schema':'oph.larmor_independent_review.v1','full_signed_matrix_dimension':108,
        'independent_matrix_error':float(error),'random_states_with_both_D_signs':3,
        'current_energy_derivative_comparisons':9,'maximum_current_derivative_error':float(max_derivative),
        'maximum_independent_density_continuity_error':max_continuity,
        'exact_fixed_state_contraction_upper':str(q),
        'scope':'The nonlinear quasistatic law is a semiclassical state-dependent construction, not an arbitrary-reference CP channel; no long-time interval error or continuum g is inferred.',
        'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'verdict':'PASS'}


if __name__=='__main__':
    value=run();(HERE/'independent_review.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(value)
