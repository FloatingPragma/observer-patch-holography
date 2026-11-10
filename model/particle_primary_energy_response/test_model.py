"""Independent finite matrix and complete-record response checks."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm, sqrtm

from . import model as m


def density(H):
    e,V=np.linalg.eigh(H);beta=.01;order=4;scale=2
    y=beta*(max(e)-e)/(2*scale)
    from math import factorial
    t=sum(y**j/factorial(j) for j in range(order+1))
    p=t**(2*scale);p/=sum(p)
    slope=beta*(1-y**order/(factorial(order)*t))
    return (V*p)@V.conj().T,p,float(min(slope))


@pytest.mark.parametrize('axis',(0,1,2))
def test_literal_full_gauss_control_current_projection(axis):
    H,A,B,phi,row=m.q6_two_cell_control(axis)
    assert row['dimension']==68
    assert np.allclose(H,H.conj().T) and np.allclose(A,A.conj().T)
    assert np.allclose(B,B.conj().T)
    assert row['one_word_current_norm_squared']==pytest.approx(4.5,abs=1e-12)
    assert np.vdot(phi,phi)==1
    exact=m.current_witness()
    assert exact['axes'][axis]['projected_commutator_norm_squared_at_a_one']==F(9,2)
    assert exact['corner_projected_commutator_norm_squared_coefficient']==F(27,2)


def test_original_susceptibility_and_full_current_measurement():
    H,A,B,phi,row=m.q6_two_cell_control(1);rho,p,slope=density(H)
    q=float(np.trace(rho@(A@(H@A-A@H)-(H@A-A@H)@A)).real)
    assert q>0
    assert q>=min(p)*slope*np.linalg.norm(H@A-A@H,'fro')**2-1e-10
    assert abs(np.trace(rho@B))<1e-12
    C=np.linalg.norm(B,2);tau=1e-5
    U=np.linalg.solve(np.eye(len(H))+1j*tau*(H+A)/2,np.eye(len(H))-1j*tau*(H+A)/2)
    out=U@rho@U.conj().T;mean=float(np.trace(out@B).real)
    assert mean<0 and mean/tau==pytest.approx(-q,rel=1e-5)
    psi=np.asarray(sqrtm(rho),complex);total=np.zeros_like(rho);weights=[]
    for sign in (-1,1):
        E=(np.eye(len(H))+sign*B/C)/2
        root=np.asarray(sqrtm(E),complex);retained=root@U@psi
        weights.append(float(np.vdot(retained,retained).real))
        total+=retained@retained.conj().T
    assert sum(weights)==pytest.approx(1,abs=1e-12)
    assert np.trace(total)==pytest.approx(1,abs=1e-12)
    assert C*(weights[1]-weights[0])==pytest.approx(mean,abs=1e-10)
    # The witness is not the prepared density; replacing it changes the law.
    assert not np.allclose(rho,np.outer(phi,phi.conj()))


def test_relative_evolution_derivative_and_current_remainder():
    H,A,B,phi,row=m.q6_two_cell_control(0);rho,p,slope=density(H)
    T=.02;eps=1e-5
    W=lambda lam:expm(1j*H*T)@expm(-1j*(H+lam*A)*T)
    # Quadrature is independent of the finite-difference perturbation.
    from numpy.polynomial.legendre import leggauss
    nodes,weights=leggauss(40);K=np.zeros_like(H)
    for x,w in zip(nodes,weights):
        t=T*(x+1)/2;U=expm(1j*H*t);K+=w*T/2*(U@A@U.conj().T)
    O=np.diag(np.linspace(-1,1,len(H)))
    fd=(W(eps).conj().T@O@W(eps)-W(-eps).conj().T@O@W(-eps))/(2*eps)
    assert np.linalg.norm(fd-1j*(K@O-O@K),2)<2e-8
    q=float(np.trace(rho@(A@(H@A-A@H)-(H@A-A@H)@A)).real)
    R=np.linalg.norm(H,2);k=np.linalg.norm(A,2);C=2*R*k;G=R+k
    t=1e-4;U=expm(-1j*(H+A)*t);b=float(np.trace(U@rho@U.conj().T@B).real)
    assert abs(b+t*q)<=2*t*t*G*G*C
    Q=np.linalg.solve(np.eye(len(H))+1j*t*(H+A)/2,np.eye(len(H))-1j*t*(H+A)/2)
    assert np.linalg.norm(Q-U,2)<=t**3*G**3/12


@pytest.mark.parametrize('level',(0,1,2,8))
def test_same_source_factored_pulse_and_adversarial_domain(level):
    member=m.contract.member(level);bound=m.bind_member(member);m.verify_bound(member,bound)
    r=bound['result'];N=r['source_density_floor']['negative_exponent']
    q=r['original_double_commutator_lower'];t=r['pulse_duration_source_clock'];d=r['complete_binary_total_variation_lower']
    assert q['negative_exponent']==t['negative_exponent']==N
    assert d['negative_exponent']==2*N
    assert F(r['projected_current_norm_squared_lower'])==F(27,2)/F(member['source']['spacing'])**4
    assert F(d['rational_multiplier'])==F(q['rational_multiplier'])*F(t['rational_multiplier'])/(4*F(r['current_norm_upper']))
    assert not r['pulse_control_is_a_proved_proper_metric']
    assert not r['Lorentz_stress_or_Einstein_dynamics_proved']
    assert r['physical_mass_table_GeV']==[]
    forged=deepcopy(bound);forged['result']['pulse_control_is_a_proved_proper_metric']=True
    forged['result_sha256']=m.contract.digest(forged['result'])
    with pytest.raises(ValueError):m.verify_bound(member,forged)


def test_controlled_energy_balance_retains_applied_work():
    H,A,B,phi,row=m.q6_two_cell_control(2);rho,p,slope=density(H)
    local=(A,H-A);f=(1.3,.7);dotf=(.2,-.1)
    energies=[f[i]*local[i] for i in range(2)];controlled=sum(energies)
    for i in range(2):
        derivative=dotf[i]*local[i]+1j*(controlled@energies[i]-energies[i]@controlled)
        flux=sum(1j*(energies[i]@energies[j]-energies[j]@energies[i]) for j in range(2))
        assert np.allclose(derivative+flux,dotf[i]*local[i],atol=1e-12)
    assert not np.allclose(sum(dotf[i]*local[i] for i in range(2)),0)
