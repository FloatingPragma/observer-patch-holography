"""Free-cap Heisenberg and normalized spectator-reduction controls."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
sys.dont_write_bytecode=True
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'common_dynamics/photon_size3'))
import providers

def run():
    pins=providers.verify()
    a=[]
    for j in range(2):
        op=s.zeros(3);op[0,j+1]=1;a.append(op)
    Omega=s.Matrix([[3,s.Rational(1,5)],[s.Rational(1,5),4]])
    H=s.eye(3)*s.trace(Omega)/2
    H[1:3,1:3]+=Omega
    for j in range(2):
        assert H*a[j]-a[j]*H==-sum((Omega[j,l]*a[l] for l in range(2)),s.zeros(3))
    # Diagonal source fields: verify exact projected first-order equations.
    freqs=[s.Integer(3),s.Integer(4)];electric_a=s.Rational(1,2)
    Hdiag=s.diag(s.Rational(7,2),s.Rational(13,2),s.Rational(15,2))
    naive=s.zeros(3)
    for j,w in enumerate(freqs):
        theta=s.sqrt(electric_a/(2*w))*(a[j]+a[j].T)
        electric=s.I*s.sqrt(w/(2*electric_a))*(a[j].T-a[j])
        assert (s.I*(Hdiag*theta-theta*Hdiag)-electric_a*electric)==s.zeros(3)
        assert (s.I*(Hdiag*electric-electric*Hdiag)+w*w/electric_a*theta)==s.zeros(3)
        naive+=w*(a[j]*a[j].T+a[j].T*a[j])/2
    naive_defect=naive-Hdiag
    assert naive_defect!=s.zeros(3)
    ccr_defect=a[0]*a[0].T-a[0].T*a[0]-s.eye(3)
    assert max(abs(x) for x in ccr_defect.diagonal())==2
    # The oriented limit's two Fourier waves satisfy both Maxwell equations.
    x,t,p=s.symbols('x t p',real=True,positive=True)
    phi=p*(x-t)
    directions=[s.Matrix([0,1,0]),s.Matrix([0,0,1])]
    ex=s.Matrix([1,0,0]);zero=s.zeros(3,1)
    for e in directions:
        A=e*s.cos(phi);E=-e*p*s.sin(phi);B=-ex.cross(e)*p*s.sin(phi)
        curl=lambda f:ex.cross(f.diff(x))
        assert s.simplify(A.diff(t)+E)==zero
        assert s.simplify(B.diff(t)+curl(E))==zero
        assert s.simplify(E.diff(t)-curl(B))==zero
        assert E[0]==B[0]==0
    # Correlated spectators and an untouched external reference.
    rng=np.random.default_rng(29021)
    state=rng.normal(size=(4,3,2))+1j*rng.normal(size=(4,3,2))
    state/=np.linalg.norm(state)
    q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    u=np.diag(np.exp(-1j*np.array([1.,2.,2.,0.])/4))
    local=np.array([[1,.1j,0,0],[0,.7,0,.1],[0,0,.8,0],[0,0,.1j,.6]],complex)/2
    reduced=lambda v:np.einsum('asr,bsq->arbq',v,v.conj()).reshape(8,8)
    worst=0.
    for weight in (1.,.37):
        v=np.sqrt(weight)*state
        source=np.einsum('ab,st,btr->asr',u,q,v)
        observed=np.einsum('ab,bsr->asr',u,v)
        worst=max(worst,float(np.linalg.norm(reduced(source)-reduced(observed))))
        assert abs(np.trace(reduced(source)).real-weight)<1e-13
        branch=np.einsum('ab,bsr->asr',local,source)
        local_branch=np.einsum('ab,bsr->asr',local,observed)
        worst=max(worst,float(np.linalg.norm(reduced(branch)-reduced(local_branch))))
        assert np.trace(reduced(branch)).real<=weight
    assert worst<2e-14
    # Replacing a correlated state with marginal products destroys its joint law.
    correlated=np.zeros((4,3,2),complex)
    correlated[0,0,0]=correlated[1,1,1]=1/np.sqrt(2)
    rr=reduced(correlated).reshape(4,2,4,2)
    marg_o=np.einsum('arbr->ab',rr);marg_r=np.einsum('aras->rs',rr)
    difference=float(np.linalg.norm(rr.reshape(8,8)-np.kron(marg_o,marg_r)))
    assert difference>.49
    assert providers.verify()==pins
    out={'schema':'oph.row29.free_source_and_spectator_contract.v1',
      'checks':{'general_two_mode_compressed_Heisenberg_equations':True,
        'exact_projected_E_and_A_dynamics':True,'naive_truncated_quadratic_rejected':True,
        'cap_CCR_defect_retained':True,'oriented_Faraday_and_Ampere_both_polarizations':True,
        'normalized_and_subnormalized_spectator_reduction':True,
        'correlated_reference_and_local_instrument_branch':True,'false_joint_product_rejected':True,
        'source_pins_unchanged':True},
      'naive_truncated_H_defect_diagonal':[str(x) for x in naive_defect.diagonal()],
      'cap_CCR_defect_norm':'2','spectator_trace_max_residual':worst,
      'false_joint_product_Frobenius_difference':difference,
      'observable_scope':'Selected radiative Fourier component of the physical field; full field first moments only when spectator first moments vanish.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'SOURCE_CONTRACT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':out['checks'],'spectator_residual':worst}),flush=True)
if __name__=='__main__':run()
