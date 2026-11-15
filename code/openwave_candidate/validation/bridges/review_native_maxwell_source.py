"""Independent Codex-worker endpoint/line-integral/augmented-exponential audit."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
import native_maxwell_source as producer

def run():
    c=producer.Construction(3,((0,1),(1,2),(2,0)),.23);model=producer.UniformFiberReduction(c,6)
    pos=np.array([[0,0,0],[1,0,0],[.5,np.sqrt(3)/2,0]],float);k=np.array([.7,.4,.3]);h=model.hamiltonian
    psi=np.array([1,2j,0,1-1j,0,0],complex)/np.sqrt(7)
    qerr=0.;jerr=0.;continuity=0.
    def charge(state):
        p=abs(state)**2;q=p[:3].copy()
        for e,(a,b) in enumerate(c.seams):q[a]+=p[3+e]/2;q[b]+=p[3+e]/2
        return q
    K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    A=np.block([[np.zeros((3,3)),1j*K],[-1j*K,np.zeros((3,3))]])
    transverse=np.eye(3)-np.outer(k,k)/(k@k);field=np.zeros(6,complex);work=0.;dt=.01
    nodes,weights=np.polynomial.legendre.leggauss(8)
    for n in range(40):
        state=expm(-1j*h*n*dt)@psi;s=producer.source(model,state,pos,k);rho=np.outer(state,state.conj());currents=[]
        for e,(a,b) in enumerate(c.seams):
            block=np.zeros_like(h);t=3+e
            for port in (a,b):block[port,t]=h[port,t];block[t,port]=h[t,port]
            velocity=np.diag(-1j*(block@rho-rho@block)).real
            currents.append((velocity[b]-velocity[a])/2)
        line=np.zeros(3,complex)
        for current,(a,b) in zip(currents,c.seams):
            dx=pos[b]-pos[a]
            integral=quad(lambda u:np.cos(k@(pos[a]+u*dx)),0,1,epsabs=1e-13)[0]-1j*quad(lambda u:np.sin(k@(pos[a]+u*dx)),0,1,epsabs=1e-13)[0]
            line+=current*dx*integral
        jerr=max(jerr,float(np.max(abs(line-s['current_fourier']))))
        eps=1e-5;fd=(charge(expm(-1j*h*eps)@state)-charge(expm(1j*h*eps)@state))/(2*eps)
        qerr=max(qerr,float(np.max(abs(fd-s['effective_charge_derivative']))))
        continuity=max(continuity,float(abs(fd@np.exp(-1j*pos@k)+1j*k@line)))
        force=np.r_[-transverse@line,np.zeros(3)];aug=np.zeros((7,7),complex);aug[:6,:6]=A;aug[:6,6]=force
        initial=np.r_[field,1]
        for node,weight in zip(nodes,weights):
            at=(node+1)*dt/2;value=(expm(aug*at)@initial)[:6]
            work+=weight*dt/2*np.vdot(value,force).real
        field=(expm(aug*dt)@initial)[:6]
    actual=producer.execute(40);reported=np.array(actual['field_real'])+1j*np.array(actual['field_imag'])
    ferr=float(np.max(abs(field-reported)));werr=float(abs(work-np.vdot(field,field).real/2))
    if max(jerr,ferr,werr)>1e-11 or max(qerr,continuity)>1e-8:raise ValueError('independent Maxwell bridge mismatch')
    paths=[Path(__file__),Path(producer.__file__),producer.PROVIDER]
    return {'schema':'oph.native.maxwell.independent_review.v1','reviewer':'Independent Codex worker; not Fable opposite-team review','outcome':'PASS_WITH_ONE_WAY_SUPPLIED_FIELD_SCOPE','endpoint_charge_central_difference_error':qerr,'independent_line_quadrature_current_error':jerr,'finite_difference_Fourier_continuity_error':continuity,'augmented_exponential_field_error':ferr,'independent_Gauss_Legendre_work_energy_error':werr,'scope':'Actual coherent native expectation current, half-transit effective charge and supplied straight-line embedding feed existing linear continuum mode. No joint quantum source/field, reciprocal backreaction, finite-volume radiation or derived field dispersion/clock. Held source error and anisotropy/embedding are declared assumptions.','full_criterion_met':False,'source_sha256':{str(p.relative_to(producer.BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
if __name__=='__main__':
    r=run();Path(__file__).with_name('NATIVE_MAXWELL_SOURCE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
