"""Independent source reconstruction and complete outcome/reference checks."""
import importlib.util,sys,json,hashlib
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[2];HERE=Path(__file__).parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.construction import Construction
from audit_car_native_algebra import phase,run as algebra_audit
spec=importlib.util.spec_from_file_location('review_target_car',HERE/'fermionic_shared_records.py');target=importlib.util.module_from_spec(spec);spec.loader.exec_module(target)


def run():
    algebra=algebra_audit();c=Construction(3,((0,1),(1,2),(2,0)),mass=.23);m=target.FermionicSharedRecords(c,7,3)
    h=np.zeros((m.n,m.n));ops={law:[{} for _ in c.seams] for law in ('C','E','MIX')}
    for col,(mask,f) in enumerate(m.states):
        u=[(mask>>i)&1 for i in range(6)]
        h[col,col]=.23*sum((1 if i<3 else -1)*u[i] for i in range(6))
        for e,(a,b) in enumerate(c.seams):
            for side,x in enumerate((a,b)):
                for src,dst,delta in ((x,3+e,1),(3+e,x,-1)):
                    if u[src] and not u[dst]:
                        word=list(f);word[2*e+side]=(word[2*e+side]+delta)%3
                        row=m.index[(mask^(1<<src)^(1<<dst),tuple(word))]
                        h[row,col]+=(1 if side==0 else -1)*phase(u,src,dst)
            d=u[a]-u[b]
            for law in ('C','E','MIX'):
                for selected in (('C','E') if law=='MIX' else (law,)):
                    for coin in (0,1):
                        j=(d-(2*coin-1)*abs(d))//2
                        newmask=mask if not j else mask^(1<<a)^(1<<b)
                        word=list(f);word[2*e]=(word[2*e]+j)%3;word[2*e+1]=(word[2*e+1]-j)%3
                        key=(selected,coin,j if selected=='C' else d)
                        k=ops[law][e].setdefault(key,np.zeros((m.n,m.n)))
                        k[m.index[(newmask,tuple(word))],col]=(phase(u,a,b) if j else 1)/(2 if law=='MIX' else np.sqrt(2))
    error=0.
    for law in ops:
        for e in range(3):
            for meta,k in m.operators[law][e]:error=max(error,float(np.max(abs(k-ops[law][e][(meta['law'],meta['coin'],meta['record_value'])]))))
    # Arbitrary coherent state: reconstruct each event's full poststate and energy.
    psi=np.arange(1,m.n+1)+1j*np.arange(m.n,0,-1);psi/=np.linalg.norm(psi);rho=np.outer(psi,psi.conj())
    total=np.zeros_like(rho);prob=[];energies=[]
    for e in range(3):
        for key,k in sorted(ops['MIX'][e].items()):
            vector=k@psi/np.sqrt(3);branch=np.outer(vector,vector.conj());total+=branch
            prob.append(float(np.trace(branch).real));energies.append(float(np.trace(h@branch).real))
    unsigned=sum((abs(k)@rho@abs(k).T/3 for e in range(3) for k in ops['MIX'][e].values()),start=np.zeros_like(rho))
    phase_control=float(np.linalg.norm(unsigned-total,ord='nuc')/2)
    output,ledger=m.commit(rho);out_error=float(np.max(abs(output-total)))
    # Complete local operation on an external reference: Kraus vectors retain
    # original off-diagonal reference blocks, no implicit measurement of it.
    entangled=np.zeros((2,m.n),complex);entangled[0,0]=1/np.sqrt(2);entangled[1,-1]=1j/np.sqrt(2)
    cq=[]
    for e in range(3):
        for k in ops['MIX'][e].values():cq.append((entangled@k.T/np.sqrt(3)).reshape(-1))
    joint=sum((np.outer(x,x.conj()) for x in cq),start=np.zeros((2*m.n,2*m.n),complex))
    ref=np.trace(joint.reshape(2,m.n,2,m.n),axis1=1,axis2=3)
    checks={'independent_actual_H':float(np.max(abs(h-m.hamiltonian)))<1e-12,'every_C_E_MIX_K':error<1e-12,
       'complete_poststates':out_error<1e-12,'complete_event_probabilities':max(abs(np.array(prob)-[r['probability'] for r in ledger['outcomes']]))<1e-12,
       'complete_event_energy':max(abs(np.array(energies)-[r['unnormalized_final_energy'] for r in ledger['outcomes']]))<1e-12,
       'external_reference_retained':bool(np.allclose(ref,np.eye(2)/2)),
       'external_reference_joint_positive':float(np.linalg.eigvalsh(joint).min())>-1e-12,
       'dropping_CAR_phases_changes_output':phase_control>1e-3}
    checks={k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    files=[Path(__file__),HERE/'fermionic_shared_records.py',HERE/'FERMIONIC_SHARED_RECORDS.json',HERE/'audit_car_native_algebra.py']
    return {'schema':'oph.c9.fermionic_shared_independent_review.v1','checks':checks,'algebra':algebra,
       'dimension':m.n,'phase_omission_trace_distance':phase_control,'Hamiltonian_error':float(np.max(abs(h-m.hamiltonian))),'every_Kraus_error':error,'complete_output_error':out_error,
       'scope':'Exact finite fixed-number CAR source/instrument attachment. Local discreteZ3 covariance; continuous integer-representative gauge rejected. Positive physical sea/species and monitoring stationarity are not established.',
       'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'full_criterion_met':False}

if __name__=='__main__':
    r=run();HERE.joinpath('FERMIONIC_SHARED_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
