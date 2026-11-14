"""Independent full-state modular field attachment audit, without producer builders."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(BASE/'next/C9'));sys.path.insert(0,str(BASE/'candidate'))
from fermionic_shared_records import FermionicSharedRecords
from oph_openwave.modular_native import LocalCoarseModular
from oph_openwave.construction import Construction
from audit_car_native_algebra import phase


def run():
    producer=HERE/'modular_field_adapter.py';receipt=HERE/'MODULAR_FIELD_ADAPTER.json';saved=json.loads(receipt.read_text())
    paths=[Path(__file__),producer,receipt,BASE/'next/C9/fermionic_shared_records.py',BASE/'candidate/oph_openwave/modular_native.py',BASE/'candidate/oph_openwave/construction.py']
    before={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};results=[]
    for car in (False,True):
        c=Construction(3,((0,1),(1,2),(2,0)),mass=.23);m=FermionicSharedRecords(c,7,6) if car else LocalCoarseModular(c,3)
        h=np.zeros((m.n,m.n));w=np.zeros_like(h);f=np.array([x[1] for x in m.states]);reps=np.where(f<=m.modulus//2,f,f-m.modulus)
        # At the even antipode the sign is immaterial to the squared energy.
        electric=np.diag(.37/2*np.sum(reps**2,axis=1))
        for col,(position,word) in enumerate(m.states):
            u=[(position>>i)&1 for i in range(6)] if car else [int(position==i) for i in range(6)]
            h[col,col]=.23*sum((1 if z<3 else -1)*u[z] for z in range(6))
            for e,(a,b) in enumerate(c.seams):
                for side,x in enumerate((a,b)):
                    for src,dst,delta in ((x,3+e,1),(3+e,x,-1)):
                        if u[src] and not u[dst]:
                            target=position^(1<<src)^(1<<dst) if car else dst
                            record=list(word);record[2*e+side]=(record[2*e+side]+delta)%m.modulus
                            h[m.index[(target,tuple(record))],col]+=(1 if side==0 else -1)*(phase(u,src,dst) if car else 1)
            target=tuple((x+d)%m.modulus for x,d in zip(word,(1,-1,1,-1,1,-1)))
            w[m.index[(position,target)],col]=1
        field=electric+.61*(2*np.eye(m.n)-w-w.T);total=h+field
        psi=np.eye(m.n,dtype=complex)[:,0];state=expm(-.7j*total)@psi;rho=np.outer(state,state.conj())
        data=saved['experiments'][int(car)];targetrho=np.array(data['final_state_real'])+1j*np.array(data['final_state_imag'])
        posterror=float(np.max(abs(rho-targetrho)))
        charges=np.array([[(position>>z)&1 if car else int(position==z) for z in range(6)] for position,_ in m.states])
        free=expm(-.7j*h)@psi
        feedback=float(np.max(abs((abs(state)**2-abs(free)**2)@charges)))
        labels=sorted({position for position,_ in m.states});uniform=np.zeros((m.n,len(labels)))
        for q,label in enumerate(labels):
            ids=[i for i,(position,_) in enumerate(m.states) if position==label];uniform[ids,q]=1/np.sqrt(len(ids))
        code=uniform@uniform.T;codebreak=float(np.linalg.norm((np.eye(m.n)-code)@electric@uniform))
        energygain=float(np.vdot(state,field@state).real-np.vdot(psi,field@psi).real)
        conservation=float(np.vdot(state,total@state).real-np.vdot(psi,total@psi).real)
        finite=expm(-1e-5j*total)@state;past=expm(1e-5j*total)@state
        derivative=float((np.vdot(finite,field@finite).real-np.vdot(past,field@past).real)/2e-5)
        commutator=float(np.vdot(state,1j*(h@field-field@h)@state).real)
        lawresults={};ledgererror=0.
        for law in ('C','E','MIX'):
            if car:operations=m.operators[law]
            else:operations=m.operators if law=='MIX' else [tuple((meta,np.sqrt(2)*k) for meta,k in ops if meta['law']==law) for ops in m.operators]
            records=[]
            for ops in operations:
                for meta,k in ops:
                    branch=k@state/np.sqrt(3)
                    records.append({'probability':float(np.vdot(branch,branch).real),'energy':float(np.vdot(branch,total@branch).real)})
            expected=data['all_native_law_work'][law]['complete_outcomes']
            ledgererror=max(ledgererror,max(abs(a['probability']-b['probability']) for a,b in zip(records,expected)),max(abs(a['energy']-b['unnormalized_total_energy']) for a,b in zip(records,expected)))
            work=sum(r['energy'] for r in records)-float(np.vdot(state,total@state).real)
            lawresults[law]={'total_probability':sum(r['probability'] for r in records),'work':work}
        periodic=1j*(electric@w-w@electric);delta=np.array([electric[np.argmax(w[:,i]),np.argmax(w[:,i])]-electric[i,i] for i in range(m.n)])
        faraday=float(np.max(abs(periodic-1j*w@np.diag(delta))))
        checks={'independent_sourceH':bool(np.allclose(h,m.hamiltonian)),'complete_state_reproduced':posterror<1e-12,'fieldgain_reproduced':abs(energygain-data['energy']['field_gain'])<1e-12,
          'field_feedback_on_source':feedback>1e-5,'uniform_compression_rejected':codebreak>1e-3,'closed_total_energy':abs(conservation)<1e-12,'actual_exchange_power':abs(derivative-commutator)<1e-8,'every_native_probability_energy_reproduced':ledgererror<1e-12,'exact_periodic_Faraday':faraday<1e-12}
        if not all(checks.values()):raise ValueError(checks)
        results.append({'domain':'CAR3N6' if car else 'oneLoadN3','checks':checks,'feedback_position_difference':feedback,'uniform_code_electric_leak':codebreak,'state_error':posterror,'ledger_error':ledgererror,'field_gain':energygain,'field_exchange_power':commutator,'power_finite_difference_error':abs(derivative-commutator),'native_work':lawresults})
    after={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if before!=after:raise ValueError('sourcechanged duringreview')
    return {'schema':'oph.c10.modular_field_independent_review.v1','checks':{'source_snapshot_unchanged':True,'both_fullsource_experiments_pass':True},'experiments':results,'source_sha256':before,
       'scope':'Actualfull18/120state reciprocal modularfield extension. Periodic half-record energy is a declared replacement of integerL05halfF²; +2g magnetic scalar included. Fullstate nativework retained. No uniform/character compression, Maxwell3D, physicalalpha/spin/species or outgoingradiation claim.','full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('MODULAR_FIELD_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['experiments']))
