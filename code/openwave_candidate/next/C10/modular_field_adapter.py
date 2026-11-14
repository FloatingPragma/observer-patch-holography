"""Full-state periodic-half-record field attachment to actual native source."""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
from scipy.linalg import expm
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'candidate'));sys.path.insert(0,str(BASE/'next/C9'))
from oph_openwave.modular_native import LocalCoarseModular
from oph_openwave.construction import Construction
from fermionic_shared_records import FermionicSharedRecords,transport
from fermionic_uniform_reference import reduction


def centered(value,modulus):return (value+modulus//2)%modulus-modulus//2


class ModularField:
    def __init__(self,model,kappa=.37,face_coefficient=.61):
        if not math.isfinite(kappa) or not math.isfinite(face_coefficient) or kappa<0 or face_coefficient<0:raise ValueError('finite nonnegative coefficients required')
        self.model=model;self.n=model.n;self.kappa=kappa;self.g=face_coefficient
        self.car=isinstance(model,FermionicSharedRecords);c=model.construction
        if c.seams!=((0,1),(1,2),(2,0)):raise ValueError('declared triangle face required')
        self.delta=np.array([1,-1]*3);self.representatives=np.array([[centered(x,model.modulus) for x in f] for _,f in model.states])
        self.electric_values=.5*kappa*(self.representatives**2).sum(axis=1)
        self.electric=np.diag(self.electric_values);self.W=np.zeros((self.n,self.n))
        for col,(position,f) in enumerate(model.states):
            record=tuple((np.array(f)+self.delta)%model.modulus)
            self.W[model.index[(position,record)],col]=1
        self.magnetic=face_coefficient*(2*np.eye(self.n)-self.W-self.W.T)
        self.field=self.electric+self.magnetic;self.source=model.hamiltonian.copy();self.hamiltonian=self.source+self.field
        self.operators=model.operators['MIX'] if self.car else model.operators
        self.operators_by_law=model.operators if self.car else {'MIX':self.operators,**{law:[tuple((meta,np.sqrt(2)*k) for meta,k in ops if meta['law']==law) for ops in self.operators] for law in ('C','E')}}
        # Existing source half-hop convention, including CAR transport phases.
        self.halfhops=[np.zeros((self.n,self.n)) for _ in range(6)]
        for col,(position,f) in enumerate(model.states):
            for e,ends in enumerate(c.seams):
                for side,port in enumerate(ends):
                    if self.car:
                        move=transport(position,port,3+e,6)
                        if move is None:continue
                        target,phase=move
                    else:
                        if position!=port:continue
                        target,phase=3+e,1
                    record=list(f);record[2*e+side]=(record[2*e+side]+1)%model.modulus
                    self.halfhops[2*e+side][model.index[(target,tuple(record))],col]=(1 if side==0 else -1)*phase

    def commit(self,rho,law='MIX'):
        out=np.zeros_like(rho);ledger=[]
        for ops in self.operators_by_law[law]:
            for meta,k in ops:
                branch=k@rho@k.T/3;out+=branch
                ledger.append({**meta,'probability':float(np.trace(branch).real),'unnormalized_total_energy':float(np.trace(self.hamiltonian@branch).real),'unnormalized_field_energy':float(np.trace(self.field@branch).real)})
        return out,ledger


def experiment(car=False):
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    model=FermionicSharedRecords(c,7,6) if car else LocalCoarseModular(c,3)
    field=ModularField(model);n=model.n;comm=lambda a,b:a@b-b@a
    omega=np.exp(2j*np.pi/model.modulus);weyl=0.;integer_error=0.
    for h,t in enumerate(field.halfhops):
        z=np.diag(np.exp(2j*np.pi*field.representatives[:,h]/model.modulus));f=np.diag(field.representatives[:,h])
        weyl=max(weyl,float(np.max(abs(z@t-omega*t@z))))
        integer_error=max(integer_error,float(np.max(abs(comm(f,t)-t))))
    # The periodic electric commutator uses the ACTUAL finite energy difference.
    difference=np.array([field.electric_values[np.argmax(field.W[:,col])]-field.electric_values[col] for col in range(n)])
    faraday=float(np.max(abs(1j*comm(field.electric,field.W)-1j*field.W@np.diag(difference))))
    wrong=1j*field.kappa*field.W@np.diag(field.representatives@field.delta+field.delta@field.delta/2)
    integer_faraday_error=float(np.max(abs(1j*comm(field.electric,field.W)-wrong)))
    reciprocal=float(np.max(abs(1j*comm(field.hamiltonian,field.field)+1j*comm(field.hamiltonian,field.source))))
    if car:_,v,_,_,_=reduction(model)
    else:
        v=np.zeros((n,6));counts=np.bincount([s[0] for s in model.states])
        for row,(pos,_) in enumerate(model.states):v[row,pos]=1/np.sqrt(counts[pos])
    p=v@v.T;leak=float(np.linalg.norm((np.eye(n)-p)@field.electric@v))
    charcomm=float(np.linalg.norm(comm(field.hamiltonian,p)))
    # Definite records are an explicit preparation distinct from the uniform code.
    psi=np.zeros(n,complex);psi[0]=1;rho=np.outer(psi,psi.conj())
    u=expm(-.7j*field.hamiltonian);final=u@rho@u.conj().T
    free=expm(-.7j*field.source)@psi
    charge=np.array([[((s[0]>>z)&1) if car else int(s[0]==z) for z in range(6)] for s in model.states])
    positions=np.diag(final).real@charge;control=(abs(free)**2)@charge
    held=expm(-.7j*field.field)@psi
    held_field_gain=float(np.vdot(held,field.field@held).real-np.trace(field.field@rho).real)
    gauss_function=model.gauss if car else model.gauss_residue
    gauss_preserved=all(gauss_function(model.states[np.argmax(field.W[:,q])])==gauss_function(state) for q,state in enumerate(model.states))
    field_charge_comm=max(float(np.max(abs(comm(field.field,np.diag(charge[:,z]))))) for z in range(6))
    es0=float(np.trace(field.source@rho).real);ef0=float(np.trace(field.field@rho).real)
    es1=float(np.trace(field.source@final).real);ef1=float(np.trace(field.field@final).real)
    monitored,ledger=field.commit(final);after=float(np.trace(field.hamiltonian@monitored).real);before=es1+ef1
    all_laws={}
    for law in ('C','E','MIX'):
        post,records=field.commit(final,law)
        all_laws[law]={'complete_outcomes':records,'trace':float(np.trace(post).real),'apparatus_work':float(np.trace(field.hamiltonian@(post-final)).real),'after_energy':float(np.trace(field.hamiltonian@post).real)}
    return {'domain':'fixedN3 CAR' if car else 'one load','modulus':model.modulus,'dimension':n,'checks':{'face_preserves_modular_Gauss':bool(gauss_preserved),'field_no_spurious_load_charge_current':bool(field_charge_comm<1e-12),'held_source_field_energy_control':bool(abs(held_field_gain)<1e-12),'Hermitian_full_H':bool(np.allclose(field.hamiltonian,field.hamiltonian.T)),'unitary_cyclic_face':bool(np.array_equal(field.W.T@field.W,np.eye(n))),'cyclic_Weyl_identity':bool(weyl<1e-12),'exact_periodic_Faraday':bool(faraday<1e-12),'integer_wrap_identity_rejected':bool(integer_error>1),'integer_Faraday_formula_rejected':bool(integer_faraday_error>1e-3),'reciprocal_energy_operator':bool(reciprocal<1e-12),'uniform_character_compression_rejected':bool(leak>1e-3 and charcomm>1e-3),'closed_unitary_total_energy':bool(abs(es1+ef1-es0-ef0)<1e-12),'actual_field_energy_gain':bool(ef1-ef0>1e-5),'reciprocal_load_feedback':bool(np.max(abs(positions-control))>1e-5),'complete_native_instrument':bool(all(abs(sum(r['probability'] for r in law['complete_outcomes'])-1)<1e-12 for law in all_laws.values())),'native_total_energy_accounting':bool(abs(sum(r['unnormalized_total_energy'] for r in ledger)-after)<1e-12)},'errors':{'Weyl':weyl,'periodic_Faraday':faraday,'wrong_integer_commutator':integer_error,'wrong_integer_Faraday':integer_faraday_error,'reciprocal_energy_operator':reciprocal,'uniform_code_electric_leak':leak,'fullH_character_projector_commutator':charcomm},'energy':{'source_before':es0,'source_after':es1,'field_before':ef0,'field_after':ef1,'field_gain':ef1-ef0,'total_unitary_error':es1+ef1-es0-ef0,'native_before':before,'native_after':after,'native_apparatus_work':after-before,'held_source_field_energy_gain':held_field_gain},'full_position_occupations':positions.tolist(),'uncoupled_position_occupations':control.tolist(),'complete_native_outcomes':ledger,'all_native_law_work':all_laws,'final_state_real':final.real.tolist(),'final_state_imag':final.imag.tolist()}


def run():
    experiments=[experiment(False),experiment(True)]
    paths=[Path(__file__),BASE/'next/C9/fermionic_shared_records.py',BASE/'candidate/oph_openwave/modular_native.py',BASE/'candidate/oph_openwave/construction.py']
    return {'schema':'oph.c10.modular_field_adapter.v1','parameters':{'mass':.23,'kappa':.37,'face_coefficient':.61,'time':.7,'electric_choice':'periodic_half_record_energy = .5*kappa*sum_h rN(Fh)^2; rN=(F+floor(N/2))modN−floor(N/2)','magnetic_choice':'g*(2I−W−Wdagger); +2g scalar relative to FableL05 E3','face_translation':[1,-1,1,-1,1,-1],'initial_state':'definite zero records, loadatport0 or CARports0,1,2 occupied','native_measurement':'same complete localcoarse C/E/MIX, separate final instruments'},'experiments':experiments,'scope':'Full finite modular sharedrecord source with reciprocal periodicrecord fieldattachment. Borrowed integer E2b convention explicitly replaced by periodicrepresentative potential; no hardwall/cyclic identity transfer or uniformcode compression. Singleface energyexchange is not distantflux, Maxwell radiation, physicalbackreaction calibration or spacetime.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
if __name__=='__main__':
 r=run();Path(__file__).with_name('MODULAR_FIELD_ADAPTER.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps([{k:v for k,v in x.items() if k in ['domain','dimension','checks','energy']} for x in r['experiments']]))
