"""Complete charged reservoir exchange on one shared-record CAR source."""
from pathlib import Path
import sys,json,hashlib,math,importlib.util
import numpy as np
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'next/C9'));sys.path.insert(0,str(BASE/'candidate'))
from fermionic_shared_records import FermionicSharedRecords,transport,FOCK
from oph_openwave.construction import Construction
from oph_openwave.fiber_preparation import FiberPreparation
_spec=importlib.util.spec_from_file_location('c10_portable_fock_provider',FOCK)
_fock=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_fock)
def annihilator(index,modes=4):return _fock.operator(modes,index,False)


class BoundarySource(FermionicSharedRecords):
    def __init__(self,modulus=3):
        c=Construction(3,((0,1),(1,2),(2,0)),.23)
        self.construction=c;self.positions=6;self.modulus=modulus;self.number=None
        prep=FiberPreparation(c,modulus);states=[]
        for mask in range(64):
            b=3-mask.bit_count();g=[int(z<3)-((mask>>z)&1)-int(z==0)*b for z in range(6)]
            for cycle in range(modulus):
                f=[0]*6;res=g.copy()
                for edge,value in zip(prep.chords,(cycle,)):
                    f[edge]=value;a,t=prep.edges[edge];res[a]-=value;res[t]+=value
                for z in reversed(prep.order[1:]):
                    parent,edge=prep.parent[z];a,t=prep.edges[edge]
                    f[edge]=(res[z] if a==z else -res[z])%modulus
                    res[parent]+=res[z];res[z]=0
                if res[0]%modulus:raise ValueError('inconsistent dressed boundary Gauss sector')
                states.append((mask,tuple(f)))
        self.states=tuple(states);self.index={s:i for i,s in enumerate(states)};self.n=len(states)
        self.hamiltonian=np.zeros((self.n,self.n))
        for col,(mask,f) in enumerate(states):
            self.hamiltonian[col,col]=c.mass*sum((1 if z<3 else -1)*((mask>>z)&1) for z in range(6))
            for e,ends in enumerate(c.seams):
                for side,port in enumerate(ends):
                    for src,dst,delta in ((port,3+e,1),(3+e,port,-1)):
                        moved=transport(mask,src,dst,6)
                        if moved is None:continue
                        record=list(f);record[2*e+side]=(record[2*e+side]+delta)%modulus
                        self.hamiltonian[self.index[(moved[0],tuple(record))],col]+=(1 if side==0 else -1)*moved[1]
        self.operators={law:[self.kraus(e,law) for e in range(3)] for law in ('C','E','MIX')}

    def gauss(self,state):
        internal=list(super().gauss(state));internal[0]=(internal[0]+3-state[0].bit_count())%self.modulus
        return tuple(internal)

    def reservoir_instrument(self):
        # Reuse exact existing CAR annihilators and the mass-instrument swap method.
        a=np.asarray(annihilator(0,7)).astype(float);r=np.asarray(annihilator(6,7)).astype(float)
        swap=np.eye(128)-a.T@a-r.T@r+a.T@r+r.T@a
        if not np.array_equal(swap.T@swap,np.eye(128)):raise ValueError('complete CAR swap unitarity failed')
        result=[]
        for initial in (0,1):
            for observed in (0,1):
                k=np.zeros((self.n,self.n))
                for col,(mask,f) in enumerate(self.states):
                    for target_mask in np.flatnonzero(swap[64*observed:64*(observed+1),mask+64*initial]):
                        k[self.index[(int(target_mask),f)],col]=swap[target_mask+64*observed,mask+64*initial]/math.sqrt(2)
                result.append(({'reservoir_initial':initial,'reservoir_final':observed,'source_charge_change':initial-observed,'boundary_record_change':observed-initial},k))
        return tuple(result)


def run():
    source=BoundarySource(3);h=source.hamiltonian;n=np.array([mask.bit_count() for mask,_ in source.states]);ids=np.flatnonzero(n==3)
    e,v=np.linalg.eigh(h[np.ix_(ids,ids)]);selected=abs(e-e[0])<1e-10
    rho=np.zeros((source.n,source.n));rho[np.ix_(ids,ids)]=v[:,selected]@v[:,selected].T/sum(selected)
    reservoir_energy=.23;ops=source.reservoir_instrument();records=[];out=np.zeros_like(rho);gauss=True;boundary=True
    # Complete final energy eigenspaces are resolved inside each actual number sector.
    energy_bases={q:np.linalg.eigh(h[np.ix_(np.flatnonzero(n==q),np.flatnonzero(n==q))]) for q in range(7)}
    for meta,k in ops:
        branch=k@rho@k.T;out+=branch;prob=float(np.trace(branch));energy=float(np.trace(h@branch))
        energyrecords=[]
        for q,(values,vectors) in energy_bases.items():
            qids=np.flatnonzero(n==q);probs=np.diag(vectors.T@branch[np.ix_(qids,qids)]@vectors).real
            groups=[]
            for value,p in zip(values,probs):
                if groups and abs(value-groups[-1]['energy'])<1e-9:groups[-1]['probability']+=float(p);groups[-1]['degeneracy']+=1
                else:groups.append({'energy':float(value),'probability':float(p),'degeneracy':1})
            energyrecords.extend({'source_number':q,**g} for g in groups)
        for row,col in zip(*np.nonzero(k)):
            initial_state=source.states[col];final_state=source.states[row]
            gauss &= source.gauss(initial_state)==source.gauss(final_state)==(1,1,1,0,0,0)
            bi=3-initial_state[0].bit_count();bf=3-final_state[0].bit_count()
            boundary &= meta['reservoir_initial']-bi==meta['reservoir_final']-bf
        work=energy+reservoir_energy*meta['reservoir_final']*prob-(float(e[0])+reservoir_energy*meta['reservoir_initial'])*prob
        records.append({**meta,'initial_source_number':3,'initial_source_energy':float(e[0]),'initial_source_reference':'complete N3 ground multiplet','reservoir_preparation_probability':.5,'status':'addition' if meta['source_charge_change']==1 else 'removal' if meta['source_charge_change']==-1 else 'no_change','probability':prob,'unnormalized_source_energy':energy,'unnormalized_reservoir_energy':reservoir_energy*meta['reservoir_final']*prob,'complete_number_energy_eigenspace_records':energyrecords,'unnormalized_apparatus_work':work,'poststate_real':branch.tolist()})
    completeness=float(np.max(abs(sum(k.T@k for _,k in ops)-np.eye(source.n))))
    files=[Path(__file__),BASE/'next/C9/fermionic_shared_records.py',BASE/'candidate/oph_openwave/fiber_preparation.py',FOCK,BASE/'validation/bridges/providers/mass_instrument_model.py']
    return {'schema':'oph.c10.charged_reservoir.v1','checks':{'complete_reservoir_instrument':bool(completeness<1e-12),'source_and_boundary_Gauss_preserved':bool(gauss and boundary),'all_preparation_probabilities_retained':bool(abs(sum(r['probability'] for r in records)-1)<1e-12),'all_final_energy_number_outcomes_retained':bool(all(abs(sum(x['probability'] for x in r['complete_number_energy_eigenspace_records'])-r['probability'])<1e-12 for r in records))},'parameters':{'graph':'same triangle CAR halfhop source','mass':.23,'modulus':3,'source_dimension':source.n,'all_source_numbers':list(range(7)),'observed_input_output_numbers':[2,3,4],'source_reference':'actual full N3 sharedrecord ground multiplet normalizedprojector','source_reference_energy':float(e[0]),'source_reference_degeneracy':int(sum(selected)),'boundary_definition':'integer B=3-Nsource attached at port0; source G=n+divF+B delta0=(1,1,1,0,0,0) mod3; reservoir G_R=nreservoir-B exactinteger','reservoir_preparation':'fair initialoccupation0/1; initial label retained','reservoir_energy':reservoir_energy,'boundary_energy':0,'work_scope':'declared source H plus .23*n_reservoir; boundary energy set to zero as an explicit apparatus choice; no complete physical boundary/preparation cost','clock':'same Hnative hopparameter for passiveevolution; reservoir swap is supplied complete apparatus pulse, duration/physicalclock not derived'},'records':records,'completeness_error':completeness,'final_source_number_probabilities':{str(q):float(np.trace(out[np.ix_(np.flatnonzero(n==q),np.flatnonzero(n==q))])) for q in range(7)},'mean_apparatus_work':sum(r['unnormalized_apparatus_work'] for r in records),'scope':'Concrete source creation/removal with explicit boundary/reference dressing and complete reservoir outcomes, on one unchanged internal sharedrecord H. No unrelated-sector ground subtraction used as measured chargegap. No physicalspecies/vacuum/coupling/battery implementation or mass claim.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}

if __name__=='__main__':
 r=run();Path(__file__).with_name('CHARGED_RESERVOIR.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
